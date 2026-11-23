#!/usr/bin/env python
"""An OBJECTIVE keep-rule for programs: no hand labelling (D199).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

User, 2026-11-22: "we don't want to have to hand label all of the programs, if there is a program
that is not biologically relevant/confounds the ranking, we need a way to filter that objectively."

Every criterion below is computed from data already on disk. None requires reading a gene list,
naming a program, or any judgement call.

THRESHOLDS ARE NOT INVENTED HERE. D37 already fixed an objective keep-rule for Step 8 axes -
"above a shuffled-data null, stable under resampling, and not technical (plate above chance < 0.2;
RNA quality and purity rho < 0.5)" - so this reuses D37's 0.5 for covariate associations and
D37's "above chance" framing for plate, rather than choosing new numbers.

FOUR TESTS, in two families.

(A) ARTEFACT - the program tracks a technical covariate:
  A1 plate      r2_plate >= that cohort's OWN plate chance + PLATE_MARGIN. The chance level ranges
                0.000-0.378 across cohorts (small groups overlap plates structurally, D37), so an
                absolute threshold over-flags high-chance cohorts and misses low-chance ones.
  A2 RNA quality  |r_length| >= COVAR_MAX. The quality index is the length-based degradation proxy
                from 57/61; the pipeline corrects FOR it, so a program still tracking it is residual
                technical structure.

(B) CONFOUND - the program is real biology but cannot help a RANKING:
  B1 proliferation  |r_proliferation| >= COVAR_MAX. This is the axis D184 removes and the one D196
                found ranks approved cytotoxics well in COAD/PRAD/LUSC - good for recovering
                chemotherapy, useless for discovery.
  B2 generic activity  |Spearman(program opposition, ||drug activity||)| >= COVAR_MAX, within
                (cohort, window). If a program orders drugs the same way a bare "did this drug do
                anything" magnitude does, it carries no mechanism-specific information and will
                drag a ranking toward promiscuous compounds. THIS is the test D37 has no analogue
                for, because D37 screened axes of tumour biology, not of drug response.

A program is KEPT only if it fails to trigger all four. The output keeps every column so a
different rule can be applied later without recomputing.

Output: results/lincs/programs/program_filter.tsv
"""
import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
OUT = RES / "programs"
ACT = RES / "regulon_contributions_v2_tf/regulator_contributions_by_window_all.tsv.gz"
MODS = ROOT / "results/step6/aracne_project_v2_tf/diagnostics/regulator_modules_v2_tf.tsv"

COVAR_MAX = float(os.environ.get("COVAR_MAX", "0.5"))        # D37's rho threshold
PLATE_MARGIN = float(os.environ.get("PLATE_MARGIN", "0.10"))  # excess over the cohort's own chance


def main() -> int:
    P = pd.read_csv(OUT / "programs.tsv", sep="\t")
    D = pd.read_csv(OUT / "program_drug.tsv.gz", sep="\t")
    # only pull columns programs.tsv does not already carry - it gained r2_plate / plate_chance
    # in the D199 fix, and re-merging them produces r2_plate_x / _y and a silent AttributeError.
    want = [c for c in ("r_length", "r_proliferation", "r2_plate") if c not in P.columns]
    M = pd.read_csv(MODS, sep="\t", usecols=["group_id", "module"] + want)
    P = P.merge(M, left_on=["cohort", "program"], right_on=["group_id", "module"], how="left")

    # ||drug activity|| per (cohort, drug, window): the bare "did this drug do anything" scalar
    print("computing global drug activity magnitude ...", flush=True)
    acc = []
    for ch in pd.read_csv(ACT, sep="\t", usecols=["cohort", "pert_id", "window", "drug_activity"],
                          chunksize=3_000_000, low_memory=False):
        ch["sq"] = ch.drug_activity ** 2
        acc.append(ch.groupby(["cohort", "pert_id", "window"], sort=False)
                   .agg(ss=("sq", "sum"), n=("sq", "size")).reset_index())
    G = (pd.concat(acc, ignore_index=True)
         .groupby(["cohort", "pert_id", "window"], as_index=False)[["ss", "n"]].sum())
    G["act_norm"] = np.sqrt(G.ss / G.n)                      # RMS activity, size-independent
    del acc

    D = D.merge(G[["cohort", "pert_id", "window", "act_norm"]], on=["cohort", "pert_id", "window"],
                how="left")
    rho = (D.dropna(subset=["opposition", "act_norm"])
           .groupby(["cohort", "program", "window"])
           .apply(lambda t: t.opposition.corr(t.act_norm, method="spearman"))
           .rename("rho_generic").reset_index())
    # a program is judged by its WORST window: one window behaving like a magnitude score is enough
    # a program is judged by its WORST window; NaN where a window has no variance to correlate
    def _worst(s):
        s = s.dropna()
        return np.nan if s.empty else s.loc[s.abs().idxmax()]
    rho = rho.groupby(["cohort", "program"]).rho_generic.apply(_worst).reset_index()
    P = P.merge(rho, on=["cohort", "program"], how="left")

    P["fail_plate"] = P.r2_plate >= P.plate_chance + PLATE_MARGIN
    P["fail_quality"] = P.r_length.abs() >= COVAR_MAX
    P["fail_proliferation"] = P.r_proliferation.abs() >= COVAR_MAX
    P["fail_generic"] = P.rho_generic.abs() >= COVAR_MAX
    # D199 AMENDED, user 2026-11-22: "we don't want to remove all generic activity and
    # proliferation as of yet". These are REPORTED as annotations, never used to drop a program.
    # The artefact flags (plate, RNA quality) stay as flags because a plate-tracking module is not
    # biology under any reading; the two CONFOUND flags become inputs to a specificity score
    # instead of an exclusion.
    fails = ["fail_plate", "fail_quality", "fail_proliferation", "fail_generic"]
    P["n_flag"] = P[fails].sum(axis=1)
    P["artefact"] = P.fail_plate | P.fail_quality

    P.drop(columns=["group_id", "module"], errors="ignore").to_csv(
        OUT / "program_filter.tsv", sep="\t", index=False, float_format="%.4g")

    print(f"\n{len(P)} programs over {P.cohort.nunique()} cohorts")
    for f in fails:
        print(f"  {f:22s} {int(P[f].sum()):4d}")
    print(f"  {'artefact (plate|quality)':22s} {int(P.artefact.sum()):4d}")
    print(f"\nflags per program: {P.n_flag.value_counts().sort_index().to_dict()}")
    print(f"\nrho_generic: median |rho| {P.rho_generic.abs().median():.3f}, "
          f">= {COVAR_MAX} in {int(P.fail_generic.sum())}")
    (OUT / "program_filter_run.json").write_text(json.dumps(dict(
        script="79aj_program_filter.py", covar_max=COVAR_MAX, plate_margin=PLATE_MARGIN,
        thresholds_from="D37 axis keep-rule", n_programs=int(len(P)), n_artefact=int(P.artefact.sum()),
        fails={f: int(P[f].sum()) for f in fails}), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
