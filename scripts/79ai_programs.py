#!/usr/bin/env python
"""Programs within a cohort: disease-defined regulator modules, projected onto the drug side (D198).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

User, 2026-11-22: "let's attempt to identify the programs within a cohort."

WHERE THE PROGRAMS COME FROM, AND WHY NOT FROM THE DRUGS. A program here is a module of regulators
that co-vary ACROSS PATIENTS within one cohort - `58_regulator_modules.R`, average-linkage on
1 - |r| of patient NES, cut at |r| >= 0.7. Defining them on the disease side matters:

  * it is tumour biology, not a property of whatever compounds LINCS happens to contain, so
    projecting drugs onto them is not circular;
  * D196 showed that factorising the DRUG x regulator matrix and scoring by cosine in that basis
    is an orthonormal rotation - it cannot change a cosine, and the measurement confirmed it
    (AUC converges to full-space at k >= 40). Disease-defined modules are NOT a rotation of the
    drug matrix, so they are not covered by that null.

WHY WITHIN COHORT ONLY. D197 measured what a regulator-level quantity costs when carried across
cohorts: with IDENTICAL cell lines but different networks (BRCA-LumA vs BRCA-Basal) the same
drug's regulator activity correlates only +0.667, and cross-tissue +0.330, because regulon CONTENT
overlap is 0.063-0.093. Programs are therefore built and scored inside one cohort and never
transferred.

WHAT IS COMPUTED, per (cohort, program):
    disease_score    mean disease NES over the module's regulators    (how the tumour runs it)
    n_regulators     module size
    hallmark         best-matching MSigDB hallmark set over the members' TARGETS, with its q
    is_proliferation |r_proliferation| >= 0.5   from 58
    is_technical     r2_plate >= that cohort's OWN plate chance level + PLATE_MARGIN. D37 warns
                     that plate R^2 sits above chance for structural reasons (small groups overlap
                     plates), and the chance level ranges 0.000-0.378 ACROSS cohorts, so an
                     absolute threshold over- and under-flags at the same time. Chance comes from
                     `quality_correction/checks.tsv` (`plate_chance`). NOTE: the pipeline never
                     regresses plate out - D8 fixed the signature model with no plate term, and
                     61_quality_correct.R corrects purity + an RNA-quality index that was built
                     because it explained the plate pattern. Plate is handled by PROXY, so a
                     residual association survives into patient NES and into these modules.
and per (cohort, program, drug, WINDOW), window-separated throughout (D183/D191):
    drug_score       mean drug regulator-activity over the module's regulators
    opposition       -disease_score * drug_score / |disease_score|, the program-level analogue of
                     rev's per-regulator term: positive = the drug pushes the program against the
                     direction the tumour runs it.

Outputs:
  results/lincs/programs/programs.tsv              one row per (cohort, program)
  results/lincs/programs/program_drug.tsv.gz       one row per (cohort, program, drug, window)
  results/lincs/programs/run.json
"""
import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
OUT = RES / "programs"
MODS = ROOT / "results/step6/aracne_project_v2_tf/diagnostics/regulator_modules_v2_tf.tsv"
EDGES = ROOT / "data/processed/lincs/trace_aracne/aracne_regulon_edges.tsv.gz"
GMT = ROOT / "data/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt"
ACT = RES / "regulon_contributions_v2_tf/regulator_contributions_by_window_all.tsv.gz"
MIN_SIZE = int(os.environ.get("MIN_SIZE", "5"))
PLATE_MARGIN = float(os.environ.get("PLATE_MARGIN", "0.10"))   # excess over the cohort's chance level


def hallmarks():
    out = {}
    for line in GMT.read_text().splitlines():
        f = line.split("\t")
        out[f[0].replace("HALLMARK_", "")] = set(f[2:])
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    M = pd.read_csv(MODS, sep="\t")
    if "members_all" not in M.columns:
        print("ABORT: regulator_modules lacks members_all - re-run 58_regulator_modules.R.", file=sys.stderr)
        return 1
    M = M[M["size"] >= MIN_SIZE].copy()
    M["members"] = M.members_all.str.split(",")

    import glob
    ck = pd.concat([pd.read_csv(f, sep="\t") for f in
                    glob.glob(str(ROOT / "results/step6/aracne_project_v2_tf/quality_correction*/checks.tsv"))])
    plate_chance = ck.groupby("group_id").plate_chance.median().to_dict()

    idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t", usecols=["gene_id", "symbol"])
    ens2sym = dict(idm.dropna().values)
    E = pd.read_csv(EDGES, sep="\t", usecols=["network", "regulator", "gene_id", "in_lincs_gene_space"])
    E = E[E.in_lincs_gene_space]
    HM = hallmarks()
    universe = {ens2sym.get(g) for g in E.gene_id.unique()} - {None}
    NU = len(universe)
    HM = {k: (v & universe) for k, v in HM.items()}

    # network per cohort, as the contribution table already resolved it
    act_cols = ["cohort", "pert_id", "cmap_name", "window", "regulator", "disease_nes", "drug_activity"]

    chunks = []
    for ch in pd.read_csv(ACT, sep="\t", usecols=lambda c: c in act_cols, chunksize=3_000_000,
                          low_memory=False):
        chunks.append(ch)
    A = pd.concat(chunks, ignore_index=True)
    del chunks

    prog_rows, pd_rows = [], []
    for coh, md in M.groupby("group_id"):
        sub = A[A.cohort == coh]
        if sub.empty:
            continue
        nes = sub.drop_duplicates("regulator").set_index("regulator").disease_nes
        piv = sub.pivot_table(index=["pert_id", "cmap_name", "window"], columns="regulator",
                              values="drug_activity")
        for r in md.itertuples():
            mem = [m for m in r.members if m in piv.columns]
            if len(mem) < MIN_SIZE:
                continue
            ds = float(nes.reindex(mem).mean())
            # hallmark enrichment over the union of the members' targets
            tg = {ens2sym.get(g) for g in E[E.regulator.isin(mem)].gene_id.unique()} - {None}
            tg &= universe
            best, bestp = "", 1.0
            for h, hs in HM.items():
                k = len(tg & hs)
                if k < 5:
                    continue
                p = hypergeom.sf(k - 1, NU, len(hs), len(tg))
                if p < bestp:
                    best, bestp = h, p
            prog_rows.append(dict(cohort=coh, program=int(r.module), n_regulators=len(mem),
                                  hub=r.hub, disease_score=round(ds, 4),
                                  hallmark=best, hallmark_q=float(min(1.0, bestp * len(HM))),
                                  is_proliferation=bool(abs(r.r_proliferation) >= 0.5),
                                  plate_chance=round(float(plate_chance.get(coh, float("nan"))), 4),
                                  r2_plate=round(float(r.r2_plate), 4),
                                  is_technical=bool(r.r2_plate >= plate_chance.get(coh, 0.0) + PLATE_MARGIN),
                                  members=",".join(mem)))
            v = piv[mem].mean(axis=1)
            opp = -ds * v / abs(ds) if ds != 0 else v * np.nan
            d = v.rename("drug_score").reset_index()
            d["cohort"], d["program"] = coh, int(r.module)
            d["opposition"] = opp.values
            pd_rows.append(d)
        print(f"{coh:22s} {len(md)} modules -> {sum(1 for x in prog_rows if x['cohort']==coh)} programs",
              flush=True)

    P = pd.DataFrame(prog_rows)
    P.to_csv(OUT / "programs.tsv", sep="\t", index=False)
    D = pd.concat(pd_rows, ignore_index=True)
    D.to_csv(OUT / "program_drug.tsv.gz", sep="\t", index=False, float_format="%.6g")

    nb = int((~P.is_proliferation & ~P.is_technical).sum())
    print(f"\n{len(P)} programs over {P.cohort.nunique()} cohorts "
          f"({nb} biological, {int(P.is_proliferation.sum())} proliferation, "
          f"{int(P.is_technical.sum())} plate-tracking)")
    print(f"{len(D):,} (cohort, program, drug, window) rows")
    (OUT / "run.json").write_text(json.dumps(dict(
        script="79ai_programs.py", min_size=MIN_SIZE, n_programs=len(P),
        n_cohorts=int(P.cohort.nunique()), n_rows=int(len(D)),
        source_modules=str(MODS.relative_to(ROOT)),
        note="programs are disease-defined (patient co-variation) and never transferred across "
             "cohorts; see D197 for why"), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
