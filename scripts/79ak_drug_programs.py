#!/usr/bin/env python
"""Programs in the DRUG profiles, per cohort and window, and how they relate to the disease ones (D200).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

User, 2026-11-22: "let's focus on identifying the programs within the disease cohort and then
within the drug profiles."

D198 built programs from the DISEASE side - regulators co-varying across PATIENTS. This builds the
other half: regulators co-varying across DRUGS, within one cohort and one timepoint window. They
answer different questions and there is no reason they should coincide:

    disease program   regulators the tumour runs together        (what is dysregulated)
    drug program      regulators a perturbation moves together   (what is pharmacologically reachable)

A drug program with no disease counterpart is biology no tumour in this cohort varies along; a
disease program with no drug counterpart is dysregulation nothing in LINCS can move. The overlap
is the part where a ranking can actually operate.

SAME CLUSTERING RULE AS 58, deliberately, so the two sides are comparable rather than merely
adjacent: average-linkage on 1 - |r|, cut at |r| >= MOD_R (0.7), modules of >= MIN_SIZE.

WINDOW-SEPARATED, never pooled (D183: ~80% of a window's reliable variance is absent from the
other; D198: aggregating to programs does NOT stabilise that, Spearman +0.337).

Outputs:
  results/lincs/programs/drug_programs.tsv        one row per (cohort, window, program)
  results/lincs/programs/program_overlap.tsv      disease x drug program Jaccard, per cohort
"""
import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
OUT = RES / "programs"
ACT = RES / "regulon_contributions_v2_tf/regulator_contributions_by_window_all.tsv.gz"

MOD_R = float(os.environ.get("MOD_R", "0.7"))      # as 58
MIN_SIZE = int(os.environ.get("MIN_SIZE", "5"))
MIN_DRUGS = int(os.environ.get("MIN_DRUGS", "40"))


def modules(A: pd.DataFrame):
    """A: drugs x regulators. Returns {module_id: [regulators]}, 58's rule."""
    R = np.corrcoef(A.values.T)
    R = np.nan_to_num(R, nan=0.0)
    np.fill_diagonal(R, 1.0)
    d = 1.0 - np.abs(R)
    np.fill_diagonal(d, 0.0)
    d = (d + d.T) / 2.0
    lab = fcluster(linkage(squareform(d, checks=False), method="average"), t=1 - MOD_R,
                   criterion="distance")
    out = {}
    for m in np.unique(lab):
        mem = A.columns[lab == m].tolist()
        if len(mem) >= MIN_SIZE:
            out[int(m)] = mem
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    cols = ["cohort", "pert_id", "window", "regulator", "drug_activity", "disease_nes"]
    chunks = [ch for ch in pd.read_csv(ACT, sep="\t", usecols=lambda c: c in cols,
                                       chunksize=3_000_000, low_memory=False)]
    A = pd.concat(chunks, ignore_index=True)
    del chunks

    dis = pd.read_csv(OUT / "programs.tsv", sep="\t")
    dis["members"] = dis.members.str.split(",")
    dmap = {c: {int(r.program): set(r.members) for r in g.itertuples()}
            for c, g in dis.groupby("cohort")}

    rows, ov = [], []
    for (coh, win), g in A.groupby(["cohort", "window"]):
        M = g.pivot_table(index="pert_id", columns="regulator", values="drug_activity")
        M = M.dropna(axis=1, how="any")
        if M.shape[0] < MIN_DRUGS or M.shape[1] < 50:
            continue
        nes = g.drop_duplicates("regulator").set_index("regulator").disease_nes
        mods = modules(M)
        for mid, mem in mods.items():
            sc = M[mem].mean(axis=1)
            # hub = member most correlated with the module score
            rr = {m: abs(np.corrcoef(M[m].values, sc.values)[0, 1]) for m in mem}
            hub = max(rr, key=rr.get)
            rows.append(dict(cohort=coh, window=win, program=mid, n_regulators=len(mem),
                             hub=hub, hub_r=round(rr[hub], 3),
                             disease_score=round(float(nes.reindex(mem).mean()), 4),
                             members=",".join(mem)))
            # best-matching DISEASE program
            best, bj = None, 0.0
            for did, dm in dmap.get(coh, {}).items():
                j = len(dm & set(mem)) / len(dm | set(mem))
                if j > bj:
                    best, bj = did, j
            ov.append(dict(cohort=coh, window=win, drug_program=mid, n_regulators=len(mem),
                           best_disease_program=best, jaccard=round(bj, 4)))
        print(f"{coh:22s} {win:5s} drugs {M.shape[0]:5d} regs {M.shape[1]:4d} -> {len(mods):3d} programs",
              flush=True)

    P = pd.DataFrame(rows)
    P.to_csv(OUT / "drug_programs.tsv", sep="\t", index=False)
    O = pd.DataFrame(ov)
    O.to_csv(OUT / "program_overlap.tsv", sep="\t", index=False)

    print(f"\n{len(P)} drug programs over {P.cohort.nunique()} cohorts x {P.window.nunique()} windows")
    print(f"  per (cohort, window): median {int(P.groupby(['cohort','window']).size().median())}")
    print(f"  size: median {int(P.n_regulators.median())}, max {int(P.n_regulators.max())}")
    print(f"\nbest Jaccard to a DISEASE program: median {O.jaccard.median():.3f}, "
          f">=0.25 in {int((O.jaccard>=0.25).sum())}/{len(O)}, >=0.50 in {int((O.jaccard>=0.50).sum())}")
    (OUT / "drug_programs_run.json").write_text(json.dumps(dict(
        script="79ak_drug_programs.py", mod_r=MOD_R, min_size=MIN_SIZE, min_drugs=MIN_DRUGS,
        n_drug_programs=int(len(P)), median_jaccard_to_disease=float(O.jaccard.median())), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
