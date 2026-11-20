#!/usr/bin/env python
"""79c2: condition-level drug profiles - keep the dose and time structure (D153, proposed).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
SPEC: docs/SPEC_condition_level_profiles.md

79c averages every active signature of a (drug, line) over doses AND times into one profile, and the
index carries only (pert_id, ModelID, n_sig). This writes the same thing keyed on
(pert_id, ModelID, pert_dose, pert_time) instead, ALONGSIDE the existing files - nothing in
production is replaced or read differently.

Two deliberate differences from 79c, both from the spec:

  * `g`, the generic drug response projected out of the specific arm, is computed per
    (ModelID, pert_time) rather than per ModelID. 79c's own header records that this projection is
    load-bearing; and the generic response at 6 h is not the one at 24 h, so pooling them would
    subtract a predominantly-24 h vector from 6 h profiles and manufacture a new artefact in exactly
    the arm D152 is about. Strata with fewer than MIN_DRUGS_FOR_G drugs are left unprojected and
    flagged `generic_removed = False`.

  * TAS_MIN is a parameter. The default reproduces production (0.4). It is NOT lowered here - the
    quality floor is a separate pre-registration (D152). This script only makes that question
    testable by keeping the conditions.

Outputs (data/processed/lincs/):
  drug_condition_profiles_raw.f32       conditions x genes, float32, profile-major
  drug_condition_profiles_specific.f32  same, generic response removed per (line, time)
  drug_condition_profiles_index.tsv     pert_id, ModelID, pert_dose, pert_time, n_sig, median_tas,
                                        generic_removed, n_drugs_in_stratum
  drug_condition_profiles_run.json      manifest
"""
import json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed/lincs"
TAS_MIN = float(os.environ.get("TAS_MIN", "0.4"))
MIN_DRUGS_FOR_G = int(os.environ.get("MIN_DRUGS_FOR_G", "30"))
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(v, "8")


def main():
    t0 = time.time()
    X = np.load(PROC / "level5_bing_hiq.npy", mmap_mode="r")
    sig = pd.read_csv(PROC / "level5_bing_hiq_sigs.tsv.gz", sep="\t")
    genes = pd.read_csv(PROC / "level5_bing_genes.tsv", sep="\t", dtype={"gene_id": str})

    # gene mapping: identical to 79c, so the gene axis matches the production artifact exactly
    idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t", dtype={"entrez_id": str})
    idm = idm[idm.primary_for_hgnc == True]
    by_entrez = idm.dropna(subset=["entrez_id"]).drop_duplicates("entrez_id").set_index("entrez_id").gene_id
    v23 = genes.gene_id.map(by_entrez)
    v23 = v23.fillna(genes.ensembl_id.where(genes.ensembl_id.isin(idm.gene_id)))
    keep = v23.notna() & ~v23.duplicated()
    gpos = np.where(keep.values)[0]
    print(f"genes kept (unique v23): {int(keep.sum())}")

    act = sig[sig.tas >= TAS_MIN]
    # pandas groupby DROPS rows whose key contains NaN. 4 active signatures have a NaN pert_dose,
    # which silently lost 3 drug x line pairs (one entirely) and broke the reproduction gate.
    # Use an explicit sentinel so "dose not recorded" is its own condition rather than vanishing.
    n_nan_dose = int(act.pert_dose.isna().sum())
    n_nan_time = int(act.pert_time.isna().sum())
    act = act.assign(pert_dose=act.pert_dose.fillna(-1.0), pert_time=act.pert_time.fillna(-1.0))
    grp = act.groupby(["pert_id", "ModelID", "pert_dose", "pert_time"]).indices
    assert sum(len(v) for v in grp.values()) == len(act), "signatures lost in the condition grouping"
    idx = pd.DataFrame(list(grp.keys()), columns=["pert_id", "ModelID", "pert_dose", "pert_time"])
    idx["n_sig"] = [len(v) for v in grp.values()]
    rows = act.index.values
    tasv = act.tas.values
    idx["median_tas"] = [float(np.median(tasv[v])) for v in grp.values()]
    idx["dose_recorded"] = idx.pert_dose >= 0
    idx["time_recorded"] = idx.pert_time >= 0
    print(f"conditions: {len(idx):,} (from {len(act):,} active signatures, TAS >= {TAS_MIN}); "
          f"{idx.groupby(['pert_id','ModelID']).ngroups:,} drug x line pairs")

    P = np.empty((len(idx), len(gpos)), dtype=np.float32)
    for i, v in enumerate(grp.values()):
        P[i] = np.asarray(X[np.sort(rows[v])][:, gpos], dtype=np.float64).mean(axis=0)
    assert not np.isnan(P).any()
    P.tofile(PROC / "drug_condition_profiles_raw.f32")
    print(f"raw written, {time.time()-t0:.0f}s")

    # generic response per (ModelID, pert_time)
    S = P.astype(np.float64)
    idx["generic_removed"] = False
    idx["n_drugs_in_stratum"] = 0
    for (model, ptime), ix in idx.groupby(["ModelID", "pert_time"]).indices.items():
        n_drugs = idx.pert_id.iloc[ix].nunique()
        idx.loc[idx.index[ix], "n_drugs_in_stratum"] = n_drugs
        if n_drugs < MIN_DRUGS_FOR_G:
            continue
        g = S[ix].mean(axis=0)
        nrm = np.linalg.norm(g)
        if nrm == 0:
            continue
        g /= nrm
        S[ix] = S[ix] - np.outer(S[ix] @ g, g)
        idx.loc[idx.index[ix], "generic_removed"] = True
    S.astype(np.float32).tofile(PROC / "drug_condition_profiles_specific.f32")
    idx.to_csv(PROC / "drug_condition_profiles_index.tsv", sep="\t", index=False)

    strata = idx.groupby(["ModelID", "pert_time"]).ngroups
    man = dict(script="79c2_condition_profiles.py", spec="docs/SPEC_condition_level_profiles.md",
               tas_min=TAS_MIN, min_drugs_for_g=MIN_DRUGS_FOR_G,
               n_conditions=int(len(idx)), n_signatures=int(len(act)),
               n_drug_line_pairs=int(idx.groupby(["pert_id", "ModelID"]).ngroups),
               n_genes=int(len(gpos)), n_strata=int(strata),
               n_conditions_generic_removed=int(idx.generic_removed.sum()),
               pct_generic_removed=round(100 * float(idx.generic_removed.mean()), 2),
               doses_per_pair_median=float(idx.groupby(["pert_id", "ModelID"]).pert_dose.nunique().median()),
               times_per_pair_median=float(idx.groupby(["pert_id", "ModelID"]).pert_time.nunique().median()),
               n_signatures_nan_dose=n_nan_dose, n_signatures_nan_time=n_nan_time,
               sentinel_for_unrecorded=-1.0,
               runtime_s=round(time.time() - t0, 1))
    (PROC / "drug_condition_profiles_run.json").write_text(json.dumps(man, indent=1))
    print(f"specific written; {strata:,} (line, time) strata, "
          f"{man['pct_generic_removed']}% of conditions had the generic response removed")
    print(f"done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    sys.exit(main())
