#!/usr/bin/env python
"""D135 stage S0, part B -- FEASIBILITY AND QC for docs/SPEC_drug_genetic_concordance.md. Counts only.

NO CONCORDANCE IS COMPUTED HERE. No drug profile is correlated with any perturbation profile, no positive
control is run, no AUROC, no null, no statistic that compares a drug with a knockout. Every quantity is an
annotation count, a presence count, or a PER-PROFILE statistic. Definitions were fixed in docs/DECISIONS.md
row D135-S0, committed and pushed before this ran.

Inputs
  results/lincs/d135/perturbation_sigs.tsv.gz   per-signature strength + target-mRNA QC (from s0a)
  results/lincs/target_role/d121/units.tsv      target annotations (rows with non-empty tid)
  results/lincs/cell_coverage.tsv               cell_iname -> ModelID
  data/processed/lincs/drug_line_profiles_{index,genes}.tsv, drug_line_profiles_bing.f32
  data/processed/lincs/level5_bing_hiq_sigs.tsv.gz        drug TAS and timepoints
  data/raw/lincs2020/siginfo_beta.txt                     coverage counts at every filter
  data/raw/depmap_24q4/CRISPRInferredCommonEssentials.csv, OmicsExpressionProteinCodingGenesTPMLogp1.csv,
  data/raw/depmap_24q4/OmicsSomaticMutationsMatrix{Damaging,Hotspot}.csv

Outputs (results/lincs/d135/)
  perturbation_coverage.tsv   per (arm, gene): signatures and lines at each filter
  perturbation_qc.tsv         per (arm, gene, line): strength and target-mRNA QC of that perturbation cell
  paired_cells.tsv            per (arm, drug, gene, line): the real testable scope, with every stratum
  attrition.tsv               counts at every step from 377 genes to clean-positive-capable cells
  strata_summary.tsv          how the paired set distributes over the spec Section 6 strata
  d135_s0_run.json            headline numbers, provenance, every definition used
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "8"

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
RAW = ROOT / "data/raw/lincs2020"
DEP = ROOT / "data/raw/depmap_24q4"
PROC = ROOT / "data/processed/lincs"
RES = ROOT / "results/lincs"
OUT = RES / "d135"
OUT.mkdir(parents=True, exist_ok=True)

N_GENES = 10112
CHUNK = 512
KD_ARMS = ["trt_xpr", "trt_sh.cgs", "trt_sh"]
NZ2_FLOORS = [0, 50, 100, 200, 400, 800]
TAS_FLOORS = [0.1, 0.2, 0.3, 0.4]
EXPR_CUT = 1.0          # log2(TPM+1); descriptive, chosen for reporting only
QC_CUTS = {"z_le_-1": -1.0, "z_le_-2": -2.0, "z_lt_0": 0.0}
ILLUS_PERT_NZ2 = 100    # illustrative only, adopted by nothing
ILLUS_DRUG_NZ2 = 100


def qstat(x, name=""):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {"n": 0, "note": f"{name}: no finite values"}
    q = np.percentile(x, [0, 5, 25, 50, 75, 95, 100])
    return {"n": int(x.size), "mean": round(float(x.mean()), 4), "min": round(float(q[0]), 4),
            "p05": round(float(q[1]), 4), "p25": round(float(q[2]), 4), "median": round(float(q[3]), 4),
            "p75": round(float(q[4]), 4), "p95": round(float(q[5]), 4), "max": round(float(q[6]), 4)}


def main():
    t0 = time.time()
    prov = {}

    # ------------------------------------------------------------------ 1. target scope (D135-S0(1))
    u_all = pd.read_csv(RES / "target_role/d121/units.tsv", sep="\t")
    u = u_all[u_all.tid.notna() & (u_all.tid.astype(str).str.strip() != "")].copy()
    rows = []
    for _, r in u.iterrows():
        for g in str(r.genes).split("|"):
            if g and g != "nan":
                rows.append(dict(cmap_name=r.cmap_name, tid=r.tid, gene=g, orient=r.orient,
                                 actions=r.actions, ttype=r.ttype))
    dg_ann = pd.DataFrame(rows).drop_duplicates()
    genes = sorted(dg_ann.gene.unique())
    drug_genes = dg_ann.groupby("cmap_name").gene.apply(lambda s: sorted(set(s)))
    drug_units = dg_ann.groupby("cmap_name").tid.nunique()
    gene_drugs = dg_ann.groupby("gene").cmap_name.apply(lambda s: sorted(set(s)))
    prov["scope"] = dict(unit_rows_total=int(len(u_all)), unit_rows_used=int(len(u)),
                         compounds=int(u.cmap_name.nunique()), units=int(u.tid.nunique()),
                         target_genes=len(genes),
                         drug_gene_pairs=int(len(dg_ann[["cmap_name", "gene"]].drop_duplicates())),
                         drugs_with_1_gene=int((drug_genes.apply(len) == 1).sum()),
                         drugs_with_ge2_genes=int((drug_genes.apply(len) >= 2).sum()),
                         genes_per_drug=qstat(drug_genes.apply(len), "genes_per_drug"))

    # ------------------------------------------------------------------ 2. drug universe (D135-S0(2))
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    gsp = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t", dtype={"gene_id": str})
    assert len(gsp) == N_GENES
    dsig = pd.read_csv(PROC / "level5_bing_hiq_sigs.tsv.gz", sep="\t")
    p2c = dsig.drop_duplicates("pert_id").set_index("pert_id").cmap_name
    idx["cmap_name"] = idx.pert_id.map(p2c)
    assert idx.cmap_name.notna().all()
    drug_lines = sorted(set(idx.ModelID))
    raw = np.memmap(PROC / "drug_line_profiles_bing.f32", dtype=np.float32, mode="r",
                    shape=(len(idx), N_GENES))
    nz2 = np.empty(len(idx), np.int32); nz4 = np.empty(len(idx), np.int32)
    mab = np.empty(len(idx), np.float32)
    for i in range(0, len(idx), CHUNK):
        j = min(i + CHUNK, len(idx))
        a = np.abs(np.asarray(raw[i:j], dtype=np.float32))
        nz2[i:j] = (a >= 2).sum(1); nz4[i:j] = (a >= 4).sum(1); mab[i:j] = a.mean(1)
    idx["nz2"], idx["nz4"], idx["mean_abs"] = nz2, nz4, mab
    act = dsig[(dsig.is_hiq == 1) & (dsig.tas >= 0.4)]
    tas = act.groupby(["pert_id", "ModelID"]).tas.max().rename("tas_max")
    idx = idx.merge(tas, left_on=["pert_id", "ModelID"], right_index=True, how="left")
    dcl = idx.groupby(["cmap_name", "ModelID"]).agg(
        drug_nz2=("nz2", "max"), drug_nz4=("nz4", "max"), drug_mean_abs=("mean_abs", "max"),
        drug_tas_max=("tas_max", "max"), n_pert_ids=("pert_id", "nunique")).reset_index()
    tp = act.groupby(["cmap_name", "ModelID"]).pert_time.agg(
        drug_n_timepoints="nunique", drug_timepoints=lambda s: "|".join(sorted({str(x) for x in s})))
    dcl = dcl.merge(tp, left_on=["cmap_name", "ModelID"], right_index=True, how="left")
    present = {c: set(g) for c, g in dcl.groupby("cmap_name").ModelID}
    ann_present = {c: v for c, v in present.items() if c in set(dg_ann.cmap_name)}
    prov["drug_universe"] = dict(
        profile_rows=int(len(idx)), compounds=int(idx.cmap_name.nunique()), lines=len(drug_lines),
        genes=N_GENES, annotated_compounds_profiled=len(ann_present),
        annotated_compound_line_profiles=int(dcl.cmap_name.isin(ann_present).sum()),
        annotated_lines=int(dcl[dcl.cmap_name.isin(ann_present)].ModelID.nunique()),
        note="already carries the 79c floor is_hiq == 1 and per-signature tas >= 0.4 (D135-S0(2))")

    # ------------------------------------------------------------------ 3. coverage at every filter (S0(3))
    cov = pd.read_csv(RES / "cell_coverage.tsv", sep="\t").dropna(subset=["ModelID"])
    cov = cov[cov.OncotreePrimaryDisease != "Non-Cancerous"]
    lmap = cov.drop_duplicates("cell_iname").set_index("cell_iname").ModelID
    si = pd.read_csv(RAW / "siginfo_beta.txt", sep="\t", low_memory=False,
                     usecols=["sig_id", "cmap_name", "pert_type", "cell_iname", "is_hiq", "tas", "pert_time"])
    si = si[si.pert_type.isin(KD_ARMS + ["trt_oe"]) & si.cmap_name.isin(genes)].copy()
    si["ModelID"] = si.cell_iname.map(lmap)
    cov_rows, cov_gene = [], []
    for arm in KD_ARMS + ["trt_oe"]:
        a = si[si.pert_type == arm]
        a_m = a[a.ModelID.notna()]
        a_d = a_m[a_m.ModelID.isin(drug_lines)]
        for lab, d, lcol in [("all_signatures", a, "cell_iname"),
                             ("on_mapped_cancer_line", a_m, "ModelID"),
                             ("on_mapped_cancer_line_hiq", a_m[a_m.is_hiq == 1], "ModelID"),
                             ("on_drug_matrix_line", a_d, "ModelID"),
                             ("on_drug_matrix_line_hiq", a_d[a_d.is_hiq == 1], "ModelID")]:
            cov_rows.append(dict(arm=arm, filter=lab, signatures=int(len(d)),
                                 genes=int(d.cmap_name.nunique()), lines=int(d[lcol].nunique())))
        for g, d in a_d.groupby("cmap_name"):
            cov_gene.append(dict(arm=arm, gene=g, signatures=int(len(d)), lines=int(d.ModelID.nunique()),
                                 signatures_hiq=int((d.is_hiq == 1).sum()),
                                 lines_hiq=int(d[d.is_hiq == 1].ModelID.nunique())))
    coverage = pd.DataFrame(cov_rows)
    pd.DataFrame(cov_gene).to_csv(OUT / "perturbation_coverage.tsv", sep="\t", index=False)
    prov["coverage_by_filter"] = coverage.to_dict("records")

    # spec Section 2 check -- which filter, if any, reproduces each published figure (D135-S0(11))
    spec = {"trt_xpr": dict(targets=208, signatures=837, lines=8),
            "trt_sh.cgs": dict(targets=209, signatures=1627, lines=10),
            "trt_sh": dict(targets=207, signatures=6118, lines=12)}
    s2 = {}
    for arm, claim in spec.items():
        hits = [r["filter"] for r in cov_rows if r["arm"] == arm and r["signatures"] == claim["signatures"]
                and r["genes"] == claim["targets"] and r["lines"] == claim["lines"]]
        s2[arm] = dict(spec_claim=claim, reproduced_by=hits,
                       measured={r["filter"]: {k: r[k] for k in ("signatures", "genes", "lines")}
                                 for r in cov_rows if r["arm"] == arm})
    prov["spec_section2_check"] = s2

    # ------------------------------------------------------------------ 4. perturbation cells (S0(7),(8))
    ps = pd.read_csv(OUT / "perturbation_sigs.tsv.gz", sep="\t")
    prov["perturbation_signatures_extracted"] = {a: int((ps.arm == a).sum()) for a in ps.arm.unique()}
    prov["perturbation_strength"] = {
        a: {k: qstat(ps.loc[ps.arm == a, k], k) for k in ["nz2", "nz4", "mean_abs", "tas"]}
        for a in ps.arm.unique()}
    prov["perturbation_pert_time"] = {a: ps.loc[ps.arm == a, "pert_time"].value_counts().to_dict()
                                      for a in ps.arm.unique()}

    # target-mRNA QC, per signature, per arm, split by feature space (D135-S0(8))
    qc = {}
    for a in ps.arm.unique():
        d = ps[(ps.arm == a) & ps.target_on_gene_space]
        up = a == "trt_oe"
        block = dict(signatures_total=int((ps.arm == a).sum()),
                     signatures_with_target_on_gene_space=int(len(d)),
                     signatures_target_absent=int((ps.arm == a).sum() - len(d)),
                     target_z=qstat(d.target_z, "target_z"))
        for lab, cut in QC_CUTS.items():
            lab2 = lab.replace("le_-", "ge_").replace("lt_0", "gt_0") if up else lab
            v = (d.target_z >= -cut) if up else (d.target_z <= cut)
            block[f"pass_{lab2}"] = int(v.sum())
            block[f"pass_rate_{lab2}"] = round(float(v.mean()), 4) if len(d) else None
        for fs in ["landmark", "best inferred"]:
            dd = d[d.target_feature_space == fs]
            sub = dict(signatures=int(len(dd)), target_z=qstat(dd.target_z, "z"),
                       genes=int(dd.cmap_name.nunique()))
            for lab, cut in QC_CUTS.items():
                lab2 = lab.replace("le_-", "ge_").replace("lt_0", "gt_0") if up else lab
                v = (dd.target_z >= -cut) if up else (dd.target_z <= cut)
                sub[f"pass_rate_{lab2}"] = round(float(v.mean()), 4) if len(dd) else None
            block[fs.replace(" ", "_")] = sub
        qc[a] = block
    prov["target_mrna_qc_per_signature"] = qc
    prov["target_mrna_qc_reading"] = (
        "Pre-stated in D135-S0(8), before the numbers were seen: a shRNA failure and a CRISPR failure do not "
        "mean the same thing. shRNA degrades the transcript, so reduced mRNA is the direct expected readout "
        "and a failure is evidence the knockdown did not work. CRISPR disrupts the GENE; a frameshift "
        "transcript escaping nonsense-mediated decay leaves mRNA near normal while the protein is gone, so a "
        "CRISPR perturbation can be fully effective with no mRNA movement. A CRISPR mRNA-QC failure is NOT a "
        "failed perturbation and is never scored as one. Only 48 of the 377 targets are L1000 landmark genes; "
        "for the rest the 'measured' mRNA is a best-inferred model output, which weakens the QC further.")

    # per (arm, gene, line) perturbation cell
    def best_z(s):
        s = s.dropna()
        return float(s.min()) if len(s) else np.nan

    cell = ps.groupby(["arm", "cmap_name", "ModelID"]).agg(
        n_sigs=("sig_id", "size"), pert_nz2=("nz2", "max"), pert_nz4=("nz4", "max"),
        pert_mean_abs=("mean_abs", "max"), pert_tas=("tas", "max"), n_hiq=("is_hiq", "sum"),
        target_z_min=("target_z", best_z), target_z_max=("target_z", "max"),
        target_on_gene_space=("target_on_gene_space", "max"),
        target_feature_space=("target_feature_space", "first"),
        pert_times=("pert_time", lambda s: "|".join(sorted({str(x) for x in s})))).reset_index()
    cell = cell.rename(columns={"cmap_name": "gene"})
    oe = cell.arm == "trt_oe"
    for lab, cut in QC_CUTS.items():
        # knockdown direction; left NaN on the overexpression arm, which has its own mirrored column
        cell[f"qc_{lab}"] = np.where(cell.target_on_gene_space & ~oe,
                                     cell.target_z_min <= cut if cut < 0 else cell.target_z_min < 0, np.nan)
        cell[f"qc_oe_{lab.replace('le_-', 'ge_').replace('lt_0', 'gt_0')}"] = np.where(
            cell.target_on_gene_space & oe,
            cell.target_z_max >= -cut if cut < 0 else cell.target_z_max > 0, np.nan)
    cell.to_csv(OUT / "perturbation_qc.tsv", sep="\t", index=False)
    prov["perturbation_cells"] = {a: dict(cells=int((cell.arm == a).sum()),
                                          genes=int(cell[cell.arm == a].gene.nunique()),
                                          lines=int(cell[cell.arm == a].ModelID.nunique()))
                                  for a in cell.arm.unique()}

    # ------------------------------------------------------------------ 5. strata (D135-S0(9))
    ess = pd.read_csv(DEP / "CRISPRInferredCommonEssentials.csv")
    ess_sym = set(ess.Essentials.astype(str).str.split(" ").str[0])
    prov["strata_sources"] = dict(common_essentials=len(ess_sym),
                                  target_genes_common_essential=int(sum(g in ess_sym for g in genes)))

    def dep_matrix(path, want):
        hdr = pd.read_csv(path, nrows=0)
        cmap = {}
        for c in hdr.columns[1:]:
            cmap.setdefault(c.split(" (")[0], c)
        use = [hdr.columns[0]] + [cmap[g] for g in want if g in cmap]
        d = pd.read_csv(path, usecols=use).rename(columns={hdr.columns[0]: "ModelID"})
        d.columns = ["ModelID"] + [c.split(" (")[0] for c in d.columns[1:]]
        return d.set_index("ModelID"), sorted(set(want) - set(cmap))

    expr, expr_missing = dep_matrix(DEP / "OmicsExpressionProteinCodingGenesTPMLogp1.csv", genes)
    dam, dam_missing = dep_matrix(DEP / "OmicsSomaticMutationsMatrixDamaging.csv", genes)
    hot, hot_missing = dep_matrix(DEP / "OmicsSomaticMutationsMatrixHotspot.csv", genes)
    prov["strata_sources"].update(
        expression_genes_missing=expr_missing, damaging_genes_missing=len(dam_missing),
        hotspot_genes_missing=len(hot_missing),
        expression_lines_missing=sorted(set(drug_lines) - set(expr.index)),
        note="a gene or line with no DepMap column/row is reported MISSING and never imputed (D135-S0(9))")

    def look(mat, g, l):
        try:
            return float(mat.at[l, g])
        except Exception:
            return np.nan

    # ------------------------------------------------------------------ 6. paired cells (D135-S0(6))
    rec = []
    for _, c in cell.iterrows():
        for d in gene_drugs.get(c.gene, []):
            if c.ModelID not in present.get(d, ()):
                continue
            rec.append(dict(arm=c.arm, drug=d, gene=c.gene, ModelID=c.ModelID,
                            n_pert_sigs=c.n_sigs, pert_nz2=c.pert_nz2, pert_nz4=c.pert_nz4,
                            pert_tas=c.pert_tas, pert_n_hiq=c.n_hiq,
                            target_on_gene_space=c.target_on_gene_space,
                            target_feature_space=c.target_feature_space,
                            target_z_min=c.target_z_min, target_z_max=c.target_z_max))
    pc = pd.DataFrame(rec)
    if len(pc):
        ann = dg_ann.groupby(["cmap_name", "gene"]).agg(
            orient=("orient", lambda s: "|".join(sorted(set(s.dropna().astype(str))))),
            actions=("actions", lambda s: "|".join(sorted(set(s.dropna().astype(str))))),
            tids=("tid", lambda s: "|".join(sorted(set(s))))).reset_index()
        pc = pc.merge(ann, left_on=["drug", "gene"], right_on=["cmap_name", "gene"], how="left").drop(
            columns=["cmap_name"])
        pc = pc.merge(dcl, left_on=["drug", "ModelID"], right_on=["cmap_name", "ModelID"],
                      how="left").drop(columns=["cmap_name"])
        pc["drug_n_target_genes"] = pc.drug.map(drug_genes.apply(len))
        pc["drug_n_units"] = pc.drug.map(drug_units)
        pc["essential"] = pc.gene.isin(ess_sym)
        pc["expr_log2tpm"] = [look(expr, g, l) for g, l in zip(pc.gene, pc.ModelID)]
        pc["expressed"] = pc.expr_log2tpm >= EXPR_CUT
        pc["expr_missing"] = pc.expr_log2tpm.isna()
        pc["mut_damaging"] = [look(dam, g, l) for g, l in zip(pc.gene, pc.ModelID)]
        pc["mut_hotspot"] = [look(hot, g, l) for g, l in zip(pc.gene, pc.ModelID)]
        oe_p = pc.arm == "trt_oe"
        for lab, cut in QC_CUTS.items():
            pc[f"qc_{lab}"] = np.where(pc.target_on_gene_space & ~oe_p,
                                       pc.target_z_min <= cut if cut < 0 else pc.target_z_min < 0, np.nan)
            pc[f"qc_oe_{lab.replace('le_-', 'ge_').replace('lt_0', 'gt_0')}"] = np.where(
                pc.target_on_gene_space & oe_p,
                pc.target_z_max >= -cut if cut < 0 else pc.target_z_max > 0, np.nan)
        pc = pc.sort_values(["arm", "gene", "ModelID", "drug"])
    pc.to_csv(OUT / "paired_cells.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ 7. attrition (D135-S0(6))
    attr = []
    ann_drug_lines = {}
    for g, ds in gene_drugs.items():
        s = set()
        for d in ds:
            s |= present.get(d, set())
        ann_drug_lines[g] = s
    for arm in KD_ARMS + ["trt_oe"]:
        a_si = si[si.pert_type == arm]
        a_m = a_si[a_si.ModelID.notna()]
        a_d = a_m[a_m.ModelID.isin(drug_lines)]
        ce = cell[cell.arm == arm]
        p = pc[pc.arm == arm] if len(pc) else pc
        base = dict(arm=arm)
        attr.append({**base, "step": "1 target genes (annotation)", "genes": len(genes),
                     "lines": None, "cells": None, "drugs": int(dg_ann.cmap_name.nunique())})
        attr.append({**base, "step": "2 gene has any signature of this arm",
                     "genes": int(a_si.cmap_name.nunique()), "lines": int(a_si.cell_iname.nunique()),
                     "cells": int(len(a_si)), "drugs": None})
        attr.append({**base, "step": "3 signature on a mapped cancer line",
                     "genes": int(a_m.cmap_name.nunique()), "lines": int(a_m.ModelID.nunique()),
                     "cells": int(len(a_m)), "drugs": None})
        attr.append({**base, "step": "4 signature on one of the 73 drug-matrix lines",
                     "genes": int(a_d.cmap_name.nunique()), "lines": int(a_d.ModelID.nunique()),
                     "cells": int(len(a_d)), "drugs": None})
        attr.append({**base, "step": "5 (gene, line) perturbation cells", "genes": int(ce.gene.nunique()),
                     "lines": int(ce.ModelID.nunique()), "cells": int(len(ce)), "drugs": None})
        ce2 = ce[[l in ann_drug_lines.get(g, set()) for g, l in zip(ce.gene, ce.ModelID)]]
        attr.append({**base, "step": "6 (gene, line) cells with >= 1 annotated drug profiled in that line",
                     "genes": int(ce2.gene.nunique()), "lines": int(ce2.ModelID.nunique()),
                     "cells": int(len(ce2)), "drugs": None})
        attr.append({**base, "step": "7 PAIRED (drug, gene, line) cells",
                     "genes": int(p.gene.nunique()) if len(p) else 0,
                     "lines": int(p.ModelID.nunique()) if len(p) else 0, "cells": int(len(p)),
                     "drugs": int(p.drug.nunique()) if len(p) else 0})
        if len(p):
            for f in NZ2_FLOORS:
                q = p[p.pert_nz2 >= f]
                attr.append({**base, "step": f"8 paired + CANDIDATE perturbation floor nz2 >= {f} (costed, "
                                             f"not chosen)", "genes": int(q.gene.nunique()),
                             "lines": int(q.ModelID.nunique()), "cells": int(len(q)),
                             "drugs": int(q.drug.nunique())})
            for f in TAS_FLOORS:
                q = p[p.pert_tas >= f]
                attr.append({**base, "step": f"9 paired + CANDIDATE perturbation floor tas >= {f} (costed, "
                                             f"not chosen)", "genes": int(q.gene.nunique()),
                             "lines": int(q.ModelID.nunique()), "cells": int(len(q)),
                             "drugs": int(q.drug.nunique())})
            q = p[p.pert_n_hiq >= 1]
            attr.append({**base, "step": "10 paired + CANDIDATE is_hiq >= 1 signature (costed, not chosen)",
                         "genes": int(q.gene.nunique()), "lines": int(q.ModelID.nunique()),
                         "cells": int(len(q)), "drugs": int(q.drug.nunique())})
    pd.DataFrame(attr).to_csv(OUT / "attrition.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ 8. strata distribution + bottom line
    srows, bottom = [], {}
    for arm in KD_ARMS:
        p = pc[pc.arm == arm] if len(pc) else pc
        if not len(p):
            bottom[arm] = dict(paired_cells=0)
            continue
        conds = {
            "L3 essential target (excused)": p.essential,
            "L8 target not expressed (log2TPM+1 < 1) (excused)": ~p.expressed & ~p.expr_missing,
            "L8 target expression MISSING": p.expr_missing,
            "L7 drug has > 1 annotated target gene (excused)": p.drug_n_target_genes > 1,
            "L7 drug has > 1 annotated target unit": p.drug_n_units > 1,
            f"L5 weak perturbation (nz2 < {ILLUS_PERT_NZ2}, illustrative)": p.pert_nz2 < ILLUS_PERT_NZ2,
            f"L6 weak drug response (nz2 < {ILLUS_DRUG_NZ2}, illustrative)": p.drug_nz2 < ILLUS_DRUG_NZ2,
            "L1 only one drug timepoint available": p.drug_n_timepoints <= 1,
            "target mRNA not on the gene space (QC impossible)": ~p.target_on_gene_space.astype(bool),
            "target mRNA QC fails at z <= -1": p["qc_z_le_-1"] == 0,
            "target gene mutated (damaging)": p.mut_damaging > 0,
            "target gene mutated (hotspot)": p.mut_hotspot > 0,
        }
        for k, v in conds.items():
            v = v.fillna(False).astype(bool) if hasattr(v, "fillna") else v
            srows.append(dict(arm=arm, stratum=k, cells=int(v.sum()),
                              frac=round(float(v.mean()), 4), drugs=int(p.loc[v, "drug"].nunique()),
                              genes=int(p.loc[v, "gene"].nunique())))
        structural = (~p.essential) & p.expressed & (p.drug_n_target_genes == 1)
        structural2 = (~p.essential) & p.expressed & (p.drug_n_target_genes <= 2)
        by_floor = {}
        for f in [0, 50, 100, 200]:
            m = structural & (p.pert_nz2 >= f) & (p.drug_nz2 >= ILLUS_DRUG_NZ2)
            m2 = structural2 & (p.pert_nz2 >= f) & (p.drug_nz2 >= ILLUS_DRUG_NZ2)
            lmq = m & (p.target_feature_space == "landmark") & (p["qc_z_le_-1"] == 1)
            by_floor[f"pert_nz2>={f}"] = dict(
                cells_1_target=int(m.sum()), genes_1_target=int(p.loc[m, "gene"].nunique()),
                drugs_1_target=int(p.loc[m, "drug"].nunique()),
                cells_le2_targets=int(m2.sum()),
                cells_landmark_and_qc_pass=int(lmq.sum()), genes_landmark_and_qc_pass=int(
                    p.loc[lmq, "gene"].nunique()))
        # cumulative attrition over the exculpatory strata, in a fixed order
        import pandas as _pd
        keep = _pd.Series(True, index=p.index)
        cum = {"0 paired cells": int(len(p))}
        for lab, v in [("1 drop common-essential targets", ~p.essential),
                       ("2 drop targets not expressed in the line", p.expressed),
                       ("3 drop drugs with > 1 annotated target gene", p.drug_n_target_genes == 1),
                       (f"4 drop weak drug profiles (nz2 < {ILLUS_DRUG_NZ2}, illustrative)",
                        p.drug_nz2 >= ILLUS_DRUG_NZ2),
                       (f"5 drop weak perturbations (nz2 < {ILLUS_PERT_NZ2}, illustrative)",
                        p.pert_nz2 >= ILLUS_PERT_NZ2)]:
            keep = keep & v.fillna(False).astype(bool)
            cum[lab] = int(keep.sum())
        clean = (~p.essential) & p.expressed & (p.drug_n_target_genes == 1) & \
                (p.pert_nz2 >= ILLUS_PERT_NZ2) & (p.drug_nz2 >= ILLUS_DRUG_NZ2)
        clean2 = (~p.essential) & p.expressed & (p.drug_n_target_genes <= 2) & \
                 (p.pert_nz2 >= ILLUS_PERT_NZ2) & (p.drug_nz2 >= ILLUS_DRUG_NZ2)
        qcf = p["qc_z_le_-1"] == 1
        bottom[arm] = dict(
            paired_cells=int(len(p)), drugs=int(p.drug.nunique()), genes=int(p.gene.nunique()),
            lines=int(p.ModelID.nunique()),
            structural_only_no_strength_floor=dict(
                cells=int(structural.sum()), drugs=int(p.loc[structural, "drug"].nunique()),
                genes=int(p.loc[structural, "gene"].nunique()),
                lines=int(p.loc[structural, "ModelID"].nunique()),
                cells_le2_targets=int(structural2.sum())),
            clean_by_perturbation_floor=by_floor, cumulative_attrition=cum,
            cells_with_pert_tas_ge_0p4=int((p.pert_tas >= 0.4).sum()),
            note_tas_asymmetry="the drug side already requires per-signature tas >= 0.4; this is how many "
                               "paired cells have a PERTURBATION that would clear the same bar",
            no_stratum_strict_1_target=dict(cells=int(clean.sum()), drugs=int(p.loc[clean, "drug"].nunique()),
                                            genes=int(p.loc[clean, "gene"].nunique()),
                                            lines=int(p.loc[clean, "ModelID"].nunique())),
            no_stratum_le2_targets=dict(cells=int(clean2.sum()),
                                        drugs=int(p.loc[clean2, "drug"].nunique()),
                                        genes=int(p.loc[clean2, "gene"].nunique()),
                                        lines=int(p.loc[clean2, "ModelID"].nunique())),
            no_stratum_strict_AND_mrna_qc=dict(cells=int((clean & qcf).sum()),
                                               genes=int(p.loc[clean & qcf, "gene"].nunique()),
                                               drugs=int(p.loc[clean & qcf, "drug"].nunique())),
            cells_per_line=p.ModelID.value_counts().to_dict(),
            cells_per_gene_top=p.gene.value_counts().head(15).to_dict(),
            orient_lower=int((p.orient.astype(str).str.contains("lower")).sum()),
            orient_raise=int((p.orient.astype(str).str.contains("raise")).sum()))
    pd.DataFrame(srows).to_csv(OUT / "strata_summary.tsv", sep="\t", index=False)
    prov["bottom_line"] = bottom
    prov["bottom_line_note"] = (
        "'no stratum' = gene not DepMap common-essential AND target expressed (log2TPM+1 >= 1) in that line "
        "AND the drug carries exactly one annotated target gene AND perturbation nz2 >= 100 AND drug nz2 "
        f">= {ILLUS_DRUG_NZ2}. The two nz2 numbers are ILLUSTRATIVE and are adopted by nothing; S1 chooses "
        "floors. The mRNA-QC variant is reported for CRISPR as a SEPARATE line and is not applied, per "
        "D135-S0(8). Limitation 2 (non-catalytic scaffolding roles) is NOT COMPUTABLE - no curated list "
        "exists in the repo - so it is MISSING from every count here and the clean-positive numbers are "
        "therefore an UPPER bound.")
    prov["limitation_2_missing"] = ("spec Section 6 limitation 2 (knockout removes scaffolding an inhibitor "
                                    "leaves intact) cannot be stratified: no curated list of targets with "
                                    "non-catalytic function is in the repo, and inventing one would be "
                                    "imputation. Reported missing, per D135-S0(9).")
    prov["constraints_held"] = ["no concordance computed", "no drug-to-perturbation comparison",
                                "no positive controls", "no AUROC", "no null", "no threshold chosen",
                                "no CTRP / PRISM / GDSC or any outcome label read",
                                "no production score changed, nothing re-ranked"]
    prov["threads"] = 8
    prov["runtime_s"] = round(time.time() - t0, 1)
    json.dump(prov, open(OUT / "d135_s0_run.json", "w"), indent=1)
    print(json.dumps({k: prov[k] for k in ["scope", "drug_universe", "spec_section2_check",
                                           "perturbation_cells", "bottom_line"]}, indent=1)[:12000])
    print(f"done in {prov['runtime_s']}s")


if __name__ == "__main__":
    main()
