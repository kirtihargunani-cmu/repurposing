#!/usr/bin/env python
"""Quick check (reported only, 2026-11-08): drugs aimed at targets every cell needs (chemotherapy). Runs in `celligner`.

Question: why does the layered score not recover where these drugs are approved, and would response biomarkers?
Drugs: anticancer FDA drugs whose D122 rows carry essential_target (a main-target gene is a DepMap common essential),
plus anticancer drugs without a protein target (alkylating / platinum), grouped into classes by main target:
TOP2 poisons, TOP1 inhibitors, TYMS antimetabolites, tubulin binders, nucleotide-synthesis / polymerase, DNA-damaging
(no protein target). Cohort-level markers (TCGA tumors): proliferation index (mean log2 median TPM of 10 proliferation
genes), SLFN11, TOP2A, TYMS, TUBB3 (median TPM), HR-loss share (truncating BRCA1 / BRCA2 / PALB2, MC3).
Tests (Hub and ChEMBL labels, excluded rows left out):
  W  within-cancer pooled AUC (79l) of layered_pct and D / A / O among essential-target anticancer drugs
  B  between-cancer: per drug, cancers where it is approved vs the other cancers where it was scored; pooled AUC over
     drugs (U statistics summed), 1000 within-drug permutations; scores: layered_pct, D, A, O, proliferation, and
     each class's expected marker (+: SLFN11, HR loss for DNA-damaging classes, TOP2A for TOP2; -: TYMS for TYMS
     drugs, TUBB3 for tubulin binders; the sign is applied before the AUC).
Outputs results/lincs/essential_check/{markers.tsv, drugs.tsv, within.tsv, between.tsv}
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results/lincs"
OUT = RES / "essential_check"
sys.path.insert(0, str(ROOT / "scripts"))
import evidence_flags  # noqa: E402
evidence_flags.STRICT = False
from fda_rows import fda_rows  # noqa: E402

_s = importlib.util.spec_from_file_location("ku", ROOT / "scripts/79l_known_use.py"); KU = importlib.util.module_from_spec(_s); _s.loader.exec_module(KU)
PROLIF = ["MKI67", "TOP2A", "CCNB1", "CDK1", "BUB1", "AURKA", "MCM2", "PCNA", "CENPF", "TPX2"]
MARK = ["SLFN11", "TOP2A", "TYMS", "TUBB3"]
SEED, N_PERM = 0, 1000


def drug_class(genes, no_target):
    g = set(genes)
    if no_target:
        return "DNA-damaging (no protein target)"
    if g & {"TOP2A", "TOP2B"}:
        return "TOP2 poison"
    if "TOP1" in g:
        return "TOP1 inhibitor"
    if "TYMS" in g:
        return "TYMS antimetabolite"
    if any(x.startswith(("TUBA", "TUBB")) for x in g):
        return "tubulin binder"
    if g & {"RRM1", "RRM2", "RRM2B"} or any(x.startswith("POL") for x in g):
        return "nucleotide / polymerase"
    return "other essential"


EXPECT = {"DNA-damaging (no protein target)": [("SLFN11", 1), ("hr_loss", 1)], "TOP2 poison": [("SLFN11", 1), ("TOP2A", 1)],
          "TOP1 inhibitor": [("SLFN11", 1), ("hr_loss", 1)], "TYMS antimetabolite": [("TYMS", -1)],
          "tubulin binder": [("TUBB3", -1)], "nucleotide / polymerase": [("SLFN11", 1)], "other essential": []}


def markers(groups):
    b = importlib.util.spec_from_file_location("b79", ROOT / "scripts/79b_target_gate.py"); b79 = importlib.util.module_from_spec(b); b.loader.exec_module(b79)
    idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t").query("primary_for_hgnc").drop_duplicates("symbol")
    s2g = dict(zip(idm.symbol, idm.gene_id)); g2s = {v: k for k, v in s2g.items()}
    ids = {s2g[x] for x in set(PROLIF + MARK) if x in s2g}
    ex = b79.group_expression(ids, groups, 1.0)
    ex["symbol"] = ex.gene_id.map(g2s)
    M = ex.pivot_table(index="group_id", columns="symbol", values="median_tpm")
    out = pd.DataFrame(index=M.index)
    out["proliferation"] = np.log2(M[[x for x in PROLIF if x in M]] + 1).mean(1)
    for x in MARK:
        out[x] = M.get(x)
    # HR loss: truncating BRCA1 / BRCA2 / PALB2 in MC3, share of the cohort's profiled tumors
    sg = pd.read_csv(ROOT / "results/step1/sample_groups_toil.tsv", sep="\t")
    sg = sg[(sg.role == "tumor") & sg.group_id.isin(groups)]
    maf = pd.read_csv(ROOT / "data/raw/tcga_pancan/mc3.v0.2.8.PUBLIC.maf.gz", sep="\t", low_memory=False,
                      usecols=["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", "FILTER"])
    maf["sample"] = maf.Tumor_Sample_Barcode.str[:15]
    prof = set(maf["sample"])
    tr = maf[maf.Hugo_Symbol.isin(["BRCA1", "BRCA2", "PALB2"]) & maf.Variant_Classification.isin(
        ["Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"])]
    hit = set(tr["sample"])
    s = sg[sg["sample"].isin(prof)]
    hr = s.assign(h=s["sample"].isin(hit)).groupby("group_id").h.agg(["mean", "size"])
    out["hr_loss"] = hr["mean"].where(hr["size"] >= 10)
    return out


def pooled(groups_of, y, s):
    """sum of within-group U / sum of pairs; groups_of: array of group keys."""
    num = den = 0.0
    for k in np.unique(groups_of):
        m = groups_of == k
        yy, ss = y[m], s[m]
        ok = ~np.isnan(ss)
        yy, ss = yy[ok], ss[ok]
        p, n = yy.sum(), (~yy).sum()
        if p == 0 or n == 0:
            continue
        r = pd.Series(ss).rank().values
        num += r[yy].sum() - p * (p + 1) / 2
        den += p * n
    return (num / den if den else np.nan), den


def main():
    rng = np.random.default_rng(SEED)
    OUT.mkdir(parents=True, exist_ok=True)
    g, _ = fda_rows()
    groups = sorted(g.group_id.unique())
    MK = markers(groups)
    MK.to_csv(OUT / "markers.tsv", sep="\t", float_format="%.4g")
    lab = pd.read_csv(RES / "known_use/labels.tsv", sep="\t")
    cols = ["group_id", "fda_ingredient", "layered_pct", "layered_D", "layered_A", "layered_O", "layered_genes", "essential_target",
            "no_protein_target", "evidence_targets", "eligible"]
    X = g[[c for c in cols if c in g]].rename(columns={"group_id": "cohort", "eligible": "eligible_d122"})
    X = X.merge(lab[["cohort", "fda_ingredient", "drug", "anticancer", "lab_hub", "lab_chembl", "excl_hub", "excl_chembl"]],
                on=["cohort", "fda_ingredient"], how="inner")
    X = X[X.anticancer.astype(bool)]
    X["no_target"] = X.no_protein_target.astype(str).eq("True")
    X["essential"] = X.essential_target.astype(str).eq("True") | X.no_target
    X["main_genes"] = X.evidence_targets.fillna("").map(lambda s: [t.split(":")[0] for t in s.split(";") if t])
    X["cls"] = [drug_class(gs, nt) for gs, nt in zip(X.main_genes, X.no_target)]
    X = X.join(MK, on="cohort")
    E = X[X.essential]
    E.drop(columns=["main_genes"]).to_csv(OUT / "drugs.tsv", sep="\t", index=False, float_format="%.4g")
    print("essential / DNA-damaging anticancer drugs:", E.drug.nunique(), "rows", len(E))
    print(E.drop_duplicates("drug").groupby("cls").drug.apply(lambda s: ", ".join(sorted(s))).to_string())

    rows = []
    for src in ("hub", "chembl"):
        y = f"lab_{src}"
        A = X[~X[f"excl_{src}"].astype(bool) & X.eligible_d122.astype(bool)]
        Ae = A[A.essential & ~A.no_target]
        for sc in ("layered_pct", "layered_D", "layered_A", "layered_O"):
            r = KU.pooled_auc(Ae, y, sc)
            rows.append(dict(source=src, test="W: essential-target drugs vs each other, within cancer", score=sc, **r))
        rows.append(dict(source=src, test="W: all anticancer drugs (reference)", score="layered_pct", **KU.pooled_auc(A, y, "layered_pct")))
    W = pd.DataFrame(rows)
    W.to_csv(OUT / "within.tsv", sep="\t", index=False)

    brow = []
    for src in ("hub", "chembl"):
        y = f"lab_{src}"
        B = E[~E[f"excl_{src}"].astype(bool)].copy()
        B["y"] = B[y].astype(bool)
        tests = [("all", "layered_pct", 1), ("all", "layered_D", 1), ("all", "layered_A", 1), ("all", "layered_O", 1),
                 ("all", "proliferation", 1)]
        for cls, ex in EXPECT.items():
            for m, sign in ex:
                tests.append((cls, m, sign))
        for cls, m, sign in tests:
            D = B if cls == "all" else B[B.cls == cls]
            if D.empty:
                continue
            s = D[m].astype(float).values * sign
            yy, gg = D.y.values, D.drug.values
            obs, pairs = pooled(gg, yy, s)
            null = []
            for _ in range(N_PERM):
                yp = yy.copy()
                for k in np.unique(gg):
                    i = np.where(gg == k)[0]
                    yp[i] = rng.permutation(yp[i])
                null.append(pooled(gg, yp, s)[0])
            null = np.array(null)
            n_drugs = int(sum(1 for k in np.unique(gg) if 0 < yy[gg == k].sum() < (gg == k).sum()))
            brow.append(dict(source=src, drugs=cls, score=m + (" (inverse)" if sign < 0 else ""), auc=round(obs, 4), pairs=int(pairs),
                             n_drugs=n_drugs, n_approved_pairs=int(yy.sum()),
                             p_perm=round((1 + np.sum(null >= obs)) / (1 + N_PERM), 4) if np.isfinite(obs) else np.nan))
    BT = pd.DataFrame(brow)
    BT.to_csv(OUT / "between.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print(W[["source", "test", "score", "auc", "n_cohorts", "n_pos", "n_neg"]].to_string(index=False))
    print(BT.to_string(index=False))
    print(MK.round(2).sort_values("proliferation").to_string())


if __name__ == "__main__":
    main()
