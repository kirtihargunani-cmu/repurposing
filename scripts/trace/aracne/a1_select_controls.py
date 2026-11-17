#!/usr/bin/env python
"""TRACE-ARACNe stage 1: choose the five negative controls of PROCEDURE Section 7 FROM DATA.

Selection reads only: drug annotation, ChEMBL units, cohort/line matching, DepMap line
expression and mutations, the OmniPath route index and the ARACNe regulons. It reads NO
drug response profile, NO D128-A2/A4 output and NO outcome label, so the controls are
fixed before anything is observed about them (PROCEDURE Section 8).

Output: results/lincs/trace_aracne/control_selection.json  (the chosen pairs + why)
"""
import json, os, re, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
OUT = ROOT / "results/lincs/trace_aracne"
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"
SEED = 20261117          # frozen: the date, not tuned
HOP_CAP_NON_TF = 2       # same frozen P6 cap the CollecTRI trace used


def net_for(cohort, override):
    """Production rule (scripts/79e_drug_reversal.py:125 + D97): config/group_networks.tsv is an
    OVERRIDE map for the subtype cohorts; every other cohort strips its source prefix."""
    if cohort in override:
        return override[cohort]
    return "COADREAD" if cohort in ("TCGA-COAD", "TCGA-READ") else cohort.split("-", 1)[1]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ovr = pd.read_csv(ROOT / "config/group_networks.tsv", sep="\t").set_index("group_id").net.to_dict()
    units = pd.read_csv(ROOT / "results/lincs/target_role/d121/units.tsv", sep="\t",
                        dtype={"tid": str, "action_flag": str})
    annot = pd.read_csv(ROOT / "results/lincs/drug_annotations.tsv", sep="\t")
    idx = pd.read_csv(ROOT / "data/processed/lincs/drug_line_profiles_index.tsv", sep="\t")
    lw = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    rank = pd.read_csv(ROOT / "results/lincs/fda_ranking_v2_consensus.tsv", sep="\t")
    genes = pd.read_csv(ROOT / "data/processed/lincs/drug_line_profiles_genes.tsv", sep="\t")
    sym2ens = dict(zip(genes.symbol, genes.gene_id))
    P = pd.read_parquet(ROOT / "results/lincs/d128/a1/target_tf_paths.parquet")
    P = P[P.network == "omnipath"]
    E = pd.read_csv(ROOT / "data/processed/lincs/trace_aracne/aracne_regulon_edges.tsv.gz", sep="\t")
    E = E[E.in_lincs_gene_space]
    regs_by_net = E.groupby("network").regulator.apply(set).to_dict()
    size_by_net = {n: d.groupby("regulator").size().to_dict() for n, d in E.groupby("network")}

    # ---- the pool: single-gene ChEMBL units whose drug is profiled in a celligner-matched line ----
    u1 = units[units.k == 1].copy()
    pid = annot[["pert_id", "cmap_name"]].drop_duplicates().merge(idx, on="pert_id")
    best = pid.sort_values("n_sig", ascending=False).drop_duplicates(["cmap_name", "ModelID"])
    lwm = lw[lw.in_key].sort_values("weight", ascending=False)
    pool = (u1.merge(best, on="cmap_name").merge(lwm[["group_id", "ModelID", "weight", "CellLineName"]],
                                                 on="ModelID"))
    pool["net"] = [net_for(c, ovr) for c in pool.group_id]
    pool = pool[pool.net.isin(regs_by_net)].copy()   # cohort must have an ARACNe network
    print(f"pool: {len(pool)} (unit x cohort x line) rows, {pool.cmap_name.nunique()} drugs, "
          f"{pool.group_id.nunique()} cohorts")

    # ---- how many ARACNe regulators each target reaches in its cohort's network ----
    def reach(target, net):
        if net not in regs_by_net:
            return []
        is_reg = target in regs_by_net[net]
        cap = 0 if is_reg else HOP_CAP_NON_TF
        pp = P[(P.target_gene == target) & (P.hops <= cap) & (P.path_status == "resolved")]
        return sorted({t for t in pp.tf if t in regs_by_net[net]})

    # ---- DepMap line expression + mutations (selection inputs, not response) ----
    ex = pd.read_csv(ROOT / "data/raw/depmap_24q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                     index_col=0, low_memory=False)
    ex.columns = [c.split(" (")[0] for c in ex.columns]
    mut = pd.read_csv(ROOT / "data/raw/depmap_24q4/OmicsSomaticMutations.csv", low_memory=False,
                      usecols=["ModelID", "HugoSymbol", "ProteinChange", "LikelyLoF", "Hotspot", "OncogeneHighImpact"])

    chosen, notes = {}, {}

    # ---------- NC1: target not expressed in the source line -> L5 must break ----------
    cand = pool[pool.genes.isin(ex.columns) & pool.ModelID.isin(ex.index)].copy()
    cand["tpm"] = [float(ex.at[r.ModelID, r.genes]) for r in cand.itertuples()]
    cand["n_reg"] = [len(reach(r.genes, r.net)) for r in cand.itertuples()]
    nc1 = cand[(cand.tpm < 0.25) & (cand.n_reg >= 3)].sort_values(
        ["n_reg", "weight"], ascending=False)
    assert len(nc1), "no unexpressed-target candidate"
    r = nc1.iloc[0]
    chosen["NC1_target_not_expressed"] = dict(
        control="target not expressed in the source line; L5 must break",
        cmap_name=r.cmap_name, target=r.genes, group_id=r.group_id, net=r.net,
        ModelID=r.ModelID, line=r.CellLineName, unit_id=f"{r.cmap_name}::{r.tid}",
        orient=r.orient, n_aracne_regulators=int(r.n_reg))
    notes["NC1_target_not_expressed"] = (
        f"{r.genes} log2(TPM+1) = {r.tpm:.3f} in {r.CellLineName} (DepMap 24Q4), i.e. not expressed, "
        f"while the target still reaches {int(r.n_reg)} ARACNe regulators in {r.net} so L1-L4 can run. "
        f"Chosen as the candidate with the most reachable regulators among {len(nc1)} rows with TPM<0.25 "
        f"and >=3 regulators -- deliberately the HARDEST such case for the procedure.")

    # ---------- NC2: RAF inhibitor in BRAF-wild-type, RAS-active ----------
    rafs = sorted(set(pool[pool.genes.isin(["BRAF", "RAF1", "ARAF"])].cmap_name))
    bm = mut[mut.HugoSymbol == "BRAF"]
    braf_mut = set(bm[bm.ProteinChange.fillna("").str.contains("V600")].ModelID)
    ras = mut[mut.HugoSymbol.isin(["KRAS", "NRAS", "HRAS"]) & (mut.Hotspot.fillna(False).astype(bool) | mut.OncogeneHighImpact.fillna(False).astype(bool))]
    ras_act = set(ras.ModelID)
    nc2c = pool[pool.genes.isin(["BRAF", "RAF1", "ARAF"]) & pool.ModelID.isin(ras_act)
                & ~pool.ModelID.isin(braf_mut)]
    nc2c = nc2c.copy(); nc2c["n_reg"] = [len(reach(r.genes, r.net)) for r in nc2c.itertuples()]
    nc2c = nc2c[nc2c.n_reg >= 1].sort_values(["n_reg", "weight"], ascending=False)
    assert len(nc2c), "no RAF-paradox candidate"
    r = nc2c.iloc[0]
    rs = ras[ras.ModelID == r.ModelID].iloc[0]
    chosen["NC2_raf_paradox"] = dict(
        control="RAF inhibitor in BRAF-wild-type, RAS-active context; the chain must not close",
        cmap_name=r.cmap_name, target=r.genes, group_id=r.group_id, net=r.net,
        ModelID=r.ModelID, line=r.CellLineName, unit_id=f"{r.cmap_name}::{r.tid}",
        orient=r.orient, n_aracne_regulators=int(r.n_reg))
    notes["NC2_raf_paradox"] = (
        f"{r.CellLineName} carries no BRAF V600 mutation and a hotspot {rs.HugoSymbol} "
        f"{rs.ProteinChange} (DepMap 24Q4 somatic calls), the paradoxical-activation context in which a "
        f"RAF inhibitor activates rather than inhibits ERK signalling. The predicted footprint sign is "
        f"therefore expected to be WRONG, not merely weak. {len(nc2c)} candidate rows; "
        f"most reachable regulators chosen.")

    # ---------- NC3: random drug-cohort pairing ----------
    rng = np.random.default_rng(SEED)
    sh = pool.drop_duplicates(["cmap_name", "genes"]).reset_index(drop=True)
    cohorts = sorted(pool.group_id.unique())
    for _ in range(500):
        i = int(rng.integers(len(sh))); c = cohorts[int(rng.integers(len(cohorts)))]
        row = sh.iloc[i]
        if row.group_id == c:
            continue
        net = net_for(c, ovr)
        ln = lwm[lwm.group_id == c]
        ln = ln[ln.ModelID.isin(set(best[best.cmap_name == row.cmap_name].ModelID))]
        if not len(ln) or len(reach(row.genes, net)) < 1:
            continue
        ln = ln.iloc[0]
        chosen["NC3_random_pairing"] = dict(
            control="drug paired with a cohort it was not matched to; L2 or L7 must break",
            cmap_name=row.cmap_name, target=row.genes, group_id=c, net=net,
            ModelID=ln.ModelID, line=ln.CellLineName, unit_id=f"{row.cmap_name}::{row.tid}",
            orient=row.orient, n_aracne_regulators=len(reach(row.genes, net)))
        notes["NC3_random_pairing"] = (
            f"seed {SEED}, first draw satisfying: cohort != the drug's own matched cohort "
            f"({row.group_id}), the drug has a profile in a celligner line of the drawn cohort, and the "
            f"target reaches >=1 ARACNe regulator there. No property of {c} was used to pick it.")
        break
    assert "NC3_random_pairing" in chosen

    # ---------- NC4: a cytotoxic, where the "mechanism" is cell death ----------
    # dna_acting in the production ranking covers only 13 alkylators, none of which has a single-gene
    # ChEMBL unit in the pool, so the class is taken from the curated MOA strings instead: the classical
    # cytotoxic chemotherapy classes, whose transcriptional response is growth arrest and death.
    dna = set(rank[rank.dna_acting == True].cmap_name)
    CYTO = re.compile(r"topoisomerase|tubulin|microtubule|DNA alkylating|antimetab|thymidylate|"
                      r"dihydrofolate|ribonucleotide reductase|DNA inhibitor", re.I)
    moa = annot[["cmap_name", "hub_moa", "lincs_moa"]].drop_duplicates("cmap_name")
    nc4c = pool.merge(moa, on="cmap_name", how="left")
    is_cyto = (nc4c.hub_moa.fillna("").str.contains(CYTO) | nc4c.lincs_moa.fillna("").str.contains(CYTO)
               | nc4c.cmap_name.isin(dna))
    # a drug that ALSO has a targeted mechanism (kinase/receptor/etc.) is not a clean cytotoxicity
    # control, because part of the chain is then supposed to close. Keep only MOAs that are cytotoxic
    # and nothing else.
    # The classical cytotoxic effectors - TOP1, TOP2A, TYMS, DHFR, GART, RRM1, TUBB - have ZERO
    # OmniPath route rows: they are not in the signalling network, so the trace cannot start on a
    # chemotherapeutic acting through its own effector. That is recorded below as a structural finding.
    # The realistic cytotoxicity failure mode is therefore a drug that KILLS but carries a targeted
    # annotation that does route to regulators: the chain can close on the annotation while the drug's
    # actual effect is death. Rule: every MOA term cytotoxic, target reaches >=1 ARACNe regulator,
    # most reachable regulators first.
    EFFECTORS = ["TOP1", "TOP2A", "TYMS", "DHFR", "GART", "RRM1", "TUBB"]
    unroutable = [g for g in EFFECTORS if not len(P[P.target_gene == g])]

    def pure_cytotoxic(row):
        terms = [t.strip() for f in (row.hub_moa, row.lincs_moa)
                 for t in str(f if pd.notna(f) else "").split("|") if t.strip()]
        return bool(terms) and all(CYTO.search(t) for t in terms)
    nc4c = nc4c[is_cyto].copy()
    nc4c["pure"] = [pure_cytotoxic(r) for r in nc4c.itertuples()]
    nc4c["n_reg"] = [len(reach(r.genes, r.net)) for r in nc4c.itertuples()]
    nc4c = nc4c[(nc4c.n_reg >= 1) & nc4c.pure].sort_values(["n_reg", "weight"], ascending=False)
    assert len(nc4c), "no cytotoxic candidate"
    r = nc4c.iloc[0]
    chosen["NC4_cytotoxic"] = dict(
        control="cytotoxic drug with a routable targeted annotation; L6 and the cytotoxicity "
                "alternative must fire rather than the chain closing on the annotation",
        cmap_name=r.cmap_name, target=r.genes, group_id=r.group_id, net=r.net,
        ModelID=r.ModelID, line=r.CellLineName, unit_id=f"{r.cmap_name}::{r.tid}",
        orient=r.orient, n_aracne_regulators=int(r.n_reg), moa=str(r.hub_moa),
        effectors_unroutable=unroutable)
    notes["NC4_cytotoxic"] = (
        f"STRUCTURAL FINDING FIRST: {', '.join(unroutable)} have zero OmniPath route rows, so a "
        f"classical chemotherapeutic acting through its own effector cannot be traced at all - it would "
        f"be 'not evaluable' at L3, which is NC5's test, not a cytotoxicity test. "
        f"The cytotoxicity control therefore uses a drug whose MOA is purely cytotoxic "
        f"('{r.hub_moa}') but whose annotated unit ({r.genes}) does route to {int(r.n_reg)} ARACNe "
        f"regulators in {r.net}: {r.cmap_name} in {r.CellLineName}. The chain is able to close here, and "
        f"if it does, the procedure cannot tell a targeted mechanism from cell death. "
        f"{len(nc4c)} such rows in the pool; most reachable regulators chosen.")

    # ---------- NC5: target with no anchored TF -> L3 must say 'not evaluable' ----------
    nc5c = pool.copy()
    nc5c["n_reg"] = [len(reach(r.genes, r.net)) for r in nc5c.itertuples()]
    # a target with ZERO OmniPath rows is trivially not evaluable; the sharp test is a target the
    # signalling network DOES know about, which still anchors no regulator in this cohort's ARACNe net.
    nc5c["n_routes"] = [int((P.target_gene == g).sum()) for g in nc5c.genes]
    nc5c = nc5c[(nc5c.n_reg == 0) & (nc5c.n_routes > 0)].sort_values(
        ["n_routes", "weight"], ascending=False)
    assert len(nc5c), "no unanchored-target candidate"
    r = nc5c.iloc[0]
    n_any = int(((P.target_gene == r.genes)).sum())
    chosen["NC5_no_anchored_tf"] = dict(
        control="target reaches no regulator; L3 must return 'not evaluable' rather than a path",
        cmap_name=r.cmap_name, target=r.genes, group_id=r.group_id, net=r.net,
        ModelID=r.ModelID, line=r.CellLineName, unit_id=f"{r.cmap_name}::{r.tid}",
        orient=r.orient, n_aracne_regulators=0)
    notes["NC5_no_anchored_tf"] = (
        f"{r.genes} HAS {n_any} OmniPath route rows, but none resolves to a regulator with an "
        f"ARACNe regulon in {r.net} within the frozen hop cap. {len(nc5c)} such rows in the pool; the "
        f"highest-celligner-weight line chosen. The correct behaviour is an empty footprint and a "
        f"'not evaluable' L3, NOT a footprint assembled from ambiguous routes.")

    payload = dict(seed=SEED, hop_cap_non_tf=HOP_CAP_NON_TF, pool_rows=int(len(pool)),
                   controls=chosen, why=notes,
                   reads_no_response_data=True, reads_no_outcome_label=True)
    (OUT / "control_selection.json").write_text(json.dumps(payload, indent=1, default=str))
    for k, v in chosen.items():
        print(f"\n{k}: {v['cmap_name']} / {v['target']} / {v['group_id']} / {v['line']} "
              f"({v['ModelID']}), {v['n_aracne_regulators']} regulators")
        print("   ", notes[k])
    print("\nwrote", OUT / "control_selection.json")


if __name__ == "__main__":
    sys.exit(main())
