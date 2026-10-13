#!/usr/bin/env python
"""Split each reversal score into a proliferation part and a mechanism part (runs in `celligner`).

Why (user 2026-10-12): proliferation reversal is shared by all tumors (caveat C4: 44 growth-signaling inhibitors are top
10% by rev in >= 20 of 41 groups) and confounds lineage regulators (C3: anti-estrogens rank higher in basal than in
luminal A). The user's point: many effective drugs work by stopping proliferation, so that part is kept as its own
score, not discarded. Diagnostic run: all scores are reported side by side; nothing replaces rev / spec yet.
Proliferation genes: HALLMARK_E2F_TARGETS + HALLMARK_G2M_CHECKPOINT (303 on the LINCS profiles).
Proliferation-linked regulators of a network (fixed before looking at drug results):
  P1  targets enriched for proliferation genes (scripts/79h: hypergeometric FDR < 0.05, fold >= 2)
  P2  drug-induced activity tracks the profile's proliferation score: |Pearson r| >= R_MIN over all drug x line
      profiles (specific version, as scored); catches regulators whose VIPER activity moves with growth arrest
      without proliferation genes among their targets
  prolif_linked = P1 or P2. FOXA1 / GATA3 in BRCA are neither (r -0.43 / -0.23), hence the adjustment below.
Proliferation adjustment: per network and regulator, the activity is regressed on the profile's proliferation score
(OLS over all profiles) and the residual kept: A_adj = A - beta * pi. Removes the growth-arrest part of every
regulator's activity (also the part that is genuine biology, e.g. ER driving proliferation in luminal tumors).
Scores per disease (41 groups) and per Step 10 stratum, same weights / direct-borrowed split / formula as 79e / 79g:
  rev               as before (all MRs, activity A)                              [copied from 79e / 79g]
  rev_prolif        proliferation-linked MRs only, activity A                    "reverses proliferation"
  rev_nonprolif     the other MRs, activity A                                    exclusion only
  rev_mech          the other MRs, activity A_adj                                exclusion + adjustment
  prolif_effect     minus the weighted proliferation score of the drug's raw profiles (z units; > 0 = lowers
                    proliferation genes); disease-independent apart from the line weights
  spec_mech         rev_mech robust z minus the mean over reference diseases (other tissue keys, >= MIN_MRS_REF
                    non-proliferation MRs), each scored the same way;  strata: spec_sib_mech vs the sibling strata
  prolif_share      |rev_prolif| * sqrt(n_p) / (|rev_prolif| * sqrt(n_p) + |rev_nonprolif| * sqrt(n_np)): share of
                    the summed MR reversal that comes from proliferation-linked MRs
  *_pct             percentile within the disease / stratum; tier_mech / tier_spec_mech as tier_rev in 79e
Outputs: results/lincs/prolif_regulators_p2.tsv (P2 correlations), drug_reversal_split.tsv,
         stratum_reversal_split.tsv, split_controls.tsv (known-drug percentiles under every score),
         split_overview.tsv (generic-drug counts, MR counts, control means per score)
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed/lincs"
RES = ROOT / "results/lincs"
MR_FDR, TOP_PCT, MIN_MRS_REF, MIN_REF, R_MIN = 0.10, 10, 20, 10, 0.5
CONTROLS = [("BRCA-LumA", "high", ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]),
            ("BRCA-LumB", "high", ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]),
            ("TCGA-PRAD", "high", ["bicalutamide", "enzalutamide", "flutamide", "nilutamide"]),
            ("TCGA-COAD", "high", ["fluorouracil", "irinotecan", "oxaliplatin"]),
            ("TCGA-LUAD", "high", ["erlotinib", "gefitinib", "afatinib", "osimertinib"]),
            ("BRCA-Basal", "low", ["fulvestrant", "tamoxifen", "raloxifene"])]
ANTI_ER = ["fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene"]
STRATUM_CONTROLS = [("BRCA-LumA-PC3-low", "ER end", ANTI_ER), ("BRCA-Her2-PC3-low", "ESR1 / XBP1 end", ANTI_ER),
                    ("BRCA-LumA-PC3-high", "proliferation end",
                     ["palbociclib", "abemaciclib", "ribociclib", "paclitaxel", "docetaxel", "vinorelbine"])]
SCORES = ["rev", "rev_prolif", "rev_nonprolif", "rev_mech", "prolif_effect", "spec", "spec_mech"]


def net_of(cohort):
    return "COADREAD" if cohort in ("TCGA-COAD", "TCGA-READ") else cohort.split("-", 1)[1]


def rz(v):
    return (v - np.median(v)) / (1.4826 * np.median(np.abs(v - np.median(v))))


def pct(v):
    return (pd.Series(v).rank(method="max", pct=True).values * 100).round(1)


def mr_rev(D, regs, m):
    c = regs.get_indexer(m.regulator)
    return -(D[:, c] * np.sign(m.nes.values)).sum(1) / np.sqrt(max(len(c), 1))


def main():
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    idx["col"] = np.arange(len(idx))
    by_drug = idx.groupby("pert_id")
    genes = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
    gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:]
           for l in open(ROOT / "data/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt")}
    pm = genes.symbol.isin(set(gmt["HALLMARK_E2F_TARGETS"]) | set(gmt["HALLMARK_G2M_CHECKPOINT"])).values
    pi_raw = np.fromfile(PROC / "drug_line_profiles_bing.f32", dtype=np.float32).reshape(len(idx), -1)[:, pm].mean(1)
    pi_spec = np.fromfile(PROC / "drug_line_profiles_bing_specific.f32",
                          dtype=np.float32).reshape(len(idx), -1)[:, pm].mean(1).astype(np.float64)
    print(f"proliferation genes on profiles: {pm.sum()}")

    grp = pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t").set_index("group_id")
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
    models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", index_col=0)
    models = models[models.OncotreePrimaryDisease != "Non-Cancerous"]
    gmr = pd.read_csv(ROOT / "results/step6/aracne_project/group_mr.tsv", sep="\t")
    smr = pd.read_csv(ROOT / "results/step10/stratum_mr.tsv", sep="\t")
    p1 = pd.read_csv(RES / "prolif_regulators.tsv", sep="\t")
    p1 = p1[p1.prolif_linked].groupby("net").regulator.apply(set)
    gmr = gmr[gmr.fdr < MR_FDR]
    smr = smr[smr.fdr < MR_FDR]

    nets, p2_rows = {}, []
    def load(net):
        if net not in nets:
            regs = pd.Index((PROC / "drug_viper_specific" / f"{net}.regulators.txt").read_text().split())
            A = np.fromfile(PROC / "drug_viper_specific" / f"{net}.f32",
                            dtype=np.float32).reshape(len(idx), len(regs)).astype(np.float64)
            pc = pi_spec - pi_spec.mean()
            Ac = A - A.mean(0)
            beta = Ac.T @ pc / (pc @ pc)
            r = beta * pc.std() / Ac.std(0)
            A_adj = A - np.outer(pi_spec, beta)
            linked = pd.Index(regs).isin(p1.get(net, set())) | (np.abs(r) >= R_MIN)
            p2_rows.append(pd.DataFrame({"net": net, "regulator": regs, "r_prolif": r.round(3),
                                         "p1": pd.Index(regs).isin(p1.get(net, set())),
                                         "p2": np.abs(r) >= R_MIN, "prolif_linked": linked}))
            nets[net] = (regs, A, A_adj, set(regs[linked]))
        return nets[net]

    def activity(drugs, sources, key_row, wg, A, A_adj):
        km = models[models.OncotreeLineage.isin(str(key_row.lineages).split(";"))]
        if isinstance(key_row.diseases, str) and key_row.diseases:
            km = km[km.OncotreePrimaryDisease.isin(key_row.diseases.split(";"))]
        D, Da = np.empty((len(drugs), A.shape[1])), np.empty((len(drugs), A.shape[1]))
        pe = np.empty(len(drugs))
        for i, (d, src) in enumerate(zip(drugs, sources)):
            q = by_drug.get_group(d)
            q = q[q.ModelID.isin(km.index)] if src == "direct" else q[~q.ModelID.isin(km.index)]
            ww = wg.reindex(q.ModelID).fillna(0).values if src == "direct" else np.zeros(len(q))
            ww = ww if ww.sum() > 0 else np.ones(len(q))
            ww = ww / ww.sum()
            D[i], Da[i], pe[i] = ww @ A[q.col.values], ww @ A_adj[q.col.values], -(ww @ pi_raw[q.col.values])
        return D, Da, pe

    def split_scores(D, Da, regs, linked, m):
        isp = m.regulator.isin(linked)
        mp, mn = m[isp], m[~isp]
        rp, rn, rm = mr_rev(D, regs, mp), mr_rev(D, regs, mn), mr_rev(Da, regs, mn)
        a, b = np.abs(rp) * np.sqrt(len(mp)), np.abs(rn) * np.sqrt(len(mn))
        return dict(rev_prolif=rp, rev_nonprolif=rn, rev_mech=rm, prolif_share=a / (a + b + 1e-12),
                    n_mrs_prolif=len(mp), n_mrs_nonprolif=len(mn)), mn

    def refs_mech(Da, regs, linked, own_key):
        out = []
        for g2, m2 in gmr.groupby("group_id"):
            if g2 not in keys.index or keys.loc[g2, "key"] == own_key:
                continue
            m2 = m2[m2.regulator.isin(regs) & ~m2.regulator.isin(linked)]
            if len(m2) >= MIN_MRS_REF:
                out.append(rz(mr_rev(Da, regs, m2)))
        return out

    def finish(o, cols):
        for c in ["rev_prolif", "rev_nonprolif", "rev_mech", "prolif_effect", "spec_mech"] + cols:
            o[c + "_pct"] = pct(o[c].values)
        for t, s, c in [("tier_mech", "rev_mech", "rev_mech_pct"), ("tier_spec_mech", "rev_mech", "spec_mech_pct")] + \
                       [(f"tier_{x}", "rev_mech", f"{x}_pct") for x in cols]:
            hit = (o[s] > 0) & (o[c] >= 100 - TOP_PCT)
            o[t] = np.where(hit & o.passes_gate & o.approved, 1, np.where(hit & o.passes_gate, 2, 3))
        o["prolif_driven"] = (o.prolif_share > 0.5) & (o.rev_pct >= 100 - TOP_PCT) & (o.rev_mech_pct < 100 - TOP_PCT)
        return o

    # 41 groups
    par = pd.read_csv(RES / "drug_reversal.tsv", sep="\t", low_memory=False)
    w = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    keep = ["group_id", "pert_id", "cmap_name", "source", "passes_gate", "approved", "multi_target", "hub_moa",
            "n_mrs", "low_mrs", "rev", "rev_pct", "spec", "spec_pct", "tier_rev", "tier_spec"]
    rows = []
    for g, gd in par.groupby("group_id"):
        regs, A, A_adj, linked = load(net_of(grp.loc[g, "cohort"]))
        D, Da, pe = activity(gd.pert_id, gd.source, keys.loc[g], w[w.group_id == g].set_index("ModelID").weight,
                             A, A_adj)
        m = gmr[(gmr.group_id == g) & gmr.regulator.isin(regs)]
        sc, mn = split_scores(D, Da, regs, linked, m)
        o = gd[keep].copy()
        for k, v in sc.items():
            o[k] = np.round(v, 3) if isinstance(v, np.ndarray) else v
        o["prolif_effect"] = pe.round(3)
        ref = refs_mech(Da, regs, linked, keys.loc[g, "key"])
        o["n_ref_mech"] = len(ref)
        o["spec_mech"] = (rz(sc["rev_mech"]) - np.mean(ref, axis=0)).round(3) \
            if len(ref) >= MIN_REF and len(mn) >= 2 else np.nan
        rows.append(finish(o, []))
        print(f"{g:<20} MRs {len(m):>4} (prolif {sc['n_mrs_prolif']:>3})  "
              f"rho(rev, rev_mech) {o.rev.corr(o.rev_mech, method='spearman'):.2f}", flush=True)
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(RES / "drug_reversal_split.tsv", sep="\t", index=False)
    pd.concat(p2_rows).to_csv(RES / "prolif_regulators_p2.tsv", sep="\t", index=False)

    # strata
    st = pd.read_csv(RES / "stratum_reversal.tsv", sep="\t", low_memory=False)
    sw = pd.read_csv(ROOT / "results/celligner/stratum_line_weights_long.tsv", sep="\t")
    srows = []
    for s, gd in st.groupby("stratum"):
        p = gd.parent_group.iloc[0]
        regs, A, A_adj, linked = load(net_of(grp.loc[p, "cohort"]))
        D, Da, pe = activity(gd.pert_id, gd.source, keys.loc[p], sw[sw.stratum == s].set_index("ModelID").weight,
                             A, A_adj)
        m = smr[(smr.group_id == s) & smr.regulator.isin(regs)]
        sc, mn = split_scores(D, Da, regs, linked, m)
        o = gd[["stratum", "parent_group", "weak", "tentative", "possible_technical", "pert_id", "cmap_name", "source",
                "passes_gate", "approved", "multi_target", "hub_moa", "n_mrs", "rev", "rev_pct", "spec", "spec_pct",
                "spec_sib", "spec_sib_pct", "tier_rev", "tier_spec", "tier_sib", "siblings"]].copy()
        for k, v in sc.items():
            o[k] = np.round(v, 3) if isinstance(v, np.ndarray) else v
        o["prolif_effect"] = pe.round(3)
        ref = refs_mech(Da, regs, linked, keys.loc[p, "key"])
        z = rz(sc["rev_mech"])
        o["spec_mech"] = (z - np.mean(ref, axis=0)).round(3) if len(ref) >= MIN_REF and len(mn) >= 2 else np.nan
        sibz = []
        for t in gd.siblings.iloc[0].split(";"):
            mt = smr[(smr.group_id == t) & smr.regulator.isin(regs) & ~smr.regulator.isin(linked)]
            if len(mt) >= 2:
                sibz.append(rz(mr_rev(Da, regs, mt)))
        o["spec_sib_mech"] = (z - np.mean(sibz, axis=0)).round(3) if sibz and len(mn) >= 2 else np.nan
        srows.append(finish(o, ["spec_sib_mech"]))
        print(f"{s:<22} MRs {len(m):>4} (prolif {sc['n_mrs_prolif']:>3})", flush=True)
    sout = pd.concat(srows, ignore_index=True)
    sout.to_csv(RES / "stratum_reversal_split.tsv", sep="\t", index=False)

    # controls: best LINCS ID per drug, percentile under every score
    ctl = []
    for g, exp, names in CONTROLS:
        x = out[out.group_id == g]
        for nm in names:
            r = x[x.cmap_name == nm]
            if len(r):
                ctl.append(dict(unit=g, expected=exp, drug=nm, **{c: r[c + "_pct"].max() for c in SCORES},
                                prolif_share=r.prolif_share.max()))
    for s, exp, names in STRATUM_CONTROLS:
        x = sout[sout.stratum == s]
        for nm in names:
            r = x[x.cmap_name == nm]
            if len(r):
                ctl.append(dict(unit=s, expected=exp, drug=nm, **{c: r[c + "_pct"].max() for c in SCORES},
                                spec_sib=r.spec_sib_pct.max(), spec_sib_mech=r.spec_sib_mech_pct.max(),
                                prolif_share=r.prolif_share.max()))
    ctl = pd.DataFrame(ctl)
    ctl.to_csv(RES / "split_controls.tsv", sep="\t", index=False)

    # overview
    ok = out[~out.low_mrs]
    ov = [("groups", out.group_id.nunique()),
          ("prolif_linked_mrs_frac_median", round((out.drop_duplicates("group_id").n_mrs_prolif
                                                  / out.drop_duplicates("group_id").n_mrs).median(), 3))]
    for c in ["rev", "rev_prolif", "rev_nonprolif", "rev_mech", "prolif_effect", "spec", "spec_mech"]:
        top = ok[ok[c + "_pct"] >= 100 - TOP_PCT].groupby("pert_id").size()
        ov.append((f"drugs_top10_{c}_in_ge20_groups", int((top >= 20).sum())))
    for c in ["rev_nonprolif", "rev_mech", "rev_prolif", "prolif_effect", "spec_mech"]:
        ov.append((f"median_rho_rev_vs_{c}", round(ok.groupby("group_id").apply(
            lambda d: d.rev.corr(d[c], method="spearman")).median(), 3)))
    for exp in ["high", "low"]:
        d = ctl[ctl.expected == exp]
        for c in SCORES:
            ov.append((f"controls_{exp}_{c}_pct_mean", round(d[c].mean(), 1)))
    ov.append(("prolif_driven_rows_frac_of_top10_rev",
               round(ok.prolif_driven.sum() / (ok.rev_pct >= 100 - TOP_PCT).sum(), 3)))
    ov = pd.DataFrame(ov, columns=["metric", "value"])
    ov.to_csv(RES / "split_overview.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(ov.to_string(index=False))
    print(ctl.to_string(index=False))


if __name__ == "__main__":
    main()
