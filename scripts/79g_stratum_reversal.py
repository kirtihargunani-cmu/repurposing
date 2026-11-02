#!/usr/bin/env python
"""Drug reversal and specificity for the Step 10 strata (runs in `celligner`).

Instructions: main:docs/notes/subgrouping_for_drug_side.md (sections 4-6); user decisions 2026-10-12.
Inputs:
  stratum MRs            results/step10/stratum_mr.tsv (disease side, scripts/95: msviper on the stratum signature,
                         parent's network and null; same columns as group_mr.tsv + parent_group)
  stratum drug lists     results/lincs/stratum_drugs.tsv (scripts/79b STRATA=1: parent's tissue key and direct /
                         borrowed split, target check judged on the stratum's own tumors)
  stratum line weights   results/celligner/stratum_line_weights_long.tsv (scripts/74b: mean of its own patients)
  drug activity          drug x line regulator activity of the parent's network (scripts/79d, PROFILE=specific)
  parent results         results/lincs/drug_reversal.tsv (41-group run, scripts/79e)
Scores (same definitions as scripts/79e; drug activity per stratum = its own line weights):
  rev, rev_raw, r_rev, rev_pct, rev_z   reversal of the stratum's MRs (FDR < MR_FDR)
  spec, spec_pct                        vs the 41-group diseases of other tissue keys (>= MIN_MRS_REF MRs)
  spec_sib, spec_sib_pct                rev_z minus the mean rev_z of the same drug activity scored against the
                                        sibling strata's MRs (COAD CMS1/3 vs CMS2/4; each READ cluster vs the other
                                        two; high vs low end of the same breast axis)
  tier_rev / tier_spec / tier_sib       1 approved + passes gate + (rev > 0 and top TOP_PCT% by rev / spec / spec_sib),
                                        2 same but not approved, 3 rest
  parent_rev_pct, parent_tier_rev       the drug in the parent group's 41-group run
Flags: weak (Her2 PC4, disease side), tentative (READ clusters), possible_technical (low ends of Basal PC3 and Her2 PC4:
housekeeping / ribosomal genes; rank last).
Outputs: results/lincs/stratum_reversal.tsv (stratum x drug), stratum_reversal_summary.tsv (per stratum, incl. overlap
         with the parent and the siblings), stratum_reversal_controls.tsv (expected sibling contrasts)
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed/lincs"
RES = ROOT / "results/lincs"
MR_FDR, TOP_PCT, MULTI_TARGET, MIN_MRS_REF, MIN_REF = 0.10, 10, 5, 20, 10
POSSIBLE_TECHNICAL = {"BRCA-Basal-PC3-low", "BRCA-Her2-PC4-low"}
# expected sibling contrasts (sanity check): drug classes that should be more specific to one end than the other
CONTROLS = [("BRCA-LumA-PC3-low", "ER end of the proliferation axis",
             ["fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene"]),
            ("BRCA-Her2-PC3-low", "ESR1 / XBP1 end",
             ["fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene"]),
            ("BRCA-LumA-PC3-high", "proliferation end", ["palbociclib", "abemaciclib", "ribociclib", "paclitaxel",
                                                         "docetaxel", "vinorelbine"])]


def net_of(cohort):
    return "COADREAD" if cohort in ("TCGA-COAD", "TCGA-READ") else cohort.split("-", 1)[1]


def rz(v):
    return (v - np.median(v)) / (1.4826 * np.median(np.abs(v - np.median(v))))


def pct(v):
    return (pd.Series(v).rank(method="max", pct=True).values * 100).round(1)


def main():
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    idx["col"] = np.arange(len(idx))
    by_drug = idx.groupby("pert_id")
    grp = pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t").set_index("group_id")
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
    models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", index_col=0)
    models = models[models.OncotreePrimaryDisease != "Non-Cancerous"]
    smr = pd.read_csv(ROOT / "results/step10/stratum_mr.tsv", sep="\t")
    gmr = pd.read_csv(ROOT / "results/step6/aracne_project/group_mr.tsv", sep="\t")
    sd = pd.read_csv(RES / "stratum_drugs.tsv", sep="\t", low_memory=False)
    sw = pd.read_csv(ROOT / "results/celligner/stratum_line_weights_long.tsv", sep="\t")
    par = pd.read_csv(RES / "drug_reversal.tsv", sep="\t", low_memory=False)
    par = par.drop_duplicates(["group_id", "pert_id"]).set_index(["group_id", "pert_id"])

    sib_set = sw.drop_duplicates("stratum").set_index("stratum").sibling_set
    ref_ok = gmr[gmr.fdr < MR_FDR].groupby("group_id").size()
    ref_ok = [g for g in ref_ok[ref_ok >= MIN_MRS_REF].index if g in keys.index]

    nets = {}
    def read(kind, net):
        f = PROC / f"drug_viper_{kind}" / f"{net}.f32"
        if not f.exists():
            return None, None
        regs = pd.Index((PROC / f"drug_viper_{kind}" / f"{net}.regulators.txt").read_text().split())
        return regs, np.fromfile(f, dtype=np.float32).reshape(len(idx), len(regs)).astype(np.float64)

    def load(net):
        if net not in nets:
            regs, A = read("specific", net)
            rr, R = read("raw", net)
            nets[net] = (regs, A, R.T[rr.get_indexer(regs)].T if R is not None else None)
        return nets[net]

    def mr_score(D, regs, m):
        m = m[(m.fdr < MR_FDR) & m.regulator.isin(regs)]
        c = regs.get_indexer(m.regulator)
        return -(D[:, c] * np.sign(m.nes.values)).sum(1) / np.sqrt(len(c)), m

    rows, summ = [], []
    for s, gd in sd.groupby("group_id"):
        p = gd.parent_group.iloc[0]
        net = net_of(grp.loc[p, "cohort"])
        regs, A, A_raw = load(net)
        k = keys.loc[p]
        km = models[models.OncotreeLineage.isin(str(k.lineages).split(";"))]
        if isinstance(k.diseases, str) and k.diseases:
            km = km[km.OncotreePrimaryDisease.isin(k.diseases.split(";"))]
        wg = sw[sw.stratum == s].set_index("ModelID").weight

        D = np.empty((len(gd), len(regs)))
        D_raw = np.full_like(D, np.nan)
        n_lines = np.zeros(len(gd), int)
        for i, (d, src) in enumerate(zip(gd.pert_id, gd.source)):
            q = by_drug.get_group(d)
            q = q[q.ModelID.isin(km.index)] if src == "direct" else q[~q.ModelID.isin(km.index)]
            ww = wg.reindex(q.ModelID).fillna(0).values if src == "direct" else np.zeros(len(q))
            ww = ww if ww.sum() > 0 else np.ones(len(q))
            D[i] = (ww[:, None] * A[q.col.values]).sum(0) / ww.sum()
            if A_raw is not None:
                D_raw[i] = (ww[:, None] * A_raw[q.col.values]).sum(0) / ww.sum()
            n_lines[i] = len(q)

        rev, mrs = mr_score(D, regs, smr[smr.group_id == s])
        n = len(mrs)
        dn = mrs.nes.values - mrs.nes.values.mean()
        Dc = D[:, regs.get_indexer(mrs.regulator)]
        Dc = Dc - Dc.mean(1, keepdims=True)
        r_rev = -(Dc @ dn) / (np.sqrt((Dc ** 2).sum(1)) * np.linalg.norm(dn) + 1e-12)
        rev_z = rz(rev)
        ref = [rz(mr_score(D, regs, gmr[gmr.group_id == g])[0]) for g in ref_ok if keys.loc[g, "key"] != k.key]
        sibs = [t for t in sib_set.index if sib_set[t] == sib_set[s] and t != s]
        sibz = [rz(mr_score(D, regs, smr[smr.group_id == t])[0]) for t in sibs]

        o = gd[["group_id", "parent_group", "weak", "pert_id", "cmap_name", "source", "passes_gate", "gate_reason",
                "n_targets", "hub_moa", "hub_clinical_phase", "chembl_max_phase", "pharmacology_flag"]].copy()
        o = o.rename(columns={"group_id": "stratum"})
        o["n_lines_used"] = n_lines
        o["n_mrs"] = n
        o["rev"] = rev.round(3)
        o["rev_raw"] = (-(D_raw[:, regs.get_indexer(mrs.regulator)] * np.sign(mrs.nes.values)).sum(1)
                        / np.sqrt(n)).round(3)
        o["r_rev"] = r_rev.round(3)
        o["rev_pct"] = pct(rev)
        o["rev_z"] = rev_z.round(3)
        o["spec"] = (rev_z - np.mean(ref, axis=0)).round(3) if len(ref) >= MIN_REF else np.nan
        o["spec_pct"] = pct(o.spec.values)
        o["n_ref_diseases"] = len(ref)
        o["spec_sib"] = (rev_z - np.mean(sibz, axis=0)).round(3)
        o["spec_sib_pct"] = pct(o.spec_sib.values)
        o["siblings"] = ";".join(sibs)
        o["approved"] = (o.chembl_max_phase == 4) | (o.hub_clinical_phase == "Launched")
        o["multi_target"] = o.n_targets > MULTI_TARGET
        for t, col in [("tier_rev", "rev_pct"), ("tier_spec", "spec_pct"), ("tier_sib", "spec_sib_pct")]:
            hit = (o.rev > 0) & (o[col] >= 100 - TOP_PCT)
            o[t] = np.where(hit & o.passes_gate & o.approved, 1, np.where(hit & o.passes_gate, 2, 3))
        pk = [(p, d) for d in o.pert_id]
        o["parent_rev_pct"] = par.rev_pct.reindex(pk).values
        o["parent_tier_rev"] = par.tier_rev.reindex(pk).values
        o["tentative"] = s.startswith("READ-")
        o["possible_technical"] = s in POSSIBLE_TECHNICAL
        o = o.sort_values(["tier_sib", "spec_sib"], ascending=[True, False])
        o["rank_sib"] = np.arange(1, len(o) + 1)
        rows.append(o)

        top = set(o.pert_id[(o.rev > 0) & (o.rev_pct >= 100 - TOP_PCT)])
        pr = par.loc[p]
        ptop = set(pr.index[(pr.rev > 0) & (pr.rev_pct >= 100 - TOP_PCT)])
        both = o.set_index("pert_id").rev.to_frame().join(pr.rev.rename("prev"), how="inner")
        t1 = o[o.tier_sib == 1]
        summ.append(dict(stratum=s, parent=p, weak=bool(o.weak.iloc[0]), tentative=bool(o.tentative.iloc[0]),
                         possible_technical=s in POSSIBLE_TECHNICAL, n_mrs=n, n_drugs=len(o),
                         n_pass_gate=int(o.passes_gate.sum()), siblings=";".join(sibs),
                         rho_rev_vs_parent=round(both.rev.corr(both.prev, method="spearman"), 3),
                         jaccard_top_rev_vs_parent=round(len(top & ptop) / len(top | ptop), 3),
                         n_top_rev_not_in_parent_top=len(top - ptop),
                         rho_rev_vs_spec_sib=round(o.rev.corr(o.spec_sib, method="spearman"), 3),
                         n_tier1_rev=int((o.tier_rev == 1).sum()), n_tier1_spec=int((o.tier_spec == 1).sum()),
                         n_tier1_sib=len(t1), n_tier1_rev_and_sib=int(((o.tier_rev == 1) & (o.tier_sib == 1)).sum()),
                         top_tier1_sib=", ".join(f"{a} {b:.1f}" for a, b in zip(t1.cmap_name.head(8),
                                                                                  t1.spec_sib.head(8)))))
        print(f"{s:<20} MRs {n:>4}  vs parent rho {summ[-1]['rho_rev_vs_parent']:.2f}  "
              f"top10 J {summ[-1]['jaccard_top_rev_vs_parent']:.2f}  tier1 sib {len(t1):>3}  "
              f"{summ[-1]['top_tier1_sib'][:80]}", flush=True)

    out = pd.concat(rows, ignore_index=True)
    out.to_csv(RES / "stratum_reversal.tsv", sep="\t", index=False)
    summ = pd.DataFrame(summ)
    # sibling similarity of rev rankings (compare: COAD~READ rho 0.96, top-10% Jaccard 0.71 in the 41-group run)
    R = out.pivot_table(index="pert_id", columns="stratum", values="rev")
    tops = {s: set(d.pert_id[(d.rev > 0) & (d.rev_pct >= 100 - TOP_PCT)]) for s, d in out.groupby("stratum")}
    tsib = {s: set(d.pert_id[d.tier_sib <= 2]) for s, d in out.groupby("stratum")}
    def sib_stats(r):
        sibs = r.siblings.split(";")
        return pd.Series(dict(
            rho_rev_vs_siblings=round(np.mean([R[r.stratum].corr(R[t], method="spearman") for t in sibs]), 3),
            jaccard_top_rev_vs_siblings=round(np.mean([len(tops[r.stratum] & tops[t]) / len(tops[r.stratum] | tops[t])
                                                       for t in sibs]), 3),
            jaccard_tier12_sib_vs_siblings=round(np.mean([len(tsib[r.stratum] & tsib[t]) / max(1, len(tsib[r.stratum]
                                                          | tsib[t])) for t in sibs]), 3)))
    summ = pd.concat([summ, summ.apply(sib_stats, axis=1)], axis=1)
    summ.to_csv(RES / "stratum_reversal_summary.tsv", sep="\t", index=False)

    ctl = []
    for s, why, names in CONTROLS:
        x = out[out.stratum == s]
        sibs = x.siblings.iloc[0].split(";")
        for nm in names:
            r = x[x.cmap_name == nm]
            if not len(r):
                ctl.append(dict(stratum=s, expected=why, drug=nm, source="not in list"))
                continue
            rs = out[out.stratum.isin(sibs) & (out.cmap_name == nm)]
            ctl.append(dict(stratum=s, expected=why, drug=nm, n_ids=len(r), rev_pct=r.rev_pct.max(),
                            sibling_rev_pct=rs.rev_pct.max(), spec_sib=r.spec_sib.max(),
                            spec_sib_pct=r.spec_sib_pct.max(), tier_sib=r.tier_sib.min()))
    ctl = pd.DataFrame(ctl)
    ctl.to_csv(RES / "stratum_reversal_controls.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(ctl.to_string(index=False))


if __name__ == "__main__":
    main()
