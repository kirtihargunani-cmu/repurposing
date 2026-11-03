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
VERSION=v2 [REGSET=tf]: v2 MRs, activity, P1 list and 79e scores (suffix _v2[_tf] on inputs and outputs); the strata
part is skipped (no v2 Step 10 strata yet).
SCORE=weighted (default for v2, D74): every regulator tested in the disease, weight = disease NES (no MR_FDR cutoff), in
  every score (rev_prolif / rev_nonprolif / rev_mech: the weighted sum over the proliferation-linked / other regulators).
REGULON=consensus (VERSION=v2, D75): ARACNe v2 TF (disease network) and CollecTRI side by side; per set: P1 from its own
  79h list (prolif_regulators_v2_tf / _collectri), P2 and the proliferation adjustment on its own activity; then, per
  disease, the per-regulator mean over the sets scored on both sides (as 79e). A regulator is proliferation-linked if it
  is linked in any of its sets (fixed 2026-10-28). Inputs: 79e drug_reversal_v2_consensus.tsv; outputs suffix _v2_consensus.
STRATA=1 (v2; user 2026-10-29): also the D79 subgroups, same rules: MRs from the <set>_strata directories, line
  weights stratum_line_weights_long_v2, 79e STRATA=1 output stratum_reversal<OUT>.tsv; spec_mech vs the 41-group
  diseases of other tissue keys, spec_sib_mech vs the sibling subgroup (proliferation-linked regulators removed, adjusted
  activity). Output stratum_reversal_split<OUT>.tsv.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed/lincs"
RES = ROOT / "results/lincs"
MR_FDR, TOP_PCT, MIN_MRS_REF, MIN_REF, R_MIN = 0.10, 10, 20, 10, 0.5
SFX = ("_v2" if os.environ.get("VERSION", "v1") == "v2" else "") + ("_tf" if os.environ.get("REGSET", "full") == "tf" else "")
assert SFX != "_tf", "REGSET=tf needs VERSION=v2"
STRATA_V2 = os.environ.get("STRATA", "0") == "1"
STRATA = SFX == "" or STRATA_V2
assert not STRATA_V2 or SFX.startswith("_v2"), "STRATA=1 needs VERSION=v2"
REGULON = os.environ.get("REGULON", "aracne_project"); CONS = REGULON == "consensus"
assert not CONS or SFX == "_v2", "consensus: VERSION=v2"
SCORE = os.environ.get("SCORE", "weighted" if SFX.startswith("_v2") else "fdr"); assert SCORE in ("weighted", "fdr")
OUT = ("_v2_consensus" if CONS else SFX)
if SCORE == "fdr" and OUT.startswith("_v2"):
    OUT += "_fdr"
# ARACNE_SET=tf (default) | tf_base | tf_sig1 | tf_sig2 (D83): which ARACNe TF-only target lists the consensus uses;
#   other than tf, outputs get the suffix _<ARACNE_SET> (e.g. drug_reversal_v2_consensus_tf_sig1.tsv).
ARACNE_SET = os.environ.get("ARACNE_SET", "tf"); assert ARACNE_SET in ("tf", "tf_base", "tf_sig1", "tf_sig2")
XS = "" if ARACNE_SET == "tf" else f"_{ARACNE_SET}"
if CONS: OUT += XS
# PROFILE=specific (default) | raw and RUN_TAG (D85, 2026-11-02): activity from drug_viper_<PROFILE><set><RUN_TAG>, the
#   proliferation adjustment and P2 on the same profiles' proliferation score (prolif_effect stays on raw profiles);
#   79e input and outputs with RUN_TAG plus "_raw" for raw
PROFILE = os.environ.get("PROFILE", "specific"); assert PROFILE in ("specific", "raw")
RUN_TAG = os.environ.get("RUN_TAG", "")
OUT += RUN_TAG + ("_raw" if PROFILE == "raw" else "")
# IDENTITY_FILE (D91, 2026-11-03): table net / regulator / identity (scripts/d91/d91_identity.R). Adds, for the 41
#   groups: rev_identity (identity-linked, non-proliferation regulators, activity A), rev_mech_id (rev_mech without
#   identity-linked regulators, except the drug's own expressed targets from 79b) and spec_mech_id (as spec_mech);
#   outputs get "_id" appended to the suffix. Default: off, outputs unchanged.
IDENT_FILE = os.environ.get("IDENTITY_FILE", "")
IN79E = OUT                                    # 79e output suffix (input here)
if IDENT_FILE:
    OUT += "_id"
CONS_SETS = {"aracne": (f"aracne_project_v2_{ARACNE_SET}", f"_v2_{ARACNE_SET}"), "collectri": ("collectri_v2", "_collectri")}
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


# NET_MAP / MR_OVERRIDE (D92, 2026-11-03): groups mapped to another network, their group_mr rows from another Step 6
#   folder (tumor ARACNe only, not consensus). Defaults: off.
NETMAP = (pd.read_csv(os.environ["NET_MAP"], sep="\t").set_index("group_id").net.to_dict() if os.environ.get("NET_MAP") else {})
MR_OVR = os.environ.get("MR_OVERRIDE", "")


def mr_override(mr):
    if not MR_OVR:
        return mr
    o = pd.read_csv(ROOT / MR_OVR / "group_mr.tsv", sep="\t")
    o = o[o.group_id.isin(NETMAP)]
    return pd.concat([mr[~mr.group_id.isin(o.group_id)], o], ignore_index=True)


def rz(v):
    return (v - np.median(v)) / (1.4826 * np.median(np.abs(v - np.median(v))))


def pct(v):
    return (pd.Series(v).rank(method="max", pct=True).values * 100).round(1)


def mr_rev(D, regs, m):
    c = regs.get_indexer(m.regulator)
    w = m.nes.values if SCORE == "weighted" else np.sign(m.nes.values)
    return -(D[:, c] * w).sum(1) / max(np.sqrt((w ** 2).sum()), 1e-12)


def wts(m):
    return m.nes.values if SCORE == "weighted" else np.sign(m.nes.values)


def rev_extra(D, regs, base, pool, extras):
    """rev over base regulators, plus per drug the pool regulators among its extras (D91 exception)."""
    w = wts(base); num = -(D[:, regs.get_indexer(base.regulator)] * w).sum(1); den = np.full(D.shape[0], (w ** 2).sum())
    if len(pool):
        for i, ex in enumerate(extras):
            if ex:
                q = pool[pool.regulator.isin(ex)]
                if len(q):
                    w2 = wts(q); num[i] -= D[i, regs.get_indexer(q.regulator)] @ w2; den[i] += (w2 ** 2).sum()
    return num / np.sqrt(np.maximum(den, 1e-12))


def consensus(Ds, regs, m, linked):
    """Set-tagged ("aracne|X") activity matrices / disease NES -> per-regulator means over the sets in m."""
    base = m.regulator.str.split("|", n=1).str[1].values
    cols = regs.get_indexer(m.regulator)
    ub, inv = np.unique(base, return_inverse=True); cnt = np.bincount(inv, minlength=len(ub))
    out = []
    for D in Ds:
        Dc = np.zeros((D.shape[0], len(ub))); np.add.at(Dc.T, inv, D[:, cols].T); out.append(Dc / cnt)
    mc = pd.DataFrame({"group_id": m.group_id.iloc[0] if len(m) else "", "regulator": ub,
                       "nes": np.bincount(inv, weights=m.nes.values, minlength=len(ub)) / cnt,
                       "fdr": m.groupby(base).fdr.min().reindex(ub).values})
    lk = set(pd.Series(base)[m.regulator.isin(linked).values])
    return out, pd.Index(ub), mc, lk


def main():
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    idx["col"] = np.arange(len(idx))
    by_drug = idx.groupby("pert_id")
    genes = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
    gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:]
           for l in open(ROOT / "data/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt")}
    pm = genes.symbol.isin(set(gmt["HALLMARK_E2F_TARGETS"]) | set(gmt["HALLMARK_G2M_CHECKPOINT"])).values
    # PROFILE_PREFIX (D90, 2026-11-03): <prefix>_raw.f32 / <prefix>_specific.f32 with the 79c layout replace the 79c profiles
    pp = os.environ.get("PROFILE_PREFIX", "")
    f_raw, f_spec = (Path(pp + "_raw.f32"), Path(pp + "_specific.f32")) if pp else (PROC / "drug_line_profiles_bing.f32", PROC / "drug_line_profiles_bing_specific.f32")
    pi_raw = np.fromfile(f_raw, dtype=np.float32).reshape(len(idx), -1)[:, pm].mean(1)
    pi_spec = np.fromfile(f_spec,
                          dtype=np.float32).reshape(len(idx), -1)[:, pm].mean(1).astype(np.float64)
    print(f"proliferation genes on profiles: {pm.sum()}")

    grp = pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t").set_index("group_id")
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
    models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", index_col=0)
    models = models[models.OncotreePrimaryDisease != "Non-Cancerous"]
    if CONS:
        gmr = pd.concat([pd.read_csv(ROOT / f"results/step6/{d}/group_mr.tsv", sep="\t").assign(regulator=lambda x, k=k: k + "|" + x.regulator)
                         for k, (d, _) in CONS_SETS.items()], ignore_index=True)
    else:
        gmr = mr_override(pd.read_csv(ROOT / f"results/step6/aracne_project{SFX}/group_mr.tsv", sep="\t"))
    if STRATA_V2:
        sets = [(k, d) for k, (d, _) in CONS_SETS.items()] if CONS else [("", f"aracne_project{SFX}")]
        smr = pd.concat([pd.read_csv(ROOT / f"results/step6/{d}_strata/group_mr.tsv", sep="\t")
                         .assign(regulator=lambda x, k=k: (k + "|" if k else "") + x.regulator) for k, d in sets], ignore_index=True)
    else:
        smr = pd.read_csv(ROOT / "results/step10/stratum_mr.tsv", sep="\t") if STRATA else None
    def p1_of(tag):
        x = pd.read_csv(RES / f"prolif_regulators{tag}.tsv", sep="\t")
        if os.environ.get("P1_EXTRA"):   # D92: P1 lists of extra networks (79h EXTRA_RDS run)
            x = pd.concat([x, pd.read_csv(os.environ["P1_EXTRA"], sep="\t")], ignore_index=True)
        return x[x.prolif_linked].groupby("net").regulator.apply(set)
    p1 = {k: p1_of(t) for k, (_, t) in CONS_SETS.items()} if CONS else p1_of(SFX)
    if SCORE == "fdr":
        gmr = gmr[gmr.fdr < MR_FDR]
    smr = smr[smr.fdr < MR_FDR] if STRATA and SCORE == "fdr" else smr

    nets, p2_rows = {}, []
    def load(net):
        if net not in nets and CONS:
            parts = [load_one(net, f"_v2_{ARACNE_SET}", p1["aracne"], "aracne|"), load_one("ALL", "_collectri", p1["collectri"], "collectri|")]
            nets[net] = (pd.Index(np.concatenate([x[0] for x in parts])), np.hstack([x[1] for x in parts]),
                         np.hstack([x[2] for x in parts]), set().union(*(x[3] for x in parts)))
        elif net not in nets:
            nets[net] = load_one(net, SFX, p1, "")
        return nets[net]

    def load_one(net, tag, p1s, prefix):
        key = (net, tag)
        if key not in nets:
            vdir = PROC / f"drug_viper_{PROFILE}{tag}{RUN_TAG}"
            regs = pd.Index((vdir / f"{net}.regulators.txt").read_text().split())
            A = np.fromfile(vdir / f"{net}.f32",
                            dtype=np.float32).reshape(len(idx), len(regs)).astype(np.float64)
            pa = pi_raw.astype(np.float64) if PROFILE == "raw" else pi_spec
            pc = pa - pa.mean()
            Ac = A - A.mean(0)
            beta = Ac.T @ pc / (pc @ pc)
            r = beta * pc.std() / Ac.std(0)
            A_adj = A - np.outer(pa, beta)
            linked = pd.Index(regs).isin(p1s.get(net, set())) | (np.abs(r) >= R_MIN)
            p2_rows.append(pd.DataFrame({"set": prefix.rstrip("|") or "aracne", "net": net, "regulator": regs, "r_prolif": r.round(3),
                                         "p1": pd.Index(regs).isin(p1s.get(net, set())),
                                         "p2": np.abs(r) >= R_MIN, "prolif_linked": linked}))
            regs = pd.Index(prefix + regs)
            nets[key] = (regs, A, A_adj, set(regs[linked]))
        return nets[key]

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
            if CONS:   # Da / regs / linked are the set-tagged ones here; average over g2's own sets
                (Da2,), regs2, m2, lk2 = consensus([Da], regs, m2[m2.regulator.isin(regs)], linked)
                m2 = m2[~m2.regulator.isin(lk2)]
                if len(m2) >= MIN_MRS_REF:
                    out.append(rz(mr_rev(Da2, regs2, m2)))
                continue
            m2 = m2[m2.regulator.isin(regs) & ~m2.regulator.isin(linked)]
            if len(m2) >= MIN_MRS_REF:
                out.append(rz(mr_rev(Da, regs, m2)))
        return out

    ident, dtargets = {}, {}
    if IDENT_FILE:
        it = pd.read_csv(IDENT_FILE, sep="\t")
        ident = it[it.identity].groupby("net").regulator.apply(set).to_dict()
        dd = pd.read_csv(RES / "disease_drugs.tsv", sep="\t", usecols=["group_id", "pert_id", "targets_expressed"], low_memory=False)
        dtargets = {(g, p): set(str(t).split("|")) for g, p, t in dd.itertuples(index=False) if isinstance(t, str)}

    def refs_mech_id(Da, regs, linked, idset, extras, own_key):
        out = []
        for g2, m2 in gmr.groupby("group_id"):
            if g2 not in keys.index or keys.loc[g2, "key"] == own_key:
                continue
            if CONS:
                (Da2,), regs2, m2, lk2 = consensus([Da], regs, m2[m2.regulator.isin(regs)], linked)
                m2 = m2[~m2.regulator.isin(lk2)]; D2, r2 = Da2, regs2
            else:
                m2 = m2[m2.regulator.isin(regs) & ~m2.regulator.isin(linked)]; D2, r2 = Da, regs
            base, pool = m2[~m2.regulator.isin(idset)], m2[m2.regulator.isin(idset)]
            if len(base) >= MIN_MRS_REF:
                out.append(rz(rev_extra(D2, r2, base, pool, extras)))
        return out

    def finish(o, cols):
        for c in ["rev_prolif", "rev_nonprolif", "rev_mech", "prolif_effect", "spec_mech"] + cols:
            o[c + "_pct"] = pct(o[c].values)
        for t, s, c in [("tier_mech", "rev_mech", "rev_mech_pct"), ("tier_spec_mech", "rev_mech", "spec_mech_pct")] + \
                       [(f"tier_{x}", "rev_mech", f"{x}_pct") for x in cols]:
            hit = (o[s] > 0) & (o[c] >= 100 - TOP_PCT)
            ok_t = o.passes_gate & ~o.unmatched_lines & ~o.low_transfer if "unmatched_lines" in o else o.passes_gate   # D76, D77
            o[t] = np.where(hit & ok_t & o.approved, 1, np.where(hit & ok_t, 2, 3))
        o["prolif_driven"] = (o.prolif_share > 0.5) & (o.rev_pct >= 100 - TOP_PCT) & (o.rev_mech_pct < 100 - TOP_PCT)
        return o

    # 41 groups
    par = pd.read_csv(RES / f"drug_reversal{IN79E}.tsv", sep="\t", low_memory=False)
    w = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    keep = ["group_id", "pert_id", "cmap_name", "source", "passes_gate", "approved", "multi_target", "hub_moa",
            *(["line_match", "unmatched_lines", "transfer_conf", "low_transfer"] if CONS or SFX.startswith("_v2") else []),
            "n_mrs", "low_mrs", "rev", "rev_pct", "spec", "spec_pct", "tier_rev", "tier_spec"]
    rows = []
    for g, gd in par.groupby("group_id"):
        regs, A, A_adj, linked = load(NETMAP.get(g, net_of(grp.loc[g, "cohort"])))
        D, Da, pe = activity(gd.pert_id, gd.source, keys.loc[g], w[w.group_id == g].set_index("ModelID").weight,
                             A, A_adj)
        m = gmr[(gmr.group_id == g) & gmr.regulator.isin(regs)]
        regs0, Da0, linked0 = regs, Da, linked
        if CONS:
            (D, Da), regs, m, linked = consensus([D, Da], regs, m, linked)
        sc, mn = split_scores(D, Da, regs, linked, m)
        o = gd[keep].copy()
        for k, v in sc.items():
            o[k] = np.round(v, 3) if isinstance(v, np.ndarray) else v
        o["prolif_effect"] = pe.round(3)
        ref = refs_mech(Da0, regs0, linked0, keys.loc[g, "key"])
        o["n_ref_mech"] = len(ref)
        o["spec_mech"] = (rz(sc["rev_mech"]) - np.mean(ref, axis=0)).round(3) \
            if len(ref) >= MIN_REF and len(mn) >= 2 else np.nan
        cols = []
        if IDENT_FILE:
            idset = ident.get(NETMAP.get(g, net_of(grp.loc[g, "cohort"])), set())
            extras = [dtargets.get((g, p), set()) for p in gd.pert_id]
            base, pool = mn[~mn.regulator.isin(idset)], mn[mn.regulator.isin(idset)]
            o["rev_mech_id"] = rev_extra(Da, regs, base, pool, extras).round(3)
            o["rev_identity"] = np.round(mr_rev(D, regs, pool), 3) if len(pool) else 0.0
            o["n_mrs_identity"] = len(pool)
            ref2 = refs_mech_id(Da0, regs0, linked0, idset, extras, keys.loc[g, "key"])
            o["spec_mech_id"] = (rz(o.rev_mech_id.values) - np.mean(ref2, axis=0)).round(3) \
                if len(ref2) >= MIN_REF and len(base) >= 2 else np.nan
            cols = ["rev_mech_id", "rev_identity", "spec_mech_id"]
        rows.append(finish(o, cols))
        print(f"{g:<20} MRs {len(m):>4} (prolif {sc['n_mrs_prolif']:>3})  "
              f"rho(rev, rev_mech) {o.rev.corr(o.rev_mech, method='spearman'):.2f}", flush=True)
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(RES / f"drug_reversal_split{OUT}.tsv", sep="\t", index=False)
    pd.concat(p2_rows).to_csv(RES / f"prolif_regulators_p2{OUT}.tsv", sep="\t", index=False)

    # strata
    if STRATA:
        if STRATA_V2:
            st = pd.read_csv(RES / f"stratum_reversal{OUT}.tsv", sep="\t", low_memory=False).rename(columns={"group_id": "stratum"})
            sw = pd.read_csv(ROOT / "results/celligner/stratum_line_weights_long_v2.tsv", sep="\t")
        else:
            st = pd.read_csv(RES / "stratum_reversal.tsv", sep="\t", low_memory=False)
            sw = pd.read_csv(ROOT / "results/celligner/stratum_line_weights_long.tsv", sep="\t")
        srows = []
        for s, gd in st.groupby("stratum"):
            p = gd.parent_group.iloc[0]
            regs, A, A_adj, linked = load(net_of(grp.loc[p, "cohort"]))
            D, Da, pe = activity(gd.pert_id, gd.source, keys.loc[p], sw[sw.stratum == s].set_index("ModelID").weight,
                                 A, A_adj)
            m = smr[(smr.group_id == s) & smr.regulator.isin(regs)]
            regs0, Da0, linked0 = regs, Da, linked
            if CONS:
                (D, Da), regs, m, linked = consensus([D, Da], regs, m, linked)
            sc, mn = split_scores(D, Da, regs, linked, m)
            o = gd[[c for c in ["stratum", "parent_group", "weak", "tentative", "possible_technical", "pert_id", "cmap_name",
                    "source", "passes_gate", "approved", "multi_target", "hub_moa", "line_match", "unmatched_lines",
                    "transfer_conf", "low_transfer", "n_mrs", "rev", "rev_pct", "spec", "spec_pct", "spec_sib",
                    "spec_sib_pct", "tier_rev", "tier_spec", "tier_sib", "siblings"] if c in gd]].copy()
            for k, v in sc.items():
                o[k] = np.round(v, 3) if isinstance(v, np.ndarray) else v
            o["prolif_effect"] = pe.round(3)
            ref = refs_mech(Da0, regs0, linked0, keys.loc[p, "key"])
            z = rz(sc["rev_mech"])
            o["spec_mech"] = (z - np.mean(ref, axis=0)).round(3) if len(ref) >= MIN_REF and len(mn) >= 2 else np.nan
            sibz = []
            for t in gd.siblings.iloc[0].split(";"):
                if CONS:   # sibling's own sets, averaged on this subgroup's drug activity
                    (Da2,), regs2, mt, lk2 = consensus([Da0], regs0, smr[(smr.group_id == t) & smr.regulator.isin(regs0)], linked0)
                    mt = mt[~mt.regulator.isin(lk2)]
                else:
                    Da2, regs2 = Da, regs
                    mt = smr[(smr.group_id == t) & smr.regulator.isin(regs) & ~smr.regulator.isin(linked)]
                if len(mt) >= 2:
                    sibz.append(rz(mr_rev(Da2, regs2, mt)))
            o["spec_sib_mech"] = (z - np.mean(sibz, axis=0)).round(3) if sibz and len(mn) >= 2 else np.nan
            srows.append(finish(o, ["spec_sib_mech"]))
            print(f"{s:<22} MRs {len(m):>4} (prolif {sc['n_mrs_prolif']:>3})", flush=True)
        sout = pd.concat(srows, ignore_index=True)
        sout.to_csv(RES / f"stratum_reversal_split{OUT}.tsv", sep="\t", index=False)

    # controls: best LINCS ID per drug, percentile under every score
    ctl = []
    for g, exp, names in CONTROLS:
        x = out[out.group_id == g]
        for nm in names:
            r = x[x.cmap_name == nm]
            if len(r):
                ctl.append(dict(unit=g, expected=exp, drug=nm, **{c: r[c + "_pct"].max() for c in SCORES},
                                prolif_share=r.prolif_share.max()))
    for s, exp, names in (STRATUM_CONTROLS if STRATA and not STRATA_V2 else []):
        x = sout[sout.stratum == s]
        for nm in names:
            r = x[x.cmap_name == nm]
            if len(r):
                ctl.append(dict(unit=s, expected=exp, drug=nm, **{c: r[c + "_pct"].max() for c in SCORES},
                                spec_sib=r.spec_sib_pct.max(), spec_sib_mech=r.spec_sib_mech_pct.max(),
                                prolif_share=r.prolif_share.max()))
    ctl = pd.DataFrame(ctl)
    ctl.to_csv(RES / f"split_controls{OUT}.tsv", sep="\t", index=False)

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
    ov.to_csv(RES / f"split_overview{OUT}.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(ov.to_string(index=False))
    print(ctl.to_string(index=False))


if __name__ == "__main__":
    main()
