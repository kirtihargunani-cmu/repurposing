#!/usr/bin/env python
"""Rank each disease's drugs by how strongly they reverse the disease's master regulators (runs in `celligner`).

Inputs: drug x line regulator activity (scripts/79d_drug_viper.R, project ARACNe network of the disease's cohort),
disease master regulators (results/step6/aracne_project/group_mr.tsv: msviper NES of tumor vs normal, BH FDR),
per-disease drug lists with the target gate (results/lincs/disease_drugs.tsv, scripts/79b), Celligner disease x line
weights (results/celligner/disease_line_weights_long.tsv).
0. Drug-specific activity: regulator activity of the line-specific profiles (generic drug response removed at the
   gene level, scripts/79c, VIPER in scripts/79d, PROFILE=specific). rev_raw = the same score on raw profiles where
   data/processed/lincs/drug_viper_raw/ has the network (for comparison only).
1. Drug activity per disease (same direct / borrowed split as scripts/79b):
     direct    mean over the drug's lines in the disease's strict tissue key, weighted by the disease's Celligner
               line weights (equal weights if all are 0)
     borrowed  unweighted mean over all the drug's (>= 3) cancer lines (distance weighting adds nothing, scripts/77)
2. Disease master regulators (MRs): FDR < MR_FDR in the disease and scored on the LINCS genes; s_i = sign of NES.
   Reversal score  rev = -sum_i s_i * a_i / sqrt(n)   (a_i = drug activity of MR i, n = number of MRs)
   positive = the drug lowers the activity of MRs that are up in the tumor and raises those that are down.
   Also r_rev = -Pearson(disease NES, drug activity) over the MRs (pattern only, ignores effect size).
3. No p-values: a null of random regulator sets (same size, MR signs) was 3-5x narrower than the spread of rev across
   drugs (sd 1.0-1.3 vs 3.1-4.8 in ACC / BLCA; MRs move together in co-regulated modules, random sets do not), so
   ~38% of drugs reached p < 0.05. Drugs are ranked instead: rev_pct = percent of the disease's drugs with a lower rev.
   null_sd_regset in the summary records that set-permutation spread (N_PERM sets) for reference.
4. Disease specificity (user 2026-10-12, option A): tumor-vs-normal MRs are dominated by what all tumors share
   (proliferation), so growth-signaling inhibitors reverse almost every disease and anti-estrogens rank higher in basal
   than in luminal A breast cancer. rev_z = robust z of rev within the disease ((rev - median) / (1.4826 MAD));
   The same drug activity (this disease's network and lines) is also scored against the MRs of each reference disease
   (other tissue keys, so sibling subtypes are not counted against each other, with >= MIN_MRS_REF MRs; MRs
   restricted to regulators of this network), each robust-z-scored over the disease's drugs.
   spec = rev_z minus the mean over the reference diseases: reverses this disease more than diseases in general.
   Works for drugs measured only in this disease's tissue (no need for the drug to be scored in other diseases).
   spec_pct = percentile of spec within the disease.
5. Tiers (user 2026-10-12: approved drugs are the top-tier output; the list should not be narrowed by much), both
   kept for comparison. tier_rev: reversing = rev > 0 and top TOP_PCT percent by rev. tier_spec: rev > 0 and top
   TOP_PCT percent by spec. Tier 1 approved + passes target gate + reversing; tier 2 not approved + passes gate +
   reversing; tier 3 everything else. Drugs with > MULTI_TARGET annotated targets are flagged, not removed.
   low_mrs: diseases with < MIN_MRS_REF MRs (scores rest on very few regulators).
ONLY_GROUPS (comma-separated) limits the run, for testing.
Line match (user 2026-10-28, D76): line_match = Celligner weight of the drug's best measured line in the tissue key / the
disease's top in-key line weight (direct drugs; NaN for borrowed). unmatched_lines = direct and line_match < LINE_MATCH_MIN
(0.1, fixed from the distribution before any ranking was looked at: bimodal, 25% of direct pairs < 0.01, 50% > 0.19).
E.g. BRCA-Basal drugs measured only in MCF7 (0.0004). Flagged drugs stay ranked but are kept out of tiers 1-2.
Borrowed drugs (v2, D77): conf and low_transfer from 79l VERSION=v2 (results/lincs/drug_transfer_confidence_v2.tsv; run
79l first); low_transfer (conf < 0.5) is kept out of tiers 1-2 the same way.
Evidence status per disease (user 2026-10-28, D78; summary columns): best_line_match = best line_match any drug could
reach (the disease's best in-key line with LINCS compound data); n_t12_matched / n_t12_borrowed = drugs in tier 1-2 (by
rev or spec) that are direct in a matched line / borrowed with transfer conf >= 0.5; evidence_status = "matched" (>= 1
matched direct drug in tiers 1-2), "borrowed only", or "no adequate compound"; evidence_note = the sentence for reports.
STRATA=1 (VERSION=v2; D79 subgroups, 2026-10-29; replaces 79g for v2): units = the Step 10 v2 subgroups. MRs from the
<set>_strata directories (50_viper.R STRATA=), drug lists results/lincs/stratum_drugs_v2.tsv (79b), line weights
results/celligner/stratum_line_weights_long_v2.tsv (74b); tissue key, network and borrowed-drug transfer confidence of the
parent group. spec: vs the 41-group diseases of other tissue keys, as for groups. spec_sib = rev_z minus the mean rev_z
of the same drug activity scored against the sibling subgroups' MRs; tier_sib as tier_spec. Outputs stratum_reversal*.
SCORE=weighted (default for VERSION=v2; user 2026-10-28, D73) | fdr (v1 default): weighted uses every regulator tested in
the disease, weight w_i = disease NES, no FDR cutoff: rev = -sum_i w_i * a_i / sqrt(sum_i w_i^2) (with w_i = sign of an
FDR < MR_FDR MR this is the formula above). Why: FDR < 0.05 / 0.10 lists jumped between null versions (COAD 425 -> 0)
while the NES ranking stayed put (rho >= 0.98; results/step6/null_adjust_test). Reference diseases for spec: every other
tissue key, scored the same way. n_mrs = regulators weighted; n_eff = (sum w^2)^2 / sum w^4 (how many carry the score);
the set-permutation null permutes the weights over the tested regulators.
VERSION=v2 [REGSET=tf]: v2 disease MRs (results/step6/aracne_project_v2[_tf]) and drug activity on the same regulons
(data/processed/lincs/drug_viper_specific_v2[_tf]); outputs get the same suffix. No raw-profile activity for v2 (rev_raw NaN).
Outputs: results/lincs/drug_reversal{SFX}.tsv (disease x drug), drug_reversal_summary{SFX}.tsv (per disease),
         drug_reversal_controls{SFX}.tsv (known-mechanism drugs, sanity check)
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed/lincs"
RES = ROOT / "results/lincs"
MR_FDR, TOP_PCT, N_PERM, MULTI_TARGET, SEED = 0.10, 10, 200, 5, 0
MIN_MRS_REF, MIN_REF = 20, 10
LINE_MATCH_MIN = 0.1
SFX = ("_v2" if os.environ.get("VERSION", "v1") == "v2" else "") + ("_tf" if os.environ.get("REGSET", "full") == "tf" else "")
assert SFX != "_tf", "REGSET=tf needs VERSION=v2"
# REGULON=collectri | dorothea_ABC (VERSION=v2; user 2026-10-28): disease MRs from results/step6/<REGULON>_v2, drug activity
# from drug_viper_specific_<REGULON>/ALL (79d; one context-free regulon for all diseases); outputs suffixed _v2_<REGULON>
# REGULON=consensus (VERSION=v2; user 2026-10-28, D75): per regulator, the mean over ARACNe v2 TF (the disease's network)
# and CollecTRI of the drug activity, and the mean of the disease NES over the same sets; a regulator's sets = those
# where it is tested in the disease and scored on the drug side, so drug and disease values average the same regulons.
# Inputs results/step6/aracne_project_v2_tf + collectri_v2 and drug_viper_specific_v2_tf + _collectri; outputs _v2_consensus
REGULON = os.environ.get("REGULON", "aracne_project")
CONS = REGULON == "consensus"
CONS_SETS = {"aracne": ("results/step6/aracne_project_v2_tf", "_v2_tf"), "collectri": ("results/step6/collectri_v2", "_collectri")}
CURATED = REGULON not in ("aracne_project", "consensus")
assert not CONS or SFX == "_v2", "consensus: VERSION=v2 (REGSET not used)"
assert not CURATED or SFX == "_v2", "curated regulons: VERSION=v2 (REGSET not used)"
MR_DIR = f"results/step6/{REGULON}_v2" if CURATED else f"results/step6/aracne_project{SFX}"
DTAG = f"_{REGULON}" if CURATED else SFX
if CURATED or CONS:
    SFX = f"_v2_{REGULON}"
SCORE = os.environ.get("SCORE", "weighted" if MR_DIR.endswith(("_v2", "_v2_tf")) or CURATED or CONS else "fdr")
assert SCORE in ("weighted", "fdr")
if SCORE == "fdr" and SFX.startswith("_v2"):
    SFX += "_fdr"
STRATA = os.environ.get("STRATA", "0") == "1"
assert not STRATA or SFX.startswith("_v2"), "STRATA=1 needs VERSION=v2"
PFX = "stratum_reversal" if STRATA else "drug_reversal"
# known mechanism -> expected to reverse the disease (sanity check, not used for scoring)
CONTROLS = [("BRCA-LumA", ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]),
            ("BRCA-LumB", ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]),
            ("TCGA-PRAD", ["bicalutamide", "enzalutamide", "flutamide", "nilutamide"]),
            ("TCGA-COAD", ["fluorouracil", "irinotecan", "oxaliplatin"]),
            ("TCGA-LUAD", ["erlotinib", "gefitinib", "afatinib", "osimertinib"])]
# ER-negative: anti-estrogens should not be specific here
NEG_CONTROLS = [("BRCA-Basal", ["fulvestrant", "tamoxifen", "raloxifene"])]


def net_of(cohort):
    return "COADREAD" if cohort in ("TCGA-COAD", "TCGA-READ") else cohort.split("-", 1)[1]


def main():
    rng = np.random.default_rng(SEED)
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    idx["col"] = np.arange(len(idx))
    grp = pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t").set_index("group_id")
    if CONS:   # set-tagged regulators ("aracne|ESR1"), combined per group below
        mr = pd.concat([pd.read_csv(ROOT / d / "group_mr.tsv", sep="\t").assign(regulator=lambda x, k=k: k + "|" + x.regulator)
                        for k, (d, _) in CONS_SETS.items()], ignore_index=True)
    else:
        mr = pd.read_csv(ROOT / MR_DIR / "group_mr.tsv", sep="\t")
    dd = pd.read_csv(RES / "disease_drugs.tsv", sep="\t", low_memory=False)
    mr_ref, par = mr, {}                      # reference diseases for spec: always the 41 groups
    if STRATA:
        dirs = [(k, d) for k, (d, _) in CONS_SETS.items()] if CONS else [("", MR_DIR)]
        mr = pd.concat([pd.read_csv(ROOT / f"{d}_strata" / "group_mr.tsv", sep="\t")
                        .assign(regulator=lambda x, k=k: (k + "|" if k else "") + x.regulator) for k, d in dirs], ignore_index=True)
        dd = pd.read_csv(RES / "stratum_drugs_v2.tsv", sep="\t", low_memory=False)
        par = dd.drop_duplicates("group_id").set_index("group_id").parent_group.to_dict()
    P = lambda g: par.get(g, g)
    tcf = RES / "drug_transfer_confidence_v2.tsv"
    tc = pd.read_csv(tcf, sep="\t") if SFX.startswith("_v2") else None
    assert tc is not None or not SFX.startswith("_v2"), "run 79l VERSION=v2 first"
    w = (pd.read_csv(RES.parent / "celligner/stratum_line_weights_long_v2.tsv", sep="\t").rename(columns={"stratum": "group_id"})
         if STRATA else pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t"))
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
    models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", index_col=0)
    models = models[models.OncotreePrimaryDisease != "Non-Cancerous"]

    nets = {}
    def read(kind, net, dtag=DTAG):
        if CONS and dtag == DTAG:   # both sets side by side, set-tagged columns
            ra, Aa = read(kind, net, CONS_SETS["aracne"][1]); rc, Ac = read(kind, "ALL", CONS_SETS["collectri"][1])
            if Aa is None or Ac is None:
                return None, None
            return pd.Index(["aracne|" + r for r in ra] + ["collectri|" + r for r in rc]), np.hstack([Aa, Ac])
        f = PROC / f"drug_viper_{kind}{dtag}" / f"{net}.f32"
        if not f.exists():
            return None, None
        regs = pd.Index((PROC / f"drug_viper_{kind}{dtag}" / f"{net}.regulators.txt").read_text().split())
        return regs, np.fromfile(f, dtype=np.float32).reshape(len(idx), len(regs)).astype(np.float64)

    def consensus(D, regs, m):
        """Set-tagged drug activity / disease NES -> per-regulator means over the sets present on both sides."""
        base = m.regulator.str.split("|", n=1).str[1].values
        cols = regs.get_indexer(m.regulator)
        ub, inv = np.unique(base, return_inverse=True)
        Dc = np.zeros((D.shape[0], len(ub))); cnt = np.bincount(inv, minlength=len(ub))
        np.add.at(Dc.T, inv, D[:, cols].T)
        mc = pd.DataFrame({"regulator": ub, "nes": np.bincount(inv, weights=m.nes.values) / cnt, "n_sets": cnt,
                           "fdr": m.groupby(base).fdr.min().reindex(ub).values, "group_id": m.group_id.iloc[0]})
        return Dc / cnt, pd.Index(ub), mc

    def load(net):
        if net not in nets:
            regs, A = read("specific", net)                                   # profiles x regulators
            rr, R = read("raw", net)
            nets[net] = (regs, A, R.T[rr.get_indexer(regs)].T if R is not None else None)
        return nets[net]

    sel = (lambda d: d) if SCORE == "weighted" else (lambda d: d[d.fdr < MR_FDR])
    wt = (lambda d: d.nes.values) if SCORE == "weighted" else (lambda d: np.sign(d.nes.values))
    ok_ref = set(sel(mr_ref).groupby("group_id").size().loc[lambda x: x >= MIN_MRS_REF].index)
    gkey = keys.key
    rz = lambda v: (v - np.median(v)) / (1.4826 * np.median(np.abs(v - np.median(v))))
    rows, summ = [], []
    only = [x for x in os.environ.get("ONLY_GROUPS", "").split(",") if x]
    if only:
        dd = dd[dd.group_id.isin(only)]
    for g, gd in dd.groupby("group_id"):
        regs, A, A_raw = load("ALL" if CURATED else net_of(grp.loc[P(g), "cohort"]))
        k = keys.loc[P(g)]
        km = models[models.OncotreeLineage.isin(str(k.lineages).split(";"))]
        if isinstance(k.diseases, str) and k.diseases:
            km = km[km.OncotreePrimaryDisease.isin(k.diseases.split(";"))]
        wg = w[w.group_id == g].set_index("ModelID").weight

        # drug activity per disease (drugs x scored regulators)
        D = np.empty((len(gd), len(regs)), dtype=np.float64)
        D_raw = np.empty_like(D)
        n_lines = np.zeros(len(gd), int)
        lmatch = np.full(len(gd), np.nan)
        top_w = wg.reindex(km.index).fillna(0).max()
        lincs_key = sorted(set(idx.ModelID) & set(km.index))
        best_lm = float(wg.reindex(lincs_key).fillna(0).max() / top_w) if top_w > 0 and lincs_key else 0.0
        best_line = models.CellLineName.get(wg.reindex(lincs_key).fillna(0).idxmax()) if lincs_key else ""
        by_drug = idx.groupby("pert_id")
        for i, (d, src) in enumerate(zip(gd.pert_id, gd.source)):
            p = by_drug.get_group(d)
            if src == "direct":
                p = p[p.ModelID.isin(km.index)]
                ww = wg.reindex(p.ModelID).fillna(0).values
                lmatch[i] = ww.max() / top_w if top_w > 0 and len(ww) else 0.0
            else:
                p = p[~p.ModelID.isin(km.index)]
                ww = np.zeros(len(p))
            ww = ww if ww.sum() > 0 else np.ones(len(p))
            D[i] = (ww[:, None] * A[p.col.values]).sum(0) / ww.sum()
            if A_raw is not None:
                D_raw[i] = (ww[:, None] * A_raw[p.col.values]).sum(0) / ww.sum()
            n_lines[i] = len(p)

        m = mr[(mr.group_id == g) & mr.regulator.isin(regs)]
        if CONS:
            regs0, D0 = regs, D
            D, regs, m = consensus(D0, regs0, m)
            D_raw = None
        mrs = sel(m)
        tested = regs.get_indexer(m.regulator)
        ci = regs.get_indexer(mrs.regulator)
        s = wt(mrs)
        n = len(ci)
        norm = np.sqrt((s ** 2).sum())
        n_eff = (s ** 2).sum() ** 2 / ((s ** 4).sum() + 1e-12)
        if n < 2:   # nothing to reverse (e.g. no MR at MR_FDR): group left out, reported
            print(f"{g:<20} MRs {n:>4}  skipped (< 2 MRs at FDR < {MR_FDR})", flush=True)
            continue
        rev = -(D[:, ci] * s).sum(1) / norm
        dn = mrs.nes.values
        Dc = D[:, ci] - D[:, ci].mean(1, keepdims=True)
        r_rev = -(Dc @ (dn - dn.mean())) / (np.sqrt((Dc ** 2).sum(1)) * np.linalg.norm(dn - dn.mean()) + 1e-12)
        null = np.empty((len(gd), N_PERM))
        for j in range(N_PERM):
            null[:, j] = -(D[:, rng.choice(tested, n, replace=False)] * rng.permutation(s)).sum(1) / norm

        o = gd[["group_id", "pert_id", "cmap_name", "source", "passes_gate", "gate_reason", "n_targets",
                "hub_moa", "hub_clinical_phase", "chembl_max_phase", "pharmacology_flag"]].copy()
        o["n_lines_used"] = n_lines
        o["line_match"] = np.round(lmatch, 4)
        o["unmatched_lines"] = o.line_match < LINE_MATCH_MIN
        if tc is not None:
            t_ = tc[tc.group_id == P(g)].set_index("pert_id")
            o["transfer_conf"] = o.pert_id.map(t_.conf).values
            o["low_transfer"] = o.pert_id.map(t_.low_transfer).fillna(False).astype(bool).values
        else:
            o["low_transfer"] = False
        o["n_mrs"] = n
        o["rev"] = rev.round(3)
        o["rev_raw"] = (-(D_raw[:, ci] * s).sum(1) / norm).round(3) if A_raw is not None else np.nan
        o["r_rev"] = r_rev.round(3)
        o["rev_pct"] = (pd.Series(rev).rank(method="max").values - 1) / len(rev) * 100
        o["rev_pct"] = o.rev_pct.round(1)
        o["approved"] = (o.chembl_max_phase == 4) | (o.hub_clinical_phase == "Launched")
        o["multi_target"] = o.n_targets > MULTI_TARGET
        reversing = (o.rev > 0) & (o.rev_pct >= 100 - TOP_PCT)
        ok_t = o.passes_gate & ~o.unmatched_lines & ~o.low_transfer
        o["tier_rev"] = np.where(reversing & ok_t & o.approved, 1, np.where(reversing & ok_t, 2, 3))
        mad = 1.4826 * np.median(np.abs(rev - np.median(rev)))
        o["rev_z"] = ((rev - np.median(rev)) / mad).round(3)
        refz = []
        for g2 in sorted(ok_ref):
            if g2 == g or gkey[g2] == gkey[P(g)]:
                continue
            m2 = sel(consensus(D0, regs0, mr_ref[(mr_ref.group_id == g2) & mr_ref.regulator.isin(regs0)])[2] if CONS
                     else mr_ref[(mr_ref.group_id == g2) & mr_ref.regulator.isin(regs)])
            m2 = m2[m2.regulator.isin(regs)]
            if len(m2) < MIN_MRS_REF:
                continue
            c2 = regs.get_indexer(m2.regulator)
            w2 = wt(m2)
            refz.append(rz(-(D[:, c2] * w2).sum(1) / np.sqrt((w2 ** 2).sum())))
        o["n_ref_diseases"] = len(refz)
        o["spec"] = (o.rev_z.values - np.mean(refz, axis=0)).round(3) if len(refz) >= MIN_REF else np.nan
        o["spec_pct"] = (o.spec.rank(method="max", pct=True) * 100).round(1)
        if STRATA:   # vs the sibling subgroups of the same parent
            sibz = []
            for g2 in sorted(set(mr.group_id) - {g}):
                if P(g2) != P(g):
                    continue
                m2 = sel(consensus(D0, regs0, mr[(mr.group_id == g2) & mr.regulator.isin(regs0)])[2] if CONS
                         else mr[(mr.group_id == g2) & mr.regulator.isin(regs)])
                m2 = m2[m2.regulator.isin(regs)]
                if len(m2) >= 2:
                    w2 = wt(m2); sibz.append(rz(-(D[:, regs.get_indexer(m2.regulator)] * w2).sum(1) / np.sqrt((w2 ** 2).sum())))
            o.insert(1, "parent_group", P(g))
            o["siblings"] = ";".join(sorted({x for x in mr.group_id if P(x) == P(g)} - {g}))
            o["spec_sib"] = (o.rev_z.values - np.mean(sibz, axis=0)).round(3) if sibz else np.nan
            o["spec_sib_pct"] = (o.spec_sib.rank(method="max", pct=True) * 100).round(1)
            sp_ = (o.rev > 0) & (o.spec_sib_pct >= 100 - TOP_PCT)
            o["tier_sib"] = np.where(sp_ & ok_t & o.approved, 1, np.where(sp_ & ok_t, 2, 3))
        rows.append(o)
        t1 = o[o.tier_rev == 1].sort_values("rev", ascending=False)
        summ.append(dict(group_id=g, parent_group=P(g), net=net_of(grp.loc[P(g), "cohort"]), n_mrs=n, n_eff=round(float(n_eff), 1), n_tested=len(tested), score=SCORE,
                         best_line_match=round(best_lm, 4), best_line=best_line, n_key_lines_lincs=len(lincs_key),
                         n_drugs=len(o), n_reversing=int(reversing.sum()), n_tier1_rev=len(t1),
                         n_tier2_rev=int((o.tier_rev == 2).sum()), n_tier1_multi_target=int(t1.multi_target.sum()),
                         null_sd_regset=round(float(null.std()), 3), rev_sd=round(float(rev.std()), 3),
                         top_tier1_rev=", ".join(f"{a} {b:.1f}" for a, b in zip(t1.cmap_name.head(8), t1.rev.head(8)))))
        print(f"{g:<20} MRs {n:>4}  drugs {len(o):>5}  reversing {int(reversing.sum()):>4}  tier1 {len(t1):>3}  "
              f"{summ[-1]['top_tier1_rev'][:90]}", flush=True)

    out = pd.concat(rows, ignore_index=True)
    summ = pd.DataFrame(summ)

    ok_ref = set(summ.loc[summ.n_mrs >= MIN_MRS_REF, "group_id"])
    sp = (out.rev > 0) & (out.spec_pct >= 100 - TOP_PCT)
    ok_t = out.passes_gate & ~out.unmatched_lines & ~out.low_transfer
    out["tier_spec"] = np.where(sp & ok_t & out.approved, 1, np.where(sp & ok_t, 2, 3))
    out["low_mrs"] = ~out.group_id.isin(ok_ref)
    out = out.sort_values(["group_id", "tier_spec", "spec"], ascending=[True, True, False])
    out["rank_spec"] = out.groupby("group_id").cumcount() + 1
    out["rank_rev"] = out.groupby("group_id").rev.rank(ascending=False, method="first").astype(int)
    out.to_csv(RES / f"{PFX}{SFX}.tsv", sep="\t", index=False)
    s2 = []
    for g, o in out.groupby("group_id"):
        t1 = o[o.tier_spec == 1]
        t12 = (o.tier_rev <= 2) | (o.tier_spec <= 2)
        n_m = int((t12 & (o.source == "direct")).sum()); n_b = int((t12 & (o.source == "borrowed")).sum())
        status = "matched" if n_m else ("borrowed only" if n_b else "no adequate compound")
        s2.append(dict(group_id=g, n_t12_matched=n_m, n_t12_borrowed=n_b, evidence_status=status,
                       n_tier1_spec=len(t1), n_tier2_spec=int((o.tier_spec == 2).sum()),
                       n_tier1_both=int(((o.tier_spec == 1) & (o.tier_rev == 1)).sum()),
                       low_mrs=bool(o.low_mrs.iloc[0]),
                       top_tier1_spec=", ".join(f"{a} {b:.1f}" for a, b in zip(t1.cmap_name.head(8), t1.spec.head(8)))))
    summ = summ.merge(pd.DataFrame(s2), on="group_id")
    def note(r):
        lm = (f"No compound was tested in a line resembling this disease (best line match {r.best_line_match:.4f}, "
              f"{r.best_line or 'no line'}; cutoff {LINE_MATCH_MIN})") if r.best_line_match < LINE_MATCH_MIN else \
             f"Compounds tested in matched lines (best {r.best_line}, line match {r.best_line_match:.2f})"
        if r.evidence_status == "matched":
            return f"{lm}; {r.n_t12_matched} tier 1-2 drug(s) rest on matched-line data."
        if r.evidence_status == "borrowed only":
            return (f"{lm}. The {r.n_t12_borrowed} candidate(s) are borrowed from other tissues and rest on the drug acting "
                    f"consistently across cell lines (transfer confidence >= 0.5), not on data from this disease's lines.")
        return f"{lm}, and no borrowed drug passes; an adequate compound cannot be reported."
    summ["evidence_note"] = summ.apply(note, axis=1)
    summ.to_csv(RES / f"{PFX}_summary{SFX}.tsv", sep="\t", index=False)

    ctl = []
    for expected, cl in [("high", CONTROLS), ("low", NEG_CONTROLS)]:
        for g, names in cl:
            x = out[out.group_id == g]
            for nm in names:
                r = x[x.cmap_name == nm]
                if len(r):   # several LINCS IDs per drug: best one
                    ctl.append(dict(group_id=g, drug=nm, expected=expected, source=r.source.iloc[0],
                                    n_ids=len(r), rev_pct=r.rev_pct.max(), spec_pct=r.spec_pct.max(),
                                    tier_rev=r.tier_rev.min(), tier_spec=r.tier_spec.min()))
                else:
                    ctl.append(dict(group_id=g, drug=nm, expected=expected, source="not in list"))
    ctl = pd.DataFrame(ctl)
    ctl.to_csv(RES / f"{PFX}_controls{SFX}.tsv", sep="\t", index=False)
    print(ctl.to_string(index=False))


if __name__ == "__main__":
    main()
