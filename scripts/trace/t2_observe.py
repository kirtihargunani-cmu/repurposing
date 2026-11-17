#!/usr/bin/env python
"""Mechanism trace, stage 2: links L5, L6, L7 and the verdict.

Procedure: docs/PROCEDURE_mechanism_trace.md (Sections 2-8, 10).

Written and run ONLY after scripts/trace/t1_predict.py and its predicted footprints were
committed (the Section 8 step 4 freeze, commit "TRACE FREEZE"). Every `expected` string
below is a literal, written from the frozen L1-L4 record before the corresponding
observation is read; none is rewritten after seeing the data.

NOT read anywhere in this file: results/lincs/d128/a4/* (the hub-dominance result that
L6 is supposed to find on its own), the D130 corrected-null tables, and any outcome label
(CTRP / PRISM / GDSC). No score is written to any production file, nothing is ranked.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
OUT = ROOT / "results/lincs/trace"
A1 = ROOT / "results/lincs/d128/a1"
A2 = ROOT / "results/lincs/d128/a2"
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"
F_DIVISOR = 2.0          # P9, frozen in scripts/d128/config.py; used only to restate f(z), not to make a new score

PAIRS = json.loads((OUT / "t1_run.json").read_text())["pairs"]

# Section 3: the underlying SOURCE each link draws on. Two links with the same source class are one link
# reported twice. Declared per pair because the O component of L2 has a different source in each.
SOURCES = {
    "fulvestrant_ESR1_BRCA_LumA": {
        "L1": ["chembl_binding"],
        # O for this pair is a VIPER NES over the same CollecTRI regulon L3/L4 use -> no new source from O
        "L2": ["depmap_dependency", "tcga_alterations", "collectri_regulon"],
        "L3": ["collectri_regulon"],                 # zero-hop anchor: no signalling prior is used at all
        "L4": ["collectri_regulon"],
        "L5": ["lincs_response"],
        "L6": [],                                     # an audit of L3-L5; adds no source
        "L7": ["tcga_gtex_contrast"],
    },
    "vemurafenib_BRAF_TCGA_SKCM": {
        "L1": ["chembl_binding"],
        # O for this pair is logFC vs GTEx -> the same source L7 uses
        "L2": ["depmap_dependency", "tcga_alterations", "tcga_gtex_contrast"],
        "L3": ["omnipath_signed_prior", "collectri_regulon"],
        "L4": ["collectri_regulon"],
        "L5": ["lincs_response"],
        "L6": [],
        "L7": ["tcga_gtex_contrast"],
    },
}
# Section 3 entanglements between links that DO each bring a new source class, and so would otherwise be
# counted twice. Both traces use a direct tissue match, which is the case Section 3 names.
ENTANGLED = [("L5", "L7",
              "direct tissue match: the cell line and the tumour cohort are the same tissue, so L5 and L7 share "
              "the tissue's biology (Section 3); with a borrowed line they would not")]

# ---- expectations, written from the frozen L1-L4 record BEFORE the observation is read --------------------------
EXPECT = {
    "L5": {
        "fulvestrant_ESR1_BRCA_LumA":
            "If the chain holds, the 478 frozen genes should move in the predicted direction in MCF7 under "
            "fulvestrant more often than not: concordance well above 50%, carried by the 405 expected-DOWN genes "
            "(ESR1-activated targets losing their activator). Two sub-predictions that can fail independently: "
            "(a) landmark (directly measured) genes should be at least as concordant as inferred genes; "
            "(b) the 8.6% of genes uniquely claimed by ESR1 should be at least as concordant as the 91.4% that "
            "other regulators also claim - if they are not, the movement is not specifically ESR1's.",
        "vemurafenib_BRAF_TCGA_SKCM":
            "If the chain holds, the 1063 frozen genes should move in the predicted direction in A375 under "
            "vemurafenib more often than not, carried by the 903 expected-DOWN genes. Same two sub-predictions: "
            "(a) landmark >= inferred; (b) the 17.0% uniquely attributable genes >= the 83.0% shared ones. A375 "
            "is BRAF-V600E, so this is the permissive case: if concordance fails here it fails everywhere.",
    },
    "L6": {
        "fulvestrant_ESR1_BRCA_LumA":
            "Structurally this link is already pinned: L3 found exactly ONE regulator, so 100% of the predicted "
            "evidence is one regulon and the only question L6 can still ask is how concentrated the OBSERVED "
            "signal is across genes. Expect the gene-level participation ratio to be a large fraction of the "
            "covered genes (no single gene should carry the link), while the regulator-level concentration is 1.00 "
            "by construction. The honest reading is that this trace has one mechanistic source whatever the genes do.",
        "vemurafenib_BRAF_TCGA_SKCM":
            "L4 already shows one regulator claiming 65.1% of the predicted genes. If the observed signal is "
            "distributed like the structure, the same regulator should carry roughly that share of the observed "
            "|contribution| too, and removing it should move the link's mean contribution substantially. If it "
            "does, 1063 genes are not 1063 pieces of evidence and L6 breaks the appearance of breadth.",
    },
    "L7": {
        "fulvestrant_ESR1_BRCA_LumA":
            "If the chain closes, the genes the drug is predicted to push DOWN should be UP in BRCA-LumA vs its "
            "normal reference, and vice versa - i.e. expected_sign == -sign(logFC) for a clear majority of the "
            "significant footprint genes. ESR1's own O component is a positive NES in this cohort, so the "
            "luminal-A tumour is predicted to be running the ESR1 programme that fulvestrant removes.",
        "vemurafenib_BRAF_TCGA_SKCM":
            "Same closure test. But L2 already recorded a warning that must be carried here: BRAF's expression "
            "logFC vs GTEx normal skin is NEGATIVE (-0.62, O percentile 35.9), so the normal reference may already "
            "run the programme, and the D128-A3.2 reference artefact is live for this pair. Expect closure to be "
            "weaker than fulvestrant's, or inverted.",
    },
}


def order_before(a, b):
    return int(a[1:]) < int(b[1:])


def participation_ratio(x):
    """Effective number of contributors: (sum|x|)^2 / sum(x^2). n for equal contributors, 1 for a single one."""
    x = np.abs(np.asarray(x, dtype=float))
    s2 = float((x ** 2).sum())
    return float(s2 and (float(x.sum()) ** 2) / s2)


def main():
    genes_tab = pd.read_csv(ROOT / "data/processed/lincs/drug_line_profiles_genes.tsv", sep="\t")
    idx = pd.read_csv(ROOT / "data/processed/lincs/drug_line_profiles_index.tsv", sep="\t")
    annot = pd.read_csv(ROOT / "results/lincs/drug_annotations.tsv", sep="\t")
    E = pd.read_parquet(A1 / "regulon_edges.parquet")
    reg_of_gene = E[E.in_lincs_gene_space].groupby("gene_id").regulator.apply(set)
    prolif = pd.read_csv(ROOT / "results/lincs/prolif_regulators_p2_v2_consensus.tsv", sep="\t")
    prolif_collectri = prolif[(prolif.net == "ALL") & (prolif.set == "collectri")].set_index("regulator")
    C = pd.read_parquet(A2 / "gene_contributions.parquet")
    SB = pd.read_parquet(A2 / "scores_by_line.parquet")

    M = np.memmap(ROOT / "data/processed/lincs/drug_line_profiles_bing_specific.f32", dtype="float32",
                  mode="r", shape=(len(idx), len(genes_tab)))
    gpos = {g: i for i, g in enumerate(genes_tab.gene_id)}

    links, summary = [], {}
    for key, p in PAIRS.items():
        fp = pd.read_csv(OUT / f"predicted_footprint_{key}.tsv", sep="\t")
        l3 = pd.read_csv(OUT / f"l3_paths_{key}.tsv", sep="\t")
        signed = fp[fp.expected_sign != 0].copy()

        # ------------------------------------------------------------------ L5: predicted -> observed
        # DISCREPANCY (reported, not worked around): results/lincs/d128/a2/gene_contributions.parquet is a SEEDED
        # SAMPLE (a2_score.py --gene-contrib sample, n=2000, gene_arm all_bing only), covering 186 of the 257
        # scored units. vemurafenib::CHEMBL5145 x ACH-000219 is in it; fulvestrant::CHEMBL206 x ACH-000019 is NOT.
        # So the observation is rebuilt here from the frozen footprint and the same profile matrix A2 used, with
        # score.py's own formula f(z) = tanh(z / 2), contribution = expected_sign * f(z); it is then checked gene
        # by gene against the A2 sample where that exists, and its mean against scores_by_line.M in every case.
        pid = sorted(set(annot[annot.cmap_name == p["cmap_name"]].pert_id))
        rows = idx[(idx.ModelID == p["ModelID"]) & (idx.pert_id.isin(pid))].sort_values("n_sig", ascending=False)
        assert len(rows) >= 1, (key, "no profile row for this drug x line")
        prow = int(rows.index[0])
        obs = C[(C.unit_id == p["unit_id"]) & (C.ModelID == p["ModelID"]) & (C.score == "M")]
        cov = signed.copy()
        cov["z"] = [float(M[prow, gpos[g]]) if g in gpos else np.nan for g in cov.gene_id]
        cov = cov[np.isfinite(cov.z)].copy()
        cov["f_z"] = np.tanh(cov.z / F_DIVISOR)
        cov["contribution"] = cov.expected_sign * cov.f_z
        # check 1: gene-by-gene against the A2 sample, when A2 sampled this unit x line
        chk = cov.merge(obs[["gene_id", "pert_id", "expected_sign", "z", "contribution"]], on="gene_id",
                        how="inner", suffixes=("", "_a2"))
        chk = chk[chk.pert_id == idx.loc[prow].pert_id]
        a2_check = (f"{len(chk)} genes checked against the A2 sample, max |z| difference "
                    f"{float((chk.z - chk.z_a2).abs().max()):.2e}, "
                    f"{int((chk.expected_sign != chk.expected_sign_a2).sum())} sign disagreements"
                    if len(chk) else "A2 did not sample this unit x line (see note above); no gene-level check")
        bad = int((chk.expected_sign != chk.expected_sign_a2).sum()) if len(chk) else 0
        # check 2: the mean against the frozen D128-A2 score for this (cohort, pert, line, unit)
        sb = SB[(SB.unit_id == p["unit_id"]) & (SB.ModelID == p["ModelID"]) &
                (SB.pert_id == idx.loc[prow].pert_id) & (SB.cohort == p["group_id"]) &
                (SB.gene_arm == "all_bing")]
        m_a2 = float(sb.M.iloc[0]) if len(sb) else np.nan
        m_here = float(cov.contribution.mean())
        m_check = (f"mean rebuilt contribution {m_here:+.6f} vs D128-A2 M {m_a2:+.6f} "
                   f"(difference {abs(m_here - m_a2):.2e}, {int(sb.M_n_genes.iloc[0]) if len(sb) else -1} genes in A2)")
        reps = "; ".join(f"{r.pert_id} ({int(r.n_sig)} sig)" for r in rows.itertuples())
        cov["concordant"] = np.sign(cov.z) == cov.expected_sign
        conc = float(cov.concordant.mean())
        by_space = cov.groupby("feature_space").concordant.agg(["mean", "size"])
        by_uniq = cov.groupby("uniquely_attributable").concordant.agg(["mean", "size"])

        # generic-response control (Section 6): the SAME frozen footprint under every other compound in this line
        line_rows = idx.index[idx.ModelID == p["ModelID"]].to_numpy()
        cols = np.array([gpos[g] for g in cov.gene_id])
        e = cov.expected_sign.to_numpy(float)
        Z = np.asarray(M[line_rows][:, cols], dtype=float)
        means = (np.tanh(Z / F_DIVISOR) * e).mean(axis=1)
        this_i = int(np.where(line_rows == prow)[0][0])
        this_mean = float(means[this_i])
        pct_generic = float((means < this_mean).mean() * 100)
        conc_all = (np.sign(Z) == e).mean(axis=1)
        pct_conc = float((conc_all < conc_all[this_i]).mean() * 100)
        # amplitude control (Section 6, D131/D133)
        ampl = np.abs(np.asarray(M[line_rows], dtype=float)).mean(axis=1)
        pct_ampl = float((ampl < ampl[this_i]).mean() * 100)
        r_ampl = float(np.corrcoef(ampl, means)[0, 1])

        links.append(dict(
            pair=key, link="L5",
            question="Did those genes move as predicted, in which cell line, matched how?",
            expected=EXPECT["L5"][key], observed=(
                f"line {p['line']} ({p['ModelID']}), pert_id {idx.loc[prow].pert_id}, "
                f"{int(idx.loc[prow].n_sig)} signature(s), specific (drug-specific) profile arm. "
                f"{len(cov)} of {len(signed)} frozen genes have a response value. "
                f"Replicate profiles for this drug in this line: {reps}; the one with most signatures is used. "
                f"Sign concordance with the frozen prediction: {100 * conc:.1f}%. "
                f"By feature space: " + ", ".join(f"{k} {100 * v['mean']:.1f}% (n={int(v['size'])})"
                                                  for k, v in by_space.iterrows()) + ". "
                f"Uniquely attributable genes {100 * by_uniq.loc[True, 'mean']:.1f}% "
                f"(n={int(by_uniq.loc[True, 'size'])}) vs shared-attribution "
                f"{100 * by_uniq.loc[False, 'mean']:.1f}% (n={int(by_uniq.loc[False, 'size'])}). "
                f"Mean expected-signed f(z) over the footprint = {this_mean:+.4f}; the same frozen footprint "
                f"under the {len(line_rows)} compounds profiled in this line puts it at the "
                f"{pct_generic:.1f}th percentile (concordance at the {pct_conc:.1f}th). "
                f"Reproduction checks: {a2_check}; {m_check}."),
            subpredictions=[
                dict(id="a", statement="landmark (directly measured) genes at least as concordant as inferred",
                     value=f"{100 * float(by_space.loc['landmark', 'mean']):.1f}% vs "
                           f"{100 * float(by_space.loc['best inferred', 'mean']):.1f}%",
                     passed=bool(by_space.loc["landmark", "mean"] >= by_space.loc["best inferred", "mean"])),
                dict(id="b", statement=("genes uniquely claimed by this unit's regulators at least as concordant "
                                        "as shared-attribution genes; failing means the movement is not "
                                        "specifically this target's (Section 5)"),
                     value=f"{100 * float(by_uniq.loc[True, 'mean']):.1f}% (n={int(by_uniq.loc[True, 'size'])}) "
                           f"vs {100 * float(by_uniq.loc[False, 'mean']):.1f}% "
                           f"(n={int(by_uniq.loc[False, 'size'])})",
                     passed=bool(by_uniq.loc[True, "mean"] >= by_uniq.loc[False, "mean"])),
            ],
            verdict=("holds" if conc >= 0.60 and pct_generic >= 90 else
                     "weak" if conc > 0.50 else "breaks"),
            provenance=("results/lincs/d128/a2/gene_contributions.parquet (D128-A2, arm 'specific'); "
                        "data/processed/lincs/drug_line_profiles_bing_specific.f32 for the generic-response and "
                        "amplitude controls; results/celligner/disease_line_weights_long.tsv for the line match"),
            independence=(
                f"{p['line']} is a DIRECT tissue match for {p['group_id']} (celligner key "
                f"{'Breast' if 'BRCA' in p['group_id'] else 'Skin'}, in_key=True), not a borrowed line, so by "
                "Section 3 L5 and L7 are NOT independent: the same tissue biology sets both the cell line's "
                "baseline and the tumour contrast. The line carries only "
                f"{'2.99' if 'BRCA' in p['group_id'] else '0.84'}% of this cohort's celligner weight, so it is "
                "one line's response standing for the cohort. L5 shares the footprint with L3/L4, so a regulon "
                "error moves the prediction and the test of it together."),
            concentration=(
                f"largest single gene's share of |contribution| over the covered footprint: "
                f"{float(cov.contribution.abs().max() / cov.contribution.abs().sum()):.4f}; "
                f"effective number of contributing genes (participation ratio) "
                f"{participation_ratio(cov.contribution.fillna(0)):.0f} of {len(cov)}."),
            alternatives=[
                (f"generic response: the same frozen footprint scores at the {pct_generic:.1f}th percentile of all "
                 f"{len(line_rows)} compounds profiled in {p['line']}. "
                 + ("Not ruled out - a large share of compounds in this line move these genes the same way."
                    if pct_generic < 90 else
                    "Largely ruled out at this link - few other compounds in this line move these genes this way.")),
                (f"amplitude (D131/D133): this drug-line's mean |z| is at the {pct_ampl:.1f}th percentile of the "
                 f"line; across the line's compounds the footprint mean correlates r={r_ampl:+.2f} with overall "
                 "response magnitude. "
                 + ("Not ruled out: in this line the footprint mean largely tracks how big the response is."
                    if abs(r_ampl) >= 0.5 else "Weakly supported as a confound at this link.")),
                ("cytotoxicity: dna_acting=False for this drug in the ranking row, and the proliferation check is "
                 "carried at L6 through the regulators, not here."),
                ("off-target: concordance cannot distinguish the annotated target from any off-target that "
                 "converges on the same regulators."),
            ]))

        # ------------------------------------------------------------------ L6: concentration
        own_tfs = [t for t in l3.tf if t in set(E.regulator)]
        share, uniq_share, pr_by_tf = {}, {}, {}
        tot = float(cov.contribution.abs().sum())
        for t in own_tfs:
            m = cov.gene_id.map(lambda g: t in reg_of_gene.get(g, set()))
            share[t] = float(cov.contribution[m].abs().sum() / tot) if tot else np.nan
            only_t = cov.gene_id.map(
                lambda g: reg_of_gene.get(g, set()) & set(own_tfs) == {t})
            uniq_share[t] = float(cov.contribution[only_t].abs().sum() / tot) if tot else np.nan
        hub = max(share, key=share.get)
        hub_mask = cov.gene_id.map(lambda g: hub in reg_of_gene.get(g, set()))
        drop_hub = cov[~hub_mask]
        mean_all = float((cov.expected_sign * cov.f_z).mean())
        mean_drop = float((drop_hub.expected_sign * drop_hub.f_z).mean()) if len(drop_hub) else np.nan
        conc_drop = float((np.sign(drop_hub.z) == drop_hub.expected_sign).mean()) if len(drop_hub) else np.nan
        pr_tf = participation_ratio(list(share.values()))
        hub_prolif = prolif_collectri.r_prolif.get(hub, np.nan), bool(prolif_collectri.prolif_linked.get(hub, False))
        summary[key] = dict(
            n_regulators=len(own_tfs), hub=hub, hub_structural_share=float(
                (signed.gene_id.map(lambda g: hub in reg_of_gene.get(g, set()))).mean()),
            hub_signal_share=share[hub], hub_unique_signal_share=uniq_share[hub],
            effective_regulators=pr_tf, effective_genes=participation_ratio(cov.contribution.fillna(0)),
            n_covered=len(cov), concordance=conc, concordance_without_hub=conc_drop,
            mean_with_hub=mean_all, mean_without_hub=mean_drop,
            generic_percentile=pct_generic, amplitude_percentile=pct_ampl, amplitude_r=r_ampl,
            hub_r_prolif=hub_prolif[0], hub_prolif_linked=hub_prolif[1])
        pd.DataFrame([dict(pair=key, regulator=t, structural_genes=int(
            signed.gene_id.map(lambda g: t in reg_of_gene.get(g, set())).sum()),
            signal_share=share[t], exclusive_signal_share=uniq_share[t],
            r_prolif=prolif_collectri.r_prolif.get(t, np.nan),
            prolif_linked=bool(prolif_collectri.prolif_linked.get(t, False))) for t in own_tfs]
        ).sort_values("signal_share", ascending=False).to_csv(
            OUT / f"l6_concentration_{key}.tsv", sep="\t", index=False)

        links.append(dict(
            pair=key, link="L6",
            question="Is the signal carried by many independent pieces of evidence, or one?",
            expected=EXPECT["L6"][key], observed=(
                f"{len(own_tfs)} regulator(s) contribute. Largest single regulator by observed |contribution|: "
                f"{hub}, carrying {100 * share[hub]:.1f}% of the link's signal "
                f"({100 * uniq_share[hub]:.1f}% exclusively, i.e. through genes no other regulator of this unit "
                f"claims); structurally it claims "
                f"{100 * summary[key]['hub_structural_share']:.1f}% of the frozen genes. Effective number of "
                f"independent regulators (participation ratio over regulator shares) = {pr_tf:.2f} of "
                f"{len(own_tfs)}. Effective number of contributing genes = "
                f"{summary[key]['effective_genes']:.0f} of {len(cov)}. Dropping {hub}'s genes entirely leaves "
                + (f"{len(drop_hub)} genes, concordance {100 * conc_drop:.1f}% (vs {100 * conc:.1f}% with them) and "
                   f"mean expected-signed f(z) {mean_drop:+.4f} (vs {mean_all:+.4f}). "
                   if len(drop_hub) else
                   "NOTHING - the hub is the unit's only regulator, so the leave-one-out test is not evaluable "
                   "and the link is, by construction, a single piece of evidence. ")
                + f"{hub} proliferation correlation in the CollecTRI set: r_prolif="
                f"{hub_prolif[0]}, prolif_linked={hub_prolif[1]}."),
            verdict=("breaks" if share[hub] >= 0.5 and len(own_tfs) > 1 else
                     "weak" if len(own_tfs) == 1 or share[hub] >= 0.5 else "holds"),
            provenance=("computed in this script from results/lincs/d128/a1/regulon_edges.parquet, the frozen "
                        "footprint and results/lincs/d128/a2/gene_contributions.parquet; proliferation link from "
                        "results/lincs/prolif_regulators_p2_v2_consensus.tsv (collectri / ALL). "
                        "No D128-A4 output was read before this record was written."),
            independence=("this link is an audit OF L3-L5, not a new source: it shares the regulon with L3/L4 and "
                          "the response with L5. It cannot add evidence, only remove the appearance of it. "
                          "The proliferation annotation is the one outside input, and it shares the CollecTRI "
                          "regulon with L4."),
            concentration=(f"{share[hub]:.3f} - {hub}'s share of the observed signal; effective regulators "
                           f"{pr_tf:.2f}."),
            alternatives=[
                (f"the hub regulator {hub} is proliferation-linked in the CollecTRI set (r_prolif "
                 f"{hub_prolif[0]}), so a general slow-down of a proliferating line reproduces most of this "
                 "link's evidence. The cytotoxicity / generic-response alternative is NOT ruled out."
                 if hub_prolif[1] else
                 f"the hub regulator {hub} is not flagged proliferation-linked in the CollecTRI set "
                 f"(r_prolif {hub_prolif[0]}), which weakens but does not exclude a proliferation explanation."),
                ("the regulator shares are not disjoint: a gene claimed by two of the unit's regulators is "
                 "counted in both, so the exclusive share is the conservative reading."),
                ("D130 / D133: within a cohort the regulator activities behave as ~10 independent terms and a "
                 "~1000-gene footprint as ~44; a gene count is not an evidence count even before the hub is "
                 "considered."),
            ]))

        # ------------------------------------------------------------------ L7: disease closure
        con = ROOT / f"data/processed/step4_v2/{p['group_id']}__gtex.tsv.gz"
        DE = pd.read_csv(con, sep="\t")[["gene_id", "symbol", "logFC", "z", "adj.P.Val"]]
        j = cov.merge(DE, on="gene_id", how="inner", suffixes=("", "_de"))
        sig = j[j["adj.P.Val"] < 0.05].copy()
        sig["reverses"] = np.sign(sig.logFC) == -sig.expected_sign
        rev = float(sig.reverses.mean())
        wrev = float((sig.reverses * sig.z.abs()).sum() / sig.z.abs().sum())
        hub_sig = sig[sig.gene_id.map(lambda g: hub in reg_of_gene.get(g, set()))]
        nonhub_sig = sig[~sig.gene_id.map(lambda g: hub in reg_of_gene.get(g, set()))]
        tgt = DE[DE.symbol == p["target"]]
        links.append(dict(
            pair=key, link="L7",
            question=("Are those same genes dysregulated in the tumour in the direction the drug reverses?"),
            expected=EXPECT["L7"][key], observed=(
                f"contrast {p['group_id']} vs gtex (P2 order: gtex first). {len(j)} frozen+covered genes present "
                f"in the contrast, {len(sig)} with adj.P.Val < 0.05 (P8). The drug's predicted direction opposes "
                f"the tumour's direction for {100 * rev:.1f}% of them ({100 * wrev:.1f}% when weighted by |z|). "
                f"Split by the L6 hub: {hub}-claimed genes {100 * float(hub_sig.reverses.mean()):.1f}% "
                f"(n={len(hub_sig)}), other genes {100 * float(nonhub_sig.reverses.mean()):.1f}% "
                f"(n={len(nonhub_sig)}). The target itself: "
                + (f"{p['target']} logFC {float(tgt.logFC.iloc[0]):+.2f}, z {float(tgt.z.iloc[0]):+.2f}, "
                   f"adj.P {float(tgt['adj.P.Val'].iloc[0]):.2e}" if len(tgt) else "not in the contrast table")),
            verdict=("holds" if rev >= 0.60 else "weak" if rev > 0.50 else "breaks"),
            provenance=(f"data/processed/step4_v2/{p['group_id']}__gtex.tsv.gz (limma + voomByGroup, P1: logFC "
                        "sign, z weight; P8: adj.P.Val < 0.05)"),
            independence=("the tumour contrast is a source neither L5 nor L4 uses, so L7 is the one genuinely "
                          "new measurement after L2. But its gene set is the frozen footprint, so a regulon error "
                          "still selects which genes are asked about; and the normal reference here is the same "
                          "GTEx reference L2's O component used"
                          + (" (for this pair O came from a VIPER NES over the same regulon instead, which is "
                             "worse: L2's O, L4 and L7's gene selection then all rest on one regulon)."
                             if key.startswith("fulvestrant") else ".")),
            concentration=(f"largest single gene's share of the |z| weight used here: "
                           f"{float(sig.z.abs().max() / sig.z.abs().sum()):.4f}; "
                           f"effective number of genes {participation_ratio(sig.z):.0f} of {len(sig)}. "
                           f"{hub}-claimed genes are {100 * len(hub_sig) / max(len(sig), 1):.1f}% of the "
                           "closure evidence."),
            alternatives=[
                ("reference artefact (D128-A3.2): the comparison is tumour vs GTEx normal tissue, and the "
                 "target's own expression is "
                 + (f"{float(tgt.logFC.iloc[0]):+.2f} logFC" if len(tgt) else "unavailable")
                 + " in this contrast. A negative value means the normal reference already runs the programme "
                   "and the closure direction can invert for reasons that have nothing to do with the drug. "
                 + ("NOT ruled out." if len(tgt) and float(tgt.logFC.iloc[0]) < 0 else
                    "Target expression is higher in tumour, so this particular inversion is not active, but the "
                    "reference remains bulk normal tissue of different composition.")),
                ("composition: GTEx normal tissue and the tumour differ in cell composition; a programme can "
                 "appear dysregulated because the mixture changed."),
                ("the closure is a cohort-level statement; only "
                 + ("1.07%" if key.startswith("fulvestrant") else "49.79%")
                 + " of patients carry an alteration of the target (L2), and the trace says nothing about which "
                   "patients these dysregulated genes belong to."),
            ]))

    # ------------------------------------------------------------------ verdict: the weakest link (Section 4)
    order = {"holds": 0, "weak": 1, "not evaluable": 2, "breaks": 3}
    verdicts = []
    pre = json.loads((OUT / "links_L1_L4.json").read_text())
    for key in PAIRS:
        ls = [l for l in pre if l["pair"] == key] + [l for l in links if l["pair"] == key]
        worst_level = max(order[l["verdict"]] for l in ls)
        worst = [l for l in ls if order[l["verdict"]] == worst_level]          # forward order preserved
        carried = [dict(link=l["link"], alternative=a) for l in ls for a in l["alternatives"]
                   if "ot ruled out" in a.lower()]
        carried += [dict(link=l["link"], alternative="sub-prediction " + sp["id"] + " FAILED: " + sp["statement"]
                         + " -> " + sp["value"])
                    for l in ls for sp in l.get("subpredictions", []) if not sp["passed"]]
        src = SOURCES[key]
        classes = sorted({c for v in src.values() for c in v})
        # Section 3: a link counts only if it brings a source class no EARLIER link already brought.
        seen, adds_new, adds_nothing = set(), [], []
        for lk in ("L1", "L2", "L3", "L4", "L5", "L6", "L7"):
            fresh = set(src[lk]) - seen
            (adds_new if fresh else adds_nothing).append(lk)
            seen |= set(src[lk])
        discounts = [e[2] for e in ENTANGLED if e[0] in adds_new and e[1] in adds_new]
        verdicts.append(dict(
            pair=key, chain=[dict(link=l["link"], verdict=l["verdict"]) for l in ls],
            break_point=worst[0]["link"], all_links_at_worst=[l["link"] for l in worst],
            verdict=worst[0]["verdict"], why=worst[0]["observed"][:500],
            alternatives_carried=carried,
            independence=dict(
                total_links=len(ls), source_classes=classes,
                links_adding_no_new_source=adds_nothing,
                links_adding_a_new_source=adds_new,
                effective_independent_links=len(adds_new) - len(discounts),
                discounts=discounts),
            summary=summary[key]))
    (OUT / "links_L5_L7.json").write_text(json.dumps(links, indent=1, default=str))
    (OUT / "verdicts.json").write_text(json.dumps(verdicts, indent=1, default=str))
    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    (OUT / "t2_run.json").write_text(json.dumps(dict(
        stage="TRACE-T2 (links L5-L7, verdict)", git_sha_at_run=sha,
        freeze_commit_subject="TRACE FREEZE: predicted footprints for both mechanism-trace pairs (L1-L4)",
        read_d128_a4=False, read_d130_nulls=False, read_outcome_labels=False,
        summary=summary), indent=1, default=str))
    for k, v in summary.items():
        print(k, json.dumps(v, indent=1, default=str))
    for v in verdicts:
        print(v["pair"], "->", v["break_point"], v["verdict"], [c["verdict"] for c in v["chain"]],
              "eff_links", v["independence"]["effective_independent_links"],
              "of", v["independence"]["total_links"])


if __name__ == "__main__":
    sys.exit(main())
