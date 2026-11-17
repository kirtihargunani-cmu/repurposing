#!/usr/bin/env python
"""TRACE-ARACNe stage 3 (= PROCEDURE stage T2): links L5, L6, L7, the verdicts, and the
procedure's FALSE-POSITIVE RATE over the five Section 7 negative controls.

Written and run only after scripts/trace/aracne/t1a_predict.py and its predicted footprints
were committed (commit "TRACE-ARACNe FREEZE"). Every `expected` string is derived
mechanically from the frozen L1-L4 record, not written after seeing the observation.

PRE-STATED before any number below was computed:
  * "the chain closes" == no link in L1..L7 has verdict "breaks" or "not evaluable".
  * the procedure's false-positive rate is the fraction of the FIVE negative controls on
    which the chain closes. It is reported whatever it turns out to be.

NO GENE SUBSETTING (PROCEDURE Section 5 as amended 2026-11-17): the whole frozen footprint
is used. Landmark-vs-inferred and unique-vs-shared splits are REPORTED as description only
and no verdict depends on them - the transcriptional response is holistic.

NOT read: results/lincs/d128/a4/*, the D130 corrected-null tables, and any outcome label
(CTRP / PRISM / GDSC).
"""
import json, os, subprocess, sys
from pathlib import Path

import numpy as np, pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
OUT = ROOT / "results/lincs/trace_aracne"
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"
F_DIVISOR = 2.0            # P9, frozen in scripts/d128/config.py
CONC_HOLDS, GENERIC_HOLDS = 0.60, 90.0     # same L5 thresholds the CollecTRI trace used
REV_HOLDS = 0.60                            # same L7 threshold


def participation_ratio(x):
    x = np.abs(np.asarray(x, dtype=float))
    s2 = float((x ** 2).sum())
    return float(s2 and (float(x.sum()) ** 2) / s2)


def main():
    run = json.loads((OUT / "t1a_run.json").read_text())
    PAIRS, COUNTS = run["pairs"], run["counts"]
    pre = json.loads((OUT / "links_L1_L4.json").read_text())

    genes_tab = pd.read_csv(ROOT / "data/processed/lincs/drug_line_profiles_genes.tsv", sep="\t")
    idx = pd.read_csv(ROOT / "data/processed/lincs/drug_line_profiles_index.tsv", sep="\t")
    annot = pd.read_csv(ROOT / "results/lincs/drug_annotations.tsv", sep="\t")
    prolif = pd.read_csv(ROOT / "results/lincs/prolif_regulators_p2_v2_consensus.tsv", sep="\t")
    E = pd.read_csv(ROOT / "data/processed/lincs/trace_aracne/aracne_regulon_edges.tsv.gz", sep="\t")
    E = E[E.in_lincs_gene_space]
    M = np.memmap(ROOT / "data/processed/lincs/drug_line_profiles_bing_specific.f32", dtype="float32",
                  mode="r", shape=(len(idx), len(genes_tab)))
    gpos = {g: i for i, g in enumerate(genes_tab.gene_id)}

    links, summary = [], {}
    for key, p in PAIRS.items():
        net = p["net"]
        En = E[E.network == net]
        members = En.groupby("regulator").gene_id.apply(set).to_dict()
        pro = prolif[(prolif.net == net) & (prolif.set == "aracne")].set_index("regulator")
        fp = pd.read_csv(OUT / f"predicted_footprint_{key}.tsv", sep="\t")
        l3 = pd.read_csv(OUT / f"l3_paths_{key}.tsv", sep="\t")
        own_tfs = sorted(l3[l3.path_status == "resolved"].tf) if len(l3) else []
        signed = fp[fp.expected_sign != 0].copy() if len(fp) else fp

        # ------------------------------------------------------------------ L5
        if not len(signed):
            links.append(dict(
                pair=key, kind=p["kind"], control=p["control"], network=net, link="L5",
                question="Did those genes move as predicted, in which cell line, matched how?",
                expected="L3 reached no regulator and L4 produced no footprint, so there is nothing to test.",
                observed="no frozen footprint; L5 cannot be evaluated.",
                verdict="not evaluable",
                provenance="results/lincs/trace_aracne/predicted_footprint_%s.tsv (empty)" % key,
                independence="not applicable", concentration="not applicable",
                alternatives=["nothing was predicted, so nothing is explained."]))
            summary[key] = dict(network=net, n_signed=0, concordance=None, mean_signed_fz=None,
                                generic_percentile=None, reversal=None)
            for lk in ("L6", "L7"):
                links.append(dict(
                    pair=key, kind=p["kind"], control=p["control"], network=net, link=lk,
                    question="(downstream of an empty footprint)",
                    expected="nothing to concentrate / nothing to close on.",
                    observed="not evaluable: the footprint is empty.", verdict="not evaluable",
                    provenance="-", independence="-", concentration="-",
                    alternatives=["-"]))
            continue

        pid = sorted(set(annot[annot.cmap_name == p["cmap_name"]].pert_id))
        rows = idx[(idx.ModelID == p["ModelID"]) & (idx.pert_id.isin(pid))].sort_values(
            "n_sig", ascending=False)
        assert len(rows) >= 1, (key, "no profile row")
        prow = int(rows.index[0])
        cov = signed.copy()
        cov["z"] = [float(M[prow, gpos[g]]) if g in gpos else np.nan for g in cov.gene_id]
        cov = cov[np.isfinite(cov.z)].copy()
        cov["f_z"] = np.tanh(cov.z / F_DIVISOR)
        cov["contribution"] = cov.expected_sign * cov.f_z
        cov["concordant"] = np.sign(cov.z) == cov.expected_sign
        conc = float(cov.concordant.mean())
        this_mean = float(cov.contribution.mean())

        # generic-response and amplitude controls over every compound profiled in this line
        line_rows = idx.index[idx.ModelID == p["ModelID"]].to_numpy()
        cols = np.array([gpos[g] for g in cov.gene_id])
        e = cov.expected_sign.to_numpy(float)
        Z = np.asarray(M[line_rows][:, cols], dtype=float)
        means = (np.tanh(Z / F_DIVISOR) * e).mean(axis=1)
        this_i = int(np.where(line_rows == prow)[0][0])
        pct_generic = float((means < this_mean).mean() * 100)
        conc_all = (np.sign(Z) == e).mean(axis=1)
        pct_conc = float((conc_all < conc_all[this_i]).mean() * 100)
        ampl = np.abs(np.asarray(M[line_rows], dtype=float)).mean(axis=1)
        pct_ampl = float((ampl < ampl[this_i]).mean() * 100)
        r_ampl = float(np.corrcoef(ampl, means)[0, 1])
        by_space = cov.groupby("feature_space").concordant.agg(["mean", "size"])

        links.append(dict(
            pair=key, kind=p["kind"], control=p["control"], network=net, link="L5",
            question="Did those genes move as predicted, in which cell line, matched how?",
            expected=(f"The {COUNTS[key]['signed']} frozen genes ({COUNTS[key]['up']} expected UP, "
                      f"{COUNTS[key]['down']} expected DOWN, from {COUNTS[key]['n_regulators']} "
                      f"{net} regulators) should move in the predicted direction in {p['line']} more "
                      f"often than not, and the whole-footprint mean should sit high against the other "
                      f"compounds profiled in that line. "
                      + ("If this is a real mechanism the chain should hold here."
                         if p["kind"] == "real" else
                         f"THIS IS A NEGATIVE CONTROL - {p['control']}")),
            observed=(f"line {p['line']} ({p['ModelID']}), pert_id {idx.loc[prow].pert_id}, "
                      f"{int(idx.loc[prow].n_sig)} signature(s), drug-specific profile arm. "
                      f"{len(cov)} of {len(signed)} frozen genes have a response value. "
                      f"Whole-footprint sign concordance {100 * conc:.1f}%. "
                      f"Mean expected-signed f(z) = {this_mean:+.4f}; against the {len(line_rows)} "
                      f"compounds profiled in this line that is the {pct_generic:.1f}th percentile "
                      f"(concordance the {pct_conc:.1f}th). "
                      f"DESCRIPTIVE ONLY, no verdict rests on it - by feature space: "
                      + ", ".join(f"{k} {100 * v['mean']:.1f}% (n={int(v['size'])})"
                                  for k, v in by_space.iterrows()) + "."),
            verdict=("holds" if conc >= CONC_HOLDS and pct_generic >= GENERIC_HOLDS else
                     "weak" if conc > 0.50 else "breaks"),
            provenance=("data/processed/lincs/drug_line_profiles_bing_specific.f32 (drug-specific arm); "
                        "frozen footprint from results/lincs/trace_aracne/"),
            independence=("L5 shares the footprint with L3/L4, so a regulon error moves the prediction "
                          "and the test of it together. The ARACNe regulon is learned from the cohort's "
                          "TUMOUR co-expression while the response is a CELL LINE perturbation, so L5 "
                          "is not self-confirming, but the expected direction still encodes a tumour "
                          "co-expression axis rather than a measured regulatory direction."),
            concentration=(f"largest single gene's share of |contribution|: "
                           f"{float(cov.contribution.abs().max() / cov.contribution.abs().sum()):.4f}; "
                           f"effective number of contributing genes "
                           f"{participation_ratio(cov.contribution):.0f} of {len(cov)}."),
            alternatives=[
                (f"generic response: the same frozen footprint sits at the {pct_generic:.1f}th percentile "
                 f"of all {len(line_rows)} compounds profiled in {p['line']}. "
                 + ("Not ruled out - many other compounds move these genes the same way."
                    if pct_generic < GENERIC_HOLDS else
                    "Largely ruled out at this link.")),
                (f"amplitude: this drug-line's mean |z| is at the {pct_ampl:.1f}th percentile of the line; "
                 f"across the line's compounds the footprint mean correlates r={r_ampl:+.2f} with overall "
                 f"response magnitude. "
                 + ("Not ruled out: the footprint mean largely tracks how big the response is."
                    if abs(r_ampl) >= 0.5 else "Weakly supported as a confound at this link.")),
                "off-target: concordance cannot separate the annotated target from any off-target that "
                "converges on the same regulators.",
            ]))

        # ------------------------------------------------------------------ L6
        shares, prs = {}, {}
        for t in own_tfs:
            mem = members.get(t, set())
            sub = cov[cov.gene_id.isin(mem)]
            shares[t] = (float(sub.contribution.abs().sum() / cov.contribution.abs().sum())
                         if len(sub) else 0.0)
            prs[t] = len(sub)
        hub = max(shares, key=shares.get) if shares else ""
        hub_share = shares.get(hub, float("nan"))
        n_eff_reg = participation_ratio(list(shares.values())) if shares else float("nan")
        hub_pro = (float(pro.r_prolif.get(hub, np.nan)), bool(pro.prolif_linked.get(hub, False)))
        pro_linked = [t for t in own_tfs if bool(pro.prolif_linked.get(t, False))]
        pro_mass = float(sum(shares.get(t, 0.0) for t in pro_linked))
        links.append(dict(
            pair=key, kind=p["kind"], control=p["control"], network=net, link="L6",
            question="Is the evidence concentrated in one regulator, and is that regulator a "
                     "proliferation proxy?",
            expected=(f"L4's structure put the largest regulator at "
                      f"{COUNTS[key]['structural_hub_share']} of the predicted genes over "
                      f"{COUNTS[key]['n_regulators']} regulators. If the observed signal is distributed "
                      f"like the structure, roughly that share of observed |contribution| should sit "
                      f"with one regulator. A footprint that is really one regulon is one piece of "
                      f"evidence, not {COUNTS[key]['signed']}."),
            observed=(f"{len(own_tfs)} contributing regulators. Largest observed share of "
                      f"|contribution|: {hub} at {hub_share:.3f} (shares are NOT disjoint - a gene in "
                      f"two regulons counts in both). Effective number of regulators (participation "
                      f"ratio over the shares) {n_eff_reg:.2f} of {len(own_tfs)}. "
                      f"{hub} proliferation correlation in the {net} ARACNe set: r_prolif="
                      f"{hub_pro[0]:.3f}, prolif_linked={hub_pro[1]}. "
                      f"{len(pro_linked)} of {len(own_tfs)} contributing regulators are "
                      f"proliferation-linked, carrying {pro_mass:.3f} of the non-disjoint share mass. "
                      f"Top regulators: "
                      + "; ".join(f"{t}({shares[t]:.3f},{prs[t]}g)"
                                  for t in sorted(shares, key=shares.get, reverse=True)[:8])),
            verdict=("breaks" if hub_share >= 0.50 or n_eff_reg < 2 else
                     "weak" if hub_share >= 0.30 or pro_mass >= 0.50 else "holds"),
            provenance=(f"data/processed/lincs/trace_aracne/aracne_regulon_edges.tsv.gz (network {net}); "
                        f"results/lincs/prolif_regulators_p2_v2_consensus.tsv ({net} / aracne)"),
            independence=("this link is an AUDIT of L3-L5, not a new source: it shares the regulon with "
                          "L3/L4 and the response with L5. The proliferation annotation is the one "
                          "outside input, and it shares the ARACNe regulon with L4."),
            concentration=f"{hub} {hub_share:.3f}; effective regulators {n_eff_reg:.2f} of {len(own_tfs)}",
            alternatives=[
                ((f"proliferation: {hub} is proliferation-linked in the {net} ARACNe set "
                  f"(r_prolif={hub_pro[0]:.3f}), so the link cannot separate the drug's mechanism from "
                  f"growth arrest. NOT ruled out.") if hub_pro[1] else
                 (f"proliferation: {hub} is not flagged proliferation-linked in {net} "
                  f"(r_prolif={hub_pro[0]:.3f}); but {len(pro_linked)} other contributing regulators are, "
                  f"carrying {pro_mass:.3f} of the share mass."
                  + (" NOT ruled out." if pro_mass >= 0.3 else ""))),
                "the regulator shares are not disjoint: a gene claimed by two of the unit's regulators "
                "is counted in both, so these are upper bounds per regulator.",
                "D130/D133: within a cohort the regulator activities behave as ~10 independent terms and "
                "a ~1000-gene footprint as ~44; a gene count is not an evidence count.",
            ]))

        # ------------------------------------------------------------------ L7
        con = ROOT / f"data/processed/step4_v2/{p['group_id']}__gtex.tsv.gz"
        DE = pd.read_csv(con, sep="\t")[["gene_id", "symbol", "logFC", "z", "adj.P.Val"]]
        j = cov.merge(DE, on="gene_id", how="inner", suffixes=("", "_de"))
        sig = j[j["adj.P.Val"] < 0.05].copy()
        rev = float(sig.reverses.mean()) if False else None
        sig["reverses"] = np.sign(sig.logFC) == -sig.expected_sign
        rev = float(sig.reverses.mean()) if len(sig) else float("nan")
        wrev = (float((sig.reverses * sig.z.abs()).sum() / sig.z.abs().sum()) if len(sig) else float("nan"))
        tgt = DE[DE.symbol == p["target"]]
        tlog = float(tgt.logFC.iloc[0]) if len(tgt) else float("nan")
        links.append(dict(
            pair=key, kind=p["kind"], control=p["control"], network=net, link="L7",
            question="Are those same genes dysregulated in the tumour in the direction the drug reverses?",
            expected=(f"If the chain closes, the genes the drug is predicted to push DOWN should be UP in "
                      f"{p['group_id']} vs its normal reference and vice versa, for a clear majority of "
                      f"the significant footprint genes. "
                      + ("" if p["kind"] == "real" else f"THIS IS A NEGATIVE CONTROL - {p['control']}")),
            observed=(f"contrast {p['group_id']} vs gtex. {len(j)} frozen+covered genes in the contrast, "
                      f"{len(sig)} at adj.P.Val < 0.05. The drug's predicted direction opposes the "
                      f"tumour's for {100 * rev:.1f}% of them ({100 * wrev:.1f}% weighted by |z|). "
                      f"The target itself: "
                      + (f"{p['target']} logFC {tlog:+.2f}, z {float(tgt.z.iloc[0]):+.2f}, "
                         f"adj.P {float(tgt['adj.P.Val'].iloc[0]):.2e}" if len(tgt)
                         else "not in the contrast table")),
            verdict=("not evaluable" if not len(sig) else
                     "holds" if rev >= REV_HOLDS else "weak" if rev > 0.50 else "breaks"),
            provenance=f"data/processed/step4_v2/{p['group_id']}__gtex.tsv.gz (limma + voomByGroup; "
                       f"P1 logFC sign, z weight; P8 adj.P.Val < 0.05)",
            independence=("the tumour contrast is a source L5 does not use. BUT the ARACNe regulon that "
                          "selected these genes was itself learned from this cohort's tumour expression, "
                          "so L4 and L7 share a data source (D133.1) - unlike the CollecTRI trace, where "
                          "the regulon was literature-curated. L7 is NOT an independent closure here."),
            concentration=(f"largest single gene's share of the |z| weight: "
                           f"{float(sig.z.abs().max() / sig.z.abs().sum()):.4f}; effective number of "
                           f"genes {participation_ratio(sig.z):.0f} of {len(sig)}" if len(sig) else "-"),
            alternatives=[
                ("reference artefact (D128-A3.2): the comparison is tumour vs GTEx normal tissue and the "
                 f"target's own logFC is {tlog:+.2f}. "
                 + ("NOT ruled out - the normal reference already runs the programme." if tlog < 0 else
                    "Target expression is higher in tumour, so this inversion is not active, but the "
                    "reference remains bulk normal tissue of different composition.")),
                "composition: GTEx normal and tumour differ in cell composition; a programme can look "
                "dysregulated because the mixture changed.",
                "ENTANGLEMENT: the gene set asked about here was selected by a regulon learned from this "
                "same tumour cohort's co-expression, so agreement is partly guaranteed by construction.",
            ]))
        summary[key] = dict(network=net, n_signed=int(len(signed)), n_covered=int(len(cov)),
                            concordance=round(conc, 4), mean_signed_fz=round(this_mean, 6),
                            generic_percentile=round(pct_generic, 1),
                            amplitude_percentile=round(pct_ampl, 1), r_amplitude=round(r_ampl, 3),
                            n_regulators=len(own_tfs), hub=hub, hub_share=round(hub_share, 4),
                            n_eff_regulators=round(n_eff_reg, 2),
                            prolif_share=round(pro_mass, 4),
                            reversal=None if rev != rev else round(rev, 4))

    # ------------------------------------------------------------------ verdicts + false-positive rate
    order = {"holds": 0, "weak": 1, "not evaluable": 2, "breaks": 3}
    verdicts = []
    for key, p in PAIRS.items():
        ls = [l for l in pre if l["pair"] == key] + [l for l in links if l["pair"] == key]
        worst_level = max(order[l["verdict"]] for l in ls)
        worst = [l for l in ls if order[l["verdict"]] == worst_level]
        carried = [dict(link=l["link"], alternative=a) for l in ls for a in l["alternatives"]
                   if "ot ruled out" in a.lower()]
        closed = not any(l["verdict"] in ("breaks", "not evaluable") for l in ls)
        verdicts.append(dict(
            pair=key, kind=p["kind"], control=p["control"], network=p["net"],
            chain=[dict(link=l["link"], verdict=l["verdict"]) for l in ls],
            break_point=worst[0]["link"], all_links_at_worst=[l["link"] for l in worst],
            verdict=worst[0]["verdict"], chain_closes=closed,
            why=worst[0]["observed"][:400], alternatives_carried=carried,
            summary=summary.get(key)))

    controls = [v for v in verdicts if v["kind"] == "control"]
    n_closed = sum(v["chain_closes"] for v in controls)
    fpr = dict(definition="a negative control 'closes' when no link in L1..L7 has verdict 'breaks' or "
                          "'not evaluable'; stated before the numbers were computed",
               n_controls=len(controls), n_closed=int(n_closed),
               false_positive_rate=round(n_closed / len(controls), 3),
               per_control={v["pair"]: dict(closes=v["chain_closes"], break_point=v["break_point"],
                                            verdict=v["verdict"]) for v in controls})

    (OUT / "links_L5_L7.json").write_text(json.dumps(links, indent=1, default=str))
    (OUT / "verdicts.json").write_text(json.dumps(verdicts, indent=1, default=str))
    (OUT / "false_positive_rate.json").write_text(json.dumps(fpr, indent=1, default=str))
    recs = pd.DataFrame(pre + links)[["pair", "kind", "network", "link", "question", "expected",
                                      "observed", "verdict", "provenance", "independence",
                                      "concentration", "alternatives"]]
    recs.to_csv(OUT / "trace_records.tsv", sep="\t", index=False)
    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    (OUT / "t2a_run.json").write_text(json.dumps(dict(
        stage="TRACE-ARACNe T2 (L5-L7, verdicts, false-positive rate)", git_sha_at_run=sha,
        freeze_commit_subject="TRACE-ARACNe FREEZE: L1-L4 predicted footprints, 2 real pairs + 5 negative controls",
        gene_subsetting_used=False, read_d128_a4=False, read_outcome_labels=False,
        thresholds=dict(L5_concordance_holds=CONC_HOLDS, L5_generic_percentile_holds=GENERIC_HOLDS,
                        L7_reversal_holds=REV_HOLDS),
        summary=summary, false_positive_rate=fpr), indent=1, default=str))

    print(f"\n{'pair':38s} {'net':10s} {'break':6s} {'verdict':14s} closes  chain")
    for v in verdicts:
        print(f"{v['pair']:38s} {v['network']:10s} {v['break_point']:6s} {v['verdict']:14s} "
              f"{str(v['chain_closes']):6s} {[c['verdict'][:4] for c in v['chain']]}")
    print(f"\nFALSE-POSITIVE RATE over {fpr['n_controls']} negative controls: "
          f"{fpr['n_closed']}/{fpr['n_controls']} = {fpr['false_positive_rate']}")
    for k, v in summary.items():
        print(" ", k, json.dumps(v, default=str))


if __name__ == "__main__":
    sys.exit(main())
