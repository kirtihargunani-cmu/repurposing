#!/usr/bin/env python
"""TRACE-ARACNe stage 2 (= PROCEDURE stage T1): links L1-L4 and the FROZEN predicted footprints,
for the two real pairs AND the five Section 7 negative controls, on cohort-matched ARACNe regulons.

Differences from scripts/trace/t1_predict.py, which ran the same procedure on CollecTRI:
  * the regulon is data/processed/lincs/trace_aracne/aracne_regulon_edges.tsv.gz, subset to the
    ARACNe network of the pair's own cohort (config/group_networks.tsv). CollecTRI is one global
    literature network; ARACNe is a per-cohort co-expression network, so the footprint is
    cohort-specific by construction.
  * tfmode is continuous in (-1, 1) rather than a curated +/-1, and every edge carries a likelihood.
    The headline footprint stays UNWEIGHTED (P7, w=1) so it is comparable with the CollecTRI trace;
    the likelihood is carried per gene for a secondary read, not used in the prediction.
  * the independence story inverts and is written into L4: an ARACNe regulon IS learned from tumour
    co-expression, the thing CollecTRI was not, so L4 and L7 share a data source here (D133.1).

Reads NO drug response profile, NO D128-A2/A4 output, NO null and NO outcome label.
Output: results/lincs/trace_aracne/{predicted_footprint_<pair>.tsv, l3_paths_<pair>.tsv,
        links_L1_L4.json, t1a_run.json}
"""
import json, os, subprocess, sys
from collections import defaultdict
from pathlib import Path

import numpy as np, pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
OUT = ROOT / "results/lincs/trace_aracne"
A1 = ROOT / "results/lincs/d128/a1"
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"
HOP_CAP_NON_TF = 2        # frozen P6, unchanged from the CollecTRI trace

REAL = {
    "fulvestrant_ESR1_BRCA_LumA": dict(
        kind="real", label="TF-anchored", cmap_name="fulvestrant", target="ESR1",
        group_id="BRCA-LumA", ModelID="ACH-000019", line="MCF7",
        control="not a control: the TF-anchored trace, re-run on ARACNe"),
    "vemurafenib_BRAF_TCGA_SKCM": dict(
        kind="real", label="Kinase", cmap_name="vemurafenib", target="BRAF",
        group_id="TCGA-SKCM", ModelID="ACH-000219", line="A375",
        control="not a control: the kinase trace, re-run on ARACNe"),
}


def net_for(cohort, override):
    """Production rule (scripts/79e_drug_reversal.py:125 + D97): config/group_networks.tsv is an
    OVERRIDE map for the subtype cohorts; every other cohort strips its source prefix."""
    if cohort in override:
        return override[cohort]
    return "COADREAD" if cohort in ("TCGA-COAD", "TCGA-READ") else cohort.split("-", 1)[1]


def build_pairs():
    sel = json.loads((OUT / "control_selection.json").read_text())
    ovr = pd.read_csv(ROOT / "config/group_networks.tsv", sep="\t").set_index("group_id").net.to_dict()
    pairs = {}
    for k, v in REAL.items():
        pairs[k] = dict(v, net=net_for(v["group_id"], ovr))
    for name, c in sel["controls"].items():
        pairs[name] = dict(kind="control", label=name, cmap_name=c["cmap_name"],
                           target=c["target"], group_id=c["group_id"], ModelID=c["ModelID"],
                           line=c["line"], net=c["net"], control=c["control"],
                           why=sel["why"][name])
    return pairs


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    PAIRS = build_pairs()
    genes = pd.read_csv(ROOT / "data/processed/lincs/drug_line_profiles_genes.tsv", sep="\t")
    sym = dict(zip(genes.gene_id, genes.symbol))
    fspace = dict(zip(genes.gene_id, genes.feature_space))
    units = pd.read_csv(ROOT / "results/lincs/target_role/d121/units.tsv", sep="\t",
                        dtype={"tid": str, "action_flag": str})
    comp = pd.read_csv(ROOT / "results/lincs/target_role/gene_components.tsv.gz", sep="\t")
    rank = pd.read_csv(ROOT / "results/lincs/fda_ranking_v2_consensus.tsv", sep="\t")
    lw = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    P = pd.read_parquet(A1 / "target_tf_paths.parquet")
    P = P[P.network == "omnipath"]
    E = pd.read_csv(ROOT / "data/processed/lincs/trace_aracne/aracne_regulon_edges.tsv.gz", sep="\t")
    E = E[E.in_lincs_gene_space]
    ENET = {n: d for n, d in E.groupby("network")}

    links, counts = [], {}
    for key, p in PAIRS.items():
        net = p["net"]
        En = ENET[net]
        regs_here = set(En.regulator)
        claims_all = En.groupby("gene_id").regulator.nunique()          # Section 5, within this network
        reg_of_gene = En.groupby("gene_id").regulator.apply(set)

        u = units[(units.cmap_name == p["cmap_name"]) & (units.genes == p["target"])]
        assert len(u) == 1, (key, len(u))
        u = u.iloc[0]
        action = 1 if u.orient == "raise" else -1

        # ----------------------------------------------------------------- L1 drug -> target
        rr = rank[(rank.group_id == p["group_id"]) & (rank.cmap_name == p["cmap_name"])]
        n_units = float(rr.layered_n_units.iloc[0]) if len(rr) else float("nan")
        links.append(dict(
            pair=key, kind=p["kind"], control=p["control"], network=net, link="L1",
            question="Does the drug bind the target, and on what evidence?",
            expected=(f"{p['cmap_name']} binds {p['target']} as a ChEMBL mechanism unit with a measured "
                      f"potency. This link is network-independent: it is the same evidence the CollecTRI "
                      f"trace used, and switching regulon cannot change it."),
            observed=(f"ttype={u.ttype}; genes={u.genes}; k={u.k}; actions={u.actions}; orient={u.orient} "
                      f"(action sign {action:+d}); pChEMBL p_u={u.p_u}; relative binding r_u={u.r_u}; "
                      f"binding_unmeasured={u.binding_unmeasured}; units for this drug in the "
                      f"{p['group_id']} ranking row={n_units:.0f}"),
            verdict="holds" if (u.ttype == "SINGLE PROTEIN" and not u.binding_unmeasured) else "weak",
            provenance="results/lincs/target_role/d121/units.tsv (D121 layer 2); ChEMBL via D105b",
            independence="sole source is the ChEMBL annotation; shared with L2's unit definition and "
                         "with L4's action sign. The trace is NOT independent of this annotation.",
            concentration="1.00 - one unit, one gene, one annotation source.",
            alternatives=[
                "off-target: a single-protein unit does not exclude unannotated off-targets.",
                (f"polypharmacology: at a 10 uM reference the exposure exceeds this unit's measured "
                 f"affinity by ~{10 ** (float(u.p_u) - 5):.0f}x, so selectivity at the profiled "
                 f"concentration is NOT established." if pd.notna(u.p_u) else
                 "polypharmacology: no pChEMBL for this unit, so selectivity is not established at all."),
            ]))

        # ----------------------------------------------------------------- L2 target -> disease
        c = comp[(comp.group_id == p["group_id"]) & (comp.gene == p["target"])]
        if len(c):
            c = c.iloc[0]
            l2_obs = (f"dep_z={c.dep_z} -> D pct {c.D}; alt_share={c.alt_share:.4f} "
                      f"({100 * c.alt_share:.2f}% of patients) -> A pct {c.A}; o_raw={c.o_raw} "
                      f"(source {c.o_source}) -> O pct {c.O}")
            if len(rr):
                r = rr.iloc[0]
                l2_obs += (f"; layered_S={r.layered_S:.2f}, layered_pct={r.layered_pct:.2f}, "
                           f"minority_alteration={r.minority_alteration}")
            else:
                l2_obs += ("; this drug has NO ranking row in this cohort - the pairing does not exist "
                           "in production, which is what a random pairing should look like")
            l2_verdict = "holds" if (c.D > 90 or c.A > 90 or c.O > 90) else "weak"
            o_source, alt_share = c.o_source, float(c.alt_share)
        else:
            l2_obs = (f"{p['target']} has NO component row in {p['group_id']}: the production pipeline "
                      f"never evaluated this target against this cohort.")
            l2_verdict = "breaks"
            o_source, alt_share = "absent", float("nan")
        lrow = lw[(lw.group_id == p["group_id"]) & (lw.ModelID == p["ModelID"])]
        links.append(dict(
            pair=key, kind=p["kind"], control=p["control"], network=net, link="L2",
            question="Is the target altered, over-expressed or depended on in this cohort?",
            expected=(f"{p['target']} is relevant to {p['group_id']} through dependency, alteration share "
                      f"or over-expression." if p["kind"] == "real" else
                      f"THIS IS A NEGATIVE CONTROL. {p['control']}"),
            observed=l2_obs, verdict=l2_verdict,
            provenance=("results/lincs/target_role/gene_components.tsv.gz (D120); "
                        "results/lincs/fda_ranking_v2_consensus.tsv (D121/D122); "
                        "alteration = TCGA pan-cancer; dependency = DepMap 24Q4"),
            independence=(("O is a VIPER NES over a regulon, so L2's O and L4 share the regulon source."
                           if o_source == "NES" else
                           "O is logFC vs GTEx, the same contrast family L7 uses, so L2's O and L7 are "
                           "entangled.") if len(c) else "no component row: nothing to entangle."),
            concentration=(f"largest component's share: {max(c.D, c.A, c.O) / (c.D + c.A + c.O):.3f}"
                           if len(c) else "not applicable"),
            alternatives=([f"only {100 * alt_share:.2f}% of patients carry an alteration."] if len(c) else
                          ["the target is not evaluated in this cohort at all; any downstream number is "
                           "computed on a pairing the pipeline does not make."]) +
                         ["the dependency component is a cell-line average (DepMap), not a patient "
                          "measurement."],
            ))

        # ----------------------------------------------------------------- L3 target -> ARACNe regulators
        pp = P[P.target_gene == p["target"]]
        is_reg = p["target"] in regs_here
        hop_cap = 0 if is_reg else HOP_CAP_NON_TF
        used = pp[(pp.hops <= hop_cap)].copy()
        used["has_aracne_regulon"] = used.tf.isin(regs_here)
        used["aracne_regulon_genes"] = used.tf.map(En.groupby("regulator").size()).fillna(0).astype(int)
        used = used[used.has_aracne_regulon]
        res = used[used.path_status == "resolved"]
        amb = used[used.path_status == "ambiguous"]
        used.sort_values(["hops", "aracne_regulon_genes"], ascending=[True, False]).to_csv(
            OUT / f"l3_paths_{key}.tsv", sep="\t", index=False)
        tot = int(res.aracne_regulon_genes.sum())
        top = res.sort_values("aracne_regulon_genes", ascending=False).head(3)
        n_wider = int(((pp.hops <= 3) & pp.tf.isin(regs_here)).sum())
        links.append(dict(
            pair=key, kind=p["kind"], control=p["control"], network=net, link="L3",
            question="What signed path reaches which regulators, in how many hops?",
            expected=(("zero-hop TF anchor: the target is itself an ARACNe regulator in this cohort"
                       if is_reg else
                       f"signalling: the target is not a regulator, so it must reach regulators within "
                       f"{HOP_CAP_NON_TF} hops") if p["kind"] == "real" else
                      f"THIS IS A NEGATIVE CONTROL. {p['control']}"),
            observed=(f"network {net}; target is an ARACNe regulator there: {is_reg}; hop cap {hop_cap}. "
                      f"{len(pp)} OmniPath route rows for this target; {len(res)} resolved-sign regulators "
                      f"with an ARACNe regulon, {len(amb)} sign-ambiguous (excluded). "
                      + ("Resolved regulators and regulon sizes: "
                         + "; ".join(f"{t.tf}(h{t.hops},{t.path_sign:+d},{t.aracne_regulon_genes}g)"
                                     for t in res.sort_values(['hops', 'aracne_regulon_genes'],
                                                              ascending=[True, False]).head(25).itertuples())
                         + (f" ... and {len(res) - 25} more" if len(res) > 25 else "")
                         if len(res) else "NO regulator is reached.")),
            verdict="holds" if len(res) >= 1 else "not evaluable",
            provenance=(f"results/lincs/d128/a1/target_tf_paths.parquet (OmniPath D126 cache) for the "
                        f"routes; data/processed/lincs/trace_aracne/aracne_regulon_edges.tsv.gz "
                        f"(network {net}) for which of them has a regulon"),
            independence=("L3 and L4 both derive from the ARACNe regulon: a regulon error moves both. "
                          "For a zero-hop anchor L3 carries no independent signalling evidence at all."
                          if hop_cap == 0 else
                          "L3 adds the OmniPath signed prior, a source L4 does not use; but every "
                          "regulator it names is scored through its ARACNe regulon, so L3 and L4 still "
                          "share that regulon."),
            concentration=(f"largest regulator's share of total regulon gene slots: "
                           f"{(res.aracne_regulon_genes.max() / tot) if tot else float('nan'):.3f} "
                           f"({top.iloc[0].tf}: {int(top.iloc[0].aracne_regulon_genes)} of {tot} slots "
                           f"over {len(res)} regulators)" if len(res) else "no regulators reached"),
            alternatives=[
                f"{len(amb)} regulators are reached by sign-ambiguous routes and are silently dropped.",
                (f"hop cap is frozen at {HOP_CAP_NON_TF}. At 3 hops this target reaches {n_wider} "
                 f"regulators with an ARACNe regulon instead of {len(used)}."),
                (f"the ARACNe network is {net}, built from that cohort's tumours; which regulators exist "
                 f"at all is a property of that cohort's co-expression, not of the signalling biology."),
            ]))

        # ----------------------------------------------------------------- L4 regulators -> genes
        acc_sign, acc_lik, acc_tf = defaultdict(list), defaultdict(list), defaultdict(set)
        for t in res.itertuples():
            sub = En[En.regulator == t.tf]
            for g, m, lk in zip(sub.gene_id, sub.tfmode, sub.likelihood):
                acc_sign[g].append(int(t.path_sign) * (1 if m >= 0 else -1))
                acc_lik[g].append(float(lk))
                acc_tf[g].add(t.tf)
        rows = []
        for g, signs in acc_sign.items():
            s = set(signs)
            rows.append(dict(gene_id=g, symbol=sym.get(g, g), feature_space=fspace.get(g, "not measured"),
                             expected_sign_pre_action=(signs[0] if len(s) == 1 else 0),
                             n_unit_regulators_claiming=len(acc_tf[g]),
                             mean_likelihood=float(np.mean(acc_lik[g]))))
        fp = pd.DataFrame(rows, columns=["gene_id", "symbol", "feature_space",
                                         "expected_sign_pre_action", "n_unit_regulators_claiming",
                                         "mean_likelihood"])
        if len(fp):
            fp["expected_sign"] = fp.expected_sign_pre_action * action
            fp["n_regulators_claiming"] = fp.gene_id.map(claims_all).fillna(0).astype(int)
            fp["uniquely_attributable"] = fp.n_regulators_claiming == 1
            fp["unit_id"] = f"{p['cmap_name']}::{u.tid}"
            fp["cmap_name"] = p["cmap_name"]
            fp = fp.sort_values(["expected_sign", "symbol"])
        signed = fp[fp.expected_sign != 0] if len(fp) else fp
        uniq_frac = float(signed.uniquely_attributable.mean()) if len(signed) else float("nan")
        hub_counts = {t: sum(1 for g in signed.gene_id if t in acc_tf[g]) for t in res.tf} if len(signed) else {}
        hub_tf = max(hub_counts, key=hub_counts.get) if hub_counts else ""
        hub_share = (hub_counts[hub_tf] / len(signed)) if hub_counts and len(signed) else float("nan")
        fp.to_csv(OUT / f"predicted_footprint_{key}.tsv", sep="\t", index=False)
        counts[key] = dict(network=net, footprint_rows=int(len(fp)), signed=int(len(signed)),
                           up=int((signed.expected_sign > 0).sum()) if len(signed) else 0,
                           down=int((signed.expected_sign < 0).sum()) if len(signed) else 0,
                           conflicted=int((fp.expected_sign == 0).sum()) if len(fp) else 0,
                           n_regulators=int(len(res)),
                           unique_fraction=None if np.isnan(uniq_frac) else round(uniq_frac, 4),
                           structural_hub_tf=hub_tf,
                           structural_hub_share=None if hub_share != hub_share else round(hub_share, 4),
                           median_regulators_claiming=(float(signed.n_regulators_claiming.median())
                                                       if len(signed) else None))
        links.append(dict(
            pair=key, kind=p["kind"], control=p["control"], network=net, link="L4",
            question="Which genes, with what weight and expected direction?",
            expected=((f"the {len(res)} regulators of L3 imply a signed gene footprint; each gene's "
                       f"expected direction is action ({u.orient}, {action:+d}) x path sign x sign(tfmode). "
                       f"Weight unweighted (P7, w=1), as in the CollecTRI trace.")
                      if len(res) else "L3 reached no regulator, so there is no footprint to predict."),
            observed=(f"{len(fp)} footprint genes in the LINCS gene space; {len(signed)} carry a signed "
                      f"expectation ({counts[key]['up']} UP, {counts[key]['down']} DOWN), "
                      f"{counts[key]['conflicted']} conflicted between this unit's own regulators and "
                      f"excluded. Attribution: median "
                      f"{counts[key]['median_regulators_claiming']} regulators in the whole {net} network "
                      f"claim each signed gene; {100 * uniq_frac:.1f}% of signed genes are claimed by "
                      f"exactly one. Structural hub {hub_tf} claims {hub_share:.3f} of the signed genes "
                      f"over {len(res)} contributing regulators."
                      if len(signed) else "empty footprint."),
            verdict=("not evaluable" if not len(signed) else "weak" if uniq_frac < 0.5 else "holds"),
            provenance=(f"data/processed/lincs/trace_aracne/aracne_regulon_edges.tsv.gz, network {net} "
                        f"(viper::aracne2regulon on the cohort's TCGA tumour matrix, "
                        f"scripts/44_aracne_regulons.R); action sign from D121 units.tsv"),
            independence=("same regulon as L3 - one source, reported twice. AND, unlike the CollecTRI "
                          "trace: an ARACNe regulon IS learned from tumour co-expression, which is the "
                          "source L7's tumour contrast uses and the source L2's O component uses when it "
                          "is a NES. The D133.1 entanglement that did NOT apply to CollecTRI DOES apply "
                          "here - L4 and L7 are not independent."),
            concentration=(f"structural: largest single regulator ({hub_tf}) claims {hub_share:.3f} of the "
                           f"{len(signed)} signed genes, over {len(res)} contributing regulators."
                           if len(signed) else "no footprint"),
            alternatives=[
                (f"attribution: {100 * (1 - uniq_frac):.1f}% of signed genes are claimed by other "
                 f"regulators in {net} too." if len(signed) else "no footprint to attribute."),
                "conflicted genes are dropped rather than counted against the prediction.",
                ("ARACNe tfmode is the SIGN OF A CORRELATION across that cohort's tumours, not a "
                 "measured direction of regulation. A gene expected DOWN is a gene that co-varies "
                 "positively with the regulator in tumours; a drug that moves the cell along the same "
                 "axis for any reason will satisfy the prediction."),
                (f"the regulon is {net}-specific: the same target in a different cohort would get a "
                 f"different footprint from the same signalling route."),
            ]))

    (OUT / "links_L1_L4.json").write_text(json.dumps(links, indent=1, default=str))
    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    (OUT / "t1a_run.json").write_text(json.dumps(dict(
        stage="TRACE-ARACNe T1 (links L1-L4, frozen prediction, 2 real + 5 negative controls)",
        git_sha_at_run=sha, procedure="docs/PROCEDURE_mechanism_trace.md",
        regulon="data/processed/lincs/trace_aracne/aracne_regulon_edges.tsv.gz",
        hop_cap_non_tf=HOP_CAP_NON_TF, weighting="unweighted (P7, w=1); likelihood carried not used",
        reads_no_response_data=True, reads_no_outcome_label=True,
        pairs=PAIRS, counts=counts), indent=1, default=str))
    for k, v in counts.items():
        print(f"{k:38s} {v['network']:10s} regs={v['n_regulators']:3d} signed={v['signed']:5d} "
              f"(up {v['up']}, down {v['down']}, conflicted {v['conflicted']}) "
              f"hub={v['structural_hub_tf']} {v['structural_hub_share']}")
    print("wrote", OUT)


if __name__ == "__main__":
    sys.exit(main())
