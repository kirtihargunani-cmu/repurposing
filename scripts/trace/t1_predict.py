#!/usr/bin/env python
"""Mechanism trace, stage 1: links L1-L4 and the FROZEN predicted footprint.

Procedure: docs/PROCEDURE_mechanism_trace.md (Sections 1-6, 8 steps 1-4, 10).

This stage reads ONLY drug annotation, disease-side components and the drug-agnostic
A1 route index. It reads NO drug response profile, NO D128-A2 contribution, NO null and
NO outcome label. Its output is committed BEFORE stage 2 (t2_observe.py) is written or run,
so the prediction is on record before the observation (Section 8 step 4).

Outputs results/lincs/trace/:
  predicted_footprint_<pair>.tsv   per-gene expected direction + attribution ambiguity
  links_L1_L4.json                 the seven-field records for L1-L4
  l3_paths_<pair>.tsv              the signed routes target -> regulator
  t1_run.json                      manifest
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
OUT = ROOT / "results/lincs/trace"
A1 = ROOT / "results/lincs/d128/a1"
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"

PAIRS = {
    "fulvestrant_ESR1_BRCA_LumA": dict(
        label="TF-anchored", cmap_name="fulvestrant", unit_id="fulvestrant::CHEMBL206",
        target="ESR1", group_id="BRCA-LumA", ModelID="ACH-000019", line="MCF7",
        expected_hops=0, expected_route="zero-hop TF anchor: ESR1 is itself a CollecTRI regulator"),
    "vemurafenib_BRAF_TCGA_SKCM": dict(
        label="Kinase", cmap_name="vemurafenib", unit_id="vemurafenib::CHEMBL5145",
        target="BRAF", group_id="TCGA-SKCM", ModelID="ACH-000219", line="A375",
        expected_hops=2, expected_route="signalling: BRAF is not a regulator, reaches TFs at 1-2 hops"),
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(ROOT / "data/processed/lincs/drug_line_profiles_genes.tsv", sep="\t")
    sym = dict(zip(genes.gene_id, genes.symbol))
    units = pd.read_csv(ROOT / "results/lincs/target_role/d121/units.tsv", sep="\t",
                        dtype={"tid": str, "action_flag": str})
    comp = pd.read_csv(ROOT / "results/lincs/target_role/gene_components.tsv.gz", sep="\t")
    rank = pd.read_csv(ROOT / "results/lincs/fda_ranking_v2_consensus.tsv", sep="\t")
    lw = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    P = pd.read_parquet(A1 / "target_tf_paths.parquet")
    F = pd.read_parquet(A1 / "target_gene_footprints.parquet")
    E = pd.read_parquet(A1 / "regulon_edges.parquet")

    # how many DISTINCT CollecTRI regulators claim each gene, over every regulator in this analysis (Section 5)
    claims = E[E.in_lincs_gene_space].groupby("gene_id").regulator.nunique()
    reg_of_gene = E[E.in_lincs_gene_space].groupby("gene_id").regulator.apply(set)

    links, manifest_counts = [], {}
    for key, p in PAIRS.items():
        u = units[(units.cmap_name == p["cmap_name"]) & (units.genes == p["target"])]
        assert len(u) == 1, (key, len(u))
        u = u.iloc[0]
        action = 1 if u.orient == "raise" else -1

        # ---------------- L1 drug -> target --------------------------------------------------------------
        links.append(dict(
            pair=key, link="L1", question="Does the drug bind the target, and on what evidence?",
            expected=(f"{p['cmap_name']} binds {p['target']} as a single-protein ChEMBL mechanism target "
                      f"with a direct, measured potency, and it is the drug's best-binding unit (r_u = 1)."),
            observed=(f"ttype={u.ttype}; genes={u.genes}; k={u.k}; actions={u.actions}; orient={u.orient} "
                      f"(action sign {action:+d}); pChEMBL p_u={u.p_u}; relative binding r_u={u.r_u}; "
                      f"binding_unmeasured={u.binding_unmeasured}; units for this drug in the ranking row="
                      f"{float(rank[(rank.group_id == p['group_id']) & (rank.cmap_name == p['cmap_name'])].layered_n_units.iloc[0]):.0f}"),
            verdict="holds",
            provenance="results/lincs/target_role/d121/units.tsv (D121 layer 2); ChEMBL mechanism records via D105b",
            independence=("sole source is the ChEMBL annotation; shared with L2's unit definition and with L4's "
                          "action sign. Section 9: the trace is NOT independent of this annotation."),
            concentration=("1.00 - one unit, one gene, one annotation source. A single ChEMBL curation error "
                           "moves the whole chain."),
            alternatives=[
                "off-target: a single-protein unit does not exclude unannotated off-targets; "
                "ChEMBL records only what was assayed.",
                ("polypharmacology: at a 10 uM reference concentration the exposure exceeds this unit's "
                 f"measured affinity (10^-{u.p_u} M) by ~{10 ** (float(u.p_u) - 5):.0f}x, so selectivity at the "
                 "concentration LINCS profiles were run at is NOT established by this link. The per-signature "
                 "dose is not read here."),
            ]))

        # ---------------- L2 target -> disease ------------------------------------------------------------
        c = comp[(comp.group_id == p["group_id"]) & (comp.gene == p["target"])]
        assert len(c) == 1
        c = c.iloc[0]
        r = rank[(rank.group_id == p["group_id"]) & (rank.cmap_name == p["cmap_name"])].iloc[0]
        lrow = lw[(lw.group_id == p["group_id"]) & (lw.ModelID == p["ModelID"])]
        links.append(dict(
            pair=key, link="L2",
            question="Is the target altered, over-expressed or depended on in this cohort, in what fraction of patients?",
            expected=(f"{p['target']} is relevant to {p['group_id']} through at least one of dependency (D), "
                      f"alteration share (A) or over-expression (O)."),
            observed=(f"dep_z={c.dep_z} -> D pct {c.D}; alt_share={c.alt_share:.4f} "
                      f"({100 * c.alt_share:.2f}% of patients) -> A pct {c.A}; o_raw={c.o_raw} "
                      f"(source {c.o_source}) -> O pct {c.O}; layered_S={r.layered_S:.2f}, "
                      f"layered_pct={r.layered_pct:.2f}; minority_alteration={r.minority_alteration}"),
            verdict="holds" if (c.D > 90 or c.A > 90 or c.O > 90) else "weak",
            provenance=("results/lincs/target_role/gene_components.tsv.gz (D120 components); "
                        "results/lincs/fda_ranking_v2_consensus.tsv (D121/D122 layered); "
                        "alteration = TCGA pan-cancer; dependency = DepMap 24Q4; "
                        f"O source for this gene = {c.o_source}"),
            independence=(("O is a VIPER NES over the SAME CollecTRI/ARACNe regulon that L3 and L4 use, so L2's "
                           "strongest component and L4 are NOT independent for this pair; D is DepMap, shared "
                           "with any essentiality reasoning at L6.")
                          if c.o_source == "NES" else
                          ("O is logFC vs GTEx, the same contrast family that L7 uses, so L2's O and L7 are "
                           "entangled; D is DepMap, shared with any essentiality reasoning at L6. A (alteration) "
                           "is independent of both the regulon and the response.")),
            concentration=(f"largest single component's share of the three percentile components: "
                           f"{max(c.D, c.A, c.O) / (c.D + c.A + c.O):.3f} "
                           f"(D {c.D / (c.D + c.A + c.O):.3f}, A {c.A / (c.D + c.A + c.O):.3f}, "
                           f"O {c.O / (c.D + c.A + c.O):.3f})"),
            alternatives=[
                f"only {100 * c.alt_share:.2f}% of patients carry an alteration - the cohort-level claim is "
                f"{'expression/dependency-driven, not genotype-driven' if c.alt_share < 0.2 else 'carried by a large but sub-half genotype subset'}.",
                "the dependency component is a cell-line average (DepMap), not a patient measurement.",
                ("reference artefact (Section 6, D128-A3.2): O vs GTEx normal skin is NEGATIVE here, "
                 "so the normal reference already expresses the target more than the tumour does."
                 if c.o_source != "NES" and c.o_raw < 0 else
                 "reference artefact not excluded: O rests on the same normal reference L7 will use."),
            ]))

        # ---------------- L3 target -> regulators ---------------------------------------------------------
        pp = P[(P.target_gene == p["target"]) & (P.network == "omnipath")]
        is_reg = p["target"] in set(E.regulator)
        hop_cap = 0 if is_reg else 2
        used = pp[(pp.hops <= hop_cap) & pp.tf_has_regulon].copy()
        res = used[used.path_status == "resolved"]
        amb = used[used.path_status == "ambiguous"]
        used.sort_values(["hops", "tf_regulon_genes"], ascending=[True, False]).to_csv(
            OUT / f"l3_paths_{key}.tsv", sep="\t", index=False)
        tot_reg_genes = int(res.tf_regulon_genes.sum())
        top = res.sort_values("tf_regulon_genes", ascending=False).head(3)
        links.append(dict(
            pair=key, link="L3", question="What signed path reaches which regulators, in how many hops?",
            expected=p["expected_route"],
            observed=(f"target is a CollecTRI regulator: {is_reg}; subset rule gives hop cap {hop_cap}. "
                      f"{len(res)} resolved-sign regulators with a regulon, {len(amb)} sign-ambiguous "
                      f"(excluded from the signed footprint). Resolved regulators and their regulon sizes: "
                      + "; ".join(f"{t.tf}(h{t.hops},{t.path_sign:+d},{t.tf_regulon_genes}g)"
                                  for t in res.sort_values(['hops', 'tf_regulon_genes'],
                                                           ascending=[True, False]).itertuples())),
            verdict="holds" if len(res) >= 1 else "not evaluable",
            provenance="results/lincs/d128/a1/target_tf_paths.parquet (D128-A1, OmniPath D126 cache + CollecTRI)",
            independence=("L3 and L4 both derive from the regulon (Section 3): a regulon error moves both. "
                          "For a zero-hop TF anchor L3 carries no independent signalling evidence at all - it is "
                          "the statement 'the target is a TF', and the whole mechanism is L4."
                          if hop_cap == 0 else
                          "L3 adds the OmniPath signed prior, which is a source L4 does not use; but every "
                          "regulator it names is scored through its CollecTRI regulon, so L3 and L4 still share "
                          "the regulon (Section 3)."),
            concentration=(f"largest regulator's share of the regulators' total regulon genes: "
                           f"{(res.tf_regulon_genes.max() / tot_reg_genes) if tot_reg_genes else float('nan'):.3f} "
                           f"({top.iloc[0].tf}: {int(top.iloc[0].tf_regulon_genes)} of {tot_reg_genes} gene slots "
                           f"over {len(res)} regulators)"),
            alternatives=[
                (f"{len(amb)} regulators are reached by sign-ambiguous routes and are silently dropped; "
                 f"if their sign were resolved the predicted directions could change."),
                ("hop cap is a frozen choice (MAX_HOPS=2, P6). At 3 hops this target reaches "
                 f"{int(((pp.hops <= 3) & pp.tf_has_regulon).sum())} regulators with a regulon instead of "
                 f"{int(((pp.hops <= hop_cap) & pp.tf_has_regulon).sum())}; "
                 "the footprint is a function of that cap, not only of the biology."),
            ]))

        # ---------------- L4 regulators -> genes ----------------------------------------------------------
        fp = F[F.unit_id == p["unit_id"]].copy()
        assert len(fp), p["unit_id"]
        fp["symbol"] = fp.gene_id.map(sym)
        fp["expected_sign"] = fp.expected_sign_pre_action * action
        fp["n_regulators_claiming"] = fp.gene_id.map(claims).fillna(0).astype(int)
        own = set(res.tf)
        fp["n_unit_regulators_claiming"] = fp.gene_id.map(
            lambda g: len(reg_of_gene.get(g, set()) & own))
        fp["uniquely_attributable"] = fp.n_regulators_claiming == 1
        signed = fp[fp.expected_sign != 0]
        uniq_frac = float(signed.uniquely_attributable.mean())
        # structural hub share: largest single regulator's share of the signed footprint
        hub = {t: len(reg_of_gene_members(reg_of_gene, signed.gene_id, t)) for t in own}
        hub_tf = max(hub, key=hub.get) if hub else ""
        hub_share = (hub[hub_tf] / len(signed)) if hub else float("nan")
        fp = fp.sort_values(["expected_sign", "symbol"])
        cols = ["unit_id", "cmap_name", "unit_kind", "subset", "gene_id", "symbol",
                "expected_sign_pre_action", "expected_sign", "status", "feature_space",
                "n_regulators_claiming", "n_unit_regulators_claiming", "uniquely_attributable"]
        fp[cols].to_csv(OUT / f"predicted_footprint_{key}.tsv", sep="\t", index=False)
        manifest_counts[key] = dict(footprint_rows=int(len(fp)), signed=int(len(signed)),
                                    up=int((signed.expected_sign > 0).sum()),
                                    down=int((signed.expected_sign < 0).sum()),
                                    conflicted=int((fp.expected_sign == 0).sum()),
                                    unique_fraction=round(uniq_frac, 4),
                                    structural_hub_tf=hub_tf, structural_hub_share=round(hub_share, 4),
                                    median_regulators_claiming=float(signed.n_regulators_claiming.median()))
        links.append(dict(
            pair=key, link="L4", question="Which genes, with what weight and expected direction?",
            expected=(f"the regulators of L3 imply a signed gene footprint; with the drug's action "
                      f"({u.orient}, {action:+d}) each gene's expected direction is "
                      f"action x path sign x tfmode. Weight is unweighted (P7, w=1)."),
            observed=(f"{len(fp)} footprint genes in the LINCS gene space; {len(signed)} carry a signed "
                      f"expectation ({int((signed.expected_sign > 0).sum())} expected UP, "
                      f"{int((signed.expected_sign < 0).sum())} expected DOWN), "
                      f"{int((fp.expected_sign == 0).sum())} conflicted/ambiguous and excluded. "
                      f"Attribution (Section 5): median {signed.n_regulators_claiming.median():.0f} "
                      f"CollecTRI regulators claim each signed gene; only {100 * uniq_frac:.1f}% of signed genes "
                      f"are claimed by exactly one regulator in the whole network. Within this unit's own "
                      f"regulator set, {100 * float((signed.n_unit_regulators_claiming > 1).mean()):.1f}% of signed "
                      f"genes are claimed by more than one of its regulators."),
            verdict="weak" if uniq_frac < 0.5 else "holds",
            provenance=("results/lincs/d128/a1/target_gene_footprints.parquet and regulon_edges.parquet "
                        "(D128-A1, CollecTRI P3); action sign from D121 units.tsv"),
            independence=("same regulon as L3 - one source, reported twice (Section 3). CollecTRI is "
                          "literature-curated, so it is NOT learned from the tumour co-expression L7 uses; "
                          "that entanglement (D133.1) would apply to an ARACNe regulon and does not apply here. "
                          "It does apply to L2's O component for this pair, which is a VIPER NES."
                          if c.o_source == "NES" else
                          "same regulon as L3 - one source, reported twice (Section 3). CollecTRI is "
                          "literature-curated, not learned from tumour co-expression, so the D133.1 "
                          "ARACNe entanglement with L7 does not apply here."),
            concentration=(f"structural: largest single regulator ({hub_tf}) claims {hub_share:.3f} of the "
                           f"{len(signed)} signed genes, over {len(own)} contributing regulators. "
                           f"A footprint whose genes nearly all come from one regulon is one piece of evidence."),
            alternatives=[
                f"attribution: {100 * (1 - uniq_frac):.1f}% of the signed genes are claimed by other regulators "
                f"too; their movement is attributable to this target only by choice, not by observation "
                "(Section 5).",
                "conflicted genes are dropped rather than counted against the prediction; the footprint is "
                "therefore the subset of the regulon that happens to agree with itself.",
                "CollecTRI tfmode signs are curation, and a wrong tfmode flips a gene's expectation without "
                "any data saying so.",
            ]))

    (OUT / "links_L1_L4.json").write_text(json.dumps(links, indent=1, default=str))
    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    (OUT / "t1_run.json").write_text(json.dumps(dict(
        stage="TRACE-T1 (links L1-L4, frozen prediction)", git_sha_at_run=sha,
        procedure="docs/PROCEDURE_mechanism_trace.md",
        reads_no_response_data=True, reads_no_outcome_label=True,
        inputs=["results/lincs/target_role/d121/units.tsv",
                "results/lincs/target_role/gene_components.tsv.gz",
                "results/lincs/fda_ranking_v2_consensus.tsv",
                "results/celligner/disease_line_weights_long.tsv",
                "results/lincs/d128/a1/target_tf_paths.parquet",
                "results/lincs/d128/a1/target_gene_footprints.parquet",
                "results/lincs/d128/a1/regulon_edges.parquet",
                "data/processed/lincs/drug_line_profiles_genes.tsv"],
        pairs=PAIRS, counts=manifest_counts), indent=1, default=str))
    for k, v in manifest_counts.items():
        print(k, v)
    print("wrote", OUT)


def reg_of_gene_members(reg_of_gene, gene_ids, tf):
    return [g for g in gene_ids if tf in reg_of_gene.get(g, set())]


if __name__ == "__main__":
    sys.exit(main())
