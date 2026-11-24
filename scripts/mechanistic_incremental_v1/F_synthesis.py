#!/usr/bin/env python
"""Experiment F (config block F, amendment 1f): per-pair mechanism evidence assembled from Stages A-E. No reranking.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
Protocol 5cc2227 + 6c29417. Reads only the study's own stage outputs and the 2026-11-23 production tables they were
built on. One row per (cancer, drug, window, target unit) for production-eligible FDA rows with a timepoint; a status
per evidence layer by the frozen rules; rule-chosen worked cases.
Outputs: results/lincs/mechanistic_incremental_v1/F_synthesis/{F_mechanism_evidence.tsv.gz, F_evidence_summary.tsv,
  F_mechanism_cases.md, run.json}  (F_decision_matrix.md and F_final_report.md are written from the stage reports)
"""
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("A", HERE / "A_resolvable.py")
A = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(A)
CFG, IN, RES = A.CFG, A.IN, A.RES
ST = A.WT / CFG["output_root"]
OUT = ST / "F_synthesis"
OUT.mkdir(parents=True, exist_ok=True)
LAYERS = ["disease", "perturbation", "temporal", "target", "network", "program"]


def main():
    t0 = time.time()
    G = pd.read_csv(RES / "lincs/gene_ranking_v2_tf.tsv", sep="\t", low_memory=False)
    G["ing"] = G.fda_ingredient.fillna(G.cmap_name)
    E = G[G.eligible & G.window.isin(["early", "late"]) & G.gene_pct.notna()].copy()
    # perturbation + split availability (A)
    H = pd.read_csv(ST / "A_resolvable/A_halves.tsv.gz", sep="\t", usecols=["unit", "pert_id", "window"])
    E["split_available"] = [k in set(zip(H.unit, H.pert_id, H.window)) for k in zip(E.unit, E.pert_id, E.window)]
    E["st_disease"] = "supported"
    E["st_perturbation"] = np.select([E.gene_pct >= 75, E.gene_pct <= 25], ["supported", "conflicting"], "unknown")
    # temporal (B, production arm)
    B = pd.read_csv(ST / "B_temporal/B_transitions.tsv.gz", sep="\t")
    B = B[B.arm == "production"].drop_duplicates(["unit", "pert_id"]).set_index(["unit", "pert_id"])
    def temporal(r):
        k = (r.unit, r.pert_id)
        if k not in B.index:
            return "not_tested", None
        b = B.loc[k]; cls = b.rev_gene_class
        z = b[f"rev_gene_{r.window}"] / b[f"rev_gene_se_{r.window}"] if b[f"rev_gene_se_{r.window}"] > 0 else 0
        if cls == "sign_reversal":
            return "conflicting", cls
        return ("supported" if cls == "persistent_same_sign" or abs(z) > 1.96 else "unknown"), cls
    tt = [temporal(r) for r in E.itertuples()]
    E["st_temporal"] = [t[0] for t in tt]; E["b_class"] = [t[1] for t in tt]
    # target (C, primary role)
    C = pd.read_csv(ST / "C_target_recovery/C_pairs.tsv.gz", sep="\t")
    Cp = C[C.role == "primary"].groupby(["unit", "pert_id", "window"]).agg(c_exact=("exact", "max"), c_down=("down", "max"),
                                                                          primary_targets=("target", lambda s: ",".join(sorted(set(s)))))
    E = E.merge(Cp.reset_index(), on=["unit", "pert_id", "window"], how="left")
    best = E[["c_exact", "c_down"]].max(axis=1)
    E["st_target"] = np.where(E.primary_targets.isna(), "not_tested",
                              np.select([best >= 0.75, best <= 0.25], ["supported", "conflicting"], "unknown"))
    # network (D flags of the primary-target regulator in pairs that include the unit's network)
    D = pd.read_csv(ST / "D_network_sensitivity/D_regulator_flags.tsv.gz", sep="\t")
    ovr = pd.read_csv(IN / "config/group_networks.tsv", sep="\t").set_index("group_id").net.to_dict()
    cohort_of = pd.read_csv(RES / "step1/disease_groups_toil.tsv", sep="\t").drop_duplicates("group_id").set_index("group_id").cohort.to_dict()
    parent_of = pd.read_csv(RES / "celligner/stratum_line_weights_long_v2.tsv", sep="\t").drop_duplicates("stratum").set_index("stratum").parent.to_dict()
    avail = {p.stem for p in A.VDIR.glob("*.f32")}
    netof = {u: A.net_for(parent_of.get(u, u), ovr, avail, cohort_of) for u in E.unit.unique()}
    Dl = pd.concat([D.rename(columns={"net_a": "net"}).drop(columns="net_b"), D.rename(columns={"net_b": "net"}).drop(columns="net_a")])
    Dg = Dl.groupby(["net", "window", "regulator"]).flag.apply(set).to_dict()
    def network(r):
        if not isinstance(r.primary_targets, str):
            return "not_tested"
        fl = set().union(*[Dg.get((netof[r.unit], r.window, t), set()) for t in r.primary_targets.split(",")])
        if not fl:
            return "not_tested"
        return "conflicting" if "rewired" in fl else ("supported" if fl == {"context_robust"} else "unknown")
    E["st_network"] = [network(r) for r in E.itertuples()]
    # program (E): linked to a reproducibly perturbable program of that cancer and window
    P = pd.read_csv(ST / "E_program_access/programs/programs.tsv", sep="\t")
    P["members_list"] = P.members.str.split(",").apply(set)
    REL = pd.read_csv(ST / "E_program_access/E_repeatability.tsv", sep="\t")
    pert = set(zip(REL[REL.perturbable == True].unit, REL[REL.perturbable == True].window, REL[REL.perturbable == True].program))
    om = pd.read_csv(IN / "data/raw/omnipath/omnipath_signed_directed.tsv.gz", sep="\t", usecols=["source_genesymbol", "target_genesymbol"])
    hop = om.groupby("source_genesymbol").target_genesymbol.apply(set).to_dict()
    progs_by_unit = {u: g for u, g in P.groupby("cohort")}
    def program(r):
        if not isinstance(r.primary_targets, str) or r.unit not in progs_by_unit:
            return "not_tested", None
        t = set(r.primary_targets.split(","))
        link = t | set().union(*[hop.get(g, set()) for g in t])
        hits = [int(p.program) for p in progs_by_unit[r.unit].itertuples() if p.members_list & link]
        if not hits:
            return "not_tested", None
        good = [h for h in hits if (r.unit, r.window, h) in pert]
        return ("supported" if good else "unknown"), ",".join(map(str, good or hits))
    pp = [program(r) for r in E.itertuples()]
    E["st_program"] = [p[0] for p in pp]; E["linked_programs"] = [p[1] for p in pp]
    lab = pd.read_csv(RES / "lincs/labels/pair_labels_v2_tf.tsv.gz", sep="\t", low_memory=False, usecols=["unit_id", "pert_id", "label_phase"])
    E = E.merge(lab.drop_duplicates(["unit_id", "pert_id"]).rename(columns={"unit_id": "unit"}), on=["unit", "pert_id"], how="left")
    E["novel"] = E.label_phase.fillna(0) <= 0
    keep = ["level", "unit", "window", "ing", "pert_id", "cmap_name", "layered_genes", "layered_T", "layered_D", "layered_A",
            "layered_O", "gene_pct", "gene_rank", "rev_gene", "n_conditions", "split_available", "b_class", "primary_targets",
            "c_exact", "c_down", "linked_programs", "label_phase", "novel"] + [f"st_{l}" for l in LAYERS]
    E[keep].to_csv(OUT / "F_mechanism_evidence.tsv.gz", sep="\t", index=False, float_format="%.4g")
    SUM = pd.concat([E.groupby(["window", f"st_{l}"]).size().rename("rows").reset_index().rename(columns={f"st_{l}": "status"}).assign(layer=l)
                     for l in LAYERS])[["layer", "window", "status", "rows"]]
    SUM["share"] = SUM.rows / SUM.groupby(["layer", "window"]).rows.transform("sum")
    SUM.to_csv(OUT / "F_evidence_summary.tsv", sep="\t", index=False, float_format="%.4f")

    # ---- rule-chosen cases ----
    inst = pd.read_csv(RES / "lincs/d127/instances.tsv", sep="\t")
    pos = set(zip(inst[inst.pair_set == "primary"].high, inst[inst.pair_set == "primary"].drug.str.lower()))
    Gm = G[G.window.isin(["early", "late"])].copy()
    Ch = pd.read_csv(RES / "lincs/gene_ranking_chemo_v2_tf.tsv", sep="\t", low_memory=False)
    Ch = Ch[Ch.window.isin(["early", "late"]) & Ch.gene_pct.notna()].assign(ing=lambda d: d.fda_ingredient.fillna(d.cmap_name))
    sec = set(C[(C.role == "secondary") & C.pchembl.notna()].pert_id)
    cats = {
        "i positive control (D127 primary pair, high cohort)": E[[(u, str(c).lower()) in pos for u, c in zip(E.unit, E.cmap_name)]],
        "ii high target evidence, low gene-level reversal": E[(E.layered_T >= 75) & (E.gene_pct <= 25)],
        "iii lower target evidence, high gene-level reversal": E[(E.layered_T < 60) & (E.gene_pct >= 90)],
        "iv confident time transition": E[E.b_class.isin(["sign_reversal", "persistent_same_sign", "early_only", "late_only"])],
        "v secondary-target drug (measured binding)": E[E.pert_id.isin(sec)],
        "vi filtered out (T < 50)": Gm[(Gm.path == "D/A/O") & ~Gm.eligible & Gm.rev_gene.notna()].assign(gene_pct=lambda d: d.groupby(["unit", "window"]).rev_gene.rank(pct=True) * 100),
        "vii chemotherapy table top": Ch,
        "viii top novel, >= 2 conditions": E[E.novel & (E.n_conditions >= 2)],
    }
    lines = ["# Worked mechanism cases (chosen by the frozen category rules, not by outcome)\n",
             "One case per category and window: the highest gene-level percentile in the category, ties by ingredient name. "
             "Statuses come from the stage outputs; nothing here re-ranks drugs.\n"]
    for cat, df in cats.items():
        lines.append(f"\n## {cat}\n")
        for w in ("early", "late"):
            d = df[df.window == w].sort_values(["gene_pct", "ing"], ascending=[False, True])
            if not len(d):
                lines.append(f"- **{w}:** no row in this category.\n"); continue
            r = d.iloc[0]
            st = " · ".join(f"{l} {r.get(f'st_{l}', 'n/a')}" for l in LAYERS if f"st_{l}" in d.columns)
            lines.append(f"- **{w}: {r.ing}** in {r.unit} (target unit {r.get('layered_genes', 'n/a')}; T {r.get('layered_T', float('nan')):.1f}; "
                         f"gene-level percentile {r.gene_pct:.0f}; {int(r.n_conditions) if pd.notna(r.get('n_conditions')) else 'n/a'} conditions)"
                         + (f". Evidence: {st}" if st else "") + ".\n")
    (OUT / "F_mechanism_cases.md").write_text("".join(lines))
    run = dict(stage="F_synthesis", protocol="5cc2227+6c29417", evidence_rows=len(E), runtime_s=round(time.time() - t0, 1),
               finished=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
    (OUT / "run.json").write_text(json.dumps(run, indent=1))
    pd.set_option("display.width", 200)
    print(SUM.to_string(index=False)); print("".join(lines))


if __name__ == "__main__":
    sys.exit(main())
