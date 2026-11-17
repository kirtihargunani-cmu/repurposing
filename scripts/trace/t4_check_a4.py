#!/usr/bin/env python
"""Mechanism trace, stage 4: AFTER-THE-FACT check of L6 against D128-A4 and of L5 against D130.

Run and committed only after scripts/trace/t2_observe.py and its L6 record were committed
(commit "TRACE: links L5-L7..."). Nothing here feeds back into the trace records; it only
reports whether the independent L6 reproduces what A4 already found, and whether the L5 verdict
agrees with the D130 corrected null.
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
OUT = ROOT / "results/lincs/trace"


def main():
    mine = json.loads((OUT / "t2_run.json").read_text())["summary"]
    rows = []

    a4 = pd.read_csv(ROOT / "results/lincs/d128/a4/braf_hub_dominance.tsv", sep="\t")
    v = a4[(a4.cmap_name == "vemurafenib") & (a4.ModelID == "ACH-000219") &
           (a4.gene_arm == "all_bing")].iloc[0]
    m = mine["vemurafenib_BRAF_TCGA_SKCM"]
    for name, theirs, ours in [
            ("hub share of footprint genes", v.share_genes_MYC, m["hub_structural_share"]),
            ("hub share of |contribution|", v.share_abs_numerator_MYC, m["hub_signal_share"]),
            ("M with hub", v.M, m["mean_with_hub"]),
            ("M without hub", v.M_without_MYC, m["mean_without_hub"]),
            ("hub gene count", v.n_genes_MYC, round(m["hub_structural_share"] * m["n_covered"]))]:
        rows.append(dict(pair="vemurafenib_BRAF_TCGA_SKCM", quantity=name, trace_L6=round(float(ours), 4),
                         d128_a4=round(float(theirs), 4), agrees=abs(float(ours) - float(theirs)) < 5e-3,
                         note="hub identified independently as MYC by both"))
    rows.append(dict(pair="vemurafenib_BRAF_TCGA_SKCM", quantity="hub share of SIGNED numerator",
                     trace_L6=float("nan"), d128_a4=round(float(v.share_signed_numerator_MYC), 4), agrees=None,
                     note="A4 also reports the signed-numerator share (0.859), a harsher reading than the "
                          "absolute share the trace used; the trace did not compute it"))

    d130 = pd.read_csv(ROOT / "results/lincs/d130/a4_controls_corrected.tsv", sep="\t")
    for pair, pid, mid in [("fulvestrant_ESR1_BRCA_LumA", "BRD-A90490067", "ACH-000019"),
                           ("vemurafenib_BRAF_TCGA_SKCM", "BRD-K56343971", "ACH-000219")]:
        src = d130
        s = src[(src.pert_id == pid) & (src.ModelID == mid) & (src.gene_arm == "all_bing") &
                (src.profile_arm == "specific")]
        if not len(s):
            e = pd.read_csv(ROOT / "results/lincs/d130/esr1_control_corrected.tsv", sep="\t")
            s = e[(e.pert_id == pid) & (e.ModelID == mid) & (e.gene_arm == "all_bing") &
                  (e.profile_arm == "specific")]
        s = s.iloc[0]
        rows.append(dict(pair=pair, quantity="corrected null (crp_sign_1000_z / cgn_1000_z)",
                         trace_L6=mine[pair]["generic_percentile"],
                         d128_a4=f"{s.crp_sign_1000_z:.2f} / {s.cgn_1000_z:.2f}", agrees=None,
                         note=("the trace's L5 reference is the same footprint under every compound in the line "
                               "(percentile in trace_L6); D130's is a correlation-preserving permutation. Both "
                               "put this pair below any conventional threshold, from different references.")))
    pd.DataFrame(rows).to_csv(OUT / "l6_agreement_with_d128_a4.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).to_string())


if __name__ == "__main__":
    main()
