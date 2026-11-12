#!/usr/bin/env python
"""Precursor (REPORTED ONLY): does reversal computed INSIDE a regulon's own target genes surface anti-estrogens
in the breast cohorts, where the aggregate regulator score does not?

No threshold, no pass/fail, no D-number; nothing in production is read for scoring or written. This tests one
premise before any spec is written. D126-B1 found that anti-estrogens are never examined on their own axis, because
ESR1 is not among a drug's 25 largest regulator shifts, and the drug-side dimension diagnostic put ESR1 at median
rank #744 on PC1. The question here is whether dropping a level - scoring on the regulon's target GENES with their
tfmode / likelihood weights kept, instead of on the collapsed regulator activity - recovers the ESR1 axis.

The quantity, parallel to production's  rev = -sum_i w_i a_i / sqrt(sum_i w_i^2)  (79e, SCORE=weighted) but one
level down, restricted to regulator i's own targets:

    w_ij   = |tfmode_ij| * likelihood_ij * d_j        target importance x the disease's gene-level z
    rev_i  = - (sum_j w_ij g_kj) / sqrt(sum_j w_ij^2)  for drug k

d_j is the Step 4 tumor-vs-GTEx z (the same signature Step 6 reads); g_kj is the drug's line-specific gene profile
(79c), consensus over the cohort's lines exactly as 79e builds drug activity: Celligner-weighted over in-key lines
for direct drugs, unweighted over out-of-key lines for borrowed ones.

Readouts, per cohort:
  1. where each anti-estrogen ranks among all the cohort's drugs under ESR1-restricted reversal;
  2. the same under production's aggregate rev, for comparison (results/lincs/drug_reversal_split_v2_consensus.tsv);
  3. where ESR1 ranks among ALL regulators by how highly it ranks the anti-estrogens - i.e. is ESR1's regulon
     special, or would any regulon do as well;
  4. the same numbers for the proliferation regulators, which is what the score currently keys on;
  5. estradiol as an anti-control: an ER agonist should not look like a reverser.

Inputs are read-only. Output: results/lincs/diagnostics/precursor_esr1_regulon/.
"""
import os

# --- config block: every constant here and nowhere else ---------------------------------------------------------
THREADS = 8
COHORT_NET = {"BRCA-Basal": "BRCA_BAS", "BRCA-Her2": "BRCA",          # 79e's map: NETMAP, else net_of(cohort)
              "BRCA-LumA": "BRCA_LUMA", "BRCA-LumB": "BRCA_LUMB"}
ARMS = ["aracne", "collectri"]          # production scores the consensus of the two; CollecTRI is the curated set
ANTI_ESTROGENS = ["fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene", "clomifene"]
ANTI_CONTROL = ["estradiol"]            # ER agonist: should NOT look like a reverser
PROLIF = ["FOXM1", "E2F1", "E2F2", "E2F7", "E2F8", "MYBL2", "CENPA", "MYC", "TEAD4"]
MIN_TARGETS = 10                        # as the pipeline prunes regulons elsewhere
WEIGHT = "abs_tfmode_x_likelihood"      # target importance; the sign of the move is already in d_j and g_kj
SIGNATURE = "z_tumor_vs_gtex.tsv"       # breast cohorts use the GTEx arm (BLCA is the adjacent one; not used here)

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = str(THREADS)

import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.sparse import csr_matrix  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PROC = ROOT / "data/processed/lincs"
REG_DIR = PROC / "precursor_regulons"
PROD = ROOT / "results/lincs/drug_reversal_split_v2_consensus.tsv"
SIG = ROOT / "results/step4_v2" / SIGNATURE
WEIGHTS = ROOT / "results/celligner/disease_line_weights_long.tsv"
OUT = ROOT / "results/lincs/diagnostics/precursor_esr1_regulon"


def pct(v):
    """Percentile of each entry among the whole vector, as 79e's rev_pct: share of drugs scoring lower."""
    r = np.argsort(np.argsort(v))
    return r / max(len(v) - 1, 1)


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    genes = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    idx["col"] = np.arange(len(idx))
    G = np.memmap(PROC / "drug_line_profiles_bing_specific.f32", dtype=np.float32,
                  mode="r", shape=(len(idx), len(genes)))
    gpos = pd.Index(genes.gene_id)
    print(f"drug x line gene profiles: {G.shape[0]} profiles x {G.shape[1]} genes")

    sig = pd.read_csv(SIG, sep="\t")
    sig["gene_id"] = sig.gene_id.astype(str).str.strip('"')
    lw = pd.read_csv(WEIGHTS, sep="\t")
    prod = pd.read_csv(PROD, sep="\t") if PROD.exists() else None

    rows, reg_rows, summary = [], [], {}
    for coh, arm in [(c, a) for c in COHORT_NET for a in ARMS]:
        net = COHORT_NET[coh] if arm == "aracne" else "collectri"
        d = sig.set_index("gene_id")[coh].reindex(gpos).to_numpy(float)
        have = np.isfinite(d)
        d = np.where(have, d, 0.0)

        drugs = pd.read_csv(PROC / f"d126_act/{coh}.drugs.tsv", sep="\t")      # the cohort's scored drugs (79e)
        w = lw[lw.group_id == coh].set_index("ModelID")
        inkey = set(w.index[w.in_key.astype(bool)])
        wg = w.weight

        # drug x gene consensus, built as 79e builds drug activity
        M = np.zeros((len(drugs), G.shape[1]), dtype=np.float64)
        by = idx.groupby("pert_id")
        keep = np.ones(len(drugs), bool)
        for i, (pid, src) in enumerate(zip(drugs.pert_id, drugs.source)):
            if pid not in by.groups:
                keep[i] = False
                continue
            p = by.get_group(pid)
            p = p[p.ModelID.isin(inkey)] if src == "direct" else p[~p.ModelID.isin(inkey)]
            if not len(p):
                keep[i] = False
                continue
            ww = wg.reindex(p.ModelID).fillna(0).to_numpy() if src == "direct" else np.zeros(len(p))
            ww = ww if ww.sum() > 0 else np.ones(len(p))
            M[i] = (ww[:, None] * np.asarray(G[p.col.values], dtype=np.float64)).sum(0) / ww.sum()
        drugs, M = drugs[keep].reset_index(drop=True), M[keep]

        # regulon-restricted reversal for every regulator at once: REV = -(M @ Wn), Wn columns unit-norm
        reg = pd.read_csv(REG_DIR / f"{net}.tsv.gz", sep="\t")
        reg = reg[reg.target.isin(gpos)]
        reg = reg[reg.groupby("regulator").regulator.transform("size") >= MIN_TARGETS]
        regs = sorted(reg.regulator.unique())
        ri = pd.Index(regs).get_indexer(reg.regulator)
        gi = gpos.get_indexer(reg.target)
        wv = reg.tfmode.abs().to_numpy() * reg.likelihood.to_numpy() * d[gi]
        W = csr_matrix((wv, (gi, ri)), shape=(len(gpos), len(regs))).toarray()
        nrm = np.linalg.norm(W, axis=0)
        ok = nrm > 0
        REV = -(M @ (W[:, ok] / nrm[ok]))                      # drugs x regulators
        regs = [r for r, k in zip(regs, ok) if k]
        P = np.apply_along_axis(pct, 0, REV)                   # percentile within the cohort's drugs, per regulator
        print(f"{coh:12s} {arm:10s} net {net:10s} {len(drugs):5d} drugs  {len(regs):5d} regulons >= {MIN_TARGETS} targets")

        # polarity check: the SAME construction over all genes must agree with production's regulator-level rev.
        # If it does not, the restricted numbers below are a bug, not a finding.
        gv = -(M @ (d / np.linalg.norm(d)))
        rho_prod = np.nan
        if prod is not None:
            pc = prod[prod.group_id == coh].set_index("cmap_name").rev
            pc = pc[~pc.index.duplicated()]
            a = pc.reindex(drugs.cmap_name.astype(str)).to_numpy(float)
            m = np.isfinite(a)
            rho_prod = float(spearmanr(gv[m], a[m]).correlation) if m.sum() > 10 else np.nan
        print(f"             polarity check: unrestricted gene-level rev vs production rev, "
              f"Spearman {rho_prod:.3f} over {int(np.isfinite(a).sum()) if prod is not None else 0} drugs")

        rpos = {r: j for j, r in enumerate(regs)}
        dpos = {n: i for i, n in enumerate(drugs.cmap_name.astype(str))}
        prod_c = prod[prod.group_id == coh].set_index("cmap_name") if prod is not None else None

        # 3. is ESR1's regulon special? score every regulon by the mean percentile it gives the anti-estrogens
        ae_i = [dpos[n] for n in ANTI_ESTROGENS if n in dpos]
        reg_score = P[ae_i].mean(0) if ae_i else np.full(len(regs), np.nan)
        order = np.argsort(-reg_score)
        esr1_rank = int(np.where(np.array(regs)[order] == "ESR1")[0][0]) + 1 if "ESR1" in rpos else None
        reg_rows += [dict(cohort=coh, arm=arm, regulator=r, mean_ae_pct=round(float(reg_score[rpos[r]]), 4),
                          rank=int(np.where(order == rpos[r])[0][0]) + 1, n_targets=int((ri == j).sum()))
                     for j, r in enumerate(regs) if r in ("ESR1", *PROLIF)]

        for name in ANTI_ESTROGENS + ANTI_CONTROL:
            if name not in dpos:
                continue
            i = dpos[name]
            row = dict(cohort=coh, arm=arm, drug=name, kind="anti_control" if name in ANTI_CONTROL else "anti_estrogen",
                       esr1_pct=round(float(P[i, rpos["ESR1"]]), 4) if "ESR1" in rpos else None,
                       esr1_rev=round(float(REV[i, rpos["ESR1"]]), 4) if "ESR1" in rpos else None)
            for p in PROLIF:
                row[f"{p.lower()}_pct"] = round(float(P[i, rpos[p]]), 4) if p in rpos else None
            row["median_regulon_pct"] = round(float(np.median(P[i])), 4)
            row["esr1_rank_among_regulons"] = (int((P[i] > P[i, rpos["ESR1"]]).sum()) + 1) if "ESR1" in rpos else None
            row["n_regulons"] = len(regs)
            if prod_c is not None and name in prod_c.index:
                r0 = prod_c.loc[name]
                r0 = r0.iloc[0] if isinstance(r0, pd.DataFrame) else r0
                row["prod_rev_pct"] = round(float(r0.rev_pct), 4) if pd.notna(r0.get("rev_pct")) else None
                row["prod_spec_pct"] = round(float(r0.spec_pct), 4) if pd.notna(r0.get("spec_pct")) else None
            rows.append(row)

        summary[f"{coh}|{arm}"] = dict(cohort=coh, arm=arm, net=net, n_drugs=int(len(drugs)), n_regulons=len(regs),
                            polarity_rho_vs_production=round(rho_prod, 3),
                            esr1_in_network="ESR1" in rpos,
                            esr1_rank_by_mean_ae_pct=esr1_rank,
                            esr1_mean_ae_pct=(round(float(reg_score[rpos["ESR1"]]), 4) if "ESR1" in rpos else None),
                            best_regulon=str(np.array(regs)[order][0]),
                            best_regulon_mean_ae_pct=round(float(reg_score[order][0]), 4))

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "anti_estrogen_rows.tsv", sep="\t", index=False)
    pd.DataFrame(reg_rows).to_csv(OUT / "regulon_ranking.tsv", sep="\t", index=False)
    (OUT / "run.json").write_text(json.dumps(dict(weight=WEIGHT, min_targets=MIN_TARGETS,
                                                  signature=SIGNATURE, cohorts=summary), indent=1))

    pd.set_option("display.width", 250)
    ae = out[out.kind == "anti_estrogen"]
    cols = ["cohort", "arm", "drug", "esr1_pct", "prod_rev_pct", "prod_spec_pct", "median_regulon_pct",
            "esr1_rank_among_regulons", "foxm1_pct", "e2f1_pct"]
    print("\n=== anti-estrogens: ESR1-restricted reversal percentile vs production ===")
    print(ae[[c for c in cols if c in ae]].to_string(index=False))
    print("\n=== estradiol (ER agonist, should not look like a reverser) ===")
    print(out[out.kind == "anti_control"][[c for c in cols if c in out]].to_string(index=False))
    print("\n=== is ESR1's regulon special? rank among all regulons by mean anti-estrogen percentile ===")
    for c, s in summary.items():
        print(f"  {c:24s} ESR1 rank {s['esr1_rank_by_mean_ae_pct']} of {s['n_regulons']}"
              f"  (ESR1 mean pct {s['esr1_mean_ae_pct']}, best regulon {s['best_regulon']}"
              f" {s['best_regulon_mean_ae_pct']})")
    for a in ARMS:
        x = ae[ae.arm == a]
        print(f"\n[{a}] ESR1 percentile over the anti-estrogen rows: median {x.esr1_pct.median():.3f}, "
              f"max {x.esr1_pct.max():.3f}")
    print(f"\nESR1 median percentile over all anti-estrogen rows: {ae.esr1_pct.median():.3f}"
          f"  |  production rev_pct median: {ae.prod_rev_pct.median():.3f}"
          if "prod_rev_pct" in ae else "")
    print(f"\nwritten to {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
