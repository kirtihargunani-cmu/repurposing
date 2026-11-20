#!/usr/bin/env python
"""D153 section 10.3: is the proliferation adjustment's beta the same at 6 h and 24 h?

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

79i removes growth arrest from every regulator's activity by OLS, A_adj = A - beta * pi, with ONE
beta per (network, regulator) fitted over all profiles pooled. pi is the profile's proliferation
score: the mean over HALLMARK_E2F_TARGETS + HALLMARK_G2M_CHECKPOINT on the RAW profiles (79i:191).

If growth arrest relates to regulator activity differently early and late - which pharmacology says
it must, since at 6 h a cytostatic drug has barely arrested anything - then a pooled beta, fitted on
data that is overwhelmingly 24 h, is applied to 6 h conditions where the slope is different. For a
drug whose mechanism is early and whose arrest is late, that subtracts a late-fitted component from
an early profile that does not yet contain one.

The test is WITHIN-CELL by construction: only (drug, line) cells holding BOTH a 6 h and a 24 h
condition contribute, so the same drugs and lines appear on both sides and the comparison cannot be
confounded by which drugs were profiled when.

Per (network, regulator): fit beta_6 over that network's 6 h conditions and beta_24 over its 24 h
conditions, both on the same cells, and compare.

Output: results/systems_validation/observability/beta_window_check.{tsv,json}
"""
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
PROC = ROOT / "data/processed/lincs"
OUT = ROOT / "results/systems_validation/observability"
VDIR = PROC / "drug_viper_specific_v2_tf_cond"
MIN_N = 40
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "6"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    idx = pd.read_csv(PROC / "drug_condition_profiles_index.tsv", sep="\t")
    genes = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
    gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:]
           for l in open(ROOT / "data/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt")}
    pm = genes.symbol.isin(set(gmt["HALLMARK_E2F_TARGETS"]) | set(gmt["HALLMARK_G2M_CHECKPOINT"])).values
    print(f"proliferation genes on the profiles: {int(pm.sum())}")
    raw = np.fromfile(PROC / "drug_condition_profiles_raw.f32", dtype=np.float32).reshape(len(idx), -1)
    idx["pi"] = raw[:, pm].mean(1).astype(np.float64)

    # within-cell set: cells holding BOTH a 6 h and a 24 h condition
    idx["window"] = np.where(idx.pert_time <= 6, "early", np.where(idx.pert_time >= 24, "late", "other"))
    w = idx[idx.window.isin(["early", "late"])]
    both = w.groupby(["pert_id", "ModelID"]).window.nunique()
    cells = both[both == 2].index
    sel = w.set_index(["pert_id", "ModelID"]).loc[cells].reset_index()
    print(f"cells with BOTH windows: {len(cells):,}  -> {len(sel):,} conditions "
          f"({int((sel.window=='early').sum()):,} early, {int((sel.window=='late').sum()):,} late)")
    print(f"pi: early mean {sel[sel.window=='early'].pi.mean():+.4f}, "
          f"late mean {sel[sel.window=='late'].pi.mean():+.4f}")

    rows = []
    for net in sorted(p.stem for p in VDIR.glob("*.f32")):
        regs = (VDIR / f"{net}.regulators.txt").read_text().split()
        A = np.fromfile(VDIR / f"{net}.f32", dtype=np.float32).reshape(len(idx), len(regs))
        for wname in ["early", "late"]:
            pass
        e = sel[sel.window == "early"].index.values
        l = sel[sel.window == "late"].index.values
        # sel index is positional into idx only after reset; rebuild positions
        key = idx.pert_id + "|" + idx.ModelID + "|" + idx.pert_dose.astype(str) + "|" + idx.pert_time.astype(str)
        skey = sel.pert_id + "|" + sel.ModelID + "|" + sel.pert_dose.astype(str) + "|" + sel.pert_time.astype(str)
        posmap = pd.Series(np.arange(len(idx)), index=key.values)
        pe = posmap[skey[sel.window == "early"].values].values
        pl = posmap[skey[sel.window == "late"].values].values
        if len(pe) < MIN_N or len(pl) < MIN_N:
            print(f"{net}: too few conditions ({len(pe)} early, {len(pl)} late)"); continue
        pie, pil = idx.pi.values[pe], idx.pi.values[pl]
        for j, r in enumerate(regs):
            ae, al = A[pe, j].astype(np.float64), A[pl, j].astype(np.float64)
            be = np.polyfit(pie, ae, 1)[0] if np.std(pie) > 1e-9 else np.nan
            bl = np.polyfit(pil, al, 1)[0] if np.std(pil) > 1e-9 else np.nan
            if not (np.isfinite(be) and np.isfinite(bl)):
                continue
            rows.append(dict(network=net, regulator=r, beta_early=be, beta_late=bl,
                             n_early=len(pe), n_late=len(pl),
                             r_early=float(np.corrcoef(pie, ae)[0, 1]),
                             r_late=float(np.corrcoef(pil, al)[0, 1])))
        print(f"{net}: {len(regs)} regulators, {len(pe)} early / {len(pl)} late conditions")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "beta_window_check.tsv", sep="\t", index=False)
    d = df.beta_late - df.beta_early
    res = dict(test="D153 section 10.3 - proliferation beta, 6 h vs 24 h, within-cell",
               n_cells_both_windows=int(len(cells)), n_regulator_rows=int(len(df)),
               beta_early=dict(median=round(float(df.beta_early.median()), 4),
                               p10=round(float(df.beta_early.quantile(.1)), 4),
                               p90=round(float(df.beta_early.quantile(.9)), 4)),
               beta_late=dict(median=round(float(df.beta_late.median()), 4),
                              p10=round(float(df.beta_late.quantile(.1)), 4),
                              p90=round(float(df.beta_late.quantile(.9)), 4)),
               spearman_early_vs_late=round(float(df.beta_early.corr(df.beta_late, method="spearman")), 4),
               median_abs_difference=round(float(d.abs().median()), 4),
               frac_sign_disagreement=round(float((np.sign(df.beta_early) != np.sign(df.beta_late)).mean()), 4),
               median_abs_beta=round(float(pd.concat([df.beta_early, df.beta_late]).abs().median()), 4))
    res["difference_relative_to_typical_beta"] = round(res["median_abs_difference"] / max(res["median_abs_beta"], 1e-9), 3)
    (OUT / "beta_window_check.json").write_text(json.dumps(res, indent=1))
    print(f"\n{len(df):,} (network, regulator) fits over {len(cells):,} both-window cells\n")
    print(f"  beta at  6h: median {res['beta_early']['median']:+.4f}  (p10 {res['beta_early']['p10']:+.4f}, p90 {res['beta_early']['p90']:+.4f})")
    print(f"  beta at 24h: median {res['beta_late']['median']:+.4f}  (p10 {res['beta_late']['p10']:+.4f}, p90 {res['beta_late']['p90']:+.4f})")
    print(f"\n  Spearman(beta_6h, beta_24h) across regulators = {res['spearman_early_vs_late']:+.4f}")
    print(f"  median |beta_24h - beta_6h| = {res['median_abs_difference']:.4f}  "
          f"= {100*res['difference_relative_to_typical_beta']:.0f}% of a typical |beta|")
    print(f"  regulators where beta CHANGES SIGN between windows: {100*res['frac_sign_disagreement']:.1f}%")


if __name__ == "__main__":
    sys.exit(main())
