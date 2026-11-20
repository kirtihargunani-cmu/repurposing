#!/usr/bin/env python
"""D157 follow-up: did the D156 TAS floor CHANGE the beta-window answer, or just BLUR it?

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

The per-window beta check (D154) ran at the uniform TAS >= 0.4 floor and found beta changing sign
between 6 h and 24 h for 25.2% of regulators, Spearman(beta_6h, beta_24h) = +0.687. Re-run on the
D156 dataset (6 h floor 0.25) the same statistics are 28.7% and +0.594 - more discordance, on 2,300
both-window cells instead of 1,373.

Two explanations, opposite in meaning:
  (a) the newly admitted 6 h signatures are weaker, so beta_6h is estimated more noisily and the
      extra discordance is measurement error;
  (b) the larger sample reveals real heterogeneity the old subset was too small and too
      selected to show.

This separates them WITHOUT rebuilding anything: stratify the both-window cells by whether their
early condition would have passed the OLD floor, and fit beta separately in each stratum.
  stratum HIGH: max median_tas over the cell's early conditions >= 0.4  (existed before D156)
  stratum NEW : max median_tas over the cell's early conditions <  0.4  (admitted by D156)
If HIGH reproduces ~25% / +0.69 and NEW is worse, explanation (a) holds and the headline number
should be reported from the HIGH stratum with the NEW one carried as a coverage gain.

Output: results/systems_validation/observability/beta_window_tas_strata.{tsv,json}
"""
import json, os
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
PROC = ROOT / "data/processed/lincs"
OUT = ROOT / "results/systems_validation/observability"
VDIR = PROC / "drug_viper_specific_v2_tf_cond"
MIN_N = 40
OLD_FLOOR = 0.4
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "6"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    idx = pd.read_csv(PROC / "drug_condition_profiles_index.tsv", sep="\t")
    genes = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
    gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:]
           for l in open(ROOT / "data/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt")}
    pm = genes.symbol.isin(set(gmt["HALLMARK_E2F_TARGETS"]) | set(gmt["HALLMARK_G2M_CHECKPOINT"])).values
    raw = np.fromfile(PROC / "drug_condition_profiles_raw.f32", dtype=np.float32).reshape(len(idx), -1)
    idx["pi"] = raw[:, pm].mean(1).astype(np.float64)
    del raw
    idx["window"] = np.where(idx.pert_time <= 6, "early", np.where(idx.pert_time >= 24, "late", "other"))
    idx["pos"] = np.arange(len(idx))

    w = idx[idx.window.isin(["early", "late"])]
    nw = w.groupby(["pert_id", "ModelID"]).window.nunique()
    cells = set(nw[nw == 2].index)
    sel = w[pd.MultiIndex.from_arrays([w.pert_id, w.ModelID]).isin(cells)].copy()

    # stratify on the cell's BEST early TAS
    et = (sel[sel.window == "early"].groupby(["pert_id", "ModelID"]).median_tas.max())
    sel["early_tas"] = pd.MultiIndex.from_arrays([sel.pert_id, sel.ModelID]).map(et)
    sel["stratum"] = np.where(sel.early_tas >= OLD_FLOOR, "HIGH", "NEW")
    print(f"both-window cells: {len(cells):,}")
    for s in ("HIGH", "NEW"):
        d = sel[sel.stratum == s]
        print(f"  {s:<5} {d.groupby(['pert_id','ModelID']).ngroups:>5,} cells | "
              f"{int((d.window=='early').sum()):>6,} early / {int((d.window=='late').sum()):>6,} late conditions | "
              f"median early TAS {d[d.window=='early'].median_tas.median():.3f}")

    rows = []
    nets = sorted(p.stem for p in VDIR.glob("*.f32"))
    for ni, net in enumerate(nets, 1):
        regs = (VDIR / f"{net}.regulators.txt").read_text().split()
        A = np.fromfile(VDIR / f"{net}.f32", dtype=np.float32).reshape(len(idx), len(regs))
        for s in ("HIGH", "NEW"):
            d = sel[sel.stratum == s]
            pe = d[d.window == "early"].pos.values
            pl = d[d.window == "late"].pos.values
            if len(pe) < MIN_N or len(pl) < MIN_N:
                continue
            pie, pil = idx.pi.values[pe], idx.pi.values[pl]
            ce, cl = pie - pie.mean(), pil - pil.mean()
            Ae = A[pe].astype(np.float64); Al = A[pl].astype(np.float64)
            be = (Ae - Ae.mean(0)).T @ ce / (ce @ ce)
            bl = (Al - Al.mean(0)).T @ cl / (cl @ cl)
            rows.append(pd.DataFrame({"net": net, "stratum": s, "regulator": regs,
                                      "beta_6h": be, "beta_24h": bl}))
        del A
        if ni % 10 == 0 or ni == len(nets):
            print(f"  ...{ni}/{len(nets)} networks", flush=True)

    out = pd.concat(rows, ignore_index=True)
    out.to_csv(OUT / "beta_window_tas_strata.tsv", sep="\t", index=False)
    from scipy.stats import spearmanr
    summ = {}
    print(f"\n{'stratum':<8}{'fits':>9}{'med b6h':>10}{'med b24h':>10}{'spearman':>10}"
          f"{'med|diff|':>11}{'% of |b|':>10}{'sign flip':>11}")
    for s in ("HIGH", "NEW"):
        d = out[(out.stratum == s)].dropna(subset=["beta_6h", "beta_24h"])
        if not len(d):
            continue
        rho = spearmanr(d.beta_6h, d.beta_24h).correlation
        md = np.median(np.abs(d.beta_24h - d.beta_6h))
        typ = np.median(np.abs(np.concatenate([d.beta_6h, d.beta_24h])))
        flip = float((np.sign(d.beta_6h) != np.sign(d.beta_24h)).mean())
        summ[s] = dict(n_fits=int(len(d)), median_beta_6h=float(d.beta_6h.median()),
                       median_beta_24h=float(d.beta_24h.median()), spearman=float(rho),
                       median_abs_diff=float(md), pct_of_typical_beta=float(100 * md / typ),
                       sign_change_rate=flip)
        print(f"{s:<8}{len(d):>9,}{d.beta_6h.median():>+10.4f}{d.beta_24h.median():>+10.4f}"
              f"{rho:>+10.4f}{md:>11.4f}{100*md/typ:>9.0f}%{100*flip:>10.1f}%")
    json.dump(summ, open(OUT / "beta_window_tas_strata.json", "w"), indent=2)
    print(f"\nreference: D154 at the uniform 0.4 floor was spearman +0.687, sign flip 25.2% "
          f"(1,373 cells); the D156 pooled run is +0.594 / 28.7% (2,300 cells)")


if __name__ == "__main__":
    main()
