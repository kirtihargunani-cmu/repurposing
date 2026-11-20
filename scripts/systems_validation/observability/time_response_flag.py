#!/usr/bin/env python
"""Time-response difference flag per (drug, line) — the user's simplification, 2026-11-19.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

Rather than building a second scoring pipeline per window, flag the drugs whose response DIFFERS
between windows and leave the existing score alone. Everything needed already exists:
  - condition-level profiles keyed on (drug, line, dose, time)   79c2, D153
  - condition-level regulator activity across all 38 networks     79d RUN_TAG=_cond
  - the window assignment (early <= 6 h, late >= 24 h)

For each (drug, line) holding BOTH windows, take the mean regulator-activity vector in each window
and correlate them. A drug whose mechanism reads the same early and late gives a high correlation;
one that reverses — as all five anti-estrogens do (D152) — gives a low or negative one.

This is SPEC_condition_level_profiles section 11's `window_agreement` field, obtained without
touching the scoring code.

Output: results/systems_validation/observability/time_response_flag.tsv.gz + run json
"""
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
PROC = ROOT / "data/processed/lincs"
VDIR = PROC / "drug_viper_specific_v2_tf_cond"
OUT = ROOT / "results/systems_validation/observability"
R_DISCORDANT = 0.0          # r < 0 = the two windows disagree in direction
R_CONCORDANT = 0.5
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    idx = pd.read_csv(PROC / "drug_condition_profiles_index.tsv", sep="\t")
    annot = pd.read_csv(ROOT / "results/lincs/drug_annotations.tsv", sep="\t")
    name = dict(annot.drop_duplicates("pert_id")[["pert_id", "cmap_name"]].values)
    idx["window"] = np.where(idx.pert_time <= 6, "early",
                             np.where(idx.pert_time >= 24, "late", "other"))
    w = idx[idx.window != "other"]
    nb = w.groupby(["pert_id", "ModelID"]).window.nunique()
    both = nb[nb == 2].index
    print(f"{len(idx):,} conditions | {len(nb):,} (drug, line) cells with a usable window | "
          f"{len(both):,} hold BOTH")

    rows = []
    for net in sorted(p.stem for p in VDIR.glob("*.f32")):
        regs = (VDIR / f"{net}.regulators.txt").read_text().split()
        A = np.fromfile(VDIR / f"{net}.f32", dtype=np.float32).reshape(len(idx), len(regs))
        sub = w.set_index(["pert_id", "ModelID"]).loc[both].reset_index()
        for (pid, model), g in sub.groupby(["pert_id", "ModelID"]):
            ge = g[g.window == "early"]; gl = g[g.window == "late"]
            if not len(ge) or not len(gl):
                continue
            ae = A[ge.index.values if False else g.index.values[g.window.values == "early"]].mean(0)
            al = A[g.index.values[g.window.values == "late"]].mean(0)
            if np.std(ae) == 0 or np.std(al) == 0:
                continue
            r = float(np.corrcoef(ae, al)[0, 1])
            rows.append(dict(network=net, pert_id=pid, cmap_name=name.get(pid), ModelID=model,
                             n_early=int((g.window == "early").sum()),
                             n_late=int((g.window == "late").sum()),
                             r_early_late=round(r, 4),
                             flag=("discordant" if r < R_DISCORDANT else
                                   "partial" if r < R_CONCORDANT else "concordant")))
        print(f"  {net}: {len(rows):,} rows so far", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "time_response_flag.tsv.gz", sep="\t", index=False)
    vc = df.flag.value_counts(normalize=True)
    res = dict(step="time-response difference flag per (drug, line, network)",
               n_rows=int(len(df)), n_cells_both_windows=int(len(both)),
               n_networks=int(df.network.nunique()),
               r_median=round(float(df.r_early_late.median()), 4),
               r_p10=round(float(df.r_early_late.quantile(.1)), 4),
               r_p90=round(float(df.r_early_late.quantile(.9)), 4),
               frac_discordant=round(float((df.r_early_late < 0).mean()), 4),
               frac_concordant=round(float((df.r_early_late >= R_CONCORDANT).mean()), 4),
               thresholds=dict(discordant_below=R_DISCORDANT, concordant_at_or_above=R_CONCORDANT))
    (OUT / "time_response_flag_run.json").write_text(json.dumps(res, indent=1))
    print(f"\n{len(df):,} (drug, line, network) rows over {len(both):,} both-window cells\n")
    print(f"  correlation between early and late regulator activity:")
    print(f"     median {res['r_median']:+.3f}   p10 {res['r_p10']:+.3f}   p90 {res['r_p90']:+.3f}")
    print(f"  DISCORDANT (r < 0, the response reverses): {100*res['frac_discordant']:.1f}%")
    print(f"  concordant (r >= 0.5):                     {100*res['frac_concordant']:.1f}%")
    ae = df[df.cmap_name.isin(["fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene"])]
    if len(ae):
        print(f"\n  anti-estrogens ({len(ae)} rows): median r {ae.r_early_late.median():+.3f}, "
              f"{100*(ae.r_early_late < 0).mean():.0f}% discordant")


if __name__ == "__main__":
    sys.exit(main())
