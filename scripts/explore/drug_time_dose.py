#!/usr/bin/env python
"""Exploratory (2026-11-03): early vs late responders and a dose approximation for LINCS compounds (runs in `celligner`).

Input: level-5 hiq signatures (landmark + best-inferred genes, MODZ z; scripts/76b), all signatures, not only active
ones, since inactive low doses and early time points are the information here. Similarity on the 978 landmark genes.
Activity = TAS >= 0.4 (the drug-side activity rule).
1 Time (drug x line with 6 h and 24 h signatures): per time point, the signature at the dose shared by both time points
  (highest shared dose; else each time point's highest dose). Classes:
    early_sustained  active at 6 h and 24 h, profiles correlate r >= 0.5   (the 24 h response is the 6 h one, persisting)
    early_evolving   active at 6 h and 24 h, r < 0.5                      (a second program develops)
    late             inactive at 6 h, active at 24 h                       (secondary / slow response)
    transient        active at 6 h, inactive at 24 h
    inactive         neither
2 Dose (drug x line at 24 h with >= 3 doses): TAS per dose; minimum active dose (MAD, lowest dose with TAS >= 0.4);
  potency proxy = log10 MAD; primary vs secondary: correlation of the profile at the lowest active dose with the one
  at the highest dose (low r = the high-dose response is a different program: secondary / off-target); share of
  the highest-dose response explained by the low-dose profile (projection).
3 Anti-estrogens in MCF7: estrogen-response vs luminal-identity genes by time point and by dose.
Outputs: results/lincs/time_dose/{time_classes.tsv, dose_response.tsv, summary.txt}
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
P = ROOT / "data/processed/lincs"; OUT = ROOT / "results/lincs/time_dose"; OUT.mkdir(parents=True, exist_ok=True)
TAS, R_SUST = 0.4, 0.5
s = pd.read_csv(P / "level5_bing_hiq_sigs.tsv.gz", sep="\t")
g = pd.read_csv(P / "level5_bing_genes.tsv", sep="\t")
X = np.load(P / "level5_bing_hiq.npy", mmap_mode="r")
lm = np.where(g.feature_space.values == "landmark")[0]
s = s[s.ModelID.notna() & (s.pert_type == "trt_cp")].copy()
s["t"] = s.pert_time.round().astype(int); s["dose"] = s.pert_dose.round(3)
ann = pd.read_csv(ROOT / "results/lincs/drug_annotations.tsv", sep="\t", usecols=["pert_id", "hub_moa", "lincs_moa"]).drop_duplicates("pert_id")
lines = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", usecols=["ModelID", "CellLineName"]).set_index("ModelID").CellLineName
log = []
def say(x):
    print(x); log.append(x)

def prof(rows):  # mean landmark profile over signatures (replicate plates)
    return np.asarray(X[np.sort(rows)][:, lm], dtype=np.float64).mean(0)

# ---- 1 time
rows = []
for (pid, mid), d in s[s.t.isin([6, 24])].groupby(["pert_id", "ModelID"]):
    if set(d.t) != {6, 24}:
        continue
    shared = set(d[d.t == 6].dose) & set(d[d.t == 24].dose)
    pick = {t: (d[(d.t == t) & (d.dose == max(shared))] if shared else d[(d.t == t) & (d.dose == d[d.t == t].dose.max())]) for t in (6, 24)}
    a6, a24 = pick[6].tas.max(), pick[24].tas.max()
    p6, p24 = prof(pick[6].index.values), prof(pick[24].index.values)
    r = float(np.corrcoef(p6, p24)[0, 1])
    cls = ("early_sustained" if r >= R_SUST else "early_evolving") if a6 >= TAS and a24 >= TAS else \
          "late" if a24 >= TAS else "transient" if a6 >= TAS else "inactive"
    rows.append((pid, mid, d.cmap_name.iloc[0], max(shared) if shared else np.nan, a6, a24, r,
                 np.linalg.norm(p6) / np.sqrt(len(lm)), np.linalg.norm(p24) / np.sqrt(len(lm)), cls))
T = pd.DataFrame(rows, columns=["pert_id", "ModelID", "cmap_name", "shared_dose", "tas_6h", "tas_24h", "r_6_24", "size_6h", "size_24h", "time_class"])
T = T.merge(ann, on="pert_id", how="left"); T["line"] = T.ModelID.map(lines)
T.to_csv(OUT / "time_classes.tsv", sep="\t", index=False, float_format="%.4g")
say(f"== time: {len(T)} drug x line pairs with 6 h and 24 h ({T.pert_id.nunique()} drugs, {T.ModelID.nunique()} lines)")
say(T.time_class.value_counts().to_string())
act = T[T.time_class != "inactive"]
say(f"active pairs {len(act)}: early (6 h active) {act.time_class.str.startswith('early').mean():.0%} + transient "
    f"{(act.time_class == 'transient').mean():.0%}; late {(act.time_class == 'late').mean():.0%}; "
    f"median r(6 h, 24 h) where both active {T[T.time_class.str.startswith('early')].r_6_24.median():.2f}; "
    f"size 24 h / 6 h median {(T.size_24h / T.size_6h)[T.time_class.str.startswith('early')].median():.2f}")
moa = T.assign(moa=T.hub_moa.fillna(T.lincs_moa)).dropna(subset=["moa"])
moa = moa[moa.time_class != "inactive"].groupby("moa").agg(pairs=("pert_id", "size"), early=("time_class", lambda x: x.str.startswith("early").mean()),
                                                           late=("time_class", lambda x: (x == "late").mean()), r=("r_6_24", "median"))
say("mechanisms with >= 8 active pairs, by share late:\n" + moa[moa.pairs >= 8].sort_values("late", ascending=False).round(2).to_string())

# ---- 2 dose (24 h)
rows = []
for (pid, mid), d in s[s.t == 24].groupby(["pert_id", "ModelID"]):
    doses = np.sort(d.dose.unique())
    if len(doses) < 3:
        continue
    tas = np.array([d[d.dose == x].tas.max() for x in doses])
    if (tas < TAS).all():
        rows.append((pid, mid, d.cmap_name.iloc[0], len(doses), doses.min(), doses.max(), np.nan, tas.max(), np.nan, np.nan, "inactive"))
        continue
    lo = doses[np.argmax(tas >= TAS)]; hi = doses[-1]
    pl, ph = prof(d[d.dose == lo].index.values), prof(d[d.dose == hi].index.values)
    r = float(np.corrcoef(pl, ph)[0, 1]) if lo != hi else 1.0
    expl = float((pl @ ph) ** 2 / ((pl @ pl) * (ph @ ph))) if lo != hi else 1.0
    kind = "only top dose" if lo == hi else ("same program" if r >= R_SUST else "program changes with dose")
    rows.append((pid, mid, d.cmap_name.iloc[0], len(doses), doses.min(), hi, lo, tas.max(), r, expl, kind))
D = pd.DataFrame(rows, columns=["pert_id", "ModelID", "cmap_name", "n_doses", "dose_min", "dose_max", "min_active_dose", "tas_max",
                                "r_low_high", "share_explained", "dose_class"])
D["log10_mad"] = np.log10(D.min_active_dose)
D = D.merge(ann, on="pert_id", how="left"); D["line"] = D.ModelID.map(lines)
D.to_csv(OUT / "dose_response.tsv", sep="\t", index=False, float_format="%.4g")
say(f"\n== dose (24 h, >= 3 doses): {len(D)} drug x line pairs ({D.pert_id.nunique()} drugs)")
say(D.dose_class.value_counts().to_string())
a = D[D.dose_class.isin(["same program", "program changes with dose"])]
say(f"minimum active dose (uM) quartiles: {np.round(a.min_active_dose.quantile([.25, .5, .75]).values, 3).tolist()}; "
    f"median r(low, high) {a.r_low_high.median():.2f}")

# ---- 3 anti-estrogens in MCF7
canon = ["GREB1", "MYB", "XBP1", "PGR", "TFF1", "CA12", "SIAH2", "PDZK1", "CCND1", "IGFBP4", "NRIP1"]
lum = ["GATA3", "FOXA1", "ERBB4", "ANXA9", "HEBP1", "MINDY1", "SCCPDH", "SPDEF", "TBX3", "AR"]
gi = {x: i for i, x in enumerate(g.gene_symbol)}
ci, li = [gi[x] for x in canon if x in gi], [gi[x] for x in lum if x in gi]
ae = s[(s.ModelID == "ACH-000019") & s.cmap_name.isin(["fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene", "clomifene", "ZK-164015", "ICI-182780"])]
say("\n== anti-estrogens in MCF7 (mean z): estrogen-response vs luminal identity")
for key, grp in [("time", ae.groupby("t")), ("dose at 24 h (uM)", ae[ae.t == 24].groupby(pd.cut(ae[ae.t == 24].dose, [0, 0.1, 0.5, 2, 5, 100])))]:
    for k, d in grp:
        if len(d) < 3:
            continue
        v = np.asarray(X[np.sort(d.index.values)], dtype=np.float64)
        say(f"{key:>18} {str(k):>12} n={len(d):3d}  estrogen-response {v[:, ci].mean():+.2f}  luminal {v[:, li].mean():+.2f}  active {np.mean(d.tas >= TAS):.0%}")
(OUT / "summary.txt").write_text("\n".join(log) + "\n")
