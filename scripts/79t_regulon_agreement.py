#!/usr/bin/env python
"""Drug-side agreement between the tumor-learned (ARACNe TF-only, v2) and curated (CollecTRI) regulons, per regulator.

Question (user 2026-11-02, after D83 / C25): how many regulators behave like ESR1, i.e. are read in opposite directions
on drug profiles by the two regulon sets? Diagnostic only; no list or score changes.
Input: data/processed/lincs/drug_viper_specific_v2_tf/<NET>.f32 and drug_viper_specific_collectri/ALL.f32 (79d), the
same 9,678 drug x line profiles; DoRothEA ABC (drug_viper_specific_dorothea_ABC) as a curated-vs-curated baseline.
Classes fixed before looking (per network x regulator scored by both):
  rho_all     Spearman of ARACNe vs CollecTRI NES over all profiles. disagree < 0; weak 0-0.3; agree >= 0.3
  ESR1-type   on profiles where CollecTRI reads the regulator strongly (|NES| >= 2, >= 30 such profiles), the ARACNe
              NES has the opposite sign in > 50%
Joined: D82 identity share of the regulon (histology labels, results/step5/identity_edges_v2_tf/regulon_summary.tsv).
Outputs: results/lincs/regulon_agreement_v2_tf.tsv (per network x regulator), regulon_agreement_v2_tf_summary.tsv
"""
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed/lincs"
OUT = ROOT / "results/lincs"
STRONG, MIN_STRONG = 2.0, 30
N_PROF = len(pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t"))


def load(d, net):
    regs = (d / f"{net}.regulators.txt").read_text().split()
    m = np.fromfile(d / f"{net}.f32", dtype=np.float32).reshape(N_PROF, len(regs))  # profile-major
    return pd.DataFrame(m, columns=regs)



def compare(a, b):
    """per shared regulator: Spearman over all profiles; sign disagreement where b is strong."""
    shared = sorted(set(a.columns) & set(b.columns))
    ra, rb = a[shared].rank(), b[shared].rank()
    rho = ra.corrwith(rb)
    rows = []
    for r in shared:
        s = np.abs(b[r].values) >= STRONG
        n = int(s.sum())
        opp = float(np.mean(np.sign(a[r].values[s]) != np.sign(b[r].values[s]))) if n else np.nan
        rows.append((r, rho[r], n, opp, float(np.median(np.abs(a[r].values)))))
    return pd.DataFrame(rows, columns=["regulator", "rho_all", "n_strong", "frac_opposite_strong", "median_abs_nes_a"])


def main():
    col = load(PROC / "drug_viper_specific_collectri", "ALL")
    dor = load(PROC / "drug_viper_specific_dorothea_ABC", "ALL")
    ad = PROC / "drug_viper_specific_v2_tf"
    nets = sorted(p.stem for p in ad.glob("*.f32"))
    res = []
    for net in nets:
        t = compare(load(ad, net), col)
        t.insert(0, "network", net)
        res.append(t)
    res = pd.concat(res, ignore_index=True)
    base = compare(dor, col)
    base.columns = ["regulator"] + [f"dorothea_{c}" for c in base.columns[1:]]
    res = res.merge(base[["regulator", "dorothea_rho_all", "dorothea_frac_opposite_strong"]], on="regulator", how="left")
    idn = pd.read_csv(ROOT / "results/step5/identity_edges_v2_tf/regulon_summary.tsv", sep="\t")
    idn = idn[idn.label_set == "histology"][["network", "regulator", "frac_identity"]]
    res = res.merge(idn, on=["network", "regulator"], how="left")

    res["class"] = pd.cut(res.rho_all, [-np.inf, 0, 0.3, np.inf], right=False, labels=["disagree", "weak", "agree"])
    res["esr1_type"] = ((res.n_strong >= MIN_STRONG) & (res.frac_opposite_strong > 0.5)).astype("boolean")
    res.loc[res.n_strong < MIN_STRONG, "esr1_type"] = pd.NA
    res.to_csv(OUT / "regulon_agreement_v2_tf.tsv", sep="\t", index=False, float_format="%.4g")

    print(f"network x regulator pairs scored by both: {len(res)} ({res.regulator.nunique()} regulators, {len(nets)} networks)")
    print("class counts:\n", res["class"].value_counts().reindex(["disagree", "weak", "agree"]).to_string())
    ev = res[res.n_strong >= MIN_STRONG]
    print(f"pairs with >= {MIN_STRONG} strong CollecTRI profiles: {len(ev)}; ESR1-type: {int(ev.esr1_type.sum())}")
    print(f"median rho ARACNe vs CollecTRI {res.rho_all.median():.2f}; DoRothEA vs CollecTRI (baseline) "
          f"{base.dorothea_rho_all.median():.2f} over {len(base)} regulators")
    per_tf = res.groupby("regulator").agg(
        n_nets=("network", "size"), n_disagree=("class", lambda c: int((c == "disagree").sum())),
        n_esr1_type=("esr1_type", lambda v: int(v.fillna(False).astype(bool).sum())),
        median_rho=("rho_all", "median"), dorothea_rho=("dorothea_rho_all", "first"),
        median_identity=("frac_identity", "median")).reset_index().sort_values("median_rho")
    per_tf.to_csv(OUT / "regulon_agreement_v2_tf_summary.tsv", sep="\t", index=False, float_format="%.4g")
    print("\nregulators by networks where they disagree (rho < 0):")
    print(per_tf.n_disagree.gt(0).sum(), "regulators disagree in >= 1 network;",
          (per_tf.n_disagree >= per_tf.n_nets / 2).sum(), "in at least half their networks")
    print("\nlowest 25 regulators by median rho:\n", per_tf.head(25).round(2).to_string(index=False))
    print("\nESR1 rows:\n", res[res.regulator == "ESR1"].round(2).to_string(index=False))
    # identity share vs agreement
    j = res.dropna(subset=["frac_identity"])
    if len(j):
        j = j.assign(id_bin=pd.cut(j.frac_identity, [-0.01, 0.1, 0.3, 1.0], labels=["<10%", "10-30%", ">=30%"]))
        print("\nagreement by D82 identity share (histology networks):\n",
              j.groupby("id_bin", observed=True).agg(n=("rho_all", "size"), median_rho=("rho_all", "median"),
                                                     share_disagree=("class", lambda c: (c == "disagree").mean())).round(2))


if __name__ == "__main__":
    main()
