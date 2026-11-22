#!/usr/bin/env python
"""Do top-ranked novel drugs share the mechanism profile of a disease's approved drugs? (D184)

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
SPEC: docs/SPEC_approved_profile_concordance.md, including the 2026-11-22 AMENDMENT (section 12) that
withdrew the pooled design. PRE-REGISTERED; every threshold, set definition, null and metric below comes
from that document and none is tuned to a result.

WINDOW-MATCHED. The unit is (cohort, window), window in {early, late}, and the two are reported separately
and never pooled. A pooled profile is the mean over whatever timepoints a drug happened to be measured at,
which is not the same quantity across drugs - for a 24 h-only drug it IS the 24 h profile - and D183 showed
~80% of a window's reliable variance is absent from the other. Timepoint composition would therefore be an
unbalanced covariate between the sets.

Sets, per (cohort, window), all requiring a profile in THAT window and `rev > 0`:
  A  approved for this cohort: label_phase >= 3 (D181 labels)            [secondary: label_phase == 4]
  T  top TOP_PCT by layered_pct, with NO label for this cohort
  N  the null: unlabelled AND not in T
Minimum sizes per spec: |A| >= 3 and |T| >= 10, else the cell is dropped and listed.

Statistic: conc(d) = cosine(x_d, centroid(A \\ {d})) on L2-normalised regulator contribution vectors,
leave-one-out so an approved drug never scores against a centroid containing itself.
Effect: median conc(T) - median conc(N), per cell; pooled across cells by the median of per-cell effects.
Null: N_PERM permutations shuffling the T/N assignment WITHIN a cell. conc does not change under the
shuffle (the centroid is fixed by A), so labels are permuted and the statistic recomputed.

Reported twice per spec section 7: raw, and after projecting out the cell's leading principal component
(the generic shared axis, 25.3% of variance in the LumA probe). A conclusion counts only if it holds in both.

Section 6 MOA control, the part that makes the test worth running: the whole thing repeated on the subset of
T sharing NO annotated target and NO MOA class with any member of A. Surviving that is mechanistic
convergence by a different route; not surviving means the test rediscovered drug classes, which target
annotation alone would have given.

Env: TOP_PCT (10) | N_PERM (10000) | SEED (0) | MIN_A (3) | MIN_T (10)
Output: results/lincs/approved_concordance/{per_drug,per_cell,acceptance}.tsv, run.json
"""
import json
import os
import re

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
SFX = "_v2_tf"
PROF = RES / f"regulon_contributions{SFX}" / "regulator_contributions_by_window_all.tsv.gz"
OUT = RES / "approved_concordance"
TOP_PCT = float(os.environ.get("TOP_PCT", "10"))
N_PERM = int(os.environ.get("N_PERM", "10000"))
SEED = int(os.environ.get("SEED", "0"))
MIN_A = int(os.environ.get("MIN_A", "3"))
MIN_T = int(os.environ.get("MIN_T", "10"))
# pre-registered acceptance (spec section 8, criterion 2 restated in section 12 for the window-matched n)
ACC_P = 0.01
ACC_BREADTH = 0.60
ACC_MOA_P = 0.05
norm_name = lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower())


def tokens(s):
    return {t for t in re.split(r"[|,;/]+", str(s or "")) if t and t.lower() not in ("nan", "none")}


def load():
    prof = pd.read_csv(PROF, sep="\t", low_memory=False,
                       usecols=["cohort", "pert_id", "cmap_name", "window", "regulator", "rev_contrib", "rev"])
    lab = pd.read_csv(RES / "labels/pair_labels_v2_tf.tsv.gz", sep="\t", low_memory=False)
    lab = lab[lab.unit_kind.eq("group")][["unit_id", "pert_id", "label_phase"]].rename(columns={"unit_id": "cohort"})
    lab = lab.groupby(["cohort", "pert_id"], as_index=False).label_phase.max()
    sp = pd.read_csv(RES / f"drug_reversal_split{SFX}.tsv", sep="\t", low_memory=False,
                     usecols=["group_id", "pert_id", "layered_pct", "passes_gate", "evidence_targets"])
    sp = sp.rename(columns={"group_id": "cohort"}).drop_duplicates(["cohort", "pert_id"])
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t", low_memory=False).drop_duplicates("pert_id").set_index("pert_id")
    return prof, lab, sp, ann


def cell_sets(meta):
    """A / T / N for one (cohort, window) cell, per spec section 3 + amendment."""
    m = meta[meta.rev > 0].copy()
    A = m[m.label_phase >= 3].pert_id.tolist()
    unl = m[m.label_phase == 0]
    if not len(unl):
        return A, [], []
    cut = np.nanpercentile(unl.layered_pct.astype(float), 100 - TOP_PCT) if unl.layered_pct.notna().any() else np.inf
    T = unl[unl.layered_pct >= cut].pert_id.tolist()
    N = unl[~unl.pert_id.isin(T)].pert_id.tolist()
    return A, T, N


def conc_vector(X, ids, A_idx):
    """cosine of each row to the centroid of A, leave-one-out for A members."""
    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)
    if len(A_idx) == 0:
        return np.full(len(ids), np.nan)
    S = Xn[A_idx].sum(0)
    out = np.empty(len(ids))
    inA = np.zeros(len(ids), bool); inA[A_idx] = True
    for i in range(len(ids)):
        c = (S - Xn[i]) if inA[i] else S
        n = len(A_idx) - (1 if inA[i] else 0)
        if n <= 0:
            out[i] = np.nan; continue
        c = c / n
        out[i] = float(Xn[i] @ c / (np.linalg.norm(c) + 1e-12))
    return out


def drop_pc1(X):
    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)
    C = Xn - Xn.mean(0)
    try:
        _, _, vt = np.linalg.svd(C, full_matrices=False)
    except np.linalg.LinAlgError:
        return X
    v = vt[0]
    return Xn - np.outer(Xn @ v, v)


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    prof, lab, sp, ann = load()
    rng = np.random.default_rng(SEED)
    rows, cells, dropped = [], [], []

    for (cohort, window), g in prof.groupby(["cohort", "window"]):
        M = g.pivot_table(index="pert_id", columns="regulator", values="rev_contrib", aggfunc="first")
        X = np.nan_to_num(M.values)
        ids = list(M.index)
        meta = (pd.DataFrame({"pert_id": ids})
                .merge(g.drop_duplicates("pert_id")[["pert_id", "cmap_name", "rev"]], on="pert_id", how="left")
                .merge(lab[lab.cohort == cohort][["pert_id", "label_phase"]], on="pert_id", how="left")
                .merge(sp[sp.cohort == cohort][["pert_id", "layered_pct", "evidence_targets"]], on="pert_id", how="left"))
        meta["label_phase"] = meta.label_phase.fillna(0)
        A, T, N = cell_sets(meta)
        if len(A) < MIN_A or len(T) < MIN_T:
            dropped.append(dict(cohort=cohort, window=window, n_A=len(A), n_T=len(T), n_N=len(N),
                                reason="|A|<MIN_A" if len(A) < MIN_A else "|T|<MIN_T"))
            continue
        pos = {p: i for i, p in enumerate(ids)}
        A_idx = [pos[p] for p in A]
        conc_raw = conc_vector(X, ids, A_idx)
        conc_pc1 = conc_vector(drop_pc1(X), ids, A_idx)

        # MOA / target disjointness against A
        tg = meta.set_index("pert_id").evidence_targets.to_dict()
        A_tg, A_moa = set(), set()
        for p in A:
            A_tg |= tokens(tg.get(p))
            A_moa |= {norm_name(x) for x in tokens(ann.hub_moa.get(p))}
        disj = {}
        for p in ids:
            dt = tokens(tg.get(p))
            dm = {norm_name(x) for x in tokens(ann.hub_moa.get(p))}
            disj[p] = (not (dt & A_tg)) and (not (dm & A_moa))

        setof = {p: ("A" if p in set(A) else "T" if p in set(T) else "N" if p in set(N) else "-") for p in ids}
        for i, p in enumerate(ids):
            rows.append(dict(cohort=cohort, window=window, pert_id=p, cmap_name=meta.cmap_name.iloc[i],
                             set=setof[p], conc=conc_raw[i], conc_pc1removed=conc_pc1[i],
                             moa_disjoint=bool(disj[p]), layered_pct=meta.layered_pct.iloc[i],
                             rev=meta.rev.iloc[i], label_phase=meta.label_phase.iloc[i]))

        ti = np.array([pos[p] for p in T]); ni = np.array([pos[p] for p in N])
        rec = dict(cohort=cohort, window=window, n_A=len(A), n_T=len(T), n_N=len(N),
                   ceiling=float(np.nanmedian(conc_raw[A_idx])))
        for tag, cv in (("raw", conc_raw), ("pc1removed", conc_pc1)):
            if len(ni) < 5:
                rec[f"effect_{tag}"] = np.nan; continue
            rec[f"medT_{tag}"] = float(np.nanmedian(cv[ti]))
            rec[f"medN_{tag}"] = float(np.nanmedian(cv[ni]))
            rec[f"effect_{tag}"] = rec[f"medT_{tag}"] - rec[f"medN_{tag}"]
            pool = np.concatenate([ti, ni]); k = len(ti)
            perm = np.array([np.nanmedian(cv[s[:k]]) - np.nanmedian(cv[s[k:]])
                             for s in (rng.permutation(pool) for _ in range(N_PERM))])
            rec[f"perm_{tag}"] = perm          # kept for the pooled test below
        # MOA-disjoint arm
        td = np.array([pos[p] for p in T if disj[p]])
        rec["n_T_moa_disjoint"] = len(td)
        if len(td) >= 5 and len(ni) >= 5:
            rec["effect_moa_disjoint"] = float(np.nanmedian(conc_raw[td]) - np.nanmedian(conc_raw[ni]))
            pool = np.concatenate([td, ni]); k = len(td)
            rec["perm_moa"] = np.array([np.nanmedian(conc_raw[s[:k]]) - np.nanmedian(conc_raw[s[k:]])
                                        for s in (rng.permutation(pool) for _ in range(N_PERM))])
        else:
            rec["effect_moa_disjoint"] = np.nan
        cells.append(rec)
        print(f"  {cohort:22s} {window:5s} |A|={len(A):>3} |T|={len(T):>4} |N|={len(N):>4} "
              f"effect {rec.get('effect_raw', float('nan')):+.4f}", flush=True)

    P = pd.DataFrame(rows)
    P.to_csv(OUT / "per_drug.tsv", sep="\t", index=False, float_format="%.6g")
    C = pd.DataFrame([{k: v for k, v in c.items() if not isinstance(v, np.ndarray)} for c in cells])
    C.to_csv(OUT / "per_cell.tsv", sep="\t", index=False, float_format="%.6g")

    # pooled permutation test, per window: median of per-cell effects vs the same statistic under permutation
    acc = []
    for window in sorted(C.window.unique()):
        sel = [c for c in cells if c["window"] == window]
        row = dict(window=window, n_cells=len(sel))
        for tag, key, pk in (("raw", "effect_raw", "perm_raw"), ("pc1removed", "effect_pc1removed", "perm_pc1removed"),
                             ("moa_disjoint", "effect_moa_disjoint", "perm_moa")):
            eff = np.array([c.get(key, np.nan) for c in sel], float)
            use = [c for c in sel if pk in c and np.isfinite(c.get(key, np.nan))]
            if not use:
                row[f"effect_{tag}"] = np.nan; row[f"p_{tag}"] = np.nan; row[f"breadth_{tag}"] = np.nan
                continue
            obs = float(np.nanmedian([c[key] for c in use]))
            null = np.median(np.vstack([c[pk] for c in use]), axis=0)
            row[f"effect_{tag}"] = round(obs, 5)
            row[f"p_{tag}"] = round(float((null >= obs).mean()), 5)
            row[f"breadth_{tag}"] = round(float(np.nanmean(eff[np.isfinite(eff)] > 0)), 4)
            row[f"n_{tag}"] = len(use)
        row["ceiling"] = round(float(np.nanmedian([c["ceiling"] for c in sel])), 4)
        row["c1_effect_p"] = bool(np.isfinite(row.get("p_raw", np.nan)) and row["effect_raw"] > 0 and row["p_raw"] < ACC_P)
        row["c2_breadth"] = bool(np.isfinite(row.get("breadth_raw", np.nan)) and row["breadth_raw"] >= ACC_BREADTH)
        row["c3_pc1"] = bool(np.isfinite(row.get("p_pc1removed", np.nan)) and row["effect_pc1removed"] > 0
                             and row["p_pc1removed"] < ACC_P)
        row["c4_moa"] = bool(np.isfinite(row.get("p_moa_disjoint", np.nan)) and row["effect_moa_disjoint"] > 0
                             and row["p_moa_disjoint"] < ACC_MOA_P)
        row["VERDICT"] = ("PASS" if all(row[k] for k in ("c1_effect_p", "c2_breadth", "c3_pc1", "c4_moa"))
                          else "CLASS_RECOVERY" if (row["c1_effect_p"] and row["c2_breadth"] and row["c3_pc1"])
                          else "NULL")
        acc.append(row)
    A_ = pd.DataFrame(acc)
    A_.to_csv(OUT / "acceptance.tsv", sep="\t", index=False)
    print("\n" + A_.to_string(index=False), flush=True)
    if dropped:
        print(f"\ndropped cells ({len(dropped)}): " + ", ".join(f"{d['cohort']}/{d['window']}({d['reason']})" for d in dropped))
    json.dump(dict(spec="docs/SPEC_approved_profile_concordance.md section 12 (window-matched)",
                   top_pct=TOP_PCT, n_perm=N_PERM, seed=SEED, min_A=MIN_A, min_T=MIN_T,
                   thresholds=dict(p=ACC_P, breadth=ACC_BREADTH, moa_p=ACC_MOA_P),
                   cells_used=len(cells), cells_dropped=dropped,
                   verdicts=A_[["window", "VERDICT"]].to_dict("records"),
                   runtime_s=round(time.time() - t0, 1)), open(OUT / "run.json", "w"), indent=1)
    print(f"\n-> {OUT}  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
