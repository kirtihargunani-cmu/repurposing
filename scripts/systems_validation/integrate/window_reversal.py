#!/usr/bin/env python
"""Per-window drug reversal scores: score_early and score_late, per drug x disease.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
SPEC: SPEC_condition_level_profiles.md section 11.

Why a separate script rather than a flag on 79e: 79e's aggregation assumes one row per
(drug, line), and that assumption is threaded through it - four attempts to retrofit a window
filter produced four distinct failures (reshape on a filtered index, get_group on drugs absent from
a window, an integer cast on NaN counts). Rather than keep patching a production scorer, this
reimplements only the scoring formula, which is small and explicit:

    rev = -(D . s) / ||s||          D = drug regulator activity, s = disease regulator NES

exactly as 79e:285, with the drug's line weighting from celligner as 79e does. 79e itself is
untouched and still produces production.

Per (cohort, drug, window): activity is the mean over that drug's conditions IN THAT WINDOW, in
cohort-matched lines, line-weighted. Windows: early <= 6 h, late >= 24 h.

Output: results/systems_validation/integrated/window_reversal.tsv.gz + run json
"""
import json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
PROC = ROOT / "data/processed/lincs"
VDIR = PROC / "drug_viper_specific_v2_tf_cond"
OUT = ROOT / "results/systems_validation/integrated"
MIN_MRS = 2
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"


def net_for(cohort, ovr, available):
    for c in ([ovr[cohort]] if cohort in ovr else []) + \
             (["COADREAD"] if cohort in ("TCGA-COAD", "TCGA-READ") else []) + \
             [cohort.split("-", 1)[1], cohort.split("-", 1)[1].split("_")[0]]:
        if c in available:
            return c
    return None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    idx = pd.read_csv(PROC / "drug_condition_profiles_index.tsv", sep="\t")
    idx["col"] = np.arange(len(idx))
    idx["window"] = np.where(idx.pert_time <= 6, "early",
                             np.where(idx.pert_time >= 24, "late", "other"))
    annot = pd.read_csv(ROOT / "results/lincs/drug_annotations.tsv", sep="\t")
    nm = dict(annot.drop_duplicates("pert_id")[["pert_id", "cmap_name"]].values)
    mr = pd.read_csv(ROOT / "results/step6/aracne_project_v2_tf/group_mr.tsv", sep="\t")
    lw = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    ovr = pd.read_csv(ROOT / "config/group_networks.tsv", sep="\t").set_index("group_id").net.to_dict()
    avail = {p.stem for p in VDIR.glob("*.f32")}
    print(f"{len(idx):,} conditions | {len(avail)} networks | "
          f"{int((idx.window=='early').sum()):,} early, {int((idx.window=='late').sum()):,} late")

    rows = []
    for cohort in sorted(mr.group_id.unique()):
        net = net_for(cohort, ovr, avail)
        if net is None:
            continue
        regs = (VDIR / f"{net}.regulators.txt").read_text().split()
        A = np.fromfile(VDIR / f"{net}.f32", dtype=np.float32).reshape(len(idx), len(regs)).astype(np.float64)
        rpos = {r: i for i, r in enumerate(regs)}
        m = mr[(mr.group_id == cohort) & mr.regulator.isin(rpos)]
        if len(m) < MIN_MRS:
            continue
        ci = np.array([rpos[r] for r in m.regulator])
        s = m.nes.values.astype(np.float64)
        norm = np.linalg.norm(s)
        wg = lw[(lw.group_id == cohort) & lw.in_key].set_index("ModelID").weight
        sub = idx[idx.ModelID.isin(wg.index) & (idx.window != "other")]
        if not len(sub):
            continue
        for (pid, win), g in sub.groupby(["pert_id", "window"]):
            ww = wg.reindex(g.ModelID).fillna(0).values
            ww = ww if ww.sum() > 0 else np.ones(len(g))
            D = (ww[:, None] * A[g.col.values][:, ci]).sum(0) / ww.sum()
            rows.append(dict(cohort=cohort, network=net, pert_id=pid, cmap_name=nm.get(pid),
                             window=win, n_conditions=len(g), n_lines=g.ModelID.nunique(),
                             n_doses=g.pert_dose.nunique(), n_mrs=len(m),
                             rev=round(float(-(D @ s) / norm), 5)))
        print(f"{cohort:22s} {net:12s} rows {len(rows):,}  [{time.time()-t0:.0f}s]", flush=True)

    df = pd.DataFrame(rows)
    # pivot to one row per (cohort, drug) carrying both windows, per spec section 11
    piv = df.pivot_table(index=["cohort", "network", "pert_id", "cmap_name"], columns="window",
                         values=["rev", "n_conditions", "n_doses"], aggfunc="first")
    piv.columns = [f"{a}_{b}" for a, b in piv.columns]
    piv = piv.reset_index()
    for w in ("early", "late"):
        for c in (f"rev_{w}", f"n_conditions_{w}", f"n_doses_{w}"):
            if c not in piv:
                piv[c] = np.nan
    piv["window_coverage"] = np.where(piv.rev_early.notna() & piv.rev_late.notna(), "both",
                             np.where(piv.rev_early.notna(), "early_only",
                             np.where(piv.rev_late.notna(), "late_only", "neither")))
    both = piv.window_coverage == "both"
    piv["window_agreement"] = np.where(~both, "not_evaluable",
                              np.where(np.sign(piv.rev_early) == np.sign(piv.rev_late),
                                       "concordant", "discordant"))
    piv.to_csv(OUT / "window_reversal.tsv.gz", sep="\t", index=False)

    b = piv[both]
    res = dict(step="per-window reversal score", formula="rev = -(D . s)/||s|| as 79e:285",
               n_rows=int(len(piv)), n_cohorts=int(piv.cohort.nunique()),
               coverage=dict(piv.window_coverage.value_counts().items()),
               n_both=int(both.sum()),
               spearman_early_late=round(float(b.rev_early.corr(b.rev_late, method="spearman")), 4),
               frac_discordant_of_both=round(float((piv[both].window_agreement == "discordant").mean()), 4),
               median_n_conditions=dict(early=float(piv.n_conditions_early.median()),
                                        late=float(piv.n_conditions_late.median())))
    (OUT / "window_reversal_run.json").write_text(json.dumps(res, indent=1, default=str))
    print(f"\n{len(piv):,} (cohort, drug) rows over {piv.cohort.nunique()} cohorts")
    print(f"  coverage: {res['coverage']}")
    print(f"  Spearman(rev_early, rev_late) over the {int(both.sum()):,} 'both' rows = {res['spearman_early_late']:+.4f}")
    print(f"  sign DISCORDANT between windows: {100*res['frac_discordant_of_both']:.1f}% of 'both' rows")
    print(f"  median conditions behind a score: early {res['median_n_conditions']['early']}, "
          f"late {res['median_n_conditions']['late']}")


if __name__ == "__main__":
    sys.exit(main())
