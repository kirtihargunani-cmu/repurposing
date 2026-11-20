#!/usr/bin/env python
"""Target-agnostic candidate mechanisms: which regulators does a drug MOVE, and is any of them
disease-relevant in the opposing direction?

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

MOTIVATION (user, 2026-11-20). The pipeline scores drugs through their ANNOTATED primary target
(ChEMBL via D121), so a drug acting through a secondary or unannotated target is invisible by
construction. Half the drugs (158 of 317) have >= 2 annotated targets and 127 secondary units are
down-weighted by relative binding (median r_u 0.24); the rest have nothing annotated beyond the
primary. The goal is discovery - a drug not already used for a disease, even at lower specificity
or through a secondary upregulated target, is the finding the pipeline exists to make.

So the target never enters the question. Per (cohort, drug, window) and per regulator R:

    contribution_R = -(activity_R  x  disease_NES_R)

the per-regulator term of 79e's reversal score. A LARGE POSITIVE contribution means the drug moves R
opposite to the disease: R is up in the tumour and the drug pushes it down, or vice versa. The top
contributors are the candidate mechanisms, annotated or not.

THE PROLIFERATION PROBLEM IS HANDLED IN THE DESIGN, NOT AFTER. Every cytostatic drug lowers growth
regulators in every tumour, which is why an earlier whole-network attempt returned PI3K/mTOR
inhibitors as the top hit in every cohort. Each candidate therefore carries `prolif_linked` from 79h
(P1 enrichment for HALLMARK_E2F_TARGETS/G2M, or P2 activity tracking the proliferation score), and
the summary reports hits with and without them separately. A proliferation-driven hit is not
suppressed - it is labelled, because for some diseases growth arrest IS the mechanism.

Each candidate is also labelled `annotated` (the regulator is one of the drug's ChEMBL target genes)
or `unannotated` - the latter being the secondary mechanisms this is built to surface.

Output: results/systems_validation/integrated/candidate_mechanisms.tsv.gz + run json
"""
import json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
PROC = ROOT / "data/processed/lincs"
VDIR = PROC / "drug_viper_specific_v2_tf_cond"
OUT = ROOT / "results/systems_validation/integrated"
TOP_N = 10
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"


def net_for(cohort, ovr, avail):
    for c in ([ovr[cohort]] if cohort in ovr else []) + \
             (["COADREAD"] if cohort in ("TCGA-COAD", "TCGA-READ") else []) + \
             [cohort.split("-", 1)[1], cohort.split("-", 1)[1].split("_")[0]]:
        if c in avail:
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
    units = pd.read_csv(ROOT / "results/lincs/target_role/d121/units.tsv", sep="\t",
                        dtype={"tid": str, "action_flag": str})
    tgt = units.groupby("cmap_name").genes.apply(
        lambda s: {g for x in s for g in str(x).split("|")}).to_dict()
    mr = pd.read_csv(ROOT / "results/step6/aracne_project_v2_tf/group_mr.tsv", sep="\t")
    lw = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    ovr = pd.read_csv(ROOT / "config/group_networks.tsv", sep="\t").set_index("group_id").net.to_dict()
    pro = pd.read_csv(ROOT / "results/lincs/prolif_regulators_p2_v2_consensus.tsv", sep="\t")
    pro = pro[pro.set == "aracne"].set_index(["net", "regulator"]).prolif_linked.to_dict()
    avail = {p.stem for p in VDIR.glob("*.f32")}

    rows = []
    for cohort in sorted(mr.group_id.unique()):
        net = net_for(cohort, ovr, avail)
        if net is None:
            continue
        regs = (VDIR / f"{net}.regulators.txt").read_text().split()
        A = np.fromfile(VDIR / f"{net}.f32", dtype=np.float32).reshape(len(idx), len(regs)).astype(np.float64)
        rpos = {r: i for i, r in enumerate(regs)}
        m = mr[(mr.group_id == cohort) & mr.regulator.isin(rpos)]
        if len(m) < 10:
            continue
        ci = np.array([rpos[r] for r in m.regulator])
        s = m.nes.values.astype(np.float64)
        rn = m.regulator.values
        wg = lw[(lw.group_id == cohort) & lw.in_key].set_index("ModelID").weight
        sub = idx[idx.ModelID.isin(wg.index) & (idx.window != "other")]
        for (pid, win), g in sub.groupby(["pert_id", "window"]):
            ww = wg.reindex(g.ModelID).fillna(0).values
            ww = ww if ww.sum() > 0 else np.ones(len(g))
            D = (ww[:, None] * A[g.col.values][:, ci]).sum(0) / ww.sum()
            contrib = -(D * s)                       # per-regulator opposition
            top = np.argsort(contrib)[::-1][:TOP_N]
            drug = nm.get(pid)
            tg = tgt.get(drug, set())
            for rank, j in enumerate(top, 1):
                r = rn[j]
                rows.append(dict(cohort=cohort, network=net, pert_id=pid, cmap_name=drug,
                                 window=win, rank=rank, regulator=r,
                                 contribution=round(float(contrib[j]), 4),
                                 disease_nes=round(float(s[j]), 3),
                                 drug_activity=round(float(D[j]), 3),
                                 disease_direction="up" if s[j] > 0 else "down",
                                 prolif_linked=bool(pro.get((net, r), False)),
                                 annotated=bool(r in tg),
                                 n_conditions=len(g)))
        print(f"{cohort:22s} {net:12s} rows {len(rows):,}  [{time.time()-t0:.0f}s]", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "candidate_mechanisms.tsv.gz", sep="\t", index=False)
    top1 = df[df["rank"] == 1]
    res = dict(step="target-agnostic candidate mechanisms", top_n=TOP_N,
               n_rows=int(len(df)), n_cohorts=int(df.cohort.nunique()),
               n_drug_cohort_window=int(len(top1)),
               top1_prolif_linked=round(float(top1.prolif_linked.mean()), 4),
               top1_annotated=round(float(top1.annotated.mean()), 4),
               top1_annotated_excl_prolif=round(float(top1[~top1.prolif_linked].annotated.mean()), 4),
               any_top10_annotated=round(float(df.groupby(["cohort","pert_id","window"]).annotated.any().mean()), 4),
               disease_direction_of_top1=dict(top1.disease_direction.value_counts().items()))
    (OUT / "candidate_mechanisms_run.json").write_text(json.dumps(res, indent=1, default=str))
    print(f"\n{len(df):,} candidate rows; {len(top1):,} (cohort, drug, window) triples\n")
    print(f"  top-1 candidate is proliferation-linked : {100*res['top1_prolif_linked']:.1f}%")
    print(f"  top-1 candidate is an ANNOTATED target  : {100*res['top1_annotated']:.1f}%")
    print(f"     ... excluding proliferation-linked   : {100*res['top1_annotated_excl_prolif']:.1f}%")
    print(f"  any of the top {TOP_N} is annotated          : {100*res['any_top10_annotated']:.1f}%")
    print(f"  disease direction of the top-1 regulator: {res['disease_direction_of_top1']}")


if __name__ == "__main__":
    sys.exit(main())
