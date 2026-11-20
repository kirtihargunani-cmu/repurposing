#!/usr/bin/env python
"""Route 1 as the known-answer set for Route 2, across four proliferation arms.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

User, 2026-11-20: "lets also run route 1, it can be used for validation of route 2. we can validate
if route 2 after removing proliferation uncovers these known secondary targets", and
"should we remove proliferation before or after".

ROUTE 2 (`candidate_mechanisms.py`) is target-agnostic: per (cohort, drug, window) it ranks every
regulator by how much it contributes to opposing the disease, `contrib = -(D_R * s_R)`, and calls the
top contributors candidate mechanisms. On the TF-only network it returned 0.0% annotated top-1
candidates - which is structural, not a bug: only 16 of 387 annotated target genes (4.1%) are
regulators there, because scored regulators are transcription factors while drug targets are kinases
and receptors. On the FULL ARACNe set 230 of 387 (59.4%) are regulators, including 111 of 126
secondary-target genes (88%).

ROUTE 1 is the known answer: where does Route 2 rank a target we already know the drug binds? With
187 drug x scoreable-secondary-target pairs over 817 cohort-weighted pairs, recovery can be measured
against the ~0.2% expected by chance from a ~4,600-regulator list.

THE FOUR ARMS answer "remove proliferation before or after", by measuring rather than arguing:
  none       no removal - the current Route 2
  universal  A - beta * mean(pi over that drug's lines in the window)   (D157: drops the drug's
             general cytostatic effect, keeps line-selective growth arrest)
  full       A - beta * pi                                             (D154/D157 `full`)
  full_filter  `full`, then drop prolif_linked regulators from the ranking
`none` -> `full` is the "subtract" axis; `full` -> `full_filter` is the "filter" axis. If recovery
RISES with removal, growth arrest was crowding out real mechanism - the user's hypothesis. If it
FALLS, what we call proliferation was carrying genuine target signal.

The adjustment is linear in pi, so no adjusted matrix is ever materialised:
    weighted_mean(A - outer(pi, beta)) = weighted_mean(A) - weighted_mean(pi) * beta
beta is fitted per (network, regulator, window) exactly as 79i does (D154: per window, because a
pooled beta changes sign for 25.7% of regulators).

Output: results/systems_validation/integrated/route1_validation.{tsv.gz,json}
"""
import json, os, time
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
PROC = ROOT / "data/processed/lincs"
VDIR = PROC / Path(os.environ.get("VDIR", "drug_viper_specific_v2_condfull"))
MRF = ROOT / os.environ.get("MR", "results/step6/aracne_project_v2/group_mr.tsv")
P1F = ROOT / os.environ.get("P1", "results/lincs/prolif_regulators_v2.tsv")
OUT = ROOT / "results/systems_validation/integrated"
TOP_N = 10
MIN_WIN_N = 40
CHUNK = 512
ARMS = ["none", "universal", "full", "full_filter"]
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "6"


def net_for(cohort, ovr, avail):
    for c in ([ovr[cohort]] if cohort in ovr else []) + \
             (["COADREAD"] if cohort in ("TCGA-COAD", "TCGA-READ") else []) + \
             [cohort.split("-", 1)[1], cohort.split("-", 1)[1].split("_")[0]]:
        if c in avail:
            return c
    return None


def window_betas(A, pi, win, regs_n):
    """beta per regulator per window, fitted exactly as 79i (D154). Chunked over regulators."""
    out = {}
    for w in ("early", "late"):
        m = (win == w).values
        if m.sum() < MIN_WIN_N:
            out[w] = None
            continue
        p = pi[m].astype(np.float64)
        pc = p - p.mean()
        den = pc @ pc
        if den <= 0:
            out[w] = None
            continue
        b = np.empty(regs_n)
        for a in range(0, regs_n, CHUNK):
            blk = A[m, a:a + CHUNK].astype(np.float64)
            b[a:a + CHUNK] = (blk - blk.mean(0)).T @ pc / den
        out[w] = b
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    idx = pd.read_csv(PROC / "drug_condition_profiles_index.tsv", sep="\t")
    idx["col"] = np.arange(len(idx))
    idx["window"] = np.where(idx.pert_time <= 6, "early",
                             np.where(idx.pert_time >= 24, "late", "other"))
    genes = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
    gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:]
           for l in open(ROOT / "data/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt")}
    pm = genes.symbol.isin(set(gmt["HALLMARK_E2F_TARGETS"]) | set(gmt["HALLMARK_G2M_CHECKPOINT"])).values
    pi = np.fromfile(PROC / "drug_condition_profiles_specific.f32",
                     dtype=np.float32).reshape(len(idx), -1)[:, pm].mean(1).astype(np.float64)
    idx["pi"] = pi
    # D157 universal component: the drug's cross-line mean of pi within the window
    idx["pi_univ"] = idx.groupby(["pert_id", "window"]).pi.transform("mean")

    annot = pd.read_csv(ROOT / "results/lincs/drug_annotations.tsv", sep="\t")
    nm = dict(annot.drop_duplicates("pert_id")[["pert_id", "cmap_name"]].values)
    units = pd.read_csv(ROOT / "results/lincs/target_role/d121/units.tsv", sep="\t",
                        dtype={"tid": str, "action_flag": str})
    units["r_u"] = pd.to_numeric(units.r_u, errors="coerce")
    units["is_primary"] = units.r_u >= units.groupby("cmap_name").r_u.transform("max") - 1e-9
    gene2role = {}
    for _, u in units.iterrows():
        for g in str(u.genes).split("|"):
            k = (u.cmap_name, g)
            # a gene is "secondary" only if it is secondary in EVERY unit naming it
            gene2role[k] = gene2role.get(k, True) and bool(u.is_primary)
    tgt = {}
    for (d, g) in gene2role:
        tgt.setdefault(d, set()).add(g)

    mr = pd.read_csv(MRF, sep="\t")
    lw = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    ovr = pd.read_csv(ROOT / "config/group_networks.tsv", sep="\t").set_index("group_id").net.to_dict()
    p1 = pd.read_csv(P1F, sep="\t")
    p1 = p1[p1.prolif_linked].groupby("net").regulator.apply(set).to_dict()
    avail = {p.stem for p in VDIR.glob("*.f32")}
    print(f"networks available in {VDIR.name}: {len(avail)} -> {sorted(avail)}", flush=True)
    if not avail:
        raise SystemExit("no VIPER output yet")

    rows = []
    for cohort in sorted(mr.group_id.unique()):
        net = net_for(cohort, ovr, avail)
        if net is None:
            continue
        regs = (VDIR / f"{net}.regulators.txt").read_text().split()
        A = np.fromfile(VDIR / f"{net}.f32", dtype=np.float32).reshape(len(idx), len(regs))
        rpos = {r: i for i, r in enumerate(regs)}
        m = mr[(mr.group_id == cohort) & mr.regulator.isin(rpos)]
        if len(m) < 10:
            del A
            continue
        ci = np.array([rpos[r] for r in m.regulator])
        s = m.nes.values.astype(np.float64)
        rn = m.regulator.values
        plink = np.array([r in p1.get(net, set()) for r in rn])
        beta = window_betas(A, pi, idx.window, len(regs))
        wg = lw[(lw.group_id == cohort) & lw.in_key].set_index("ModelID").weight
        sub = idx[idx.ModelID.isin(wg.index) & (idx.window != "other")]
        nreg = len(rn)
        for (pid, win), g in sub.groupby(["pert_id", "window"]):
            drug = nm.get(pid)
            tg = {t for t in tgt.get(drug, set()) if t in set(rn)}
            if not tg:
                continue                      # Route 1 needs a scoreable annotated target
            ww = wg.reindex(g.ModelID).fillna(0).values
            ww = ww if ww.sum() > 0 else np.ones(len(g))
            ww = ww / ww.sum()
            D0 = (ww[:, None] * A[g.col.values][:, ci]).sum(0).astype(np.float64)
            b = beta.get(win)
            bci = b[ci] if b is not None else np.zeros(nreg)
            pbar_full = float(ww @ g.pi.values)
            pbar_univ = float(ww @ g.pi_univ.values)
            Ds = {"none": D0, "full": D0 - pbar_full * bci, "universal": D0 - pbar_univ * bci}
            Ds["full_filter"] = Ds["full"]
            for arm in ARMS:
                contrib = -(Ds[arm] * s)
                keep = ~plink if arm == "full_filter" else np.ones(nreg, bool)
                cc = np.where(keep, contrib, -np.inf)
                order = np.argsort(cc)[::-1]
                rank_of = np.empty(nreg, int); rank_of[order] = np.arange(1, nreg + 1)
                nk = int(keep.sum())
                for t in sorted(tg):
                    j = int(np.where(rn == t)[0][0])
                    if not keep[j]:
                        continue
                    rows.append(dict(cohort=cohort, network=net, pert_id=pid, cmap_name=drug,
                                     window=win, arm=arm, target=t,
                                     is_primary=bool(gene2role.get((drug, t), True)),
                                     rank=int(rank_of[j]), n_regulators=nk,
                                     contribution=round(float(contrib[j]), 4),
                                     disease_nes=round(float(s[j]), 3),
                                     prolif_linked=bool(plink[j]),
                                     n_conditions=len(g)))
        del A
        print(f"{cohort:22s} {net:12s} rows {len(rows):,}  [{time.time()-t0:.0f}s]", flush=True)

    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit("no Route 1 pairs scored")
    df.to_csv(OUT / "route1_validation.tsv.gz", sep="\t", index=False)

    res = {"n_rows": int(len(df)), "n_pairs": int(df.groupby(["cohort","pert_id","window","target"]).ngroups),
           "n_drugs": int(df.cmap_name.nunique()), "n_targets": int(df.target.nunique()),
           "n_cohorts": int(df.cohort.nunique()), "arms": {}}
    print(f"\n{df.cmap_name.nunique()} drugs x {df.target.nunique()} targets over "
          f"{df.cohort.nunique()} cohorts; {res['n_pairs']:,} (cohort, drug, window, target) pairs\n")
    hdr = f"{'arm':<13}{'n':>7}{'chance@10':>11}{'rec@10':>9}{'lift':>7}{'rec@50':>9}{'med rank':>10}{'med pct':>9}"
    print(hdr); print("-" * len(hdr))
    for arm in ARMS:
        for sub_name, d in (("all", df[df.arm == arm]),
                            ("secondary", df[(df.arm == arm) & ~df.is_primary])):
            if not len(d):
                continue
            ch = float((TOP_N / d.n_regulators).mean())
            r10 = float((d["rank"] <= TOP_N).mean()); r50 = float((d["rank"] <= 50).mean())
            pct = float((d["rank"] / d.n_regulators).median())
            key = arm if sub_name == "all" else f"{arm}:secondary"
            res["arms"][key] = dict(n=int(len(d)), chance_at_10=round(ch, 5),
                                    recovery_at_10=round(r10, 5), lift=round(r10 / ch, 2) if ch else None,
                                    recovery_at_50=round(r50, 5), median_rank=float(d["rank"].median()),
                                    median_rank_pct=round(pct, 4))
            lab = arm if sub_name == "all" else f"  `-secondary"
            print(f"{lab:<13}{len(d):>7,}{100*ch:>10.2f}%{100*r10:>8.2f}%{r10/ch if ch else 0:>7.1f}"
                  f"{100*r50:>8.2f}%{d['rank'].median():>10.0f}{100*pct:>8.1f}%")
    json.dump(res, open(OUT / "route1_validation.json", "w"), indent=2)
    print(f"\n[{time.time()-t0:.0f}s] -> {OUT/'route1_validation.tsv.gz'}")


if __name__ == "__main__":
    main()
