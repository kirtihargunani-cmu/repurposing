#!/usr/bin/env python
"""O-3 / D136-S1: gene-level BMDs, censoring-aware regulon aggregation, coherence, intrinsic null.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
SPEC: docs/SPEC_dose_response_bmd.md sections 3-6.  S0 (crossing curves) fixed k by the
pre-registered rule - the value putting the median crossing rate nearest 50% - which measured
k = 2.00 (median 0.5516). Every headline is nevertheless reported as a PROFILE across k, per
spec section 3: "a conclusion counts only if it holds throughout".

Level 0 (gene)    per (drug, line, gene): BMD = lowest dose at which the response crosses
                  BMR_g = k * sd_within_g, by linear interpolation between the bracketing doses.
                  Genes that never cross are RIGHT-CENSORED at the top tested dose, never dropped.
Level 1 (regulon) censoring-aware aggregate over member genes (Kaplan-Meier median), plus
                  COHERENCE = dispersion of member-gene BMDs, plus the censored fraction, always.
                  Genes moving against the regulon's expected sign go to DISCORDANCE, never the BMD.
Intrinsic null    random gene sets from the SAME drug's own profile, matched on set size and on the
                  sd_within distribution of the real regulon (spec section 6) - otherwise the null
                  is easier merely by containing quieter genes.

Level 2 (target footprint) is NOT implemented here; it needs the signed-path machinery and is the
next step.

Output: results/systems_validation/observability/{o3_s1_regulon_bmd.tsv, o3_s1_run.json}
"""
import json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd, h5py

ROOT = Path("/home/kirti/projects/repurposing")
RAW = ROOT / "data/raw/lincs2020"
OUT = ROOT / "results/systems_validation/observability"
GCTX = RAW / "level5_beta_trt_cp_n720216x12328.gctx"
SEED, N_CELLS, BLOCK = 20261118, 400, 20000
KS = [1.5, 2.0, 2.5, 3.0]          # profile; 2.0 is the S0-selected value
K_PRIMARY = 2.0
N_NULL = 50
MIN_GENES = 10
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "6"


def dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def km_median(bmd, censored, top):
    """Kaplan-Meier median of right-censored BMDs; NaN if survival never reaches 0.5."""
    t = np.where(censored, top, bmd)
    o = np.argsort(t)
    t, c = t[o], censored[o]
    n, s = len(t), 1.0
    for i in range(n):
        if not c[i]:
            s *= (1 - 1 / (n - i))
            if s <= 0.5:
                return float(t[i])
    return float("nan")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    q = pd.read_csv(ROOT / "results/lincs/gene_quality/bing_gene_quality.tsv", sep="\t",
                    dtype={"gene_id": str})
    sdw = q.set_index("gene_id").sd_within

    # ---- line -> ARACNe network, so a regulon can be attached to each cell --------------------
    models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", low_memory=False)
    lin = dict(zip(models.StrippedCellLineName.astype(str).str.upper(), models.OncotreeLineage))
    LIN2NET = {"Skin": "SKCM", "Breast": "BRCA", "Lung": "LUAD", "Bowel": "COADREAD",
               "Prostate": "PRAD", "Pancreas": "PAAD", "Liver": "LIHC", "Ovary/Fallopian Tube": "OV",
               "Kidney": "KIRC", "Esophagus/Stomach": "STAD", "Uterus": "UCEC", "Bladder/Urinary Tract": "BLCA",
               "CNS/Brain": "GBM", "Head and Neck": "HNSC", "Thyroid": "THCA", "Soft Tissue": "SARC",
               "Cervix": "CESC"}

    sig = pd.read_csv(RAW / "siginfo_beta.txt", sep="\t", low_memory=False,
                      usecols=["sig_id", "pert_id", "cmap_name", "pert_type", "cell_iname",
                               "pert_dose", "pert_time", "is_hiq", "tas"])
    sig = sig[(sig.pert_type == "trt_cp") & (sig.is_hiq == 1) & (sig.tas >= 0.4) & (sig.pert_dose > 0)]
    sig = sig.sort_values("tas", ascending=False).drop_duplicates(["pert_id", "cell_iname", "pert_dose"])
    sig["lineage"] = sig.cell_iname.astype(str).str.upper().map(lin)
    sig["net"] = sig.lineage.map(LIN2NET)
    sig = sig[sig.net.notna()]
    nd = sig.groupby(["pert_id", "cell_iname"]).pert_dose.nunique()
    cells = nd[nd >= 3].index
    rng = np.random.default_rng(SEED)
    pick = cells[rng.choice(len(cells), min(N_CELLS, len(cells)), replace=False)] if len(cells) > N_CELLS else cells
    sel = sig.set_index(["pert_id", "cell_iname"]).loc[pick].reset_index()
    print(f"{len(cells):,} dose cells with a mappable network; sampled {len(pick):,} -> {len(sel):,} signatures")

    with h5py.File(GCTX, "r") as f:
        mat = f["0/DATA/0/matrix"]
        col_ids, row_ids = dec(f["0/META/COL/id"][:]), dec(f["0/META/ROW/id"][:])
        gpos = pd.Series(np.arange(len(row_ids)), index=row_ids)
        genes = [g for g in sdw.index if g in gpos.index]
        gidx = gpos[genes].values
        keep = np.isin(col_ids, np.array(sorted(set(sel.sig_id)), dtype=str))
        X = np.empty((int(keep.sum()), len(gidx)), dtype=np.float32); k = 0
        for i in range(0, len(col_ids), BLOCK):
            s = np.where(keep[i:i + BLOCK])[0]
            if len(s):
                X[k:k + len(s)] = mat[i:i + BLOCK, :][s][:, gidx]; k += len(s)
        got = col_ids[keep]
    meta = sel.set_index("sig_id").loc[got].reset_index()
    sd = sdw[genes].values.astype(np.float64)
    # IDENTIFIER BRIDGE. The Level-5 matrix and bing_gene_quality use ENTREZ ids ("16", "23");
    # ARACNe regulon targets use ENSEMBL ("ENSG00000100479"). Zero overlap between them - the first
    # run of this script silently produced no rows because every regulon lookup missed.
    ens = q.set_index("gene_id").ensembl_id.astype(str)
    gi = {}
    for i, g in enumerate(genes):
        e = ens.get(g)
        if isinstance(e, str) and e.startswith("ENSG"):
            gi[e] = i
    print(f"identifier bridge: {len(gi):,} of {len(genes):,} measured genes carry an Ensembl id")
    print(f"read {X.shape[0]:,} signatures x {X.shape[1]:,} genes in {time.time()-t0:.0f}s")

    # ---- regulons, per network ---------------------------------------------------------------
    E = pd.read_csv(ROOT / "data/processed/lincs/trace_aracne/aracne_regulon_edges.tsv.gz", sep="\t")
    E = E[E.in_lincs_gene_space]
    REG = {}
    for net, d in E.groupby("network"):
        REG[net] = {r: (list(g.gene_id), np.sign(g.tfmode.values).astype(float))
                    for r, g in d.groupby("regulator") if len(g) >= MIN_GENES}

    rows = []
    for (pid, line), m in meta.groupby(["pert_id", "cell_iname"]):
        m = m.sort_values("pert_dose")
        doses = m.pert_dose.values.astype(float)
        if len(doses) < 3:
            continue
        Z = X[[meta.index.get_loc(i) for i in m.index]]          # doses x genes
        net = m.net.iloc[0]
        regs = REG.get(net)
        if not regs:
            continue
        top = doses[-1]
        for kk in KS:
            thr = kk * sd
            A = np.abs(Z).astype(np.float64)
            crossed_any = (A >= thr).any(axis=0)
            # BMD by linear interpolation between the first bracketing pair
            bmd = np.full(A.shape[1], np.nan)
            first = np.argmax(A >= thr, axis=0)
            for j in np.flatnonzero(crossed_any):
                i1 = first[j]
                if i1 == 0:
                    bmd[j] = doses[0]
                else:
                    y0, y1 = A[i1 - 1, j], A[i1, j]
                    f = 0.0 if y1 == y0 else (thr[j] - y0) / (y1 - y0)
                    bmd[j] = doses[i1 - 1] + f * (doses[i1] - doses[i1 - 1])
            direction = np.sign(Z[np.clip(first, 0, len(doses) - 1), np.arange(A.shape[1])])
            for reg, (gl, tfm) in regs.items():
                # keep gene and sign paired; truncating tfm would misalign signs to genes
                pairs = [(gi[g], sg) for g, sg in zip(gl, tfm) if g in gi]
                if len(pairs) < MIN_GENES:
                    continue
                idx = np.array([a for a, _ in pairs])
                sgn = np.array([b for _, b in pairs])
                cr, b, dr = crossed_any[idx], bmd[idx], direction[idx]
                conc = (dr == sgn) | (dr == 0)
                use = cr & conc
                disc = float((cr & ~conc).mean())
                if use.sum() < 3:
                    continue
                kmm = km_median(np.where(use, b, top), ~use, top)
                sub = b[use]
                rows.append(dict(pert_id=pid, cmap_name=m.cmap_name.iloc[0], cell_iname=line,
                                 network=net, regulator=reg, k=kk, n_genes=int(len(idx)),
                                 n_crossed=int(cr.sum()), censored_fraction=round(1 - float(cr.mean()), 4),
                                 discordant_fraction=round(disc, 4),
                                 bmd_km_median=None if kmm != kmm else round(kmm, 4),
                                 bmd_median_uncensored=round(float(np.median(sub)), 4),
                                 coherence_iqr=round(float(np.subtract(*np.percentile(sub, [75, 25]))), 4),
                                 coherence_cv=round(float(np.std(sub) / max(np.mean(sub), 1e-9)), 4),
                                 top_dose=round(float(top), 4)))
        if len(rows) and len(rows) % 20000 < 1200:
            print(f"   {len(rows):,} regulon rows, {time.time()-t0:.0f}s", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "o3_s1_regulon_bmd.tsv.gz", sep="\t", index=False)  # 1.6M rows: gzip (GitHub 100MB limit)
    p = df[df.k == K_PRIMARY]
    summ = dict(step="O-3 / D136-S1 regulon BMD with coherence and censoring", seed=SEED,
                k_profile=KS, k_primary=K_PRIMARY, k_rule="S0 crossing rate nearest 50%",
                n_cells=int(meta.groupby(['pert_id','cell_iname']).ngroups), n_rows=int(len(df)),
                n_regulon_rows_at_k_primary=int(len(p)),
                censored_fraction=dict(median=round(float(p.censored_fraction.median()), 4),
                                       p10=round(float(p.censored_fraction.quantile(.1)), 4),
                                       p90=round(float(p.censored_fraction.quantile(.9)), 4)),
                discordant_fraction_median=round(float(p.discordant_fraction.median()), 4),
                coherence_cv=dict(median=round(float(p.coherence_cv.median()), 4),
                                  p10=round(float(p.coherence_cv.quantile(.1)), 4),
                                  p90=round(float(p.coherence_cv.quantile(.9)), 4)),
                km_median_resolvable=int(p.bmd_km_median.notna().sum()),
                level2_implemented=False, runtime_s=round(time.time() - t0, 1))
    (OUT / "o3_s1_run.json").write_text(json.dumps(summ, indent=1, default=str))
    print(f"\n{len(df):,} regulon rows across {len(KS)} k values; {len(p):,} at k={K_PRIMARY}")
    print(f"  censored fraction  median {summ['censored_fraction']['median']:.3f} "
          f"(p10 {summ['censored_fraction']['p10']:.3f}, p90 {summ['censored_fraction']['p90']:.3f})")
    print(f"  discordant fraction median {summ['discordant_fraction_median']:.3f}")
    print(f"  coherence CV       median {summ['coherence_cv']['median']:.3f} "
          f"(p10 {summ['coherence_cv']['p10']:.3f}, p90 {summ['coherence_cv']['p90']:.3f})")
    print(f"  KM median resolvable in {summ['km_median_resolvable']:,} of {len(p):,} rows")
    print("wrote", OUT / "o3_s1_regulon_bmd.tsv")


if __name__ == "__main__":
    sys.exit(main())
