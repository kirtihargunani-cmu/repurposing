#!/usr/bin/env python
"""O-3 / D136-S2: the INTRINSIC NULL for regulon dose-response coherence.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
SPEC: docs/SPEC_dose_response_bmd.md section 6 -
  "Random gene sets drawn from the same drug's own profiles, matched on set size AND on the
   sd_within distribution of the real regulon - otherwise the null is easier simply by containing
   quieter genes."

Why this is required. D148 measured coherence CV (coefficient of variation of member-gene benchmark
doses) at median 0.955, but a dispersion number means nothing without a baseline: scattered relative
to what? The null answers that with gene sets that share the real regulon's size and noise profile
but none of its biology.

Matching procedure: genes are binned into deciles of sd_within (replicate standard deviation). For a
real regulon, the per-decile counts are recorded, and every null set draws exactly that many genes
from each decile. Size and noise profile are therefore matched by construction, and the only thing
that differs is whether the genes share a regulator.

Reported per (compound, line, regulon): observed coherence, the null median, and an empirical
percentile - the share of null sets whose coherence is TIGHTER (lower CV) than the observed one. A
real regulatory event should sit in the LOW tail: tighter than most random sets of the same size and
noise.

Output: results/systems_validation/observability/{o3_s2_null.tsv.gz, o3_s2_run.json}
"""
import json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd, h5py

ROOT = Path("/home/kirti/projects/repurposing")
RAW = ROOT / "data/raw/lincs2020"
OUT = ROOT / "results/systems_validation/observability"
GCTX = RAW / "level5_beta_trt_cp_n720216x12328.gctx"
SEED, N_CELLS, BLOCK = 20261118, 400, 20000      # same seed/sample as S1, so rows correspond
K = 2.0                                           # the S0-selected value
N_NULL, MIN_GENES, N_DECILES = 40, 10, 10
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "6"


def dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def main():
    t0 = time.time()
    q = pd.read_csv(ROOT / "results/lincs/gene_quality/bing_gene_quality.tsv", sep="\t",
                    dtype={"gene_id": str})
    sdw = q.set_index("gene_id").sd_within
    models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", low_memory=False)
    lin = dict(zip(models.StrippedCellLineName.astype(str).str.upper(), models.OncotreeLineage))
    LIN2NET = {"Skin": "SKCM", "Breast": "BRCA", "Lung": "LUAD", "Bowel": "COADREAD",
               "Prostate": "PRAD", "Pancreas": "PAAD", "Liver": "LIHC", "Ovary/Fallopian Tube": "OV",
               "Kidney": "KIRC", "Esophagus/Stomach": "STAD", "Uterus": "UCEC",
               "Bladder/Urinary Tract": "BLCA", "CNS/Brain": "GBM", "Head and Neck": "HNSC",
               "Thyroid": "THCA", "Soft Tissue": "SARC", "Cervix": "CESC"}

    sig = pd.read_csv(RAW / "siginfo_beta.txt", sep="\t", low_memory=False,
                      usecols=["sig_id", "pert_id", "cmap_name", "pert_type", "cell_iname",
                               "pert_dose", "pert_time", "is_hiq", "tas"])
    sig = sig[(sig.pert_type == "trt_cp") & (sig.is_hiq == 1) & (sig.tas >= 0.4) & (sig.pert_dose > 0)]
    sig = sig.sort_values("tas", ascending=False).drop_duplicates(["pert_id", "cell_iname", "pert_dose"])
    sig["net"] = sig.cell_iname.astype(str).str.upper().map(lin).map(LIN2NET)
    sig = sig[sig.net.notna()]
    nd = sig.groupby(["pert_id", "cell_iname"]).pert_dose.nunique()
    cells = nd[nd >= 3].index
    rng = np.random.default_rng(SEED)
    pick = cells[rng.choice(len(cells), min(N_CELLS, len(cells)), replace=False)] if len(cells) > N_CELLS else cells
    sel = sig.set_index(["pert_id", "cell_iname"]).loc[pick].reset_index()

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
    ens = q.set_index("gene_id").ensembl_id.astype(str)
    gi = {}
    for i, g in enumerate(genes):
        e = ens.get(g)
        if isinstance(e, str) and e.startswith("ENSG"):
            gi[e] = i
    # sd_within deciles: the matching strata
    dq = np.quantile(sd, np.linspace(0, 1, N_DECILES + 1)[1:-1])
    strata = np.digitize(sd, dq)
    by_stratum = [np.flatnonzero(strata == s) for s in range(N_DECILES)]
    print(f"read {X.shape[0]:,} signatures; {len(gi):,} genes with Ensembl ids; "
          f"sd_within deciles sized {[len(b) for b in by_stratum]}")

    E = pd.read_csv(ROOT / "data/processed/lincs/trace_aracne/aracne_regulon_edges.tsv.gz", sep="\t")
    E = E[E.in_lincs_gene_space]
    REG = {}
    for net, d in E.groupby("network"):
        REG[net] = {r: (list(g.gene_id), np.sign(g.tfmode.values).astype(float))
                    for r, g in d.groupby("regulator") if len(g) >= MIN_GENES}

    def coherence(bmd_sub):
        return float(np.std(bmd_sub) / max(np.mean(bmd_sub), 1e-9))

    rows = []
    nrng = np.random.default_rng(SEED + 1)
    for (pid, line), m in meta.groupby(["pert_id", "cell_iname"]):
        m = m.sort_values("pert_dose")
        doses = m.pert_dose.values.astype(float)
        if len(doses) < 3:
            continue
        regs = REG.get(m.net.iloc[0])
        if not regs:
            continue
        Z = X[m.index.values]
        A = np.abs(Z).astype(np.float64)
        thr = K * sd
        crossed = (A >= thr).any(axis=0)
        first = np.argmax(A >= thr, axis=0)
        bmd = np.full(A.shape[1], np.nan)
        for j in np.flatnonzero(crossed):
            i1 = first[j]
            if i1 == 0:
                bmd[j] = doses[0]
            else:
                y0, y1 = A[i1 - 1, j], A[i1, j]
                f = 0.0 if y1 == y0 else (thr[j] - y0) / (y1 - y0)
                bmd[j] = doses[i1 - 1] + f * (doses[i1] - doses[i1 - 1])
        direction = np.sign(Z[np.clip(first, 0, len(doses) - 1), np.arange(A.shape[1])])

        for reg, (gl, tfm) in regs.items():
            pairs = [(gi[g], sg) for g, sg in zip(gl, tfm) if g in gi]
            if len(pairs) < MIN_GENES:
                continue
            idx = np.array([a for a, _ in pairs]); sgn = np.array([b for _, b in pairs])
            conc = (direction[idx] == sgn) | (direction[idx] == 0)
            use = crossed[idx] & conc
            if use.sum() < 5:
                continue
            obs = coherence(bmd[idx][use])
            obs_cens = 1 - float(crossed[idx].mean())
            # null: same size, same sd_within decile composition, no shared regulator
            counts = np.bincount(strata[idx], minlength=N_DECILES)
            nulls, ncens = np.empty(N_NULL), np.empty(N_NULL)
            for b in range(N_NULL):
                draw = np.concatenate([nrng.choice(by_stratum[s], counts[s], replace=False)
                                       for s in range(N_DECILES) if counts[s]])
                cu = crossed[draw]
                ncens[b] = 1 - float(cu.mean())
                nulls[b] = coherence(bmd[draw][cu]) if cu.sum() >= 5 else np.nan
            v = nulls[~np.isnan(nulls)]
            if len(v) < 10:
                continue
            rows.append(dict(pert_id=pid, cmap_name=m.cmap_name.iloc[0], cell_iname=line,
                             regulator=reg, n_genes=int(len(idx)), n_used=int(use.sum()),
                             obs_coherence_cv=round(obs, 4),
                             null_coherence_cv_median=round(float(np.median(v)), 4),
                             pct_null_tighter=round(float((v < obs).mean()), 4),
                             obs_censored=round(obs_cens, 4),
                             null_censored_median=round(float(np.median(ncens)), 4),
                             n_null_valid=int(len(v))))
        if len(rows) % 20000 < 1000 and rows:
            print(f"   {len(rows):,} rows, {time.time()-t0:.0f}s", flush=True)

    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "o3_s2_null.tsv.gz", sep="\t", index=False)
    tighter = float((df.pct_null_tighter < 0.05).mean())
    summ = dict(step="O-3 / D136-S2 intrinsic null for regulon coherence", spec="section 6",
                seed=SEED, k=K, n_null=N_NULL, matched_on=["set size", "sd_within decile counts"],
                n_rows=int(len(df)),
                obs_coherence_median=round(float(df.obs_coherence_cv.median()), 4),
                null_coherence_median=round(float(df.null_coherence_cv_median.median()), 4),
                median_pct_null_tighter=round(float(df.pct_null_tighter.median()), 4),
                frac_rows_better_than_95pct_of_null=round(tighter, 4),
                obs_censored_median=round(float(df.obs_censored.median()), 4),
                null_censored_median=round(float(df.null_censored_median.median()), 4),
                runtime_s=round(time.time() - t0, 1))
    (OUT / "o3_s2_run.json").write_text(json.dumps(summ, indent=1))
    print(f"\n{len(df):,} regulon rows with a valid matched null (N={N_NULL} draws each)")
    print(f"  observed coherence CV   median {summ['obs_coherence_median']}")
    print(f"  NULL coherence CV       median {summ['null_coherence_median']}")
    print(f"  share of null sets tighter than observed: median {summ['median_pct_null_tighter']}")
    print(f"  rows beating 95% of their own null:       {100*tighter:.2f}%")
    print(f"  censored  observed {summ['obs_censored_median']} vs null {summ['null_censored_median']}")
    print("wrote", OUT / "o3_s2_null.tsv.gz")


if __name__ == "__main__":
    sys.exit(main())
