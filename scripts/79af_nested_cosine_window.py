#!/usr/bin/env python
"""D150's nested regulon cosine, re-run WINDOW-SEPARATED with the network mapping fixed (D191).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
SPEC: docs/SPEC_nested_cosine_window.md  (pre-registered; read it before changing anything here)

The cosine is UNCHANGED from `systems_validation/integrate/nested_cosine.py`. For each regulator r
of the cohort's own ARACNe network:

    v_d[r] = tumour-vs-normal z over r's regulon genes
    v_x[r] = drug response z over THE SAME genes
    cos[r] = cosine(v_d[r], v_x[r])        negative = the drug pushes those genes the other way

Two things differ, and only two:

(a) THE NETWORK MAPPING. nested_cosine.py:37 hardcodes BRCA-Her2 -> BRCA_BAS. There is no
    Her2-specific ARACNe network, and production resolves an unlisted group through its TCGA
    cohort (TCGA-BRCA -> BRCA), so Her2 was scored against BASAL biology. Rather than correct the
    literal, this imports `net_for()` from window_reversal.py - the function that already encodes
    the production rule (and that D171 had to fix for exactly this class of bug) - so the mapping
    has ONE source of truth.

(b) THE WINDOW. D150 used drug_line_profiles_bing_specific.f32, which 79c averages over doses AND
    times. This uses the condition-level profiles and splits early (<= 6 h) from late (>= 24 h);
    `other` and the unrecorded sentinel are DROPPED, never pooled (D180). Aggregation mirrors
    window_reversal.py - conditions, in-key lines only, line-weighted - because D182 established
    that following 79e instead silently disagrees with the window table.

Registered caveat (SPEC s2): the contribution-vector reparameterisation is NOT a live hypothesis -
it tracks the raw-z cosine at Spearman 0.942, 98.5% sign agreement. Only the time axis is.
"""
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")          # D165: the budget is < half of 64 threads

import numpy as np
import pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
PROC = ROOT / "data/processed/lincs"
OUT = ROOT / ("results/systems_validation/integrated/nested_cosine_window"
               + ("_all" if os.environ.get("ALL") == "1" else ""))

MIN_GENES = 15                                           # unchanged from D150
TOP_REGS = int(os.environ.get("TOP_REGS", "25"))          # D193 granular table cap; 0 = all
# The SPEC registers these four, for a like-for-like re-run of D150. ALL=1 instead sweeps every
# group that has a tumour-vs-normal contrast and writes to a SEPARATE directory, so the registered
# four-cohort result is never overwritten by the descriptive sweep. The ALL run is DESCRIPTIVE:
# section 4's acceptance applies to the registered cohorts only and is not re-evaluated over 41.
SPEC_COHORTS = ["BRCA-LumA", "BRCA-Her2", "TCGA-SKCM", "TCGA-LUAD"]
ALL = os.environ.get("ALL") == "1"
ANTI_ESTROGEN = {"fulvestrant", "tamoxifen", "raloxifene", "clomifene", "bazedoxifene", "toremifene"}
BRAF = {"vemurafenib", "dabrafenib"}


def _load_net_for():
    """Import net_for() from window_reversal.py - one source of truth for the mapping."""
    p = ROOT / "scripts/systems_validation/integrate/window_reversal.py"
    spec = importlib.util.spec_from_file_location("window_reversal", p)
    m = importlib.util.module_from_spec(spec)
    sys.modules["window_reversal"] = m
    spec.loader.exec_module(m)
    return m.net_for


def window_of(t):
    """D180 windows. The unrecorded sentinel is -1.0 and must NEVER read as early."""
    if t is None or (isinstance(t, float) and np.isnan(t)) or t <= 0:
        return None
    if t <= 6:
        return "early"
    if t >= 24:
        return "late"
    return None                                          # `other` is dropped, not pooled


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    net_for = _load_net_for()

    genes = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
    idx = pd.read_csv(PROC / "drug_condition_profiles_index.tsv", sep="\t")
    M = np.memmap(PROC / "drug_condition_profiles_specific.f32", dtype="float32", mode="r",
                  shape=(len(idx), len(genes)))
    gpos = {g: i for i, g in enumerate(genes.gene_id)}

    annot = pd.read_csv(ROOT / "results/lincs/drug_annotations.tsv", sep="\t")
    pid2name = dict(annot.drop_duplicates("pert_id")[["pert_id", "cmap_name"]].values)

    E = pd.read_csv(PROC / "trace_aracne/aracne_regulon_edges.tsv.gz", sep="\t")
    E = E[E.in_lincs_gene_space]
    avail = set(E.network.unique())
    mr = pd.read_csv(ROOT / "results/step6/aracne_project_v2_tf/group_mr.tsv", sep="\t")
    lw = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    ovr = pd.read_csv(ROOT / "config/group_networks.tsv", sep="\t").set_index("group_id").net.to_dict()
    cohort_of = (pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t")
                 .drop_duplicates("group_id").set_index("group_id").cohort.to_dict())

    idx["window"] = idx.pert_time.map(window_of)

    cohorts = SPEC_COHORTS
    if ALL:
        cohorts = sorted(c for c in cohort_of
                         if (ROOT / f"data/processed/step4_v2/{c}__gtex.tsv.gz").exists())
        print(f"ALL=1: descriptive sweep over {len(cohorts)} cohorts with a contrast "
              f"(the SPEC's acceptance covers only {SPEC_COHORTS})", flush=True)

    per_drug, summary, nets_used, skipped, granular = [], [], {}, [], []
    for cohort in cohorts:
        net = net_for(cohort, ovr, avail, cohort_of)
        nets_used[cohort] = net
        con = ROOT / f"data/processed/step4_v2/{cohort}__gtex.tsv.gz"
        if net is None or not con.exists():
            skipped.append(dict(cohort=cohort, reason="no network" if net is None else "no contrast"))
            print(f"{cohort}: SKIPPED ({'no network' if net is None else 'no contrast'})")
            continue

        dz = pd.read_csv(con, sep="\t")[["gene_id", "z"]].dropna().set_index("gene_id").z
        # D193: keep each edge's OWN enrichment factor w = tfmode * likelihood, aligned to the
        # gene list, so the cosine can be taken over the regulon's weighted contributions rather
        # than over bare z-scores. tfmode carries the SIGN of the edge (activating / repressing)
        # and likelihood its strength - exactly what "the regulon data point's own enrichment
        # factor" means, and what D150 threw away before the vectors ever met.
        regs = {r: (list(d.gene_id), (d.tfmode * d.likelihood).to_numpy(dtype=np.float64))
                for r, d in E[E.network == net].groupby("regulator")}
        nes = mr[mr.group_id == cohort].set_index("regulator").nes

        # in-key matched lines and their weights, exactly as window_reversal.py uses them
        L = lw[(lw.group_id == cohort) & lw.in_key][["ModelID", "weight"]].drop_duplicates("ModelID")
        wmap = dict(L.values)
        if not wmap:
            skipped.append(dict(cohort=cohort, reason="insufficient matched evidence: no in-key line"))
            print(f"{cohort}: SKIPPED (insufficient matched evidence - no in-key line)")
            continue

        sub = idx[idx.ModelID.isin(wmap) & idx.window.notna()]
        if sub.empty:
            skipped.append(dict(cohort=cohort, reason="insufficient matched evidence: no windowed condition"))
            print(f"{cohort}: SKIPPED (insufficient matched evidence - no windowed condition)")
            continue

        # line-weighted mean profile per (drug, window), over that drug's conditions in in-key lines
        keys, vecs = [], []
        for (pid, win), g in sub.groupby(["pert_id", "window"], sort=False):
            w = g.ModelID.map(wmap).to_numpy(dtype=np.float64)
            if w.sum() <= 0:
                continue
            P = np.asarray(M[g.index.to_numpy()], dtype=np.float64)
            keys.append((pid, win))
            vecs.append((P * w[:, None]).sum(0) / w.sum())
        if not keys:
            skipped.append(dict(cohort=cohort, reason="insufficient matched evidence: no profile"))
            continue
        X = np.vstack(vecs)

        rows = []
        for r, (gl0, w0) in regs.items():
            keep = [i for i, g in enumerate(gl0) if g in gpos and g in dz.index]
            if len(keep) < MIN_GENES:
                continue
            gl = [gl0[i] for i in keep]
            wj = w0[keep]                                   # the regulon's own enrichment factors
            vd = dz[gl].to_numpy(dtype=np.float64)
            S = X[:, np.array([gpos[g] for g in gl])]

            def _cos(vd_, S_):
                nd_ = np.linalg.norm(vd_)
                nx_ = np.linalg.norm(S_, axis=1)
                out = np.full(len(S_), np.nan)
                ok_ = (nx_ > 0) & (nd_ > 0)
                if nd_ > 0:
                    out[ok_] = (S_[ok_] @ vd_) / (nx_[ok_] * nd_)
                return out

            cos_z = _cos(vd, S)                             # D150, unchanged: bare z-scores
            cos_e = _cos(wj * vd, S * wj)                   # enrichment-weighted contributions
            if np.all(np.isnan(cos_z)) and np.all(np.isnan(cos_e)):
                continue
            rows.append(pd.DataFrame(dict(row=np.arange(len(S)), regulator=r,
                                          cos=cos_z, cos_enr=cos_e, n_genes=len(gl),
                                          nes=float(nes.get(r, np.nan)))))
        if not rows:
            skipped.append(dict(cohort=cohort, reason="no usable regulon"))
            continue

        df = pd.concat(rows, ignore_index=True).dropna(subset=["cos", "cos_enr"], how="all")
        df["w"] = df.nes.abs().fillna(0)
        kk = pd.DataFrame(keys, columns=["pert_id", "window"])
        df["pert_id"] = kk.pert_id.to_numpy()[df.row.to_numpy()]
        df["window"] = kk.window.to_numpy()[df.row.to_numpy()]

        # vectorised rather than groupby.apply: this pandas has no `include_groups`, and the
        # weighted mean is a ratio of two sums, so it needs no per-group Python at all.
        df["cw"] = df.cos * df.w
        df["cwe"] = df.cos_enr * df.w
        df["neg"] = (df.cos < 0).astype(float)
        df["neg_e"] = (df.cos_enr < 0).astype(float)
        agg = df.groupby(["pert_id", "window"]).agg(
            n_regulators=("cos", "size"), cos_mean=("cos", "mean"),
            cos_enr_mean=("cos_enr", "mean"),
            _cw=("cw", "sum"), _cwe=("cwe", "sum"), _w=("w", "sum"),
            frac_opposing=("neg", "mean"), frac_opposing_enr=("neg_e", "mean")).reset_index()
        agg["cos_weighted"] = agg._cw / agg._w.clip(lower=1e-9)
        agg["cos_enr_weighted"] = agg._cwe / agg._w.clip(lower=1e-9)
        agg = agg.drop(columns=["_cw", "_cwe", "_w"])

        # D193 granular view: one row per (cohort, window, drug, regulator). Capped to the
        # TOP_REGS regulators by |disease NES| so the file stays usable; TOP_REGS=0 keeps all.
        gr = df.copy()
        gr["cohort"], gr["network"] = cohort, net
        gr["cmap_name"] = gr.pert_id.map(pid2name)
        if TOP_REGS:
            keepr = (nes.abs().sort_values(ascending=False).head(TOP_REGS).index)
            gr = gr[gr.regulator.isin(set(keepr))]
        granular.append(gr[["cohort", "network", "window", "pert_id", "cmap_name", "regulator",
                            "n_genes", "nes", "cos", "cos_enr"]])
        agg["cohort"], agg["network"] = cohort, net
        agg["cmap_name"] = agg.pert_id.map(pid2name)
        # opposing percentile: within (cohort, window), rank of -cos_weighted
        agg["opposing_pct"] = agg.groupby("window").cos_weighted.rank(pct=True, ascending=False) * 100
        agg["opposing_pct_enr"] = agg.groupby("window").cos_enr_weighted.rank(pct=True, ascending=False) * 100
        per_drug.append(agg)

        for win, g in agg.groupby("window"):
            summary.append(dict(cohort=cohort, network=net, window=win, n_drugs=int(len(g)),
                                n_regulons_usable=int(df.regulator.nunique()),
                                cos_weighted_median=round(float(g.cos_weighted.median()), 5),
                                cos_enr_weighted_median=round(float(g.cos_enr_weighted.median()), 5),
                                frac_opposing_median=round(float(g.frac_opposing.median()), 4),
                                frac_opposing_enr_median=round(float(g.frac_opposing_enr.median()), 4)))
            print(f"{cohort:11s} {net:10s} {win:5s} n={len(g):5d} "
                  f"cos_w median {g.cos_weighted.median():+.4f}  [{time.time()-t0:.0f}s]", flush=True)

    if not per_drug:
        print("nothing scored", file=sys.stderr)
        return 1
    D = pd.concat(per_drug, ignore_index=True)
    D.to_csv(OUT / "per_drug.tsv.gz", sep="\t", index=False, float_format="%.6g")
    if granular:
        G = pd.concat(granular, ignore_index=True)
        G.to_csv(OUT / "per_regulator.tsv.gz", sep="\t", index=False, float_format="%.6g")
        print(f"granular: {len(G):,} (cohort, window, drug, regulator) rows -> per_regulator.tsv.gz")
    S = pd.DataFrame(summary)
    S.to_csv(OUT / "summary.tsv", sep="\t", index=False)

    # ---- controls, named in SPEC s3 before the run
    ctl = []
    for label, names, coh in (("anti_estrogen", ANTI_ESTROGEN, "BRCA-LumA"),
                              ("anti_estrogen", ANTI_ESTROGEN, "BRCA-Her2"),
                              ("braf", BRAF, "TCGA-SKCM")):
        for win in ("early", "late"):
            g = D[(D.cohort == coh) & (D.window == win) & D.cmap_name.isin(names)]
            ctl.append(dict(control=label, cohort=coh, window=win, n=int(len(g)),
                            median_opposing_pct=round(float(g.opposing_pct.median()), 2) if len(g) else None,
                            median_opposing_pct_enr=round(float(g.opposing_pct_enr.median()), 2) if len(g) else None,
                            median_cos_weighted=round(float(g.cos_weighted.median()), 5) if len(g) else None,
                            median_cos_enr_weighted=round(float(g.cos_enr_weighted.median()), 5) if len(g) else None,
                            drugs=",".join(sorted(g.cmap_name.dropna().unique()))))
    C = pd.DataFrame(ctl)
    C.to_csv(OUT / "controls.tsv", sep="\t", index=False)
    print("\n" + C.to_string(index=False))

    # ---- the single pre-registered acceptance criterion (SPEC s4.1)
    ae = {w: C[(C.control == "anti_estrogen") & (C.cohort == "BRCA-LumA") & (C.window == w)]
          for w in ("early", "late")}
    e = ae["early"].median_opposing_pct.iloc[0] if len(ae["early"]) else None
    l = ae["late"].median_opposing_pct.iloc[0] if len(ae["late"]) else None
    verdict = dict(criterion="anti-estrogens in BRCA-LumA: median opposing pct higher EARLY than LATE, and early > 50",
                   early=e, late=l,
                   passed=bool(e is not None and l is not None and e > l and e > 50))
    print(f"\nACCEPTANCE: early={e} late={l} -> {'PASS' if verdict['passed'] else 'FAIL'}")

    (OUT / "run.json").write_text(json.dumps(dict(
        script="79af_nested_cosine_window.py", spec="docs/SPEC_nested_cosine_window.md",
        min_genes=MIN_GENES, cohorts=cohorts, mode="all" if ALL else "spec", networks_resolved=nets_used,
        d150_mapping_bug="BRCA-Her2 was hardcoded to BRCA_BAS; net_for() now resolves it",
        skipped=skipped, acceptance=verdict, n_rows=int(len(D)),
        runtime_s=round(time.time() - t0, 1)), indent=1))
    print("wrote", OUT / "run.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
