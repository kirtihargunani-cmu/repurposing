#!/usr/bin/env python
"""D128 A0: label-blind feasibility audit (docs/IMPLEMENTATION_target_mechanism_scoring.md Section 4 A0; D128-A0).

Builds one row per (cohort, pert_id, source line, target unit) of the production lead population and labels it
eligible / excluded / unscorable with a reason, separately for M and R. Counts only: no M or R value is computed, no
outcome label (PRISM / GDSC values, D127 outcomes, known use) is read. Section 0 values come from config.py and are
NOT frozen here.

Outputs results/lincs/d128/a0/:
  coverage.tsv                 per target unit: anchored TFs by hop, footprint genes (landmark / inferred), sign
                               resolution, hub fraction, floors at x0.5 / x1 / x2
  eligibility.tsv.gz           per (cohort, pert_id, line, unit): M_status, R_status, status, reasons (large; not in git)
  eligibility_summary.tsv      counts by level x status x reason
  cohort_feasibility.tsv       per cohort: population, direct / borrowed, contrast, eligible pairs
  source_line_context.tsv      per LINCS line: lineage, profiles, genotype panel, ERBB2 / RB1 copy number, RB1 status
  control_availability.tsv     plan control compounds: pert_ids, lines, in-context lines, unit coverage
  wave_feasibility.tsv         Section 5 waves: expected status from this audit
  patient_strata.tsv           TCGA genotype strata counts (MC3, GISTIC)
  profile_value_audit.tsv      |z| distribution of the 79c matrices (P9) and orientation checks
  floor_admission.tsv          rows admitted at each P10 floor multiplier
  pert_multiplicity.tsv        compounds with several pert_ids (P12)
  unmatched_ids.tsv            every ID that failed to join, by kind
  external_screen_prior_use.tsv  scripts / decisions that read PRISM or GDSC (grep only)
  run_manifest.json            git SHA, input checksums, config, counts, timings
"""
import argparse
import json
import os
import pickle
import re
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import config as C  # noqa: E402
import d128io as IO  # noqa: E402

A0 = C.OUT / "a0"
T0 = time.time()
LOG = []


def log(msg):
    line = f"[{time.time() - T0:7.1f}s] {msg}"
    LOG.append(line); print(line, flush=True)


# ------------------------------------------------------------------------------------------------------------------
# footprint: target unit -> anchored TFs (signed, by hop) -> CollecTRI genes
# ------------------------------------------------------------------------------------------------------------------
# Signs are bitmasks, not sets: 1 = {+1}, 2 = {-1}, 3 = both (ambiguous). The first version used Python sets inside
# the BFS and a pandas filter per reached TF, which made the loop single-threaded and compute-bound (see D128-A0
# result: 43 min without finishing). Masks, one cached BFS per gene and pre-filtered regulon gene lists replace it.
PLUS, MINUS, BOTH = 1, 2, 3


def mul(m1, m2):
    """Product of two sign masks."""
    if m1 == BOTH or m2 == BOTH:
        return BOTH
    return PLUS if m1 == m2 else MINUS


def build_graph(edges):
    g = {}
    for s, d, sign, disp in edges[["src", "dst", "sign", "sign_disputed"]].itertuples(index=False):
        g.setdefault(s, []).append((d, BOTH if bool(disp) else (PLUS if int(sign) == 1 else MINUS)))
    return g


def reach_regulators(graph, start, max_hops, regs):
    """BFS from start, keeping only regulators: {reg: (shortest hop, sign mask over all shortest paths)}.

    The mask is over *shortest* paths, so a walk at the largest hop arm also answers every smaller arm (filter on hop).
    """
    dist, mask = {start: 0}, {start: PLUS}
    out = {}
    if start in regs:
        out[start] = (0, PLUS)
    frontier = [start]
    for h in range(1, max_hops + 1):
        nxt = {}
        for u in frontier:
            mu = mask[u]
            for v, me in graph.get(u, ()):
                if v in dist and dist[v] < h:
                    continue
                s = mul(mu, me)
                nxt[v] = nxt.get(v, 0) | s
        for v, s in nxt.items():
            if v not in dist:
                dist[v], mask[v] = h, s
                if v in regs:
                    out[v] = (h, s)
            elif dist[v] == h:                      # reached again at the same (shortest) hop
                mask[v] |= s
                if v in regs:
                    out[v] = (h, mask[v])
        frontier = [v for v in nxt if dist.get(v) == h]
    return out


class Reach:
    """Per-gene reachable-regulator cache, persisted per (network, source checksum, max hops) as the plan asks."""

    def __init__(self, name, graph, regs, max_hops, checksum, cache_dir):
        self.name, self.graph, self.regs, self.max_hops = name, graph, regs, max_hops
        self.path = Path(cache_dir) / f"reach_{name}_{checksum[:12]}_h{max_hops}.pkl"
        self.mem, self.hits, self.misses = {}, 0, 0
        if self.path.exists():
            self.mem = pickle.loads(self.path.read_bytes())

    def __call__(self, gene):
        r = self.mem.get(gene)
        if r is None:
            self.misses += 1
            r = self.mem[gene] = reach_regulators(self.graph, gene, self.max_hops, self.regs)
        else:
            self.hits += 1
        return r

    def save(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_bytes(pickle.dumps(self.mem, protocol=4))
        return self.path


def unit_footprint(genes, reach, reg_genes, max_hops, _cache={}):
    """Counts for one target unit. Expected gene sign up to the drug's action = path sign x tfmode.

    `reg_genes[tf] = (gene ids, tfmodes)`, already restricted to the LINCS gene space, so no pandas in this loop.
    Memoised on (network, unit gene set, hop arm): units share target genes across drugs.
    """
    key = (reach.name, frozenset(genes), max_hops)
    if key in _cache:
        return _cache[key]
    tf_hop, tf_mask = {}, {}
    for g in genes:
        for t, (h, m) in reach(g).items():
            if h > max_hops:
                continue
            if t not in tf_hop or h < tf_hop[t]:
                tf_hop[t], tf_mask[t] = h, m
            elif h == tf_hop[t]:
                tf_mask[t] |= m
    gene_exp, gene_any, per_tf = {}, set(), {}
    for t, h in tf_hop.items():
        gl = reg_genes.get(t)
        if gl is None:
            continue
        ids, modes = gl
        gene_any.update(ids)
        if tf_mask[t] != BOTH:
            p = 1 if tf_mask[t] == PLUS else -1
            for gid, mode in zip(ids, modes):
                e = p * mode
                prev = gene_exp.get(gid)
                gene_exp[gid] = e if prev is None else (prev if prev == e else 0)   # 0 = conflicting expectations
            per_tf[t] = set(ids)
    resolved = {g for g, s in gene_exp.items() if s}
    conflicted = {g for g, s in gene_exp.items() if not s}
    hub = max((len(v & resolved) for v in per_tf.values()), default=0)
    hop0 = [t for t, h in tf_hop.items() if h == 0 and t in per_tf]
    hop0_resolved = set().union(*[per_tf[t] for t in hop0]) & resolved if hop0 else set()
    out = dict(
        tf_hop=tf_hop, n_tfs={h: sum(1 for x in tf_hop.values() if x <= h) for h in C.HOP_ARMS},
        n_anchored_tfs=len(tf_hop), n_resolved_tfs=sum(1 for m in tf_mask.values() if m != BOTH),
        n_ambiguous_tfs=sum(1 for m in tf_mask.values() if m == BOTH),
        genes_any=gene_any, genes_resolved=resolved, genes_conflicted=conflicted, hop0_resolved=hop0_resolved,
        hub_fraction=(hub / len(resolved)) if resolved else np.nan,
        hub_tf=max(per_tf, key=lambda t: len(per_tf[t] & resolved)) if per_tf and resolved else "")
    _cache[key] = out
    return out


def floors(kind, n_res, n_res_tfs, n_hop0, mult):
    fg, ft = C.FLOOR_GENES * mult, C.FLOOR_TFS * mult
    if kind == "tf":
        return n_hop0 >= fg
    return n_res >= fg and n_res_tfs >= ft


# ------------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=C.THREADS)
    a = ap.parse_args()
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = str(a.threads)
    if A0.resolve().is_relative_to((C.ROOT / "results/lincs").resolve()) and A0.parent.name != "d128":
        sys.exit("output root must be results/lincs/d128 (no writes to production names)")
    A0.mkdir(parents=True, exist_ok=True)
    unmatched = []

    # --- inputs ----------------------------------------------------------------------------------------------------
    idx = IO.profile_index(); genes = IO.profile_genes()
    log(f"profile index {len(idx)} rows ({idx.pert_id.nunique()} pert_ids, {idx.ModelID.nunique()} lines); genes {len(genes)}")
    gene_space = set(genes.gene_id)
    land = set(genes.gene_id[genes.feature_space == "landmark"])
    sym2id = dict(zip(genes.symbol, genes.gene_id))
    ann = pd.read_csv(C.DRUG_ANNOT, sep="\t", usecols=["pert_id", "cmap_name"]).drop_duplicates("pert_id")
    idx = idx.merge(ann, on="pert_id", how="left")
    for p in idx.pert_id[idx.cmap_name.isna()].unique():
        unmatched.append(dict(kind="pert_id_without_cmap_name", id=p, detail="profile index"))

    # --- A0 action 5 / P9: profile values and orientation ----------------------------------------------------------
    spec = IO.profiles(C.PROFILES_SPECIFIC, len(idx), len(genes))
    raw = IO.profiles(C.PROFILES_RAW, len(idx), len(genes))
    rows = []
    lm_cols = np.flatnonzero(genes.feature_space.values == "landmark")
    inf_cols = np.flatnonzero(genes.feature_space.values != "landmark")
    for arm, M in (("specific", spec), ("raw", raw)):
        X = np.asarray(M)
        for space, cols in (("all", None), ("landmark", lm_cols), ("inferred", inf_cols)):
            v = np.abs(X if cols is None else X[:, cols]).ravel()
            fin = np.isfinite(v)
            r = dict(arm=arm, feature_space=space, n_values=int(v.size), n_nonfinite=int((~fin).sum()))
            v = v[fin]
            for q, x in zip(C.Z_QUANTILES, np.quantile(v, C.Z_QUANTILES)):
                r[f"abs_z_q{q}"] = round(float(x), 4)
            for t in C.Z_THRESHOLDS:
                r[f"share_abs_z_gt{t:g}"] = round(float((v > t).mean()), 5)
            f = np.tanh(v / C.F_DIVISOR)
            r["share_f_gt0.9"] = round(float((f > 0.9).mean()), 5)
            r["median_f"] = round(float(np.median(f)), 4)
            rows.append(r)
        sd = X.std(axis=1)
        rows.append(dict(arm=arm, feature_space="per_profile_sd", n_values=len(sd), n_nonfinite=int((~np.isfinite(sd)).sum()),
                         **{f"abs_z_q{q}": round(float(x), 4) for q, x in zip(C.Z_QUANTILES, np.quantile(sd, C.Z_QUANTILES))},
                         n_zero_sd=int((sd == 0).sum())))
        del X
    rng = np.random.default_rng(C.SEED)
    pick = rng.choice(len(idx), 500, replace=False)
    cors = [np.corrcoef(np.asarray(spec[i]), np.asarray(raw[i]))[0, 1] for i in pick]
    cors_shift = [np.corrcoef(np.asarray(spec[i]), np.asarray(raw[(i + 1) % len(idx)]))[0, 1] for i in pick]
    rows.append(dict(arm="orientation", feature_space="raw_vs_specific_same_row_r",
                     **{f"abs_z_q{q}": round(float(x), 4) for q, x in zip(C.Z_QUANTILES, np.quantile(cors, C.Z_QUANTILES))},
                     n_values=len(cors)))
    rows.append(dict(arm="orientation", feature_space="raw_vs_specific_shifted_row_r",
                     **{f"abs_z_q{q}": round(float(x), 4) for q, x in zip(C.Z_QUANTILES, np.quantile(cors_shift, C.Z_QUANTILES))},
                     n_values=len(cors_shift)))
    pva = pd.DataFrame(rows)
    pva.to_csv(A0 / "profile_value_audit.tsv", sep="\t", index=False)
    med_same, med_shift = float(np.median(cors)), float(np.median(cors_shift))
    if not med_same > med_shift + 0.2:
        sys.exit(f"orientation check failed: raw-vs-specific same-row r {med_same:.3f} vs shifted {med_shift:.3f}")
    prof_sd = np.asarray(spec).std(axis=1)
    log(f"profile audit done; same-row r median {med_same:.3f} vs shifted {med_shift:.3f}")

    # --- regulons and signaling prior ------------------------------------------------------------------------------
    col = IO.collectri(C.COLLECTRI)
    regs = set(col.regulator)
    n_col_targets_off = int((~col.target.isin(gene_space)).sum())
    unmatched.append(dict(kind="collectri_edges_target_not_in_lincs_genes", id=str(n_col_targets_off),
                          detail=f"of {len(col)} edges ({col.target.nunique()} targets)"))
    ck = col[col.target.isin(gene_space)]            # restrict once: the inner loop must not touch pandas
    reg_genes = {r: (tuple(d.target), tuple(d.tfmode.astype(int))) for r, d in ck.groupby("regulator")}
    om = IO.signaling(C.OMNIPATH); sg = IO.signaling(C.SIGNOR)
    MAXH = max(C.HOP_ARMS)
    REACH = {"omnipath": Reach("omnipath", build_graph(om), regs, MAXH, IO.sha256(C.OMNIPATH), A0 / "cache"),
             "signor": Reach("signor", build_graph(sg), regs, MAXH, IO.sha256(C.SIGNOR), A0 / "cache")}
    G_om = REACH["omnipath"].graph
    log(f"CollecTRI {len(regs)} regulators, {len(col)} edges ({len(ck)} in the LINCS gene space, {len(reg_genes)} regulators); "
        f"OmniPath {len(om)} edges; SIGNOR {len(sg)} edges; reach cache "
        + ", ".join(f"{k}: {len(v.mem)} genes preloaded" for k, v in REACH.items()))

    # --- target units ----------------------------------------------------------------------------------------------
    U = pd.read_csv(C.UNITS, sep="\t", dtype={"action_flag": str, "tid": str}).fillna({"action_flag": ""})
    U["genes_str"] = U.genes
    U["genes"] = U.genes.str.split("|")
    # 79t2 leaves tid empty for units that come from an annotation rather than a ChEMBL target record, so the id
    # falls back to the unit's genes; a residual collision gets a suffix and is recorded, never silently merged.
    U["unit_id"] = U.cmap_name + "::" + U.tid.fillna("annot:" + U.genes_str)
    n_annot = int(U.tid.isna().sum())
    unmatched.append(dict(kind="unit_without_chembl_tid_id_from_genes", id=str(n_annot), detail=f"of {len(U)} D121 units"))
    dup = U.unit_id.duplicated(keep=False)
    if dup.any():
        for u in sorted(U.unit_id[dup].unique()):
            unmatched.append(dict(kind="duplicate_unit_id_suffixed", id=u, detail=f"{int((U.unit_id == u).sum())} rows"))
        U.loc[dup, "unit_id"] = U.unit_id[dup] + "#" + U[dup].groupby("unit_id").cumcount().add(1).astype(str)
    assert U.unit_id.is_unique, "unit_id still not unique"
    cov_rows, fp = [], {}
    t_fp = time.time()
    for n_done, r in enumerate(U.itertuples(index=False), 1):
        out = {arm: unit_footprint(r.genes, REACH[arm], reg_genes, MAXH) for arm in ("omnipath", "signor")}
        prim = unit_footprint(r.genes, REACH["omnipath"], reg_genes, C.MAX_HOPS)
        if n_done % 25 == 0 or n_done == len(U):
            el = time.time() - t_fp
            log(f"footprints {n_done}/{len(U)} units, {el:.0f}s elapsed, ETA {el / n_done * (len(U) - n_done):.0f}s "
                f"(reach cache hits {sum(v.hits for v in REACH.values())}, misses {sum(v.misses for v in REACH.values())})")
        fp[r.unit_id] = prim
        kind = "tf" if any(g in regs for g in r.genes) else "signaling"
        res_l = prim["genes_resolved"] & land
        d = dict(unit_id=r.unit_id, cmap_name=r.cmap_name, tid=r.tid, ttype=r.ttype, genes="|".join(r.genes), k=r.k,
                 orient=r.orient, action_flag=r.action_flag, actions=r.actions, kind=kind,
                 genes_in_omnipath=sum(1 for g in r.genes if g in G_om),
                 n_anchored_tfs=prim["n_anchored_tfs"], n_resolved_tfs=prim["n_resolved_tfs"],
                 n_ambiguous_tfs=prim["n_ambiguous_tfs"], ambiguous_fraction=round(prim["n_ambiguous_tfs"] / prim["n_anchored_tfs"], 3) if prim["n_anchored_tfs"] else np.nan,
                 n_genes_total=len(prim["genes_any"]), n_genes_total_landmark=len(prim["genes_any"] & land),
                 n_sign_resolved=len(prim["genes_resolved"]), n_resolved_landmark=len(res_l),
                 n_resolved_inferred=len(prim["genes_resolved"]) - len(res_l), n_conflicted=len(prim["genes_conflicted"]),
                 n_hop0_resolved=len(prim["hop0_resolved"]), n_hop0_resolved_landmark=len(prim["hop0_resolved"] & land),
                 hub_fraction=round(prim["hub_fraction"], 3) if prim["hub_fraction"] == prim["hub_fraction"] else np.nan,
                 hub_tf=prim["hub_tf"],
                 **{f"tfs_le{h}hop_omnipath": out["omnipath"]["n_tfs"][h] for h in C.HOP_ARMS},
                 **{f"tfs_le{h}hop_signor": out["signor"]["n_tfs"][h] for h in C.HOP_ARMS},
                 n_sign_resolved_signor=len(out["signor"]["genes_resolved"]))
        for m in C.FLOOR_MULTIPLIERS:
            d[f"floor_x{m:g}"] = floors(kind, d["n_sign_resolved"], d["n_resolved_tfs"], d["n_hop0_resolved"], m)
            d[f"floor_x{m:g}_landmark"] = floors(kind, d["n_resolved_landmark"],
                                                d["n_resolved_tfs"], d["n_hop0_resolved_landmark"], m)
        d["M_unit_status"] = ("ambiguous_action" if r.action_flag in C.M_AMBIGUOUS_FLAGS else
                              "no_anchored_tf" if not prim["n_anchored_tfs"] else
                              "below_floor" if not d["floor_x1"] else "eligible")
        cov_rows.append(d)
        for g in r.genes:
            if g not in sym2id and g not in regs and g not in G_om:
                unmatched.append(dict(kind="unit_gene_not_in_lincs_collectri_omnipath", id=g, detail=r.unit_id))
            elif g not in G_om and g not in regs:
                unmatched.append(dict(kind="unit_gene_not_in_omnipath", id=g, detail=r.unit_id))
    cov = pd.DataFrame(cov_rows)
    for v in REACH.values():
        log(f"reach cache {v.name}: {len(v.mem)} genes, {v.hits} hits / {v.misses} misses -> {v.save().name}")
    cov.to_csv(A0 / "coverage.tsv", sep="\t", index=False)
    log(f"coverage: {len(cov)} units; M-eligible at x1 {int((cov.M_unit_status == 'eligible').sum())}")

    # --- population and source lines --------------------------------------------------------------------------------
    L = pd.read_csv(C.LAYERED, sep="\t", usecols=["unit", "pert_id", "cmap_name", "level", "layered_S", "no_protein_target"])
    L = L.rename(columns={"unit": "cohort"})
    sg_ = pd.read_csv(C.SPLIT_GROUPS, sep="\t", usecols=["group_id", "pert_id", "source", "passes_gate", "line_match", "transfer_conf"],
                      low_memory=False).rename(columns={"group_id": "cohort"})
    ss_ = pd.read_csv(C.SPLIT_STRATA, sep="\t", usecols=["stratum", "parent_group", "pert_id", "source", "passes_gate", "line_match",
                                                          "transfer_conf"], low_memory=False).rename(columns={"stratum": "cohort"})
    split = pd.concat([sg_.assign(parent_group=sg_.cohort), ss_])
    P = L.merge(split, on=["cohort", "pert_id"], how="left", validate="one_to_one")
    for r in P[P.source.isna()].itertuples(index=False):
        unmatched.append(dict(kind="population_row_not_in_split_table", id=f"{r.cohort}|{r.pert_id}", detail=r.cmap_name))
    wg = pd.read_csv(C.LINE_W_GROUPS, sep="\t", usecols=["group_id", "ModelID", "weight", "in_key"]).rename(columns={"group_id": "cohort"})
    ws = pd.read_csv(C.LINE_W_STRATA, sep="\t", usecols=["stratum", "ModelID", "weight", "in_key"]).rename(columns={"stratum": "cohort"})
    W = pd.concat([wg, ws]); W = W[W.in_key]
    key_lines = {c: dict(zip(d.ModelID, d.weight)) for c, d in W.groupby("cohort")}
    drug_lines = {p: list(d.ModelID) for p, d in idx.groupby("pert_id")}
    nsig = dict(zip(zip(idx.pert_id, idx.ModelID), idx.n_sig))
    prow = dict(zip(zip(idx.pert_id, idx.ModelID), range(len(idx))))
    units_by_name = {n: list(d.unit_id) for n, d in U.groupby("cmap_name")}

    # --- disease gene sets for R (P1, P2, P8) ------------------------------------------------------------------------
    de_sig, de_info = {}, {}
    for level, (proc, summ) in (("group_id", C.DE_GROUPS), ("stratum", C.DE_STRATA)):
        S = pd.read_csv(summ, sep="\t")
        for c in P.cohort[P.level == level].unique():
            ok = S[(S.group_id == c) & (S.status == "ok")].contrast.tolist() if "group_id" in S.columns else []
            con = next((x for x in C.CONTRAST_ORDER if x in ok), None)
            d = IO.de_table(proc, c, con) if con else None
            if d is None:
                de_info[c] = dict(contrast=con or "none", n_de_genes=0)
                unmatched.append(dict(kind="cohort_without_de_contrast", id=c, detail=f"{summ.name}: ok={ok}"))
                continue
            sig = set(d.gene_id[d[C.DE_P_COL] < C.R_P_MAX]) & gene_space
            de_sig[c] = sig
            de_info[c] = dict(contrast=con, n_de_genes=len(sig), n_de_landmark=len(sig & land))
    log(f"DE gene sets: {len(de_sig)} of {P.cohort.nunique()} cohorts")

    # --- eligibility rows ---------------------------------------------------------------------------------------------
    covd = cov.set_index("unit_id")
    out, src_check = [], []
    for r in P.itertuples(index=False):
        lines = drug_lines.get(r.pert_id, [])
        kl = key_lines.get(r.cohort, {})
        direct = [l for l in lines if l in kl]
        src = "direct" if direct else ("borrowed" if len(lines) >= 3 else "none")
        src_check.append((r.source, src))
        use = direct if direct else (lines if src == "borrowed" else [])
        wts = np.array([kl.get(l, 0.0) for l in use]) if direct else np.ones(len(use))
        wts = wts / wts.sum() if wts.sum() > 0 else np.full(len(use), 1 / max(len(use), 1))
        uids = units_by_name.get(r.cmap_name, [])
        base = dict(cohort=r.cohort, level=r.level, parent_group=r.parent_group, pert_id=r.pert_id, cmap_name=r.cmap_name,
                    source_prod=r.source, source_a0=src, passes_gate=r.passes_gate, line_match=r.line_match,
                    transfer_conf=r.transfer_conf)
        if not uids:
            reason = "no_protein_target" if r.no_protein_target is True or r.no_protein_target == "True" else "no_target_unit"
            for l, w in zip(use or [""], wts if len(use) else [np.nan]):
                out.append(dict(base, ModelID=l, line_weight=w, unit_id="", status="unscorable", M_status="unscorable",
                                R_status="unscorable", reason=reason))
            continue
        for u in uids:
            cu = covd.loc[u]
            fpu = fp[u]
            sig = de_sig.get(r.cohort)
            nR = len(fpu["genes_any"] & sig) if sig is not None else 0
            nR_l = len(fpu["genes_any"] & sig & land) if sig is not None else 0
            R_ok = sig is not None and nR >= C.FLOOR_GENES
            R_reason = "" if R_ok else ("no_de_contrast" if sig is None else "R_below_floor")
            for l, w in zip(use or [""], wts if len(use) else [np.nan]):
                i = prow.get((r.pert_id, l))
                null = i is None or not np.isfinite(prof_sd[i]) or prof_sd[i] == 0
                M_s = "eligible" if cu.M_unit_status == "eligible" and not null else "unscorable"
                R_s = "eligible" if R_ok and not null else "unscorable"
                reasons = [x for x in ([] if M_s == "eligible" else [f"M:{cu.M_unit_status if not null else 'null_profile'}"]) +
                           ([] if R_s == "eligible" else [f"R:{R_reason if not null else 'null_profile'}"])]
                if not r.passes_gate:
                    status, reasons = "excluded", ["fails_production_target_gate"] + reasons
                elif src == "none":
                    status, reasons = "unscorable", ["no_source_line"] + reasons
                else:
                    status = "eligible" if (M_s == "eligible" or R_s == "eligible") else "unscorable"
                out.append(dict(base, ModelID=l, line_weight=round(float(w), 4), n_sig=nsig.get((r.pert_id, l)),
                                unit_id=u, unit_kind=cu.kind, action_flag=cu.action_flag, n_sign_resolved=cu.n_sign_resolved,
                                n_anchored_tfs=cu.n_anchored_tfs, n_R_genes=nR, n_R_genes_landmark=nR_l,
                                M_status=M_s, R_status=R_s, status=status, reason=";".join(reasons)))
    E = pd.DataFrame(out)
    E.to_csv(A0 / "eligibility.tsv.gz", sep="\t", index=False, compression="gzip")
    sc = pd.Series([f"{a}->{b}" for a, b in src_check]).value_counts()
    log(f"eligibility rows {len(E)}; status {E.status.value_counts().to_dict()}; source prod->a0 {sc.to_dict()}")

    # summary
    E["reason_main"] = E.reason.fillna("").str.split(";").str[0].replace("", "eligible")
    summ = E.groupby(["level", "status", "M_status", "R_status", "reason_main"]).size().rename("n_rows").reset_index()
    pairs = E.groupby(["level", "cohort", "pert_id"]).status.agg(lambda s: "eligible" if (s == "eligible").any() else
                                                                 ("excluded" if (s == "excluded").all() else "unscorable"))
    ps = pairs.groupby(level=0).value_counts().rename("n_pairs").reset_index()
    summ.to_csv(A0 / "eligibility_summary.tsv", sep="\t", index=False)
    ps.to_csv(A0 / "eligibility_pairs_summary.tsv", sep="\t", index=False)

    # floors (P10)
    fa = []
    for m in C.FLOOR_MULTIPLIERS:
        for space in ("all", "landmark"):
            okm = set(cov.unit_id[cov[f"floor_x{m:g}" + ("" if space == "all" else "_landmark")] &
                                  ~cov.action_flag.isin(C.M_AMBIGUOUS_FLAGS)])
            Eg = E[E.passes_gate.astype(bool) & (E.unit_id != "")]
            rok = Eg.n_R_genes >= C.FLOOR_GENES * m if space == "all" else Eg.n_R_genes_landmark >= C.FLOOR_GENES * m
            fa.append(dict(multiplier=m, feature_space=space, units_M_admitted=len(okm), units_total=len(cov),
                           rows_M_admitted=int(Eg.unit_id.isin(okm).sum()), rows_R_admitted=int(rok.sum()),
                           rows_gate_with_unit=len(Eg),
                           pairs_M_admitted=int(Eg[Eg.unit_id.isin(okm)][["cohort", "pert_id"]].drop_duplicates().shape[0]),
                           pairs_R_admitted=int(Eg[rok][["cohort", "pert_id"]].drop_duplicates().shape[0])))
    pd.DataFrame(fa).to_csv(A0 / "floor_admission.tsv", sep="\t", index=False)

    # cohort feasibility
    cf = []
    for c, d in E.groupby("cohort"):
        pr = d[["pert_id", "source_a0", "passes_gate"]].drop_duplicates("pert_id")
        el = d[d.status == "eligible"]
        cf.append(dict(cohort=c, level=d.level.iloc[0], parent_group=d.parent_group.iloc[0], n_pairs=len(pr),
                       n_gate=int(pr.passes_gate.astype(bool).sum()), n_direct=int((pr.source_a0 == "direct").sum()),
                       n_borrowed=int((pr.source_a0 == "borrowed").sum()), n_in_key_lines=len(key_lines.get(c, {})),
                       best_line_match=d.line_match.max(), **de_info.get(c, {}),
                       pairs_M_eligible=el[el.M_status == "eligible"].pert_id.nunique(),
                       pairs_R_eligible=el[el.R_status == "eligible"].pert_id.nunique(),
                       pairs_any_eligible=el.pert_id.nunique()))
    pd.DataFrame(cf).to_csv(A0 / "cohort_feasibility.tsv", sep="\t", index=False)
    log("cohort feasibility written")

    # --- source-line context (DepMap) ----------------------------------------------------------------------------------
    models = pd.read_csv(C.DEPMAP / "Model.csv", usecols=["ModelID", "StrippedCellLineName", "OncotreeLineage", "OncotreePrimaryDisease"])
    lines = idx.groupby("ModelID").agg(n_profiles=("pert_id", "size"), n_sig=("n_sig", "sum")).reset_index().merge(models, on="ModelID", how="left")
    hot = IO.depmap_matrix(C.DEPMAP / "OmicsSomaticMutationsMatrixHotspot.csv", set(C.PANEL))
    dam = IO.depmap_matrix(C.DEPMAP / "OmicsSomaticMutationsMatrixDamaging.csv", set(C.PANEL))
    cn = IO.depmap_matrix(C.DEPMAP / "OmicsCNGene.csv", {"ERBB2", "RB1", "CDKN2A", "ESR1"})
    ex = IO.depmap_matrix(C.DEPMAP / "OmicsExpressionProteinCodingGenesTPMLogp1.csv", {"ESR1", "ERBB2", "RB1", "AR"})
    mut = pd.read_csv(C.DEPMAP / "OmicsSomaticMutations.csv", usecols=["ModelID", "HugoSymbol", "ProteinChange"], low_memory=False)
    mut = mut[mut.HugoSymbol.isin(C.PANEL) & mut.ModelID.isin(lines.ModelID)]
    pchg = mut.groupby(["ModelID", "HugoSymbol"]).ProteinChange.agg(lambda s: ",".join(sorted(set(map(str, s)))))
    for m in lines.ModelID:
        for tab, nm in ((hot, "hotspot"), (cn, "cn")):
            if m not in tab.index:
                unmatched.append(dict(kind=f"line_not_in_depmap_{nm}", id=m, detail=""))
    for gname in C.PANEL:
        lines[f"{gname}_hot"] = lines.ModelID.map(hot[gname]) if gname in hot else np.nan
        lines[f"{gname}_dam"] = lines.ModelID.map(dam[gname]) if gname in dam else np.nan
        lines[f"{gname}_protein"] = [pchg.get((m, gname), "") for m in lines.ModelID]
    for gname in ("ERBB2", "RB1", "CDKN2A", "ESR1"):
        lines[f"{gname}_cn"] = lines.ModelID.map(cn[gname]).round(3)
        allv = cn[gname].dropna()
        lines[f"{gname}_cn_pct_all_depmap"] = lines[f"{gname}_cn"].map(lambda x: round(float((allv < x).mean() * 100), 1) if x == x else np.nan)
    for gname in ("ESR1", "ERBB2", "RB1", "AR"):
        lines[f"{gname}_tpm_log1p"] = lines.ModelID.map(ex[gname]).round(2)
    lines["RB1_null_candidate"] = (lines.RB1_dam.fillna(0) > 0) | (lines.RB1_tpm_log1p < 1)
    # Absolute relative-copy-number cut, not a percentile: the DepMap ERBB2 distribution is bimodal above it
    # (q0.95 1.57, q0.98 3.61), and a 99th-percentile rule missed SKBR3 at CN 11.0 / percentile 98.7. Changed before
    # any score was computed; the percentile stays as a reported column (D128-A0 result).
    lines["ERBB2_amp_candidate"] = lines.ERBB2_cn >= 2.0
    egfr_sens = re.compile(r"p\.(L858R|L861Q|G719[A-Z]|E746_A750del|.*74[5-9].*del|.*75[0-9].*del|.*7[67][0-9].*ins)")
    lines["EGFR_sensitizing"] = lines.EGFR_protein.map(lambda s: bool(egfr_sens.search(s)) if s else False)
    lines["EGFR_T790M"] = lines.EGFR_protein.str.contains("T790M")
    lines["BRAF_V600"] = lines.BRAF_protein.str.contains("V600")
    lines["KRAS_NRAS_hotspot"] = lines.KRAS_protein.str.contains(r"G12|G13|Q61") | lines.NRAS_protein.str.contains(r"G12|G13|Q61")
    lines["BRCA1_2_damaging"] = (lines.BRCA1_dam.fillna(0) > 0) | (lines.BRCA2_dam.fillna(0) > 0)
    lines.sort_values("n_profiles", ascending=False).to_csv(A0 / "source_line_context.tsv", sep="\t", index=False)
    log(f"source-line context: {len(lines)} lines; RB1-null candidates {lines.RB1_null_candidate.sum()}; "
        f"ERBB2-amp candidates {lines.ERBB2_amp_candidate.sum()}; EGFR sensitizing {lines.EGFR_sensitizing.sum()}")

    # --- controls (availability only) ------------------------------------------------------------------------------------
    lc = lines.set_index("ModelID")
    ctx = {"ESR1": lambda m: lc.loc[m, "ESR1_tpm_log1p"] >= 5,
           "ERBB2": lambda m: bool(lc.loc[m, "ERBB2_amp_candidate"]),
           "BRAF": lambda m: bool(lc.loc[m, "BRAF_V600"]), "MEK": lambda m: bool(lc.loc[m, "BRAF_V600"] or lc.loc[m, "KRAS_NRAS_hotspot"]),
           "EGFR": lambda m: bool(lc.loc[m, "EGFR_sensitizing"]) and not bool(lc.loc[m, "EGFR_T790M"]),
           "CDK4/6": lambda m: not bool(lc.loc[m, "RB1_null_candidate"]) and lc.loc[m, "OncotreeLineage"] == "Breast",
           "PARP": lambda m: bool(lc.loc[m, "BRCA1_2_damaging"])}
    ctl = []
    inpop = set(P.cmap_name)
    for cls, names in C.CONTROLS.items():
        for n in names:
            d = idx[idx.cmap_name.str.lower() == n.lower()]
            us = cov[cov.cmap_name.str.lower() == n.lower()]
            in_ctx = [m for m in d.ModelID.unique() if m in lc.index and ctx[cls](m)]
            ctl.append(dict(cls=cls, compound=n, n_pert_ids=d.pert_id.nunique(), pert_ids="|".join(sorted(d.pert_id.unique())),
                            n_lines=d.ModelID.nunique(), n_sig_total=int(d.n_sig.sum()),
                            in_context_lines="|".join(sorted(f"{lc.loc[m, 'StrippedCellLineName']}({int(d[d.ModelID == m].n_sig.sum())})" for m in in_ctx)),
                            n_sig_in_context=int(d[d.ModelID.isin(in_ctx)].n_sig.sum()),
                            in_lead_population=n in inpop or n.lower() in {x.lower() for x in inpop if isinstance(x, str)},
                            units="|".join(f"{x.genes}:{x.orient}{'/' + x.action_flag if x.action_flag else ''}" for x in us.itertuples()),
                            unit_kind="|".join(us.kind), n_sign_resolved="|".join(map(str, us.n_sign_resolved)),
                            n_resolved_landmark="|".join(map(str, us.n_resolved_landmark)),
                            n_anchored_tfs="|".join(map(str, us.n_anchored_tfs)), M_unit_status="|".join(us.M_unit_status)))
    CT = pd.DataFrame(ctl)
    CT.to_csv(A0 / "control_availability.tsv", sep="\t", index=False)

    # --- P12: pert_id multiplicity ------------------------------------------------------------------------------------
    pm = idx.groupby("cmap_name").agg(n_pert_ids=("pert_id", "nunique"), n_lines=("ModelID", "nunique"), n_sig=("n_sig", "sum"))
    pm = pm[pm.n_pert_ids > 1].reset_index()
    pm["in_lead_population"] = pm.cmap_name.isin(inpop)
    pop_mult = P.groupby(["cohort", "cmap_name"]).pert_id.nunique()
    pm["max_pert_ids_per_cohort_in_population"] = pm.cmap_name.map(pop_mult.groupby(level=1).max())
    pm.sort_values("n_pert_ids", ascending=False).to_csv(A0 / "pert_multiplicity.tsv", sep="\t", index=False)

    # --- TCGA patient strata ----------------------------------------------------------------------------------------------
    surv = pd.read_csv(C.TCGA / "Survival_SupplementalTable_S1_20171025_xena_sp", sep="\t", usecols=["_PATIENT", "cancer type abbreviation"])
    ctype = dict(zip(surv._PATIENT, surv["cancer type abbreviation"]))
    maf = pd.read_csv(C.TCGA / "mc3.v0.2.8.PUBLIC.maf.gz", sep="\t", usecols=["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode",
                      "HGVSp_Short", "Protein_position", "FILTER"], low_memory=False)
    maf["patient"] = maf.Tumor_Sample_Barcode.str[:12]
    maf["ctype"] = maf.patient.map(ctype)
    seq = maf.groupby("ctype").patient.unique().map(set)
    m = maf[(maf.FILTER == "PASS") & maf.Hugo_Symbol.isin(["BRAF", "NRAS", "KRAS", "HRAS", "EGFR", "RB1", "ERBB2"])]
    pos = pd.to_numeric(m.Protein_position.astype(str).str.split("/").str[0].str.split("-").str[0], errors="coerce")
    hs = m.HGVSp_Short.fillna("")
    sets = {
        "BRAF_V600": set(m.patient[(m.Hugo_Symbol == "BRAF") & hs.str.startswith("p.V600")]),
        "RAS_hotspot": set(m.patient[m.Hugo_Symbol.isin(["NRAS", "KRAS", "HRAS"]) & hs.str.contains(r"p\.(?:G12|G13|Q61)")]),
        "NRAS_hotspot": set(m.patient[(m.Hugo_Symbol == "NRAS") & hs.str.contains(r"p\.(?:G12|G13|Q61)")]),
        "KRAS_hotspot": set(m.patient[(m.Hugo_Symbol == "KRAS") & hs.str.contains(r"p\.(?:G12|G13|Q61)")]),
        "EGFR_activating": set(m.patient[(m.Hugo_Symbol == "EGFR") & (
            hs.str.contains(r"p\.(?:L858R|L861Q|G719)") |
            ((m.Variant_Classification == "In_Frame_Del") & pos.between(745, 759)) |
            ((m.Variant_Classification == "In_Frame_Ins") & pos.between(762, 775)))]),
        "EGFR_T790M": set(m.patient[(m.Hugo_Symbol == "EGFR") & hs.str.contains("T790M")]),
        "RB1_truncating": set(m.patient[(m.Hugo_Symbol == "RB1") & m.Variant_Classification.isin(
            ["Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"])]),
    }
    gis = pd.read_csv(C.TCGA / "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz", sep="\t", index_col=0)
    gis = gis.loc[["ERBB2", "RB1"]]
    gpat = {c[:12]: c for c in gis.columns if c[13:15] in ("01", "03", "06")}
    erbb2_amp = {p for p, c in gpat.items() if gis.at["ERBB2", c] == 2}
    rb1_del = {p for p, c in gpat.items() if gis.at["RB1", c] == -2}
    cn_pat = set(gpat)
    sub = pd.read_csv(C.TCGA / "TCGASubtype.20170308.tsv.gz", sep="\t", usecols=["sampleID", "Subtype_mRNA", "Subtype_Selected"])
    sub["patient"] = sub.sampleID.str[:12]
    brca_sub = dict(zip(sub.patient[sub.Subtype_Selected.astype(str).str.startswith("BRCA")], sub.Subtype_mRNA))
    st = []
    def add(cohort, stratum, pats, base, note=""):
        st.append(dict(cohort=cohort, stratum=stratum, n_patients=len(pats & base), n_base=len(base), note=note))
    sk = seq.get("SKCM", set())
    add("SKCM", "BRAF_V600", sets["BRAF_V600"], sk)
    add("SKCM", "NRAS_hotspot_BRAF_wt", sets["NRAS_hotspot"] - sets["BRAF_V600"], sk)
    add("SKCM", "BRAF_wt_RAS_wt", sk - sets["BRAF_V600"] - sets["RAS_hotspot"], sk, "neither BRAF V600 nor RAS hotspot")
    lu = seq.get("LUAD", set())
    add("LUAD", "EGFR_activating", sets["EGFR_activating"], lu)
    add("LUAD", "EGFR_T790M", sets["EGFR_T790M"], lu)
    add("LUAD", "KRAS_hotspot_EGFR_wt", sets["KRAS_hotspot"] - sets["EGFR_activating"], lu)
    add("LUAD", "EGFR_wt_KRAS_wt", lu - sets["EGFR_activating"] - sets["KRAS_hotspot"], lu)
    th = seq.get("THCA", set())
    add("THCA", "BRAF_V600", sets["BRAF_V600"], th)
    add("THCA", "RAS_hotspot", sets["RAS_hotspot"] - sets["BRAF_V600"], th)
    br_cn = {p for p in cn_pat if ctype.get(p) == "BRCA"}
    add("BRCA", "ERBB2_amp_GISTIC", erbb2_amp, br_cn)
    add("BRCA", "RB1_loss(deep_del_or_truncating)", rb1_del | sets["RB1_truncating"], br_cn | seq.get("BRCA", set()))
    for sname in ("LumA", "LumB", "Her2", "Basal", "Normal"):
        sp = {p for p, s in brca_sub.items() if s == sname}
        add("BRCA", f"PAM50_{sname}:ERBB2_amp", erbb2_amp & sp, br_cn & sp)
        add("BRCA", f"PAM50_{sname}:RB1_loss", (rb1_del | sets["RB1_truncating"]) & sp, (br_cn | seq.get("BRCA", set())) & sp)
    pd.DataFrame(st).to_csv(A0 / "patient_strata.tsv", sep="\t", index=False)
    log("patient strata written")

    # --- waves (Section 5) ---------------------------------------------------------------------------------------------------
    def ctx_sig(cls):
        return int(CT[CT.cls == cls].n_sig_in_context.sum())
    def ctx_lines(cls):
        return "; ".join(f"{r.compound}: {r.in_context_lines}" for r in CT[(CT.cls == cls) & (CT.in_context_lines != "")].itertuples())
    def unit_ok(cls):
        return "; ".join(f"{r.compound}: {r.M_unit_status}" for r in CT[(CT.cls == cls) & (CT.M_unit_status != "")].itertuples())
    stmap = {(r["cohort"], r["stratum"]): r["n_patients"] for r in st}
    waves = [
        dict(wave="A1", cohorts="BRCA-LumA/LumB vs BRCA-Basal", cls="ESR1", in_context_signatures=ctx_sig("ESR1"),
             in_context=ctx_lines("ESR1"), unit_status=unit_ok("ESR1"), patients="existing cohorts"),
        dict(wave="A2", cohorts="BRCA-Her2 vs BRCA-Basal", cls="ERBB2", in_context_signatures=ctx_sig("ERBB2"),
             in_context=ctx_lines("ERBB2"), unit_status=unit_ok("ERBB2"), patients=f"ERBB2 amp (GISTIC) {stmap[('BRCA', 'ERBB2_amp_GISTIC')]}"),
        dict(wave="A3", cohorts="SKCM BRAF V600 vs WT/NRAS", cls="BRAF", in_context_signatures=ctx_sig("BRAF"),
             in_context=ctx_lines("BRAF"), unit_status=unit_ok("BRAF"),
             patients=f"V600 {stmap[('SKCM', 'BRAF_V600')]}, NRAS {stmap[('SKCM', 'NRAS_hotspot_BRAF_wt')]}, WT {stmap[('SKCM', 'BRAF_wt_RAS_wt')]}"),
        dict(wave="A4", cohorts="LUAD EGFR activating vs WT/KRAS", cls="EGFR", in_context_signatures=ctx_sig("EGFR"),
             in_context=ctx_lines("EGFR"), unit_status=unit_ok("EGFR"),
             patients=f"EGFR act {stmap[('LUAD', 'EGFR_activating')]}, KRAS {stmap[('LUAD', 'KRAS_hotspot_EGFR_wt')]}, WT {stmap[('LUAD', 'EGFR_wt_KRAS_wt')]}"),
        dict(wave="B1", cohorts="THCA BRAF-like vs RAS-like", cls="BRAF|MEK", in_context_signatures=ctx_sig("BRAF") + ctx_sig("MEK"),
             in_context=ctx_lines("MEK"), unit_status=unit_ok("MEK"),
             patients=f"BRAF V600 {stmap[('THCA', 'BRAF_V600')]}, RAS {stmap[('THCA', 'RAS_hotspot')]}"),
        dict(wave="B2", cohorts="BRCA RB1-intact vs RB1-deficient", cls="CDK4/6", in_context_signatures=ctx_sig("CDK4/6"),
             in_context=ctx_lines("CDK4/6"), unit_status=unit_ok("CDK4/6"),
             patients=f"RB1 loss {stmap[('BRCA', 'RB1_loss(deep_del_or_truncating)')]}"),
        dict(wave="B3", cohorts="OV/BRCA HRD vs HR-proficient", cls="PARP", in_context_signatures=ctx_sig("PARP"),
             in_context=ctx_lines("PARP"), unit_status=unit_ok("PARP"), patients="HRD status not audited (no HRD table on disk)"),
    ]
    pd.DataFrame(waves).to_csv(A0 / "wave_feasibility.tsv", sep="\t", index=False)

    # --- A0 action 7: external screens, prior use (grep only) -------------------------------------------------------------
    ext = []
    pat = re.compile(r"prism|gdsc|ctrp", re.I)
    for f in sorted((C.ROOT / "scripts").rglob("*")):
        if f.is_file() and f.suffix in (".py", ".R", ".sh", ".html") and "d128" not in f.parts:
            txt = f.read_text(errors="ignore")
            hits = [l.strip() for l in txt.splitlines() if pat.search(l)]
            if hits:
                doc = txt.split('"""')[1].strip().splitlines()[0] if '"""' in txt else ""
                ext.append(dict(kind="script", where=str(f.relative_to(C.ROOT)), n_mentions=len(hits),
                                screens="|".join(sorted({m.lower() for h in hits for m in pat.findall(h)})),
                                first_line_of_docstring=doc[:200], example=hits[0][:200]))
    for l in (C.ROOT / "docs/DECISIONS.md").read_text().splitlines():
        if l.startswith("|") and pat.search(l):
            did = l.split("|")[1].strip()
            if did.startswith("D128"):
                continue
            ms = [mm.start() for mm in pat.finditer(l)]
            ext.append(dict(kind="decision", where=did, n_mentions=len(ms),
                            screens="|".join(sorted({x.lower() for x in pat.findall(l)})),
                            first_line_of_docstring="", example=l[max(0, ms[0] - 120): ms[0] + 120].replace("\n", " ")))
    pd.DataFrame(ext).to_csv(A0 / "external_screen_prior_use.tsv", sep="\t", index=False)

    # --- unmatched, manifest --------------------------------------------------------------------------------------------------
    for c in sorted(set(P.cohort) - set(key_lines)):
        unmatched.append(dict(kind="cohort_without_in_key_lines", id=c, detail=""))
    for n in sorted(set(P.cmap_name.dropna()) - set(U.cmap_name)):
        unmatched.append(dict(kind="population_compound_without_d121_unit", id=n, detail=""))
    um = pd.DataFrame(unmatched).drop_duplicates()
    um.to_csv(A0 / "unmatched_ids.tsv", sep="\t", index=False)

    inputs = [C.PROFILES_SPECIFIC, C.PROFILES_RAW, C.PROFILE_INDEX, C.PROFILE_GENES, C.DRUG_ANNOT, C.LAYERED, C.SPLIT_GROUPS,
              C.SPLIT_STRATA, C.UNITS, C.LINE_W_GROUPS, C.LINE_W_STRATA, C.COLLECTRI, C.OMNIPATH, C.SIGNOR]
    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=C.ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "scripts/d128"], cwd=C.ROOT, capture_output=True, text=True).stdout.strip()
    man = dict(stage="D128-A0", git_sha=sha, scripts_dirty=bool(dirty), frozen=C.FROZEN, threads=a.threads,
               inputs={str(p.relative_to(C.ROOT)): dict(sha256=IO.sha256(p), bytes=p.stat().st_size) for p in inputs},
               config={k: (v if isinstance(v, (int, float, str, bool, list, tuple, dict)) else str(v))
                       for k, v in vars(C).items() if k.isupper() and not isinstance(v, Path)},
               counts=dict(population_pairs=len(P), eligibility_rows=len(E), status=E.status.value_counts().to_dict(),
                           units=len(cov), source_prod_vs_a0=sc.to_dict()),
               log=LOG, seconds=round(time.time() - T0, 1))
    (A0 / "run_manifest.json").write_text(json.dumps(man, indent=1, default=str))
    log("done")


if __name__ == "__main__":
    main()
