#!/usr/bin/env python
"""D126 run-order B1: signed path consistency, all cohorts (spec B1/B3; choices in D126-B1).

For each (FDA drug, cohort) that is D122-eligible, has a curated main target (D105b `evidence_targets`) and a known
action direction: take the TOP_K regulators the drug moves hardest in that cohort's consensus profile (the vector
79e scores, dumped by its D126_ACT_DIR option), ask whether a signed path runs from the drug's main targets to each
of them in the B0 prior, and compare the propagated sign -- drug action x product of edge signs -- with the sign the
drug actually produced.

    mechanism_coverage        = reachable / TOP_K          (reachable = reached with a resolved sign)
    mechanism_consistent_frac = consistent / reachable
    mechanism_inconsistent    = coverage >= FLAG_COVERAGE_MIN and consistent_frac < FLAG_CONSISTENT_MAX

Regulators reached only by shortest paths that disagree on net sign are `sign_ambiguous`: excluded from both the
numerator and the denominator, and counted. A regulator that *is* a main target has distance 0 and sign +1, so the
predicted sign is the drug's own action -- which is the case B5.1 turns on.

Reported, never ranked or gated (D120 / D121 convention, spec B5): per-drug and per-cohort tables, the path-length
distribution, the single-main-target subset, and the same run on the SIGNOR prior as the sensitivity arm.

Outputs: results/lincs/d126/{b1_rows,b1_drugs,b1_paths,b3_criteria}.tsv, b1_run.json.
"""
import json
import re
import sys
from collections import deque
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fda_rows as F  # noqa: E402
from b1_config import (ACT_DIR, ANTI_ESTROGEN_TARGETS, B51_CONTROL_MIN, B52_RHO_MAX, B53_COVERAGE_MIN,  # noqa: E402
                       B53_SHARE_MIN, FLAG_CONSISTENT_MAX, FLAG_COVERAGE_MIN, MCF7_LINE, OUT, PRIOR, PRIOR_SIGNOR,
                       TOP_K)


# ---- the prior ------------------------------------------------------------------------------------------------
def load_prior(path):
    e = pd.read_csv(path, sep="\t")
    adj = {}
    for s, t, sg in zip(e.source_genesymbol, e.target_genesymbol, e.sign):
        adj.setdefault(str(s), []).append((str(t), int(sg)))
    out_deg = {n: len(v) for n, v in adj.items()}
    nodes = set(adj) | {t for v in adj.values() for t, _ in v}
    tot_deg = dict(out_deg)
    for v in adj.values():
        for t, _ in v:
            tot_deg[t] = tot_deg.get(t, 0) + 1
    return adj, out_deg, tot_deg, nodes


def reach(adj, targets):
    """Multi-source BFS with sign propagation. Returns {node: (distance, frozenset of net signs at that distance)}.

    A target reaches itself at distance 0 with sign +1 (the empty path). Signs are collected only at the shortest
    distance, so a longer path never overrides or joins a shorter one (D126-B1 (4))."""
    seen = {t: (0, {1}) for t in targets}
    q = deque((t, 0) for t in targets)
    while q:
        u, d = q.popleft()
        if seen[u][0] != d:
            continue
        for v, sg in adj.get(u, ()):
            ns = {x * sg for x in seen[u][1]}
            if v not in seen:
                seen[v] = (d + 1, set(ns))
                q.append((v, d + 1))
            elif seen[v][0] == d + 1:
                seen[v][1].update(ns)
    return {k: (d, frozenset(s)) for k, (d, s) in seen.items()}


# ---- the drugs ------------------------------------------------------------------------------------------------
def target_genes(x):
    """D105b evidence targets, written 'GENE:value/value', separated by | or ; (as 79t / d127 read them)."""
    if not isinstance(x, str):
        return []
    return [p.split(":")[0].strip() for p in re.split(r"[|;]", x) if p.strip()]


def action_sign(moa):
    """D126-B1 (3): +1 agonist / activator, -1 antagonist / inhibitor / degrader / blocker, as 79v ligands() reads it."""
    moa = (moa or "").lower()
    if not moa:
        return None
    ago = bool(re.search(r"agonist|activator", moa.replace("antagonist", "")))
    ant = bool(re.search(r"antagonist|inhibitor|degrader|blocker", moa))
    return None if ago == ant else (1 if ago else -1)


def drug_table():
    g, s = F.fda_rows()
    a = pd.concat([g.assign(kind="group"), s.assign(kind="subgroup")], ignore_index=True)
    a = a[F.eligible(a)].copy()
    ann = pd.read_csv(ROOT / "results/lincs/drug_annotations.tsv", sep="\t", dtype=str).drop_duplicates("cmap_name")
    for c in ["hub_moa", "lincs_moa"]:
        if c not in ann:
            ann[c] = None
    moa = {n: " ".join(x for x in [h, l] if isinstance(x, str))
           for n, h, l in zip(ann.cmap_name, ann.hub_moa, ann.lincs_moa)}
    a["targets"] = a.evidence_targets.map(target_genes)
    a["action"] = a.cmap_name.map(lambda n: action_sign(moa.get(n)))
    return a


# ---- B1 -----------------------------------------------------------------------------------------------------
def run_arm(a, adj, out_deg, nodes, label):
    rows, paths = [], []
    cache = {}
    for gid, gd in a.groupby("group_id"):
        f = ACT_DIR / f"{gid}.f32"
        if not f.exists():
            continue
        regs = np.array((ACT_DIR / f"{gid}.regulators.txt").read_text().split())
        drugs = pd.read_csv(ACT_DIR / f"{gid}.drugs.tsv", sep="\t")
        D = np.fromfile(f, dtype=np.float32).reshape(len(drugs), len(regs))
        pos = {p: i for i, p in enumerate(drugs.pert_id)}
        for r in gd.itertuples():
            if r.pert_id not in pos or r.action is None or not r.targets:
                continue
            tset = tuple(sorted({t for t in r.targets if t in nodes}))
            if not tset:
                continue
            if tset not in cache:
                cache[tset] = reach(adj, [t for t in tset if t in adj or t in nodes])
            R = cache[tset]
            v = D[pos[r.pert_id]]
            top = np.argsort(-np.abs(v))[:TOP_K]
            n_cons = n_reach = n_amb = 0
            for i in top:
                reg, obs = str(regs[i]), int(np.sign(v[i]))
                hit = R.get(reg)
                if hit is None or obs == 0:
                    continue
                d, signs = hit
                if len(signs) != 1:
                    n_amb += 1
                    continue
                n_reach += 1
                pred = r.action * next(iter(signs))
                n_cons += int(pred == obs)
                paths.append((label, gid, r.cmap_name, reg, d, pred, obs, int(pred == obs)))
            cov = n_reach / TOP_K
            frac = n_cons / n_reach if n_reach else np.nan
            rows.append(dict(arm=label, group_id=gid, kind=r.kind, drug=r.cmap_name, pert_id=r.pert_id,
                             n_targets=len(r.targets), n_targets_in_prior=len(tset), action=r.action,
                             mechanism_coverage=round(cov, 4), n_reachable=n_reach, n_consistent=n_cons,
                             n_sign_ambiguous=n_amb,
                             mechanism_consistent_frac=round(frac, 4) if n_reach else np.nan,
                             mechanism_inconsistent=bool(cov >= FLAG_COVERAGE_MIN and n_reach
                                                         and frac < FLAG_CONSISTENT_MAX),
                             max_out_degree=max(out_deg.get(t, 0) for t in tset)))
    return pd.DataFrame(rows), pd.DataFrame(paths, columns=["arm", "group_id", "drug", "regulator", "hops",
                                                            "predicted", "observed", "consistent"])


# ---- B3: the criteria -----------------------------------------------------------------------------------------
def mcf7_breast_groups():
    """Breast cohorts whose in-key LINCS lines include MCF7 (D126-B1 (7))."""
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#")
    models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv")
    mid = set(models.loc[models.CellLineName.astype(str).str.upper().str.replace("-", "") == MCF7_LINE, "ModelID"])
    if not mid:
        mid = set(models.loc[models.StrippedCellLineName.astype(str).str.upper() == MCF7_LINE, "ModelID"]) \
            if "StrippedCellLineName" in models else set()
    lin = models.set_index("ModelID").OncotreeLineage.to_dict()
    lineages = {lin.get(m) for m in mid}
    idx = pd.read_csv(ROOT / "data/processed/lincs/drug_line_profiles_index.tsv", sep="\t")
    has_lincs = bool(mid & set(idx.ModelID))
    out = []
    for r in keys.itertuples():
        ls = set(str(r.lineages).split(";")) if isinstance(r.lineages, str) else set()
        if ls & lineages and has_lincs:
            out.append(r.group_id)
    return set(out), sorted(mid)


def criteria(rows, a, out_deg):
    g = rows[rows.arm == "omnipath"]
    res = []

    # B5.1 known-failure control
    ae = set(a[a.targets.map(lambda L: bool(set(L) & ANTI_ESTROGEN_TARGETS))].cmap_name)
    breast, mid = mcf7_breast_groups()
    ctl = g[g.drug.isin(ae) & g.group_id.isin(breast)]
    share = float(ctl.mechanism_inconsistent.mean()) if len(ctl) else np.nan
    res.append(dict(criterion="B5.1 known-failure control (anti-estrogens, MCF7 breast cohorts)",
                    value=f"flag fires on {int(ctl.mechanism_inconsistent.sum())} of {len(ctl)} rows"
                          f" ({share:.3f})" if len(ctl) else "no scorable rows",
                    threshold=f"> {B51_CONTROL_MIN}", met=bool(len(ctl) and share > B51_CONTROL_MIN)))

    # B5.2 degree confound
    per = g.groupby("drug").agg(flag=("mechanism_inconsistent", "mean"), deg=("max_out_degree", "max")).dropna()
    rho = float(spearmanr(per.flag, per.deg).correlation) if len(per) > 2 else np.nan
    res.append(dict(criterion="B5.2 degree confound (flag rate vs main-target out-degree)",
                    value=f"Spearman rho {rho:.3f} over {len(per)} drugs",
                    threshold=f"|rho| < {B52_RHO_MAX}", met=bool(abs(rho) < B52_RHO_MAX)))

    # B5.3 coverage floor
    pc = g.groupby("drug").mechanism_coverage.median()
    sh = float((pc >= B53_COVERAGE_MIN).mean()) if len(pc) else np.nan
    res.append(dict(criterion="B5.3 coverage floor",
                    value=f"median coverage >= {B53_COVERAGE_MIN} for {int((pc >= B53_COVERAGE_MIN).sum())}"
                          f" of {len(pc)} scorable drugs ({sh:.3f})",
                    threshold=f">= {B53_SHARE_MIN}", met=bool(sh >= B53_SHARE_MIN)))
    return pd.DataFrame(res), dict(control_rows=int(len(ctl)), control_share=share, rho=rho, coverage_share=sh,
                                   mcf7_models=mid, breast_groups=sorted(breast))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    a = drug_table()
    print(f"eligible FDA rows {len(a)} over {a.group_id.nunique()} cohorts; "
          f"{a[a.action.notna() & a.targets.map(bool)].cmap_name.nunique()} drugs with a target and a direction")

    arms, allrows, allpaths = {}, [], []
    for label, path in [("omnipath", PRIOR), ("signor", PRIOR_SIGNOR)]:
        adj, out_deg, tot_deg, nodes = load_prior(path)
        r, p = run_arm(a, adj, out_deg, nodes, label)
        arms[label] = (r, out_deg)
        allrows.append(r)
        allpaths.append(p)
        print(f"[{label}] {len(r)} (drug, cohort) rows, {r.drug.nunique()} drugs; "
              f"median coverage {r.mechanism_coverage.median():.3f}, "
              f"median consistency {r.mechanism_consistent_frac.median():.3f}, "
              f"flag fires on {int(r.mechanism_inconsistent.sum())} ({r.mechanism_inconsistent.mean():.3f})")

    rows = pd.concat(allrows, ignore_index=True)
    paths = pd.concat(allpaths, ignore_index=True)
    rows.to_csv(OUT / "b1_rows.tsv", sep="\t", index=False)
    paths.to_csv(OUT / "b1_paths.tsv", sep="\t", index=False)

    per = (rows[rows.arm == "omnipath"].groupby("drug")
           .agg(cohorts=("group_id", "nunique"), n_targets=("n_targets", "max"),
                coverage=("mechanism_coverage", "median"), consistency=("mechanism_consistent_frac", "median"),
                flag_rate=("mechanism_inconsistent", "mean"), out_degree=("max_out_degree", "max"))
           .reset_index().sort_values("flag_rate", ascending=False))
    per.to_csv(OUT / "b1_drugs.tsv", sep="\t", index=False, float_format="%.4g")

    crit, detail = criteria(rows, a, arms["omnipath"][1])
    crit.to_csv(OUT / "b3_criteria.tsv", sep="\t", index=False)

    single = rows[(rows.arm == "omnipath") & (rows.n_targets == 1)]
    multi = rows[(rows.arm == "omnipath") & (rows.n_targets > 1)]
    run = dict(top_k=TOP_K, flag=dict(coverage_min=FLAG_COVERAGE_MIN, consistent_max=FLAG_CONSISTENT_MAX),
               criteria=crit.to_dict("records"), detail=detail,
               hops=paths[paths.arm == "omnipath"].hops.value_counts().sort_index().to_dict(),
               single_target=dict(n=len(single), coverage=float(single.mechanism_coverage.median()),
                                  consistency=float(single.mechanism_consistent_frac.median()),
                                  flag=float(single.mechanism_inconsistent.mean())),
               multi_target=dict(n=len(multi), coverage=float(multi.mechanism_coverage.median()),
                                 consistency=float(multi.mechanism_consistent_frac.median()),
                                 flag=float(multi.mechanism_inconsistent.mean())))
    (OUT / "b1_run.json").write_text(json.dumps(run, indent=1, default=float))

    pd.set_option("display.width", 220)
    print("\nB3 criteria:")
    print(crit.to_string(index=False))
    print("\nsingle vs multi main target:", json.dumps({k: run[k] for k in ["single_target", "multi_target"]}, indent=1, default=float))
    print("\npath lengths (hops):", run["hops"])
    print("\nmost-flagged drugs:")
    print(per.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
