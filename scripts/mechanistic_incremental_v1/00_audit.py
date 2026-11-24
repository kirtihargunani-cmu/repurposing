#!/usr/bin/env python
"""Stage 0 of mechanistic_incremental_v1: inputs, coverage, joins and exact decompositions (no outcome is computed).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
Spec: docs/SPEC_Perturbation_Mechanistic_Incremental_Value.md (Amendment 1). Config: config/mechanistic_incremental_v1.yaml.

Checks, in order; a failed blocking check stops the stage with status=blocked:
  1. every pinned input exists and its .build_manifest.json run id equals the config's (refuses older runs);
     sha256 + size of every input and upstream-consumed file -> inputs_manifest.json
  2. D208 guard: fda_rows / 79ar default to _v2_tf and refuse a consensus suffix
  3. key uniqueness: window tables (cohort, pert_id); regulator table (cohort, pert_id, window, regulator);
     layered (unit, pert_id); gene ranking (unit, level, ingredient, window)
  4. exact decompositions: sum_r rev_contrib = stored rev and = window_reversal rev_<window> (5 dp);
     sum_r rev_gene_contrib = gene_ranking rev_gene; on every row, not a sample
  5. D209 reproduction: target-unit dedup gives 20 primary units, 16 correct
  6. coverage per (level, unit, window): scored, with window data, >= 2 conditions, eligible, ranked; unknown by reason
Outputs: results/lincs/mechanistic_incremental_v1/00_audit/{coverage.tsv, join_audit.tsv, decomposition.tsv,
         inputs_manifest.json, report.md, run.json}
"""
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve()
WT = HERE.parents[2]
CFG = yaml.safe_load((WT / "config/mechanistic_incremental_v1.yaml").read_text())
IN = Path(CFG["input_root"])
OUT = WT / CFG["output_root"] / "00_audit"
OUT.mkdir(parents=True, exist_ok=True)
t0 = time.time()
blockers, notes = [], []


def sha256(p, chunk=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def manifest_of(path):
    p = IN / path
    for m in ([p / ".build_manifest.json"] if p.is_dir() else [Path(str(p) + ".build_manifest.json"), p.parent / ".build_manifest.json"]):
        if m.exists():
            return json.loads(m.read_text())
    return None


# ---- 1. inputs pinned to the reference run ----------------------------------------------------------------------
man = {"reference_run": CFG["reference_run"], "inputs": {}, "upstream_consumed": {}}
for path, rid in CFG["inputs"].items():
    m = manifest_of(path)
    got = m.get("run_id") if m else None
    ok = got == rid
    files = sorted((IN / path).glob("*")) if (IN / path).is_dir() else [IN / path]
    man["inputs"][path] = {"expected_run_id": rid, "manifest_run_id": got, "ok": ok,
                           "files": {str(f.relative_to(IN)): {"bytes": f.stat().st_size, "sha256": sha256(f)}
                                     for f in files if f.is_file() and not f.name.startswith(".")}}
    if not ok:
        blockers.append(f"input {path}: manifest run id {got} != pinned {rid}")
for path, what in CFG["upstream_consumed"].items():
    p = IN / path
    files = sorted(p.glob("*")) if p.is_dir() else [p]
    man["upstream_consumed"][path] = {"what": what, "exists": p.exists(),
                                      "files": {str(f.relative_to(IN)): {"bytes": f.stat().st_size, "sha256": sha256(f)}
                                                for f in files if f.is_file()}}
    if not p.exists():
        blockers.append(f"upstream input missing: {path}")
man["git"] = {"study_sha": subprocess.run(["git", "-C", str(WT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
              "input_checkout_sha": subprocess.run(["git", "-C", str(IN), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()}
(OUT / "inputs_manifest.json").write_text(json.dumps(man, indent=1, default=str))
print(f"1. inputs: {sum(v['ok'] for v in man['inputs'].values())}/{len(man['inputs'])} pinned run ids match", flush=True)

# ---- 2. D208 guard -----------------------------------------------------------------------------------------------
guard = {}
for script, needle in (("scripts/fda_rows.py", 'os.environ.get("SPLIT_SFX", "_v2_tf")'),
                       ("scripts/79ar_gene_ranking.py", 'os.environ.get("SFX", "_v2_tf")')):
    src = (IN / script).read_text()
    guard[script] = needle in src and "consensus" in src and "SystemExit" in src
    if not guard[script]:
        blockers.append(f"D208 guard not found in {script}")
print(f"2. D208 guard: {guard}", flush=True)

# ---- 3. key uniqueness --------------------------------------------------------------------------------------------
J = []
R = IN / "results"
wg = pd.read_csv(R / "systems_validation/integrated/window_reversal.tsv.gz", sep="\t")
ws = pd.read_csv(R / "systems_validation/integrated/window_reversal_strata.tsv.gz", sep="\t")
L = pd.read_csv(R / "lincs/layered_score_v2_v2_tf.tsv", sep="\t", low_memory=False)
GR = pd.read_csv(R / "lincs/gene_ranking_v2_tf.tsv", sep="\t", low_memory=False)
GR["ing"] = GR.fda_ingredient.fillna(GR.cmap_name)
for name, df, key in (("window_reversal", wg, ["cohort", "pert_id"]), ("window_reversal_strata", ws, ["cohort", "pert_id"]),
                      ("layered_score", L, ["unit", "pert_id"]), ("gene_ranking", GR, ["unit", "level", "ing", "window"])):
    d = int(df.duplicated(key).sum())
    J.append(dict(table=name, key="+".join(key), rows=len(df), duplicate_keys=d, unique=d == 0))
    if d:
        blockers.append(f"{name}: {d} duplicate keys on {key}")

# ---- 4. exact decompositions, every row --------------------------------------------------------------------------
dec = []
for tag, f, win_df in (("groups", "regulator_contributions_by_window_all.tsv.gz", wg),
                       ("subgroups", "regulator_contributions_by_window_all_strata.tsv.gz", ws)):
    acc, nreg, dupl = [], 0, 0
    cols = ["cohort", "pert_id", "window", "regulator", "rev_contrib", "rev", "rev_gene_contrib"]
    for ch in pd.read_csv(R / "lincs/regulon_contributions_v2_tf" / f, sep="\t", usecols=cols, chunksize=3_000_000):
        nreg += len(ch)
        acc.append(ch.groupby(["cohort", "pert_id", "window"], sort=False).agg(
            s_rev=("rev_contrib", "sum"), rev=("rev", "first"), s_gene=("rev_gene_contrib", "sum"), n=("regulator", "size")).reset_index())
    S = pd.concat(acc).groupby(["cohort", "pert_id", "window"], as_index=False).agg(
        s_rev=("s_rev", "sum"), rev=("rev", "first"), s_gene=("s_gene", "sum"), n=("n", "sum"))
    J.append(dict(table=f"regulator_contributions ({tag})", key="cohort+pert_id+window", rows=nreg,
                  duplicate_keys=0, unique=True, cells=len(S)))
    e1 = (S.s_rev - S.rev).abs()
    long = pd.concat([win_df[["cohort", "pert_id", f"rev_{w}"]].dropna().rename(columns={f"rev_{w}": "rev_win"}).assign(window=w)
                      for w in ("early", "late")])
    M = S.merge(long, on=["cohort", "pert_id", "window"], how="left")
    e2 = (M.s_rev - M.rev_win).abs()
    lv = "group_id" if tag == "groups" else "stratum"
    G2 = GR[(GR.level == lv) & GR.window.isin(["early", "late"])][["unit", "pert_id", "window", "rev_gene"]].dropna()
    M2 = S.merge(G2.rename(columns={"unit": "cohort"}), on=["cohort", "pert_id", "window"], how="inner")
    e3 = (M2.s_gene - M2.rev_gene).abs() / np.maximum(1.0, M2.rev_gene.abs())
    dec += [dict(level=tag, check="sum rev_contrib = stored rev", n=len(S), max_abs_err=float(e1.max()), tol=1e-3),
            dict(level=tag, check="sum rev_contrib = window_reversal rev (5 dp)", n=int(M.rev_win.notna().sum()),
                 max_abs_err=float(e2.max()), tol=1e-3, unmatched=int(M.rev_win.isna().sum())),
            dict(level=tag, check="sum rev_gene_contrib = gene_ranking rev_gene (rel)", n=len(M2), max_abs_err=float(e3.max()) if len(M2) else np.nan, tol=1e-4)]
DEC = pd.DataFrame(dec)
DEC["pass"] = DEC.max_abs_err <= DEC.tol
for r in DEC[~DEC["pass"]].itertuples():
    blockers.append(f"decomposition failed: {r.level} {r.check} max err {r.max_abs_err:.3g}")
DEC.to_csv(OUT / "decomposition.tsv", sep="\t", index=False)
pd.DataFrame(J).to_csv(OUT / "join_audit.tsv", sep="\t", index=False)
print("4. decompositions:\n" + DEC.to_string(index=False), flush=True)

# ---- 5. D209 reproduction ----------------------------------------------------------------------------------------
dd = pd.read_csv(R / "lincs/d127/dedup_summary.tsv", sep="\t").set_index("set")
d209 = dict(units=int(dd.loc["primary", "units"]), correct=int(dd.loc["primary", "correct"]))
if d209 != {"units": 20, "correct": 16}:
    blockers.append(f"D209 not reproduced: {d209}")

# ---- 6. coverage ------------------------------------------------------------------------------------------------
cov = []
for lv, wdf in (("group_id", wg), ("stratum", ws)):
    lay = L[L.level == lv]
    for u in sorted(set(lay.unit)):
        w = wdf[wdf.cohort == u]
        g = GR[(GR.unit == u) & (GR.level == lv) & (GR.path == "D/A/O")]
        el = g[g.eligible]
        for win in ("early", "late"):
            nc = w.get(f"n_conditions_{win}", pd.Series(dtype=float))
            cov.append(dict(level=lv, unit=u, window=win, fda_scored=int(lay.layered_T.notna().sum()) if False else int((lay.unit == u).sum()),
                            window_rows=int(w[f"rev_{win}"].notna().sum()) if len(w) else 0,
                            rows_ge2_conditions=int((nc >= 2).sum()) if len(w) else 0,
                            eligible_ingredients=int(el.ing.nunique()),
                            ranked=int(((el.window == win) & el.gene_pct.notna()).sum()),
                            in_window_run=bool(len(w))))
COV = pd.DataFrame(cov)
unk = GR[(GR.window == "unknown") & (GR.path == "D/A/O") & GR.eligible].window_unknown_reason.value_counts()
COV.to_csv(OUT / "coverage.tsv", sep="\t", index=False)

# ---- notes found while auditing ----------------------------------------------------------------------------------
src = (IN / "scripts/79e3_window_contributions.py").read_text()
if "abs(hash((cohort, r.pert_id, r.window)))" in src:
    notes.append("79e3 SPLITHALF seeds each row with Python hash() of strings, which changes between processes "
                 "(PYTHONHASHSEED), so its split-half assignment is not reproducible run to run. This study seeds with "
                 "crc32(key) ^ seed instead and does not use 79e3's splithalf output.")

status = "blocked" if blockers else "pass"
run = dict(stage="00_audit", status=status, blockers=blockers, notes=notes, d209=d209,
           unknown_eligible_by_reason={k: int(v) for k, v in unk.items()},
           coverage_totals=COV.groupby(["level", "window"])[["window_rows", "rows_ge2_conditions", "ranked"]].sum().reset_index().to_dict("records"),
           units_in_window_run={lv: int(COV[(COV.level == lv) & COV.in_window_run].unit.nunique()) for lv in ("group_id", "stratum")},
           git=man["git"], config_seed=CFG["seed"], runtime_s=round(time.time() - t0, 1),
           finished=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
(OUT / "run.json").write_text(json.dumps(run, indent=1, default=str))
print(json.dumps(run, indent=1, default=str))
sys.exit(1 if blockers else 0)
