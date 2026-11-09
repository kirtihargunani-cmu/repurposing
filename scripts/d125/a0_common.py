"""Shared harness for D125 section A0 (docs/SPEC_graph_metric_scoring.md; D125-A0 in docs/DECISIONS.md).

A scoring pass replaces the profiles of a set of drug x line pairs (columns of the 79c index) and runs the production
path on them: profiles (79c rules, each line's production generic response g_l) -> VIPER (79d, ARACNe v2 TF + CollecTRI,
changed columns only, spliced into copies of the production matrices) -> 79e (groups, subgroups) -> 79i (STRATA=1),
all unchanged, in a sandbox under data/processed/lincs/d125_a0/. Nothing at the top of results/lincs is written: the
sandbox's results/lincs holds links to the inputs only, and every pass output carries a RUN_TAG.
"""
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a0_config import (TAS_MIN, PROFILE_IDENTITY_TOL, VIPER_IDENTITY_TOL, R_POOLED_MIN_DRUGS, CI_LEVEL)  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PROC = ROOT / "data/processed/lincs"
RES = ROOT / "results/lincs"
OUT = RES / "d125/a0"                       # A0g outputs
WORK = PROC / "d125_a0"                     # sandbox, pass scores, timings (not tracked)
SBX = WORK / "sbx"
PASSES = WORK / "passes"
PY = os.environ.get("PY", str(Path.home() / "miniforge3/envs/celligner/bin/python"))
RSCRIPT = os.environ.get("RSCRIPT", str(Path.home() / "miniforge3/envs/repurposing/bin/Rscript"))
VIPER_SETS = {  # 79e / 79i consensus inputs: production dir name -> 79d environment
    "drug_viper_specific_v2_tf": {"VERSION": "v2", "REGSET": "tf"},
    "drug_viper_specific_collectri": {"REGULON": "collectri"},
}
SANDBOX_INPUTS = ["disease_drugs.tsv", "stratum_drugs_v2.tsv", "drug_transfer_confidence_v2.tsv",
                  "prolif_regulators_v2_tf.tsv", "prolif_regulators_collectri.tsv"]
GROUP_TABLE, STRATUM_TABLE = "drug_reversal_split_v2_consensus", "stratum_reversal_split_v2_consensus"


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


# ---------------------------------------------------------------------------------------------------------------------
# signatures and profiles
# ---------------------------------------------------------------------------------------------------------------------
def production_index():
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    idx["col"] = np.arange(len(idx))
    return idx


def active_signatures():
    """Active level-5 signatures (is_hiq == 1, TAS >= TAS_MIN), with their row in level5_bing_hiq.npy (= table index, as
    79c), the 79c column of their pair and the level-4 replicate count from siginfo (distil_ids)."""
    sig = pd.read_csv(PROC / "level5_bing_hiq_sigs.tsv.gz", sep="\t")
    assert (sig.is_hiq == 1).all()
    act = sig[sig.tas >= TAS_MIN].copy()
    act["row"] = act.index.values
    si = pd.read_csv(ROOT / "data/raw/lincs2020/siginfo_beta.txt", sep="\t", usecols=["sig_id", "distil_ids"],
                     low_memory=False)
    act = act.merge(si, on="sig_id", how="left")
    assert act.distil_ids.notna().all()
    act["n_l4"] = act.distil_ids.str.count(r"\|") + 1
    idx = production_index()
    act = act.merge(idx[["pert_id", "ModelID", "col"]], on=["pert_id", "ModelID"], how="left")
    assert act.col.notna().all(), "an active signature has no 79c profile"
    act["col"] = act.col.astype(int)
    n_sig = act.groupby("col").size().reindex(idx.col).values
    assert (n_sig == idx.n_sig.values).all(), "active signatures differ from the 79c index (n_sig)"
    return act


class Profiles:
    """79c gene mapping and each line's production g_l; rebuilds raw / specific profiles from signature rows."""

    def __init__(self):
        genes = pd.read_csv(PROC / "level5_bing_genes.tsv", sep="\t", dtype={"gene_id": str})
        idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t", dtype={"entrez_id": str})
        idm = idm[idm.primary_for_hgnc == True]  # noqa: E712  (as 79c)
        by_entrez = idm.dropna(subset=["entrez_id"]).drop_duplicates("entrez_id").set_index("entrez_id").gene_id
        v23 = genes.gene_id.map(by_entrez)
        v23 = v23.fillna(genes.ensembl_id.where(genes.ensembl_id.isin(idm.gene_id)))
        keep = v23.notna() & ~v23.duplicated()
        self.gpos = np.where(keep.values)[0]
        pg = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
        assert (pg.gene_id.values == v23[keep].values).all(), "79c gene mapping not reproduced"
        self.X = np.load(PROC / "level5_bing_hiq.npy", mmap_mode="r")
        self.idx = production_index()
        self.n, self.g = len(self.idx), len(self.gpos)
        self.raw = np.fromfile(PROC / "drug_line_profiles_bing.f32", dtype=np.float32).reshape(self.n, self.g)
        self.spec = np.fromfile(PROC / "drug_line_profiles_bing_specific.f32", dtype=np.float32).reshape(self.n, self.g)
        self.line_of = self.idx.ModelID.values
        self.gl = {}
        for m, ix in self.idx.groupby("ModelID").indices.items():
            v = self.raw[ix].astype(np.float64).mean(0)
            self.gl[m] = v / np.linalg.norm(v)

    def build(self, col, rows):
        """raw and specific profile of 79c column `col` from level-5 rows `rows` (mean, then g_l projected out)."""
        p = np.asarray(self.X[np.sort(np.asarray(rows))][:, self.gpos], dtype=np.float64).mean(0)
        g = self.gl[self.line_of[col]]
        return p, p - (p @ g) * g

    def identity_check(self, act, cols):
        """Rebuild the production profiles of `cols` from all their active signatures; must match 79c's files."""
        rows = act.groupby("col").row.apply(list)
        d_raw = d_spec = 0.0
        for c in cols:
            r, s = self.build(c, rows[c])
            d_raw = max(d_raw, float(np.abs(r - self.raw[c]).max()))
            d_spec = max(d_spec, float(np.abs(s - self.spec[c]).max()))
        ok = d_raw <= PROFILE_IDENTITY_TOL and d_spec <= PROFILE_IDENTITY_TOL
        log(f"profile identity on {len(cols)} columns: max |diff| raw {d_raw:.2e}, specific {d_spec:.2e} -> {'ok' if ok else 'FAILED'}")
        if not ok:
            raise SystemExit("profile identity check failed: the rebuilt profiles are not the production profiles")
        return dict(n_cols=len(cols), max_diff_raw=d_raw, max_diff_specific=d_spec)


# ---------------------------------------------------------------------------------------------------------------------
# sandbox
# ---------------------------------------------------------------------------------------------------------------------
def build_sandbox():
    """SBX: copies of 79d / 79e / 79i; config, R and results/* (not lincs) linked; data mirrored with links, except
    data/processed/lincs, a real folder of links into which pass folders (RUN_TAG names) are written."""
    if SBX.exists():
        shutil.rmtree(SBX)
    (SBX / "scripts").mkdir(parents=True)
    for s in ["79d_drug_viper.R", "79e_drug_reversal.py", "79i_split_reversal.py"]:
        shutil.copy2(ROOT / "scripts" / s, SBX / "scripts" / s)
    for d in ["config", "R"]:
        (SBX / d).symlink_to(ROOT / d)
    (SBX / "results/lincs").mkdir(parents=True)
    for d in (ROOT / "results").iterdir():
        if d.name != "lincs":
            (SBX / "results" / d.name).symlink_to(d)
    for f in SANDBOX_INPUTS:
        (SBX / "results/lincs" / f).symlink_to(RES / f)
    for d in (ROOT / "data").iterdir():
        if d.name != "processed":
            (SBX / "data" / d.name).parent.mkdir(parents=True, exist_ok=True)
            (SBX / "data" / d.name).symlink_to(d)
    for d in (ROOT / "data/processed").iterdir():
        if d.name != "lincs":
            (SBX / "data/processed").mkdir(parents=True, exist_ok=True)
            (SBX / "data/processed" / d.name).symlink_to(d)
    (SBX / "data/processed/lincs").mkdir(parents=True)
    for d in PROC.iterdir():
        if d.name != WORK.name:
            (SBX / "data/processed/lincs" / d.name).symlink_to(d)
    log(f"sandbox built: {SBX}")


def _env(threads, **kw):
    e = dict(os.environ)
    e.update({k: str(threads) for k in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]})
    e.update({k: str(v) for k, v in kw.items()})
    return e


def _run(cmd, env, logf):
    with open(logf, "w") as fh:
        r = subprocess.run(cmd, cwd=SBX, env=env, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode:
        raise SystemExit(f"failed ({r.returncode}): {' '.join(cmd)}; log {logf}")


def viper_subset(tag, cols, spec_rows, threads):
    """79d (both sets) on the specific profiles `spec_rows` (len(cols) x genes) of 79c columns `cols`.
    Returns {set dir: {net: (regulators, matrix len(cols) x regs)}}; the subset folders are removed."""
    pdir = SBX / f"sub{tag}"
    pdir.mkdir(exist_ok=True)
    idx = production_index().iloc[cols]
    idx[["pert_id", "ModelID", "n_sig"]].to_csv(pdir / "index.tsv", sep="\t", index=False)
    np.ascontiguousarray(spec_rows, dtype=np.float32).tofile(pdir / "specific.f32")
    out = {}
    for vset, kw in VIPER_SETS.items():
        rt = f"{tag}_sub"
        env = _env(1, PROFILE="specific", PROFILE_FILE=pdir / "specific.f32", INDEX_FILE=pdir / "index.tsv",
                   RUN_TAG=rt, N_CORES=threads, **kw)
        _run([RSCRIPT, "scripts/79d_drug_viper.R"], env, WORK / "logs" / f"{tag}_79d_{vset}.log")
        d = SBX / "data/processed/lincs" / f"{vset}{rt}"
        out[vset] = {}
        for f in sorted(d.glob("*.f32")):
            regs = (d / f"{f.stem}.regulators.txt").read_text().split()
            out[vset][f.stem] = (regs, np.fromfile(f, dtype=np.float32).reshape(len(cols), len(regs)))
        shutil.rmtree(d)
    shutil.rmtree(pdir)
    return out


def check_viper_identity(sub, cols):
    """VIPER of the rebuilt production profiles must equal the production VIPER columns."""
    worst, stale = 0.0, {}
    used = networks_used()
    for vset, nets in sub.items():
        prod = {p.stem for p in (PROC / vset).glob("*.f32")}
        if set(nets) - prod:
            raise SystemExit(f"{vset}: 79d produces networks production lacks ({sorted(set(nets) - prod)})")
        old = prod - set(nets)   # production-only files (e.g. ESCA / LGG of 10-25, before the subtype networks)
        if vset.endswith("_v2_tf") and old & used or vset.endswith("collectri") and old:
            raise SystemExit(f"{vset}: production networks that 79d no longer builds are used by a cohort ({sorted(old)})")
        stale[vset] = sorted(old)
        for net, (regs, M) in nets.items():
            preg = (PROC / vset / f"{net}.regulators.txt").read_text().split()
            if regs != preg:
                raise SystemExit(f"{vset}/{net}: regulator list differs from production")
            P = np.fromfile(PROC / vset / f"{net}.f32", dtype=np.float32).reshape(-1, len(preg))
            worst = max(worst, float(np.abs(P[cols] - M).max()))
    ok = worst <= VIPER_IDENTITY_TOL
    log(f"VIPER identity on {len(cols)} columns: max |diff| {worst:.2e} -> {'ok' if ok else 'FAILED'}; "
        f"unused production-only networks {stale}")
    if not ok:
        raise SystemExit("VIPER identity check failed: 79d no longer reproduces the production activity")
    return dict(n_cols=len(cols), max_diff_viper=worst, unused_production_only_networks=stale)


def networks_used():
    """ARACNe networks read by the 41 groups (79e / 79i: config/group_networks.tsv, else the cohort); subgroups use
    their parent's network."""
    g = pd.read_csv(RES / f"{GROUP_TABLE}.tsv", sep="\t", usecols=["group_id"]).group_id.unique()
    grp = pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t").set_index("group_id")
    nm = pd.read_csv(ROOT / "config/group_networks.tsv", sep="\t").set_index("group_id").net.to_dict()
    net = lambda c: "COADREAD" if c in ("TCGA-COAD", "TCGA-READ") else c.split("-", 1)[1]  # noqa: E731
    return {nm.get(x, net(grp.loc[x, "cohort"])) for x in g}


def stage_pass(tag, cols, raw_rows, spec_rows, sub):
    """Write the pass inputs into the sandbox: full profile files and spliced VIPER folders (RUN_TAG = tag)."""
    pre = SBX / f"prof{tag}"
    for kind, base, new in [("raw", PROF.raw, raw_rows), ("specific", PROF.spec, spec_rows)]:
        a = np.array(base, dtype=np.float32)
        a[cols] = new
        a.tofile(f"{pre}_{kind}.f32")
    for vset, nets in sub.items():
        d = SBX / "data/processed/lincs" / f"{vset}{tag}"
        d.mkdir()
        for net, (regs, M) in nets.items():
            preg = (PROC / vset / f"{net}.regulators.txt").read_text().split()
            assert regs == preg
            P = np.fromfile(PROC / vset / f"{net}.f32", dtype=np.float32).reshape(-1, len(preg))
            P[cols] = M
            P.tofile(d / f"{net}.f32")
            shutil.copy2(PROC / vset / f"{net}.regulators.txt", d / f"{net}.regulators.txt")
    return pre


def stage_baseline(tag):
    """Production profiles and VIPER, linked (no copy): the sandbox scoring of unchanged inputs."""
    pre = SBX / f"prof{tag}"
    Path(f"{pre}_raw.f32").symlink_to(PROC / "drug_line_profiles_bing.f32")
    Path(f"{pre}_specific.f32").symlink_to(PROC / "drug_line_profiles_bing_specific.f32")
    for vset in VIPER_SETS:
        (SBX / "data/processed/lincs" / f"{vset}{tag}").symlink_to(PROC / vset)
    return pre


def score_pass(tag, pre, threads):
    """79e (groups, then subgroups) and 79i (STRATA=1) in the sandbox; returns the compact spec_mech table and removes
    every file of the pass."""
    base = dict(VERSION="v2", REGULON="consensus", RUN_TAG=tag)
    t0 = time.time()
    _run([PY, "scripts/79e_drug_reversal.py"], _env(threads, **base), WORK / "logs" / f"{tag}_79e.log")
    _run([PY, "scripts/79e_drug_reversal.py"], _env(threads, STRATA=1, **base), WORK / "logs" / f"{tag}_79e_sub.log")
    _run([PY, "scripts/79i_split_reversal.py"], _env(threads, STRATA=1, PROFILE_PREFIX=pre, **base),
         WORK / "logs" / f"{tag}_79i.log")
    rl = SBX / "results/lincs"
    g = pd.read_csv(rl / f"{GROUP_TABLE}{tag}.tsv", sep="\t", low_memory=False,
                    usecols=["group_id", "pert_id", "source", "passes_gate", "spec_mech", "spec_mech_pct"])
    s = pd.read_csv(rl / f"{STRATUM_TABLE}{tag}.tsv", sep="\t", low_memory=False,
                    usecols=["stratum", "parent_group", "pert_id", "source", "passes_gate", "spec_mech", "spec_mech_pct"])
    out = pd.concat([g.rename(columns={"group_id": "cohort"}).assign(kind="group", parent_group=lambda x: x.cohort),
                     s.rename(columns={"stratum": "cohort"}).assign(kind="subgroup")], ignore_index=True)
    assert not out.duplicated(["kind", "cohort", "pert_id"]).any()
    secs = time.time() - t0
    cleanup(tag)
    return out, secs


def cleanup(tag):
    for p in list((SBX / "results/lincs").glob(f"*{tag}*")) + list(SBX.glob(f"prof{tag}_*")):
        p.unlink()
    for vset in VIPER_SETS:
        d = SBX / "data/processed/lincs" / f"{vset}{tag}"
        if d.is_symlink():
            d.unlink()
        elif d.exists():
            shutil.rmtree(d)


def save_pass(tag, df):
    PASSES.mkdir(parents=True, exist_ok=True)
    df.to_csv(PASSES / f"{tag}.tsv.gz", sep="\t", index=False)


def load_pass(tag):
    f = PASSES / f"{tag}.tsv.gz"
    return pd.read_csv(f, sep="\t") if f.exists() else None


def record_time(stage, n_cols, secs):
    f = WORK / "timings.tsv"
    row = pd.DataFrame([dict(stage=stage, n_cols=n_cols, seconds=round(secs, 1), when=time.strftime("%Y-%m-%d %H:%M"))])
    row.to_csv(f, sep="\t", index=False, mode="a", header=not f.exists())


def init_work():
    for d in [WORK, WORK / "logs", PASSES, OUT]:
        d.mkdir(parents=True, exist_ok=True)


PROF = None


def profiles():
    global PROF
    if PROF is None:
        PROF = Profiles()
    return PROF


# ---------------------------------------------------------------------------------------------------------------------
# which drugs draw on which pairs (79e line selection)
# ---------------------------------------------------------------------------------------------------------------------
def pairs_used(rows):
    """rows: kind, cohort, parent_group, pert_id, source. Returns rows + list of 79c columns each score draws on: direct
    -> the drug's lines in the cohort's tissue key with Celligner weight > 0 (all of them if every weight is 0);
    borrowed -> the drug's lines outside the key (79e / 79i activity())."""
    idx = production_index()
    by = idx.groupby("pert_id")
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
    models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", index_col=0)
    models = models[models.OncotreePrimaryDisease != "Non-Cancerous"]
    wg_ = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    ws_ = pd.read_csv(ROOT / "results/celligner/stratum_line_weights_long_v2.tsv", sep="\t").rename(columns={"stratum": "group_id"})
    used, inkey = [], {}
    for (kind, c), d in rows.groupby(["kind", "cohort"], sort=False):
        p = d.parent_group.iloc[0]
        k = keys.loc[p]
        km = models[models.OncotreeLineage.isin(str(k.lineages).split(";"))]
        if isinstance(k.diseases, str) and k.diseases:
            km = km[km.OncotreePrimaryDisease.isin(k.diseases.split(";"))]
        w = (wg_ if kind == "group" else ws_)
        wg = w[w.group_id == c].set_index("ModelID").weight
        inkey[(kind, c)] = idx.col[idx.ModelID.isin(km.index)].values
        for pid, src in zip(d.pert_id, d.source):
            q = by.get_group(pid)
            q = q[q.ModelID.isin(km.index)] if src == "direct" else q[~q.ModelID.isin(km.index)]
            ww = wg.reindex(q.ModelID).fillna(0).values if src == "direct" else np.zeros(len(q))
            used.append(q.col.values[ww > 0] if ww.sum() > 0 else q.col.values)
    out = rows.copy()
    out["cols_used"] = used
    return out, inkey


def cohort_inkey_lines():
    """(kind, cohort) -> set of the cohort's in-key DepMap lines with LINCS profiles (79e direct lines), for the 41
    groups and the 14 subgroups (tissue key of the parent group)."""
    idx = production_index()
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
    models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", index_col=0)
    models = models[models.OncotreePrimaryDisease != "Non-Cancerous"]
    g = pd.read_csv(RES / f"{GROUP_TABLE}.tsv", sep="\t", usecols=["group_id"]).drop_duplicates()
    s = pd.read_csv(RES / f"{STRATUM_TABLE}.tsv", sep="\t", usecols=["stratum", "parent_group"]).drop_duplicates()
    out = {}
    for kind, c, p in [("group", x, x) for x in g.group_id] + [("subgroup", a, b) for a, b in zip(s.stratum, s.parent_group)]:
        k = keys.loc[p]
        km = models[models.OncotreeLineage.isin(str(k.lineages).split(";"))]
        if isinstance(k.diseases, str) and k.diseases:
            km = km[km.OncotreePrimaryDisease.isin(k.diseases.split(";"))]
        out[(kind, c)] = set(idx.ModelID[idx.ModelID.isin(km.index)])
    return out


# ---------------------------------------------------------------------------------------------------------------------
# reliability
# ---------------------------------------------------------------------------------------------------------------------
def spearman_ci(r, n):
    """Fisher z CI with the Bonett-Wright variance for Spearman's rho."""
    from scipy.stats import norm
    if n <= 3 or not np.isfinite(r) or abs(r) >= 1:
        return np.nan, np.nan
    se = np.sqrt((1 + r ** 2 / 2) / (n - 3))
    z = norm.ppf(0.5 + CI_LEVEL / 2)
    return float(np.tanh(np.arctanh(r) - z * se)), float(np.tanh(np.arctanh(r) + z * se))


def r_by_cohort(a, b, changed_cols, used):
    """Per cohort: Spearman of spec_mech_pct between passes a and b over gate-passing drugs whose score draws on >= 1
    changed column (D125-A0 (3)). `used` = pairs_used() output for the cohort rows."""
    ch = set(int(c) for c in changed_cols)
    u = used.assign(touched=[bool(ch.intersection(map(int, cs))) for cs in used.cols_used])
    m = a.merge(b[["kind", "cohort", "pert_id", "spec_mech_pct"]], on=["kind", "cohort", "pert_id"], suffixes=("_a", "_b"))
    m = m.merge(u[["kind", "cohort", "pert_id", "touched"]], on=["kind", "cohort", "pert_id"], how="left")
    rows = []
    for (kind, c), d in m.groupby(["kind", "cohort"], sort=True):
        gate = d[d.passes_gate.astype(bool)]
        x = gate[gate.touched.fillna(False) & gate.spec_mech_pct_a.notna() & gate.spec_mech_pct_b.notna()]
        r = x.spec_mech_pct_a.corr(x.spec_mech_pct_b, method="spearman") if len(x) >= 3 else np.nan
        rows.append(dict(kind=kind, cohort=c, n_gate=len(gate), n_drugs=len(x), r_c=r))
    return pd.DataFrame(rows)


def finish_reliability(t, estimator):
    """CI (if absent), r_pooled and the median used for pooled cohorts."""
    t = t.copy()
    if "ci_lo" not in t:
        ci = [spearman_ci(r, n) for r, n in zip(t.r_c, t.n_drugs)]
        t["ci_lo"], t["ci_hi"] = [c[0] for c in ci], [c[1] for c in ci]
    t["r_pooled"] = t.n_drugs < R_POOLED_MIN_DRUGS
    med = float(t.loc[~t.r_pooled, "r_c"].median())
    t["r_used"] = np.where(t.r_pooled, med, t.r_c)
    t["estimator"] = estimator
    return t, med


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x)))
