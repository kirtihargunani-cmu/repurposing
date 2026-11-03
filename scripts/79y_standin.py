#!/usr/bin/env python
"""Closest-compound stand-ins between FDA drugs and LINCS compounds (user 2026-11-02; rules fixed in config/drug_gate.yaml:
standin, before looking at results).

Many FDA drugs were never profiled in a cohort's own lines (borrowed) or only in lines unlike it (unmatched_lines /
low_transfer), while some other LINCS compound with the same target and a similar effect was. This step pairs them.
Similarity of two compounds = mean over their shared profiled lines of the Pearson correlation of their line-specific
regulator activity (CollecTRI, one context-free regulon for all lines: data/processed/lincs/drug_viper_specific_collectri;
line average removed in 79c/79d). Pair accepted if >= min_shared_lines shared lines, similarity >= min_similarity and
>= 1 shared annotated target (LINCS / Repurposing Hub, ADME genes of the target gate excluded); with
shared_target_evidence = chembl (user 2026-11-02) both compounds must also bind that target in ChEMBL 37 by the receptor
gate's rule (scripts/chembl_lib.py; cache data/raw/chembl/ChEMBL_37/standin/). LINCS ids of the same
FDA ingredient or the same name never pair.
Per cohort (41 diseases, 14 subgroups):
  matched compound = passes the target gate, tested in the cohort's matched lines (direct), no unmatched_lines flag
  - FDA drug needing a stand-in (borrowed or weak evidence): its accepted matched compounds, best first (standin_rank)
  - matched compound: its accepted FDA drugs in the cohort's FDA list, best first (closest_rank)
Validation (diseases decide adoption; the 14 subgroups are reported separately on spec_sib_mech, similarity outside
their parent's key lines; user 2026-11-02): FDA drugs that are themselves matched in a cohort get
their best accepted matched stand-in, with similarity measured only on lines outside the cohort's tissue key (as in
use: a borrowed drug was never profiled there; otherwise the score's own profiles would enter the similarity); Spearman of the stand-in's spec_mech percentile with the drug's own, with a
bootstrap 95% CI, against random shared-target pairs (no similarity rule; 200 draws). Adopted if Spearman >= min and
its CI lies above the random pairs' median. Also by number of shared lines and similarity band.
Reports show a stand-in next to the drug's own score; it never replaces it.
Inputs: drug_line_profiles_index.tsv, drug_viper_specific_collectri/ALL.f32, drug_annotations.tsv,
        drug_reversal_split_v2_consensus.tsv, stratum_reversal_split_v2_consensus.tsv, fda_matches.tsv (scripts/fda_rows.py)
Outputs: results/lincs/standin_pairs_v2.tsv (accepted pairs with both sides' scores), standin_validation_v2.tsv,
         standin_run.json
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
RES, PROC = ROOT / "results/lincs", ROOT / "data/processed/lincs"
sys.path.insert(0, str(ROOT / "scripts"))
from chembl_lib import RELEASE, Binding  # noqa: E402
from fda_rows import fda_rows, truthy  # noqa: E402

CFG = yaml.safe_load(open(ROOT / "config/drug_gate.yaml"))
SC = CFG["standin"]
MIN_LINES, MIN_R, NEED_T = SC["min_shared_lines"], SC["min_similarity"], SC["require_shared_target"]
SEED = 0


def adme_genes():
    hg = pd.read_csv(sorted((ROOT / "data/raw/hgnc").glob("hgnc_complete_set_*.txt"))[-1], sep="\t", dtype=str,
                     usecols=["symbol", "gene_group_id"]).dropna()
    ids = {str(i) for i in CFG["adme"]["hgnc_group_ids"]}
    return {r.symbol for r in hg.itertuples() if set(r.gene_group_id.split("|")) & ids} | set(CFG["adme"]["genes"])


def similarity():
    """(perts, S, N): S[i, j] = mean per-line Pearson r of compounds i, j over shared lines, N = shared line count"""
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    d = PROC / "drug_viper_specific_collectri"
    regs = (d / "ALL.regulators.txt").read_text().split()
    A = np.fromfile(d / "ALL.f32", dtype=np.float32).reshape(len(idx), len(regs)).astype(np.float64)
    A = np.nan_to_num(A)
    A = A - A.mean(1, keepdims=True)
    A /= np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-12)
    perts = sorted(idx.pert_id.unique())
    pi = {p: i for i, p in enumerate(perts)}
    S = np.zeros((len(perts), len(perts)), dtype=np.float32)
    N = np.zeros((len(perts), len(perts)), dtype=np.int16)
    for _, x in idx.reset_index().groupby("ModelID"):
        rows, ii = x["index"].values, np.array([pi[p] for p in x.pert_id])
        C = A[rows] @ A[rows].T
        S[np.ix_(ii, ii)] += C.astype(np.float32)
        N[np.ix_(ii, ii)] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        S = np.where(N > 0, S / np.maximum(N, 1), np.nan).astype(np.float32)
    np.fill_diagonal(N, 0)
    vec = {(p, m): A[k] for k, (p, m) in enumerate(zip(idx.pert_id, idx.ModelID))}
    lines = idx.groupby("pert_id").ModelID.apply(set).to_dict()
    return perts, pi, S, N, vec, lines


def main():
    t0 = time.time()
    perts, pi, S, N, vec, lines = similarity()
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
    models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", index_col=0)
    models = models[models.OncotreePrimaryDisease != "Non-Cancerous"]

    def key_lines(g):
        k = keys.loc[g]
        km = models[models.OncotreeLineage.isin(str(k.lineages).split(";"))]
        if isinstance(k.diseases, str) and k.diseases:
            km = km[km.OncotreePrimaryDisease.isin(k.diseases.split(";"))]
        return set(km.index)

    def sim_excl(a, b, excl):
        """similarity over shared lines outside excl (validation: the drug's own matched lines are hidden)"""
        sh = (lines.get(a, set()) & lines.get(b, set())) - excl
        if not sh:
            return np.nan, 0
        return float(np.mean([vec[(a, m)] @ vec[(b, m)] for m in sh])), len(sh)
    adme = adme_genes()
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t").drop_duplicates("pert_id").set_index("pert_id")
    tg = {p: {t for c in ("lincs_target", "hub_target") for t in str(r[c]).split("|") if t not in ("nan", "") and t not in adme}
          for p, r in ann.iterrows()}
    nm = ann.cmap_name.to_dict()
    # ---- shared targets: annotated on both; with shared_target_evidence = chembl, ChEMBL-backed binding by both too
    G, Sg = fda_rows()
    uni = set(G.pert_id) | set(Sg.pert_id)
    for f in ("drug_reversal_split_v2_consensus.tsv", "stratum_reversal_split_v2_consensus.tsv"):
        x = pd.read_csv(RES / f, sep="\t", low_memory=False, usecols=["pert_id", "source", "unmatched_lines", "passes_gate"])
        uni |= set(x[(x.source == "direct") & ~x.unmatched_lines.map(truthy) & x.passes_gate.map(truthy)].pert_id)
    uni = {p for p in uni if tg.get(p)}
    EVID = SC.get("shared_target_evidence", "annotation")
    B = Binding(uni, ROOT / f"data/raw/chembl/{RELEASE}/standin") if EVID == "chembl" else None
    _sh = {}

    def shared(a, b):
        k = (a, b) if a < b else (b, a)
        if k not in _sh:
            s_ = tg.get(a, set()) & tg.get(b, set())
            _sh[k] = {g for g in s_ if B.binds(a, g) and B.binds(b, g)} if B else s_
        return _sh[k]
    fm = pd.read_csv(RES / "fda_matches.tsv", sep="\t").drop_duplicates("pert_id").set_index("pert_id").fda_ingredient.to_dict()
    fda_by = {"disease": G, "subgroup": Sg}
    cols = ["group_id", "pert_id", "cmap_name", "source", "unmatched_lines", "low_transfer", "line_match", "rev_pct",
            "spec_pct", "spec_mech", "spec_mech_pct"]
    full = {"disease": pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False),
            "subgroup": pd.read_csv(RES / "stratum_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
            .rename(columns={"stratum": "group_id"})}
    for k, df in full.items():
        df = df[df.passes_gate.map(truthy)].copy()
        for f in ("unmatched_lines", "low_transfer"):
            df[f] = df[f].map(truthy)
        if k == "subgroup":
            cols_k = cols + ["spec_sib_mech_pct", "parent_group"]
        else:
            df["spec_sib_mech_pct"], df["parent_group"] = np.nan, ""
            cols_k = cols + ["spec_sib_mech_pct", "parent_group"]
        full[k] = df[cols_k]

    def accept(a, b):
        """(ok, r, n, shared targets) for compounds a, b"""
        if a not in pi or b not in pi or a == b:
            return False, np.nan, 0, ""
        if fm.get(a) and fm.get(a) == fm.get(b) or nm.get(a) == nm.get(b):
            return False, np.nan, 0, ""
        i, j = pi[a], pi[b]
        n, r = int(N[i, j]), float(S[i, j])
        sh = shared(a, b)
        ok = n >= MIN_LINES and r >= MIN_R and (bool(sh) or not NEED_T)
        return ok, r, n, "|".join(sorted(sh))

    pairs = []
    for lvl, df in full.items():
        F = fda_by[lvl]
        for g, x in df.groupby("group_id"):
            matched = x[(x.source == "direct") & ~x.unmatched_lines]
            fx = x[x.pert_id.isin(F[F.group_id == g].pert_id)]
            for fr in fx.itertuples(index=False):
                fi = pi.get(fr.pert_id)
                if fi is None:
                    continue
                for cr in matched.itertuples(index=False):
                    if cr.pert_id == fr.pert_id:
                        continue
                    ok, r, n, sh = accept(fr.pert_id, cr.pert_id)
                    if ok:
                        pairs.append(dict(level=lvl, group_id=g, parent_group=fr.parent_group,
                                          fda_pert=fr.pert_id, fda_drug=F.loc[F.pert_id == fr.pert_id, "drug"].iloc[0],
                                          fda_source=fr.source, fda_unmatched=fr.unmatched_lines, fda_low_transfer=fr.low_transfer,
                                          fda_rev_pct=fr.rev_pct, fda_spec_pct=fr.spec_pct, fda_spec_mech_pct=fr.spec_mech_pct,
                                          fda_spec_sib_mech_pct=fr.spec_sib_mech_pct,
                                          fda_needs_standin=(fr.source == "borrowed") or fr.unmatched_lines or fr.low_transfer,
                                          cmp_pert=cr.pert_id, cmp_name=cr.cmap_name, cmp_is_fda=cr.pert_id in fm,
                                          cmp_line_match=cr.line_match, cmp_rev_pct=cr.rev_pct, cmp_spec_pct=cr.spec_pct,
                                          cmp_spec_mech_pct=cr.spec_mech_pct, cmp_spec_sib_mech_pct=cr.spec_sib_mech_pct,
                                          similarity=round(r, 3), n_shared_lines=n, shared_targets=sh))
    P = pd.DataFrame(pairs)
    P["standin_rank"] = P.groupby(["level", "group_id", "fda_pert"]).similarity.rank(ascending=False, method="first").astype(int)
    P["closest_rank"] = P.groupby(["level", "group_id", "cmp_pert"]).similarity.rank(ascending=False, method="first").astype(int)
    P.sort_values(["level", "group_id", "fda_drug", "standin_rank"]).to_csv(RES / "standin_pairs_v2.tsv", sep="\t", index=False)

    # ---------------- validation: FDA drugs matched themselves; diseases (decides adoption), then subgroups (reported)
    rng = np.random.default_rng(SEED)

    def boot(a, b, n_boot):
        a, b = np.asarray(a), np.asarray(b)
        st = [spearmanr(a[k], b[k])[0] for k in (rng.integers(0, len(a), len(a)) for _ in range(n_boot))]
        return np.nanpercentile(st, [2.5, 97.5])

    def validate(level, D, F, col, pre):
        """D: cohort table, F: FDA rows, col: lead percentile column; similarity outside the (parent's) key lines"""
        V, rand_pool = [], []
        for g, x in D.groupby("group_id"):
            matched = x[(x.source == "direct") & ~x.unmatched_lines]
            mset = matched.set_index("pert_id")
            par = x.parent_group.iloc[0] if level == "subgroup" else g
            excl = key_lines(par)
            fdam = matched[matched.pert_id.isin(F[F.group_id == g].pert_id)]
            for fr in fdam.itertuples(index=False):
                best, cand_t = None, []
                for cr in matched.itertuples(index=False):
                    if cr.pert_id == fr.pert_id or (fm.get(cr.pert_id) and fm.get(cr.pert_id) == fm.get(fr.pert_id)) or cr.cmap_name == fr.cmap_name:
                        continue
                    sh = shared(fr.pert_id, cr.pert_id)
                    if sh:
                        cand_t.append(cr.pert_id)
                    r, n = sim_excl(fr.pert_id, cr.pert_id, excl)
                    ok = n >= MIN_LINES and r >= MIN_R and (bool(sh) or not NEED_T)
                    if ok and (best is None or r > best[1]):
                        best = (cr.pert_id, r, n)
                own = getattr(fr, col)
                if best:
                    V.append(dict(level=level, group_id=g, fda_pert=fr.pert_id, own=own, cmp=best[0],
                                  standin=mset.loc[best[0], col], similarity=best[1], n_shared_lines=best[2]))
                if cand_t:
                    rand_pool.append((own, [mset.loc[c, col] for c in cand_t]))
        V = pd.DataFrame(V)
        out = []
        rho = spearmanr(V.own, V.standin)[0] if len(V) > 2 else np.nan
        lo, hi = boot(V.own, V.standin, SC["validation"]["n_boot"]) if len(V) > 2 else (np.nan, np.nan)
        rr = [spearmanr([o for o, _ in rand_pool], [c[rng.integers(len(c))] for _, c in rand_pool])[0] for _ in range(200)]
        rand_med = float(np.median(rr))
        out.append(dict(set=pre + "accepted stand-ins", n=len(V), spearman=round(rho, 3), ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                        mae_pct=round(float((V.own - V.standin).abs().mean()), 1) if len(V) else np.nan))
        out.append(dict(set=pre + "random shared-target pairs (median of 200)", n=len(rand_pool), spearman=round(rand_med, 3),
                        ci_lo=round(float(np.percentile(rr, 2.5)), 3), ci_hi=round(float(np.percentile(rr, 97.5)), 3), mae_pct=np.nan))
        for lab, m in [("shared lines = 2", V.n_shared_lines == 2), ("shared lines >= 3", V.n_shared_lines >= 3),
                       ("similarity 0.4-0.6", V.similarity < 0.6), ("similarity >= 0.6", V.similarity >= 0.6)]:
            x = V[m]
            out.append(dict(set=pre + lab, n=len(x), spearman=round(spearmanr(x.own, x.standin)[0], 3) if len(x) > 2 else np.nan,
                            ci_lo=np.nan, ci_hi=np.nan, mae_pct=round(float((x.own - x.standin).abs().mean()), 1) if len(x) else np.nan))
        return V, out, bool(rho >= SC["validation"]["min_spearman"] and lo > rand_med)

    Vd, Od, adopted = validate("disease", full["disease"], G, "spec_mech_pct", "")
    Vs, Os, adopted_sub = validate("subgroup", full["subgroup"], Sg, "spec_sib_mech_pct", "subgroup: ")
    O = pd.DataFrame(Od + Os)
    O.to_csv(RES / "standin_validation_v2.tsv", sep="\t", index=False)
    pd.concat([Vd, Vs], ignore_index=True).to_csv(RES / "standin_validation_pairs_v2.tsv", sep="\t", index=False)

    need = P[P.fda_needs_standin & (P.standin_rank == 1)]
    allneed = sum(((fda_by[l].source == "borrowed") | fda_by[l].unmatched_lines.map(truthy) | fda_by[l].low_transfer.map(truthy)).sum()
                  for l in fda_by)
    run = dict(date=time.strftime("%Y-%m-%d"), config=SC, chembl=B.n if B else None, adopted=adopted, adopted_subgroups_same_rule=adopted_sub,
               validation=O.to_dict("records"),
               pairs=len(P), fda_rows_needing=int(allneed),
               fda_rows_with_standin=int(len(need)),
               fda_rows_with_standin_by_level=need.level.value_counts().to_dict(),
               cohorts_with_standins=int(need.group_id.nunique()),
               matched_compounds_with_close_fda=int(P[P.closest_rank == 1].drop_duplicates(["level", "group_id", "cmp_pert"]).shape[0]),
               seconds=round(time.time() - t0))
    json.dump(run, open(RES / "standin_run.json", "w"), indent=1, default=str)
    print(O.to_string(index=False))
    print(json.dumps({k: v for k, v in run.items() if k not in ("validation", "config")}, indent=1))


if __name__ == "__main__":
    main()
