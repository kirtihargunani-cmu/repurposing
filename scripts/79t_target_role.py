#!/usr/bin/env python
"""D120: target-role score. Is the drug's main target needed, altered or turned up in this cancer? (runs in `celligner`)

Rule: docs/DECISIONS.md D120 (fixed before running). Reported beside the ranking; never ranks, gates or tiers.
Gene components per group c x gene g (within-group percentiles over genes, ties averaged, 0-100):
  D dependency      79n gene_dep (weighted mean of -z Chronos over in-key lines), percentile over the group's genes
  A alteration      f = share of the group's TCGA tumors with an MC3 PASS missense / in-frame indel at a recurrent
                    position (>= 3 MC3 tumors pan-cancer) or GISTIC +2; f < 5% -> 0; A = f - median f over groups
  O over-activity   regulator in the group's MR tables (ARACNe TF v2, CollecTRI v2): mean NES, percentile over the
                    group's regulators; else logFC tumor vs GTEx (adjacent if skipped), percentile over genes
Direction (ChEMBL action_type): lowering -> as is; raising -> D, O = 100 - pct, A not scored; other -> lowering, flagged.
Gene score T = mean of available oriented components; drug row = main target with the highest T; drug level = mean
over gate-passing compounds of the FDA ingredient (79l drug_table rows).
Evaluation and decomposition (a)-(i): see D120. (j) drug-level vs cohort-specific part and (k) number of main targets were
added after seeing (a)-(i) (reported, not part of the rule). Outputs results/lincs/target_role/.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
OUT = RES / "target_role"
sys.path.insert(0, str(ROOT / "scripts"))
from fda_rows import truthy  # noqa: E402

MIN_F, MIN_REC, SEED, N_PERM, N_BOOT = 0.05, 3, 0, 1000, 1000
LOWER = {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR", "NEGATIVE MODULATOR", "ALLOSTERIC ANTAGONIST",
         "INVERSE AGONIST", "DEGRADER", "DISRUPTING AGENT", "STABILISER"}
RAISE = {"AGONIST", "PARTIAL AGONIST", "POSITIVE ALLOSTERIC MODULATOR", "POSITIVE MODULATOR", "ACTIVATOR", "OPENER", "RELEASING AGENT"}
PROT_SRC = {"chembl_mechanism", "confirmed_annotation", "annotation"}
COMP = ["D", "A", "O"]
POS_CTL = [(g, d) for g in ["BRCA-LumA", "BRCA-LumB"] for d in ["FULVESTRANT", "TAMOXIFEN CITRATE", "RALOXIFENE HYDROCHLORIDE", "PALBOCICLIB"]]
NEG_CTL = [("BRCA-Basal", d) for d in ["FULVESTRANT", "TAMOXIFEN CITRATE", "RALOXIFENE HYDROCHLORIDE"]]

_s = importlib.util.spec_from_file_location("ku", ROOT / "scripts/79l_known_use.py"); KU = importlib.util.module_from_spec(_s); _s.loader.exec_module(KU)


def pct(s):
    return s.rank(pct=True, method="average") * 100


# ------------------------------------------------------------------ gene components
def comp_dependency(groups):
    G = pd.read_csv(RES / "dependency/gene_dep.tsv.gz", sep="\t", index_col=0)
    return {c: G[c].dropna() for c in groups if c in G}, {c: pct(G[c].dropna()) for c in groups if c in G}


def alteration_pairs():
    """(gene, sample) pairs with an activating alteration in TCGA (MC3 recurrent missense / in-frame, GISTIC +2), and
    the set of profiled samples (in MC3 or GISTIC). Shared by groups (79t) and subgroups (79t3)."""
    maf = pd.read_csv(ROOT / "data/raw/tcga_pancan/mc3.v0.2.8.PUBLIC.maf.gz", sep="\t", low_memory=False,
                      usecols=["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", "Protein_position", "FILTER"])
    maf["sample"] = maf.Tumor_Sample_Barcode.str[:15]
    prof_mut = set(maf["sample"])
    m = maf[(maf.FILTER == "PASS") & maf.Variant_Classification.isin(["Missense_Mutation", "In_Frame_Del", "In_Frame_Ins"])].copy()
    m["pos"] = m.Protein_position.astype(str).str.split("/").str[0]
    m = m[m.pos.str.match(r"^\d")]
    rec = m.drop_duplicates(["Hugo_Symbol", "pos", "sample"]).groupby(["Hugo_Symbol", "pos"])["sample"].nunique()
    rec = rec[rec >= MIN_REC].reset_index()[["Hugo_Symbol", "pos"]]
    mut = m.merge(rec, on=["Hugo_Symbol", "pos"])[["Hugo_Symbol", "sample"]].rename(columns={"Hugo_Symbol": "gene"})
    amp, prof_cn = [], set()
    for ch in pd.read_csv(ROOT / "data/raw/tcga_pancan/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz", sep="\t",
                          index_col=0, chunksize=2000):
        prof_cn |= set(ch.columns)
        v = ch.values
        r, k = np.nonzero(v == 2)
        amp.append(pd.DataFrame({"gene": ch.index.values[r], "sample": np.array(ch.columns)[k]}))
    amp = pd.concat(amp)
    alt = pd.concat([mut, amp]).drop_duplicates()
    return alt, prof_mut | prof_cn, dict(recurrent_positions=int(len(rec)), mutation_pairs=int(len(mut)), amp_pairs=int(len(amp)))


def alteration_shares(alt, prof, members):
    """members {unit: set(samples)} -> DataFrame genes x units of alteration shares (units with >= 10 profiled tumors)."""
    s2u = [(smp, u) for u, ss in members.items() for smp in ss if smp in prof]
    S = pd.DataFrame(s2u, columns=["sample", "unit"])
    den = S.groupby("unit")["sample"].nunique()
    den = den[den >= 10]
    cnt = alt.merge(S, on="sample").groupby(["unit", "gene"])["sample"].nunique().unstack(0).fillna(0)
    return cnt.reindex(columns=den.index).fillna(0) / den, den


def comp_alteration(groups):
    sg = pd.read_csv(ROOT / "results/step1/sample_groups_toil.tsv", sep="\t")
    sg = sg[(sg.role == "tumor") & sg.group_id.isin(groups)]
    alt, prof, info = alteration_pairs()
    F, den = alteration_shares(alt, prof, {g: set(d["sample"]) for g, d in sg.groupby("group_id")})
    genes_all = set(pd.read_csv(RES / "dependency/gene_dep.tsv.gz", sep="\t", usecols=[0]).iloc[:, 0]) | set(F.index)
    F = F.reindex(sorted(genes_all)).fillna(0)
    F0 = F.where(F >= MIN_F, 0)
    A = F0.sub(F0.median(1), axis=0)
    info["tumors_profiled"] = den.to_dict()
    return {c: F[c] for c in F}, {c: pct(A[c]) for c in A}, info


def comp_overactivity(groups, mr_dirs=("aracne_project_v2_tf", "collectri_v2"), summary="results/step4_v2/signature_summary.tsv",
                      proc="data/processed/step4_v2"):
    nes = []
    for d in mr_dirs:
        x = pd.read_csv(ROOT / f"results/step6/{d}/group_mr.tsv", sep="\t", usecols=["group_id", "regulator", "nes"])
        nes.append(x)
    N = pd.concat(nes).groupby(["group_id", "regulator"]).nes.mean()
    raw, P, src = {}, {}, {}
    ss = pd.read_csv(ROOT / summary, sep="\t")
    for c in groups:
        p, r, s = pd.Series(dtype=float), pd.Series(dtype=float), pd.Series(dtype=object)
        ok = ss[(ss.group_id == c) & (ss.status == "ok")].contrast.tolist()
        con = "gtex" if "gtex" in ok else ("adjacent" if "adjacent" in ok else None)
        if con:
            e = pd.read_csv(ROOT / f"{proc}/{c}__{con}.tsv.gz", sep="\t", usecols=["symbol", "logFC"])
            e = e.dropna().drop_duplicates("symbol").set_index("symbol").logFC
            p, r, s = pct(e), e, pd.Series(f"logFC_vs_{con}", index=e.index)
        if c in N.index.get_level_values(0):
            n = N.loc[c]
            p = pd.concat([p.drop(n.index, errors="ignore"), pct(n)])
            r = pd.concat([r.drop(n.index, errors="ignore"), n])
            s = pd.concat([s.drop(n.index, errors="ignore"), pd.Series("NES", index=n.index)])
        raw[c], P[c], src[c] = r, p, s
    return raw, P, src


# ------------------------------------------------------------------ drug targets and action
def drug_actions(names):
    """cmap_name -> {gene: action class}, from the ChEMBL mechanisms of the D105b lookup."""
    spec = importlib.util.spec_from_file_location("i2", ROOT / "scripts/79i2_target_presence.py")
    i2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(i2)
    w = i2.w
    w.CACHE = w.RAW / "chembl" / w.RELEASE / "d105_all"
    d2m, _, _, _, _ = w.chembl(sorted(names))
    mech = json.load(open(w.CACHE / "mechanisms.json")) + json.load(open(w.CACHE / "mechanisms_by_parent.json"))
    T = {t["target_chembl_id"]: t for t in json.load(open(w.CACHE / "targets.json")) + json.load(open(w.CACHE / "targets_by_parent.json"))}
    out = {}
    for d in names:
        ms = d2m.get(d, set())
        acts = {}
        for m in mech:
            if m["direct_interaction"] != 1 or not m["target_chembl_id"] or not ({m["molecule_chembl_id"], m["parent_molecule_chembl_id"]} & ms):
                continue
            t = T.get(m["target_chembl_id"])
            if not t:
                continue
            try:
                genes = i2.target_genes(t)
            except SystemExit:
                continue
            for g in genes:
                acts.setdefault(g, set()).add(m["action_type"])
        out[d] = acts
    return out


def orient(acts):
    """set of action types -> (orientation, flag)."""
    lo, hi = bool(acts & LOWER), bool(acts & RAISE)
    if lo and hi:
        return "lower", "mixed"
    if hi:
        return "raise", ""
    if lo:
        return "lower", ""
    return "lower", "direction_assumed"


def gene_score(c, g, orient_, Pd, Pa, Po, signed=True, drop=None):
    v = {}
    raise_ = signed and orient_ == "raise"
    if c in Pd and g in Pd[c].index:
        v["D"] = 100 - Pd[c][g] if raise_ else Pd[c][g]
    if c in Pa and g in Pa[c].index and not raise_:
        v["A"] = Pa[c][g]
    if c in Po and g in Po[c].index:
        v["O"] = 100 - Po[c][g] if raise_ else Po[c][g]
    if drop:
        v.pop(drop, None)
    return v


def main():
    rng = np.random.default_rng(SEED)
    OUT.mkdir(parents=True, exist_ok=True)
    sp = pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
    groups = sorted(sp.group_id.unique())
    Rd, Pd = comp_dependency(groups)
    Ra, Pa, ainfo = comp_alteration(groups)
    Ro, Po, Osrc = comp_overactivity(groups)
    print("components:", {k: len(v) for k, v in (("D", Pd), ("A", Pa), ("O", Po))}, "groups", flush=True)
    # gene components (long) for the record
    rows = []
    for c in groups:
        idx = set(Pd.get(c, pd.Series(dtype=float)).index) | set(Po.get(c, pd.Series(dtype=float)).index)
        if c in Ra:
            idx |= set(Ra[c][Ra[c] > 0].index)
        for g in idx:
            rows.append((c, g, Rd.get(c, {}).get(g, np.nan), Pd.get(c, {}).get(g, np.nan), Ra.get(c, {}).get(g, np.nan),
                         Pa.get(c, {}).get(g, np.nan), Ro.get(c, {}).get(g, np.nan), Po.get(c, {}).get(g, np.nan), Osrc.get(c, {}).get(g, "")))
    GC = pd.DataFrame(rows, columns=["group_id", "gene", "dep_z", "D", "alt_share", "A", "o_raw", "O", "o_source"])
    GC.to_csv(OUT / "gene_components.tsv.gz", sep="\t", index=False, float_format="%.4g")

    fm = KU.fda_map()
    x = sp[sp.pert_id.isin(fm.index) & sp.passes_gate.map(truthy)].copy()
    x["fda_ingredient"] = x.pert_id.map(fm)
    x["main_targets"] = x.evidence_targets.fillna("").map(lambda s: [t.split(":")[0] for t in s.split(";") if t])
    x.loc[~x.target_source.isin(PROT_SRC), "main_targets"] = x.loc[~x.target_source.isin(PROT_SRC), "main_targets"].map(lambda _: [])
    ACT = drug_actions(set(x.cmap_name))

    def best(c, cmap, tg, signed=True, drop=None):
        bestv = None
        for g in tg:
            o, flag = orient(ACT.get(cmap, {}).get(g, set()))
            v = gene_score(c, g, o, Pd, Pa, Po, signed, drop)
            if not v:
                continue
            t = float(np.mean(list(v.values())))
            if bestv is None or t > bestv[0]:
                bestv = (t, g, o, flag, v)
        return bestv

    VARS = {"T": dict(), "T_noD": dict(drop="D"), "T_noA": dict(drop="A"), "T_noO": dict(drop="O"), "T_unsigned": dict(signed=False)}
    rec = []
    for r in x[["group_id", "pert_id", "cmap_name", "fda_ingredient", "main_targets", "target_source"]].itertuples(index=False):
        o = dict(group_id=r.group_id, pert_id=r.pert_id, cmap_name=r.cmap_name, fda_ingredient=r.fda_ingredient,
                 target_source=r.target_source, n_main_targets=len(r.main_targets))
        for name, kw in VARS.items():
            b = best(r.group_id, r.cmap_name, r.main_targets, **kw)
            o[name] = b[0] if b else np.nan
            if name == "T" and b:
                t, g, ori, flag, v = b
                o.update(gene=g, action=ori, action_flag=flag, D=v.get("D", np.nan), A=v.get("A", np.nan), O=v.get("O", np.nan),
                         n_components=len(v), dep_z=Rd.get(r.group_id, {}).get(g, np.nan), alt_share=Ra.get(r.group_id, {}).get(g, np.nan),
                         o_raw=Ro.get(r.group_id, {}).get(g, np.nan), o_source=Osrc.get(r.group_id, {}).get(g, ""),
                         actions="|".join(sorted(a for a in ACT.get(r.cmap_name, {}).get(g, set()) if a)))
        # each component alone: best gene on that component only
        for cpt in COMP:
            vals = []
            for g in r.main_targets:
                ori, _ = orient(ACT.get(r.cmap_name, {}).get(g, set()))
                v = gene_score(r.group_id, g, ori, Pd, Pa, Po)
                if cpt in v:
                    vals.append(v[cpt])
            o[f"{cpt}_alone"] = max(vals) if vals else np.nan
        rec.append(o)
    C = pd.DataFrame(rec)
    C.to_csv(OUT / "compound_target_role.tsv", sep="\t", index=False, float_format="%.3f")
    num = ["T", "T_noD", "T_noA", "T_noO", "T_unsigned", "D", "A", "O", "D_alone", "A_alone", "O_alone", "dep_z", "alt_share", "o_raw", "n_components", "n_main_targets"]
    txt = ["gene", "action", "action_flag", "o_source", "actions", "target_source"]
    agg = {c: (c, "mean") for c in num}
    agg.update({c: (c, lambda s: "|".join(sorted({str(v) for v in s if isinstance(v, str) and v}))) for c in txt})
    DR = C.groupby(["group_id", "fda_ingredient"]).agg(**agg).reset_index().rename(columns={"group_id": "cohort"})
    DR["target_role_pct"] = DR.groupby("cohort")["T"].transform(lambda s: s.rank(pct=True, method="max") * 100)
    lab = pd.read_csv(RES / "known_use/labels.tsv", sep="\t")
    L = lab.merge(DR, on=["cohort", "fda_ingredient"], how="left")
    L["essential_target"] = L.gene.map(lambda g: False)  # filled below
    ess = set(pd.read_csv(ROOT / "data/raw/depmap_24q4/CRISPRInferredCommonEssentials.csv").iloc[:, 0].str.split(" ").str[0])
    L["essential_target"] = L.gene.fillna("").map(lambda s: any(g in ess for g in s.split("|") if g))
    L["T_plus_spec"] = (L["T"].groupby(L.cohort).rank(pct=True) * 100 + L.spec_mech_pct_pregate.groupby(L.cohort).rank(pct=True) * 100) / 2
    keep = ["cohort", "fda_ingredient", "drug", "anticancer", "lab_hub", "lab_chembl", "excl_hub", "excl_chembl", "eligible",
            "spec_mech_pct", "spec_mech_pct_pregate", "target_role_pct"] + [c for c in num + txt] + ["essential_target", "T_plus_spec"]
    L[keep].to_csv(OUT / "drug_target_role.tsv", sep="\t", index=False, float_format="%.3f")

    # ---------------------------------------------------------------- evaluation (reported only)
    ev = []
    def run(src, popname, X, score, label, perm=False):
        a = KU.pooled_auc(X, f"lab_{src}", score, rng if perm else None, N_PERM if perm else 0, N_BOOT if perm else 0)
        ev.append(dict(source=src, rows=popname, score=label, **a))
    for src in ("hub", "chembl"):
        base = L[L.anticancer & ~L[f"excl_{src}"]]
        P = base[base["T"].notna()]
        run(src, "primary: anticancer, gate-passing, T present", P, "T", "T (target role, primary)", True)
        run(src, "primary: anticancer, gate-passing, T present", P, "spec_mech_pct_pregate", "spec_mech_pct_pregate (same rows)", True)
        E = P[P.eligible.map(truthy)]
        run(src, "D118 primary rows (eligible), T present", E, "T", "T", True)
        run(src, "D118 primary rows (eligible), T present", E, "spec_mech_pct", "spec_mech_pct (same rows)", True)
        for v, lbl in (("D_alone", "(a) D alone"), ("A_alone", "(a) A alone"), ("O_alone", "(a) O alone"),
                       ("T_noD", "(b) T without D"), ("T_noA", "(b) T without A"), ("T_noO", "(b) T without O"),
                       ("T_unsigned", "(c) T unsigned"), ("T_plus_spec", "(h) mean of T and spec_mech_pct_pregate pct")):
            run(src, "primary rows", P, v, lbl, True)
        for v, lbl in (("D_alone", "(a) D alone"), ("A_alone", "(a) A alone"), ("O_alone", "(a) O alone")):
            Q = base[base[v].notna()]
            run(src, f"anticancer, gate-passing, {v} present", Q, v, lbl + ", own rows")
            run(src, f"anticancer, gate-passing, {v} present", Q, "spec_mech_pct_pregate", "spec_mech_pct_pregate, same rows")
        for f, nm in ((True, "common-essential main target"), (False, "other main target")):
            run(src, f"(f) primary rows, {nm}", P[P.essential_target == f], "T", "T")
        cs = sorted(P.cohort.unique())
        for h, sel in (("odd", cs[0::2]), ("even", cs[1::2])):
            run(src, f"(i) split half: groups {h} (sorted)", P[P.cohort.isin(sel)], "T", "T")
        # (j)-(k) added after seeing (a)-(i), not part of the D120 rule: drug-level vs cohort-specific signal; best-of-N targets
        for v in ("T", "D_alone", "A_alone", "O_alone"):
            m = L.groupby("fda_ingredient")[v].transform("mean")
            L[f"{v}_drugmean"], L[f"{v}_dev"] = m, L[v] - m
        P = L[L.anticancer & ~L[f"excl_{src}"] & L["T"].notna()]
        for v, lbl in (("T_drugmean", "(j) drug's mean T over all groups (no cohort information)"),
                       ("T_dev", "(j) T minus the drug's mean T (cohort-specific part)"),
                       ("D_alone_dev", "(j) D minus drug mean"), ("A_alone_dev", "(j) A minus drug mean"), ("O_alone_dev", "(j) O minus drug mean"),
                       ("n_main_targets", "(k) number of main targets as the score")):
            run(src, "primary rows (post hoc)", P, v, lbl, True)
        for k, sel in (("1", P.n_main_targets <= 1), ("2-3", P.n_main_targets.between(1.01, 3)), (">3", P.n_main_targets > 3)):
            run(src, f"(k) primary rows, {k} main targets (post hoc)", P[sel], "T", "T")
    EV = pd.DataFrame(ev)
    EV.to_csv(OUT / "eval_known_use.tsv", sep="\t", index=False)

    # (d) coverage, correlations, variance share
    S = L[L["T"].notna()]
    dec = [dict(item="FDA drug rows (gate-passing, labels table)", value=len(L)),
           dict(item="rows with a protein main target", value=int(L.n_main_targets.gt(0).sum())),
           dict(item="rows with T", value=len(S))]
    for cpt in COMP:
        dec.append(dict(item=f"rows with {cpt} on the best target", value=int(S[cpt].notna().sum())))
        dec.append(dict(item=f"rows with {cpt} on any main target", value=int(L[f"{cpt}_alone"].notna().sum())))
    for n in (1, 2, 3):
        dec.append(dict(item=f"rows with {n} component(s)", value=int(S.n_components.round().eq(n).sum())))
    for a in ("lower", "raise"):
        dec.append(dict(item=f"rows scored as {a}", value=int(S.action.str.contains(a).sum())))
    dec.append(dict(item="rows flagged direction_assumed", value=int(S.action_flag.str.contains("direction_assumed").sum())))
    dec.append(dict(item="rows flagged mixed", value=int(S.action_flag.str.contains("mixed").sum())))
    for a, b in [("D", "A"), ("D", "O"), ("A", "O"), ("D", "T"), ("A", "T"), ("O", "T"), ("T", "spec_mech_pct_pregate"), ("T", "spec_mech_pct")]:
        rr = []
        for _, d in S.groupby("cohort"):
            d = d.dropna(subset=[a, b])
            if len(d) >= 10:
                rr.append(d[a].corr(d[b], method="spearman"))
        dec.append(dict(item=f"median within-group Spearman {a} vs {b}", value=round(float(np.nanmedian(rr)), 3)))
        if b == "T" and a in COMP:
            dec.append(dict(item=f"(d) {a}: median share of T's within-group variance (rho^2)", value=round(float(np.nanmedian(np.square(rr))), 3)))
    # (e) top-10 kept when each component dropped
    for v in ("T_noD", "T_noA", "T_noO", "T_unsigned"):
        kept = []
        for _, d in S.groupby("cohort"):
            if len(d) < 20:
                continue
            a = set(d.nlargest(10, "T").fda_ingredient); b = set(d.dropna(subset=[v]).nlargest(10, v).fda_ingredient)
            kept.append(len(a & b))
        dec.append(dict(item=f"(e) top-10 by T kept by {v} (median of groups, of 10)", value=float(np.median(kept))))
    # how often each component is the largest contributor in the top 10
    top = S.groupby("cohort", group_keys=False).apply(lambda d: d.nlargest(10, "T"))
    lead = top[COMP].idxmax(axis=1)
    for cpt in COMP:
        dec.append(dict(item=f"top-10 rows where {cpt} is the highest component", value=int((lead == cpt).sum())))
    dec.append(dict(item="top-10 rows with A > 90 and alteration share < 50% (minority-carrier driven)",
                    value=int(((top.A > 90) & (top.alt_share < 0.5)).sum())))
    dec.append(dict(item="top-10 rows with a common-essential main target", value=int(top.essential_target.sum())))
    dec.append(dict(item="top-10 rows total", value=len(top)))
    rr = [d["T"].corr(d.n_main_targets, method="spearman") for _, d in S.groupby("cohort")]
    dec.append(dict(item="(k) median within-group Spearman T vs number of main targets (post hoc)", value=round(float(np.nanmedian(rr)), 3)))
    DEC = pd.DataFrame(dec)
    DEC.to_csv(OUT / "decomposition.tsv", sep="\t", index=False)
    top[["cohort", "drug", "gene", "action", "T", "D", "A", "O", "dep_z", "alt_share", "o_raw", "o_source", "spec_mech_pct_pregate",
         "lab_hub", "lab_chembl", "anticancer", "essential_target"]].to_csv(OUT / "top10.tsv", sep="\t", index=False, float_format="%.3f")
    # (g) controls: yardstick genes and breast anti-estrogens, decoded
    ctl = []
    for coh, gene in KU_YARD:
        r = GC[(GC.group_id == coh) & (GC.gene == gene)]
        ctl.append(dict(kind="yardstick gene", cohort=coh, item=gene, **(r[["dep_z", "D", "alt_share", "A", "o_raw", "O", "o_source"]].iloc[0].to_dict() if len(r) else {})))
    for kind, pairs in (("positive control drug", POS_CTL), ("negative control drug", NEG_CTL)):
        for coh, ing in pairs:
            r = L[(L.cohort == coh) & (L.fda_ingredient == ing)]
            ctl.append(dict(kind=kind, cohort=coh, item=ing, **(r[["gene", "action", "T", "target_role_pct", "D", "A", "O", "dep_z", "alt_share", "o_raw",
                                                                   "o_source", "spec_mech_pct_pregate"]].iloc[0].to_dict() if len(r) else {})))
    pd.DataFrame(ctl).to_csv(OUT / "controls.tsv", sep="\t", index=False, float_format="%.3f")
    json.dump(dict(min_alt_share=MIN_F, min_recurrent_tumors=MIN_REC, seed=SEED, n_perm=N_PERM, n_boot=N_BOOT,
                   groups=len(groups), alteration=ainfo), open(OUT / "run.json", "w"), indent=1, default=str)
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 80)
    print(EV[["source", "rows", "score", "auc", "n_cohorts", "n_pos", "n_neg", "p_perm", "ci_lo", "ci_hi"]].to_string(index=False))
    print(DEC.to_string(index=False))
    print(pd.read_csv(OUT / "controls.tsv", sep="\t").to_string(index=False))


KU_YARD = [("TCGA-LUAD", "EGFR"), ("BRCA-Her2", "ERBB2"), ("BRCA-LumA", "ESR1"), ("BRCA-LumB", "ESR1"),
           ("TCGA-PRAD", "AR"), ("TCGA-SKCM", "BRAF"), ("TCGA-STAD", "KIT"), ("TCGA-LIHC", "CTNNB1")]   # 79n YARDSTICK

if __name__ == "__main__":
    main()
