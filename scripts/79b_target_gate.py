#!/usr/bin/env python
"""Per-disease drug lists with a target-presence gate and pharmacology flags (user, 2026-10-09: items 2 and 3).

Drugs per disease (41 Step 4 signature groups): active (is_hiq, TAS >= 0.4) in a LINCS line of the group's strict
tissue key ("direct"), or not direct but active in >= MIN_LINES cancer lines elsewhere ("borrowed"; cross-line
consensus, leave-one-line-out r ~0.46 at 3-4 lines).
Target gate (rules and thresholds in config/drug_gate.yaml; gene classes from HGNC gene groups of the pinned HGNC
snapshot): targets = LINCS compoundinfo + Repurposing Hub (scripts/79a_drug_annotations.py), symbols mapped to
GENCODE v23 Ensembl through HGNC (current, previous and alias symbols; unique hits only) and results/step5/id_map.tsv.
A gene is expressed in a disease if >= min_fraction of the group's primary_per_patient tumors have TPM >= tpm_min.
ADME genes (drug-metabolism / transport groups) are not used as targets. gate_reason, first rule that applies:
  target_expressed         >= 1 non-ADME target expressed                                 -> pass
  family_member_expressed  no target expressed, but another member of one of its specific -> pass
                           HGNC groups (size within family limits) is
  adme_only                only ADME genes annotated                                      -> pass
  no_target                no annotated / mappable target                                 -> pass
  target_not_expressed     none of the above                                              -> fail
Pharmacology flag (locally acting: topical-only route or local ATC group) carried for review, not used to exclude.
STRATA=1 (user 2026-10-12): the same lists for the Step 10 strata (results/step10/stratum_patients.tsv): tissue key,
direct / borrowed split and rules of the parent group; "expressed" judged on the stratum's own tumors. Outputs get the
prefix stratum_ instead of disease_ (target expression: target_expression_strata.tsv.gz); adds parent_group and weak.
Outputs: results/lincs/disease_drugs.tsv (one row per disease x drug), disease_drugs_summary.tsv,
         data/processed/lincs/target_expression.tsv.gz (target gene x group: median TPM, fraction of tumors TPM >= 1)
"""
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
RES = ROOT / "results/lincs"
PROC = ROOT / "data/processed/lincs"
TAS_MIN, MIN_LINES = 0.4, 3
CONFIG = ROOT / "config/drug_gate.yaml"


def hgnc_file():
    return sorted((RAW / "hgnc").glob("hgnc_complete_set_*.txt"))[-1]


def gene_groups(h2g):
    """HGNC gene group id -> set of v23 gene_ids, and gene_id -> group ids"""
    hg = pd.read_csv(hgnc_file(), sep="\t", dtype=str, usecols=["hgnc_id", "symbol", "gene_group_id"])
    hg = hg.dropna(subset=["gene_group_id"]).assign(gid=lambda d: d.gene_group_id.str.split("|")).explode("gid")
    hg["gene_id"] = hg.hgnc_id.map(h2g)
    hg = hg.dropna(subset=["gene_id"])
    return hg.groupby("gid").gene_id.apply(set).to_dict(), hg.groupby("gene_id").gid.apply(set).to_dict(), hg


def symbol_map():
    """symbol -> v23 gene_id via HGNC current / previous / alias symbols (ambiguous names dropped)"""
    hg = pd.read_csv(hgnc_file(), sep="\t", dtype=str, usecols=["hgnc_id", "symbol", "prev_symbol", "alias_symbol"])
    idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t")
    h2g = idm[idm.primary_for_hgnc].drop_duplicates("hgnc_id").set_index("hgnc_id").gene_id
    rows = [(s, r.hgnc_id, 0) for r in hg.itertuples() for s in [r.symbol]]
    for col, rank in [("prev_symbol", 1), ("alias_symbol", 2)]:
        for r in hg.dropna(subset=[col]).itertuples():
            rows += [(s.strip().strip('"'), r.hgnc_id, rank) for s in getattr(r, col).split("|")]
    t = pd.DataFrame(rows, columns=["name", "hgnc_id", "rank"])
    t = t[t["rank"] == t.groupby("name")["rank"].transform("min")]       # current symbol beats previous beats alias
    t = t[~t.duplicated("name", keep=False)]                               # drop names still ambiguous
    return t.set_index("name").hgnc_id.map(h2g).dropna(), h2g


def group_expression(genes, groups, tpm_min, tum=None):
    """tum: optional group_id x sample table (strata); default = the groups' primary_per_patient tumors"""
    if tum is None:
        cov = pd.read_csv(ROOT / "results/step3/covariates_toil.tsv", sep="\t", low_memory=False)
        tum = cov[(cov.role == "tumor") & (cov.primary_per_patient == True) & cov.group_id.isin(groups)]
        tum = tum.drop_duplicates(["group_id", "sample"])
    f = RAW / "solid_toil/TcgaTargetGtex_rsem_gene_tpm.gz"
    header = pd.read_csv(f, sep="\t", nrows=0).columns
    use = sorted(set(tum["sample"]) & set(header))
    parts = []
    for ch in pd.read_csv(f, sep="\t", usecols=[header[0]] + use, chunksize=5000, index_col=0):
        ch.index = ch.index.str.split(".").str[0]
        parts.append(ch[ch.index.isin(genes)])
    x = pd.concat(parts)
    tpm = np.power(2.0, x) - 0.001                                        # log2(TPM+0.001) -> TPM
    out = []
    for g, d in tum[tum["sample"].isin(use)].groupby("group_id"):
        t = tpm[d["sample"]]
        out.append(pd.DataFrame({"gene_id": t.index, "group_id": g, "median_tpm": t.median(axis=1).values,
                                 "frac_expressed": (t >= tpm_min).mean(axis=1).values}))
    return pd.concat(out)


def main():
    PROC.mkdir(parents=True, exist_ok=True)
    import yaml
    cfg = yaml.safe_load(open(CONFIG))
    tpm_min, min_frac = cfg["expression"]["tpm_min"], cfg["expression"]["min_fraction"]
    fmin, fmax = cfg["family"]["min_group_size"], cfg["family"]["max_group_size"]
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t")
    sig = pd.read_csv(ROOT / "results/step4/signature_summary.tsv", sep="\t")
    groups = sorted(set(sig.loc[sig.status == "ok", "group_id"]))
    strata = os.environ.get("STRATA", "0") == "1"
    parent, weak, tum, pfx = {g: g for g in groups}, {}, None, "disease"
    if strata:
        sp = pd.read_csv(ROOT / "results/step10/stratum_patients.tsv", sep="\t")
        groups = sorted(sp.stratum.unique())
        parent = sp.drop_duplicates("stratum").set_index("stratum").group_id.to_dict()
        weak = sp.drop_duplicates("stratum").set_index("stratum").weak.to_dict()
        tum, pfx = sp.rename(columns={"group_id": "parent", "stratum": "group_id"})[["group_id", "sample"]], "stratum"
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
    models = pd.read_csv(RAW / "depmap_24q4/Model.csv", index_col=0)
    cancer_models = models[models.OncotreePrimaryDisease != "Non-Cancerous"]

    # active drug x line (cancer lines)
    cov = pd.read_csv(RES / "cell_coverage.tsv", sep="\t").dropna(subset=["ModelID"])
    s = pd.read_csv(RAW / "lincs2020/siginfo_beta.txt", sep="\t", low_memory=False,
                    usecols=["pert_id", "pert_type", "cell_iname", "is_hiq", "tas"])
    s = s[(s.pert_type == "trt_cp") & (s.is_hiq == 1) & (s.tas >= TAS_MIN)].merge(cov[["cell_iname", "ModelID"]])
    act = s[s.ModelID.isin(cancer_models.index)][["pert_id", "ModelID"]].drop_duplicates()

    # targets -> genes
    smap, h2g = symbol_map()
    grp_genes, gene_grps, hgg = gene_groups(h2g)
    adme = set().union(*(grp_genes.get(str(i), set()) for i in cfg["adme"]["hgnc_group_ids"]))
    adme |= set(hgg.loc[hgg.symbol.isin(cfg["adme"]["genes"]), "gene_id"])
    fam_ok = {g: m for g, m in grp_genes.items() if fmin <= len(m) <= fmax} if cfg["family"].get("enabled", True) else {}
    tg = ann[["pert_id", "lincs_target", "hub_target"]].copy()
    tg["symbols"] = [sorted({x.strip() for v in (a, b) if isinstance(v, str) for x in v.split("|") if x.strip()})
                     for a, b in zip(tg.lincs_target, tg.hub_target)]
    tg = tg.explode("symbols").dropna(subset=["symbols"])
    tg["gene_id"] = tg.symbols.map(smap)
    print(f"drug-target pairs: {len(tg)}; symbols mapped to v23 genes: {tg.gene_id.notna().mean():.1%}; "
          f"drugs with >= 1 mapped target: {tg.dropna(subset=['gene_id']).pert_id.nunique()}")
    targets = set(tg.gene_id.dropna()) - adme
    family = {t: set().union(*(fam_ok[g] for g in gene_grps.get(t, set()) if g in fam_ok)) - {t} for t in targets}
    expr = group_expression(targets | set().union(*family.values()), groups, tpm_min, tum)
    expr.round(4).to_csv(PROC / ("target_expression_strata.tsv.gz" if strata else "target_expression.tsv.gz"),
                         sep="\t", index=False)
    expressed = expr[expr.frac_expressed >= min_frac].groupby("group_id").gene_id.apply(set).to_dict()
    gid2sym = tg.dropna(subset=["gene_id"]).drop_duplicates("gene_id").set_index("gene_id").symbols
    drug_targets = tg.dropna(subset=["gene_id"]).groupby("pert_id").gene_id.apply(set).to_dict()
    print(f"ADME genes excluded as targets: {len(adme)} ({len(set(tg.gene_id.dropna()) & adme)} occur as drug targets); "
          f"targets with a family (groups of {fmin}-{fmax}): {sum(bool(v) for v in family.values())} of {len(targets)}")

    keep = ["cmap_name", "n_active_lines", "hub_clinical_phase", "chembl_max_phase", "hub_moa", "organ_systems",
            "approved_indications", "route_topical_only", "atc_local", "pharmacology_flag"]
    a = ann.set_index("pert_id")[keep]
    rows, summ = [], []
    for g in groups:
        k = keys.loc[parent[g]]
        km = cancer_models[cancer_models.OncotreeLineage.isin(str(k.lineages).split(";"))]
        if isinstance(k.diseases, str) and k.diseases:
            km = km[km.OncotreePrimaryDisease.isin(k.diseases.split(";"))]
        direct = set(act.loc[act.ModelID.isin(km.index), "pert_id"])
        ne = act[~act.ModelID.isin(km.index)].groupby("pert_id").ModelID.nunique()
        borrowed = set(ne[ne >= MIN_LINES].index) - direct
        ex = expressed.get(g, set())
        for d in sorted(direct | borrowed):
            t_all = drug_targets.get(d, set())
            t = t_all - adme
            te = t & ex
            fam = {x for x in t - te if family.get(x, set()) & ex}
            if te:
                status = "target_expressed"
            elif fam:
                status = "family_member_expressed"
            elif t_all and not t:
                status = "adme_only"
            elif not t_all:
                status = "no_target"
            else:
                status = "target_not_expressed"
            rows.append(dict(group_id=g, key=k.key, pert_id=d, source="direct" if d in direct else "borrowed",
                             n_lines_elsewhere=int(ne.get(d, 0)), gate_reason=status, n_targets=len(t),
                             targets_expressed="|".join(sorted(gid2sym.get(x, x) for x in te)),
                             targets_family_rescued="|".join(sorted(gid2sym.get(x, x) for x in fam)),
                             targets_not_expressed="|".join(sorted(gid2sym.get(x, x) for x in t - te - fam)),
                             targets_adme="|".join(sorted(gid2sym.get(x, x) for x in t_all & adme))))
    out = pd.DataFrame(rows).join(a, on="pert_id")
    out["passes_gate"] = out.gate_reason != "target_not_expressed"
    if strata:
        out.insert(1, "parent_group", out.group_id.map(parent))
        out.insert(2, "weak", out.group_id.map(weak))
    out.to_csv(RES / f"{pfx}_drugs.tsv", sep="\t", index=False)

    for g, d in out.groupby("group_id"):
        summ.append(dict(group_id=g, key=d.key.iloc[0], n_direct=int((d.source == "direct").sum()),
                         n_borrowed=int((d.source == "borrowed").sum()), n_total=len(d),
                         **{f"n_{r}": int((d.gate_reason == r).sum()) for r in
                            ["target_expressed", "family_member_expressed", "adme_only", "no_target", "target_not_expressed"]},
                         n_pass=int(d.passes_gate.sum()),
                         n_pass_flagged_local=int((d.passes_gate & d.pharmacology_flag.eq(True)).sum()),
                         n_pass_approved=int((d.passes_gate & (pd.to_numeric(d.chembl_max_phase, errors="coerce") >= 4)).sum())))
    sm = pd.DataFrame(summ).sort_values("n_pass", ascending=False)
    sm.to_csv(RES / f"{pfx}_drugs_summary.tsv", sep="\t", index=False)
    sha = lambda f: hashlib.sha256(open(f, "rb").read()).hexdigest()[:16]
    run = {"config": str(CONFIG.relative_to(ROOT)), "config_sha256_16": sha(CONFIG), "rules": cfg,
           "hgnc_file": hgnc_file().name, "hgnc_sha256_16": sha(hgnc_file()), "tas_min": TAS_MIN,
           "min_lines_borrowed": MIN_LINES, "n_adme_genes": len(adme), "n_rows": len(out)}
    json.dump(run, open(RES / f"{pfx}_drugs_run.json", "w"), indent=1)
    if strata:
        print(sm.to_string(index=False))
        return
    pd.set_option("display.width", 250)
    print(sm.to_string(index=False))
    # example: drugs removed in kidney because their targets are not expressed there
    chk = out[out.cmap_name.isin(["simeprevir", "thapsigargin", "BIBR-1532", "doxycycline", "abiraterone", "letrozole"])]
    print(chk[chk.group_id == "TCGA-KIRC"][["cmap_name", "gate_reason", "targets_expressed", "targets_family_rescued",
                                             "targets_adme", "targets_not_expressed"]].to_string(index=False))
    x = out[(out.group_id == "TCGA-KIRC") & (out.gate_reason == "target_not_expressed")]
    print(x[["cmap_name", "hub_moa", "targets_not_expressed"]].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
