#!/usr/bin/env python
"""Literature and trial links for the top v2 drug-disease pairs (user 2026-10-29; runs in `celligner`).

Pairs: per disease (41 groups), the union of the top TOP_N drugs by spec_mech (reports lead with it, D80) and by rev,
among drugs that pass the target gate, are not flagged (unmatched_lines / low_transfer) and have rev > 0
(results/lincs/drug_reversal_split_v2_consensus.tsv). Unnamed LINCS compounds (BRD-...) are skipped.
Searches (same for every pair; no manual curation):
  PubMed       esearch, drug name in title/abstract (also without hyphen for codes such as AZD-9291) AND the disease's
               terms (DISEASE below) in title/abstract; count and the RETMAX most relevant PMIDs with titles
  ClinicalTrials.gov  API v2, intervention = drug, condition = the disease's trial condition; count and up to RETMAX
               studies (NCT id, phase, status, title), most relevant first
A hit means the pair has been studied, not that it worked; counts are a pointer for reading, not evidence.
ClinicalTrials.gov's relevance search matches loosely (related drugs, e.g. daunorubicin -> anthracycline-regimen trials;
'Wilms tumor' -> WT1 trials), so trial_count is an upper bound; read the linked trials. PubMed (title / abstract) is stricter.
Output: results/lincs/literature_v2_consensus.tsv (one row per pair, links as URLs), literature_v2_summary.tsv
MODE=fda (user 2026-10-29): the top FDA_N FDA-approved drugs of every disease and subgroup (79p rankings), searched by
  their FDA ingredient name (salt removed) OR their LINCS name when it is a real name (not BRD-... or a code); subgroups
  use their parent disease's terms. Also drug_review_url: PubMed reviews about the drug itself. Output
  results/lincs/literature_fda_top_v2.tsv (unit, level, parent, fda_rank, drug, ...).
MODE=fda_all (user 2026-11-02): every FDA drug row of the rankings page (scripts/fda_rows.py: target-gate pass, Drugs@FDA
  match, one per ingredient; flagged drugs included), same searches as MODE=fda. A subgroup uses its parent's terms, so
  each disease x drug is searched once and shared with its subgroups. PubMed and ClinicalTrials.gov run in two threads
  (PubMed <= 3 requests / s, NCBI's limit without a key); answers are cached line by line in
  data/interim/lincs/literature_fda_all_{pubmed,trials}.jsonl, so an interrupted run resumes. Output
  results/lincs/literature_fda_all_v2.tsv (unit, level, parent_group, pert_id, drug, ...).
"""
import os
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
TOP_N, RETMAX, FDA_N = 10, 3, 3
MODE = os.environ.get("MODE", "pairs")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# disease -> (PubMed title/abstract terms, ClinicalTrials.gov condition)
DISEASE = {
    "BRCA-Basal": (['triple-negative breast', 'triple negative breast', 'basal-like breast'], "triple negative breast cancer"),
    "BRCA-Her2": (['HER2-positive breast', 'HER2 positive breast', 'HER2-enriched'], "HER2-positive breast cancer"),
    "BRCA-LumA": (['luminal A', 'hormone receptor-positive breast', 'ER-positive breast', 'estrogen receptor-positive breast'], "hormone receptor positive breast cancer"),
    "BRCA-LumB": (['luminal B', 'hormone receptor-positive breast', 'ER-positive breast', 'estrogen receptor-positive breast'], "hormone receptor positive breast cancer"),
    "CESC-AD": (['cervical adenocarcinoma', 'cervical cancer'], "cervical cancer"),
    "CESC-SCC": (['cervical squamous', 'cervical cancer'], "cervical cancer"),
    "ESCA-AD": (['esophageal adenocarcinoma', 'oesophageal adenocarcinoma', 'gastroesophageal junction'], "esophageal adenocarcinoma"),
    "ESCA-SCC": (['esophageal squamous', 'oesophageal squamous'], "esophageal squamous cell carcinoma"),
    "LGG-IDHmut-codel": (['oligodendroglioma', 'lower-grade glioma', 'low-grade glioma', 'IDH-mutant glioma'], "oligodendroglioma"),
    "LGG-IDHmut-noncodel": (['astrocytoma', 'lower-grade glioma', 'low-grade glioma', 'IDH-mutant glioma'], "astrocytoma"),
    "LGG-IDHwt": (['IDH-wildtype glioma', 'lower-grade glioma', 'low-grade glioma', 'glioma'], "glioma"),
    "PCPG-PGL": (['paraganglioma'], "paraganglioma"),
    "PCPG-PHEO": (['pheochromocytoma'], "pheochromocytoma"),
    "SARC-DDLPS": (['dedifferentiated liposarcoma', 'liposarcoma'], "liposarcoma"),
    "SARC-LMS": (['leiomyosarcoma'], "leiomyosarcoma"),
    "SARC-UPS": (['undifferentiated pleomorphic sarcoma', 'malignant fibrous histiocytoma', 'myxofibrosarcoma'], "undifferentiated pleomorphic sarcoma"),
    "TARGET-NBL": (['neuroblastoma'], "neuroblastoma"),
    "TARGET-WT": (['Wilms tumor', 'Wilms tumour', 'nephroblastoma'], "Wilms tumor"),
    "TCGA-ACC": (['adrenocortical carcinoma', 'adrenocortical cancer'], "adrenocortical carcinoma"),
    "TCGA-BLCA": (['bladder cancer', 'urothelial carcinoma', 'bladder carcinoma'], "bladder cancer"),
    "TCGA-COAD": (['colon cancer', 'colorectal cancer', 'colon adenocarcinoma'], "colorectal cancer"),
    "TCGA-GBM": (['glioblastoma'], "glioblastoma"),
    "TCGA-HNSC": (['head and neck squamous', 'head and neck cancer', 'oral squamous'], "head and neck squamous cell carcinoma"),
    "TCGA-KICH": (['chromophobe renal', 'chromophobe'], "chromophobe renal cell carcinoma"),
    "TCGA-KIRC": (['clear cell renal', 'renal cell carcinoma'], "clear cell renal cell carcinoma"),
    "TCGA-KIRP": (['papillary renal', 'renal cell carcinoma'], "papillary renal cell carcinoma"),
    "TCGA-LIHC": (['hepatocellular carcinoma', 'liver cancer'], "hepatocellular carcinoma"),
    "TCGA-LUAD": (['lung adenocarcinoma', 'non-small cell lung', 'non-small-cell lung'], "non-small cell lung cancer"),
    "TCGA-LUSC": (['lung squamous', 'squamous non-small cell lung', 'squamous cell lung'], "squamous non-small cell lung cancer"),
    "TCGA-OV": (['ovarian cancer', 'ovarian carcinoma', 'high-grade serous'], "ovarian cancer"),
    "TCGA-PAAD": (['pancreatic cancer', 'pancreatic adenocarcinoma', 'pancreatic ductal'], "pancreatic cancer"),
    "TCGA-PRAD": (['prostate cancer', 'prostate adenocarcinoma'], "prostate cancer"),
    "TCGA-READ": (['rectal cancer', 'colorectal cancer'], "colorectal cancer"),
    "TCGA-SKCM": (['melanoma'], "melanoma"),
    "TCGA-STAD": (['gastric cancer', 'stomach cancer', 'gastric adenocarcinoma'], "gastric cancer"),
    "TCGA-THCA": (['thyroid cancer', 'thyroid carcinoma', 'papillary thyroid'], "thyroid cancer"),
    "TCGA-UCS": (['uterine carcinosarcoma', 'carcinosarcoma'], "uterine carcinosarcoma"),
    "TGCT-NSGCT": (['nonseminoma', 'non-seminomatous', 'germ cell tumor', 'germ cell tumour'], "germ cell tumor"),
    "TGCT-SEM": (['seminoma', 'germ cell tumor', 'germ cell tumour'], "seminoma"),
    "UCEC-EEC": (['endometrial cancer', 'endometrioid', 'endometrial carcinoma'], "endometrial cancer"),
    "UCEC-SER": (['serous endometrial', 'uterine serous', 'endometrial cancer'], "endometrial cancer"),
}


def get(url, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.loads(r.read().decode())
        except Exception:
            time.sleep(2 * (k + 1))
    return None


def names(drug):
    v = {drug}
    if re.match(r"^[A-Za-z]+-\d", drug):
        v.add(drug.replace("-", ""))
    v.add(drug.replace("-", " ") if not re.match(r"^[A-Za-z]+-\d", drug) else drug)
    return sorted(v)


def pubmed(drug, terms, extra=()):
    q = "(" + " OR ".join(f'"{n}"[tiab]' for n in sorted(set(names(drug)) | set(extra))) + ") AND (" + " OR ".join(f'"{t}"[tiab]' for t in terms) + ")"
    r = get(f"{EUTILS}/esearch.fcgi?db=pubmed&retmode=json&sort=relevance&retmax={RETMAX}&term=" + urllib.parse.quote(q))
    time.sleep(0.35)
    if not r:
        return q, None, []
    return q, int(r["esearchresult"]["count"]), r["esearchresult"]["idlist"]


def trials(drug, cond):
    url = ("https://clinicaltrials.gov/api/v2/studies?countTotal=true&pageSize=%d&fields=NCTId,BriefTitle,Phase,OverallStatus"
           "&query.intr=%s&query.cond=%s" % (RETMAX, urllib.parse.quote(drug), urllib.parse.quote(cond)))
    r = get(url)
    if not r:
        return None, []
    out = []
    for s in r.get("studies", []):
        p = s["protocolSection"]
        out.append((p["identificationModule"]["nctId"], "/".join(p.get("designModule", {}).get("phases", []) or ["NA"]),
                    p["statusModule"]["overallStatus"], p["identificationModule"]["briefTitle"]))
    return r.get("totalCount"), out


def main():
    o = pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
    ok = o[o.passes_gate & ~o.unmatched_lines & ~o.low_transfer & (o.rev > 0) & ~o.cmap_name.str.startswith("BRD-")]
    a = ok.sort_values("spec_mech", ascending=False).groupby("group_id").head(TOP_N).assign(top_spec_mech=True)
    b = ok.sort_values("rev", ascending=False).groupby("group_id").head(TOP_N).assign(top_rev=True)
    sel = pd.concat([a, b]).groupby(["group_id", "cmap_name"], as_index=False).agg(
        source=("source", "first"), approved=("approved", "first"), hub_moa=("hub_moa", "first"),
        spec_mech=("spec_mech", "max"), spec_mech_pct=("spec_mech_pct", "max"), rev=("rev", "max"), rev_pct=("rev_pct", "max"),
        transfer_conf=("transfer_conf", "first"), top_spec_mech=("top_spec_mech", "any"), top_rev=("top_rev", "any"))
    sel[["top_spec_mech", "top_rev"]] = sel[["top_spec_mech", "top_rev"]].fillna(False)
    rows, pmids = [], set()
    t0 = time.time()
    for i, r in enumerate(sel.itertuples(index=False)):
        terms, cond = DISEASE[r.group_id]
        q, n_pm, ids = pubmed(r.cmap_name, terms)
        n_ct, st = trials(r.cmap_name, cond)
        pmids.update(ids)
        rows.append(dict(r._asdict(), pubmed_count=n_pm, pmids=ids, trial_count=n_ct, trials=st, pubmed_query=q))
        if i % 50 == 0:
            print(f"{i}/{len(sel)} pairs ({time.time() - t0:.0f}s)", flush=True)
    titles = {}
    ids = sorted(pmids)
    for k in range(0, len(ids), 150):
        r = get(f"{EUTILS}/esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(ids[k:k + 150]))
        time.sleep(0.35)
        if r:
            for pid in r["result"].get("uids", []):
                d = r["result"][pid]
                titles[pid] = f'{d.get("title", "").rstrip(".")} ({d.get("source", "")} {d.get("pubdate", "")[:4]})'
    out = []
    for r in rows:
        d = {k: v for k, v in r.items() if k not in ("pmids", "trials")}
        d["pubmed_search_url"] = "https://pubmed.ncbi.nlm.nih.gov/?term=" + urllib.parse.quote(r["pubmed_query"])
        for j in range(RETMAX):
            pid = r["pmids"][j] if j < len(r["pmids"]) else ""
            d[f"paper{j + 1}"] = f"https://pubmed.ncbi.nlm.nih.gov/{pid}/" if pid else ""
            d[f"paper{j + 1}_title"] = titles.get(pid, "")
        d["trials_search_url"] = ("https://clinicaltrials.gov/search?intr=" + urllib.parse.quote(r["cmap_name"])
                                  + "&cond=" + urllib.parse.quote(DISEASE[r["group_id"]][1]))
        for j in range(RETMAX):
            t = r["trials"][j] if j < len(r["trials"]) else None
            d[f"trial{j + 1}"] = f"https://clinicaltrials.gov/study/{t[0]}" if t else ""
            d[f"trial{j + 1}_info"] = f"{t[1]}, {t[2]}: {t[3]}" if t else ""
        out.append(d)
    out = pd.DataFrame(out).sort_values(["group_id", "spec_mech"], ascending=[True, False])
    out.to_csv(RES / "literature_v2_consensus.tsv", sep="\t", index=False)
    sm = out.groupby("group_id").agg(pairs=("cmap_name", "size"), with_pubmed=("pubmed_count", lambda v: int((v > 0).sum())),
                                     with_trials=("trial_count", lambda v: int((v > 0).sum())),
                                     top_spec_mech_with_pubmed=("pubmed_count", lambda v: int((v[out.loc[v.index, "top_spec_mech"]] > 0).sum())))
    sm.to_csv(RES / "literature_v2_summary.tsv", sep="\t")
    print(sm.to_string()); print(f"pairs {len(out)}; with PubMed hits {int((out.pubmed_count > 0).sum())}; with trials {int((out.trial_count > 0).sum())}")


def is_code(n):
    return n.startswith("BRD-") or bool(re.match(r"^[A-Za-z]{1,5}-?\d", n))


def main_fda():
    g = pd.read_csv(RES / "fda_ranking_v2_consensus.tsv", sep="\t").rename(columns={"group_id": "unit"}).assign(level="disease", parent_group=None)
    s = pd.read_csv(RES / "fda_ranking_strata_v2_consensus.tsv", sep="\t").rename(columns={"stratum": "unit"}).assign(level="subgroup")
    sel = pd.concat([g, s], ignore_index=True)
    sel = sel[sel.fda_rank <= FDA_N].copy()
    sel["disease_key"] = sel.parent_group.fillna(sel.unit)
    rows, pmids = [], set()
    for r in sel.itertuples(index=False):
        terms, cond = DISEASE[r.disease_key]
        extra = () if is_code(r.cmap_name) else tuple(names(r.cmap_name))
        q, n_pm, ids = pubmed(r.drug, terms, extra)
        n_ct, st = trials(r.drug, cond)
        pmids.update(ids)
        rows.append(dict(unit=r.unit, level=r.level, parent_group=r.parent_group, fda_rank=r.fda_rank, drug=r.drug,
                         cmap_name=r.cmap_name, pubmed_count=n_pm, pmids=ids, trial_count=n_ct, trials=st, pubmed_query=q,
                         trial_cond=cond))
    titles = {}
    ids = sorted(pmids)
    for k in range(0, len(ids), 150):
        r = get(f"{EUTILS}/esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(ids[k:k + 150]))
        time.sleep(0.35)
        if r:
            for pid in r["result"].get("uids", []):
                d = r["result"][pid]
                titles[pid] = f'{d.get("title", "").rstrip(".")} ({d.get("source", "")} {d.get("pubdate", "")[:4]})'
    out = []
    for r in rows:
        d = {k: v for k, v in r.items() if k not in ("pmids", "trials", "trial_cond")}
        d["pubmed_search_url"] = "https://pubmed.ncbi.nlm.nih.gov/?term=" + urllib.parse.quote(r["pubmed_query"])
        d["drug_review_url"] = ("https://pubmed.ncbi.nlm.nih.gov/?sort=relevance&term="
                                + urllib.parse.quote(f'"{r["drug"]}"[tiab] AND review[pt]'))
        for j in range(RETMAX):
            pid = r["pmids"][j] if j < len(r["pmids"]) else ""
            d[f"paper{j + 1}"] = f"https://pubmed.ncbi.nlm.nih.gov/{pid}/" if pid else ""
            d[f"paper{j + 1}_title"] = titles.get(pid, "")
        d["trials_search_url"] = ("https://clinicaltrials.gov/search?intr=" + urllib.parse.quote(r["drug"])
                                  + "&cond=" + urllib.parse.quote(r["trial_cond"]))
        for j in range(RETMAX):
            t = r["trials"][j] if j < len(r["trials"]) else None
            d[f"trial{j + 1}"] = f"https://clinicaltrials.gov/study/{t[0]}" if t else ""
            d[f"trial{j + 1}_info"] = f"{t[1]}, {t[2]}: {t[3]}" if t else ""
        out.append(d)
    out = pd.DataFrame(out)
    out.to_csv(RES / "literature_fda_top_v2.tsv", sep="\t", index=False)
    print(f"{len(out)} cohort x drug rows ({out.unit.nunique()} cohorts); PubMed hits {int((out.pubmed_count > 0).sum())}; "
          f"trials {int((out.trial_count > 0).sum())}")
    print(out[["unit", "fda_rank", "drug", "pubmed_count", "trial_count"]].head(30).to_string(index=False))


def main_fda_all():
    import sys
    import threading
    sys.path.insert(0, str(ROOT / "scripts"))
    from fda_rows import fda_rows
    g, s = fda_rows()
    g = g.assign(level="disease", parent_group=None)
    s = s.assign(level="subgroup")
    sel = pd.concat([g, s], ignore_index=True)
    sel["disease_key"] = sel.parent_group.fillna(sel.group_id)
    sel["extra"] = sel.cmap_name.map(lambda n: "" if is_code(n) else "|".join(names(n)))
    q = sel.drop_duplicates(["disease_key", "drug", "extra"])[["disease_key", "drug", "extra"]]
    cache = ROOT / "data/interim/lincs"
    cache.mkdir(parents=True, exist_ok=True)
    print(f"{len(sel)} cohort x drug rows -> {len(q)} distinct searches", flush=True)

    def load(f):
        d = {}
        if f.exists():
            for line in open(f):
                r = json.loads(line)
                d[tuple(r["k"])] = r["v"]
        return d

    def worker(kind):
        f = cache / f"literature_fda_all_{kind}.jsonl"
        done = load(f)
        todo = [r for r in q.itertuples(index=False) if (r.disease_key, r.drug, r.extra) not in done]
        t0 = time.time()
        with open(f, "a") as fh:
            for i, r in enumerate(todo):
                terms, cond = DISEASE[r.disease_key]
                if kind == "pubmed":
                    qq, n, ids = pubmed(r.drug, terms, tuple(r.extra.split("|")) if r.extra else ())
                    v = dict(q=qq, n=n, ids=ids)
                else:
                    n, st = trials(r.drug, cond)
                    v = dict(n=n, st=st, cond=cond)
                if v["n"] is None:      # failed request: not cached, retried on the next run
                    continue
                fh.write(json.dumps(dict(k=[r.disease_key, r.drug, r.extra], v=v)) + "\n")
                fh.flush()
                if i % 250 == 0:
                    print(f"{kind}: {i}/{len(todo)} ({time.time() - t0:.0f}s)", flush=True)
    th = [threading.Thread(target=worker, args=(k,)) for k in ("pubmed", "trials")]
    [t.start() for t in th]
    [t.join() for t in th]
    pm, ct = load(cache / "literature_fda_all_pubmed.jsonl"), load(cache / "literature_fda_all_trials.jsonl")
    ids = sorted({i for v in pm.values() for i in v["ids"]})
    tf = cache / "literature_fda_all_titles.json"
    titles = json.load(open(tf)) if tf.exists() else {}
    need = [i for i in ids if i not in titles]
    for k in range(0, len(need), 150):
        r = get(f"{EUTILS}/esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(need[k:k + 150]))
        time.sleep(0.35)
        if r:
            for pid in r["result"].get("uids", []):
                d = r["result"][pid]
                titles[pid] = f'{d.get("title", "").rstrip(".")} ({d.get("source", "")} {d.get("pubdate", "")[:4]})'
    json.dump(titles, open(tf, "w"))
    out = []
    for r in sel.itertuples(index=False):
        k = (r.disease_key, r.drug, r.extra)
        P, C = pm.get(k), ct.get(k)
        d = dict(unit=r.group_id, level=r.level, parent_group=r.parent_group, pert_id=r.pert_id, cmap_name=r.cmap_name,
                 drug=r.drug, pubmed_count=P["n"] if P else None, trial_count=C["n"] if C else None,
                 pubmed_query=P["q"] if P else "")
        d["pubmed_search_url"] = "https://pubmed.ncbi.nlm.nih.gov/?term=" + urllib.parse.quote(d["pubmed_query"])
        d["drug_review_url"] = ("https://pubmed.ncbi.nlm.nih.gov/?sort=relevance&term="
                                + urllib.parse.quote(f'"{r.drug}"[tiab] AND review[pt]'))
        pids = P["ids"] if P else []
        for j in range(RETMAX):
            pid = pids[j] if j < len(pids) else ""
            d[f"paper{j + 1}"] = f"https://pubmed.ncbi.nlm.nih.gov/{pid}/" if pid else ""
            d[f"paper{j + 1}_title"] = titles.get(pid, "")
        d["trials_search_url"] = ("https://clinicaltrials.gov/search?intr=" + urllib.parse.quote(r.drug)
                                  + "&cond=" + urllib.parse.quote(DISEASE[r.disease_key][1]))
        st = C["st"] if C else []
        for j in range(RETMAX):
            t = st[j] if j < len(st) else None
            d[f"trial{j + 1}"] = f"https://clinicaltrials.gov/study/{t[0]}" if t else ""
            d[f"trial{j + 1}_info"] = f"{t[1]}, {t[2]}: {t[3]}" if t else ""
        out.append(d)
    out = pd.DataFrame(out)
    out.to_csv(RES / "literature_fda_all_v2.tsv", sep="\t", index=False)
    print(f"{len(out)} rows; searched {int(out.pubmed_count.notna().sum())}; PubMed hits {int((out.pubmed_count > 0).sum())}; "
          f"trials {int((out.trial_count > 0).sum())}; failed searches pubmed {int(out.pubmed_count.isna().sum())} "
          f"trials {int(out.trial_count.isna().sum())}")


if __name__ == "__main__":
    {"fda": main_fda, "fda_all": main_fda_all}.get(MODE, main)()
