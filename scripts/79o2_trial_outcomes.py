#!/usr/bin/env python
"""Why trials failed: stopped-trial reasons and completed-trial outcomes for the top FDA drugs (D108, user 2026-11-06).

Pairs: fda_rank <= TOP_N in results/lincs/fda_ranking_v2_consensus.tsv (diseases) and fda_ranking_strata_v2_consensus.tsv
(subgroups, parent's trial condition). Trial condition per disease from 79o (DISEASE).
STAGE=fetch
  ClinicalTrials.gov API v2 (intervention = FDA drug name, condition), all pages, cached in
  data/interim/lincs/trial_outcomes_ctgov.jsonl. A trial is kept only if the drug is a listed intervention (any of its
  names - FDA, LINCS, ChEMBL, Hub - salt words removed, whole word, case-insensitive, in an intervention name or
  other-name): the search itself matches loosely (other drugs, regimens).
  Stopped trials (TERMINATED, WITHDRAWN, SUSPENDED): whyStopped sorted by the D108 keyword rules (stop_class).
  Completed trials: PMIDs of linked papers (reference type RESULT or DERIVED), abstracts from PubMed (cached in
  data/interim/lincs/trial_outcomes_abstracts.json); one reading-queue item per trial with >= 1 abstract:
  results/lincs/d108/reading_queue.jsonl.
  Output: results/lincs/d108/trials.tsv (pair x trial).
STAGE=merge
  Reads the model answers results/lincs/d108/readings/*.jsonl ({nct, outcome, reason_class, reason, quote, pmid,
  drug_role}; prompt results/lincs/d108/reading_prompt.md); trials where the drug is background in every arm are not
  counted as the drug failing;
  an answer is accepted only if its quote occurs verbatim (whitespace / case normalized) in that trial's abstracts.
  Outputs: results/lincs/d108/trials_read.tsv (trials.tsv + model columns), results/lincs/trial_outcomes_v2.tsv
  (one row per pair: counts and the failed / stopped trials with reasons).
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
OUT = RES / "d108"
CACHE = ROOT / "data/interim/lincs"
TOP_N = int(os.environ.get("TOP_N", 10))
STAGE = os.environ.get("STAGE", "fetch")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
sys.path.insert(0, str(ROOT / "scripts"))
import importlib.util  # noqa: E402
_s = importlib.util.spec_from_file_location("o79", ROOT / "scripts/79o_literature.py"); o79 = importlib.util.module_from_spec(_s); _s.loader.exec_module(o79)
_s = importlib.util.spec_from_file_location("p79", ROOT / "scripts/79p_fda_ranking.py"); p79 = importlib.util.module_from_spec(_s); _s.loader.exec_module(p79)
DISEASE = o79.DISEASE
STOPPED = {"TERMINATED", "WITHDRAWN", "SUSPENDED"}
FIELDS = ("NCTId,BriefTitle,OverallStatus,WhyStopped,Phase,HasResults,InterventionName,InterventionOtherName,"
          "Condition,ReferencePMID,ReferenceType,StartDate,CompletionDate,EnrollmentCount")

# D108 stop-reason rules, first match wins
RULES = [
    ("efficacy", r"futil|(lack|absence|insufficient|no|limited|poor|low)\s+(of\s+)?(clinical\s+)?(efficacy|benefit|activity|response|effect)|"
                 r"interim analys|did not meet|not meet (the|its)|ineffective|negative (result|outcome|data)|efficacy (was |is )?(not|insufficient)"),
    ("safety", r"toxicit|adverse (event|effect|reaction)|safety (concern|issue|signal|reason)|side effect|death|not tolera[tb]|intolerab|unacceptable risk"),
    ("accrual", r"accru|enrol|recruit|participant|subjects? (were|was) (not|difficult)|patient availability"),
    ("business", r"sponsor|fund|financ|business|compan|commercial|strateg|portfolio|drug (supply|availab)|supply|investigator (left|leaving|relocat)|PI (left|leaving)"),
]
NO_SAFETY = re.compile(r"no (safety|toxicity) (issue|concern|signal)|not (due|related) to safety|without safety", re.I)


# D119: negated clauses ("not due to any safety concern", "Not related to safety concerns or lack of efficacy") and the death
# of a study person are removed before the efficacy / safety rules
NEG = re.compile(r"\b(?:not|no|non|without|unrelated)\b(?:\s+\w+){0,2}?\s+(?:due|related|relat\w*|based|driven|a consequence|consequence|"
                 r"because|any|major|new|efficacy analys\w*)\b[^.;,]*", re.I)   # negation introducing a cause ("not due to ...")
STAFF_DEATH = re.compile(r"death of (the |a |our |an )?(co-?)?(principal |study |lead )?(investigator|PI|physician|staff|sponsor|coordinator)\w*", re.I)


def stop_class(why):
    if not isinstance(why, str) or not why.strip():
        return "not_given"
    clean = NEG.sub(" ", STAFF_DEATH.sub(" ", why))
    for name, pat in RULES:
        text = clean if name in ("efficacy", "safety") else why
        if name == "safety" and NO_SAFETY.search(why):
            continue
        if re.search(pat, text, re.I):
            return name
    return "other"


def get(url):
    for k in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "repurposing-pipeline"}), timeout=60) as r:
                return r.read()
        except Exception:
            time.sleep(2 * (k + 1))
    return None


def ctgov(drug, cond):
    out, token = [], None
    while True:
        u = ("https://clinicaltrials.gov/api/v2/studies?pageSize=200&fields=" + FIELDS + "&query.intr=" + urllib.parse.quote(drug)
             + "&query.cond=" + urllib.parse.quote(cond) + (f"&pageToken={token}" if token else ""))
        b = get(u)
        if b is None:
            return None
        r = json.loads(b)
        for s in r.get("studies", []):
            p = s["protocolSection"]
            iv = p.get("armsInterventionsModule", {}).get("interventions", [])
            out.append(dict(
                nct=p["identificationModule"]["nctId"], title=p["identificationModule"].get("briefTitle", ""),
                status=p["statusModule"].get("overallStatus", ""), why_stopped=p["statusModule"].get("whyStopped", ""),
                phase="/".join(p.get("designModule", {}).get("phases", []) or ["NA"]), has_results=bool(s.get("hasResults")),
                interventions=[i.get("name", "") for i in iv], other_names=[n for i in iv for n in i.get("otherNames", []) or []],
                conditions=p.get("conditionsModule", {}).get("conditions", []),
                refs=[(x.get("pmid"), x.get("type")) for x in p.get("referencesModule", {}).get("references", []) if x.get("pmid")],
                start=p["statusModule"].get("startDateStruct", {}).get("date", ""),
                completion=p["statusModule"].get("completionDateStruct", {}).get("date", ""),
                enrollment=p.get("designModule", {}).get("enrollmentInfo", {}).get("count")))
        token = r.get("nextPageToken")
        if not token:
            return out
        time.sleep(0.3)


def pairs():
    g = pd.read_csv(RES / "fda_ranking_v2_consensus.tsv", sep="\t").rename(columns={"group_id": "unit"}).assign(level="disease", parent_group=None)
    s = pd.read_csv(RES / "fda_ranking_strata_v2_consensus.tsv", sep="\t").rename(columns={"stratum": "unit"}).assign(level="subgroup")
    x = pd.concat([g, s], ignore_index=True)
    x = x[x.fda_rank <= TOP_N].copy()
    x["disease_key"] = x.parent_group.fillna(x.unit)
    x["condition"] = x.disease_key.map(lambda k: DISEASE[k][1])
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t", low_memory=False).drop_duplicates("cmap_name").set_index("cmap_name")

    def names(r):
        v = {r.drug, r.fda_ingredient, r.cmap_name}
        if r.cmap_name in ann.index:
            v |= {ann.loc[r.cmap_name, "chembl_pref_name"], ann.loc[r.cmap_name, "hub_name"]}
        out = set()
        for n in v:
            if isinstance(n, str) and n.strip() and not o79.is_code(n):
                w = " ".join(t for t in re.split(r"[\s\-]+", n.lower()) if t.upper() not in p79.SALTS)
                if len(w) >= 4:
                    out.add(w)
        return sorted(out)
    x["names"] = [names(r) for r in x.itertuples(index=False)]
    return x


def is_listed(trial, names):
    text = " | ".join(trial["interventions"] + trial["other_names"]).lower().replace("-", " ")
    return any(re.search(r"\b" + re.escape(n) + r"\b", text) for n in names)


def abstracts(pmids):
    f = CACHE / "trial_outcomes_abstracts.json"
    have = json.load(open(f)) if f.exists() else {}
    need = [p for p in sorted(set(pmids)) if p not in have]
    for k in range(0, len(need), 150):
        b = get(f"{EUTILS}/efetch.fcgi?db=pubmed&retmode=xml&id=" + ",".join(need[k:k + 150]))
        time.sleep(0.4)
        if b is None:
            continue
        for art in ET.fromstring(b).iter("PubmedArticle"):
            pid = art.findtext(".//PMID")
            title = "".join(art.find(".//ArticleTitle").itertext()) if art.find(".//ArticleTitle") is not None else ""
            ab = " ".join(("".join(a.itertext()) if a.get("Label") is None else f'{a.get("Label")}: ' + "".join(a.itertext()))
                          for a in art.iter("AbstractText"))
            jr = art.findtext(".//Journal/ISOAbbreviation") or ""
            yr = art.findtext(".//PubDate/Year") or (art.findtext(".//PubDate/MedlineDate") or "")[:4]
            have[pid] = dict(title=title, abstract=ab, journal=jr, year=yr)
        json.dump(have, open(f, "w"))
    return have


def fetch():
    OUT.mkdir(parents=True, exist_ok=True)
    x = pairs()
    cf = CACHE / "trial_outcomes_ctgov.jsonl"
    cache = {}
    if cf.exists():
        for line in open(cf):
            r = json.loads(line)
            cache[tuple(r["k"])] = r["v"]
    q = x.drop_duplicates(["drug", "condition"])[["drug", "condition"]]
    print(f"{len(x)} pairs ({(x.level == 'disease').sum()} disease, {(x.level == 'subgroup').sum()} subgroup); {len(q)} searches", flush=True)
    with open(cf, "a") as fh:
        for i, r in enumerate(q.itertuples(index=False)):
            if (r.drug, r.condition) in cache:
                continue
            v = ctgov(r.drug, r.condition)
            if v is None:
                print("failed:", r.drug, r.condition, flush=True)
                continue
            cache[(r.drug, r.condition)] = v
            fh.write(json.dumps(dict(k=[r.drug, r.condition], v=v)) + "\n"); fh.flush()
            if i % 50 == 0:
                print(f"  {i}/{len(q)}", flush=True)
    rows, n_raw = [], 0
    for r in x.itertuples(index=False):
        for t in cache.get((r.drug, r.condition), []):
            n_raw += 1
            if not is_listed(t, r.names):
                continue
            rows.append(dict(unit=r.unit, level=r.level, parent_group=r.parent_group, fda_rank=r.fda_rank, drug=r.drug,
                             condition=r.condition, nct=t["nct"], title=t["title"], phase=t["phase"], status=t["status"],
                             why_stopped=t["why_stopped"], stop_class=stop_class(t["why_stopped"]) if t["status"] in STOPPED else "",
                             has_results=t["has_results"], enrollment=t["enrollment"], start=t["start"], completion=t["completion"],
                             interventions="; ".join(t["interventions"]),
                             result_pmids=";".join(p for p, ty in t["refs"] if ty in ("RESULT", "DERIVED"))))
    tr = pd.DataFrame(rows)
    tr.to_csv(OUT / "trials.tsv", sep="\t", index=False)
    print(f"trial hits {n_raw}; drug listed as intervention {len(tr)} ({tr.nct.nunique()} trials) in {tr.groupby(['unit','drug']).ngroups} pairs")
    print(tr.drop_duplicates("nct").status.value_counts().to_dict())
    st = tr[tr.status.isin(STOPPED)].drop_duplicates("nct")
    print("stopped trials by class:", st.stop_class.value_counts().to_dict())
    comp = tr[(tr.status == "COMPLETED") & tr.result_pmids.fillna("").ne("")].drop_duplicates("nct")
    ab = abstracts([p for v in comp.result_pmids for p in v.split(";")])
    with open(OUT / "reading_queue.jsonl", "w") as fh:
        n = 0
        for t in comp.itertuples(index=False):
            papers = [dict(pmid=p, **ab[p]) for p in t.result_pmids.split(";") if p in ab and ab[p]["abstract"].strip()]
            if papers:
                fh.write(json.dumps(dict(nct=t.nct, title=t.title, drug=t.drug, condition=t.condition, phase=t.phase,
                                         interventions=t.interventions, papers=papers)) + "\n")
                n += 1
    print(f"completed trials {tr[tr.status == 'COMPLETED'].nct.nunique()}; with a linked result paper {len(comp)}; queued with abstracts {n}")


def norm(s):
    return re.sub(r"\s+", " ", str(s)).strip().lower()


def merge():
    tr = pd.read_csv(OUT / "trials.tsv", sep="\t", dtype={"result_pmids": str})
    q = {json.loads(l)["nct"]: json.loads(l) for l in open(OUT / "reading_queue.jsonl")}
    ans = {}
    for f in sorted((OUT / "readings").glob("*.jsonl")):
        for line in open(f):
            if line.strip():
                a = json.loads(line)
                ans[a["nct"]] = a
    acc = rej = 0
    rec = {}
    for nct, a in ans.items():
        if nct not in q:
            continue
        text = norm(" ".join(p["title"] + " " + p["abstract"] for p in q[nct]["papers"]))
        ok = bool(a.get("quote")) and norm(a["quote"]) in text
        if a.get("outcome") in ("not_reported", "not_this_trial"):
            ok = True
        acc += ok; rej += not ok
        rec[nct] = dict(model_drug_role=a.get("drug_role", ""), model_outcome=a.get("outcome") if ok else "not_reported", model_reason_class=a.get("reason_class") if ok else "",
                        model_reason=a.get("reason") if ok else "", model_quote=a.get("quote") if ok else "",
                        model_pmid=a.get("pmid") if ok else "", quote_check="ok" if ok else "rejected")
    m = tr.join(pd.DataFrame.from_dict(rec, orient="index"), on="nct")
    m.to_csv(OUT / "trials_read.tsv", sep="\t", index=False)
    print(f"answers {len(ans)} of {len(q)} queued; accepted {acc}, rejected (quote not in abstract) {rej}")
    # a negative trial in which the drug is part of the backbone in every arm does not test it (reported, not counted)
    role = OUT / "roles"   # role-only pass (D108 amendment): overrides drug_role
    if role.exists():
        rr = {}
        for f in sorted(role.glob("*.jsonl")):
            for line in open(f):
                if line.strip():
                    a = json.loads(line); rr[a["nct"]] = a["drug_role"]
        m["model_drug_role"] = m.nct.map(rr).fillna(m.model_drug_role)
    m["background"] = ~m.model_drug_role.isin(["single_agent", "combination"]) & m.model_drug_role.notna() & m.model_drug_role.ne("")
    # third amendment: stopped efficacy / safety trials get a role per (nct, drug) from roles/stopped.jsonl
    sf = OUT / "roles/stopped.jsonl"
    stop_role = {}
    if sf.exists():
        for line in open(sf):
            if line.strip():
                a = json.loads(line); stop_role[(a["nct"], a["drug"])] = (a["drug_role"], a.get("stop_attributable_to_drug", "yes"))
    is_stop = m.stop_class.isin(["efficacy", "safety"])
    m.loc[is_stop, "model_drug_role"] = [stop_role.get((n, d), ("", ""))[0] for n, d in zip(m.nct[is_stop], m.drug[is_stop])]
    m["stop_attributable"] = ""
    m.loc[is_stop, "stop_attributable"] = [stop_role.get((n, d), ("", ""))[1] for n, d in zip(m.nct[is_stop], m.drug[is_stop])]
    is_stop = is_stop & m.stop_attributable.ne("no")   # fourth amendment: a stop blamed on another drug does not count
    m["background"] = ~m.model_drug_role.isin(["single_agent", "combination"]) & m.model_drug_role.notna() & m.model_drug_role.ne("")
    m["failed"] = (is_stop | m.model_outcome.eq("negative")) & ~m.background
    m["fail_eff"] = m.failed & (m.stop_class.eq("efficacy") | m.model_reason_class.isin(["lack_of_efficacy", "both"]))
    m["fail_saf"] = m.failed & (m.stop_class.eq("safety") | m.model_reason_class.isin(["toxicity", "both"]))
    m["stopped_other"] = m.stop_class.isin(["accrual", "business", "other", "not_given"])
    m.to_csv(OUT / "trials_read.tsv", sep="\t", index=False)

    def summ(d):
        f = d[d.failed]
        why = [f"{r.nct} ({r.phase}, {r.status.lower()}; {r.model_drug_role or 'role not judged'}): " + (f"stopped for {r.stop_class}: {r.why_stopped}" if r.stop_class in ("efficacy", "safety")
               else f"negative [{r.model_reason_class}]: {r.model_reason} \"{r.model_quote}\" (PMID {r.model_pmid}; model-extracted)")
               for r in f.itertuples(index=False)]
        return pd.Series(dict(n_trials=d.nct.nunique(), n_completed=int((d.status == "COMPLETED").sum()),
                              n_positive=int((d.model_outcome.eq("positive") & ~d.background).sum()), n_mixed=int((d.model_outcome.eq("mixed") & ~d.background).sum()),
                              n_failed=int(d.failed.sum()), n_failed_efficacy=int(d.fail_eff.sum()), n_failed_safety=int(d.fail_saf.sum()),
                              n_not_testing_drug=int(d.background.sum()), n_stopped_not_failed=int(d.stopped_other.sum()), failures=" || ".join(why)))
    s = m.groupby(["unit", "level", "parent_group", "fda_rank", "drug"], dropna=False).apply(summ).reset_index()
    s.to_csv(RES / "trial_outcomes_v2.tsv", sep="\t", index=False)
    print(f"pairs with a listed trial {len(s)}; with a failed trial {(s.n_failed > 0).sum()} "
          f"(efficacy {(s.n_failed_efficacy > 0).sum()}, safety {(s.n_failed_safety > 0).sum()}); with a positive completed trial {(s.n_positive > 0).sum()}")


if __name__ == "__main__":
    {"fetch": fetch, "merge": merge}[STAGE]()
