#!/usr/bin/env python
"""Evidence flags (D109, user 2026-11-06): dose_inconsistent, inactive_matched_line, generic_mechanism, low_exposure.

Signatures: LINCS 2020 siginfo_beta (all, any quality); active = is_hiq and TAS >= 0.4 (as 79c). Lines behind a score
as 79e: direct = the drug's active lines in the cohort's tissue key (in_key in the Celligner weights), borrowed = its
active lines outside the key. Matched line = in-key line with weight >= 0.1 x the cohort's top in-key weight (D76).
(1) dose_inconsistent: per drug x line and time point with an active signature, signatures at doses >= 2x the highest
    active dose (any quality) have median TAS < 0.5 x the median TAS of the active ones; profiles without a higher dose
    are not judged. Row flag: >= 50% of the judged lines behind the score inconsistent.
(2) inactive_matched_line: borrowed row; the drug has >= 1 is_hiq signature in a matched line, none active. D116: over the
    drug's sibling compounds (same FDA ingredient, else cmap_name); per-compound value kept as inactive_matched_line_own.
    inconclusive_matched_line: only low-quality signatures there (reported, not a flag).
(3) generic_mechanism: ChEMBL curated mechanism (molecule or parent) in GENERIC.
(4) low_exposure (FDA drugs, STAGE=exposure): lowest active dose in the lines behind the score > 10 x label Cmax (uM).
STAGE=flags   -> results/lincs/evidence_flags_v2.tsv (groups and subgroups; low_exposure filled by STAGE=exposure),
                 results/lincs/d109/profile_dose.tsv, and the label queue results/lincs/d109/label_queue.jsonl
                 (FDA drugs of the 79p eligible set: openFDA label text, cached in data/interim/lincs/d109_labels.json)
STAGE=exposure -> reads results/lincs/d109/cmax/*.jsonl ({drug, cmax_value, cmax_unit, dose, quote}), checks each quote is
                 verbatim in the label text, converts to uM with the ChEMBL parent full_mwt, fills low_exposure.
"""
import json
import os
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
OUT = RES / "d109"
CACHE = ROOT / "data/interim/lincs"
STAGE = os.environ.get("STAGE", "flags")
TAS_MIN, HI_DOSE, TAS_RATIO, ROW_SHARE, LM_MIN, EXPO = 0.4, 2.0, 0.5, 0.5, 0.1, 10.0
GENERIC = {"80S Ribosome inhibitor", "Protein synthesis inhibitor", "Sodium/potassium-transporting ATPase inhibitor",
           "DNA-directed RNA polymerase inhibitor"}


def get(url):
    for k in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "repurposing-pipeline"}), timeout=60) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            time.sleep(2 * (k + 1))
        except Exception:
            time.sleep(2 * (k + 1))
    return None


def load_sigs():
    s = pd.read_csv(ROOT / "data/raw/lincs2020/siginfo_beta.txt", sep="\t", low_memory=False,
                    usecols=["pert_id", "cell_iname", "pert_dose", "pert_dose_unit", "pert_time", "tas", "is_hiq", "pert_type"])
    s = s[s.pert_type == "trt_cp"].copy()
    s["dose"] = pd.to_numeric(s.pert_dose, errors="coerce")
    s.loc[s.pert_dose_unit.astype(str).str.lower().isin(["nm"]), "dose"] /= 1000
    s = s[s.pert_dose_unit.astype(str).str.lower().isin(["um", "nm"]) & s.dose.notna()]
    s["active"] = (s.is_hiq == 1) & (s.tas >= TAS_MIN)
    cl = pd.read_csv(ROOT / "data/processed/lincs/compound_line.tsv.gz", sep="\t").dropna(subset=["ModelID"])
    m = cl.drop_duplicates("cell_iname").set_index("cell_iname").ModelID
    s["ModelID"] = s.cell_iname.map(m)
    return s.dropna(subset=["ModelID"])


def profile_dose(s):
    """per drug x line: judged (bool), inconsistent (bool), lowest active dose"""
    out = []
    act = s[s.active]
    keys = set(zip(act.pert_id, act.ModelID))
    sub = s[[k in keys for k in zip(s.pert_id, s.ModelID)]]
    for (p, mdl), d in sub.groupby(["pert_id", "ModelID"]):
        judged = incons = False
        for t, dt in d.groupby("pert_time"):
            a = dt[dt.active]
            if a.empty:
                continue
            hi = dt[dt.dose >= HI_DOSE * a.dose.max()]
            if hi.empty:
                continue
            judged = True
            if hi.tas.median() < TAS_RATIO * a.tas.median():
                incons = True
        out.append(dict(pert_id=p, ModelID=mdl, judged=judged, inconsistent=incons, min_active_dose=d[d.active].dose.min()))
    return pd.DataFrame(out)


def flags():
    OUT.mkdir(parents=True, exist_ok=True)
    s = load_sigs()
    pdz = profile_dose(s)
    pdz.to_csv(OUT / "profile_dose.tsv", sep="\t", index=False)
    print(f"active drug x line profiles {len(pdz)}; judged {int(pdz.judged.sum())}; inconsistent {int(pdz.inconsistent.sum())}")
    pdz = pdz.set_index(["pert_id", "ModelID"])
    idx = pd.read_csv(ROOT / "data/processed/lincs/active_profiles_index.tsv.gz", sep="\t")
    lines_of = idx.groupby("pert_id").ModelID.apply(set).to_dict()
    hiq = s[s.is_hiq == 1].groupby("pert_id").ModelID.apply(set).to_dict()
    anyq = s.groupby("pert_id").ModelID.apply(set).to_dict()
    # generic mechanism via ChEMBL (molecule or parent id)
    mech = json.load(open(ROOT / "data/raw/chembl/ChEMBL_37/d105_all/mechanisms.json"))
    gen_ids = {i for x in mech if x.get("mechanism_of_action") in GENERIC for i in (x["molecule_chembl_id"], x.get("parent_molecule_chembl_id"))}
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t", low_memory=False).drop_duplicates("pert_id").set_index("pert_id")
    gen = ann.chembl_id.isin(gen_ids)
    # D116: sibling compounds of a drug (FDA ingredient for FDA drugs, else cmap_name)
    fmz = pd.read_csv(RES / "fda_matches.tsv", sep="\t").drop_duplicates("pert_id").set_index("pert_id").fda_ingredient
    dkey = {p: ("fda:" + fmz[p]) if p in fmz.index else ("cmap:" + str(c)) for p, c in ann.cmap_name.items()}
    sib = {}
    for p, k in dkey.items():
        sib.setdefault(k, set()).add(p)
    sibs = lambda p: sib.get(dkey.get(p), set()) | {p}
    wg = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    ws = pd.read_csv(ROOT / "results/celligner/stratum_line_weights_long_v2.tsv", sep="\t").rename(columns={"stratum": "group_id"})
    res = []
    for split, w in (("drug_reversal_split_v2_consensus.tsv", wg), ("stratum_reversal_split_v2_consensus.tsv", ws)):
        df = pd.read_csv(RES / split, sep="\t", low_memory=False, usecols=lambda c: c in ("group_id", "stratum", "pert_id", "source"))
        df = df.rename(columns={"stratum": "group_id"})
        for g, d in df.groupby("group_id"):
            wk = w[(w.group_id == g) & w.in_key.astype(str).str.lower().eq("true")]
            key = set(wk.ModelID)
            top = wk.weight.max() if len(wk) else 0
            matched = set(wk.ModelID[wk.weight >= LM_MIN * top]) if top > 0 else set()
            for p, src in zip(d.pert_id, d.source):
                L = lines_of.get(p, set())
                L = (L & key) if src == "direct" else (L - key)
                j = [pdz.loc[(p, m)] for m in L if (p, m) in pdz.index]
                jj = [r for r in j if r.judged]
                inc = sum(r.inconsistent for r in jj)
                min_dose = min([r.min_active_dose for r in j], default=np.nan)
                inact = incon = own = False
                if src == "borrowed":
                    mh = hiq.get(p, set()) & matched
                    own = bool(mh)                       # hiq signatures in a matched line, none active (else it would be direct)
                    S = sibs(p)                          # D116: judged over the drug's sibling compounds
                    act_m = any(lines_of.get(q, set()) & matched for q in S)
                    inact = (not act_m) and any(hiq.get(q, set()) & matched for q in S)
                    incon = (not act_m) and (not inact) and any(anyq.get(q, set()) & matched for q in S)
                res.append(dict(group_id=g, pert_id=p, source=src, n_lines=len(L), n_dose_judged=len(jj), n_dose_inconsistent=inc,
                                dose_inconsistent=bool(jj) and inc / len(jj) >= ROW_SHARE, inactive_matched_line=inact, inactive_matched_line_own=own,
                                inconclusive_matched_line=incon, generic_mechanism=bool(gen.get(p, False)),
                                min_active_dose_uM=min_dose))
    ef = pd.DataFrame(res)
    ef["low_exposure"] = pd.NA
    ef.to_csv(RES / "evidence_flags_v2.tsv", sep="\t", index=False)
    for c in ("dose_inconsistent", "inactive_matched_line", "inconclusive_matched_line", "generic_mechanism"):
        print(f"{c}: {int(ef[c].sum())} of {len(ef)} rows ({ef[ef[c]].pert_id.nunique()} drugs)")
    label_queue(ann)


def label_queue(ann):
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    import evidence_flags; evidence_flags.STRICT = False   # D117: this script builds the flags
    from fda_rows import fda_rows
    g, s = fda_rows()
    drugs = pd.concat([g, s]).drop_duplicates("fda_ingredient")[["drug", "fda_ingredient", "pert_id"]]
    cf = CACHE / "d109_labels.json"
    lab = json.load(open(cf)) if cf.exists() else {}
    for r in drugs.itertuples(index=False):
        if r.drug in lab:
            continue
        txt = ""
        for field in ("openfda.generic_name", "openfda.substance_name"):
            for nm in dict.fromkeys([r.drug, r.fda_ingredient]):
                b = get(f"https://api.fda.gov/drug/label.json?limit=5&search={field}:%22{urllib.parse.quote(nm)}%22")
                time.sleep(0.3)
                if not b:
                    continue
                best = ""
                for x in json.loads(b).get("results", []):
                    t = "\n".join(" ".join(x.get(k, [])) for k in ("pharmacokinetics", "clinical_pharmacology", "dosage_and_administration"))
                    if "cmax" in t.lower() or "peak" in t.lower():
                        best = max(best, t, key=len) if best else t
                if best:
                    txt = best
                    break
            if txt:
                break
        lab[r.drug] = txt
        json.dump(lab, open(cf, "w"))
    with open(OUT / "label_queue.jsonl", "w") as fh:
        for r in drugs.itertuples(index=False):
            if lab.get(r.drug):
                fh.write(json.dumps(dict(drug=r.drug, label_text=lab[r.drug])) + "\n")
    print(f"FDA drugs {len(drugs)}; label text with Cmax / peak found {sum(bool(lab.get(d)) for d in drugs.drug)}")


def norm(t):
    return re.sub(r"\s+", " ", str(t)).strip().lower()


UNIT = {"ng/ml": 1e-3, "mcg/ml": 1.0, "ug/ml": 1.0, "µg/ml": 1.0, "mg/l": 1.0, "mg/ml": 1e3, "pg/ml": 1e-6, "ng/l": 1e-6,
        "um": None, "µm": None, "nm": None, "mcg/l": 1e-3, "ug/l": 1e-3, "µg/l": 1e-3}


def exposure():
    ef = pd.read_csv(RES / "evidence_flags_v2.tsv", sep="\t")
    lab = json.load(open(CACHE / "d109_labels.json"))
    ans = {}
    for f in sorted((OUT / "cmax").glob("*.jsonl")):
        for line in open(f):
            if line.strip():
                a = json.loads(line)
                ans[a["drug"]] = a
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    import evidence_flags; evidence_flags.STRICT = False   # D117: this script builds the flags
    from fda_rows import fda_rows, drug_name
    g, s = fda_rows()
    fd = pd.concat([g, s])[["group_id", "pert_id", "drug"]]
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t", low_memory=False).drop_duplicates("pert_id").set_index("pert_id")
    # D111: Cmax belongs to the ingredient -> every LINCS compound of an FDA ingredient (fda_matches, same drug name)
    fm = pd.read_csv(RES / "fda_matches.tsv", sep="\t").drop_duplicates("pert_id")
    p2d = dict(zip(fm.pert_id, fm.fda_ingredient.map(drug_name)))
    p2c = dict(zip(fm.pert_id, fm.cmap_name))
    nl = json.load(open(ROOT / "data/raw/chembl/ChEMBL_37/d105_all/name_lookup.json"))   # cmap name -> ChEMBL id (D105 cache)
    d2p = {}
    for p_, d_ in sorted(p2d.items()):
        d2p.setdefault(d_, []).append(p_)
    mwf = CACHE / "d109_mw.json"
    mw = json.load(open(mwf)) if mwf.exists() else {}
    def molw(cid):
        if isinstance(cid, str) and cid not in mw:
            b = get(f"https://www.ebi.ac.uk/chembl/api/data/molecule/{cid}.json")
            if b:
                j = json.loads(b)
                par = (j.get("molecule_hierarchy") or {}).get("parent_chembl_id") or cid
                if par != cid:
                    b2 = get(f"https://www.ebi.ac.uk/chembl/api/data/molecule/{par}.json")
                    j = json.loads(b2) if b2 else j
                mw[cid] = float((j.get("molecule_properties") or {}).get("full_mwt") or "nan")
                json.dump(mw, open(mwf, "w"))
        return mw.get(cid, np.nan) if isinstance(cid, str) else np.nan

    def judge(drug, a, text, source):
        ok = bool(a.get("quote")) and norm(a["quote"]) in norm(text) and a.get("cmax_value") not in (None, "")
        # D111: weight from any compound of the ingredient (fda_rows compound first, then the others in pert_id order)
        pids = list(dict.fromkeys(list(fd[fd.drug == drug].pert_id[:1]) + d2p.get(drug, [])))
        m_ = np.nan
        if ok:
            for cid in [ann.chembl_id.get(pid) for pid in pids] + [nl.get(p2c.get(pid)) for pid in pids]:   # D111: annotation id, else ChEMBL name lookup
                m_ = molw(cid)
                if m_ == m_:
                    break
        u = str(a.get("cmax_unit", "")).lower()
        u = re.sub(r"\bmicrograms?\b", "ug", re.sub(r"\bnanograms?\b", "ng", re.sub(r"\s+per\s+", "/", u)))   # D111 spellings
        u = re.sub(r"\s+", "", u).replace("μ", "µ").replace("mcmol", "µmol").replace("umol/l", "um").replace("µmol/l", "µm").replace("nmol/l", "nm")
        v = pd.to_numeric(str(a.get("cmax_value", "")).replace(",", ""), errors="coerce")
        um = np.nan
        if ok and v == v:
            if u in ("um", "µm"):
                um = v
            elif u == "nm":
                um = v / 1000
            elif u in UNIT and UNIT[u] is not None and m_ == m_:
                um = v * UNIT[u] * 1000 / m_       # ug/mL -> mg/L -> umol/L
        return dict(drug=drug, source=source, url=a.get("url", ""), cmax_value=a.get("cmax_value"), cmax_unit=a.get("cmax_unit"),
                    dose=a.get("dose"), quote=a.get("quote"), note=a.get("note", ""), quote_check="ok" if ok else "rejected", mw=m_, cmax_uM=um)

    rows = [judge(drug, a, lab.get(drug, ""), "openfda_label") for drug, a in ans.items()]
    # D109b: archived labels and reviews (precedence label > review), quote checked against the document text
    have = {r["drug"] for r in rows if r["cmax_uM"] == r["cmax_uM"]}
    dq = OUT / "doc_queue.jsonl"
    dtext = {}
    if dq.exists():
        for line in open(dq):
            q = json.loads(line)
            tf = CACHE / "d109_docs" / (__import__("hashlib").md5(q["url"].encode()).hexdigest() + ".txt")
            dtext[(q["drug"], q["url"])] = (q["doc_type"], tf.read_text(errors="ignore") if tf.exists() else " ".join(q["passages"]))
    dans = []
    for f in sorted((OUT / "cmax_docs").glob("*.jsonl")) if (OUT / "cmax_docs").exists() else []:
        for line in open(f):
            if line.strip():
                dans.append(json.loads(line))
    for typ in ("label", "review"):
        for a in dans:
            k = (a["drug"], a.get("url", ""))
            if a["drug"] in have or k not in dtext or dtext[k][0] != typ:
                continue
            r = judge(a["drug"], a, dtext[k][1], "fda_" + typ)
            rows.append(r)
            if r["cmax_uM"] == r["cmax_uM"]:
                have.add(a["drug"])
    cm = pd.DataFrame(rows)
    cm.to_csv(OUT / "cmax.tsv", sep="\t", index=False)
    print(f"Cmax answers {len(cm)}; quote ok {int((cm.quote_check == 'ok').sum())}; drugs with a Cmax in uM {cm[cm.cmax_uM.notna()].drug.nunique()} "
          f"(by source: {cm[cm.cmax_uM.notna()].drop_duplicates('drug').source.value_counts().to_dict()})")
    cm.to_csv(OUT / "cmax.tsv", sep="\t", index=False)
    best = cm[cm.cmax_uM.notna()].drop_duplicates("drug")[["drug", "cmax_uM", "source"]].rename(columns={"source": "cmax_source"})
    ef = ef.drop(columns=[c for c in ("low_exposure", "drug", "cmax_uM", "cmax_source", "exposure_ratio") if c in ef])
    ef["drug"] = ef.pert_id.map(p2d)   # D111: every compound of an FDA ingredient, not only the fda_rows pick
    ef = ef.merge(best, on="drug", how="left")
    ef["exposure_ratio"] = ef.min_active_dose_uM / ef.cmax_uM
    ef["low_exposure"] = np.where(ef.exposure_ratio.notna(), ef.exposure_ratio > EXPO, pd.NA)
    ef.to_csv(RES / "evidence_flags_v2.tsv", sep="\t", index=False)
    x = ef[ef.drug.notna()]
    print(f"FDA rows {len(x)}; exposure assessed {int(x.exposure_ratio.notna().sum())}; low_exposure {int((x.low_exposure == True).sum())} "
          f"({x[x.low_exposure == True].drug.nunique()} drugs: {sorted(x[x.low_exposure == True].drug.unique())[:40]})")


if __name__ == "__main__":
    {"flags": flags, "exposure": exposure}[STAGE]()
