#!/usr/bin/env python
"""D123: free fraction in plasma (fu) and free-drug occupancy (rule: docs/DECISIONS.md D123). Runs in `celligner`.

STAGE=fetch  openFDA label text for the 307 layered-score FDA ingredients (cache data/interim/lincs/d123_labels.json),
             protein-binding passages (+ Cmax passages for drugs without a D109 Cmax), fallback D109b document texts;
             writes results/lincs/d123/queue.jsonl and chunk_{1..N}.jsonl for the subagent readers.
STAGE=merge  reads results/lincs/d123/readings/chunk_*.jsonl, checks quotes, writes free_fraction.tsv (fu, Cmax) and
             occupancy.tsv (per drug x unit: total / free occupancy; drug level = strongest-binding unit).
"""
import hashlib
import json
import os
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results/lincs"
OUT = RES / "d123"
CACHE = ROOT / "data/interim/lincs"
DOCS = CACHE / "d109_docs"
N_CHUNKS = int(os.environ.get("N_CHUNKS", 4))
BIND = re.compile(r"protein[- ]bind|bound to (?:human )?(?:plasma|serum)|(?:plasma|serum) protein|unbound|free fraction|fraction unbound", re.I)
CMAX = re.compile(r"\bC ?max\b|peak (?:plasma |serum )?(?:concentration|level)", re.I)
SECTIONS = ("pharmacokinetics", "clinical_pharmacology", "description")


def get(url):
    for _ in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "repurposing-research"}), timeout=60) as r:
                return r.read().decode()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return ""
            time.sleep(2)
        except Exception:
            time.sleep(2)
    return ""


def sentences(t):
    return [s for s in re.split(r"(?<=[.;])\s+(?=[A-Z(])", re.sub(r"\s+", " ", t)) if s]


def passages(text, pat, width=1, maxn=8):
    ss = sentences(text)
    hit = [i for i, s in enumerate(ss) if pat.search(s)]
    out, used = [], set()
    for i in hit:
        rng = [j for j in range(max(0, i - width), min(len(ss), i + width + 1)) if j not in used]
        if rng:
            used.update(rng)
            out.append(" ".join(ss[j] for j in rng))
        if len(out) >= maxn:
            break
    return out


def drugs():
    D = pd.read_csv(RES / "target_role/d121/drug_decoded.tsv", sep="\t").drop_duplicates("fda_ingredient")[["drug", "fda_ingredient"]]
    c = pd.read_csv(RES / "d109/cmax.tsv", sep="\t")
    have = set(c.dropna(subset=["cmax_uM"]).drug.str.lower())
    D["need_cmax"] = ~D.drug.str.lower().isin(have)
    return D.reset_index(drop=True)


def fetch():
    OUT.mkdir(parents=True, exist_ok=True)
    D = drugs()
    cf = CACHE / "d123_labels.json"
    lab = json.load(open(cf)) if cf.exists() else {}
    for r in D.itertuples(index=False):
        if r.drug in lab:
            continue
        cands = []
        for field in ("openfda.generic_name", "openfda.substance_name"):
            for nm in dict.fromkeys([r.drug, r.fda_ingredient]):
                b = get(f"https://api.fda.gov/drug/label.json?limit=10&search={field}:%22{urllib.parse.quote(nm)}%22")
                time.sleep(0.3)
                if b:
                    for x in json.loads(b).get("results", []):
                        cands.append("\n".join(" ".join(x.get(k, [])) for k in SECTIONS))
            if cands:
                break
        withb = [t for t in cands if BIND.search(t)]
        lab[r.drug] = max(withb or cands or [""], key=len)
        json.dump(lab, open(cf, "w"))
    dq = [json.loads(l) for l in open(RES / "d109/doc_queue.jsonl")]
    urls = {}
    for q in dq:
        urls.setdefault(q["drug"], []).append(q["url"])
    rows = []
    for r in D.itertuples(index=False):
        t = lab.get(r.drug, "")
        pb, src, url = passages(t, BIND), "openfda_label", ""
        if not pb:
            for u in urls.get(r.drug, []):
                tf = DOCS / (hashlib.md5(u.encode()).hexdigest() + ".txt")
                if tf.exists():
                    pb = passages(tf.read_text(errors="ignore"), BIND)
                    if pb:
                        src, url = "fda_document", u
                        break
        pc = passages(t, CMAX, maxn=6) if r.need_cmax else []
        if pb or pc:
            rows.append(dict(drug=r.drug, fda_ingredient=r.fda_ingredient, source=src, url=url, need_cmax=bool(r.need_cmax),
                             binding_passages=pb, cmax_passages=pc))
    with open(OUT / "queue.jsonl", "w") as fh:
        for x in rows:
            fh.write(json.dumps(x) + "\n")
    for k in range(N_CHUNKS):
        with open(OUT / f"chunk_{k + 1}.jsonl", "w") as fh:
            for x in rows[k::N_CHUNKS]:
                fh.write(json.dumps(x) + "\n")
    print(f"drugs {len(D)}; label text {sum(bool(lab.get(d)) for d in D.drug)}; with binding passages "
          f"{sum(bool(x['binding_passages']) for x in rows)} (documents {sum(x['source'] == 'fda_document' for x in rows)}); "
          f"need Cmax {int(D.need_cmax.sum())}, with Cmax passages {sum(bool(x['cmax_passages']) for x in rows)}; queued {len(rows)}")


def norm(t):
    return re.sub(r"\s+", " ", str(t)).strip().lower()


UNIT = {"ng/ml": 1e-3, "µg/ml": 1.0, "mg/l": 1.0, "mg/ml": 1e3, "pg/ml": 1e-6, "ng/l": 1e-6,
        "mcg/l": 1e-3, "ug/l": 1e-3, "µg/l": 1e-3}


def merge():
    Q = {json.loads(l)["drug"]: json.loads(l) for l in open(OUT / "queue.jsonl")}
    A = [json.loads(l) for f in sorted((OUT / "readings").glob("chunk_*.jsonl")) for l in open(f) if l.strip()]
    mwf = CACHE / "d109_mw.json"
    mw = json.load(open(mwf)) if mwf.exists() else {}
    c109 = pd.read_csv(RES / "d109/cmax.tsv", sep="\t")
    mw_by_drug = {k: v for k, v in zip(c109.drug, c109.mw) if v == v}
    # drug -> ChEMBL ids of its LINCS compounds (fda_matches + annotations), parent molecular weight as D111 (cached)
    fmx = pd.read_csv(RES / "fda_matches.tsv", sep="\t")
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t", low_memory=False)[["cmap_name", "chembl_id"]].dropna()
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from fda_rows import drug_name
    fmx["drug"] = fmx.fda_ingredient.map(drug_name)
    cids = fmx.merge(ann, on="cmap_name").groupby("drug").chembl_id.apply(lambda v: sorted(set(v))).to_dict()

    def molw(cid):
        if cid not in mw:
            b = get(f"https://www.ebi.ac.uk/chembl/api/data/molecule/{cid}.json")
            if b:
                j = json.loads(b)
                par = (j.get("molecule_hierarchy") or {}).get("parent_chembl_id") or cid
                if par != cid:
                    b2 = get(f"https://www.ebi.ac.uk/chembl/api/data/molecule/{par}.json")
                    j = json.loads(b2) if b2 else j
                mw[cid] = float((j.get("molecule_properties") or {}).get("full_mwt") or "nan")
                json.dump(mw, open(mwf, "w"))
        return mw.get(cid, float("nan"))

    def mw_of(drug):
        if drug in mw_by_drug:
            return mw_by_drug[drug]
        for c in cids.get(drug, []):
            v = molw(c)
            if v == v:
                return v
        return float("nan")
    rows = []
    for a in A:
        q = Q.get(a["drug"], {})
        tb, tc = norm(" ".join(q.get("binding_passages", []))), norm(" ".join(q.get("cmax_passages", [])))
        fu = a.get("fu")
        fu_ok = fu not in (None, "") and bool(a.get("fu_quote")) and norm(a["fu_quote"]) in tb
        cm = a.get("cmax_value")
        cm_ok = cm not in (None, "") and bool(a.get("cmax_quote")) and norm(a["cmax_quote"]) in tc
        um = np.nan
        if cm_ok:
            u = norm(a.get("cmax_unit", "")).replace("μ", "µ").replace("mcg", "µg").replace("ug/", "µg/")
            m_ = mw_of(a["drug"])
            if u in ("um", "µm"):
                um = float(cm)
            elif u == "nm":
                um = float(cm) / 1e3
            elif u in UNIT and m_ and m_ == m_:
                um = float(cm) * UNIT[u] / float(m_) * 1e3
        rows.append(dict(drug=a["drug"], fu=float(fu) if fu_ok else np.nan, fu_text=a.get("fu_text", ""), fu_quote=a.get("fu_quote", ""),
                         fu_flags=a.get("fu_flags", ""), fu_check="ok" if fu_ok else ("rejected" if fu not in (None, "") else "none"),
                         cmax_value=cm, cmax_unit=a.get("cmax_unit", ""), cmax_dose=a.get("cmax_dose", ""), cmax_quote=a.get("cmax_quote", ""),
                         cmax_check="ok" if cm_ok else ("rejected" if cm not in (None, "") else "none"), cmax_uM_d123=um,
                         source=q.get("source", ""), url=q.get("url", ""), note=a.get("note", "")))
    F = pd.DataFrame(rows)
    D = drugs()
    F = D.merge(F, on="drug", how="left")
    old = c109.dropna(subset=["cmax_uM"]).drop_duplicates("drug").set_index("drug").cmax_uM
    F["cmax_uM"] = F.drug.map(old).fillna(F.cmax_uM_d123)
    F["cmax_source"] = np.where(F.drug.isin(old.index), "d109", np.where(F.cmax_uM_d123.notna(), "d123", ""))
    F.to_csv(OUT / "free_fraction.tsv", sep="\t", index=False)
    # occupancy per drug x unit (units from D121)
    U = pd.read_csv(RES / "target_role/d121/units.tsv", sep="\t")
    dd = pd.read_csv(RES / "target_role/d121/drug_decoded.tsv", sep="\t").drop_duplicates("fda_ingredient")
    fm = pd.read_csv(RES / "fda_matches.tsv", sep="\t").drop_duplicates("cmap_name")
    U["fda_ingredient"] = U.cmap_name.map(dict(zip(fm.cmap_name, fm.fda_ingredient)))
    U["drug"] = U.fda_ingredient.map(dict(zip(dd.fda_ingredient, dd.drug)))
    U = U.merge(F[["drug", "fu", "cmax_uM"]], on="drug", how="left")
    K = 10 ** (-U.p_u) * 1e6   # µM
    U["occ_total"] = U.cmax_uM / (U.cmax_uM + K)
    U["occ_free"] = U.cmax_uM * U.fu / (U.cmax_uM * U.fu + K)
    U.to_csv(OUT / "occupancy_units.tsv", sep="\t", index=False, float_format="%.4g")
    best = U.dropna(subset=["p_u"]).sort_values("p_u", ascending=False).drop_duplicates("drug")
    O = best[["drug", "fda_ingredient", "genes", "p_u", "fu", "cmax_uM", "occ_total", "occ_free"]].copy()
    O["low_free_occupancy"] = O.occ_free < 0.5
    O.to_csv(OUT / "occupancy.tsv", sep="\t", index=False, float_format="%.4g")
    print(f"drugs {len(D)}; fu read {int(F.fu.notna().sum())} (rejected quotes {int((F.fu_check == 'rejected').sum())}); "
          f"Cmax: D109 {int((F.cmax_source == 'd109').sum())}, new D123 {int((F.cmax_source == 'd123').sum())} "
          f"(rejected {int((F.cmax_check == 'rejected').sum())}); free occupancy for {int(O.occ_free.notna().sum())} drugs, "
          f"low_free_occupancy {int(O.low_free_occupancy.sum())}; total occupancy for {int(O.occ_total.notna().sum())}")


if __name__ == "__main__":
    {"fetch": fetch, "merge": merge}[os.environ.get("STAGE", "fetch")]()
