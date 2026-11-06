#!/usr/bin/env python
"""D109b: Cmax passages from Drugs@FDA archived labels and reviews for FDA drugs without a verified label Cmax.

Drugs: 79p eligible set (fda_ranking_v2_consensus / strata) minus drugs with cmax_uM in results/lincs/d109/cmax.tsv.
Documents: ApplicationDocs of the drug's NDA / BLA applications (fda_matches.tsv fda_appl_nda_bla), type 2 (Label,
most recent) and 3 (Review; table-of-contents pages followed to clinical pharmacology / biopharmaceutics PDFs).
Text by pdftotext, cached in data/interim/lincs/d109_docs/. Queue: results/lincs/d109/doc_queue.jsonl
({drug, doc_type, url, passages}); passages = +-700 characters around Cmax / peak / maximum concentration mentions.
"""
import hashlib
import json
import re
import subprocess
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
FDA = ROOT / "data/raw/drugsatfda/extracted"
CACHE = ROOT / "data/interim/lincs/d109_docs"
CACHE.mkdir(parents=True, exist_ok=True)
PAT = re.compile(r"C\s?max|peak (plasma|serum|blood)? ?(concentration|level)|maximum (plasma |serum )?concentration", re.I)


def fetch(url):
    h = hashlib.md5(url.encode()).hexdigest()
    f = CACHE / h
    if f.exists():
        return f.read_bytes()
    for k in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (research pipeline)"}), timeout=90) as r:
                b = r.read()
            f.write_bytes(b)
            time.sleep(0.5)
            return b
        except Exception:
            time.sleep(3 * (k + 1))
    return None


def text_of(url):
    tf = CACHE / (hashlib.md5(url.encode()).hexdigest() + ".txt")
    if tf.exists():
        return tf.read_text()
    b = fetch(url)
    if not b or not b[:5] == b"%PDF-":
        return ""
    pf = CACHE / (hashlib.md5(url.encode()).hexdigest() + ".pdf")
    pf.write_bytes(b)
    subprocess.run(["pdftotext", "-layout", str(pf), str(tf)], capture_output=True, timeout=300)
    return tf.read_text(errors="ignore") if tf.exists() else ""


def review_pdfs(url):
    if url.lower().endswith(".pdf"):
        return [url]
    b = fetch(url)
    if not b:
        return []
    html = b.decode("latin-1")
    links = [urllib.parse.urljoin(url, h) for h in re.findall(r'href="([^"]+\.pdf)"', html, re.I)]
    keep = [l for l in links if re.search(r"clinpharm|biopharm|ClinPharmR|BioPharmR|pharmr", l, re.I)]
    return keep or links[:3]


def passages(t):
    out, last = [], -1
    for m in PAT.finditer(t):
        a, b = max(0, m.start() - 700), m.end() + 700
        if a <= last:
            continue
        seg = t[a:b]
        if re.search(r"\d", seg):
            out.append(seg); last = b
        if len(out) >= 12:
            break
    return out


def main():
    g = pd.read_csv(RES / "fda_ranking_v2_consensus.tsv", sep="\t")
    s = pd.read_csv(RES / "fda_ranking_strata_v2_consensus.tsv", sep="\t")
    drugs = pd.concat([g, s]).drop_duplicates("drug")[["drug", "cmap_name", "fda_ingredient"]]
    cm = pd.read_csv(RES / "d109/cmax.tsv", sep="\t")
    have = set(cm.drug[cm.cmax_uM.notna()])
    drugs = drugs[~drugs.drug.isin(have)]
    fm = pd.read_csv(RES / "fda_matches.tsv", sep="\t").drop_duplicates("fda_ingredient").set_index("fda_ingredient").fda_appl_nda_bla
    docs = pd.read_csv(FDA / "ApplicationDocs.txt", sep="\t", dtype=str, encoding="latin-1", on_bad_lines="skip")
    print(f"drugs needing Cmax {len(drugs)}", flush=True)
    n = 0
    with open(RES / "d109/doc_queue.jsonl", "w") as fh:
        for r in drugs.itertuples(index=False):
            appl = [a for a in str(fm.get(r.fda_ingredient, "")).split(";") if a and a != "nan"]
            d = docs[docs.ApplNo.isin(appl)].sort_values("ApplicationDocsDate", ascending=False)
            cands = [("label", u) for u in d[d.ApplicationDocsTypeID == "2"].ApplicationDocsURL.head(2)]
            for u in d[d.ApplicationDocsTypeID == "3"].ApplicationDocsURL.tail(3):   # oldest = original review
                cands += [("review", p) for p in review_pdfs(u)[:4]]
            got = False
            for typ, u in cands:
                t = text_of(u)
                ps = passages(t)
                if ps:
                    fh.write(json.dumps(dict(drug=r.drug, doc_type=typ, url=u, passages=ps)) + "\n"); fh.flush()
                    got = True
                    if typ == "label":
                        break
            n += got
            print(f"  {r.drug}: {len(cands)} docs, passages {'yes' if got else 'no'}", flush=True)
    print(f"drugs with document passages {n} of {len(drugs)}")


if __name__ == "__main__":
    main()
