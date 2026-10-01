#!/usr/bin/env python
"""Step 1: disease grounding for the Toil (TCGA/TARGET) cohorts.

Reads the hand-curated map in config/toil_disease_curation.tsv, validates every
MONDO / UBERON / OncoTree / GTEx-tissue identifier against the pinned raw files,
and writes:
  results/step1/disease_grounding_toil.tsv   one row per cohort
  results/step1/mondo_descendants.tsv        long format: mondo_id, descendant_id, label
Exits non-zero if any identifier fails validation.
"""
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
OUT = ROOT / "results/step1"
XREF_PREFIXES = ["EFO", "DOID", "NCIT", "UMLS", "OMIM", "Orphanet", "ICDO"]
TUMOR_TYPES = re.compile(r"Tumor|Cancer|Metastatic|New Primary")


def parse_obo(path, prefix):
    """Return {id: dict(name, obsolete, parents, xrefs)} for terms with the given prefix."""
    terms = {}
    for block in path.read_text().split("\n[Term]\n")[1:]:
        block = block.split("\n[Typedef]")[0]
        m = re.search(r"^id: (\S+)", block, re.M)
        if not m or not m.group(1).startswith(prefix):
            continue
        terms[m.group(1)] = dict(
            name=re.search(r"^name: (.+)", block, re.M).group(1),
            obsolete="is_obsolete: true" in block,
            parents=re.findall(rf"^is_a: ({prefix}\S+)", block, re.M),
            xrefs=re.findall(r"^xref: (\S+)", block, re.M),
        )
    return terms


def descendants(terms):
    children = defaultdict(set)
    for tid, t in terms.items():
        if not t["obsolete"]:
            for p in t["parents"]:
                children[p].add(tid)

    def walk(tid):
        seen, stack = set(), [tid]
        while stack:
            for c in children[stack.pop()]:
                if c not in seen:
                    seen.add(c)
                    stack.append(c)
        return seen

    return walk


def main():
    cur = pd.read_csv(ROOT / "config/toil_disease_curation.tsv", sep="\t", dtype=str).fillna("")
    mondo = parse_obo(RAW / "mondo/mondo.obo", "MONDO:")
    uberon = parse_obo(RAW / "uberon/uberon-basic.obo", "UBERON:")
    oncotree = {x["code"]: x for x in json.loads((RAW / "oncotree/tumorTypes.json").read_text())}
    with gzip.open(RAW / "solid_toil/TcgaTargetGTEX_phenotype.txt.gz", "rt", encoding="latin-1") as fh:
        pheno = pd.read_csv(fh, sep="\t", dtype=str).fillna("")
    pheno.columns = ["sample", "category", "disease", "site", "sample_type", "sex", "study"]
    gtex = pheno[(pheno.study == "GTEX") & (pheno.sample_type == "Normal Tissue")]
    gtex_counts = gtex.category.value_counts()
    walk = descendants(mondo)

    errors, rows, desc_rows = [], [], []
    for r in cur.itertuples(index=False):
        m = mondo.get(r.mondo_id)
        if m is None or m["obsolete"]:
            errors.append(f"{r.cohort}: MONDO {r.mondo_id} missing or obsolete")
            continue
        u = uberon.get(r.uberon)
        if u is None or u["obsolete"]:
            errors.append(f"{r.cohort}: UBERON {r.uberon} missing or obsolete")
        codes = [c for c in r.oncotree.split(";") if c]
        for c in codes:
            if c not in oncotree:
                errors.append(f"{r.cohort}: OncoTree {c} not in pinned release")
        tissues = [t for t in r.gtex_tissues.split(";") if t]
        for t in tissues:
            if t not in gtex_counts:
                errors.append(f"{r.cohort}: GTEx tissue '{t}' not in Toil phenotype")

        coh = pheno[(pheno.study == r.toil_study) & (pheno.category == r.toil_category)]
        xrefs = {p: ";".join(sorted(x.split(":", 1)[1] for x in m["xrefs"] if x.startswith(p + ":")))
                 for p in XREF_PREFIXES}
        # cross-check: does the OncoTree code's NCIt match one of MONDO's NCIt xrefs?
        onc_ncit = {n for c in codes for n in oncotree.get(c, {}).get("externalReferences", {}).get("NCI", [])}
        ncit_agree = "n/a" if not onc_ncit else ("yes" if onc_ncit & set(xrefs["NCIT"].split(";")) else "no")
        desc = sorted(walk(r.mondo_id))
        desc_rows += [(r.mondo_id, d, mondo[d]["name"]) for d in desc]

        rows.append(dict(
            cohort=r.cohort,
            toil_study=r.toil_study,
            toil_category=r.toil_category,
            mondo_id=r.mondo_id,
            mondo_label=m["name"],
            n_mondo_descendants=len(desc),
            **{f"{p.lower()}_xref": v for p, v in xrefs.items()},
            oncotree=r.oncotree,
            oncotree_name=";".join(oncotree[c]["name"] for c in codes if c in oncotree),
            oncotree_ncit_matches_mondo=ncit_agree,
            uberon=r.uberon,
            uberon_label=u["name"] if u else "",
            gtex_tissues=r.gtex_tissues,
            gtex_match=r.gtex_match,
            depmap_lineage="",  # DepMap not pulled yet
            arm=r.arm,
            mixed_histology=r.mixed_histology,
            n_tumor=int(coh.sample_type.str.contains(TUMOR_TYPES).sum()),
            n_adjacent_normal=int((coh.sample_type == "Solid Tissue Normal").sum()),
            n_gtex_normal=int(sum(gtex_counts.get(t, 0) for t in tissues)),
            notes=r.notes,
        ))

    if errors:
        print("VALIDATION FAILED:\n  " + "\n  ".join(errors), file=sys.stderr)
        sys.exit(1)

    OUT.mkdir(parents=True, exist_ok=True)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "disease_grounding_toil.tsv", sep="\t", index=False)
    pd.DataFrame(desc_rows, columns=["mondo_id", "descendant_id", "descendant_label"]).to_csv(
        OUT / "mondo_descendants.tsv", sep="\t", index=False)
    print(f"{len(out)} cohorts validated; {out.n_tumor.sum()} tumors; "
          f"{(out.n_tumor > 0).sum()} cohorts with tumors")


if __name__ == "__main__":
    main()
