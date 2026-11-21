#!/usr/bin/env python
"""Fetch ChEMBL drug_indication: per-(molecule, indication) highest clinical phase reached (D178 input).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python   (run with -I: the output dir is untrusted)

WHY. D176 showed our ranking predicts how OFTEN a drug is trialled (11.73 / 3.65 / 1.71 per candidate
pair, p = 5.1e-12) but NOT whether those trials succeed (0.484 / 0.513 / 0.475, p = 0.787). That test
was limited to the 312 (unit, drug) pairs with a readable abstract. `max_phase_for_ind` gives a
COMPLETE, machine-readable outcome label for every drug-indication pair ChEMBL knows about, which is
the standard endpoint in the probability-of-success literature: success = advancing to the next
phase, phase 3 success = approval.

This script ONLY fetches. No joining, no scoring, no interpretation - those go in D178 behind a
pre-registered spec, because the temptation to tune is highest right after a null.

Output: data/raw/chembl_indications/drug_indication.tsv (+ _run.json)
"""
import json
import sys
import time
import urllib.request
from pathlib import Path

BASE = "https://www.ebi.ac.uk/chembl/api/data/drug_indication?format=json&limit=1000&offset={}"
UA = "repurposing-pipeline/1.0 (academic research; python-urllib)"   # no personal details, ever


def main(out_dir):
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    rows, offset, t0 = [], 0, time.time()
    while True:
        req = urllib.request.Request(BASE.format(offset), headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=120) as r:
            d = json.loads(r.read().decode())
        page = d.get("drug_indications", [])
        if not page:
            break
        for x in page:
            refs = x.get("indication_refs") or []
            ncts = sorted({n.strip() for rf in refs if rf.get("ref_type") == "ClinicalTrials"
                           for n in str(rf.get("ref_id", "")).split(",") if n.strip().startswith("NCT")})
            rows.append(dict(molecule_chembl_id=x.get("molecule_chembl_id"),
                             parent_molecule_chembl_id=x.get("parent_molecule_chembl_id"),
                             max_phase_for_ind=x.get("max_phase_for_ind"),
                             efo_id=x.get("efo_id"), efo_term=x.get("efo_term"),
                             mesh_id=x.get("mesh_id"), mesh_heading=x.get("mesh_heading"),
                             n_nct=len(ncts), ncts="|".join(ncts)))
        offset += len(page)
        total = d.get("page_meta", {}).get("total_count")
        print(f"  {offset:,} / {total:,}", flush=True)
        if total and offset >= total:
            break
        time.sleep(0.25)                      # be polite to EBI
    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(out / "drug_indication.tsv", sep="\t", index=False)
    (out / "drug_indication_run.json").write_text(json.dumps(
        dict(source=BASE.split("?")[0], n_rows=int(len(df)),
             n_molecules=int(df.molecule_chembl_id.nunique()),
             n_indications=int(df.efo_id.nunique()),
             fetched=time.strftime("%Y-%m-%dT%H:%M:%S"), runtime_s=round(time.time() - t0, 1)), indent=1))
    print(f"\n{len(df):,} rows, {df.molecule_chembl_id.nunique():,} molecules, "
          f"{df.efo_id.nunique():,} indications in {time.time() - t0:.0f}s")
    print(df.max_phase_for_ind.value_counts(dropna=False).sort_index().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
