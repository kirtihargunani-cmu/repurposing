#!/usr/bin/env bash
# Step 9 replication cohorts (PI decision 2026-10-08: two per cancer). All open access.
#   breast      METABRIC (Illumina HT-12 microarray, ~1,980 tumors; cBioPortal datahub, pinned commit)
#               CPTAC-BRCA (GDC CPTAC-2 STAR counts, 133 tumors; clinical/PAM50 + proteome/phosphoproteome from
#               cBioPortal brca_cptac_2020)
#   colorectal  GSE39582 (Affymetrix U133 Plus 2, 585 samples incl. 19 normals; GEO series matrix + GPL570 annotation)
#               CPTAC-COAD (GDC CPTAC-2 STAR counts, 105 colon + 1 rectum tumors; clinical + proteome/phosphoproteome
#               from cBioPortal coad_cptac_2019)
# CPTAC RNA: tumors only (no normal RNA-seq) -> patient-level replication, not tumor-vs-normal.
# Output: data/raw/replication/<cohort>/; GDC file list in gdc_manifest.tsv; sha256 appended to data/raw/checksums.sha256
set -euo pipefail
cd "$(dirname "$0")/.."
R=data/raw/replication
DATAHUB_SHA=18763f8eee86d4886d7210068f745f1e509498ee   # cBioPortal/datahub master, pinned 2026-10-08
DH=https://media.githubusercontent.com/media/cBioPortal/datahub/$DATAHUB_SHA/public
get() {  # skip files already complete (resuming a finished file returns HTTP 416)
  mkdir -p "$1"; local out="$1/${3:-$(basename "$2")}"
  [ -s "$out" ] && [ "$(curl -fsSIL "$2" | grep -i '^content-length' | tail -1 | tr -dc 0-9)" = "$(stat -c %s "$out")" ] && return 0
  curl -fsSL --retry 5 -C - -o "$out" "$2"
}

# METABRIC
for f in data_mrna_illumina_microarray.txt data_clinical_patient.txt data_clinical_sample.txt; do
  get $R/metabric $DH/brca_metabric/$f
done
echo "$DATAHUB_SHA" > $R/metabric/DATAHUB_COMMIT

# CPTAC clinical / subtype / protein tables (cBioPortal)
for s in brca_cptac_2020:cptac_brca coad_cptac_2019:cptac_coad; do
  study=${s%%:*}; dir=$R/${s##*:}
  for f in data_clinical_patient.txt data_clinical_sample.txt data_protein_quantification.txt data_phosphoprotein_quantification.txt; do
    get $dir/cbioportal $DH/$study/$f
  done
  echo "$DATAHUB_SHA" > $dir/cbioportal/DATAHUB_COMMIT
done

# GSE39582
get $R/gse39582 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39582/matrix/GSE39582_series_matrix.txt.gz
get $R/gse39582 https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz

# CPTAC-2 RNA-seq STAR counts (GDC open access), split by primary site
python3 - "$R" <<'PY'
import json, sys, tarfile, io, time, urllib.request
from pathlib import Path
R = Path(sys.argv[1])
q = {"filters": {"op": "and", "content": [
        {"op": "=", "content": {"field": "cases.project.project_id", "value": "CPTAC-2"}},
        {"op": "=", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
        {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}},
        {"op": "in", "content": {"field": "cases.primary_site", "value": ["breast", "colon", "rectum"]}}]},
     "fields": "file_id,file_name,md5sum,cases.submitter_id,cases.primary_site,cases.samples.submitter_id,cases.samples.sample_type",
     "format": "JSON", "size": 1000}
req = urllib.request.Request("https://api.gdc.cancer.gov/files", json.dumps(q).encode(), {"Content-Type": "application/json"})
hits = json.load(urllib.request.urlopen(req))["data"]["hits"]
rows = []
for h in hits:
    c = h["cases"][0]; site = c["primary_site"].lower()
    cohort = "cptac_brca" if site == "breast" else "cptac_coad"
    rows.append((cohort, h["file_id"], h["file_name"], h["md5sum"], c["submitter_id"],
                 ";".join(s["submitter_id"] for s in c["samples"]), ";".join(s["sample_type"] for s in c["samples"]), site))
for cohort in ("cptac_brca", "cptac_coad"):
    d = R / cohort / "star_counts"; d.mkdir(parents=True, exist_ok=True)
    sel = [r for r in rows if r[0] == cohort]
    with open(R / cohort / "gdc_manifest.tsv", "w") as fh:
        fh.write("file_id\tfile_name\tmd5\tcase\tsamples\tsample_types\tprimary_site\n")
        for r in sel: fh.write("\t".join(r[1:]) + "\n")
    todo = [r for r in sel if not (d / r[2]).exists()]
    for i in range(0, len(todo), 20):
        ids = [r[1] for r in todo[i:i + 20]]
        req = urllib.request.Request("https://api.gdc.cancer.gov/data", json.dumps({"ids": ids}).encode(),
                                     {"Content-Type": "application/json"})
        for attempt in range(5):  # GDC resets long transfers now and then
            try:
                data = urllib.request.urlopen(req, timeout=600).read(); break
            except (ConnectionError, OSError) as e:
                print(f"{cohort} batch {i // 20}: {e}; retry {attempt + 1}", flush=True); time.sleep(15)
        else:
            sys.exit(f"{cohort}: GDC download failed after 5 attempts")
        if len(ids) == 1:
            (d / todo[i][2]).write_bytes(data)
        else:
            with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tf:
                for m in tf.getmembers():
                    if m.isfile() and m.name.endswith(".tsv"):
                        (d / Path(m.name).name).write_bytes(tf.extractfile(m).read())
    import hashlib
    bad = [r[2] for r in sel if hashlib.md5((d / r[2]).read_bytes()).hexdigest() != r[3]]
    print(f"{cohort}: {len(sel)} files, md5 mismatches: {len(bad)}")
    if bad: sys.exit(f"md5 mismatch: {bad[:5]}")
PY

(cd data/raw && find replication -type f ! -path '*/star_counts/*' | sort | xargs sha256sum | sed 's#  #  ./#') > data/raw/.new.sha256
cat data/raw/checksums.sha256 data/raw/.new.sha256 | sort -u -k2 -o data/raw/checksums.sha256 && rm data/raw/.new.sha256
echo "REPLICATION DOWNLOAD DONE"
