#!/usr/bin/env bash
# Drug annotations: Broad Drug Repurposing Hub (2020-03-24 release; clinical phase, MOA, targets, disease area,
# indication; samples table links Broad IDs = LINCS pert_id and InChIKeys). md5 = S3 ETag. ChEMBL (routes, ATC,
# indications) is pulled through its REST API by scripts/79a_drug_annotations.py, raw JSON kept in data/raw/chembl.
set -euo pipefail
cd "$(dirname "$0")/.."
D=data/raw/repurposing_hub
B=https://s3.amazonaws.com/data.clue.io/repurposing/downloads
mkdir -p $D
while read -r name md5; do
  [ -s "$D/$name" ] || curl -fsSL --retry 5 -o "$D/$name" "$B/$name"
  echo "$md5  $D/$name" | md5sum -c --quiet -
done <<'LIST'
repurposing_drugs_20200324.txt b4b3233dad28d20a522c0200cc4123d1
repurposing_samples_20200324.txt 44067d6ea0e63c70f51e0469d60232ec
LIST
(cd data/raw && find repurposing_hub -type f | sort | xargs sha256sum | sed 's#  #  ./#') > data/raw/.new.sha256
cat data/raw/checksums.sha256 data/raw/.new.sha256 | sort -u -k2 -o data/raw/checksums.sha256 && rm data/raw/.new.sha256
echo "DRUG ANNOTATIONS DONE"
