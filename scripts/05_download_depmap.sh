#!/usr/bin/env bash
# Celligner / drug-side sources: DepMap 24Q4 Public from Figshare (article 27993248, CC BY 4.0). The DepMap portal
# download API sits behind a browser challenge; 24Q4 is the newest full release on Figshare. Verifies Figshare md5,
# appends sha256 to data/raw/checksums.sha256.
set -euo pipefail
cd "$(dirname "$0")/.."
D=data/raw/depmap_24q4
mkdir -p $D
# name  figshare-file-id  md5
while read -r name id md5; do
  [ -s "$D/$name" ] || curl -fsSL --retry 5 -C - -o "$D/$name" "https://ndownloader.figshare.com/files/$id"
  echo "$md5  $D/$name" | md5sum -c --quiet -
done <<'EOF'
README.txt 51065795 54628c15a10b9ff6536db245cfed5231
Model.csv 51065297 675210d17675f3517b0ce39a3c274f16
OmicsProfiles.csv 51065723 7195088a5fc41eddb3f896c6e984a992
OmicsDefaultModelProfiles.csv 51065339 4d342a03701c0a4c0742fedb56c1de40
OmicsExpressionProteinCodingGenesTPMLogp1.csv 51065489 71794802b750ce77c422dad0720a40af
OmicsSomaticMutations.csv 51065732 7bdba347a1602fe96d5654a74d6e52f1
OmicsSomaticMutationsMatrixHotspot.csv 51065750 a5aeb1deef897ead3c955c372148d840
OmicsSomaticMutationsMatrixDamaging.csv 51065747 cb20fdbe1cf3b9b0d8ed4f53e1f399b6
OmicsFusionFiltered.csv 51065693 16ebb6c5cab1e066d904ac6d7488f82e
OmicsCNGene.csv 51065324 4696565fa74dc805983bc278b4a41802
CRISPRGeneEffect.csv 51064667 6edf7ade09b9b34199210b559d4745d3
CRISPRInferredCommonEssentials.csv 51064916 545edef0d377f252db6e5e675201bd77
EOF
(cd data/raw && find depmap_24q4 -type f | sort | xargs sha256sum | sed 's#  #  ./#') > data/raw/.new.sha256
cat data/raw/checksums.sha256 data/raw/.new.sha256 | sort -u -k2 -o data/raw/checksums.sha256 && rm data/raw/.new.sha256
echo "DEPMAP DONE"
