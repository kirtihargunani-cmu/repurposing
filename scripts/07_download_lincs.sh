#!/usr/bin/env bash
# LINCS 2020 (CMap / clue.io "beta" release) metadata from the public S3 build: cell lines, signatures (which
# perturbagen in which line, dose, time, QC), compounds (targets, MOA) and genes. Signature matrices are not pulled
# here. Single-part S3 ETags are md5 and verified; siginfo is a multipart upload (ETag not an md5), so sha256 only.
# Appends sha256 to data/raw/checksums.sha256.
set -euo pipefail
cd "$(dirname "$0")/.."
D=data/raw/lincs2020
B=https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020
mkdir -p $D
while read -r name md5; do
  [ -s "$D/$name" ] || curl -fsSL --retry 5 -C - -o "$D/$name" "$B/$name"
  [ "$md5" = "-" ] || echo "$md5  $D/$name" | md5sum -c --quiet -
done <<'LIST'
cellinfo_beta.txt c4686b4bcd2bad8fa64e229932c8d486
compoundinfo_beta.txt bf8e3a15ad026b47903c98d625195d24
geneinfo_beta.txt 45c725d17ce6c377f1e7de07b821a5f0
siginfo_beta.txt -
LIST
(cd data/raw && find lincs2020 -type f | sort | xargs sha256sum | sed 's#  #  ./#') > data/raw/.new.sha256
cat data/raw/checksums.sha256 data/raw/.new.sha256 | sort -u -k2 -o data/raw/checksums.sha256 && rm data/raw/.new.sha256
echo "LINCS DONE"
