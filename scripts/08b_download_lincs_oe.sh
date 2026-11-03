#!/usr/bin/env bash
# LINCS 2020 level 5 overexpression signatures (trt_oe, MODZ z; 34,171 signatures x 12,328 genes), 1.7 GB GCTX from
# the public S3 build, for the D84 perturbation benchmark. Resumable (curl -C -); size checked against S3, then
# sha256 appended to data/raw/checksums.sha256. Subsetting: scripts/79u_overexpression_subset.py.
set -euo pipefail
cd "$(dirname "$0")/.."
D=data/raw/lincs2020
F=level5_beta_trt_oe_n34171x12328.gctx
URL=https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/$F
SIZE=1685209064
mkdir -p $D
until [ -f "$D/$F" ] && [ "$(stat -c%s "$D/$F")" = "$SIZE" ]; do
  curl -fsSL --retry 10 --retry-delay 10 -C - -o "$D/$F" "$URL" || { echo "curl exited $?; resuming"; sleep 10; }
done
echo "size OK: $SIZE"
(cd data/raw && sha256sum lincs2020/$F | sed 's#  #  ./#') > data/raw/.new.sha256
cat data/raw/checksums.sha256 data/raw/.new.sha256 | sort -u -k2 -o data/raw/checksums.sha256 && rm data/raw/.new.sha256
echo "LINCS LEVEL5 OE DONE"
