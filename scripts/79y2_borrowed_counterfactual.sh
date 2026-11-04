#!/usr/bin/env bash
# Borrowed counterfactual for the stand-in comparison (79y2; user 2026-11-02).
# Builds a sandbox (default logs/sbx_borrowed_cf) whose results/lincs holds only links to 79e / 79i inputs plus a copy of
# disease_drugs.tsv in which every FDA drug tested in its cohort's lines (direct) with >= 3 lines elsewhere is switched to
# borrowed (same for stratum_drugs_v2.tsv, the 14 subgroups); reruns 79e (diseases, then STRATA=1) and 79i (STRATA=1:
# diseases and subgroups) there (consensus, v2), then 79y2 compares stand-ins with these borrowed scores.
# Production results are never written: every sandbox output is a new file in the sandbox.
# Env: PY (celligner python), SBX (sandbox dir), THREADS (BLAS threads, default 8).
set -euo pipefail
cd "$(dirname "$0")/.."
R=$(pwd -P)
PY=${PY:-$HOME/miniforge3/envs/celligner/bin/python}
SBX=${SBX:-$R/logs/sbx_borrowed_cf}
export OMP_NUM_THREADS=${THREADS:-8} OPENBLAS_NUM_THREADS=${THREADS:-8} MKL_NUM_THREADS=${THREADS:-8}
rm -rf "$SBX"; mkdir -p "$SBX/scripts" "$SBX/results/lincs"
cp scripts/79e_drug_reversal.py scripts/79i_split_reversal.py scripts/79i2_target_presence.py scripts/79w_receptor_gate.py "$SBX/scripts/"   # 79i2: D102c gate (reads 79w helpers)
for d in config data R; do ln -s "$R/$d" "$SBX/$d"; done
for d in "$R"/results/*; do n=$(basename "$d"); [ "$n" != lincs ] && ln -s "$d" "$SBX/results/$n"; done
for f in drug_transfer_confidence_v2.tsv prolif_regulators_v2_tf.tsv prolif_regulators_collectri.tsv drug_annotations.tsv; do   # drug_annotations: 79i2 (via 79w)
  ln -s "$R/results/lincs/$f" "$SBX/results/lincs/$f"; done
"$PY" - "$R" "$SBX" <<'EOF'
import sys, pandas as pd
R, B = sys.argv[1], sys.argv[2]
fm = set(pd.read_csv(f"{R}/results/lincs/fda_matches.tsv", sep="\t").pert_id)
for f in ("disease_drugs.tsv", "stratum_drugs_v2.tsv"):
    d = pd.read_csv(f"{R}/results/lincs/{f}", sep="\t", low_memory=False)
    m = (d.source == "direct") & d.pert_id.isin(fm) & (d.n_lines_elsewhere >= 3)
    d.loc[m, "source"] = "borrowed"
    d.to_csv(f"{B}/results/lincs/{f}", sep="\t", index=False)
    print(f, "switched direct FDA rows to borrowed:", int(m.sum()))
EOF
(cd "$SBX" && VERSION=v2 REGULON=consensus "$PY" scripts/79e_drug_reversal.py > e.log 2>&1)
(cd "$SBX" && VERSION=v2 REGULON=consensus STRATA=1 "$PY" scripts/79e_drug_reversal.py > e_sub.log 2>&1)
(cd "$SBX" && VERSION=v2 REGULON=consensus STRATA=1 "$PY" scripts/79i_split_reversal.py > i.log 2>&1)
(cd "$SBX" && "$PY" scripts/79i2_target_presence.py > i2.log 2>&1)
SBX="$SBX" "$PY" scripts/79y2_standin_vs_borrowed.py
