#!/usr/bin/env bash
# D97 production switch (user decision 2026-11-03, docs/DECISIONS.md). <= 24 threads. Backups made beforehand (_pre_d97).
set -euo pipefail
cd "$(dirname "$0")/../.."
RS=~/miniforge3/envs/repurposing/bin/Rscript; PY=~/miniforge3/envs/celligner/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
step() { echo "$(date +%T) $1"; }
step "1 regulons: add the 12 networks to the production file"
$RS -e 'p <- "data/processed/step5/regulons_aracne_project_v2_tf.rds"; a <- readRDS(p); b <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf_d96all.rds")
  m <- read.delim("config/group_networks.tsv"); n <- setdiff(unique(m$net), names(a)); stopifnot(all(n %in% names(b)))
  saveRDS(c(a, b[n]), p); cat("added:", n, "\n")'
sha256sum data/processed/step5/regulons_aracne_project_v2_tf.rds > /tmp/claude-1014/d97_sha.txt
$PY - <<'P'
from pathlib import Path
f = Path("results/step5/regulon_checksums.sha256"); new = open("/tmp/claude-1014/d97_sha.txt").read().split()
lines = [l for l in f.read_text().splitlines() if not l.endswith(" " + new[1])] + [f"{new[0]}  {new[1]}"]
f.write_text("\n".join(lines) + "\n")
P
G=$(tail -n +2 config/group_networks.tsv | cut -f1 | paste -sd,)
step "2 Step 6: switched groups (group + patient VIPER) and subgroups, merged into production"
REGULON=aracne_project VERSION=v2 REGSET=tf RES_TAG=_d97run ONLY_GROUPS=$G SEED_GROUP=1 N_CORES=12 $RS scripts/50_viper.R > logs/d97_viper.log 2>&1
REGULON=aracne_project VERSION=v2 REGSET=tf STRATA=results/step10_v2/stratum_patients.tsv RES_TAG=_d97run \
  ONLY_GROUPS=LGG-noncodel-LGr2,LGG-noncodel-LGr3,PHEO-kinase-signaling,PHEO-pseudohypoxia-Wnt SEED_GROUP=1 N_CORES=4 $RS scripts/50_viper.R > logs/d97_viper_strata.log 2>&1
$PY scripts/d97/d97_merge_step6.py > logs/d97_merge.log 2>&1
step "3 P1 lists and drug activity"
VERSION=v2 REGSET=tf $RS scripts/79h_prolif_regulators.R > logs/d97_79h.log 2>&1
VERSION=v2 REGSET=tf PROFILE=specific N_CORES=20 $RS scripts/79d_drug_viper.R > logs/d97_79d.log 2>&1
step "4 rankings: consensus and TF arm, groups and subgroups"
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
( VERSION=v2 REGULON=consensus $PY scripts/79e_drug_reversal.py > logs/d97_79e_cons.log 2>&1
  VERSION=v2 REGULON=consensus STRATA=1 $PY scripts/79e_drug_reversal.py > logs/d97_79e_cons_strata.log 2>&1
  VERSION=v2 REGULON=consensus STRATA=1 $PY scripts/79i_split_reversal.py > logs/d97_79i_cons.log 2>&1 ) &
( VERSION=v2 REGSET=tf $PY scripts/79e_drug_reversal.py > logs/d97_79e_tf.log 2>&1
  VERSION=v2 REGSET=tf STRATA=1 $PY scripts/79e_drug_reversal.py > logs/d97_79e_tf_strata.log 2>&1
  VERSION=v2 REGSET=tf STRATA=1 $PY scripts/79i_split_reversal.py > logs/d97_79i_tf.log 2>&1 ) &
wait
step "5 downstream"
for s in 79w_receptor_gate 79p_fda_ranking 79o_literature 79y_standin; do step "  $s"; $PY scripts/$s.py > logs/d97_$s.log 2>&1; done
step "  79y2 counterfactual"; bash scripts/79y2_borrowed_counterfactual.sh > logs/d97_79y2_sh.log 2>&1; $PY scripts/79y2_standin_vs_borrowed.py > logs/d97_79y2.log 2>&1
for s in 79l_known_use 79n_dependency 79m_gate_sensitivity 79x_rankings_page 79z_standin_page 79q_fda_page; do step "  $s"; $PY scripts/$s.py > logs/d97_$s.log 2>&1; done
step "D97 SWITCH DONE"
