#!/usr/bin/env bash
# Conda env "aracne" + ARACNe-AP for scripts/43_run_aracne.sh: OpenJDK 17 (17.0.18 as used for the v1 and v2 builds),
# ARACNe-AP (califano-lab/ARACNe-AP commit 30967ea) built with Ant into $ARACNE_DIR/dist/aracne.jar.
# Env vars: ARACNE_ENV (default aracne), ARACNE_DIR (default ~/tools/ARACNe-AP).
set -euo pipefail
cd "$(dirname "$0")/.."
MAMBA=${MAMBA:-~/miniforge3/bin/mamba}
ENV=${ARACNE_ENV:-aracne}
T=${ARACNE_DIR:-~/tools/ARACNe-AP}
$MAMBA env list | grep -qE "^$ENV " || $MAMBA create -y -n $ENV -c conda-forge openjdk=17.0.18
$MAMBA install -y -n $ENV -c conda-forge ant   # build only
[ -d $T ] || git clone -q https://github.com/califano-lab/ARACNe-AP.git $T
git -C $T checkout -q 30967ea
PREFIX=$($MAMBA env list | awk -v e=$ENV '$1 == e {print $NF}')
(cd $T && PATH=$PREFIX/bin:$PATH JAVA_HOME=$PREFIX/lib/jvm ant -q main)
$PREFIX/bin/java -jar $T/dist/aracne.jar 2>&1 | head -3 || true
echo "ARACNE SETUP DONE: $PREFIX/bin/java, $T/dist/aracne.jar"
