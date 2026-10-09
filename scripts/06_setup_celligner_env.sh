#!/usr/bin/env bash
# Conda env "celligner" for scripts/70_celligner.py: Celligner Python port (broadinstitute/celligner d9c9246, the 2023
# refactor with makeNewReference / computeMetricsForOutput; PyPI 1.1.0 is the older API) + its mnnpy fork
# (DeKegel/mnnpy bf1dd73; the submodule URL is SSH, so it is cloned over HTTPS). Versions follow the repo's
# requirements.txt where conda has them; umap-learn 0.5.3 / setuptools < 70 for sklearn 1.1 compatibility.
# rpy2 needs R on PATH: run with PATH=~/miniforge3/envs/celligner/bin:$PATH.
set -euo pipefail
cd "$(dirname "$0")/.."
MAMBA=~/miniforge3/bin/mamba
$MAMBA create -y -n celligner -c conda-forge -c bioconda python=3.10 numpy=1.23.3 pandas=1.4.4 scikit-learn=1.1.2 \
  scipy=1.9.1 scanpy=1.9.1 anndata=0.8.0 louvain=0.7.1 python-igraph=0.9.11 "rpy2=3.5" cython=0.29.32 \
  "matplotlib<3.7" umap-learn=0.5.3 pynndescent=0.5.8 "setuptools<70" numba r-base bioconductor-limma pyyaml gxx_linux-64
T=~/tools/celligner
[ -d $T ] || git clone -q https://github.com/broadinstitute/celligner.git $T
git -C $T checkout -q d9c9246
[ -f $T/mnnpy/setup.py ] || { rmdir $T/mnnpy 2>/dev/null || true; git clone -q https://github.com/DeKegel/mnnpy.git $T/mnnpy; }
git -C $T/mnnpy checkout -q bf1dd736f5b8a6b4b571e1ed572e73451162316a
PY=~/miniforge3/envs/celligner/bin/python
(cd $T/mnnpy && $PY -m pip install --no-deps --no-build-isolation -q .)
(cd $T && $PY -m pip install --no-deps -e . -q)
PATH=~/miniforge3/envs/celligner/bin:$PATH $PY -c 'from celligner import Celligner; print("celligner OK")'
$MAMBA env export -n celligner > environment.celligner.lock.yml
echo "CELLIGNER ENV DONE"
