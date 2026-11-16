#!/usr/bin/env bash
# Create (or update) the conda env from environment.yml and install the package from src/.
# Needs conda.anaconda.org + pypi.org (see `submit.sh check`). Run by `submit.sh env`, or directly.
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/config.sh"
set +u; source "$CONDA_BASE/etc/profile.d/conda.sh"; set -u
cd "$REPO"

solver=conda
command -v mamba > /dev/null 2>&1 && solver=mamba
echo "$solver $($solver --version | head -1)"
if conda env list | awk '{print $1}' | grep -qx "$CONDA_ENV"; then
    "$solver" env update -n "$CONDA_ENV" -f environment.yml --prune
else
    "$solver" env create -n "$CONDA_ENV" -f environment.yml
fi
set +u; conda activate "$CONDA_ENV"; set -u
pip install -e . --no-deps
for t in blastn mafft minimap2 nhmmer Rscript quarto; do printf '%-10s %s\n' "$t" "$(command -v "$t" || echo MISSING)"; done
Rscript -e 'cat("DNAshapeR:", requireNamespace("DNAshapeR", quietly = TRUE), "\n")'
python -m pytest -q
