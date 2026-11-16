#!/usr/bin/env bash
# What a node has: CPUs, memory, disk, conda + the env's tools, which download hosts it reaches.
# Read-only. `submit.sh check` runs it on every node in NODES.
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/config.sh"

echo "== node"; hostname; nproc; free -g 2>/dev/null | head -2
echo "== disk at $REPO"; df -h "$REPO" | tail -1
echo "== conda ($CONDA_BASE)"
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    set +u; source "$CONDA_BASE/etc/profile.d/conda.sh"; set -u
    conda --version
    if conda activate "$CONDA_ENV" 2> /dev/null; then
        echo "env $CONDA_ENV: active"
        for t in python blastn makeblastdb tblastn mafft minimap2 hmmbuild nhmmer Rscript quarto; do
            printf '  %-12s %s\n' "$t" "$(command -v "$t" || echo MISSING)"
        done
        python -c "import centromere" 2> /dev/null && echo "  package centromere: ok" \
            || echo "  package centromere: MISSING (pip install -e . --no-deps)"
    else
        echo "env $CONDA_ENV: not created (submit.sh env)"
    fi
else
    echo "no conda.sh in $CONDA_BASE"
fi

echo "== hosts (any HTTP code = reachable, 000 = blocked)"
for h in https://zenodo.org https://static-content.springer.com https://hgdownload.soe.ucsc.edu \
         http://sgd-archive.yeastgenome.org https://conda.anaconda.org https://pypi.org \
         https://github.com; do
    code=$(curl -sS -o /dev/null -w '%{http_code}' --max-time 20 "$h/" 2> /dev/null) || code=000
    printf '  %-40s %s\n' "$h" "$code"
done
