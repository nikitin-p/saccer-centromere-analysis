# Settings for the Slurm scripts in scripts/cluster/. Any variable can be overridden from the
# environment: PARTITION=long bash scripts/cluster/submit.sh run scrap142

REPO=${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}

# --- Slurm ---
PARTITION=${PARTITION:-main}
NODES=${NODES:-"vavilov0 vavilov1"}   # `submit.sh check` runs once on each
# BLAST and minimap2 run as job arrays of SHARDS tasks, at most MAX_PARALLEL at a time
# (user limit on the cluster: 64 CPU / 512 GB at once, over all jobs)
SHARDS=${SHARDS:-16}
MAX_PARALLEL=${MAX_PARALLEL:-8}

# --- conda env from environment.yml ---
CONDA_ENV=${CONDA_ENV:-centromere}
if [ -z "${CONDA_BASE:-}" ]; then
    if command -v conda > /dev/null 2>&1; then CONDA_BASE=$(conda info --base 2> /dev/null)
    else CONDA_BASE=$HOME/miniconda3; fi
fi

# --- datasets: name -> config (output goes to output/<name>/) ---
DATASETS=${DATASETS:-"scrap142 loegler1086"}
config_of() { echo "config.$1.yaml"; }
