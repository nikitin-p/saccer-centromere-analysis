#!/usr/bin/env bash
# Body of every Slurm job from submit.sh: activate the conda env, cd to the repo, run "$@".
#   bash scripts/cluster/job.sh bash scripts/run_all.sh config.scrap142.yaml --only hmm
# In a job array (JOB_SHARDS set by submit.sh) the task gets SHARD=<task id> NSHARDS=JOB_SHARDS,
# which blast.sh / run_minimap2.sh / run_all.sh use to take every NSHARDS-th genome.
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/config.sh"

if [ -z "${NO_CONDA:-}" ]; then
    [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ] || { echo "no conda in $CONDA_BASE (set CONDA_BASE)" >&2; exit 1; }
    set +u; source "$CONDA_BASE/etc/profile.d/conda.sh"; conda activate "$CONDA_ENV"; set -u
fi
cd "$REPO"
export THREADS=${SLURM_CPUS_PER_TASK:-${THREADS:-4}}
unset SHARD NSHARDS
if [ -n "${JOB_SHARDS:-}" ]; then
    export SHARD=${SLURM_ARRAY_TASK_ID:?JOB_SHARDS set outside a job array} NSHARDS=$JOB_SHARDS
fi
echo "[$(date '+%F %T')] $(hostname -s) cpus=$THREADS${NSHARDS:+ shard=$SHARD/$NSHARDS}: $*"
"$@"
echo "[$(date '+%F %T')] done"
