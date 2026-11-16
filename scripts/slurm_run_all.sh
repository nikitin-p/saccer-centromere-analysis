#!/usr/bin/env bash
# Example SLURM job for scripts/run_all.sh. Adjust partition/account and resources to the cluster.
#
#   sbatch scripts/slurm_run_all.sh config.scrap142.yaml
#   sbatch --time=3-00:00:00 --mem=64G scripts/slurm_run_all.sh config.loegler1086.yaml
#   sbatch scripts/slurm_run_all.sh config.loegler1086.yaml --from trees     # resume
#
# Extra arguments after the config go to run_all.sh. The job is restartable: finished steps
# are skipped (output/<set>/.done/), BLAST / minimap2 skip finished genomes.
#SBATCH --job-name=centromere
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/slurm-%x-%j.out
set -euo pipefail

# conda env (environment.yml); `conda activate` needs the shell hook in batch jobs
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV:-centromere}"

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
export THREADS="${SLURM_CPUS_PER_TASK:-4}"
bash scripts/run_all.sh "$@"
