#!/usr/bin/env bash
# Run full pipeline on cluster with ~/sacer/ directory structure

set -e

CLUSTER_BASE="${CLUSTER_BASE:-$HOME/sacer}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"

# BLAST environment variables
export QUERY="$CLUSTER_BASE/sacer_centromere.fasta"
export GENOME_DIR="$CLUSTER_BASE/genome_assemblies"
export DB_DIR="$CLUSTER_BASE/blast_db_genome_assemblies"
export RESULTS_DIR="$CLUSTER_BASE/blast_results"
export COMBINED_OUTPUT="$CLUSTER_BASE/centromeres.tsv"

run_blast() {
    echo "Running BLAST..."
    echo "  QUERY: $QUERY"
    echo "  GENOME_DIR: $GENOME_DIR"
    echo "  OUTPUT: $COMBINED_OUTPUT"
    bash "$SCRIPT_DIR/blast.sh"
}

run_analysis() {
    echo "Running analysis pipeline..."
    python3 "$SCRIPT_DIR/run_pipeline.py" --config "$REPO_DIR/config.cluster.yaml" "$@"
}

case "${1:-all}" in
    blast)
        run_blast
        ;;
    analysis)
        shift
        run_analysis "$@"
        ;;
    all)
        run_blast
        run_analysis
        ;;
    *)
        echo "Usage: $0 [blast|analysis|all] [--start-from N]"
        exit 1
        ;;
esac
