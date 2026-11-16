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

# L2 locator: 2-kb S288C CEN flanks (data/s288c_cen_flanks.fasta) against the same DBs.
# QUERY / COMBINED_OUTPUT are exported above for the CEN mode, so they are overridden here.
run_flanks() {
    echo "Running flank BLAST (L2 locator)..."
    MODE=flanks QUERY="$REPO_DIR/data/s288c_cen_flanks.fasta" \
        COMBINED_OUTPUT="$CLUSTER_BASE/cen_flanks.tsv" bash "$SCRIPT_DIR/blast.sh"
}

run_analysis() {
    echo "Running analysis pipeline..."
    python3 "$SCRIPT_DIR/run_pipeline.py" --config "$REPO_DIR/config.cluster.yaml" "$@"
}

case "${1:-all}" in
    blast)
        run_blast
        ;;
    flanks)
        run_flanks
        ;;
    analysis)
        shift
        run_analysis "$@"
        ;;
    all)
        run_blast
        run_flanks
        run_analysis
        ;;
    *)
        echo "Usage: $0 [blast|flanks|analysis|all] [--start-from N]"
        exit 1
        ;;
esac
