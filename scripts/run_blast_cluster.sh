#!/usr/bin/env bash
# Run BLAST pipeline on cluster with ~/sacer/ directory structure

set -e

CLUSTER_BASE="${CLUSTER_BASE:-$HOME/sacer}"

export QUERY="$CLUSTER_BASE/sacer_centromere.fasta"
export GENOME_DIR="$CLUSTER_BASE/genome_assemblies"
export DB_DIR="$CLUSTER_BASE/blast_db_genome_assemblies"
export RESULTS_DIR="$CLUSTER_BASE/blast_results"
export COMBINED_OUTPUT="$CLUSTER_BASE/centromeres.tsv"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Running BLAST with cluster paths:"
echo "  QUERY: $QUERY"
echo "  GENOME_DIR: $GENOME_DIR"
echo "  DB_DIR: $DB_DIR"
echo "  RESULTS_DIR: $RESULTS_DIR"
echo "  COMBINED_OUTPUT: $COMBINED_OUTPUT"

bash "$SCRIPT_DIR/blast.sh"
