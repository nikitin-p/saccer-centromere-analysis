#!/usr/bin/env bash
# Run PCAn v1.0 (Helsen et al. 2026; Zenodo 10.5281/zenodo.17293587) on every
# assembly, Saccharomyces motifs, no BLAST synteny check.
#
# Usage: bash scripts/run_pcan.sh PCAN_DIR GENOME_DIR OUT_DIR [PYTHON]
#   PCAN_DIR    folder holding AutomatedCENretrieval.py (unzipped PCAn/)
#   GENOME_DIR  *.fasta.gz / *.fasta / *.fa assemblies
#   OUT_DIR     gets <sample>/CENsequences.txt per assembly
#   PYTHON      interpreter of the PCAn conda env (needs fimo on PATH next to it);
#               default ~/miniconda3/envs/pcan/bin/python
#
# PCAn asks three questions on stdin: genome path, genus index (0 = Saccharomyces),
# synteny BLAST (No). It writes into the current directory, so each run gets its own.
set -euo pipefail

PCAN_DIR=$(cd "$1" && pwd)
GENOME_DIR=$(cd "$2" && pwd)
mkdir -p "$3"
OUT_DIR=$(cd "$3" && pwd)
PYTHON=${4:-$HOME/miniconda3/envs/pcan/bin/python}
export PATH="$(dirname "$PYTHON"):$PATH"   # fimo lives in the env's bin/

command -v fimo >/dev/null || { echo "fimo not found next to $PYTHON" >&2; exit 1; }

for genome in "$GENOME_DIR"/*.fasta.gz "$GENOME_DIR"/*.fasta "$GENOME_DIR"/*.fa; do
    [ -e "$genome" ] || continue
    sample=$(basename "$genome"); sample=${sample%.gz}; sample=${sample%.*}
    run="$OUT_DIR/$sample"
    if [ -s "$run/CENsequences.txt" ]; then
        echo "skip $sample (done)"; continue
    fi
    rm -rf "$run"; mkdir -p "$run"
    fasta="$genome"
    if [[ "$genome" == *.gz ]]; then
        fasta="$run/genome.fasta"
        gunzip -c "$genome" > "$fasta"
    fi
    echo "PCAn $sample"
    (cd "$run" && printf '%s\n0\nNo\n' "$fasta" | "$PYTHON" "$PCAN_DIR/AutomatedCENretrieval.py" > pcan.log 2>&1)
    rm -rf "$run/FIMOresults"
    [[ "$fasta" == "$run/genome.fasta" ]] && rm -f "$fasta"
done
