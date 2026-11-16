#!/usr/bin/env bash
# BLAST reference centromeres against every genome assembly in GENOME_DIR.
#
# Accepts *.fa, *.fasta, *.fna, *.re.fa, each optionally gzipped.
# Sample name = file name without those extensions (e.g. GCA_949124585.1).
# Output: one 9-column TSV per sample, concatenated into COMBINED_OUTPUT:
#   sample qseqid qstart qend sseqid sstart send sseq evalue
#
# MODE=flanks (L2 locator): query = 2-kb S288C CEN flanks (scripts/make_cen_flanks.py),
# megablast against the same DBs; per-sample *_flanks.tsv and COMBINED_OUTPUT (default
# cen_flanks.tsv) with 12 columns:
#   sample qseqid qstart qend qlen sseqid sstart send pident length evalue bitscore

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

MODE="${MODE:-cen}"
GENOME_DIR="${GENOME_DIR:-genome_assemblies}"
DB_DIR="${DB_DIR:-blast_db_genome_assemblies}"
RESULTS_DIR="${RESULTS_DIR:-blast_results}"
THREADS="${THREADS:-1}"

case "$MODE" in
    cen)
        QUERY="${QUERY:-$SCRIPT_DIR/data/sacer_centromere.fasta}"
        COMBINED_OUTPUT="${COMBINED_OUTPUT:-centromeres.tsv}"
        SUFFIX="_db_results.tsv"
        BLAST_ARGS=(-evalue 1e-5 -word_size 7
                    -outfmt "6 qseqid qstart qend sseqid sstart send sseq evalue")
        ;;
    flanks)
        QUERY="${QUERY:-$SCRIPT_DIR/data/s288c_cen_flanks.fasta}"
        COMBINED_OUTPUT="${COMBINED_OUTPUT:-cen_flanks.tsv}"
        SUFFIX="_flanks.tsv"
        BLAST_ARGS=(-task megablast -evalue 1e-20
                    -outfmt "6 qseqid qstart qend qlen sseqid sstart send pident length evalue bitscore")
        ;;
    *)
        echo "Unknown MODE=$MODE (cen | flanks)" >&2
        exit 1
        ;;
esac

mkdir -p "$DB_DIR" "$RESULTS_DIR"

shopt -s nullglob
FILES=("$GENOME_DIR"/*.{fa,fasta,fna}{,.gz})
if [ ${#FILES[@]} -eq 0 ]; then
    echo "No assemblies (*.fa/*.fasta/*.fna[.gz]) found in $GENOME_DIR" >&2
    exit 1
fi

for FILE in "${FILES[@]}"; do
    BASENAME=$(basename "$FILE")
    BASENAME=${BASENAME%.gz}
    BASENAME=${BASENAME%.*}       # .fa / .fasta / .fna
    BASENAME=${BASENAME%.re}      # Loegler et al. *.re.fa

    DB_PATH="$DB_DIR/${BASENAME}_db"
    RESULT_PATH="$RESULTS_DIR/${BASENAME}${SUFFIX}"

    # Skip finished samples so an interrupted run over ~1000 genomes can resume
    if [ -s "$RESULT_PATH" ]; then
        echo "Skipping $BASENAME (results exist)"
        continue
    fi

    if [ ! -f "$DB_PATH.nsq" ] && [ ! -f "$DB_PATH.nal" ]; then
        if [[ "$FILE" == *.gz ]]; then
            gzip -cd "$FILE" | makeblastdb -in - -title "$BASENAME" -dbtype nucl -out "$DB_PATH" > /dev/null
        else
            makeblastdb -in "$FILE" -dbtype nucl -out "$DB_PATH" > /dev/null
        fi
    fi

    # Write to a temp file first so a crash never leaves a partial result that
    # the resume check above would mistake for a finished sample
    blastn \
        -query "$QUERY" \
        -db "$DB_PATH" \
        "${BLAST_ARGS[@]}" \
        -num_threads "$THREADS" | \
        awk -v name="$BASENAME" '{print name "\t" $0}' > "$RESULT_PATH.tmp"
    mv "$RESULT_PATH.tmp" "$RESULT_PATH"

    echo "Finished processing $FILE ($(wc -l < "$RESULT_PATH") hits)"
done

echo "Combining results..."
cat "$RESULTS_DIR"/*"$SUFFIX" > "$COMBINED_OUTPUT"
echo "Combined output written to $COMBINED_OUTPUT ($(wc -l < "$COMBINED_OUTPUT") hits)"
