#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

QUERY="${QUERY:-$SCRIPT_DIR/data/sacer_centromere.fasta}"
GENOME_DIR="${GENOME_DIR:-genome_assemblies}"
DB_DIR="${DB_DIR:-blast_db_genome_assemblies}"
RESULTS_DIR="${RESULTS_DIR:-blast_results}"
COMBINED_OUTPUT="${COMBINED_OUTPUT:-centromeres.tsv}"

mkdir -p "$DB_DIR"
mkdir -p "$RESULTS_DIR"

for FILE in "$GENOME_DIR"/*.re.fa; do
    BASENAME=$(basename "$FILE" .re.fa)

    DB_PATH="$DB_DIR/${BASENAME}_db"
    RESULT_PATH="$RESULTS_DIR/${BASENAME}_db_results.tsv"

    makeblastdb -in "$FILE" -dbtype nucl -out "$DB_PATH"

    blastn \
        -query "$QUERY" \
        -db "$DB_PATH" \
        -evalue 1e-5 \
        -outfmt "6 qseqid qstart qend sseqid sstart send sseq evalue" \
        -word_size 7 | \
        awk -v name="$BASENAME" '{print name "\t" $0}' > "$RESULT_PATH"

    echo "Finished processing $FILE"
done

echo "Combining results..."
cat "$RESULTS_DIR"/*_db_results.tsv > "$COMBINED_OUTPUT"
echo "Combined output written to $COMBINED_OUTPUT"
