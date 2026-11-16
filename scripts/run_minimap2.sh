#!/usr/bin/env bash
# L3 step 0.1: align every assembly to S288C R64 (nuclear chromosomes only) with minimap2.
#
# Usage: bash scripts/run_minimap2.sh CONFIG [PYTHON]
#   CONFIG  e.g. config.scrap.yaml (sections input, cde, kinetochore, locate)
#   PYTHON  interpreter with PyYAML; default .venv/bin/python
#
# Writes <locate.paf_dir>/<accession>.paf; an existing non-empty PAF is skipped
# (same resume behaviour as blast.sh). The reference without chrM is written once
# to <paf_dir>/sacCer3.nuclear.fa.
# SHARD/NSHARDS (env): align only every NSHARDS-th assembly starting at SHARD (Slurm array task).
set -euo pipefail
SHARD=${SHARD:-0}
NSHARDS=${NSHARDS:-1}

CONFIG=$1
PYTHON=${2:-.venv/bin/python}

read -r MM2 PRESET THREADS PAF_DIR REF GDIR MAP < <("$PYTHON" - "$CONFIG" <<'PY'
import os, sys, yaml
c = yaml.safe_load(open(sys.argv[1]))
l = c["locate"]
x = lambda p: os.path.expanduser(str(p))
print(x(l.get("minimap2") or "minimap2"), l.get("preset", "asm5"), l.get("threads", 4), x(l["paf_dir"]),
      x(c["kinetochore"]["s288c_genome"]), x(c["cde"]["genome_dir"]), x(c["input"]["sample_map"]))
PY
)

mkdir -p "$PAF_DIR"
NUC="$PAF_DIR/sacCer3.nuclear.fa"
if [ ! -s "$NUC" ]; then   # temp + mv: array tasks may start together
    awk '/^>/{keep = ($1 != ">chrM")} keep' "$REF" > "$NUC.tmp.$SHARD"
    mv "$NUC.tmp.$SHARD" "$NUC"
fi

i=-1
tail -n +2 "$MAP" | cut -f1 | while read -r acc; do
    i=$((i + 1))
    [ $((i % NSHARDS)) -eq "$SHARD" ] || continue
    out="$PAF_DIR/$acc.paf"
    if [ -s "$out" ]; then echo "skip $acc (done)"; continue; fi
    echo "minimap2 -x $PRESET $acc"
    "$MM2" -x "$PRESET" -c --cs -t "$THREADS" "$NUC" "$GDIR/$acc.fasta.gz" > "$out.tmp" 2> "$out.log"
    mv "$out.tmp" "$out"
done
