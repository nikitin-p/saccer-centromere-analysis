#!/usr/bin/env bash
# Run every assembly-based step for one dataset, from BLAST to the rendered report.
#
#   bash scripts/run_all.sh CONFIG [--from STEP] [--only STEP] [--with STEP,...] [--list]
#
# CONFIG    config.scrap142.yaml | config.loegler1086.yaml | ...
# --from    start at STEP (earlier steps are assumed done); --only: run just STEP
# --with    also run optional steps: validate (A: Helsen/PCAn), peter748 (G), phenotypes
# --list    print the steps and exit
#
# Run from anywhere inside the conda env (`conda activate centromere`): every tool comes from
# PATH. A finished step leaves <output.dir>/.done/<step>; a rerun skips it (delete the marker
# or use --from/--only to redo it). Each step logs to <output.dir>/logs/<step>.log.
# THREADS (default 4) is passed to BLAST; the configs set align.workers / locate.threads.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1   # step logs show progress while a step runs (was empty for 18 h)
PY="${PYTHON:-python}"
THREADS="${THREADS:-4}"

STEPS=(blast_cen pipeline flanks_blast l2_flanks minimap2 l3_diversity hmm reconcile
       cen_like synthetic cdeii_steps kinetochore kinetochore_domains kinetochore_assoc
       transposons ltr_cde diversity trees grammar report)
OPTIONAL=(validate peter748 phenotypes)

usage() { sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
[ $# -ge 1 ] || usage 1
CONFIG="$1"; shift
FROM="" ONLY="" WITH=""
while [ $# -gt 0 ]; do
    case "$1" in
        --from) FROM="$2"; shift 2 ;;
        --only) ONLY="$2"; shift 2 ;;
        --with) WITH="$2"; shift 2 ;;
        --list) printf '%s\n' "${STEPS[@]}"; echo "optional: ${OPTIONAL[*]}"; exit 0 ;;
        -h|--help) usage ;;
        *) echo "unknown argument: $1" >&2; usage 1 ;;
    esac
done
[ -f "$CONFIG" ] || { echo "config not found: $CONFIG" >&2; exit 1; }

# Paths from the config
read -r OUT INPUT GDIR FLANK_HITS < <("$PY" - "$CONFIG" <<'PY'
import sys, yaml
c = yaml.safe_load(open(sys.argv[1]))
print(c["output"]["dir"], c["input"]["path"], c["cde"]["genome_dir"], c["locate"]["flank_hits"])
PY
)
mkdir -p "$OUT/.done" "$OUT/logs"

for tool in blastn makeblastdb tblastn mafft minimap2 hmmbuild nhmmer Rscript quarto; do
    command -v "$tool" > /dev/null || { echo "missing on PATH: $tool (conda activate centromere?)" >&2; exit 1; }
done
ls "$GDIR"/*.fasta.gz > /dev/null 2>&1 || { echo "no assemblies in $GDIR (run scripts/download_data.py)" >&2; exit 1; }
[ -s data/external/s288c/sacCer3.fa ] || { echo "no S288C genome: python scripts/download_data.py --set reference" >&2; exit 1; }

blast_env() {   # shared BLAST DBs per dataset
    export GENOME_DIR="$GDIR" DB_DIR="$OUT/blastdb" THREADS
}

run_step() {
    local step=$1
    case "$step" in
        blast_cen)
            blast_env
            RESULTS_DIR="$OUT/blast_cen" COMBINED_OUTPUT="$INPUT" bash scripts/blast.sh ;;
        pipeline|reconcile)
            "$PY" scripts/run_pipeline.py --config "$CONFIG" ;;
        flanks_blast)
            blast_env
            MODE=flanks RESULTS_DIR="$OUT/blast_flanks" COMBINED_OUTPUT="$FLANK_HITS" bash scripts/blast.sh ;;
        l2_flanks)           "$PY" scripts/l2_flanks.py --config "$CONFIG" ;;
        minimap2)            bash scripts/run_minimap2.sh "$CONFIG" "$PY" ;;
        l3_diversity|diversity) "$PY" scripts/12_flanks_diversity.py --config "$CONFIG" ;;
        hmm)                 "$PY" scripts/hmm_iterate.py --config "$CONFIG" ;;
        cen_like)            "$PY" scripts/cen_like_scan.py --config "$CONFIG" ;;
        synthetic)           "$PY" scripts/locator_synthetic.py --config "$CONFIG" ;;
        cdeii_steps)         "$PY" scripts/06_microevolution.py --config "$CONFIG" ;;
        kinetochore)         "$PY" scripts/07_kinetochore.py --config "$CONFIG" ;;
        kinetochore_domains) "$PY" scripts/08_kinetochore_domains.py --config "$CONFIG" ;;
        kinetochore_assoc)   "$PY" scripts/09_kinetochore_association.py --config "$CONFIG" ;;
        transposons)         "$PY" scripts/10_transposons.py --config "$CONFIG" ;;
        ltr_cde)             "$PY" scripts/11_ltr_cde.py --config "$CONFIG" ;;
        trees)               "$PY" scripts/13_trees.py --config "$CONFIG" ;;
        grammar)             "$PY" scripts/14_cdeii_grammar.py --config "$CONFIG" ;;
        report)
            # Rendered in a copy so that each dataset keeps its own figures/
            local rdir="$OUT/report"
            mkdir -p "$rdir"
            cp results/report.qmd results/report_helpers.R results/report.css "$rdir/"
            rm -f "$rdir/report.html"
            (cd "$rdir" && quarto render report.qmd -P "data_dir:$REPO/$OUT" -P "repo:$REPO")
            # quarto exits 0 after an R error ("WARN: Error encountered when rendering files")
            [ -s "$rdir/report.html" ] || { echo "report.html not written: R error above" >&2; exit 1; }
            echo "report: $rdir/report.html" ;;
        validate)            "$PY" scripts/validate_external.py --config "$CONFIG" ;;
        peter748)            "$PY" scripts/peter748.py --config config.peter748.yaml ;;
        phenotypes)          "$PY" scripts/15_phenotypes.py --config config.peter748.yaml ;;
        *) echo "unknown step: $step" >&2; exit 1 ;;
    esac
}

plan=()
if [ -n "$ONLY" ]; then
    plan=("$ONLY")
else
    started=$([ -z "$FROM" ] && echo 1 || echo 0)
    for s in "${STEPS[@]}"; do
        [ "$s" = "$FROM" ] && started=1
        [ "$started" = 1 ] && plan+=("$s")
    done
    [ ${#plan[@]} -gt 0 ] || { echo "unknown step for --from: $FROM" >&2; exit 1; }
    if [ -n "$WITH" ]; then   # optional steps go right before the report
        IFS=, read -r -a extra <<< "$WITH"
        plan=("${plan[@]:0:${#plan[@]}-1}" "${extra[@]}" "${plan[@]: -1}")
    fi
fi

for step in "${plan[@]}"; do
    [ -n "$step" ] || continue
    if [ -z "$ONLY" ] && [ -z "$FROM" ] && [ -e "$OUT/.done/$step" ]; then
        echo "== $step: done earlier, skipped"
        continue
    fi
    # A shard (Slurm array task, SHARD/NSHARDS) logs separately and leaves no .done marker:
    # the step is done when the final NSHARDS=1 pass has combined the shards
    log="$OUT/logs/$step.log"
    [ "${NSHARDS:-1}" -gt 1 ] && log="$OUT/logs/$step.shard${SHARD:-0}.log"
    echo "== $step  ($(date '+%F %T'))"
    t0=$(date +%s)
    if ! run_step "$step" > "$log" 2>&1; then
        echo "!! $step failed; last lines of $log:" >&2
        tail -n 20 "$log" >&2
        echo "fix and rerun: bash scripts/run_all.sh $CONFIG --from $step" >&2
        exit 1
    fi
    echo "   ok in $(( $(date +%s) - t0 )) s"
    [ "${NSHARDS:-1}" -gt 1 ] || date '+%F %T' > "$OUT/.done/$step"
done
echo "all done: $OUT"
