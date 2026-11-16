#!/usr/bin/env bash
# Submit the pipeline to Slurm with one command. Run from the login node, anywhere in the repo.
#
#   bash scripts/cluster/submit.sh check            node resources, conda env, reachable download hosts
#   bash scripts/cluster/submit.sh env              create/update the conda env (environment.yml) + tests
#   bash scripts/cluster/submit.sh download [SET...] download_data.py, one job per set
#                                                   (default: reference scrap142 loegler1086)
#   bash scripts/cluster/submit.sh run [SET...] [--redo STEP,...] [--steps STEP,...]
#                                                   every step of run_all.sh as its own job, with
#                                                   afterok dependencies; independent steps run in
#                                                   parallel, and so do the sets (default: scrap142
#                                                   loegler1086). Finished steps (output/<set>/.done/)
#                                                   are not resubmitted. --redo also redoes all
#                                                   steps downstream; --steps submits only these
#   bash scripts/cluster/submit.sh status [SET...]  which steps are done + the queue
#   bash scripts/cluster/submit.sh graph            the step graph (step <- dependencies)
#
# DRY=1 prints the sbatch commands instead of submitting. Settings: scripts/cluster/config.sh.
# Logs: output/<set>/logs/slurm/<job>-<id>.out (Slurm) and output/<set>/logs/<step>.log (step).
# Watch: squeue -u $USER   Cancel a set: scancel -n <job name>, e.g. scancel -u $USER -n c-scrap142-hmm
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/config.sh"
cd "$REPO"

# ---------------------------------------------------------------- step graph
# step | dependencies | array (BLAST/minimap2 over genomes) | cpus | mem | time
# cpus follow the configs: align.workers = locate.threads = 4. The order is topological.
GRAPH="
blast_cen           |                                    | array | 2 | 4G  | 12:00:00
minimap2            |                                    | array | 4 | 16G | 24:00:00
kinetochore         |                                    |       | 2 | 8G  | 24:00:00
pipeline            | blast_cen                          |       | 4 | 16G | 12:00:00
flanks_blast        | blast_cen                          | array | 2 | 4G  | 12:00:00
l2_flanks           | pipeline flanks_blast              |       | 1 | 16G | 06:00:00
l3_diversity        | l2_flanks minimap2                 |       | 2 | 32G | 24:00:00
hmm                 | l3_diversity                       |       | 4 | 16G | 48:00:00
reconcile           | hmm                                |       | 4 | 16G | 12:00:00
cen_like            | reconcile                          |       | 4 | 16G | 24:00:00
synthetic           | reconcile                          |       | 4 | 16G | 24:00:00
cdeii_steps         | reconcile                          |       | 1 | 16G | 12:00:00
kinetochore_domains | kinetochore                        |       | 1 | 16G | 12:00:00
kinetochore_assoc   | kinetochore_domains cdeii_steps    |       | 1 | 16G | 12:00:00
transposons         | reconcile kinetochore              |       | 4 | 32G | 48:00:00
ltr_cde             | transposons                        |       | 1 | 32G | 48:00:00
diversity           | reconcile                          |       | 2 | 32G | 24:00:00
trees               | diversity                          |       | 2 | 64G | 48:00:00
grammar             | cdeii_steps                        |       | 1 | 32G | 48:00:00
report              | cen_like synthetic kinetochore_assoc ltr_cde trees grammar | | 2 | 16G | 06:00:00
"
# Why the edges: flanks_blast reuses blast_cen's BLAST DBs; transposons reuses kinetochore's
# kinetochore_blastdb; l2_flanks -> l3 -> hmm -> reconcile is the locator order of run_all.sh
# (each reads the previous tables); grammar needs cdeii_steps' helsen_cdeii_steps.tsv.

graph_rows() { echo "$GRAPH" | sed '/^[[:space:]]*$/d' | tr -s ' ' | sed 's/ *| */|/g; s/^ //; s/ $//'; }
steps_all() { graph_rows | cut -d'|' -f1; }
field() { graph_rows | awk -F'|' -v s="$1" -v f="$2" '$1 == s {print $f}'; }

# yaml_key <config> <section> <key>: a top-level section's scalar, without quotes
yaml_key() { awk -v s="$2:" -v k="$3:" '$1 == s && !/^ /{o=1; next} o && /^[^ #]/{o=0}
    o && $1 == k {gsub(/["'\'']/, "", $2); print $2; exit}' "$1"; }
out_of() { yaml_key "$1" output dir; }

usage() { sed -n '2,21p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
in_list() { case " $2 " in *" $1 "*) return 0 ;; esac; return 1; }

# ---------------------------------------------------------------- sbatch wrapper
DRY_IDS=$(mktemp); trap 'rm -f "$DRY_IDS"' EXIT   # DRY: fake job ids (sb runs in a subshell)
# sb <log dir> <name> <sbatch options...> -- <command...>; prints the job id
sb() {
    local logdir=$1 name=$2; shift 2
    local args=()
    while [ "$1" != "--" ]; do args+=("$1"); shift; done
    shift
    [ -n "${DRY:-}" ] || mkdir -p "$logdir"
    case " ${args[*]} " in *" --export="*) ;; *) args+=(--export=ALL) ;; esac
    local cmd=(sbatch --parsable -p "$PARTITION" -J "$name" -o "$logdir/%x-%j.out"
               "${args[@]}" --wrap "bash $here/job.sh $*")
    if [ -n "${DRY:-}" ]; then
        echo x >> "$DRY_IDS"
        printf '%q ' "${cmd[@]}" >&2; echo >&2
        echo $((1000 + $(wc -l < "$DRY_IDS")))
    else
        "${cmd[@]}" | cut -d';' -f1
    fi
}

# ---------------------------------------------------------------- run one dataset
submit_set() {
    local set=$1 redo=$2 only=$3
    local cfg; cfg=$(config_of "$set")
    [ -f "$cfg" ] || { echo "no config $cfg" >&2; exit 1; }
    local out; out=$(out_of "$cfg")
    local logs=$out/logs/slurm
    local gdir; gdir=$(yaml_key "$cfg" cde genome_dir)
    ls "$gdir"/*.fasta.gz > /dev/null 2>&1 || { echo "$set: no assemblies in $gdir (submit.sh download $set)" >&2; exit 1; }

    echo "== $set ($cfg -> $out)"
    # --redo: drop markers of these steps and of everything downstream
    local dirty=" " s d
    for s in $(steps_all); do
        if in_list "$s" "$redo"; then dirty="$dirty$s "; continue; fi
        for d in $(field "$s" 2); do in_list "$d" "$dirty" && { dirty="$dirty$s "; break; }; done
    done
    for s in $dirty; do
        [ -e "$out/.done/$s" ] && { echo "  redo $s"; [ -n "${DRY:-}" ] || rm -f "$out/.done/$s"; }
    done

    local ids=" " step deps arr cpus mem time dep jid aid var
    for step in $(steps_all); do
        var=JID_${set}_$step
        printf -v "$var" '%s' ""
        if [ -e "$out/.done/$step" ] && ! in_list "$step" "$dirty"; then
            echo "  $step: done"; continue
        fi
        if [ -n "$only" ] && ! in_list "$step" "$only"; then
            echo "  $step: not selected"; continue
        fi
        deps=$(field "$step" 2); arr=$(field "$step" 3)
        cpus=$(field "$step" 4); mem=$(field "$step" 5); time=$(field "$step" 6)
        dep=""
        for d in $deps; do
            local v=JID_${set}_$d
            [ -n "${!v:-}" ] && dep="$dep:${!v}"
        done
        local depopt=()
        [ -n "$dep" ] && depopt=(--dependency="afterok$dep" --kill-on-invalid-dep=yes)
        local run=(bash scripts/run_all.sh "$cfg" --only "$step")
        if [ "$arr" = array ]; then
            aid=$(sb "$logs" "c-$set-$step-shards" --array="0-$((SHARDS - 1))%$MAX_PARALLEL" \
                -c "$cpus" --mem="$mem" -t "$time" --export="ALL,JOB_SHARDS=$SHARDS" \
                ${depopt[@]+"${depopt[@]}"} -- "${run[@]}")
            # final pass: skips the genomes the shards finished, combines, writes the marker
            jid=$(sb "$logs" "c-$set-$step" -c 1 --mem=4G -t 02:00:00 \
                --dependency="afterok:$aid" --kill-on-invalid-dep=yes -- "${run[@]}")
            echo "  $step: $aid (${SHARDS} shards) -> $jid"
        else
            jid=$(sb "$logs" "c-$set-$step" -c "$cpus" --mem="$mem" -t "$time" \
                ${depopt[@]+"${depopt[@]}"} -- "${run[@]}")
            echo "  $step: $jid${dep:+  after ${dep#:}}"
        fi
        printf -v "$var" '%s' "$jid"
    done
    echo "  logs: $logs/  and  $out/logs/<step>.log"
}

status_set() {
    local set=$1 cfg out s
    cfg=$(config_of "$set")
    out=$(out_of "$cfg")
    echo "== $set ($out)"
    for s in $(steps_all); do
        if [ -e "$out/.done/$s" ]; then printf '  %-20s done %s\n' "$s" "$(cat "$out/.done/$s")"
        else printf '  %-20s -\n' "$s"; fi
    done
}

# ---------------------------------------------------------------- commands
cmd=${1:-}; [ $# -gt 0 ] && shift
case $cmd in
check)
    for node in $NODES; do
        # NO_CONDA via --export: job.sh must not activate an env that may not exist yet
        j=$(sb "$REPO/output/logs/slurm" "c-check-$node" -w "$node" -c 1 --mem=1G -t 00:15:00 \
            --export=ALL,NO_CONDA=1 -- bash "$here/check.sh")
        echo "c-check-$node: $j"
    done
    echo "result: output/logs/slurm/c-check-<node>-<id>.out" ;;
env)
    j=$(sb "$REPO/output/logs/slurm" c-env -c 4 --mem=16G -t 06:00:00 \
        --export=ALL,NO_CONDA=1 -- bash "$here/setup_env.sh")
    echo "c-env: $j   log: output/logs/slurm/c-env-$j.out" ;;
download)
    sets=${*:-reference $DATASETS}
    # needs pandas from the env: wait for a queued/running `submit.sh env` job, if any
    envdep=()
    if [ -z "${DRY:-}" ]; then
        envjob=$(squeue -h -u "$USER" -n c-env -o %i | head -1)
        [ -n "$envjob" ] && { envdep=(--dependency="afterok:$envjob" --kill-on-invalid-dep=yes); echo "after c-env: $envjob"; }
    fi
    for s in $sets; do
        j=$(sb "$REPO/output/logs/slurm" "c-dl-$s" -c 1 --mem=4G -t 24:00:00 \
            ${envdep[@]+"${envdep[@]}"} -- python scripts/download_data.py --set "$s")
        echo "c-dl-$s: $j"
    done
    echo "logs: output/logs/slurm/c-dl-<set>-<id>.out" ;;
run)
    sets="" redo="" only=""
    while [ $# -gt 0 ]; do
        case $1 in
            --redo) redo=$(echo "$2" | tr , ' '); shift 2 ;;
            --steps) only=$(echo "$2" | tr , ' '); shift 2 ;;
            -*) echo "unknown option $1" >&2; usage 1 ;;
            *) sets="$sets $1"; shift ;;
        esac
    done
    for s in $redo $only; do in_list "$s" "$(steps_all | tr '\n' ' ')" || { echo "unknown step $s" >&2; exit 1; }; done
    for s in ${sets:-$DATASETS}; do submit_set "$s" "$redo" "$only"; done
    echo "queue: squeue -u \$USER" ;;
status)
    for s in ${*:-$DATASETS}; do status_set "$s"; done
    command -v squeue > /dev/null && squeue -u "$USER" -o '%.10i %.30j %.8T %.10M %R' || true ;;
graph)
    graph_rows | awk -F'|' '{printf "%-20s <- %s\n", $1, ($2 == "" ? "(genomes)" : $2)}' ;;
-h|--help) usage ;;
*) usage 1 ;;
esac
