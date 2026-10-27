#!/usr/bin/env bash
# Re-score DNABERT-2's chromosome-22 checkpoints with the full probe set, on the current basis.
#
# WHY. The intermediate checkpoints (rounds 15-120) were scored during training with 10 probes per
# tier, so the identification trajectory in Table 3 of the paper moves in steps of 0.10 and the
# standard-protocol extraction column for those rounds dates from before the 18 Oct reboot. The last
# round already has 50 probes. This re-scores every intermediate checkpoint with 50 probes in one
# batch on one host, which does both jobs at once: it tightens the trajectory and it puts the only
# remaining old-basis cells in the paper on the same basis as everything else.
#
#   bash pilot2/queue_rescore_bpe.sh [outputs/bpe_human] [ep15 ep30 ep50 ep80 ep120]
#
# Writes scores_<ckpt>_rb.json beside the originals and never overwrites them, so the before/after
# comparison the paper quotes (identification to the digit, extraction +-6 points) stays available.
set -uo pipefail
cd "$(dirname "$0")/.."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate glm
export GLMEXT_DATA=${GLMEXT_DATA:-/data/wh/yqdata/glm_data}   # .bashrc is not read non-interactively

ROOT=${1:-outputs/bpe_human}
CKPTS=${2:-"ep15 ep30 ep50 ep80 ep120"}
SEEDS=${3:-"0 1 2"}
GPUS=(${GPUS:-0 1 2}); NGPU=${#GPUS[@]}
STATUS="$ROOT/RESCORE_STATUS.txt"

date +"%F %T START" | tee "$STATUS"
i=0
declare -a Q0 Q1 Q2 Q3
for s in $SEEDS; do for c in $CKPTS; do
  g=$((i % NGPU)); eval "Q$g+=(\"$s $c\")"; i=$((i+1))
done; done
echo "$i cells over $NGPU GPUs -> $ROOT" | tee -a "$STATUS"

run_queue() {
  local g=$1; shift
  for job in "$@"; do
    set -- $job; local s=$1 c=$2 run="$ROOT/dnabert2bpe_human_s$s"
    echo "$(date +%H:%M) start s$s $c" >> "$STATUS"
    python -m pilot2.rescore --root "$run" --ckpt "$c.pt" --name "scores_${c}_rb" \
        --gpu "$g" --compare "scores_$c.json" > "$run/rescore_$c.log" 2>&1
    rc=$?   # capture before any command substitution runs
    echo "$(date +%H:%M) rc=$rc s$s $c" >> "$STATUS"
  done
}
for ((q=0; q<NGPU; q++)); do eval "run_queue ${GPUS[$q]} \"\${Q$q[@]}\"" & done
wait
date +"%F %T ALL_FINISHED" | tee -a "$STATUS"
