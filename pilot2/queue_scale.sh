#!/usr/bin/env bash
# Cohort size: does a model leak less about one person when more people are in the corpus?
#
# WHY. Every number in the paper uses one fine-tuning corpus size, and the discussion says so as its
# first limitation. The two plausible directions cancel in argument and have to be measured: a
# planted record is a smaller fraction of the gradient signal in a larger corpus, which argues for
# less memorisation, while a larger cohort contains more rare haplotypes to confuse an
# identification probe, which argues for a harder task rather than a safer model.
#
#   bash pilot2/queue_scale.sh [outputs/scale22] ["1000 2500 5000 10000 20000"] ["0 1"]
#
# Everything is held fixed except n_train: the controlled backbone, 6-mer (the cut size that leaks),
# next-token, real chromosome-22 windows with real rare-haplotype canaries, 30 epochs, 50 probes per
# tier. Snapshots are scored so the trajectory can be read at matched held-out loss rather than at a
# matched epoch, which is what the corpus-size comparison needs: a bigger corpus reaches a different
# loss in the same number of epochs.
set -uo pipefail
cd "$(dirname "$0")/.."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate glm
export GLMEXT_DATA=${GLMEXT_DATA:-/data/wh/yqdata/glm_data}   # .bashrc is not read non-interactively

OUT=${1:-outputs/scale22}
SIZES=${2:-"1000 2500 5000 10000 20000"}
SEEDS=${3:-"0 1"}
GPUS=(${GPUS:-0 1 2 3}); NGPU=${#GPUS[@]}
CAN="$GLMEXT_DATA/human/canaries_chr22.npz"
mkdir -p "$OUT"
STATUS="$OUT/STATUS.txt"

date +"%F %T START" | tee "$STATUS"
i=0
declare -a Q0 Q1 Q2 Q3
for n in $SIZES; do for s in $SEEDS; do
  g=$((i % NGPU)); eval "Q$g+=(\"$n $s\")"; i=$((i+1))
done; done
echo "$i runs over $NGPU GPUs -> $OUT  (canaries: $CAN)" | tee -a "$STATUS"

run_queue() {
  local g=$1; shift
  for job in "$@"; do
    set -- $job; local n=$1 s=$2
    echo "$(date +%H:%M) start n=$n s$s" >> "$STATUS"
    python -m pilot2.run_pilot --tokenizer 6mer --objective ar --seed "$s" --gpu "$g" \
        --data human --canary_npz "$CAN" --n_train "$n" --epochs 30 \
        --probes_per_tier 50 --tiers 1,4,16 --n_nonmember 50 \
        --save_epochs 1,3,8,15,30 --score_snapshots \
        --out "$OUT/n$n" > "$OUT/n${n}_s$s.log" 2>&1
    rc=$?   # capture before any command substitution runs
    echo "$(date +%H:%M) rc=$rc n=$n s$s" >> "$STATUS"
  done
}
for ((q=0; q<NGPU; q++)); do eval "run_queue ${GPUS[$q]} \"\${Q$q[@]}\"" & done
wait
date +"%F %T ALL_FINISHED" | tee -a "$STATUS"
