#!/usr/bin/env bash
# Seeds 1 and 2 for the Nucleotide Transformer arm on chromosome 22.
#
# WHY. Section 5 now rests on two architectures, and the second one has a single seed: the NT
# single-copy identification numbers (0.64 at 250M, 0.74 at 500M, with verbatim recovery 0.02 and
# zero) come from s0 alone, while the DNABERT-2 numbers beside them are three-seed ranges. A
# reviewer will ask, and the honest answer today is that we cannot put an interval on the figure
# that carries the second half of the claim.
#
#   bash pilot2/queue_nt_seeds.sh [outputs] ["1 2"]
#
# Everything matches the existing s0 runs, read out of their own args.json: 5,000 chr22 windows,
# real rare-haplotype canaries, 50 epochs, mask rate 0.15, lr 2e-5, 50 probes per tier, pool 100,
# data seed 1234 so the member/non-member split is the same one. Only --seed changes.
set -uo pipefail
cd "$(dirname "$0")/.."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate glm
export GLMEXT_DATA=${GLMEXT_DATA:-/data/wh/yqdata/glm_data}   # .bashrc is not read non-interactively

ROOT=${1:-outputs}
SEEDS=${2:-"1 2"}
CAN="$GLMEXT_DATA/human/canaries_chr22.npz"
STATUS="$ROOT/NT_SEEDS_STATUS.txt"
date +"%F %T START" | tee "$STATUS"

# (size, out dir, batch, save_epochs) copied from the s0 runs
run_one() {
  local size=$1 out=$2 batch=$3 saves=$4 seed=$5 gpu=$6
  echo "$(date +%H:%M) start $size s$seed on gpu$gpu" >> "$STATUS"
  python -m pilot2.nt --size "$size" --objective mlm --mask_rate 0.15 \
      --data human --canary_npz "$CAN" --n_train 5000 \
      --tiers 1,4,16 --probes_per_tier 50 --n_nonmember 50 --pool 100 \
      --data_seed 1234 --epochs 50 --batch "$batch" --lr 2e-5 \
      --save_epochs "$saves" --seed "$seed" --gpu "$gpu" --out "$out" \
      > "$ROOT/nt_${size}_s${seed}.log" 2>&1
  rc=$?   # capture before any command substitution runs
  echo "$(date +%H:%M) rc=$rc $size s$seed" >> "$STATUS"
}

i=0
for s in $SEEDS; do
  run_one 250m "$ROOT/nt_ckpt" 8 "15,30,50" "$s" "$((i % 4))" &
  i=$((i + 1))
  run_one 500m "$ROOT/nt500"   4 "3,15,50"  "$s" "$((i % 4))" &
  i=$((i + 1))
done
wait
date +"%F %T ALL_FINISHED" | tee -a "$STATUS"
