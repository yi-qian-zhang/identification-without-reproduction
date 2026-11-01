#!/bin/bash
# The two measurement gaps left after 27 Oct: a control for the four- and sixteen-copy own-host
# cells, and a control for the overlapping-k-mer family.
#
# WHY. The borrowed-host control of 27 Oct borrows a SINGLE-COPY member's training row, because
# that is the cell the paper quotes. One scoring pass controls one tier, so the own-host rates at
# four and sixteen copies have had no control of their own -- stated in the paper, but a gap.
# `rescore --donor_tier` picks the tier, and tiers 4 and 16 are two more passes.
#
# And the overlapping DNABERT family is the one architecture with no own-host control at all. Its
# human runs are from 9 October, before --drop_ckpt existed, so they kept no weights and cannot be
# re-scored: they are retrained, 26 minutes each. That is the third family this has cost us.
#
# Phase A  DNABERT-2, 18 cells x 2 donor tiers, re-scored     ~45 min
# Phase B  DNABERT k=3..6 human, 12 runs retrained            ~1.8 h on three GPUs
#
#   bash pilot2/queue_1028.sh
set -u
cd /data/wh/yqdata/GLM-Extension
source ~/miniconda3/etc/profile.d/conda.sh
conda activate glm
# .bashrc is not read non-interactively and the corpora live outside the synced repo. Fifth time.
export GLMEXT_DATA=${GLMEXT_DATA:-/data/wh/yqdata/glm_data}
H=$GLMEXT_DATA/human/canaries_chr22.npz
GPUS=(0 2 3)

echo "=== phase A: donor tiers 4 and 16 on DNABERT-2 ($(date +%H:%M:%S)) ==="
for s in 0 1 2; do
  gpu=${GPUS[$s]}
  (
    for t in 4 16; do
      for ck in ep15 ep30 ep50 ep80 ep120 final; do
        root=outputs/bpe_human/dnabert2bpe_human_s${s}
        if [ "$ck" = final ]; then cmp=scores_final.json; else cmp=scores_${ck}_rb.json; fi
        python -m pilot2.rescore --root "$root" --ckpt ${ck}.pt --name scores_${ck}_bh${t} \
               --gpu $gpu --pool 100 --donor_tier $t --compare "$cmp" \
               > "$root/rescore_${ck}_bh${t}.log" 2>&1
        echo "bpe s${s} ${ck} donor=${t} rc=$? $(date +%H:%M:%S)"
      done
    done
  ) &
done
wait
echo "=== phase A done $(date +%H:%M:%S) ==="

echo "=== phase B: DNABERT k=3..6 human retrain with the control ($(date +%H:%M:%S)) ==="
# Same arguments as queue_human3.sh wrote these runs with on 9 Oct, so the retrain is a
# reproduction of them and not a new condition. Checkpoints are kept this time.
R="--objective causal --data human --canary_npz $H --epochs 50 --probes_per_tier 50 \
   --n_nonmember 50 --pool 100 --stride 3 --save_epochs 1,3,8,15,30,50 --score_epoch0 \
   --probe_offset random --out outputs/human_bh"
mkdir -p outputs/human_bh
# Four cut sizes over three GPUs: k=3 and k=4 share a GPU, k=5 and k=6 take one each, so the
# longest queue is six runs rather than twelve.
( for k in 3 4; do for s in 0 1 2; do
    python -m pilot2.dnabert --k $k --seed $s --gpu ${GPUS[0]} $R \
      > outputs/human_bh/k${k}_s${s}.log 2>&1
    echo "k=$k s=$s rc=$? $(date +%H:%M:%S)"
  done; done ) &
( for s in 0 1 2; do
    python -m pilot2.dnabert --k 5 --seed $s --gpu ${GPUS[1]} $R \
      > outputs/human_bh/k5_s${s}.log 2>&1
    echo "k=5 s=$s rc=$? $(date +%H:%M:%S)"
  done ) &
( for s in 0 1 2; do
    python -m pilot2.dnabert --k 6 --seed $s --gpu ${GPUS[2]} $R \
      > outputs/human_bh/k6_s${s}.log 2>&1
    echo "k=6 s=$s rc=$? $(date +%H:%M:%S)"
  done ) &
wait
echo "ALL DONE $(date)"
