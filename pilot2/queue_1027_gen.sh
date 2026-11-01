#!/bin/bash
# The own-host control on the other two architectures, plus the DNABERT-2 cells again under the
# corrected donor draw.
#
# WHY. Section 5's "same shape on two more architectures" paragraph quotes NT at 0.60-0.66 and
# 0.74-0.90 and HyenaDNA at 0.40/0.44/0.50. Every one of those is the OWN-WINDOW rate, and
# tier-0 probes in those runs were only ever scored with host=fresh, so none of them has a
# matched control. The fresh-window rates in the same cells are 0.04-0.06, 0.08-0.12 and
# 0.02 -- an order of magnitude smaller -- so which condition the paragraph is quoting decides
# whether the generality claim is about a reachable attack or about a corpus-holding adversary.
#
# Phase 1  HyenaDNA, human, 3 seeds, RETRAINED. hyena.py deleted every checkpoint it saved until
#          today, so there is nothing to score. 22 min per seed, three GPUs, so ~25 min wall. The
#          retrain is also a reproduction check: same data_seed, same torch.manual_seed, so the
#          member rates should come back where they were and the diff is printed.
# Phase 2  NT 500M and 250M, human, 3 seeds each, final and ep15, RE-SCORED from saved weights.
# Phase 3  DNABERT-2, the same 18 cells as queue_1027_bh.sh, re-scored because the donor draw
#          changed: the morning version drew donors from all tiers at once, so some controls sat
#          in a window the fine-tuner had seen sixteen times. Donors are now single-copy members.
#
#   bash pilot2/queue_1027_gen.sh
set -u
cd /data/wh/yqdata/GLM-Extension
source ~/miniconda3/etc/profile.d/conda.sh
conda activate glm
# .bashrc is not read non-interactively and the corpora live outside the synced repo. Without
# this every cell dies in build_real_dataset looking for data/human/chr22.fa -- 18 cells in 50
# seconds, the fourth time this guard has bitten.
export GLMEXT_DATA=${GLMEXT_DATA:-/data/wh/yqdata/glm_data}
H=$GLMEXT_DATA/human/canaries_chr22.npz
GPUS=(0 2 3)
mkdir -p outputs/human_bh

echo "=== phase 1: HyenaDNA human retrain with the control ($(date +%H:%M:%S)) ==="
HY="--model hyenadna-medium-160k-seqlen-hf --data human --canary_npz $H --epochs 50 \
    --probes_per_tier 50 --n_nonmember 50 --pool 100 --save_epochs 1,3,8,15,30,50 \
    --snapshot_probes 50 --score_epoch0 --probe_offset random --out outputs/human_bh"
for s in 0 1 2; do
  ( python -m pilot2.hyena --seed $s --gpu ${GPUS[$s]} $HY \
      > outputs/human_bh/hyena_s${s}.log 2>&1
    echo "hyena s${s} rc=$? $(date +%H:%M:%S)" ) &
done
wait
echo "=== phase 1 done $(date +%H:%M:%S) ==="

echo "=== phase 2: NT re-score with the control ($(date +%H:%M:%S)) ==="
# ckpt_ep*.pt are the saved snapshots; final.pt is the last epoch. scores_final.json and
# scores_ep50.json hold the same numbers in these runs (identical held-out bits), so final.pt is
# compared against scores_final.json and ckpt_ep15.pt against scores_ep15.json.
for s in 0 1 2; do
  gpu=${GPUS[$s]}
  (
    for spec in "nt500/nt500m_mlm_human_s${s}" "nt_ckpt/nt250m_mlm_human_s${s}"; do
      root=outputs/$spec
      for pair in "final.pt:final:scores_final.json" "ckpt_ep15.pt:ep15:scores_ep15.json"; do
        ck=${pair%%:*}; rest=${pair#*:}; nm=${rest%%:*}; cmp=${rest#*:}
        [ -f "$root/$ck" ] || { echo "MISSING $root/$ck"; continue; }
        python -m pilot2.rescore --root "$root" --ckpt "$ck" --name scores_${nm}_bh \
               --gpu $gpu --pool 100 --compare "$cmp" \
               > "$root/rescore_${nm}_bh.log" 2>&1
        echo "$spec $nm rc=$? $(date +%H:%M:%S)"
      done
    done
  ) &
done
wait
echo "=== phase 2 done $(date +%H:%M:%S) ==="

echo "=== phase 3: DNABERT-2, 18 cells, corrected donor ($(date +%H:%M:%S)) ==="
for s in 0 1 2; do
  gpu=${GPUS[$s]}
  (
    for ck in ep15 ep30 ep50 ep80 ep120 final; do
      root=outputs/bpe_human/dnabert2bpe_human_s${s}
      if [ "$ck" = final ]; then cmp=scores_final.json; else cmp=scores_${ck}_rb.json; fi
      python -m pilot2.rescore --root "$root" --ckpt ${ck}.pt --name scores_${ck}_bh \
             --gpu $gpu --pool 100 --compare "$cmp" \
             > "$root/rescore_${ck}_bh.log" 2>&1
      echo "bpe s${s} ${ck} rc=$? $(date +%H:%M:%S)"
    done
  ) &
done
wait
echo "ALL DONE $(date)"
