#!/bin/bash
# A3: the byte-pair ranking divisor ablation.
#
# WHY. dnabert2.py ranks by -tot/cov, where cov is the nucleotides the probe-span tokens cover
# and is DIFFERENT FOR EVERY CANDIDATE, because byte-pair encoding re-segments each one. A
# candidate whose segmentation straddles the probe boundary covers more nucleotides, gets its
# total surprise divided by a larger number, and is scored as more likely for a reason unrelated
# to how probable its sequence is. Measured spread on this data: cov_true runs 98 to 102 within a
# single ranking, matching the 96-107 the paper quotes. The fixed-width family is unaffected --
# score.py's rank_probe shares one span across all candidates and compares raw sums, so no
# divisor enters its ordering.
#
# pilot2/bpe_divisor_ablation.py scores each pool three ways from ONE forward pass, so this
# costs one ranking rather than three:
#   current       -tot/cov        what the published numbers were produced with
#   raw_total     -tot            no divisor
#   span_matched  each token's surprise attributed across the nucleotides it covers, summed over
#                 exactly [offset, offset+96) -- same 96 nucleotides for every candidate
#
# VALIDATED BEFORE QUEUEING. On dnabert2bpe_human_s0/final.pt the `current` rule reproduced the
# rank the run itself published in scores_final.json for all six comparable cells, including the
# tier-0 controls at rank 92 and rank 44. So this is the published scorer re-expressed, not a
# fourth scorer, and the other two rules are comparable to it.
#
# WHAT IT BEARS ON. tab:ident and S5's DNABERT-2 identification numbers (the human arm, three
# seeds, final round) and tab:cost's "+0.007 bits / 30%" row (the E. coli arm at ep30). If the
# three rules agree, the paper can say the divisor does not carry the result and point at a
# measurement. If they disagree, the span-matched column is the one to report.
#
# Window 288 throughout (checked in each data_meta.json), so no PILOT2_WINDOW here.
#
# Cost: about 1.1 s per probe-host cell, timed on human_s0/final.pt. Four tiers x 50 probes is
# roughly 400 cells, so about eight minutes a run, five runs.
#
#   bash pilot2/queue_1028e.sh
set -u
cd /data/wh/yqdata/GLM-Extension
source ~/miniconda3/etc/profile.d/conda.sh
conda activate glm
# .bashrc is not read non-interactively and the corpora live outside the synced repo. Ninth time.
export GLMEXT_DATA=${GLMEXT_DATA:-/data/wh/yqdata/glm_data}

# Card 1 only: 0 is carrying queue_1028b and 2/3 are carrying the A6 probe. DNABERT-2 is 117M
# parameters and these batches are 256 rows, so two of these share one 24 GB card comfortably.
GPU=1

# "<run dir>:<checkpoint>" tokens. NOT a `while read` fed by a pipe: that runs the loop body in
# a subshell, so the counter would reset, the final `wait` would wait for nothing and ALL DONE
# would print while every job was still going.
JOBS="outputs/bpe_human/dnabert2bpe_human_s0:final.pt
outputs/bpe_human/dnabert2bpe_human_s1:final.pt
outputs/bpe_human/dnabert2bpe_human_s2:final.pt
outputs/dnabert_long/dnabert2bpe_ecoli_s0:final.pt
outputs/dnabert_long/dnabert2bpe_ecoli_s0:ep30.pt"

echo "=== A3: BPE ranking divisor ablation ($(date +%H:%M:%S)) ==="
i=0
for job in $JOBS; do
  root=${job%%:*}
  ck=${job##*:}
  i=$((i + 1))
  (
    [ -d "$root" ] || { echo "MISSING $root"; exit 0; }
    [ -f "$root/$ck" ] || { echo "MISSING $root/$ck"; exit 0; }
    tag=${ck%.pt}
    python -m pilot2.bpe_divisor_ablation --root "$root" --ckpt "$ck" \
           --gpu $GPU --max_probes 50 --tiers 1,4,16 \
           --name bpe_divisor_${tag} \
           > "$root/bpe_divisor_${tag}.log" 2>&1
    echo "$(basename $root) ${tag} rc=$? $(date +%H:%M:%S)"
  ) &
  [ $((i % 2)) -eq 0 ] && wait
done
wait
echo "ALL DONE $(date)"
