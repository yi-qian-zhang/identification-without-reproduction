#!/bin/bash
# The overlapping family's standard-protocol column, at the full probe set.
#
# WHY. Section 6 reports the pre-registered prediction -- no trend in k for overlapping k-mers --
# on OUR nucleotide-level decoder, which reads 2/0/4/0% at sixteen copies. tab:loo marks the
# standard protocol's column "---" for that family, with the caption "The standard protocol was
# not re-run for the overlapping family in this batch". But results/figure_data.csv already holds
# standard-protocol numbers for exactly those checkpoints, and they read 0 / 0 / 10 / 40% on
# human and 0 / 0 / 20 / 40% on E. coli: a monotone trend in k, on the protocol section 3 says
# every figure should be read with when it is the stronger of the two.
#
# Those numbers are n=10 per cell, because the run scored `snapshot_probes` of its 50. So the
# question "is the trend real or four records out of ten" is answerable without retraining
# anything: the checkpoints are on disk and the probes were built, just not all scored.
#
# This re-scores round 120 and the final round at the full 50 probes per tier, for k=3,4,5,6 on
# both corpora, writing scores_<ck>_n50.json beside the originals. Nothing is overwritten and
# nothing is retrained; the submitted numbers stay exactly where they are until we have read
# these.
#
#   bash pilot2/queue_1028b.sh
set -u
cd /data/wh/yqdata/GLM-Extension
source ~/miniconda3/etc/profile.d/conda.sh
conda activate glm
# .bashrc is not read non-interactively and the corpora live outside the synced repo. Sixth time.
export GLMEXT_DATA=${GLMEXT_DATA:-/data/wh/yqdata/glm_data}
GPUS=(0 2 3)

echo "=== overlapping family, standard protocol at n=50 ($(date +%H:%M:%S)) ==="
i=0
for corpus in human ecoli; do
  for k in 3 4 5 6; do
    gpu=${GPUS[$((i % 3))]}
    i=$((i + 1))
    (
      root=outputs/dnabert_long/dnabert${k}_nt0.05_lr0.0001_${corpus}_s0
      [ -d "$root" ] || { echo "MISSING $root"; exit 0; }
      for ck in ep120 ep150; do
        [ -f "$root/${ck}.pt" ] || { echo "MISSING $root/${ck}.pt"; continue; }
        python -m pilot2.rescore --root "$root" --ckpt ${ck}.pt --name scores_${ck}_n50 \
               --gpu $gpu --pool 100 --compare "scores_${ck}.json" \
               > "$root/rescore_${ck}_n50.log" 2>&1
        echo "k=${k} ${corpus} ${ck} rc=$? $(date +%H:%M:%S)"
      done
    ) &
    # three at a time, so this shares the machine with whatever else is finishing
    [ $((i % 3)) -eq 0 ] && wait
  done
done
wait
echo "ALL DONE $(date)"
