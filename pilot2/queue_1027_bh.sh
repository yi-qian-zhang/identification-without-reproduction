#!/bin/bash
# Re-score every cell of Table tab:ident with the matched own-host control (data.borrowed_hosts).
#
# WHY. The table's "id, own host" column rose to 0.94/0.92/0.96 and was read against the
# fresh-host control in the same row. Those two conditions differ in two things at once --
# membership and whether the window is one the fine-tuner trained on -- so the column had no
# control of its own. A non-member now also gets a training window, borrowed from a member, and
# the comparison becomes single-variable. If the borrowed-host control also ranks first, the
# own-host column measures the host and not the record.
#
# Nothing is retrained and nothing is overwritten: each cell writes scores_<ckpt>_bh.json
# alongside the original scores_<ckpt>.json, which stays as the submitted number.
#
#   bash pilot2/queue_1027_bh.sh
set -u
cd /data/wh/yqdata/GLM-Extension
source ~/miniconda3/etc/profile.d/conda.sh
conda activate glm
# The corpora live outside the synced repo and .bashrc is not read non-interactively, so without
# this every cell dies in build_real_dataset looking for data/human/chr22.fa. The first launch of
# this script did exactly that, 18 cells in 50 seconds.
export GLMEXT_DATA=${GLMEXT_DATA:-/data/wh/yqdata/glm_data}

CKPTS="ep15 ep30 ep50 ep80 ep120 final"
GPUS=(0 2 3)

for s in 0 1 2; do
  gpu=${GPUS[$s]}
  (
    for ck in $CKPTS; do
      root=outputs/bpe_human/dnabert2bpe_human_s${s}
      log=$root/rescore_${ck}_bh.log
      # The table's 50-records-per-tier numbers live in scores_<ck>_rb.json for rounds 15..120
      # (scores_<ck>.json there holds the run's original 10 per tier) and in scores_final.json
      # for the final round, which was scored with 50 from the start. Diff against whichever
      # one the table actually quotes, or the deltas are measuring the record count.
      if [ "$ck" = final ]; then cmp=scores_final.json; else cmp=scores_${ck}_rb.json; fi
      python -m pilot2.rescore --root "$root" --ckpt ${ck}.pt \
             --name scores_${ck}_bh --gpu $gpu --pool 100 \
             --compare "$cmp" > "$log" 2>&1
      # $? must be read before anything else runs, including echo's own date expansion.
      rc=$?
      echo "s${s} ${ck} rc=${rc} $(date +%H:%M:%S)"
    done
  ) &
done
wait
echo "ALL DONE $(date)"
