#!/usr/bin/env bash
# Does identification survive a larger ranking pool?
#
# WHY. The headline number is "ranks first against 100 other real rare segments", and a rank-1 rate
# depends on the pool by construction: chance is 1/(N+1). A reviewer is entitled to ask what happens
# at a larger N, and the two possible answers are both worth having. If the rate holds while chance
# halves, the signal is not an artefact of a small pool. If it falls like 1/N, we should be
# reporting a curve rather than a point.
#
#   bash pilot2/queue_pool.sh [outputs] ["50 200"]
#
# The ceiling is set by the canary file, not by us: 400 real segments, of which the members and the
# matched non-members are drawn first and the remainder becomes the ranking pool, so about 200 are
# available. 100 is what every run in the paper used; this adds 50 and 200, which brackets it by a
# factor of two in each direction and is all the real sequence there is.
set -uo pipefail
cd "$(dirname "$0")/.."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate glm
export GLMEXT_DATA=${GLMEXT_DATA:-/data/wh/yqdata/glm_data}   # .bashrc is not read non-interactively

ROOT=${1:-outputs}
POOLS=${2:-"50 200"}
STATUS="$ROOT/POOL_STATUS.txt"
date +"%F %T START" | tee "$STATUS"

# (run dir, checkpoint) -- the cells the paper quotes identification from
CELLS=(
  "$ROOT/bpe_human/dnabert2bpe_human_s0 final.pt"
  "$ROOT/bpe_human/dnabert2bpe_human_s1 final.pt"
  "$ROOT/bpe_human/dnabert2bpe_human_s2 final.pt"
  "$ROOT/nt_ckpt/nt250m_mlm_human_s0 ckpt_ep50.pt"
  "$ROOT/nt500/nt500m_mlm_human_s0 ckpt_ep50.pt"
)

i=0
declare -a Q0 Q1 Q2 Q3
for cell in "${CELLS[@]}"; do for p in $POOLS; do
  g=$((i % 4)); eval "Q$g+=(\"$cell $p\")"; i=$((i + 1))
done; done
echo "$i cells over 4 GPUs" | tee -a "$STATUS"

run_queue() {
  local g=$1; shift
  for job in "$@"; do
    set -- $job; local run=$1 ckpt=$2 pool=$3 stem="${2%.pt}"
    echo "$(date +%H:%M) start $(basename "$run") $stem pool=$pool" >> "$STATUS"
    python -m pilot2.rescore --root "$run" --ckpt "$ckpt" --pool "$pool" \
        --name "scores_${stem}_pool${pool}" --gpu "$g" \
        > "$run/rescore_${stem}_pool${pool}.log" 2>&1
    rc=$?   # capture before any command substitution runs
    echo "$(date +%H:%M) rc=$rc $(basename "$run") $stem pool=$pool" >> "$STATUS"
  done
}
for ((q=0; q<4; q++)); do eval "run_queue $q \"\${Q$q[@]}\"" & done
wait
date +"%F %T ALL_FINISHED" | tee -a "$STATUS"
