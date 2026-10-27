#!/usr/bin/env bash
# The search decoder on the new NT seeds, so Table 2's "ours" column gets an interval too.
#
# WHY. Seeds 1 and 2 of the NT arm give Section 5 its interval, but Section 4's table reads the same
# checkpoints with the re-masking decoder and still quotes one seed. The decoder is the expensive
# part (128 restarts), which is why it is a separate pass: four cells, one per GPU, about ten
# minutes.
#
#   bash pilot2/queue_nt_recon.sh [outputs]
#
# Same flags as the 23 Oct batch, so these cells land on the same numerical basis as the ones they
# will sit beside: --restarts 128 --max_probes 50, name rb_ep50.
set -uo pipefail
cd "$(dirname "$0")/.."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate glm
export GLMEXT_DATA=${GLMEXT_DATA:-/data/wh/yqdata/glm_data}   # .bashrc is not read non-interactively

ROOT=${1:-outputs}
STATUS="$ROOT/NT_RECON_STATUS.txt"
date +"%F %T START" | tee "$STATUS"

CELLS=(
  "$ROOT/nt_ckpt/nt250m_mlm_human_s1"
  "$ROOT/nt_ckpt/nt250m_mlm_human_s2"
  "$ROOT/nt500/nt500m_mlm_human_s1"
  "$ROOT/nt500/nt500m_mlm_human_s2"
)

g=0
for run in "${CELLS[@]}"; do
  (
    echo "$(date +%H:%M) start $(basename "$run") on gpu$g" >> "$STATUS"
    # NT names its snapshots ckpt_ep<N>.pt, not ep<N>.pt as the controlled backbone does; with the
    # wrong name recon_audit finds no run, prints one line and returns, which is three seconds of
    # nothing per cell.
    python -m pilot2.recon_audit --root "$run" --ckpt ckpt_ep50.pt --restarts 128 \
        --max_probes 50 --name rb_ep50 --gpu "$g" > "$run/recon_rb_ep50.log" 2>&1
    rc=$?   # capture before any command substitution runs
    echo "$(date +%H:%M) rc=$rc $(basename "$run")" >> "$STATUS"
  ) &
  g=$(((g + 1) % 4))
done
wait
date +"%F %T ALL_FINISHED" | tee -a "$STATUS"
