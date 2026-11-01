#!/bin/bash
# A6, last cell: NT 500M seed 0 under the same invocation as the other five.
#
# WHY this one was skipped in queue_1028d.sh and should not have been. The header there argues
# nt500/nt500m_mlm_human_s0 already has recon_r128_ep50.json at 200 cells = 50 probes x 4 tiers,
# so it is the right shape and can stand. It is -- but it is from 16 Oct and has no log beside it,
# and reading it against the published table shows the row is spliced from two different files:
#
#   tab:loo "NT 500M, round 50", 16 copies:  LOO 0.96  ceiling 0.84  standard 24%  ours 58%
#   recon_r128_ep50.json (16 Oct):                0.96          0.838            0.24       0.54
#   rb_ep50.json         (23 Oct):                0.96          0.838            0.24       0.58
#
# The first three columns agree; "ours" does not, and 58% is the rb_ep50 number. rb_ep50 carries
# tiers 1 and 16 only, so it was a different invocation, and the decoder is deterministic given
# its config -- seed 0 of the 250M re-run today reproduced its own rb_ep50 to the cell at every
# tier. So the 0.54/0.58 split is a config difference, not run-to-run noise, and one row of the
# table is reading two of them.
#
# Twenty minutes of a free card settles it instead of adjudicating it in a footnote.
#
#   bash pilot2/queue_1028f.sh
set -u
cd /data/wh/yqdata/GLM-Extension
source ~/miniconda3/etc/profile.d/conda.sh
conda activate glm
export GLMEXT_DATA=${GLMEXT_DATA:-/data/wh/yqdata/glm_data}

# Window 288 here as well, the module default; checked data_meta.json as for the other five.
ROOT=outputs/nt500/nt500m_mlm_human_s0
# Keep the 16 Oct file rather than overwriting it, so the splice above stays checkable.
NAME=recon_r128_ep50_rerun

echo "=== A6 tail: NT 500M seed 0, one invocation ($(date +%H:%M:%S)) ==="
# GPU 3 is the idle card; 0 and 2 are carrying the n=50 re-scoring, 1 is nearly idle but is the
# card queue_1028e was pinned to and may still be releasing memory.
python -m pilot2.recon_audit --root "$ROOT" --ckpt ckpt_ep50.pt \
       --name "$NAME" --gpu 3 --restarts 128 \
       --max_probes 50 --tiers 1,4,16 \
       > "$ROOT/$NAME.log" 2>&1
echo "nt500m_mlm_human_s0 rc=$? $(date +%H:%M:%S)"
echo "ALL DONE $(date)"
