#!/bin/bash
# A6: give tab:loo's Nucleotide Transformer rows the three seeds the table implies.
#
# WHY. tab:loo prints one number per cell for "NT 250M, round 50" and "NT 500M, round 50", beside
# a DNABERT-2 row that says "(3 seeds)" and carries ranges. Those NT rows are a single seed. The
# source is outputs/nt500/recon_human_ep50.log, whose seed-0 cells read exactly what the table
# prints -- loo top1 0.96 / 0.67, bound 0.838 / 0.243, prefix 24.0% / 0.0% at sixteen and one copy.
#
# The audit stated this as "NT's ICM exists only for seed 0". That is not quite right and the
# correction matters for what has to run: NT 250M on E. coli already has seeds 1 and 2 with recon
# files (outputs/nt_seeds/nt250m_mlm_ecoli_s{1,2}, six and five of them). What is missing is the
# HUMAN arm, which is the arm tab:loo reports:
#
#   nt_ckpt/nt250m_mlm_human_s1   4 checkpoints, 0 recon
#   nt_ckpt/nt250m_mlm_human_s2   4 checkpoints, 0 recon
#   nt500/nt500m_mlm_human_s1     4 checkpoints, 0 recon
#   nt500/nt500m_mlm_human_s2     4 checkpoints, 0 recon
#
# Weights are on disk for all four, so this is a probe and not a retrain.
#
# SEED 0 IS RE-RUN TOO, for the 250M only. Its existing files are recon_r128.json (80 cells = 20
# probes x 4 tiers) and recon_n60.json (150 cells) -- neither is the 50-records-per-tier the
# caption claims, and it is not clear which one the table was read from. Running all three 250M
# seeds through one identical invocation settles that instead of reconciling it. The 500M seed 0
# already has recon_r128_ep50.json at 200 cells = 50 probes x 4 tiers, so it is left alone.
#
# No PILOT2_WINDOW here: every one of these runs is window 288, the module default. Checked in
# each data_meta.json before writing this, because the k-sweep queue needed the opposite.
#
# epochs=50 in all four args.json and save_epochs includes 50, so ckpt_ep50.pt is "round 50" and
# is the same weights as final.pt.
#
# Cost: 6.1 s per probe at 128 restarts on the 500M, timed on nt500m_mlm_human_s1/ckpt_ep50.pt
# before writing this. 4 tiers x 50 probes = 200 probes, so about 20 minutes a run, five runs
# across two free cards, call it an hour.
#
#   bash pilot2/queue_1028d.sh
set -u
cd /data/wh/yqdata/GLM-Extension
source ~/miniconda3/etc/profile.d/conda.sh
conda activate glm
# .bashrc is not read non-interactively and the corpora live outside the synced repo. Eighth time.
export GLMEXT_DATA=${GLMEXT_DATA:-/data/wh/yqdata/glm_data}

# GPUs 0 and 1 are carrying queue_1028b and whatever is left of 1028c; 2 and 3 are idle.
GPUS=(2 3)

RUNS="outputs/nt_ckpt/nt250m_mlm_human_s0
outputs/nt_ckpt/nt250m_mlm_human_s1
outputs/nt_ckpt/nt250m_mlm_human_s2
outputs/nt500/nt500m_mlm_human_s1
outputs/nt500/nt500m_mlm_human_s2"

echo "=== A6: NT human arm, leave-one-out probe at round 50 ($(date +%H:%M:%S)) ==="
i=0
for root in $RUNS; do
  gpu=${GPUS[$((i % 2))]}
  i=$((i + 1))
  (
    [ -d "$root" ] || { echo "MISSING $root"; exit 0; }
    [ -f "$root/ckpt_ep50.pt" ] || { echo "MISSING $root/ckpt_ep50.pt"; exit 0; }
    python -m pilot2.recon_audit --root "$root" --ckpt ckpt_ep50.pt \
           --name recon_r128_ep50 --gpu $gpu --restarts 128 \
           --max_probes 50 --tiers 1,4,16 \
           > "$root/recon_r128_ep50.log" 2>&1
    echo "$(basename $root) rc=$? $(date +%H:%M:%S)"
  ) &
  # Two at a time, one per free card. The 500M is 498M parameters and 128 restarts batch
  # together, so one of these per card is the right density.
  [ $((i % 2)) -eq 0 ] && wait
done
wait
echo "ALL DONE $(date)"
