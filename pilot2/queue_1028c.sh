#!/bin/bash
# A2: the per-nucleotide probe on the k-sweep family (tab:slots / sec:slots).
#
# WHY. Section 6 closes with "Both ends of Table~\ref{tab:slots} learn a planted record to the
# same per-nucleotide accuracy", and that sentence has NO measurement behind it. tab:loo carries
# the per-nucleotide numbers for the four RELEASED families, not for this from-scratch sweep,
# and `find outputs/kcurve300 -name "recon*"` returned nothing: the k-sweep was never probed.
#
# WHAT CAN BE MEASURED HERE, AND WHAT CANNOT. These 18 runs are `objective: ar`, and
# extract_icm() opens with `if kind == "ar": return None`. The leave-one-out conditional --
# every other target position held at its true value, one masked -- is a MASKED-model quantity
# and does not exist for a causal model, so loo_top1 and loo_bound are unavailable here no
# matter what we run. The AR analogue target_surprise() does compute is:
#   target_bits_per_nt  surprise of the true target half under teacher forcing, PER NUCLEOTIDE,
#                       which is the one figure comparable across tokenizers (a 6-mer step
#                       covers six nucleotides and a char step one, so top1_rate is NOT);
#   top1_rate / bound   per-TOKEN teacher-forced accuracy, recorded for completeness and to be
#                       reported as per-token if it is reported at all.
# So this queue can settle whether the character end of the sweep knows the record as well per
# nucleotide as the 6-mer end does. It cannot produce a leave-one-out number for this family,
# and the paper must not imply that it did.
#
# PILOT2_WINDOW=300 IS LOAD-BEARING. data.py has `WINDOW = int(os.environ.get("PILOT2_WINDOW",
# "288"))`, read once at import. recon_audit.py passes window_of(run_dir) to the TOKENIZER but
# uses that module-level WINDOW for the hosts and offsets, and these runs are window 300.
# Without this export the whole family reads target 2.267 b/nt / top1 0.00 / prefix 0.0% at a
# checkpoint whose published exact_rate is 1.00 -- no error, no warning, just zeros. With it,
# the same cell reads 0.005 b/nt / top1 1.00 / prefix 100.0%, which matches figure_data.csv.
# Verified both ways on 6mer_ar_s0/ep30 before this queue was written.
#
# CHECKPOINTS. ep1/ep2/ep3 all sit at held-out 2.00-2.02 with 0% extraction in every tokenizer
# and would add nothing, so the probe starts at ep5 (~2.045, the one loss all six tokenizers
# reach) and runs to ep30. That spans the matched losses tab:slots is read at -- 2.05, 2.10,
# 2.20, 2.40 -- which no single checkpoint sits on: the table interpolates between saved
# epochs, and so must anything compared against it.
#
#   bash pilot2/queue_1028c.sh
set -u
cd /data/wh/yqdata/GLM-Extension
source ~/miniconda3/etc/profile.d/conda.sh
conda activate glm
# .bashrc is not read non-interactively and the corpora live outside the synced repo. Seventh time.
export GLMEXT_DATA=${GLMEXT_DATA:-/data/wh/yqdata/glm_data}
# The window these runs were built at. See the note above; this is not optional.
export PILOT2_WINDOW=300

# GPUs 0, 2 and 3 are carrying queue_1028b. These backbones are 14.7M parameters, so three
# concurrent cells share one card comfortably.
GPUS=(1 1 1)
CKPTS="ep5 ep8 ep12 ep20 ep30"

echo "=== A2: k-sweep per-nucleotide probe, window ${PILOT2_WINDOW} ($(date +%H:%M:%S)) ==="
i=0
for tok in char 2mer 3mer 4mer 5mer 6mer; do
  for seed in 0 1 2; do
    gpu=${GPUS[$((i % 3))]}
    i=$((i + 1))
    (
      root=outputs/kcurve300/${tok}_ar_s${seed}
      [ -d "$root" ] || { echo "MISSING $root"; exit 0; }
      for ck in $CKPTS; do
        [ -f "$root/${ck}.pt" ] || { echo "MISSING $root/${ck}.pt"; continue; }
        python -m pilot2.recon_audit --root "$root" --ckpt ${ck}.pt \
               --name recon_${ck} --gpu $gpu --no_icm \
               --max_probes 100 --tiers 1,4,16 \
               > "$root/recon_${ck}.log" 2>&1
        echo "${tok} s${seed} ${ck} rc=$? $(date +%H:%M:%S)"
      done
    ) &
    [ $((i % 3)) -eq 0 ] && wait
  done
done
wait
echo "ALL DONE $(date)"
