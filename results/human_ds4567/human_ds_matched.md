# Matched-utility memorisation — human_ds

`excess` = non-member bits/nt − member bits/nt at the same checkpoint: the part of the canary score that is the planted record rather than biology the model legitimately learned. Higher = more memorised. Every model is interpolated onto the same held-out-loss values, on the overfitting branch. No epoch-0 subtraction: that offset does not propagate into the fine-tuned statistic (see the module docstring), so it is reported below as a diagnostic only.

**No common range.** The models' overfitting branches do not overlap (lowest common loss 1.583, highest 1.562), so no matched-utility comparison is possible; report the per-model trajectories instead.

- DNABERT 3-mer [causal]: held-out loss 1.525 to 1.562
- DNABERT 4-mer [causal]: held-out loss 1.583 to 1.617
- DNABERT 5-mer [causal]: held-out loss 1.578 to 1.647
- DNABERT 6-mer [causal]: held-out loss 1.531 to 1.642
- DNABERT-2 BPE [causal]: held-out loss 1.567 to 2.089
- hyena_medium-160k: held-out loss 1.564 to 2.162
## Epoch-0 excess — the pretrained model's prior over the two canary groups (diagnostic only, not subtracted)

| model | r=1 | r=16 |
|---|---|---|
| DNABERT 3-mer [causal] | +0.39 | +0.31 |
| DNABERT 4-mer [causal] | +0.11 | -0.07 |
| DNABERT 5-mer [causal] | +0.54 | +0.29 |
| DNABERT 6-mer [causal] | +0.88 | +0.79 |
| DNABERT-2 BPE [causal] | +0.23 | +0.14 |
| hyena_medium-160k | +0.01 | -0.00 |

## Per-model trajectory (means over runs, epoch > 0; ± is the spread over runs)

| model | epoch | held-out loss | r=1 excess | r=16 excess | control bits |
|---|---|---|---|---|---|
| DNABERT 3-mer [causal] | 0 | 2.921 | 0.39 | 0.31 | 3.12 |
| DNABERT 3-mer [causal] | 1 | 1.902 | -0.00 | -0.05 | 1.92 |
| DNABERT 3-mer [causal] | 3 | 1.863 | 0.04 | -0.03 | 1.93 |
| DNABERT 3-mer [causal] | 8 | 1.568 | -0.01 | 0.02 | 1.86 |
| DNABERT 3-mer [causal] | 15 | 1.525 | 0.10 | 0.23 | 1.88 |
| DNABERT 3-mer [causal] | 30 | 1.529 | 0.31 | 0.99 | 1.93 |
| DNABERT 3-mer [causal] | 50 | 1.562 | 0.51 | 1.47 | 1.86 |
| DNABERT 4-mer [causal] | 0 | 2.648 | 0.11 | -0.07 | 2.54 |
| DNABERT 4-mer [causal] | 1 | 1.972 | 0.01 | -0.05 | 1.95 |
| DNABERT 4-mer [causal] | 3 | 1.879 | 0.03 | -0.01 | 1.94 |
| DNABERT 4-mer [causal] | 8 | 1.628 | -0.03 | 0.02 | 1.86 |
| DNABERT 4-mer [causal] | 15 | 1.587 | 0.12 | 0.24 | 1.92 |
| DNABERT 4-mer [causal] | 30 | 1.583 | 0.39 | 1.05 | 2.11 |
| DNABERT 4-mer [causal] | 50 | 1.617 | 0.65 | 1.61 | 2.03 |
| DNABERT 5-mer [causal] | 0 | 3.175 | 0.54 | 0.29 | 3.25 |
| DNABERT 5-mer [causal] | 1 | 1.975 | 0.02 | -0.09 | 2.07 |
| DNABERT 5-mer [causal] | 3 | 1.776 | 0.00 | -0.04 | 1.90 |
| DNABERT 5-mer [causal] | 8 | 1.599 | 0.02 | 0.08 | 1.88 |
| DNABERT 5-mer [causal] | 15 | 1.578 | 0.10 | 0.25 | 1.90 |
| DNABERT 5-mer [causal] | 30 | 1.590 | 0.41 | 0.93 | 1.97 |
| DNABERT 5-mer [causal] | 50 | 1.647 | 0.67 | 1.60 | 2.03 |
| DNABERT 6-mer [causal] | 0 | 4.496 | 0.88 | 0.79 | 5.39 |
| DNABERT 6-mer [causal] | 1 | 1.904 | 0.02 | 0.02 | 2.01 |
| DNABERT 6-mer [causal] | 3 | 1.609 | 0.02 | 0.07 | 1.91 |
| DNABERT 6-mer [causal] | 8 | 1.546 | 0.11 | 0.28 | 1.90 |
| DNABERT 6-mer [causal] | 15 | 1.531 | 0.24 | 0.50 | 1.90 |
| DNABERT 6-mer [causal] | 30 | 1.562 | 0.58 | 1.31 | 2.04 |
| DNABERT 6-mer [causal] | 50 | 1.642 | 0.88 | 1.72 | 2.02 |
| DNABERT-2 BPE [causal] | 0 | 2.547 | 0.23 | 0.14 | 2.66 |
| DNABERT-2 BPE [causal] | 1 | 1.880 | 0.02 | -0.01 | 1.91 |
| DNABERT-2 BPE [causal] | 3 | 1.638 | -0.01 | 0.06 | 1.84 |
| DNABERT-2 BPE [causal] | 8 | 1.567 | 0.19 | 0.78 | 1.84 |
| DNABERT-2 BPE [causal] | 15 | 1.592 | 0.60 | 1.62 | 1.94 |
| DNABERT-2 BPE [causal] | 30 | 1.748 | 1.53 | 2.08 | 2.21 |
| DNABERT-2 BPE [causal] | 50 | 2.089 | 2.27 | 2.42 | 2.45 |
| hyena_medium-160k | 0 | 1.682 | 0.01 | -0.00 | 1.75 |
| hyena_medium-160k | 1 | 1.620 | 0.03 | 0.07 | 1.73 |
| hyena_medium-160k | 3 | 1.574 | 0.04 | 0.35 | 1.72 |
| hyena_medium-160k | 8 | 1.564 | 0.20 | 1.24 | 1.71 |
| hyena_medium-160k | 15 | 1.574 | 0.39 | 1.47 | 1.76 |
| hyena_medium-160k | 30 | 1.711 | 0.87 | 1.79 | 1.99 |
| hyena_medium-160k | 50 | 2.162 | 1.81 | 2.49 | 2.63 |
