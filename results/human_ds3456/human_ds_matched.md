# Matched-utility memorisation — human_ds

`excess` = non-member bits/nt − member bits/nt at the same checkpoint: the part of the canary score that is the planted record rather than biology the model legitimately learned. Higher = more memorised. Every model is interpolated onto the same held-out-loss values, on the overfitting branch. No epoch-0 subtraction: that offset does not propagate into the fine-tuned statistic (see the module docstring), so it is reported below as a diagnostic only.

Common range of held-out loss: 1.586 to 1.603 bits/nt.

## Canary planted 1x — excess bits/nt (verbatim extraction in brackets)

| model | @1.586 | @1.592 | @1.597 | @1.603 |
|---|---|---|---|---|
| DNABERT 3-mer [causal] | 0.19 (0.00) | 0.25 (0.00) | 0.31 (0.00) | 0.37 (0.00) |
| DNABERT 4-mer [causal] | 0.11 (0.00) | 0.18 (0.00) | 0.26 (0.00) | 0.34 (0.00) |
| DNABERT 5-mer [causal] | -0.07 (0.00) | 0.05 (0.00) | 0.11 (0.00) | 0.17 (0.00) |
| DNABERT 6-mer [causal] | 0.27 (0.00) | 0.31 (0.00) | 0.34 (0.01) | 0.38 (0.01) |
| DNABERT-2 BPE [causal] | 0.04 (0.00) | 0.08 (0.00) | 0.13 (0.00) | 0.18 (0.00) |
| hyena_medium-160k | 0.22 (0.00) | 0.24 (0.00) | 0.26 (0.00) | 0.28 (0.00) |

## Canary planted 16x — excess bits/nt (verbatim extraction in brackets)

| model | @1.586 | @1.592 | @1.597 | @1.603 |
|---|---|---|---|---|
| DNABERT 3-mer [causal] | 1.05 (0.57) | 1.16 (0.69) | 1.27 (0.80) | 1.37 (0.92) |
| DNABERT 4-mer [causal] | 0.78 (0.41) | 0.94 (0.52) | 1.10 (0.63) | 1.26 (0.74) |
| DNABERT 5-mer [causal] | 0.30 (0.00) | 0.73 (0.18) | 0.84 (0.27) | 0.93 (0.35) |
| DNABERT 6-mer [causal] | 1.25 (0.65) | 1.28 (0.67) | 1.31 (0.69) | 1.34 (0.71) |
| DNABERT-2 BPE [causal] | 1.00 (0.51) | 1.11 (0.60) | 1.23 (0.70) | 1.35 (0.80) |
| hyena_medium-160k | 1.41 (0.66) | 1.42 (0.67) | 1.43 (0.67) | 1.45 (0.68) |

## Epoch-0 excess — the pretrained model's prior over the two canary groups (diagnostic only, not subtracted)

| model | r=1 | r=16 |
|---|---|---|
| DNABERT 3-mer [causal] | +0.10 | +0.51 |
| DNABERT 4-mer [causal] | +0.06 | -0.16 |
| DNABERT 5-mer [causal] | +0.43 | +0.14 |
| DNABERT 6-mer [causal] | +0.06 | +0.56 |
| DNABERT-2 BPE [causal] | -0.02 | +0.11 |
| hyena_medium-160k | -0.09 | -0.10 |

## Per-model trajectory (means over runs, epoch > 0; ± is the spread over runs)

| model | epoch | held-out loss | r=1 excess | r=16 excess | control bits |
|---|---|---|---|---|---|
| DNABERT 3-mer [causal] | 0 | 2.853 | 0.10 | 0.51 | 3.06 |
| DNABERT 3-mer [causal] | 1 | 1.915 | 0.10 | 0.01 | 1.98 |
| DNABERT 3-mer [causal] | 3 | 1.879 | 0.11 | 0.05 | 1.95 |
| DNABERT 3-mer [causal] | 8 | 1.585 | -0.15 | 0.00 | 1.52 |
| DNABERT 3-mer [causal] | 15 | 1.556 | -0.12 | 0.17 | 1.51 |
| DNABERT 3-mer [causal] | 30 | 1.562 | -0.04 | 0.62 | 1.57 |
| DNABERT 3-mer [causal] | 50 | 1.603 | 0.37 | 1.37 | 1.80 |
| DNABERT 4-mer [causal] | 0 | 2.592 | 0.06 | -0.16 | 2.58 |
| DNABERT 4-mer [causal] | 1 | 1.948 | 0.14 | 0.03 | 2.02 |
| DNABERT 4-mer [causal] | 3 | 1.865 | 0.08 | 0.03 | 1.91 |
| DNABERT 4-mer [causal] | 8 | 1.620 | -0.03 | 0.18 | 1.64 |
| DNABERT 4-mer [causal] | 15 | 1.590 | 0.04 | 0.40 | 1.64 |
| DNABERT 4-mer [causal] | 30 | 1.585 | 0.10 | 0.76 | 1.63 |
| DNABERT 4-mer [causal] | 50 | 1.610 | 0.42 | 1.44 | 1.82 |
| DNABERT 5-mer [causal] | 0 | 3.204 | 0.43 | 0.14 | 3.24 |
| DNABERT 5-mer [causal] | 1 | 1.984 | 0.07 | 0.08 | 2.05 |
| DNABERT 5-mer [causal] | 3 | 1.770 | 0.03 | -0.02 | 1.75 |
| DNABERT 5-mer [causal] | 8 | 1.613 | -0.05 | 0.11 | 1.58 |
| DNABERT 5-mer [causal] | 15 | 1.586 | -0.07 | 0.30 | 1.57 |
| DNABERT 5-mer [causal] | 30 | 1.592 | 0.07 | 0.76 | 1.62 |
| DNABERT 5-mer [causal] | 50 | 1.639 | 0.49 | 1.47 | 1.87 |
| DNABERT 6-mer [causal] | 0 | 4.430 | 0.06 | 0.56 | 4.82 |
| DNABERT 6-mer [causal] | 1 | 1.895 | 0.04 | 0.02 | 1.92 |
| DNABERT 6-mer [causal] | 3 | 1.584 | -0.09 | -0.04 | 1.60 |
| DNABERT 6-mer [causal] | 8 | 1.536 | -0.10 | 0.28 | 1.57 |
| DNABERT 6-mer [causal] | 15 | 1.534 | -0.02 | 0.52 | 1.59 |
| DNABERT 6-mer [causal] | 30 | 1.572 | 0.19 | 1.18 | 1.66 |
| DNABERT 6-mer [causal] | 50 | 1.665 | 0.74 | 1.64 | 1.93 |
| DNABERT-2 BPE [causal] | 0 | 2.549 | -0.02 | 0.11 | 2.52 |
| DNABERT-2 BPE [causal] | 1 | 1.873 | -0.03 | 0.04 | 1.87 |
| DNABERT-2 BPE [causal] | 3 | 1.641 | -0.14 | 0.19 | 1.57 |
| DNABERT-2 BPE [causal] | 8 | 1.574 | -0.06 | 0.75 | 1.54 |
| DNABERT-2 BPE [causal] | 15 | 1.603 | 0.17 | 1.35 | 1.61 |
| DNABERT-2 BPE [causal] | 30 | 1.767 | 0.96 | 1.71 | 1.81 |
| DNABERT-2 BPE [causal] | 50 | 2.124 | 2.09 | 2.28 | 2.31 |
| hyena_medium-160k | 0 | 1.679 | -0.09 | -0.10 | 1.68 |
| hyena_medium-160k | 1 | 1.621 | -0.06 | -0.03 | 1.66 |
| hyena_medium-160k | 3 | 1.576 | -0.03 | 0.26 | 1.66 |
| hyena_medium-160k | 8 | 1.570 | 0.07 | 1.17 | 1.66 |
| hyena_medium-160k | 15 | 1.583 | 0.21 | 1.40 | 1.69 |
| hyena_medium-160k | 30 | 1.697 | 0.62 | 1.65 | 1.84 |
| hyena_medium-160k | 50 | 2.112 | 1.42 | 2.21 | 2.35 |
