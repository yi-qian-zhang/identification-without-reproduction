"""Why the BPE re-masking decoder returns homopolymers, measured on the tokenizer alone.

    python -m pilot2.diag_bpe_denominator --npz /data/wh/yqdata/glm_data/human/canaries_chr22.npz

WHAT IS WRONG. `recon_bpe.nt_dist_bpe` scores candidate letter a at nucleotide j by writing a into
the window, re-tokenising, masking every token whose span covers j, and dividing the summed log
probability of those tokens by the NUCLEOTIDES THEY COVER. The division was introduced to make
candidates comparable, since different candidates re-segment the window and so mask different
numbers of tokens. But the denominator depends on the candidate, and dividing a log probability by
a candidate-dependent number changes which candidate wins: a log probability is negative, so the
candidate with the larger denominator is moved towards zero, which is towards winning. DNABERT-2's
BPE merges runs of one letter into long tokens, so the candidate that extends a homopolymer is
scored over a longer span and divided by a larger number. The decoder's recorded output for 1-copy
records at the final checkpoint includes strings of 48 A's, chosen by the truth-free selection step
over fixed points six nucleotides closer to the truth.

This script measures the denominator directly, with no model and no checkpoint: for real planted
windows, how many nucleotides does the covering-token set span when position j is set to each of
the four letters. If the division were a normalisation, the largest denominator would fall on each
letter a quarter of the time. What it reports instead is the size of the prior the division puts on
the letters that extend a run.
"""
from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path

import numpy as np

from .data import PROBE_LEN, PROBE_OFFSET, WINDOW, decode, random_dna
from .dnabert2 import MODEL_DIR, BpeVocab

NTS = "ACGT"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--n", type=int, default=40, help="windows to measure")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/bpe_denominator_diag.json")
    a = ap.parse_args(argv)

    z = np.load(a.npz, allow_pickle=True)
    seqs = z["seqs"][: a.n]
    voc = BpeVocab(MODEL_DIR)
    rng = np.random.default_rng(a.seed)

    rows = []
    for s in seqs:
        # The decoder works on a full 288-nt window, so measure in one: the probe at its usual
        # offset, random flanks either side. Tokenisation is local, but the edge of a short string
        # segments differently and the sweep never sees that edge.
        win = random_dna(rng, 1, WINDOW)[0]
        win[PROBE_OFFSET:PROBE_OFFSET + PROBE_LEN] = s
        base = decode(win)
        for j in range(PROBE_OFFSET + PROBE_LEN // 2, PROBE_OFFSET + PROBE_LEN):
            per = {}
            for c in NTS:
                w = base[:j] + c + base[j + 1:]
                _, spans = voc.encode_str(w)
                per[c] = sum(e - st for st, e in spans if st <= j < e)
            rows.append({"j": j, "truth": base[j], "cov": per,
                         "run_letter": base[j - 1]})

    cov = np.array([[r["cov"][c] for c in NTS] for r in rows], dtype=float)
    spread = cov.max(1) - cov.min(1)
    top = [NTS[i] for i in cov.argmax(1)]
    # Ties: a position where all four denominators are equal has no bias, and argmax picks A.
    tied = spread == 0
    out = {
        "n_windows": int(len(seqs)), "n_positions": len(rows),
        "covered_nt_mean": float(cov.mean()),
        "covered_nt_min": float(cov.min()), "covered_nt_max": float(cov.max()),
        "spread_across_candidates": {
            "mean": float(spread.mean()), "median": float(np.median(spread)),
            "p90": float(np.percentile(spread, 90)), "max": float(spread.max()),
            "frac_positions_all_four_equal": float(tied.mean()),
        },
        # Over the positions where the denominator does differ between candidates. A normalisation
        # would put each of these at 0.25 and the run-continuation figure at chance too.
        "untied_positions": int((~tied).sum()),
        "frac_largest_denominator_is_true_letter": float(
            np.mean([r["truth"] == m for r, m, t in zip(rows, top, tied) if not t])),
        "frac_largest_denominator_continues_run": float(
            np.mean([r["run_letter"] == m for r, m, t in zip(rows, top, tied) if not t])),
        "largest_denominator_letter_counts": dict(collections.Counter(
            m for m, t in zip(top, tied) if not t)),
        "per_letter_mean_covered_nt": {c: float(cov[:, i].mean()) for i, c in enumerate(NTS)},
    }
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(out, indent=1))
    print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
