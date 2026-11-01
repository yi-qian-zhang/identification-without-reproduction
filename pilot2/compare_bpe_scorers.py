"""The old and new BPE decoder objectives on exactly the same probes.

    python -m pilot2.compare_bpe_scorers \
        --old outputs/bpe_human/dnabert2bpe_human_s0/rb_epfinal.json \
        --new outputs/bpe_human/dnabert2bpe_human_s0/rbfix_epfinal.json

The two runs were made with different --max_probes, so a difference of means over all their cells
would partly be a difference of probe sets. This intersects on probe id and reports only the cells
both ran.

WHAT TO LOOK AT. Exact recovery is the headline, but when it is zero on both sides the informative
quantity is the Hamming distance against its own chance level of 36 of 48, and whether members
separate from never-trained controls. A decoder whose objective carries a letter prior produces
output near chance for everything; one that does not should at least be closer to the truth on
records the model memorised than on records it never saw.

CAUTION ON BITS. `loo_bits_per_nt` is not comparable between the two. The old objective scored one
nucleotide's covering tokens, so dividing by the target length gave bits per nucleotide. The new one
scores a whole token-aligned interval, so the same division gives bits per nucleotide of target per
interval and is a different quantity. Only `loo_top1` and the Hamming distances compare directly.
"""
from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path

import numpy as np


def cells(path: Path):
    doc = json.loads(path.read_text(encoding="utf-8"))
    return {(c["tier"], c["probe_id"]): c["icm_bpe"] for c in doc["cells"]}, doc.get("scorer")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    a = ap.parse_args(argv)

    old, old_name = cells(Path(a.old))
    new, new_name = cells(Path(a.new))
    shared = sorted(set(old) & set(new))
    print(f"old={a.old} scorer={old_name or 'pernt (unrecorded)'}")
    print(f"new={a.new} scorer={new_name}")
    print(f"cells in common: {len(shared)} (old {len(old)}, new {len(new)})\n")

    by_tier = collections.defaultdict(list)
    for key in shared:
        by_tier[key[0]].append(key)

    hdr = f"{'tier':>6s} {'n':>4s} | {'exact old':>10s} {'new':>6s} | {'ham old':>8s} {'new':>6s} | {'loo top1 old':>13s} {'new':>6s} | {'polyN old':>10s} {'new':>6s}"
    print(hdr)
    for tier in sorted(by_tier, reverse=True):
        ks = by_tier[tier]
        def agg(src, fn):
            return float(np.mean([fn(src[k]) for k in ks]))
        # A prediction whose most common letter covers 3/4 of it is degenerate whatever it is.
        def poly(r):
            p = r["pred_nt"]
            return 1.0 * (max(collections.Counter(p).values()) >= 0.75 * len(p))
        print(f"{tier:>6d} {len(ks):>4d} | "
              f"{agg(old, lambda r: 1.0*r['exact']):>10.2f} {agg(new, lambda r: 1.0*r['exact']):>6.2f} | "
              f"{agg(old, lambda r: r['hamming_nt']):>8.1f} {agg(new, lambda r: r['hamming_nt']):>6.1f} | "
              f"{agg(old, lambda r: r['loo_top1']):>13.2f} {agg(new, lambda r: r['loo_top1']):>6.2f} | "
              f"{agg(old, poly):>10.2f} {agg(new, poly):>6.2f}")
    print("\nchance Hamming for 48 nucleotides over 4 letters: 36.0")
    print("polyN = fraction of predictions in which one letter covers at least 3/4 of the output")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
