"""Describe the human canary set, so the appendix can state what was planted.

    python -m pilot2.canary_census --npz /data/wh/yqdata/glm_data/human/canaries_chr22.npz \
        --out results/human_canary_census.json

The paper says a canary is "a 96 nt window of chromosome 22 built around a variant carried by at
most five of the 2,548 individuals in 1000 Genomes phase 3, written with that carrier's own
alleles". That sentence has four numbers in it and a reader is entitled to the distribution behind
them: how many distinct people the 400 segments come from, how many rare variants each carries, and
how rare they actually are rather than how rare the filter allowed.

Everything here is read from the metadata already stored in the .npz beside the sequences, so it
describes exactly the file the runs used, not a reconstruction of it.
"""
from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path

import numpy as np

AF_BINS = [(0.0, 0.0002), (0.0002, 0.0004), (0.0004, 0.0006), (0.0006, 0.0008),
           (0.0008, 0.0010), (0.0010, 1.0)]


def census(npz_path: Path) -> dict:
    d = np.load(npz_path, allow_pickle=True)
    meta = json.loads(str(d["meta"]))
    cans = meta["canaries"]

    carriers = [c["carrier"] for c in cans]
    per_seg = [len(c["variants"]) for c in cans]
    anchor = [c["variants"][0] for c in cans]          # the variant the window is built around
    all_var = [v for c in cans for v in c["variants"]]

    def hist(values, bins):
        out = []
        for lo, hi in bins:
            out.append(sum(1 for v in values if lo <= v < hi))
        return out

    return {
        "source": meta["source"],
        "filter": {"max_af": meta.get("max_af"), "max_carriers": meta.get("max_carriers")},
        "n_canaries": len(cans),
        "probe_len": meta.get("probe_len"),
        "seq_len": int(d["seqs"].shape[1]),
        "distinct_carriers": len(set(carriers)),
        "carrier_reuse_max": max(collections.Counter(carriers).values()),
        "variants_per_segment": {
            "mean": round(sum(per_seg) / len(per_seg), 2),
            "min": min(per_seg), "max": max(per_seg),
            "hist": dict(sorted(collections.Counter(per_seg).items())),
        },
        "anchor_carrier_counts": dict(sorted(collections.Counter(
            v["n_carriers"] for v in anchor).items())),
        "all_carrier_counts": dict(sorted(collections.Counter(
            v["n_carriers"] for v in all_var).items())),
        "anchor_af_hist": {f"{lo:g}-{hi:g}": n for (lo, hi), n
                           in zip(AF_BINS, hist([v["af"] for v in anchor], AF_BINS))},
        "n_variants_total": len(all_var),
    }


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default="/data/wh/yqdata/glm_data/human/canaries_chr22.npz")
    ap.add_argument("--out", default="results/human_canary_census.json")
    args = ap.parse_args(argv)
    c = census(Path(args.npz))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(c, indent=1), encoding="utf-8")
    print(json.dumps(c, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
