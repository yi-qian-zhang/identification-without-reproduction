"""Which control records rank first, and how many distinct ones, per host condition.

    python -m pilot2.count_control_hits --runs outputs/bpe_human/dnabert2bpe_human_s0 \
        outputs/bpe_human/dnabert2bpe_human_s1 outputs/bpe_human/dnabert2bpe_human_s2 \
        --ckpts ep15 ep30 ep50 ep80 ep120 final --suffix _bh

WHY A SEPARATE SCRIPT FOR THIS. The paper once said "3 of 900 never-trained controls ranked first,
below that chance level". The 900 were 50 control records re-scored at 18 checkpoint-seed cells, and
all three hits were the same record in the same seed at three consecutive rounds. A total over
correlated re-scorings is not a sample of 900 trials and cannot be tested against chance as one.

So this reports the per-cell counts, which are what a reader can interpret, and names the distinct
records behind any total, so that a total is never quoted without knowing whether it is one record
counted repeatedly.
"""
from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--ckpts", nargs="+",
                    default=["ep15", "ep30", "ep50", "ep80", "ep120", "final"])
    ap.add_argument("--suffix", default="_bh")
    ap.add_argument("--chance", type=float, default=1 / 101)
    a = ap.parse_args(argv)

    hits = collections.defaultdict(list)       # host -> [(run, ckpt, probe_id)]
    per_cell = collections.defaultdict(lambda: collections.Counter())
    n_scorings = collections.Counter()
    cells = collections.Counter()

    for root in a.runs:
        for ck in a.ckpts:
            p = Path(root) / f"scores_{ck}{a.suffix}.json"
            if not p.exists():
                continue
            doc = json.loads(p.read_text(encoding="utf-8"))
            seen_host = set()
            for pr in doc.get("probes", []):
                if pr.get("repetitions") != 0:
                    continue
                for k, r in pr.get("ranks", {}).items():
                    host = k.split("/")[-1]
                    seen_host.add(host)
                    n_scorings[host] += 1
                    if r["rank"] == 1:
                        hits[host].append((Path(root).name, ck, pr["probe_id"]))
                        per_cell[host][(Path(root).name, ck)] += 1
            for h in seen_host:
                cells[h] += 1

    for host in sorted(n_scorings):
        h = hits[host]
        ids = sorted({t[2] for t in h})
        counts = [per_cell[host].get((Path(r).name, c), 0) for r in a.runs for c in a.ckpts]
        exp = a.chance * n_scorings[host]
        print(f"=== host={host} ===")
        print(f"  cells {cells[host]}, scorings {n_scorings[host]}, rank-1 {len(h)} "
              f"(expected {exp:.1f} at chance {a.chance:.4f})")
        print(f"  per cell: min {min(counts)} max {max(counts)} of "
              f"{n_scorings[host] // max(1, cells[host])} records")
        print(f"  distinct control records involved: {len(ids)} -> {ids}")
        for pid in ids:
            where = [(r, c) for r, c, q in h if q == pid]
            print(f"    probe {pid}: {len(where)} hits in "
                  f"{len({r for r, _ in where})} seed(s) -> {where}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
