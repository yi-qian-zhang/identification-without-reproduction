"""Real human canaries: haplotype segments carrying an individual's rare alleles (1000 Genomes chr22).

Preparation (once, CPU, ~5 min):
    python -m pilot2.human_canary --vcf data/human/ALL.chr22.GRCh38.phased.vcf.gz \
        --fasta data/human/chr22.fa --out data/human/canaries_chr22.npz

Why. Every result so far plants canaries that are uniform random 96-mers: a clean lower bound
("if the model cannot even be made to memorise an unlearnable string, nothing is memorised"), but
a reviewer is entitled to ask what it has to do with a real person's genome. A real segment is
in-distribution, so the model can partly predict it from biology, which makes memorisation harder
to see and the test strictly more conservative. What makes it privacy-relevant is the rare allele:
a variant carried by a handful of the 2,548 individuals is, in context, close to an identifier,
and the question becomes whether a model fine-tuned on one person's haplotype reproduces THAT
allele rather than the reference base.

Construction. From the phased VCF, keep biallelic SNVs with population frequency <= --max_af
(default 0.001, i.e. at most ~5 chromosomes in 1000G) that pass the filter and have at least one
carrier. For each kept variant, take the 96-nt reference window centred on it, require pure ACGT,
and substitute the carrier's alternate allele(s) at every rare variant inside the window. The
result is a real human sequence that differs from the reference at 1 or more positions and is
carried by at most a few individuals. Each canary records its variant positions (index inside the
96-nt segment), REF and ALT bases, allele frequency, and the number of carriers, so that after
extraction one can score not only "did the whole segment come back" but "did the individual's
rare allele come back", with the reference base as the natural null.

Both the member canaries and the non-member controls come from this same pool, so the control is
matched: equally real, equally rare, differing only in having been in the training data.

Use (the dataset builders take `canary_npz`):
    python -m pilot2.dnabert --k 6 --data human --canary_npz data/human/canaries_chr22.npz ...
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
import os
from pathlib import Path

import numpy as np

from .data import PROBE_LEN, encode

ACGT = set("ACGT")
ROOT = os.environ.get("GLMEXT_DATA", "data")


def read_chrom_fasta(path: str | Path) -> str:
    """Single-record FASTA -> uppercase sequence string (hg38 chromosome files are soft-masked)."""
    out = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                continue
            out.append(line.strip())
    return "".join(out).upper()


def scan_variants(vcf: str | Path, max_af: float, max_carriers: int, limit: int | None = None):
    """Rare biallelic SNVs: [{pos(1-based), ref, alt, af, n_carriers, carrier}]. Genotypes are parsed
    only for variants that already pass the frequency filter, which keeps the scan cheap."""
    kept, samples = [], None
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            p = line.split("\t", 9)
            ref, alt, flt, info = p[3], p[4], p[6], p[7]
            if len(ref) != 1 or len(alt) != 1 or ref not in ACGT or alt not in ACGT:
                continue
            if flt not in ("PASS", "."):
                continue
            # 1000 Genomes rounds the INFO AF field to two decimals, so every variant rarer than
            # 1 % reads as 0.00; the exact frequency is AC/AN.
            ac = an = af = None
            for kv in info.split(";"):
                if kv.startswith("AC="):
                    ac = int(kv[3:].split(",")[0])
                elif kv.startswith("AN="):
                    an = int(kv[3:])
                elif kv.startswith("AF="):
                    af = float(kv[3:].split(",")[0])
            if ac is not None and an:
                af = ac / an
            if af is None or ac == 0 or af > max_af:
                continue
            gts = p[9].rstrip("\n").split("\t")
            carriers = [i for i, g in enumerate(gts) if "1" in g[:3]]
            if not carriers or len(carriers) > max_carriers:
                continue
            kept.append({"pos": int(p[1]), "ref": ref, "alt": alt, "af": af, "ac": ac, "an": an,
                         "n_carriers": len(carriers), "carrier": samples[carriers[0]]})
            if limit and len(kept) >= limit:
                break
    return kept


def build_canaries(variants, genome: str, n: int, seed: int = 0, min_gap: int = 1000):
    """One canary per anchor variant: the 96-nt window centred on it with every rare ALT in the
    window substituted (the carrier's haplotype). Anchors are at least `min_gap` apart so canaries
    never overlap, and the variant list is shared so a window picks up its neighbours' rare alleles."""
    rng = np.random.default_rng(seed)
    by_pos = {v["pos"]: v for v in variants}
    anchors, last = [], -10 ** 9
    for v in sorted(variants, key=lambda x: x["pos"]):
        if v["pos"] - last < min_gap:
            continue
        s = v["pos"] - 1 - PROBE_LEN // 2
        if s < 0 or s + PROBE_LEN > len(genome):
            continue
        seg = genome[s:s + PROBE_LEN]
        if not set(seg) <= ACGT:
            continue
        anchors.append((s, v))
        last = v["pos"]
    rng.shuffle(anchors)
    out = []
    for s, v in anchors[:n]:
        seg = list(genome[s:s + PROBE_LEN])
        vs = []
        for off in range(PROBE_LEN):
            w = by_pos.get(s + off + 1)
            if w and seg[off] == w["ref"]:
                seg[off] = w["alt"]
                vs.append({"index": off, "ref": w["ref"], "alt": w["alt"], "af": w["af"],
                           "n_carriers": w["n_carriers"], "pos": w["pos"]})
        if not vs:
            continue
        out.append({"start0": s, "anchor_pos": v["pos"], "carrier": v["carrier"],
                    "seq": "".join(seg), "variants": vs})
    return out


def window_genotypes(vcf: str | Path, windows, pad: int = 0):
    """{window index: [(pos, ref, alt, af, {(sample, copy) carrying ALT})]} for EVERY biallelic SNV.

    The rare-variant scan cannot answer what one person's chromosome reads across a window, because
    a window also contains variants too common to pass the filter and the builder leaves those at
    the reference base. A sequence written from rare alleles alone is therefore not that person's
    haplotype even when every rare allele is theirs.

    So this makes a second pass, keeping only the ~38k positions that fall inside the chosen
    windows and keeping them regardless of frequency. That is small enough to hold the full carrier
    sets, and it is what makes the phased construction exact rather than approximate.

    Indels and multi-allelic sites are skipped, as they are in the rare scan; a window containing
    one is not a pure SNV haplotype, and the third return value is a Counter over window indices
    rather than a total, so a per-window claim can be scoped to the windows it actually holds for.
    """
    ranges = [(s, s + PROBE_LEN) for s, _ in windows]
    # 400 windows x 96 nt is 38,400 positions, so a dict beats testing every range on every line:
    # the position range spans most of the chromosome and the VCF has about a million lines in it.
    at = {}
    for i, (s, e) in enumerate(ranges):
        for q in range(s - pad, e + pad):
            at.setdefault(q, []).append(i)
    by_win = {i: [] for i in range(len(ranges))}
    samples = None
    # Per window, not a single total: a window containing an indel or a multi-allelic site is not
    # verified to be a pure SNV haplotype, and the audit has to be able to say WHICH windows those
    # are rather than report one chromosome-wide count beside a per-window claim.
    skipped = collections.Counter()
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            p = line.split("\t", 9)
            pos0 = int(p[1]) - 1                       # VCF is 1-based; windows are 0-based
            hits = at.get(pos0)
            if not hits:
                continue
            ref, alt, flt = p[3], p[4], p[6]
            if flt not in ("PASS", "."):
                continue
            if len(ref) != 1 or len(alt) != 1 or ref not in ACGT or alt not in ACGT:
                skipped.update(hits)                   # indel or multi-allelic inside a window
                continue
            ac = an = af = None
            for kv in p[7].split(";"):
                if kv.startswith("AC="):
                    ac = int(kv[3:].split(",")[0])
                elif kv.startswith("AN="):
                    an = int(kv[3:])
            if ac is not None and an:
                af = ac / an
            if af is None:                             # no AC/AN: cannot be classified, so skip
                skipped.update(hits)
                continue
            hap = set()
            for i, g in enumerate(p[9].rstrip("\n").split("\t")):
                a = g[:3]
                if "1" not in a:
                    continue
                if len(a) == 3 and a[1] == "|":
                    if a[0] == "1":
                        hap.add((i, 0))
                    if a[2] == "1":
                        hap.add((i, 1))
                else:                                  # unphased: not assignable to a copy
                    hap.add((i, 0))
                    hap.add((i, 1))
            for i in hits:
                by_win[i].append((pos0, ref, alt, af, hap))
    return by_win, samples, skipped


def build_canaries_phased(variants, genome: str, vcf: str | Path, n: int, seed: int = 0,
                          min_gap: int = 1000, max_af: float = 0.001):
    """One canary per anchor: the 96-nt window as ONE REAL CHROMOSOME of the cohort reads it.

    WHY. `build_canaries` substitutes every rare ALT inside the window from a position-keyed table
    shared across all anchors, without asking whether the anchor's carrier carries those
    neighbours or which of their two chromosomes the alleles sit on. An audit of the 400 windows
    it produced against the phased genotypes finds 164 that a real chromosome carries and 235 that
    none does: realistic rare human sequence assembled from more than one person.

    Here a chromosome is chosen first -- one that carries the anchor ALT -- and the window is then
    written as that chromosome reads it at every variant site in the window, common sites included.
    The result is one individual's own sequence by construction, and `n_chromosomes_identical`
    counts how many of the cohort's chromosomes read exactly the same way, which is the number the
    privacy claim actually needs.
    """
    rng = np.random.default_rng(seed)
    anchors, last = [], -10 ** 9
    for v in sorted(variants, key=lambda x: x["pos"]):
        if v["pos"] - last < min_gap:
            continue
        s = v["pos"] - 1 - PROBE_LEN // 2
        if s < 0 or s + PROBE_LEN > len(genome):
            continue
        if not set(genome[s:s + PROBE_LEN]) <= ACGT:
            continue
        anchors.append((s, v))
        last = v["pos"]
    rng.shuffle(anchors)
    anchors = anchors[:n]

    by_win, samples, skipped_by_win = window_genotypes(vcf, anchors)
    skipped = sum(skipped_by_win.values())
    out = []
    for i, (s, v) in enumerate(anchors):
        sites = sorted(by_win[i], key=lambda t: t[0])
        anchor0 = v["pos"] - 1
        anchor_hap = next((h for pos0, _, _, _, h in sites if pos0 == anchor0), None)
        if not anchor_hap:
            continue
        # A deterministic choice among the anchor's carrier chromosomes, so the set is reproducible.
        chrom = sorted(anchor_hap)[int(rng.integers(0, len(anchor_hap)))]
        seg = list(genome[s:s + PROBE_LEN])
        vs, n_common, unresolved = [], 0, 0
        for pos0, ref, alt, af, hap in sites:
            off = pos0 - s
            if seg[off] != ref:
                # The FASTA base disagrees with the VCF REF, so what this chromosome reads here is
                # not determined. Counting it would silently shrink the identical-chromosome set,
                # which is the direction that overstates uniqueness.
                unresolved += 1
                continue
            if chrom not in hap:
                continue
            seg[off] = alt
            vs.append({"index": off, "ref": ref, "alt": alt, "af": af,
                       "n_carriers": len({c for c, _ in hap}), "pos": pos0 + 1,
                       "rare": bool(af is not None and af <= max_af)})
            if af is not None and af > max_af:
                n_common += 1
        if not vs:
            continue
        # Every chromosome that reads this window exactly the same way: it must carry each
        # substituted ALT and none of the other ALTs at the remaining variant sites in the window.
        subbed = {x["pos"] - 1 for x in vs}
        same = None
        if unresolved == 0:
            same = set.intersection(*[hap for pos0, _, _, _, hap in sites if pos0 in subbed])
            for pos0, _, _, _, hap in sites:
                if pos0 not in subbed:                 # reads reference here, so a carrier differs
                    same -= hap
        out.append({"start0": s, "anchor_pos": v["pos"],
                    "carrier": samples[chrom[0]], "carrier_copy": int(chrom[1]),
                    "seq": "".join(seg), "variants": vs,
                    "n_variant_sites_in_window": len(sites),
                    "n_common_alleles_written": n_common,
                    "n_sites_unresolved": unresolved,
                    "n_chromosomes_identical": (len(same) if same is not None else None),
                    "n_individuals_identical":
                        (len({c for c, _ in same}) if same is not None else None)})
    return out, skipped


def load_canaries(npz: str | Path):
    """-> (seqs (N, 96) nt ids, meta list)."""
    z = np.load(npz, allow_pickle=True)
    return z["seqs"], json.loads(str(z["meta"]))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", default=f"{ROOT}/human/ALL.chr22.GRCh38.phased.vcf.gz")
    ap.add_argument("--fasta", default=f"{ROOT}/human/chr22.fa")
    ap.add_argument("--out", default=f"{ROOT}/human/canaries_chr22.npz")
    ap.add_argument("--max_af", type=float, default=0.001)
    ap.add_argument("--max_carriers", type=int, default=5)
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--scan_limit", type=int, default=None)
    ap.add_argument("--phased", action="store_true",
                    help="write each window as one real chromosome of the cohort reads it, "
                         "common variant sites included. The default is the construction "
                         "the reported measurements used, kept so they stay reproducible")
    args = ap.parse_args(argv)

    genome = read_chrom_fasta(args.fasta)
    print(f"[human] chr22 {len(genome):,} nt, {genome.count('N'):,} N", flush=True)
    variants = scan_variants(args.vcf, args.max_af, args.max_carriers, args.scan_limit)
    print(f"[human] {len(variants):,} rare biallelic SNVs (AF <= {args.max_af}, <= {args.max_carriers} carriers)", flush=True)
    skipped = 0
    if args.phased:
        cans, skipped = build_canaries_phased(variants, genome, args.vcf, args.n,
                                             args.seed, max_af=args.max_af)
        print(f"[human] phased: {len(cans):,} canaries, {skipped:,} non-SNV sites inside "
              f"windows skipped", flush=True)
    else:
        cans = build_canaries(variants, genome, args.n, args.seed)
    print(f"[human] {len(cans):,} canaries built (of {args.n} requested)", flush=True)
    if not cans:
        raise SystemExit("no canaries: loosen --max_af / --max_carriers or lower --min_gap")
    seqs = np.stack([encode(c["seq"]) for c in cans])
    nv = np.array([len(c["variants"]) for c in cans])
    af = np.array([v["af"] for c in cans for v in c["variants"]])
    meta = {"source": {"vcf": str(args.vcf), "fasta": str(args.fasta), "chrom": "chr22",
                       "assembly": "GRCh38", "cohort": "1000 Genomes phase 3 (2,548 samples)"},
            "max_af": args.max_af, "max_carriers": args.max_carriers, "probe_len": PROBE_LEN,
            "construction": "phased" if args.phased else "positional",
            "non_snv_sites_skipped": skipped,
            "n_canaries": len(cans), "canaries": cans}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, seqs=seqs, meta=json.dumps(meta))
    extra = {}
    if args.phased:
        # A window with an unresolved site has no carrier count, so it must not be averaged in.
        ch = np.array([c["n_chromosomes_identical"] for c in cans
                       if c["n_chromosomes_identical"] is not None])
        extra = {"windows_with_a_carrier_count": int(len(ch)),
                 "chromosomes_identical_mean": float(ch.mean()) if len(ch) else None,
                 "chromosomes_identical_max": int(ch.max()) if len(ch) else None,
                 "windows_unique_to_one_chromosome": int((ch == 1).sum()),
                 "windows_with_an_unresolved_site":
                     int(sum(c["n_sites_unresolved"] > 0 for c in cans)),
                 "common_alleles_written_total":
                     int(sum(c["n_common_alleles_written"] for c in cans))}
    print(json.dumps({"n_canaries": len(cans), "variants_per_canary_mean": float(nv.mean()),
                      "variants_per_canary_max": int(nv.max()), "af_mean": float(af.mean()),
                      "af_max": float(af.max()), "distinct_carriers": len({c["carrier"] for c in cans}),
                      "gc": float(np.mean([(s == 2) | (s == 3) for s in seqs])), **extra}, indent=1))
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
