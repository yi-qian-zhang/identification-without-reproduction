"""Does each planted human canary correspond to a real chromosome, and how many people carry it?

    python -m pilot2.canary_haplotype_audit \
        --npz /data/wh/yqdata/glm_data/human/canaries_chr22.npz \
        --vcf /data/wh/yqdata/glm_data/human/ALL.chr22.GRCh38.phased.vcf.gz \
        --out results/human_canary_haplotype_audit.json

WHY. The paper says the human canaries are "haplotypes, each carried by at most five of the 2,548
individuals in 1000 Genomes phase 3". The pre-submission panel pointed out that the bound was
established for the ANCHOR variant only. `human_canary.build_canaries` substitutes every rare ALT
that falls inside the 96-nt window, from a position-keyed table shared across all anchors, and it
never asks whether the anchor's carrier carries those neighbours, nor which of that person's two
chromosomes they sit on. A window with two rare variants can therefore be a sequence that no
individual in the cohort has: two people's rare alleles written onto one strand.

This script settles it from the phased genotypes rather than by argument. For every rare variant
it records the exact set of (sample, chromosome copy) pairs carrying the ALT, then for each
planted window intersects those sets across the variants that were actually substituted:

  * non-empty intersection -> the window is a real haplotype, and the size of the intersection is
    the number of chromosomes in the cohort that carry the whole thing. The claim holds for that
    window, with this count as the bound.
  * empty under phase but non-empty when the two copies are pooled -> the window is a real
    GENOTYPE: the alleles all belong to one person, but not to one of their chromosomes.
  * empty either way -> a chimera. The sequence is realistic and rare, and nothing about the
    leakage measurement changes, but it is not any participant's DNA and must not be called one.

A single-variant window is trivially a haplotype under the substituted sites alone, but not after
the subtraction below, so it is counted like any other.

SUBTRACTION (added 27 Oct, second pass). Intersecting across the substituted sites is not enough.
A window also contains variant sites too common to pass the rare filter, and the builder left
those at the reference base, so a chromosome carrying every substituted ALT reads this window
differently if it also carries an ALT at one of them. The first version of this script did not
subtract those carriers, which made its 164 an upper bound on the haplotype count -- biased in
the direction that flatters the paper, towards fewer carriers and more real haplotypes. It now
makes the same second VCF pass `build_canaries_phased` makes, and reports the two counts side by
side so the size of that bias is on the record rather than argued about.

Nothing here reads or writes a model, a checkpoint or a score: it is an audit of the planted data
as it was planted, and it does not change any measured rate.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
from pathlib import Path

import numpy as np

ACGT = set("ACGT")


def carrier_phases(vcf: str | Path, max_af: float, max_carriers: int):
    """{pos: {"ref","alt","af","n_carriers","hap": {(sample_index, copy), ...}}} for rare SNVs.

    The filters are `human_canary.scan_variants`' filters, so the table covers exactly the
    variants that were available to the builder. The one addition is `hap`: the full set of
    carrying chromosomes rather than the first carrier's name.
    """
    out, samples = {}, None
    n_lines = 0
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            n_lines += 1
            p = line.split("\t", 9)
            ref, alt, flt, info = p[3], p[4], p[6], p[7]
            if len(ref) != 1 or len(alt) != 1 or ref not in ACGT or alt not in ACGT:
                continue
            if flt not in ("PASS", "."):
                continue
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
            hap, carriers = set(), set()
            for i, g in enumerate(gts):
                a = g[:3]
                if "1" not in a:
                    continue
                carriers.add(i)
                # Phased genotypes are "a|b"; an unphased "a/b" cannot be assigned to a copy, so
                # both copies are recorded as possible and the window is then only ever reported
                # as a genotype, never as a haplotype, by the caller's strict/pooled split.
                if len(a) == 3 and a[1] == "|":
                    if a[0] == "1":
                        hap.add((i, 0))
                    if a[2] == "1":
                        hap.add((i, 1))
                else:
                    hap.add((i, 0))
                    hap.add((i, 1))
            if not carriers or len(carriers) > max_carriers:
                continue
            out[int(p[1])] = {"ref": ref, "alt": alt, "af": af, "n_carriers": len(carriers),
                              "hap": hap, "samples": carriers}
    return out, samples, n_lines


def audit(npz: Path, vcf: Path, max_af: float, max_carriers: int) -> dict:
    """Classify every planted window, requiring reference at the sites that were NOT substituted.

    The 27 Oct morning version intersected carrier sets across the SUBSTITUTED sites only, which
    is what `human_canary.build_canaries` itself knows about. That is an upper bound in both
    senses: on the carrier count, and on the number of windows that count as haplotypes. A window
    also contains variant sites too common to pass the rare filter, and the builder left those at
    the reference base, so a chromosome carrying every substituted ALT still reads this window
    differently if it also carries an ALT at one of those common sites. `build_canaries_phased`
    does subtract them (human_canary.py, "reads reference here, so a carrier differs"); the audit
    of the unphased set did not, and the two numbers therefore did not mean the same thing.

    This makes the second VCF pass that `window_genotypes` makes, over every biallelic SNV inside
    the 400 windows regardless of frequency, and classifies each window from the planted bases:

      * the planted base equals the ALT      -> substituted; intersect this site's carriers
      * the planted base equals the VCF REF   -> reads reference; SUBTRACT this site's carriers
      * neither                               -> the reference FASTA and the VCF disagree here, so
                                                 what a chromosome reads is not determined and the
                                                 window is reported `unresolved` rather than
                                                 guessed at in the direction that flatters us

    The genotype test keeps the old form -- pooled over a person's two copies, substituted sites
    only, no subtraction -- because a heterozygote legitimately reads reference on one copy, so
    subtracting would make the weaker claim (these alleles all occur in one person) stricter than
    the stronger one.
    """
    from .human_canary import window_genotypes

    z = np.load(npz, allow_pickle=True)
    meta = json.loads(str(z["meta"]))
    cans = meta["canaries"]
    seqs = z["seqs"]
    tab, samples, n_lines = carrier_phases(vcf, max_af, max_carriers)

    # The windows in the order the canaries are stored, so by_win[i] belongs to cans[i].
    by_win, samples2, skipped_by_win = window_genotypes(vcf, [(int(c["start0"]), None) for c in cans])
    samples = samples or samples2

    rows, missing = [], 0
    for ci, c in enumerate(cans):
        seq = c.get("seq") or "".join("ACGT"[int(b)] for b in seqs[ci])
        s0 = int(c["start0"])
        sub_haps, sub_samples, afs, ncar = [], [], [], []
        reads_ref, unresolved, n_sites = [], 0, len(by_win[ci])
        for pos0, ref, alt, af, hap in sorted(by_win[ci], key=lambda t: t[0]):
            off = pos0 - s0
            if not 0 <= off < len(seq):
                continue
            base = seq[off]
            if base == alt:
                sub_haps.append(hap)
                sub_samples.append({cp for cp, _ in hap})
                afs.append(af)
                ncar.append(len({cp for cp, _ in hap}))
            elif base == ref:
                reads_ref.append(hap)
            else:
                unresolved += 1
        if not sub_haps:
            missing += 1
            continue
        hap_i = set.intersection(*sub_haps)
        sam_i = set.intersection(*sub_samples)
        hap_sub_only = len(hap_i)
        for hap in reads_ref:                          # a carrier here reads a different window
            hap_i = hap_i - hap
        status = ("unresolved" if unresolved else
                  "haplotype" if hap_i else "genotype" if sam_i else "chimera")
        rows.append({
            "canary_index": ci, "anchor_pos": c["anchor_pos"],
            "n_variants": len(sub_haps),
            "n_variant_sites_in_window": n_sites,
            "n_sites_reading_reference": len(reads_ref),
            "n_sites_unresolved": unresolved,
            # An indel or multi-allelic site inside the window that neither pass can classify. A
            # chromosome carrying one reads this window differently, so a window with any of
            # these is not VERIFIED to be a haplotype even when it is classified as one, and the
            # count below lets the claim be made over the windows that have none.
            "n_non_snv_sites_in_window": int(skipped_by_win.get(ci, 0)),
            "n_haplotypes_substituted_only": hap_sub_only,
            "n_haplotypes_carrying_all": len(hap_i),
            "n_individuals_carrying_all": len(sam_i),
            "anchor_n_carriers": ncar[0] if ncar else None,
            "max_af": max([a for a in afs if a is not None], default=None),
            "status": status,
            # The carrying samples' identifiers are deliberately NOT recorded. They were here
            # while this was being written and they do not belong in a published result file:
            # the paper's ethics statement says we publish no identifier that would link a
            # sequence to a sample, and a row holding `anchor_pos` beside the names of the five
            # people who carry it is exactly that, in one place and ready to use. The counts
            # above carry everything the claim needs.
        })

    by_status = collections.Counter(r["status"] for r in rows)
    comp = [r for r in rows if r["n_variants"] > 1]
    single = [r for r in rows if r["n_variants"] == 1]
    hap_rows = [r for r in rows if r["status"] == "haplotype"]
    moved = [r for r in rows
             if r["n_haplotypes_substituted_only"] > 0 and r["n_haplotypes_carrying_all"] == 0]
    return {
        "npz": str(npz), "vcf": str(vcf), "max_af": max_af, "max_carriers": max_carriers,
        "n_samples": len(samples), "n_rare_variants_in_table": len(tab),
        "n_canaries": len(rows), "windows_with_no_substituted_site": missing,
        "non_snv_sites_skipped_in_windows": int(sum(skipped_by_win.values())),
        "windows_with_a_non_snv_site": sum(1 for r in rows if r["n_non_snv_sites_in_window"]),
        "pure_snv_by_status": dict(collections.Counter(
            r["status"] for r in rows if not r["n_non_snv_sites_in_window"])),
        "by_status": dict(by_status),
        "reclassified_by_subtraction": {
            "n": len(moved),
            "note": "windows the substituted-sites-only intersection called haplotypes and that "
                    "no chromosome actually carries once the other variant sites in the window "
                    "are required to read reference",
        },
        "single_variant": {"n": len(single),
                           "by_status": dict(collections.Counter(r["status"] for r in single))},
        "composite": {"n": len(comp),
                      "by_status": dict(collections.Counter(r["status"] for r in comp))},
        "haplotypes_carrying_all": {
            "n": len(hap_rows),
            "hist": dict(sorted(collections.Counter(
                r["n_haplotypes_carrying_all"] for r in hap_rows).items())),
            "max": max([r["n_haplotypes_carrying_all"] for r in hap_rows], default=None),
        },
        "individuals_carrying_all_max": max(
            [r["n_individuals_carrying_all"] for r in rows if r["n_individuals_carrying_all"]],
            default=None),
        "variants_per_canary_hist": dict(sorted(collections.Counter(
            r["n_variants"] for r in rows).items())),
        "variant_sites_per_window_hist": dict(sorted(collections.Counter(
            r["n_variant_sites_in_window"] for r in rows).items())),
        "rows": rows,
    }

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--max_af", type=float, default=0.001)
    ap.add_argument("--max_carriers", type=int, default=5)
    ap.add_argument("--out", default="results/human_canary_haplotype_audit.json")
    a = ap.parse_args(argv)

    res = audit(Path(a.npz), Path(a.vcf), a.max_af, a.max_carriers)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(res, indent=1), encoding="utf-8")
    brief = {k: v for k, v in res.items() if k != "rows"}
    print(json.dumps(brief, indent=1))
    print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
