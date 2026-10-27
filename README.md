# Identification without Reproduction

Code and measurements for *Identification without Reproduction: Privacy Leakage from Fine-Tuned
Genomic Language Models*.

A genomic language model can stop reproducing a training record and go on identifying the person it
came from. We measure both on eight released checkpoints from four families — Nucleotide Transformer
(250M, 500M), the DNABERT *k*-mer family (3/4/5/6-mer), DNABERT-2 and HyenaDNA — fine-tuned on
*E. coli*, yeast, and human chromosome 22 with planted records, and on one small backbone we train
ourselves so that tokenisation can be varied with everything else held fixed.

Three things it reports:

- **Released checkpoints leak at a held-out loss a practitioner would accept.** The cheapest costs
  0.007 bits/nt above the run's own optimum and returns 30% of planted records.
- **Reproduction is the wrong alarm.** At checkpoints where verbatim recovery reads exactly zero in
  every seed, the same model still ranks one individual's chromosome-22 haplotype first against 100
  real rare segments 0.92–0.96 of the time, against a chance level of 0.0099.
- **Whether a record can be reached at all is set by how many token slots it occupies**, not by how
  much the model knows: models that answer any single hidden nucleotide identically return 92% or
  0% of whole records. This is why no single extraction test is correct across tokenizers.

Every input is public: 1000 Genomes phase 3 (chr22, GRCh38), the GUE benchmark, and the reference
genomes of *E. coli* K-12 MG1655 and *S. cerevisiae* S288C.

## Start here

Every figure and table in the paper is drawn from one published file, so the quickest check needs no
GPU, no model download and no genome:

```bash
pip install matplotlib     # the only dependency; the csv is read with the stdlib
python -m pilot2.figures_privacy     # results/ -> paper/figures
python -m pilot2.appendix_tables     # results/ -> paper/tables
```

`results/figure_data.csv` is 27,106 rows, one per (run, snapshot, tier, host); see
[`results/README.md`](results/README.md) for the columns and for which directory holds which cell.
The measured JSON behind each row is in `results/<sweep>/<run>/`, next to the `args.json` the run was
launched with, so any number in the paper can be traced to the file it came from without rerunning
anything.

## What is where

| | |
|---|---|
| `pilot2/` | All code: experiments, the four decoders, scoring, figures, tables |
| `results/` | The measured output behind every table and figure, plus `figure_data.csv` |
| `notes/` | Dated working logs — every experiment, every bug, and what each one changed |
| `pilot/` | A frozen first pilot, kept because the week-1 log cites it |

## What is deliberately not here

**The fine-tuned checkpoints**, because they contain the planted records.

**The decoders' predicted and target nucleotide strings.** For the human arms a target is a 48 nt
window of chromosome 22 written with a real carrier's rare alleles, so publishing one would publish a
reconstructed haplotype. They are stripped from every result file rather than from the human ones
only: classifying by directory name is a mistake this project has already made once, and the cost of
making it here is a real person's sequence. Everything the paper measures is kept — `exact`,
`hamming_nt`, the per-restart Hamming vector and the leave-one-out bounds — and the generated tables
are byte-identical with the strings removed.

**The canary files.** `pilot2/human_canary.py` rebuilds them from the public VCF.

## Reproducing the measurements

These need GPUs. Every number in the paper came from 4×RTX 4090 (24 GB); a single 24 GB card is
enough for any one cell. Set `GLMEXT_DATA` to a directory holding the corpora — it is read from the
environment and not from a shell profile, because a non-interactive SSH does not read one.

```bash
pip install -r requirements.txt
export GLMEXT_DATA=/path/to/data
```

**Build the human canaries.** Each is a 96 nt window of chromosome 22 built around a variant carried
by at most five of the 2,548 individuals in 1000 Genomes phase 3, written with that carrier's own
alleles. Both commands read from `$GLMEXT_DATA/human/`, so neither needs an argument:

```bash
python -m pilot2.human_canary    # -> $GLMEXT_DATA/human/canaries_chr22.npz
python -m pilot2.canary_census   # -> results/human_canary_census.json
```

**Fine-tune a released checkpoint with those records planted**, saving every snapshot so a decoder
can be asked about any round afterwards:

```bash
python -m pilot2.nt --size 250m --data human --epochs 50 --mask_rate 0.15 \
    --canary_npz $GLMEXT_DATA/human/canaries_chr22.npz \
    --probes_per_tier 50 --pool 100 --seed 0 --out outputs/nt_human_s0
```

`pilot2/dnabert.py`, `pilot2/dnabert2.py` and `pilot2/hyena.py` take the same shape for the other
families; `pilot2/run_pilot.py` runs one (tokenizer, objective, seed) cell of the controlled
backbone. The `pilot2/queue_*.sh` scripts are the exact invocations each arm was launched with,
dated, with the reasoning for the run in the header.

**Decode.** Which decoder is correct depends on the model's vocabulary, and this is the point of the
paper rather than an implementation detail:

| decoder | for | |
|---|---|---|
| `pilot2/recon_audit.py` | fixed-width non-overlapping tokens | re-masking ICM, 128 restarts |
| `pilot2/recon_overlap.py` | overlapping *k*-mers (DNABERT) | decodes at nucleotide level |
| `pilot2/recon_bpe.py` | variable-length vocabularies (DNABERT-2) | span-aware |
| — | prefix-conditioned | the protocol the field uses, run inline by each training script |

Each reports a leave-one-out ceiling alongside its own result: per-position top-1 accuracy raised to
the slot count, which bounds what *any* decoder of that kind could recover. A zero below a high
ceiling is a limit of the attack, not a property of the model, and the paper says so where that
happens.

**Score an existing checkpoint again** without retraining it — this is how the single-basis
re-measurement in the appendix was produced:

```bash
python -m pilot2.recon_audit --root outputs/nt_human_s0 --ckpt ckpt_ep50.pt --restarts 128
python -m pilot2.rescore --root outputs/nt_human_s0 --ckpt ckpt_ep50.pt --name scores_ep50_rb
python -m pilot2.utility_at_leak --root outputs/nt_human_s0   # downstream score at the same points
```

Note the checkpoint naming: Nucleotide Transformer writes `ckpt_ep<N>.pt` where the controlled
backbone writes `ep<N>.pt`. Passing the wrong one finds no run and returns quietly.

**Rebuild the derived table** after new runs:

```bash
python -m pilot2.export_figure_data --root outputs --out results/figure_data.csv
```

## Two cautions about the numbers

**Extraction rates carry more uncertainty than the probe count suggests.** The re-masking decode is a
discrete iteration and is sensitive to the model's last-bit arithmetic: re-running an identical cell
after the host's floating-point behaviour changed moved its extraction rate by six points, comparable
to the binomial standard error on 50 probes. Identification does not share this sensitivity. Expect
agreement on identification to the digit and on extraction to within several points.

**Read figures against held-out loss, not epoch count.** Tokenizers converge at very different rates,
so comparing them at a matched epoch compares two different amounts of training.

## Credit

This work began as an extension of the **GLM-Memorization** framework by Alexander Nemecek, which
implements the released-model backends, the canary machinery and the baseline attacks. The code here
is a separate implementation under `pilot2/` and imports nothing from that framework, so it is not
redistributed here. It follows [Genome-Factory](http://arxiv.org/abs/2509.12266) for the model
backends and tokenizer conventions. The predecessor paper is Nemecek, Li, Jiang, Vaidya and Ayday,
*Quantifying memorization and privacy risks in genomic language models*, ACM-BCB 2026,
[10.1145/3807503.3819470](https://doi.org/10.1145/3807503.3819470).

Apache License 2.0, as in [LICENSE](LICENSE).
