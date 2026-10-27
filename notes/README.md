# Working logs

The dated record of how the measurements in the paper were produced — every experiment, every
bug, and what each one changed. The paper's availability statement points here.

| | |
|---|---|
| `07_gate1_raw_reports.md` | The first gate: the reports that decided whether the full grid was worth running |
| `08_experiment_log_week1_EN.md` | Week 1 |
| `09_experiment_log_week2_EN.md` | Week 2 |
| `17_experiment_log_week3_EN.md` | Week 3, and the longest — the released-model arms, the four decoders, the human chromosome-22 arm, and the single-basis re-measurement |
| `figures/` | The generated per-sweep reports and trajectory plots the logs quote |

Entries are numbered `§<week>.<n>` and dated. A claim in the paper can be traced through a log
entry to the run directory under [`../results/`](../results/) that produced it, and from there to
the `args.json` the run was launched with.

## What these are, and are not

They are a laboratory notebook, written as the work happened and not rewritten afterwards. They
record results that did not survive, numbers that were wrong before they were right, and the
reasons for each change. Several entries exist only because a measurement was found to be wrong —
the single-basis re-measurement of 19–23 October is the largest of them, and the paper reports the
re-measured values throughout.

A recurring theme, recorded because it cost the most time: **the failures in this project are
silent.** A `.bashrc` that a non-interactive SSH never reads; a bare `%` eating a table row; a
dispatch keyed on a path; an exporter blind to re-scored files; a status line reporting the exit
code of `date`; a `\ref` that lost its backslash and printed as ordinary text through a clean
compile. None of them raised an error. The practice that catches them is the same every time:
count the output files, and read the output rather than the log.

Project planning, correspondence and talk material are not part of the record and are not
published here.
