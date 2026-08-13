# Agent rules for SAMRust

Follow [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md). Project name is **SAMRust**
(crates/packages: `samrust*`).

## Rule 1

Execute **one milestone at a time**. Do not implement later milestones early.

## Rule 2

Before starting a milestone:

1. Read `DEVELOPMENT_PLAN.md`
2. Read current code and tests
3. Run baseline tests

## Rule 3

Each milestone must ship **implementation + tests + documentation** together.

## Rule 4

Performance code without benchmarks must not merge.

## Rule 5

Any parallel implementation must first prove:

```text
1-thread output == N-thread output
```

## Rule 6

Any pysam-compatible method must have a pysam oracle test.

Critical BAM analytics (`fetch` / `count` / `count_coverage` / depth /
`pileup_counts` / candidate recount) must also be compared against
[rubam](https://github.com/victormar1/rubam) for **output**, **runtime**, and
**resource usage** (RSS/CPU), with results under `benchmark/results/`. See
`DEVELOPMENT_PLAN.md` §1.2 / §15 / §19.6. Features and results track **pysam**;
runtime tracks **rubam**. If rubam lacks an equivalent API or result, keep the
workload row in the 16T three-way table and fill rubam cells with **NA**
(`wall_s` / RSS / digest / vs-rubam); document why — do not silently skip.

## Rule 7

Do not change coordinate semantics. Python API is always **0-based half-open**.
All `+1`/`-1` conversions belong in `coords.rs`.

## Rule 8

Do not implement hot loops in Python. Hot paths belong in Rust.

## Rule 9

Do not copy entire BAM / genome results to every worker for multithreading.

## Rule 10

Every new dependency needs a written rationale.

## Rule 11

Never modify user resequencing source data. Real-data tests are **read-only**.
Outputs go to `benchmark/results` or scratch.

## Rule 12

If real-data paths are missing: stop real-data benchmarks, continue synthetic/unit
work. Never fabricate benchmark results.

## Rule 13

Heavy compute (fungal BAM, full-genome / multi-thread scaling, M7–M8 recount
gates, Tier 2–4) must be submitted to the HPC Slurm partition **`qcpu_18i`**.
Do not run these on the login node. See `DEVELOPMENT_PLAN.md` §19.5.
Login node is OK only for Tier-0 fixtures, short unit/pytest, and compile/clippy.

## M0–M8 notes

- **M0** must not implement BAM algorithms, pileup, VCF, or parallel runtime.
- **M1** builds fixtures, discovery, and oracle baselines only — still no BAM engine.
- **M2** implements linear BAM decode only — no indexed fetch, pileup, or Python AlignmentFile.
- **M3–M6** deliver Python surface, indexed fetch, parallel runtime, and count/depth/coverage.
- **M7** delivers parallel `pileup_counts` with pysam-normalized parity and 1T==NT bit-exact.
- **M8** delivers `samrust recount` (benchmark utility only) and fungal `ALT_COUNT >= 10` gate.
- **M9** delivers Python `VariantFile` read path (VCF / VCF.gz / BCF, indexed fetch, pysam + bcftools oracle). No writer/caller. 16T three-way table: rubam = **NA**.
- **M10** delivers profiling-driven BAM efficiency vs rubam (stats: 1 chunk/thread; clipped CIGAR inner loops). Formal 16T three-way table on `qcpu_18i`.
- **M11** delivers CRAM evaluation: noodles sequential + CRAI fetch vs pysam; stats stay BAM-only; optional HTSlib backend designed not implemented. Does not block BAM/VCF v0.1.
- **M12** delivers v0.1 release (docs gate, Linux wheel, GitHub). Further work is out of the current stop point unless requested.
