# Changelog

## 0.1.1 (2026-08-13)

Correctness / performance / compatibility patch from the v0.1 code review
([REVIEW.md](REVIEW.md)). No coordinate-semantics or stats-algorithm changes.

### Fixed

- `iter_batches`: 1-thread mode returned only the first batch then EOF;
  multi-thread mode hash-deduplicated away identical records (cat-bam inputs),
  ignored `batch_size`, and materialized the whole file. Now: 1T streams
  batch-by-batch; MT processes bounded 1 Mb waves, honors `batch_size`, keeps
  exact duplicates, and both modes include the unmapped-unplaced tail via
  BAI pseudo-bin metadata. 1T == NT record-for-record.
- Re-iterating an `AlignmentFile` mid-stream (`iter(af)`) no longer silently
  drops records held in the prefetch buffer; iteration continues from the
  logical position (pysam semantics).
- `VariantFile.fetch(contig)` returned empty when the header lacked the
  contig length (`stop.unwrap_or(0)`); unbounded fetch now falls back to a
  sequential scan.
- Accessing `header` / `references` / `lengths` / `nreferences` after
  `close()` raised `PanicException`; now raises `ValueError`.
- `parallel_fetch` replaced hash-based dedup (2 String allocations per record,
  collision risk, duplicate loss) with region-merge + start-ownership
  filtering: overlapping regions are merged first, each record is emitted
  exactly once per region union, exact duplicates are preserved.
- Negative `start` / `stop` raised `OverflowError`; now `ValueError` like
  pysam. `stop` beyond the contig length is clamped (pysam behavior);
  out-of-range `start` yields empty results.

### Added

- `AlignedSegment.reference_end` (0-based exclusive; `None` for unmapped /
  empty CIGAR, matching pysam).
- `AlignmentFile.header` returns a pysam-style dict (`HD` / `SQ` / `RG` /
  `PG`) instead of a single-line string.
- B-array tags (`B,c` … `B,i`, `B,f`) round-trip as `array.array` with the
  original typecode — element- and type-equal to pysam (was: Rust debug text).
- `tests/parity/test_random_regions.py`: seeded random-region differential
  test vs pysam (count / fetch / coverage / depth / pileup at threads 1 & 4);
  fixture mode always runs, real-data mode behind `SAMRUST_REAL_DATA=1`.
- `tests/parity/test_regression_v011.py`: regression coverage for every fix
  above plus VCF no-length-contig fetch.
- `pytest --require-oracles`: turns oracle-missing skips into failures for
  CI / release gates.

### Changed

- `recount` / `parallel_map_regions` reuse one `IndexedAlignmentReader` per
  thread (`rayon::map_init`) instead of reopening BAM+index per chunk / site
  (16T recount was ~5.6× slower than necessary).
- `__next__` batch decode releases the GIL (`py.allow_threads`); batch and
  VCF iterators consume `std::vec::IntoIter` instead of cloning each record.
- Three-way benchmark (`scripts/run_three_way_benchmark.py`) measures each
  tool in an isolated subprocess (independent RSS high-water marks) and
  digests canonical outputs only, so digests match iff results match.
- CI pytest matrix now covers Python 3.10–3.13 (matches release wheels).
- Removed dead code: `config.rs`, `BoundedChannel`, `compat.py`; merged the
  duplicated `BASE_BUCKET` LUT into `samrust-core/src/base.rs`.

## 0.1.0 (2026-08-13)

First public tag (`v0.1.0`). Linux x86_64 wheels via GitHub Releases (not PyPI).

### Added

- M0 repository scaffold: Cargo workspace, PyO3/maturin package, CLI stub, CI.
- M1 test infrastructure: Tier-0 fixtures, real-data discovery, pysam/samtools baselines, metadata collector.
- M2 Rust BAM core: `AlignmentReader` linear iteration, `Record`/`Cigar`/`Tags`/`Header`, `samrust dump-records`, pysam parity gate.
- M3 Python `AlignmentFile` / `AlignedSegment` with batch iteration and GIL release on indexed paths.
- M4 Indexed `fetch` (BAI/CSI via noodles).
- M5 Parallel scheduler / `iter_batches` / `parallel_fetch` (rayon + crossbeam-channel).
- M6 `count`, `count_coverage`, `depth_blocks`, `depth_numpy` (serial + parallel, NumPy when available).
- M7 parallel `pileup_counts` (BQ/MAPQ/flag filters; indel-aware skip of del/refskip; serial==pysam normalized; 1T==NT).
- M8 `samrust recount` candidate-site utility + fungal `ALT_COUNT >= 10` validation gate / thread scaling.
- HPC helpers: `scripts/submit_m7_m8_validation.sh` / `run_m7_m8_heavy.sh` (Slurm partition `qcpu_18i`).
- M9 Python `VariantFile` / `VariantRecord` read path: VCF, VCF.gz+TBI, BCF+CSI; header samples/contigs; sequential iteration; indexed `fetch` (0-based half-open). `pos` is 1-based like pysam. Writer / caller not in scope. 16T three-way table: rubam cells **NA** (no VariantFile API).
- M11 CRAM evaluation: noodles `cram`/`fasta` features; `AlignmentFile` sequential iteration + indexed `fetch` on `.cram`+`.crai`+FASTA (mode `rb`/`rc`); pysam `"rc"` oracle. Stats / parallel APIs stay BAM-only (`NotImplemented` on CRAM). Optional `htslib-backend` designed, not implemented; does not block v0.1.
- M12 v0.1 release: GitHub distribution, manylinux wheel workflow, INSTALL/DEVELOPMENT/API_COMPATIBILITY/NOTICE, Interval `proptest` gate, published 16T table.

### Changed

- Stats hot path (`count` / `parallel_count` / `depth_profile` / `coverage_profile_with_filter` / `pileup_counts`) streams noodles BAM records via `IndexedAlignmentReader::for_each_raw` (scratch reuse). Accumulators walk CIGAR in place and do not materialize owned `Record` (qname / seq `String` / tags / CIGAR `Vec`). Python `fetch` / `AlignedSegment` still uses `fetch_records` / `Record::from_noodles`.
- M10 BAM efficiency vs rubam: stats paths (`count` / depth / coverage / `pileup_counts`) use `Scheduler::stats()` (1 chunk/thread, matching rubam `n.div_ceil(nt)`). Fetch still uses 4 chunks/thread. CIGAR M/=/X inner loops clip to the interval (`overlap_span`); depth range-increments the overlap slice; coverage/pileup classify bases with a 256-byte LUT. Formal 16T table: job 2312423 (`BENCHMARKS.md`).
- `DEVELOPMENT_PLAN.md` §1.2 / §15 / §19.6 / §21 / Rule 6: durable benchmark rule — **features and results track pysam; runtime tracks rubam** (pysam is the slow baseline, not the speed target). If rubam has no equivalent API or result, the 16T three-way table keeps the row and fills rubam cells with **NA**.
- Login-node 1T fungal 100kb timing after this change: `benchmark/results/hotpath_for_each_raw.1t_login_100kb.json` (not a replacement three-way table; host differs from prior Slurm `bnode17` numbers).

### Fixed

- `count` 1-thread vs N-thread fork: serial path dropped placed-unmapped mates (`BAM_FUNMAP` with POS set) that pysam `count`/`fetch` include. Serial and parallel now share `parallel_count` (fetch hits + start ownership) so 1T == NT == pysam `nofilter`.
- `depth_numpy` / `depth_blocks`: no longer count CIGAR `D`/`N`. Semantics match samtools depth and rubam `get_depths` (aligned M/=/X, including ambiguous bases). Three-way pysam oracle updated off `count_coverage` A+C+G+T.
