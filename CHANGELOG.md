# Changelog

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
