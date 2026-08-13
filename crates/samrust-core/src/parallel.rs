//! Parallel region scheduler and merge utilities (M5).
//!
//! `rayon`: data-parallel worker pool for region chunks.
//! `crossbeam-channel`: reserved for bounded producer/consumer pipelines (v0.2).

use std::path::Path;

use rayon::prelude::*;

use crate::coords::Interval;
use crate::error::{Result, SamRustError};
use crate::header::Header;
use crate::indexed::IndexedAlignmentReader;
use crate::record::Record;

/// A genomic sub-region assigned to one worker.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegionChunk {
    pub contig: String,
    pub interval: Interval,
    pub chunk_id: usize,
}

/// Adaptive region chunking.
///
/// Fetch / batch engines keep ~4 chunks per thread (plan §8.3) so long
/// contigs can rebalance. Count / depth / coverage / pileup use
/// [`Scheduler::stats`]: one chunk per thread, matching rubam `n.div_ceil(nt)`.
#[derive(Debug, Clone)]
pub struct Scheduler {
    pub target_chunks_per_thread: usize,
}

impl Default for Scheduler {
    fn default() -> Self {
        Self {
            target_chunks_per_thread: 4,
        }
    }
}

impl Scheduler {
    pub fn new(target_chunks_per_thread: usize) -> Self {
        Self {
            target_chunks_per_thread: target_chunks_per_thread.max(1),
        }
    }

    /// Stats paths: one indexed query per worker.
    ///
    /// M10 evidence (fungal 100 kb, job 2312392): `target_chunks_per_thread=4`
    /// opened 64 BAI queries at 16T and depth *regressed* 8T→16T
    /// (0.091 s → 0.102 s). rubam fast mode uses 1 chunk/thread and plateaued
    /// (~0.058 s). Fetch still uses [`Scheduler::default`].
    pub fn stats() -> Self {
        Self {
            target_chunks_per_thread: 1,
        }
    }

    /// Split `[start, stop)` into `threads * target_chunks_per_thread` chunks (minimum 1).
    pub fn chunk_interval(
        &self,
        contig: impl Into<String>,
        start: u64,
        stop: u64,
        threads: usize,
    ) -> Result<Vec<RegionChunk>> {
        let interval = Interval::new(start, stop)?;
        if interval.is_empty() {
            return Ok(Vec::new());
        }
        let contig = contig.into();
        let threads = threads.max(1);
        let n_chunks = (threads * self.target_chunks_per_thread).max(1);
        let len = interval.len();
        let chunk_size = len.div_ceil(n_chunks as u64);
        let mut chunks = Vec::new();
        let mut pos = interval.start.0;
        let mut chunk_id = 0usize;
        while pos < interval.stop.0 {
            let end = (pos + chunk_size).min(interval.stop.0);
            chunks.push(RegionChunk {
                contig: contig.clone(),
                interval: Interval::new(pos, end)?,
                chunk_id,
            });
            chunk_id += 1;
            pos = end;
        }
        Ok(chunks)
    }
}

/// Ordered merge of per-chunk results keyed by `chunk_id`.
pub fn ordered_merge<T: Send>(mut chunks: Vec<(usize, T)>, ordered: bool) -> Vec<T> {
    if ordered {
        chunks.sort_by_key(|(id, _)| *id);
    }
    chunks.into_iter().map(|(_, v)| v).collect()
}

/// Exactly-once ownership test shared by parallel count / fetch paths.
///
/// A record belongs to `chunk` iff its 0-based alignment start (clamped up to
/// `parent.start`, so reads hanging left of the parent region are owned by the
/// first chunk) falls inside `chunk`. Records without an alignment start
/// (unmapped, no POS) are never owned: they do not appear in indexed region
/// queries.
pub(crate) fn start_owned_by_interval(start: i64, parent: &Interval, chunk: &Interval) -> bool {
    if start < 0 {
        return false;
    }
    let s = (start as u64).max(parent.start.0);
    chunk.contains(s)
}

/// Parallel map over region chunks.
///
/// Each rayon worker opens its indexed reader once (`map_init`) and reuses it
/// across all chunks it processes — previously every chunk reopened the BAM
/// and its index.
pub fn parallel_map_regions<T, F>(
    bam_path: &Path,
    chunks: Vec<RegionChunk>,
    threads: usize,
    ordered: bool,
    f: F,
) -> Result<Vec<T>>
where
    T: Send,
    F: Fn(&Header, &mut IndexedAlignmentReader, &RegionChunk) -> Result<T> + Send + Sync,
{
    if chunks.is_empty() {
        return Ok(Vec::new());
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .map_err(|e| SamRustError::InvalidArgument(e.to_string()))?;

    let results: Result<Vec<(usize, T)>> = pool.install(|| {
        chunks
            .par_iter()
            .map_init(
                || IndexedAlignmentReader::open(bam_path),
                |slot, chunk| {
                    let reader = slot.as_mut().map_err(|e| {
                        SamRustError::InvalidArgument(format!("reopen {}: {e}", bam_path.display()))
                    })?;
                    let header = reader.header().clone();
                    let value = f(&header, reader, chunk)?;
                    Ok((chunk.chunk_id, value))
                },
            )
            .collect()
    });
    let pairs = results?;
    Ok(ordered_merge(pairs, ordered))
}

/// A fetch unit: `chunk` partitions `parent`, and ownership is decided by
/// clamping record starts to `parent` (see [`start_owned_by_chunk`]).
#[derive(Debug, Clone)]
pub struct FetchWindow {
    pub parent: RegionChunk,
    pub chunk: RegionChunk,
}

/// Chunk length for whole-file parallel iteration.
///
/// Bounds per-wave memory of `iter_batches(threads>1)`: each in-flight chunk
/// holds at most the records overlapping ~1 Mb of reference.
pub const WHOLE_FILE_CHUNK_LEN: u64 = 1_000_000;

/// Build windows covering every placed record of the file, in file order
/// (contigs in header order, then by position). Unmapped records without a
/// position are not indexed and must be collected separately via
/// [`IndexedAlignmentReader::unmapped_tail_records`].
pub fn whole_file_windows(header: &Header) -> Vec<FetchWindow> {
    let mut windows = Vec::new();
    let mut next_id = 0usize;
    for (contig, &len) in header.references().iter().zip(header.lengths()) {
        if len == 0 {
            continue;
        }
        let parent = RegionChunk {
            contig: contig.clone(),
            interval: Interval::new(0, len).expect("nonzero contig length"),
            chunk_id: 0,
        };
        let mut pos = 0u64;
        while pos < len {
            let end = (pos + WHOLE_FILE_CHUNK_LEN).min(len);
            windows.push(FetchWindow {
                parent: parent.clone(),
                chunk: RegionChunk {
                    contig: contig.clone(),
                    interval: Interval::new(pos, end).expect("pos < end"),
                    chunk_id: next_id,
                },
            });
            next_id += 1;
            pos = end;
        }
    }
    windows
}

/// Fetch one wave of windows in parallel; per-window records come back in
/// wave order, each record exactly once (ownership-filtered).
pub fn parallel_fetch_wave(
    bam_path: &Path,
    wave: &[FetchWindow],
    threads: usize,
) -> Result<Vec<Vec<Record>>> {
    if wave.is_empty() {
        return Ok(Vec::new());
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .map_err(|e| SamRustError::InvalidArgument(e.to_string()))?;
    let results: Vec<Result<Vec<Record>>> = pool.install(|| {
        wave.par_iter()
            .map_init(
                || IndexedAlignmentReader::open(bam_path),
                |slot, w| {
                    let reader = slot.as_mut().map_err(|e| {
                        SamRustError::InvalidArgument(format!("reopen {}: {e}", bam_path.display()))
                    })?;
                    let mut records = reader.fetch_records(&w.chunk.contig, w.chunk.interval)?;
                    records.retain(|r| {
                        start_owned_by_interval(
                            r.reference_start(),
                            &w.parent.interval,
                            &w.chunk.interval,
                        )
                    });
                    Ok(records)
                },
            )
            .collect()
    });
    results.into_iter().collect()
}

/// Parallel fetch over arbitrary regions.
///
/// Overlapping / adjacent regions on the same contig are merged before
/// chunking, and each record is emitted exactly once via positional ownership
/// (no hashing, no false-positive dedup of identical records). Output order is
/// genomic: contigs in header order, then by position.
pub fn parallel_fetch_regions(
    bam_path: &Path,
    regions: &[(String, Interval)],
    threads: usize,
) -> Result<Vec<Record>> {
    if regions.is_empty() {
        return Ok(Vec::new());
    }
    let header = {
        let reader = IndexedAlignmentReader::open(bam_path)?;
        reader.header().clone()
    };

    // Group by contig, clamp to contig length, sort, merge overlapping/adjacent.
    let n_contigs = header.nreferences();
    let mut by_contig: Vec<Vec<Interval>> = vec![Vec::new(); n_contigs];
    for (contig, interval) in regions {
        let id = header.reference_id(contig)? as usize;
        let len = header.lengths()[id];
        let start = interval.start.0.min(len);
        let stop = interval.stop.0.min(len);
        if start < stop {
            by_contig[id].push(Interval::new(start, stop)?);
        }
    }

    let scheduler = Scheduler::default();
    let mut windows = Vec::new();
    let mut next_id = 0usize;
    for (id, mut ivs) in by_contig.into_iter().enumerate() {
        if ivs.is_empty() {
            continue;
        }
        ivs.sort_by_key(|iv| iv.start.0);
        let contig = header.references()[id].clone();
        let mut merged: Vec<Interval> = Vec::with_capacity(ivs.len());
        for iv in ivs {
            if let Some(last) = merged.last_mut() {
                if iv.start.0 <= last.stop.0 {
                    if iv.stop.0 > last.stop.0 {
                        *last = Interval::new(last.start.0, iv.stop.0)?;
                    }
                    continue;
                }
            }
            merged.push(iv);
        }
        for m in merged {
            let parent = RegionChunk {
                contig: contig.clone(),
                interval: m,
                chunk_id: 0,
            };
            for c in scheduler.chunk_interval(&contig, m.start.0, m.stop.0, threads)? {
                windows.push(FetchWindow {
                    parent: parent.clone(),
                    chunk: RegionChunk {
                        chunk_id: next_id,
                        ..c
                    },
                });
                next_id += 1;
            }
        }
    }

    let mut out = Vec::new();
    for records in parallel_fetch_wave(bam_path, &windows, threads)? {
        out.extend(records);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scheduler_splits_interval() {
        let sched = Scheduler::default();
        let chunks = sched.chunk_interval("chr1", 0, 1000, 2).unwrap();
        assert!(chunks.len() >= 2);
        assert_eq!(chunks.first().unwrap().interval.start.0, 0);
        assert_eq!(chunks.last().unwrap().interval.stop.0, 1000);
        assert_eq!(chunks.len(), 8, "fetch default is 4 chunks/thread");
    }

    #[test]
    fn stats_scheduler_is_one_chunk_per_thread() {
        let chunks = Scheduler::stats()
            .chunk_interval("chr1", 0, 100_000, 16)
            .unwrap();
        assert_eq!(chunks.len(), 16);
        assert_eq!(chunks[0].interval.start.0, 0);
        assert_eq!(chunks[15].interval.stop.0, 100_000);
    }

    #[test]
    fn ordered_merge_sorts_by_chunk_id() {
        let merged = ordered_merge(vec![(2, "c"), (0, "a"), (1, "b")], true);
        assert_eq!(merged, vec!["a", "b", "c"]);
    }

    #[test]
    fn ownership_clamps_to_parent_start() {
        let parent = Interval::new(100, 300).unwrap();
        let first = Interval::new(100, 200).unwrap();
        let second = Interval::new(200, 300).unwrap();
        // Read starting left of the parent belongs to the first chunk.
        assert!(start_owned_by_interval(50, &parent, &first));
        assert!(!start_owned_by_interval(50, &parent, &second));
        assert!(start_owned_by_interval(250, &parent, &second));
        assert!(!start_owned_by_interval(250, &parent, &first));
        // Unmapped (no position) is never owned.
        assert!(!start_owned_by_interval(-1, &parent, &first));
    }
}
