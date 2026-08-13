//! Parallel region scheduler and merge utilities (M5).
//!
//! `rayon`: data-parallel worker pool for region chunks.
//! `crossbeam-channel`: bounded channels for backpressure between producers/consumers.

use std::path::Path;

use crossbeam_channel::{bounded, Receiver, Sender};
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

/// Bounded channel wrapper (default capacity `2 * workers`).
#[derive(Debug)]
pub struct BoundedChannel<T> {
    pub tx: Sender<T>,
    pub rx: Receiver<T>,
}

impl<T> BoundedChannel<T> {
    pub fn new(capacity: usize) -> Self {
        let (tx, rx) = bounded(capacity);
        Self { tx, rx }
    }

    pub fn with_workers(workers: usize) -> Self {
        Self::new(workers.saturating_mul(2).max(2))
    }
}

/// Ordered merge of per-chunk results keyed by `chunk_id`.
pub fn ordered_merge<T: Send>(mut chunks: Vec<(usize, T)>, ordered: bool) -> Vec<T> {
    if ordered {
        chunks.sort_by_key(|(id, _)| *id);
    }
    chunks.into_iter().map(|(_, v)| v).collect()
}

/// Batch of decoded records (parallel batch engine).
#[derive(Debug, Clone, Default)]
pub struct RecordBatch {
    pub chunk_id: usize,
    pub records: Vec<Record>,
}

/// Parallel map over region chunks; each worker opens its own indexed reader.
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
            .map(|chunk| {
                let mut reader = IndexedAlignmentReader::open(bam_path)?;
                let header = reader.header().clone();
                let value = f(&header, &mut reader, chunk)?;
                Ok((chunk.chunk_id, value))
            })
            .collect()
    });
    let pairs = results?;
    Ok(ordered_merge(pairs, ordered))
}

/// Parallel fetch with deduplication at chunk boundaries.
pub fn parallel_fetch_records(
    bam_path: &Path,
    chunks: Vec<RegionChunk>,
    threads: usize,
    ordered: bool,
) -> Result<Vec<Record>> {
    let batches: Vec<RecordBatch> = parallel_map_regions(
        bam_path,
        chunks,
        threads,
        ordered,
        |_header, reader, chunk| {
            let records = reader.fetch_records(&chunk.contig, chunk.interval)?;
            Ok(RecordBatch {
                chunk_id: chunk.chunk_id,
                records,
            })
        },
    )?;
    let mut seen = std::collections::HashSet::new();
    let mut out = Vec::new();
    for batch in batches {
        for record in batch.records {
            let key = record_dedup_key(&record);
            if seen.insert(key) {
                out.push(record);
            }
        }
    }
    Ok(out)
}

fn record_dedup_key(record: &Record) -> (String, u16, i32, i64, Option<String>) {
    (
        record.query_name().to_string(),
        record.flag(),
        record.reference_id(),
        record.reference_start(),
        record.cigarstring(),
    )
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
}
