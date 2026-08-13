//! Pileup base counts (M7).

use std::path::Path;

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::Flags;

use crate::base::BASE_BUCKET;
use crate::coords::Interval;
use crate::error::{Result, SamRustError};
use crate::indexed::{raw_alignment_start_0based, IndexedAlignmentReader};
use crate::parallel::{parallel_map_regions, Scheduler};
use crate::record::Record;

/// Per-position pileup counts (A/C/G/T/N + total depth).
#[derive(Debug, Clone, Default)]
pub struct PileupCounts {
    pub start: u64,
    pub a: Vec<u32>,
    pub c: Vec<u32>,
    pub g: Vec<u32>,
    pub t: Vec<u32>,
    pub n: Vec<u32>,
    pub depth: Vec<u32>,
}

impl PileupCounts {
    pub fn new(start: u64, len: usize) -> Self {
        Self {
            start,
            a: vec![0; len],
            c: vec![0; len],
            g: vec![0; len],
            t: vec![0; len],
            n: vec![0; len],
            depth: vec![0; len],
        }
    }

    pub fn len(&self) -> usize {
        self.depth.len()
    }

    pub fn is_empty(&self) -> bool {
        self.depth.is_empty()
    }

    pub fn merge(&mut self, other: &Self) {
        assert_eq!(self.start, other.start);
        assert_eq!(self.depth.len(), other.depth.len());
        for i in 0..self.depth.len() {
            self.a[i] += other.a[i];
            self.c[i] += other.c[i];
            self.g[i] += other.g[i];
            self.t[i] += other.t[i];
            self.n[i] += other.n[i];
            self.depth[i] += other.depth[i];
        }
    }
}

/// Flag / quality filters for pileup (defaults exclude secondary/supplementary/qcfail/duplicate/unmapped).
#[derive(Debug, Clone, Copy)]
pub struct PileupFilter {
    pub min_base_quality: u8,
    pub min_mapping_quality: u8,
    pub exclude_unmapped: bool,
    pub exclude_secondary: bool,
    pub exclude_supplementary: bool,
    pub exclude_qcfail: bool,
    pub exclude_duplicate: bool,
}

impl Default for PileupFilter {
    fn default() -> Self {
        Self {
            min_base_quality: 0,
            min_mapping_quality: 0,
            exclude_unmapped: true,
            exclude_secondary: true,
            exclude_supplementary: true,
            exclude_qcfail: true,
            exclude_duplicate: true,
        }
    }
}

impl PileupFilter {
    pub fn passes(&self, record: &Record) -> bool {
        self.passes_flags_mapq(Flags::from(record.flag()), record.mapping_quality())
    }

    pub fn passes_raw(&self, raw: &bam::Record) -> bool {
        // Missing MAPQ → 255 (pysam / `Record::from_noodles`), not rubam's 0.
        let mapq = raw.mapping_quality().map(u8::from).unwrap_or(255);
        self.passes_flags_mapq(raw.flags(), mapq)
    }

    fn passes_flags_mapq(&self, flags: Flags, mapq: u8) -> bool {
        if self.exclude_unmapped && flags.is_unmapped() {
            return false;
        }
        if self.exclude_secondary && flags.is_secondary() {
            return false;
        }
        if self.exclude_supplementary && flags.is_supplementary() {
            return false;
        }
        if self.exclude_qcfail && flags.is_qc_fail() {
            return false;
        }
        if self.exclude_duplicate && flags.is_duplicate() {
            return false;
        }
        mapq >= self.min_mapping_quality
    }
}

/// Accumulate pileup counts for one read into `counts`.
pub fn add_record_pileup(
    raw: &bam::Record,
    interval: Interval,
    counts: &mut PileupCounts,
    filter: PileupFilter,
) -> Result<()> {
    if !filter.passes_raw(raw) {
        return Ok(());
    }
    if raw.flags().is_unmapped() {
        return Ok(());
    }
    let start = raw_alignment_start_0based(raw)?;
    if start < 0 {
        return Ok(());
    }

    let seq = raw.sequence();
    let quals = raw.quality_scores().as_bytes();
    let mut ref_pos = start as u64;
    let mut query_pos = 0usize;

    for item in raw.cigar().iter() {
        if ref_pos >= interval.stop.0 {
            break;
        }
        let op = item.map_err(SamRustError::from)?;
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if let Some((lo, hi)) = interval.overlap_span(ref_pos, len) {
                    let skip = (lo - ref_pos) as usize;
                    let n = (hi - lo) as usize;
                    let idx0 = (lo - interval.start.0) as usize;
                    for j in 0..n {
                        let idx = idx0 + j;
                        if idx >= counts.depth.len() {
                            break;
                        }
                        let qidx = query_pos + skip + j;
                        if qidx < quals.len() && quals[qidx] >= filter.min_base_quality {
                            counts.depth[idx] += 1;
                            if let Some(base) = seq.get(qidx) {
                                increment_pileup_base(counts, idx, base);
                            }
                        }
                    }
                }
                ref_pos += len as u64;
                query_pos += len;
            }
            Kind::Insertion | Kind::SoftClip => query_pos += len,
            // Deletion / ref-skip: advance reference only. Base-depth pileup skips these
            // (matches pysam pileup `is_del` / `is_refskip` exclusion for base counts).
            Kind::Deletion | Kind::Skip => {
                ref_pos += len as u64;
            }
            _ => {}
        }
    }
    Ok(())
}

fn increment_pileup_base(counts: &mut PileupCounts, idx: usize, base: u8) {
    match BASE_BUCKET[base as usize] {
        0 => counts.a[idx] += 1,
        1 => counts.c[idx] += 1,
        2 => counts.g[idx] += 1,
        3 => counts.t[idx] += 1,
        _ => counts.n[idx] += 1,
    }
}

/// 0=A, 1=C, 2=G, 3=T, 4=other.
/// Serial pileup counts for a region.
pub fn pileup_counts(
    reader: &mut IndexedAlignmentReader,
    contig: &str,
    interval: Interval,
    filter: PileupFilter,
) -> Result<PileupCounts> {
    let len = interval.len() as usize;
    let mut counts = PileupCounts::new(interval.start.0, len);
    reader.for_each_raw(contig, interval, |raw| {
        add_record_pileup(raw, interval, &mut counts, filter)
    })?;
    Ok(counts)
}

/// Parallel pileup counts; stitches owned-interval chunk results into the parent region.
pub fn parallel_pileup_counts(
    bam_path: &Path,
    contig: &str,
    interval: Interval,
    filter: PileupFilter,
    threads: usize,
) -> Result<PileupCounts> {
    let chunks =
        Scheduler::stats().chunk_interval(contig, interval.start.0, interval.stop.0, threads)?;
    let partial: Vec<PileupCounts> =
        parallel_map_regions(bam_path, chunks, threads, true, |_header, reader, chunk| {
            pileup_counts(reader, &chunk.contig, chunk.interval, filter)
        })?;
    Ok(stitch_pileup_counts(interval, partial))
}

fn stitch_pileup_counts(parent: Interval, profiles: Vec<PileupCounts>) -> PileupCounts {
    let mut out = PileupCounts::new(parent.start.0, parent.len() as usize);
    for p in profiles {
        for i in 0..p.len() {
            let abs = p.start + i as u64;
            if abs >= parent.start.0 && abs < parent.stop.0 {
                let idx = (abs - parent.start.0) as usize;
                out.a[idx] += p.a[i];
                out.c[idx] += p.c[i];
                out.g[idx] += p.g[i];
                out.t[idx] += p.t[i];
                out.n[idx] += p.n[i];
                out.depth[idx] += p.depth[i];
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::IndexedAlignmentReader;

    fn fixture_bam() -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/small.bam")
    }

    #[test]
    fn pileup_indel_adjacent_nonzero_and_parallel_bit_exact() {
        let bam = fixture_bam();
        let iv = Interval::new(180, 230).unwrap();
        let filter = PileupFilter::default();
        let mut reader = IndexedAlignmentReader::open(&bam).unwrap();
        let serial = pileup_counts(&mut reader, "chr1", iv, filter).unwrap();
        assert!(serial.depth.iter().any(|&d| d > 0));

        let parallel = parallel_pileup_counts(&bam, "chr1", iv, filter, 4).unwrap();
        assert_eq!(serial.a, parallel.a);
        assert_eq!(serial.c, parallel.c);
        assert_eq!(serial.g, parallel.g);
        assert_eq!(serial.t, parallel.t);
        assert_eq!(serial.n, parallel.n);
        assert_eq!(serial.depth, parallel.depth);
    }

    #[test]
    fn pileup_bq_filter_reduces_or_keeps_depth() {
        let bam = fixture_bam();
        let iv = Interval::new(0, 100).unwrap();
        let mut reader = IndexedAlignmentReader::open(&bam).unwrap();
        let loose = pileup_counts(
            &mut reader,
            "chr1",
            iv,
            PileupFilter {
                min_base_quality: 0,
                ..PileupFilter::default()
            },
        )
        .unwrap();
        let mut reader = IndexedAlignmentReader::open(&bam).unwrap();
        let strict = pileup_counts(
            &mut reader,
            "chr1",
            iv,
            PileupFilter {
                min_base_quality: 30,
                ..PileupFilter::default()
            },
        )
        .unwrap();
        let sum = |v: &[u32]| v.iter().map(|&x| u64::from(x)).sum::<u64>();
        assert!(sum(&strict.depth) <= sum(&loose.depth));
    }

    #[test]
    fn pileup_mapq_filter_excludes_low_mapq_region_reads() {
        let bam = fixture_bam();
        // Fixture includes low-MAPQ / secondary neighborhood near 590+.
        let iv = Interval::new(590, 650).unwrap();
        let mut reader = IndexedAlignmentReader::open(&bam).unwrap();
        let all_mapq = pileup_counts(
            &mut reader,
            "chr1",
            iv,
            PileupFilter {
                min_mapping_quality: 0,
                exclude_secondary: false,
                exclude_supplementary: false,
                exclude_duplicate: false,
                exclude_qcfail: false,
                ..PileupFilter::default()
            },
        )
        .unwrap();
        let mut reader = IndexedAlignmentReader::open(&bam).unwrap();
        let high_mapq = pileup_counts(
            &mut reader,
            "chr1",
            iv,
            PileupFilter {
                min_mapping_quality: 60,
                exclude_secondary: false,
                exclude_supplementary: false,
                exclude_duplicate: false,
                exclude_qcfail: false,
                ..PileupFilter::default()
            },
        )
        .unwrap();
        let sum = |v: &[u32]| v.iter().map(|&x| u64::from(x)).sum::<u64>();
        assert!(sum(&high_mapq.depth) <= sum(&all_mapq.depth));
    }

    #[test]
    fn pileup_filter_passes_raw_matches_owned_and_missing_mapq_is_255() {
        let bam = fixture_bam();
        let iv = Interval::new(0, 1000).unwrap();
        let mut reader = IndexedAlignmentReader::open(&bam).unwrap();
        let owned = reader.fetch_records("chr1", iv).unwrap();
        let filter = PileupFilter {
            min_mapping_quality: 60,
            ..PileupFilter::default()
        };
        let mut n = 0usize;
        reader
            .for_each_raw("chr1", iv, |raw| {
                assert_eq!(filter.passes_raw(raw), filter.passes(&owned[n]));
                n += 1;
                Ok(())
            })
            .unwrap();
        assert_eq!(n, owned.len());
        // Missing MAPQ → 255 (pysam / Record::from_noodles), so min_mapq=255 still passes.
        let missing = bam::Record::default();
        assert!(missing.mapping_quality().is_none());
        let keep_missing = PileupFilter {
            min_mapping_quality: 255,
            exclude_unmapped: false,
            ..PileupFilter::default()
        };
        assert!(keep_missing.passes_raw(&missing));
    }
}
