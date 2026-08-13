//! Read counting, depth, and coverage (M6).

use std::path::Path;

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::Flags;

use crate::coords::Interval;
use crate::error::{Result, SamRustError};
use crate::indexed::{raw_alignment_start_0based, IndexedAlignmentReader};
use crate::parallel::{parallel_map_regions, Scheduler};
use crate::record::Record;

/// Read filter presets mirroring common pysam `read_callback` names.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ReadFilter {
    /// Count all overlapping reads (pysam default / `nofilter`).
    #[default]
    NoFilter,
    /// Exclude unmapped, secondary, QC-fail, duplicate (pysam `all`).
    All,
}

impl ReadFilter {
    pub fn passes(&self, record: &Record) -> bool {
        self.passes_flags(Flags::from(record.flag()))
    }

    pub fn passes_raw(&self, raw: &bam::Record) -> bool {
        self.passes_flags(raw.flags())
    }

    fn passes_flags(&self, flags: Flags) -> bool {
        match self {
            Self::NoFilter => true,
            Self::All => {
                !(flags.is_unmapped()
                    || flags.is_secondary()
                    || flags.is_qc_fail()
                    || flags.is_duplicate())
            }
        }
    }
}

/// Count records overlapping a region (pysam `AlignmentFile.count` semantics).
///
/// Indexed fetch already applies BAM overlap (including **placed unmapped**
/// mates that still have POS set). `nofilter` therefore counts every fetch hit;
/// do not re-filter with [`crate::indexed::record_overlaps_interval`], which
/// drops `BAM_FUNMAP` and under-counts relative to pysam.
pub fn count(
    reader: &mut IndexedAlignmentReader,
    contig: &str,
    interval: Interval,
    filter: ReadFilter,
) -> Result<u64> {
    let mut n = 0u64;
    reader.for_each_raw(contig, interval, |raw| {
        if filter.passes_raw(raw) {
            n += 1;
        }
        Ok(())
    })?;
    Ok(n)
}

/// Parallel `count` with the same result as [`count`] for every thread count.
///
/// Each overlapping record is owned by one chunk via alignment start (clamped
/// into the parent interval) so spanning reads are not double-counted.
pub fn parallel_count(
    bam_path: &Path,
    contig: &str,
    interval: Interval,
    filter: ReadFilter,
    threads: usize,
) -> Result<u64> {
    if threads <= 1 || interval.is_empty() {
        let mut reader = IndexedAlignmentReader::open(bam_path)?;
        return count(&mut reader, contig, interval, filter);
    }
    let chunks =
        Scheduler::stats().chunk_interval(contig, interval.start.0, interval.stop.0, threads)?;
    let partial: Vec<u64> =
        parallel_map_regions(bam_path, chunks, threads, true, |_header, reader, chunk| {
            let mut n = 0u64;
            reader.for_each_raw(&chunk.contig, chunk.interval, |raw| {
                if count_owned_by_chunk(raw, interval, chunk.interval, filter)? {
                    n += 1;
                }
                Ok(())
            })?;
            Ok(n)
        })?;
    Ok(partial.into_iter().sum())
}

/// True when `raw` passes `filter` and its start-owned position falls in `chunk`.
fn count_owned_by_chunk(
    raw: &bam::Record,
    parent: Interval,
    chunk: Interval,
    filter: ReadFilter,
) -> Result<bool> {
    if !filter.passes_raw(raw) {
        return Ok(false);
    }
    let start = raw_alignment_start_0based(raw)?;
    if start < 0 {
        return Ok(false);
    }
    let start = start as u64;
    let owner = if start < parent.start.0 {
        parent.start.0
    } else {
        start
    };
    Ok(owner >= chunk.start.0 && owner < chunk.stop.0)
}

/// Per-position depth accumulator for `[start, stop)`.
#[derive(Debug, Clone, Default)]
pub struct DepthProfile {
    pub start: u64,
    pub depth: Vec<u32>,
}

impl DepthProfile {
    pub fn new(start: u64, len: usize) -> Self {
        Self {
            start,
            depth: vec![0; len],
        }
    }

    pub fn len(&self) -> usize {
        self.depth.len()
    }

    pub fn is_empty(&self) -> bool {
        self.depth.is_empty()
    }

    /// Merge another profile for the same interval (parallel reduce).
    pub fn merge(&mut self, other: &Self) {
        assert_eq!(self.start, other.start);
        assert_eq!(self.depth.len(), other.depth.len());
        for (a, b) in self.depth.iter_mut().zip(&other.depth) {
            *a += b;
        }
    }
}

/// Per-position base counts (A/C/G/T) for coverage.
#[derive(Debug, Clone, Default)]
pub struct CoverageProfile {
    pub start: u64,
    pub a: Vec<u32>,
    pub c: Vec<u32>,
    pub g: Vec<u32>,
    pub t: Vec<u32>,
}

impl CoverageProfile {
    pub fn new(start: u64, len: usize) -> Self {
        Self {
            start,
            a: vec![0; len],
            c: vec![0; len],
            g: vec![0; len],
            t: vec![0; len],
        }
    }

    pub fn len(&self) -> usize {
        self.a.len()
    }

    pub fn is_empty(&self) -> bool {
        self.a.is_empty()
    }

    pub fn merge(&mut self, other: &Self) {
        assert_eq!(self.start, other.start);
        assert_eq!(self.a.len(), other.a.len());
        for i in 0..self.a.len() {
            self.a[i] += other.a[i];
            self.c[i] += other.c[i];
            self.g[i] += other.g[i];
            self.t[i] += other.t[i];
        }
    }
}

/// Add base counts for coverage with optional quality threshold.
pub fn add_record_coverage(
    raw: &bam::Record,
    interval: Interval,
    profile: &mut CoverageProfile,
    quality_threshold: u8,
) -> Result<()> {
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
                        if idx >= profile.len() {
                            break;
                        }
                        let qidx = query_pos + skip + j;
                        if qidx < quals.len() && quals[qidx] >= quality_threshold {
                            if let Some(base) = seq.get(qidx) {
                                increment_base(profile, idx, base);
                            }
                        }
                    }
                }
                ref_pos += len as u64;
                query_pos += len;
            }
            Kind::Insertion | Kind::SoftClip => query_pos += len,
            Kind::Deletion | Kind::Skip => ref_pos += len as u64,
            _ => {}
        }
    }
    Ok(())
}

/// Add aligned bases (CIGAR M/=/X) to depth for positions in `interval`.
///
/// Matches **samtools depth** / rubam `get_depths`: deletions (`D`) and
/// ref-skips (`N`) do not contribute. Ambiguous query bases still increment
/// depth (unlike `count_coverage`, which only tallies A/C/G/T).
/// Filter matches pysam `read_callback="all"` (exclude unmapped / secondary /
/// qcfail / duplicate; **keep** supplementary). Sequence is not decoded.
pub fn add_record_depth(
    raw: &bam::Record,
    interval: Interval,
    profile: &mut DepthProfile,
) -> Result<()> {
    if !ReadFilter::All.passes_raw(raw) {
        return Ok(());
    }
    let start = raw_alignment_start_0based(raw)?;
    if start < 0 {
        return Ok(());
    }

    let mut ref_pos = start as u64;
    for item in raw.cigar().iter() {
        if ref_pos >= interval.stop.0 {
            break;
        }
        let op = item.map_err(SamRustError::from)?;
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if let Some((lo, hi)) = interval.overlap_span(ref_pos, len) {
                    let i0 = (lo - interval.start.0) as usize;
                    let i1 = ((hi - interval.start.0) as usize).min(profile.depth.len());
                    for d in &mut profile.depth[i0..i1] {
                        *d += 1;
                    }
                }
                ref_pos += len as u64;
            }
            Kind::Deletion | Kind::Skip => ref_pos += len as u64,
            _ => {}
        }
    }
    Ok(())
}

/// 0=A, 1=C, 2=G, 3=T, 4=other. Avoids `to_ascii_uppercase` in the inner loop.
const BASE_BUCKET: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
};

fn increment_base(cov: &mut CoverageProfile, idx: usize, base: u8) {
    match BASE_BUCKET[base as usize] {
        0 => cov.a[idx] += 1,
        1 => cov.c[idx] += 1,
        2 => cov.g[idx] += 1,
        3 => cov.t[idx] += 1,
        _ => {}
    }
}

/// Compute depth profile for a region (serial).
pub fn depth_profile(
    reader: &mut IndexedAlignmentReader,
    contig: &str,
    interval: Interval,
) -> Result<DepthProfile> {
    let len = interval.len() as usize;
    let mut profile = DepthProfile::new(interval.start.0, len);
    reader.for_each_raw(contig, interval, |raw| {
        add_record_depth(raw, interval, &mut profile)
    })?;
    Ok(profile)
}

/// Compute coverage (A/C/G/T) for a region (serial).
pub fn coverage_profile(
    reader: &mut IndexedAlignmentReader,
    contig: &str,
    interval: Interval,
    quality_threshold: u8,
) -> Result<CoverageProfile> {
    coverage_profile_with_filter(reader, contig, interval, quality_threshold, ReadFilter::All)
}

/// Coverage with explicit read filter.
pub fn coverage_profile_with_filter(
    reader: &mut IndexedAlignmentReader,
    contig: &str,
    interval: Interval,
    quality_threshold: u8,
    filter: ReadFilter,
) -> Result<CoverageProfile> {
    let len = interval.len() as usize;
    let mut profile = CoverageProfile::new(interval.start.0, len);
    reader.for_each_raw(contig, interval, |raw| {
        if filter.passes_raw(raw) {
            add_record_coverage(raw, interval, &mut profile, quality_threshold)?;
        }
        Ok(())
    })?;
    Ok(profile)
}

/// Run-length encoded depth blocks `(start, length, depth)`.
pub fn depth_blocks(profile: &DepthProfile) -> Vec<(u64, u64, u32)> {
    if profile.depth.is_empty() {
        return Vec::new();
    }
    let mut blocks = Vec::new();
    let mut run_start = profile.start;
    let mut run_depth = profile.depth[0];
    let mut run_len = 1u64;
    for (i, &d) in profile.depth.iter().enumerate().skip(1) {
        if d == run_depth {
            run_len += 1;
        } else {
            blocks.push((run_start, run_len, run_depth));
            run_start = profile.start + i as u64;
            run_depth = d;
            run_len = 1;
        }
    }
    blocks.push((run_start, run_len, run_depth));
    blocks
}

/// Parallel depth over chunked intervals; merges partial profiles by absolute position.
pub fn parallel_depth_profile(
    bam_path: &Path,
    contig: &str,
    interval: Interval,
    threads: usize,
) -> Result<DepthProfile> {
    let chunks =
        Scheduler::stats().chunk_interval(contig, interval.start.0, interval.stop.0, threads)?;
    let profiles: Vec<DepthProfile> =
        parallel_map_regions(bam_path, chunks, threads, true, |_header, reader, chunk| {
            depth_profile(reader, &chunk.contig, chunk.interval)
        })?;
    Ok(stitch_depth_profiles(interval, profiles))
}

/// Parallel coverage profile.
pub fn parallel_coverage_profile(
    bam_path: &Path,
    contig: &str,
    interval: Interval,
    quality_threshold: u8,
    threads: usize,
) -> Result<CoverageProfile> {
    parallel_coverage_profile_with_filter(
        bam_path,
        contig,
        interval,
        quality_threshold,
        ReadFilter::All,
        threads,
    )
}

pub fn parallel_coverage_profile_with_filter(
    bam_path: &Path,
    contig: &str,
    interval: Interval,
    quality_threshold: u8,
    filter: ReadFilter,
    threads: usize,
) -> Result<CoverageProfile> {
    let chunks =
        Scheduler::stats().chunk_interval(contig, interval.start.0, interval.stop.0, threads)?;
    let profiles: Vec<CoverageProfile> =
        parallel_map_regions(bam_path, chunks, threads, true, |_header, reader, chunk| {
            coverage_profile_with_filter(
                reader,
                &chunk.contig,
                chunk.interval,
                quality_threshold,
                filter,
            )
        })?;
    Ok(stitch_coverage_profiles(interval, profiles))
}

fn stitch_depth_profiles(parent: Interval, profiles: Vec<DepthProfile>) -> DepthProfile {
    let mut out = DepthProfile::new(parent.start.0, parent.len() as usize);
    for p in profiles {
        for (i, &d) in p.depth.iter().enumerate() {
            let abs = p.start + i as u64;
            if abs >= parent.start.0 && abs < parent.stop.0 {
                let idx = (abs - parent.start.0) as usize;
                out.depth[idx] += d;
            }
        }
    }
    out
}

fn stitch_coverage_profiles(parent: Interval, profiles: Vec<CoverageProfile>) -> CoverageProfile {
    let mut out = CoverageProfile::new(parent.start.0, parent.len() as usize);
    for p in profiles {
        for i in 0..p.len() {
            let abs = p.start + i as u64;
            if abs >= parent.start.0 && abs < parent.stop.0 {
                let idx = (abs - parent.start.0) as usize;
                out.a[idx] += p.a[i];
                out.c[idx] += p.c[i];
                out.g[idx] += p.g[i];
                out.t[idx] += p.t[i];
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
    fn count_nofilter_includes_placed_unmapped() {
        let bam = fixture_bam();
        let iv = Interval::new(0, 100).unwrap();
        let mut reader = IndexedAlignmentReader::open(&bam).unwrap();
        let nofilter = count(&mut reader, "chr1", iv, ReadFilter::NoFilter).unwrap();
        let mut reader = IndexedAlignmentReader::open(&bam).unwrap();
        let all = count(&mut reader, "chr1", iv, ReadFilter::All).unwrap();
        assert!(
            nofilter > all,
            "placed unmapped mate must be counted under nofilter only"
        );
        assert_eq!(
            nofilter,
            parallel_count(&bam, "chr1", iv, ReadFilter::NoFilter, 4).unwrap()
        );
    }

    #[test]
    fn count_nofilter_vs_all() {
        let bam = fixture_bam();
        let mut reader = IndexedAlignmentReader::open(&bam).unwrap();
        let iv = Interval::new(0, 1000).unwrap();
        let all = count(&mut reader, "chr1", iv, ReadFilter::NoFilter).unwrap();
        let mut reader = IndexedAlignmentReader::open(&bam).unwrap();
        let filtered = count(&mut reader, "chr1", iv, ReadFilter::All).unwrap();
        assert!(all >= filtered);
    }

    #[test]
    fn count_one_thread_matches_n_threads() {
        let bam = fixture_bam();
        let iv = Interval::new(0, 1000).unwrap();
        let serial = parallel_count(&bam, "chr1", iv, ReadFilter::NoFilter, 1).unwrap();
        for threads in [2, 4, 8] {
            let parallel = parallel_count(&bam, "chr1", iv, ReadFilter::NoFilter, threads).unwrap();
            assert_eq!(serial, parallel, "1T={serial} {threads}T={parallel}");
        }
        let mut reader = IndexedAlignmentReader::open(&bam).unwrap();
        assert_eq!(
            serial,
            count(&mut reader, "chr1", iv, ReadFilter::NoFilter).unwrap()
        );
    }

    #[test]
    fn depth_profile_nonzero_on_mapped_reads() {
        let bam = fixture_bam();
        let mut reader = IndexedAlignmentReader::open(&bam).unwrap();
        let iv = Interval::new(10, 60).unwrap();
        let profile = depth_profile(&mut reader, "chr1", iv).unwrap();
        assert!(profile.depth.iter().any(|&d| d > 0));
    }

    #[test]
    fn depth_skips_deletion_and_refskip_but_counts_mismatch_bases() {
        let bam = fixture_bam();
        let mut reader = IndexedAlignmentReader::open(&bam).unwrap();
        let iv = Interval::new(0, 1000).unwrap();
        let profile = depth_profile(&mut reader, "chr1", iv).unwrap();
        // indel1 10M2I8M3D10M @ 200 → D occupies 218..221
        assert_eq!(profile.depth[218], 0);
        assert_eq!(profile.depth[219], 0);
        assert_eq!(profile.depth[220], 0);
        assert!(profile.depth[217] > 0);
        assert!(profile.depth[221] > 0);
        // splice1 15M50N15M @ 400 → N occupies 415..465
        assert_eq!(profile.depth[415], 0);
        assert_eq!(profile.depth[464], 0);
        assert!(profile.depth[414] > 0);
        assert!(profile.depth[465] > 0);
        // eqx1 20=2X18= @ 500, query bases at the 2X are N — still depth
        assert!(profile.depth[520] > 0);
        assert!(profile.depth[521] > 0);
        let parallel = parallel_depth_profile(&bam, "chr1", iv, 4).unwrap();
        assert_eq!(profile.depth, parallel.depth);
    }

    #[test]
    fn coverage_one_thread_matches_n_threads() {
        let bam = fixture_bam();
        let iv = Interval::new(0, 1000).unwrap();
        let mut reader = IndexedAlignmentReader::open(&bam).unwrap();
        let serial =
            coverage_profile_with_filter(&mut reader, "chr1", iv, 0, ReadFilter::All).unwrap();
        for threads in [2, 4, 8] {
            let parallel = parallel_coverage_profile_with_filter(
                &bam,
                "chr1",
                iv,
                0,
                ReadFilter::All,
                threads,
            )
            .unwrap();
            assert_eq!(serial.a, parallel.a, "A 1T vs {threads}T");
            assert_eq!(serial.c, parallel.c, "C 1T vs {threads}T");
            assert_eq!(serial.g, parallel.g, "G 1T vs {threads}T");
            assert_eq!(serial.t, parallel.t, "T 1T vs {threads}T");
        }
    }
}
