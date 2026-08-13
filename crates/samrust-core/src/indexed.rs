//! Indexed BAM reader and region fetch (M4).

use std::fs::File;
use std::path::{Path, PathBuf};

use noodles::bam;
use noodles::bam::io::IndexedReader;
use noodles::bgzf;

use crate::coords::Interval;
use crate::error::{Result, SamRustError};
use crate::header::Header;
use crate::record::Record;

type RawIndexedReader = IndexedReader<bgzf::io::Reader<File>>;

/// BAM reader with an associated BAI/CSI index.
pub struct IndexedAlignmentReader {
    path: PathBuf,
    reader: RawIndexedReader,
    header: Header,
    raw_header: noodles::sam::Header,
}

impl IndexedAlignmentReader {
    /// Open a BAM and load `<path>.bai` or `<path>.csi`.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(&path)
            .map_err(|e| SamRustError::from_index_io(&path.display().to_string(), e))?;
        let mut reader = reader;
        let raw_header = reader.read_header().map_err(SamRustError::from)?;
        let header = Header::from_noodles(&raw_header);
        Ok(Self {
            path,
            reader,
            header,
            raw_header,
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn raw_header(&self) -> &noodles::sam::Header {
        &self.raw_header
    }

    /// Fetch records overlapping `contig` on the 0-based half-open interval.
    pub fn fetch(&mut self, contig: &str, interval: Interval) -> Result<FetchIter<'_>> {
        self.header.reference_id(contig)?;
        if interval.is_empty() {
            return Ok(FetchIter::empty(&self.header));
        }
        let region = interval
            .to_noodles_region(contig)?
            .ok_or_else(|| SamRustError::InvalidArgument("empty fetch interval".into()))?;
        let query = self
            .reader
            .query(&self.raw_header, &region)
            .map_err(SamRustError::from)?;
        Ok(FetchIter {
            query: Some(query),
            header: &self.header,
            scratch: bam::Record::default(),
        })
    }

    /// Collect all records for a fetch (Python `fetch` / tests).
    ///
    /// Stats paths should use [`Self::for_each_raw`] instead of materializing
    /// owned [`Record`] values (qname, sequence `String`, tags, CIGAR `Vec`).
    pub fn fetch_records(&mut self, contig: &str, interval: Interval) -> Result<Vec<Record>> {
        self.fetch(contig, interval)?.collect::<Result<Vec<_>>>()
    }

    /// Read the unmapped tail: records without a reference id (unmapped, no
    /// POS), which are not covered by any index bin and are therefore
    /// invisible to region queries.
    ///
    /// Seeks via the index pseudo-bin offset when present (O(tail)); otherwise
    /// falls back to a full scan. Placed-unmapped records (FUNMAP with POS)
    /// are indexed under their contig and are NOT returned here.
    pub fn unmapped_tail_records(&mut self) -> Result<Vec<Record>> {
        let header = self.header.clone();
        let iter = self.reader.query_unmapped().map_err(SamRustError::from)?;
        let mut out = Vec::new();
        for result in iter {
            let raw = result.map_err(SamRustError::from)?;
            if raw.reference_sequence_id().is_none() {
                out.push(Record::from_noodles(&raw, &header)?);
            }
        }
        Ok(out)
    }

    /// Visit each indexed-fetch hit as a noodles BAM record, reusing one scratch buffer.
    ///
    /// Does **not** build an owned [`Record`]. Python `fetch` still goes through
    /// [`Self::fetch`] / [`Self::fetch_records`].
    pub fn for_each_raw<F>(&mut self, contig: &str, interval: Interval, mut visit: F) -> Result<()>
    where
        F: FnMut(&bam::Record) -> Result<()>,
    {
        self.header.reference_id(contig)?;
        if interval.is_empty() {
            return Ok(());
        }
        let region = interval
            .to_noodles_region(contig)?
            .ok_or_else(|| SamRustError::InvalidArgument("empty fetch interval".into()))?;
        let mut query = self
            .reader
            .query(&self.raw_header, &region)
            .map_err(SamRustError::from)?;
        let mut scratch = bam::Record::default();
        loop {
            match query.read_record(&mut scratch) {
                Ok(0) => break,
                Ok(_) => visit(&scratch)?,
                Err(e) => return Err(SamRustError::from(e)),
            }
        }
        Ok(())
    }
}

/// 0-based alignment start, or `-1` when POS is absent (pysam unmapped).
pub(crate) fn raw_alignment_start_0based(raw: &bam::Record) -> Result<i64> {
    match raw.alignment_start() {
        Some(Ok(pos)) => Ok(i64::try_from(pos.get()).unwrap_or(1) - 1),
        Some(Err(e)) => Err(SamRustError::from(e)),
        None => Ok(-1),
    }
}

/// Iterator over records from an indexed region query.
pub struct FetchIter<'a> {
    query: Option<bam::io::reader::Query<'a, bgzf::io::Reader<File>>>,
    header: &'a Header,
    scratch: bam::Record,
}

impl FetchIter<'_> {
    fn empty(header: &Header) -> FetchIter<'_> {
        FetchIter {
            query: None,
            header,
            scratch: bam::Record::default(),
        }
    }
}

impl Iterator for FetchIter<'_> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let query = self.query.as_mut()?;
        match query.read_record(&mut self.scratch) {
            Ok(0) => {
                self.query = None;
                None
            }
            Ok(_) => Some(Record::from_noodles(&self.scratch, self.header)),
            Err(e) => Some(Err(SamRustError::from(e))),
        }
    }
}

/// Whether a mapped record overlaps a 0-based half-open interval on the reference.
///
/// Unmapped records return `false` even when POS is set. Do **not** use this for
/// pysam `count`/`fetch` parity — placed unmapped mates are valid fetch hits.
pub fn record_overlaps_interval(record: &Record, interval: Interval) -> bool {
    if record.is_unmapped() || record.reference_start() < 0 {
        return false;
    }
    let start = record.reference_start() as u64;
    let end = start.saturating_add(record.cigar().reference_length());
    start < interval.stop.0 && end > interval.start.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AlignmentReader;

    fn fixture_bam() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/small.bam")
    }

    #[test]
    fn fetch_full_contig_matches_linear_scan() {
        let bam = fixture_bam();
        let mut linear = AlignmentReader::open(&bam).unwrap();
        let linear_recs: Vec<_> = linear
            .records()
            .map(|r| r.unwrap())
            .filter(|r| r.reference_name() == Some("chr1"))
            .collect();

        let mut indexed = IndexedAlignmentReader::open(&bam).unwrap();
        let iv = Interval::new(0, 1000).unwrap();
        let fetched = indexed.fetch_records("chr1", iv).unwrap();
        assert_eq!(fetched.len(), linear_recs.len());
        assert_eq!(fetched, linear_recs);
    }

    #[test]
    fn fetch_subregion_and_empty_interval() {
        let bam = fixture_bam();
        let mut indexed = IndexedAlignmentReader::open(&bam).unwrap();
        let iv = Interval::new(10, 60).unwrap();
        let recs = indexed.fetch_records("chr1", iv).unwrap();
        assert!(!recs.is_empty());
        assert!(
            recs.iter().any(|r| r.query_name() == "unmap_placed"),
            "placed unmapped mate at POS=50 must be a fetch hit"
        );
        for r in &recs {
            if r.is_unmapped() {
                assert!(r.reference_start() >= iv.start.0 as i64);
                assert!(r.reference_start() < iv.stop.0 as i64);
            } else {
                assert!(record_overlaps_interval(r, iv));
            }
        }

        let empty = indexed.fetch_records("chr1", Interval::new(100, 100).unwrap());
        assert_eq!(empty.unwrap().len(), 0);
    }

    #[test]
    fn for_each_raw_matches_fetch_records_len() {
        let bam = fixture_bam();
        let mut indexed = IndexedAlignmentReader::open(&bam).unwrap();
        let iv = Interval::new(10, 60).unwrap();
        let owned = indexed.fetch_records("chr1", iv).unwrap();
        let mut n = 0u64;
        let mut saw_placed_unmapped = false;
        indexed
            .for_each_raw("chr1", iv, |raw| {
                n += 1;
                if raw.flags().is_unmapped() {
                    saw_placed_unmapped = true;
                }
                Ok(())
            })
            .unwrap();
        assert_eq!(n as usize, owned.len());
        assert!(saw_placed_unmapped);
    }

    #[test]
    fn invalid_contig_errors() {
        let bam = fixture_bam();
        let mut indexed = IndexedAlignmentReader::open(&bam).unwrap();
        assert!(indexed
            .fetch("missing", Interval::new(0, 10).unwrap())
            .is_err());
    }

    #[test]
    fn unmapped_tail_finds_unplaced_record() {
        // The Tier-0 fixture ends with `unmap1` (FUNMAP, no POS): invisible to
        // region queries, so it must come from the tail scan.
        let bam = fixture_bam();
        let mut indexed = IndexedAlignmentReader::open(&bam).unwrap();
        let tail = indexed.unmapped_tail_records().unwrap();
        assert_eq!(tail.len(), 1);
        assert_eq!(tail[0].query_name(), "unmap1");
        assert!(tail[0].is_unmapped());
        assert_eq!(tail[0].reference_id(), -1);
    }
}
