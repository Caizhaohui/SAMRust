//! CRAM sequential + indexed fetch (M11 evaluation).
//!
//! Stats (`count` / depth / coverage / `pileup_counts`) stay BAM-only: the hot
//! path walks `noodles::bam::Record`. CRAM decode uses `sam::alignment::RecordBuf`.

use std::fs::File;
use std::path::{Path, PathBuf};

use noodles::cram;
use noodles::fasta;
use noodles::fasta::repository::adapters::IndexedReader as FastaIndexedReader;

use crate::coords::Interval;
use crate::error::{Result, SamRustError};
use crate::header::Header;
use crate::record::Record;

type RawCramReader = cram::io::indexed_reader::IndexedReader<File>;

/// Indexed CRAM reader with an associated FASTA repository and `.crai`.
pub struct CramAlignmentReader {
    path: PathBuf,
    fasta: PathBuf,
    reader: RawCramReader,
    header: Header,
    raw_header: noodles::sam::Header,
    sequential: Vec<Record>,
    sequential_pos: usize,
    sequential_loaded: bool,
}

impl CramAlignmentReader {
    /// Open `path.cram` + `path.cram.crai` using `fasta` (must have `.fai`).
    pub fn open(path: impl AsRef<Path>, fasta: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let fasta = fasta.as_ref().to_path_buf();
        let repository = fasta_repository(&fasta)?;
        let mut reader = cram::io::indexed_reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(&path)
            .map_err(|e| SamRustError::from_index_io(&format!("{}.crai", path.display()), e))?;
        let raw_header = reader.read_header().map_err(SamRustError::from)?;
        let header = Header::from_noodles(&raw_header);
        Ok(Self {
            path,
            fasta,
            reader,
            header,
            raw_header,
            sequential: Vec::new(),
            sequential_pos: 0,
            sequential_loaded: false,
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn fasta_path(&self) -> &Path {
        &self.fasta
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Sequential records. M11 materializes the file once (evaluation / fixture-scale).
    pub fn read_batch(&mut self, batch_size: usize) -> Result<Vec<Record>> {
        if batch_size == 0 {
            return Ok(Vec::new());
        }
        self.ensure_sequential()?;
        let end = self
            .sequential_pos
            .saturating_add(batch_size)
            .min(self.sequential.len());
        let out = self.sequential[self.sequential_pos..end].to_vec();
        self.sequential_pos = end;
        Ok(out)
    }

    /// Collect remaining sequential records.
    pub fn records_all(&mut self) -> Result<Vec<Record>> {
        self.ensure_sequential()?;
        let out = self.sequential[self.sequential_pos..].to_vec();
        self.sequential_pos = self.sequential.len();
        Ok(out)
    }

    fn ensure_sequential(&mut self) -> Result<()> {
        if self.sequential_loaded {
            return Ok(());
        }
        let sam_header = self.raw_header.clone();
        let header = self.header.clone();
        let mut all = Vec::new();
        for result in self.reader.records(&sam_header) {
            let raw = result.map_err(SamRustError::from)?;
            all.push(Record::from_alignment(&raw, &header, &sam_header)?);
        }
        self.sequential = all;
        self.sequential_pos = 0;
        self.sequential_loaded = true;
        Ok(())
    }

    /// Indexed fetch on a 0-based half-open interval (requires `.crai`).
    ///
    /// Opens a fresh reader so sequential iteration state is not disturbed.
    pub fn fetch_records(&mut self, contig: &str, interval: Interval) -> Result<Vec<Record>> {
        self.header.reference_id(contig)?;
        if interval.is_empty() {
            return Ok(Vec::new());
        }
        let mut fresh = Self::open(&self.path, &self.fasta)?;
        let region = interval
            .to_noodles_region(contig)?
            .ok_or_else(|| SamRustError::InvalidArgument("empty fetch interval".into()))?;
        let sam_header = fresh.raw_header.clone();
        let header = fresh.header.clone();
        let query = fresh
            .reader
            .query(&sam_header, &region)
            .map_err(SamRustError::from)?;
        let mut out = Vec::new();
        for result in query.records() {
            let raw = result.map_err(SamRustError::from)?;
            out.push(Record::from_alignment(&raw, &header, &sam_header)?);
        }
        Ok(out)
    }
}

fn fasta_repository(fasta: &Path) -> Result<fasta::Repository> {
    let reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(fasta)
        .map_err(|e| {
            SamRustError::Io(std::io::Error::new(
                e.kind(),
                format!("failed to open FASTA {}: {e}", fasta.display()),
            ))
        })?;
    Ok(fasta::Repository::new(FastaIndexedReader::new(reader)))
}

/// True when `path` looks like a CRAM file.
pub fn is_cram_path(path: &Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("cram"))
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AlignmentReader;

    fn fixture_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures")
    }

    #[test]
    fn cram_sequential_matches_bam_qnames() {
        let dir = fixture_dir();
        let cram = dir.join("small.cram");
        let fasta = dir.join("small.fa");
        if !cram.is_file() {
            return;
        }
        let mut bam = AlignmentReader::open(dir.join("small.bam")).unwrap();
        let bam_names: Vec<_> = bam
            .records()
            .map(|r| r.unwrap().query_name().to_string())
            .collect();
        let mut reader = CramAlignmentReader::open(&cram, &fasta).unwrap();
        let cram_names: Vec<_> = reader
            .records_all()
            .unwrap()
            .into_iter()
            .map(|r| r.query_name().to_string())
            .collect();
        assert_eq!(cram_names, bam_names);
    }

    #[test]
    fn cram_empty_interval_fetch_is_empty() {
        let dir = fixture_dir();
        let cram = dir.join("small.cram");
        let fasta = dir.join("small.fa");
        if !cram.is_file() {
            return;
        }
        let mut reader = CramAlignmentReader::open(&cram, &fasta).unwrap();
        let recs = reader
            .fetch_records("chr1", Interval::new(100, 100).unwrap())
            .unwrap();
        assert!(recs.is_empty());
    }

    #[test]
    fn cram_fetch_subregion_nonempty() {
        let dir = fixture_dir();
        let cram = dir.join("small.cram");
        let fasta = dir.join("small.fa");
        if !cram.is_file() {
            return;
        }
        let mut reader = CramAlignmentReader::open(&cram, &fasta).unwrap();
        let recs = reader
            .fetch_records("chr1", Interval::new(10, 60).unwrap())
            .unwrap();
        assert!(!recs.is_empty());
    }
}
