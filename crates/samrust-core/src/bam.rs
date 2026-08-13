//! Linear BAM reader (M2). Indexed fetch in [`crate::indexed`] (M4).

use std::fs::File;
use std::path::{Path, PathBuf};

use noodles::bam;
use noodles::bgzf;

use crate::error::{Result, SamRustError};
use crate::header::Header;
use crate::record::Record;

type RawReader = bam::io::Reader<bgzf::io::Reader<File>>;

/// Sequential BAM reader with cached header.
pub struct AlignmentReader {
    path: PathBuf,
    reader: RawReader,
    header: Header,
    raw_header: noodles::sam::Header,
}

impl AlignmentReader {
    /// Open a BAM file and read its header.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let file = File::open(&path).map_err(|e| {
            SamRustError::Io(std::io::Error::new(
                e.kind(),
                format!("failed to open {}: {e}", path.display()),
            ))
        })?;
        let mut reader = bam::io::Reader::new(file);
        let raw_header = reader.read_header().map_err(SamRustError::from)?;
        let header = Header::from_noodles(&raw_header);
        Ok(Self {
            path,
            reader,
            header,
            raw_header,
        })
    }

    /// Path used to open this reader.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Reference dictionary / header view.
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Underlying noodles SAM header (for advanced use / later milestones).
    pub fn raw_header(&self) -> &noodles::sam::Header {
        &self.raw_header
    }

    /// Reopen the BAM from the same path (pysam `reset`).
    pub fn reset(&mut self) -> Result<()> {
        *self = Self::open(&self.path)?;
        Ok(())
    }

    /// Decode up to `batch_size` records in one call (GIL-free batch path for M3).
    pub fn read_batch(&mut self, batch_size: usize) -> Result<Vec<Record>> {
        if batch_size == 0 {
            return Ok(Vec::new());
        }
        let mut batch = Vec::with_capacity(batch_size);
        let mut iter = self.reader.records();
        let header = &self.header;
        for _ in 0..batch_size {
            match iter.next() {
                Some(Ok(raw)) => batch.push(Record::from_noodles(&raw, header)?),
                Some(Err(e)) => return Err(SamRustError::from(e)),
                None => break,
            }
        }
        Ok(batch)
    }

    /// Iterate records in file order (including unmapped).
    pub fn records(&mut self) -> RecordIter<'_> {
        RecordIter {
            header: &self.header,
            inner: self.reader.records(),
        }
    }
}

/// Iterator yielding owned [`Record`] values.
pub struct RecordIter<'a> {
    header: &'a Header,
    inner: bam::io::reader::Records<'a, bgzf::io::Reader<File>>,
}

impl Iterator for RecordIter<'_> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.next()? {
            Ok(raw) => Some(Record::from_noodles(&raw, self.header)),
            Err(e) => Some(Err(SamRustError::from(e))),
        }
    }
}
