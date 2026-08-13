//! Alignment record (pysam AlignedSegment field subset).

use noodles::bam;
use noodles::sam::alignment::record::Flags;

use crate::cigar::Cigar;
use crate::error::{Result, SamRustError};
use crate::header::Header;
use crate::tags::Tags;

/// Owned alignment record with pysam-oriented field semantics.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    query_name: String,
    flag: u16,
    reference_id: i32,
    reference_name: Option<String>,
    /// 0-based start; `-1` when unmapped (pysam).
    reference_start: i64,
    mapping_quality: u8,
    cigar: Cigar,
    query_sequence: String,
    /// Phred-scaled base qualities (same length as sequence when present).
    query_qualities: Vec<u8>,
    mate_reference_id: i32,
    mate_reference_start: i64,
    template_length: i32,
    tags: Tags,
}

impl Record {
    /// Convert a noodles BAM record into an owned SAMRust record.
    pub fn from_noodles(raw: &bam::Record, header: &Header) -> Result<Self> {
        let flags = raw.flags();
        let flag = u16::from(flags);

        let reference_id = match raw.reference_sequence_id() {
            Some(Ok(id)) => i32::try_from(id).unwrap_or(-1),
            Some(Err(e)) => return Err(SamRustError::from(e)),
            None => -1,
        };
        let reference_name = header.reference_name(reference_id).map(str::to_owned);

        let reference_start = match raw.alignment_start() {
            Some(Ok(pos)) => i64::try_from(pos.get()).unwrap_or(1) - 1, // 1-based → 0-based
            Some(Err(e)) => return Err(SamRustError::from(e)),
            None => -1,
        };

        let mapping_quality = raw.mapping_quality().map(u8::from).unwrap_or(255);

        let cigar = if flags.is_unmapped() {
            Cigar::empty()
        } else {
            Cigar::from_ops(raw.cigar().iter())?
        };

        let query_sequence: String = raw.sequence().iter().map(char::from).collect();
        let query_qualities: Vec<u8> = raw.quality_scores().iter().collect();

        let mate_reference_id = match raw.mate_reference_sequence_id() {
            Some(Ok(id)) => i32::try_from(id).unwrap_or(-1),
            Some(Err(e)) => return Err(SamRustError::from(e)),
            None => -1,
        };
        let mate_reference_start = match raw.mate_alignment_start() {
            Some(Ok(pos)) => i64::try_from(pos.get()).unwrap_or(1) - 1,
            Some(Err(e)) => return Err(SamRustError::from(e)),
            None => -1,
        };

        let query_name = raw
            .name()
            .map(|n| n.to_string())
            .unwrap_or_else(|| "*".to_string());

        let tags = Tags::from_noodles_data(raw.data().iter())?;

        Ok(Self {
            query_name,
            flag,
            reference_id,
            reference_name,
            reference_start,
            mapping_quality,
            cigar,
            query_sequence,
            query_qualities,
            mate_reference_id,
            mate_reference_start,
            template_length: raw.template_length(),
            tags,
        })
    }

    /// Convert any noodles SAM alignment record (BAM or CRAM `RecordBuf`).
    pub fn from_alignment(
        raw: &impl noodles::sam::alignment::Record,
        header: &Header,
        sam_header: &noodles::sam::Header,
    ) -> Result<Self> {
        let flags = raw.flags().map_err(SamRustError::from)?;
        let flag = u16::from(flags);

        let reference_id = match raw.reference_sequence_id(sam_header) {
            Some(Ok(id)) => i32::try_from(id).unwrap_or(-1),
            Some(Err(e)) => return Err(SamRustError::from(e)),
            None => -1,
        };
        let reference_name = header.reference_name(reference_id).map(str::to_owned);

        let reference_start = match raw.alignment_start() {
            Some(Ok(pos)) => i64::try_from(pos.get()).unwrap_or(1) - 1,
            Some(Err(e)) => return Err(SamRustError::from(e)),
            None => -1,
        };

        let mapping_quality = match raw.mapping_quality() {
            Some(Ok(mq)) => u8::from(mq),
            Some(Err(e)) => return Err(SamRustError::from(e)),
            None => 255,
        };

        let cigar = if flags.is_unmapped() {
            Cigar::empty()
        } else {
            Cigar::from_ops(raw.cigar().iter())?
        };

        let query_sequence: String = raw.sequence().iter().map(char::from).collect();
        let mut query_qualities = Vec::new();
        for q in raw.quality_scores().iter() {
            query_qualities.push(q.map_err(SamRustError::from)?);
        }

        let mate_reference_id = match raw.mate_reference_sequence_id(sam_header) {
            Some(Ok(id)) => i32::try_from(id).unwrap_or(-1),
            Some(Err(e)) => return Err(SamRustError::from(e)),
            None => -1,
        };
        let mate_reference_start = match raw.mate_alignment_start() {
            Some(Ok(pos)) => i64::try_from(pos.get()).unwrap_or(1) - 1,
            Some(Err(e)) => return Err(SamRustError::from(e)),
            None => -1,
        };

        let query_name = raw
            .name()
            .map(|n| n.to_string())
            .unwrap_or_else(|| "*".to_string());

        let tags = Tags::from_noodles_data(raw.data().iter())?;

        Ok(Self {
            query_name,
            flag,
            reference_id,
            reference_name,
            reference_start,
            mapping_quality,
            cigar,
            query_sequence,
            query_qualities,
            mate_reference_id,
            mate_reference_start,
            template_length: raw.template_length().map_err(SamRustError::from)?,
            tags,
        })
    }

    pub fn query_name(&self) -> &str {
        &self.query_name
    }
    pub fn flag(&self) -> u16 {
        self.flag
    }
    pub fn reference_id(&self) -> i32 {
        self.reference_id
    }
    pub fn reference_name(&self) -> Option<&str> {
        self.reference_name.as_deref()
    }
    pub fn reference_start(&self) -> i64 {
        self.reference_start
    }
    /// One past the last aligned reference base (0-based), pysam
    /// `reference_end` semantics: `None` for unmapped reads or empty CIGAR.
    pub fn reference_end(&self) -> Option<i64> {
        if self.is_unmapped() || self.reference_start < 0 || self.cigar.is_empty() {
            return None;
        }
        Some(self.reference_start + self.cigar.reference_length() as i64)
    }
    pub fn mapping_quality(&self) -> u8 {
        self.mapping_quality
    }
    pub fn cigar(&self) -> &Cigar {
        &self.cigar
    }
    pub fn cigarstring(&self) -> Option<String> {
        self.cigar.cigarstring()
    }
    pub fn query_sequence(&self) -> &str {
        &self.query_sequence
    }
    pub fn query_length(&self) -> usize {
        self.query_sequence.len()
    }
    pub fn query_qualities(&self) -> &[u8] {
        &self.query_qualities
    }
    pub fn mate_reference_id(&self) -> i32 {
        self.mate_reference_id
    }
    pub fn mate_reference_start(&self) -> i64 {
        self.mate_reference_start
    }
    pub fn template_length(&self) -> i32 {
        self.template_length
    }
    pub fn tags(&self) -> &Tags {
        &self.tags
    }

    fn flags(&self) -> Flags {
        Flags::from(self.flag)
    }

    pub fn is_paired(&self) -> bool {
        self.flags().is_segmented()
    }
    pub fn is_proper_pair(&self) -> bool {
        self.flags().is_properly_segmented()
    }
    pub fn is_unmapped(&self) -> bool {
        self.flags().is_unmapped()
    }
    pub fn mate_is_unmapped(&self) -> bool {
        self.flags().is_mate_unmapped()
    }
    pub fn is_reverse(&self) -> bool {
        self.flags().is_reverse_complemented()
    }
    pub fn mate_is_reverse(&self) -> bool {
        self.flags().is_mate_reverse_complemented()
    }
    pub fn is_read1(&self) -> bool {
        self.flags().is_first_segment()
    }
    pub fn is_read2(&self) -> bool {
        self.flags().is_last_segment()
    }
    pub fn is_secondary(&self) -> bool {
        self.flags().is_secondary()
    }
    pub fn is_qcfail(&self) -> bool {
        self.flags().is_qc_fail()
    }
    pub fn is_duplicate(&self) -> bool {
        self.flags().is_duplicate()
    }
    pub fn is_supplementary(&self) -> bool {
        self.flags().is_supplementary()
    }
}
