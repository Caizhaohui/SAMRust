//! VCF / BCF / VCF.gz reading (M9). Python fetch is 0-based half-open.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use noodles::bcf;
use noodles::vcf;
use noodles::vcf::variant::record::samples::keys::key as format_key;
use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
use noodles::vcf::variant::record_buf::samples::sample::Value as SampleValue;
use noodles::vcf::variant::Record as _;
use noodles::vcf::variant::RecordBuf;

use crate::coords::Interval;
use crate::error::{Result, SamRustError};

/// Header fields needed for the Python `VariantFile` surface.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariantHeader {
    pub samples: Vec<String>,
    pub contigs: Vec<String>,
    pub lengths: Vec<u64>,
}

impl VariantHeader {
    fn from_noodles(header: &vcf::Header) -> Self {
        let samples: Vec<String> = header.sample_names().iter().cloned().collect();
        let mut contigs = Vec::new();
        let mut lengths = Vec::new();
        for (name, contig) in header.contigs() {
            contigs.push(name.clone());
            lengths.push(contig.length().unwrap_or(0) as u64);
        }
        Self {
            samples,
            contigs,
            lengths,
        }
    }

    pub fn contig_index(&self, name: &str) -> Result<usize> {
        self.contigs
            .iter()
            .position(|c| c == name)
            .ok_or_else(|| SamRustError::InvalidArgument(format!("unknown contig: {name}")))
    }

    pub fn contig_length(&self, name: &str) -> Result<u64> {
        let i = self.contig_index(name)?;
        Ok(self.lengths[i])
    }
}

/// Owned INFO value (pysam-oriented).
#[derive(Debug, Clone, PartialEq)]
pub enum VariantInfoValue {
    Integer(i32),
    Float(f32),
    Flag,
    Character(char),
    String(String),
    IntegerArray(Vec<Option<i32>>),
    FloatArray(Vec<Option<f32>>),
    StringArray(Vec<Option<String>>),
}

/// One sample's GT / DP / AD (M9 priority FORMAT fields).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct VariantSample {
    /// Allele indices; `None` is a missing allele (`.`).
    pub gt: Option<Vec<Option<i32>>>,
    pub dp: Option<i32>,
    pub ad: Option<Vec<Option<i32>>>,
}

/// Owned variant record with pysam-oriented coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct VariantRecord {
    pub chrom: String,
    /// 0-based start.
    pub start: u64,
    /// 0-based exclusive stop.
    pub stop: u64,
    pub id: Option<String>,
    pub ref_allele: String,
    pub alts: Vec<String>,
    pub qual: Option<f32>,
    pub filter: Vec<String>,
    pub info: BTreeMap<String, VariantInfoValue>,
    pub format: Vec<String>,
    pub samples: Vec<VariantSample>,
}

impl VariantRecord {
    /// pysam `VariantRecord.pos`: 1-based.
    pub fn pos(&self) -> u64 {
        self.start.saturating_add(1)
    }
}

/// Reader for VCF, bgzipped VCF, or BCF.
pub struct VariantReader {
    path: PathBuf,
    kind: VariantKind,
    header: VariantHeader,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VariantKind {
    Vcf,
    Bcf,
}

impl VariantReader {
    /// Open a VCF / VCF.gz / BCF and read the header.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let kind = detect_kind(&path)?;
        let header = match kind {
            VariantKind::Vcf => {
                let mut reader = vcf::io::reader::Builder::default()
                    .build_from_path(&path)
                    .map_err(SamRustError::from)?;
                let raw = reader.read_header().map_err(SamRustError::from)?;
                VariantHeader::from_noodles(&raw)
            }
            VariantKind::Bcf => {
                let mut reader = bcf::io::reader::Builder::default()
                    .build_from_path(&path)
                    .map_err(SamRustError::from)?;
                let raw = reader.read_header().map_err(SamRustError::from)?;
                VariantHeader::from_noodles(&raw)
            }
        };
        Ok(Self { path, kind, header })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn header(&self) -> &VariantHeader {
        &self.header
    }

    /// Sequential records (file order).
    pub fn records(&self) -> Result<Vec<VariantRecord>> {
        match self.kind {
            VariantKind::Vcf => read_vcf_records(&self.path),
            VariantKind::Bcf => read_bcf_records(&self.path),
        }
    }

    /// Indexed region query; coordinates are 0-based half-open.
    ///
    /// Falls back to a sequential overlap filter when no index is present
    /// (plain VCF). Unknown contig is an error.
    pub fn fetch(&self, contig: &str, interval: Interval) -> Result<Vec<VariantRecord>> {
        self.header.contig_index(contig)?;
        if interval.is_empty() {
            return Ok(Vec::new());
        }
        match self.try_indexed_fetch(contig, interval) {
            Ok(recs) => Ok(recs),
            Err(SamRustError::MissingIndex(_)) => Ok(self
                .records()?
                .into_iter()
                .filter(|r| record_overlaps(r, contig, interval))
                .collect()),
            Err(e) => Err(e),
        }
    }

    /// Fetch `[start, ∞)` on a contig whose length is unknown from the header.
    ///
    /// An unbounded indexed query is not expressible for tabix/CSI (the binning
    /// scheme caps the end bound), so this always uses a sequential scan.
    pub fn fetch_from(&self, contig: &str, start: u64) -> Result<Vec<VariantRecord>> {
        self.header.contig_index(contig)?;
        Ok(self
            .records()?
            .into_iter()
            .filter(|r| r.chrom == contig && r.stop > start)
            .collect())
    }

    fn try_indexed_fetch(&self, contig: &str, interval: Interval) -> Result<Vec<VariantRecord>> {
        let region = interval
            .to_noodles_region(contig)?
            .ok_or_else(|| SamRustError::InvalidArgument("empty fetch interval".into()))?;
        match self.kind {
            VariantKind::Vcf => {
                let mut reader = vcf::io::indexed_reader::Builder::default()
                    .build_from_path(&self.path)
                    .map_err(|e| {
                        SamRustError::from_index_io(&self.path.display().to_string(), e)
                    })?;
                let header = reader.read_header().map_err(SamRustError::from)?;
                let query = reader.query(&header, &region).map_err(SamRustError::from)?;
                let mut out = Vec::new();
                for result in query.records() {
                    let rec = result.map_err(SamRustError::from)?;
                    let buf = RecordBuf::try_from_variant_record(&header, &rec)
                        .map_err(SamRustError::from)?;
                    let converted = variant_from_buf(&header, &buf)?;
                    if record_overlaps(&converted, contig, interval) {
                        out.push(converted);
                    }
                }
                Ok(out)
            }
            VariantKind::Bcf => {
                let mut reader = bcf::io::indexed_reader::Builder::default()
                    .build_from_path(&self.path)
                    .map_err(|e| {
                        SamRustError::from_index_io(&self.path.display().to_string(), e)
                    })?;
                let header = reader.read_header().map_err(SamRustError::from)?;
                let query = reader.query(&header, &region).map_err(SamRustError::from)?;
                let mut out = Vec::new();
                for result in query.records() {
                    let rec = result.map_err(SamRustError::from)?;
                    let buf = RecordBuf::try_from_variant_record(&header, &rec)
                        .map_err(SamRustError::from)?;
                    let converted = variant_from_buf(&header, &buf)?;
                    if record_overlaps(&converted, contig, interval) {
                        out.push(converted);
                    }
                }
                Ok(out)
            }
        }
    }
}

fn detect_kind(path: &Path) -> Result<VariantKind> {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    if name.ends_with(".bcf") {
        Ok(VariantKind::Bcf)
    } else if name.ends_with(".vcf") || name.ends_with(".vcf.gz") || name.ends_with(".vcf.bgz") {
        Ok(VariantKind::Vcf)
    } else {
        Err(SamRustError::InvalidArgument(format!(
            "unsupported variant file: {}",
            path.display()
        )))
    }
}

fn record_overlaps(record: &VariantRecord, contig: &str, interval: Interval) -> bool {
    record.chrom == contig && record.start < interval.stop.0 && record.stop > interval.start.0
}

fn read_vcf_records(path: &Path) -> Result<Vec<VariantRecord>> {
    let mut reader = vcf::io::reader::Builder::default()
        .build_from_path(path)
        .map_err(SamRustError::from)?;
    let header = reader.read_header().map_err(SamRustError::from)?;
    let mut out = Vec::new();
    for result in reader.record_bufs(&header) {
        let buf = result.map_err(SamRustError::from)?;
        out.push(variant_from_buf(&header, &buf)?);
    }
    Ok(out)
}

fn read_bcf_records(path: &Path) -> Result<Vec<VariantRecord>> {
    let mut reader = bcf::io::reader::Builder::default()
        .build_from_path(path)
        .map_err(SamRustError::from)?;
    let header = reader.read_header().map_err(SamRustError::from)?;
    let mut out = Vec::new();
    for result in reader.record_bufs(&header) {
        let buf = result.map_err(SamRustError::from)?;
        out.push(variant_from_buf(&header, &buf)?);
    }
    Ok(out)
}

fn variant_from_buf(header: &vcf::Header, rec: &RecordBuf) -> Result<VariantRecord> {
    let chrom = rec.reference_sequence_name().to_string();
    let start1 = rec.variant_start().map(|p| p.get() as u64).unwrap_or(1);
    let end1 = rec.variant_end(header).map_err(SamRustError::from)?.get() as u64;
    let interval = Interval::from_1based_inclusive(start1, end1)?;

    let ids = rec.ids().as_ref();
    let id = if ids.is_empty() {
        None
    } else {
        Some(ids.iter().cloned().collect::<Vec<_>>().join(";"))
    };

    let alts = rec.alternate_bases().as_ref().to_vec();
    let filter: Vec<String> = rec.filters().as_ref().iter().cloned().collect();
    let format: Vec<String> = rec.format().as_ref().iter().cloned().collect();

    let mut info = BTreeMap::new();
    for (key, value) in rec.info().as_ref() {
        if let Some(v) = value.as_ref() {
            info.insert(key.clone(), convert_info(v));
        }
    }

    let n_samples = header.sample_names().len();
    let mut samples = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let Some(sample) = rec.samples().get_index(i) else {
            samples.push(VariantSample::default());
            continue;
        };
        samples.push(convert_sample(&sample));
    }

    Ok(VariantRecord {
        chrom,
        start: interval.start.0,
        stop: interval.stop.0,
        id,
        ref_allele: rec.reference_bases().to_string(),
        alts,
        qual: rec.quality_score(),
        filter,
        info,
        format,
        samples,
    })
}

fn convert_info(value: &InfoValue) -> VariantInfoValue {
    match value {
        InfoValue::Integer(n) => VariantInfoValue::Integer(*n),
        InfoValue::Float(n) => VariantInfoValue::Float(*n),
        InfoValue::Flag => VariantInfoValue::Flag,
        InfoValue::Character(c) => VariantInfoValue::Character(*c),
        InfoValue::String(s) => VariantInfoValue::String(s.clone()),
        InfoValue::Array(arr) => match arr {
            noodles::vcf::variant::record_buf::info::field::value::Array::Integer(v) => {
                VariantInfoValue::IntegerArray(v.clone())
            }
            noodles::vcf::variant::record_buf::info::field::value::Array::Float(v) => {
                VariantInfoValue::FloatArray(v.clone())
            }
            noodles::vcf::variant::record_buf::info::field::value::Array::Character(v) => {
                VariantInfoValue::StringArray(
                    v.iter().map(|c| c.map(|ch| ch.to_string())).collect(),
                )
            }
            noodles::vcf::variant::record_buf::info::field::value::Array::String(v) => {
                VariantInfoValue::StringArray(v.clone())
            }
        },
    }
}

fn convert_sample(sample: &vcf::variant::record_buf::samples::Sample<'_>) -> VariantSample {
    let gt = sample
        .get(format_key::GENOTYPE)
        .flatten()
        .and_then(parse_gt);
    let dp = sample.get("DP").flatten().and_then(|v| match v {
        SampleValue::Integer(n) => Some(*n),
        _ => None,
    });
    let ad = sample.get("AD").flatten().and_then(|v| match v {
        SampleValue::Array(
            noodles::vcf::variant::record_buf::samples::sample::value::Array::Integer(v),
        ) => Some(v.clone()),
        SampleValue::Integer(n) => Some(vec![Some(*n)]),
        _ => None,
    });
    VariantSample { gt, dp, ad }
}

fn parse_gt(value: &SampleValue) -> Option<Vec<Option<i32>>> {
    match value {
        SampleValue::Genotype(gt) => Some(
            gt.as_ref()
                .iter()
                .map(|allele| {
                    allele
                        .position()
                        .map(|p| i32::try_from(p).unwrap_or(i32::MAX))
                })
                .collect(),
        ),
        SampleValue::String(s) => s
            .parse::<vcf::variant::record_buf::samples::sample::value::Genotype>()
            .ok()
            .map(|gt| {
                gt.as_ref()
                    .iter()
                    .map(|allele| {
                        allele
                            .position()
                            .map(|p| i32::try_from(p).unwrap_or(i32::MAX))
                    })
                    .collect()
            }),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_vcf() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/small.vcf.gz")
    }

    #[test]
    fn header_samples_and_contigs() {
        let reader = VariantReader::open(fixture_vcf()).unwrap();
        assert_eq!(reader.header().samples, ["sample1"]);
        assert_eq!(reader.header().contigs, ["chr1", "chr2"]);
    }

    #[test]
    fn sequential_len_and_snp_coords() {
        let recs = VariantReader::open(fixture_vcf())
            .unwrap()
            .records()
            .unwrap();
        assert_eq!(recs.len(), 5);
        assert_eq!(recs[0].chrom, "chr1");
        assert_eq!(recs[0].start, 14);
        assert_eq!(recs[0].stop, 15);
        assert_eq!(recs[0].pos(), 15);
        assert_eq!(recs[0].alts, ["G"]);
        assert_eq!(recs[2].start, 259);
        assert_eq!(recs[2].stop, 262);
    }

    #[test]
    fn fetch_half_open_matches_expected_ids() {
        let reader = VariantReader::open(fixture_vcf()).unwrap();
        let hit = reader
            .fetch("chr1", Interval::new(14, 15).unwrap())
            .unwrap();
        assert_eq!(hit.len(), 1);
        assert_eq!(hit[0].pos(), 15);
        let miss = reader.fetch("chr1", Interval::new(0, 14).unwrap()).unwrap();
        assert!(miss.is_empty());
        let empty = reader
            .fetch("chr1", Interval::new(100, 100).unwrap())
            .unwrap();
        assert!(empty.is_empty());
    }
}
