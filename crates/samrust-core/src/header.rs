//! BAM/SAM header view (reference dictionary).

use noodles::sam;

use crate::error::{Result, SamRustError};

/// Parsed BAM/SAM header focused on the reference dictionary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Header {
    references: Vec<String>,
    lengths: Vec<u64>,
}

impl Header {
    /// Build from a noodles SAM header.
    pub fn from_noodles(header: &sam::Header) -> Self {
        let mut references = Vec::new();
        let mut lengths = Vec::new();
        for (name, map) in header.reference_sequences() {
            references.push(name.to_string());
            lengths.push(u64::try_from(map.length().get()).unwrap_or(u64::MAX));
        }
        Self {
            references,
            lengths,
        }
    }

    /// Number of reference sequences.
    pub fn nreferences(&self) -> usize {
        self.references.len()
    }

    /// Reference names in dictionary order.
    pub fn references(&self) -> &[String] {
        &self.references
    }

    /// Reference lengths in dictionary order.
    pub fn lengths(&self) -> &[u64] {
        &self.lengths
    }

    /// Resolve a reference id to a name.
    pub fn reference_name(&self, id: i32) -> Option<&str> {
        if id < 0 {
            return None;
        }
        self.references.get(id as usize).map(String::as_str)
    }

    /// Resolve a reference name to an id.
    pub fn reference_id(&self, name: &str) -> Result<i32> {
        self.references
            .iter()
            .position(|n| n == name)
            .map(|i| i as i32)
            .ok_or_else(|| SamRustError::InvalidArgument(format!("unknown reference: {name}")))
    }
}

/// One SQ record: `SN`, `LN`, plus any nonstandard fields in file order.
#[derive(Debug, Clone, PartialEq)]
pub struct SqEntry {
    pub sn: String,
    pub ln: u64,
    pub extra: Vec<(String, String)>,
}

/// pysam-compatible header dictionary content (`{'HD':.., 'SQ':[..], 'RG':[..], 'PG':[..]}`).
///
/// Field values are strings except `LN` (u64). Nonstandard fields are carried
/// through verbatim from the SAM header records.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct HeaderDict {
    /// HD record fields in file order (VN first when present).
    pub hd: Vec<(String, String)>,
    /// SQ records per reference, in file order.
    pub sq: Vec<SqEntry>,
    /// RG records: `(ID, fields)` in file order.
    pub rg: Vec<(String, Vec<(String, String)>)>,
    /// PG records: `(ID, fields)` in file order.
    pub pg: Vec<(String, Vec<(String, String)>)>,
}

fn bstr_to_string(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

impl HeaderDict {
    /// Extract a pysam-style view from a noodles SAM header.
    pub fn from_noodles(header: &sam::Header) -> Self {
        let mut out = Self::default();

        if let Some(hd) = header.header() {
            let v = hd.version();
            out.hd
                .push(("VN".into(), format!("{}.{}", v.major(), v.minor())));
            for (tag, value) in hd.other_fields() {
                out.hd
                    .push((bstr_to_string(tag.as_ref()), bstr_to_string(value.as_ref())));
            }
        }

        for (name, map) in header.reference_sequences() {
            let extras = map
                .other_fields()
                .iter()
                .map(|(tag, value)| (bstr_to_string(tag.as_ref()), bstr_to_string(value.as_ref())))
                .collect();
            out.sq.push(SqEntry {
                sn: bstr_to_string(name.as_ref()),
                ln: map.length().get() as u64,
                extra: extras,
            });
        }

        for (id, map) in header.read_groups() {
            let fields = map
                .other_fields()
                .iter()
                .map(|(tag, value)| (bstr_to_string(tag.as_ref()), bstr_to_string(value.as_ref())))
                .collect();
            out.rg.push((bstr_to_string(id.as_ref()), fields));
        }

        for (id, map) in header.programs().as_ref() {
            let fields = map
                .other_fields()
                .iter()
                .map(|(tag, value)| (bstr_to_string(tag.as_ref()), bstr_to_string(value.as_ref())))
                .collect();
            out.pg.push((bstr_to_string(id.as_ref()), fields));
        }

        out
    }
}
