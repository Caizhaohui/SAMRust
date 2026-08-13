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
