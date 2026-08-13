//! CIGAR helpers (pysam-compatible cigarstring).

use std::fmt;

use noodles::sam::alignment::record::cigar::op::{Kind, Op};

use crate::error::{Result, SamRustError};

/// Owned CIGAR as a list of (op char, length) pairs.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Cigar {
    ops: Vec<(char, u32)>,
}

impl Cigar {
    /// Empty CIGAR (unmapped / unavailable).
    pub fn empty() -> Self {
        Self { ops: Vec::new() }
    }

    /// Whether there are no operations.
    pub fn is_empty(&self) -> bool {
        self.ops.is_empty()
    }

    /// Number of CIGAR operations.
    pub fn len(&self) -> usize {
        self.ops.len()
    }

    /// Iterate `(op, length)` pairs.
    pub fn iter(&self) -> impl Iterator<Item = (char, u32)> + '_ {
        self.ops.iter().copied()
    }

    /// pysam-style cigartuples: `(op_code, length)` with op 0=M, 1=I, …
    pub fn cigartuples(&self) -> Vec<(u8, u32)> {
        self.ops
            .iter()
            .map(|&(kind, len)| (kind_to_op(kind), len))
            .collect()
    }

    /// Reference span consumed on reference (excludes clipping/insertions).
    pub fn reference_length(&self) -> u64 {
        self.ops
            .iter()
            .filter(|&&(kind, _)| matches!(kind, 'M' | 'D' | 'N' | '=' | 'X'))
            .map(|&(_, len)| u64::from(len))
            .sum()
    }

    /// pysam-style cigarstring (`"10M2I8M"`), or `None` when empty.
    pub fn cigarstring(&self) -> Option<String> {
        if self.ops.is_empty() {
            None
        } else {
            Some(self.to_string())
        }
    }

    /// Build from noodles CIGAR ops.
    pub fn from_ops<I>(ops: I) -> Result<Self>
    where
        I: IntoIterator<Item = std::io::Result<Op>>,
    {
        let mut out = Vec::new();
        for item in ops {
            let op = item.map_err(SamRustError::from)?;
            out.push((
                kind_char(op.kind()),
                u32::try_from(op.len()).unwrap_or(u32::MAX),
            ));
        }
        Ok(Self { ops: out })
    }
}

impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (kind, len) in &self.ops {
            write!(f, "{len}{kind}")?;
        }
        Ok(())
    }
}

fn kind_to_op(kind: char) -> u8 {
    match kind {
        'M' => 0,
        'I' => 1,
        'D' => 2,
        'N' => 3,
        'S' => 4,
        'H' => 5,
        'P' => 6,
        '=' => 7,
        'X' => 8,
        _ => 0,
    }
}

fn kind_char(kind: Kind) -> char {
    match kind {
        Kind::Match => 'M',
        Kind::Insertion => 'I',
        Kind::Deletion => 'D',
        Kind::Skip => 'N',
        Kind::SoftClip => 'S',
        Kind::HardClip => 'H',
        Kind::Pad => 'P',
        Kind::SequenceMatch => '=',
        Kind::SequenceMismatch => 'X',
    }
}
