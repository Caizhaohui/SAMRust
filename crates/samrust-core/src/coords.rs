//! Coordinate system helpers.
//!
//! Python-facing APIs are always **0-based, half-open** `[start, stop)`.
//! All `+1` / `-1` conversions must live in this module — nowhere else.

use crate::error::{Result, SamRustError};

/// 0-based reference position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Position(pub u64);

/// Half-open interval `[start, stop)` on a single contig.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Interval {
    /// Inclusive start (0-based).
    pub start: Position,
    /// Exclusive stop (0-based).
    pub stop: Position,
}

impl Interval {
    /// Create a validated half-open interval.
    pub fn new(start: u64, stop: u64) -> Result<Self> {
        if start > stop {
            return Err(SamRustError::InvalidArgument(format!(
                "interval start ({start}) must be <= stop ({stop})"
            )));
        }
        Ok(Self {
            start: Position(start),
            stop: Position(stop),
        })
    }

    /// Length in bases (`stop - start`).
    pub fn len(&self) -> u64 {
        self.stop.0 - self.start.0
    }

    /// Overlap of `[ref_pos, ref_pos + op_len)` with this 0-based half-open interval.
    ///
    /// Returns `[lo, hi)` in the same coordinate space, or `None` if disjoint.
    /// This is interval arithmetic, not a 0-based/1-based conversion.
    #[inline]
    pub fn overlap_span(&self, ref_pos: u64, op_len: usize) -> Option<(u64, u64)> {
        let op_end = ref_pos.saturating_add(op_len as u64);
        let lo = ref_pos.max(self.start.0);
        let hi = op_end.min(self.stop.0);
        (lo < hi).then_some((lo, hi))
    }

    /// Whether the interval is empty (`start == stop`).
    pub fn is_empty(&self) -> bool {
        self.start == self.stop
    }

    /// Convert to 1-based inclusive coordinates for tools that require them.
    ///
    /// Empty intervals map to `(start+1, start)` which callers must treat carefully.
    pub fn to_1based_inclusive(&self) -> (u64, u64) {
        if self.is_empty() {
            (self.start.0.saturating_add(1), self.start.0)
        } else {
            (self.start.0 + 1, self.stop.0)
        }
    }

    /// Build a noodles [`Region`] for indexed BAM query (1-based inclusive).
    ///
    /// Returns `None` for empty intervals (caller should yield no records).
    pub fn to_noodles_region(&self, contig: &str) -> Result<Option<noodles::core::Region>> {
        if self.is_empty() {
            return Ok(None);
        }
        let (start, end) = self.to_1based_inclusive();
        let start_pos = noodles::core::Position::try_from(start as usize).map_err(|_| {
            SamRustError::InvalidArgument(format!("region start out of range: {start}"))
        })?;
        let end_pos = noodles::core::Position::try_from(end as usize).map_err(|_| {
            SamRustError::InvalidArgument(format!("region end out of range: {end}"))
        })?;
        Ok(Some(noodles::core::Region::new(
            contig,
            start_pos..=end_pos,
        )))
    }

    /// Build from 1-based inclusive coordinates.
    pub fn from_1based_inclusive(start: u64, end: u64) -> Result<Self> {
        if start == 0 {
            return Err(SamRustError::InvalidArgument(
                "1-based start must be >= 1".into(),
            ));
        }
        if end < start {
            return Err(SamRustError::InvalidArgument(format!(
                "1-based inclusive end ({end}) must be >= start ({start})"
            )));
        }
        Self::new(start - 1, end)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_half_open() {
        let iv = Interval::new(100, 200).unwrap();
        assert_eq!(iv.len(), 100);
        let (s, e) = iv.to_1based_inclusive();
        assert_eq!((s, e), (101, 200));
        let back = Interval::from_1based_inclusive(s, e).unwrap();
        assert_eq!(back, iv);
    }

    #[test]
    fn rejects_inverted_interval() {
        assert!(Interval::new(10, 5).is_err());
    }

    #[test]
    fn empty_interval_ok() {
        let iv = Interval::new(42, 42).unwrap();
        assert!(iv.is_empty());
        assert_eq!(iv.len(), 0);
    }

    #[test]
    fn overlap_span_clips_to_interval() {
        let iv = Interval::new(10, 20).unwrap();
        assert_eq!(iv.overlap_span(0, 5), None);
        assert_eq!(iv.overlap_span(20, 5), None);
        assert_eq!(iv.overlap_span(8, 4), Some((10, 12)));
        assert_eq!(iv.overlap_span(15, 10), Some((15, 20)));
        assert_eq!(iv.overlap_span(10, 10), Some((10, 20)));
    }
}

#[cfg(test)]
mod property_tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn rejects_inverted_start_stop(start in 0u64..1_000_000, delta in 1u64..10_000) {
            prop_assert!(Interval::new(start + delta, start).is_err());
        }

        #[test]
        fn empty_has_no_noodles_region(start in 0u64..1_000_000) {
            let iv = Interval::new(start, start).unwrap();
            prop_assert!(iv.is_empty());
            prop_assert_eq!(iv.len(), 0);
            prop_assert!(iv.to_noodles_region("chr1").unwrap().is_none());
        }

        #[test]
        fn overlap_span_is_subset_of_interval(
            start in 0u64..500,
            len in 1u64..200,
            pos in 0u64..700,
            op in 1usize..80,
        ) {
            let iv = Interval::new(start, start + len).unwrap();
            if let Some((lo, hi)) = iv.overlap_span(pos, op) {
                prop_assert!(lo >= iv.start.0);
                prop_assert!(hi <= iv.stop.0);
                prop_assert!(lo < hi);
            }
        }

        #[test]
        fn half_open_roundtrip(start in 0u64..50_000, len in 1u64..5_000) {
            let iv = Interval::new(start, start + len).unwrap();
            let (s, e) = iv.to_1based_inclusive();
            let back = Interval::from_1based_inclusive(s, e).unwrap();
            prop_assert_eq!(back, iv);
        }
    }
}
