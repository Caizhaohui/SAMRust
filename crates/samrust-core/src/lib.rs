//! SAMRust core library.
//!
//! M2: linear BAM access. M4+: indexed fetch, parallel runtime, depth, pileup.
//! M9: VCF/BCF VariantFile read path. M11: CRAM sequential + indexed fetch.

#![deny(unsafe_code)]

pub mod bam;
pub(crate) mod base;
pub mod cigar;
pub mod coords;
pub mod cram;
pub mod depth;
pub mod error;
pub mod header;
pub mod indexed;
pub mod parallel;
pub mod pileup;
pub mod record;
pub mod recount;
pub mod tags;
pub mod vcf;

pub use bam::{AlignmentReader, RecordIter};
pub use cigar::Cigar;
pub use coords::{Interval, Position};
pub use cram::{is_cram_path, CramAlignmentReader};
pub use depth::{
    count, coverage_profile, coverage_profile_with_filter, depth_blocks, depth_profile,
    parallel_count, parallel_coverage_profile, parallel_coverage_profile_with_filter,
    parallel_depth_profile, CoverageProfile, DepthProfile, ReadFilter,
};
pub use error::{Result, SamRustError};
pub use header::{Header, HeaderDict, SqEntry};
pub use indexed::{FetchIter, IndexedAlignmentReader};
pub use parallel::{
    ordered_merge, parallel_fetch_regions, parallel_fetch_wave, parallel_map_regions,
    whole_file_windows, FetchWindow, RegionChunk, Scheduler,
};
pub use pileup::{parallel_pileup_counts, pileup_counts, PileupCounts, PileupFilter};
pub use record::Record;
pub use recount::{filter_alt_ge, load_sites, recount_bam, rows_to_tsv, RecountRow, RecountSite};
pub use tags::{TagValue, Tags};
pub use vcf::{VariantHeader, VariantInfoValue, VariantReader, VariantRecord, VariantSample};

/// Library version string (kept in sync with Cargo package version).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
