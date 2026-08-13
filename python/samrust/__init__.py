"""SAMRust: Rust-native, multi-threaded, pysam-compatible HTS processing."""

from __future__ import annotations

from samrust._samrust import (
    AlignedSegment,
    AlignmentFile,
    VariantFile,
    VariantRecord,
    __version__,
    version,
)

__all__ = [
    "AlignedSegment",
    "AlignmentFile",
    "VariantFile",
    "VariantRecord",
    "__version__",
    "version",
]
