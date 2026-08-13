"""Typing stubs for SAMRust (v0.1.1 surface)."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from typing import Any

__version__: str

def version() -> str: ...

class AlignedSegment:
    query_name: str
    flag: int
    reference_id: int
    reference_name: str | None
    reference_start: int
    reference_end: int | None
    mapping_quality: int
    cigarstring: str | None
    cigartuples: list[tuple[int, int]]
    query_sequence: str
    query_length: int
    query_qualities: list[int]
    next_reference_id: int
    next_reference_start: int
    template_length: int
    is_paired: bool
    is_proper_pair: bool
    is_unmapped: bool
    mate_is_unmapped: bool
    is_reverse: bool
    mate_is_reverse: bool
    is_read1: bool
    is_read2: bool
    is_secondary: bool
    is_qcfail: bool
    is_duplicate: bool
    is_supplementary: bool
    def has_tag(self, tag: str) -> bool: ...
    def get_tag(self, tag: str) -> Any: ...

class BatchIterator:
    """Iterator of `list[AlignedSegment]` batches (AlignmentFile.iter_batches)."""

    def __iter__(self) -> BatchIterator: ...
    def __next__(self) -> list[AlignedSegment]: ...

class AlignmentFile:
    filename: str
    mode: str
    references: list[str]
    lengths: list[int]
    nreferences: int
    header: dict[str, Any]
    def __init__(
        self,
        filename: str,
        mode: str = "rb",
        reference_filename: str | None = None,
    ) -> None: ...
    def close(self) -> None: ...
    def reset(self) -> None: ...
    def __enter__(self) -> AlignmentFile: ...
    def __exit__(self, *args: object) -> bool: ...
    def __iter__(self) -> AlignmentFile: ...
    def __next__(self) -> AlignedSegment: ...
    def fetch(
        self, contig: str, start: int | None = None, stop: int | None = None
    ) -> Iterator[AlignedSegment]: ...
    def count(
        self,
        contig: str,
        start: int | None = None,
        stop: int | None = None,
        read_callback: str = "nofilter",
        threads: int = 1,
    ) -> int: ...
    def count_coverage(
        self,
        contig: str,
        start: int | None = None,
        stop: int | None = None,
        quality_threshold: int = 15,
        read_callback: str = "all",
        threads: int = 1,
    ) -> tuple[Any, Any, Any, Any]: ...
    def depth_blocks(
        self,
        contig: str,
        start: int | None = None,
        stop: int | None = None,
        threads: int = 1,
    ) -> list[tuple[int, int, int]]: ...
    def depth_numpy(
        self,
        contig: str,
        start: int | None = None,
        stop: int | None = None,
        threads: int = 1,
    ) -> Any: ...
    def pileup_counts(
        self,
        contig: str,
        start: int | None = None,
        stop: int | None = None,
        min_base_quality: int = 0,
        min_mapping_quality: int = 0,
        threads: int = 1,
    ) -> dict[str, Any]: ...
    def iter_batches(
        self, batch_size: int = 256, threads: int = 1, ordered: bool = True
    ) -> BatchIterator: ...
    def parallel_fetch(
        self,
        regions: Iterable[Sequence[object]],
        threads: int = 1,
        ordered: bool = True,
    ) -> list[AlignedSegment]: ...

class VariantHeader:
    samples: list[str]
    contigs: list[str]

class _VariantSample:
    def __getitem__(self, key: str) -> Any: ...  # "GT" / "DP" / "AD"

class _VariantRecordSamples:
    def __len__(self) -> int: ...
    def __getitem__(self, key: int | str) -> _VariantSample: ...

class VariantRecord:
    chrom: str
    contig: str
    pos: int  # 1-based, like pysam
    start: int  # 0-based
    stop: int  # 0-based exclusive
    id: str | None
    ref: str
    alts: Any
    alleles: Any
    qual: float | None
    filter: list[str]
    format: list[str]
    info: Any
    samples: _VariantRecordSamples

class VariantFile:
    filename: str
    header: VariantHeader
    samples: list[str]
    def __init__(self, filename: str, mode: str = "r") -> None: ...
    def close(self) -> None: ...
    def __enter__(self) -> VariantFile: ...
    def __exit__(self, *args: object) -> bool: ...
    def __iter__(self) -> VariantFile: ...
    def __next__(self) -> VariantRecord: ...
    def fetch(
        self,
        contig: str | None = None,
        start: int | None = None,
        stop: int | None = None,
    ) -> Iterator[VariantRecord]: ...
