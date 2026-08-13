"""Typing stubs for SAMRust (M3–M6 surface)."""

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

class AlignmentFile:
    filename: str
    mode: str
    references: list[str]
    lengths: list[int]
    nreferences: int
    header: dict[str, Any]
    def __init__(self, filename: str, mode: str = "rb") -> None: ...
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
    def iter_batches(
        self, batch_size: int = 256, threads: int = 1, ordered: bool = True
    ) -> list[list[AlignedSegment]]: ...
    def parallel_fetch(
        self,
        regions: Iterable[Sequence[object]],
        threads: int = 1,
        ordered: bool = True,
    ) -> list[AlignedSegment]: ...
