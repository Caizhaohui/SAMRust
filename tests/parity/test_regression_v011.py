"""v0.1.1 regression tests (REVIEW.md P0–P2).

P0a iter_batches streaming/exactly-once, P0b VariantFile.fetch without contig
length, P0c closed-file header access, P0d re-iteration buffer, P1c
parallel_fetch ownership semantics, P2a reference_end / header dict,
P2b B-array tags, P2c coordinate validation / clamping.
"""

from __future__ import annotations

import array
from pathlib import Path

import pytest

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "small.bam"

pysam = pytest.importorskip("pysam")
samrust = pytest.importorskip("samrust")


@pytest.fixture(scope="module")
def fixture_names():
    with pysam.AlignmentFile(str(FIXTURE), "rb") as af:
        return [r.query_name for r in af.fetch(until_eof=True)]


# --- P0a: iter_batches -------------------------------------------------------


def test_iter_batches_1t_streams_all_records(fixture_names) -> None:
    with samrust.AlignmentFile(str(FIXTURE), "rb") as sr:
        batches = list(sr.iter_batches(batch_size=5, threads=1))
    assert [len(b) for b in batches] == [5, 5, 4]
    assert [r.query_name for b in batches for r in b] == fixture_names


def test_iter_batches_mt_matches_1t_including_unmapped_tail(fixture_names) -> None:
    # The fixture ends with unmap1 (FUNMAP, no POS): invisible to region
    # queries, so MT must pick it up from the unmapped tail scan.
    with samrust.AlignmentFile(str(FIXTURE), "rb") as sr:
        mt = [r.query_name for b in sr.iter_batches(batch_size=3, threads=4) for r in b]
    assert mt == fixture_names


def test_iter_batches_mt_respects_batch_size() -> None:
    with samrust.AlignmentFile(str(FIXTURE), "rb") as sr:
        sizes = [len(b) for b in sr.iter_batches(batch_size=4, threads=4)]
    assert sizes == [4, 4, 4, 2]


def test_iter_batches_1t_mt_various_batch_sizes(fixture_names) -> None:
    with samrust.AlignmentFile(str(FIXTURE), "rb") as sr:
        for bs in (1, 2, 7, 256):
            got_1t = [r.query_name for b in sr.iter_batches(batch_size=bs, threads=1) for r in b]
            got_mt = [r.query_name for b in sr.iter_batches(batch_size=bs, threads=8) for r in b]
            assert got_1t == fixture_names, f"1T bs={bs}"
            assert got_mt == fixture_names, f"MT bs={bs}"


# --- P0d: re-iteration must not drop prefetched records ----------------------


def test_reiteration_continues_from_logical_position(fixture_names) -> None:
    sr = samrust.AlignmentFile(str(FIXTURE), "rb")
    it = iter(sr)
    next(it)
    rest = [r.query_name for r in it]
    assert rest == fixture_names[1:]
    sr.close()


def test_reset_rewinds(fixture_names) -> None:
    sr = samrust.AlignmentFile(str(FIXTURE), "rb")
    it = iter(sr)
    next(it)
    sr.reset()
    assert [r.query_name for r in sr] == fixture_names
    sr.close()


# --- P0c: closed file raises ValueError, not PanicException ------------------


def test_closed_file_header_access_raises_value_error() -> None:
    sr = samrust.AlignmentFile(str(FIXTURE), "rb")
    sr.close()
    for attr in ("references", "lengths", "nreferences", "header"):
        with pytest.raises(ValueError):
            getattr(sr, attr)


# --- P1c: parallel_fetch semantics -------------------------------------------


def test_parallel_fetch_overlapping_regions_return_union() -> None:
    with pysam.AlignmentFile(str(FIXTURE), "rb") as af:
        expected = sorted(r.query_name for r in af.fetch("chr1", 0, 1000))
    with samrust.AlignmentFile(str(FIXTURE), "rb") as sr:
        got = sorted(
            r.query_name
            for r in sr.parallel_fetch([("chr1", 0, 1000), ("chr1", 0, 1000)], threads=4)
        )
    assert got == expected


def test_parallel_fetch_1t_matches_nt_ordered() -> None:
    with pysam.AlignmentFile(str(FIXTURE), "rb") as af:
        expected = [r.query_name for r in af.fetch("chr1", 0, 1000)]
    with samrust.AlignmentFile(str(FIXTURE), "rb") as sr:
        for threads in (1, 2, 8):
            got = [r.query_name for r in sr.parallel_fetch([("chr1", 0, 1000)], threads=threads)]
            assert got == expected, f"threads={threads}"


def test_parallel_fetch_preserves_exact_duplicate_records(tmp_path) -> None:
    # cat-bam scenario: two byte-identical records must both survive (the old
    # hash-based dedup dropped one).
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    unsorted = tmp_path / "dup.unsorted.bam"
    bam_path = tmp_path / "dup.bam"

    def mkrec(qname, start):
        a = pysam.AlignedSegment()
        a.query_name = qname
        a.query_sequence = "ACGT"
        a.flag = 0
        a.reference_id = 0
        a.reference_start = start
        a.mapping_quality = 20
        a.cigar = ((0, 4),)
        return a

    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        out.write(mkrec("dupA", 10))
        out.write(mkrec("dupA", 10))
        out.write(mkrec("dupB", 10))
        out.write(mkrec("solo", 100))
    pysam.sort("-o", str(bam_path), str(unsorted))
    pysam.index(str(bam_path))

    with samrust.AlignmentFile(str(bam_path), "rb") as sr:
        got = [r.query_name for r in sr.parallel_fetch([("chr1", 0, 1000)], threads=4)]
    assert got == ["dupA", "dupA", "dupB", "solo"]


# --- P2a: reference_end / header dict ----------------------------------------


def test_reference_end_matches_pysam() -> None:
    with pysam.AlignmentFile(str(FIXTURE), "rb") as py, samrust.AlignmentFile(
        str(FIXTURE), "rb"
    ) as sr:
        for a, b in zip(py.fetch(until_eof=True), sr.fetch("chr1", 0, 1000)):
            assert b.reference_end == a.reference_end, a.query_name
        # unmapped placed read -> None
        unmap = [r for r in sr.fetch("chr1", 0, 100) if r.is_unmapped]
        assert unmap and unmap[0].reference_end is None


def test_header_dict_matches_pysam() -> None:
    with pysam.AlignmentFile(str(FIXTURE), "rb") as py:
        expected = dict(py.header)
    with samrust.AlignmentFile(str(FIXTURE), "rb") as sr:
        got = sr.header
    assert got == expected


# --- P2b: B-array tags --------------------------------------------------------


def test_b_array_tags_match_pysam(tmp_path) -> None:
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    bam_path = tmp_path / "arr.bam"
    a = pysam.AlignedSegment()
    a.query_name = "arr1"
    a.query_sequence = "ACGT"
    a.flag = 0
    a.reference_id = 0
    a.reference_start = 10
    a.mapping_quality = 20
    a.cigar = ((0, 4),)
    tags = [
        ("Ta", array.array("b", [1, 2])),
        ("Tb", array.array("B", [1, 2])),
        ("Tc", array.array("h", [1, 2])),
        ("Td", array.array("H", [1, 2])),
        ("Te", array.array("i", [1, 2])),
        ("Tf", array.array("I", [1, 2])),
        ("Tg", array.array("f", [1.5, 2.5])),
    ]
    for name, val in tags:
        a.set_tag(name, val)
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        out.write(a)
    pysam.index(str(bam_path))

    with pysam.AlignmentFile(str(bam_path), "rb") as py:
        pr = next(py.fetch("chr1", 0, 100))
    with samrust.AlignmentFile(str(bam_path), "rb") as sr:
        srec = next(sr.fetch("chr1", 0, 100))
    for name, _ in tags:
        sv, pv = srec.get_tag(name), pr.get_tag(name)
        assert sv == pv, f"{name}: {sv!r} != {pv!r}"
        assert sv.typecode == pv.typecode


# --- P2c: coordinate validation / clamping ------------------------------------


def test_negative_coordinates_raise_value_error() -> None:
    with samrust.AlignmentFile(str(FIXTURE), "rb") as sr:
        with pytest.raises(ValueError):
            sr.count("chr1", -1, 100)
        with pytest.raises(ValueError):
            list(sr.fetch("chr1", -5, 100))
        with pytest.raises(ValueError):
            sr.depth_numpy("chr1", 0, -1)
        with pytest.raises(ValueError):
            sr.pileup_counts("chr1", -1, 100)


def test_stop_beyond_contig_is_clamped() -> None:
    with pysam.AlignmentFile(str(FIXTURE), "rb") as py:
        expected_count = py.count("chr1", 0, 5000)
        expected_fetch = len(list(py.fetch("chr1", 0, 5000)))
    with samrust.AlignmentFile(str(FIXTURE), "rb") as sr:
        assert sr.count("chr1", 0, 5000) == expected_count
        assert len(list(sr.fetch("chr1", 0, 5000))) == expected_fetch
        # start beyond contig end -> empty (pysam behavior)
        assert sr.count("chr1", 1500, 5000) == 0
        assert list(sr.fetch("chr1", 1500, 5000)) == []


def test_unknown_contig_raises_value_error() -> None:
    with samrust.AlignmentFile(str(FIXTURE), "rb") as sr:
        with pytest.raises(ValueError):
            sr.count("chrMissing", 0, 100)
