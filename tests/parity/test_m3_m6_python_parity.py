"""M3–M6 Python parity vs pysam on Tier-0 fixture."""

from __future__ import annotations

from pathlib import Path

import pytest

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "small.bam"


@pytest.fixture(scope="module")
def bam_pair():
    pysam = pytest.importorskip("pysam")
    samrust = pytest.importorskip("samrust")
    if not FIXTURE.is_file():
        pytest.skip("missing fixture")
    return (
        pysam.AlignmentFile(str(FIXTURE), "rb"),
        samrust.AlignmentFile(str(FIXTURE), "rb"),
    )


def test_m3_header_and_sequential_fields(bam_pair) -> None:
    py, sr = bam_pair
    assert list(py.references) == sr.references
    assert list(py.lengths) == sr.lengths
    assert py.nreferences == sr.nreferences

    for a, b in zip(py.fetch(until_eof=True), sr):
        assert a.query_name == b.query_name
        assert a.flag == b.flag
        assert a.reference_id == b.reference_id
        assert a.reference_start == b.reference_start
        assert a.mapping_quality == b.mapping_quality
        assert a.cigarstring == b.cigarstring
        assert a.query_length == b.query_length
        if a.has_tag("NM"):
            assert b.has_tag("NM")
            assert a.get_tag("NM") == b.get_tag("NM")


def test_m4_fetch_regions(bam_pair) -> None:
    py, sr = bam_pair
    regions = [
        ("chr1", 0, 100),
        ("chr1", 200, 250),
        ("chr1", 590, 650),
        ("chr1", 100, 100),  # empty
        ("chr2", 0, 100),
    ]
    for contig, start, stop in regions:
        py_ids = [
            (r.query_name, r.flag, r.reference_start, r.cigarstring)
            for r in py.fetch(contig, start, stop)
        ]
        sr_ids = [
            (r.query_name, r.flag, r.reference_start, r.cigarstring)
            for r in sr.fetch(contig, start, stop)
        ]
        assert sr_ids == py_ids, f"mismatch in {contig}:{start}-{stop}"


def test_m5_parallel_fetch_matches_serial(bam_pair) -> None:
    _, sr = bam_pair
    regions = [("chr1", 0, 700), ("chr2", 0, 100)]
    serial = []
    for contig, start, stop in regions:
        serial.extend(
            (r.query_name, r.flag, r.reference_start)
            for r in sr.fetch(contig, start, stop)
        )
    parallel = [
        (r.query_name, r.flag, r.reference_start)
        for r in sr.parallel_fetch(regions, threads=4, ordered=True)
    ]
    assert sorted(parallel) == sorted(serial)


def test_m5_thread_counts_identical_for_count(bam_pair) -> None:
    py, sr = bam_pair
    values = [sr.count("chr1", 0, 1000, threads=t) for t in (1, 2, 4, 8)]
    assert len(set(values)) == 1
    assert values[0] == py.count("chr1", 0, 1000)


def test_m6_count_and_coverage(bam_pair) -> None:
    py, sr = bam_pair
    assert sr.count("chr1", 0, 1000) == py.count("chr1", 0, 1000)
    assert sr.count("chr1", 590, 650) == py.count("chr1", 590, 650)
    # Placed unmapped mate (unmap_placed @ POS=50) is in pysam nofilter count.
    assert sr.count("chr1", 0, 100) == py.count("chr1", 0, 100)
    assert sr.count("chr1", 0, 100, read_callback="all") == py.count(
        "chr1", 0, 100, read_callback="all"
    )

    pa, pc, pg, pt = py.count_coverage("chr1", 0, 200, quality_threshold=15)
    sa, sc, sg, st = sr.count_coverage("chr1", 0, 200, quality_threshold=15)
    assert list(sa) == list(pa)
    assert list(sc) == list(pc)
    assert list(sg) == list(pg)
    assert list(st) == list(pt)


def test_m6_depth_parallel_matches_serial(bam_pair) -> None:
    _, sr = bam_pair
    d1 = list(sr.depth_numpy("chr1", 0, 200, threads=1))
    d4 = list(sr.depth_numpy("chr1", 0, 200, threads=4))
    assert d1 == d4


def test_m6_coverage_parallel_matches_serial(bam_pair) -> None:
    _, sr = bam_pair
    c1 = sr.count_coverage("chr1", 0, 200, quality_threshold=0, threads=1)
    for threads in (2, 4, 8):
        c_n = sr.count_coverage("chr1", 0, 200, quality_threshold=0, threads=threads)
        for a, b in zip(c1, c_n):
            assert list(a) == list(b), f"coverage mismatch threads={threads}"


def _pysam_aligned_depth(af, contig: str, start: int, stop: int) -> list[int]:
    """samtools depth / rubam get_depths: M/=/X only; N bases count; keep supplementary."""
    depth = [0] * (stop - start)
    for read in af.fetch(contig, start, stop):
        if (
            read.is_unmapped
            or read.is_secondary
            or read.is_qcfail
            or read.is_duplicate
            or read.reference_start is None
            or read.reference_start < 0
        ):
            continue
        ref_pos = int(read.reference_start)
        for op, op_len in read.cigartuples or []:
            if op in (0, 7, 8):
                for i in range(op_len):
                    pos = ref_pos + i
                    if start <= pos < stop:
                        depth[pos - start] += 1
                ref_pos += op_len
            elif op in (2, 3):
                ref_pos += op_len
    return depth


def test_m6_depth_matches_pysam_aligned_oracle(bam_pair) -> None:
    py, sr = bam_pair
    regions = [("chr1", 0, 1000), ("chr1", 200, 250), ("chr1", 390, 480)]
    for contig, start, stop in regions:
        oracle = _pysam_aligned_depth(py, contig, start, stop)
        got = [int(x) for x in sr.depth_numpy(contig, start, stop, threads=1)]
        assert got == oracle, f"depth mismatch in {contig}:{start}-{stop}"
        got4 = [int(x) for x in sr.depth_numpy(contig, start, stop, threads=4)]
        assert got4 == oracle

    # count_coverage omits non-ACGT; depth still counts eqx1 2X (NN) at 520-521.
    a, c, g, t = py.count_coverage("chr1", 520, 522, quality_threshold=0)
    cov = int(a[0] + c[0] + g[0] + t[0])
    dep = int(sr.depth_numpy("chr1", 520, 522)[0])
    assert cov == 0
    assert dep == 1
