"""M11 CRAM evaluation: sequential + indexed fetch vs pysam; stats stay BAM-only."""

from __future__ import annotations

from pathlib import Path

import pytest

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures"
CRAM = FIXTURE / "small.cram"
FASTA = FIXTURE / "small.fa"
BAM = FIXTURE / "small.bam"


@pytest.fixture(scope="module")
def cram_pair():
    pysam = pytest.importorskip("pysam")
    samrust = pytest.importorskip("samrust")
    if not CRAM.is_file() or not FASTA.is_file():
        pytest.skip("missing small.cram fixture")
    return (
        pysam.AlignmentFile(str(CRAM), "rc", reference_filename=str(FASTA)),
        samrust.AlignmentFile(str(CRAM), "rc", reference_filename=str(FASTA)),
    )


def test_m11_cram_header_and_sequential(cram_pair) -> None:
    py, sr = cram_pair
    assert list(py.references) == sr.references
    assert list(py.lengths) == sr.lengths
    py_ids = [(r.query_name, r.flag, r.reference_start, r.cigarstring) for r in py]
    sr_ids = [(r.query_name, r.flag, r.reference_start, r.cigarstring) for r in sr]
    assert sr_ids == py_ids


def test_m11_cram_fetch_matches_pysam(cram_pair) -> None:
    py, sr = cram_pair
    regions = [
        ("chr1", 0, 100),
        ("chr1", 200, 250),
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

    # Empty half-open [100, 100) is empty in SAMRust (and pysam BAM).
    # pysam CRAM/HTSlib still emits overlaps at that point — do not copy that.
    assert [
        (r.query_name, r.flag, r.reference_start, r.cigarstring)
        for r in sr.fetch("chr1", 100, 100)
    ] == []


def test_m11_cram_stats_not_implemented(cram_pair) -> None:
    _, sr = cram_pair
    with pytest.raises(RuntimeError, match="CRAM"):
        sr.count("chr1", 0, 100)
    with pytest.raises(RuntimeError, match="CRAM"):
        sr.count_coverage("chr1", 0, 100)
    with pytest.raises(RuntimeError, match="CRAM"):
        sr.depth_numpy("chr1", 0, 100)
    with pytest.raises(RuntimeError, match="CRAM"):
        sr.pileup_counts("chr1", 0, 100)


def test_m11_bam_still_opens_without_reference() -> None:
    samrust = pytest.importorskip("samrust")
    if not BAM.is_file():
        pytest.skip("missing fixture")
    af = samrust.AlignmentFile(str(BAM), "rb")
    assert af.nreferences == 2


def test_m11_cram_rb_mode_and_sibling_fasta() -> None:
    samrust = pytest.importorskip("samrust")
    if not CRAM.is_file() or not FASTA.is_file():
        pytest.skip("missing small.cram fixture")
    af = samrust.AlignmentFile(str(CRAM), "rb", reference_filename=str(FASTA))
    assert af.nreferences == 2
    # same directory as small.fa
    af2 = samrust.AlignmentFile(str(CRAM), "rc")
    assert af2.references == af.references
