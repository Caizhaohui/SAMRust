"""M7 pileup parity: pysam-normalized serial + parallel bit-exact."""

from __future__ import annotations

from pathlib import Path

import pytest

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "small.bam"

# Match samrust PileupFilter::default() exclude flags.
# BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP|BAM_FSUPPLEMENTARY
_FLAG_FILTER = 0x4 | 0x100 | 0x200 | 0x400 | 0x800


def _pysam_base_counts(bam, contig: str, start: int, stop: int, min_bq: int = 0):
    """Normalized A/C/G/T/N + depth, skipping del/refskip (same as SAMRust base pileup)."""
    length = stop - start
    a = [0] * length
    c = [0] * length
    g = [0] * length
    t = [0] * length
    n = [0] * length
    depth = [0] * length
    for col in bam.pileup(
        contig,
        start,
        stop,
        truncate=True,
        min_base_quality=min_bq,
        flag_filter=_FLAG_FILTER,
        stepper="all",
    ):
        pos = col.reference_pos
        if pos < start or pos >= stop:
            continue
        idx = pos - start
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip:
                continue
            q = pr.alignment.query_qualities
            qpos = pr.query_position
            if qpos is None:
                continue
            if q is not None and q[qpos] < min_bq:
                continue
            base = (pr.alignment.query_sequence or "N")[qpos].upper()
            depth[idx] += 1
            if base == "A":
                a[idx] += 1
            elif base == "C":
                c[idx] += 1
            elif base == "G":
                g[idx] += 1
            elif base == "T":
                t[idx] += 1
            else:
                n[idx] += 1
    return a, c, g, t, n, depth


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


@pytest.mark.parametrize(
    "region",
    [
        ("chr1", 0, 150),
        ("chr1", 180, 230),  # indel-adjacent
        ("chr1", 590, 650),  # secondary/dup neighborhood
        ("chr1", 0, 700),
        ("chr2", 0, 50),
    ],
)
def test_m7_serial_pileup_matches_pysam_normalized(bam_pair, region) -> None:
    py, sr = bam_pair
    contig, start, stop = region
    pa, pc, pg, pt, pn, pd = _pysam_base_counts(py, contig, start, stop, min_bq=0)
    counts = sr.pileup_counts(contig, start, stop, min_base_quality=0, threads=1)
    assert list(counts["A"]) == pa
    assert list(counts["C"]) == pc
    assert list(counts["G"]) == pg
    assert list(counts["T"]) == pt
    assert list(counts["N"]) == pn
    assert list(counts["depth"]) == pd


@pytest.mark.parametrize("min_bq", [0, 10, 20])
def test_m7_bq_filter_matches_pysam(bam_pair, min_bq) -> None:
    py, sr = bam_pair
    contig, start, stop = "chr1", 0, 200
    pa, pc, pg, pt, pn, pd = _pysam_base_counts(py, contig, start, stop, min_bq=min_bq)
    counts = sr.pileup_counts(contig, start, stop, min_base_quality=min_bq, threads=1)
    assert list(counts["A"]) == pa
    assert list(counts["C"]) == pc
    assert list(counts["G"]) == pg
    assert list(counts["T"]) == pt
    assert list(counts["N"]) == pn
    assert list(counts["depth"]) == pd


def test_m7_parallel_pileup_bit_exact_vs_serial(bam_pair) -> None:
    _, sr = bam_pair
    contig, start, stop = "chr1", 0, 700
    serial = sr.pileup_counts(contig, start, stop, min_base_quality=0, threads=1)
    for threads in (2, 4, 8):
        parallel = sr.pileup_counts(contig, start, stop, min_base_quality=0, threads=threads)
        for key in ("A", "C", "G", "T", "N", "depth"):
            assert list(parallel[key]) == list(serial[key]), f"{key} mismatch threads={threads}"
