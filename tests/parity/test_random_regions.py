"""Random-region differential testing vs pysam (DEVELOPMENT_PLAN §17).

Fixture mode always runs (Tier-0, ~100 regions). Real-data mode activates when
``benchmark/real_data_manifest.tsv`` points at an existing BAM + index and
``SAMRUST_REAL_DATA=1`` is set (Rule 12: never fabricate; Rule 13: heavy real
data belongs on qcpu_18i — this mode is for scheduled/benchmark jobs).

For every region we compare, at threads 1 and 4:

- ``count`` (nofilter) vs ``pysam.count``
- ``fetch`` record identity vs ``pysam.fetch``
- ``count_coverage`` (quality_threshold=0) vs ``pysam.count_coverage``
- ``depth_numpy`` vs the samtools-depth oracle (CIGAR M/=/X)
- ``pileup_counts`` vs normalized ``pysam.pileup``
"""

from __future__ import annotations

import os
import random
from pathlib import Path

import pytest

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "small.bam"
MANIFEST = Path(__file__).resolve().parents[2] / "benchmark" / "real_data_manifest.tsv"

pysam = pytest.importorskip("pysam")
samrust = pytest.importorskip("samrust")

_FLAG_FILTER = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
SEED = 20260813


def _rand_regions(rng: random.Random, contigs: list[tuple[str, int]], n: int, sizes):
    regions = []
    for _ in range(n):
        contig, clen = rng.choice(contigs)
        size = min(rng.choice(sizes), clen)
        start = rng.randint(0, max(0, clen - size))
        regions.append((contig, start, start + size))
    return regions


def _pysam_aligned_depth(af, contig: str, start: int, stop: int) -> list[int]:
    depth = [0] * max(0, stop - start)
    if stop <= start:
        return depth
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


def _pysam_pileup_counts(af, contig: str, start: int, stop: int) -> dict[str, list[int]]:
    length = max(0, stop - start)
    out = {k: [0] * length for k in ("A", "C", "G", "T", "N", "depth")}
    if length == 0:
        return out
    for col in af.pileup(
        contig,
        start,
        stop,
        truncate=True,
        min_base_quality=0,
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
            qpos = pr.query_position
            if qpos is None:
                continue
            base = (pr.alignment.query_sequence or "N")[qpos].upper()
            out["depth"][idx] += 1
            out[base if base in "ACGT" else "N"][idx] += 1
    return out


def _check_region(py, sr, contig: str, start: int, stop: int, with_fetch: bool) -> None:
    where = f"{contig}:{start}-{stop}"
    for threads in (1, 4):
        assert sr.count(contig, start, stop, threads=threads) == py.count(
            contig, start, stop
        ), f"count {where} T={threads}"

        pa, pc, pg, pt = py.count_coverage(contig, start, stop, quality_threshold=0)
        sa, sc, sg, st = sr.count_coverage(contig, start, stop, quality_threshold=0, threads=threads)
        assert (list(sa), list(sc), list(sg), list(st)) == (
            list(pa),
            list(pc),
            list(pg),
            list(pt),
        ), f"coverage {where} T={threads}"

        oracle_depth = _pysam_aligned_depth(py, contig, start, stop)
        got_depth = [int(x) for x in sr.depth_numpy(contig, start, stop, threads=threads)]
        assert got_depth == oracle_depth, f"depth {where} T={threads}"

        oracle_pu = _pysam_pileup_counts(py, contig, start, stop)
        got_pu = sr.pileup_counts(contig, start, stop, threads=threads)
        for k in ("A", "C", "G", "T", "N", "depth"):
            assert [int(x) for x in got_pu[k]] == oracle_pu[k], f"pileup.{k} {where} T={threads}"

    if with_fetch:
        expected = [
            (r.query_name, r.flag, r.reference_start, r.cigarstring)
            for r in py.fetch(contig, start, stop)
        ]
        got = [
            (r.query_name, r.flag, r.reference_start, r.cigarstring)
            for r in sr.fetch(contig, start, stop)
        ]
        assert got == expected, f"fetch {where}"


def test_random_regions_fixture() -> None:
    if not FIXTURE.is_file():
        pytest.skip("missing fixture")
    rng = random.Random(SEED)
    with pysam.AlignmentFile(str(FIXTURE), "rb") as py:
        contigs = list(zip(py.references, py.lengths))
        regions = _rand_regions(rng, contigs, 100, sizes=[1, 7, 100, 500, 1000])
        with samrust.AlignmentFile(str(FIXTURE), "rb") as sr:
            for contig, start, stop in regions:
                _check_region(py, sr, contig, start, stop, with_fetch=True)


@pytest.mark.skipif(
    os.environ.get("SAMRUST_REAL_DATA") != "1",
    reason="set SAMRUST_REAL_DATA=1 to run real-data random regions (qcpu_18i jobs)",
)
def test_random_regions_real_data() -> None:
    if not MANIFEST.is_file():
        pytest.skip("real-data manifest missing (Rule 12)")
    lines = [ln for ln in MANIFEST.read_text().splitlines()[1:] if ln.strip()]
    bam = None
    for ln in lines:
        cols = ln.split("\t")
        if len(cols) >= 3 and Path(cols[1]).is_file() and Path(cols[2]).is_file():
            bam = cols[1]
            break
    if bam is None:
        pytest.skip("no usable real BAM in manifest (Rule 12)")

    rng = random.Random(SEED)
    with pysam.AlignmentFile(bam, "rb") as py:
        contigs = [(c, l) for c, l in zip(py.references, py.lengths) if l >= 1000]
        regions = _rand_regions(rng, contigs, 200, sizes=[1000, 10_000, 100_000])
        with samrust.AlignmentFile(bam, "rb") as sr:
            for contig, start, stop in regions:
                # bound fetch memory on large windows
                _check_region(py, sr, contig, start, stop, with_fetch=(stop - start) <= 10_000)
