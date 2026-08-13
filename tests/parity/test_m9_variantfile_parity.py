"""M9 VariantFile parity vs pysam (and bcftools view -H counts).

Coordinates: fetch / start / stop are 0-based half-open. ``pos`` is 1-based
like pysam.VariantRecord.pos.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
VCF_GZ = FIXTURES / "small.vcf.gz"
VCF_PLAIN = FIXTURES / "small.vcf"
BCF = FIXTURES / "small.bcf"

REGIONS = [
    ("chr1", 0, 100),
    ("chr1", 14, 15),
    ("chr1", 0, 14),
    ("chr1", 249, 261),
    ("chr1", 100, 100),
    ("chr2", 0, 50),
]


def _require_vcf() -> Path:
    if not VCF_GZ.is_file():
        pytest.skip("missing fixture; run: python scripts/prepare_fixture.py")
    return VCF_GZ


def _approx_info(value):
    if isinstance(value, float):
        return pytest.approx(value, rel=1e-5, abs=1e-5)
    if isinstance(value, (list, tuple)):
        return type(value)(_approx_info(v) for v in value)
    return value


def _record_snapshot(rec) -> dict:
    sample = rec.samples[0]
    return {
        "chrom": rec.chrom,
        "pos": rec.pos,
        "start": rec.start,
        "stop": rec.stop,
        "id": rec.id,
        "ref": rec.ref,
        "alts": tuple(rec.alts or ()),
        "qual": rec.qual,
        "filter": list(rec.filter),
        "info": {k: _approx_info(v) for k, v in dict(rec.info).items()},
        "gt": tuple(sample["GT"]),
        "dp": sample["DP"],
        "ad": tuple(sample["AD"]),
    }


def test_m9_header_and_sequential_fields() -> None:
    pysam = pytest.importorskip("pysam")
    samrust = pytest.importorskip("samrust")
    path = _require_vcf()

    py = pysam.VariantFile(str(path))
    sr = samrust.VariantFile(str(path))

    assert list(py.header.samples) == list(sr.header.samples)
    assert list(py.header.contigs) == list(sr.header.contigs)
    assert list(sr.samples) == ["sample1"]

    py_recs = [_record_snapshot(r) for r in py]
    sr_recs = [_record_snapshot(r) for r in sr]
    assert sr_recs == py_recs
    assert len(sr_recs) == 5
    py.close()
    sr.close()


def test_m9_fetch_zero_based_half_open() -> None:
    pysam = pytest.importorskip("pysam")
    samrust = pytest.importorskip("samrust")
    path = _require_vcf()

    py = pysam.VariantFile(str(path))
    sr = samrust.VariantFile(str(path))
    for contig, start, stop in REGIONS:
        py_ids = [(r.chrom, r.pos, r.id, r.start, r.stop) for r in py.fetch(contig, start, stop)]
        sr_ids = [(r.chrom, r.pos, r.id, r.start, r.stop) for r in sr.fetch(contig, start, stop)]
        assert sr_ids == py_ids, f"mismatch in {contig}:{start}-{stop}"
    py.close()
    sr.close()


def test_m9_bcftools_view_count() -> None:
    pysam = pytest.importorskip("pysam")
    samrust = pytest.importorskip("samrust")
    path = _require_vcf()
    bcftools = shutil.which("bcftools")
    if not bcftools:
        pytest.skip("bcftools not on PATH")

    out = subprocess.check_output([bcftools, "view", "-H", str(path)], text=True)
    n_cli = len([ln for ln in out.splitlines() if ln.strip()])
    n_py = sum(1 for _ in pysam.VariantFile(str(path)))
    n_sr = sum(1 for _ in samrust.VariantFile(str(path)))
    assert n_sr == n_py == n_cli == 5


@pytest.mark.skipif(not VCF_PLAIN.is_file(), reason="plain VCF fixture missing")
def test_m9_plain_vcf_iterates() -> None:
    pytest.importorskip("pysam")
    samrust = pytest.importorskip("samrust")
    recs = list(samrust.VariantFile(str(VCF_PLAIN)))
    assert len(recs) == 5
    assert recs[0].chrom == "chr1"
    assert recs[0].start == 14


@pytest.mark.skipif(not BCF.is_file(), reason="BCF fixture missing")
def test_m9_bcf_matches_vcf_gz() -> None:
    samrust = pytest.importorskip("samrust")
    gz = [_record_snapshot(r) for r in samrust.VariantFile(str(_require_vcf()))]
    bcf = [_record_snapshot(r) for r in samrust.VariantFile(str(BCF))]
    assert bcf == gz


def test_m9_context_manager_and_missing_contig() -> None:
    samrust = pytest.importorskip("samrust")
    path = _require_vcf()
    with samrust.VariantFile(str(path)) as vf:
        assert vf.header.samples == ["sample1"]
    with pytest.raises(Exception):
        samrust.VariantFile(str(path)).fetch("missing_contig", 0, 10)


def test_m9_fetch_contig_without_length_in_header(tmp_path) -> None:
    """v0.1.1 P0b: fetch(contig) on a header without contig length must not
    silently return empty; unbounded fetch falls back to a sequential scan."""
    pysam = pytest.importorskip("pysam")
    samrust = pytest.importorskip("samrust")

    vcf = tmp_path / "nolength.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
        "chr1\t20\t.\tC\tT\t60\tPASS\t.\tGT\t1/1\n"
    )
    gz = tmp_path / "nolength.vcf.gz"
    pysam.tabix_compress(str(vcf), str(gz), force=True)
    pysam.tabix_index(str(gz), preset="vcf", force=True)

    with pysam.VariantFile(str(gz)) as py:
        expected = [r.pos for r in py.fetch("chr1")]
    with samrust.VariantFile(str(gz)) as sr:
        assert [r.pos for r in sr.fetch("chr1")] == expected
        # explicit stop still uses the index
        assert [r.pos for r in sr.fetch("chr1", 0, 15)] == [10]
        # negative coordinates -> ValueError (pysam semantics)
        with pytest.raises(ValueError):
            list(sr.fetch("chr1", -5, 15))
