"""Tier-0 fixture integrity tests (no SAMRust BAM engine yet)."""

from __future__ import annotations

from pathlib import Path

import pytest

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"


@pytest.fixture(scope="module")
def fixture_paths() -> dict[str, Path]:
    required = {
        "fa": FIXTURES / "small.fa",
        "fai": FIXTURES / "small.fa.fai",
        "bam": FIXTURES / "small.bam",
        "bai": FIXTURES / "small.bam.bai",
        "vcf": FIXTURES / "small.vcf.gz",
        "tbi": FIXTURES / "small.vcf.gz.tbi",
    }
    missing = [k for k, p in required.items() if not p.is_file()]
    if missing:
        pytest.skip(
            "fixtures missing; run: python scripts/prepare_fixture.py "
            f"(missing: {', '.join(missing)})"
        )
    return required


def test_fixture_files_exist(fixture_paths: dict[str, Path]) -> None:
    assert all(p.is_file() for p in fixture_paths.values())


def test_bam_has_expected_edge_cases(fixture_paths: dict[str, Path]) -> None:
    pysam = pytest.importorskip("pysam")
    with pysam.AlignmentFile(str(fixture_paths["bam"]), "rb") as bam:
        flags = {r.query_name: r.flag for r in bam.fetch(until_eof=True)}
    assert "pair1" in flags
    assert flags.get("dup1") == 1024
    assert flags.get("sec1") == 256
    assert flags.get("sup1") == 2048
    assert flags.get("qcfail1") == 512
    assert flags.get("unmap1") == 4
    assert flags.get("unmap_placed") == 133


def test_vcf_multiallelic(fixture_paths: dict[str, Path]) -> None:
    pysam = pytest.importorskip("pysam")
    with pysam.VariantFile(str(fixture_paths["vcf"])) as vf:
        recs = list(vf.fetch())
    assert len(recs) >= 5
    multi = [r for r in recs if len(r.alts or []) > 1]
    assert multi, "expected multi-allelic record"
