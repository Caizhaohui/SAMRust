"""M8 recount CLI smoke on fixture-derived candidate sites."""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import pytest

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "small.bam"
ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def recount_bin():
    bin_path = ROOT / "target" / "debug" / "samrust"
    if not bin_path.is_file():
        pytest.skip("samrust binary not built")
    return bin_path


def test_m8_recount_cli_smoke(recount_bin: Path) -> None:
    if not FIXTURE.is_file():
        pytest.skip("missing fixture")
    with tempfile.TemporaryDirectory() as td:
        sites = Path(td) / "sites.bed"
        out = Path(td) / "out.tsv"
        sites.write_text("chr1\t50\t51\tA>C\nchr1\t100\t101\tG>T\n")
        subprocess.run(
            [
                str(recount_bin),
                "recount",
                "--bam",
                str(FIXTURE),
                "--sites",
                str(sites),
                "--sample",
                "fixture",
                "--threads",
                "2",
                "--output",
                str(out),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        text = out.read_text().strip().splitlines()
        assert text[0].startswith("sample\tchrom\tpos\tref\t")
        assert len(text) == 3  # header + 2 sites
