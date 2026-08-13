"""M2 parity gate: samrust vs pysam on Tier-0 fixture."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests" / "fixtures" / "small.bam"
SCRIPT = ROOT / "scripts" / "parity_bam_core.py"


@pytest.mark.skipif(not FIXTURE.is_file(), reason="missing small.bam fixture")
def test_fixture_parity_zero_mismatches() -> None:
    pytest.importorskip("pysam")
    if not (ROOT / "target" / "debug" / "samrust").is_file():
        subprocess.run(
            [
                "cargo",
                "--config",
                'source.crates-io.replace-with="ustc"',
                "--config",
                'source.ustc.registry="sparse+https://mirrors.ustc.edu.cn/crates.io-index/"',
                "build",
                "-p",
                "samrust-cli",
                "-q",
            ],
            cwd=ROOT,
            check=True,
        )
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--bam", str(FIXTURE), "--repo-root", str(ROOT)],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        pytest.fail(proc.stderr or proc.stdout or f"exit {proc.returncode}")
    assert "0 mismatches" in proc.stdout
