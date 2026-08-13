"""Smoke tests for discovery helpers (uses local config when present)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts"
LOCAL_CONFIG = ROOT / "benchmark" / "configs" / "myceliophthora.local.yaml"


def _load(name: str):
    path = SCRIPTS / name
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_discover_module_imports() -> None:
    pytest.importorskip("yaml")
    mod = _load("discover_real_data.py")
    assert hasattr(mod, "discover")


@pytest.mark.skipif(not LOCAL_CONFIG.is_file(), reason="no local real-data config")
def test_discover_finds_markdup_bams() -> None:
    pytest.importorskip("yaml")
    mod = _load("discover_real_data.py")
    cfg = mod.load_config(LOCAL_CONFIG)
    rows = mod.discover(cfg)
    assert rows, "expected at least one discovered BAM"
    assert all(Path(r.bam).is_file() for r in rows)
    indexed = [r for r in rows if r.index]
    assert indexed, "expected indexed BAMs"
