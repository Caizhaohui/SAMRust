#!/usr/bin/env python3
"""Collect host / toolchain / dataset metadata for benchmark JSON payloads."""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _run(cmd: list[str]) -> str | None:
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
        return out.strip().splitlines()[0] if out.strip() else ""
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None


def _git_commit(repo: Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None


def _mem_total_bytes() -> int | None:
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemTotal:"):
                    kb = int(line.split()[1])
                    return kb * 1024
    except OSError:
        return None
    return None


def _cpu_model() -> str | None:
    try:
        with open("/proc/cpuinfo") as fh:
            for line in fh:
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        return None
    return None


def collect(repo_root: Path, extra: dict[str, Any] | None = None) -> dict[str, Any]:
    meta: dict[str, Any] = {
        "date_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": _git_commit(repo_root),
        "hostname": platform.node(),
        "os": platform.platform(),
        "linux_kernel": platform.release(),
        "cpu_model": _cpu_model(),
        "cpu_count": os.cpu_count(),
        "ram_bytes": _mem_total_bytes(),
        "filesystem": None,
        "storage_type": "unknown",
        "rustc": _run(["rustc", "--version"]),
        "cargo": _run(["cargo", "--version"]),
        "python": sys.version.split()[0],
        "pysam": None,
        "rubam": None,
        "samtools": _run(["samtools", "--version"]),
        "bcftools": _run(["bcftools", "--version"]),
        "samrust": None,
        "tmpdir": os.environ.get("TMPDIR"),
    }
    try:
        import pysam  # type: ignore

        meta["pysam"] = getattr(pysam, "__version__", "unknown")
    except Exception:
        meta["pysam"] = None
    try:
        import rubam  # type: ignore

        meta["rubam"] = getattr(rubam, "__version__", "unknown")
    except Exception:
        meta["rubam"] = None
    try:
        import samrust  # type: ignore

        meta["samrust"] = getattr(samrust, "__version__", "unknown")
    except Exception:
        meta["samrust"] = None

    # Best-effort filesystem type for repo path
    try:
        st = os.statvfs(repo_root)
        meta["filesystem"] = {
            "f_frsize": st.f_frsize,
            "f_blocks": st.f_blocks,
            "f_bavail": st.f_bavail,
        }
    except OSError:
        pass

    if extra:
        meta.update(extra)
    return meta


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmark/results/host_metadata.json"),
    )
    args = parser.parse_args()
    payload = collect(args.repo_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {args.output}")
    which = {k: shutil.which(k) for k in ("samtools", "bcftools", "python")}
    print("PATH tools:", which)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
