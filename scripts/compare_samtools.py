#!/usr/bin/env python3
"""Generate samtools oracle baseline JSON for a BAM (Tier-0 fixture or Tier-1 region)."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
from collect_benchmark_metadata import collect  # noqa: E402


def _run(cmd: list[str]) -> tuple[str, float]:
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return proc.stdout, time.perf_counter() - t0


def parse_idxstats(text: str) -> list[dict[str, Any]]:
    rows = []
    for line in text.splitlines():
        if not line.strip():
            continue
        name, length, mapped, unmapped = line.split("\t")
        rows.append(
            {
                "name": name,
                "length": int(length),
                "mapped": int(mapped),
                "unmapped": int(unmapped),
            }
        )
    return rows


def parse_depth(text: str) -> list[dict[str, Any]]:
    rows = []
    for line in text.splitlines():
        if not line.strip():
            continue
        contig, pos, depth = line.split("\t")[:3]
        rows.append({"contig": contig, "pos": int(pos), "depth": int(depth)})
    return rows


def baseline_bam(
    bam_path: Path,
    *,
    contig: str,
    start: int,
    stop: int,
) -> dict[str, Any]:
    # samtools region is 1-based inclusive
    region = f"{contig}:{start + 1}-{stop}"
    view_c, t_view = _run(["samtools", "view", "-c", str(bam_path), region])
    idxstats, t_idx = _run(["samtools", "idxstats", str(bam_path)])
    # depth/mpileup: keep small windows to avoid multi-MB JSON dumps
    depth_stop = min(stop, start + 1000)
    depth_region = f"{contig}:{start + 1}-{depth_stop}"
    depth, t_depth = _run(
        ["samtools", "depth", "-a", "-r", depth_region, str(bam_path)]
    )
    pile_stop = min(stop, start + 50)
    pile_region = f"{contig}:{start + 1}-{pile_stop}"
    mpileup, t_mpileup = _run(
        ["samtools", "mpileup", "-r", pile_region, "-o", "-", str(bam_path)]
    )
    flagstat, t_flag = _run(["samtools", "flagstat", str(bam_path)])
    depth_rows = parse_depth(depth)

    return {
        "tool": "samtools",
        "samtools_version": subprocess.check_output(
            ["samtools", "--version"], text=True
        )
        .splitlines()[0]
        .strip(),
        "bam": str(bam_path.resolve()),
        "region": {
            "contig": contig,
            "start_0based": start,
            "stop_0based_exclusive": stop,
            "samtools_region": region,
            "depth_region": depth_region,
        },
        "metrics": {
            "view_count": int(view_c.strip() or "0"),
            "idxstats": parse_idxstats(idxstats),
            "depth": depth_rows,
            "depth_sum": sum(r["depth"] for r in depth_rows),
            "mpileup_lines": [ln for ln in mpileup.splitlines() if ln.strip()],
            "flagstat": flagstat.strip().splitlines(),
        },
        "runtime_seconds": {
            "view_c": t_view,
            "idxstats": t_idx,
            "depth": t_depth,
            "mpileup": t_mpileup,
            "flagstat": t_flag,
        },
    }


def infer_contig(bam_path: Path) -> tuple[str, int]:
    out, _ = _run(["samtools", "idxstats", str(bam_path)])
    for row in parse_idxstats(out):
        if row["name"] != "*" and row["length"] > 0:
            return row["name"], row["length"]
    raise RuntimeError(f"no contigs in idxstats for {bam_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--contig", type=str, default=None)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int, default=200)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmark/results/samtools_baseline.json"),
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    args = parser.parse_args()

    if shutil.which("samtools") is None:
        print("ERROR: samtools not found on PATH", file=sys.stderr)
        return 2
    if not args.bam.is_file():
        print(f"ERROR: BAM not found: {args.bam}", file=sys.stderr)
        return 2

    contig = args.contig
    stop = args.stop
    if contig is None:
        contig, length = infer_contig(args.bam)
        stop = min(length, args.stop if args.stop > args.start else 200)

    payload = {
        "metadata": collect(args.repo_root, {"benchmark_kind": "samtools_baseline"}),
        "baseline": baseline_bam(
            args.bam, contig=contig, start=args.start, stop=stop
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
