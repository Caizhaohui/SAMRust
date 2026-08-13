#!/usr/bin/env python3
"""Generate pysam oracle baseline JSON for a BAM (Tier-0 fixture or Tier-1 region)."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

import pysam

# Allow `python scripts/compare_pysam.py` to import sibling helpers.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from collect_benchmark_metadata import collect  # noqa: E402


def record_identity(read: pysam.AlignedSegment) -> dict[str, Any]:
    return {
        "qname": read.query_name,
        "flag": int(read.flag),
        "reference_id": int(read.reference_id),
        "reference_name": None if read.is_unmapped else read.reference_name,
        "reference_start": int(read.reference_start),
        "mapping_quality": int(read.mapping_quality),
        "cigar": read.cigarstring,
        "query_length": int(read.query_length) if read.query_sequence is not None else 0,
        "is_duplicate": bool(read.is_duplicate),
        "is_secondary": bool(read.is_secondary),
        "is_supplementary": bool(read.is_supplementary),
        "is_qcfail": bool(read.is_qcfail),
        "is_unmapped": bool(read.is_unmapped),
        "tags": {k: _jsonable(v) for k, v in (read.get_tags() or [])},
    }


def _jsonable(value: Any) -> Any:
    if isinstance(value, (bytes, bytearray)):
        return value.decode("ascii", errors="replace")
    if isinstance(value, (int, float, str, bool)) or value is None:
        return value
    return str(value)


def pileup_snapshot(
    bam: pysam.AlignmentFile, contig: str, start: int, stop: int
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for col in bam.pileup(contig, start, stop, truncate=True, min_base_quality=0):
        if col.reference_pos < start or col.reference_pos >= stop:
            continue
        bases = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip:
                continue
            base = (pr.alignment.query_sequence[pr.query_position] or "N").upper()
            if base not in bases:
                base = "N"
            bases[base] += 1
        rows.append(
            {
                "pos": int(col.reference_pos),
                "nsegments": int(col.nsegments),
                **bases,
            }
        )
    return rows


def baseline_bam(
    bam_path: Path,
    *,
    contig: str | None,
    start: int,
    stop: int,
    max_records: int,
) -> dict[str, Any]:
    t0 = time.perf_counter()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        refs = list(bam.references)
        lengths = [int(x) for x in bam.lengths]
        if contig is None:
            contig = refs[0]
            start = 0
            stop = min(lengths[0], max(stop - start, 1) if stop > start else 200)

        identities: list[dict[str, Any]] = []
        for i, read in enumerate(bam.fetch(until_eof=True)):
            if i >= max_records:
                break
            identities.append(record_identity(read))

        count = bam.count(contig, start, stop)
        coverage = bam.count_coverage(contig, start, stop, quality_threshold=0)
        cov_lists = [list(map(int, arr)) for arr in coverage]
        region_len = max(0, stop - start)
        # Keep full arrays only for small windows (fixtures / tiny Tier-1 probes).
        store_full_coverage = region_len <= 1000
        pile = pileup_snapshot(bam, contig, start, min(stop, start + 50))

        mapped = None
        unmapped = None
        try:
            mapped = int(bam.mapped)
            unmapped = int(bam.unmapped)
        except ValueError:
            # Some BAMs without index stats
            pass

    elapsed = time.perf_counter() - t0
    metrics: dict[str, Any] = {
        "count": int(count),
        "count_coverage_sums": [sum(x) for x in cov_lists],
        "pileup_head": pile,
    }
    if store_full_coverage:
        metrics["count_coverage"] = cov_lists
    else:
        metrics["count_coverage"] = "omitted_large_region"
        metrics["count_coverage_first_100"] = [arr[:100] for arr in cov_lists]

    return {
        "tool": "pysam",
        "pysam_version": pysam.__version__,
        "bam": str(bam_path.resolve()),
        "region": {"contig": contig, "start": start, "stop": stop},
        "header": {
            "nreferences": len(refs),
            "references": refs,
            "lengths": lengths,
            "mapped": mapped,
            "unmapped": unmapped,
        },
        "records_head": identities,
        "metrics": metrics,
        "runtime_seconds": elapsed,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--contig", type=str, default=None)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int, default=200)
    parser.add_argument("--max-records", type=int, default=100)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmark/results/pysam_baseline.json"),
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    args = parser.parse_args()

    if not args.bam.is_file():
        print(f"ERROR: BAM not found: {args.bam}", file=sys.stderr)
        return 2

    payload = {
        "metadata": collect(args.repo_root, {"benchmark_kind": "pysam_baseline"}),
        "baseline": baseline_bam(
            args.bam,
            contig=args.contig,
            start=args.start,
            stop=args.stop,
            max_records=args.max_records,
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
