#!/usr/bin/env python3
"""Generate rubam baseline JSON for a BAM (Tier-0 fixture or Tier-1 region).

Coordinate notes (rubam 0.3.x):
- ``AlignmentFile.fetch/count/count_coverage/pileup``: 0-based half-open (pysam-like)
- Free functions ``get_depths`` / ``pileup_bases`` / ``count_reads``: **1-based inclusive**

See DEVELOPMENT_PLAN.md §15.2 / §19.6.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

try:
    import rubam
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "rubam is required: pip install 'rubam>=0.3.13'\n"
        f"ImportError: {exc}"
    ) from exc

sys.path.insert(0, str(Path(__file__).resolve().parent))
from collect_benchmark_metadata import collect  # noqa: E402


def _rubam_version() -> str:
    return str(getattr(rubam, "__version__", getattr(rubam, "version", "unknown")))


def _safe_int(value: Any, default: int | None = None) -> int | None:
    if value is None:
        return default
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def record_identity(read: Any) -> dict[str, Any]:
    qseq = getattr(read, "query_sequence", None)
    qlen = 0
    if qseq is not None:
        qlen = _safe_int(getattr(read, "query_length", None), default=len(qseq)) or 0
    return {
        "qname": read.query_name,
        "flag": _safe_int(read.flag, 0) or 0,
        "reference_id": _safe_int(getattr(read, "reference_id", None)),
        "reference_name": None
        if getattr(read, "is_unmapped", False)
        else getattr(read, "reference_name", None),
        "reference_start": _safe_int(getattr(read, "reference_start", None), -1),
        "mapping_quality": _safe_int(getattr(read, "mapping_quality", None), 0) or 0,
        "cigar": getattr(read, "cigarstring", None),
        "query_length": qlen,
        "is_duplicate": bool(getattr(read, "is_duplicate", False)),
        "is_secondary": bool(getattr(read, "is_secondary", False)),
        "is_supplementary": bool(getattr(read, "is_supplementary", False)),
        "is_qcfail": bool(getattr(read, "is_qcfail", False)),
        "is_unmapped": bool(getattr(read, "is_unmapped", False)),
    }


def pileup_via_alignment_file(
    bam: Any, contig: str, start: int, stop: int, *, min_bq: int = 0, threads: int = 1
) -> list[dict[str, Any]]:
    """0-based half-open region via AlignmentFile.pileup."""
    rows: list[dict[str, Any]] = []
    for col in bam.pileup(
        contig,
        start,
        stop,
        min_mapq=0,
        min_bq=min_bq,
        truncate=True,
        num_threads=max(1, threads),
    ):
        pos = int(col.reference_pos)
        if pos < start or pos >= stop:
            continue
        rows.append(
            {
                "pos": pos,
                "nsegments": int(getattr(col, "nsegments", 0) or 0),
                "A": int(getattr(col, "a", 0) or 0),
                "C": int(getattr(col, "c", 0) or 0),
                "G": int(getattr(col, "g", 0) or 0),
                "T": int(getattr(col, "t", 0) or 0),
                "N": int(getattr(col, "n", 0) or 0),
                "depth": int(getattr(col, "depth", 0) or 0),
            }
        )
    return rows


def pileup_bases_1based(
    bam_path: Path, contig: str, start0: int, stop0: int, *, min_bq: int = 0, threads: int = 1
) -> dict[str, Any]:
    """Convert 0-based half-open → rubam 1-based inclusive ``pileup_bases``."""
    if stop0 <= start0:
        return {"positions": [], "A": [], "C": [], "G": [], "T": [], "N": [], "depth": []}
    start1 = start0 + 1
    end1 = stop0  # inclusive end of last base in [start0, stop0)
    positions, a, c, g, t, n, depth = rubam.pileup_bases(
        str(bam_path),
        contig,
        start1,
        end1,
        step=1,
        min_mapq=0,
        min_bq=min_bq,
        max_depth=1_000_000,
        num_threads=max(1, threads),
        flag_filter=0,
    )
    return {
        "positions_1based": list(map(int, positions)),
        "A": list(map(int, a)),
        "C": list(map(int, c)),
        "G": list(map(int, g)),
        "T": list(map(int, t)),
        "N": list(map(int, n)),
        "depth": list(map(int, depth)),
    }


def baseline_bam(
    bam_path: Path,
    *,
    contig: str | None,
    start: int,
    stop: int,
    max_records: int,
    threads: int,
    min_bq: int,
) -> dict[str, Any]:
    t0 = time.perf_counter()
    with rubam.AlignmentFile(str(bam_path), "rb") as bam:
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

        count = int(bam.count(contig, start, stop))
        coverage = bam.count_coverage(contig, start, stop, quality_threshold=min_bq)
        cov_lists = [list(map(int, arr)) for arr in coverage]
        region_len = max(0, stop - start)
        store_full_coverage = region_len <= 1000
        pile_stop = min(stop, start + 50)
        pile_af = pileup_via_alignment_file(
            bam, contig, start, pile_stop, min_bq=min_bq, threads=threads
        )

    # Free-function path (1-based inclusive)
    count_reads = int(
        rubam.count_reads(
            str(bam_path),
            contig,
            start + 1,
            stop,  # inclusive
            min_mapq=0,
            flag_required=0,
            flag_filtered=0,
        )
    )
    depths = rubam.get_depths(
        str(bam_path),
        contig,
        start + 1,
        stop if stop > start else start + 1,
        step=1,
        min_mapq=0,
        min_bq=min_bq,
        max_depth=1_000_000,
        num_threads=max(1, threads),
    )
    depth_sum = int(sum(map(int, depths[1]))) if stop > start else 0
    pile_fn = pileup_bases_1based(
        bam_path, contig, start, pile_stop, min_bq=min_bq, threads=threads
    )

    elapsed = time.perf_counter() - t0
    metrics: dict[str, Any] = {
        "count": count,
        "count_reads_1based_nofilter": count_reads,
        "count_coverage_sums": [sum(x) for x in cov_lists],
        "get_depths_sum": depth_sum,
        "pileup_alignmentfile_head": pile_af,
        "pileup_bases_head": {
            k: (v[:50] if isinstance(v, list) else v) for k, v in pile_fn.items()
        },
    }
    if store_full_coverage:
        metrics["count_coverage"] = cov_lists
    else:
        metrics["count_coverage"] = "omitted_large_region"
        metrics["count_coverage_first_100"] = [arr[:100] for arr in cov_lists]

    return {
        "tool": "rubam",
        "rubam_version": _rubam_version(),
        "bam": str(bam_path.resolve()),
        "region": {"contig": contig, "start": start, "stop": stop},
        "threads": threads,
        "min_bq": min_bq,
        "header": {
            "nreferences": len(refs),
            "references": refs,
            "lengths": lengths,
        },
        "records_head": identities,
        "metrics": metrics,
        "runtime_seconds": elapsed,
        "notes": {
            "alignmentfile_coords": "0-based half-open",
            "free_function_coords": "1-based inclusive",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--contig", type=str, default=None)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int, default=200)
    parser.add_argument("--max-records", type=int, default=100)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--min-bq", type=int, default=0)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmark/results/rubam_baseline.json"),
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
        "metadata": collect(args.repo_root, {"benchmark_kind": "rubam_baseline"}),
        "baseline": baseline_bam(
            args.bam,
            contig=args.contig,
            start=args.start,
            stop=args.stop,
            max_records=args.max_records,
            threads=args.threads,
            min_bq=args.min_bq,
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
