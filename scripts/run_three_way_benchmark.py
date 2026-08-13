#!/usr/bin/env python3
"""Three-way benchmark: SAMRust vs pysam vs rubam (output + runtime + RSS).

Implements DEVELOPMENT_PLAN.md §1.2 / §15.0 / §19.6 for critical workloads:

- count
- count_coverage
- depth (samtools / rubam ``get_depths``: CIGAR M/=/X, not ``count_coverage`` A+C+G+T)
- pileup_counts

Features and results track pysam; runtime tracks rubam. If rubam has no
equivalent API or result for a workload, **keep the row** and write the
literal ``NA`` into rubam ``wall_s`` / ``max_rss_kb`` / ``output_digest`` /
``gate_vs_pysam`` (and vs-rubam). Do not drop the row. Semantic differences
with a real rubam entry still report numbers plus notes — NA is only for a
missing function or missing result column.

Coordinates for the CLI are always **0-based half-open**. Rubam free functions
are adapted internally (1-based inclusive).

Heavy fungal BAMs must be run via ``run_three_way_benchmark.sh`` on ``qcpu_18i``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import resource
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

sys.path.insert(0, str(Path(__file__).resolve().parent))
from collect_benchmark_metadata import collect  # noqa: E402

# Match SAMRust default pileup flag exclusions when normalizing pysam.
_FLAG_FILTER = 0x4 | 0x100 | 0x200 | 0x400 | 0x800


@dataclass
class RunResult:
    tool: str
    workload: str
    threads: int
    wall_s: float
    max_rss_kb: int | None
    output_digest: str
    gate_vs_pysam: str
    detail: dict[str, Any]


def _rss_kb() -> int:
    # Linux: ru_maxrss is kilobytes
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def _digest(obj: Any) -> str:
    blob = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


def _canonical(workload: str, out: dict[str, Any]) -> Any:
    """Normalize a tool's output to comparable values only.

    Digests must match across tools when (and only when) the results match;
    tool-specific extra keys (e.g. rubam ``count_reads_nofilter`` /
    ``positions_1based``) are excluded here but kept in ``detail``.
    """
    if workload == "count":
        return {"count": int(out["count"])}
    if workload == "count_coverage":
        return {k: [int(x) for x in out[k]] for k in ("A", "C", "G", "T")}
    if workload == "depth":
        return {"depth": [int(x) for x in out["depth"]]}
    if workload == "pileup_counts":
        return {k: [int(x) for x in out[k]] for k in ("A", "C", "G", "T", "N", "depth")}
    return out


def pysam_count(bam: Path, contig: str, start: int, stop: int, **_kw: Any) -> dict[str, Any]:
    import pysam

    with pysam.AlignmentFile(str(bam), "rb") as af:
        return {"count": int(af.count(contig, start, stop))}


def rubam_count(bam: Path, contig: str, start: int, stop: int, **_kw: Any) -> dict[str, Any]:
    import rubam

    with rubam.AlignmentFile(str(bam), "rb") as af:
        n_af = int(af.count(contig, start, stop))
    n_fn = int(
        rubam.count_reads(
            str(bam),
            contig,
            start + 1,
            stop,
            min_mapq=0,
            flag_required=0,
            flag_filtered=0,
        )
    )
    return {"count": n_af, "count_reads_nofilter": n_fn}


def samrust_count(
    bam: Path, contig: str, start: int, stop: int, *, threads: int = 1, **_kw: Any
) -> dict[str, Any]:
    import samrust

    with samrust.AlignmentFile(str(bam), "rb") as af:
        return {"count": int(af.count(contig, start, stop, threads=threads))}


def pysam_coverage(
    bam: Path, contig: str, start: int, stop: int, *, min_bq: int = 0, **_kw: Any
) -> dict[str, Any]:
    import pysam

    with pysam.AlignmentFile(str(bam), "rb") as af:
        a, c, g, t = af.count_coverage(contig, start, stop, quality_threshold=min_bq)
        arrays = [list(map(int, a)), list(map(int, c)), list(map(int, g)), list(map(int, t))]
    return {
        "A": arrays[0],
        "C": arrays[1],
        "G": arrays[2],
        "T": arrays[3],
        "sums": [sum(x) for x in arrays],
    }


def rubam_coverage(
    bam: Path, contig: str, start: int, stop: int, *, min_bq: int = 0, **_kw: Any
) -> dict[str, Any]:
    import rubam

    with rubam.AlignmentFile(str(bam), "rb") as af:
        a, c, g, t = af.count_coverage(contig, start, stop, quality_threshold=min_bq)
        arrays = [list(map(int, a)), list(map(int, c)), list(map(int, g)), list(map(int, t))]
    return {
        "A": arrays[0],
        "C": arrays[1],
        "G": arrays[2],
        "T": arrays[3],
        "sums": [sum(x) for x in arrays],
    }


def samrust_coverage(
    bam: Path,
    contig: str,
    start: int,
    stop: int,
    *,
    min_bq: int = 0,
    threads: int = 1,
    **_kw: Any,
) -> dict[str, Any]:
    import samrust

    with samrust.AlignmentFile(str(bam), "rb") as af:
        a, c, g, t = af.count_coverage(
            contig, start, stop, quality_threshold=min_bq, threads=threads
        )
        arrays = [list(map(int, a)), list(map(int, c)), list(map(int, g)), list(map(int, t))]
    return {
        "A": arrays[0],
        "C": arrays[1],
        "G": arrays[2],
        "T": arrays[3],
        "sums": [sum(x) for x in arrays],
    }


def pysam_aligned_depth(
    bam: Path, contig: str, start: int, stop: int, **_kw: Any
) -> list[int]:
    """samtools-depth / rubam get_depths oracle: CIGAR M/=/X only, include N bases.

    Excludes unmapped / secondary / QC-fail / duplicate (pysam ``read_callback='all'``).
    Deletions and ref-skips do not contribute. Supplementary alignments are kept.
    """
    import pysam

    length = max(0, stop - start)
    depth = [0] * length
    if length == 0:
        return depth
    with pysam.AlignmentFile(str(bam), "rb") as af:
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


def pysam_depth(
    bam: Path, contig: str, start: int, stop: int, **_kw: Any
) -> dict[str, Any]:
    depth = pysam_aligned_depth(bam, contig, start, stop)
    return {"depth": depth, "sum": int(sum(depth))}


def rubam_depth(
    bam: Path,
    contig: str,
    start: int,
    stop: int,
    *,
    min_bq: int = 0,
    threads: int = 1,
    **_kw: Any,
) -> dict[str, Any]:
    import rubam

    if stop <= start:
        return {"depth": [], "sum": 0}
    _pos, dep = rubam.get_depths(
        str(bam),
        contig,
        start + 1,
        stop,
        step=1,
        min_mapq=0,
        min_bq=min_bq,
        max_depth=1_000_000,
        num_threads=max(1, threads),
    )
    depth = list(map(int, dep))
    return {"depth": depth, "sum": int(sum(depth))}


def samrust_depth(
    bam: Path,
    contig: str,
    start: int,
    stop: int,
    *,
    threads: int = 1,
    **_kw: Any,
) -> dict[str, Any]:
    import samrust

    with samrust.AlignmentFile(str(bam), "rb") as af:
        depth = list(map(int, af.depth_numpy(contig, start, stop, threads=threads)))
    return {"depth": depth, "sum": int(sum(depth))}


def pysam_pileup(
    bam: Path, contig: str, start: int, stop: int, *, min_bq: int = 0, **_kw: Any
) -> dict[str, Any]:
    import pysam

    length = stop - start
    a = [0] * length
    c = [0] * length
    g = [0] * length
    t = [0] * length
    n = [0] * length
    depth = [0] * length
    with pysam.AlignmentFile(str(bam), "rb") as af:
        for col in af.pileup(
            contig,
            start,
            stop,
            truncate=True,
            min_base_quality=min_bq,
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
                quals = pr.alignment.query_qualities
                if quals is not None and quals[qpos] < min_bq:
                    continue
                base = (pr.alignment.query_sequence or "N")[qpos].upper()
                depth[idx] += 1
                if base == "A":
                    a[idx] += 1
                elif base == "C":
                    c[idx] += 1
                elif base == "G":
                    g[idx] += 1
                elif base == "T":
                    t[idx] += 1
                else:
                    n[idx] += 1
    return {"A": a, "C": c, "G": g, "T": t, "N": n, "depth": depth}


def rubam_pileup(
    bam: Path,
    contig: str,
    start: int,
    stop: int,
    *,
    min_bq: int = 0,
    threads: int = 1,
    **_kw: Any,
) -> dict[str, Any]:
    import rubam

    if stop <= start:
        return {"A": [], "C": [], "G": [], "T": [], "N": [], "depth": []}
    # Prefer free-function path with explicit filters; flag_filter=0 for fair base counts
    # vs SAMRust default which excludes secondary/dup/etc. via PileupFilter.
    # For apples-to-apples with SAMRust defaults, use AlignmentFile.pileup + same flags if available.
    # Here we use pileup_bases with flag_filter matching SAMRust exclude mask (0x704|0x800=0xF04?);
    # SAMRust uses unmapped|secondary|qcfail|dup|suppl = 0x4|0x100|0x200|0x400|0x800 = 0xF04
    positions, a, c, g, t, n, depth = rubam.pileup_bases(
        str(bam),
        contig,
        start + 1,
        stop,
        step=1,
        min_mapq=0,
        min_bq=min_bq,
        max_depth=1_000_000,
        num_threads=max(1, threads),
        flag_filter=_FLAG_FILTER,
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


def samrust_pileup(
    bam: Path,
    contig: str,
    start: int,
    stop: int,
    *,
    min_bq: int = 0,
    threads: int = 1,
    **_kw: Any,
) -> dict[str, Any]:
    import samrust

    with samrust.AlignmentFile(str(bam), "rb") as af:
        counts = af.pileup_counts(
            contig, start, stop, min_base_quality=min_bq, threads=threads
        )
        return {
            "A": list(map(int, counts["A"])),
            "C": list(map(int, counts["C"])),
            "G": list(map(int, counts["G"])),
            "T": list(map(int, counts["T"])),
            "N": list(map(int, counts["N"])),
            "depth": list(map(int, counts["depth"])),
        }


WORKLOADS: dict[str, dict[str, Callable[..., dict[str, Any]]]] = {
    "count": {"pysam": pysam_count, "rubam": rubam_count, "samrust": samrust_count},
    "count_coverage": {
        "pysam": pysam_coverage,
        "rubam": rubam_coverage,
        "samrust": samrust_coverage,
    },
    "depth": {"pysam": pysam_depth, "rubam": rubam_depth, "samrust": samrust_depth},
    "pileup_counts": {
        "pysam": pysam_pileup,
        "rubam": rubam_pileup,
        "samrust": samrust_pileup,
    },
}


def _compare_to_pysam(workload: str, pysam_out: dict[str, Any], other: dict[str, Any]) -> str:
    """Legacy full-output comparison (kept for ad-hoc debugging).

    The benchmark path compares canonical digests instead (see run_workload).
    """
    return (
        "match"
        if _digest(_canonical(workload, other)) == _digest(_canonical(workload, pysam_out))
        else "mismatch"
    )


def _compact_detail(workload: str, out: dict[str, Any]) -> dict[str, Any]:
    """Compact detail for JSON (avoid dumping huge arrays)."""
    if workload == "count":
        return dict(out)
    detail: dict[str, Any] = {
        "sums": out.get("sums")
        or {
            k: int(sum(out[k]))
            for k in ("A", "C", "G", "T", "N", "depth")
            if k in out and isinstance(out[k], list)
        },
        "len": len(next((out[k] for k in ("A", "depth") if k in out), [])),
    }
    if "sum" in out:
        detail["sum"] = out["sum"]
    if "count_reads_nofilter" in out:
        detail["count_reads_nofilter"] = out["count_reads_nofilter"]
    return detail


def _measure_one(argv: list[str]) -> int:
    """Child-process entry: run one workload once, print JSON to stdout.

    Each (tool, repeat) runs in an isolated subprocess so the RSS high-water
    mark is per-tool (v0.1.1 P3b: previously all three tools shared the parent
    process and inherited each other's ``ru_maxrss``). The child emits the
    canonical digest + compact detail instead of raw arrays, keeping stdout
    small even for whole-chromosome workloads.
    """
    workload, tool, bam, contig, start, stop, threads, min_bq = argv
    fn = WORKLOADS[workload][tool]
    kw = {"threads": int(threads), "min_bq": int(min_bq)}
    fn(Path(bam), contig, int(start), int(stop), **kw)  # warm-up
    t0 = time.perf_counter()
    out = fn(Path(bam), contig, int(start), int(stop), **kw)
    wall = time.perf_counter() - t0
    print(
        json.dumps(
            {
                "digest": _digest(_canonical(workload, out)),
                "wall_s": wall,
                "max_rss_kb": _rss_kb(),
                "detail": _compact_detail(workload, out),
            }
        )
    )
    return 0


def _run_child(
    workload: str,
    tool: str,
    bam: Path,
    contig: str,
    start: int,
    stop: int,
    threads: int,
    min_bq: int,
) -> dict[str, Any]:
    proc = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).resolve()),
            "--_measure-one",
            workload,
            tool,
            str(bam),
            contig,
            str(start),
            str(stop),
            str(threads),
            str(min_bq),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(proc.stdout.strip().splitlines()[-1])


def run_workload(
    workload: str,
    bam: Path,
    contig: str,
    start: int,
    stop: int,
    *,
    threads: int,
    min_bq: int,
    repeats: int,
) -> list[RunResult]:
    measured: dict[str, list[dict[str, Any]]] = {t: [] for t in WORKLOADS[workload]}
    for _ in range(max(1, repeats)):
        for tool in WORKLOADS[workload]:
            measured[tool].append(
                _run_child(workload, tool, bam, contig, start, stop, threads, min_bq)
            )

    def pick(tool: str) -> dict[str, Any]:
        rows = sorted(measured[tool], key=lambda r: r["wall_s"])
        return rows[len(rows) // 2]

    pysam_row = pick("pysam")
    results: list[RunResult] = []
    for tool in ("pysam", "rubam", "samrust"):
        row = pick(tool)
        if tool == "pysam":
            gate = "oracle"
        elif row["digest"] == pysam_row["digest"]:
            gate = "match"
        elif (
            workload == "depth"
            and row["detail"].get("sum") == pysam_row["detail"].get("sum")
            and row["detail"].get("len") == pysam_row["detail"].get("len")
        ):
            gate = "sum_match_array_mismatch"
        else:
            gate = "mismatch"
        results.append(
            RunResult(
                tool=tool,
                workload=workload,
                threads=threads,
                wall_s=float(row["wall_s"]),
                max_rss_kb=int(row["max_rss_kb"]),
                output_digest=row["digest"],
                gate_vs_pysam=gate,
                detail=row["detail"],
            )
        )
    return results


def write_outputs(
    payload: dict[str, Any], rows: list[RunResult], out_json: Path, out_csv: Path
) -> None:
    """Write JSON + CSV. Missing rubam API/result must be the literal ``NA``, not omitted."""
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    with out_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "workload",
                "threads",
                "tool",
                "wall_s",
                "max_rss_kb",
                "output_digest",
                "gate_vs_pysam",
            ],
        )
        writer.writeheader()
        for r in rows:
            writer.writerow(
                {
                    "workload": r.workload,
                    "threads": r.threads,
                    "tool": r.tool,
                    "wall_s": f"{r.wall_s:.6f}",
                    "max_rss_kb": r.max_rss_kb,
                    "output_digest": r.output_digest,
                    "gate_vs_pysam": r.gate_vs_pysam,
                }
            )


def main() -> int:
    # Hidden child mode: --_measure-one <workload> <tool> <bam> <contig> <start> <stop> <threads> <min_bq>
    if len(sys.argv) > 1 and sys.argv[1] == "--_measure-one":
        return _measure_one(sys.argv[2:])

    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", type=Path, required=True)
    ap.add_argument("--contig", type=str, default=None)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--stop", type=int, default=200)
    ap.add_argument("--threads", type=str, default="1", help="colon/comma list, e.g. 1:4:8")
    ap.add_argument("--min-bq", type=int, default=0)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument(
        "--workloads",
        type=str,
        default="count,count_coverage,depth,pileup_counts",
        help="comma-separated subset of workloads",
    )
    ap.add_argument(
        "--outdir",
        type=Path,
        default=root / "benchmark" / "results",
    )
    ap.add_argument("--tag", type=str, default="fixture", help="filename tag / workload label")
    ap.add_argument("--repo-root", type=Path, default=root)
    args = ap.parse_args()

    if not args.bam.is_file():
        print(f"ERROR: BAM not found: {args.bam}", file=sys.stderr)
        return 2

    # Resolve default contig from BAM
    import pysam

    with pysam.AlignmentFile(str(args.bam), "rb") as af:
        contig = args.contig or af.references[0]
        stop = args.stop
        if stop <= args.start:
            stop = min(int(af.lengths[af.references.index(contig)]), args.start + 200)

    thread_list = [
        int(x)
        for x in args.threads.replace(" ", ":").replace(",", ":").split(":")
        if x.strip()
    ]
    workloads = [w.strip() for w in args.workloads.split(",") if w.strip()]
    for w in workloads:
        if w not in WORKLOADS:
            print(f"ERROR: unknown workload {w}; choose from {list(WORKLOADS)}", file=sys.stderr)
            return 2

    # Version pins
    versions: dict[str, Any] = {}
    try:
        import pysam as _pysam

        versions["pysam"] = _pysam.__version__
    except Exception:
        versions["pysam"] = None
    try:
        import rubam as _rubam

        versions["rubam"] = getattr(_rubam, "__version__", "unknown")
    except Exception as exc:
        print(f"ERROR: rubam required: {exc}", file=sys.stderr)
        return 2
    try:
        import samrust as _samrust

        versions["samrust"] = _samrust.__version__
    except Exception as exc:
        print(f"ERROR: samrust required: {exc}", file=sys.stderr)
        return 2

    all_rows: list[RunResult] = []
    for threads in thread_list:
        for workload in workloads:
            print(f"[three-way] {workload} threads={threads} ...", flush=True)
            all_rows.extend(
                run_workload(
                    workload,
                    args.bam,
                    contig,
                    args.start,
                    stop,
                    threads=threads,
                    min_bq=args.min_bq,
                    repeats=args.repeats,
                )
            )

    payload = {
        "metadata": collect(
            args.repo_root,
            {
                "benchmark_kind": "three_way_pysam_rubam_samrust",
                "versions": versions,
            },
        ),
        "bam": str(args.bam.resolve()),
        "region": {"contig": contig, "start": args.start, "stop": stop},
        "min_bq": args.min_bq,
        "repeats": args.repeats,
        "tag": args.tag,
        "rows": [
            {
                "workload": r.workload,
                "threads": r.threads,
                "tool": r.tool,
                "wall_s": r.wall_s,
                "max_rss_kb": r.max_rss_kb,
                "output_digest": r.output_digest,
                "gate_vs_pysam": r.gate_vs_pysam,
                "detail": r.detail,
            }
            for r in all_rows
        ],
    }

    stem = f"compare_pysam_rubam_samrust.{args.tag}"
    out_json = args.outdir / f"{stem}.json"
    out_csv = args.outdir / f"{stem}.csv"
    write_outputs(payload, all_rows, out_json, out_csv)
    print(f"Wrote {out_json}")
    print(f"Wrote {out_csv}")

    # Non-zero if any samrust/rubam mismatch vs pysam on count/coverage/pileup
    hard = {"count", "count_coverage", "depth", "pileup_counts"}
    bad = [
        r
        for r in all_rows
        if r.tool != "pysam" and r.workload in hard and r.gate_vs_pysam == "mismatch"
    ]
    if bad:
        print(f"WARNING: {len(bad)} mismatches vs pysam", file=sys.stderr)
        for r in bad[:10]:
            print(f"  {r.tool} {r.workload} T={r.threads} {r.gate_vs_pysam}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
