#!/usr/bin/env python3
"""M7 fungal BAM: random-region pileup serial==pysam and parallel==serial."""

from __future__ import annotations

import argparse
import json
import random
import time
from pathlib import Path

_FLAG_FILTER = 0x4 | 0x100 | 0x200 | 0x400 | 0x800


def pysam_counts(bam, contig, start, stop, min_bq=0):
    length = stop - start
    a = [0] * length
    c = [0] * length
    g = [0] * length
    t = [0] * length
    n = [0] * length
    depth = [0] * length
    for col in bam.pileup(
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
    return a, c, g, t, n, depth


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--bam",
        type=Path,
        default=Path(
            "/hpcfs/fhome/caizhh/18_WJX_work/02_Genetic_variation_diversity/"
            "05_align/Mt-35-15d-1/Mt-35-15d-1.markdup.bam"
        ),
    )
    ap.add_argument("--n-regions", type=int, default=20)
    ap.add_argument("--width", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-bq", type=int, default=0)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("benchmark/results/m7_fungal_pileup.json"),
    )
    args = ap.parse_args()

    import pysam
    import samrust

    if not args.bam.is_file():
        raise SystemExit(f"missing BAM: {args.bam}")

    py = pysam.AlignmentFile(str(args.bam), "rb")
    sr = samrust.AlignmentFile(str(args.bam), "rb")
    rng = random.Random(args.seed)
    refs = list(py.references)
    lens = list(py.lengths)

    mismatches = []
    parallel_mismatches = []
    regions = []
    t0 = time.perf_counter()
    for _ in range(args.n_regions):
        i = rng.randrange(len(refs))
        contig = refs[i]
        length = int(lens[i])
        if length < args.width:
            start, stop = 0, length
        else:
            start = rng.randrange(0, length - args.width)
            stop = start + args.width
        regions.append((contig, start, stop))

        pa, pc, pg, pt, pn, pd = pysam_counts(py, contig, start, stop, args.min_bq)
        serial = sr.pileup_counts(
            contig, start, stop, min_base_quality=args.min_bq, threads=1
        )
        for key, expected in (
            ("A", pa),
            ("C", pc),
            ("G", pg),
            ("T", pt),
            ("N", pn),
            ("depth", pd),
        ):
            got = list(serial[key])
            if got != expected:
                diffs = sum(1 for x, y in zip(got, expected) if x != y)
                mismatches.append(
                    {"region": [contig, start, stop], "channel": key, "diffs": diffs}
                )

        parallel = sr.pileup_counts(
            contig, start, stop, min_base_quality=args.min_bq, threads=args.threads
        )
        for key in ("A", "C", "G", "T", "N", "depth"):
            if list(parallel[key]) != list(serial[key]):
                parallel_mismatches.append(
                    {"region": [contig, start, stop], "channel": key}
                )

    elapsed = time.perf_counter() - t0
    payload = {
        "bam": str(args.bam),
        "n_regions": args.n_regions,
        "width": args.width,
        "threads": args.threads,
        "elapsed_s": elapsed,
        "pysam_mismatches": mismatches[:20],
        "parallel_mismatches": parallel_mismatches[:20],
        "gate_ok": not mismatches and not parallel_mismatches,
        "regions_sample": regions[:5],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({k: payload[k] for k in ("gate_ok", "elapsed_s", "n_regions")}, indent=2))
    print(f"wrote {args.out}")
    return 0 if payload["gate_ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
