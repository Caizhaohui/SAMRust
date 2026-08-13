#!/usr/bin/env python3
"""M2 parity: compare samrust dump-records vs pysam (0 mismatches gate).

Usage:
  python scripts/parity_bam_core.py --bam tests/fixtures/small.bam
  python scripts/parity_bam_core.py --bam /path/to/real.bam --limit 100000
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pysam

SELECTED_TAGS = ("NM", "RG", "MD", "AS", "XS")


def find_samrust_bin(repo: Path) -> Path:
    candidates = [
        repo / "target" / "debug" / "samrust",
        repo / "target" / "release" / "samrust",
    ]
    for path in candidates:
        if path.is_file():
            return path
    # build debug binary
    cmd = [
        "cargo",
        "--config",
        'source.crates-io.replace-with="ustc"',
        "--config",
        'source.ustc.registry="sparse+https://mirrors.ustc.edu.cn/crates.io-index/"',
        "build",
        "-p",
        "samrust-cli",
        "-q",
    ]
    subprocess.run(cmd, cwd=repo, check=True)
    path = repo / "target" / "debug" / "samrust"
    if not path.is_file():
        raise FileNotFoundError("samrust binary not found after build")
    return path


def load_samrust(bin_path: Path, bam: Path, limit: int) -> tuple[dict, list[dict]]:
    cmd = [str(bin_path), "dump-records", "--bam", str(bam)]
    if limit:
        cmd.extend(["--limit", str(limit)])
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    header = None
    records: list[dict] = []
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        obj = json.loads(line)
        if obj.get("type") == "header":
            header = obj
        elif obj.get("type") == "record":
            records.append(obj)
    if header is None:
        raise RuntimeError("samrust dump missing header line")
    return header, records


def pysam_records(bam: Path, limit: int) -> tuple[dict, list[dict]]:
    with pysam.AlignmentFile(str(bam), "rb") as af:
        header = {
            "type": "header",
            "nreferences": af.nreferences,
            "references": list(af.references),
            "lengths": [int(x) for x in af.lengths],
        }
        records: list[dict] = []
        for i, r in enumerate(af.fetch(until_eof=True)):
            if limit and i >= limit:
                break
            tags = {}
            for key in SELECTED_TAGS:
                if r.has_tag(key):
                    tags[key] = r.get_tag(key)
            records.append(
                {
                    "type": "record",
                    "qname": r.query_name,
                    "flag": int(r.flag),
                    "reference_id": int(r.reference_id),
                    "reference_start": int(r.reference_start),
                    "mapping_quality": int(r.mapping_quality),
                    "cigar": r.cigarstring,
                    "query_length": int(r.query_length),
                    "tags": tags,
                }
            )
    return header, records


def normalize_tag(value):
    # pysam may return numpy ints
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    return value


def compare(samrust_h, samrust_recs, pysam_h, pysam_recs) -> list[str]:
    mismatches: list[str] = []
    if samrust_h["nreferences"] != pysam_h["nreferences"]:
        mismatches.append(
            f"header.nreferences samrust={samrust_h['nreferences']} pysam={pysam_h['nreferences']}"
        )
    if samrust_h["references"] != pysam_h["references"]:
        mismatches.append("header.references differ")
    if samrust_h["lengths"] != pysam_h["lengths"]:
        mismatches.append("header.lengths differ")

    if len(samrust_recs) != len(pysam_recs):
        mismatches.append(
            f"record_count samrust={len(samrust_recs)} pysam={len(pysam_recs)}"
        )
        return mismatches

    fields = (
        "qname",
        "flag",
        "reference_id",
        "reference_start",
        "mapping_quality",
        "cigar",
        "query_length",
    )
    for i, (a, b) in enumerate(zip(samrust_recs, pysam_recs)):
        for field in fields:
            if a.get(field) != b.get(field):
                mismatches.append(
                    f"record[{i}].{field} samrust={a.get(field)!r} pysam={b.get(field)!r} qname={b.get('qname')}"
                )
        # selected tags present on either side
        keys = set(a.get("tags", {})) | set(b.get("tags", {}))
        for key in keys:
            av = a.get("tags", {}).get(key)
            bv = normalize_tag(b.get("tags", {}).get(key))
            if av != bv:
                mismatches.append(
                    f"record[{i}].tags.{key} samrust={av!r} pysam={bv!r} qname={b.get('qname')}"
                )
        if len(mismatches) >= 50:
            mismatches.append("... truncated after 50 mismatches")
            break
    return mismatches


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    args = parser.parse_args()

    if not args.bam.is_file():
        print(f"ERROR: BAM not found: {args.bam}", file=sys.stderr)
        return 2
    if shutil.which("cargo") is None and not (args.repo_root / "target/debug/samrust").is_file():
        print("ERROR: cargo/samrust binary unavailable", file=sys.stderr)
        return 2

    bin_path = find_samrust_bin(args.repo_root)
    sr_h, sr_recs = load_samrust(bin_path, args.bam, args.limit)
    py_h, py_recs = pysam_records(args.bam, args.limit)
    mismatches = compare(sr_h, sr_recs, py_h, py_recs)

    if mismatches:
        print(f"FAIL: {len(mismatches)} mismatches", file=sys.stderr)
        for line in mismatches:
            print(line, file=sys.stderr)
        return 1

    print(
        f"OK: 0 mismatches ({len(sr_recs)} records) bam={args.bam} limit={args.limit or 'all'}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
