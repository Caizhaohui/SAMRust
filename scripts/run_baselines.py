#!/usr/bin/env python3
"""Orchestrate M1 baseline generation (synthetic required; real Tier-1 optional)."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import yaml


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def fill_tier1_contig(tiers_path: Path) -> dict:
    if not tiers_path.is_file():
        return {}
    with tiers_path.open() as fh:
        tiers = yaml.safe_load(fh) or {}
    tier1 = tiers.get("tier1") or {}
    bam = tier1.get("bam")
    if not bam:
        return tiers
    import pysam

    with pysam.AlignmentFile(bam, "rb") as af:
        contig = af.references[0]
        length = int(af.lengths[0])
    region_len = int(tier1.get("region_length_bp") or 1_000_000)
    tier1["region_contig"] = contig
    tier1["region_start"] = 0
    tier1["region_stop"] = min(length, region_len)
    tiers["tier1"] = tier1
    with tiers_path.open("w") as fh:
        yaml.safe_dump(tiers, fh, sort_keys=False)
    return tiers


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=root)
    parser.add_argument(
        "--fixture-bam",
        type=Path,
        default=None,
        help="Synthetic BAM (default: tests/fixtures/small.bam)",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Optional real-data YAML; if set, also run discovery + Tier-1 baselines",
    )
    parser.add_argument(
        "--skip-real",
        action="store_true",
        help="Only generate synthetic baselines",
    )
    args = parser.parse_args()
    repo = args.repo_root
    scripts = repo / "scripts"
    fixture_bam = args.fixture_bam or (repo / "tests" / "fixtures" / "small.bam")
    results = repo / "benchmark" / "results"
    results.mkdir(parents=True, exist_ok=True)

    # Ensure fixtures exist
    if not fixture_bam.is_file():
        run([sys.executable, str(scripts / "prepare_fixture.py"), "--out-dir", str(fixture_bam.parent)])

    run([sys.executable, str(scripts / "collect_benchmark_metadata.py"), "--repo-root", str(repo), "--output", str(results / "host_metadata.json")])
    run(
        [
            sys.executable,
            str(scripts / "compare_pysam.py"),
            "--bam",
            str(fixture_bam),
            "--contig",
            "chr1",
            "--start",
            "0",
            "--stop",
            "200",
            "--output",
            str(results / "pysam_baseline.json"),
            "--repo-root",
            str(repo),
        ]
    )
    run(
        [
            sys.executable,
            str(scripts / "compare_samtools.py"),
            "--bam",
            str(fixture_bam),
            "--contig",
            "chr1",
            "--start",
            "0",
            "--stop",
            "200",
            "--output",
            str(results / "samtools_baseline.json"),
            "--repo-root",
            str(repo),
        ]
    )

    if args.skip_real or args.config is None:
        print("Synthetic baselines ready (real-data skipped).")
        return 0

    config = args.config
    if not config.is_file():
        print(f"WARNING: config not found ({config}); stopping real-data portion.", file=sys.stderr)
        return 0

    manifest = repo / "benchmark" / "real_data_manifest.tsv"
    tiers_out = repo / "benchmark" / "configs" / "tiers.generated.yaml"
    try:
        run(
            [
                sys.executable,
                str(scripts / "discover_real_data.py"),
                "--config",
                str(config),
                "--manifest",
                str(manifest),
                "--tiers-out",
                str(tiers_out),
            ]
        )
    except subprocess.CalledProcessError as exc:
        print(f"STOP real-data discovery (exit {exc.returncode}). Synthetic baselines kept.", file=sys.stderr)
        return 0

    tiers = fill_tier1_contig(tiers_out)
    tier1 = tiers.get("tier1") or {}
    bam = tier1.get("bam")
    contig = tier1.get("region_contig")
    start = int(tier1.get("region_start") or 0)
    stop = int(tier1.get("region_stop") or 0)
    if not bam or not contig or stop <= start:
        print("WARNING: Tier-1 incomplete; skip real baselines.", file=sys.stderr)
        return 0

    # Keep real Tier-1 windows modest for M1 (<1 min goal when possible)
    stop = min(stop, start + 50_000)
    run(
        [
            sys.executable,
            str(scripts / "compare_pysam.py"),
            "--bam",
            str(bam),
            "--contig",
            str(contig),
            "--start",
            str(start),
            "--stop",
            str(stop),
            "--max-records",
            "1000",
            "--output",
            str(results / "pysam_baseline_tier1.json"),
            "--repo-root",
            str(repo),
        ]
    )
    run(
        [
            sys.executable,
            str(scripts / "compare_samtools.py"),
            "--bam",
            str(bam),
            "--contig",
            str(contig),
            "--start",
            str(start),
            "--stop",
            str(stop),
            "--output",
            str(results / "samtools_baseline_tier1.json"),
            "--repo-root",
            str(repo),
        ]
    )

    summary = {
        "synthetic": {
            "pysam": str(results / "pysam_baseline.json"),
            "samtools": str(results / "samtools_baseline.json"),
        },
        "tier1": {
            "pysam": str(results / "pysam_baseline_tier1.json"),
            "samtools": str(results / "samtools_baseline_tier1.json"),
            "region": {"contig": contig, "start": start, "stop": stop, "bam": bam},
        },
    }
    (results / "baseline_index.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("Wrote real Tier-1 baselines + baseline_index.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
