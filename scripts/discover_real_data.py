#!/usr/bin/env python3
"""Discover Myceliophthora resequencing BAM/VCF paths (read-only).

Reads a YAML config (CLI-overridable). Writes:
  benchmark/real_data_manifest.tsv
  benchmark/configs/tiers.generated.yaml

Does not modify source sequencing data.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


MANIFEST_COLUMNS = [
    "sample_id",
    "bam",
    "index",
    "vcf",
    "condition",
    "timepoint",
    "temperature_c",
    "tier_hint",
]


@dataclass
class SampleRow:
    sample_id: str
    temperature_c: str
    timepoint: str
    condition: str
    bam: str
    index: str
    vcf: str
    tier_hint: str


def load_config(path: Path) -> dict[str, Any]:
    with path.open() as fh:
        data = yaml.safe_load(fh) or {}
    if not isinstance(data, dict):
        raise ValueError(f"config must be a mapping: {path}")
    return data


def find_bam(analysis_dir: Path, sample_id: str) -> tuple[Path | None, Path | None]:
    """Prefer primary markdup BAM under 05_align; never pick evidence/pytest trees."""
    align_root = analysis_dir / "05_align" / sample_id
    candidates = [
        align_root / f"{sample_id}.markdup.bam",
        align_root / f"{sample_id}.bam",
    ]
    for bam in candidates:
        if bam.is_file():
            for idx in (Path(str(bam) + ".bai"), Path(str(bam) + ".csi")):
                if idx.is_file():
                    return bam, idx
            return bam, None
    return None, None


def find_vcf(analysis_dir: Path, sample_id: str) -> Path | None:
    var_root = analysis_dir / "06_variants_call" / sample_id
    for name in (
        f"{sample_id}.freebayes.vcf.gz",
        f"{sample_id}.vcf.gz",
        f"{sample_id}.bcf",
    ):
        path = var_root / name
        if path.is_file():
            return path
    return None


def assign_tier_hint(temperature: str, timepoint: str, sample_id: str) -> str:
    # Heuristic labels only — tiers.generated.yaml selects concrete Tier1/2/3 sets.
    if sample_id.endswith("-1") and timepoint in {"3", "3d"}:
        return "tier2_candidate"
    if temperature in {"35", "45", "50", "52"}:
        return "tier3_pool"
    return "discovered"


def discover(config: dict[str, Any]) -> list[SampleRow]:
    meta_path = Path(config["sample_metadata"])
    analysis_dir = Path(config["analysis_dir"])
    if not meta_path.is_file():
        raise FileNotFoundError(f"sample_metadata not found: {meta_path}")
    if not analysis_dir.is_dir():
        raise FileNotFoundError(f"analysis_dir not found: {analysis_dir}")

    rows: list[SampleRow] = []
    with meta_path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for rec in reader:
            sample_id = (rec.get("sampleID_canonical") or rec.get("sampleID_raw") or "").strip()
            if not sample_id:
                continue
            bam, index = find_bam(analysis_dir, sample_id)
            if bam is None:
                continue
            vcf = find_vcf(analysis_dir, sample_id)
            temp = str(rec.get("culture_temperature_c", "")).strip()
            days = str(rec.get("passage_duration_days", "")).strip()
            condition = f"T{temp}C_{days}d" if temp or days else ""
            rows.append(
                SampleRow(
                    sample_id=sample_id,
                    temperature_c=temp,
                    timepoint=days,
                    condition=condition,
                    bam=str(bam.resolve()),
                    index=str(index.resolve()) if index else "",
                    vcf=str(vcf.resolve()) if vcf else "",
                    tier_hint=assign_tier_hint(temp, days, sample_id),
                )
            )
    return rows


def write_manifest(path: Path, rows: list[SampleRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow(
                {
                    "sample_id": r.sample_id,
                    "bam": r.bam,
                    "index": r.index,
                    "vcf": r.vcf,
                    "condition": r.condition,
                    "timepoint": r.timepoint,
                    "temperature_c": r.temperature_c,
                    "tier_hint": r.tier_hint,
                }
            )


def choose_tiers(rows: list[SampleRow]) -> dict[str, Any]:
    """Build Tier1/2/3 selectors without hard-coding sample IDs when possible."""
    indexed = [r for r in rows if r.index]
    if not indexed:
        return {"error": "no indexed BAMs discovered"}

    # Prefer mid coverage heuristic later; for now pick first indexed as Tier2.
    tier2 = indexed[0]

    # Diversify Tier3 by temperature then timepoint.
    by_temp: dict[str, list[SampleRow]] = {}
    for r in indexed:
        by_temp.setdefault(r.temperature_c or "NA", []).append(r)
    tier3: list[SampleRow] = []
    for _temp, group in sorted(by_temp.items()):
        # pick earliest and latest timepoint in group when available
        group_sorted = sorted(group, key=lambda x: (x.timepoint, x.sample_id))
        for candidate in (group_sorted[0], group_sorted[-1]):
            if candidate.sample_id not in {x.sample_id for x in tier3}:
                tier3.append(candidate)
        if len(tier3) >= 6:
            break

    # Contig for Tier1: filled by baselines after reading BAM header; placeholder here.
    return {
        "tier1": {
            "sample_id": tier2.sample_id,
            "bam": tier2.bam,
            "index": tier2.index,
            "region_contig": null_placeholder(),
            "region_start": 0,
            "region_length_bp": 1_000_000,
            "notes": "contig filled by run_baselines / metadata collector from BAM header",
        },
        "tier2": {
            "sample_id": tier2.sample_id,
            "bam": tier2.bam,
            "index": tier2.index,
            "notes": "full-genome sample for correctness + scaling",
        },
        "tier3": [
            {
                "sample_id": r.sample_id,
                "bam": r.bam,
                "index": r.index,
                "temperature_c": r.temperature_c,
                "timepoint": r.timepoint,
                "condition": r.condition,
            }
            for r in tier3
        ],
    }


def null_placeholder() -> None:
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="YAML config with analysis_dir / sample_metadata (and optional path overrides)",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("benchmark/real_data_manifest.tsv"),
        help="Output manifest TSV (gitignored by default)",
    )
    parser.add_argument(
        "--tiers-out",
        type=Path,
        default=Path("benchmark/configs/tiers.generated.yaml"),
        help="Generated tier selection YAML (gitignored)",
    )
    parser.add_argument(
        "--analysis-dir",
        type=Path,
        default=None,
        help="Override config analysis_dir",
    )
    parser.add_argument(
        "--sample-metadata",
        type=Path,
        default=None,
        help="Override config sample_metadata",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    if args.analysis_dir:
        config["analysis_dir"] = str(args.analysis_dir)
    if args.sample_metadata:
        config["sample_metadata"] = str(args.sample_metadata)

    for key in ("analysis_dir", "sample_metadata"):
        if key not in config:
            print(f"ERROR: config missing '{key}'", file=sys.stderr)
            return 2

    try:
        rows = discover(config)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        print("STOP real-data discovery. Prepare coordinate-sorted indexed BAM first.", file=sys.stderr)
        return 3

    if not rows:
        print(
            "ERROR: no BAM files discovered for metadata samples under analysis_dir.",
            file=sys.stderr,
        )
        print(
            "Expected e.g. <analysis_dir>/05_align/<sample>/<sample>.markdup.bam (+ .bai).",
            file=sys.stderr,
        )
        return 4

    missing_index = [r.sample_id for r in rows if not r.index]
    write_manifest(args.manifest, rows)
    tiers = choose_tiers(rows)
    args.tiers_out.parent.mkdir(parents=True, exist_ok=True)
    with args.tiers_out.open("w") as fh:
        yaml.safe_dump(tiers, fh, sort_keys=False)

    print(f"Discovered {len(rows)} samples with BAM")
    print(f"  manifest: {args.manifest}")
    print(f"  tiers:    {args.tiers_out}")
    if missing_index:
        print(f"WARNING: {len(missing_index)} BAMs missing index: {', '.join(missing_index[:8])}...")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
