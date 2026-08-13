#!/usr/bin/env python3
"""M8 fungal domain validation: ALT_COUNT>=10 gate + thread scaling.

Reads real BAM + candidate BED (read-only). Writes JSON under benchmark/results/.
Does not modify original resequencing data.
"""

from __future__ import annotations

import argparse
import json
import resource
import time
from collections import defaultdict
from pathlib import Path

# Match samrust PileupFilter::default()
_FLAG_FILTER = 0x4 | 0x100 | 0x200 | 0x400 | 0x800


def load_sites_bed(path: Path, limit: int | None = None):
    sites = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            chrom, start, _stop, allele = cols[0], int(cols[1]), cols[2], cols[3]
            ref, alt = allele.split(">", 1)
            sites.append((chrom, start, ref, alt))  # start 0-based
            if limit is not None and len(sites) >= limit:
                break
    return sites


def pysam_recount(bam_path: Path, sample: str, sites, min_bq: int = 0):
    import pysam

    bam = pysam.AlignmentFile(str(bam_path), "rb")
    # Group by chrom for fewer seek patterns
    by_chrom: dict[str, list[tuple[int, str, str, int]]] = defaultdict(list)
    for i, (chrom, pos0, ref, alt) in enumerate(sites):
        by_chrom[chrom].append((pos0, ref, alt, i))

    rows = [None] * len(sites)
    for chrom, items in by_chrom.items():
        items.sort(key=lambda x: x[0])
        for pos0, ref, alt, idx in items:
            a = c = g = t = n = depth = 0
            for col in bam.pileup(
                chrom,
                pos0,
                pos0 + 1,
                truncate=True,
                min_base_quality=min_bq,
                flag_filter=_FLAG_FILTER,
                stepper="all",
            ):
                if col.reference_pos != pos0:
                    continue
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
                    depth += 1
                    if base == "A":
                        a += 1
                    elif base == "C":
                        c += 1
                    elif base == "G":
                        g += 1
                    elif base == "T":
                        t += 1
                    else:
                        n += 1
            alt_u = alt.upper()
            alt_count = {"A": a, "C": c, "G": g, "T": t}.get(alt_u, n if len(alt) == 1 else 0)
            rows[idx] = {
                "sample": sample,
                "chrom": chrom,
                "pos": pos0 + 1,
                "ref": ref,
                "alt": alt,
                "A": a,
                "C": c,
                "G": g,
                "T": t,
                "N": n,
                "DP": depth,
                "ALT_COUNT": alt_count,
                "AF": (alt_count / depth) if depth else 0.0,
            }
    bam.close()
    return rows


def load_recount_tsv(path: Path):
    rows = []
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            rec = dict(zip(header, cols))
            rows.append(
                {
                    "sample": rec["sample"],
                    "chrom": rec["chrom"],
                    "pos": int(rec["pos"]),
                    "ref": rec["ref"],
                    "alt": rec.get("alt", ""),
                    "A": int(rec["A"]),
                    "C": int(rec["C"]),
                    "G": int(rec["G"]),
                    "T": int(rec["T"]),
                    "N": int(rec["N"]),
                    "DP": int(rec["DP"]),
                    "ALT_COUNT": int(rec["ALT_COUNT"]),
                    "AF": float(rec["AF"]),
                }
            )
    return rows


def samrust_recount(bam_path: Path, sample: str, sites_path: Path, threads: int, min_bq: int = 0):
    import subprocess
    import tempfile

    out = tempfile.NamedTemporaryFile(prefix="samrust_recount_", suffix=".tsv", delete=False)
    out_path = Path(out.name)
    out.close()
    cmd = [
        "samrust",
        "recount",
        "--bam",
        str(bam_path),
        "--sites",
        str(sites_path),
        "--sample",
        sample,
        "--threads",
        str(threads),
        "--min-base-quality",
        str(min_bq),
        "--output",
        str(out_path),
    ]
    t0 = time.perf_counter()
    rss0 = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    rss1 = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    rows = load_recount_tsv(out_path)
    out_path.unlink(missing_ok=True)
    return rows, {
        "elapsed_s": elapsed,
        "ru_maxrss_kb": max(rss0, rss1),
        "stderr": proc.stderr.strip(),
    }


def alt_ge_set(rows, threshold: int = 10):
    return {
        (r["chrom"], r["pos"], r["ref"], r["alt"])
        for r in rows
        if r["ALT_COUNT"] >= threshold
    }


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
    ap.add_argument(
        "--sites",
        type=Path,
        default=Path(
            "/hpcfs/fhome/caizhh/18_WJX_work/02_Genetic_variation_diversity/"
            "07_candidates_union/Mt-35.candidates.bed"
        ),
    )
    ap.add_argument("--sample", default="Mt-35-15d-1")
    ap.add_argument("--limit", type=int, default=None, help="optional site limit for smoke runs")
    ap.add_argument("--min-bq", type=int, default=0)
    ap.add_argument("--threshold", type=int, default=10)
    ap.add_argument(
        "--threads",
        default="1,2,4,8,16,32",
        help="comma-separated thread counts for scaling",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("benchmark/results/m8_fungal_validation.json"),
    )
    ap.add_argument("--skip-scaling", action="store_true")
    ap.add_argument(
        "--samrust-tsv",
        type=Path,
        default=None,
        help="Reuse an existing samrust recount TSV for the ALT>=threshold gate (skip threads=1 recount)",
    )
    ap.add_argument(
        "--scaling-outdir",
        type=Path,
        default=Path("benchmark/results"),
        help="Directory for per-thread recount TSV outputs during scaling",
    )
    args = ap.parse_args()

    if not args.bam.is_file():
        raise SystemExit(f"missing BAM: {args.bam}")
    if not args.sites.is_file():
        raise SystemExit(f"missing sites: {args.sites}")

    sites_path = args.sites
    work_sites = sites_path
    if args.limit is not None:
        import tempfile

        subset = load_sites_bed(sites_path, limit=args.limit)
        tmp = tempfile.NamedTemporaryFile(
            prefix="m8_sites_", suffix=".bed", delete=False, mode="w"
        )
        for chrom, pos0, ref, alt in subset:
            tmp.write(f"{chrom}\t{pos0}\t{pos0 + 1}\t{ref}>{alt}\n")
        tmp.close()
        work_sites = Path(tmp.name)

    print(f"Loading pysam reference recount for {work_sites} ...", flush=True)
    sites = load_sites_bed(work_sites)
    t0 = time.perf_counter()
    py_rows = pysam_recount(args.bam, args.sample, sites, min_bq=args.min_bq)
    py_elapsed = time.perf_counter() - t0
    py_set = alt_ge_set(py_rows, args.threshold)
    print(
        f"pysam: {len(sites)} sites, ALT>={args.threshold}: {len(py_set)} in {py_elapsed:.1f}s",
        flush=True,
    )

    if args.samrust_tsv is not None:
        print(f"Reusing SAMRust TSV: {args.samrust_tsv}", flush=True)
        sr_rows = load_recount_tsv(args.samrust_tsv)
        meta1 = {
            "elapsed_s": None,
            "ru_maxrss_kb": None,
            "stderr": f"reused {args.samrust_tsv}",
        }
    else:
        print("SAMRust recount threads=1 ...", flush=True)
        sr_rows, meta1 = samrust_recount(
            args.bam, args.sample, work_sites, threads=1, min_bq=args.min_bq
        )
    sr_set = alt_ge_set(sr_rows, args.threshold)
    gate_ok = sr_set == py_set
    only_sr = sorted(sr_set - py_set)[:20]
    only_py = sorted(py_set - sr_set)[:20]
    print(
        f"samrust: ALT>={args.threshold}: {len(sr_set)}; gate_ok={gate_ok}; "
        f"elapsed={meta1.get('elapsed_s')}",
        flush=True,
    )
    if not gate_ok:
        print(f"only_samrust sample: {only_sr}", flush=True)
        print(f"only_pysam sample: {only_py}", flush=True)

    scaling = []
    if not args.skip_scaling and gate_ok:
        args.scaling_outdir.mkdir(parents=True, exist_ok=True)
        for t in [int(x) for x in args.threads.replace(":", ",").replace(" ", ",").split(",") if x.strip()]:
            print(f"scaling threads={t} ...", flush=True)
            out_tsv = args.scaling_outdir / f"m8_recount_t{t}.tsv"
            import subprocess

            cmd = [
                "samrust",
                "recount",
                "--bam",
                str(args.bam),
                "--sites",
                str(work_sites),
                "--sample",
                args.sample,
                "--threads",
                str(t),
                "--min-base-quality",
                str(args.min_bq),
                "--output",
                str(out_tsv),
            ]
            t0 = time.perf_counter()
            rss0 = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
            proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
            elapsed = time.perf_counter() - t0
            rss1 = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
            meta = {
                "threads": t,
                "elapsed_s": elapsed,
                "ru_maxrss_kb": max(rss0, rss1),
                "output": str(out_tsv),
                "stderr": proc.stderr.strip(),
            }
            scaling.append(meta)
            print(f"  {elapsed:.2f}s rss_kb={meta['ru_maxrss_kb']} -> {out_tsv}", flush=True)

    payload = {
        "bam": str(args.bam),
        "sites": str(args.sites),
        "sample": args.sample,
        "n_sites": len(sites),
        "threshold": args.threshold,
        "pysam_alt_ge": len(py_set),
        "samrust_alt_ge": len(sr_set),
        "gate_ok": gate_ok,
        "pysam_elapsed_s": py_elapsed,
        "samrust_threads1": meta1,
        "scaling": scaling,
        "only_samrust_sample": only_sr,
        "only_pysam_sample": only_py,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {args.out}", flush=True)

    if args.limit is not None and work_sites != sites_path:
        work_sites.unlink(missing_ok=True)

    return 0 if gate_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
