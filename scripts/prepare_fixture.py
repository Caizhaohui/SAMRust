#!/usr/bin/env python3
"""Build Tier-0 synthetic BAM/VCF fixtures for SAMRust CI and parity tests.

Requires: pysam, and samtools/bcftools (bgzip/tabix) on PATH.
Does not touch real resequencing data.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pysam


CHR1_LEN = 1000
CHR2_LEN = 500


def _make_seq(length: int, motif: str) -> str:
    reps = (length // len(motif)) + 1
    return (motif * reps)[:length]


def write_fasta(path: Path) -> None:
    chr1 = _make_seq(CHR1_LEN, "ACGT")
    chr2 = _make_seq(CHR2_LEN, "GCTA")
    path.write_text(f">chr1\n{chr1}\n>chr2\n{chr2}\n")
    pysam.faidx(str(path))


def _aligned_segment(
    header: pysam.AlignmentHeader,
    *,
    qname: str,
    flag: int,
    ref_name: str | None,
    ref_start: int,
    mapq: int,
    cigar: str | None,
    seq: str,
    qual: str,
    next_ref_name: str | None = None,
    next_ref_start: int = -1,
    tlen: int = 0,
    tags: list[tuple] | None = None,
) -> pysam.AlignedSegment:
    a = pysam.AlignedSegment(header)
    a.query_name = qname
    a.flag = flag
    if ref_name is None:
        a.reference_id = -1
        a.reference_start = -1
    else:
        a.reference_name = ref_name
        a.reference_start = ref_start
    a.mapping_quality = mapq
    if cigar is None:
        a.cigar = None
    else:
        a.cigarstring = cigar
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array(qual)
    if next_ref_name is not None:
        a.next_reference_name = next_ref_name
        a.next_reference_start = next_ref_start
    a.template_length = tlen
    if tags:
        for tag in tags:
            a.set_tag(*tag)
    return a


def build_bam(bam_path: Path) -> None:
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": "chr1", "LN": CHR1_LEN},
            {"SN": "chr2", "LN": CHR2_LEN},
        ],
        "RG": [{"ID": "synth", "SM": "small", "PL": "ILLUMINA"}],
    }
    h = pysam.AlignmentHeader.from_dict(header)
    records: list[pysam.AlignedSegment] = [
        _aligned_segment(
            h,
            qname="pair1",
            flag=99,
            ref_name="chr1",
            ref_start=10,
            mapq=60,
            cigar="50M",
            seq="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC",
            qual="I" * 50,
            next_ref_name="chr1",
            next_ref_start=80,
            tlen=120,
            tags=[("RG", "synth"), ("NM", 0, "i")],
        ),
        _aligned_segment(
            h,
            qname="pair1",
            flag=147,
            ref_name="chr1",
            ref_start=80,
            mapq=60,
            cigar="50M",
            seq="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC",
            qual="I" * 50,
            next_ref_name="chr1",
            next_ref_start=10,
            tlen=-120,
            tags=[("RG", "synth"), ("NM", 0, "i")],
        ),
        _aligned_segment(
            h,
            qname="indel1",
            flag=0,
            ref_name="chr1",
            ref_start=200,
            mapq=40,
            cigar="10M2I8M3D10M",
            seq="ACGTACGTACXXACGTACGTACGTACGTAC",
            qual="I" * 30,
            tags=[("NM", 5, "i")],
        ),
        _aligned_segment(
            h,
            qname="clip1",
            flag=0,
            ref_name="chr1",
            ref_start=300,
            mapq=30,
            cigar="5S20M5H",
            seq="NNNNNACGTACGTACGTACGTACGT",
            qual="#####" + ("I" * 20),
            tags=[("NM", 0, "i")],
        ),
        _aligned_segment(
            h,
            qname="splice1",
            flag=0,
            ref_name="chr1",
            ref_start=400,
            mapq=50,
            cigar="15M50N15M",
            seq="ACGTACGTACGTACGACGTACGTACGTACG",
            qual="I" * 30,
        ),
        _aligned_segment(
            h,
            qname="eqx1",
            flag=0,
            ref_name="chr1",
            ref_start=500,
            mapq=55,
            cigar="20=2X18=",
            seq="ACGTACGTACGTACGTACGTNNACGTACGTACGTACGTAC",
            qual="I" * 40,
        ),
        _aligned_segment(
            h,
            qname="dup1",
            flag=1024,
            ref_name="chr1",
            ref_start=600,
            mapq=20,
            cigar="25M",
            seq="ACGTACGTACGTACGTACGTACGTA",
            qual="I" * 25,
        ),
        _aligned_segment(
            h,
            qname="sec1",
            flag=256,
            ref_name="chr1",
            ref_start=610,
            mapq=10,
            cigar="25M",
            seq="ACGTACGTACGTACGTACGTACGTA",
            qual="I" * 25,
        ),
        _aligned_segment(
            h,
            qname="sup1",
            flag=2048,
            ref_name="chr1",
            ref_start=620,
            mapq=15,
            cigar="25M",
            seq="ACGTACGTACGTACGTACGTACGTA",
            qual="I" * 25,
        ),
        _aligned_segment(
            h,
            qname="qcfail1",
            flag=512,
            ref_name="chr1",
            ref_start=630,
            mapq=5,
            cigar="25M",
            seq="ACGTACGTACGTACGTACGTACGTA",
            qual="!" * 25,
        ),
        _aligned_segment(
            h,
            qname="unmap1",
            flag=4,
            ref_name=None,
            ref_start=-1,
            mapq=0,
            cigar=None,
            seq="ACGTACGTAC",
            qual="!!!!!!!!!!",
        ),
        # Placed unmapped mate: FLAG 0x4 with RNAME/POS set (pysam count/fetch include these).
        _aligned_segment(
            h,
            qname="unmap_placed",
            flag=133,
            ref_name="chr1",
            ref_start=50,
            mapq=0,
            cigar=None,
            seq="ACGTACGTAC",
            qual="!!!!!!!!!!",
            next_ref_name="chr1",
            next_ref_start=10,
        ),
        _aligned_segment(
            h,
            qname="mapq0",
            flag=0,
            ref_name="chr2",
            ref_start=10,
            mapq=0,
            cigar="30M",
            seq="GCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
            qual="I" * 30,
        ),
        _aligned_segment(
            h,
            qname="mapq255",
            flag=0,
            ref_name="chr2",
            ref_start=50,
            mapq=255,
            cigar="30M",
            seq="GCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
            qual="I" * 30,
        ),
    ]

    mapped = [r for r in records if not r.is_unmapped]
    unmapped = [r for r in records if r.is_unmapped]
    mapped.sort(key=lambda r: (r.reference_id, r.reference_start, r.query_name, r.flag))
    ordered = mapped + unmapped

    tmp = bam_path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(tmp), "wb", header=header) as out:
        for r in ordered:
            out.write(r)

    pysam.sort("-o", str(bam_path), str(tmp))
    tmp.unlink(missing_ok=True)
    pysam.index(str(bam_path))


def build_vcf(vcf_path: Path) -> None:
    plain = Path(str(vcf_path).removesuffix(".gz"))
    plain.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=500>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t15\t.\tA\tG\t60\tPASS\tDP=20;AF=0.5\tGT:DP:AD\t0/1:20:10,10
chr1\t250\tins1\tA\tATT\t40\tPASS\tDP=12;AF=0.25\tGT:DP:AD\t0/1:12:9,3
chr1\t260\tdel1\tACG\tA\t35\tPASS\tDP=15;AF=0.2\tGT:DP:AD\t0/1:15:12,3
chr1\t510\tmulti1\tC\tG,T\t80\tPASS\tDP=30;AF=0.3,0.2\tGT:DP:AD\t1/2:30:5,15,10
chr2\t20\t.\tG\tA\t10\tLowQual\tDP=5;AF=0.1\tGT:DP:AD\t0/1:5:4,1
"""
    )
    pysam.tabix_compress(str(plain), str(vcf_path), force=True)
    pysam.tabix_index(str(vcf_path), preset="vcf", force=True)
    # Keep uncompressed VCF for M9 sequential-read tests (`plain` is small.vcf).


def build_bcf(vcf_gz: Path, bcf_path: Path) -> None:
    bcftools = shutil.which("bcftools")
    if not bcftools:
        print("WARN: bcftools not on PATH; skip small.bcf", file=sys.stderr)
        return
    subprocess.run(
        [bcftools, "view", "-Ob", "-o", str(bcf_path), str(vcf_gz)],
        check=True,
    )
    subprocess.run([bcftools, "index", "-f", str(bcf_path)], check=True)


def build_cram(bam_path: Path, fasta_path: Path, cram_path: Path) -> None:
    samtools = shutil.which("samtools")
    if not samtools:
        print("WARN: samtools not on PATH; skip small.cram", file=sys.stderr)
        return
    subprocess.run(
        [
            samtools,
            "view",
            "-C",
            "-T",
            str(fasta_path),
            "-o",
            str(cram_path),
            str(bam_path),
        ],
        check=True,
    )
    subprocess.run([samtools, "index", str(cram_path)], check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "tests" / "fixtures",
        help="Output directory for Tier-0 fixtures",
    )
    args = parser.parse_args()
    out: Path = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    fasta = out / "small.fa"
    bam = out / "small.bam"
    vcf = out / "small.vcf.gz"
    bcf = out / "small.bcf"

    write_fasta(fasta)
    build_bam(bam)
    build_vcf(vcf)
    build_bcf(vcf, bcf)
    cram = out / "small.cram"
    build_cram(bam, fasta, cram)

    expected = [
        fasta,
        Path(str(fasta) + ".fai"),
        bam,
        Path(str(bam) + ".bai"),
        vcf,
        Path(str(vcf) + ".tbi"),
        Path(str(vcf).removesuffix(".gz")),
    ]
    if bcf.is_file():
        expected.extend([bcf, Path(str(bcf) + ".csi")])
    if cram.is_file():
        expected.extend([cram, Path(str(cram) + ".crai")])
    missing = [p for p in expected if not p.exists()]
    if missing:
        print("ERROR: missing outputs:", ", ".join(str(p) for p in missing), file=sys.stderr)
        return 1

    with pysam.AlignmentFile(str(bam), "rb") as af:
        n = sum(1 for _ in af.fetch(until_eof=True))
    print(f"Wrote fixtures to {out}")
    print(f"  BAM records (until_eof): {n}")
    print(f"  VCF: {vcf}")
    if bcf.is_file():
        print(f"  BCF: {bcf}")
    if cram.is_file():
        print(f"  CRAM: {cram}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
