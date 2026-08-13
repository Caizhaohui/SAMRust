# SAMRust

Rust-native, multi-threaded, **pysam-compatible** HTS access for Linux/HPC.

Correctness tracks **pysam**. Performance tracks **[rubam](https://github.com/victormar1/rubam)**. If rubam has no equivalent API, the 16T table keeps the row and fills rubam cells with **NA**.

## Why it is faster than pysam

pysam is excellent, but region scans, depth, and pileup often spend time in Python-bound loops. SAMRust keeps pysam coordinates and filters, and runs the hot path in Rust (`noodles` decode, rayon region chunks, NumPy when present).

On a 100 kb fungal BAM window (16 threads, job 2312423): `count` **3.6×** vs pysam, `depth` **25×**, `count_coverage` **59×**, `pileup_counts` **562×**. Versus rubam on the same node: **1.03–4.1×** (see [BENCHMARKS.md](BENCHMARKS.md)).

## pysam compatible?

v0.1 covers the AlignmentFile APIs used most in resequencing (`fetch`, `count`, `count_coverage`, depth, `pileup_counts`) plus a **read-only** `VariantFile`. Python coordinates are always **0-based half-open**. See [COMPATIBILITY.md](COMPATIBILITY.md) / [API_COMPATIBILITY.md](API_COMPATIBILITY.md).

Not a drop-in for the entire pysam surface (no writer/caller, no `FastaFile`/`TabixFile` APIs, CRAM stats not implemented).

## Install

Linux x86_64, Python 3.10+.

```bash
# GitHub Release wheel (v0.1) — pick the CPython tag (cp310–cp313) from
# https://github.com/Caizhaohui/SAMRust/releases
pip install ./samrust-0.1.0-cp312-*.whl
```

From source: see [INSTALL.md](INSTALL.md).

## Minimal example

```python
import samrust

bam = samrust.AlignmentFile("sample.bam", "rb")
print(samrust.__version__, bam.references)

n = bam.count("chr1", 0, 1000, threads=8)
a, c, g, t = bam.count_coverage("chr1", 0, 1000, threads=8)
depth = bam.depth_numpy("chr1", 0, 1000, threads=8)
pu = bam.pileup_counts("chr1", 0, 1000, threads=8)

for rec in bam.fetch("chr1", 0, 1000):  # [start, stop)
    print(rec.query_name, rec.reference_start, rec.cigarstring)

vf = samrust.VariantFile("sample.vcf.gz")
for rec in vf.fetch("chr1", 0, 1_000_000):
    print(rec.chrom, rec.pos, rec.ref, rec.alts)
```

`rec.pos` on `VariantRecord` is **1-based**, like pysam. `fetch` start/stop stay 0-based half-open.

## Project layout

| Path | Role |
|------|------|
| `crates/samrust-core` | Pure Rust core (no Python) |
| `crates/samrust-python` | PyO3 bindings |
| `crates/samrust-cli` | CLI utilities (`samrust recount`, …) |
| `python/samrust` | Python package surface |
| `tests/` | unit / parity / fixtures |
| `benchmark/` | configs, scripts, published 16T table |

## Develop

See [DEVELOPMENT.md](DEVELOPMENT.md) and [AGENTS.md](AGENTS.md).

```bash
cargo fmt --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
maturin develop
python -m pytest
```

## License

MIT — [LICENSE](LICENSE). Third-party behavior references: [NOTICE](NOTICE).

## Status

**v0.1.0** (M12) — BAM analytics + VariantFile read path; CRAM sequential/fetch evaluation. Post-v0.1 work (PyPI, HTSlib backend, FastaFile) is out of scope unless requested.
