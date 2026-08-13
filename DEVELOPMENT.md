# Developing SAMRust

Follow [AGENTS.md](AGENTS.md) and [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md). One milestone at a time.

## Toolchain

- Rust **1.82+** (workspace `rust-version`)
- Python **3.10+**
- `maturin`, `pytest`
- Optional oracles: `pysam`, `samtools`, `bcftools`, NumPy
- Example HPC env: conda `mt-var`

```bash
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
maturin develop
python -m pytest tests -q
```

Regenerate Tier-0 fixtures (committed under `tests/fixtures/`):

```bash
python scripts/prepare_fixture.py
```

BAM core dump parity:

```bash
cargo build -p samrust-cli
python scripts/parity_bam_core.py --bam tests/fixtures/small.bam
```

## Tests that belong where

| Kind | Where | CI |
|------|-------|----|
| Rust unit + property (`Interval`) | `cargo test --workspace` | GitHub |
| Python parity vs pysam | `tests/parity/` | GitHub (Tier-0 fixtures) |
| Real fungal BAM / 16T / recount | `scripts/submit_*.sh` | **not** GitHub; Slurm `qcpu_18i` |

Do not run fungal full-genome / 16T scaling on the login node.

```bash
bash scripts/submit_m7_m8_validation.sh
bash scripts/submit_three_way_benchmark.sh
```

`qcpu_18i` nodes have about 24 CPUs; do not request 32+ `cpus-per-task`. That is why the published table is **1/4/8/16T**, not 32T.

## Layout

| Path | Role |
|------|------|
| `crates/samrust-core` | BAM/CRAM/VCF core (`#![deny(unsafe_code)]`) |
| `crates/samrust-python` | PyO3 `AlignmentFile` / `VariantFile` |
| `crates/samrust-cli` | `samrust recount` and dump helpers |
| `python/samrust` | import surface |
| `tests/fixtures` | Tier-0 BAM/CRAM/VCF (CI) |
| `benchmark/results/` | local HPC outputs (gitignored except the published 16T files) |

## Release

Tags matching `v*` run [`.github/workflows/release.yml`](.github/workflows/release.yml): manylinux x86_64 wheels uploaded to the GitHub Release. See [INSTALL.md](INSTALL.md).
