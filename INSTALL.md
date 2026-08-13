# Install SAMRust

v0.1 ships **Linux x86_64** wheels. Python **3.10–3.13**.

## From a GitHub Release wheel

1. Open [Releases](https://github.com/Caizhaohui/SAMRust/releases) and download the `manylinux` wheel for your CPython ABI (`cp310`, `cp311`, `cp312`, or `cp313`).
2. Install:

```bash
pip install samrust-0.1.0-*.whl
python -c "import samrust; print(samrust.__version__)"
```

Use the asset name from the [Releases](https://github.com/Caizhaohui/SAMRust/releases) page (`cp310`–`cp313`). PyPI `pip install samrust` is **not** wired in v0.1 (GitHub Release is the distribution channel).

## From source (editable)

Needs a Rust toolchain (MSRV **1.82**) and [maturin](https://www.maturin.rs/).

```bash
git clone https://github.com/Caizhaohui/SAMRust.git
cd SAMRust
pip install maturin
maturin develop          # debug
# maturin develop --release
python -c "import samrust; print(samrust.__version__)"
```

`pip install -e .` also works (pyproject build-backend is maturin).

### Local wheel

```bash
maturin build --release -o dist
pip install dist/samrust-*.whl
```

A wheel built on an HPC login node is often `linux_x86_64` (not manylinux). Prefer the CI manylinux artifact for other machines.

### crates.io mirror (HPC)

If crates.io is blocked, point Cargo at a mirror **locally** (do not commit this file; `.cargo/config.toml` is gitignored):

```toml
# .cargo/config.toml
[source.crates-io]
replace-with = "ustc"

[source.ustc]
registry = "sparse+https://mirrors.ustc.edu.cn/crates.io-index/"
```

## Optional runtime extras

| Extra | Why |
|-------|-----|
| NumPy | `depth_numpy` / `count_coverage` return `ndarray` when importable |
| A FASTA + `.fai` | CRAM `AlignmentFile` (`mode="rb"` or `"rc"`) |
| pysam / samtools / bcftools | parity tests and fixture regeneration only |

```bash
pip install numpy
python scripts/prepare_fixture.py   # needs pysam + samtools + bcftools
```

## HPC notes

Heavy fungal BAM jobs belong on Slurm partition **`qcpu_18i`**, not the login node. See [DEVELOPMENT.md](DEVELOPMENT.md).
