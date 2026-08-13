# SAMRust

**面向 Linux / HPC 的 Rust-native、多线程、pysam 兼容 HTS 读写库。**

> Correctness tracks **pysam**. Performance tracks **[rubam](https://github.com/victormar1/rubam)**.
>
> 功能与结果对齐 pysam；运行效率对齐 rubam。rubam 没有的 API 在 16T 三方表填 **NA**，不静默跳过。

当前版本 **v0.1.0**（2026-08-13）。Python 3.10–3.13，Linux x86_64。许可证 [MIT](LICENSE)。

- 仓库：https://github.com/Caizhaohui/SAMRust
- Release / wheel：https://github.com/Caizhaohui/SAMRust/releases/tag/v0.1.0

---

## 目录

- [定位](#定位)
- [功能一览](#功能一览)
- [安装](#安装)
- [快速开始](#快速开始)
- [坐标系统](#坐标系统)
- [Python API](#python-api)
  - [AlignmentFile（BAM）](#alignmentfilebam)
  - [AlignedSegment](#alignedsegment)
  - [统计：count / coverage / depth / pileup](#统计count--coverage--depth--pileup)
  - [并行 fetch](#并行-fetch)
  - [CRAM](#cram)
  - [VariantFile](#variantfile)
- [CLI](#cli)
- [与 pysam / rubam 的关系](#与-pysam--rubam-的关系)
- [明确不做](#明确不做)
- [性能](#性能)
- [架构](#架构)
- [仓库结构](#仓库结构)
- [开发与测试](#开发与测试)
- [文档索引](#文档索引)
- [许可与致谢](#许可与致谢)

---

## 定位

SAMRust 把重测序里最常用的 pysam 路径（区域 `fetch`、`count`、`count_coverage`、depth、pileup、只读 `VariantFile`）做到 **Python 语义不变、热点在 Rust**。解码走 [noodles](https://github.com/zaeleus/noodles)，区域并行走 rayon，统计路径不把整条 BAM 记录物化成 Python 对象。

适合：

- 已有 pysam 脚本，只想把 `count` / coverage / depth / pileup 换成更快实现
- HPC 上对嗜热毁丝霉等真菌 BAM 做窗口统计、候选位点 recount
- 需要 **1 线程输出 == N 线程输出** 的确定性并行

不适合：当作完整 pysam / samtools / bcftools 替代品（见 [明确不做](#明确不做)）。

---

## 功能一览

| 能力 | 状态 | 说明 |
|------|------|------|
| BAM 顺序迭代 + 索引 `fetch`（BAI/CSI） | ✅ | Python 0-based half-open |
| `count` / `count_coverage` / `depth_*` / `pileup_counts` | ✅ | 串行 + `threads=`；1T==NT |
| `parallel_fetch` / `iter_batches` | ✅ | BAM only |
| `VariantFile` 读路径 | ✅ | VCF / VCF.gz+TBI / BCF+CSI；无 writer |
| CRAM 顺序迭代 + CRAI `fetch` | 评估级 | 需 FASTA；**统计 API 未实现** |
| `samrust recount` | ✅ | 候选位点 A/C/G/T/N/DP/ALT；**不是 caller** |
| 写 BAM/VCF、variant caller、mapper | ❌ | 不在 v0.1 |
| `FastaFile` / `TabixFile` 独立 API | planned | post-v0.1 |
| PyPI `pip install samrust` | ❌ | v0.1 走 GitHub Release wheel |
| Windows | ❌ | 以 Linux/HPC 为先 |

完整矩阵：[COMPATIBILITY.md](COMPATIBILITY.md)。

---

## 安装

### 环境

- **OS**：Linux x86_64（manylinux2014 / glibc 2.17+）
- **Python**：3.10 / 3.11 / 3.12 / 3.13
- 可选：**NumPy**（`depth_numpy` / `count_coverage` / `pileup_counts` 返回 `ndarray`）
- CRAM：参考 FASTA + `.fai`

v0.1 **未发布到 PyPI**。

### 从 GitHub Release 安装 wheel

打开 [Releases](https://github.com/Caizhaohui/SAMRust/releases/tag/v0.1.0)，按 CPython ABI 选择文件：

| ABI | 文件名（节选） |
|-----|----------------|
| 3.10 | `samrust-0.1.0-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl` |
| 3.11 | `…-cp311-…` |
| 3.12 | `…-cp312-…` |
| 3.13 | `…-cp313-…` |

```bash
pip install https://github.com/Caizhaohui/SAMRust/releases/download/v0.1.0/samrust-0.1.0-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
python -c "import samrust; print(samrust.__version__)"
```

### 从源码（开发）

需要 Rust **1.82+** 与 [maturin](https://www.maturin.rs/)。

```bash
git clone https://github.com/Caizhaohui/SAMRust.git
cd SAMRust
pip install maturin
maturin develop --release
python -c "import samrust; print(samrust.__version__)"
```

本机打 wheel：`maturin build --release -o dist`。登录节点打出的包往往不是 manylinux，换机器请用 Release 产物。

HPC 若访问 crates.io 失败，在**本地未提交**的 `.cargo/config.toml` 配置镜像，见 [INSTALL.md](INSTALL.md)。

---

## 快速开始

```python
import samrust

bam = samrust.AlignmentFile("sample.bam", "rb")
print(samrust.__version__, bam.references, bam.lengths)

# 区域均为 0-based half-open [start, stop)
n = bam.count("chr1", 0, 1000, threads=8)
a, c, g, t = bam.count_coverage("chr1", 0, 1000, threads=8)
depth = bam.depth_numpy("chr1", 0, 1000, threads=8)
pu = bam.pileup_counts("chr1", 0, 1000, threads=8)

for rec in bam.fetch("chr1", 0, 1000):
    print(rec.query_name, rec.reference_start, rec.cigarstring, rec.get_tag("NM"))

vf = samrust.VariantFile("sample.vcf.gz")
for rec in vf.fetch("chr1", 0, 1_000_000):
    print(rec.chrom, rec.pos, rec.ref, rec.alts)  # pos 为 1-based
```

从 pysam 迁移时，把 `import pysam` 换成 `import samrust as pysam` 通常只需核对 [与 pysam / rubam 的关系](#与-pysam--rubam-的关系) 中的差异。

---

## 坐标系统

**所有 Python 区间参数（`fetch` / `count` / coverage / depth / pileup / `VariantFile.fetch`）都是 0-based half-open `[start, stop)`。**

```text
bam.fetch("chr1", 100, 200)   # 覆盖参考位置 100, 101, …, 199
```

- 空区间 `[start, start)` 返回空（与 pysam **BAM** 一致）。pysam **CRAM** 在 `start==stop` 时仍可能吐出覆盖该点的记录，SAMRust **不复制**该行为。
- `AlignedSegment.reference_start`：0-based，与 pysam 相同。
- `VariantRecord.pos`：**1-based**（与 pysam 相同）；`start` / `stop` 仍为 0-based half-open。
- Rust 内部与 noodles 的 1-based 转换只允许出现在 `coords.rs`。

---

## Python API

包入口：`AlignmentFile`、`AlignedSegment`、`VariantFile`、`VariantRecord`、`__version__`。

### AlignmentFile（BAM）

```python
bam = samrust.AlignmentFile("sample.bam", "rb")
# 属性：filename, mode, references, lengths, nreferences, header
# header 为 dict：nreferences / references / lengths
# 上下文管理器、close()、reset()、顺序 for rec in bam
```

- 只支持读：`mode="rb"`（CRAM 另支持 `"rc"`）。
- 需要索引的操作（`fetch`、统计、并行）要求同目录 `.bai` 或 CSI。
- `count` 的 `read_callback`：`"nofilter"`（默认，与 pysam `count` 默认一致，含 placed-unmapped）或 `"all"`（排除 unmapped / secondary / QC-fail / duplicate）。
- `count_coverage` 的 `read_callback` 默认是 `"all"`（与 pysam 一致），另有 `quality_threshold=15`。

### AlignedSegment

常用字段与 pysam 对齐：

| 属性 / 方法 | 含义 |
|-------------|------|
| `query_name`, `flag`, `mapping_quality` | QNAME / FLAG / MAPQ（缺省 MAPQ 为 255） |
| `reference_id`, `reference_name`, `reference_start` | 参考；start 为 0-based |
| `cigarstring`, `cigartuples` | CIGAR |
| `query_sequence`, `query_length`, `query_qualities` | 序列与碱基质量 |
| `next_reference_*`, `template_length` | 配对信息 |
| `is_paired` / `is_unmapped` / `is_reverse` / `is_secondary` / `is_duplicate` / … | FLAG 布尔属性 |
| `has_tag(tag)`, `get_tag(tag)` | 辅助标签 |

`fetch` 仍会物化 `AlignedSegment`。`count` / coverage / depth / pileup 在 Rust 内流式扫描 noodles 记录，不经过这一层。

### 统计：count / coverage / depth / pileup

```python
n = bam.count("chr1", 0, 1000, read_callback="nofilter", threads=8)

a, c, g, t = bam.count_coverage(
    "chr1", 0, 1000, quality_threshold=15, read_callback="all", threads=8
)

blocks = bam.depth_blocks("chr1", 0, 1000, threads=8)
# [(start, length, depth), ...]  连续同深度 run-length

depth = bam.depth_numpy("chr1", 0, 1000, threads=8)
# 有 NumPy 时为 ndarray；否则为 list[int]

pu = bam.pileup_counts(
    "chr1", 0, 1000, min_base_quality=0, min_mapping_quality=0, threads=8
)
# dict：A / C / G / T / N / depth（等长数组）
```

语义要点：

- **并行**：`threads>1` 时统计路径每线程 **1 个索引块**（与 rubam fast mode 对齐）。必须满足 **1 线程输出 == N 线程输出**。
- **depth**：只计 CIGAR `M/=/X`，**不计** `D`/`N`；含模糊碱基 N。与 samtools depth、rubam `get_depths` 对齐，**不等于** `count_coverage` 的 A+C+G+T（后者受 BQ 与 `read_callback` 约束）。
- **pileup_counts**：默认排除 unmapped / secondary / supplementary / qcfail / duplicate；跳过 del/refskip 上的碱基。这是 **计数数组**，不是 pysam 列式 `pileup()` 迭代器。
- **CRAM** 上上述统计方法抛 `RuntimeError`（M11 范围外）。

### 并行 fetch

```python
# 多区域；ordered=True 时按块序合并
recs = bam.parallel_fetch([("chr1", 0, 1000), ("chr2", 0, 500)], threads=8, ordered=True)

batches = bam.iter_batches(batch_size=256, threads=8, ordered=True)
```

fetch 调度仍为每线程约 4 块（与统计路径的 1 块不同）。CRAM 上这些 API 同样 `NotImplemented`。

### CRAM

```python
cram = samrust.AlignmentFile(
    "sample.cram", "rc", reference_filename="ref.fa"
)
# 若同目录有 sample.fa / .fasta / .fna，可省略 reference_filename
# 需要 sample.cram.crai 与 FASTA .fai

for rec in cram: ...
for rec in cram.fetch("chr1", 0, 1000): ...

cram.count("chr1", 0, 1000)  # RuntimeError：统计仍为 BAM-only
```

noodles CRAM 3.x 评估级读路径；加密 / 非 3.x / HTSlib-only codec 不在 v0.1。optional `--features htslib-backend` **仅设计、未实现**。

### VariantFile

只读。支持明文 VCF、`VCF.gz`+`.tbi`、BCF+`.csi`。

```python
vf = samrust.VariantFile("sample.vcf.gz", "r")  # 或 "rb"
print(vf.samples, vf.header.contigs)

for rec in vf:  # 顺序
    print(rec.chrom, rec.pos, rec.ref, rec.alts, rec.qual, rec.filter, rec.info)
    print(rec.samples["NA12878"]["GT"], rec.samples["NA12878"]["DP"], rec.samples["NA12878"]["AD"])

for rec in vf.fetch("chr1", 0, 1_000_000):
    ...
```

| 字段 | 坐标 / 语义 |
|------|-------------|
| `fetch(contig, start, stop)` | 0-based half-open |
| `rec.pos` | 1-based（pysam） |
| `rec.start` / `rec.stop` | 0-based half-open |
| `rec.samples[i 或 name]["GT"/"DP"/"AD"]` | 读路径常用 FORMAT |

无 writer、无 mutation、无 variant caller。16T 三方表 **VariantFile 行 rubam = NA**（rubam 无等价 API）。

---

## CLI

```bash
cargo build -p samrust-cli --release
# 或安装后使用 target/release/samrust

samrust version
samrust dump-records --bam sample.bam --limit 1000
samrust recount \
  --bam sample.bam \
  --sites candidates.bed \
  --sample SAMPLE \
  --threads 8 \
  --min-alt 10 \
  --output recount.tsv
```

`recount` 是 **benchmark / 校验工具**，不是 SNP caller。位点文件：BED（`chrom start stop REF>ALT`）或 TSV（`chrom pos ref alt`）。真菌数据上门禁：`ALT_COUNT >= 10` 的位点集合与 pysam 一致（15453 候选 → 6895 位点，见 [BENCHMARKS.md](BENCHMARKS.md)）。

---

## 与 pysam / rubam 的关系

| 维度 | 规范 |
|------|------|
| Python 坐标、FLAG 过滤、count/coverage 默认 callback | **pysam** |
| BAM 统计墙钟 / 资源 | **rubam**（同条件 wall / RSS） |
| rubam 没有的功能（如 `VariantFile`） | 三方表保留该行，rubam 填 **NA** |

不是 100% pysam：

- 无 `AlignmentFile` 写入、无 SAM 文本写、无 `pileup()` 列式迭代器（有 `pileup_counts`）
- 无独立 `FastaFile` / `TabixFile`
- CRAM 无 count/depth/coverage/pileup
- `header` 目前是精简 dict，不是完整 SAM header 对象
- `count_coverage` 的 pysam 默认 `read_callback="all"`；`count` 默认 `"nofilter"`

深度语义与 **samtools depth / rubam `get_depths`** 对齐（M/=/X），不要用 `count_coverage` 之和当 depth oracle。

---

## 明确不做

v0.1 冻结：read mapper、SNP/Indel caller、de novo assembler、完整 samtools/bcftools CLI、GUI、Windows 优先、一次性覆盖全部 pysam API、把 tokio 引入 BAM 热路径。

真实重测序 BAM **只读**；输出只进 `benchmark/results` 或 scratch。

---

## 性能

正式表：作业 **2312423**，Slurm `qcpu_18i` / **bnode2** / Xeon Silver 4116，`Mt-35-15d-1.markdup.bam`，`NC_016472.1:0-100000`，release，median of 3。原始文件：`benchmark/results/compare_pysam_rubam_samrust.fungal_mt35.{csv,json}`。

16T 墙钟（秒）；加速比 = 对照 / SAMRust（>1 表示 SAMRust 更快）：

| 负载 | SAMRust | pysam | rubam | vs pysam | vs rubam |
|------|---------|-------|-------|----------|----------|
| count | 0.033 | 0.118 | 0.135 | 3.59× | **4.11×** |
| count_coverage | 0.076 | 4.472 | 0.106 | 59× | **1.41×** |
| depth | 0.041 | 1.022 | 0.048 | 25× | **1.17×** |
| pileup_counts | 0.090 | 50.7 | 0.093 | 562× | **1.03×** |
| VariantFile | — | — | **NA** | — | **NA** |

SAMRust 1→8T 约 1.9–2.3×；16T 在 100 kb 窗口上相对 8T 持平。未跑 **32T**（`qcpu_18i` 节点约 24 核）。coverage 的 1T vs-rubam 不公平（rubam 未传 `threads`）。叙事与 scaling / RSS：[BENCHMARKS.md](BENCHMARKS.md)。

实现要点（M10）：统计路径每线程 1 次 BAM+BAI 查询；CIGAR 用区间裁剪；A/C/G/T 用 256 项 LUT。

---

## 架构

```text
Python  samrust.AlignmentFile / VariantFile
            │  PyO3
            ▼
crates/samrust-python     绑定、GIL release、NumPy
            │
            ▼
crates/samrust-core       #![deny(unsafe_code)]
  bam / indexed / cram / vcf
  coords.rs               唯一允许 +1/-1 的地方
  depth / pileup          for_each_raw 流式热点
  parallel                rayon + 有界 channel；stats=1 chunk/thread
            │
            ▼
noodles 0.115             BAM / CRAM / VCF / BCF / FASTA
```

热循环禁止写在 Python。并行禁止把整份 BAM/基因组结果拷到每个 worker。

---

## 仓库结构

| 路径 | 作用 |
|------|------|
| `crates/samrust-core` | 纯 Rust 核心 |
| `crates/samrust-python` | PyO3 扩展 `samrust._samrust` |
| `crates/samrust-cli` | `samrust` CLI |
| `python/samrust` | 包表面 |
| `tests/fixtures` | Tier-0 `small.bam/.cram/.vcf/.bcf`（CI） |
| `tests/parity` | pysam / bcftools oracle |
| `benchmark/` | 三方对比脚本；16T 表已入库 |
| `scripts/` | fixture、parity、Slurm 提交 |
| `.github/workflows` | Rust / Python CI；`v*` tag 打 manylinux wheel |

---

## 开发与测试

```bash
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace          # 含 Interval proptest
maturin develop --release
python -m pytest tests -q
python scripts/prepare_fixture.py
```

| 测试 | 在哪跑 |
|------|--------|
| 单元 / property / Tier-0 pytest | 登录节点或 GitHub Actions |
| 真菌 BAM、16T、M7/M8 recount | **仅** Slurm 分区 **`qcpu_18i`**，禁止登录节点 |

```bash
bash scripts/submit_three_way_benchmark.sh
bash scripts/submit_m7_m8_validation.sh
```

Agent 规则：[AGENTS.md](AGENTS.md)。开发流程：[DEVELOPMENT.md](DEVELOPMENT.md)、[DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md)。

---

## 文档索引

| 文件 | 内容 |
|------|------|
| [INSTALL.md](INSTALL.md) | wheel / 源码 / Cargo 镜像 |
| [COMPATIBILITY.md](COMPATIBILITY.md) | pysam API 矩阵 |
| [BENCHMARKS.md](BENCHMARKS.md) | 16T 三方表、RSS、M8 recount |
| [CHANGELOG.md](CHANGELOG.md) | v0.1.0 变更 |
| [DEVELOPMENT.md](DEVELOPMENT.md) | 本地开发 |
| [NOTICE](NOTICE) | 行为参考（pysam / rubam / noodles / HTSlib），无大段抄代码 |
| [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md) | 里程碑 M0–M12 与设计约束 |

---

## 许可与致谢

MIT，见 [LICENSE](LICENSE)。

实现为独立代码。API 与 oracle 参考了 [pysam](https://github.com/pysam-developers/pysam)、[rubam](https://github.com/victormar1/rubam)、[noodles](https://github.com/zaeleus/noodles)、samtools / bcftools / HTSlib，详见 [NOTICE](NOTICE)。

---

## 状态

**v0.1.0** 完成 BAM 分析路径、VariantFile 读路径、CRAM 评估读路径，以及 GitHub manylinux 发布。后续（PyPI、HTSlib backend、`FastaFile`）需另行提出。
