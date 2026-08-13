# SAMRust

**面向 Linux / HPC 的 Rust-native、多线程、pysam 兼容 HTS 读写库。**

[![Python CI](https://github.com/Caizhaohui/SAMRust/actions/workflows/python.yml/badge.svg)](https://github.com/Caizhaohui/SAMRust/actions/workflows/python.yml)
[![Rust CI](https://github.com/Caizhaohui/SAMRust/actions/workflows/rust.yml/badge.svg)](https://github.com/Caizhaohui/SAMRust/actions/workflows/rust.yml)
[![Release](https://img.shields.io/github/v/release/Caizhaohui/SAMRust)](https://github.com/Caizhaohui/SAMRust/releases)
[![Python](https://img.shields.io/badge/python-3.10--3.13-blue)](https://github.com/Caizhaohui/SAMRust)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

当前版本 **v0.1.1**（2026-08-13）。Python 3.10–3.13，Linux x86_64。许可证 [MIT](LICENSE)。

- 仓库：https://github.com/Caizhaohui/SAMRust
- Release / wheel：https://github.com/Caizhaohui/SAMRust/releases/tag/v0.1.1
- v0.1.1 修复与兼容性变更：[CHANGELOG](CHANGELOG.md) / [REVIEW](REVIEW.md) / [COMPATIBILITY](COMPATIBILITY.md)

---

## 目录

- [定位](#定位)
- [功能一览](#功能一览)
- [安装](#安装)
- [快速开始](#快速开始)
- [坐标系统](#坐标系统)
- [Python API 参考](#python-api-参考)
  - [AlignmentFile（BAM / CRAM）](#alignmentfilebam--cram)
  - [AlignedSegment](#alignedsegment)
  - [统计：count / coverage / depth / pileup](#统计count--coverage--depth--pileup)
  - [并行：parallel_fetch / iter_batches](#并行parallel_fetch--iter_batches)
  - [CRAM](#cram)
  - [VariantFile（VCF / BCF）](#variantfilevcf--bcf)
  - [错误处理](#错误处理)
- [CLI](#cli)
- [与 pysam / rubam 的关系](#与-pysam--rubam-的关系)
- [明确不做](#明确不做)
- [性能](#性能)
- [架构](#架构)
- [仓库结构](#仓库结构)
- [开发与测试](#开发与测试)
- [故障排查](#故障排查)
- [路线图](#路线图)
- [文档索引](#文档索引)
- [许可与致谢](#许可与致谢)

---

## 定位

SAMRust 把重测序里最常用的 pysam 路径（区域 `fetch`、`count`、`count_coverage`、depth、pileup、只读 `VariantFile`）做到 **Python 语义不变、热点在 Rust**。解码走 [noodles](https://github.com/zaeleus/noodles)，区域并行走 rayon，统计路径不把整条 BAM 记录物化成 Python 对象。

适合：

- 已有 pysam 脚本，只想把 `count` / coverage / depth / pileup 换成更快实现
- HPC 上对真菌等重测序 BAM 做窗口统计、候选位点 recount
- 需要 **1 线程输出 == N 线程输出** 的确定性并行

不适合：当作完整 pysam / samtools / bcftools 替代品（见 [明确不做](#明确不做)）。

设计基线：**功能和结果对齐 pysam；运行效率对齐 rubam**（[DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md) §1.2）。

---

## 功能一览

| 能力 | 状态 | 说明 |
|------|------|------|
| BAM 顺序迭代 + 索引 `fetch`（BAI/CSI） | ✅ | Python 0-based half-open |
| `count` / `count_coverage` / `depth_*` / `pileup_counts` | ✅ | 串行 + `threads=`；1T==NT bit-exact |
| `parallel_fetch` / `iter_batches` | ✅ | BAM only；区域合并 + 位置归属，无哈希去重 |
| `VariantFile` 读路径 | ✅ | VCF / VCF.gz+TBI / BCF+CSI；无 writer |
| CRAM 顺序迭代 + CRAI `fetch` | 评估级 | 需参考 FASTA；**统计 API 未实现** |
| `samrust recount` | ✅ | 候选位点 A/C/G/T/N/DP/ALT；**不是 caller** |
| 写 BAM/VCF、variant caller、mapper | ❌ | 不在 v0.1 |
| `FastaFile` / `TabixFile` 独立 API | planned | post-v0.1 |
| PyPI `pip install samrust` | ❌ | v0.1 走 GitHub Release wheel |
| Windows | ❌ | 以 Linux/HPC 为先 |

完整矩阵：[COMPATIBILITY.md](COMPATIBILITY.md)。

---

## 安装

### 环境要求

- **OS**：Linux x86_64（manylinux2014 / glibc 2.17+）
- **Python**：3.10 / 3.11 / 3.12 / 3.13
- 可选：**NumPy**（`depth_numpy` / `count_coverage` / `pileup_counts` 返回 `ndarray`；无 NumPy 时退化为 `list`）
- CRAM：参考 FASTA + `.fai`

v0.1 **未发布到 PyPI**。

### 从 GitHub Release 安装 wheel（推荐）

打开 [Releases](https://github.com/Caizhaohui/SAMRust/releases/tag/v0.1.1)，按 CPython ABI 选择文件：

| ABI | 文件名（节选） |
|-----|----------------|
| 3.10 | `samrust-0.1.1-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl` |
| 3.11 | `…-cp311-…` |
| 3.12 | `…-cp312-…` |
| 3.13 | `…-cp313-…` |

```bash
pip install https://github.com/Caizhaohui/SAMRust/releases/download/v0.1.1/samrust-0.1.1-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
python -c "import samrust; print(samrust.__version__)"
```

### 从源码安装（开发）

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
    print(rec.query_name, rec.reference_start, rec.reference_end,
          rec.cigarstring, rec.get_tag("NM"))

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

- `start` / `stop` 省略时为 `None`（= 整条 contig）；负数抛 `ValueError`（与 pysam 一致）。
- `stop` 超出 contig 长度时静默截断到 contig 末尾（pysam 行为）；`start` 越界返回空结果。
- 空区间 `[start, start)` 返回空（与 pysam **BAM** 一致）。pysam **CRAM** 在 `start==stop` 时仍可能吐出覆盖该点的记录，SAMRust **不复制**该行为。
- `AlignedSegment.reference_start`：0-based；`reference_end`：0-based exclusive，均与 pysam 相同。
- `VariantRecord.pos`：**1-based**（与 pysam 相同）；`start` / `stop` 仍为 0-based half-open。
- Rust 内部与 noodles 的 1-based 转换只允许出现在 `coords.rs`。

---

## Python API 参考

包入口：`AlignmentFile`、`AlignedSegment`、`VariantFile`、`VariantRecord`、`version()`、`__version__`。

### AlignmentFile（BAM / CRAM）

```python
bam = samrust.AlignmentFile("sample.bam", mode="rb", reference_filename=None)
```

| 参数 | 默认 | 说明 |
|------|------|------|
| `filename` | — | BAM / CRAM 路径（按扩展名识别 `.cram`） |
| `mode` | `"rb"` | 只读；CRAM 可用 `"rc"` |
| `reference_filename` | `None` | CRAM 参考 FASTA；省略时尝试同目录 `sample.fa/.fasta/.fna` |

属性与方法：

| 成员 | 类型 / 签名 | 说明 |
|------|-------------|------|
| `filename` / `mode` | `str` | 路径与打开模式 |
| `references` | `list[str]` | 参考序列名 |
| `lengths` | `list[int]` | 参考序列长度 |
| `nreferences` | `int` | 参考序列条数 |
| `header` | `dict` | pysam 风格：`HD` / `SQ`（`SN`/`LN`）/ `RG` / `PG`（v0.1.1 起） |
| `close()` / `reset()` | — | 关闭 / 回到文件开头重新迭代 |
| `for rec in bam` | 迭代器 | 顺序遍历（含 unmapped 尾记录）；中途 `iter(bam)` 从逻辑位置继续 |
| `fetch(contig, start=None, stop=None)` | `Iterator[AlignedSegment]` | 索引区域查询，需 `.bai` / CSI |
| `count(...)` | `int` | 见[统计](#统计count--coverage--depth--pileup) |
| `count_coverage(...)` | `(A, C, G, T)` | 同上 |
| `depth_blocks(...)` | `list[(start, len, depth)]` | 同上 |
| `depth_numpy(...)` | `ndarray` 或 `list[int]` | 同上 |
| `pileup_counts(...)` | `dict[str, array]` | 同上 |
| `parallel_fetch(regions, threads=1, ordered=True)` | `list[AlignedSegment]` | 见[并行](#并行parallel_fetch--iter_batches) |
| `iter_batches(batch_size=256, threads=1, ordered=True)` | `BatchIterator` | 同上 |

- 需要索引的操作（`fetch`、统计、并行）要求同目录 `.bai` 或 CSI。
- 支持上下文管理器：`with samrust.AlignmentFile(...) as bam: ...`。
- `close()` 后再访问 `header` / `references` / `lengths` / `nreferences` 抛 `ValueError`（v0.1.1 起，不再是 `PanicException`）。

### AlignedSegment

常用字段与 pysam 对齐：

| 属性 / 方法 | 含义 |
|-------------|------|
| `query_name`, `flag`, `mapping_quality` | QNAME / FLAG / MAPQ（缺省 MAPQ 为 255） |
| `reference_id`, `reference_name`, `reference_start` | 参考；start 为 0-based |
| `reference_end` | 0-based exclusive 末端；unmapped（含 placed-unmapped）或无 CIGAR 时为 `None`（v0.1.1 起） |
| `cigarstring`, `cigartuples` | CIGAR（tuples 为 `(op, len)`，op 编码同 pysam） |
| `query_sequence`, `query_length`, `query_qualities` | 序列与碱基质量 |
| `next_reference_id`, `next_reference_start`, `template_length` | 配对信息 |
| `is_paired` / `is_proper_pair` / `is_unmapped` / `mate_is_unmapped` / `is_reverse` / `mate_is_reverse` / `is_read1` / `is_read2` / `is_secondary` / `is_qcfail` / `is_duplicate` / `is_supplementary` | FLAG 布尔属性 |
| `has_tag(tag)`, `get_tag(tag)` | 辅助标签 |

`get_tag` 返回类型（v0.1.1 起与 pysam 逐类型一致）：

| BAM tag 类型 | Python 返回 |
|--------------|-------------|
| `A` / `Z` / `H` | `str` |
| `c` / `s` / `i` | `int` |
| `C` / `S` / `I` | `int` |
| `f` | `float` |
| `B,c/B,C/B,s/B,S/B,i/B,I` | `array.array`（保留原 typecode） |
| `B,f` | `array.array('f')` |

`fetch` 会物化 `AlignedSegment`；`count` / coverage / depth / pileup 在 Rust 内流式扫描 noodles 记录，不经过这一层。

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

- **并行确定性**：`threads>1` 时统计路径每线程 **1 个索引块**（与 rubam fast mode 对齐）。保证 **1 线程输出 == N 线程输出**（bit-exact）。
- **`read_callback`**：`count` 默认 `"nofilter"`（与 pysam `count` 默认一致，含 placed-unmapped）；`count_coverage` 默认 `"all"`（排除 unmapped / secondary / QC-fail / duplicate，与 pysam 一致）。
- **depth**：只计 CIGAR `M/=/X`，**不计** `D`/`N`；含模糊碱基 N。与 samtools depth、rubam `get_depths` 对齐，**不等于** `count_coverage` 的 A+C+G+T（后者受 BQ 与 `read_callback` 约束）。
- **pileup_counts**：默认排除 unmapped / secondary / supplementary / qcfail / duplicate；跳过 del/refskip 上的碱基。这是 **计数数组**，不是 pysam 列式 `pileup()` 迭代器（`pileup()` 推迟到 v0.2）。
- **CRAM** 上上述统计方法抛 `RuntimeError`（M11 范围外）。

### 并行：parallel_fetch / iter_batches

```python
# 多区域并行 fetch；重叠区域先合并再切块，每条记录恰好出现一次
# （完全相同的重复记录，如 cat-bam 场景，会被保留）
recs = bam.parallel_fetch([("chr1", 0, 1000), ("chr2", 0, 500)], threads=8, ordered=True)

# 全文件分批迭代：1T 为流式；MT 按 ~1 Mb 波次处理，内存有界
# 两种模式都包含文件末尾的 unmapped-unplaced 记录，且 1T == NT 逐条一致
it = bam.iter_batches(batch_size=256, threads=8, ordered=True)
for batch in it:        # batch: list[AlignedSegment]，长度 <= batch_size
    ...
```

fetch 调度为每线程约 4 块（与统计路径的 1 块不同）。CRAM 上这些 API 抛 `NotImplementedError`。

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

### VariantFile（VCF / BCF）

只读。支持明文 VCF、`VCF.gz`+`.tbi`、BCF+`.csi`。

```python
vf = samrust.VariantFile("sample.vcf.gz", mode="r")  # 或 "rb"
print(vf.samples)            # ["NA12878", ...]
print(vf.header.samples)     # VariantHeader: samples / contigs
print(vf.header.contigs)

for rec in vf:               # 顺序迭代
    print(rec.chrom, rec.pos, rec.ref, rec.alts, rec.qual, rec.filter, rec.info)
    print(rec.samples["NA12878"]["GT"], rec.samples["NA12878"]["DP"])

for rec in vf.fetch("chr1", 0, 1_000_000):   # 索引区域查询
    ...
```

`VariantRecord` 字段：

| 字段 | 坐标 / 语义 |
|------|-------------|
| `fetch(contig=None, start=None, stop=None)` | 0-based half-open；`contig=None` 为全文件 |
| `rec.chrom` / `rec.contig` | 染色体名 |
| `rec.pos` | 1-based（pysam） |
| `rec.start` / `rec.stop` | 0-based half-open |
| `rec.id` / `rec.ref` / `rec.alts` / `rec.alleles` | ID / REF / ALT / 全部等位 |
| `rec.qual` / `rec.filter` / `rec.info` | QUAL / FILTER / INFO |
| `rec.format` | 本记录出现的 FORMAT 键 |
| `rec.samples[i 或 name]["GT"/"DP"/"AD"]` | 读路径常用 FORMAT |

边界行为（v0.1.1）：

- header 无 contig length 时，`fetch(contig)`（`stop=None`）回退顺序扫描，不再返回空；显式 `stop` 仍走索引。
- 无 `.tbi` / `.csi` 的 VCF.gz：`fetch` 抛 `ValueError`（请先 `pysam.tabix_index` / `bcftools index`）；顺序迭代不受影响。
- 无 writer、无 mutation、无 variant caller。16T 三方表 **VariantFile 行 rubam = NA**（rubam 无等价 API）。

### 错误处理

| 情况 | 异常 |
|------|------|
| 负坐标、未知 contig、`start > stop`、无索引 VCF fetch、close 后访问 header | `ValueError` |
| 文件不存在 / 无法打开 / 索引缺失 | `OSError`（`IOError`） |
| CRAM 上调用统计 / 并行 API | `RuntimeError` / `NotImplementedError` |
| BAM / VCF 解析失败 | `ValueError`（含 noodles 错误上下文） |

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
- `header` 是 pysam 风格 dict（HD/SQ/RG/PG），不是可写的完整 SAM header 对象
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

实现要点（M10 + v0.1.1）：统计路径每线程 1 次 BAM+BAI 查询；CIGAR 用区间裁剪；A/C/G/T 用 256 项 LUT；并行统计 / recount 每线程复用一个 reader（`map_init`），不再每块重开 BAM+索引。

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
  parallel                rayon；stats=1 chunk/thread；fetch≈4 chunks/thread
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
| `python/samrust` | 包表面（`__init__.py` / `typing.pyi`） |
| `tests/fixtures` | Tier-0 `small.bam/.cram/.vcf/.bcf`（CI） |
| `tests/parity` | pysam / bcftools oracle、随机区域差分、v0.1.1 回归 |
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
python -m pytest tests -q --require-oracles   # CI 门禁：oracle 缺失即失败
SAMRUST_REAL_DATA=1 python -m pytest tests/parity/test_random_regions.py -q  # 真实数据差分
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

## 故障排查

| 症状 | 原因与处理 |
|------|-----------|
| `fetch` / 统计报 "index not found" | BAM 需要同目录 `.bai` 或 CSI：`samtools index sample.bam` |
| `VariantFile.fetch` 抛 `ValueError` | VCF.gz 缺 `.tbi` / `.csi`：`pysam.tabix_index("x.vcf.gz", preset="vcf")` |
| CRAM 打开失败 | 缺参考 FASTA 或 `.crai`：传 `reference_filename=` 并 `samtools faidx ref.fa` |
| CRAM 上 `count` 抛 `RuntimeError` | 统计 API 是 BAM-only（M11 范围外），先转 BAM |
| `depth_numpy` 返回 list 而非 ndarray | 未安装 NumPy：`pip install numpy` |
| 换机器后 wheel 无法安装 | 登录节点本地构建非 manylinux；改用 GitHub Release 产物 |
| 行为与 pysam 不一致 | 先查 [COMPATIBILITY.md](COMPATIBILITY.md)；未记录则提 issue |

---

## 路线图

v0.1.x：缺陷修复与 pysam 兼容性补齐（见 [CHANGELOG](CHANGELOG.md)）。

post-v0.1（需另行提出，不承诺排期）：

- `pileup()` 列式迭代器、流式 `fetch`（自持 reader）
- PyPI 发布
- `FastaFile` / `TabixFile` 独立 API
- HTSlib backend（`--features htslib-backend`，已设计未实现）
- CRAM 统计 / 并行 API

---

## 文档索引

| 文件 | 内容 |
|------|------|
| [INSTALL.md](INSTALL.md) | wheel / 源码 / Cargo 镜像 |
| [COMPATIBILITY.md](COMPATIBILITY.md) | pysam API 矩阵 + v0.1.1 语义说明 |
| [BENCHMARKS.md](BENCHMARKS.md) | 16T 三方表、RSS、M8 recount |
| [CHANGELOG.md](CHANGELOG.md) | v0.1.0 / v0.1.1 变更 |
| [REVIEW.md](REVIEW.md) | v0.1 代码审查与 v0.1.1 整改清单 |
| [DEVELOPMENT.md](DEVELOPMENT.md) | 本地开发 |
| [NOTICE](NOTICE) | 行为参考（pysam / rubam / noodles / HTSlib），无大段抄代码 |
| [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md) | 里程碑 M0–M12 与设计约束 |

---

## 许可与致谢

MIT，见 [LICENSE](LICENSE)。

实现为独立代码。API 与 oracle 参考了 [pysam](https://github.com/pysam-developers/pysam)、[rubam](https://github.com/victormar1/rubam)、[noodles](https://github.com/zaeleus/noodles)、samtools / bcftools / HTSlib，详见 [NOTICE](NOTICE)。

---

## 状态

**v0.1.1** 在 v0.1.0（BAM 分析路径、VariantFile 读路径、CRAM 评估读路径、GitHub manylinux 发布）之上完成审查整改：`iter_batches` / `parallel_fetch` 正确性重写、统计并行 reader 复用、`reference_end` 与 pysam 风格 header、B-array tag 完整支持、坐标校验对齐 pysam。后续（PyPI、HTSlib backend、`FastaFile`、`pileup()` 迭代器）需另行提出。
