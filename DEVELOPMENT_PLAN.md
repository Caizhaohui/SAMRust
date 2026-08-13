# SAMRust 详细系统开发计划

> 工作名：**SAMRust**（可后续改名）  
> 定位：**面向 Linux/HPC 的 Rust-native、多线程、高性能 pysam 兼容实现**  
> 参考项目：`pysam-developers/pysam`、[`victormar1/rubam`](https://github.com/victormar1/rubam)、`zaeleus/noodles`、HTSlib/samtools  
> **对比基线（强制）**：`pysam` + `rubam` — 关键功能须比较输出结果、运行时间、计算资源（见 §1.2 / §15 / §19.6）  
> **Benchmark 总原则**：**功能和结果对齐 pysam；运行效率对齐 rubam。** rubam 没有的功能或结果，在 16T 三方表对应单元格填 **NA**（禁止静默跳过）。pysam 是正确性 oracle 与慢对照，不是速度目标。  
> 主要真实测试数据：**嗜热毁丝霉（Myceliophthora thermophila）基因组重测序数据**  
> HPC 重负载测试默认递交队列：Slurm 分区 **`qcpu_18i`**（详见 §19.5）

---

# 1. 项目目标

## 1.1 核心目标

开发一个以 Rust 为核心、通过 PyO3 暴露 Python API 的高性能 HTS 数据访问与统计库，重点兼容 pysam 中最常用、最影响性能的接口，并通过 Rust 多线程、批处理、零拷贝/低拷贝数据通路提升 BAM/VCF 大规模处理效率。

项目不是简单“把 pysam 翻译成 Rust”，而是：

1. **保持 pysam 用户熟悉的 Python API 与语义**；
2. **把性能敏感的循环从 Python 移到 Rust**；
3. **提供真正可扩展的多线程 region/batch/pileup/depth 引擎**；
4. **默认保持确定性输出与严格正确性**；
5. **优先服务 Linux/HPC、WGS/WES/RNA-seq/微生物与真菌重测序场景**；
6. 对已有 pysam 脚本尽量做到低成本迁移；
7. **以 pysam 与 [rubam](https://github.com/victormar1/rubam) 为双基线**：关键 API 必须对比输出正确性、wall-clock 运行时间与计算资源（CPU / RSS / 线程效率）。

## 1.2 Benchmark 总原则（强制，不可省略）

**SAMRust 功能和结果对齐 pysam；运行效率对齐 rubam。**

这是贯穿验收、优化与发布的耐久规则，不是某次 benchmark 的临时注释。

| 维度 | 对齐对象 | 含义 |
|------|----------|------|
| **功能 / 正确性 / 输出** | **pysam**（oracle） | Python API 覆盖与语义、0-based half-open 坐标、默认过滤、计数结果以 pysam 为准。未解释的 mismatch = **正确性 gate 失败**。 |
| **性能 / 资源** | **[rubam](https://github.com/victormar1/rubam)** | 同一 workload、同一线程数下，wall-clock（及 RSS/CPU）以 rubam 为效率目标。**pysam 是慢对照基线，不是速度目标。** |

### 1.2.1 rubam 缺失时：三方表填 NA

rubam **没有**的功能，或 **没有**可对齐的结果字段时：

1. **不得**为对齐 rubam 而砍掉 pysam 已有功能，也 **不得**为对齐 rubam 输出而偏离 pysam oracle；
2. **不得**静默跳过该行；16T（及 1/4/8/16）三方表必须保留该 workload；
3. rubam 侧单元格一律写字面量 **`NA`**（wall_s / max_rss_kb / output_digest / gate / vs-rubam 加速比）；
4. notes / COMPATIBILITY / `benchmark/results/` 写明原因（例如「rubam 0.3.x 无 VariantFile」）；
5. SAMRust 对该功能的正确性 gate **只对 pysam**（及计划中的 samtools/bcftools）；**不得**因 rubam 为 NA 而标「性能验收完成」。

CSV 约定（`compare_pysam_rubam_samrust.*.csv`）：

```text
workload,threads,tool,wall_s,max_rss_kb,output_digest,gate_vs_pysam
variantfile,16,pysam,<sec>,<kb>,<digest>,oracle
variantfile,16,samrust,<sec>,<kb>,<digest>,match
variantfile,16,rubam,NA,NA,NA,NA
```

展示表约定：

| 负载 | SAMRust 16T | pysam 16T | rubam 16T | vs pysam | vs rubam |
|------|-------------|-----------|-----------|----------|----------|
| count / coverage / depth / pileup_counts | 实测 | 实测 | 实测 | 加速比 | 加速比 |
| VariantFile（例） | 实测（对 pysam） | 实测 | **NA** | 匹配/否 | **NA** |

已记录的「有入口但语义不同」仍用数字 + notes，不用 NA（例如缺失 MAPQ：pysam/SAMRust=255，rubam 可能为 0——报实测墙钟，输出差异写入 notes）。**NA 只用于 rubam 完全没有该功能或没有该结果列。**

已知语义例外必须写进计划 / `COMPATIBILITY.md` / `benchmark/results/`：

- **depth**（`depth_numpy` / `depth_blocks` / `get_depths`）：与 **samtools depth / rubam `get_depths`** 对齐，计 CIGAR `M/=/X`（含模糊碱基），**不计** `D`/`N`；这 **不等于** pysam `count_coverage` 的 A+C+G+T 之和。
- 缺失 MAPQ：SAMRust / pysam 为 **255**；rubam 可能为 0。不得为追平 rubam 而改 pysam 语义。

三方对比（SAMRust / pysam / rubam）对关键功能（§15.0）**强制**：有 rubam 入口则报输出 + runtime + RSS/CPU；无入口则该侧 **NA** + 原因。结果写入 `benchmark/results/`。

## 1.3 建议用户体验

```python
import samrust as pysam

bam = pysam.AlignmentFile("sample.bam", "rb")

for read in bam.fetch("chr1", 1000, 2000):
    ...
```

同时增加高性能扩展 API：

```python
bam = samrust.AlignmentFile("sample.bam", "rb", threads=16)

result = bam.count_coverage(
    "chr1",
    0,
    1_000_000,
    threads=16,
)

for batch in bam.iter_batches(batch_size=8192, threads=16):
    ...
```

---

# 2. 明确范围

## 2.1 第一阶段必须实现

### BAM/SAM

- `AlignmentFile`
- header 读取
- `AlignedSegment`
- 顺序 records iteration
- BAI/CSI indexed fetch
- `fetch()`
- `count()`
- `count_coverage()`
- `pileup()`
- depth/base-count 核心
- flags/CIGAR/sequence/quality/tags
- BAM index 检测与错误提示

### Python

- PyO3 bindings
- maturin 打包
- Linux wheels
- Python 3.10+ 优先
- GIL 释放
- NumPy 高效返回路径

### 多线程

- region scatter
- adaptive chunking
- bounded worker pool
- batch processing
- ordered result merge
- memory backpressure

### VCF/BCF

第一阶段后半程实现：

- `VariantFile`
- `VariantRecord`
- header/sample access
- sequential iteration
- indexed fetch

---

## 2.2 第一阶段明确不做

以下内容必须冻结，禁止 Codex/Grok 在早期里程碑擅自扩展：

- 不开发 read mapper
- 不开发 SNP/Indel caller
- 不开发 de novo assembler
- 不重写完整 samtools CLI
- 不重写完整 bcftools CLI
- 不开发 GUI
- 不支持 Windows 优先优化
- 不追求一次性 100% pysam API 覆盖
- 不在 M0-M7 阶段投入 CRAM 完整兼容
- 不把 tokio/async 引入本地 BAM 热路径，除非 benchmark 证明有收益

项目核心是：

> **pysam-compatible high-performance data access + parallel analytics**

而不是：

> samtools/bcftools replacement suite

---

# 3. 技术路线

## 3.1 默认后端

首选：

```text
Python
  │
  ▼
PyO3
  │
  ▼
SAMRust Rust API
  │
  ├── parallel scheduler
  ├── batch engine
  ├── pileup/depth engine
  │
  ▼
noodles
  │
  ▼
BAM / SAM / VCF / BCF / BGZF / FASTA
```

推荐依赖：

```toml
noodles = { version = "PINNED_VERSION", features = [
  "bam", "sam", "bgzf", "csi", "fasta", "vcf", "bcf", "tabix"
] }
rayon = "PINNED_VERSION"
crossbeam-channel = "PINNED_VERSION"
pyo3 = { version = "PINNED_VERSION", features = ["extension-module"] }
numpy = "PINNED_VERSION"
thiserror = "PINNED_VERSION"
smallvec = "PINNED_VERSION"
clap = { version = "PINNED_VERSION", features = ["derive"] }
```

要求：

- 创建项目时解析最新兼容版本后锁定 `Cargo.lock`；
- 不使用 `*` 或过宽依赖版本；
- 每次依赖升级单独 PR；
- 性能相关依赖升级必须重新 benchmark。

---

## 3.2 可选 HTSlib compatibility backend

不要在 M0–M11 实现。M11 只写设计；**不得**作为 BAM/VCF v0.1（M12）blocker。

后期若 noodles 无法覆盖真实 CRAM（非 3.0/3.1、缺参考、HTSlib-only codec、加密容器），可以增加 Cargo feature：

```text
--features htslib-backend
```

建议形态（未落地）：

- 新可选 crate 或 `cfg(feature)` 模块包装 `rust-htslib`（C HTSlib）。与 `samrust-core` 的 `#![deny(unsafe_code)]` 隔离，不把 `unsafe` 扩到 BAM stats 热路径。
- **默认 backend 仍为 noodles**（BAM / VCF / 本里程碑的 CRAM 读路径）。
- Feature 打开时：仅当 noodles 打不开该 CRAM 时 fallback；或用于与 HTSlib 行为对拍。
- 用途：CRAM compatibility fallback、与 HTSlib 行为差异验证、极端格式 edge cases。
- 不把 HTSlib 接到 `count` / depth / coverage / `pileup_counts`（那些路径继续走 `noodles::bam::Record` + `for_each_raw`）。

M11 评估：Tier-0 `small.cram`（samtools CRAM 3.x）上 noodles 顺序迭代 + CRAI fetch 已与 pysam 对拍，因此 **v0.1 不需要** 先实现该 feature。

---

# 4. 核心设计原则

## 4.1 pysam 语义优先

Python API 的语义以 pysam 为最终规范，而不是以内部 Rust 实现方便为准。

特别是：

### 坐标系统

所有 Python API：

```text
0-based
half-open
[start, stop)
```

例如：

```python
bam.fetch("chr1", 100, 200)
```

表示参考序列：

```text
100 <= position < 200
```

禁止在 Python 层同时出现 1-based inclusive API。

内部转换必须集中在：

```text
coords.rs
```

严禁各模块自行写 `+1/-1`。

---

## 4.2 Correctness before performance

任何性能优化都必须满足：

```text
correctness test pass（相对 pysam oracle）
↓
benchmark（相对 rubam 效率目标 + pysam 慢对照）
↓
optimization
↓
correctness re-test
↓
accept
```

禁止：

```text
先优化 → 再想办法解释结果差异
为追平 rubam 速度而改变 pysam 语义
```

效率目标见 §1.2 / §21：有 rubam 等价入口时尽量与 rubam 同量级；慢于 rubam 须 profiling，而不是改以 pysam 为速度上限。rubam 无入口则三方表填 NA，不做 vs-rubam 比值 gate。

---

## 4.3 Streaming first

避免把全基因组逐位点结果全部保存在内存。

默认设计：

```text
BAM
 ↓
chunk
 ↓
worker
 ↓
result block
 ↓
ordered merge
 ↓
Python/CLI consumer
```

提供两类接口：

### streaming

```python
for block in bam.depth_blocks(...):
    ...
```

### materialized

```python
positions, depth = bam.depth_numpy(...)
```

全量 materialization 是显式选择，不是默认行为。

---

## 4.4 Bounded memory

所有 producer/consumer channel 必须设置容量。

禁止：

```text
unbounded queue
```

建议初始：

```text
queue_capacity = 2 × worker_count
```

根据 benchmark 调整。

---

## 4.5 Deterministic output

默认：

```text
ordered=True
```

即使 region 并行处理，结果仍按：

```text
contig order
position order
input region order
```

返回。

可额外提供：

```python
ordered=False
```

用于最大吞吐场景。

---

# 5. 推荐仓库结构

```text
samrust/
├── Cargo.toml
├── Cargo.lock
├── pyproject.toml
├── README.md
├── LICENSE
├── CHANGELOG.md
├── DEVELOPMENT_PLAN.md
├── COMPATIBILITY.md
├── BENCHMARKS.md
│
├── crates/
│   ├── samrust-core/
│   │   └── src/
│   │       ├── lib.rs
│   │       ├── error.rs
│   │       ├── coords.rs
│   │       ├── config.rs
│   │       ├── header.rs
│   │       ├── record.rs
│   │       ├── cigar.rs
│   │       ├── tags.rs
│   │       ├── bam.rs
│   │       ├── index.rs
│   │       ├── fetch.rs
│   │       ├── batch.rs
│   │       ├── scheduler.rs
│   │       ├── count.rs
│   │       ├── coverage.rs
│   │       ├── pileup.rs
│   │       ├── variant.rs
│   │       └── fasta.rs
│   │
│   ├── samrust-python/
│   │   └── src/
│   │       ├── lib.rs
│   │       ├── alignment_file.rs
│   │       ├── aligned_segment.rs
│   │       ├── pileup.rs
│   │       ├── variant_file.rs
│   │       └── numpy_bridge.rs
│   │
│   └── samrust-cli/
│       └── src/main.rs
│
├── python/
│   └── samrust/
│       ├── __init__.py
│       ├── compat.py
│       └── typing.pyi
│
├── tests/
│   ├── unit/
│   ├── integration/
│   ├── parity/
│   ├── property/
│   ├── regression/
│   └── fixtures/
│
├── benchmark/
│   ├── scripts/
│   ├── configs/
│   ├── results/
│   └── reports/
│
├── scripts/
│   ├── discover_real_data.py
│   ├── prepare_fixture.py
│   ├── compare_pysam.py
│   ├── compare_rubam.py
│   ├── compare_samtools.py
│   ├── run_three_way_benchmark.sh
│   └── run_scaling_benchmark.sh
│
└── .github/workflows/
    ├── rust.yml
    ├── python.yml
    └── release.yml
```

注意：

- 早期不要拆更多 crates；
- `core` 必须完全不依赖 Python；
- Python binding 只负责对象转换与异常映射；
- 算法不得写在 PyO3 class method 内。

---

# 6. 核心 Rust API

## 6.1 AlignmentReader

初始设计：

```rust
pub struct AlignmentReader {
    // backend reader
    // parsed header
    // optional index
    // configuration
}
```

核心 API：

```rust
impl AlignmentReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self>;

    pub fn header(&self) -> &Header;

    pub fn records(&mut self) -> Result<RecordIter<'_>>;

    pub fn fetch(
        &self,
        region: &Region,
    ) -> Result<FetchIter>;

    pub fn count(
        &self,
        region: &Region,
        options: &CountOptions,
    ) -> Result<u64>;
}
```

---

## 6.2 Record

至少支持：

```text
query_name
flag
reference_id
reference_name
reference_start
reference_end
mapping_quality
cigar
query_sequence
query_qualities
mate_reference_id
mate_reference_start
template_length
aux tags
```

flag helper：

```text
is_paired
is_proper_pair
is_unmapped
mate_is_unmapped
is_reverse
mate_is_reverse
is_read1
is_read2
is_secondary
is_qcfail
is_duplicate
is_supplementary
```

---

# 7. Python pysam compatibility surface

## 7.1 AlignmentFile

目标：

```python
samrust.AlignmentFile(
    filename,
    mode="rb",
    reference_filename=None,
    threads=0,
)
```

第一阶段属性：

```text
filename
mode
header
references
lengths
nreferences
mapped
unmapped
```

第一阶段方法：

```text
fetch
count
count_coverage
pileup
close
reset
__iter__
__enter__
__exit__
```

---

## 7.2 AlignedSegment

要求尽量保持 pysam 属性命名：

```text
query_name
query_sequence
query_qualities
flag
reference_id
reference_name
reference_start
reference_end
mapping_quality
cigar
cigarstring
next_reference_id
next_reference_start
template_length
```

以及 flags helper。

---

## 7.3 性能扩展 API

这些不是 pysam parity API，而是 SAMRust extension：

```python
bam.iter_batches(
    batch_size=8192,
    threads=8,
    ordered=True,
)
```

```python
bam.parallel_fetch(
    regions,
    threads=16,
    ordered=True,
)
```

```python
bam.depth_blocks(
    contig,
    start,
    stop,
    block_size=1_000_000,
    threads=16,
)
```

```python
bam.pileup_counts(
    contig,
    start,
    stop,
    threads=16,
    min_base_quality=13,
    min_mapping_quality=0,
)
```

---

# 8. 多线程架构

# 8.1 不采用“多个线程共享一个 BAM reader”

初期原则：

> 每个 indexed worker 使用独立 reader handle。

原因：

- 简化线程安全；
- 避免 reader 内部锁竞争；
- region query 天然适合独立 handle；
- 更容易控制生命周期。

---

## 8.2 Region scatter

输入：

```text
contig:start-stop
```

转换为：

```text
Chunk 0
Chunk 1
Chunk 2
...
```

每个 chunk：

```rust
struct RegionChunk {
    id: usize,
    contig_id: usize,
    start: u64,
    stop: u64,
}
```

worker 输出：

```rust
struct ChunkResult<T> {
    id: usize,
    value: T,
}
```

merge 根据 `id` 重建顺序。

---

## 8.3 Adaptive chunking

不能固定按“线程数等分基因组”。

需要支持：

```text
min_chunk_size
max_chunk_size
target_chunks_per_thread
```

初始默认（**fetch / batch**）：

```text
target_chunks_per_thread = 4
```

即：

```text
8 threads → ~32 chunks
```

这样可降低不同区域 coverage 不均造成的负载不平衡。

**count / depth / coverage / pileup（M10）** 使用 `Scheduler::stats()`：

```text
target_chunks_per_thread = 1
```

证据：真菌 100 kb、16T、`chunks=4` 时 depth 从 8T 的 0.091 s **回退**到 16T 的 0.102 s（64 次 BAI query）。rubam fast mode 为 `n.div_ceil(nt)`（每线程 1 块）并在 ~0.058 s 平台。短窗口上额外分块的开销大于负载均衡收益。长 contig 的 fetch 仍用 4。

---

## 8.4 Chunk boundary correctness

这是 pileup/depth 最关键的工程问题之一。

read 可能跨 chunk 边界。

禁止简单：

```text
chunk1 0-1Mb
chunk2 1Mb-2Mb
query exactly chunk boundaries
```

然后直接拼接而不考虑 overlapping alignment。

正确策略：

- indexed query 允许返回与 region overlap 的 read；
- worker 只对属于其 `owned interval` 的 reference positions 输出；
- read 可跨边界，但 position ownership 不重叠。

定义：

```text
query interval = owned interval
output interval = owned interval
```

依赖 BAM indexed query overlap semantics。

所有 chunk-boundary 结果必须与单线程结果 bit-exact。

---

## 8.5 Producer / worker / merge

推荐结构：

```text
                     ┌──────── worker 0 ────────┐
                     │                          │
region scheduler ────┼──────── worker 1 ────────┼── ordered merge ── consumer
                     │                          │
                     └──────── worker N ────────┘
```

对于 batch sequential scan：

```text
BAM reader
   │
   ▼
record batch producer
   │
   ▼
bounded channel
   │
   ├── worker
   ├── worker
   ├── worker
   └── worker
```

---

## 8.6 GIL

所有长时间 Rust 计算：

```text
BAM decode
pileup
coverage
count
parallel fetch
```

必须在不持有 Python GIL 的状态下运行。

PyO3 binding 层仅在：

```text
input parse
Python object creation
exception conversion
result handoff
```

阶段持有 GIL。

CI 中增加 GIL concurrency regression test。

---

# 9. Batch Engine

## 9.1 目标

降低：

- Python/Rust FFI 调用次数；
- Python object creation；
- allocator pressure；
- cache miss。

---

## 9.2 Batch 类型

```rust
pub struct RecordBatch {
    records: Vec<Record>,
}
```

后续优化可改为 columnar：

```text
flags: Vec<u16>
positions: Vec<i64>
mapq: Vec<u8>
...
```

但 M0-M6 禁止过早 columnar rewrite。

先建立性能 baseline。

---

## 9.3 Python batch API

MVP 返回轻量 batch object：

```python
for batch in bam.iter_batches(batch_size=8192):
    flags = batch.flags
    positions = batch.reference_start
```

后期支持：

```text
NumPy arrays
PyArrow RecordBatch
Polars DataFrame
```

Arrow/Polars 不进入 v0.1 必选范围。

---

# 10. count / depth / coverage

## 10.1 count

第一原则：

必须逐项复制 pysam 默认过滤语义。

测试内容：

```text
secondary
supplementary
duplicate
QC fail
unmapped
read_callback
```

不要凭经验自行定义默认 filter。

---

## 10.2 count_coverage

目标 API：

```python
A, C, G, T = bam.count_coverage(
    contig,
    start,
    stop,
    quality_threshold=15,
    threads=16,
)
```

返回优先：

```text
NumPy uint32 arrays
```

兼容模式再转换成 pysam-like arrays/tuples。

---

## 10.3 depth

内部单独实现高性能 depth primitive：

```rust
DepthEngine
```

不要通过：

```text
pileup → 再求 depth
```

强制共用全部复杂 pileup 状态。

可复用 read/CIGAR traversal，但允许 depth 使用更轻的数据结构。

---

# 11. Pileup Engine

pileup 是项目核心技术壁垒之一。

## 11.1 第一阶段输出

```text
reference_pos
A
C
G
T
N
depth
```

## 11.2 第二阶段兼容 pysam PileupColumn

```text
reference_pos
nsegments
pileups
```

## 11.3 必须正确处理 CIGAR

```text
M
I
D
N
S
H
P
=
X
```

特别测试：

- insertion
- deletion
- splice/refskip
- soft clip
- supplementary alignment
- long indels

---

## 11.4 base quality

过滤次序与 pysam/samtools parity test 必须固定。

至少支持：

```text
min_base_quality
min_mapping_quality
max_depth
flag_filter
```

---

# 12. VCF/BCF 模块

M9 开始。

目标 API：

```python
vf = samrust.VariantFile("sample.vcf.gz")

for rec in vf.fetch("chr1", 0, 1_000_000):
    ...
```

优先支持：

```text
CHROM
POS
ID
REF
ALT
QUAL
FILTER
INFO
FORMAT
samples
GT
DP
AD
```

后期再增加 mutation/writer parity。

---

# 13. 嗜热毁丝霉真实测试数据设计

## 13.1 数据目录

项目测试脚本默认读取配置文件：

```yaml
reference_dir: /hpcfs/fhome/caizhh/18_WJX_work/Myceliophthora_thermophila_ATCC42464
raw_data_dir: /hpcfs/fhome/caizhh/18_WJX_work/X101SC260610181-Z01-J001
sample_metadata: /hpcfs/fhome/caizhh/18_WJX_work/sampleID_meta.csv
analysis_dir: /hpcfs/fhome/caizhh/18_WJX_work/02_Genetic_variation_diversity
```

保存为：

```text
benchmark/configs/myceliophthora.yaml
```

注意：

- 仓库不得硬编码 HPC 路径；
- 配置文件模板可提交；
- 本地真实路径加入 `.gitignore`；
- 测试代码必须允许命令行覆盖。

---

## 13.2 BAM 自动发现

创建：

```text
scripts/discover_real_data.py
```

行为：

1. 读取 `sampleID_meta.csv`；
2. 搜索 analysis_dir 中：

```text
*.bam
*.bam.bai
*.cram
*.crai
*.vcf.gz
*.bcf
```

3. 建立：

```text
benchmark/real_data_manifest.tsv
```

格式：

```text
sample_id	bam	index	vcf	condition	timepoint
```

4. 只发现，不修改原始数据；
5. 不自动重新比对全量 FASTQ；
6. 如果不存在 BAM，明确退出并提示需要先准备 coordinate-sorted indexed BAM。

---

# 14. 测试数据分级

## Tier 0：synthetic fixture

目标：

- CI
- edge case
- millisecond/second-level tests

包含：

```text
small.fa
small.bam
small.bam.bai
small.vcf.gz
small.vcf.gz.tbi
```

人工覆盖：

```text
M/I/D/N/S/H/=/X
paired/unpaired
duplicate
secondary
supplementary
QC fail
MAPQ 0/255
low BQ
indel
multi-allelic VCF
```

---

## Tier 1：真实 BAM 小区域

从嗜热毁丝霉任一代表样本选：

```text
1 contig
或
0.5-2 Mb region
```

目标：

- 快速 parity
- 开发阶段频繁运行

耗时目标：

```text
< 1 min
```

---

## Tier 2：单个样本全基因组

目标：

- full genome correctness
- thread scaling
- memory benchmark

测试：

```text
1 / 2 / 4 / 8 / 16 / 32 threads
```

**计算资源**：必须递交到 HPC 队列 **`qcpu_18i`**（见 §19.5），不要在登录节点执行。

---

## Tier 3：多个代表样本

从 metadata 自动选择：

- 不同培养温度；
- 不同时间点；
- 不同测序深度；
- 一个低 coverage；
- 一个中等 coverage；
- 一个较高 coverage。

不要写死 sample ID。

---

## Tier 4：全部重测序样本

只用于：

- release benchmark
- throughput
- HPC scaling
- stress test

不进入普通 CI。

**计算资源**：一律通过 **`qcpu_18i`** 队列递交（见 §19.5）。

---

# 15. Correctness Oracle 与对比基线

SAMRust 的正确性与性能验收采用 **双主基线 + 补充工具**（总原则见 §1.2）：

| 角色 | 工具 | 用途 |
|------|------|------|
| **语义 / API oracle（主）** | [pysam](https://github.com/pysam-developers/pysam) | **正确性对齐对象**：Python API 行为、坐标语义、默认过滤、bit-exact / normalized 输出 |
| **同品类竞品基线（主）** | [rubam](https://github.com/victormar1/rubam) | **效率对齐对象**：同为 Rust-native / PyO3 / noodles 路线；对比输出、runtime、资源。pysam 不作速度目标。 |
| 补充 CLI oracle | samtools / bcftools | `depth` / `mpileup` / `view -c` 等可比操作 |
| 可选诊断 | bam-readcount | 候选位点 allele support（非 bit-exact gate） |

> **规则**：凡列为“关键功能”的 API（见 §15.0），发布与里程碑 gate 必须同时给出相对 **pysam** 与 **rubam** 的：
> 1. **输出结果**（相对 pysam：一致 / 已记录的语义差；相对 rubam：可解释或已归档）；
> 2. **运行时间**（wall-clock，含 1T 与多线程；**效率目标 = 同条件 rubam**）；
> 3. **计算资源**（峰值 RSS、CPU 利用率或等效 `/usr/bin/time -v` 指标）。
>
> rubam **没有**等价功能或结果时：该功能不得标为“已完成性能验收”；可先完成 pysam 正确性 gate；**16T 三方表必须保留该行**，rubam 侧 `wall_s` / `max_rss_kb` / `output_digest` / `gate` / vs-rubam 一律填字面量 **`NA`**（见 §1.2.1），notes 写明原因。**禁止静默跳过。**

---

## 15.0 关键功能对比清单（强制）

以下功能必须纳入三方（SAMRust / pysam / rubam）对比。rubam 若某版本尚无等价 API：先尽量用最接近入口（如 `get_depths` / `pileup_bases` / `AlignmentFile.count`）做 apples-to-apples；**完全没有入口或没有可对齐结果列**时，16T 三方表 rubam 单元格填 **`NA`**（§1.2.1），不得删行。

| 功能 | pysam | rubam（典型入口） | 输出比较 | 时间 / 资源 |
|------|-------|-------------------|----------|-------------|
| indexed `fetch` / region iterate | `AlignmentFile.fetch` | `AlignmentFile.fetch` | record 身份集合 | ✅ |
| `count` | `AlignmentFile.count` | `count` / `count_reads` | 整数相等（统一 `read_callback`/flag） | ✅ |
| `count_coverage` | `count_coverage` | `count_coverage` / 等价 base 计数 | A/C/G/T 数组 | ✅ |
| depth / coverage profile | `count_coverage` 或 depth 脚本 | `get_depths` / `get_depths_numpy` | 每位点 depth | ✅ |
| pileup base counts | `pileup` 规范化 A/C/G/T/N/DP | `pileup` / `pileup_bases` | 规范化计数 | ✅ |
| candidate-site recount | pysam pileup 脚本 | rubam 等价 pileup/count | `ALT_COUNT>=10` 位点集合 | ✅ |
| parallel scaling（同功能） | 通常受 GIL 限制（作对照） | `num_threads` / 多线程 API | 1T==NT（SAMRust）+ 三方 runtime/RSS | ✅ |
| `VariantFile`（读路径） | `VariantFile` | **无等价 API** | pysam + bcftools | rubam 列 **NA** |

参数对齐要求：

- 坐标一律 **0-based half-open**（若 rubam 某 API 为 1-based inclusive，转换只允许在适配层）；
- 统一 `min_bq` / `min_mapq` / flag 过滤（或显式记录差异）；
- 同一 BAM、同一 region 列表、同一节点、同一文件系统层级（shared vs scratch 分开报）。

---

## 15.1 pysam

主要验证：

```text
AlignmentFile
AlignedSegment
fetch
count
count_coverage
pileup
VariantFile
```

角色：

- **Python API 语义最终规范**（正确性对齐对象，见 §1.2）；
- 默认过滤 / BQ / MAPQ 行为以 pysam 文档与实测为准；
- SAMRust 与 pysam 的 unexplained mismatch = **正确性 gate 失败**；
- pysam 的 runtime 只作为慢对照记录，**不是** SAMRust 的速度目标。

---

## 15.2 rubam（https://github.com/victormar1/rubam）

`rubam` 是与 SAMRust 最接近的开源竞品：纯 Rust + noodles、PyO3 绑定、多线程 depth/pileup/count，并宣称相对 pysam/samtools 的 bit-exact 校验。

主要对比面：

```text
AlignmentFile / AlignedSegment
fetch / count / count_coverage / pileup
get_depths / get_depths_numpy / pileup_bases
count_reads / flag_stats
VariantFile — 无等价 API → 16T 三方表 rubam 列填 NA
```

规则：

1. **输出**：在参数对齐后，关键计数类结果应与 pysam **同时**可解释；若 SAMRust≡pysam 但与 rubam 不一致，必须归因（过滤默认、坐标、del/refskip、quality 阈值、缺失 MAPQ 等），写入 `benchmark/results/` 差异报告，不得静默忽略。**不得为对齐 rubam 而改 pysam 语义。**
2. **运行时间**：**效率对齐对象是 rubam**（仅当 rubam **有**等价入口时）。同一 workload 记录 `SAMRust` / `pysam` / `rubam` 的 wall-clock（1T 与目标线程数）。目标是尽量与 rubam 同量级；pysam 仅作慢对照。
3. **计算资源**：至少记录峰值 RSS（`/usr/bin/time -v` 或等价）；多线程任务另记 CPU% / 效率 `speedup/n`。资源目标同样对照 rubam。
4. **版本钉扎**：benchmark 元数据必须记录 `rubam` 版本（pip/git commit）。
5. 安装与脚本：`scripts/compare_rubam.py` + `scripts/run_three_way_benchmark.py` / `.sh` + `scripts/submit_three_way_benchmark.sh`；重负载递交 `qcpu_18i`（§19.5）。
6. **无等价功能或结果**：保留该 workload 行；rubam 侧 `wall_s` / `max_rss_kb` / `output_digest` / `gate` / vs-rubam 填字面量 **`NA`**；notes 写明原因（例如「rubam 0.3.x 无 VariantFile」）。禁止静默跳过、禁止删行。有入口但语义不同仍报实测数字 + notes，不用 NA。

---

## 15.3 samtools

主要验证：

```text
view
view -c
depth
mpileup
idxstats
```

只比较双方明确具有相同语义的结果。

---

## 15.4 bcftools

验证：

```text
VCF/BCF parse
indexed region records
core fields
sample genotype
```

---

## 15.5 bam-readcount

如果环境存在：

用于真实候选位点 A/C/G/T/indel support 的补充验证。

注意：

bam-readcount 与 pysam pileup 默认参数/过滤语义并不完全相同。

因此：

- 必须显式统一参数；
- 未统一时只做诊断，不作为 bit-exact gate。

---

# 16. 嗜热毁丝霉实际研究场景验收

增加一个 domain-specific benchmark：

```text
candidate-site recount
```

输入：

```text
candidate VCF/BED
+
BAM
```

输出：

```text
sample
chrom
pos
ref
A
C
G
T
N
indel_support
depth
alt_support
allele_frequency
```

实际验收：

```text
ALT supporting reads >= 10
```

得到的有效位点集合必须与 **pysam** reference implementation 完全一致；在参数可对齐时，同时与 **rubam** 等价 pileup/count 结果交叉验证，并记录三方 runtime / RSS（见 §15.0 / §19.6）。

这不是 v0.1 variant caller，而是 pileup/count API 的真实应用测试。

---

# 17. Differential Testing

创建：

```text
tests/parity/test_random_regions.py
```

每个真实 BAM：

1. 读取 contig length；
2. 固定 random seed；
3. 随机生成 1000 个 region；
4. region length 覆盖：

```text
1 bp
10 bp
100 bp
1 kb
10 kb
100 kb
```

比较：

```text
fetch record identities
count
coverage arrays
pileup counts
```

工具：

```text
pysam
rubam（参数对齐后的等价 API）
samrust
```

输出 mismatch report：

```text
sample
region
metric
pysam
rubam
samrust
```

release gate：

```text
vs pysam: 0 unexplained mismatches
vs rubam: 0 unexplained mismatches（或已归档的语义差说明）
```

---

# 18. Property Testing

使用 `proptest`。

重点测试：

### coordinate conversion

```text
0-based half-open ↔ internal coordinate
```

### chunk split

任意：

```text
start < stop
threads >= 1
```

必须满足：

```text
union(chunks) == original interval
intersection(chunk_i, chunk_j) == empty
```

### merge

```text
parallel result == serial result
```

---

# 19. Benchmark Framework

## 19.1 工具

Rust microbenchmark：

```text
criterion
```

end-to-end：

```text
hyperfine
/usr/bin/time -v
perf stat
```

可选：

```text
flamegraph
heaptrack
```

---

## 19.2 所有 benchmark 必须记录

```text
date
git commit
CPU model
CPU count
RAM
Linux kernel
filesystem
storage type
Rust version
Python version
pysam version
rubam version（pip 或 git commit）
samtools version
samrust version
sample
BAM size
read count
thread count
tool（samrust|pysam|rubam）
wall_clock_s
max_rss_kb
cpu_percent（如可测）
```

---

## 19.6 三方性能与资源对比（pysam + rubam + SAMRust）

总原则（§1.2）：**功能和结果对齐 pysam；运行效率对齐 rubam。** pysam 是正确性 oracle 与慢对照，不是速度目标。rubam 没有的功能或结果，16T 三方表对应单元格填 **NA**。

关键功能（§15.0）的每个正式 benchmark 必须产出可比表格，写入：

```text
benchmark/results/compare_pysam_rubam_samrust.<workload>.json
benchmark/results/compare_pysam_rubam_samrust.<workload>.csv
```

最小列：

```text
workload, region_or_sites, threads,
tool, wall_s, max_rss_kb, cpu_percent,
output_digest_or_gate, notes
```

要求：

1. **同一节点 / 同一队列作业**内连续跑完三方（避免跨天机器噪声）；重负载用 `qcpu_18i`。
2. **输出**：先做正确性比对，再报时间；输出不一致时不得只报“更快”。
3. **时间**：至少 warm-up 1 + measured 3（正式 release 按 §19.4 用 5）；报告 median。
4. **资源**：峰值 RSS 必报；多线程另报 speedup 与 parallel_efficiency。
5. **公平性**：
   - pysam 多线程若受 GIL 限制，仍作为对照如实记录；
   - rubam / SAMRust 使用各自推荐的多线程参数，线程数集合对齐（如 1/2/4/8/16）；
   - **rubam `AlignmentFile.count_coverage` 不传 `threads`**：该入口使用 rubam 默认内部线程，三方脚本保持这一调用方式，避免把“默认多线程 rubam”当成 SAMRust 1T 的不公平对照；depth / pileup 等显式 `num_threads` 入口则按请求线程数对齐；
   - 禁止把不同过滤默认当成同输入硬比速度；
   - depth 与 `count_coverage` 不得混用 oracle（depth = samtools/rubam M/=/X，不是 A+C+G+T）。
6. 脚本入口（规划 / 维护）：
   - `scripts/compare_pysam.py`
   - `scripts/compare_rubam.py`
   - `scripts/run_three_way_benchmark.sh`（递交 `qcpu_18i`）
7. rubam 无等价功能或结果时：**不得删行、不得静默跳过**。CSV / 展示表 rubam 侧填字面量 **`NA`**（`wall_s`、`max_rss_kb`、`output_digest`、`gate_vs_pysam`、vs-rubam）；notes 写原因。SAMRust 正确性仍只对 pysam（及 samtools/bcftools）gate。效率对齐 **仅** 适用于 rubam 有等价入口的负载；NA 行不得标“性能验收完成”。

---

## 19.3 HPC 文件系统噪声

必须区分：

### shared filesystem benchmark

实际生产场景。

### local scratch benchmark

用于判断 CPU/算法性能。

如果节点有 `$TMPDIR` / local NVMe：

复制代表 BAM fixture 到 local scratch 后测试。

禁止只报告 shared filesystem 单次结果。

---

## 19.4 重复次数

至少：

```text
warm-up = 1
measured runs = 5
```

报告：

```text
median
min
max
```

不要只报告 best-of-one。

---

## 19.5 计算资源与作业队列（HPC）

开发与测试过程中，**真实 BAM / 全基因组 / 多线程 scaling / 真菌 recount** 等任务需要大量 CPU 与长时间运行，**禁止在登录节点上直接跑重负载**（易被系统 SIGKILL，也会干扰他人）。

### 默认队列

| 项 | 值 |
|----|----|
| 调度系统 | Slurm |
| 分区 / 队列 | **`qcpu_18i`** |
| 适用任务 | Tier 1–4 真实数据、M7 fungal pileup、M8 recount/`ALT_COUNT>=10` gate、1–32 线程 scaling、release benchmark |

登录节点仅允许：

- Tier 0 fixture / 单元测试 / 短 pytest
- `cargo check` / `clippy` / 小范围编译
- 作业脚本编写、`sbatch`/`squeue` 管理

### 提交方式

使用 Slurm 分区 `qcpu_18i`（本集群示例）：

```bash
sbatch --partition=qcpu_18i \
  --job-name=samrust-bench \
  --cpus-per-task=16 \
  --mem=32G \
  --time=08:00:00 \
  --output=benchmark/results/slurm_%j.out \
  --error=benchmark/results/slurm_%j.err \
  scripts/run_m7_m8_heavy.sh
```

或在脚本头写入：

```bash
#!/usr/bin/env bash
#SBATCH --partition=qcpu_18i
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=benchmark/results/slurm_%j.out
#SBATCH --error=benchmark/results/slurm_%j.err
```

注意：`qcpu_18i` 节点一般为 **24 CPU / 64G**；不要申请超过单节点可满足的 `cpus-per-task`（例如 32 会被拒绝）。thread scaling 默认测到 16；若需 32 线程，改用可提供 ≥32 CPU 的分区，或跨作业对比。

资源申请按任务调整：短区域 parity 可用较少核；thread scaling 至少申请覆盖本次最大线程数的 CPU。

### Agent / 开发约定

1. 识别到将跑真菌全量 BAM、候选位点 recount、多线程 scaling 时，**先写可递交脚本，再 `sbatch -p qcpu_18i`**，不要在交互会话里硬跑。
2. 作业内环境需自行激活 conda / 设置 `PATH`（含 `samrust` release 二进制与 `mt-var` 等依赖）。
3. 结果仍写入 `benchmark/results/` 或节点 `$TMPDIR`/scratch；**只读**原始重测序数据，禁止改写。
4. 若队列不可用或路径缺失：停止该真实数据 benchmark，继续 synthetic/unit，**不得伪造结果**。

---

# 20. Thread Scaling Benchmark

运行：

```text
threads = 1, 2, 4, 8, 16, 32
```

计算：

```text
speedup(n) = T1 / Tn
parallel_efficiency(n) = speedup(n) / n
```

分别 benchmark（每项均含 pysam / rubam / SAMRust，见 §19.6）：

```text
count
count_coverage
depth
pileup_counts
parallel_fetch
candidate recount
```

输出：

```text
benchmark/results/scaling.csv
benchmark/results/compare_pysam_rubam_samrust.*.csv
```

---

# 21. 性能目标

这些是工程目标，不允许为了达标牺牲正确性。总原则见 §1.2：**功能和结果对齐 pysam，效率对齐 rubam**（rubam 无入口则 16T 表填 NA）。

## v0.1 target

### Correctness

```text
pysam parity: 0 unexplained mismatch
rubam comparable ops: 0 unexplained mismatch（或已文档化的语义差）
samtools comparable operations: 0 unexplained mismatch
```

### single thread

目标：

```text
效率：尽量与 rubam 同量级（同 workload / 同线程）
pysam：仅作慢对照，不是速度目标；SAMRust 不应慢于 pysam
相对 rubam：记录比值；若慢于 rubam > 1.5×，必须 profiling 并说明原因
```

初始警戒线：

```text
SAMRust 1T > 1.25 × pysam runtime  → 必须 profiling（已慢于慢对照，不可接受）
SAMRust 1T > 1.50 × rubam runtime → 必须 profiling / 优化 backlog（未达到效率目标）
```

### multi-thread

在 Tier 2 fungal BAM 上：

```text
8T vs SAMRust 1T
目标 speedup >= 2×
```

对 pileup/depth/count_coverage：

```text
8T SAMRust vs pysam reference
目标 >= 2×（pysam 为慢对照）
8T SAMRust vs rubam（同线程数）
必须报告比值；**效率目标 = 尽量与 rubam 同量级**（警戒：> 1.5× 需解释）
rubam 无等价入口：该行 vs-rubam / wall / RSS 填 NA，不套用 1.5× 警戒
```

如果文件系统 I/O 导致未达到，必须同时给 local scratch 结果。

### memory

必须报告 RSS（SAMRust / pysam / rubam 三方同表）。

禁止以数倍内存交换速度而不说明。

---

# 22. Profiling Strategy

优化前运行：

```bash
perf stat ...
perf record ...
```

分析：

```text
BGZF decode
BAM parsing
CIGAR traversal
allocation
Python conversion
result merge
lock contention
filesystem wait
```

每次性能 PR 必须回答：

1. bottleneck 是什么？
2. flamegraph/perf 证据是什么？
3. 修改了什么？
4. correctness 是否保持？
5. 1/2/4/8/16T 是否改善？
6. RSS 是否增加？

---

# 23. 里程碑

# M0 — Repository Initialization / Scope Freeze

## 目标

建立可维护仓库骨架。

## 任务

- Cargo workspace
- samrust-core
- samrust-python
- samrust-cli
- rustfmt
- clippy
- pytest
- maturin
- GitHub Actions Linux CI
- README
- LICENSE
- DEVELOPMENT_PLAN
- COMPATIBILITY matrix skeleton

## 验收

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --workspace
python -m pytest
```

全部通过。

## 范围冻结

M0 禁止实现：

```text
BAM algorithms
pileup
VCF
parallel processing
```

---

# M1 — Test Infrastructure & Real Data Discovery

## 目标

先建立 oracle 与 benchmark 系统。

## 任务

- synthetic BAM fixture
- synthetic VCF fixture
- `discover_real_data.py`
- real_data_manifest.tsv generation
- pysam baseline scripts
- samtools baseline scripts
- benchmark metadata collector

## 真实数据

发现嗜热毁丝霉 BAM/BAI。

自动建立 Tier 1 / Tier 2 数据集。

## 验收

能生成：

```text
benchmark/results/pysam_baseline.json
benchmark/results/samtools_baseline.json
```

---

# M2 — Rust BAM Core

## 目标

实现 Rust 线性 BAM 访问。

## 功能

- BAM open
- header
- records iterator
- Record
- flags
- sequence
- quality
- CIGAR
- aux tags

## 测试

synthetic fixture + fungal BAM first 100k reads。

对比 pysam：

```text
qname
flag
reference id
position
mapq
cigar
sequence length
selected tags
```

## Gate

```text
0 mismatches
```

---

# M3 — Python AlignmentFile / AlignedSegment

## 目标

建立最小 pysam-compatible Python surface。

## API

```python
AlignmentFile
AlignedSegment
records iteration
header
references
lengths
```

## GIL

线性 Rust decode loop 需要提供 batch path 并释放 GIL。

## 验收

已有 Python parity suite 可以同时：

```python
import pysam
import samrust
```

比较输出。

---

# M4 — Indexed Fetch

## 目标

支持 BAI/CSI indexed query。

## API

```python
bam.fetch(contig, start, stop)
```

## 必须处理

```text
missing index
invalid contig
empty interval
start == stop
region at contig start
region at contig end
read crossing interval boundary
secondary/supplementary records
```

## 测试

fungal BAM random 1000 regions。

## Gate

```text
record identity parity == 100%
```

---

# M5 — Parallel Runtime & Batch Engine

## 目标

建立后续所有性能功能共用的并行基础设施。

## 实现

```text
RegionChunk
Scheduler
WorkerPool
BoundedChannel
OrderedMerger
RecordBatch
```

## API

Rust：

```text
parallel_map_regions
```

Python：

```text
iter_batches
parallel_fetch
```

## 验收

所有 thread count：

```text
1/2/4/8/16
```

结果与 serial identical。

线程数变化不得改变输出。

---

# M6 — count / depth / count_coverage

## 目标

先实现最容易准确 benchmark 的核心热点。

## API

```text
count
count_coverage
depth_blocks
depth_numpy
```

## 实现重点

- adaptive chunking
- streaming output
- NumPy output
- max memory bound

## 真实数据 benchmark

Tier 1 + Tier 2。

比较：

```text
pysam
samtools depth
SAMRust 1/2/4/8/16/32T
```

---

# M7 — Parallel Pileup

## 目标

实现项目第一核心竞争功能。

## 功能

```text
A/C/G/T/N
DP
BQ filter
MAPQ filter
flag filter
indel-aware traversal foundation
```

## 测试

synthetic CIGAR edge cases。

真实 fungal BAM：

```text
random regions
full genome
candidate sites
```

## Gate

单线程与 **pysam** / samtools normalized output 无差异。

在参数对齐后，与 **rubam** 等价 pileup/base-count 输出交叉验证（不一致须归档语义差）。

多线程与 SAMRust 单线程 bit-exact。

正式性能报告须含 pysam + rubam + SAMRust 的 runtime / RSS（§19.6）。

---

# M8 — Fungal Resequencing Domain Validation

## 目标

使用真实研究问题验证工具。

## 输入

- real BAMs
- candidate VCF/BED

## 实现 benchmark utility

```text
samrust-recount
```

注意：该 CLI 仅作为 benchmark/reference utility，不扩展成 variant caller。

## 输出

```text
sample
chrom
pos
ref
A
C
G
T
N
DP
ALT_COUNT
AF
```

## Gate

`ALT_COUNT >= 10` 的有效候选集合与 **pysam** reference 完全一致。

在可对齐时，与 **rubam** 等价计数交叉验证；并记录三方 runtime / RSS / 线程 scaling（§19.6）。

然后对代表样本运行：

```text
1/2/4/8/16/32 threads
```

记录 runtime/RSS/CPU（SAMRust；并在同作业中采集 pysam / rubam 对照）。

---

# M9 — VariantFile

## 目标

实现 pysam VariantFile 常用**读取**接口。

## 功能

- VCF
- VCF.gz（+ tabix `.tbi`）
- BCF（+ CSI）
- header（samples / contigs）
- sample names
- record fields：CHROM / POS / ID / REF / ALT / QUAL / FILTER / INFO / FORMAT / samples（GT、DP、AD）
- indexed fetch（Python **0-based half-open**；`VariantRecord.pos` 为 pysam 式 **1-based**）

不包含 writer、mutation、variant caller。

## Oracle

```text
pysam.VariantFile
bcftools view -H
```

pytest：`tests/parity/test_m9_variantfile_parity.py`。

**16T 三方表：rubam = NA。** rubam 0.3.x 无 `VariantFile` 等价入口（plan §15.2）；不得静默跳过、不得删行。正确性只对 pysam + bcftools。M9 不是性能里程碑：SAMRust / pysam 墙钟可后续补测，rubam 列在有等价 API 之前保持 **NA**。

## 验收（2026-08-13）

- 顺序迭代与 `fetch` 区域集合对 pysam 一致（Tier-0 `small.vcf.gz` / `small.vcf` / `small.bcf`）。
- 记录数与 `bcftools view -H` 一致。
- 未知 contig 报错；空区间 `fetch` 为空。

---

# M10 — Performance Optimization Wave

## 目标

只基于 profiling 做优化。BAM **功能和结果对齐 pysam；运行效率对齐 rubam**。

候选：

- allocation reuse
- SmallVec CIGAR
- pre-sized buffers
- pooled worker buffers
- faster ASCII/base classification
- branch reduction
- cache-friendly count arrays
- result block compression
- Python zero-copy arrays

禁止无 benchmark 的“感觉型优化”。

## 本波已落地（2026-08-13）

Profiling 证据（作业 2312392，Mt-35 100 kb，release，`for_each_raw` 之后）：

- 1T：count / depth / pileup 已快于 rubam；coverage 的 rubam 行内部默认 ~4 线程，不可直接当 1T。
- 16T：depth 0.102 s vs rubam 0.058 s（**1.78× 慢**，超过 §21 的 1.5× 警戒）；pileup 1.56×；depth 从 8T 0.091 s **回退**到 16T。
- 根因：统计路径 `target_chunks_per_thread=4` → 16T 打开 **64** 次 BAM+BAI；rubam fast mode 为每线程 1 块。次因：CIGAR M/=/X 内层对每个碱基做区间判断。

实施（不改 pysam 语义）：

1. `Scheduler::stats()`：count / depth / coverage / pileup 每线程 1 块；fetch / `iter_batches` 仍为 4（§8.3）。
2. `Interval::overlap_span` 裁剪 CIGAR 操作；depth 对重叠切片 `+= 1`；coverage / pileup 只扫重叠碱基。
3. A/C/G/T 分类改 256 项 LUT，去掉热循环里的 `to_ascii_uppercase`。

未做（无新证据不引入依赖）：SmallVec CIGAR（统计路径已 `cigar().iter()`）、整 BAM worker 拷贝、改 rubam `count_coverage` 传 threads。

## 验收（2026-08-13）

- cargo test / clippy `-D warnings` / pytest 31 passed（含 coverage 1T==NT）。
- 作业 **2312423**（`qcpu_18i` / bnode2 / release）：16T 四项 vs pysam **match**。
- 16T vs-rubam：count 4.11×、coverage 1.41×、depth 1.17×、pileup 1.03×（均 ≥ 1；不再慢于 rubam 1.5×）。
- VariantFile 行 rubam = **NA**。结果：`BENCHMARKS.md` 与 `benchmark/results/compare_pysam_rubam_samrust.fungal_mt35.*`。

---

# M11 — CRAM Evaluation

## 目标

评估：

```text
noodles CRAM readiness
```

如果真实兼容性不足：

设计 optional HTSlib backend。

CRAM 不得拖延 BAM/VCF v0.1 release。

## 评估结论（2026-08-13）

noodles **0.115** + `noodles-cram` **0.98.0** + `noodles-fasta` **0.66.0**（现有 `noodles` crate 的 `cram`/`fasta` feature，不是新 crate）：

| 能力 | 就绪？ | 说明 |
|------|--------|------|
| 顺序迭代 | 是 | `cram::io::indexed_reader` + `records(&sam_header)` |
| CRAI 区域查询 | 是 | `IndexedReader::query` + `Query::records()` |
| 参考序列 | 是 | `fasta::io::indexed_reader` + `fasta::Repository`（`.fai`） |
| CRAM 3.0 / 3.1 | 是 | noodles-cram 文档；fixture 由 samtools 1.23 `view -C` 生成 |
| 与 BAM stats 热路径共用 | **否** | CRAM 记录是 `sam::alignment::RecordBuf`，不是 `noodles::bam::Record`；`for_each_raw` 不能直接走 |
| 完整 HTSlib CRAM 边角 | 未证 | 加密、非 3.x、缺参考、HTSlib-only codec 不在本里程碑 |

因此 M11 **产品边界**是评估级读路径，不是完整 CRAM 引擎：

- Python `AlignmentFile`：`.cram` + `.crai` + FASTA；mode `rb` 或 `rc`；`reference_filename=` 或同目录 `{stem}.fa/.fasta/.fna`
- 顺序迭代与 indexed `fetch`（Python **0-based half-open**）对 pysam `AlignmentFile(..., "rc")` 对拍
- 空区间 `[start, start)`：SAMRust 返回空（与 pysam **BAM** 一致）。pysam CRAM/HTSlib 仍可能吐出覆盖该点的记录；**不复制该行为**
- **不做** `count` / `count_coverage` / `depth_*` / `pileup_counts` / `iter_batches` / `parallel_fetch`（明确 `NotImplemented`）
- 顺序路径物化全部记录（fixture 规模）；fetch 另开 reader，避免打乱迭代状态
- 16T 三方表：本里程碑不加 CRAM 行。若后续补测，rubam 无 CRAM API 则 rubam 列 **NA**

## Optional HTSlib backend（仅设计，不实现）

见 §3.2。M11 **不**引入 `rust-htslib`、**不**打开 `--features htslib-backend`。BAM/VCF v0.1（M12）不以 CRAM 完整兼容为门。

## Oracle

```text
pysam.AlignmentFile(..., "rc", reference_filename=...)
```

pytest：`tests/parity/test_m11_cram_parity.py`。Rust：`cram::tests` vs `small.bam` qname。

## 验收（2026-08-13）

- cargo test --workspace / clippy `-D warnings`
- pytest 含 M11 CRAM 顺序字段、fetch 区域、stats `RuntimeError`、BAM 仍可不传 reference
- fixture：`scripts/prepare_fixture.py` → `tests/fixtures/small.cram[.crai]`

---

# M12 — v0.1 Release

## Release gate

### Tests

```text
cargo test
pytest
property tests
parity tests
real fungal Tier 2 test
```

### Correctness

```text
0 unexplained parity mismatches
```

### Performance report

必须生成：

```text
BENCHMARKS.md
```

包含：

```text
pysam
samtools
SAMRust 1/2/4/8/16/32T
runtime
RSS
speedup
CPU utilization
```

### Package

```bash
pip install samrust
```

Linux x86_64 wheel。

## 验收（2026-08-13）

- 文档门：`README.md`（第一屏）、`INSTALL.md`、`API_COMPATIBILITY.md` → `COMPATIBILITY.md`、`BENCHMARKS.md`、`DEVELOPMENT.md`、`NOTICE`。
- Property tests：`proptest` on `Interval`（0-based half-open；dev-only）。
- Parity：Tier-0 pytest（含 M11 CRAM）+ 既有 M8 真菌 `ALT_COUNT >= 10` gate（6895/6895）。
- 16T 三方表：作业 2312423（`BENCHMARKS.md` + `benchmark/results/compare_pysam_rubam_samrust.fungal_mt35.*`）。线程 **1/4/8/16**；**无 32T**（`qcpu_18i` ~24 核）；**无独立 CPU% 采样**（记录 wall + max RSS）。M8 recount 含 2T。
- Package：GitHub Release **manylinux x86_64** wheels（`.github/workflows/release.yml` on `v*`）。v0.1 **不**发布 PyPI；安装见 `INSTALL.md`。仓库：https://github.com/Caizhaohui/SAMRust

---

# 24. CI 设计

普通 CI 只运行 Tier 0。

## Rust CI

```text
fmt
clippy
unit tests
property tests
```

## Python CI

```text
maturin build
install wheel
pytest parity fixture
```

## 不进入 GitHub CI

```text
real fungal BAM
large BAM
HPC benchmark
```

真实测试通过：

```text
scripts/run_real_data_validation.sh
```

手动执行。

---

# 25. Compatibility Matrix

创建：

```text
COMPATIBILITY.md
```

格式：

| pysam API | SAMRust | parity | multithread | notes |
|---|---|---|---|---|
| AlignmentFile | ✅ | tested | n/a | |
| fetch | ✅ | tested | ✅ | |
| count | ✅ | tested | ✅ | |
| count_coverage | ✅ | tested | ✅ | |
| pileup | partial | tested subset | ✅ | |
| VariantFile | ✅ | pysam + bcftools | n/a M9 | 读路径；16T 三方表 rubam 列 **NA**（无等价 API） |
| FastaFile | planned | | | |
| TabixFile | planned | | | |

任何新 API 必须更新此表。

---

# 26. Error Handling

Rust：

```rust
#[derive(thiserror::Error, Debug)]
pub enum SamRustError {
    ...
}
```

Python 映射：

```text
FileNotFoundError
ValueError
OSError
IndexError / custom SAMRustError where appropriate
```

要求：

- Rust panic 不允许穿过 FFI；
- malformed BAM 不得导致 Python process abort；
- worker panic 必须捕获并转换为错误；
- 一个 worker 错误后其余 worker 应尽快取消。

---

# 27. Logging

Rust core 使用：

```text
tracing
```

但默认静默。

调试：

```bash
SEQHTS_LOG=debug
```

性能 benchmark 禁止 debug log。

---

# 28. API 稳定性

v0.1 前：

```text
experimental
```

v0.1 后：

pysam-compatible API 尽量稳定。

SAMRust-specific parallel API 可以在：

```text
samrust.experimental
```

先迭代。

---

# 29. 文档

必须包含：

```text
README.md
INSTALL.md
API_COMPATIBILITY.md
BENCHMARKS.md
DEVELOPMENT.md
```

README 第一屏只回答：

1. SAMRust 是什么？
2. 为什么比 pysam 快？
3. 是否兼容 pysam？
4. 如何安装？
5. 最小示例。

---

# 30. 许可证与参考代码规则

参考项目包括 pysam 与 rubam。

开发原则：

- 优先参考 API、行为、测试思路和架构；
- 不应让 AI agent 大段复制第三方源码；
- 如果直接移植 substantial code，必须核对许可证并保留 attribution；
- 建议项目 LICENSE 选用 MIT 或 Apache-2.0/MIT dual license 前进行最终确认；
- `NOTICE` 中记录直接复用代码来源。

AI Agent 每次参考 rubam/pysam 代码实现时，在 PR 描述中注明：

```text
Reference inspected:
- repository
- file
- commit
- whether code was copied or only behavior/architecture referenced
```

---

# 31. AI Coding Agent 执行规则

以下规则必须放入：

```text
AGENTS.md
```

## Rule 1

严格按 milestone 执行。

不得提前做后续 milestone。

---

## Rule 2

每个 milestone 开始前：

```text
read DEVELOPMENT_PLAN.md
read current code
read current tests
run baseline tests
```

---

## Rule 3

每个 milestone 必须：

```text
implementation
+
tests
+
documentation
```

三者同时完成。

---

## Rule 4

性能代码没有 benchmark 不允许 merge。

---

## Rule 5

任何 parallel implementation 必须首先证明：

```text
1-thread output == N-thread output
```

---

## Rule 6

任何 pysam-compatible method 必须有：

```text
pysam oracle test
```

关键 BAM analytics（`fetch` / `count` / `count_coverage` / depth / `pileup_counts` / candidate recount）还必须与 [rubam](https://github.com/victormar1/rubam) 比较 **输出、runtime、RSS/CPU**，结果写入 `benchmark/results/`（§15 / §19.6）。

Benchmark 总原则（§1.2）：**功能和结果对齐 pysam；运行效率对齐 rubam。** pysam 不是速度目标。rubam 无等价功能或结果时，16T 三方表 rubam 单元格填 **NA** 并写明原因，禁止静默跳过。

---

## Rule 7

不要擅自改变坐标语义。

Python API 永远：

```text
0-based half-open
```

---

## Rule 8

不要在 Python loop 中实现热点逻辑。

热点必须放到 Rust。

---

## Rule 9

不要为了多线程把整个 BAM/全基因组结果复制到每个 worker。

---

## Rule 10

任何新 dependency 必须说明理由。

---

## Rule 11

不要修改用户原始重测序数据。

real-data test：

```text
read-only
```

输出进入项目 benchmark/results 或单独 scratch。

---

## Rule 12

如果真实数据路径不存在：

```text
STOP real-data benchmark
```

但继续完成 synthetic/unit work。

不得伪造真实 benchmark 结果。

---

# 32. 每个 Milestone 的 Codex/Grok 提交模板

```markdown
## Milestone
M5 — Parallel Runtime

## Scope
- ...

## Implemented
- ...

## Tests Added
- ...

## Correctness
- pysam comparison: ...
- serial vs parallel: ...

## Benchmark
- dataset: ...
- threads: ...
- runtime: ...
- RSS: ...

## Known Limitations
- ...

## Files Changed
- ...

## Next Milestone
M6
```

---

# 33. 推荐首次 Codex/Grok 指令

项目初始化时不要直接把整份计划一次性要求 agent 全做完。

首条任务建议：

```text
执行 DEVELOPMENT_PLAN.md 中 M0 — Repository Initialization / Scope Freeze。

严格限制在 M0：
- 初始化 Cargo workspace
- 创建 samrust-core / samrust-python / samrust-cli
- 配置 rustfmt/clippy/test
- 配置 maturin/pytest
- 配置 Linux GitHub Actions
- 创建 README、LICENSE、COMPATIBILITY.md、AGENTS.md

禁止实现：
- BAM parser
- fetch
- pileup
- depth
- VCF
- parallel runtime

完成后运行：
- cargo fmt --check
- cargo clippy --all-targets --all-features -- -D warnings
- cargo test --workspace
- python -m pytest

最后报告：
1. 创建的目录结构
2. 测试结果
3. 当前依赖版本
4. 是否存在 blocker
5. 不要开始 M1
```

M0 验收后，再单独发送 M1。

---

# 34. 建议开发顺序总结

```text
M0  repository
 ↓
M1  tests + real-data discovery + baselines
 ↓
M2  Rust BAM core
 ↓
M3  Python AlignmentFile
 ↓
M4  indexed fetch
 ↓
M5  parallel runtime + batch engine
 ↓
M6  count/depth/count_coverage
 ↓
M7  parallel pileup
 ↓
M8  fungal resequencing validation
 ↓
M9  VariantFile
 ↓
M10 profiling-driven optimization
 ↓
M11 CRAM evaluation
 ↓
M12 v0.1
```

---

# 35. v0.1 成功定义

v0.1 不需要“替代整个 pysam”。

如果达到以下目标即可认为项目成功：

1. Linux 上 `pip install` 可使用；
2. 提供稳定的 `AlignmentFile` / `AlignedSegment`；
3. `fetch/count/count_coverage/pileup` 与 **pysam** 核心语义一致（正确性），并与 **rubam** 完成可比输出交叉验证与 **效率对齐**（同条件 wall-clock / RSS）；rubam 无等价入口的功能（如 `VariantFile`）在 16T 三方表填 **NA**；
4. 嗜热毁丝霉真实 BAM 上没有未解释的 correctness mismatch（相对 pysam；相对 rubam 的差异已说明）；
5. pileup/depth/count_coverage 能有效利用 4-16 个 CPU cores；
6. 线程扩展不明显增加错误或内存失控；
7. candidate-site recount 中 `ALT_COUNT >= 10` 的位点集合与 pysam reference 一致，并完成 rubam 交叉验证；
8. 真实数据 benchmark 清楚展示 **SAMRust / pysam / rubam** 的 runtime、RSS、CPU utilization 与 scaling（功能和结果对齐 pysam；效率对齐 rubam）。**效率对齐仅适用于 rubam 有等价入口的负载**；无入口或无可对齐结果时，16T 三方表 rubam 列填 **NA**；
9. Python 用户可以用很小修改把性能关键代码从 pysam 迁移到 SAMRust。

最终项目定位：

> **SAMRust is a Rust-native, multi-threaded, pysam-compatible HTS processing library optimized for large-scale sequencing analysis on Linux/HPC systems. Correctness tracks pysam; performance tracks rubam.**

