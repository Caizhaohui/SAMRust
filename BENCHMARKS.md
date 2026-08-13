# BENCHMARKS

功能和结果对齐 **pysam**；运行效率对齐 **rubam**。rubam 没有的功能或结果，16T 三方表填 **NA**（`DEVELOPMENT_PLAN.md` §1.2.1）。

Published artifact (v0.1): `benchmark/results/compare_pysam_rubam_samrust.fungal_mt35.{csv,json}` (job **2312423**).

CPU utilization was not sampled as a separate counter; the recorded resource columns are **wall-clock** and **max RSS**. Thread set is **1/4/8/16**. **32T was not run**: `qcpu_18i` nodes have ~24 cores (do not request 32+ `cpus-per-task`). **2T** was not in this three-way campaign (M8 recount scaling did include 2T; see below).

## M10 16T 三方表（作业 2312423，现行）

- 分区 `qcpu_18i` / **bnode2** / Xeon Silver 4116
- `maturin develop --release`；M10：`Scheduler::stats()`（每线程 1 块）+ CIGAR `overlap_span` + base LUT
- `Mt-35-15d-1.markdup.bam`，`NC_016472.1:0-100000`，min_bq=0
- threads 1/4/8/16，repeats=3 median
- 四项 vs pysam **match**（count=111076）

16T 墙钟（秒）；加速比 = 对照 / SAMRust（>1 表示 SAMRust 更快）：

| 负载 | SAMRust | pysam | rubam | vs pysam | vs rubam |
|------|---------|-------|-------|----------|----------|
| count | 0.033 | 0.118 | 0.135 | 3.59× | **4.11×** |
| count_coverage | 0.076 | 4.472 | 0.106 | 59.1× | **1.41×** |
| depth | 0.041 | 1.022 | 0.048 | 25.0× | **1.17×** |
| pileup_counts | 0.090 | 50.729 | 0.093 | 562× | **1.03×** |
| VariantFile | — | — | **NA** | — | **NA** |

SAMRust 墙钟 scaling（秒）：

| 负载 | 1T | 4T | 8T | 16T |
|------|----|----|----|-----|
| count | 0.064 | 0.035 | 0.027 | 0.033 |
| count_coverage | 0.150 | 0.088 | 0.074 | 0.076 |
| depth | 0.076 | 0.045 | 0.038 | 0.041 |
| pileup_counts | 0.171 | 0.103 | 0.088 | 0.090 |

1→8T SAMRust：count 2.34×、coverage 2.02×、depth 1.98×、pileup 1.94×。16T 相对 8T 持平（100 kb 窗口）。coverage 的 1T 行 vs-rubam=0.71× **不是**公平单核对照（rubam `count_coverage` 不传 `threads`）；4T 为 1.20×。

16T max RSS（KB，benchmark 进程 `max_rss_kb`）：

| 负载 | SAMRust | pysam | rubam |
|------|---------|-------|-------|
| count | 134988 | 134988 | 134988 |
| count_coverage | 134988 | 134988 | 134988 |
| depth | 134988 | 134988 | 134988 |
| pileup_counts | 142952 | 134988 | 134988 |
| VariantFile | — | — | **NA** |

Identical RSS on several rows is the harness process sample, not a per-library allocator trace.

## M8 fungal recount gate

`ALT_COUNT >= 10` site set vs pysam: **6895 / 6895**, `gate_ok=true` (15 453 candidate sites). SAMRust recount wall / RSS:

| threads | wall_s | max RSS KB |
|---------|--------|------------|
| 1 | 282.7 | 52252 |
| 2 | 150.0 | 65684 |
| 4 | 82.3 | 116616 |
| 8 | 47.5 | 191568 |
| 16 | 27.5 | 347948 |

## Pre-M10 baseline（作业 2312392，bnode29）

同一 BAM / 区域 / release / `for_each_raw`，但统计路径仍为每线程 4 块。主机不同，不能当严格 A/B。

| 负载 | SAMRust 16T | vs rubam |
|------|-------------|----------|
| count | 0.084 | 2.00× |
| count_coverage | 0.135 | 0.89× |
| depth | 0.102 | 0.56×（8T→16T 回退） |
| pileup_counts | 0.158 | 0.64× |
| VariantFile | — | **NA** |

M10 证据：16T 打开 64 次 BAI query；rubam fast mode 为 `n.div_ceil(nt)`。改为每线程 1 块后，16T 四项 vs-rubam 均 ≥ 1.03×，且不再出现 depth 8T→16T 回退。
