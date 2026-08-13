# SAMRust v0.1 代码审查与整改计划（2026-08-13）

本文档总结对 v0.1 代码库的全面审查结果，并跟踪整改状态。整改以 **v0.1.1**
发布，分四级优先级：P0 正确性、P1 性能、P2 兼容性、P3 流程。

审查方法：通读 `samrust-core` / `samrust-python` / `samrust-cli` 全部源码，
对照 [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md) 的语义规范，并用探针脚本在
Tier-0 fixture 上逐条验证（pysam 0.23.3 为行为 oracle）。

## P0 — 正确性（已确认 bug）

| # | 问题 | 位置 | 状态 |
|---|------|------|------|
| P0a | `iter_batches` 1T 只返回第一批后 EOF；MT 路径哈希 dedup 误删同 qname 记录、忽略 `batch_size`、全量物化 | `alignment.rs` | 已修 |
| P0b | `VariantFile.fetch(contig)` 在 header 无 contig length 时返回空（`stop.unwrap_or(0)`） | `variant.rs` / `vcf.rs` | 已修 |
| P0c | `close()` 后访问 `references` / `lengths` / `header` 触发 `PanicException`（`expect`） | `alignment.rs` | 已修 |
| P0d | 顺序迭代中途重新 `iter(bam)` 静默丢弃预取缓冲中的记录 | `alignment.rs` | 已修 |

## P1 — 性能

| # | 问题 | 位置 | 状态 |
|---|------|------|------|
| P1a | recount / `parallel_map_regions` 每个 chunk / 每个位点重开 BAM+索引（16T 下 recount 慢 5.6×） | `parallel.rs` / `recount.rs` | 已修（`map_init` 每线程一个 reader） |
| P1b | `__next__` 批解码持 GIL；`__next__` / fetch 迭代器逐条 `clone()` Record | `alignment.rs` / `variant.rs` | 已修（`allow_threads` + `IntoIter`） |
| P1c | `parallel_fetch` 哈希 dedup 每条记录分配 2 个 String 且有碰撞风险 | `parallel.rs` | 已修（区域合并 + 位置归属过滤，无哈希） |
| — | fetch / VCF 全量物化（流式 fetch 需自持 reader 的自引用结构） | — | 推迟到 v0.2，内存特性已文档化 |

## P2 — pysam 兼容性

| # | 问题 | 状态 |
|---|------|------|
| P2a | 缺 `AlignedSegment.reference_end`；`header` 不是 pysam 结构（HD/SQ/RG/PG） | 已修 |
| P2b | B 型 array tag 被 stringify 成 Rust Debug 文本 | 已修（`list[int]` / `list[float]`，与 pysam `array.array` 逐元素相等） |
| P2c | 负坐标抛 `OverflowError`（应为 `ValueError`）；`stop` 超出 contig 长度时行为与 pysam 不一致 | 已修（统一 `ValueError`；stop clamp 到 contig 长度，start 越界返回空） |
| P2d | `pileup()` 列式迭代器推迟、无索引 VCF fetch 回退、depth 语义差异未记录 | 已记录于 COMPATIBILITY.md |

### 探针证伪项（不改）

- **`query_name`**：pysam 对 QNAME=`*` 返回字符串 `"*"` 而非 `None`，SAMRust 现状一致。
- **stop clamp**：pysam `fetch("chr1", 0, 5000)`（contig 长 1000）静默截断返回 11 条；
  `count("chr1", 1500, 5000)` 返回 0。SAMRust 对齐为 clamp + 空结果。
- **`reference_end`**：pysam 对 unmapped（含 placed-unmapped）返回 `None`。

## P3 — 流程 / 基础设施

| # | 问题 | 状态 |
|---|------|------|
| P3a | §17 随机区域差分测试缺失（M7 曾靠它抓到过 1T/MT 不一致） | 已补 `tests/parity/test_random_regions.py`（fixture 常驻 + 真实数据可选） |
| P3b | 三方 benchmark 同进程测 RSS（三工具共享高水位）、digest 依赖各自序列化格式 | 已改为子进程独立采 RSS + 规范化 digest |
| P3c | 根目录 `-` 垃圾文件；`SAMRust_DEVELOPMENT_PLAN.md` 旧计划；`config.rs` / `BoundedChannel` / `compat.py` 死代码；`BASE_BUCKET` 重复定义；CI 只测 py3.10/3.12；pytest 缺 oracle 时静默 skip | 已清理 |

## 边界说明

- 坐标语义不变：Python API 永远 0-based half-open（Rule 7）。
- 统计路径（count/coverage/depth/pileup）算法未动，仅修调度与接口层。
- `parallel_fetch` 语义微调：重叠区域现在先合并再切块，结果仍是区域并集、
  每记录恰好一次，但完全相同的重复记录（cat-bam 场景）不再被误删。
