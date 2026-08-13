# SAMRust Compatibility Matrix

Status legend:

- `planned` — not started
- `partial` — subset implemented
- `✅` — implemented for the stated milestone scope
- `parity` — compared against pysam / rubam / samtools / bcftools where noted

| pysam API | SAMRust | parity | multithread | notes |
|-----------|---------|--------|-------------|-------|
| AlignmentFile | ✅ | pysam | n/a | M3；rubam 交叉验证见 plan §15 |
| AlignmentFile CRAM | partial | pysam `"rc"` | n/a | M11；顺序迭代 + indexed `fetch`（`.cram`+`.crai`+FASTA）；mode `rb`/`rc`；stats/parallel **NotImplemented**；空区间与 pysam BAM 一致（pysam CRAM 可能非空）；16T 表无 CRAM 行（rubam 无 CRAM 则日后 **NA**） |
| AlignedSegment | ✅ | pysam | n/a | M3 |
| header / references / lengths | ✅ | pysam | n/a | M2–M3 |
| AlignedSegment fields (Rust `Record`) | ✅ | pysam | n/a | M2 |
| fetch | ✅ | pysam | ✅ `parallel_fetch` | M4 / M5；Python 仍物化 `AlignedSegment`；rubam 输出+runtime 对比 required |
| count | ✅ | pysam | ✅ | M6；1T==NT（含 placed unmapped）；`for_each_raw`；M10：`Scheduler::stats()` 1 chunk/thread；rubam 对比 required |
| count_coverage | ✅ | pysam | ✅ | M6；`for_each_raw`；M10：区间裁剪 + base LUT；rubam 对比 required |
| depth_blocks / depth_numpy | ✅ | pysam CIGAR oracle / samtools / rubam `get_depths` | ✅ | M6；M/=/X only（不含 D/N）；含 N 碱基；≠ `count_coverage` 之和；`for_each_raw`；M10：切片 `+= 1` |
| pileup_counts | ✅ | pysam-normalized | ✅ | M7；`for_each_raw`；M10：区间裁剪 + base LUT；vs rubam pileup/pileup_bases required |
| `samrust recount` | ✅ | ALT≥10 vs pysam | ✅ site-parallel | M8；rubam 交叉验证 + 三方 runtime/RSS required |
| VariantFile | ✅ | pysam + bcftools `view -H` | n/a (M9 read path) | M9；VCF / VCF.gz / BCF；`fetch` 0-based half-open；`pos` 1-based；无 writer；16T 三方表 rubam 列 **NA**（无等价 API） |
| FastaFile | planned | | | post-v0.1 |
| TabixFile | planned | | | post-v0.1 |
| `samrust.__version__` | ✅ | n/a | n/a | M0 |
| Tier-0 fixtures + oracle baselines | ✅ | pysam/samtools scripts | n/a | M1；扩展 rubam baselines per §15 |

正式性能验收须产出 `benchmark/results/compare_pysam_rubam_samrust.*`（见 `DEVELOPMENT_PLAN.md` §1.2 / §19.6）：**功能和结果对齐 pysam；运行效率对齐 rubam。** rubam 没有的功能或结果，16T 三方表对应单元格填 **NA**（禁止静默跳过）。Python `fetch` / `AlignedSegment` 仍物化 owned record；`count` / coverage / depth / `pileup_counts` 在 Rust 内流式访问 noodles 记录。统计路径并行为每线程 1 个索引查询（`Scheduler::stats`）；fetch 仍为每线程约 4 块。

## v0.1.1 语义说明

- **坐标校验**：`fetch` / `count` / `count_coverage` / `depth_*` / `pileup_counts` 的
  `start` / `stop` 为负时抛 `ValueError`（与 pysam 一致，不再是 `OverflowError`）；
  `stop` 超出 contig 长度时静默截断到 contig 末尾（pysam 行为）；`start` 越界返回空结果。
- **`header`**：`AlignmentFile.header` 返回 pysam 风格 dict（`HD` / `SQ` / `RG` / `PG`），
  不再是单行文本。`SQ` 条目含 `SN` / `LN`。
- **`reference_end`**：`AlignedSegment.reference_end` 为 0-based exclusive 末端；
  unmapped（含 placed-unmapped）或无 CIGAR 时返回 `None`，与 pysam 一致。
- **B 型 array tag**：`get_tag` 对 `B,c/B,C/B,s/B,S/B,i/B,I` 返回 `array.array`
  （typecode 保留），对 `B,f` 返回 `array.array('f')`；与 pysam 逐元素、逐类型相等。
- **`iter_batches`**：1T 为真流式（不再只返回首批）；MT 按 1 Mb 波次处理且
  尊重 `batch_size`；两种模式都包含文件末尾的 unmapped-unplaced 记录；
  1T 与 NT 输出逐条一致。完全相同的重复记录（cat-bam 场景）不再被去重误删。
- **重迭代**：顺序迭代中途重新 `iter(af)` 从逻辑位置继续（pysam 语义），
  不再静默丢弃预取缓冲中的记录。
- **关闭后访问**：`close()` 后访问 `header` / `references` / `lengths` /
  `nreferences` 抛 `ValueError`（不再是 `PanicException`）。
- **VCF 无 length contig**：`VariantFile.fetch(contig)`（`stop=None`）在 header
  缺少 contig length 时回退顺序扫描，不再返回空；显式 `stop` 仍走索引。
- **无索引 VCF**：`VariantFile` 打开无 `.tbi`/`.csi` 的 VCF.gz 时 `fetch` 不可用
  （抛 `ValueError`）；请先用 `pysam.tabix_index` / `bcftools index` 建索引。
  纯文本 VCF 与 BCF+CSI 不受影响。
- **`pileup()` 列式迭代器**：仍未实现（`pileup_counts` 为数组式替代），
  推迟到 v0.2。
- **depth 语义**：`depth_numpy` / `depth_blocks` 统计 CIGAR M/=/X（含 N 碱基，
  不含 D/N），与 samtools depth / rubam `get_depths` 一致；**不等于**
  `count_coverage` 的 A+C+G+T（后者还受 base-quality 阈值影响）。
- **内存特性**：`fetch` / `VariantFile.fetch` 全量物化结果（流式 fetch 需
  自持 reader 的自引用结构，推迟到 v0.2）；大区域请用 `iter_batches` 或分窗。
