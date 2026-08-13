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
