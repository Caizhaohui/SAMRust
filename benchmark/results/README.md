# Benchmark results

Most files in this directory are **gitignored** (local HPC / login-node outputs).

v0.1 publishes the M10 16T three-way table:

- `compare_pysam_rubam_samrust.fungal_mt35.csv`
- `compare_pysam_rubam_samrust.fungal_mt35.json`

Narrative: [`BENCHMARKS.md`](../../BENCHMARKS.md).

Regenerate locally:

```bash
python scripts/prepare_fixture.py
python scripts/run_baselines.py --skip-real
# fungal / 16T: submit to qcpu_18i, do not use the login node
bash scripts/submit_three_way_benchmark.sh
```
