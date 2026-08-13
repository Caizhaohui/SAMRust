#!/usr/bin/env bash
# Submit three-way pysam/rubam/SAMRust benchmark to Slurm partition qcpu_18i.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTDIR="${ROOT}/benchmark/results"
mkdir -p "$OUTDIR" "${ROOT}/benchmark/logs"

MODE="${MODE:-fungal}"
THREADS="${THREADS:-1:4:8:16}"
REPEATS="${REPEATS:-3}"

JOB_ID=$(sbatch --parsable \
  --partition=qcpu_18i \
  --job-name=samrust-3way \
  --cpus-per-task=16 \
  --mem=32G \
  --time=04:00:00 \
  --chdir="${ROOT}" \
  --output="${ROOT}/benchmark/logs/threeway_%j.out" \
  --error="${ROOT}/benchmark/logs/threeway_%j.err" \
  --export=ALL,SAMRUST_ROOT="${ROOT}",OUTDIR="${OUTDIR}",MODE="${MODE}",THREADS="${THREADS}",REPEATS="${REPEATS}" \
  "${ROOT}/scripts/run_three_way_benchmark.sh")

echo "Submitted job ${JOB_ID} on partition qcpu_18i (MODE=${MODE})"
echo "  stdout: ${ROOT}/benchmark/logs/threeway_${JOB_ID}.out"
echo "  stderr: ${ROOT}/benchmark/logs/threeway_${JOB_ID}.err"
echo "Monitor: squeue -j ${JOB_ID}"
