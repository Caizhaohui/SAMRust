#!/usr/bin/env bash
# Submit M7/M8 heavy validation to Slurm partition qcpu_18i.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTDIR="${ROOT}/benchmark/results"
mkdir -p "$OUTDIR" "${ROOT}/benchmark/logs"

JOB_ID=$(sbatch --parsable \
  --partition=qcpu_18i \
  --job-name=samrust-m7m8 \
  --cpus-per-task=16 \
  --mem=32G \
  --time=08:00:00 \
  --chdir="${ROOT}" \
  --output="${ROOT}/benchmark/logs/m7m8_%j.out" \
  --error="${ROOT}/benchmark/logs/m7m8_%j.err" \
  --export=ALL,SAMRUST_ROOT="${ROOT}",OUTDIR="${OUTDIR}",M8_THREADS="${M8_THREADS:-1:2:4:8:16}" \
  "${ROOT}/scripts/run_m7_m8_heavy.sh")

echo "Submitted job ${JOB_ID} on partition qcpu_18i"
echo "  stdout: ${ROOT}/benchmark/logs/m7m8_${JOB_ID}.out"
echo "  stderr: ${ROOT}/benchmark/logs/m7m8_${JOB_ID}.err"
echo "Monitor: squeue -j ${JOB_ID} ; tail -f ${ROOT}/benchmark/logs/m7m8_${JOB_ID}.out"
