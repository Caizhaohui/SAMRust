#!/usr/bin/env bash
# Three-way benchmark entrypoint: pysam + rubam + SAMRust.
# Login node: Tier-0 fixture only.
# Real fungal BAM / large regions: submit via submit_three_way_benchmark.sh (qcpu_18i).
set -euo pipefail

ROOT="${SAMRUST_ROOT:-/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/02_SAMRust}"
cd "$ROOT"
if [[ ! -f "${ROOT}/Cargo.toml" ]]; then
  echo "SAMRust root not found at $ROOT. Export SAMRUST_ROOT." >&2
  exit 2
fi

export PATH="/hpcfs/fhome/caizhh/.conda/envs/mt-var/bin:${ROOT}/target/release:${ROOT}/target/debug:${PATH}"
export PYTHONPATH="${ROOT}/python${PYTHONPATH:+:${PYTHONPATH}}"

OUTDIR="${OUTDIR:-${ROOT}/benchmark/results}"
mkdir -p "$OUTDIR" "${ROOT}/benchmark/logs"

MODE="${MODE:-fixture}"   # fixture | fungal
THREADS="${THREADS:-1:4:8}"
REPEATS="${REPEATS:-3}"
MIN_BQ="${MIN_BQ:-0}"

echo "[$(date -Is)] host=$(hostname) mode=${MODE} threads=${THREADS} cpus=${SLURM_CPUS_PER_TASK:-?}"

if ! python -c "import rubam" 2>/dev/null; then
  echo "Installing rubam ..."
  pip install --user 'rubam>=0.3.13'
fi

# Three-way timings are meaningless in debug. Always rebuild release so
# fungal 16T tables pick up the current hot path (for_each_raw, etc.).
echo "[$(date -Is)] Building samrust Python extension (release) ..."
maturin develop --release -m crates/samrust-python/Cargo.toml

CONTIG_ARGS=()
if [[ "$MODE" == "fixture" ]]; then
  BAM="${BAM:-${ROOT}/tests/fixtures/small.bam}"
  CONTIG_ARGS=(--contig "${CONTIG:-chr1}")
  START="${START:-0}"
  STOP="${STOP:-200}"
  TAG="${TAG:-fixture}"
elif [[ "$MODE" == "fungal" ]]; then
  BAM="${BAM:-/hpcfs/fhome/caizhh/18_WJX_work/02_Genetic_variation_diversity/05_align/Mt-35-15d-1/Mt-35-15d-1.markdup.bam}"
  if [[ -n "${CONTIG:-}" ]]; then
    CONTIG_ARGS=(--contig "$CONTIG")
  fi
  START="${START:-0}"
  STOP="${STOP:-100000}"
  TAG="${TAG:-fungal_mt35}"
else
  echo "Unknown MODE=$MODE (use fixture|fungal)" >&2
  exit 2
fi

if [[ ! -f "$BAM" ]]; then
  echo "Missing BAM: $BAM" >&2
  exit 2
fi

echo "[$(date -Is)] rubam baseline ..."
python -u scripts/compare_rubam.py \
  --bam "$BAM" \
  "${CONTIG_ARGS[@]}" \
  --start "$START" \
  --stop "$STOP" \
  --threads 1 \
  --min-bq "$MIN_BQ" \
  --output "${OUTDIR}/rubam_baseline.${TAG}.json" \
  --repo-root "$ROOT"

echo "[$(date -Is)] pysam baseline ..."
python -u scripts/compare_pysam.py \
  --bam "$BAM" \
  "${CONTIG_ARGS[@]}" \
  --start "$START" \
  --stop "$STOP" \
  --output "${OUTDIR}/pysam_baseline.${TAG}.json" \
  --repo-root "$ROOT"

echo "[$(date -Is)] three-way benchmark ..."
python -u scripts/run_three_way_benchmark.py \
  --bam "$BAM" \
  "${CONTIG_ARGS[@]}" \
  --start "$START" \
  --stop "$STOP" \
  --threads "$THREADS" \
  --min-bq "$MIN_BQ" \
  --repeats "$REPEATS" \
  --tag "$TAG" \
  --outdir "$OUTDIR" \
  --repo-root "$ROOT"

echo "[$(date -Is)] three-way benchmark finished OK"
ls -la "${OUTDIR}/compare_pysam_rubam_samrust.${TAG}".* "${OUTDIR}/rubam_baseline.${TAG}.json"
