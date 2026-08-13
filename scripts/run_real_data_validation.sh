#!/usr/bin/env bash
# Manual real-data validation entrypoint (not for GitHub CI).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CONFIG="${1:-benchmark/configs/myceliophthora.local.yaml}"
if [[ ! -f "$CONFIG" ]]; then
  echo "Missing config: $CONFIG" >&2
  echo "Copy benchmark/configs/myceliophthora.yaml to a local file and fill paths." >&2
  exit 2
fi

# Prefer mt-var if present (pysam + samtools + bcftools).
if [[ -x /hpcfs/fhome/caizhh/.conda/envs/mt-var/bin/python ]]; then
  export PATH="/hpcfs/fhome/caizhh/.conda/envs/mt-var/bin:$PATH"
fi

python scripts/run_baselines.py --config "$CONFIG"
echo "Real-data validation scaffolding complete."
