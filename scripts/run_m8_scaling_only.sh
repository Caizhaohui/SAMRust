#!/usr/bin/env bash
# M8 thread scaling only (gate already passed). Run on qcpu_18i.
set -euo pipefail

ROOT="${SAMRUST_ROOT:-/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/02_SAMRust}"
cd "$ROOT"
export PATH="/hpcfs/fhome/caizhh/.conda/envs/mt-var/bin:${ROOT}/target/release:${PATH}"

BAM="${BAM:-/hpcfs/fhome/caizhh/18_WJX_work/02_Genetic_variation_diversity/05_align/Mt-35-15d-1/Mt-35-15d-1.markdup.bam}"
SITES="${SITES:-/hpcfs/fhome/caizhh/18_WJX_work/02_Genetic_variation_diversity/07_candidates_union/Mt-35.candidates.bed}"
SAMPLE="${SAMPLE:-Mt-35-15d-1}"
OUTDIR="${OUTDIR:-${ROOT}/benchmark/results}"
mkdir -p "$OUTDIR"

THREADS_SPEC="${M8_THREADS:-1:2:4:8:16}"
IFS=':,' read -r -a THREADS <<< "${THREADS_SPEC// /:}"

echo "[$(date -Is)] host=$(hostname) cpus=${SLURM_CPUS_PER_TASK:-?} scaling=${THREADS_SPEC}"

python -u - <<PY
import json, resource, subprocess, time
from pathlib import Path

bam = Path("${BAM}")
sites = Path("${SITES}")
sample = "${SAMPLE}"
outdir = Path("${OUTDIR}")
threads = [int(x) for x in "${THREADS_SPEC}".replace(" ", ":").replace(",", ":").split(":") if x]
scaling = []
for t in threads:
    out_tsv = outdir / f"m8_recount_t{t}.tsv"
    print(f"scaling threads={t} ...", flush=True)
    cmd = [
        "samrust", "recount",
        "--bam", str(bam),
        "--sites", str(sites),
        "--sample", sample,
        "--threads", str(t),
        "--output", str(out_tsv),
    ]
    t0 = time.perf_counter()
    rss0 = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    rss1 = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    row = {
        "threads": t,
        "elapsed_s": elapsed,
        "ru_maxrss_kb": max(rss0, rss1),
        "output": str(out_tsv),
        "stderr": proc.stderr.strip(),
    }
    scaling.append(row)
    print(f"  {elapsed:.2f}s rss_kb={row['ru_maxrss_kb']}", flush=True)

gate_path = outdir / "m8_fungal_validation.json"
payload = {}
if gate_path.is_file():
    payload = json.loads(gate_path.read_text())
payload["scaling"] = scaling
payload["scaling_note"] = "updated by run_m8_scaling_only.sh"
gate_path.write_text(json.dumps(payload, indent=2) + "\n")
(outdir / "m8_scaling.json").write_text(json.dumps({"scaling": scaling}, indent=2) + "\n")
print(f"wrote {gate_path} and m8_scaling.json", flush=True)
PY

echo "[$(date -Is)] M8 scaling finished OK"
