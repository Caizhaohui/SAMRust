#!/usr/bin/env bash
# Heavy M7/M8 real-data validation — run on compute node (qcpu_18i), not login.
set -euo pipefail

# Slurm may copy this script into /var/spool; never derive ROOT from BASH_SOURCE alone.
ROOT="${SAMRUST_ROOT:-/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/02_SAMRust}"
cd "$ROOT"
if [[ ! -f "${ROOT}/Cargo.toml" ]]; then
  echo "SAMRust root not found at $ROOT (pwd=$PWD). Export SAMRUST_ROOT." >&2
  exit 2
fi

export PATH="/hpcfs/fhome/caizhh/.conda/envs/mt-var/bin:${ROOT}/target/release:${PATH}"
export PYTHONPATH="${ROOT}/python${PYTHONPATH:+:${PYTHONPATH}}"

BAM="${BAM:-/hpcfs/fhome/caizhh/18_WJX_work/02_Genetic_variation_diversity/05_align/Mt-35-15d-1/Mt-35-15d-1.markdup.bam}"
SITES="${SITES:-/hpcfs/fhome/caizhh/18_WJX_work/02_Genetic_variation_diversity/07_candidates_union/Mt-35.candidates.bed}"
SAMPLE="${SAMPLE:-Mt-35-15d-1}"
OUTDIR="${OUTDIR:-${ROOT}/benchmark/results}"
export BAM SITES SAMPLE OUTDIR
mkdir -p "$OUTDIR"

echo "[$(date -Is)] host=$(hostname) cpus=${SLURM_CPUS_PER_TASK:-?} pwd=$ROOT"

# Ensure release binary + Python extension are available on this node.
if [[ ! -x "${ROOT}/target/release/samrust" ]]; then
  echo "Building release samrust ..."
  cargo --config 'source.crates-io.replace-with="ustc"' \
    --config 'source.ustc.registry="sparse+https://mirrors.ustc.edu.cn/crates.io-index/"' \
    build -p samrust-cli --release
fi
if [[ ! -f "${ROOT}/python/samrust/_samrust.cpython-312-x86_64-linux-gnu.so" ]]; then
  echo "Building Python extension ..."
  maturin develop -m crates/samrust-python/Cargo.toml
fi

echo "[$(date -Is)] === M7 fungal random-region pileup ==="
python -u scripts/run_m7_fungal_pileup.py \
  --bam "$BAM" \
  --n-regions "${M7_N_REGIONS:-40}" \
  --width "${M7_WIDTH:-10000}" \
  --threads "${M7_THREADS:-8}" \
  --out "${OUTDIR}/m7_fungal_pileup.json"

echo "[$(date -Is)] === M7 one-contig wide window (full-contig sample) ==="
python -u - <<'PY'
import json, time
from pathlib import Path
import pysam, samrust

bam = Path(__import__("os").environ["BAM"])
outdir = Path(__import__("os").environ["OUTDIR"])
FLAG = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
py = pysam.AlignmentFile(str(bam), "rb")
sr = samrust.AlignmentFile(str(bam), "rb")
# First contig, middle 50kb window as "full-genome style" stress without scanning all contigs
contig = py.references[0]
length = int(py.lengths[0])
start = max(0, length // 2 - 25000)
stop = min(length, start + 50000)

def pysam_counts(contig, start, stop):
    n = stop - start
    a, c, g, t, nn, d = ([0] * n for _ in range(6))
    for col in py.pileup(contig, start, stop, truncate=True, min_base_quality=0,
                         flag_filter=FLAG, stepper="all"):
        pos = col.reference_pos
        if pos < start or pos >= stop:
            continue
        idx = pos - start
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip:
                continue
            qpos = pr.query_position
            if qpos is None:
                continue
            base = (pr.alignment.query_sequence or "N")[qpos].upper()
            d[idx] += 1
            if base == "A":
                a[idx] += 1
            elif base == "C":
                c[idx] += 1
            elif base == "G":
                g[idx] += 1
            elif base == "T":
                t[idx] += 1
            else:
                nn[idx] += 1
    return a, c, g, t, nn, d

t0 = time.perf_counter()
pa,pc,pg,pt,pn,pd = pysam_counts(contig, start, stop)
serial = sr.pileup_counts(contig, start, stop, min_base_quality=0, threads=1)
parallel = sr.pileup_counts(contig, start, stop, min_base_quality=0, threads=8)
elapsed = time.perf_counter() - t0
mism = []
for key, exp in (("A",pa),("C",pc),("G",pg),("T",pt),("N",pn),("depth",pd)):
    if list(serial[key]) != exp:
        mism.append(key)
par_mism = [k for k in ("A","C","G","T","N","depth") if list(parallel[k]) != list(serial[k])]
payload = {
    "contig": contig, "start": start, "stop": stop, "elapsed_s": elapsed,
    "pysam_channel_mismatches": mism, "parallel_mismatches": par_mism,
    "gate_ok": not mism and not par_mism,
}
out = outdir / "m7_fungal_wide_window.json"
out.write_text(json.dumps(payload, indent=2) + "\n")
print(json.dumps(payload))
raise SystemExit(0 if payload["gate_ok"] else 1)
PY

echo "[$(date -Is)] === M8 ALT>=10 gate + scaling ==="
SR_TSV_ARGS=()
if [[ -f "${OUTDIR}/m8_recount_t1.tsv" ]]; then
  SR_TSV_ARGS=(--samrust-tsv "${OUTDIR}/m8_recount_t1.tsv")
  echo "Reusing existing ${OUTDIR}/m8_recount_t1.tsv"
fi

python -u scripts/run_m8_validation.py \
  --bam "$BAM" \
  --sites "$SITES" \
  --sample "$SAMPLE" \
  --threads "${M8_THREADS:-1:2:4:8:16}" \
  --scaling-outdir "$OUTDIR" \
  --out "${OUTDIR}/m8_fungal_validation.json" \
  "${SR_TSV_ARGS[@]}"

echo "[$(date -Is)] M7/M8 heavy validation finished OK"
