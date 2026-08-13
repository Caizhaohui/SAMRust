# API compatibility

The living matrix is **[COMPATIBILITY.md](COMPATIBILITY.md)** (plan §25 name: `API_COMPATIBILITY.md`).

Coordinates: Python `fetch` / `count` / coverage / depth / pileup / VariantFile `fetch` are **0-based half-open**. `VariantRecord.pos` is **1-based**, like pysam.
