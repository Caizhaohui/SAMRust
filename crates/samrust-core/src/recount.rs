//! Candidate-site recount utility (M8 benchmark / validation — not a caller).

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rayon::prelude::*;

use crate::coords::Interval;
use crate::error::{Result, SamRustError};
use crate::indexed::IndexedAlignmentReader;
use crate::pileup::{pileup_counts, PileupCounts, PileupFilter};

/// One recount site (0-based reference position).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecountSite {
    pub chrom: String,
    /// 0-based position.
    pub pos0: u64,
    pub ref_allele: String,
    pub alt_allele: String,
}

/// One sample×site recount row.
#[derive(Debug, Clone, PartialEq)]
pub struct RecountRow {
    pub sample: String,
    pub chrom: String,
    /// 1-based position for output tables.
    pub pos_1based: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub a: u32,
    pub c: u32,
    pub g: u32,
    pub t: u32,
    pub n: u32,
    pub depth: u32,
    pub alt_count: u32,
    pub allele_frequency: f64,
}

/// Load sites from a BED (`chrom start stop REF>ALT`) or TSV (`chrom pos ref alt`, 1-based pos).
pub fn load_sites(path: &Path) -> Result<Vec<RecountSite>> {
    let file = File::open(path).map_err(SamRustError::from)?;
    let reader = BufReader::new(file);
    let mut sites = Vec::new();
    for (lineno, line) in reader.lines().enumerate() {
        let line = line.map_err(SamRustError::from)?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        // Skip TSV header
        if lineno == 0 && line.to_ascii_lowercase().starts_with("chrom") {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() >= 4 && cols[3].contains('>') {
            // BED: chrom start stop REF>ALT (0-based half-open)
            let start: u64 = cols[1]
                .parse()
                .map_err(|e| SamRustError::InvalidArgument(format!("bad BED start: {e}")))?;
            let allele = cols[3];
            let (reff, alt) = allele.split_once('>').ok_or_else(|| {
                SamRustError::InvalidArgument(format!("expected REF>ALT in column 4, got {allele}"))
            })?;
            sites.push(RecountSite {
                chrom: cols[0].to_string(),
                pos0: start,
                ref_allele: reff.to_string(),
                alt_allele: alt.to_string(),
            });
        } else if cols.len() >= 4 {
            // TSV: chrom pos(1-based) ref alt
            let pos1: u64 = cols[1]
                .parse()
                .map_err(|e| SamRustError::InvalidArgument(format!("bad TSV pos: {e}")))?;
            if pos1 == 0 {
                return Err(SamRustError::InvalidArgument(
                    "1-based TSV position must be >= 1".into(),
                ));
            }
            sites.push(RecountSite {
                chrom: cols[0].to_string(),
                pos0: pos1 - 1,
                ref_allele: cols[2].to_string(),
                alt_allele: cols[3].to_string(),
            });
        } else {
            return Err(SamRustError::InvalidArgument(format!(
                "unrecognized site line {}: {line}",
                lineno + 1
            )));
        }
    }
    Ok(sites)
}

fn base_count(counts: &PileupCounts, idx: usize, base: u8) -> u32 {
    match base.to_ascii_uppercase() {
        b'A' => counts.a[idx],
        b'C' => counts.c[idx],
        b'G' => counts.g[idx],
        b'T' => counts.t[idx],
        _ => counts.n[idx],
    }
}

fn row_from_counts(sample: &str, site: &RecountSite, counts: &PileupCounts) -> RecountRow {
    let idx = 0usize; // single-base window
    let a = counts.a[idx];
    let c = counts.c[idx];
    let g = counts.g[idx];
    let t = counts.t[idx];
    let n = counts.n[idx];
    let depth = counts.depth[idx];
    let alt_base = site.alt_allele.as_bytes().first().copied().unwrap_or(b'N');
    let alt_count = if site.alt_allele.len() == 1 {
        base_count(counts, idx, alt_base)
    } else {
        // Multi-base ALT: not fully modeled in base pileup; report 0 for M8 SNP gate.
        0
    };
    let allele_frequency = if depth > 0 {
        f64::from(alt_count) / f64::from(depth)
    } else {
        0.0
    };
    RecountRow {
        sample: sample.to_string(),
        chrom: site.chrom.clone(),
        pos_1based: site.pos0 + 1,
        ref_allele: site.ref_allele.clone(),
        alt_allele: site.alt_allele.clone(),
        a,
        c,
        g,
        t,
        n,
        depth,
        alt_count,
        allele_frequency,
    }
}

/// Recount all sites for one BAM.
///
/// `threads > 1` parallelizes across sites. Each rayon worker opens its
/// indexed reader once (`map_init`) and reuses it — previously every site
/// reopened the BAM and its index (the dominant 16T recount cost).
pub fn recount_bam(
    bam_path: &Path,
    sample: &str,
    sites: &[RecountSite],
    filter: PileupFilter,
    threads: usize,
) -> Result<Vec<RecountRow>> {
    if threads <= 1 {
        let mut reader = IndexedAlignmentReader::open(bam_path)?;
        let mut rows = Vec::with_capacity(sites.len());
        for site in sites {
            let interval = Interval::new(site.pos0, site.pos0 + 1)?;
            let counts = pileup_counts(&mut reader, &site.chrom, interval, filter)?;
            rows.push(row_from_counts(sample, site, &counts));
        }
        return Ok(rows);
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| SamRustError::InvalidArgument(format!("rayon pool: {e}")))?;

    pool.install(|| {
        sites
            .par_iter()
            .map_init(
                || IndexedAlignmentReader::open(bam_path),
                |slot, site| {
                    let reader = slot.as_mut().map_err(|e| {
                        SamRustError::InvalidArgument(format!("reopen {}: {e}", bam_path.display()))
                    })?;
                    let interval = Interval::new(site.pos0, site.pos0 + 1)?;
                    let counts = pileup_counts(reader, &site.chrom, interval, filter)?;
                    Ok(row_from_counts(sample, site, &counts))
                },
            )
            .collect()
    })
}

/// Format recount rows as TSV (with header).
pub fn rows_to_tsv(rows: &[RecountRow]) -> String {
    let mut out = String::from("sample\tchrom\tpos\tref\tA\tC\tG\tT\tN\tDP\tALT_COUNT\tAF\talt\n");
    for r in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\n",
            r.sample,
            r.chrom,
            r.pos_1based,
            r.ref_allele,
            r.a,
            r.c,
            r.g,
            r.t,
            r.n,
            r.depth,
            r.alt_count,
            r.allele_frequency,
            r.alt_allele
        ));
    }
    out
}

/// Sites with `ALT_COUNT >= threshold`.
pub fn filter_alt_ge(rows: &[RecountRow], threshold: u32) -> Vec<RecountRow> {
    rows.iter()
        .filter(|r| r.alt_count >= threshold)
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn load_bed_and_tsv_sites() {
        let dir = std::env::temp_dir().join("samrust_recount_sites");
        let _ = std::fs::create_dir_all(&dir);
        let bed = dir.join("sites.bed");
        {
            let mut f = File::create(&bed).unwrap();
            writeln!(f, "chr1\t10\t11\tA>G").unwrap();
            writeln!(f, "chr1\t20\t21\tC>T").unwrap();
        }
        let bed_sites = load_sites(&bed).unwrap();
        assert_eq!(bed_sites.len(), 2);
        assert_eq!(bed_sites[0].pos0, 10);
        assert_eq!(bed_sites[0].ref_allele, "A");
        assert_eq!(bed_sites[0].alt_allele, "G");

        let tsv = dir.join("sites.tsv");
        {
            let mut f = File::create(&tsv).unwrap();
            writeln!(f, "chrom\tpos\tref\talt\tvartype").unwrap();
            writeln!(f, "chr1\t11\tA\tG\tSNP").unwrap();
        }
        let tsv_sites = load_sites(&tsv).unwrap();
        assert_eq!(tsv_sites[0].pos0, 10);
        assert_eq!(tsv_sites[0].alt_allele, "G");
    }
}
