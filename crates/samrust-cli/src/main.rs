//! SAMRust CLI (M2 dump-records; M8 recount benchmark utility).

use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::{Parser, Subcommand};
use samrust_core::{
    load_sites, recount_bam, rows_to_tsv, AlignmentReader, PileupFilter, TagValue, VERSION,
};

#[derive(Debug, Parser)]
#[command(name = "samrust", version = VERSION, about = "SAMRust HTS utilities")]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Print version information.
    Version,
    /// Dump alignment records as JSONL for pysam parity checks (M2).
    DumpRecords {
        /// Input BAM path.
        #[arg(long)]
        bam: PathBuf,
        /// Maximum records to emit (0 = all).
        #[arg(long, default_value_t = 0)]
        limit: usize,
    },
    /// Recount A/C/G/T/N/DP/ALT at candidate sites (M8 benchmark — not a caller).
    Recount {
        /// Input BAM path (indexed).
        #[arg(long)]
        bam: PathBuf,
        /// Candidate sites: BED (`chrom start stop REF>ALT`) or TSV (`chrom pos ref alt`).
        #[arg(long)]
        sites: PathBuf,
        /// Sample name for output rows.
        #[arg(long)]
        sample: String,
        /// Worker threads (site-parallel when > 1).
        #[arg(long, default_value_t = 1)]
        threads: usize,
        /// Minimum base quality.
        #[arg(long, default_value_t = 0)]
        min_base_quality: u8,
        /// Minimum mapping quality.
        #[arg(long, default_value_t = 0)]
        min_mapping_quality: u8,
        /// Write TSV to this path (default: stdout).
        #[arg(long)]
        output: Option<PathBuf>,
        /// Only emit rows with ALT_COUNT >= N (0 = all rows).
        #[arg(long, default_value_t = 0)]
        min_alt: u32,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    match cli.command {
        Some(Commands::Version) | None => {
            println!("samrust {VERSION}");
            Ok(())
        }
        Some(Commands::DumpRecords { bam, limit }) => dump_records(&bam, limit),
        Some(Commands::Recount {
            bam,
            sites,
            sample,
            threads,
            min_base_quality,
            min_mapping_quality,
            output,
            min_alt,
        }) => recount(
            &bam,
            &sites,
            &sample,
            threads,
            min_base_quality,
            min_mapping_quality,
            output.as_ref(),
            min_alt,
        ),
    }
}

fn dump_records(bam: &Path, limit: usize) -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = AlignmentReader::open(bam)?;
    let header = reader.header();
    let mut out = io::stdout().lock();

    // First line: header summary
    writeln!(
        out,
        "{{\"type\":\"header\",\"nreferences\":{},\"references\":[{}],\"lengths\":[{}]}}",
        header.nreferences(),
        header
            .references()
            .iter()
            .map(|r| format!("\"{}\"", escape(r)))
            .collect::<Vec<_>>()
            .join(","),
        header
            .lengths()
            .iter()
            .map(|l| l.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )?;

    let mut n = 0usize;
    for result in reader.records() {
        let rec = result?;
        let cigar = match rec.cigarstring() {
            Some(s) => format!("\"{}\"", escape(&s)),
            None => "null".to_string(),
        };
        let tags = selected_tags_json(rec.tags());
        writeln!(
            out,
            "{{\"type\":\"record\",\"qname\":\"{}\",\"flag\":{},\"reference_id\":{},\"reference_start\":{},\"mapping_quality\":{},\"cigar\":{},\"query_length\":{},\"tags\":{{{}}}}}",
            escape(rec.query_name()),
            rec.flag(),
            rec.reference_id(),
            rec.reference_start(),
            rec.mapping_quality(),
            cigar,
            rec.query_length(),
            tags
        )?;
        n += 1;
        if limit != 0 && n >= limit {
            break;
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn recount(
    bam: &Path,
    sites_path: &Path,
    sample: &str,
    threads: usize,
    min_base_quality: u8,
    min_mapping_quality: u8,
    output: Option<&PathBuf>,
    min_alt: u32,
) -> Result<(), Box<dyn std::error::Error>> {
    let sites = load_sites(sites_path)?;
    let filter = PileupFilter {
        min_base_quality,
        min_mapping_quality,
        ..PileupFilter::default()
    };
    let t0 = Instant::now();
    let mut rows = recount_bam(bam, sample, &sites, filter, threads.max(1))?;
    if min_alt > 0 {
        rows.retain(|r| r.alt_count >= min_alt);
    }
    let elapsed = t0.elapsed();
    let tsv = rows_to_tsv(&rows);
    match output {
        Some(path) => {
            let mut f = File::create(path)?;
            f.write_all(tsv.as_bytes())?;
        }
        None => {
            let mut out = io::stdout().lock();
            out.write_all(tsv.as_bytes())?;
        }
    }
    eprintln!(
        "samrust recount: {} sites -> {} rows in {:.3}s (threads={})",
        sites.len(),
        rows.len(),
        elapsed.as_secs_f64(),
        threads.max(1)
    );
    Ok(())
}

fn selected_tags_json(tags: &samrust_core::Tags) -> String {
    const KEYS: &[&str] = &["NM", "RG", "MD", "AS", "XS"];
    let mut parts = Vec::new();
    for key in KEYS {
        if let Some(val) = tags.get(key) {
            parts.push(format!("\"{key}\":{}", tag_json(val)));
        }
    }
    parts.join(",")
}

fn tag_json(val: &TagValue) -> String {
    match val {
        TagValue::Int(n) => n.to_string(),
        TagValue::Float(f) => f.to_string(),
        TagValue::Char(c) => format!("\"{c}\""),
        TagValue::Str(s) | TagValue::Other(s) => format!("\"{}\"", escape(s)),
        TagValue::IntArray(_, v) => format!(
            "[{}]",
            v.iter()
                .map(|n| n.to_string())
                .collect::<Vec<_>>()
                .join(",")
        ),
        TagValue::FloatArray(v) => format!(
            "[{}]",
            v.iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(",")
        ),
    }
}

fn escape(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}
