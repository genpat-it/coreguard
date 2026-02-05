// Suppress warnings for legacy code modules
#![allow(dead_code, unused_variables, unused_mut)]

//! alignment-qc: Pre-alignment QC framework for identifying problematic samples
//!
//! Detects samples that would introduce gaps or inconsistencies in multi-sample
//! SNP analysis, allowing users to review them before running SNP pipelines.
//!
//! Uses two alignment strategies:
//! 1. Pairwise alignment (minimap2) - for cluster detection and outlier identification
//! 2. Reference alignment (minimap2) - for gap contagion prediction

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use log::{info, warn};
use std::path::PathBuf;

mod cluster;
mod empty_column;
mod fasta;
mod gaps;
mod gff;
mod impact;
mod metrics;
mod output;
mod pairwise;
mod reads;
mod reference;
mod scoring;
mod bam_validate;
mod snp_compare;
mod vcf;
mod vcf_compare;
mod config;
mod compare;
mod parsers;
mod pileup;
#[cfg(feature = "serve")]
mod serve;

use crate::fasta::Sample;
use crate::metrics::SampleMetrics;
use crate::output::Report;
use crate::reference::ReferenceResults;
use crate::scoring::score_samples;

/// Pre-alignment QC for identifying problematic samples in SNP analysis
#[derive(Parser, Debug)]
#[command(name = "coreguard")]
#[command(author = "A. De Ruvo <andrea.deruvo@gssi.it>")]
#[command(version)]
#[command(about = "Guard your core genome - identify problematic samples before SNP pipeline analysis")]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,

    #[command(flatten)]
    args: Args,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Compare SNPs across multiple pipelines (VCF/BAM → compact JSON for visualization)
    Compare(CompareArgs),

    /// Convert pipeline output to CoreGuard TSV format
    Convert(ConvertArgs),
}

/// Arguments for the compare subcommand
#[derive(Parser, Debug)]
struct CompareArgs {
    /// YAML configuration file defining reference, samples, and pipelines
    #[arg(short, long)]
    config: PathBuf,

    /// Output JSON file for WASM visualization
    #[arg(short, long, default_value = "report.json")]
    output: PathBuf,

    /// Use compact JSON (no pretty-printing, smaller file size)
    #[arg(long)]
    compact: bool,

    /// Compress output with gzip (.json.gz)
    #[arg(long)]
    gzip: bool,

    /// Use binary format (bincode) instead of JSON (~10x faster to parse)
    #[arg(long)]
    binary: bool,

    /// Start web server to view the report in browser
    #[cfg(feature = "serve")]
    #[arg(long)]
    serve: bool,

    /// Port for web server (default: 8765)
    #[cfg(feature = "serve")]
    #[arg(long, default_value = "8765")]
    port: u16,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

/// Arguments for the convert subcommand
#[derive(Parser, Debug)]
struct ConvertArgs {
    /// Source format: cfsan
    #[arg(value_enum)]
    format: ConvertFormat,

    /// Output file (CoreGuard TSV format)
    #[arg(short, long)]
    output: PathBuf,

    /// CFSAN snplist.txt file (genomic positions)
    #[arg(long)]
    snplist: Option<PathBuf>,

    /// CFSAN snpma.fasta file (sample alleles)
    #[arg(long)]
    snpma: Option<PathBuf>,

    /// CFSAN referenceSNP.fasta file (reference alleles)
    #[arg(long)]
    reference_snp: Option<PathBuf>,

    /// Chromosome/contig name (default: auto-detect from snplist.txt)
    #[arg(long)]
    chrom: Option<String>,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

#[derive(clap::ValueEnum, Clone, Debug)]
enum ConvertFormat {
    /// CFSAN SNP Pipeline output (snplist.txt + snpma.fasta + referenceSNP.fasta)
    Cfsan,
}

/// Legacy args for the default QC analysis command
#[derive(Parser, Debug)]
struct Args {
    /// Input directory containing FASTA files
    #[arg(short, long)]
    input: Option<PathBuf>,

    /// Individual FASTA files to analyze
    #[arg(short, long, num_args = 1..)]
    samples: Option<Vec<PathBuf>>,

    /// Input directory containing FASTQ files (alternative to FASTA input)
    #[arg(long)]
    fastq_input: Option<PathBuf>,

    /// Individual FASTQ files to analyze (supports .fq, .fastq, .fq.gz, .fastq.gz)
    /// For paired-end, provide R1 files only - R2 will be auto-detected
    #[arg(long, num_args = 1..)]
    fastq_samples: Option<Vec<PathBuf>>,

    /// Minimum read depth to consider a position "covered" (for FASTQ input)
    /// Positions with depth < this value are considered gaps
    #[arg(long, default_value_t = 5)]
    min_depth: u32,

    /// Read type for FASTQ input: 'short' (Illumina) or 'long' (ONT/PacBio)
    #[arg(long, default_value = "short")]
    read_type: String,

    /// Reference genome (REQUIRED for FASTQ input, recommended for FASTA)
    #[arg(short, long)]
    reference: Option<PathBuf>,

    /// GFF/GTF annotation file for gene zone analysis
    /// Identifies which genes/CDS are affected by alignment gaps
    #[arg(long)]
    gff: Option<PathBuf>,

    /// Snippy results directory for SNP pipeline comparison
    /// Contains sample subdirectories (e.g., TE15676_snippy/) with aligned.fa files
    #[arg(long)]
    snippy_dir: Option<PathBuf>,

    /// Snippy-core results directory for comparison
    /// Should contain snippycore.txt and snippycore.tab files
    #[arg(long)]
    snippycore_dir: Option<PathBuf>,

    /// CFSAN SNP Pipeline results directory for comparison
    /// Should contain metrics.tsv and snplist.txt files
    #[arg(long)]
    cfsan_dir: Option<PathBuf>,

    /// Directory containing Snippy VCF files for variant comparison
    /// Can be the parent directory of *_snippy folders
    #[arg(long)]
    vcf_snippy: Option<PathBuf>,

    /// Directory containing CFSAN VCF files for variant comparison
    /// Should contain samples/ subdirectory with var.flt.vcf files
    #[arg(long)]
    vcf_cfsan: Option<PathBuf>,

    /// Enable VCF comparison between pipelines
    /// Requires --vcf-snippy and --vcf-cfsan
    #[arg(long)]
    compare_variants: bool,

    /// Output file (JSON or TSV based on extension)
    #[arg(short, long, default_value = "report.json")]
    output: PathBuf,

    /// Number of threads for parallel processing
    #[arg(short, long, default_value_t = num_cpus::get())]
    threads: usize,

    /// Minimum identity threshold (flag samples below this)
    #[arg(long, default_value_t = 0.90)]
    min_identity: f64,

    /// Path to minimap2 binary (default: minimap2 in PATH)
    #[arg(long)]
    minimap2_path: Option<PathBuf>,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,

    /// Output BED files for IGV visualization
    #[arg(long)]
    bed_output: Option<PathBuf>,

    /// Fast mode: skip pairwise analysis, use reference-only for gap prediction
    /// Much faster for large datasets (N alignments instead of N²)
    #[arg(long)]
    fast: bool,

    /// Force pairwise analysis even with --fast (for identity matrix)
    #[arg(long)]
    with_pairwise: bool,

    /// Disable deduplication (by default, identical sequences are analyzed once)
    #[arg(long)]
    no_dedup: bool,

    /// Start a web server with interactive IGV.js genome browser
    /// Opens browser automatically with the report and genome visualization
    #[cfg(feature = "serve")]
    #[arg(long)]
    serve: bool,

    /// Port for the web server (default: 8765)
    #[cfg(feature = "serve")]
    #[arg(long, default_value_t = 8765)]
    port: u16,

    /// Serve an existing report without re-running analysis
    /// Provide path to a previously generated JSON report
    #[cfg(feature = "serve")]
    #[arg(long, value_name = "REPORT.json")]
    serve_only: Option<PathBuf>,

    /// Directory containing BAM files for IGV visualization
    /// BAM files should be indexed (.bai) and named like {sample}.bam or reads.sorted.bam
    #[cfg(feature = "serve")]
    #[arg(long)]
    bam_dir: Option<PathBuf>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Handle subcommands first
    if let Some(cmd) = cli.command {
        return match cmd {
            Commands::Compare(compare_args) => run_compare(compare_args),
            Commands::Convert(convert_args) => run_convert(convert_args),
        };
    }

    // Default: run the legacy QC analysis
    run_qc_analysis(cli.args)
}

/// Run the compare subcommand: YAML config → JSON report
fn run_compare(args: CompareArgs) -> Result<()> {
    // Initialize logging
    let log_level = if args.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level)).init();

    info!("coreguard compare v{}", env!("CARGO_PKG_VERSION"));
    info!("Loading configuration: {}", args.config.display());

    // Load and validate configuration
    let config = config::Config::from_yaml(&args.config)?;

    info!(
        "Configuration loaded: {} samples, {} pipelines",
        config.all_sample_ids().len(),
        config.all_pipeline_ids().len()
    );

    // Generate report
    info!("Generating comparison report...");
    let config_yaml = std::fs::read_to_string(&args.config).ok().map(|yaml| {
        // Sanitize absolute paths: make them relative to the config file directory
        let config_dir = args.config.parent()
            .and_then(|p| p.canonicalize().ok())
            .unwrap_or_default();
        let home_dir = std::env::var("HOME").unwrap_or_default();
        let mut sanitized = yaml;
        // Replace config dir prefix with relative paths
        if !config_dir.as_os_str().is_empty() {
            let prefix = format!("{}/", config_dir.display());
            sanitized = sanitized.replace(&prefix, "");
        }
        // Replace home directory with ~
        if !home_dir.is_empty() {
            sanitized = sanitized.replace(&home_dir, "~");
        }
        sanitized
    });
    let report = compare::CompareReport::from_config_with_yaml(&config, config_yaml)?;

    // Auto-detect gzip from output extension
    let use_gzip = args.gzip || args.output.extension().map(|e| e == "gz").unwrap_or(false);

    // Save report
    if args.binary {
        // Binary format (bincode)
        let output_path = if use_gzip {
            if args.output.extension().map(|e| e == "gz").unwrap_or(false) {
                args.output.clone()
            } else {
                args.output.with_extension("bin.gz")
            }
        } else {
            if args.output.extension().map(|e| e == "bin").unwrap_or(false) {
                args.output.clone()
            } else {
                args.output.with_extension("bin")
            }
        };

        if use_gzip {
            report.save_binary_gzip(&output_path)?;
            info!("Gzipped binary report saved to: {}", output_path.display());
        } else {
            report.save_binary(&output_path)?;
            info!("Binary report saved to: {}", output_path.display());
        }
    } else if use_gzip {
        // Add .gz extension if not present
        let output_path = if args.output.extension().map(|e| e == "gz").unwrap_or(false) {
            args.output.clone()
        } else {
            args.output.with_extension("json.gz")
        };
        report.save_gzip(&output_path, args.compact)?;
        info!("Gzipped report saved to: {}", output_path.display());
    } else if args.compact {
        report.save_compact(&args.output)?;
        info!("Compact report saved to: {}", args.output.display());
    } else {
        report.save(&args.output)?;
        info!("Report saved to: {}", args.output.display());
    }

    info!(
        "Done! {} samples × {} pipelines processed",
        report.summary.total_samples,
        report.summary.total_pipelines
    );

    // Start web server if requested
    #[cfg(feature = "serve")]
    if args.serve {
        let output_path = if args.binary {
            if use_gzip {
                if args.output.extension().map(|e| e == "gz").unwrap_or(false) {
                    args.output.clone()
                } else {
                    args.output.with_extension("bin.gz")
                }
            } else if args.output.extension().map(|e| e == "bin").unwrap_or(false) {
                args.output.clone()
            } else {
                args.output.with_extension("bin")
            }
        } else if use_gzip {
            if args.output.extension().map(|e| e == "gz").unwrap_or(false) {
                args.output.clone()
            } else {
                args.output.with_extension("json.gz")
            }
        } else {
            args.output.clone()
        };
        serve_compare_viewer(&output_path, args.port)?;
    }

    Ok(())
}

/// Serve the compare viewer with embedded assets
#[cfg(feature = "serve")]
fn serve_compare_viewer(report_path: &PathBuf, port: u16) -> Result<()> {
    use std::sync::Arc;
    use tiny_http::{Server, Response, Header};

    // Embed viewer.html at compile time
    const VIEWER_HTML: &str = include_str!("../viewer.html");

    // Embed WASM files at compile time
    const WASM_JS: &str = include_str!("../wasm/pkg/coreguard_wasm.js");
    const WASM_BIN: &[u8] = include_bytes!("../wasm/pkg/coreguard_wasm_bg.wasm");

    // Read the report file
    let report_content = std::fs::read(report_path)
        .with_context(|| format!("Failed to read report: {}", report_path.display()))?;
    let is_gzip = report_path.extension().map(|e| e == "gz").unwrap_or(false);

    let report_content = Arc::new(report_content);

    // Modify viewer.html to point to the correct WASM path
    let viewer_html = VIEWER_HTML.replace(
        "from './wasm/pkg/coreguard_wasm.js'",
        "from '/wasm/pkg/coreguard_wasm.js'"
    );

    let addr = format!("0.0.0.0:{}", port);
    let server = Server::http(&addr)
        .map_err(|e| anyhow::anyhow!("Failed to start server: {}", e))?;

    let url = format!("http://localhost:{}", port);
    info!("Starting viewer at {}", url);
    info!("Press Ctrl+C to stop the server");

    // Try to open browser
    if webbrowser::open(&url).is_err() {
        warn!("Could not open browser automatically. Please visit: {}", url);
    }

    for request in server.incoming_requests() {
        let path = request.url().to_string();

        let response = match path.as_str() {
            "/" | "/index.html" | "/viewer.html" => {
                Response::from_string(&viewer_html)
                    .with_header(Header::from_bytes("Content-Type", "text/html; charset=utf-8").unwrap())
            }
            "/report.json" => {
                if is_gzip {
                    // Decompress for non-gzip request
                    use flate2::read::GzDecoder;
                    use std::io::Read;
                    let mut decoder = GzDecoder::new(&report_content[..]);
                    let mut decompressed = String::new();
                    decoder.read_to_string(&mut decompressed)?;
                    Response::from_string(decompressed)
                        .with_header(Header::from_bytes("Content-Type", "application/json").unwrap())
                } else {
                    Response::from_data(report_content.as_ref().clone())
                        .with_header(Header::from_bytes("Content-Type", "application/json").unwrap())
                }
            }
            "/report.json.gz" => {
                if is_gzip {
                    Response::from_data(report_content.as_ref().clone())
                        .with_header(Header::from_bytes("Content-Type", "application/gzip").unwrap())
                        .with_header(Header::from_bytes("Content-Encoding", "gzip").unwrap())
                } else {
                    // Compress on the fly
                    use flate2::write::GzEncoder;
                    use flate2::Compression;
                    use std::io::Write;
                    let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
                    encoder.write_all(&report_content)?;
                    let compressed = encoder.finish()?;
                    Response::from_data(compressed)
                        .with_header(Header::from_bytes("Content-Type", "application/gzip").unwrap())
                }
            }
            "/wasm/pkg/coreguard_wasm.js" => {
                Response::from_string(WASM_JS)
                    .with_header(Header::from_bytes("Content-Type", "application/javascript").unwrap())
            }
            "/wasm/pkg/coreguard_wasm_bg.wasm" => {
                Response::from_data(WASM_BIN.to_vec())
                    .with_header(Header::from_bytes("Content-Type", "application/wasm").unwrap())
            }
            // Pako.js CDN fallback - redirect to CDN
            path if path.contains("pako") => {
                Response::from_string("")
                    .with_status_code(302)
                    .with_header(Header::from_bytes("Location", "https://cdnjs.cloudflare.com/ajax/libs/pako/2.1.0/pako.min.js").unwrap())
            }
            _ => {
                Response::from_string("Not found").with_status_code(404)
            }
        };

        let _ = request.respond(response);
    }

    Ok(())
}

/// Run the convert subcommand: convert pipeline output to CoreGuard TSV format
fn run_convert(args: ConvertArgs) -> Result<()> {
    // Initialize logging
    let log_level = if args.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level)).init();

    info!("coreguard convert v{}", env!("CARGO_PKG_VERSION"));

    match args.format {
        ConvertFormat::Cfsan => convert_cfsan(&args)?,
    }

    Ok(())
}

/// Convert CFSAN output to CoreGuard TSV format
fn convert_cfsan(args: &ConvertArgs) -> Result<()> {
    use std::io::{BufRead, BufReader, Write};
    use std::collections::HashMap;

    let snplist = args.snplist.as_ref()
        .context("--snplist is required for CFSAN conversion")?;
    let snpma = args.snpma.as_ref()
        .context("--snpma is required for CFSAN conversion")?;
    let reference_snp = args.reference_snp.as_ref()
        .context("--reference-snp is required for CFSAN conversion")?;

    info!("Reading CFSAN files:");
    info!("  snplist.txt: {}", snplist.display());
    info!("  snpma.fasta: {}", snpma.display());
    info!("  referenceSNP.fasta: {}", reference_snp.display());

    // 1. Read genomic positions and chromosome from snplist.txt
    let snplist_file = std::fs::File::open(snplist)
        .with_context(|| format!("Failed to open {}", snplist.display()))?;
    let reader = BufReader::new(snplist_file);

    let mut positions: Vec<(String, usize)> = Vec::new(); // (chrom, pos)
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 { continue; }
        let chrom = parts[0].to_string();
        let pos: usize = parts[1].parse().unwrap_or(0);
        if pos == 0 { continue; }
        positions.push((chrom, pos)); // Keep 1-based for output
    }
    info!("  Loaded {} positions from snplist.txt", positions.len());

    // 2. Read reference alleles from referenceSNP.fasta
    let ref_file = std::fs::File::open(reference_snp)
        .with_context(|| format!("Failed to open {}", reference_snp.display()))?;
    let reader = BufReader::new(ref_file);
    let mut ref_seq = String::new();
    for line in reader.lines() {
        let line = line?;
        if !line.starts_with('>') {
            ref_seq.push_str(line.trim());
        }
    }
    info!("  Loaded {} reference alleles", ref_seq.len());

    if ref_seq.len() != positions.len() {
        anyhow::bail!(
            "Length mismatch: referenceSNP.fasta has {} alleles, snplist.txt has {} positions",
            ref_seq.len(), positions.len()
        );
    }

    // 3. Read sample alleles from snpma.fasta
    let snpma_file = std::fs::File::open(snpma)
        .with_context(|| format!("Failed to open {}", snpma.display()))?;
    let reader = BufReader::new(snpma_file);
    let mut sample_seqs: HashMap<String, String> = HashMap::new();
    let mut sample_order: Vec<String> = Vec::new();
    let mut current: Option<String> = None;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            let name = line[1..].trim().to_string();
            if !sample_seqs.contains_key(&name) {
                sample_order.push(name.clone());
            }
            current = Some(name.clone());
            sample_seqs.entry(name).or_default();
        } else if let Some(ref name) = current {
            sample_seqs.get_mut(name).unwrap().push_str(line.trim());
        }
    }
    info!("  Loaded {} samples from snpma.fasta", sample_order.len());

    // Verify all sequences have same length
    for (name, seq) in &sample_seqs {
        if seq.len() != positions.len() {
            anyhow::bail!(
                "Sample {} has {} alleles but expected {} (from snplist.txt)",
                name, seq.len(), positions.len()
            );
        }
    }

    // 4. Write CoreGuard TSV format
    let mut output = std::fs::File::create(&args.output)
        .with_context(|| format!("Failed to create {}", args.output.display()))?;

    // Header: CHR, POS, REF, sample1, sample2, ...
    write!(output, "CHR\tPOS\tREF")?;
    for sample in &sample_order {
        write!(output, "\t{}", sample)?;
    }
    writeln!(output)?;

    // Data rows
    let ref_bytes = ref_seq.as_bytes();
    for (i, (chrom, pos)) in positions.iter().enumerate() {
        let ref_allele = ref_bytes[i] as char;
        write!(output, "{}\t{}\t{}", chrom, pos, ref_allele)?;

        for sample in &sample_order {
            let seq = sample_seqs.get(sample).unwrap();
            let allele = seq.as_bytes()[i] as char;
            write!(output, "\t{}", allele)?;
        }
        writeln!(output)?;
    }

    info!("Wrote {} positions × {} samples to {}",
          positions.len(), sample_order.len(), args.output.display());
    info!("CoreGuard TSV format ready for use with 'coreguard compare'");

    Ok(())
}

/// Run the legacy QC analysis (default command)
fn run_qc_analysis(args: Args) -> Result<()> {
    let start_time = std::time::Instant::now();
    let command_line = std::env::args().collect::<Vec<_>>().join(" ");
    let timestamp = chrono::Local::now().format("%Y-%m-%d %H:%M:%S").to_string();

    // Initialize logging
    let log_level = if args.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level)).init();

    // Set thread pool size
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

    info!("coreguard v{}", env!("CARGO_PKG_VERSION"));

    // === Serve-only mode: load existing report and serve ===
    #[cfg(feature = "serve")]
    if let Some(ref report_path) = args.serve_only {
        return serve_existing_report(
            report_path,
            args.reference.as_deref(),
            args.bed_output.as_deref(),
            args.port,
            args.bam_dir.as_deref(),
        );
    }

    info!("Using {} threads", args.threads);

    // Detect input mode: FASTQ or FASTA
    let fastq_mode = args.fastq_input.is_some() || args.fastq_samples.is_some();

    // FASTQ mode requires a reference
    if fastq_mode && args.reference.is_none() {
        anyhow::bail!("--reference is REQUIRED for FASTQ input mode (reads must be mapped to reference)");
    }

    // Load samples (FASTA mode) or virtual samples (FASTQ mode)
    let (all_samples, _fastq_results): (Vec<Sample>, Option<Vec<reads::ReadMappingResult>>) = if fastq_mode {
        // FASTQ mode: load reference first, then map reads
        let reference = fasta::load_sample(args.reference.as_ref().unwrap())?;
        info!("FASTQ mode: mapping reads to reference ({} bp)", reference.total_length);

        let (samples, results) = load_fastq_samples_and_analyze(&args, &reference)?;
        (samples, Some(results))
    } else {
        // FASTA mode: load assemblies
        let samples = load_samples(&args)?;
        (samples, None)
    };

    info!("Loaded {} samples", all_samples.len());

    if all_samples.len() < 2 {
        anyhow::bail!("Need at least 2 samples for comparison");
    }

    // Deduplication: group samples by sequence hash
    let (samples, _dedup_groups) = if args.no_dedup {
        info!("Deduplication disabled");
        let groups: std::collections::HashMap<String, Vec<String>> = all_samples
            .iter()
            .map(|s| (s.name.clone(), vec![s.name.clone()]))
            .collect();
        (all_samples, groups)
    } else {
        use std::collections::HashMap;
        let mut hash_to_samples: HashMap<String, Vec<String>> = HashMap::new();
        for sample in &all_samples {
            hash_to_samples
                .entry(sample.sequence_hash.clone())
                .or_default()
                .push(sample.name.clone());
        }

        let unique_count = hash_to_samples.len();
        let duplicate_count = all_samples.len() - unique_count;

        if duplicate_count > 0 {
            info!("DEDUP: Found {} unique sequences ({} duplicates removed)",
                  unique_count, duplicate_count);

            // Log groups with duplicates
            for (hash, names) in &hash_to_samples {
                if names.len() > 1 {
                    log::debug!("  Hash {}...: {} samples ({}, ...)",
                               &hash[..8], names.len(), names[0]);
                }
            }
        } else {
            info!("DEDUP: All {} samples are unique", all_samples.len());
        }

        // Keep only one representative per hash
        let unique_samples: Vec<Sample> = all_samples
            .into_iter()
            .filter(|s| {
                hash_to_samples.get(&s.sequence_hash)
                    .map(|names| names[0] == s.name)
                    .unwrap_or(false)
            })
            .collect();

        // Build groups: representative -> all samples with same hash
        let dedup_groups: HashMap<String, Vec<String>> = hash_to_samples
            .into_iter()
            .map(|(_, names)| (names[0].clone(), names))
            .collect();

        (unique_samples, dedup_groups)
    };

    info!("Analyzing {} unique samples", samples.len());

    // Load reference if provided
    let reference: Option<Sample> = if let Some(ref ref_path) = args.reference {
        info!("Loading reference: {}", ref_path.display());
        Some(fasta::load_sample(ref_path)?)
    } else {
        None
    };

    // Reference is now REQUIRED for gap analysis
    if reference.is_none() {
        anyhow::bail!("--reference is REQUIRED for gap analysis. Please provide a reference genome.");
    }

    if fastq_mode {
        info!("FASTQ MODE: Using read coverage for gap detection (min depth: {})", args.min_depth);
    }

    // === Analysis 1: Reference alignment (minimap2) ===
    // All samples are aligned to reference to get gap regions
    // This is the PRIMARY analysis - all other metrics derive from this
    let reference_results: ReferenceResults = if fastq_mode {
        // FASTQ mode: build reference results from coverage analysis
        info!("Building reference results from FASTQ coverage analysis...");
        let fastq_results = _fastq_results.as_ref().unwrap();
        reference::build_from_fastq_coverage(
            reference.as_ref().unwrap(),
            &samples,
            fastq_results,
        )
    } else {
        info!("Running reference alignment analysis...");
        reference::analyze_vs_reference(
            reference.as_ref().unwrap(),
            &samples,
            args.minimap2_path.as_deref(),
            args.threads,
        )?
    };

    // === Analysis 2: Pairwise gap analysis ===
    // For each pair of samples, calculate:
    // - Gap union (positions that are gaps in A OR B)
    // - Gap intersection (positions that are gaps in A AND B)
    // - Quality score: intersection / union
    info!("Computing pairwise gap analysis...");
    let pairwise_results = pairwise::analyze_pairwise_gaps(&reference_results);

    // === Cluster analysis ===
    info!("Analyzing sample clusters and outliers...");
    let sample_names: Vec<String> = samples.iter().map(|s| s.name.clone()).collect();
    let cluster_analysis = cluster::analyze_clusters(
        &sample_names,
        &pairwise_results,
        &reference_results,
    );

    // Log cluster warnings
    for warning in &cluster_analysis.warnings {
        log::warn!("{}", warning);
    }
    for rec in &cluster_analysis.recommendations {
        log::info!("Recommendation: {}", rec);
    }

    // Calculate metrics for each sample
    info!("Calculating sample metrics...");
    let metrics: Vec<SampleMetrics> = samples
        .iter()
        .map(|s| {
            metrics::calculate_metrics(
                s,
                &pairwise_results,
                &reference_results,
            )
        })
        .collect();

    // Analyze impact on core alignment (gap contagion)
    info!("Analyzing core alignment impact (gap contagion)...");
    let ref_len = reference.as_ref().unwrap().total_length;
    let impact_analysis = impact::analyze_impact(
        &sample_names,
        &pairwise_results,
        &reference_results,
        ref_len,
    );

    // Log impact warnings
    for sample_impact in &impact_analysis.samples {
        if sample_impact.risk_level == impact::RiskLevel::Critical {
            log::warn!(
                "CRITICAL: {} has {} unique gap bases ({:.1}% core reduction)",
                sample_impact.sample,
                sample_impact.unique_gap_bases,
                sample_impact.estimated_core_reduction * 100.0
            );
        } else if sample_impact.risk_level == impact::RiskLevel::High {
            log::warn!(
                "HIGH RISK: {} has {} unique gap bases ({:.1}% core reduction)",
                sample_impact.sample,
                sample_impact.unique_gap_bases,
                sample_impact.estimated_core_reduction * 100.0
            );
        }
    }

    // Score samples and generate recommendations
    info!("Scoring samples...");
    let scored_samples = score_samples(&samples, &metrics, args.min_identity);

    // Print core genome estimate (before moving impact_analysis into report)
    info!(
        "Estimated core with all samples: {} bp ({:.1}% of reference)",
        impact_analysis.estimated_core_all,
        (1.0 - impact_analysis.total_core_reduction) * 100.0
    );

    // Generate report
    let mut report = Report::new(
        scored_samples,
        &pairwise_results,
        &reference_results,
        impact_analysis,
        cluster_analysis,
    );

    // Collect sample gaps for analysis (used by gene zone and empty column)
    let sample_gaps: Vec<(String, Vec<crate::gaps::Region>)> = reference_results
        .alignments
        .iter()
        .map(|a| (a.sample_name.clone(), a.reference_uncovered_regions.clone()))
        .collect();

    // Gene zone analysis if GFF provided
    if let Some(ref gff_path) = args.gff {
        info!("Loading GFF annotation: {}", gff_path.display());
        match gff::Annotation::from_file(gff_path) {
            Ok(annotation) => {
                info!(
                    "Loaded {} genes, {} CDS from GFF",
                    annotation.total_genes, annotation.total_cds
                );

                // Run gene zone analysis
                let gene_analyses = gff::analyze_gene_zones(&annotation, &sample_gaps);

                // Calculate summary stats
                let mut all_affected_genes = std::collections::HashSet::new();
                let mut all_removed_genes = std::collections::HashSet::new();
                for analysis in &gene_analyses {
                    for af in &analysis.affected_features {
                        if af.feature.feature_type == "gene" {
                            all_affected_genes.insert(af.feature.display_name());
                            if af.is_completely_removed() {
                                all_removed_genes.insert(af.feature.display_name());
                            }
                        }
                    }
                }

                let gene_summary = output::GeneZoneSummary {
                    gff_path: gff_path.display().to_string(),
                    total_genes: annotation.total_genes,
                    total_cds: annotation.total_cds,
                    samples: gene_analyses,
                    total_affected_genes: all_affected_genes.len(),
                    genes_removed_any_sample: all_removed_genes.len(),
                };

                info!(
                    "Gene zone analysis: {} genes affected, {} removed in at least one sample",
                    gene_summary.total_affected_genes, gene_summary.genes_removed_any_sample
                );

                report = report.with_gene_zone_analysis(gene_summary);
            }
            Err(e) => {
                log::warn!("Failed to parse GFF file: {}", e);
            }
        }
    }

    // Empty column analysis (always run)
    info!("Analyzing coverage thresholds...");
    let empty_column_analysis = empty_column::analyze_empty_columns(
        reference_results.reference_length,
        &sample_gaps,
        &empty_column::default_thresholds(),
    );
    info!(
        "Coverage analysis: {} fragile regions, core at 100% threshold: {:.1}%",
        empty_column_analysis.fragile_regions.len(),
        empty_column_analysis
            .threshold_results
            .first()
            .map(|t| t.core_pct)
            .unwrap_or(0.0)
    );
    report = report.with_empty_column_analysis(empty_column_analysis);

    // SNP pipeline comparison if Snippy directory provided
    if let Some(ref snippy_dir) = args.snippy_dir {
        info!("Loading Snippy results from: {}", snippy_dir.display());
        let sample_names: Vec<String> = reference_results
            .alignments
            .iter()
            .map(|a| a.sample_name.clone())
            .collect();

        match snp_compare::load_snippy_results(snippy_dir, &sample_names, Some(&reference_results.reference_name)) {
            Ok(snippy_filtered) => {
                if snippy_filtered.is_empty() {
                    log::warn!("No Snippy results found for any samples");
                } else {
                    let comparison = snp_compare::compare_with_snp_pipeline(
                        &sample_gaps,
                        &snippy_filtered,
                        "Snippy",
                        reference_results.reference_length,
                    );
                    info!(
                        "SNP comparison: avg concordance {:.1}%, {} concordant, {} coreguard-only, {} snippy-only",
                        comparison.summary.avg_concordance_pct,
                        comparison.summary.total_concordant_bases,
                        comparison.summary.total_coreguard_only_bases,
                        comparison.summary.total_snp_only_bases
                    );
                    report = report.with_snp_comparison(comparison);
                }
            }
            Err(e) => {
                log::warn!("Failed to load Snippy results: {}", e);
            }
        }
    }

    // Snippy-core comparison if directory provided
    if let Some(ref snippycore_dir) = args.snippycore_dir {
        info!("Loading snippy-core results from: {}", snippycore_dir.display());

        let txt_path = snippycore_dir.join("snippycore.txt");
        let tab_path = snippycore_dir.join("snippycore.tab");

        match (snp_compare::parse_snippycore_txt(&txt_path), snp_compare::parse_snippycore_tab(&tab_path)) {
            (Ok(stats), Ok(snps)) => {
                info!("Loaded snippy-core stats for {} samples, {} SNP positions", stats.len(), snps.len());

                let mut comparison = snp_compare::compare_with_snippycore(
                    &sample_gaps,
                    &stats,
                    &snps,
                    reference_results.reference_length,
                );

                // Try to load distance matrix
                let distance_patterns = [
                    "hamming_distances_vanilla.tsv",
                    "hamming_distances_filtered.tsv",
                    "core.dist.tsv",
                ];
                if let Some(matrix) = snp_compare::find_distance_matrix(snippycore_dir, &distance_patterns) {
                    info!("Loaded Snippy distance matrix: {} samples", matrix.samples.len());
                    comparison.distance_matrix = Some(matrix);
                }

                info!(
                    "Snippy-core comparison: avg core diff {:.2}%, avg SNP retention {:.1}% ({} SNPs, {} in gaps)",
                    comparison.summary.avg_core_diff_pct,
                    comparison.summary.avg_snp_retention_pct,
                    comparison.summary.total_snps,
                    comparison.summary.total_snps_in_gaps
                );

                report = report.with_snippycore_comparison(comparison);
            }
            (Err(e), _) => {
                log::warn!("Failed to load snippycore.txt: {}", e);
            }
            (_, Err(e)) => {
                log::warn!("Failed to load snippycore.tab: {}", e);
            }
        }
    }

    // CFSAN comparison if directory provided
    if let Some(ref cfsan_dir) = args.cfsan_dir {
        info!("Loading CFSAN results from: {}", cfsan_dir.display());

        let metrics_path = cfsan_dir.join("metrics.tsv");
        let snplist_path = cfsan_dir.join("snplist.txt");

        match (snp_compare::parse_cfsan_metrics(&metrics_path), snp_compare::parse_cfsan_snplist_full(&snplist_path)) {
            (Ok(metrics), Ok((total_snps, snp_map))) => {
                info!("Loaded CFSAN metrics for {} samples, {} total SNP positions", metrics.len(), total_snps);

                let mut comparison = snp_compare::compare_with_cfsan(
                    &sample_gaps,
                    &metrics,
                    &snp_map,
                    total_snps,
                    reference_results.reference_length,
                );

                // Try to load distance matrix
                let distance_patterns = [
                    "snp_distance_matrix.tsv",
                    "snp_distance_matrix_preserved.tsv",
                ];
                if let Some(matrix) = snp_compare::find_distance_matrix(cfsan_dir, &distance_patterns) {
                    info!("Loaded CFSAN distance matrix: {} samples", matrix.samples.len());
                    comparison.distance_matrix = Some(matrix);
                }

                info!(
                    "CFSAN comparison: avg SNP retention {:.1}% ({} SNPs, {} in gaps)",
                    comparison.summary.avg_snp_retention_pct,
                    comparison.summary.total_snps,
                    comparison.summary.total_snps_in_gaps
                );

                report = report.with_cfsan_comparison(comparison);
            }
            (Err(e), _) => {
                log::warn!("Failed to load CFSAN metrics.tsv: {}", e);
            }
            (_, Err(e)) => {
                log::warn!("Failed to load CFSAN snplist.txt: {}", e);
            }
        }
    }

    // VCF comparison between pipelines
    // Store separately for server (full data with all positions)
    let mut vcf_comparison_full: Option<vcf_compare::VcfComparison> = None;

    if args.compare_variants {
        if let (Some(ref snippy_vcf_dir), Some(ref cfsan_vcf_dir)) = (&args.vcf_snippy, &args.vcf_cfsan) {
            info!("Loading VCF files for comparison...");

            // Load Snippy VCFs
            let snippy_vcfs = match vcf::load_vcfs_from_dir(snippy_vcf_dir, "snippy") {
                Ok(vcfs) => {
                    info!("Loaded {} Snippy VCF files", vcfs.len());
                    vcfs
                }
                Err(e) => {
                    log::warn!("Failed to load Snippy VCFs: {}", e);
                    Vec::new()
                }
            };

            // Load CFSAN VCFs
            let cfsan_vcfs = match vcf::load_vcfs_from_dir(cfsan_vcf_dir, "cfsan") {
                Ok(vcfs) => {
                    info!("Loaded {} CFSAN VCF files", vcfs.len());
                    vcfs
                }
                Err(e) => {
                    log::warn!("Failed to load CFSAN VCFs: {}", e);
                    Vec::new()
                }
            };

            if !snippy_vcfs.is_empty() && !cfsan_vcfs.is_empty() {
                // Build gap regions map for correlation
                let gap_regions: std::collections::HashMap<String, Vec<gaps::Region>> = sample_gaps.iter()
                    .map(|(name, gaps)| (name.clone(), gaps.clone()))
                    .collect();

                // Perform comparison
                let mut vcf_comparison = vcf_compare::compare_pipelines(
                    &snippy_vcfs,
                    &cfsan_vcfs,
                    "Snippy",
                    "CFSAN",
                    Some(&gap_regions),
                    None,
                );

                info!(
                    "VCF comparison: {:.1}% concordance, {} Snippy-only, {} CFSAN-only",
                    vcf_comparison.summary.concordance_rate,
                    vcf_comparison.summary.only_pipeline_a,
                    vcf_comparison.summary.only_pipeline_b
                );

                // BAM validation to detect Snippy artifacts
                // Uses the same directory as VCF files since BAMs are alongside VCFs
                let snippy_only_positions = vcf_comparison.get_snippy_only_positions();
                if !snippy_only_positions.is_empty() {
                    info!("Validating {} Snippy-only positions against BAM files...", snippy_only_positions.len());
                    match bam_validate::validate_snippy_only_positions(&snippy_only_positions, snippy_vcf_dir) {
                        Ok(validation) => {
                            if validation.artifacts > 0 {
                                info!(
                                    "BAM validation: {} artifacts detected ({:.1}%), corrected Snippy-only: {} -> {}",
                                    validation.artifacts,
                                    validation.artifact_rate,
                                    validation.original_snippy_only,
                                    validation.corrected_snippy_only
                                );
                            } else {
                                info!("BAM validation: no artifacts detected, all {} Snippy-only positions confirmed", validation.real_variants);
                            }
                            vcf_comparison.set_bam_validation(validation);
                        }
                        Err(e) => {
                            log::warn!("BAM validation failed: {}", e);
                        }
                    }
                }

                // Store full comparison for server (with all positions for pagination)
                vcf_comparison_full = Some(vcf_comparison.clone());
                // Add summary to report (without all positions - uses serde skip)
                report = report.with_vcf_comparison(vcf_comparison);
            }
        } else {
            log::warn!("VCF comparison requires both --vcf-snippy and --vcf-cfsan");
        }
    }

    // Add execution metadata
    let duration_secs = start_time.elapsed().as_secs_f64();
    let input_mode = if fastq_mode { "FASTQ" } else { "FASTA" }.to_string();
    let input_path = args.input.as_ref()
        .map(|p| p.display().to_string())
        .or_else(|| args.fastq_input.as_ref().map(|p| p.display().to_string()))
        .unwrap_or_else(|| "individual files".to_string());
    let reference_path = args.reference.as_ref()
        .map(|p| p.display().to_string())
        .unwrap_or_default();

    let report = report.with_metadata(
        command_line,
        timestamp,
        duration_secs,
        input_mode,
        input_path,
        reference_path,
    );

    // Write output
    output::write_report(&report, &args.output)?;
    info!("Report written to: {}", args.output.display());

    // Write BED files if requested
    if let Some(ref bed_dir) = args.bed_output {
        output::write_bed_files(&report, bed_dir)?;
    }

    // Print summary
    let flagged = report.samples.iter().filter(|s| s.flagged).count();
    info!(
        "Summary: {}/{} samples flagged for review",
        flagged,
        report.samples.len()
    );

    // Start web server if requested
    #[cfg(feature = "serve")]
    if args.serve {
        // Load Snippy gaps from aligned FASTAs (N positions)
        let snippy_gaps = if let Some(ref snippy_vcf_dir) = args.vcf_snippy {
            let gaps = gaps::load_snippy_gaps(snippy_vcf_dir);
            if !gaps.is_empty() {
                info!("Loaded Snippy gaps for {} samples", gaps.len());
            }
            gaps
        } else {
            std::collections::HashMap::new()
        };

        // Load CFSAN gaps from pileup files (zero coverage positions)
        let cfsan_gaps = if let Some(ref cfsan_vcf_dir) = args.vcf_cfsan {
            let gaps = gaps::load_cfsan_gaps(cfsan_vcf_dir);
            if !gaps.is_empty() {
                info!("Loaded CFSAN gaps for {} samples", gaps.len());
            }
            gaps
        } else {
            std::collections::HashMap::new()
        };

        info!("Starting web server with IGV.js genome browser...");
        serve::start_server(
            &report,
            reference.as_ref(),
            args.bed_output.as_deref(),
            args.port,
            &pairwise_results,
            vcf_comparison_full.as_ref(),
            &snippy_gaps,
            &cfsan_gaps,
            args.bam_dir.as_deref(),
        )?;
    }

    Ok(())
}

/// Serve an existing report without re-running analysis
#[cfg(feature = "serve")]
fn serve_existing_report(
    report_path: &std::path::Path,
    reference_path: Option<&std::path::Path>,
    bed_dir: Option<&std::path::Path>,
    port: u16,
    bam_dir: Option<&std::path::Path>,
) -> Result<()> {
    use std::fs::File;
    use std::io::BufReader;

    info!("Loading existing report: {}", report_path.display());

    // Load the report JSON
    let file = File::open(report_path)
        .with_context(|| format!("Failed to open report: {}", report_path.display()))?;
    let reader = BufReader::new(file);
    let report: Report = serde_json::from_reader(reader)
        .with_context(|| format!("Failed to parse report JSON: {}", report_path.display()))?;

    info!("Report loaded: {} samples", report.samples.len());

    // Optionally load reference
    let reference: Option<Sample> = if let Some(ref_path) = reference_path {
        info!("Loading reference: {}", ref_path.display());
        Some(fasta::load_sample(ref_path)?)
    } else {
        // Try to infer reference from report directory
        let report_dir = report_path.parent();
        if let Some(dir) = report_dir {
            // Look for common reference file names
            for name in &["reference.fa", "reference.fasta", "ref.fa", "ref.fasta"] {
                let ref_path = dir.join(name);
                if ref_path.exists() {
                    info!("Found reference: {}", ref_path.display());
                    match fasta::load_sample(&ref_path) {
                        Ok(sample) => return serve_with_reference(report, Some(sample), bed_dir, port, bam_dir),
                        Err(e) => log::warn!("Failed to load {}: {}", ref_path.display(), e),
                    }
                }
            }
        }
        log::warn!("No reference provided - IGV.js genome browser may not work correctly");
        None
    };

    serve_with_reference(report, reference, bed_dir, port, bam_dir)
}

#[cfg(feature = "serve")]
fn serve_with_reference(
    report: Report,
    reference: Option<Sample>,
    bed_dir: Option<&std::path::Path>,
    port: u16,
    bam_dir: Option<&std::path::Path>,
) -> Result<()> {
    info!("Starting web server with IGV.js genome browser...");
    info!("Report summary: {} total samples, {} flagged",
          report.summary.total_samples,
          report.summary.flagged_samples);

    if !report.impact_analysis.samples.is_empty() {
        info!("Impact analysis available for {} samples", report.impact_analysis.samples.len());
    }

    // Create empty pairwise results for serving existing report
    let empty_pairwise = pairwise::PairwiseResults::empty();
    let empty_snippy_gaps = std::collections::HashMap::new();
    let empty_cfsan_gaps = std::collections::HashMap::new();

    // Note: VCF comparison full data not available when loading existing report
    // (the full positions are skipped in JSON serialization)
    serve::start_server(
        &report,
        reference.as_ref(),
        bed_dir,
        port,
        &empty_pairwise,
        None,
        &empty_snippy_gaps,
        &empty_cfsan_gaps,
        bam_dir,
    )
}

fn load_samples(args: &Args) -> Result<Vec<Sample>> {
    let paths: Vec<PathBuf> = if let Some(ref input_dir) = args.input {
        // Load all FASTA files from directory (including gzipped)
        std::fs::read_dir(input_dir)?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| {
                let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("");
                // Match .fa, .fasta, .fna and their gzipped versions
                name.ends_with(".fa")
                    || name.ends_with(".fasta")
                    || name.ends_with(".fna")
                    || name.ends_with(".fa.gz")
                    || name.ends_with(".fasta.gz")
                    || name.ends_with(".fna.gz")
            })
            .collect()
    } else if let Some(ref sample_paths) = args.samples {
        sample_paths.clone()
    } else {
        anyhow::bail!("Either --input or --samples must be specified");
    };

    info!("Loading {} FASTA files...", paths.len());

    paths.iter().map(|p| fasta::load_sample(p)).collect()
}

/// Load FASTQ samples and analyze read coverage against reference
fn load_fastq_samples_and_analyze(
    args: &Args,
    reference: &Sample,
) -> Result<(Vec<Sample>, Vec<reads::ReadMappingResult>)> {
    // Load FASTQ samples
    let read_samples = if let Some(ref fastq_dir) = args.fastq_input {
        reads::load_fastq_samples(fastq_dir, &args.read_type)?
    } else if let Some(ref fastq_paths) = args.fastq_samples {
        // Load from individual file paths
        fastq_paths
            .iter()
            .filter_map(|p| {
                if reads::is_fastq_r1(p) {
                    reads::load_sample_from_file(p, &args.read_type).ok().flatten()
                } else {
                    None
                }
            })
            .collect()
    } else {
        anyhow::bail!("Either --fastq-input or --fastq-samples must be specified for FASTQ mode");
    };

    info!("Found {} FASTQ samples", read_samples.len());

    // Analyze each sample's coverage
    let mut mapping_results = Vec::new();
    let mut virtual_samples = Vec::new();

    for sample in &read_samples {
        info!("Analyzing coverage for {}...", sample.name);

        let result = reads::analyze_read_coverage(
            sample,
            reference,
            args.min_depth,
            args.threads,
            args.minimap2_path.as_deref(),
        )?;

        info!(
            "  {} - Depth: {:.1}x, Coverage: {:.1}%, Gaps: {} regions ({} bp)",
            result.name,
            result.mean_depth,
            result.coverage_percent,
            result.gap_regions.len(),
            result.total_gap_bases
        );

        // Create virtual sample for compatibility with existing analysis
        let virtual_sample = reads::read_result_to_sample(&result, reference);
        virtual_samples.push(virtual_sample);
        mapping_results.push(result);
    }

    Ok((virtual_samples, mapping_results))
}
