//! Compare module - pipeline-agnostic SNP comparison
//!
//! Processes VCF/BAM files from multiple pipelines and generates
//! a compact JSON report for visualization.

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;
use std::process::Command;

use crate::config::Config;
use crate::pileup::{self, BaseCall, BaseCallWithStats};

/// Compact report structure for WASM visualization
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CompareReport {
    /// Schema version
    #[serde(rename = "_version")]
    pub version: String,

    /// Reference genome info
    pub reference: ReferenceInfo,

    /// Sample metadata
    pub samples: HashMap<String, SampleInfo>,

    /// Pipeline metadata
    pub pipelines: HashMap<String, PipelineInfo>,

    /// Data: sample -> pipeline -> gaps/snps
    pub data: HashMap<String, HashMap<String, PipelineData>>,

    /// Polymorphic sites for distance matrix calculation (per pipeline)
    /// Key: pipeline_id -> position (0-based) -> allele data per sample
    #[serde(skip_serializing_if = "HashMap::is_empty", default)]
    pub polymorphic_sites: HashMap<String, HashMap<u32, PolymorphicSite>>,

    /// Summary statistics
    pub summary: Summary,
}

/// Allele information at a polymorphic site
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolymorphicSite {
    /// Reference allele at this position
    #[serde(rename = "ref")]
    pub ref_allele: char,
    /// Allele for each sample (key: sample_id, value: allele info)
    pub alleles: HashMap<String, SampleAllele>,
}

/// Allele call for a single sample at a polymorphic site
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleAllele {
    /// The base called (A, C, G, T, or N for ambiguous)
    pub base: char,
    /// Source of the call: "vcf", "bam", or "gap"
    pub source: String,
    /// Depth at this position (if available)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub depth: Option<u32>,
    /// Quality score from VCF (if available)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub qual: Option<f64>,
    /// Consensus percentage from BAM (0.0 - 1.0, if available)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub consensus: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceInfo {
    pub name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    pub length: usize,
    pub sequence: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleInfo {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineInfo {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    pub has_vcf: bool,
    pub has_bam: bool,
    /// Mark this pipeline as ground truth (baseline for comparison)
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub ground_truth: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PipelineData {
    /// Gap regions [start, end)
    #[serde(skip_serializing_if = "Vec::is_empty", default)]
    pub gaps: Vec<[usize; 2]>,

    /// SNP list
    #[serde(skip_serializing_if = "Vec::is_empty", default)]
    pub snps: Vec<Snp>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Snp {
    pub pos: usize,
    #[serde(rename = "ref")]
    pub ref_allele: String,
    pub alt: String,
    pub qual: f64,
    pub dp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Summary {
    pub total_samples: usize,
    pub total_pipelines: usize,
    pub generated_at: String,
    pub coreguard_version: String,
    /// Warnings (e.g., MNP decomposition notices)
    #[serde(skip_serializing_if = "Vec::is_empty", default)]
    pub warnings: Vec<String>,
}

impl CompareReport {
    /// Generate report from configuration
    pub fn from_config(config: &Config) -> Result<Self> {
        let sample_ids = config.all_sample_ids();
        let pipeline_ids = config.all_pipeline_ids();

        // Load reference
        log::info!("Loading reference: {}", config.reference.path);
        let ref_content = std::fs::read_to_string(&config.reference.path)
            .with_context(|| format!("Failed to read reference: {}", config.reference.path))?;

        let (ref_name, ref_seq) = parse_fasta(&ref_content);
        let ref_length = ref_seq.len();

        log::info!(
            "Reference: {} ({} bp)",
            ref_name,
            ref_length
        );

        // Build sample info
        let mut samples = HashMap::new();
        for sample_id in &sample_ids {
            samples.insert(
                sample_id.clone(),
                SampleInfo {
                    label: config.samples.get(sample_id).and_then(|s| s.label.clone()),
                },
            );
        }

        // Build pipeline info
        let mut pipelines = HashMap::new();
        for pipeline_id in &pipeline_ids {
            if let Some(pipeline) = config.pipelines.get(pipeline_id) {
                let has_vcf = pipeline.samples.values().any(|f| f.vcf.is_some());
                let has_bam = pipeline.samples.values().any(|f| f.bam.is_some());
                pipelines.insert(
                    pipeline_id.clone(),
                    PipelineInfo {
                        label: pipeline.label.clone(),
                        has_vcf,
                        has_bam,
                        ground_truth: pipeline.ground_truth,
                    },
                );
            }
        }

        // Process data for each sample and pipeline
        let mut data: HashMap<String, HashMap<String, PipelineData>> = HashMap::new();
        let mut total_mnps_found = 0;
        let mut total_snps_from_mnps = 0;

        for sample_id in &sample_ids {
            let mut sample_data: HashMap<String, PipelineData> = HashMap::new();

            for pipeline_id in &pipeline_ids {
                if let Some(pipeline) = config.pipelines.get(pipeline_id) {
                    if let Some(files) = pipeline.samples.get(sample_id) {
                        let mut pipeline_data = PipelineData::default();

                        // Load gaps from BAM
                        if let Some(bam_path) = &files.bam {
                            log::info!(
                                "Loading gaps for {}/{} from {}",
                                sample_id,
                                pipeline_id,
                                bam_path
                            );
                            pipeline_data.gaps =
                                load_gaps_from_bam(bam_path, ref_length, config.options.min_depth)?;
                            log::info!(
                                "  Found {} gap regions",
                                pipeline_data.gaps.len()
                            );
                        }

                        // Load SNPs from VCF
                        if let Some(vcf_path) = &files.vcf {
                            log::info!(
                                "Loading SNPs for {}/{} from {}",
                                sample_id,
                                pipeline_id,
                                vcf_path
                            );
                            let result = load_snps_from_vcf(
                                vcf_path,
                                config.options.min_qual,
                                config.options.include_indels,
                            )?;
                            pipeline_data.snps = result.snps;
                            total_mnps_found += result.mnps_found;
                            total_snps_from_mnps += result.snps_from_mnps;
                            log::info!("  Found {} SNPs", pipeline_data.snps.len());
                        }

                        sample_data.insert(pipeline_id.clone(), pipeline_data);
                    }
                }
            }

            data.insert(sample_id.clone(), sample_data);
        }

        // Build warnings
        let mut warnings = Vec::new();
        if total_mnps_found > 0 {
            warnings.push(format!(
                "Detected {} MNP(s) (multi-nucleotide polymorphisms) which were automatically decomposed into {} individual SNPs",
                total_mnps_found, total_snps_from_mnps
            ));
        }

        // Build polymorphic sites for distance matrix calculation
        log::info!("Building polymorphic sites for distance matrix...");
        let polymorphic_sites = build_polymorphic_sites(
            &data,
            &sample_ids,
            &pipeline_ids,
            config,
            &ref_name,
            &ref_seq,
        )?;
        let total_sites: usize = polymorphic_sites.values().map(|v| v.len()).max().unwrap_or(0);
        log::info!("Built polymorphic sites for {} pipelines ({} positions each)", polymorphic_sites.len(), total_sites);

        // Build summary
        let summary = Summary {
            total_samples: sample_ids.len(),
            total_pipelines: pipeline_ids.len(),
            generated_at: chrono::Utc::now().to_rfc3339(),
            coreguard_version: env!("CARGO_PKG_VERSION").to_string(),
            warnings,
        };

        Ok(CompareReport {
            version: "1.0".to_string(),
            reference: ReferenceInfo {
                name: ref_name.clone(),
                label: config.reference.label.clone(),
                length: ref_length,
                sequence: ref_seq,
            },
            samples,
            pipelines,
            data,
            polymorphic_sites,
            summary,
        })
    }

    /// Save report to JSON file
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let json = serde_json::to_string_pretty(self)?;
        std::fs::write(path.as_ref(), json)?;
        Ok(())
    }

    /// Save report to compact JSON (no pretty print, smaller file)
    pub fn save_compact<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let json = serde_json::to_string(self)?;
        std::fs::write(path.as_ref(), json)?;
        Ok(())
    }

    /// Save report to gzipped JSON file
    pub fn save_gzip<P: AsRef<Path>>(&self, path: P, compact: bool) -> Result<()> {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;

        let file = std::fs::File::create(path.as_ref())?;
        let mut encoder = GzEncoder::new(file, Compression::default());

        let json = if compact {
            serde_json::to_string(self)?
        } else {
            serde_json::to_string_pretty(self)?
        };

        encoder.write_all(json.as_bytes())?;
        encoder.finish()?;
        Ok(())
    }
}

/// Parse FASTA file, return (name, sequence)
fn parse_fasta(content: &str) -> (String, String) {
    let mut name = String::new();
    let mut seq = String::new();

    for line in content.lines() {
        if line.starts_with('>') {
            name = line[1..].split_whitespace().next().unwrap_or("").to_string();
        } else {
            seq.push_str(line.trim());
        }
    }

    (name, seq)
}

/// Load gaps from BAM using samtools depth
fn load_gaps_from_bam(bam_path: &str, ref_len: usize, min_depth: usize) -> Result<Vec<[usize; 2]>> {
    let output = Command::new("samtools")
        .args(["depth", "-a", bam_path])
        .output()
        .with_context(|| format!("Failed to run samtools depth on {}", bam_path))?;

    if !output.status.success() {
        anyhow::bail!(
            "samtools depth failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let mut gaps = Vec::new();
    let mut in_gap = false;
    let mut gap_start = 0;

    for line in stdout.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 3 {
            let pos: usize = fields[1].parse().unwrap_or(0);
            let depth: usize = fields[2].parse().unwrap_or(0);

            if depth < min_depth {
                if !in_gap {
                    gap_start = pos;
                    in_gap = true;
                }
            } else if in_gap {
                gaps.push([gap_start, pos]);
                in_gap = false;
            }
        }
    }

    if in_gap {
        gaps.push([gap_start, ref_len + 1]);
    }

    Ok(gaps)
}

/// Result of loading SNPs including MNP statistics
pub struct LoadSnpsResult {
    pub snps: Vec<Snp>,
    pub mnps_found: usize,
    pub snps_from_mnps: usize,
}

/// Load SNPs from VCF file
fn load_snps_from_vcf(vcf_path: &str, min_qual: f64, include_indels: bool) -> Result<LoadSnpsResult> {
    let content = std::fs::read_to_string(vcf_path)
        .with_context(|| format!("Failed to read VCF: {}", vcf_path))?;

    let mut snps = Vec::new();
    let mut mnps_found = 0;
    let mut snps_from_mnps = 0;

    for line in content.lines() {
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            continue;
        }

        let pos: usize = fields[1].parse().unwrap_or(0);
        let ref_allele = fields[3];
        let alt_allele = fields[4];

        // QUAL can be "." (not available) - treat as passing filter
        let qual_str = fields[5];
        let qual: f64 = if qual_str == "." {
            f64::MAX  // Pass quality filter when not available
        } else {
            qual_str.parse().unwrap_or(0.0)
        };

        // Skip low quality (but "." passes through)
        if qual < min_qual {
            continue;
        }

        // For storage, use a reasonable default when QUAL was "."
        let qual_for_storage = if qual == f64::MAX { 0.0 } else { qual };

        // Extract depth from INFO or FORMAT/SAMPLE
        let info = fields[7];
        let format = fields.get(8).copied();
        let sample = fields.get(9).copied();
        let dp = extract_depth(info, format, sample);

        // Handle MNPs (multi-nucleotide polymorphisms) - decompose into individual SNPs
        // MNPs have same length ref and alt, both > 1 (e.g., TTGGCG → CCGGCT)
        if ref_allele.len() == alt_allele.len() && ref_allele.len() > 1 {
            mnps_found += 1;
            let mut decomposed_count = 0;
            for (i, (r, a)) in ref_allele.chars().zip(alt_allele.chars()).enumerate() {
                if r != a {
                    snps.push(Snp {
                        pos: pos + i,
                        ref_allele: r.to_string(),
                        alt: a.to_string(),
                        qual: qual_for_storage,
                        dp,
                    });
                    decomposed_count += 1;
                }
            }
            snps_from_mnps += decomposed_count;
            log::debug!(
                "Decomposed MNP {}→{} at pos {} into {} SNPs",
                ref_allele, alt_allele, pos, decomposed_count
            );
        } else if ref_allele.len() != alt_allele.len() {
            // Indel (different lengths) - skip if not requested
            if !include_indels {
                continue;
            }
            snps.push(Snp {
                pos,
                ref_allele: ref_allele.to_string(),
                alt: alt_allele.to_string(),
                qual: qual_for_storage,
                dp,
            });
        } else {
            // Regular SNP (single base)
            snps.push(Snp {
                pos,
                ref_allele: ref_allele.to_string(),
                alt: alt_allele.to_string(),
                qual: qual_for_storage,
                dp,
            });
        }
    }

    if mnps_found > 0 {
        log::info!(
            "VCF {}: Decomposed {} MNPs into {} individual SNPs",
            vcf_path, mnps_found, snps_from_mnps
        );
    }

    Ok(LoadSnpsResult {
        snps,
        mnps_found,
        snps_from_mnps,
    })
}

/// Extract depth from INFO field or FORMAT/SAMPLE columns
fn extract_depth(info: &str, format: Option<&str>, sample: Option<&str>) -> usize {
    // Try INFO field first (DP=XX or ADP=XX for VarScan)
    for part in info.split(';') {
        if part.starts_with("DP=") {
            if let Ok(dp) = part[3..].parse() {
                return dp;
            }
        }
        if part.starts_with("ADP=") {
            if let Ok(dp) = part[4..].parse() {
                return dp;
            }
        }
    }

    // Try FORMAT/SAMPLE (DP in format, value in sample)
    if let (Some(fmt), Some(smp)) = (format, sample) {
        let format_fields: Vec<&str> = fmt.split(':').collect();
        let sample_fields: Vec<&str> = smp.split(':').collect();

        // Look for DP field
        for (i, field) in format_fields.iter().enumerate() {
            if *field == "DP" || *field == "SDP" {
                if let Some(value) = sample_fields.get(i) {
                    if let Ok(dp) = value.parse() {
                        return dp;
                    }
                }
            }
        }
    }

    0
}

/// Build polymorphic sites map for distance matrix calculation
///
/// Build polymorphic sites for each pipeline separately.
/// VCF SNP info for distance matrix calculation
#[derive(Debug, Clone)]
struct VcfSnpInfo {
    base: char,
    depth: Option<u32>,
    qual: Option<f64>,
}

/// Only includes "core" polymorphic positions where:
/// 1. All samples have coverage (no gaps)
/// 2. At least one sample differs from another (truly polymorphic)
/// Returns: pipeline_id -> position -> PolymorphicSite
fn build_polymorphic_sites(
    data: &HashMap<String, HashMap<String, PipelineData>>,
    sample_ids: &[String],
    pipeline_ids: &[String],
    config: &Config,
    ref_name: &str,
    ref_seq: &str,
) -> Result<HashMap<String, HashMap<u32, PolymorphicSite>>> {
    use std::collections::HashSet;

    // Step 1: Collect SNP positions PER PIPELINE (not across all pipelines)
    // pipeline -> position -> (sample -> VcfSnpInfo)
    let mut pipeline_snp_map: HashMap<String, HashMap<u32, HashMap<String, VcfSnpInfo>>> = HashMap::new();

    for (sample_id, pipelines_data) in data {
        for (pipeline_id, pipeline_data) in pipelines_data {
            let pipeline_positions = pipeline_snp_map
                .entry(pipeline_id.clone())
                .or_insert_with(HashMap::new);

            for snp in &pipeline_data.snps {
                let pos = (snp.pos - 1) as u32;  // Convert to 0-based
                let alt_char = snp.alt.chars().next().unwrap_or('N');

                pipeline_positions
                    .entry(pos)
                    .or_insert_with(HashMap::new)
                    .insert(sample_id.clone(), VcfSnpInfo {
                        base: alt_char,
                        depth: if snp.dp > 0 { Some(snp.dp as u32) } else { None },
                        qual: if snp.qual > 0.0 { Some(snp.qual) } else { None },
                    });
            }
        }
    }

    // Step 2: Build polymorphic sites for EACH pipeline using only THAT pipeline's positions
    let mut result: HashMap<String, HashMap<u32, PolymorphicSite>> = HashMap::new();

    for pipeline_id in pipeline_ids {
        // Get positions for THIS pipeline only
        let pipeline_snps = match pipeline_snp_map.get(pipeline_id) {
            Some(snps) => snps,
            None => {
                result.insert(pipeline_id.clone(), HashMap::new());
                continue;
            }
        };

        let mut positions: Vec<u32> = pipeline_snps.keys().copied().collect();
        positions.sort();

        log::info!("Building polymorphic sites for pipeline {}: {} candidate positions", pipeline_id, positions.len());

        // Collect gap positions for each sample (for faster lookup)
        let mut sample_gap_positions: HashMap<String, HashSet<u32>> = HashMap::new();
        for sample_id in sample_ids {
            let gap_regions: Vec<[usize; 2]> = data.get(sample_id)
                .and_then(|p| p.get(pipeline_id))
                .map(|d| d.gaps.clone())
                .unwrap_or_default();

            let mut gap_set = HashSet::new();
            for [start, end] in gap_regions {
                for p in start..end {
                    gap_set.insert(p as u32);
                }
            }
            sample_gap_positions.insert(sample_id.clone(), gap_set);
        }

        // Pre-read BAM bases with stats for all samples
        let mut sample_bam_bases: HashMap<String, HashMap<u32, BaseCallWithStats>> = HashMap::new();
        for sample_id in sample_ids {
            let bam_path: Option<String> = config.pipelines.get(pipeline_id)
                .and_then(|p| p.samples.get(sample_id))
                .and_then(|f| f.bam.clone());

            let bam_bases: HashMap<u32, BaseCallWithStats> = if let Some(ref path) = bam_path {
                match pileup::get_bases_with_stats(
                    Path::new(path),
                    ref_name,
                    &positions,
                    config.options.min_depth as u32,
                ) {
                    Ok(bases) => bases,
                    Err(e) => {
                        log::warn!("Failed to read BAM {}: {}", path, e);
                        HashMap::new()
                    }
                }
            } else {
                HashMap::new()
            };
            sample_bam_bases.insert(sample_id.clone(), bam_bases);
        }

        // Build sites and filter to core polymorphic positions
        let mut pipeline_sites: HashMap<u32, PolymorphicSite> = HashMap::new();
        let mut core_count = 0u32;

        for &pos in &positions {
            let ref_char = ref_seq.chars().nth(pos as usize).unwrap_or('N');
            let mut alleles: HashMap<String, SampleAllele> = HashMap::new();
            let mut all_have_coverage = true;

            // Collect alleles for all samples
            for sample_id in sample_ids {
                // Check if position is in a gap
                let in_gap = sample_gap_positions
                    .get(sample_id)
                    .map(|s| s.contains(&pos))
                    .unwrap_or(false);

                if in_gap {
                    all_have_coverage = false;
                    alleles.insert(sample_id.clone(), SampleAllele {
                        base: '-',
                        source: "gap".to_string(),
                        depth: Some(0),
                        qual: None,
                        consensus: None,
                    });
                    continue;
                }

                // Check if we have a SNP from this pipeline's VCF for this sample
                let snp_info: Option<&VcfSnpInfo> = pipeline_snps
                    .get(&pos)
                    .and_then(|sample_alts| sample_alts.get(sample_id));

                if let Some(info) = snp_info {
                    // Have SNP from VCF
                    alleles.insert(sample_id.clone(), SampleAllele {
                        base: info.base,
                        source: "vcf".to_string(),
                        depth: info.depth,
                        qual: info.qual,
                        consensus: None,  // VCF doesn't provide consensus
                    });
                } else if let Some(bam_stats) = sample_bam_bases.get(sample_id).and_then(|m| m.get(&pos)) {
                    // No SNP - get from BAM with stats
                    match &bam_stats.call {
                        BaseCall::Base(b) => {
                            alleles.insert(sample_id.clone(), SampleAllele {
                                base: *b,
                                source: "bam".to_string(),
                                depth: Some(bam_stats.depth),
                                qual: None,
                                consensus: Some(bam_stats.consensus),
                            });
                        }
                        BaseCall::Gap => {
                            all_have_coverage = false;
                            alleles.insert(sample_id.clone(), SampleAllele {
                                base: '-',
                                source: "gap".to_string(),
                                depth: Some(bam_stats.depth),
                                qual: None,
                                consensus: Some(0.0),
                            });
                        }
                        BaseCall::Ambiguous => {
                            alleles.insert(sample_id.clone(), SampleAllele {
                                base: 'N',
                                source: "ambiguous".to_string(),
                                depth: Some(bam_stats.depth),
                                qual: None,
                                consensus: Some(bam_stats.consensus),
                            });
                        }
                    }
                } else {
                    // No BAM data - assume reference (for vcf_ref mode this is what we need)
                    alleles.insert(sample_id.clone(), SampleAllele {
                        base: ref_char,
                        source: "inferred".to_string(),
                        depth: None,
                        qual: None,
                        consensus: None,
                    });
                }
            }

            // Check if this is a "core" position (all samples have coverage)
            if !all_have_coverage {
                continue;
            }
            core_count += 1;

            // Check if truly polymorphic using ACTUAL bases (BAM-validated)
            // This is pipeline-agnostic: uses real sequenced bases, not VCF+Ref assumptions
            let sample_bases: Vec<char> = sample_ids.iter()
                .filter_map(|sid| {
                    let allele = alleles.get(sid)?;
                    if allele.base == '-' || allele.base == 'N' {
                        return None;
                    }
                    // Use actual base (from VCF or BAM) - NOT reference assumption
                    Some(allele.base)
                })
                .collect();

            // Only include if all samples have valid bases AND at least one differs
            if sample_bases.len() != sample_ids.len() {
                continue;
            }

            let is_polymorphic = {
                let first = sample_bases[0];
                sample_bases.iter().any(|&b| b != first)
            };

            if is_polymorphic {
                pipeline_sites.insert(pos, PolymorphicSite {
                    ref_allele: ref_char,
                    alleles,
                });
            }
        }

        log::info!("Pipeline {}: {} core positions, {} polymorphic",
            pipeline_id, core_count, pipeline_sites.len());

        result.insert(pipeline_id.clone(), pipeline_sites);
    }

    Ok(result)
}
