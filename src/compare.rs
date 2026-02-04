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
    /// Key: pipeline_id -> position (as string for JSON compatibility) -> allele data per sample
    #[serde(default)]
    pub polymorphic_sites: HashMap<String, HashMap<String, PolymorphicSite>>,

    /// Summary statistics
    pub summary: Summary,

    /// Description (markdown content)
    #[serde(default)]
    pub description: Option<String>,

    /// Pre-computed distance matrices from pipelines (pipeline_id -> matrix data)
    #[serde(default)]
    pub pipeline_distance_matrices: HashMap<String, PipelineDistanceMatrix>,

    /// Pre-computed GT discriminating SNPs vs pipeline core SNP results
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gt_disc_vs_pipelines: Option<Vec<GtDiscVsPipelineResult>>,
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
    #[serde(default)]
    pub depth: Option<u32>,
    /// Quality score from VCF (if available)
    #[serde(default)]
    pub qual: Option<f64>,
    /// Consensus percentage from BAM (0.0 - 1.0, if available)
    #[serde(default)]
    pub consensus: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceInfo {
    pub name: String,
    #[serde(default)]
    pub label: Option<String>,
    pub length: usize,
    pub sequence: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleInfo {
    #[serde(default)]
    pub label: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineInfo {
    #[serde(default)]
    pub label: Option<String>,
    #[serde(default)]
    pub command: Option<String>,
    pub has_vcf: bool,
    pub has_bam: bool,
    /// Mark this pipeline as ground truth (baseline for comparison)
    #[serde(default)]
    pub ground_truth: bool,
    /// True if SNPs were derived from BAM pileup (no variant calling filters)
    #[serde(default)]
    pub from_bam_pileup: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PipelineData {
    /// Gap regions [start, end)
    #[serde(default)]
    pub gaps: Vec<[usize; 2]>,

    /// SNP list
    #[serde(default)]
    pub snps: Vec<Snp>,

    /// Source VCF file path (for reproducibility)
    #[serde(default)]
    pub vcf_path: Option<String>,

    /// Source BAM file path (for reproducibility)
    #[serde(default)]
    pub bam_path: Option<String>,
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
    #[serde(default)]
    pub warnings: Vec<String>,
    /// SNPs in ground truth gaps statistics (pipeline_id -> stats) - DEPRECATED, use snps_in_gaps
    #[serde(default)]
    pub snps_in_gt_gaps: Option<HashMap<String, SnpsInGapsStats>>,
    /// SNPs in gaps for ALL pipeline pairs (gap_pipeline -> snp_pipeline -> stats)
    #[serde(default)]
    pub snps_in_gaps: Option<HashMap<String, HashMap<String, SnpsInGapsStats>>>,
    /// Ground truth pileup SNP statistics (raw count from BAM without variant calling)
    #[serde(default)]
    pub ground_truth_pileup: Option<GroundTruthPileupStats>,
    /// Per-pipeline MNP decomposition statistics
    #[serde(default)]
    pub mnp_stats: Option<HashMap<String, MnpStats>>,
}

/// Ground truth SNP count from BAM pileup (without variant calling)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GroundTruthPileupStats {
    /// Total SNPs counted from ground truth BAM pileup (all samples combined)
    pub total_snps: usize,
    /// SNP count per sample (sample_id -> count)
    pub per_sample: HashMap<String, usize>,
    /// Covered positions (where we could make a call)
    pub covered_positions: usize,
    /// Comparison with VCF pipelines (pipeline_id -> comparison)
    pub pipeline_comparison: HashMap<String, PipelineVsGroundTruth>,
}

/// Comparison of a pipeline's SNP count vs ground truth pileup
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineVsGroundTruth {
    /// Pipeline SNP count (from VCF)
    pub pipeline_snps: usize,
    /// Ground truth SNP count (from BAM pileup)
    pub ground_truth_snps: usize,
    /// Numerical difference (pipeline - ground_truth)
    pub difference: i64,
    /// Percentage difference ((pipeline - ground_truth) / ground_truth * 100)
    pub percentage_diff: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnpsInGapsStats {
    pub total_snps: usize,
    pub snps_in_gaps: usize,
    pub percentage: f64,
}

/// MNP (Multi-Nucleotide Polymorphism) statistics per pipeline
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MnpStats {
    /// Number of MNPs found and decomposed
    pub mnps_found: usize,
    /// Total individual SNPs resulting from MNP decomposition
    pub snps_from_mnps: usize,
}

/// Pre-computed distance matrix from a pipeline
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineDistanceMatrix {
    /// Sample IDs in order (matches row/column order in matrix)
    pub samples: Vec<String>,
    /// Distance matrix (row-major, samples[i] vs samples[j])
    pub matrix: Vec<Vec<i64>>,
}

// ============================================================================
// GT Disc vs Pipelines (pre-computed from core SNP files)
// ============================================================================

/// Core SNP data from a pipeline's native output file (snippycore.tab or snplist.txt)
#[derive(Debug, Clone)]
pub struct CoreSnpData {
    /// Pipeline positions with per-sample alleles
    pub positions: Vec<CoreSnpPosition>,
    /// True for snippycore.tab (has alleles), false for snplist.txt (positions only)
    pub has_alleles: bool,
}

/// A single position from a core SNP file
#[derive(Debug, Clone)]
pub struct CoreSnpPosition {
    /// 0-based position
    pub pos: usize,
    /// Reference allele (if available)
    pub ref_allele: Option<String>,
    /// Sample → allele (empty if !has_alleles)
    pub alleles: HashMap<String, String>,
    /// Samples that have a SNP at this position
    pub samples_with_snp: Vec<String>,
}

/// Parse a core SNP file, auto-detecting format (snippycore.tab vs snplist.txt)
pub fn parse_core_snps(path: &str) -> anyhow::Result<CoreSnpData> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("Failed to read core SNPs file: {}", path))?;

    let first_line = content.lines().next().unwrap_or("");

    if first_line.starts_with("CHR\t") || first_line.starts_with("CHR ") {
        // snippycore.tab format
        parse_snippycore_tab_for_core(path)
    } else {
        // snplist.txt format (CFSAN)
        parse_cfsan_snplist_for_core(path)
    }
}

/// Parse snippycore.tab format: CHR\tPOS\tREF\tsample1\tsample2\t...
fn parse_snippycore_tab_for_core(path: &str) -> anyhow::Result<CoreSnpData> {
    use std::io::{BufRead, BufReader};

    let file = std::fs::File::open(path)
        .with_context(|| format!("Failed to open core SNPs file: {}", path))?;
    let reader = BufReader::new(file);

    let mut positions = Vec::new();
    let mut sample_names: Vec<String> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();

        if line.starts_with("CHR\t") {
            // Header line — extract sample names
            sample_names = parts.iter()
                .skip(3) // Skip CHR, POS, REF
                .map(|s| s.trim_end_matches("_snippy").to_string())
                .collect();
            continue;
        }

        if parts.len() < 4 || sample_names.is_empty() {
            continue;
        }

        let pos: usize = parts[1].parse().unwrap_or(0);
        if pos == 0 { continue; }
        let pos = pos - 1; // Convert to 0-based

        let ref_allele = parts[2].to_string();

        let mut alleles = HashMap::new();
        let mut samples_with_snp = Vec::new();

        for (i, sample) in sample_names.iter().enumerate() {
            if let Some(allele) = parts.get(3 + i) {
                let allele = allele.to_string();
                // A sample has a SNP if its allele differs from ref and isn't N or -
                if allele != ref_allele && allele != "N" && allele != "-" {
                    samples_with_snp.push(sample.clone());
                }
                alleles.insert(sample.clone(), allele);
            }
        }

        positions.push(CoreSnpPosition {
            pos,
            ref_allele: Some(ref_allele),
            alleles,
            samples_with_snp,
        });
    }

    log::info!("Parsed {} core SNP positions from snippycore.tab (with alleles)", positions.len());

    Ok(CoreSnpData {
        positions,
        has_alleles: true,
    })
}

/// Parse CFSAN snplist.txt format: CHROM\tPOS\tCOUNT\tsample1\tsample2\t...
fn parse_cfsan_snplist_for_core(path: &str) -> anyhow::Result<CoreSnpData> {
    use std::io::{BufRead, BufReader};

    let file = std::fs::File::open(path)
        .with_context(|| format!("Failed to open core SNPs file: {}", path))?;
    let reader = BufReader::new(file);

    let mut positions = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 {
            continue;
        }

        let pos: usize = parts[1].parse().unwrap_or(0);
        if pos == 0 { continue; }
        let pos = pos - 1; // Convert to 0-based

        let samples_with_snp: Vec<String> = parts.iter()
            .skip(3)
            .map(|s| s.to_string())
            .collect();

        positions.push(CoreSnpPosition {
            pos,
            ref_allele: None,
            alleles: HashMap::new(),
            samples_with_snp,
        });
    }

    log::info!("Parsed {} core SNP positions from snplist.txt (no alleles)", positions.len());

    Ok(CoreSnpData {
        positions,
        has_alleles: false,
    })
}

/// Result of GT disc vs pipeline comparison
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GtDiscVsPipelineResult {
    pub pipeline_id: String,
    /// Total core SNP positions in the pipeline's core_snps file
    pub pl_total_core_snps: u32,
    pub gap_intersect: GapStrategyResult,
    pub gap_union: GapStrategyResult,
    pub pairwise: PairwiseGtDiscResult,
}

/// Results for a single gap strategy (intersect or union)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GapStrategyResult {
    pub gt_disc: u32,
    pub same_pos: u32,
    /// None if pipeline has no allele data (CFSAN snplist.txt)
    pub concordant: Option<u32>,
    pub pl_snps_in_gt_gaps: u32,
}

/// Pairwise results (averages across all sample pairs)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairwiseGtDiscResult {
    pub gt_disc_avg: f64,
    pub same_pos_avg: f64,
    pub concordant_avg: Option<f64>,
    pub num_pairs: u32,
}

/// Compute GT disc vs pipeline metrics using core SNP data.
///
/// For each non-GT pipeline with `core_snps`, computes:
/// - Gap-Intersect and Gap-Union strategies for GT discriminating positions
/// - same_pos: GT disc positions present in pipeline's core SNP file
/// - concordant: among same_pos, alleles match GT (only if has_alleles)
/// - pl_snps_in_gt_gaps: pipeline core SNP positions falling in GT gaps
pub fn compute_gt_disc_vs_pipelines(
    data: &HashMap<String, HashMap<String, PipelineData>>,
    sample_ids: &[String],
    gt_pipeline_id: &str,
    core_snp_data: &HashMap<String, CoreSnpData>,
    ref_seq: &str,
) -> Vec<GtDiscVsPipelineResult> {
    use std::collections::HashSet;

    let n = sample_ids.len();
    if n < 2 {
        return Vec::new();
    }

    let ref_bytes = ref_seq.as_bytes();

    // Build per-sample GT gap sets
    let mut gt_gap_sets: Vec<HashSet<usize>> = Vec::new();
    for sample_id in sample_ids {
        let mut gaps = HashSet::new();
        if let Some(sample_data) = data.get(sample_id) {
            if let Some(gt_data) = sample_data.get(gt_pipeline_id) {
                for gap in &gt_data.gaps {
                    for pos in gap[0]..gap[1] {
                        gaps.insert(pos);
                    }
                }
            }
        }
        gt_gap_sets.push(gaps);
    }

    // GT gap union (any sample has GT gap)
    let mut gt_gap_union: HashSet<usize> = HashSet::new();
    for gs in &gt_gap_sets {
        gt_gap_union.extend(gs);
    }

    // GT gap intersection (all samples have GT gap)
    let gt_gap_intersection = if gt_gap_sets.len() > 1 {
        let mut isect = gt_gap_sets[0].clone();
        for gs in &gt_gap_sets[1..] {
            isect = isect.intersection(gs).cloned().collect();
        }
        isect
    } else if gt_gap_sets.len() == 1 {
        gt_gap_sets[0].clone()
    } else {
        HashSet::new()
    };

    // Build per-sample GT SNP maps: pos -> alt allele
    let mut gt_snp_maps: Vec<HashMap<usize, u8>> = Vec::new();
    for sample_id in sample_ids {
        let mut snp_map = HashMap::new();
        if let Some(sample_data) = data.get(sample_id) {
            if let Some(gt_data) = sample_data.get(gt_pipeline_id) {
                for snp in &gt_data.snps {
                    let pos = snp.pos - 1; // VCF is 1-based, convert to 0-based
                    let alt = snp.alt.as_bytes().first().copied().unwrap_or(b'N');
                    // Filter out bogus SNPs where alt == reference
                    if pos < ref_bytes.len() && alt != ref_bytes[pos] {
                        snp_map.insert(pos, alt);
                    }
                }
            }
        }
        gt_snp_maps.push(snp_map);
    }

    // All GT SNP positions
    let mut all_gt_snp_positions: HashSet<usize> = HashSet::new();
    for snp_map in &gt_snp_maps {
        all_gt_snp_positions.extend(snp_map.keys());
    }

    // --- Gap-Union GT discriminating positions ---
    let mut gt_disc_union: Vec<usize> = Vec::new();
    for &pos in &all_gt_snp_positions {
        if gt_gap_union.contains(&pos) { continue; }
        let alleles: Vec<Option<u8>> = (0..n).map(|idx| gt_snp_maps[idx].get(&pos).copied()).collect();
        let first = alleles[0];
        if alleles.iter().any(|a| *a != first) {
            gt_disc_union.push(pos);
        }
    }

    // --- Gap-Intersect GT discriminating positions ---
    let mut gt_disc_intersect: Vec<usize> = Vec::new();
    for &pos in &all_gt_snp_positions {
        if gt_gap_intersection.contains(&pos) { continue; }
        let alleles: Vec<Option<u8>> = (0..n)
            .filter(|&idx| !gt_gap_sets[idx].contains(&pos))
            .map(|idx| gt_snp_maps[idx].get(&pos).copied())
            .collect();
        if alleles.len() < 2 { continue; }
        let first = alleles[0];
        if alleles.iter().any(|a| *a != first) {
            gt_disc_intersect.push(pos);
        }
    }

    // For quick lookup: convert core SNP positions to HashSet per pipeline
    let mut results = Vec::new();

    for (pipeline_id, core_data) in core_snp_data {
        if pipeline_id == gt_pipeline_id { continue; }

        let core_positions: HashSet<usize> = core_data.positions.iter().map(|p| p.pos).collect();

        // Build lookup: pos -> CoreSnpPosition
        let core_pos_map: HashMap<usize, &CoreSnpPosition> = core_data.positions.iter()
            .map(|p| (p.pos, p))
            .collect();

        // --- Helper: check allele concordance at a position ---
        let check_concordance = |pos: usize, active_indices: &[usize]| -> bool {
            if !core_data.has_alleles { return false; }
            let core_pos = match core_pos_map.get(&pos) {
                Some(p) => p,
                None => return false,
            };
            active_indices.iter().all(|&idx| {
                let sample_id = &sample_ids[idx];
                let gt_allele = gt_snp_maps[idx].get(&pos).copied()
                    .unwrap_or_else(|| if pos < ref_bytes.len() { ref_bytes[pos] } else { b'N' });
                let pl_allele = core_pos.alleles.get(sample_id)
                    .and_then(|a| a.as_bytes().first().copied())
                    .unwrap_or_else(|| if pos < ref_bytes.len() { ref_bytes[pos] } else { b'N' });
                gt_allele == pl_allele
            })
        };

        // --- Gap-Union: same_pos, concordant ---
        let mut gu_same_pos: u32 = 0;
        let mut gu_concordant: u32 = 0;
        for &pos in &gt_disc_union {
            if !core_positions.contains(&pos) { continue; }
            gu_same_pos += 1;
            if core_data.has_alleles {
                let all_indices: Vec<usize> = (0..n).collect();
                if check_concordance(pos, &all_indices) {
                    gu_concordant += 1;
                }
            }
        }
        let gu_pl_in_gaps: u32 = core_positions.iter()
            .filter(|p| gt_gap_union.contains(p))
            .count() as u32;

        // --- Gap-Intersect: same_pos, concordant ---
        let mut gi_same_pos: u32 = 0;
        let mut gi_concordant: u32 = 0;
        for &pos in &gt_disc_intersect {
            if !core_positions.contains(&pos) { continue; }
            let active: Vec<usize> = (0..n).filter(|&idx| !gt_gap_sets[idx].contains(&pos)).collect();
            if active.len() < 2 { continue; }
            gi_same_pos += 1;
            if core_data.has_alleles && check_concordance(pos, &active) {
                gi_concordant += 1;
            }
        }
        let gi_pl_in_gaps: u32 = core_positions.iter()
            .filter(|p| gt_gap_intersection.contains(p))
            .count() as u32;

        // --- Pairwise ---
        let mut pw_total_disc: f64 = 0.0;
        let mut pw_total_same_pos: f64 = 0.0;
        let mut pw_total_concordant: f64 = 0.0;
        let mut num_pairs: u32 = 0;

        for i in 0..n {
            for j in (i + 1)..n {
                let pair_gap_union: HashSet<usize> = gt_gap_sets[i].union(&gt_gap_sets[j]).cloned().collect();

                let mut pair_disc: Vec<usize> = Vec::new();
                for &pos in &all_gt_snp_positions {
                    if pair_gap_union.contains(&pos) { continue; }
                    let a_i = gt_snp_maps[i].get(&pos).copied();
                    let a_j = gt_snp_maps[j].get(&pos).copied();
                    if a_i != a_j {
                        pair_disc.push(pos);
                    }
                }

                let mut p_same: u32 = 0;
                let mut p_conc: u32 = 0;
                let pair_indices = vec![i, j];
                for &pos in &pair_disc {
                    if !core_positions.contains(&pos) { continue; }
                    p_same += 1;
                    if core_data.has_alleles && check_concordance(pos, &pair_indices) {
                        p_conc += 1;
                    }
                }

                pw_total_disc += pair_disc.len() as f64;
                pw_total_same_pos += p_same as f64;
                pw_total_concordant += p_conc as f64;
                num_pairs += 1;
            }
        }

        let np = if num_pairs > 0 { num_pairs as f64 } else { 1.0 };

        results.push(GtDiscVsPipelineResult {
            pipeline_id: pipeline_id.clone(),
            pl_total_core_snps: core_data.positions.len() as u32,
            gap_intersect: GapStrategyResult {
                gt_disc: gt_disc_intersect.len() as u32,
                same_pos: gi_same_pos,
                concordant: if core_data.has_alleles { Some(gi_concordant) } else { None },
                pl_snps_in_gt_gaps: gi_pl_in_gaps,
            },
            gap_union: GapStrategyResult {
                gt_disc: gt_disc_union.len() as u32,
                same_pos: gu_same_pos,
                concordant: if core_data.has_alleles { Some(gu_concordant) } else { None },
                pl_snps_in_gt_gaps: gu_pl_in_gaps,
            },
            pairwise: PairwiseGtDiscResult {
                gt_disc_avg: pw_total_disc / np,
                same_pos_avg: pw_total_same_pos / np,
                concordant_avg: if core_data.has_alleles { Some(pw_total_concordant / np) } else { None },
                num_pairs,
            },
        });
    }

    results
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
                let has_vcf_files = pipeline.samples.values().any(|f| f.vcf.is_some());
                let has_bam = pipeline.samples.values().any(|f| f.bam.is_some());
                // Ground truth with BAM will have SNPs from pileup, so treat as having VCF for comparisons
                let has_vcf = has_vcf_files || (pipeline.ground_truth && has_bam);
                // Track if SNPs come from BAM pileup (no variant calling)
                let from_bam_pileup = pipeline.ground_truth && !has_vcf_files && has_bam;
                pipelines.insert(
                    pipeline_id.clone(),
                    PipelineInfo {
                        label: pipeline.label.clone(),
                        command: pipeline.command.clone(),
                        has_vcf,
                        has_bam,
                        ground_truth: pipeline.ground_truth,
                        from_bam_pileup,
                    },
                );
            }
        }

        // Process data for each sample and pipeline
        let mut data: HashMap<String, HashMap<String, PipelineData>> = HashMap::new();
        let mut total_mnps_found = 0;
        let mut total_snps_from_mnps = 0;
        // Per-pipeline MNP stats
        let mut pipeline_mnp_stats: HashMap<String, MnpStats> = HashMap::new();

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
                            pipeline_data.bam_path = Some(bam_path.clone());
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
                            pipeline_data.vcf_path = Some(vcf_path.clone());
                            total_mnps_found += result.mnps_found;
                            total_snps_from_mnps += result.snps_from_mnps;
                            // Track per-pipeline MNP stats
                            let entry = pipeline_mnp_stats.entry(pipeline_id.clone()).or_insert(MnpStats {
                                mnps_found: 0,
                                snps_from_mnps: 0,
                            });
                            entry.mnps_found += result.mnps_found;
                            entry.snps_from_mnps += result.snps_from_mnps;
                            log::info!("  Found {} SNPs", pipeline_data.snps.len());
                        } else if pipeline.ground_truth && files.bam.is_some() {
                            // For ground truth without VCF, load SNPs from BAM pileup
                            let bam_path = files.bam.as_ref().unwrap();
                            log::info!(
                                "Loading SNPs for {}/{} from BAM pileup {}",
                                sample_id,
                                pipeline_id,
                                bam_path
                            );
                            let pileup_snps = pileup::get_snps_from_pileup(
                                Path::new(bam_path),
                                &ref_seq,
                                &ref_name,
                                config.options.min_depth as u32,
                                config.options.min_consensus,
                            )?;
                            pipeline_data.snps = pileup_snps.iter().map(|ps| Snp {
                                pos: ps.pos,
                                ref_allele: ps.ref_base.to_string(),
                                alt: ps.alt_base.to_string(),
                                qual: 0.0,  // No quality score from pileup
                                dp: ps.depth as usize,
                            }).collect();
                            log::info!("  Found {} SNPs from pileup", pipeline_data.snps.len());
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

        // Helper function to check if position is in any gap
        fn pos_in_gaps(pos: usize, gaps: &[[usize; 2]]) -> bool {
            gaps.iter().any(|[start, end]| pos >= *start && pos < *end)
        }

        // Calculate SNPs in gaps for ALL pipeline pairs (gap_pipeline -> snp_pipeline -> stats)
        // This replaces the old snps_in_gt_gaps which only worked for ground truth
        let snps_in_gaps: Option<HashMap<String, HashMap<String, SnpsInGapsStats>>> = {
            let mut all_stats: HashMap<String, HashMap<String, SnpsInGapsStats>> = HashMap::new();

            // For each pipeline that has BAM data (source of gaps)
            for gap_pipeline_id in &pipeline_ids {
                if let Some(gap_pipeline_info) = pipelines.get(gap_pipeline_id) {
                    if !gap_pipeline_info.has_bam {
                        continue;  // No BAM = no gap information
                    }
                }

                let mut gap_pipeline_stats: HashMap<String, SnpsInGapsStats> = HashMap::new();

                // For each OTHER pipeline that has SNPs
                for snp_pipeline_id in &pipeline_ids {
                    if snp_pipeline_id == gap_pipeline_id {
                        continue;  // Don't compare pipeline with itself
                    }
                    if let Some(snp_pipeline_info) = pipelines.get(snp_pipeline_id) {
                        if !snp_pipeline_info.has_vcf {
                            continue;  // No VCF = no SNP information
                        }
                    }

                    let mut total_snps = 0usize;
                    let mut snps_in_gap = 0usize;

                    for sample_id in &sample_ids {
                        // Get gaps from gap_pipeline for this sample
                        let gaps = data
                            .get(sample_id)
                            .and_then(|s| s.get(gap_pipeline_id))
                            .map(|d| d.gaps.as_slice())
                            .unwrap_or(&[]);

                        // Get SNPs from snp_pipeline for this sample
                        if let Some(sample_data) = data.get(sample_id) {
                            if let Some(pipeline_data) = sample_data.get(snp_pipeline_id) {
                                for snp in &pipeline_data.snps {
                                    total_snps += 1;
                                    if pos_in_gaps(snp.pos, gaps) {
                                        snps_in_gap += 1;
                                    }
                                }
                            }
                        }
                    }

                    let percentage = if total_snps > 0 {
                        (snps_in_gap as f64 / total_snps as f64) * 100.0
                    } else {
                        0.0
                    };

                    gap_pipeline_stats.insert(snp_pipeline_id.clone(), SnpsInGapsStats {
                        total_snps,
                        snps_in_gaps: snps_in_gap,
                        percentage,
                    });
                }

                if !gap_pipeline_stats.is_empty() {
                    all_stats.insert(gap_pipeline_id.clone(), gap_pipeline_stats);
                }
            }

            if !all_stats.is_empty() {
                Some(all_stats)
            } else {
                None
            }
        };

        // Calculate SNPs in ground truth gaps (DEPRECATED - kept for backwards compatibility)
        let snps_in_gt_gaps = if let Some(gt_pipeline) = config.ground_truth_pipeline() {
            // Use the new snps_in_gaps data if available
            snps_in_gaps.as_ref().and_then(|all| all.get(&gt_pipeline).cloned())
        } else {
            None
        };

        // Calculate ground truth pileup SNPs (raw count from BAM without variant calling)
        let ground_truth_pileup = if let Some(gt_pipeline) = config.ground_truth_pipeline() {
            log::info!("Calculating ground truth pileup SNP count...");

            let mut per_sample: HashMap<String, usize> = HashMap::new();
            let mut total_snps = 0usize;
            let mut total_covered = 0usize;

            // For each sample, count SNPs from the ground truth BAM
            for sample_id in &sample_ids {
                let bam_path: Option<String> = config.pipelines.get(&gt_pipeline)
                    .and_then(|p| p.samples.get(sample_id))
                    .and_then(|f| f.bam.clone());

                if let Some(ref path) = bam_path {
                    match pileup::count_snps_from_pileup(
                        Path::new(path),
                        &ref_seq,
                        &ref_name,
                        config.options.min_depth as u32,
                        config.options.min_consensus,
                    ) {
                        Ok((snp_count, covered_positions)) => {
                            log::info!("  {}: {} SNPs from pileup ({} covered positions)",
                                sample_id, snp_count, covered_positions);
                            per_sample.insert(sample_id.clone(), snp_count);
                            total_snps += snp_count;
                            total_covered = total_covered.max(covered_positions);
                        }
                        Err(e) => {
                            log::warn!("Failed to count pileup SNPs for {}: {}", sample_id, e);
                        }
                    }
                }
            }

            // Compare with VCF pipelines
            let mut pipeline_comparison: HashMap<String, PipelineVsGroundTruth> = HashMap::new();

            for pipeline_id in &pipeline_ids {
                if pipeline_id == &gt_pipeline {
                    continue;
                }
                if let Some(pipeline_info) = pipelines.get(pipeline_id) {
                    if !pipeline_info.has_vcf {
                        continue;
                    }
                }

                // Count total SNPs for this pipeline across all samples
                let mut pipeline_snps = 0usize;
                for sample_id in &sample_ids {
                    if let Some(sample_data) = data.get(sample_id) {
                        if let Some(pipeline_data) = sample_data.get(pipeline_id) {
                            pipeline_snps += pipeline_data.snps.len();
                        }
                    }
                }

                let difference = pipeline_snps as i64 - total_snps as i64;
                let percentage_diff = if total_snps > 0 {
                    (difference as f64 / total_snps as f64) * 100.0
                } else {
                    0.0
                };

                pipeline_comparison.insert(pipeline_id.clone(), PipelineVsGroundTruth {
                    pipeline_snps,
                    ground_truth_snps: total_snps,
                    difference,
                    percentage_diff,
                });

                log::info!("  {} vs GT: {} vs {} SNPs (diff: {:+}, {:+.2}%)",
                    pipeline_id, pipeline_snps, total_snps, difference, percentage_diff);
            }

            if !per_sample.is_empty() {
                Some(GroundTruthPileupStats {
                    total_snps,
                    per_sample,
                    covered_positions: total_covered,
                    pipeline_comparison,
                })
            } else {
                None
            }
        } else {
            None
        };

        // Build summary
        // Only include MNP stats if any MNPs were found
        let mnp_stats = if pipeline_mnp_stats.values().any(|s| s.mnps_found > 0) {
            Some(pipeline_mnp_stats)
        } else {
            None
        };

        let summary = Summary {
            total_samples: sample_ids.len(),
            total_pipelines: pipeline_ids.len(),
            generated_at: chrono::Utc::now().to_rfc3339(),
            coreguard_version: env!("CARGO_PKG_VERSION").to_string(),
            warnings,
            snps_in_gt_gaps,
            snps_in_gaps,
            ground_truth_pileup,
            mnp_stats,
        };

        // Get description content (from file or inline)
        let description = config.get_description_content();

        // Load pre-computed distance matrices from pipelines
        let mut pipeline_distance_matrices: HashMap<String, PipelineDistanceMatrix> = HashMap::new();
        for (pipeline_id, pipeline_config) in &config.pipelines {
            if let Some(matrix_path) = &pipeline_config.distance_matrix {
                match parse_distance_matrix_tsv(matrix_path) {
                    Ok(matrix) => {
                        log::info!("Loaded distance matrix for pipeline '{}' from {}", pipeline_id, matrix_path);
                        pipeline_distance_matrices.insert(pipeline_id.clone(), matrix);
                    }
                    Err(e) => {
                        log::warn!("Failed to load distance matrix for pipeline '{}': {}", pipeline_id, e);
                    }
                }
            }
        }

        // Compute GT disc vs pipelines using core SNP files
        let gt_disc_vs_pipelines = if let Some(gt_pipeline_id) = config.ground_truth_pipeline() {
            let mut core_snp_data: HashMap<String, CoreSnpData> = HashMap::new();
            for (pipeline_id, pipeline_config) in &config.pipelines {
                if let Some(core_snps_path) = &pipeline_config.core_snps {
                    match parse_core_snps(core_snps_path) {
                        Ok(core_data) => {
                            log::info!(
                                "Loaded core SNPs for pipeline '{}': {} positions (has_alleles: {})",
                                pipeline_id, core_data.positions.len(), core_data.has_alleles
                            );
                            core_snp_data.insert(pipeline_id.clone(), core_data);
                        }
                        Err(e) => {
                            log::warn!("Failed to load core SNPs for pipeline '{}': {}", pipeline_id, e);
                        }
                    }
                }
            }

            if !core_snp_data.is_empty() {
                let results = compute_gt_disc_vs_pipelines(
                    &data,
                    &sample_ids,
                    &gt_pipeline_id,
                    &core_snp_data,
                    &ref_seq,
                );
                if !results.is_empty() {
                    log::info!("Computed GT disc vs pipelines for {} pipeline(s)", results.len());
                    Some(results)
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
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
            description,
            pipeline_distance_matrices,
            gt_disc_vs_pipelines,
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

    /// Save report to binary format (bincode) - ~10x faster to parse
    pub fn save_binary<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let binary = bincode::serialize(self)?;
        std::fs::write(path.as_ref(), binary)?;
        Ok(())
    }

    /// Save report to gzipped binary format (bincode)
    pub fn save_binary_gzip<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;

        let file = std::fs::File::create(path.as_ref())?;
        let mut encoder = GzEncoder::new(file, Compression::default());

        let binary = bincode::serialize(self)?;
        encoder.write_all(&binary)?;
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

/// Parse TSV distance matrix file
/// Format: header row with sample names, then matrix rows
/// Example:
///     sample1 sample2 sample3
/// sample1 0   5   10
/// sample2 5   0   8
/// sample3 10  8   0
fn parse_distance_matrix_tsv(path: &str) -> Result<PipelineDistanceMatrix> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("Failed to read distance matrix file: {}", path))?;

    let mut lines = content.lines().peekable();

    // Parse header to get sample order
    let header = lines.next()
        .ok_or_else(|| anyhow::anyhow!("Empty distance matrix file"))?;

    let samples: Vec<String> = header
        .split('\t')
        .skip(1) // Skip first column (row labels)
        .filter(|s| !s.is_empty())
        .map(|s| s.trim().to_string())
        .collect();

    if samples.is_empty() {
        anyhow::bail!("No sample names found in distance matrix header");
    }

    // Parse matrix rows
    let mut matrix: Vec<Vec<i64>> = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue;
        }

        // Skip row label, parse distances
        let row: Vec<i64> = parts[1..]
            .iter()
            .filter(|s| !s.is_empty())
            .map(|s| s.trim().parse::<i64>().unwrap_or(0))
            .collect();

        if row.len() == samples.len() {
            matrix.push(row);
        }
    }

    if matrix.len() != samples.len() {
        anyhow::bail!(
            "Matrix dimensions mismatch: {} samples but {} rows",
            samples.len(),
            matrix.len()
        );
    }

    Ok(PipelineDistanceMatrix { samples, matrix })
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
) -> Result<HashMap<String, HashMap<String, PolymorphicSite>>> {
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
    let mut result: HashMap<String, HashMap<String, PolymorphicSite>> = HashMap::new();

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
                    config.options.min_consensus,
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
        let mut pipeline_sites: HashMap<String, PolymorphicSite> = HashMap::new();
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
                pipeline_sites.insert(pos.to_string(), PolymorphicSite {
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
