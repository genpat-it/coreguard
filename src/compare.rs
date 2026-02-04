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
use crate::parsers::{self, CoreSnpData, CoreSnpPosition};
use crate::pileup;

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
    #[serde(default)]
    pub data: HashMap<String, HashMap<String, PipelineData>>,

    /// Pre-computed per-pipeline KPIs (replaces WASM get_kpis)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kpis: Option<PreComputedKpis>,

    /// Pre-computed per-pipeline statistics (replaces WASM get_global_stats_for_pipeline / get_pairwise_usable_stats_for_pipeline)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pipeline_stats: Option<HashMap<String, PipelineStats>>,

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

    /// Raw YAML configuration used to generate this report
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub config_yaml: Option<String>,
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
    /// Pileup options used for GT computation
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pileup_options: Option<PileupOptions>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PileupOptions {
    pub min_depth: usize,
    pub min_qual: f64,
    pub min_consensus: f64,
    pub include_indels: bool,
}

/// Pre-computed KPIs for all pipelines (replaces WASM get_kpis)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PreComputedKpis {
    /// Per-pipeline basic KPIs
    pub pipelines: HashMap<String, PipelineKpi>,
}

/// Basic KPIs for a single pipeline
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineKpi {
    /// Number of SNP positions (any sample has SNP)
    pub total_snps: u32,
    /// Number of gap regions across all samples
    pub total_gap_positions: u32,
    /// Core SNPs: positions where at least 2 samples differ
    pub core_snps: u32,
}

/// Pre-computed per-pipeline statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineStats {
    /// Global stats (gap-intersect and gap-union strategies)
    pub global: GlobalPipelineStats,
    /// Pairwise stats
    pub pairwise: PairwisePipelineStats,
}

/// Global statistics for a pipeline
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GlobalPipelineStats {
    pub ref_length: u32,
    /// Gap-intersect: positions where ALL samples have gap
    pub gap_intersect_count: u32,
    /// Gap-union: positions where ANY sample has gap
    pub gap_union_count: u32,
    /// Usable space (ref_length - gap_union)
    pub usable_intersect: u32,
    pub usable_union: u32,
    /// Total SNPs in usable space
    pub total_snps_intersect: u32,
    pub total_snps_union: u32,
    /// Consensus SNPs (all samples same alt)
    pub consensus_snps_intersect: u32,
    pub consensus_snps_union: u32,
    /// Discriminating SNPs (samples differ)
    pub disc_snps_intersect: u32,
    pub disc_snps_union: u32,
    /// Missing VCF calls (some samples have call, others don't) in usable space
    pub missing_vcf_intersect: u32,
    pub missing_vcf_union: u32,
    /// Discriminating SNP breakdown (VCF pipelines only)
    pub disc_breakdown: Option<DiscBreakdown>,
}

/// Discriminating SNP breakdown for VCF pipelines
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiscBreakdown {
    /// Gap-affected: at least one sample skipped due to gap
    pub gap_affected: u32,
    /// GT-consensus: GT pileup shows all samples agree
    pub gt_consensus: u32,
    /// Majority-rule: all but one sample agree
    pub majority_rule: u32,
    /// Confirmed: genuine disagreement
    pub confirmed: u32,
}

/// Pairwise statistics for a pipeline
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairwisePipelineStats {
    pub num_pairs: u32,
    /// Min/median/max/avg discriminating SNPs across pairs (gap-union)
    pub disc_snps_min: u32,
    pub disc_snps_median: f64,
    pub disc_snps_max: u32,
    pub disc_snps_avg: f64,
    /// Average usable space across pairs
    pub usable_space_avg: f64,
    /// Per-sample averages: sample_id -> (avg_usable_space, avg_disc_snps)
    pub per_sample: HashMap<String, SamplePairwiseStats>,
}

/// Per-sample pairwise stats (individual pair values + averages)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SamplePairwiseStats {
    pub avg_usable_space: f64,
    pub avg_disc_snps: f64,
    /// Per-pair breakdown: other_sample_id -> (usable_space, disc_snps)
    #[serde(default)]
    pub pairs: HashMap<String, PairDetail>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairDetail {
    pub usable_space: u32,
    pub disc_snps: u32,
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

/// Result of GT disc vs pipeline comparison
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GtDiscVsPipelineResult {
    pub pipeline_id: String,
    /// Total core SNP positions in the pipeline's core_snps file
    pub pl_total_core_snps: u32,
    /// Discriminating core SNPs: positions where at least 2 samples differ (contribute to Hamming distance).
    /// Exact for Snippy (allele comparison), conservative underestimate for CFSAN (subset check).
    pub pl_discriminating_core_snps: u32,
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
    /// Per-position details for drill-down (optional, may be large)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub position_details: Option<Vec<PositionDetail>>,
}

/// Detail for a single position in the GT disc vs pipeline comparison
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PositionDetail {
    pub pos: usize,
    #[serde(rename = "ref")]
    pub ref_allele: String,
    /// Sample -> GT allele ("-" if in gap)
    pub gt: HashMap<String, String>,
    /// Sample -> Pipeline allele ("-" if in gap)
    pub pl: HashMap<String, String>,
    /// concordant, discordant, not_disc_in_pl, lost, in_gt_gap, in_pl_gap
    pub status: String,
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

        // Build lookup: pos -> CoreSnpPosition
        let core_pos_map: HashMap<usize, &CoreSnpPosition> = core_data.positions.iter()
            .map(|p| (p.pos, p))
            .collect();

        // Build per-sample pipeline gap sets (from VCF/BAM data)
        let pl_gap_sets: Vec<HashSet<usize>> = sample_ids.iter().map(|sid| {
            let mut gaps = HashSet::new();
            if let Some(sample_data) = data.get(sid) {
                if let Some(pl_data) = sample_data.get(pipeline_id.as_str()) {
                    for gap in &pl_data.gaps {
                        for pos in gap[0]..gap[1] { gaps.insert(pos); }
                    }
                }
            }
            gaps
        }).collect();

        // Helper: check if a core SNP position is discriminating
        let is_discriminating = |p: &CoreSnpPosition| -> bool {
            if core_data.has_alleles {
                let valid: Vec<&str> = p.alleles.values()
                    .map(|a| a.as_str())
                    .filter(|a| *a != "-" && *a != "N")
                    .collect();
                if valid.len() < 2 { return false; }
                let first = valid[0];
                valid.iter().any(|a| *a != first)
            } else {
                let active_samples: Vec<usize> = (0..n)
                    .filter(|&idx| !pl_gap_sets[idx].contains(&p.pos))
                    .collect();
                if active_samples.len() < 2 { return false; }
                let snp_set: HashSet<&str> = p.samples_with_snp.iter().map(|s| s.as_str()).collect();
                let has_alt = active_samples.iter().any(|&idx| snp_set.contains(sample_ids[idx].as_str()));
                let has_ref = active_samples.iter().any(|&idx| !snp_set.contains(sample_ids[idx].as_str()));
                has_alt && has_ref
            }
        };

        // Only use discriminating positions for comparison with GT
        let disc_positions: HashSet<usize> = core_data.positions.iter()
            .filter(|p| is_discriminating(p))
            .map(|p| p.pos)
            .collect();
        let pl_discriminating_core_snps = disc_positions.len() as u32;

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

        // --- Helper: build a PositionDetail for a GT disc position ---
        let build_detail = |pos: usize, status: &str, active_indices: &[usize]| -> PositionDetail {
            let ref_allele = if pos < ref_bytes.len() { (ref_bytes[pos] as char).to_string() } else { "N".to_string() };
            let mut gt = HashMap::new();
            let mut pl = HashMap::new();
            for &idx in active_indices {
                let sid = sample_ids[idx].clone();
                // GT allele
                let gt_a = if gt_gap_sets[idx].contains(&pos) {
                    "-".to_string()
                } else {
                    let a = gt_snp_maps[idx].get(&pos).copied()
                        .unwrap_or_else(|| if pos < ref_bytes.len() { ref_bytes[pos] } else { b'N' });
                    (a as char).to_string()
                };
                // PL allele: use core_snps alleles if available, else gap/ref
                let pl_a = if pl_gap_sets[idx].contains(&pos) {
                    "-".to_string()
                } else if core_data.has_alleles {
                    core_pos_map.get(&pos)
                        .and_then(|cp| cp.alleles.get(&sid))
                        .map(|a| if a == "-" || a == "N" { ref_allele.clone() } else { a.chars().next().unwrap_or('?').to_string() })
                        .unwrap_or_else(|| ref_allele.clone())
                } else {
                    "?".to_string()
                };
                gt.insert(sid.clone(), gt_a);
                pl.insert(sid, pl_a);
            }
            PositionDetail { pos, ref_allele, gt, pl, status: status.to_string() }
        };

        let all_indices: Vec<usize> = (0..n).collect();

        // --- Gap-Union: same_pos, concordant, position details ---
        let mut gu_same_pos: u32 = 0;
        let mut gu_concordant: u32 = 0;
        let mut gu_details: Vec<PositionDetail> = Vec::new();
        for &pos in &gt_disc_union {
            if !disc_positions.contains(&pos) {
                // GT disc pos not discriminating in pipeline
                let has_any_pl_snp = core_pos_map.contains_key(&pos);
                if has_any_pl_snp {
                    gu_details.push(build_detail(pos, "not_disc_in_pl", &all_indices));
                } else if pl_gap_sets.iter().any(|gs| gs.contains(&pos)) {
                    gu_details.push(build_detail(pos, "in_pl_gap", &all_indices));
                } else {
                    gu_details.push(build_detail(pos, "lost", &all_indices));
                }
                continue;
            }
            gu_same_pos += 1;
            if core_data.has_alleles {
                if check_concordance(pos, &all_indices) {
                    gu_concordant += 1;
                    gu_details.push(build_detail(pos, "concordant", &all_indices));
                } else {
                    gu_details.push(build_detail(pos, "discordant", &all_indices));
                }
            } else {
                gu_details.push(build_detail(pos, "concordant", &all_indices));
            }
        }
        let gu_pl_in_gaps: u32 = disc_positions.iter()
            .filter(|p| gt_gap_union.contains(p))
            .count() as u32;
        // Add in_gt_gap positions
        for &pos in disc_positions.iter().filter(|p| gt_gap_union.contains(p)) {
            gu_details.push(build_detail(pos, "in_gt_gap", &all_indices));
        }
        // Add pl_only positions (pipeline disc, not in GT disc, not in GT gap = table "Lost")
        for &pos in disc_positions.iter().filter(|p| !gt_disc_union.contains(p) && !gt_gap_union.contains(p)) {
            gu_details.push(build_detail(pos, "pl_only", &all_indices));
        }
        gu_details.sort_by_key(|d| d.pos);

        // --- Gap-Intersect: same_pos, concordant, position details ---
        let mut gi_same_pos: u32 = 0;
        let mut gi_concordant: u32 = 0;
        let mut gi_details: Vec<PositionDetail> = Vec::new();
        for &pos in &gt_disc_intersect {
            let active: Vec<usize> = (0..n).filter(|&idx| !gt_gap_sets[idx].contains(&pos)).collect();
            if active.len() < 2 { continue; }
            if !disc_positions.contains(&pos) {
                let has_any_pl_snp = core_pos_map.contains_key(&pos);
                if has_any_pl_snp {
                    gi_details.push(build_detail(pos, "not_disc_in_pl", &all_indices));
                } else if pl_gap_sets.iter().any(|gs| gs.contains(&pos)) {
                    gi_details.push(build_detail(pos, "in_pl_gap", &all_indices));
                } else {
                    gi_details.push(build_detail(pos, "lost", &all_indices));
                }
                continue;
            }
            gi_same_pos += 1;
            if core_data.has_alleles && check_concordance(pos, &active) {
                gi_concordant += 1;
                gi_details.push(build_detail(pos, "concordant", &all_indices));
            } else if core_data.has_alleles {
                gi_details.push(build_detail(pos, "discordant", &all_indices));
            } else {
                gi_details.push(build_detail(pos, "concordant", &all_indices));
            }
        }
        let gi_pl_in_gaps: u32 = disc_positions.iter()
            .filter(|p| gt_gap_intersection.contains(p))
            .count() as u32;
        // Add in_gt_gap positions for gap-intersect
        for &pos in disc_positions.iter().filter(|p| gt_gap_intersection.contains(p)) {
            gi_details.push(build_detail(pos, "in_gt_gap", &all_indices));
        }
        // Add pl_only positions for gap-intersect (pipeline disc, not in GT disc, not in GT gap = table "Lost")
        for &pos in disc_positions.iter().filter(|p| !gt_disc_intersect.contains(p) && !gt_gap_intersection.contains(p)) {
            gi_details.push(build_detail(pos, "pl_only", &all_indices));
        }
        gi_details.sort_by_key(|d| d.pos);

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
                    if !disc_positions.contains(&pos) { continue; }
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
            pl_discriminating_core_snps,
            gap_intersect: GapStrategyResult {
                gt_disc: gt_disc_intersect.len() as u32,
                same_pos: gi_same_pos,
                concordant: if core_data.has_alleles { Some(gi_concordant) } else { None },
                pl_snps_in_gt_gaps: gi_pl_in_gaps,
                position_details: Some(gi_details),
            },
            gap_union: GapStrategyResult {
                gt_disc: gt_disc_union.len() as u32,
                same_pos: gu_same_pos,
                concordant: if core_data.has_alleles { Some(gu_concordant) } else { None },
                pl_snps_in_gt_gaps: gu_pl_in_gaps,
                position_details: Some(gu_details),
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

/// Compute GT disc vs pipeline metrics using VCF data (no core_snps needed).
///
/// Each non-GT pipeline's own VCF-derived SNPs and gaps are used for comparison.
/// Concordance is always available since we have per-sample alleles from VCFs.
pub fn compute_gt_disc_vs_pipelines_from_vcf(
    data: &HashMap<String, HashMap<String, PipelineData>>,
    sample_ids: &[String],
    pipeline_ids: &[String],
    gt_pipeline_id: &str,
    ref_seq: &str,
) -> Vec<GtDiscVsPipelineResult> {
    use std::collections::HashSet;

    let n = sample_ids.len();
    if n < 2 {
        return Vec::new();
    }

    let ref_bytes = ref_seq.as_bytes();

    // Build per-sample GT gap sets
    let gt_gap_sets: Vec<HashSet<usize>> = sample_ids.iter().map(|sid| {
        let mut gaps = HashSet::new();
        if let Some(sample_data) = data.get(sid) {
            if let Some(gt_data) = sample_data.get(gt_pipeline_id) {
                for gap in &gt_data.gaps {
                    for pos in gap[0]..gap[1] {
                        gaps.insert(pos);
                    }
                }
            }
        }
        gaps
    }).collect();

    // GT gap union / intersection
    let mut gt_gap_union: HashSet<usize> = HashSet::new();
    for gs in &gt_gap_sets { gt_gap_union.extend(gs); }

    let gt_gap_intersection = if gt_gap_sets.len() > 1 {
        let mut isect = gt_gap_sets[0].clone();
        for gs in &gt_gap_sets[1..] { isect = isect.intersection(gs).cloned().collect(); }
        isect
    } else if gt_gap_sets.len() == 1 {
        gt_gap_sets[0].clone()
    } else {
        HashSet::new()
    };

    // Build per-sample GT SNP maps: pos -> alt allele (0-based)
    let gt_snp_maps: Vec<HashMap<usize, u8>> = sample_ids.iter().map(|sid| {
        let mut snp_map = HashMap::new();
        if let Some(sample_data) = data.get(sid) {
            if let Some(gt_data) = sample_data.get(gt_pipeline_id) {
                for snp in &gt_data.snps {
                    let pos = snp.pos - 1;
                    let alt = snp.alt.as_bytes().first().copied().unwrap_or(b'N');
                    if pos < ref_bytes.len() && alt != ref_bytes[pos] {
                        snp_map.insert(pos, alt);
                    }
                }
            }
        }
        snp_map
    }).collect();

    // All GT SNP positions
    let mut all_gt_snp_positions: HashSet<usize> = HashSet::new();
    for snp_map in &gt_snp_maps { all_gt_snp_positions.extend(snp_map.keys()); }

    // Gap-Union GT discriminating positions
    let gt_disc_union: Vec<usize> = all_gt_snp_positions.iter().filter(|&&pos| {
        if gt_gap_union.contains(&pos) { return false; }
        let alleles: Vec<Option<u8>> = (0..n).map(|idx| gt_snp_maps[idx].get(&pos).copied()).collect();
        let first = alleles[0];
        alleles.iter().any(|a| *a != first)
    }).copied().collect();

    // Gap-Intersect GT discriminating positions
    let gt_disc_intersect: Vec<usize> = all_gt_snp_positions.iter().filter(|&&pos| {
        if gt_gap_intersection.contains(&pos) { return false; }
        let alleles: Vec<Option<u8>> = (0..n)
            .filter(|&idx| !gt_gap_sets[idx].contains(&pos))
            .map(|idx| gt_snp_maps[idx].get(&pos).copied())
            .collect();
        if alleles.len() < 2 { return false; }
        let first = alleles[0];
        alleles.iter().any(|a| *a != first)
    }).copied().collect();

    let mut results = Vec::new();

    for pipeline_id in pipeline_ids {
        if pipeline_id == gt_pipeline_id { continue; }

        // Build per-sample pipeline gap sets and SNP maps
        let pl_gap_sets: Vec<HashSet<usize>> = sample_ids.iter().map(|sid| {
            let mut gaps = HashSet::new();
            if let Some(sample_data) = data.get(sid) {
                if let Some(pl_data) = sample_data.get(pipeline_id.as_str()) {
                    for gap in &pl_data.gaps {
                        for pos in gap[0]..gap[1] { gaps.insert(pos); }
                    }
                }
            }
            gaps
        }).collect();

        let pl_snp_maps: Vec<HashMap<usize, u8>> = sample_ids.iter().map(|sid| {
            let mut snp_map = HashMap::new();
            if let Some(sample_data) = data.get(sid) {
                if let Some(pl_data) = sample_data.get(pipeline_id.as_str()) {
                    for snp in &pl_data.snps {
                        let pos = snp.pos - 1;
                        let alt = snp.alt.as_bytes().first().copied().unwrap_or(b'N');
                        if pos < ref_bytes.len() && alt != ref_bytes[pos] {
                            snp_map.insert(pos, alt);
                        }
                    }
                }
            }
            snp_map
        }).collect();

        let mut all_pl_snp_pos: HashSet<usize> = HashSet::new();
        for m in &pl_snp_maps { all_pl_snp_pos.extend(m.keys()); }
        let pl_total_core_snps = all_pl_snp_pos.len() as u32;

        // Discriminating: for each position, remove samples with gaps,
        // then check if remaining samples have discordant alleles
        let pl_disc = all_pl_snp_pos.iter().filter(|&&pos| {
            let alleles: Vec<Option<u8>> = (0..n)
                .filter(|&idx| !pl_gap_sets[idx].contains(&pos))
                .map(|idx| pl_snp_maps[idx].get(&pos).copied())
                .collect();
            if alleles.len() < 2 { return false; }
            let first = alleles[0];
            alleles.iter().any(|a| *a != first)
        }).count() as u32;

        // Helper: count same_pos and concordant for a set of GT disc positions
        let count_same_concordant = |disc: &[usize], active: &[usize]| -> (u32, u32) {
            let mut same_pos: u32 = 0;
            let mut concordant: u32 = 0;
            for &pos in disc {
                let pl_has_snp = active.iter().any(|&idx| pl_snp_maps[idx].contains_key(&pos));
                if !pl_has_snp { continue; }
                same_pos += 1;
                let any_pl_gap = active.iter().any(|&idx| pl_gap_sets[idx].contains(&pos));
                if any_pl_gap { continue; }
                let all_match = active.iter().all(|&idx| {
                    let gt_allele = gt_snp_maps[idx].get(&pos).copied()
                        .unwrap_or_else(|| if pos < ref_bytes.len() { ref_bytes[pos] } else { b'N' });
                    let pl_allele = pl_snp_maps[idx].get(&pos).copied()
                        .unwrap_or_else(|| if pos < ref_bytes.len() { ref_bytes[pos] } else { b'N' });
                    gt_allele == pl_allele
                });
                if all_match { concordant += 1; }
            }
            (same_pos, concordant)
        };

        let count_pl_in_gt_gaps = |gt_gaps: &HashSet<usize>| -> u32 {
            all_pl_snp_pos.iter().filter(|p| gt_gaps.contains(p)).count() as u32
        };

        let all_indices: Vec<usize> = (0..n).collect();

        // Gap-Union
        let (gu_same_pos, gu_concordant) = count_same_concordant(&gt_disc_union, &all_indices);
        let gu_pl_in_gaps = count_pl_in_gt_gaps(&gt_gap_union);

        // Gap-Intersect
        let mut gi_same_pos: u32 = 0;
        let mut gi_concordant: u32 = 0;
        for &pos in &gt_disc_intersect {
            let active: Vec<usize> = (0..n).filter(|&idx| !gt_gap_sets[idx].contains(&pos)).collect();
            if active.len() < 2 { continue; }
            let pl_has_snp = active.iter().any(|&idx| pl_snp_maps[idx].contains_key(&pos));
            if !pl_has_snp { continue; }
            gi_same_pos += 1;
            let any_pl_gap = active.iter().any(|&idx| pl_gap_sets[idx].contains(&pos));
            if any_pl_gap { continue; }
            let all_match = active.iter().all(|&idx| {
                let gt_allele = gt_snp_maps[idx].get(&pos).copied()
                    .unwrap_or_else(|| if pos < ref_bytes.len() { ref_bytes[pos] } else { b'N' });
                let pl_allele = pl_snp_maps[idx].get(&pos).copied()
                    .unwrap_or_else(|| if pos < ref_bytes.len() { ref_bytes[pos] } else { b'N' });
                gt_allele == pl_allele
            });
            if all_match { gi_concordant += 1; }
        }
        let gi_pl_in_gaps = count_pl_in_gt_gaps(&gt_gap_intersection);

        // Pairwise
        let mut pw_total_disc: f64 = 0.0;
        let mut pw_total_same_pos: f64 = 0.0;
        let mut pw_total_concordant: f64 = 0.0;
        let mut num_pairs: u32 = 0;

        for i in 0..n {
            for j in (i + 1)..n {
                let pair_gap_union: HashSet<usize> = gt_gap_sets[i].union(&gt_gap_sets[j]).cloned().collect();
                let pair_disc: Vec<usize> = all_gt_snp_positions.iter().filter(|&&pos| {
                    if pair_gap_union.contains(&pos) { return false; }
                    gt_snp_maps[i].get(&pos).copied() != gt_snp_maps[j].get(&pos).copied()
                }).copied().collect();

                let pair_indices = vec![i, j];
                let (p_same, p_conc) = count_same_concordant(&pair_disc, &pair_indices);

                pw_total_disc += pair_disc.len() as f64;
                pw_total_same_pos += p_same as f64;
                pw_total_concordant += p_conc as f64;
                num_pairs += 1;
            }
        }

        let np = if num_pairs > 0 { num_pairs as f64 } else { 1.0 };

        results.push(GtDiscVsPipelineResult {
            pipeline_id: pipeline_id.clone(),
            pl_total_core_snps: pl_total_core_snps,
            pl_discriminating_core_snps: pl_disc,
            gap_intersect: GapStrategyResult {
                gt_disc: gt_disc_intersect.len() as u32,
                same_pos: gi_same_pos,
                concordant: Some(gi_concordant),
                pl_snps_in_gt_gaps: gi_pl_in_gaps,
                position_details: None, // VCF fallback path doesn't collect details
            },
            gap_union: GapStrategyResult {
                gt_disc: gt_disc_union.len() as u32,
                same_pos: gu_same_pos,
                concordant: Some(gu_concordant),
                pl_snps_in_gt_gaps: gu_pl_in_gaps,
                position_details: None,
            },
            pairwise: PairwiseGtDiscResult {
                gt_disc_avg: pw_total_disc / np,
                same_pos_avg: pw_total_same_pos / np,
                concordant_avg: Some(pw_total_concordant / np),
                num_pairs,
            },
        });
    }

    results
}

/// Compute pre-computed KPIs and per-pipeline statistics from the data
fn compute_all_stats(
    data: &HashMap<String, HashMap<String, PipelineData>>,
    sample_ids: &[String],
    pipeline_ids: &[String],
    pipelines: &HashMap<String, PipelineInfo>,
    ref_length: usize,
    ref_seq: &str,
) -> (PreComputedKpis, HashMap<String, PipelineStats>) {
    use std::collections::HashSet;

    let ref_bytes = ref_seq.as_bytes();
    let n = sample_ids.len();
    let mut kpi_map = HashMap::new();
    let mut stats_map = HashMap::new();

    // Find GT pipeline
    let gt_pipeline_id: Option<&str> = pipelines.iter()
        .find(|(_, info)| info.ground_truth)
        .map(|(id, _)| id.as_str());

    for pipeline_id in pipeline_ids {
        // Build per-sample gap sets and SNP maps for this pipeline
        let mut gap_sets: Vec<HashSet<usize>> = Vec::with_capacity(n);
        let mut snp_maps: Vec<HashMap<usize, u8>> = Vec::with_capacity(n);

        for sample_id in sample_ids {
            let mut gaps = HashSet::new();
            let mut snps = HashMap::new();

            if let Some(sample_data) = data.get(sample_id) {
                if let Some(pd) = sample_data.get(pipeline_id) {
                    for gap in &pd.gaps {
                        for pos in gap[0]..gap[1] {
                            gaps.insert(pos);
                        }
                    }
                    for snp in &pd.snps {
                        let pos = snp.pos - 1; // 1-based to 0-based
                        let alt = snp.alt.as_bytes().first().copied().unwrap_or(b'N');
                        if pos < ref_bytes.len() && alt != ref_bytes[pos] {
                            snps.insert(pos, alt);
                        }
                    }
                }
            }

            gap_sets.push(gaps);
            snp_maps.push(snps);
        }

        // Gap-intersect: ALL samples have gap
        let gap_intersect: HashSet<usize> = if n > 0 {
            let mut isect = gap_sets[0].clone();
            for gs in &gap_sets[1..] {
                isect = isect.intersection(gs).cloned().collect();
            }
            isect
        } else {
            HashSet::new()
        };

        // Gap-union: ANY sample has gap
        let mut gap_union: HashSet<usize> = HashSet::new();
        for gs in &gap_sets {
            gap_union.extend(gs);
        }

        // All SNP positions
        let mut all_snp_positions: HashSet<usize> = HashSet::new();
        for sm in &snp_maps {
            all_snp_positions.extend(sm.keys());
        }

        // Count gap positions
        let total_gap_positions = gap_union.len() as u32;

        // Classify SNPs for each gap strategy
        let classify = |excluded_gaps: &HashSet<usize>| -> (u32, u32, u32, u32) {
            let mut total = 0u32;
            let mut consensus = 0u32;
            let mut disc = 0u32;
            let mut missing_vcf = 0u32;

            for &pos in &all_snp_positions {
                if excluded_gaps.contains(&pos) { continue; }
                total += 1;

                let alleles: Vec<Option<u8>> = (0..n).map(|idx| {
                    if gap_sets[idx].contains(&pos) {
                        None // gapped
                    } else {
                        Some(snp_maps[idx].get(&pos).copied()
                            .unwrap_or_else(|| ref_bytes.get(pos).copied().unwrap_or(b'N')))
                    }
                }).collect();

                let present: Vec<u8> = alleles.iter().filter_map(|a| *a).collect();
                if present.is_empty() { continue; }

                // Check if any sample is missing (gapped)
                let has_missing = alleles.iter().any(|a| a.is_none());
                if has_missing {
                    missing_vcf += 1;
                }

                if present.len() < 2 { continue; }
                let first = present[0];
                if present.iter().all(|&a| a == first) {
                    // Check if all are ref or all are same alt
                    if first != ref_bytes.get(pos).copied().unwrap_or(b'N') {
                        consensus += 1;
                    }
                    // If all are ref, it shouldn't be in all_snp_positions, but could happen with gap filtering
                } else {
                    disc += 1;
                }
            }

            (total, consensus, disc, missing_vcf)
        };

        let (total_i, consensus_i, disc_i, missing_i) = classify(&gap_intersect);
        let (total_u, consensus_u, disc_u, missing_u) = classify(&gap_union);

        // Discriminating SNP breakdown (for VCF pipelines with GT available)
        let disc_breakdown = if let Some(gt_id) = gt_pipeline_id {
            if pipeline_id != gt_id && pipelines.get(pipeline_id).map(|p| p.has_vcf).unwrap_or(false) {
                // Build GT SNP maps for breakdown analysis
                let mut gt_snp_maps_for_breakdown: Vec<HashMap<usize, u8>> = Vec::with_capacity(n);
                for sample_id in sample_ids {
                    let mut snps = HashMap::new();
                    if let Some(sample_data) = data.get(sample_id) {
                        if let Some(gt_data) = sample_data.get(gt_id) {
                            for snp in &gt_data.snps {
                                let pos = snp.pos - 1;
                                let alt = snp.alt.as_bytes().first().copied().unwrap_or(b'N');
                                if pos < ref_bytes.len() && alt != ref_bytes[pos] {
                                    snps.insert(pos, alt);
                                }
                            }
                        }
                    }
                    gt_snp_maps_for_breakdown.push(snps);
                }

                let mut gap_affected = 0u32;
                let mut gt_consensus_count = 0u32;
                let mut majority_rule = 0u32;
                let mut confirmed = 0u32;

                for &pos in &all_snp_positions {
                    if gap_union.contains(&pos) { continue; }

                    // Check if this is a discriminating position
                    let alleles: Vec<Option<u8>> = (0..n).map(|idx| {
                        if gap_sets[idx].contains(&pos) {
                            None
                        } else {
                            Some(snp_maps[idx].get(&pos).copied()
                                .unwrap_or_else(|| ref_bytes.get(pos).copied().unwrap_or(b'N')))
                        }
                    }).collect();
                    let present: Vec<u8> = alleles.iter().filter_map(|a| *a).collect();
                    if present.len() < 2 { continue; }
                    let first = present[0];
                    if present.iter().all(|&a| a == first) { continue; }

                    // It's discriminating - classify
                    let has_gap = alleles.iter().any(|a| a.is_none());
                    if has_gap {
                        gap_affected += 1;
                        continue;
                    }

                    // Check GT consensus at this position
                    let gt_alleles: Vec<u8> = (0..n).map(|idx| {
                        gt_snp_maps_for_breakdown[idx].get(&pos).copied()
                            .unwrap_or_else(|| ref_bytes.get(pos).copied().unwrap_or(b'N'))
                    }).collect();
                    let gt_first = gt_alleles[0];
                    if gt_alleles.iter().all(|&a| a == gt_first) {
                        gt_consensus_count += 1;
                        continue;
                    }

                    // Check majority rule: all but one agree
                    let mut counts: HashMap<u8, usize> = HashMap::new();
                    for &a in &present {
                        *counts.entry(a).or_default() += 1;
                    }
                    let max_count = counts.values().max().copied().unwrap_or(0);
                    if max_count == present.len() - 1 {
                        majority_rule += 1;
                    } else {
                        confirmed += 1;
                    }
                }

                Some(DiscBreakdown {
                    gap_affected,
                    gt_consensus: gt_consensus_count,
                    majority_rule,
                    confirmed,
                })
            } else {
                None
            }
        } else {
            None
        };

        let global = GlobalPipelineStats {
            ref_length: ref_length as u32,
            gap_intersect_count: gap_intersect.len() as u32,
            gap_union_count: gap_union.len() as u32,
            usable_intersect: (ref_length - gap_intersect.len()) as u32,
            usable_union: (ref_length - gap_union.len()) as u32,
            total_snps_intersect: total_i,
            total_snps_union: total_u,
            consensus_snps_intersect: consensus_i,
            consensus_snps_union: consensus_u,
            disc_snps_intersect: disc_i,
            disc_snps_union: disc_u,
            missing_vcf_intersect: missing_i,
            missing_vcf_union: missing_u,
            disc_breakdown,
        };

        // Pairwise stats (gap-union per pair)
        let mut pair_disc_snps: Vec<u32> = Vec::new();
        let mut pair_usable: Vec<u32> = Vec::new();
        let mut sample_disc_totals: HashMap<String, (f64, f64, u32)> = HashMap::new(); // (usable_sum, disc_sum, count)
        let mut sample_pairs: HashMap<String, HashMap<String, PairDetail>> = HashMap::new();

        for i in 0..n {
            for j in (i+1)..n {
                let pair_gap_union: HashSet<usize> = gap_sets[i].union(&gap_sets[j]).cloned().collect();
                let usable = (ref_length - pair_gap_union.len()) as u32;

                let mut disc = 0u32;
                for &pos in &all_snp_positions {
                    if pair_gap_union.contains(&pos) { continue; }
                    let a_i = snp_maps[i].get(&pos).copied()
                        .unwrap_or_else(|| ref_bytes.get(pos).copied().unwrap_or(b'N'));
                    let a_j = snp_maps[j].get(&pos).copied()
                        .unwrap_or_else(|| ref_bytes.get(pos).copied().unwrap_or(b'N'));
                    if a_i != a_j {
                        disc += 1;
                    }
                }

                pair_disc_snps.push(disc);
                pair_usable.push(usable);

                let detail = PairDetail { usable_space: usable, disc_snps: disc };

                // Accumulate per-sample stats + store pair details
                for &idx in &[i, j] {
                    let other_idx = if idx == i { j } else { i };
                    let entry = sample_disc_totals.entry(sample_ids[idx].clone()).or_default();
                    entry.0 += usable as f64;
                    entry.1 += disc as f64;
                    entry.2 += 1;
                    sample_pairs.entry(sample_ids[idx].clone())
                        .or_default()
                        .insert(sample_ids[other_idx].clone(), detail.clone());
                }
            }
        }

        let num_pairs = pair_disc_snps.len() as u32;
        let (disc_min, disc_max, disc_avg, disc_median, usable_avg) = if num_pairs > 0 {
            let mut sorted = pair_disc_snps.clone();
            sorted.sort();
            let min = sorted[0];
            let max = *sorted.last().unwrap();
            let avg = sorted.iter().map(|&x| x as f64).sum::<f64>() / num_pairs as f64;
            let median = if sorted.len() % 2 == 0 {
                (sorted[sorted.len()/2 - 1] as f64 + sorted[sorted.len()/2] as f64) / 2.0
            } else {
                sorted[sorted.len()/2] as f64
            };
            let usable_avg = pair_usable.iter().map(|&x| x as f64).sum::<f64>() / num_pairs as f64;
            (min, max, avg, median, usable_avg)
        } else {
            (0, 0, 0.0, 0.0, 0.0)
        };

        let per_sample: HashMap<String, SamplePairwiseStats> = sample_disc_totals.into_iter()
            .map(|(sample, (us, ds, cnt))| {
                let cnt = cnt as f64;
                let pairs = sample_pairs.remove(&sample).unwrap_or_default();
                (sample, SamplePairwiseStats {
                    avg_usable_space: if cnt > 0.0 { us / cnt } else { 0.0 },
                    avg_disc_snps: if cnt > 0.0 { ds / cnt } else { 0.0 },
                    pairs,
                })
            })
            .collect();

        let pairwise = PairwisePipelineStats {
            num_pairs,
            disc_snps_min: disc_min,
            disc_snps_median: disc_median,
            disc_snps_max: disc_max,
            disc_snps_avg: disc_avg,
            usable_space_avg: usable_avg,
            per_sample,
        };

        // Basic KPIs
        let core_snps = disc_u; // disc_snps using gap-union is the core SNPs count
        kpi_map.insert(pipeline_id.clone(), PipelineKpi {
            total_snps: all_snp_positions.len() as u32,
            total_gap_positions,
            core_snps,
        });

        stats_map.insert(pipeline_id.clone(), PipelineStats {
            global,
            pairwise,
        });
    }

    (PreComputedKpis { pipelines: kpi_map }, stats_map)
}

impl CompareReport {
    /// Generate report from configuration
    pub fn from_config_with_yaml(config: &Config, config_yaml: Option<String>) -> Result<Self> {
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

        let summary = Summary {
            total_samples: sample_ids.len(),
            total_pipelines: pipeline_ids.len(),
            generated_at: chrono::Utc::now().to_rfc3339(),
            coreguard_version: env!("CARGO_PKG_VERSION").to_string(),
            warnings,
            pileup_options: Some(PileupOptions {
                min_depth: config.options.min_depth,
                min_qual: config.options.min_qual,
                min_consensus: config.options.min_consensus,
                include_indels: config.options.include_indels,
            }),
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
                    match parsers::parse_core_snps(core_snps_path) {
                        Ok(core_data) => {
                            log::info!(
                                "Loaded core SNPs for pipeline '{}': {} positions ({} discriminating, has_alleles: {})",
                                pipeline_id, core_data.positions.len(), core_data.discriminating_count, core_data.has_alleles
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
                    log::info!("Computed GT disc vs pipelines (core_snps) for {} pipeline(s)", results.len());
                    Some(results)
                } else {
                    None
                }
            } else {
                // Fallback: compute from VCF data when no core_snps configured
                log::info!("No core_snps configured, computing GT disc vs pipelines from VCF data");
                let results = compute_gt_disc_vs_pipelines_from_vcf(
                    &data,
                    &sample_ids,
                    &pipeline_ids,
                    &gt_pipeline_id,
                    &ref_seq,
                );
                if !results.is_empty() {
                    log::info!("Computed GT disc vs pipelines (VCF) for {} pipeline(s)", results.len());
                    Some(results)
                } else {
                    None
                }
            }
        } else {
            None
        };

        // Pre-compute KPIs and pipeline stats
        log::info!("Pre-computing pipeline statistics...");
        let (kpis, pipeline_stats) = compute_all_stats(
            &data,
            &sample_ids,
            &pipeline_ids,
            &pipelines,
            ref_length,
            &ref_seq,
        );
        log::info!("Pre-computed stats for {} pipeline(s)", pipeline_stats.len());

        Ok(CompareReport {
            version: "2.0".to_string(),
            reference: ReferenceInfo {
                name: ref_name.clone(),
                label: config.reference.label.clone(),
                length: ref_length,
                sequence: ref_seq,
            },
            samples,
            pipelines,
            data,
            kpis: Some(kpis),
            pipeline_stats: Some(pipeline_stats),
            summary,
            description,
            pipeline_distance_matrices,
            gt_disc_vs_pipelines,
            config_yaml,
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

