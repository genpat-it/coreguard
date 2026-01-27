//! VCF comparison between SNP pipelines
//!
//! Compares variant calls between Snippy, CFSAN, and other pipelines
//! to identify discordant positions and potential causes.

use crate::bam_validate::BamValidationSummary;
use crate::vcf::{VcfFile, VariantCall};
use crate::gaps::Region;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

/// Comparison result for a single position
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PositionComparison {
    /// Chromosome/contig
    pub chrom: String,
    /// Position (1-based)
    pub pos: usize,
    /// Sample name
    pub sample: String,
    /// Variant in pipeline A (e.g., Snippy)
    pub pipeline_a: Option<VariantInfo>,
    /// Variant in pipeline B (e.g., CFSAN)
    pub pipeline_b: Option<VariantInfo>,
    /// Concordance status
    pub status: ConcordanceStatus,
    /// Coreguard coverage at this position (if available)
    pub coreguard_coverage: Option<u32>,
    /// Whether this position is in a coreguard gap region
    pub in_gap_region: bool,
    /// Likely cause of discordance
    pub likely_cause: Option<String>,
}

/// Simplified variant info for comparison
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariantInfo {
    /// Reference allele
    pub ref_allele: String,
    /// Alternate allele
    pub alt_allele: String,
    /// Quality score
    pub qual: f64,
    /// Read depth
    pub depth: Option<u32>,
    /// Allele frequency
    pub allele_freq: Option<f64>,
    /// Filter status
    pub filter: String,
    /// Genotype
    pub genotype: Option<String>,
    // === Additional quality metrics ===
    /// Reference allele observations
    pub ref_obs: Option<u32>,
    /// Alternate allele observations
    pub alt_obs: Option<u32>,
    /// Quality of alternate allele reads
    pub alt_qual: Option<u32>,
    /// Genotype quality
    pub genotype_qual: Option<u32>,
    /// Average base quality
    pub avg_base_qual: Option<u32>,
    /// P-value
    pub pvalue: Option<f64>,
}

impl From<&VariantCall> for VariantInfo {
    fn from(v: &VariantCall) -> Self {
        VariantInfo {
            ref_allele: v.ref_allele.clone(),
            alt_allele: v.alt_allele.clone(),
            qual: v.qual,
            depth: v.depth,
            allele_freq: v.allele_freq,
            filter: v.filter.clone(),
            genotype: v.genotype.clone(),
            ref_obs: v.ref_obs,
            alt_obs: v.alt_obs,
            alt_qual: v.alt_qual,
            genotype_qual: v.genotype_qual,
            avg_base_qual: v.avg_base_qual,
            pvalue: v.pvalue,
        }
    }
}

/// Concordance status between two pipelines
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConcordanceStatus {
    /// Both pipelines called the same variant
    Concordant,
    /// Both called a variant but different alleles
    DiscordantAllele,
    /// Only pipeline A called a variant
    OnlyPipelineA,
    /// Only pipeline B called a variant
    OnlyPipelineB,
    /// Both called no variant (reference)
    BothReference,
}

impl ConcordanceStatus {
    pub fn as_str(&self) -> &'static str {
        match self {
            ConcordanceStatus::Concordant => "Concordant",
            ConcordanceStatus::DiscordantAllele => "DiscordantAllele",
            ConcordanceStatus::OnlyPipelineA => "OnlyPipelineA",
            ConcordanceStatus::OnlyPipelineB => "OnlyPipelineB",
            ConcordanceStatus::BothReference => "BothReference",
        }
    }
}

/// Summary of comparison between two pipelines
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineComparisonSummary {
    /// Pipeline A name
    pub pipeline_a: String,
    /// Pipeline B name
    pub pipeline_b: String,
    /// Number of samples compared
    pub num_samples: usize,
    /// Total positions compared
    pub total_positions: usize,
    /// Concordant positions
    pub concordant: usize,
    /// Discordant alleles
    pub discordant_allele: usize,
    /// Only in pipeline A
    pub only_pipeline_a: usize,
    /// Only in pipeline B
    pub only_pipeline_b: usize,
    /// Concordance rate (%)
    pub concordance_rate: f64,
    /// Positions unique to A that are in coreguard gap regions
    pub a_only_in_gaps: usize,
    /// Positions unique to B that are in coreguard gap regions
    pub b_only_in_gaps: usize,
    /// Interpretation
    pub interpretation: String,
    /// BAM validation results (to detect Snippy artifacts)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bam_validation: Option<BamValidationSummary>,
}

/// Per-sample comparison results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleComparison {
    /// Sample name
    pub sample: String,
    /// Total SNPs in pipeline A
    pub snps_a: usize,
    /// Total SNPs in pipeline B
    pub snps_b: usize,
    /// Concordant SNPs
    pub concordant: usize,
    /// Only in A
    pub only_a: usize,
    /// Only in B
    pub only_b: usize,
    /// Concordance rate
    pub concordance_rate: f64,
    /// Discordant positions - skipped in JSON to avoid huge files
    #[serde(skip)]
    pub discordant_positions: Vec<PositionComparison>,
}

/// Complete VCF comparison results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VcfComparison {
    /// Pipeline A name
    pub pipeline_a: String,
    /// Pipeline B name
    pub pipeline_b: String,
    /// Summary statistics
    pub summary: PipelineComparisonSummary,
    /// Per-sample comparisons
    pub sample_comparisons: Vec<SampleComparison>,
    /// All discordant positions (for detailed analysis) - skipped in JSON to avoid huge files
    #[serde(skip)]
    pub discordant_positions: Vec<PositionComparison>,
    /// Top missed SNPs by pipeline B (with good quality in A)
    pub top_missed_by_b: Vec<PositionComparison>,
    /// Top missed SNPs by pipeline A (with good quality in B)
    pub top_missed_by_a: Vec<PositionComparison>,
}

/// Paginated response for discordant positions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PaginatedPositions {
    /// Total number of positions matching the filter
    pub total: usize,
    /// Current page (0-indexed)
    pub page: usize,
    /// Page size
    pub page_size: usize,
    /// Total pages
    pub total_pages: usize,
    /// Available samples for filtering
    pub available_samples: Vec<String>,
    /// Available statuses for filtering
    pub available_statuses: Vec<String>,
    /// Positions on this page
    pub positions: Vec<PositionComparison>,
}

impl VcfComparison {
    /// Get paginated discordant positions with optional filtering
    pub fn get_paginated_positions(
        &self,
        page: usize,
        page_size: usize,
        sample_filter: Option<&str>,
        status_filter: Option<&str>,
    ) -> PaginatedPositions {
        // Use discordant_positions which contains all positions from all samples
        // (already aggregated during comparison)
        let all_positions: Vec<&PositionComparison> = self.discordant_positions.iter().collect();

        // Get available samples and statuses
        let mut available_samples: Vec<String> = all_positions.iter()
            .map(|p| p.sample.clone())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        available_samples.sort();

        let available_statuses = vec![
            "Concordant".to_string(),
            "OnlyPipelineA".to_string(),
            "OnlyPipelineB".to_string(),
            "DiscordantAllele".to_string(),
        ];

        // Apply filters
        let filtered: Vec<&PositionComparison> = all_positions.into_iter()
            .filter(|p| {
                let sample_match = sample_filter
                    .map(|s| s == "all" || p.sample == s)
                    .unwrap_or(true);
                let status_match = status_filter
                    .map(|s| s == "all" || p.status.as_str() == s)
                    .unwrap_or(true);
                sample_match && status_match
            })
            .collect();

        let total = filtered.len();
        let total_pages = (total + page_size - 1) / page_size;
        let start = page * page_size;
        let end = std::cmp::min(start + page_size, total);

        let positions = if start < total {
            filtered[start..end].iter().cloned().cloned().collect()
        } else {
            Vec::new()
        };

        PaginatedPositions {
            total,
            page,
            page_size,
            total_pages,
            available_samples,
            available_statuses,
            positions,
        }
    }

    /// Set BAM validation results to detect Snippy artifacts
    pub fn set_bam_validation(&mut self, validation: BamValidationSummary) {
        self.summary.bam_validation = Some(validation);
    }

    /// Get Snippy-only positions for BAM validation
    /// Returns: Vec<(chrom, pos, sample, alt_allele)>
    pub fn get_snippy_only_positions(&self) -> Vec<(String, usize, String, String)> {
        self.discordant_positions
            .iter()
            .filter(|p| matches!(p.status, ConcordanceStatus::OnlyPipelineA))
            .filter_map(|p| {
                let alt = p.pipeline_a.as_ref()?.alt_allele.clone();
                Some((p.chrom.clone(), p.pos, p.sample.clone(), alt))
            })
            .collect()
    }
}

/// Compare variants between two pipelines
pub fn compare_pipelines(
    vcfs_a: &[VcfFile],
    vcfs_b: &[VcfFile],
    pipeline_a_name: &str,
    pipeline_b_name: &str,
    gap_regions: Option<&HashMap<String, Vec<Region>>>,
    coverage_data: Option<&HashMap<String, HashMap<usize, u32>>>,
) -> VcfComparison {
    let mut sample_comparisons = Vec::new();
    let mut all_discordant = Vec::new();
    let mut top_missed_by_b = Vec::new();
    let mut top_missed_by_a = Vec::new();

    let mut total_concordant = 0usize;
    let mut total_discordant_allele = 0usize;
    let mut total_only_a = 0usize;
    let mut total_only_b = 0usize;
    let mut total_positions = 0usize;
    let mut a_only_in_gaps = 0usize;
    let mut b_only_in_gaps = 0usize;

    // Index VCFs by sample name
    let vcfs_a_by_sample: HashMap<String, &VcfFile> = vcfs_a.iter()
        .map(|v| (v.sample.clone(), v))
        .collect();
    let vcfs_b_by_sample: HashMap<String, &VcfFile> = vcfs_b.iter()
        .map(|v| (v.sample.clone(), v))
        .collect();

    // Find common samples
    let samples_a: HashSet<&String> = vcfs_a_by_sample.keys().collect();
    let samples_b: HashSet<&String> = vcfs_b_by_sample.keys().collect();
    let common_samples: Vec<&String> = samples_a.intersection(&samples_b).cloned().collect();

    for sample in &common_samples {
        let vcf_a = vcfs_a_by_sample.get(*sample).unwrap();
        let vcf_b = vcfs_b_by_sample.get(*sample).unwrap();

        let sample_gaps = gap_regions.and_then(|g| g.get(*sample));
        let sample_coverage = coverage_data.and_then(|c| c.get(*sample));

        let comparison = compare_sample_vcfs(
            vcf_a,
            vcf_b,
            sample_gaps,
            sample_coverage,
            pipeline_a_name,
            pipeline_b_name,
        );

        total_concordant += comparison.concordant;
        total_only_a += comparison.only_a;
        total_only_b += comparison.only_b;
        total_positions += comparison.snps_a.max(comparison.snps_b);

        // Count gap-related discordances
        for pos in &comparison.discordant_positions {
            if pos.in_gap_region {
                match pos.status {
                    ConcordanceStatus::OnlyPipelineA => a_only_in_gaps += 1,
                    ConcordanceStatus::OnlyPipelineB => b_only_in_gaps += 1,
                    _ => {}
                }
            }
            all_discordant.push(pos.clone());
        }

        // Collect top missed SNPs
        for pos in &comparison.discordant_positions {
            match pos.status {
                ConcordanceStatus::OnlyPipelineA => {
                    if let Some(ref info) = pos.pipeline_a {
                        if info.qual > 100.0 && info.depth.unwrap_or(0) > 10 {
                            top_missed_by_b.push(pos.clone());
                        }
                    }
                }
                ConcordanceStatus::OnlyPipelineB => {
                    if let Some(ref info) = pos.pipeline_b {
                        if info.qual > 100.0 && info.depth.unwrap_or(0) > 10 {
                            top_missed_by_a.push(pos.clone());
                        }
                    }
                }
                _ => {}
            }
        }

        sample_comparisons.push(comparison);
    }

    // Sort top missed by quality
    top_missed_by_b.sort_by(|a, b| {
        let qual_a = a.pipeline_a.as_ref().map(|v| v.qual).unwrap_or(0.0);
        let qual_b = b.pipeline_a.as_ref().map(|v| v.qual).unwrap_or(0.0);
        qual_b.partial_cmp(&qual_a).unwrap_or(std::cmp::Ordering::Equal)
    });
    top_missed_by_b.truncate(50);

    top_missed_by_a.sort_by(|a, b| {
        let qual_a = a.pipeline_b.as_ref().map(|v| v.qual).unwrap_or(0.0);
        let qual_b = b.pipeline_b.as_ref().map(|v| v.qual).unwrap_or(0.0);
        qual_b.partial_cmp(&qual_a).unwrap_or(std::cmp::Ordering::Equal)
    });
    top_missed_by_a.truncate(50);

    let concordance_rate = if total_positions > 0 {
        (total_concordant as f64 / total_positions as f64) * 100.0
    } else {
        0.0
    };

    let interpretation = generate_interpretation(
        total_concordant,
        total_only_a,
        total_only_b,
        a_only_in_gaps,
        b_only_in_gaps,
        pipeline_a_name,
        pipeline_b_name,
    );

    let summary = PipelineComparisonSummary {
        pipeline_a: pipeline_a_name.to_string(),
        pipeline_b: pipeline_b_name.to_string(),
        num_samples: common_samples.len(),
        total_positions,
        concordant: total_concordant,
        discordant_allele: total_discordant_allele,
        only_pipeline_a: total_only_a,
        only_pipeline_b: total_only_b,
        concordance_rate,
        a_only_in_gaps,
        b_only_in_gaps,
        interpretation,
        bam_validation: None,
    };

    VcfComparison {
        pipeline_a: pipeline_a_name.to_string(),
        pipeline_b: pipeline_b_name.to_string(),
        summary,
        sample_comparisons,
        discordant_positions: all_discordant,
        top_missed_by_b,
        top_missed_by_a,
    }
}

/// Compare VCFs for a single sample
fn compare_sample_vcfs(
    vcf_a: &VcfFile,
    vcf_b: &VcfFile,
    gap_regions: Option<&Vec<Region>>,
    _coverage_data: Option<&HashMap<usize, u32>>,
    pipeline_a_name: &str,
    pipeline_b_name: &str,
) -> SampleComparison {
    let snps_a = vcf_a.passed_snps();
    let snps_b = vcf_b.passed_snps();

    // Collect all positions
    let mut all_positions: HashSet<String> = HashSet::new();
    for snp in &snps_a {
        all_positions.insert(snp.position_key());
    }
    for snp in &snps_b {
        all_positions.insert(snp.position_key());
    }

    let mut concordant = 0usize;
    let mut only_a = 0usize;
    let mut only_b = 0usize;
    let mut discordant_allele = 0usize;
    let mut discordant_positions = Vec::new();

    for pos_key in &all_positions {
        // Only consider true SNPs (ref and alt both length 1)
        // This filters out complex variants (MNPs, indels)
        let var_a = vcf_a.by_position.get(pos_key).filter(|v| v.is_snp() && v.is_pass());
        let var_b = vcf_b.by_position.get(pos_key).filter(|v| v.is_snp() && v.is_pass());

        // Parse position
        let parts: Vec<&str> = pos_key.split(':').collect();
        let chrom = parts[0].to_string();
        let pos: usize = parts.get(1).and_then(|p| p.parse().ok()).unwrap_or(0);

        // Check if in gap region
        let in_gap = gap_regions.map(|gaps| {
            gaps.iter().any(|g| pos >= g.start && pos <= g.end)
        }).unwrap_or(false);

        let (status, likely_cause) = match (var_a, var_b) {
            (Some(a), Some(b)) => {
                if a.alt_allele == b.alt_allele {
                    (ConcordanceStatus::Concordant, None)
                } else {
                    (ConcordanceStatus::DiscordantAllele, Some("Different alleles called".to_string()))
                }
            }
            (Some(a), None) => {
                let cause = infer_cause_for_missing(a, pipeline_b_name, in_gap);
                (ConcordanceStatus::OnlyPipelineA, Some(cause))
            }
            (None, Some(b)) => {
                let cause = infer_cause_for_missing(b, pipeline_a_name, in_gap);
                (ConcordanceStatus::OnlyPipelineB, Some(cause))
            }
            (None, None) => (ConcordanceStatus::BothReference, None),
        };

        match status {
            ConcordanceStatus::Concordant => concordant += 1,
            ConcordanceStatus::OnlyPipelineA => only_a += 1,
            ConcordanceStatus::OnlyPipelineB => only_b += 1,
            ConcordanceStatus::DiscordantAllele => discordant_allele += 1,
            ConcordanceStatus::BothReference => {}
        }

        // Record ALL positions (including concordant) for SNP Tracks visualization
        if status != ConcordanceStatus::BothReference {
            discordant_positions.push(PositionComparison {
                chrom,
                pos,
                sample: vcf_a.sample.clone(),
                pipeline_a: var_a.map(VariantInfo::from),
                pipeline_b: var_b.map(VariantInfo::from),
                status,
                coreguard_coverage: None,
                in_gap_region: in_gap,
                likely_cause,
            });
        }
    }

    let concordance_rate = if !all_positions.is_empty() {
        (concordant as f64 / all_positions.len() as f64) * 100.0
    } else {
        100.0
    };

    SampleComparison {
        sample: vcf_a.sample.clone(),
        snps_a: snps_a.len(),
        snps_b: snps_b.len(),
        concordant,
        only_a,
        only_b,
        concordance_rate,
        discordant_positions,
    }
}

/// Infer likely cause for a missing variant
fn infer_cause_for_missing(present_var: &VariantCall, missing_pipeline: &str, in_gap: bool) -> String {
    let mut causes = Vec::new();

    if in_gap {
        causes.push("In coreguard gap region (low coverage)");
    }

    if let Some(dp) = present_var.depth {
        if dp < 10 {
            causes.push("Low read depth (<10)");
        }
    }

    if present_var.qual < 100.0 {
        causes.push("Low quality score (<100)");
    }

    if let Some(af) = present_var.allele_freq {
        if af < 0.9 {
            causes.push("Allele frequency <90%");
        }
    }

    if causes.is_empty() {
        // High quality variant missed
        if missing_pipeline.to_lowercase().contains("cfsan") {
            causes.push("Likely VarScan2 sensitivity issue");
        } else if missing_pipeline.to_lowercase().contains("snippy") {
            causes.push("Likely FreeBayes sensitivity issue");
        } else {
            causes.push("Pipeline-specific filter or algorithm difference");
        }
    }

    causes.join("; ")
}

/// Generate interpretation text
fn generate_interpretation(
    concordant: usize,
    only_a: usize,
    only_b: usize,
    a_only_in_gaps: usize,
    b_only_in_gaps: usize,
    pipeline_a: &str,
    pipeline_b: &str,
) -> String {
    let total = concordant + only_a + only_b;
    if total == 0 {
        return "No variants to compare".to_string();
    }

    let concordance = (concordant as f64 / total as f64) * 100.0;

    let mut interpretation = if concordance > 95.0 {
        format!("Excellent agreement ({:.1}% concordance)", concordance)
    } else if concordance > 80.0 {
        format!("Good agreement ({:.1}% concordance)", concordance)
    } else if concordance > 50.0 {
        format!("Moderate agreement ({:.1}% concordance) - review discordant positions", concordance)
    } else {
        format!("Poor agreement ({:.1}% concordance) - significant pipeline differences", concordance)
    };

    // Add specific insights
    if only_a > only_b * 2 {
        interpretation.push_str(&format!(
            ". {} found {}x more unique SNPs than {} - possible {} sensitivity issue",
            pipeline_a,
            only_a / only_b.max(1),
            pipeline_b,
            pipeline_b
        ));
    } else if only_b > only_a * 2 {
        interpretation.push_str(&format!(
            ". {} found {}x more unique SNPs than {} - possible {} sensitivity issue",
            pipeline_b,
            only_b / only_a.max(1),
            pipeline_a,
            pipeline_a
        ));
    }

    if a_only_in_gaps > 0 || b_only_in_gaps > 0 {
        interpretation.push_str(&format!(
            ". {} SNPs unique to {} and {} unique to {} are in low-coverage regions",
            a_only_in_gaps, pipeline_a, b_only_in_gaps, pipeline_b
        ));
    }

    interpretation
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_concordance_status() {
        assert_eq!(ConcordanceStatus::Concordant.as_str(), "Concordant");
        assert_eq!(ConcordanceStatus::OnlyPipelineA.as_str(), "OnlyPipelineA");
    }
}
