//! SNP Pipeline Comparison Module
//!
//! Compares coreguard gap predictions with actual filtering done by SNP pipelines
//! like Snippy and CFSAN.
//!
//! Key metrics:
//! - Core size comparison: coreguard estimated core vs snippy-core actual aligned bases
//! - SNP impact: how many real SNPs fall within coreguard's predicted gap regions
//! - Low coverage correlation: LOWCOV positions from snippy-core vs coreguard gaps

use crate::gaps::Region;
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

// ============================================================================
// Snippy-core statistics from snippycore.txt
// ============================================================================

/// Per-sample statistics from snippy-core
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnippyCoreStats {
    /// Sample name
    pub sample: String,
    /// Reference length
    pub length: usize,
    /// Aligned bases (core for this sample)
    pub aligned: usize,
    /// Unaligned bases
    pub unaligned: usize,
    /// Variant positions
    pub variant: usize,
    /// Heterozygous calls
    pub het: usize,
    /// Masked positions
    pub masked: usize,
    /// Low coverage positions
    pub lowcov: usize,
}

/// Parse snippycore.txt to get per-sample statistics
pub fn parse_snippycore_txt(path: &Path) -> Result<Vec<SnippyCoreStats>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open snippycore.txt: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut stats = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with("ID\t") || line.trim().is_empty() {
            continue; // Skip header
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 8 {
            // Clean sample name: remove "_snippy" suffix if present
            let sample = parts[0].trim_end_matches("_snippy").to_string();

            // Skip Reference line
            if sample == "Reference" {
                continue;
            }

            stats.push(SnippyCoreStats {
                sample,
                length: parts[1].parse().unwrap_or(0),
                aligned: parts[2].parse().unwrap_or(0),
                unaligned: parts[3].parse().unwrap_or(0),
                variant: parts[4].parse().unwrap_or(0),
                het: parts[5].parse().unwrap_or(0),
                masked: parts[6].parse().unwrap_or(0),
                lowcov: parts[7].parse().unwrap_or(0),
            });
        }
    }

    Ok(stats)
}

// ============================================================================
// SNP positions from snippycore.tab
// ============================================================================

/// SNP position with per-sample alleles
#[derive(Debug, Clone)]
pub struct SnpPosition {
    /// Chromosome/contig
    pub chrom: String,
    /// Position (0-based)
    pub pos: usize,
    /// Reference allele
    pub ref_allele: String,
    /// Sample alleles (sample_name -> allele)
    pub alleles: HashMap<String, String>,
}

/// Parse snippycore.tab to get SNP positions
pub fn parse_snippycore_tab(path: &Path) -> Result<Vec<SnpPosition>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open snippycore.tab: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut snps = Vec::new();
    let mut sample_names: Vec<String> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();

        if line.starts_with("CHR\t") {
            // Header line - extract sample names
            sample_names = parts.iter()
                .skip(3) // Skip CHR, POS, REF
                .map(|s| s.trim_end_matches("_snippy").to_string())
                .collect();
            continue;
        }

        if parts.len() >= 4 && !sample_names.is_empty() {
            let pos: usize = parts[1].parse().unwrap_or(0);
            if pos == 0 {
                continue;
            }

            let mut alleles = HashMap::new();
            for (i, sample) in sample_names.iter().enumerate() {
                if let Some(allele) = parts.get(3 + i) {
                    alleles.insert(sample.clone(), allele.to_string());
                }
            }

            snps.push(SnpPosition {
                chrom: parts[0].to_string(),
                pos: pos - 1, // Convert to 0-based
                ref_allele: parts[2].to_string(),
                alleles,
            });
        }
    }

    log::info!("Parsed {} SNP positions from snippycore.tab", snps.len());
    Ok(snps)
}

// ============================================================================
// Redesigned comparison metrics
// ============================================================================

/// Core-level comparison for a single sample
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoreComparison {
    /// Sample name
    pub sample: String,

    // Coreguard metrics
    /// Coreguard estimated core size (non-gap positions)
    pub cg_core_size: usize,
    /// Coreguard gap bases
    pub cg_gap_bases: usize,
    /// Coreguard core percentage
    pub cg_core_pct: f64,

    // Snippy-core metrics
    /// Snippy-core aligned bases (actual core)
    pub sc_aligned: usize,
    /// Snippy-core low coverage bases
    pub sc_lowcov: usize,
    /// Snippy-core aligned percentage
    pub sc_aligned_pct: f64,

    // Comparison
    /// Core size difference (coreguard - snippy-core)
    pub core_diff: i64,
    /// Core size difference percentage
    pub core_diff_pct: f64,
    /// Low coverage overlap: positions flagged by both
    pub lowcov_overlap: usize,
    /// Coreguard gaps that are NOT low-cov in snippy
    pub cg_gaps_not_lowcov: usize,
    /// Snippy low-cov that are NOT in coreguard gaps
    pub sc_lowcov_not_in_gaps: usize,
}

/// SNP impact analysis for a single sample
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnpImpact {
    /// Sample name
    pub sample: String,
    /// Total SNPs for this sample
    pub total_snps: usize,
    /// SNPs falling in coreguard gap regions
    pub snps_in_gaps: usize,
    /// SNPs in non-gap regions (would be retained)
    pub snps_retained: usize,
    /// Percentage of SNPs that would be retained
    pub retention_pct: f64,
    /// Positions of SNPs that fall in gaps
    pub lost_snp_positions: Vec<usize>,
}

/// Complete snippy-core comparison
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnippyCoreComparison {
    /// Reference length
    pub reference_length: usize,
    /// Per-sample core comparisons
    pub core_comparisons: Vec<CoreComparison>,
    /// Per-sample SNP impact
    pub snp_impacts: Vec<SnpImpact>,
    /// Summary
    pub summary: SnippyCoreSummary,
    /// SNP distance matrix (if available)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub distance_matrix: Option<SnpDistanceMatrix>,
}

/// Summary of snippy-core comparison
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnippyCoreSummary {
    /// Number of samples compared
    pub total_samples: usize,
    /// Average core size difference (%)
    pub avg_core_diff_pct: f64,
    /// Average SNP retention (%)
    pub avg_snp_retention_pct: f64,
    /// Total SNPs across all samples
    pub total_snps: usize,
    /// Total SNPs that would be lost to gaps
    pub total_snps_in_gaps: usize,
    /// Interpretation
    pub interpretation: String,
}

/// Compare coreguard gaps with snippy-core results
pub fn compare_with_snippycore(
    coreguard_gaps: &[(String, Vec<Region>)],
    snippycore_stats: &[SnippyCoreStats],
    snp_positions: &[SnpPosition],
    reference_length: usize,
) -> SnippyCoreComparison {
    let mut core_comparisons = Vec::new();
    let mut snp_impacts = Vec::new();

    // Build maps for quick lookup
    let stats_map: HashMap<&str, &SnippyCoreStats> = snippycore_stats
        .iter()
        .map(|s| (s.sample.as_str(), s))
        .collect();

    for (sample_name, cg_regions) in coreguard_gaps {
        // Get snippy-core stats for this sample
        let sc_stats = match stats_map.get(sample_name.as_str()) {
            Some(s) => *s,
            None => {
                log::debug!("No snippy-core stats for sample: {}", sample_name);
                continue;
            }
        };

        // Calculate coreguard gap positions
        let cg_gap_positions: HashSet<usize> = cg_regions
            .iter()
            .flat_map(|r| r.start..r.end)
            .collect();
        let cg_gap_bases = cg_gap_positions.len();
        let cg_core_size = reference_length.saturating_sub(cg_gap_bases);
        let cg_core_pct = (cg_core_size as f64 / reference_length as f64) * 100.0;

        // Snippy-core metrics
        let sc_aligned_pct = (sc_stats.aligned as f64 / reference_length as f64) * 100.0;

        // Core difference
        let core_diff = cg_core_size as i64 - sc_stats.aligned as i64;
        let core_diff_pct = (core_diff as f64 / reference_length as f64) * 100.0;

        // TODO: For lowcov overlap, we'd need the actual positions from snippy
        // For now, we use the count as an approximation
        let lowcov_overlap = 0; // Would need position-level data
        let cg_gaps_not_lowcov = cg_gap_bases.saturating_sub(sc_stats.lowcov);
        let sc_lowcov_not_in_gaps = sc_stats.lowcov.saturating_sub(cg_gap_bases);

        core_comparisons.push(CoreComparison {
            sample: sample_name.clone(),
            cg_core_size,
            cg_gap_bases,
            cg_core_pct,
            sc_aligned: sc_stats.aligned,
            sc_lowcov: sc_stats.lowcov,
            sc_aligned_pct,
            core_diff,
            core_diff_pct,
            lowcov_overlap,
            cg_gaps_not_lowcov,
            sc_lowcov_not_in_gaps,
        });

        // SNP impact analysis
        let sample_snps: Vec<usize> = snp_positions
            .iter()
            .filter(|snp| {
                // Check if this sample has a variant (different from ref)
                snp.alleles.get(sample_name)
                    .map(|a| a != &snp.ref_allele && a != "-" && a != "N")
                    .unwrap_or(false)
            })
            .map(|snp| snp.pos)
            .collect();

        let total_snps = sample_snps.len();
        let snps_in_gaps: Vec<usize> = sample_snps
            .iter()
            .filter(|pos| cg_gap_positions.contains(pos))
            .copied()
            .collect();
        let snps_in_gaps_count = snps_in_gaps.len();
        let snps_retained = total_snps - snps_in_gaps_count;
        let retention_pct = if total_snps > 0 {
            (snps_retained as f64 / total_snps as f64) * 100.0
        } else {
            100.0
        };

        snp_impacts.push(SnpImpact {
            sample: sample_name.clone(),
            total_snps,
            snps_in_gaps: snps_in_gaps_count,
            snps_retained,
            retention_pct,
            lost_snp_positions: snps_in_gaps,
        });
    }

    // Calculate summary
    let total_samples = core_comparisons.len();
    let avg_core_diff_pct = if total_samples > 0 {
        core_comparisons.iter().map(|c| c.core_diff_pct.abs()).sum::<f64>() / total_samples as f64
    } else {
        0.0
    };
    let avg_snp_retention_pct = if !snp_impacts.is_empty() {
        snp_impacts.iter().map(|s| s.retention_pct).sum::<f64>() / snp_impacts.len() as f64
    } else {
        100.0
    };
    let total_snps: usize = snp_impacts.iter().map(|s| s.total_snps).sum();
    let total_snps_in_gaps: usize = snp_impacts.iter().map(|s| s.snps_in_gaps).sum();

    let interpretation = if avg_snp_retention_pct > 99.0 && avg_core_diff_pct < 1.0 {
        "Excellent: coreguard closely matches snippy-core with minimal SNP loss".to_string()
    } else if avg_snp_retention_pct > 95.0 {
        "Good: coreguard preserves most SNPs with slight core size differences".to_string()
    } else if avg_snp_retention_pct > 90.0 {
        "Moderate: some SNPs fall in coreguard gap regions - review gap thresholds".to_string()
    } else {
        "Warning: significant SNP loss in coreguard gaps - consider adjusting parameters".to_string()
    };

    SnippyCoreComparison {
        reference_length,
        core_comparisons,
        snp_impacts,
        summary: SnippyCoreSummary {
            total_samples,
            avg_core_diff_pct,
            avg_snp_retention_pct,
            total_snps,
            total_snps_in_gaps,
            interpretation,
        },
        distance_matrix: None, // Will be set by caller if available
    }
}

// ============================================================================
// CFSAN Comparison
// ============================================================================

/// Per-sample CFSAN metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CfsanSampleMetrics {
    /// Sample name
    pub sample: String,
    /// Percent of reads mapped
    pub percent_reads_mapped: f64,
    /// Average pileup depth
    pub average_depth: f64,
    /// Phase 2 SNPs (after quality filtering)
    pub phase2_snps: usize,
    /// Missing SNP matrix positions (filtered out for this sample)
    pub missing_positions: usize,
}

/// CFSAN comparison with coreguard
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CfsanComparison {
    /// Reference length
    pub reference_length: usize,
    /// CFSAN metrics per sample
    pub cfsan_metrics: Vec<CfsanSampleMetrics>,
    /// Total SNP positions in CFSAN analysis
    pub total_cfsan_snps: usize,
    /// Per-sample comparison
    pub sample_comparisons: Vec<CfsanSampleComparison>,
    /// Summary
    pub summary: CfsanSummary,
    /// SNP distance matrix (if available)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub distance_matrix: Option<SnpDistanceMatrix>,
}

/// Per-sample CFSAN comparison
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CfsanSampleComparison {
    /// Sample name
    pub sample: String,
    /// coreguard gap bases
    pub cg_gap_bases: usize,
    /// CFSAN SNPs for this sample
    pub cfsan_snps: usize,
    /// SNPs in coreguard gaps (would be lost)
    pub snps_in_gaps: usize,
    /// SNP retention percentage
    pub snp_retention_pct: f64,
    /// CFSAN missing positions that overlap with coreguard gaps
    pub missing_overlap: usize,
}

/// CFSAN comparison summary
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CfsanSummary {
    /// Total samples
    pub total_samples: usize,
    /// Average SNP retention
    pub avg_snp_retention_pct: f64,
    /// Total SNPs across samples
    pub total_snps: usize,
    /// Total SNPs in gaps
    pub total_snps_in_gaps: usize,
    /// Interpretation
    pub interpretation: String,
}

/// Parse CFSAN metrics.tsv
pub fn parse_cfsan_metrics(path: &Path) -> Result<Vec<CfsanSampleMetrics>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open CFSAN metrics.tsv: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut metrics = Vec::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i == 0 {
            continue; // Skip header
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 15 {
            // Remove quotes from sample name
            let sample = parts[0].trim_matches('"').to_string();

            metrics.push(CfsanSampleMetrics {
                sample,
                percent_reads_mapped: parts[7].parse().unwrap_or(0.0),
                average_depth: parts[10].parse().unwrap_or(0.0),
                phase2_snps: parts[13].parse().unwrap_or(0),
                missing_positions: parts[15].parse().unwrap_or(0),
            });
        }
    }

    Ok(metrics)
}

/// Parse CFSAN snplist.txt to get SNP positions per sample
pub fn parse_cfsan_snplist_full(path: &Path) -> Result<(usize, HashMap<String, Vec<usize>>)> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open CFSAN snplist: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut sample_snps: HashMap<String, Vec<usize>> = HashMap::new();
    let mut total_positions = 0usize;

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 4 {
            let pos: usize = parts[1].parse().unwrap_or(0);
            if pos == 0 {
                continue;
            }
            let pos = pos - 1; // Convert to 0-based

            total_positions += 1;

            // Samples are listed after the count
            for sample in parts.iter().skip(3) {
                sample_snps
                    .entry(sample.to_string())
                    .or_default()
                    .push(pos);
            }
        }
    }

    log::info!("Parsed {} SNP positions from CFSAN snplist", total_positions);
    Ok((total_positions, sample_snps))
}

/// Compare coreguard gaps with CFSAN results
pub fn compare_with_cfsan(
    coreguard_gaps: &[(String, Vec<Region>)],
    cfsan_metrics: &[CfsanSampleMetrics],
    cfsan_snps: &HashMap<String, Vec<usize>>,
    total_cfsan_snps: usize,
    reference_length: usize,
) -> CfsanComparison {
    let mut sample_comparisons = Vec::new();

    // Build metrics lookup
    let metrics_map: HashMap<&str, &CfsanSampleMetrics> = cfsan_metrics
        .iter()
        .map(|m| (m.sample.as_str(), m))
        .collect();

    for (sample_name, cg_regions) in coreguard_gaps {
        // Get CFSAN SNPs for this sample
        let snp_positions: Vec<usize> = cfsan_snps
            .get(sample_name)
            .cloned()
            .unwrap_or_default();

        if snp_positions.is_empty() {
            log::debug!("No CFSAN SNPs for sample: {}", sample_name);
            continue;
        }

        // Calculate coreguard gap positions
        let cg_gap_positions: HashSet<usize> = cg_regions
            .iter()
            .flat_map(|r| r.start..r.end)
            .collect();
        let cg_gap_bases = cg_gap_positions.len();

        // Count SNPs in gaps
        let snps_in_gaps = snp_positions
            .iter()
            .filter(|pos| cg_gap_positions.contains(pos))
            .count();
        let snps_retained = snp_positions.len() - snps_in_gaps;
        let snp_retention_pct = if !snp_positions.is_empty() {
            (snps_retained as f64 / snp_positions.len() as f64) * 100.0
        } else {
            100.0
        };

        // Check overlap with CFSAN missing positions (if metrics available)
        let missing_overlap = metrics_map
            .get(sample_name.as_str())
            .map(|m| {
                // Approximate: assume missing positions correlate with gaps
                std::cmp::min(m.missing_positions, cg_gap_bases)
            })
            .unwrap_or(0);

        sample_comparisons.push(CfsanSampleComparison {
            sample: sample_name.clone(),
            cg_gap_bases,
            cfsan_snps: snp_positions.len(),
            snps_in_gaps,
            snp_retention_pct,
            missing_overlap,
        });
    }

    // Summary
    let total_samples = sample_comparisons.len();
    let avg_snp_retention_pct = if total_samples > 0 {
        sample_comparisons.iter().map(|s| s.snp_retention_pct).sum::<f64>() / total_samples as f64
    } else {
        100.0
    };
    let total_snps: usize = sample_comparisons.iter().map(|s| s.cfsan_snps).sum();
    let total_snps_in_gaps: usize = sample_comparisons.iter().map(|s| s.snps_in_gaps).sum();

    let interpretation = if avg_snp_retention_pct > 99.0 {
        "Excellent: coreguard gaps do not affect CFSAN variant calling".to_string()
    } else if avg_snp_retention_pct > 95.0 {
        "Good: most CFSAN SNPs are preserved outside coreguard gap regions".to_string()
    } else if avg_snp_retention_pct > 90.0 {
        "Moderate: some CFSAN SNPs fall in coreguard gaps - review thresholds".to_string()
    } else {
        "Warning: significant SNP loss in coreguard gaps - consider adjusting parameters".to_string()
    };

    CfsanComparison {
        reference_length,
        cfsan_metrics: cfsan_metrics.to_vec(),
        total_cfsan_snps,
        sample_comparisons,
        summary: CfsanSummary {
            total_samples,
            avg_snp_retention_pct,
            total_snps,
            total_snps_in_gaps,
            interpretation,
        },
        distance_matrix: None, // Will be set by caller if available
    }
}

// ============================================================================
// Original comparison (for backwards compatibility with --snippy-dir)
// ============================================================================

/// Comparison result for a single sample
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleComparison {
    /// Sample name
    pub sample: String,
    /// Total coreguard gap bases
    pub coreguard_gap_bases: usize,
    /// Total SNP pipeline filtered bases
    pub snp_filtered_bases: usize,
    /// Bases filtered by both (concordant)
    pub concordant_bases: usize,
    /// Bases filtered only by coreguard (CG predicts gap, SNP pipeline keeps)
    pub coreguard_only_bases: usize,
    /// Bases filtered only by SNP pipeline (SNP filters, CG keeps)
    pub snp_only_bases: usize,
    /// Concordance percentage (concordant / union)
    pub concordance_pct: f64,
    /// Regions unique to coreguard
    pub coreguard_only_regions: Vec<Region>,
    /// Regions unique to SNP pipeline
    pub snp_only_regions: Vec<Region>,
}

/// Complete SNP comparison analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnpComparison {
    /// Source pipeline name (e.g., "Snippy", "CFSAN")
    pub pipeline_name: String,
    /// Reference length
    pub reference_length: usize,
    /// Per-sample comparisons
    pub samples: Vec<SampleComparison>,
    /// Summary statistics
    pub summary: ComparisonSummary,
}

/// Summary of comparison across all samples
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComparisonSummary {
    /// Total samples compared
    pub total_samples: usize,
    /// Average concordance percentage
    pub avg_concordance_pct: f64,
    /// Total concordant bases (all samples)
    pub total_concordant_bases: usize,
    /// Total coreguard-only bases
    pub total_coreguard_only_bases: usize,
    /// Total SNP-only bases
    pub total_snp_only_bases: usize,
    /// Interpretation
    pub interpretation: String,
}

/// Parse Snippy aligned.fa file and extract N (masked) positions
pub fn parse_snippy_aligned(path: &Path) -> Result<HashSet<usize>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open Snippy aligned file: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut sequence = String::new();
    for line in reader.lines() {
        let line = line?;
        if !line.starts_with('>') {
            sequence.push_str(line.trim());
        }
    }

    // Find N positions
    let n_positions: HashSet<usize> = sequence
        .char_indices()
        .filter(|(_, c)| *c == 'N' || *c == 'n')
        .map(|(i, _)| i)
        .collect();

    Ok(n_positions)
}

/// Parse Snippy directory and extract masked positions for the sample
pub fn parse_snippy_dir(dir: &Path, sample_name: &str) -> Result<HashSet<usize>> {
    // Try different naming patterns
    let patterns = [
        format!("{}_snippy.aligned.fa", sample_name),
        format!("{}.aligned.fa", sample_name),
        "snps.aligned.fa".to_string(),
    ];

    for pattern in &patterns {
        let aligned_path = dir.join(pattern);
        if aligned_path.exists() {
            return parse_snippy_aligned(&aligned_path);
        }
    }

    // If directory contains the aligned.fa directly
    let direct_path = dir.join(format!("{}_snippy.aligned.fa", sample_name));
    if direct_path.exists() {
        return parse_snippy_aligned(&direct_path);
    }

    anyhow::bail!(
        "Could not find Snippy aligned.fa in {} for sample {}",
        dir.display(),
        sample_name
    );
}

/// Parse CFSAN snplist.txt to understand which positions are included
/// Returns positions that are NOT in the SNP list (i.e., filtered out)
pub fn parse_cfsan_snplist(path: &Path, reference_length: usize) -> Result<HashMap<String, HashSet<usize>>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open CFSAN snplist: {}", path.display()))?;
    let reader = BufReader::new(file);

    // Track which samples have which positions
    let mut sample_positions: HashMap<String, HashSet<usize>> = HashMap::new();
    let mut all_samples: HashSet<String> = HashSet::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 3 {
            // Format: chrom, position, num_samples, sample1, sample2, ...
            let position: usize = parts[1].parse().unwrap_or(0);
            if position > 0 {
                let position = position - 1; // Convert to 0-based
                for sample in parts.iter().skip(3) {
                    all_samples.insert(sample.to_string());
                    sample_positions
                        .entry(sample.to_string())
                        .or_default()
                        .insert(position);
                }
            }
        }
    }

    // For each sample, positions NOT in snplist are "filtered"
    // But this isn't quite right - snplist only has variant positions
    // We need a different approach for CFSAN

    Ok(sample_positions)
}

/// Parse BAM/SAM file to extract low-coverage positions
/// This is a more generic approach that works with any aligner
pub fn parse_bam_coverage(
    _path: &Path,
    _min_depth: usize,
    _reference_length: usize,
) -> Result<HashSet<usize>> {
    // This would require a BAM parser like rust-htslib
    // For now, we'll use samtools via command line
    anyhow::bail!("BAM parsing not yet implemented - use Snippy aligned.fa instead");
}

/// Convert a set of positions to Region objects
fn positions_to_regions(positions: &HashSet<usize>, max_gap: usize) -> Vec<Region> {
    if positions.is_empty() {
        return vec![];
    }

    let mut sorted: Vec<usize> = positions.iter().copied().collect();
    sorted.sort_unstable();

    let mut regions = Vec::new();
    let mut start = sorted[0];
    let mut end = sorted[0] + 1;

    for &pos in sorted.iter().skip(1) {
        if pos <= end + max_gap {
            end = pos + 1;
        } else {
            regions.push(Region::new(start, end));
            start = pos;
            end = pos + 1;
        }
    }
    regions.push(Region::new(start, end));

    regions
}

/// Compare coreguard gaps with SNP pipeline filtered positions
pub fn compare_with_snp_pipeline(
    coreguard_gaps: &[(String, Vec<Region>)],
    snp_filtered: &HashMap<String, HashSet<usize>>,
    pipeline_name: &str,
    reference_length: usize,
) -> SnpComparison {
    let mut samples = Vec::new();
    let mut total_concordant = 0usize;
    let mut total_cg_only = 0usize;
    let mut total_snp_only = 0usize;

    for (sample_name, cg_regions) in coreguard_gaps {
        // Skip samples without SNP pipeline results
        let snp_positions = match snp_filtered.get(sample_name) {
            Some(pos) => pos.clone(),
            None => {
                log::debug!("Skipping {} - no SNP pipeline results", sample_name);
                continue;
            }
        };

        // Convert coreguard regions to position set
        let cg_positions: HashSet<usize> = cg_regions
            .iter()
            .flat_map(|r| r.start..r.end)
            .collect();

        // Calculate overlaps
        let concordant: HashSet<usize> = cg_positions.intersection(&snp_positions).copied().collect();
        let cg_only: HashSet<usize> = cg_positions.difference(&snp_positions).copied().collect();
        let snp_only: HashSet<usize> = snp_positions.difference(&cg_positions).copied().collect();

        let concordant_bases = concordant.len();
        let cg_only_bases = cg_only.len();
        let snp_only_bases = snp_only.len();

        // Calculate concordance
        let union_size = cg_positions.len() + snp_positions.len() - concordant_bases;
        let concordance_pct = if union_size > 0 {
            (concordant_bases as f64 / union_size as f64) * 100.0
        } else {
            100.0
        };

        // Convert to regions (merge nearby positions)
        let cg_only_regions = positions_to_regions(&cg_only, 10);
        let snp_only_regions = positions_to_regions(&snp_only, 10);

        total_concordant += concordant_bases;
        total_cg_only += cg_only_bases;
        total_snp_only += snp_only_bases;

        samples.push(SampleComparison {
            sample: sample_name.clone(),
            coreguard_gap_bases: cg_positions.len(),
            snp_filtered_bases: snp_positions.len(),
            concordant_bases,
            coreguard_only_bases: cg_only_bases,
            snp_only_bases,
            concordance_pct,
            coreguard_only_regions: cg_only_regions,
            snp_only_regions,
        });
    }

    // Calculate summary
    let avg_concordance = if !samples.is_empty() {
        samples.iter().map(|s| s.concordance_pct).sum::<f64>() / samples.len() as f64
    } else {
        0.0
    };

    let interpretation = if avg_concordance > 80.0 {
        "High concordance: coreguard predictions align well with SNP pipeline filtering".to_string()
    } else if avg_concordance > 50.0 {
        "Moderate concordance: some differences in filtering approach".to_string()
    } else if total_cg_only > total_snp_only {
        "Low concordance: coreguard is more conservative (flags more regions)".to_string()
    } else {
        "Low concordance: SNP pipeline applies additional quality filters".to_string()
    };

    let num_samples = samples.len();

    SnpComparison {
        pipeline_name: pipeline_name.to_string(),
        reference_length,
        samples,
        summary: ComparisonSummary {
            total_samples: num_samples,
            avg_concordance_pct: avg_concordance,
            total_concordant_bases: total_concordant,
            total_coreguard_only_bases: total_cg_only,
            total_snp_only_bases: total_snp_only,
            interpretation,
        },
    }
}

/// Load Snippy results from multiple sample directories
/// If reference_name is provided, prioritize directories ending with _referenceref
pub fn load_snippy_results(
    results_dir: &Path,
    sample_names: &[String],
    reference_name: Option<&str>,
) -> Result<HashMap<String, HashSet<usize>>> {
    let mut results = HashMap::new();

    // Priority: 1) sample_snippy_refnameref, 2) sample_snippy, 3) sample
    // Store all candidates, prioritize by type
    let mut sample_dirs_priority1: HashMap<String, std::path::PathBuf> = HashMap::new();
    let mut sample_dirs_priority2: HashMap<String, std::path::PathBuf> = HashMap::new();
    let mut sample_dirs_priority3: HashMap<String, std::path::PathBuf> = HashMap::new();

    // Build expected suffix for reference-specific directories (e.g., "TE15676ref")
    let ref_suffix = reference_name.map(|r| format!("{}ref", r));

    if let Ok(entries) = std::fs::read_dir(results_dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                let dir_name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");

                for sample in sample_names {
                    let prefix = format!("{}_snippy_", sample);

                    // Priority 1: sample_snippy_refnameref (exact reference match)
                    if let Some(ref suffix) = ref_suffix {
                        if dir_name == format!("{}{}", prefix, suffix) {
                            sample_dirs_priority1.insert(sample.clone(), path.clone());
                            continue;
                        }
                    }

                    // Priority 2: sample_snippy
                    if dir_name == format!("{}_snippy", sample) {
                        sample_dirs_priority2.entry(sample.clone()).or_insert_with(|| path.clone());
                    }
                    // Priority 3: exact sample name
                    else if dir_name == sample {
                        sample_dirs_priority3.entry(sample.clone()).or_insert_with(|| path.clone());
                    }
                }
            }
        }
    }

    // Now load results from the discovered directories, preferring higher priority
    for sample in sample_names {
        let sample_dir = sample_dirs_priority1.get(sample)
            .or_else(|| sample_dirs_priority2.get(sample))
            .or_else(|| sample_dirs_priority3.get(sample));

        if let Some(sample_dir) = sample_dir {
            match parse_snippy_dir(sample_dir, sample) {
                Ok(positions) => {
                    log::info!("Loaded Snippy results for {} from {}: {} masked positions",
                        sample, sample_dir.display(), positions.len());
                    results.insert(sample.clone(), positions);
                }
                Err(e) => {
                    log::debug!("Could not load Snippy from {}: {}", sample_dir.display(), e);
                }
            }
        }

        if !results.contains_key(sample) {
            log::warn!("Could not find Snippy results for sample: {}", sample);
        }
    }

    Ok(results)
}

// ============================================================================
// SNP Distance Matrix
// ============================================================================

/// SNP distance matrix (pairwise distances between samples)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnpDistanceMatrix {
    /// Sample names (row/column headers)
    pub samples: Vec<String>,
    /// Distance values as 2D array (row-major)
    pub distances: Vec<Vec<i64>>,
}

/// Parse a TSV distance matrix file
/// Format: header row with sample names, then rows with sample name + distances
pub fn parse_distance_matrix(path: &Path) -> Result<SnpDistanceMatrix> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open distance matrix: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut lines = reader.lines();

    // Parse header to get sample names
    let header = lines.next()
        .ok_or_else(|| anyhow::anyhow!("Empty distance matrix file"))??;

    let samples: Vec<String> = header.split('\t')
        .skip(1) // Skip empty first column
        .map(|s| s.trim().trim_end_matches("_snippy").to_string())
        .filter(|s| !s.is_empty())
        .collect();

    if samples.is_empty() {
        anyhow::bail!("No sample names found in distance matrix header");
    }

    // Parse distance rows
    let mut distances: Vec<Vec<i64>> = Vec::new();

    for line in lines {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue;
        }

        // Skip first column (sample name), parse distances
        let row: Vec<i64> = parts[1..]
            .iter()
            .take(samples.len())
            .map(|s| s.trim().parse().unwrap_or(0))
            .collect();

        if row.len() == samples.len() {
            distances.push(row);
        }
    }

    Ok(SnpDistanceMatrix { samples, distances })
}

/// Try to find and parse a distance matrix from a directory
pub fn find_distance_matrix(dir: &Path, patterns: &[&str]) -> Option<SnpDistanceMatrix> {
    for pattern in patterns {
        let path = dir.join(pattern);
        if path.exists() {
            match parse_distance_matrix(&path) {
                Ok(matrix) => {
                    log::info!("Loaded distance matrix from: {}", path.display());
                    return Some(matrix);
                }
                Err(e) => {
                    log::warn!("Failed to parse {}: {}", path.display(), e);
                }
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_positions_to_regions() {
        let positions: HashSet<usize> = [1, 2, 3, 10, 11, 12, 20].into_iter().collect();
        let regions = positions_to_regions(&positions, 5);

        // Gap between end of [1,2,3] (exclusive end=4) and 10 is 6, exceeds max_gap=5
        assert_eq!(regions.len(), 3);
        assert_eq!(regions[0].start, 1);
        assert_eq!(regions[0].end, 4);
        assert_eq!(regions[1].start, 10);
        assert_eq!(regions[1].end, 13);
        assert_eq!(regions[2].start, 20);
        assert_eq!(regions[2].end, 21);

        // With max_gap=6, [1-3] and [10-12] merge (10 <= 4+6) but not [20] (20 > 13+6)
        let regions = positions_to_regions(&positions, 6);
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].start, 1);
        assert_eq!(regions[0].end, 13);
        assert_eq!(regions[1].start, 20);
        assert_eq!(regions[1].end, 21);
    }

    #[test]
    fn test_compare_empty() {
        let cg_gaps: Vec<(String, Vec<Region>)> = vec![];
        let snp_filtered: HashMap<String, HashSet<usize>> = HashMap::new();

        let result = compare_with_snp_pipeline(&cg_gaps, &snp_filtered, "Test", 1000);
        assert_eq!(result.samples.len(), 0);
    }
}
