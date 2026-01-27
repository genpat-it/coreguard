//! Sample quality metrics calculation

use crate::fasta::Sample;
use crate::pairwise::PairwiseResults;
use crate::reference::ReferenceResults;
use serde::{Deserialize, Serialize};

/// Quality metrics for a single sample
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleMetrics {
    /// Sample name
    pub name: String,

    // === Assembly statistics ===
    pub total_length: usize,
    pub num_contigs: usize,
    pub n50: usize,

    // === Reference alignment metrics ===
    /// Identity vs reference (alignment quality)
    pub ref_identity: f64,
    /// Coverage of sample aligned to reference
    pub ref_sample_coverage: f64,
    /// Coverage of reference by sample
    pub ref_reference_coverage: f64,
    /// Unaligned bases in sample (sample-specific sequence)
    pub ref_unaligned_bases: usize,
    /// Uncovered bases in reference (gaps)
    pub ref_uncovered_bases: usize,

    // === Pairwise gap metrics ===
    /// Average gap quality score with other samples
    pub avg_gap_quality: f64,
    /// Minimum gap quality (worst compatibility)
    pub min_gap_quality: f64,
    /// Total unique gap bases this sample introduces
    pub total_unique_gaps: usize,
    /// Average pairwise core size
    pub avg_pairwise_core: usize,

    // === Derived scores ===
    /// Z-score for reference coverage (negative = low coverage outlier)
    pub coverage_zscore: f64,
    /// Z-score for gap quality (negative = incompatible gaps)
    pub gap_quality_zscore: f64,
}

/// Calculate metrics for a single sample
pub fn calculate_metrics(
    sample: &Sample,
    pairwise: &PairwiseResults,
    reference: &ReferenceResults,
) -> SampleMetrics {
    // === Reference metrics ===
    let (ref_identity, ref_sample_coverage, ref_reference_coverage, ref_unaligned_bases, ref_uncovered_bases) =
        if let Some(alignment) = reference.get(&sample.name) {
            (
                alignment.identity,
                alignment.sample_coverage,
                alignment.reference_coverage,
                alignment.sample_unaligned,
                alignment.reference_uncovered,
            )
        } else {
            (0.0, 0.0, 0.0, sample.total_length, reference.reference_length)
        };

    // === Pairwise gap metrics ===
    let pair_results = pairwise.for_sample(&sample.name);

    let avg_gap_quality = pairwise.avg_quality(&sample.name);
    let min_gap_quality = pair_results
        .iter()
        .map(|r| r.quality_score)
        .fold(f64::INFINITY, f64::min);
    let min_gap_quality = if min_gap_quality.is_infinite() { 1.0 } else { min_gap_quality };

    let total_unique_gaps = pairwise.total_unique_gaps(&sample.name);
    let avg_pairwise_core = pairwise.avg_pairwise_core(&sample.name);

    SampleMetrics {
        name: sample.name.clone(),
        total_length: sample.total_length,
        num_contigs: sample.num_contigs,
        n50: sample.n50,
        ref_identity,
        ref_sample_coverage,
        ref_reference_coverage,
        ref_unaligned_bases,
        ref_uncovered_bases,
        avg_gap_quality,
        min_gap_quality,
        total_unique_gaps,
        avg_pairwise_core,
        coverage_zscore: 0.0,
        gap_quality_zscore: 0.0,
    }
}

/// Calculate z-scores across all samples
pub fn calculate_zscores(metrics: &mut [SampleMetrics]) {
    if metrics.is_empty() {
        return;
    }

    // Coverage z-scores
    let coverage_mean: f64 = metrics.iter().map(|m| m.ref_reference_coverage).sum::<f64>() / metrics.len() as f64;
    let coverage_var: f64 = metrics
        .iter()
        .map(|m| (m.ref_reference_coverage - coverage_mean).powi(2))
        .sum::<f64>()
        / metrics.len() as f64;
    let coverage_std = coverage_var.sqrt();

    // Gap quality z-scores
    let quality_mean: f64 = metrics.iter().map(|m| m.avg_gap_quality).sum::<f64>() / metrics.len() as f64;
    let quality_var: f64 = metrics
        .iter()
        .map(|m| (m.avg_gap_quality - quality_mean).powi(2))
        .sum::<f64>()
        / metrics.len() as f64;
    let quality_std = quality_var.sqrt();

    for m in metrics.iter_mut() {
        m.coverage_zscore = if coverage_std > 0.0 {
            (m.ref_reference_coverage - coverage_mean) / coverage_std
        } else {
            0.0
        };

        m.gap_quality_zscore = if quality_std > 0.0 {
            (m.avg_gap_quality - quality_mean) / quality_std
        } else {
            0.0
        };
    }
}
