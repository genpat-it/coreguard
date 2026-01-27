//! Empty Column Analysis
//!
//! Analyzes which reference positions become "empty" (insufficient coverage)
//! at different sample inclusion thresholds.

use crate::gaps::Region;
use serde::{Deserialize, Serialize};

/// Coverage statistics for a reference position
#[derive(Debug, Clone)]
struct PositionCoverage {
    /// Number of samples that cover this position
    covered_by: usize,
}

/// Threshold analysis result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThresholdResult {
    /// Threshold percentage (0-100)
    pub threshold_pct: f64,
    /// Minimum samples required (computed from threshold)
    pub min_samples: usize,
    /// Number of positions in core at this threshold
    pub core_positions: usize,
    /// Core size as percentage of reference
    pub core_pct: f64,
    /// Number of empty positions at this threshold
    pub empty_positions: usize,
    /// Empty regions at this threshold
    pub empty_regions: Vec<Region>,
}

/// Per-position coverage histogram entry
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoverageHistogramEntry {
    /// Number of samples covering
    pub sample_count: usize,
    /// Number of positions with this coverage
    pub position_count: usize,
    /// Percentage of reference
    pub position_pct: f64,
}

/// Complete empty column analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EmptyColumnAnalysis {
    /// Reference length
    pub reference_length: usize,
    /// Total number of samples
    pub total_samples: usize,
    /// Results at different thresholds
    pub threshold_results: Vec<ThresholdResult>,
    /// Coverage histogram (how many positions have X samples covering)
    pub coverage_histogram: Vec<CoverageHistogramEntry>,
    /// Fragile regions (covered by <50% of samples)
    pub fragile_regions: Vec<FragileRegion>,
}

/// A region that's fragile (low coverage)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FragileRegion {
    /// Region coordinates
    pub region: Region,
    /// Maximum coverage in this region (samples)
    pub max_coverage: usize,
    /// Maximum coverage percentage
    pub max_coverage_pct: f64,
    /// Samples that DON'T cover this region
    pub missing_samples: Vec<String>,
}

/// Analyze empty columns at various thresholds
pub fn analyze_empty_columns(
    reference_length: usize,
    sample_gaps: &[(String, Vec<Region>)],
    thresholds: &[f64],
) -> EmptyColumnAnalysis {
    let total_samples = sample_gaps.len();

    if total_samples == 0 || reference_length == 0 {
        return EmptyColumnAnalysis {
            reference_length,
            total_samples,
            threshold_results: vec![],
            coverage_histogram: vec![],
            fragile_regions: vec![],
        };
    }

    // Build coverage array: for each position, count how many samples cover it
    // We use a compact representation: track gap intervals rather than per-position
    let mut coverage: Vec<usize> = vec![total_samples; reference_length];

    // Subtract 1 for each sample that has a gap at each position
    for (_sample_name, gaps) in sample_gaps {
        for gap in gaps {
            let start = gap.start.min(reference_length);
            let end = gap.end.min(reference_length);
            for pos in start..end {
                coverage[pos] = coverage[pos].saturating_sub(1);
            }
        }
    }

    // Build coverage histogram
    let mut histogram_counts: Vec<usize> = vec![0; total_samples + 1];
    for &cov in &coverage {
        histogram_counts[cov] += 1;
    }

    let coverage_histogram: Vec<CoverageHistogramEntry> = histogram_counts
        .iter()
        .enumerate()
        .map(|(count, &positions)| CoverageHistogramEntry {
            sample_count: count,
            position_count: positions,
            position_pct: (positions as f64 / reference_length as f64) * 100.0,
        })
        .collect();

    // Analyze at each threshold
    let threshold_results: Vec<ThresholdResult> = thresholds
        .iter()
        .map(|&threshold_pct| {
            let min_samples = ((threshold_pct / 100.0) * total_samples as f64).ceil() as usize;
            let min_samples = min_samples.max(1); // At least 1 sample required

            // Count core positions (positions with >= min_samples coverage)
            let core_positions = coverage.iter().filter(|&&c| c >= min_samples).count();
            let empty_positions = reference_length - core_positions;

            // Find empty regions at this threshold
            let empty_regions = find_empty_regions(&coverage, min_samples, reference_length);

            ThresholdResult {
                threshold_pct,
                min_samples,
                core_positions,
                core_pct: (core_positions as f64 / reference_length as f64) * 100.0,
                empty_positions,
                empty_regions,
            }
        })
        .collect();

    // Find fragile regions (covered by <50% of samples)
    let fragile_threshold = (total_samples as f64 * 0.5).ceil() as usize;
    let fragile_regions = find_fragile_regions(&coverage, fragile_threshold, sample_gaps, reference_length);

    EmptyColumnAnalysis {
        reference_length,
        total_samples,
        threshold_results,
        coverage_histogram,
        fragile_regions,
    }
}

/// Find contiguous regions where coverage < min_samples
fn find_empty_regions(coverage: &[usize], min_samples: usize, reference_length: usize) -> Vec<Region> {
    let mut regions = Vec::new();
    let mut in_empty = false;
    let mut start = 0;

    for (pos, &cov) in coverage.iter().enumerate() {
        if cov < min_samples {
            if !in_empty {
                start = pos;
                in_empty = true;
            }
        } else if in_empty {
            regions.push(Region::new(start, pos));
            in_empty = false;
        }
    }

    // Close last region if still in empty
    if in_empty {
        regions.push(Region::new(start, reference_length));
    }

    regions
}

/// Find fragile regions with detailed sample information
fn find_fragile_regions(
    coverage: &[usize],
    fragile_threshold: usize,
    sample_gaps: &[(String, Vec<Region>)],
    reference_length: usize,
) -> Vec<FragileRegion> {
    // First find the raw fragile regions
    let raw_regions = find_empty_regions(coverage, fragile_threshold, reference_length);
    let total_samples = sample_gaps.len();

    // For each region, find which samples are missing
    raw_regions
        .into_iter()
        .filter(|r| r.len() >= 100) // Only report regions >= 100bp
        .map(|region| {
            // Find max coverage in this region
            let max_cov = coverage[region.start..region.end]
                .iter()
                .copied()
                .max()
                .unwrap_or(0);

            // Find samples that have gaps overlapping this region
            let missing_samples: Vec<String> = sample_gaps
                .iter()
                .filter(|(_, gaps)| {
                    gaps.iter().any(|g| g.start < region.end && g.end > region.start)
                })
                .map(|(name, _)| name.clone())
                .collect();

            FragileRegion {
                region,
                max_coverage: max_cov,
                max_coverage_pct: (max_cov as f64 / total_samples as f64) * 100.0,
                missing_samples,
            }
        })
        .collect()
}

/// Default thresholds to analyze
pub fn default_thresholds() -> Vec<f64> {
    vec![100.0, 95.0, 90.0, 80.0, 70.0, 50.0, 25.0]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_column_analysis() {
        // 4 samples, 100bp reference
        // Sample1: gap at 20-30
        // Sample2: gap at 25-35
        // Sample3: gap at 80-100
        // Sample4: no gaps
        let sample_gaps = vec![
            ("Sample1".to_string(), vec![Region::new(20, 30)]),
            ("Sample2".to_string(), vec![Region::new(25, 35)]),
            ("Sample3".to_string(), vec![Region::new(80, 100)]),
            ("Sample4".to_string(), vec![]),
        ];

        let analysis = analyze_empty_columns(100, &sample_gaps, &[100.0, 75.0, 50.0]);

        // At 100% threshold (all 4 samples must cover):
        // - Positions 0-20: 4/4 covered
        // - Positions 20-25: 3/4 (Sample1 has gap) - NOT in core
        // - Positions 25-30: 2/4 (Sample1,2 have gap) - NOT in core
        // - Positions 30-35: 3/4 (Sample2 has gap) - NOT in core
        // - Positions 35-80: 4/4 covered
        // - Positions 80-100: 3/4 (Sample3 has gap) - NOT in core
        // Core at 100%: 0-20 + 35-80 = 20 + 45 = 65 positions

        let t100 = &analysis.threshold_results[0];
        assert_eq!(t100.threshold_pct, 100.0);
        assert_eq!(t100.min_samples, 4);
        assert_eq!(t100.core_positions, 65);
        assert_eq!(t100.empty_positions, 35);
    }

    #[test]
    fn test_coverage_histogram() {
        let sample_gaps = vec![
            ("S1".to_string(), vec![Region::new(0, 50)]),
            ("S2".to_string(), vec![Region::new(0, 50)]),
            ("S3".to_string(), vec![]),
            ("S4".to_string(), vec![]),
        ];

        let analysis = analyze_empty_columns(100, &sample_gaps, &[100.0]);

        // Positions 0-50: covered by 2 samples
        // Positions 50-100: covered by 4 samples
        assert_eq!(analysis.coverage_histogram[2].position_count, 50); // 2 samples cover 50 positions
        assert_eq!(analysis.coverage_histogram[4].position_count, 50); // 4 samples cover 50 positions
    }
}
