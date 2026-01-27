//! Sample clustering and outlier detection
//!
//! Analyzes reference alignment to detect:
//! - Clonal groups (samples with high identity to reference)
//! - Outliers (samples significantly different from reference)
//! - Gap quality patterns using pairwise gap analysis

use crate::pairwise::PairwiseResults;
use crate::reference::ReferenceResults;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Cluster analysis results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClusterAnalysis {
    /// Detected clusters of similar samples (based on reference coverage)
    pub clusters: Vec<SampleCluster>,
    /// Samples identified as outliers
    pub outliers: Vec<OutlierSample>,
    /// Overall pattern detected
    pub pattern: ClusterPattern,
    /// Warnings for SNP analysis
    pub warnings: Vec<String>,
    /// Recommendations
    pub recommendations: Vec<String>,
    /// Gap quality statistics
    pub gap_quality_stats: GapQualityStats,
}

/// A cluster of similar samples
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleCluster {
    /// Cluster ID
    pub id: usize,
    /// Samples in this cluster
    pub samples: Vec<String>,
    /// Minimum reference coverage within cluster
    pub min_coverage: f64,
    /// Average reference coverage within cluster
    pub avg_coverage: f64,
    /// Average pairwise gap quality within cluster
    pub avg_gap_quality: f64,
    /// Cluster type
    pub cluster_type: ClusterType,
}

/// Type of cluster
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum ClusterType {
    /// Nearly identical coverage (>99% ref coverage) - likely same strain
    HighCoverage,
    /// Good coverage (>95% ref coverage)
    GoodCoverage,
    /// Moderate coverage (>90% ref coverage)
    ModerateCoverage,
    /// Low coverage (<90% ref coverage)
    LowCoverage,
}

/// An outlier sample
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutlierSample {
    /// Sample name
    pub sample: String,
    /// Reference coverage
    pub reference_coverage: f64,
    /// Average gap quality with other samples
    pub avg_gap_quality: f64,
    /// Severity of outlier status
    pub severity: OutlierSeverity,
    /// Reason for being flagged
    pub reason: String,
}

/// Severity of outlier
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum OutlierSeverity {
    /// Minor outlier - slightly different
    Minor,
    /// Moderate outlier - noticeably different
    Moderate,
    /// Major outlier - very different, may cause issues
    Major,
}

/// Overall pattern in the dataset
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum ClusterPattern {
    /// All samples have high reference coverage (likely outbreak)
    AllHighCoverage,
    /// High coverage samples plus outliers
    HighCoverageWithOutliers,
    /// Variable coverage across samples
    VariableCoverage,
    /// High diversity in gap patterns
    HighGapDiversity,
    /// Single sample or insufficient data
    Insufficient,
}

impl std::fmt::Display for ClusterPattern {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ClusterPattern::AllHighCoverage => write!(f, "ALL_HIGH_COVERAGE"),
            ClusterPattern::HighCoverageWithOutliers => write!(f, "HIGH_COVERAGE_WITH_OUTLIERS"),
            ClusterPattern::VariableCoverage => write!(f, "VARIABLE_COVERAGE"),
            ClusterPattern::HighGapDiversity => write!(f, "HIGH_GAP_DIVERSITY"),
            ClusterPattern::Insufficient => write!(f, "INSUFFICIENT_DATA"),
        }
    }
}

/// Gap quality statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GapQualityStats {
    /// Average pairwise gap quality score
    pub mean_quality: f64,
    /// Minimum pairwise gap quality score
    pub min_quality: f64,
    /// Maximum pairwise gap quality score
    pub max_quality: f64,
    /// Dataset-level quality score
    pub dataset_quality: f64,
    /// Number of high-quality pairs (quality > 0.8)
    pub high_quality_pairs: usize,
    /// Number of low-quality pairs (quality < 0.2)
    pub low_quality_pairs: usize,
}

/// Perform cluster analysis on samples
pub fn analyze_clusters(
    sample_names: &[String],
    pairwise: &PairwiseResults,
    reference: &ReferenceResults,
) -> ClusterAnalysis {
    if sample_names.len() < 2 {
        return ClusterAnalysis {
            clusters: vec![],
            outliers: vec![],
            pattern: ClusterPattern::Insufficient,
            warnings: vec!["Insufficient samples for cluster analysis".to_string()],
            recommendations: vec![],
            gap_quality_stats: GapQualityStats {
                mean_quality: 1.0,
                min_quality: 1.0,
                max_quality: 1.0,
                dataset_quality: 1.0,
                high_quality_pairs: 0,
                low_quality_pairs: 0,
            },
        };
    }

    // Calculate gap quality statistics from pairwise results
    let gap_quality_stats = calculate_gap_quality_stats(pairwise);

    // Build reference coverage map
    let coverage_map: HashMap<String, f64> = reference.alignments
        .iter()
        .map(|a| (a.sample_name.clone(), a.reference_coverage))
        .collect();

    // Build per-sample gap quality map
    let mut sample_gap_quality: HashMap<String, f64> = HashMap::new();
    for name in sample_names {
        sample_gap_quality.insert(name.clone(), pairwise.avg_quality(name));
    }

    // Detect outliers based on coverage and gap quality
    let outliers = detect_outliers(sample_names, &coverage_map, &sample_gap_quality, &gap_quality_stats);

    // Detect clusters based on reference coverage
    let clusters = detect_clusters(sample_names, &coverage_map, &sample_gap_quality);

    // Determine overall pattern
    let pattern = determine_pattern(&clusters, &outliers, &coverage_map, &gap_quality_stats);

    // Generate warnings and recommendations
    let (warnings, recommendations) = generate_warnings_and_recommendations(
        &pattern,
        &clusters,
        &outliers,
        &gap_quality_stats,
        reference,
    );

    ClusterAnalysis {
        clusters,
        outliers,
        pattern,
        warnings,
        recommendations,
        gap_quality_stats,
    }
}

fn calculate_gap_quality_stats(pairwise: &PairwiseResults) -> GapQualityStats {
    if pairwise.results.is_empty() {
        return GapQualityStats {
            mean_quality: 1.0,
            min_quality: 1.0,
            max_quality: 1.0,
            dataset_quality: 1.0,
            high_quality_pairs: 0,
            low_quality_pairs: 0,
        };
    }

    let qualities: Vec<f64> = pairwise.results.iter().map(|r| r.quality_score).collect();

    let mean_quality = qualities.iter().sum::<f64>() / qualities.len() as f64;
    let min_quality = qualities.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_quality = qualities.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    let high_quality_pairs = qualities.iter().filter(|&&q| q > 0.8).count();
    let low_quality_pairs = qualities.iter().filter(|&&q| q < 0.2).count();

    GapQualityStats {
        mean_quality,
        min_quality,
        max_quality,
        dataset_quality: pairwise.dataset_quality(),
        high_quality_pairs,
        low_quality_pairs,
    }
}

fn detect_outliers(
    sample_names: &[String],
    coverage_map: &HashMap<String, f64>,
    gap_quality_map: &HashMap<String, f64>,
    stats: &GapQualityStats,
) -> Vec<OutlierSample> {
    let mut outliers = Vec::new();

    for name in sample_names {
        let coverage = *coverage_map.get(name).unwrap_or(&0.0);
        let gap_quality = *gap_quality_map.get(name).unwrap_or(&1.0);

        // Determine if outlier based on coverage and gap quality
        let (severity, reason) = if coverage < 0.90 {
            (Some(OutlierSeverity::Major), format!("Low reference coverage: {:.1}%", coverage * 100.0))
        } else if coverage < 0.95 {
            (Some(OutlierSeverity::Moderate), format!("Moderate reference coverage: {:.1}%", coverage * 100.0))
        } else if gap_quality < 0.2 && stats.mean_quality > 0.5 {
            (Some(OutlierSeverity::Major), format!("Low gap quality: {:.2} (mean: {:.2})", gap_quality, stats.mean_quality))
        } else if gap_quality < stats.mean_quality - 0.3 {
            (Some(OutlierSeverity::Moderate), format!("Below average gap quality: {:.2}", gap_quality))
        } else {
            (None, String::new())
        };

        if let Some(sev) = severity {
            outliers.push(OutlierSample {
                sample: name.clone(),
                reference_coverage: coverage,
                avg_gap_quality: gap_quality,
                severity: sev,
                reason,
            });
        }
    }

    outliers
}

fn detect_clusters(
    sample_names: &[String],
    coverage_map: &HashMap<String, f64>,
    gap_quality_map: &HashMap<String, f64>,
) -> Vec<SampleCluster> {
    // Simple clustering by coverage threshold
    let mut high_coverage: Vec<String> = Vec::new();
    let mut good_coverage: Vec<String> = Vec::new();
    let mut moderate_coverage: Vec<String> = Vec::new();
    let mut low_coverage: Vec<String> = Vec::new();

    for name in sample_names {
        let coverage = *coverage_map.get(name).unwrap_or(&0.0);
        if coverage >= 0.99 {
            high_coverage.push(name.clone());
        } else if coverage >= 0.95 {
            good_coverage.push(name.clone());
        } else if coverage >= 0.90 {
            moderate_coverage.push(name.clone());
        } else {
            low_coverage.push(name.clone());
        }
    }

    let mut clusters = Vec::new();
    let mut cluster_id = 0;

    // Create clusters for each coverage tier (only if > 1 sample)
    for (samples, cluster_type) in [
        (high_coverage, ClusterType::HighCoverage),
        (good_coverage, ClusterType::GoodCoverage),
        (moderate_coverage, ClusterType::ModerateCoverage),
        (low_coverage, ClusterType::LowCoverage),
    ] {
        if samples.len() > 1 {
            let coverages: Vec<f64> = samples.iter()
                .map(|s| *coverage_map.get(s).unwrap_or(&0.0))
                .collect();
            let gap_qualities: Vec<f64> = samples.iter()
                .map(|s| *gap_quality_map.get(s).unwrap_or(&1.0))
                .collect();

            clusters.push(SampleCluster {
                id: cluster_id,
                samples,
                min_coverage: coverages.iter().cloned().fold(f64::INFINITY, f64::min),
                avg_coverage: coverages.iter().sum::<f64>() / coverages.len() as f64,
                avg_gap_quality: gap_qualities.iter().sum::<f64>() / gap_qualities.len() as f64,
                cluster_type,
            });
            cluster_id += 1;
        }
    }

    clusters
}

fn determine_pattern(
    clusters: &[SampleCluster],
    outliers: &[OutlierSample],
    coverage_map: &HashMap<String, f64>,
    gap_stats: &GapQualityStats,
) -> ClusterPattern {
    let total_samples = coverage_map.len();

    if total_samples < 2 {
        return ClusterPattern::Insufficient;
    }

    let high_coverage_count = coverage_map.values().filter(|&&c| c >= 0.95).count();
    let major_outliers = outliers.iter().filter(|o| o.severity == OutlierSeverity::Major).count();

    if high_coverage_count == total_samples && outliers.is_empty() {
        ClusterPattern::AllHighCoverage
    } else if high_coverage_count > 0 && major_outliers > 0 {
        ClusterPattern::HighCoverageWithOutliers
    } else if gap_stats.low_quality_pairs > gap_stats.high_quality_pairs {
        ClusterPattern::HighGapDiversity
    } else {
        ClusterPattern::VariableCoverage
    }
}

fn generate_warnings_and_recommendations(
    pattern: &ClusterPattern,
    clusters: &[SampleCluster],
    outliers: &[OutlierSample],
    gap_stats: &GapQualityStats,
    reference: &ReferenceResults,
) -> (Vec<String>, Vec<String>) {
    let mut warnings = Vec::new();
    let mut recommendations = Vec::new();

    // Reference-based warnings
    for alignment in &reference.alignments {
        if alignment.reference_coverage < 0.95 {
            warnings.push(format!(
                "REFERENCE: {} covers only {:.1}% of reference - missing regions",
                alignment.sample_name,
                alignment.reference_coverage * 100.0
            ));
        }
    }

    // Gap quality warnings
    if gap_stats.dataset_quality < 0.5 {
        warnings.push(format!(
            "LOW GAP QUALITY: Dataset quality score {:.2} - samples have very different gap patterns",
            gap_stats.dataset_quality
        ));
        recommendations.push("Consider removing samples with unique gaps to increase core genome".to_string());
    }

    if gap_stats.low_quality_pairs > 0 {
        warnings.push(format!(
            "GAP DIVERSITY: {} sample pairs have low gap overlap (quality < 0.2)",
            gap_stats.low_quality_pairs
        ));
    }

    // Outlier-specific warnings
    for outlier in outliers {
        if outlier.severity == OutlierSeverity::Major {
            warnings.push(format!(
                "MAJOR OUTLIER: {} - {}",
                outlier.sample,
                outlier.reason
            ));
            recommendations.push(format!(
                "Consider excluding {} from analysis",
                outlier.sample
            ));
        }
    }

    // Pattern-specific recommendations
    match pattern {
        ClusterPattern::AllHighCoverage => {
            warnings.push("All samples have high reference coverage - OUTBREAK pattern".to_string());
            recommendations.push("Use CLOSE reference (>99% ANI) to avoid alignment artifacts".to_string());
            recommendations.push("Strict variant callers may miss real SNPs - use sensitive settings".to_string());
            recommendations.push("Verify any 'polymorphic' SNPs with BAM-level validation".to_string());
        }
        ClusterPattern::HighCoverageWithOutliers => {
            let outlier_names: Vec<String> = outliers.iter().map(|o| o.sample.clone()).collect();
            warnings.push(format!(
                "HIGH COVERAGE + OUTLIERS: Outlier(s) [{}] may reduce core genome significantly",
                outlier_names.join(", ")
            ));
            recommendations.push("Consider analyzing high-coverage samples separately".to_string());
        }
        ClusterPattern::HighGapDiversity => {
            warnings.push(format!(
                "HIGH GAP DIVERSITY: Gap quality range {:.2}-{:.2}",
                gap_stats.min_quality,
                gap_stats.max_quality
            ));
            recommendations.push("Large core reduction expected - verify this is intended".to_string());
            recommendations.push("Consider using soft-core (95%) instead of strict core (100%)".to_string());
        }
        ClusterPattern::VariableCoverage => {
            recommendations.push("Variable coverage - check assembly/sequencing quality".to_string());
        }
        ClusterPattern::Insufficient => {
            warnings.push("Insufficient samples for meaningful cluster analysis".to_string());
        }
    }

    (warnings, recommendations)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gap_quality_stats() {
        let stats = GapQualityStats {
            mean_quality: 0.5,
            min_quality: 0.1,
            max_quality: 0.9,
            dataset_quality: 0.5,
            high_quality_pairs: 3,
            low_quality_pairs: 2,
        };

        assert!(stats.mean_quality > 0.0);
        assert!(stats.min_quality < stats.max_quality);
    }
}
