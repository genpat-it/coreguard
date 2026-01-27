//! Core alignment impact prediction
//!
//! Predicts how much each sample would reduce the core alignment
//! if included in SNP analysis (due to unique gaps).
//!
//! Uses REFERENCE coordinates as the common coordinate system.
//! A position causes "gap contagion" if it's uncovered in ONE sample
//! but covered in ALL others - forcing SNP pipelines to exclude it.

use crate::gaps::Region;
use crate::pairwise::PairwiseResults;
use crate::reference::ReferenceResults;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

/// Impact prediction for a single sample
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleImpact {
    /// Sample name
    pub sample: String,

    /// Gaps unique to this sample (not shared with others)
    pub unique_gap_bases: usize,

    /// Gaps shared with at least one other sample
    pub shared_gap_bases: usize,

    /// Estimated core reduction if this sample is included (0.0 - 1.0)
    /// Higher = more damaging to core alignment
    pub estimated_core_reduction: f64,

    /// Number of unique gap regions
    pub unique_gap_regions: usize,

    /// Coordinates of unique gaps (sample-specific problems)
    pub unique_gap_coords: Vec<Region>,

    /// Average pairwise gap quality (from new pairwise analysis)
    pub avg_gap_quality: f64,

    /// Risk level
    pub risk_level: RiskLevel,

    /// Recommendation
    pub recommendation: String,
}

/// Risk level for including sample in SNP analysis
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum RiskLevel {
    /// Safe to include, minimal impact on core
    Low,
    /// Some unique gaps, review recommended
    Medium,
    /// Many unique gaps, likely to significantly reduce core
    High,
    /// Extreme unique gaps, will severely damage core alignment
    Critical,
}

impl std::fmt::Display for RiskLevel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RiskLevel::Low => write!(f, "LOW"),
            RiskLevel::Medium => write!(f, "MEDIUM"),
            RiskLevel::High => write!(f, "HIGH"),
            RiskLevel::Critical => write!(f, "CRITICAL"),
        }
    }
}

/// Complete impact analysis results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ImpactAnalysis {
    /// Per-sample impact predictions
    pub samples: Vec<SampleImpact>,

    /// Estimated core size if all samples included
    pub estimated_core_all: usize,

    /// Reference/genome length (baseline)
    pub reference_length: usize,

    /// Core reduction percentage with all samples
    pub total_core_reduction: f64,

    /// Samples recommended for exclusion
    pub exclude_recommendations: Vec<String>,

    /// Samples recommended for review
    pub review_recommendations: Vec<String>,

    /// Dataset gap quality score (from pairwise analysis)
    pub dataset_gap_quality: f64,
}

/// Analyze impact of each sample on core alignment
pub fn analyze_impact(
    sample_names: &[String],
    pairwise: &PairwiseResults,
    reference: &ReferenceResults,
    reference_length: usize,
) -> ImpactAnalysis {
    // Collect positions on reference not covered by each sample
    let mut sample_ref_gaps: HashMap<String, HashSet<usize>> = HashMap::new();

    for name in sample_names {
        sample_ref_gaps.insert(name.clone(), HashSet::new());
    }

    // Get gap positions from reference alignment
    for alignment in &reference.alignments {
        if let Some(gaps) = sample_ref_gaps.get_mut(&alignment.sample_name) {
            for region in &alignment.reference_uncovered_regions {
                for pos in region.start..region.end {
                    gaps.insert(pos);
                }
            }
        }
    }

    // Get pairwise gap quality for each sample
    let sample_gap_quality: HashMap<String, f64> = sample_names
        .iter()
        .map(|name| (name.clone(), pairwise.avg_quality(name)))
        .collect();

    // Calculate unique gaps per sample (in REFERENCE coordinates)
    let mut sample_impacts: Vec<SampleImpact> = Vec::new();

    for name in sample_names {
        let my_gaps = sample_ref_gaps.get(name).unwrap();

        // Find gaps unique to this sample (not in any other sample)
        let mut unique_positions: HashSet<usize> = my_gaps.clone();
        for (other_name, other_gaps) in &sample_ref_gaps {
            if other_name != name {
                unique_positions.retain(|pos| !other_gaps.contains(pos));
            }
        }

        let gap_quality = *sample_gap_quality.get(name).unwrap_or(&1.0);

        // Convert unique positions to regions
        let unique_gap_coords = positions_to_regions(&unique_positions);
        let unique_gap_bases = unique_positions.len();
        let shared_gap_bases = my_gaps.len() - unique_gap_bases;

        // Estimate core reduction
        let estimated_core_reduction = if reference_length > 0 {
            unique_gap_bases as f64 / reference_length as f64
        } else {
            0.0
        };

        // Determine risk level based on unique gap percentage and gap quality
        let risk_level = if estimated_core_reduction > 0.05 {
            RiskLevel::Critical  // >5% core reduction
        } else if estimated_core_reduction > 0.02 {
            RiskLevel::High      // 2-5% core reduction
        } else if estimated_core_reduction > 0.005 || gap_quality < 0.3 {
            RiskLevel::Medium    // 0.5-2% core reduction OR low gap quality
        } else {
            RiskLevel::Low       // <0.5% core reduction
        };

        let recommendation = match risk_level {
            RiskLevel::Critical => format!(
                "EXCLUDE: Would reduce core by {:.1}% ({} unique gap bases). Severely damages SNP resolution.",
                estimated_core_reduction * 100.0, unique_gap_bases
            ),
            RiskLevel::High => format!(
                "REVIEW: Would reduce core by {:.1}% ({} unique gap bases). Consider excluding.",
                estimated_core_reduction * 100.0, unique_gap_bases
            ),
            RiskLevel::Medium => format!(
                "CAUTION: Core reduction {:.2}% ({} unique gaps). Gap quality: {:.2}",
                estimated_core_reduction * 100.0, unique_gap_bases, gap_quality
            ),
            RiskLevel::Low => format!(
                "OK: Minimal impact ({:.3}%, {} unique gaps). Gap quality: {:.2}",
                estimated_core_reduction * 100.0, unique_gap_bases, gap_quality
            ),
        };

        sample_impacts.push(SampleImpact {
            sample: name.clone(),
            unique_gap_bases,
            shared_gap_bases,
            estimated_core_reduction,
            unique_gap_regions: unique_gap_coords.len(),
            unique_gap_coords,
            avg_gap_quality: gap_quality,
            risk_level,
            recommendation,
        });
    }

    // Calculate total core with all samples (union of all gaps)
    let all_gaps: HashSet<usize> = sample_ref_gaps.values()
        .flat_map(|g| g.iter().cloned())
        .collect();
    let estimated_core_all = reference_length.saturating_sub(all_gaps.len());
    let total_core_reduction = if reference_length > 0 {
        all_gaps.len() as f64 / reference_length as f64
    } else {
        0.0
    };

    // Recommendations
    let exclude_recommendations: Vec<String> = sample_impacts.iter()
        .filter(|s| s.risk_level == RiskLevel::Critical)
        .map(|s| s.sample.clone())
        .collect();

    let review_recommendations: Vec<String> = sample_impacts.iter()
        .filter(|s| s.risk_level == RiskLevel::High || s.risk_level == RiskLevel::Medium)
        .map(|s| s.sample.clone())
        .collect();

    ImpactAnalysis {
        samples: sample_impacts,
        estimated_core_all,
        reference_length,
        total_core_reduction,
        exclude_recommendations,
        review_recommendations,
        dataset_gap_quality: pairwise.dataset_quality(),
    }
}

/// Convert a set of positions to merged regions
fn positions_to_regions(positions: &HashSet<usize>) -> Vec<Region> {
    if positions.is_empty() {
        return Vec::new();
    }

    let mut sorted: Vec<usize> = positions.iter().cloned().collect();
    sorted.sort_unstable();

    let mut regions = Vec::new();
    let mut start = sorted[0];
    let mut end = sorted[0] + 1;

    for &pos in sorted.iter().skip(1) {
        if pos == end {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_positions_to_regions() {
        let mut positions = HashSet::new();
        positions.insert(1);
        positions.insert(2);
        positions.insert(3);
        positions.insert(10);
        positions.insert(11);

        let regions = positions_to_regions(&positions);

        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].start, 1);
        assert_eq!(regions[0].end, 4);
        assert_eq!(regions[1].start, 10);
        assert_eq!(regions[1].end, 12);
    }
}
