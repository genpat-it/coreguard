//! Sample problematicness scoring and recommendations

use crate::fasta::Sample;
use crate::metrics::{calculate_zscores, SampleMetrics};
use serde::{Deserialize, Serialize};

/// Flags indicating why a sample is problematic
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum ProblemFlag {
    /// Identity significantly below mean
    LowIdentity,
    /// Coverage significantly below expected
    LowCoverage,
    /// High fraction of unaligned sequence
    HighUnaligned,
    /// Assembly is highly fragmented
    HighFragmentation,
    /// Sample is an outlier based on identity z-score
    IdentityOutlier,
    /// Sample has unusually high unaligned regions
    UnalignedOutlier,
    /// Low identity to reference
    LowRefIdentity,
}

impl std::fmt::Display for ProblemFlag {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ProblemFlag::LowIdentity => write!(f, "low_identity"),
            ProblemFlag::LowCoverage => write!(f, "low_coverage"),
            ProblemFlag::HighUnaligned => write!(f, "high_unaligned"),
            ProblemFlag::HighFragmentation => write!(f, "high_fragmentation"),
            ProblemFlag::IdentityOutlier => write!(f, "identity_outlier"),
            ProblemFlag::UnalignedOutlier => write!(f, "unaligned_outlier"),
            ProblemFlag::LowRefIdentity => write!(f, "low_ref_identity"),
        }
    }
}

/// Recommendation for how to handle the sample
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum Recommendation {
    /// Sample looks fine
    Ok,
    /// Sample should be reviewed manually
    Review,
    /// Sample should likely be excluded
    Exclude,
}

impl std::fmt::Display for Recommendation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Recommendation::Ok => write!(f, "OK"),
            Recommendation::Review => write!(f, "REVIEW"),
            Recommendation::Exclude => write!(f, "EXCLUDE"),
        }
    }
}

/// A scored sample with flags and recommendations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScoredSample {
    /// Sample name
    pub name: String,
    /// Overall problematicness score (0.0 = good, 1.0 = very problematic)
    pub score: f64,
    /// Whether sample is flagged for attention
    pub flagged: bool,
    /// Specific issues found
    pub flags: Vec<ProblemFlag>,
    /// Recommendation
    pub recommendation: Recommendation,
    /// Underlying metrics
    pub metrics: SampleMetrics,
}

/// Score all samples and generate recommendations
pub fn score_samples(
    samples: &[Sample],
    metrics: &[SampleMetrics],
    min_identity: f64,
) -> Vec<ScoredSample> {
    // Calculate z-scores
    let mut metrics = metrics.to_vec();
    calculate_zscores(&mut metrics);

    // Calculate median N50 for fragmentation comparison
    let mut n50s: Vec<usize> = metrics.iter().map(|m| m.n50).collect();
    n50s.sort_unstable();
    let median_n50 = if n50s.is_empty() {
        0
    } else {
        n50s[n50s.len() / 2]
    };

    // Score each sample
    samples
        .iter()
        .zip(metrics.iter())
        .map(|(sample, m)| {
            score_sample(
                sample,
                m,
                min_identity,
                median_n50,
            )
        })
        .collect()
}

fn score_sample(
    _sample: &Sample,
    metrics: &SampleMetrics,
    min_identity: f64,
    median_n50: usize,
) -> ScoredSample {
    let mut flags = Vec::new();
    let mut score: f64 = 0.0;

    // === Reference alignment checks ===

    // Check reference identity threshold
    if metrics.ref_identity < min_identity {
        flags.push(ProblemFlag::LowIdentity);
        score += 0.3;
    }

    // Check reference coverage (below 90% is concerning)
    if metrics.ref_reference_coverage < 0.90 {
        flags.push(ProblemFlag::LowCoverage);
        score += 0.2;
    }

    // Check gap quality (low quality means gaps in different places)
    if metrics.avg_gap_quality < 0.3 {
        flags.push(ProblemFlag::HighUnaligned);
        score += 0.25;
    }

    // Check z-score outliers (more than 2 std devs from mean)
    if metrics.coverage_zscore < -2.0 {
        flags.push(ProblemFlag::IdentityOutlier);
        score += 0.2;
    }

    if metrics.gap_quality_zscore < -2.0 {
        flags.push(ProblemFlag::UnalignedOutlier);
        score += 0.2;
    }

    // Check assembly fragmentation (N50 less than half median)
    if median_n50 > 0 && metrics.n50 < median_n50 / 2 {
        flags.push(ProblemFlag::HighFragmentation);
        score += 0.15;
    }

    // === Large unique gap penalty ===
    // Samples with many unique gaps will severely reduce core genome
    let unique_gap_ratio = if metrics.total_length > 0 {
        metrics.total_unique_gaps as f64 / metrics.total_length as f64
    } else {
        0.0
    };

    if unique_gap_ratio > 0.02 {
        flags.push(ProblemFlag::LowRefIdentity);  // Reusing flag for "problematic"
        score += 0.3;
    } else if unique_gap_ratio > 0.005 {
        score += 0.1;
    }

    // Clamp score to [0, 1]
    score = score.min(1.0);

    // Determine recommendation
    let recommendation = if score >= 0.5 {
        Recommendation::Exclude
    } else if score >= 0.2 || !flags.is_empty() {
        Recommendation::Review
    } else {
        Recommendation::Ok
    };

    let flagged = !flags.is_empty();

    ScoredSample {
        name: metrics.name.clone(),
        score,
        flagged,
        flags,
        recommendation,
        metrics: metrics.clone(),
    }
}
