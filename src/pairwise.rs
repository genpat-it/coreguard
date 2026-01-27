//! Pairwise gap analysis based on reference alignment
//!
//! For each pair of samples (A, B), calculates:
//! - Gap union: positions that are gaps in A OR B (relative to reference)
//! - Gap intersection: positions that are gaps in A AND B
//! - Quality score: intersection / union (1.0 = identical gaps, 0.0 = no overlap)
//!
//! This approach uses the reference as coordinate system, making gaps comparable.

use crate::gaps::Region;
use crate::reference::ReferenceResults;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Result of pairwise gap analysis between two samples
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairwiseGapResult {
    /// Name of first sample
    pub sample_a: String,
    /// Name of second sample
    pub sample_b: String,
    /// Total bases in gap union (gaps in A OR B)
    pub gap_union_bases: usize,
    /// Total bases in gap intersection (gaps in A AND B)
    pub gap_intersection_bases: usize,
    /// Quality score: intersection / union (1.0 = identical gaps, 0.0 = no shared gaps)
    /// High score = gaps are in same places = less core reduction when combined
    pub quality_score: f64,
    /// Unique gaps in A (not in B)
    pub unique_gap_a_bases: usize,
    /// Unique gaps in B (not in A)
    pub unique_gap_b_bases: usize,
    /// Core genome size if only these two samples are included
    /// (reference_length - gap_union_bases)
    pub pairwise_core_size: usize,
}

/// Container for all pairwise gap results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairwiseResults {
    /// All pairwise comparisons
    pub results: Vec<PairwiseGapResult>,
    /// Quick lookup by sample pair
    pub by_pair: HashMap<(String, String), usize>,
    /// Reference length used for calculations
    pub reference_length: usize,
    /// Reference name
    pub reference_name: String,
}

impl PairwiseResults {
    /// Create empty results (for modes without reference)
    pub fn empty() -> Self {
        PairwiseResults {
            results: Vec::new(),
            by_pair: HashMap::new(),
            reference_length: 0,
            reference_name: String::new(),
        }
    }

    /// Check if results are empty
    pub fn is_empty(&self) -> bool {
        self.results.is_empty()
    }

    /// Get result for a specific pair
    pub fn get(&self, a: &str, b: &str) -> Option<&PairwiseGapResult> {
        self.by_pair
            .get(&(a.to_string(), b.to_string()))
            .or_else(|| self.by_pair.get(&(b.to_string(), a.to_string())))
            .map(|&idx| &self.results[idx])
    }

    /// Get all results involving a specific sample
    pub fn for_sample(&self, name: &str) -> Vec<&PairwiseGapResult> {
        self.results
            .iter()
            .filter(|r| r.sample_a == name || r.sample_b == name)
            .collect()
    }

    /// Average quality score for a sample across all its pairs
    pub fn avg_quality(&self, name: &str) -> f64 {
        let results = self.for_sample(name);
        if results.is_empty() {
            return 1.0; // No pairs = perfect quality
        }
        results.iter().map(|r| r.quality_score).sum::<f64>() / results.len() as f64
    }

    /// Total unique gap bases this sample introduces across all pairs
    pub fn total_unique_gaps(&self, name: &str) -> usize {
        self.results
            .iter()
            .filter_map(|r| {
                if r.sample_a == name {
                    Some(r.unique_gap_a_bases)
                } else if r.sample_b == name {
                    Some(r.unique_gap_b_bases)
                } else {
                    None
                }
            })
            .sum()
    }

    /// Average pairwise core size for a sample
    pub fn avg_pairwise_core(&self, name: &str) -> usize {
        let results = self.for_sample(name);
        if results.is_empty() {
            return self.reference_length;
        }
        let sum: usize = results.iter().map(|r| r.pairwise_core_size).sum();
        sum / results.len()
    }

    /// Get the worst pair (lowest quality score)
    pub fn worst_pair(&self) -> Option<&PairwiseGapResult> {
        self.results.iter().min_by(|a, b| {
            a.quality_score.partial_cmp(&b.quality_score).unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Get the best pair (highest quality score)
    pub fn best_pair(&self) -> Option<&PairwiseGapResult> {
        self.results.iter().max_by(|a, b| {
            a.quality_score.partial_cmp(&b.quality_score).unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Dataset-level quality score (average of all pairwise scores)
    pub fn dataset_quality(&self) -> f64 {
        if self.results.is_empty() {
            return 1.0;
        }
        self.results.iter().map(|r| r.quality_score).sum::<f64>() / self.results.len() as f64
    }
}

/// Analyze pairwise gap relationships based on reference alignment results
///
/// For each pair of samples:
/// 1. Get their gap regions vs reference
/// 2. Calculate union and intersection
/// 3. Compute quality score
pub fn analyze_pairwise_gaps(reference_results: &ReferenceResults) -> PairwiseResults {
    let alignments = &reference_results.alignments;
    let ref_len = reference_results.reference_length;

    if alignments.len() < 2 {
        return PairwiseResults {
            results: Vec::new(),
            by_pair: HashMap::new(),
            reference_length: ref_len,
            reference_name: reference_results.reference_name.clone(),
        };
    }

    log::info!(
        "Computing pairwise gap analysis for {} samples ({} pairs)...",
        alignments.len(),
        alignments.len() * (alignments.len() - 1) / 2
    );

    // Pre-compute gap bitmaps for each sample (for fast set operations)
    let gap_bitmaps: HashMap<String, Vec<bool>> = alignments
        .iter()
        .map(|a| {
            let bitmap = regions_to_bitmap(&a.reference_uncovered_regions, ref_len);
            (a.sample_name.clone(), bitmap)
        })
        .collect();

    let mut results = Vec::new();

    // Compute all pairs
    for i in 0..alignments.len() {
        for j in (i + 1)..alignments.len() {
            let sample_a = &alignments[i].sample_name;
            let sample_b = &alignments[j].sample_name;

            let bitmap_a = &gap_bitmaps[sample_a];
            let bitmap_b = &gap_bitmaps[sample_b];

            // Calculate union and intersection
            let (union_bases, intersection_bases) = compute_set_operations(bitmap_a, bitmap_b);

            // Quality score: intersection / union
            let quality_score = if union_bases > 0 {
                intersection_bases as f64 / union_bases as f64
            } else {
                1.0 // No gaps in either = perfect quality
            };

            // Unique gaps
            let unique_a = bitmap_a.iter().zip(bitmap_b.iter())
                .filter(|(&a, &b)| a && !b)
                .count();
            let unique_b = bitmap_a.iter().zip(bitmap_b.iter())
                .filter(|(&a, &b)| !a && b)
                .count();

            results.push(PairwiseGapResult {
                sample_a: sample_a.clone(),
                sample_b: sample_b.clone(),
                gap_union_bases: union_bases,
                gap_intersection_bases: intersection_bases,
                quality_score,
                unique_gap_a_bases: unique_a,
                unique_gap_b_bases: unique_b,
                pairwise_core_size: ref_len.saturating_sub(union_bases),
            });
        }
    }

    // Build lookup index
    let by_pair: HashMap<(String, String), usize> = results
        .iter()
        .enumerate()
        .map(|(idx, r)| ((r.sample_a.clone(), r.sample_b.clone()), idx))
        .collect();

    log::info!(
        "Pairwise analysis complete: {} pairs, avg quality: {:.3}",
        results.len(),
        if results.is_empty() { 1.0 } else {
            results.iter().map(|r| r.quality_score).sum::<f64>() / results.len() as f64
        }
    );

    PairwiseResults {
        results,
        by_pair,
        reference_length: ref_len,
        reference_name: reference_results.reference_name.clone(),
    }
}

/// Convert gap regions to a boolean bitmap for fast set operations
fn regions_to_bitmap(regions: &[Region], length: usize) -> Vec<bool> {
    let mut bitmap = vec![false; length];
    for region in regions {
        let start = region.start.min(length);
        let end = region.end.min(length);
        for pos in start..end {
            bitmap[pos] = true;
        }
    }
    bitmap
}

/// Compute union and intersection sizes from two bitmaps
fn compute_set_operations(a: &[bool], b: &[bool]) -> (usize, usize) {
    let mut union = 0;
    let mut intersection = 0;

    for (&val_a, &val_b) in a.iter().zip(b.iter()) {
        if val_a || val_b {
            union += 1;
        }
        if val_a && val_b {
            intersection += 1;
        }
    }

    (union, intersection)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_regions_to_bitmap() {
        let regions = vec![
            Region::new(10, 20),
            Region::new(30, 40),
        ];
        let bitmap = regions_to_bitmap(&regions, 50);

        assert!(!bitmap[5]);
        assert!(bitmap[10]);
        assert!(bitmap[15]);
        assert!(!bitmap[20]);
        assert!(!bitmap[25]);
        assert!(bitmap[30]);
        assert!(bitmap[35]);
        assert!(!bitmap[40]);
    }

    #[test]
    fn test_set_operations() {
        // A: gaps at 10-20
        // B: gaps at 15-25
        // Union: 10-25 = 15 positions
        // Intersection: 15-20 = 5 positions
        let mut a = vec![false; 50];
        let mut b = vec![false; 50];

        for i in 10..20 { a[i] = true; }
        for i in 15..25 { b[i] = true; }

        let (union, intersection) = compute_set_operations(&a, &b);
        assert_eq!(union, 15);
        assert_eq!(intersection, 5);
    }

    #[test]
    fn test_quality_score() {
        // Identical gaps = quality 1.0
        let a = vec![false, false, true, true, false];
        let b = vec![false, false, true, true, false];
        let (union, intersection) = compute_set_operations(&a, &b);
        let quality = intersection as f64 / union as f64;
        assert!((quality - 1.0).abs() < 0.001);

        // No overlap = quality 0.0
        let a = vec![true, true, false, false, false];
        let b = vec![false, false, false, true, true];
        let (union, intersection) = compute_set_operations(&a, &b);
        let quality = if union > 0 { intersection as f64 / union as f64 } else { 1.0 };
        assert!((quality - 0.0).abs() < 0.001);
    }
}
