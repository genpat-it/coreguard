//! Gap region tracking and analysis
//!
//! Tracks genomic regions with gaps/unaligned sequences at each analysis level.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A genomic region (0-based, half-open interval)
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Region {
    /// Start position (0-based, inclusive)
    pub start: usize,
    /// End position (0-based, exclusive)
    pub end: usize,
}

#[allow(dead_code)]
impl Region {
    pub fn new(start: usize, end: usize) -> Self {
        Region { start, end }
    }

    pub fn len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }

    pub fn is_empty(&self) -> bool {
        self.end <= self.start
    }
}

/// Gap information for a single sample at one analysis level
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleGaps {
    /// Sample name
    pub sample: String,
    /// Contig/chromosome name (for multi-contig assemblies)
    pub contig: Option<String>,
    /// Unaligned/gap regions
    pub regions: Vec<Region>,
    /// Total gap bases
    pub total_bases: usize,
}

#[allow(dead_code)]
impl SampleGaps {
    pub fn new(sample: &str) -> Self {
        SampleGaps {
            sample: sample.to_string(),
            contig: None,
            regions: Vec::new(),
            total_bases: 0,
        }
    }

    pub fn add_region(&mut self, start: usize, end: usize) {
        if end > start {
            self.total_bases += end - start;
            self.regions.push(Region::new(start, end));
        }
    }

    /// Merge overlapping/adjacent regions and sort
    pub fn merge_regions(&mut self) {
        if self.regions.is_empty() {
            return;
        }

        self.regions.sort_by_key(|r| r.start);

        let mut merged = Vec::new();
        let mut current = self.regions[0].clone();

        for region in self.regions.iter().skip(1) {
            if region.start <= current.end {
                // Overlapping or adjacent - extend
                current.end = current.end.max(region.end);
            } else {
                // Gap - save current and start new
                merged.push(current);
                current = region.clone();
            }
        }
        merged.push(current);

        self.regions = merged;
        self.total_bases = self.regions.iter().map(|r| r.len()).sum();
    }
}

/// Gaps from pairwise alignments
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairwiseGaps {
    /// Sample A name
    pub sample_a: String,
    /// Sample B name
    pub sample_b: String,
    /// Regions in sample A not aligned to B
    pub unaligned_in_a: Vec<Region>,
    /// Regions in sample B not aligned to A
    pub unaligned_in_b: Vec<Region>,
}

/// Gaps from reference alignment
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceGaps {
    /// Sample name
    pub sample: String,
    /// Regions in sample not aligned to reference
    pub sample_unaligned: Vec<Region>,
    /// Regions in reference not covered by this sample
    pub reference_uncovered: Vec<Region>,
}

/// Complete gap report across all analysis levels
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GapReport {
    /// Pairwise alignment gaps (per sample pair)
    pub pairwise: Vec<PairwiseGaps>,
    /// Reference alignment gaps (per sample)
    pub reference: Vec<ReferenceGaps>,
    /// Summary statistics
    pub summary: GapSummary,
}

/// Summary of gaps across all levels
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GapSummary {
    /// Per-sample total unaligned bases (pairwise)
    pub pairwise_unaligned: HashMap<String, usize>,
    /// Per-sample total unaligned bases (vs reference)
    pub reference_unaligned: HashMap<String, usize>,
    /// Samples with unusual gap patterns
    pub outliers: Vec<String>,
}

/// Calculate unaligned regions from aligned intervals
pub fn calculate_unaligned_regions(
    total_length: usize,
    aligned_intervals: &[(usize, usize)],
) -> Vec<Region> {
    if aligned_intervals.is_empty() {
        return vec![Region::new(0, total_length)];
    }

    // Sort and merge aligned intervals
    let mut intervals: Vec<(usize, usize)> = aligned_intervals.to_vec();
    intervals.sort_by_key(|&(start, _)| start);

    let mut merged = Vec::new();
    let mut current = intervals[0];

    for &(start, end) in intervals.iter().skip(1) {
        if start <= current.1 {
            current.1 = current.1.max(end);
        } else {
            merged.push(current);
            current = (start, end);
        }
    }
    merged.push(current);

    // Find gaps between aligned regions
    let mut unaligned = Vec::new();

    // Gap before first alignment
    if merged[0].0 > 0 {
        unaligned.push(Region::new(0, merged[0].0));
    }

    // Gaps between alignments
    for i in 0..merged.len() - 1 {
        let gap_start = merged[i].1;
        let gap_end = merged[i + 1].0;
        if gap_end > gap_start {
            unaligned.push(Region::new(gap_start, gap_end));
        }
    }

    // Gap after last alignment
    let last_end = merged.last().unwrap().1;
    if last_end < total_length {
        unaligned.push(Region::new(last_end, total_length));
    }

    unaligned
}

/// Load masked/gap regions from Snippy aligned FASTA files (N positions)
/// Returns a HashMap of sample name -> Vec<Region>
pub fn load_snippy_gaps(snippy_dir: &std::path::Path) -> std::collections::HashMap<String, Vec<Region>> {
    use std::fs;
    use std::io::{BufRead, BufReader};

    let mut result = std::collections::HashMap::new();

    // Find all *_snippy directories
    let entries = match fs::read_dir(snippy_dir) {
        Ok(e) => e,
        Err(_) => return result,
    };

    for entry in entries.flatten() {
        let path = entry.path();
        if !path.is_dir() {
            continue;
        }

        let dir_name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if !dir_name.ends_with("_snippy") {
            continue;
        }

        // Extract sample name (remove _snippy suffix)
        let sample_name = dir_name.trim_end_matches("_snippy");

        // Look for aligned.fa file
        let aligned_fa = path.join(format!("{}.aligned.fa", dir_name));
        if !aligned_fa.exists() {
            continue;
        }

        // Parse FASTA and extract N regions
        let file = match fs::File::open(&aligned_fa) {
            Ok(f) => f,
            Err(_) => continue,
        };

        let reader = BufReader::new(file);
        let mut seq = String::new();

        for line in reader.lines().flatten() {
            if line.starts_with('>') {
                continue;
            }
            seq.push_str(line.trim());
        }

        // Find N regions (1-based positions for consistency)
        let mut regions = Vec::new();
        let mut in_n_region = false;
        let mut region_start = 0;

        for (i, c) in seq.chars().enumerate() {
            let is_n = c == 'N' || c == 'n';
            if is_n && !in_n_region {
                // Start of N region
                region_start = i + 1; // 1-based
                in_n_region = true;
            } else if !is_n && in_n_region {
                // End of N region
                regions.push(Region::new(region_start, i + 1)); // end is exclusive
                in_n_region = false;
            }
        }

        // Handle N region at end of sequence
        if in_n_region {
            regions.push(Region::new(region_start, seq.len() + 1));
        }

        if !regions.is_empty() {
            log::debug!("Loaded {} Snippy gap regions for {}", regions.len(), sample_name);
            result.insert(sample_name.to_string(), regions);
        }
    }

    result
}

/// Load zero-coverage gap regions from CFSAN pileup files
/// Returns a HashMap of sample name -> Vec<Region>
pub fn load_cfsan_gaps(cfsan_dir: &std::path::Path) -> std::collections::HashMap<String, Vec<Region>> {
    use std::fs;
    use std::io::{BufRead, BufReader};

    let mut result = std::collections::HashMap::new();

    // CFSAN structure: cfsan_dir/samples/SAMPLE_NAME/reads.all.pileup
    let samples_dir = cfsan_dir.join("samples");
    if !samples_dir.exists() {
        return result;
    }

    let entries = match fs::read_dir(&samples_dir) {
        Ok(e) => e,
        Err(_) => return result,
    };

    for entry in entries.flatten() {
        let path = entry.path();
        if !path.is_dir() {
            continue;
        }

        let sample_name = match path.file_name().and_then(|n| n.to_str()) {
            Some(n) => n.to_string(),
            None => continue,
        };

        // Look for pileup file
        let pileup_file = path.join("reads.all.pileup");
        if !pileup_file.exists() {
            continue;
        }

        // Parse pileup and extract zero coverage positions
        let file = match fs::File::open(&pileup_file) {
            Ok(f) => f,
            Err(_) => continue,
        };

        let reader = BufReader::new(file);
        let mut zero_positions: Vec<usize> = Vec::new();

        for line in reader.lines().flatten() {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 4 {
                // Column 4 is depth
                if let (Ok(pos), Ok(depth)) = (fields[1].parse::<usize>(), fields[3].parse::<u32>()) {
                    if depth == 0 {
                        zero_positions.push(pos);
                    }
                }
            }
        }

        // Convert positions to regions (merge consecutive)
        if zero_positions.is_empty() {
            continue;
        }

        zero_positions.sort();
        let mut regions = Vec::new();
        let mut region_start = zero_positions[0];
        let mut region_end = zero_positions[0];

        for &pos in zero_positions.iter().skip(1) {
            if pos == region_end + 1 {
                // Consecutive - extend region
                region_end = pos;
            } else {
                // Gap - save current region and start new one
                regions.push(Region::new(region_start, region_end + 1));
                region_start = pos;
                region_end = pos;
            }
        }
        // Save last region
        regions.push(Region::new(region_start, region_end + 1));

        if !regions.is_empty() {
            log::debug!("Loaded {} CFSAN gap regions for {} ({} bp)",
                regions.len(), sample_name,
                regions.iter().map(|r| r.len()).sum::<usize>());
            result.insert(sample_name, regions);
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_unaligned() {
        // Total length 100, aligned 10-30 and 50-80
        let aligned = vec![(10, 30), (50, 80)];
        let unaligned = calculate_unaligned_regions(100, &aligned);

        assert_eq!(unaligned.len(), 3);
        assert_eq!(unaligned[0], Region::new(0, 10));   // Before first
        assert_eq!(unaligned[1], Region::new(30, 50));  // Between
        assert_eq!(unaligned[2], Region::new(80, 100)); // After last
    }

    #[test]
    fn test_merge_regions() {
        let mut gaps = SampleGaps::new("test");
        gaps.add_region(10, 20);
        gaps.add_region(15, 25);  // Overlaps
        gaps.add_region(30, 40);  // Separate

        gaps.merge_regions();

        assert_eq!(gaps.regions.len(), 2);
        assert_eq!(gaps.regions[0], Region::new(10, 25));
        assert_eq!(gaps.regions[1], Region::new(30, 40));
    }
}
