//! BAM validation for detecting Snippy artifacts
//!
//! Validates VCF-reported polymorphisms against raw BAM pileup data.
//!
//! The Snippy bug: When Snippy reports a position as polymorphic (samples differ),
//! but examining ALL BAM files shows that all samples have the SAME consensus base,
//! it's an artifact caused by complex variant decomposition.

use anyhow::Result;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::process::Command;

/// Pairwise distance between two samples
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairwiseDistance {
    pub sample_a: String,
    pub sample_b: String,
    pub distance: usize,
}

/// Distance matrix with sample order
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DistanceMatrix {
    pub samples: Vec<String>,
    /// Matrix[i][j] = distance between samples[i] and samples[j]
    pub matrix: Vec<Vec<usize>>,
    /// List of pairwise distances for easier access
    pub pairwise: Vec<PairwiseDistance>,
}

/// Summary of BAM validation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BamValidationSummary {
    /// Total positions checked
    pub total_checked: usize,
    /// Positions confirmed as real variants (samples truly differ)
    pub real_variants: usize,
    /// Positions detected as artifacts (all samples same in BAM)
    pub artifacts: usize,
    /// Positions that couldn't be validated (insufficient BAM data)
    pub no_data: usize,
    /// Artifact rate as percentage
    pub artifact_rate: f64,
    /// Corrected Snippy-only count
    pub corrected_snippy_only: usize,
    /// Original Snippy-only count
    pub original_snippy_only: usize,
    /// Corrected distance matrix (based on BAM consensus at real variant positions)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub corrected_distance_matrix: Option<DistanceMatrix>,
}

/// Get BAM consensus for a single sample at multiple positions
fn get_sample_consensus_batch(
    bam_path: &Path,
    positions: &[(String, usize)], // (chrom, pos)
) -> HashMap<String, char> {
    let mut results = HashMap::new();

    if !bam_path.exists() || positions.is_empty() {
        return results;
    }

    // Create a temporary BED file with all positions
    let bed_content: String = positions
        .iter()
        .map(|(chrom, pos)| format!("{}\t{}\t{}", chrom, pos - 1, pos)) // BED is 0-based
        .collect::<Vec<_>>()
        .join("\n");

    // Write BED to temp file (use thread ID to avoid conflicts in parallel execution)
    let thread_id = std::thread::current().id();
    let temp_bed = std::env::temp_dir().join(format!(
        "coreguard_positions_{:?}_{}.bed",
        thread_id,
        std::process::id()
    ));

    if std::fs::write(&temp_bed, &bed_content).is_err() {
        return results;
    }

    // Run samtools mpileup with positions file
    let output = Command::new("samtools")
        .args([
            "mpileup",
            "-l", temp_bed.to_str().unwrap_or(""),
            "-Q", "20", // Min base quality
            bam_path.to_str().unwrap_or(""),
        ])
        .output();

    // Clean up temp file
    let _ = std::fs::remove_file(&temp_bed);

    let output = match output {
        Ok(o) if o.status.success() => o,
        _ => return results,
    };

    let stdout = String::from_utf8_lossy(&output.stdout);

    // Parse mpileup output
    for line in stdout.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }

        let chrom = fields[0];
        let pos: usize = match fields[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        let ref_base = fields[2].chars().next().unwrap_or('N');
        let bases = fields[4];

        // Count bases (including reference matches)
        let mut counts: HashMap<char, usize> = HashMap::new();
        let mut skip_next = 0;
        let mut in_indel_len = false;
        let mut indel_len_str = String::new();

        for c in bases.chars() {
            if skip_next > 0 {
                skip_next -= 1;
                continue;
            }

            // Handle indel length parsing
            if in_indel_len {
                if c.is_ascii_digit() {
                    indel_len_str.push(c);
                    continue;
                } else {
                    // End of length, skip the indel bases
                    if let Ok(len) = indel_len_str.parse::<usize>() {
                        skip_next = len;
                    }
                    in_indel_len = false;
                    indel_len_str.clear();
                    continue;
                }
            }

            match c {
                '+' | '-' => {
                    in_indel_len = true;
                    indel_len_str.clear();
                    continue;
                }
                '^' => {
                    skip_next = 1; // Skip the next char (mapping quality)
                    continue;
                }
                '$' | '*' | '<' | '>' => continue,
                '.' | ',' => {
                    // Reference match
                    let base = ref_base.to_ascii_uppercase();
                    if "ATGC".contains(base) {
                        *counts.entry(base).or_insert(0) += 1;
                    }
                }
                _ => {
                    let base = c.to_ascii_uppercase();
                    if "ATGC".contains(base) {
                        *counts.entry(base).or_insert(0) += 1;
                    }
                }
            }
        }

        // Get most common base (require at least 3 reads)
        if let Some((base, _)) = counts
            .into_iter()
            .filter(|(_, count)| *count >= 3)
            .max_by_key(|(_, count)| *count)
        {
            let key = format!("{}:{}", chrom, pos);
            results.insert(key, base);
        }
    }

    results
}

/// Find BAM file for a sample
fn find_bam_file(sample_name: &str, bam_dir: &Path) -> Option<std::path::PathBuf> {
    let bam_patterns = [
        format!("{}_snippy/{}_snippy.bam", sample_name, sample_name),
        format!("{}/{}.bam", sample_name, sample_name),
        format!("{}.bam", sample_name),
    ];

    for pattern in &bam_patterns {
        let path = bam_dir.join(pattern);
        if path.exists() {
            return Some(path);
        }
    }
    None
}

/// Calculate pairwise distance matrix from BAM consensus data
fn calculate_distance_matrix(
    samples: &[String],
    position_sample_consensus: &HashMap<String, HashMap<String, char>>,
    real_variant_positions: &HashSet<String>,
) -> DistanceMatrix {
    let n = samples.len();
    let mut matrix = vec![vec![0usize; n]; n];
    let mut pairwise = Vec::new();

    // For each pair of samples
    for i in 0..n {
        for j in (i + 1)..n {
            let sample_a = &samples[i];
            let sample_b = &samples[j];
            let mut distance = 0usize;

            // Count differences at real variant positions
            for pos_key in real_variant_positions {
                if let Some(sample_bases) = position_sample_consensus.get(pos_key) {
                    let base_a = sample_bases.get(sample_a);
                    let base_b = sample_bases.get(sample_b);

                    if let (Some(&a), Some(&b)) = (base_a, base_b) {
                        if a != b {
                            distance += 1;
                        }
                    }
                }
            }

            matrix[i][j] = distance;
            matrix[j][i] = distance;

            pairwise.push(PairwiseDistance {
                sample_a: sample_a.clone(),
                sample_b: sample_b.clone(),
                distance,
            });
        }
    }

    DistanceMatrix {
        samples: samples.to_vec(),
        matrix,
        pairwise,
    }
}

/// Validate Snippy-only positions by checking if ALL samples have the same BAM consensus
///
/// The Snippy bug: position reported as polymorphic but all BAMs show same base
pub fn validate_snippy_only_positions(
    snippy_only_positions: &[(String, usize, String, String)], // (chrom, pos, sample, alt)
    bam_dir: &Path,
) -> Result<BamValidationSummary> {
    log::info!(
        "Validating {} Snippy-only positions for artifact detection...",
        snippy_only_positions.len()
    );

    // Get unique positions (a position might be called in multiple samples)
    let mut unique_positions: HashSet<(String, usize)> = HashSet::new();
    let mut position_to_samples: HashMap<String, HashSet<String>> = HashMap::new();
    let mut all_samples: HashSet<String> = HashSet::new();

    for (chrom, pos, sample, _alt) in snippy_only_positions {
        unique_positions.insert((chrom.clone(), *pos));
        let key = format!("{}:{}", chrom, pos);
        position_to_samples
            .entry(key)
            .or_default()
            .insert(sample.clone());
        all_samples.insert(sample.clone());
    }

    let positions_vec: Vec<(String, usize)> = unique_positions.into_iter().collect();
    let mut samples_vec: Vec<String> = all_samples.into_iter().collect();
    samples_vec.sort(); // Sort for consistent ordering

    log::info!(
        "Checking {} unique positions across {} samples...",
        positions_vec.len(),
        samples_vec.len()
    );

    // For each sample, get BAM consensus at all positions (in parallel)
    let sample_consensus: Vec<(String, HashMap<String, char>)> = samples_vec
        .par_iter()
        .filter_map(|sample| {
            let bam_path = find_bam_file(sample, bam_dir)?;
            let consensus = get_sample_consensus_batch(&bam_path, &positions_vec);
            Some((sample.clone(), consensus))
        })
        .collect();

    // Build: position -> sample -> consensus_base
    let mut position_sample_consensus: HashMap<String, HashMap<String, char>> = HashMap::new();
    for (sample, consensus_map) in &sample_consensus {
        for (pos_key, &base) in consensus_map {
            position_sample_consensus
                .entry(pos_key.clone())
                .or_default()
                .insert(sample.clone(), base);
        }
    }

    log::info!("Analyzing consensus across samples for artifact detection...");

    // Now check each position: if all samples with BAM data have the SAME consensus, it's an artifact
    let mut artifacts = 0usize;
    let mut real_variants = 0usize;
    let mut no_data = 0usize;
    let mut real_variant_positions: HashSet<String> = HashSet::new();

    for (chrom, pos) in &positions_vec {
        let pos_key = format!("{}:{}", chrom, pos);

        match position_sample_consensus.get(&pos_key) {
            Some(sample_bases) if sample_bases.len() >= 2 => {
                // We have consensus data from at least 2 samples
                let unique_bases: HashSet<char> = sample_bases.values().cloned().collect();

                if unique_bases.len() == 1 {
                    // ALL samples have the SAME base in BAM → ARTIFACT!
                    // Snippy reported this as polymorphic but BAMs show no difference
                    artifacts += 1;
                } else {
                    // Samples have DIFFERENT bases in BAM → real polymorphism
                    real_variants += 1;
                    real_variant_positions.insert(pos_key);
                }
            }
            Some(sample_bases) if sample_bases.len() == 1 => {
                // Only one sample has BAM data - can't determine if artifact
                // Count as real for now (conservative)
                real_variants += 1;
            }
            _ => {
                // No BAM data for this position
                no_data += 1;
            }
        }
    }

    let total_checked = positions_vec.len();
    let artifact_rate = if total_checked > no_data {
        100.0 * artifacts as f64 / (total_checked - no_data) as f64
    } else {
        0.0
    };

    // Original count is the number of sample-position pairs (not unique positions)
    let original_snippy_only = snippy_only_positions.len();

    // Corrected count: remove artifacts (each artifact position might affect multiple samples)
    // Count how many sample-position pairs are at artifact positions
    let artifact_positions: HashSet<String> = positions_vec
        .iter()
        .filter(|(chrom, pos)| {
            let pos_key = format!("{}:{}", chrom, pos);
            if let Some(sample_bases) = position_sample_consensus.get(&pos_key) {
                if sample_bases.len() >= 2 {
                    let unique_bases: HashSet<char> = sample_bases.values().cloned().collect();
                    return unique_bases.len() == 1; // artifact
                }
            }
            false
        })
        .map(|(chrom, pos)| format!("{}:{}", chrom, pos))
        .collect();

    let artifact_sample_pairs: usize = snippy_only_positions
        .iter()
        .filter(|(chrom, pos, _, _)| artifact_positions.contains(&format!("{}:{}", chrom, pos)))
        .count();

    let corrected_snippy_only = original_snippy_only - artifact_sample_pairs;

    // Calculate corrected distance matrix based on real variant positions only
    log::info!(
        "Calculating corrected distance matrix from {} real variant positions...",
        real_variant_positions.len()
    );

    let corrected_distance_matrix = if samples_vec.len() >= 2 && !real_variant_positions.is_empty() {
        Some(calculate_distance_matrix(
            &samples_vec,
            &position_sample_consensus,
            &real_variant_positions,
        ))
    } else {
        None
    };

    if let Some(ref dm) = corrected_distance_matrix {
        log::info!("Corrected distance matrix ({} samples):", dm.samples.len());
        for pw in &dm.pairwise {
            log::info!("  {} - {}: {} SNPs", pw.sample_a, pw.sample_b, pw.distance);
        }
    }

    log::info!(
        "BAM validation complete: {} artifacts ({:.1}%), {} real, {} no data",
        artifacts,
        artifact_rate,
        real_variants,
        no_data
    );
    log::info!(
        "Corrected Snippy-only: {} -> {} ({} artifact pairs removed)",
        original_snippy_only,
        corrected_snippy_only,
        artifact_sample_pairs
    );

    Ok(BamValidationSummary {
        total_checked,
        real_variants,
        artifacts,
        no_data,
        artifact_rate,
        corrected_snippy_only,
        original_snippy_only,
        corrected_distance_matrix,
    })
}
