//! FASTA file parsing and sample representation

use anyhow::{Context, Result};
use md5;
use needletail::parse_fastx_file;
use std::path::Path;

/// A single contig from an assembly
#[derive(Debug, Clone)]
pub struct Contig {
    /// Contig name/header
    pub name: String,
    /// Contig sequence
    pub sequence: Vec<u8>,
}

/// Represents a loaded sample with its sequences
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct Sample {
    /// Sample name (derived from filename)
    pub name: String,
    /// Individual contigs (for proper multi-contig handling)
    pub contigs: Vec<Contig>,
    /// Concatenated sequence (all contigs joined) - for legacy compatibility
    pub sequence: Vec<u8>,
    /// Total sequence length
    pub total_length: usize,
    /// Number of contigs
    pub num_contigs: usize,
    /// N50 statistic
    pub n50: usize,
    /// MD5 hash of sequence (for deduplication)
    pub sequence_hash: String,
}

impl Sample {
    /// Calculate N50 from contig lengths
    fn calculate_n50(lengths: &[usize]) -> usize {
        if lengths.is_empty() {
            return 0;
        }

        let total: usize = lengths.iter().sum();
        let half = total / 2;

        let mut sorted = lengths.to_vec();
        sorted.sort_unstable_by(|a, b| b.cmp(a)); // Sort descending

        let mut cumsum = 0;
        for len in sorted {
            cumsum += len;
            if cumsum >= half {
                return len;
            }
        }

        lengths[0]
    }
}

/// Load a sample from a FASTA file (supports gzip compression)
pub fn load_sample(path: &Path) -> Result<Sample> {
    // Handle gzipped files: sample.fa.gz -> sample, sample.fasta.gz -> sample
    let name = {
        let stem = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");

        // If the stem ends with .fa, .fasta, or .fna, strip that too
        // (handles double extensions like .fa.gz)
        let name = stem
            .strip_suffix(".fa")
            .or_else(|| stem.strip_suffix(".fasta"))
            .or_else(|| stem.strip_suffix(".fna"))
            .unwrap_or(stem);

        name.to_string()
    };

    let mut contigs = Vec::new();
    let mut sequence = Vec::new();

    let mut reader =
        parse_fastx_file(path).with_context(|| format!("Failed to open FASTA: {}", path.display()))?;

    while let Some(record) = reader.next() {
        let record = record.with_context(|| format!("Failed to parse record in {}", path.display()))?;
        let seq = record.seq();
        let contig_name = String::from_utf8_lossy(record.id()).to_string();

        // Store contig separately
        contigs.push(Contig {
            name: contig_name,
            sequence: seq.to_vec(),
        });

        // Also build concatenated sequence for legacy compatibility
        sequence.extend_from_slice(&seq);
    }

    let contig_lengths: Vec<usize> = contigs.iter().map(|c| c.sequence.len()).collect();
    let total_length = sequence.len();
    let num_contigs = contigs.len();
    let n50 = Sample::calculate_n50(&contig_lengths);

    // Compute MD5 hash for deduplication
    // Normalize: uppercase and sort contigs by sequence to handle different ordering
    let mut sorted_seqs: Vec<Vec<u8>> = contigs
        .iter()
        .map(|c| c.sequence.to_ascii_uppercase())
        .collect();
    sorted_seqs.sort();
    let combined: Vec<u8> = sorted_seqs.into_iter().flatten().collect();
    let sequence_hash = format!("{:x}", md5::compute(&combined));

    Ok(Sample {
        name,
        contigs,
        sequence,
        total_length,
        num_contigs,
        n50,
        sequence_hash,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_n50_calculation() {
        // Example: contigs of length 10, 20, 30, 40
        // Total = 100, half = 50
        // Sorted desc: 40, 30, 20, 10
        // Cumsum: 40, 70, 90, 100
        // First to exceed 50 is 30 (cumsum 70)
        let lengths = vec![10, 20, 30, 40];
        assert_eq!(Sample::calculate_n50(&lengths), 30);
    }

    #[test]
    fn test_n50_single_contig() {
        let lengths = vec![1000];
        assert_eq!(Sample::calculate_n50(&lengths), 1000);
    }

    #[test]
    fn test_n50_empty() {
        let lengths: Vec<usize> = vec![];
        assert_eq!(Sample::calculate_n50(&lengths), 0);
    }
}
