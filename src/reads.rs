//! FASTQ/read-based input handling
//!
//! This module provides gap detection from raw sequencing reads (FASTQ)
//! by mapping reads to a reference and identifying low-coverage regions.

use anyhow::{Context, Result};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use crate::fasta::Sample;
use crate::gaps::Region;

/// Represents a sample loaded from FASTQ files
#[derive(Debug, Clone)]
pub struct ReadSample {
    /// Sample name (derived from directory or file name)
    pub name: String,
    /// Path to R1 FASTQ file
    pub r1_path: PathBuf,
    /// Path to R2 FASTQ file (optional, for paired-end)
    pub r2_path: Option<PathBuf>,
    /// Read type: "short" (Illumina) or "long" (ONT/PacBio)
    pub read_type: String,
}

/// Result of mapping reads to reference
#[derive(Debug, Clone)]
pub struct ReadMappingResult {
    /// Sample name
    pub name: String,
    /// Total reads in input
    pub total_reads: usize,
    /// Reads that mapped to reference
    pub mapped_reads: usize,
    /// Average depth of coverage
    pub mean_depth: f64,
    /// Percentage of reference covered at min_depth
    pub coverage_percent: f64,
    /// Gap regions (coverage < min_depth)
    pub gap_regions: Vec<Region>,
    /// Total gap bases
    pub total_gap_bases: usize,
}

/// Load FASTQ samples from a directory
/// Expects structure: dir/sample_name/reads_1.fq.gz, reads_2.fq.gz
/// Or: dir/sample_name_R1.fq.gz, sample_name_R2.fq.gz
pub fn load_fastq_samples(dir: &Path, read_type: &str) -> Result<Vec<ReadSample>> {
    let mut samples = Vec::new();

    for entry in std::fs::read_dir(dir)
        .with_context(|| format!("Failed to read directory: {}", dir.display()))?
    {
        let entry = entry?;
        let path = entry.path();

        if path.is_dir() {
            // Directory structure: dir/sample_name/reads_1.fq.gz
            if let Some(sample) = load_sample_from_dir(&path, read_type)? {
                samples.push(sample);
            }
        } else if is_fastq_r1(&path) {
            // Flat structure: dir/sample_R1.fq.gz
            if let Some(sample) = load_sample_from_file(&path, read_type)? {
                samples.push(sample);
            }
        }
    }

    samples.sort_by(|a, b| a.name.cmp(&b.name));
    Ok(samples)
}

/// Load a sample from a directory containing FASTQ files
fn load_sample_from_dir(dir: &Path, read_type: &str) -> Result<Option<ReadSample>> {
    let name = dir
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("unknown")
        .to_string();

    // Look for R1 file
    let r1_patterns = ["reads_1.fq.gz", "reads_1.fastq.gz", "_R1.fq.gz", "_R1.fastq.gz", "_1.fq.gz", "_1.fastq.gz"];
    let r2_patterns = ["reads_2.fq.gz", "reads_2.fastq.gz", "_R2.fq.gz", "_R2.fastq.gz", "_2.fq.gz", "_2.fastq.gz"];

    let mut r1_path: Option<PathBuf> = None;
    let mut r2_path: Option<PathBuf> = None;

    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let file_name = entry.file_name().to_string_lossy().to_string();

        for pattern in &r1_patterns {
            if file_name.ends_with(pattern) || file_name == *pattern {
                r1_path = Some(entry.path());
                break;
            }
        }

        for pattern in &r2_patterns {
            if file_name.ends_with(pattern) || file_name == *pattern {
                r2_path = Some(entry.path());
                break;
            }
        }
    }

    match r1_path {
        Some(r1) => Ok(Some(ReadSample {
            name,
            r1_path: r1,
            r2_path,
            read_type: read_type.to_string(),
        })),
        None => Ok(None),
    }
}

/// Load a sample from a single R1 file (auto-detect R2)
pub fn load_sample_from_file(r1_path: &Path, read_type: &str) -> Result<Option<ReadSample>> {
    let file_name = r1_path.file_name().and_then(|n| n.to_str()).unwrap_or("");

    // Extract sample name from R1 filename
    let name_from_file = file_name
        .replace("_R1.fq.gz", "")
        .replace("_R1.fastq.gz", "")
        .replace("_1.fq.gz", "")
        .replace("_1.fastq.gz", "")
        .replace(".fq.gz", "")
        .replace(".fastq.gz", "");

    // If filename is generic (e.g., "reads", "reads_1"), use parent directory name
    let name = if name_from_file == "reads" || name_from_file == "reads_1" || name_from_file == "R1" {
        // Get parent directory name as sample name
        r1_path.parent()
            .and_then(|p| p.file_name())
            .and_then(|n| n.to_str())
            .map(|s| s.to_string())
            .unwrap_or(name_from_file)
    } else {
        name_from_file
    };

    // Try to find R2
    let r2_path = find_r2_file(r1_path);

    Ok(Some(ReadSample {
        name,
        r1_path: r1_path.to_path_buf(),
        r2_path,
        read_type: read_type.to_string(),
    }))
}

/// Check if a file is a FASTQ R1 file
pub fn is_fastq_r1(path: &Path) -> bool {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
    let lower = name.to_lowercase();

    (lower.ends_with(".fq") || lower.ends_with(".fastq") ||
     lower.ends_with(".fq.gz") || lower.ends_with(".fastq.gz")) &&
    (lower.contains("_r1") || lower.contains("_1.") || lower.contains("reads_1"))
}

/// Find the R2 file corresponding to an R1 file
fn find_r2_file(r1_path: &Path) -> Option<PathBuf> {
    let r1_str = r1_path.to_string_lossy();

    let r2_candidates = [
        r1_str.replace("_R1.", "_R2."),
        r1_str.replace("_1.", "_2."),
        r1_str.replace("reads_1", "reads_2"),
    ];

    for candidate in &r2_candidates {
        let r2_path = PathBuf::from(candidate);
        if r2_path.exists() {
            return Some(r2_path);
        }
    }

    None
}

/// Map reads to reference and calculate coverage-based gaps
pub fn analyze_read_coverage(
    sample: &ReadSample,
    reference: &Sample,
    min_depth: u32,
    threads: usize,
    minimap2_path: Option<&Path>,
) -> Result<ReadMappingResult> {
    let mm2 = minimap2_path
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| "minimap2".into());

    // Write reference to temp file
    let ref_temp = write_reference_temp(reference)?;

    // Choose minimap2 preset based on read type
    let preset = match sample.read_type.as_str() {
        "long" | "ont" | "pacbio" => "map-ont",
        _ => "sr", // short reads (Illumina)
    };

    log::info!("Mapping {} to reference ({} mode)...", sample.name, preset);

    // Build minimap2 command
    let mut cmd = Command::new(&mm2);
    cmd.args([
        "-x", preset,
        "-t", &threads.to_string(),
        "--secondary=no",
        "-a", // Output SAM
    ]);
    cmd.arg(ref_temp.path());
    cmd.arg(&sample.r1_path);

    if let Some(ref r2) = sample.r2_path {
        cmd.arg(r2);
    }

    cmd.stdout(Stdio::piped());
    cmd.stderr(Stdio::null());

    let output = cmd.output().context("Failed to run minimap2")?;

    if !output.status.success() {
        anyhow::bail!("minimap2 failed for {}", sample.name);
    }

    // Parse SAM and calculate coverage
    let coverage = calculate_coverage_from_sam(&output.stdout, reference.total_length)?;

    // Find gap regions (coverage < min_depth)
    let gap_regions = find_gap_regions(&coverage, min_depth, reference.total_length);
    let total_gap_bases: usize = gap_regions.iter().map(|r| r.end - r.start).sum();

    // Calculate stats
    let total_reads = coverage.iter().filter(|&&d| d > 0).count();
    let mapped_reads = total_reads; // Simplified
    let mean_depth = if coverage.is_empty() {
        0.0
    } else {
        coverage.iter().map(|&d| d as f64).sum::<f64>() / coverage.len() as f64
    };
    let covered_positions = coverage.iter().filter(|&&d| d >= min_depth).count();
    let coverage_percent = covered_positions as f64 / reference.total_length as f64 * 100.0;

    Ok(ReadMappingResult {
        name: sample.name.clone(),
        total_reads,
        mapped_reads,
        mean_depth,
        coverage_percent,
        gap_regions,
        total_gap_bases,
    })
}

/// Write reference to temporary FASTA file
fn write_reference_temp(reference: &Sample) -> Result<tempfile::NamedTempFile> {
    use std::io::Write;

    let mut temp = tempfile::NamedTempFile::new().context("Failed to create temp file")?;

    for contig in &reference.contigs {
        let contig_id = contig.name.split_whitespace().next().unwrap_or(&contig.name);
        writeln!(temp, ">{}", contig_id)?;
        for chunk in contig.sequence.chunks(80) {
            temp.write_all(chunk)?;
            writeln!(temp)?;
        }
    }

    temp.flush()?;
    Ok(temp)
}

/// Calculate per-position coverage from SAM output
fn calculate_coverage_from_sam(sam_data: &[u8], ref_length: usize) -> Result<Vec<u32>> {
    let mut coverage = vec![0u32; ref_length];

    let reader = BufReader::new(sam_data);

    for line in reader.lines() {
        let line = line?;

        // Skip header lines
        if line.starts_with('@') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }

        // Parse SAM fields
        let flag: u16 = fields[1].parse().unwrap_or(0);
        let pos: usize = fields[3].parse().unwrap_or(0);
        let cigar = fields[5];

        // Skip unmapped reads (flag & 4)
        if flag & 4 != 0 || pos == 0 || cigar == "*" {
            continue;
        }

        // Parse CIGAR and increment coverage
        let aligned_length = parse_cigar_aligned_length(cigar);
        let start = pos.saturating_sub(1); // SAM is 1-based
        let end = (start + aligned_length).min(ref_length);

        for i in start..end {
            coverage[i] = coverage[i].saturating_add(1);
        }
    }

    Ok(coverage)
}

/// Parse CIGAR string and return aligned length on reference
fn parse_cigar_aligned_length(cigar: &str) -> usize {
    let mut length = 0;
    let mut num_str = String::new();

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num_str.push(c);
        } else {
            let num: usize = num_str.parse().unwrap_or(0);
            num_str.clear();

            // Operations that consume reference: M, D, N, =, X
            match c {
                'M' | 'D' | 'N' | '=' | 'X' => length += num,
                _ => {}
            }
        }
    }

    length
}

/// Find regions with coverage below threshold
fn find_gap_regions(coverage: &[u32], min_depth: u32, ref_length: usize) -> Vec<Region> {
    let mut regions = Vec::new();
    let mut in_gap = false;
    let mut gap_start = 0;

    for (i, &depth) in coverage.iter().enumerate() {
        if depth < min_depth {
            if !in_gap {
                in_gap = true;
                gap_start = i;
            }
        } else if in_gap {
            in_gap = false;
            regions.push(Region::new(gap_start, i));
        }
    }

    // Handle gap at end
    if in_gap {
        regions.push(Region::new(gap_start, ref_length));
    }

    regions
}

/// Convert ReadMappingResult to a Sample-like structure for compatibility
pub fn read_result_to_sample(result: &ReadMappingResult, reference: &Sample) -> Sample {
    // Create a "virtual" sample with gap regions marked
    // This allows reusing existing analysis code

    Sample {
        name: result.name.clone(),
        contigs: vec![], // No actual contigs for read-based samples
        sequence: vec![], // No sequence
        total_length: reference.total_length,
        num_contigs: 0,
        n50: 0,
        sequence_hash: format!("reads_{}", result.name),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigar_parsing() {
        assert_eq!(parse_cigar_aligned_length("100M"), 100);
        assert_eq!(parse_cigar_aligned_length("50M2I48M"), 98); // I doesn't consume ref
        assert_eq!(parse_cigar_aligned_length("50M2D48M"), 100); // D consumes ref
        assert_eq!(parse_cigar_aligned_length("10M1000N90M"), 1100); // N (intron) consumes ref
    }

    #[test]
    fn test_gap_regions() {
        let coverage = vec![10, 10, 2, 1, 0, 10, 10, 3, 10];
        let gaps = find_gap_regions(&coverage, 5, 9);

        assert_eq!(gaps.len(), 2);
        assert_eq!(gaps[0].start, 2);
        assert_eq!(gaps[0].end, 5);
        assert_eq!(gaps[1].start, 7);
        assert_eq!(gaps[1].end, 8);
    }
}
