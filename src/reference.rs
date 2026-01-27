//! Reference-based alignment analysis using minimap2
//!
//! This module performs O(N) alignment of samples to a reference genome.
//! Much faster than pairwise O(N²) - the primary analysis for --fast mode.

use crate::fasta::Sample;
use crate::gaps::{calculate_unaligned_regions, ReferenceGaps, Region};
use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};
use tempfile::NamedTempFile;

/// Result of aligning a sample to the reference
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ReferenceAlignment {
    /// Sample name
    pub sample_name: String,
    /// Reference name
    pub reference_name: String,
    /// Sequence identity (0.0 - 1.0)
    pub identity: f64,
    /// Fraction of sample aligned to reference
    pub sample_coverage: f64,
    /// Fraction of reference covered by sample
    pub reference_coverage: f64,
    /// Number of alignment blocks
    pub num_alignments: usize,
    /// Unaligned bases in sample (sample-specific sequence)
    pub sample_unaligned: usize,
    /// Uncovered bases in reference (missing in sample)
    pub reference_uncovered: usize,
    /// Total gap bases in alignments
    pub gap_bases: usize,
    /// Regions of sample not aligned to reference (coordinates)
    pub sample_unaligned_regions: Vec<Region>,
    /// Regions of reference not covered by sample (coordinates)
    pub reference_uncovered_regions: Vec<Region>,
}

/// Container for all reference alignment results
#[derive(Debug, Clone)]
pub struct ReferenceResults {
    /// Reference sample info
    pub reference_name: String,
    pub reference_length: usize,
    /// All alignments to reference
    pub alignments: Vec<ReferenceAlignment>,
}

impl ReferenceResults {
    /// Get alignment for a specific sample
    pub fn get(&self, sample_name: &str) -> Option<&ReferenceAlignment> {
        self.alignments.iter().find(|a| a.sample_name == sample_name)
    }

    /// Average identity across all samples
    pub fn avg_identity(&self) -> f64 {
        if self.alignments.is_empty() {
            return 0.0;
        }
        self.alignments.iter().map(|a| a.identity).sum::<f64>() / self.alignments.len() as f64
    }

    /// Average sample coverage
    pub fn avg_sample_coverage(&self) -> f64 {
        if self.alignments.is_empty() {
            return 0.0;
        }
        self.alignments.iter().map(|a| a.sample_coverage).sum::<f64>() / self.alignments.len() as f64
    }

    /// Average reference coverage
    pub fn avg_reference_coverage(&self) -> f64 {
        if self.alignments.is_empty() {
            return 0.0;
        }
        self.alignments.iter().map(|a| a.reference_coverage).sum::<f64>() / self.alignments.len() as f64
    }

    /// Get all gap regions for export
    pub fn get_all_gaps(&self) -> Vec<ReferenceGaps> {
        self.alignments
            .iter()
            .map(|a| ReferenceGaps {
                sample: a.sample_name.clone(),
                sample_unaligned: a.sample_unaligned_regions.clone(),
                reference_uncovered: a.reference_uncovered_regions.clone(),
            })
            .collect()
    }
}

/// PAF record for parsing minimap2 output
#[derive(Debug)]
#[allow(dead_code)]
struct PafRecord {
    query_name: String,
    query_len: usize,
    query_start: usize,
    query_end: usize,
    target_name: String,
    target_len: usize,
    target_start: usize,
    target_end: usize,
    num_matches: usize,
    block_len: usize,
}

impl PafRecord {
    fn parse(line: &str) -> Option<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            return None;
        }

        Some(PafRecord {
            query_name: fields[0].to_string(),
            query_len: fields[1].parse().ok()?,
            query_start: fields[2].parse().ok()?,
            query_end: fields[3].parse().ok()?,
            target_name: fields[5].to_string(),
            target_len: fields[6].parse().ok()?,
            target_start: fields[7].parse().ok()?,
            target_end: fields[8].parse().ok()?,
            num_matches: fields[9].parse().ok()?,
            block_len: fields[10].parse().ok()?,
        })
    }

    fn query_aligned(&self) -> usize {
        self.query_end - self.query_start
    }

    fn target_aligned(&self) -> usize {
        self.target_end - self.target_start
    }

    fn gap_bases(&self) -> usize {
        self.block_len.saturating_sub(self.num_matches)
    }
}

/// Align all samples to a reference genome using minimap2
/// This is O(N) complexity - much faster than pairwise O(N²)
pub fn analyze_vs_reference(
    reference: &Sample,
    samples: &[Sample],
    minimap2_path: Option<&Path>,
    _threads: usize,  // Reserved for future batch optimization
) -> Result<ReferenceResults> {
    let mm2 = minimap2_path
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| "minimap2".into());

    log::info!(
        "Aligning {} samples to reference {} ({} bp)...",
        samples.len(),
        reference.name,
        reference.total_length
    );

    // Write reference to temp file
    let ref_temp = write_sample_to_temp(reference)?;

    // Write all samples to temp files
    let sample_temps: Vec<NamedTempFile> = samples
        .iter()
        .map(|s| write_sample_to_temp(s))
        .collect::<Result<Vec<_>>>()?;

    // Setup progress bar
    let pb = ProgressBar::new(samples.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} samples ({percent}%) ETA: {eta}")
            .unwrap()
            .progress_chars("=>-"),
    );
    pb.set_message("Reference alignments");

    let counter = AtomicUsize::new(0);

    // Align each sample to reference in parallel
    let alignments: Vec<ReferenceAlignment> = samples
        .par_iter()
        .zip(sample_temps.par_iter())
        .map(|(sample, sample_temp)| {
            let result = run_reference_alignment(&mm2, reference, sample, ref_temp.path(), sample_temp.path());
            let done = counter.fetch_add(1, Ordering::Relaxed) + 1;
            pb.set_position(done as u64);
            result
        })
        .collect::<Result<Vec<_>>>()?;

    pb.finish_with_message("Reference alignments complete");

    Ok(ReferenceResults {
        reference_name: reference.name.clone(),
        reference_length: reference.total_length,
        alignments,
    })
}

/// Build contig offset map: contig_id -> (offset_in_sample, contig_len)
fn build_contig_offsets(sample: &Sample) -> std::collections::HashMap<String, (usize, usize)> {
    let mut offsets = std::collections::HashMap::new();
    let mut offset = 0usize;
    for contig in &sample.contigs {
        let contig_id = contig.name.split_whitespace().next().unwrap_or(&contig.name);
        offsets.insert(contig_id.to_string(), (offset, contig.sequence.len()));
        offset += contig.sequence.len();
    }
    offsets
}

/// Write sample contigs to temporary FASTA file (each contig as separate record)
fn write_sample_to_temp(sample: &Sample) -> Result<NamedTempFile> {
    use std::io::Write;

    let mut temp = NamedTempFile::new().context("Failed to create temp file")?;

    // Write each contig as a separate FASTA record
    // Use only the contig ID (before first space) to match minimap2 behavior
    for contig in &sample.contigs {
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

/// Run minimap2 alignment of sample to reference
fn run_reference_alignment(
    mm2_path: &Path,
    reference: &Sample,
    sample: &Sample,
    ref_file: &Path,
    sample_file: &Path,
) -> Result<ReferenceAlignment> {
    let output = Command::new(mm2_path)
        .args([
            "-x", "asm5",
            "-c",
            "--secondary=no",
            "-t", "1",
        ])
        .arg(ref_file)
        .arg(sample_file)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .output()
        .context("Failed to run minimap2")?;

    if !output.status.success() {
        anyhow::bail!(
            "minimap2 failed for {} vs reference {}",
            sample.name,
            reference.name
        );
    }

    // Build contig offset maps for coordinate conversion
    let ref_offsets = build_contig_offsets(reference);
    let sample_offsets = build_contig_offsets(sample);

    // Parse PAF output
    let reader = BufReader::new(&output.stdout[..]);
    let records: Vec<PafRecord> = reader
        .lines()
        .filter_map(|l| l.ok())
        .filter_map(|l| PafRecord::parse(&l))
        .collect();

    let num_alignments = records.len();

    if num_alignments == 0 {
        return Ok(ReferenceAlignment {
            sample_name: sample.name.clone(),
            reference_name: reference.name.clone(),
            identity: 0.0,
            sample_coverage: 0.0,
            reference_coverage: 0.0,
            num_alignments: 0,
            sample_unaligned: sample.total_length,
            reference_uncovered: reference.total_length,
            gap_bases: 0,
            sample_unaligned_regions: vec![Region::new(0, sample.total_length)],
            reference_uncovered_regions: vec![Region::new(0, reference.total_length)],
        });
    }

    // Calculate metrics
    let total_matches: usize = records.iter().map(|r| r.num_matches).sum();
    let total_block_len: usize = records.iter().map(|r| r.block_len).sum();
    let identity = if total_block_len > 0 {
        total_matches as f64 / total_block_len as f64
    } else {
        0.0
    };

    // Collect aligned intervals - convert contig coordinates to concatenated coordinates
    let mut sample_intervals: Vec<(usize, usize)> = Vec::new();
    let mut reference_intervals: Vec<(usize, usize)> = Vec::new();

    for record in &records {
        // Convert sample (query) coordinates
        if let Some(&(sample_offset, _)) = sample_offsets.get(&record.query_name) {
            sample_intervals.push((
                sample_offset + record.query_start,
                sample_offset + record.query_end,
            ));
        }

        // Convert reference (target) coordinates
        if let Some(&(ref_offset, _)) = ref_offsets.get(&record.target_name) {
            reference_intervals.push((
                ref_offset + record.target_start,
                ref_offset + record.target_end,
            ));
        }
    }

    let sample_aligned: usize = records.iter().map(|r| r.query_aligned()).sum();
    let reference_aligned: usize = records.iter().map(|r| r.target_aligned()).sum();

    let sample_coverage = (sample_aligned as f64 / sample.total_length as f64).min(1.0);
    let reference_coverage = (reference_aligned as f64 / reference.total_length as f64).min(1.0);

    let gap_bases: usize = records.iter().map(|r| r.gap_bases()).sum();

    let sample_unaligned = sample.total_length.saturating_sub(sample_aligned);
    let reference_uncovered = reference.total_length.saturating_sub(reference_aligned);

    // Calculate gap regions (coordinates)
    let sample_unaligned_regions = calculate_unaligned_regions(sample.total_length, &sample_intervals);
    let reference_uncovered_regions = calculate_unaligned_regions(reference.total_length, &reference_intervals);

    Ok(ReferenceAlignment {
        sample_name: sample.name.clone(),
        reference_name: reference.name.clone(),
        identity,
        sample_coverage,
        reference_coverage,
        num_alignments,
        sample_unaligned,
        reference_uncovered,
        gap_bases,
        sample_unaligned_regions,
        reference_uncovered_regions,
    })
}

/// Build ReferenceResults from FASTQ coverage analysis
/// This is used in FASTQ mode where we have coverage-based gaps instead of assembly alignments
pub fn build_from_fastq_coverage(
    reference: &Sample,
    samples: &[Sample],
    fastq_results: &[crate::reads::ReadMappingResult],
) -> ReferenceResults {
    // Build alignments from FASTQ results
    let alignments: Vec<ReferenceAlignment> = samples
        .iter()
        .zip(fastq_results.iter())
        .map(|(sample, result)| {
            // Calculate coverage from gap regions
            let covered_bases = reference.total_length.saturating_sub(result.total_gap_bases);
            let reference_coverage = covered_bases as f64 / reference.total_length as f64;

            ReferenceAlignment {
                sample_name: sample.name.clone(),
                reference_name: reference.name.clone(),
                identity: 1.0, // Not applicable for read mapping
                sample_coverage: result.coverage_percent / 100.0,
                reference_coverage,
                num_alignments: 1, // Treat as single alignment
                sample_unaligned: 0, // Not applicable for read mapping
                reference_uncovered: result.total_gap_bases,
                gap_bases: result.total_gap_bases,
                sample_unaligned_regions: vec![], // Reads don't have "unaligned" in the same sense
                reference_uncovered_regions: result.gap_regions.clone(),
            }
        })
        .collect();

    log::info!(
        "Built reference results from {} FASTQ samples",
        alignments.len()
    );

    for (sample, result) in samples.iter().zip(fastq_results.iter()) {
        log::debug!(
            "  {} - Coverage: {:.1}%, Depth: {:.1}x, Gap regions: {}",
            sample.name,
            result.coverage_percent,
            result.mean_depth,
            result.gap_regions.len()
        );
    }

    ReferenceResults {
        reference_name: reference.name.clone(),
        reference_length: reference.total_length,
        alignments,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_paf_parse() {
        let line = "sample\t3000000\t100\t2999000\t+\tref\t3000000\t0\t2999000\t2990000\t3000000\t60";
        let record = PafRecord::parse(line).unwrap();

        assert_eq!(record.query_len, 3000000);
        assert_eq!(record.num_matches, 2990000);
        assert_eq!(record.query_aligned(), 2998900);
    }
}
