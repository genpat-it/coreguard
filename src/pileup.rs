//! BAM pileup module for extracting consensus bases at specific positions
//!
//! Uses noodles to read BAM files and determine the consensus base
//! at each position, enabling accurate SNP distance calculations.

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use noodles::bam;
use noodles::bam::bai;
use noodles::sam::alignment::Record;

/// Result of reading a base from BAM at a specific position
#[derive(Debug, Clone, PartialEq)]
pub enum BaseCall {
    /// Clear base call with sufficient coverage
    Base(char),
    /// No coverage at this position (gap)
    Gap,
    /// Ambiguous call (mixed bases)
    Ambiguous,
}

impl BaseCall {
    /// Convert to character representation
    pub fn to_char(&self) -> char {
        match self {
            BaseCall::Base(c) => *c,
            BaseCall::Gap => '-',
            BaseCall::Ambiguous => 'N',
        }
    }
}

/// Result of reading a base from BAM with additional statistics
#[derive(Debug, Clone)]
pub struct BaseCallWithStats {
    /// The base call result
    pub call: BaseCall,
    /// Read depth at this position
    pub depth: u32,
    /// Consensus percentage (0.0 - 1.0): fraction of reads agreeing on the called base
    pub consensus: f64,
}

/// Batch read bases at multiple positions from a BAM file
///
/// This is the main function to use - it efficiently reads all positions in one pass.
/// Positions should be 0-based.
///
/// Returns a map of position -> BaseCall
pub fn get_bases_at_positions(
    bam_path: &Path,
    chrom: &str,
    positions: &[u32],  // 0-based positions
    min_depth: u32,
) -> Result<HashMap<u32, BaseCall>> {
    if positions.is_empty() {
        return Ok(HashMap::new());
    }

    // Create set of positions we care about for fast lookup
    let position_set: std::collections::HashSet<u32> = positions.iter().copied().collect();

    // Find the range we need to scan
    let min_pos = *positions.iter().min().unwrap();
    let max_pos = *positions.iter().max().unwrap();

    // Open BAM file
    let mut reader = File::open(bam_path)
        .map(BufReader::new)
        .map(bam::io::Reader::new)
        .with_context(|| format!("Failed to open BAM: {}", bam_path.display()))?;

    // Read header
    let header = reader.read_header()?;

    // Try to open index for random access
    let index_path = bam_path.with_extension("bam.bai");
    let index_path_alt = {
        let mut p = bam_path.as_os_str().to_owned();
        p.push(".bai");
        PathBuf::from(p)
    };

    // Count bases at each position
    let mut position_counts: HashMap<u32, HashMap<char, u32>> = HashMap::new();

    // Try indexed access first
    if index_path.exists() || index_path_alt.exists() {
        let idx_path = if index_path.exists() { &index_path } else { &index_path_alt };

        let index = bai::read(idx_path)
            .with_context(|| format!("Failed to read BAM index: {}", idx_path.display()))?;

        // Get reference sequence ID
        let header_refs = header.reference_sequences();
        let ref_id = header_refs.get_index_of(chrom.as_bytes())
            .ok_or_else(|| anyhow::anyhow!("Chromosome {} not found in BAM header", chrom))?;

        // Create region to query
        let start = noodles::core::Position::try_from((min_pos + 1) as usize)?; // 1-based
        let end = noodles::core::Position::try_from((max_pos + 2) as usize)?; // 1-based, inclusive

        let region = noodles::core::Region::new(chrom, start..=end);

        // Query the region
        let mut query = reader.query(&header, &index, &region)?;

        while let Some(result) = query.next() {
            let record = result?;
            process_record(&record, &header, &position_set, &mut position_counts)?;
        }
    } else {
        // Fall back to full scan (slower)
        log::warn!("No BAM index found for {}, doing full scan", bam_path.display());

        for result in reader.records() {
            let record = result?;

            // Check if this read overlaps our region of interest
            if let (Some(ref_name), Some(start)) = (record.reference_sequence(&header), record.alignment_start()) {
                let (name, _) = ref_name?;
                let name_bytes: &[u8] = name.as_ref();
                if name_bytes == chrom.as_bytes() {
                    let start_pos = start?.get() as u32 - 1; // Convert to 0-based
                    let read_len = record.sequence().len() as u32;
                    let end_pos = start_pos + read_len;

                    // Check if read overlaps our position range
                    if end_pos >= min_pos && start_pos <= max_pos {
                        process_record(&record, &header, &position_set, &mut position_counts)?;
                    }
                }
            }
        }
    }

    // Convert counts to base calls
    let mut results = HashMap::new();

    for &pos in positions {
        let call = if let Some(counts) = position_counts.get(&pos) {
            let total_depth: u32 = counts.values().sum();

            if total_depth < min_depth {
                BaseCall::Gap
            } else {
                let (best_base, best_count) = counts.iter()
                    .max_by_key(|(_, count)| *count)
                    .map(|(b, c)| (*b, *c))
                    .unwrap_or(('N', 0));

                if (best_count as f64 / total_depth as f64) < 0.8 {
                    BaseCall::Ambiguous
                } else {
                    BaseCall::Base(best_base)
                }
            }
        } else {
            BaseCall::Gap  // No coverage at this position
        };

        results.insert(pos, call);
    }

    Ok(results)
}

/// Batch read bases at multiple positions from a BAM file, returning depth and consensus stats
///
/// This is similar to `get_bases_at_positions` but returns additional statistics.
/// Positions should be 0-based.
///
/// Returns a map of position -> BaseCallWithStats
pub fn get_bases_with_stats(
    bam_path: &Path,
    chrom: &str,
    positions: &[u32],  // 0-based positions
    min_depth: u32,
) -> Result<HashMap<u32, BaseCallWithStats>> {
    if positions.is_empty() {
        return Ok(HashMap::new());
    }

    // Create set of positions we care about for fast lookup
    let position_set: std::collections::HashSet<u32> = positions.iter().copied().collect();

    // Find the range we need to scan
    let min_pos = *positions.iter().min().unwrap();
    let max_pos = *positions.iter().max().unwrap();

    // Open BAM file
    let mut reader = File::open(bam_path)
        .map(BufReader::new)
        .map(bam::io::Reader::new)
        .with_context(|| format!("Failed to open BAM: {}", bam_path.display()))?;

    // Read header
    let header = reader.read_header()?;

    // Try to open index for random access
    let index_path = bam_path.with_extension("bam.bai");
    let index_path_alt = {
        let mut p = bam_path.as_os_str().to_owned();
        p.push(".bai");
        PathBuf::from(p)
    };

    // Count bases at each position
    let mut position_counts: HashMap<u32, HashMap<char, u32>> = HashMap::new();

    // Try indexed access first
    if index_path.exists() || index_path_alt.exists() {
        let idx_path = if index_path.exists() { &index_path } else { &index_path_alt };

        let index = bai::read(idx_path)
            .with_context(|| format!("Failed to read BAM index: {}", idx_path.display()))?;

        // Get reference sequence ID
        let header_refs = header.reference_sequences();
        let _ref_id = header_refs.get_index_of(chrom.as_bytes())
            .ok_or_else(|| anyhow::anyhow!("Chromosome {} not found in BAM header", chrom))?;

        // Create region to query
        let start = noodles::core::Position::try_from((min_pos + 1) as usize)?; // 1-based
        let end = noodles::core::Position::try_from((max_pos + 2) as usize)?; // 1-based, inclusive

        let region = noodles::core::Region::new(chrom, start..=end);

        // Query the region
        let mut query = reader.query(&header, &index, &region)?;

        while let Some(result) = query.next() {
            let record = result?;
            process_record(&record, &header, &position_set, &mut position_counts)?;
        }
    } else {
        // Fall back to full scan (slower)
        log::warn!("No BAM index found for {}, doing full scan", bam_path.display());

        for result in reader.records() {
            let record = result?;

            // Check if this read overlaps our region of interest
            if let (Some(ref_name), Some(start)) = (record.reference_sequence(&header), record.alignment_start()) {
                let (name, _) = ref_name?;
                let name_bytes: &[u8] = name.as_ref();
                if name_bytes == chrom.as_bytes() {
                    let start_pos = start?.get() as u32 - 1; // Convert to 0-based
                    let read_len = record.sequence().len() as u32;
                    let end_pos = start_pos + read_len;

                    // Check if read overlaps our position range
                    if end_pos >= min_pos && start_pos <= max_pos {
                        process_record(&record, &header, &position_set, &mut position_counts)?;
                    }
                }
            }
        }
    }

    // Convert counts to base calls with stats
    let mut results = HashMap::new();

    for &pos in positions {
        let stats = if let Some(counts) = position_counts.get(&pos) {
            let total_depth: u32 = counts.values().sum();

            if total_depth < min_depth {
                BaseCallWithStats {
                    call: BaseCall::Gap,
                    depth: total_depth,
                    consensus: 0.0,
                }
            } else {
                let (best_base, best_count) = counts.iter()
                    .max_by_key(|(_, count)| *count)
                    .map(|(b, c)| (*b, *c))
                    .unwrap_or(('N', 0));

                let consensus = best_count as f64 / total_depth as f64;

                if consensus < 0.8 {
                    BaseCallWithStats {
                        call: BaseCall::Ambiguous,
                        depth: total_depth,
                        consensus,
                    }
                } else {
                    BaseCallWithStats {
                        call: BaseCall::Base(best_base),
                        depth: total_depth,
                        consensus,
                    }
                }
            }
        } else {
            BaseCallWithStats {
                call: BaseCall::Gap,
                depth: 0,
                consensus: 0.0,
            }
        };

        results.insert(pos, stats);
    }

    Ok(results)
}

/// Process a single BAM record, extracting bases at positions of interest
fn process_record<R: Record>(
    record: &R,
    header: &noodles::sam::Header,
    positions: &std::collections::HashSet<u32>,
    counts: &mut HashMap<u32, HashMap<char, u32>>,
) -> Result<()> {
    // Get alignment start (1-based in SAM/BAM)
    let Some(start) = record.alignment_start() else {
        return Ok(());
    };
    let start_pos = start?.get() as u32 - 1; // Convert to 0-based

    // Get sequence
    let seq = record.sequence();
    let seq_len = seq.len();

    // Get CIGAR to properly map reference positions to query positions
    let cigar = record.cigar();

    // Walk through the alignment
    let mut ref_pos = start_pos;
    let mut query_pos: usize = 0;

    for op in cigar.iter() {
        let op = op?;
        let len = op.len();

        use noodles::sam::alignment::record::cigar::op::Kind;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                // These consume both reference and query
                for i in 0..len {
                    let rp = ref_pos + i as u32;
                    let qp = query_pos + i;

                    if positions.contains(&rp) && qp < seq_len {
                        let base = seq.get(qp).map(|b| b as char).unwrap_or('N');
                        if base != 'N' {
                            let entry = counts.entry(rp).or_insert_with(HashMap::new);
                            *entry.entry(base).or_insert(0) += 1;
                        }
                    }
                }
                ref_pos += len as u32;
                query_pos += len;
            }
            Kind::Insertion | Kind::SoftClip => {
                // Consumes query only
                query_pos += len;
            }
            Kind::Deletion | Kind::Skip => {
                // Consumes reference only
                ref_pos += len as u32;
            }
            Kind::HardClip | Kind::Pad => {
                // Consumes neither
            }
        }
    }

    Ok(())
}

use std::path::PathBuf;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_call_to_char() {
        assert_eq!(BaseCall::Base('A').to_char(), 'A');
        assert_eq!(BaseCall::Gap.to_char(), '-');
        assert_eq!(BaseCall::Ambiguous.to_char(), 'N');
    }

    #[test]
    fn test_empty_positions() {
        let result = get_bases_at_positions(
            Path::new("/nonexistent.bam"),
            "chr1",
            &[],
            10,
        );
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }
}
