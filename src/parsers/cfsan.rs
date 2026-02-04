//! Parser for CFSAN SNP Pipeline output.
//!
//! Primary: reads `snpma.fasta` (SNP matrix alignment) + `referenceSNP.fasta` (reference alleles)
//! + `snplist.txt` (genomic positions) from the same directory.
//!
//! `snpma.fasta` provides per-sample alleles at each SNP position, with `-` for gaps.
//! This gives us full allele information and accurate gap handling.
//!
//! Fallback: if `snpma.fasta` is not found, falls back to `snplist.txt`-only parsing
//! (position-only, no alleles, no gap info).

use anyhow::{Context, bail};
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use super::{CoreSnpData, CoreSnpParser, CoreSnpPosition};

pub struct CfsanSnplistParser;

impl CoreSnpParser for CfsanSnplistParser {
    fn format_name(&self) -> &str {
        "cfsan"
    }

    fn can_parse(&self, path: &Path) -> bool {
        // Accept snplist.txt, snpma.fasta, or any file in a CFSAN output directory
        path.exists()
    }

    fn parse(&self, path: &Path) -> anyhow::Result<CoreSnpData> {
        let dir = path.parent().unwrap_or(Path::new("."));

        // Try to find all CFSAN files in the same directory
        let snplist_path = if path.file_name().map_or(false, |f| f == "snplist.txt") {
            path.to_path_buf()
        } else {
            dir.join("snplist.txt")
        };
        let snpma_path = dir.join("snpma.fasta");
        let ref_snp_path = dir.join("referenceSNP.fasta");

        if snpma_path.exists() && ref_snp_path.exists() && snplist_path.exists() {
            self.parse_with_snpma(&snplist_path, &snpma_path, &ref_snp_path)
        } else if snplist_path.exists() {
            log::warn!("snpma.fasta or referenceSNP.fasta not found in {}, falling back to snplist.txt-only parsing (no alleles, no gap info)", dir.display());
            self.parse_snplist_only(&snplist_path)
        } else {
            bail!("No CFSAN output files found in {}", dir.display())
        }
    }
}

impl CfsanSnplistParser {
    /// Full parsing using snpma.fasta + referenceSNP.fasta + snplist.txt
    fn parse_with_snpma(
        &self,
        snplist_path: &Path,
        snpma_path: &Path,
        ref_snp_path: &Path,
    ) -> anyhow::Result<CoreSnpData> {
        // 1. Read genomic positions from snplist.txt
        let genomic_positions = Self::read_snplist_positions(snplist_path)?;

        // 2. Read reference alleles from referenceSNP.fasta
        let ref_seq = Self::read_single_fasta(ref_snp_path)
            .with_context(|| format!("Failed to read {}", ref_snp_path.display()))?;

        // 3. Read sample sequences from snpma.fasta
        let sample_seqs = Self::read_multi_fasta(snpma_path)
            .with_context(|| format!("Failed to read {}", snpma_path.display()))?;

        if ref_seq.len() != genomic_positions.len() {
            bail!(
                "Length mismatch: referenceSNP.fasta has {} positions, snplist.txt has {}",
                ref_seq.len(), genomic_positions.len()
            );
        }

        let ref_bytes = ref_seq.as_bytes();
        let sample_names: Vec<&String> = sample_seqs.keys().collect();

        let mut positions = Vec::new();
        let mut discriminating_count = 0;

        for (col_idx, &genomic_pos) in genomic_positions.iter().enumerate() {
            let ref_allele = (ref_bytes[col_idx] as char).to_string();
            let mut alleles = HashMap::new();
            let mut samples_with_snp = Vec::new();
            let mut non_gap_alleles: Vec<u8> = Vec::new();

            for &sample in &sample_names {
                let seq_bytes = sample_seqs[sample].as_bytes();
                if col_idx >= seq_bytes.len() { continue; }

                let base = seq_bytes[col_idx];
                if base == b'-' {
                    // Gap — record as "-" in alleles but don't count as SNP
                    alleles.insert(sample.clone(), "-".to_string());
                    continue;
                }

                let allele = (base as char).to_string();
                alleles.insert(sample.clone(), allele);

                if base != ref_bytes[col_idx] {
                    samples_with_snp.push(sample.clone());
                }

                non_gap_alleles.push(base);
            }

            // Discriminating: at least 2 non-gap samples with different alleles
            if non_gap_alleles.len() >= 2 {
                let first = non_gap_alleles[0];
                if non_gap_alleles.iter().any(|&a| a != first) {
                    discriminating_count += 1;
                }
            }

            positions.push(CoreSnpPosition {
                pos: genomic_pos,
                ref_allele: Some(ref_allele),
                alleles,
                samples_with_snp,
            });
        }

        log::info!(
            "Parsed {} core SNP positions from CFSAN snpma.fasta (with alleles, {} discriminating, {} samples)",
            positions.len(), discriminating_count, sample_names.len()
        );

        Ok(CoreSnpData {
            positions,
            has_alleles: true,
            discriminating_count,
        })
    }

    /// Fallback: parse snplist.txt only (no alleles, no gap info)
    fn parse_snplist_only(&self, path: &Path) -> anyhow::Result<CoreSnpData> {
        let file = std::fs::File::open(path)
            .with_context(|| format!("Failed to open core SNPs file: {}", path.display()))?;
        let reader = BufReader::new(file);

        let mut positions = Vec::new();

        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 4 { continue; }

            let pos: usize = parts[1].parse().unwrap_or(0);
            if pos == 0 { continue; }
            let pos = pos - 1;

            let samples_with_snp: Vec<String> = parts.iter()
                .skip(3)
                .map(|s| s.to_string())
                .collect();

            positions.push(CoreSnpPosition {
                pos,
                ref_allele: None,
                alleles: HashMap::new(),
                samples_with_snp,
            });
        }

        let all_samples: std::collections::HashSet<&str> = positions.iter()
            .flat_map(|p| p.samples_with_snp.iter().map(|s| s.as_str()))
            .collect();
        let total_samples = all_samples.len();

        let discriminating_count = positions.iter().filter(|p| {
            let n = p.samples_with_snp.len();
            n > 0 && n < total_samples
        }).count();

        log::info!(
            "Parsed {} core SNP positions from snplist.txt (no alleles, {} discriminating, {} total samples)",
            positions.len(), discriminating_count, total_samples
        );

        Ok(CoreSnpData {
            positions,
            has_alleles: false,
            discriminating_count,
        })
    }

    /// Read genomic positions from snplist.txt (column 2, 1-based → 0-based)
    fn read_snplist_positions(path: &Path) -> anyhow::Result<Vec<usize>> {
        let file = std::fs::File::open(path)
            .with_context(|| format!("Failed to open {}", path.display()))?;
        let reader = BufReader::new(file);
        let mut positions = Vec::new();

        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 2 { continue; }
            let pos: usize = parts[1].parse().unwrap_or(0);
            if pos == 0 { continue; }
            positions.push(pos - 1);
        }
        Ok(positions)
    }

    /// Read a single-sequence FASTA file and return the sequence
    fn read_single_fasta(path: &Path) -> anyhow::Result<String> {
        let file = std::fs::File::open(path)?;
        let reader = BufReader::new(file);
        let mut seq = String::new();

        for line in reader.lines() {
            let line = line?;
            if !line.starts_with('>') {
                seq.push_str(line.trim());
            }
        }
        Ok(seq)
    }

    /// Read a multi-sequence FASTA file and return sample_name -> sequence
    fn read_multi_fasta(path: &Path) -> anyhow::Result<HashMap<String, String>> {
        let file = std::fs::File::open(path)?;
        let reader = BufReader::new(file);
        let mut seqs: HashMap<String, String> = HashMap::new();
        let mut current: Option<String> = None;

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('>') {
                current = Some(line[1..].trim().to_string());
                seqs.entry(current.clone().unwrap()).or_default();
            } else if let Some(ref name) = current {
                seqs.get_mut(name).unwrap().push_str(line.trim());
            }
        }
        Ok(seqs)
    }
}
