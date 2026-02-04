//! Parser for CFSAN SNP Pipeline `snplist.txt` format.
//!
//! Format: TSV with columns `CHROM\tPOS\tCOUNT\tsample1\tsample2\t...`
//! No per-sample alleles — only which samples have a SNP at each position.
//! This is also used as a fallback parser for unknown formats.

use anyhow::Context;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use super::{CoreSnpData, CoreSnpParser, CoreSnpPosition};

pub struct CfsanSnplistParser;

impl CoreSnpParser for CfsanSnplistParser {
    fn format_name(&self) -> &str {
        "snplist.txt"
    }

    fn can_parse(&self, path: &Path) -> bool {
        // Fallback parser — accepts any readable file
        path.exists()
    }

    fn parse(&self, path: &Path) -> anyhow::Result<CoreSnpData> {
        let file = std::fs::File::open(path)
            .with_context(|| format!("Failed to open core SNPs file: {}", path.display()))?;
        let reader = BufReader::new(file);

        let mut positions = Vec::new();

        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 4 {
                continue;
            }

            let pos: usize = parts[1].parse().unwrap_or(0);
            if pos == 0 { continue; }
            let pos = pos - 1; // Convert to 0-based

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

        // Collect all unique sample names to determine total sample count
        let all_samples: std::collections::HashSet<&str> = positions.iter()
            .flat_map(|p| p.samples_with_snp.iter().map(|s| s.as_str()))
            .collect();
        let total_samples = all_samples.len();

        // Discriminating: positions where only a subset of samples has the SNP.
        // If all samples appear -> all mutated the same way vs ref -> not discriminating.
        // Conservative underestimate: two samples could have different ALT alleles
        // but we count them as non-discriminating since we lack allele info.
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
}
