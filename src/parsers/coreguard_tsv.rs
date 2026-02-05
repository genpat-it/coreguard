//! Parser for CoreGuard TSV format.
//!
//! Format: TSV with header `CHR\tPOS\tREF\tsample1\tsample2\t...`
//! Each row has per-sample alleles; samples with allele != REF (and != N/-) have a SNP.
//!
//! This format is identical to Snippy's `core.tab` output.
//! For other pipelines (CFSAN, GATK, etc.), use `coreguard convert` to generate this format.

use anyhow::Context;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use super::{CoreSnpData, CoreSnpParser, CoreSnpPosition};

pub struct CoreGuardTsvParser;

impl CoreSnpParser for CoreGuardTsvParser {
    fn format_name(&self) -> &str {
        "coreguard_tsv"
    }

    fn can_parse(&self, path: &Path) -> bool {
        let Ok(content) = std::fs::read_to_string(path) else { return false };
        let first_line = content.lines().next().unwrap_or("");
        first_line.starts_with("CHR\t") || first_line.starts_with("CHR ")
    }

    fn parse(&self, path: &Path) -> anyhow::Result<CoreSnpData> {
        let file = std::fs::File::open(path)
            .with_context(|| format!("Failed to open core SNPs file: {}", path.display()))?;
        let reader = BufReader::new(file);

        let mut positions = Vec::new();
        let mut sample_names: Vec<String> = Vec::new();

        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split('\t').collect();

            if line.starts_with("CHR\t") {
                // Header line — extract sample names
                sample_names = parts.iter()
                    .skip(3) // Skip CHR, POS, REF
                    .map(|s| s.trim_end_matches("_snippy").to_string())
                    .collect();
                continue;
            }

            if parts.len() < 4 || sample_names.is_empty() {
                continue;
            }

            let pos: usize = parts[1].parse().unwrap_or(0);
            if pos == 0 { continue; }
            let pos = pos - 1; // Convert to 0-based

            let ref_allele = parts[2].to_string();

            let mut alleles = HashMap::new();
            let mut samples_with_snp = Vec::new();

            for (i, sample) in sample_names.iter().enumerate() {
                if let Some(allele) = parts.get(3 + i) {
                    let allele = allele.to_string();
                    // A sample has a SNP if its allele differs from ref and isn't N or -
                    if allele != ref_allele && allele != "N" && allele != "-" {
                        samples_with_snp.push(sample.clone());
                    }
                    alleles.insert(sample.clone(), allele);
                }
            }

            positions.push(CoreSnpPosition {
                pos,
                ref_allele: Some(ref_allele),
                alleles,
                samples_with_snp,
            });
        }

        // Count discriminating positions: at least 2 samples with different valid alleles
        let discriminating_count = positions.iter().filter(|p| {
            let valid: Vec<&str> = p.alleles.values()
                .map(|a| a.as_str())
                .filter(|a| *a != "N" && *a != "-")
                .collect();
            if valid.len() < 2 { return false; }
            let first = valid[0];
            valid.iter().any(|a| *a != first)
        }).count();

        log::info!(
            "Parsed {} core SNP positions from CoreGuard TSV (with alleles, {} discriminating)",
            positions.len(), discriminating_count
        );

        Ok(CoreSnpData {
            positions,
            has_alleles: true,
            discriminating_count,
        })
    }
}
