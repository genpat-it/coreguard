//! VCF file parsing for variant comparison
//!
//! Parses VCF files from Snippy, CFSAN, and other SNP pipelines
//! to enable direct comparison of variant calls.

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::path::Path;

/// A single variant call from a VCF file
#[derive(Debug, Clone)]
pub struct VariantCall {
    /// Chromosome/contig name
    pub chrom: String,
    /// 1-based position
    pub pos: usize,
    /// Reference allele
    pub ref_allele: String,
    /// Alternate allele(s)
    pub alt_allele: String,
    /// Quality score (QUAL column)
    pub qual: f64,
    /// Filter status (PASS, or filter name)
    pub filter: String,
    /// Read depth (DP from INFO or FORMAT)
    pub depth: Option<u32>,
    /// Allele frequency (AF from INFO or FORMAT)
    pub allele_freq: Option<f64>,
    /// Genotype (GT from FORMAT)
    pub genotype: Option<String>,
    /// Raw INFO field for additional parsing
    pub info: String,
    /// Sample name (extracted from filename or header)
    pub sample: String,
    // === Additional quality metrics ===
    /// Reference allele observations (RO - Snippy)
    pub ref_obs: Option<u32>,
    /// Alternate allele observations (AO - Snippy)
    pub alt_obs: Option<u32>,
    /// Quality of alternate allele reads (QA - Snippy)
    pub alt_qual: Option<u32>,
    /// Genotype quality (GQ - CFSAN)
    pub genotype_qual: Option<u32>,
    /// Average base quality (ABQ - CFSAN)
    pub avg_base_qual: Option<u32>,
    /// P-value (PVAL - CFSAN)
    pub pvalue: Option<f64>,
}

impl VariantCall {
    /// Check if this is a SNP (single nucleotide polymorphism)
    pub fn is_snp(&self) -> bool {
        self.ref_allele.len() == 1
            && self.alt_allele.len() == 1
            && self.ref_allele != "."
            && self.alt_allele != "."
            && self.ref_allele != "*"
            && self.alt_allele != "*"
    }

    /// Check if this is an MNP (multi-nucleotide polymorphism)
    /// MNPs have REF and ALT of equal length > 1, e.g., AT→GC
    pub fn is_mnp(&self) -> bool {
        let ref_len = self.ref_allele.len();
        let alt_len = self.alt_allele.len();
        ref_len == alt_len
            && ref_len > 1
            && !self.ref_allele.contains(',')
            && !self.alt_allele.contains(',')
            && self.ref_allele != "."
            && self.alt_allele != "."
    }

    /// Check if this is an indel
    pub fn is_indel(&self) -> bool {
        self.ref_allele.len() != self.alt_allele.len()
            || self.ref_allele.contains(',')
            || self.alt_allele.contains(',')
    }

    /// Decompose an MNP into individual SNPs
    /// e.g., AT→GC at pos 100 becomes: A→G at 100, T→C at 101
    /// Returns empty vec if not an MNP or if REF==ALT at all positions
    pub fn decompose_mnp(&self) -> Vec<VariantCall> {
        if !self.is_mnp() {
            return vec![];
        }

        let ref_bytes = self.ref_allele.as_bytes();
        let alt_bytes = self.alt_allele.as_bytes();
        let mut snps = Vec::new();

        for (i, (r, a)) in ref_bytes.iter().zip(alt_bytes.iter()).enumerate() {
            // Only create SNP if the bases differ
            if r != a {
                snps.push(VariantCall {
                    chrom: self.chrom.clone(),
                    pos: self.pos + i,  // Adjust position for each base
                    ref_allele: (*r as char).to_string(),
                    alt_allele: (*a as char).to_string(),
                    qual: self.qual,
                    filter: self.filter.clone(),
                    depth: self.depth,
                    allele_freq: self.allele_freq,
                    genotype: self.genotype.clone(),
                    info: format!("{};DECOMPOSED_FROM_MNP={}", self.info, self.ref_allele),
                    sample: self.sample.clone(),
                    ref_obs: self.ref_obs,
                    alt_obs: self.alt_obs,
                    alt_qual: self.alt_qual,
                    genotype_qual: self.genotype_qual,
                    avg_base_qual: self.avg_base_qual,
                    pvalue: self.pvalue,
                });
            }
        }

        snps
    }

    /// Check if this variant passed filters
    pub fn is_pass(&self) -> bool {
        self.filter == "PASS" || self.filter == "." || self.filter.is_empty()
    }

    /// Get a unique key for this position (chrom:pos)
    pub fn position_key(&self) -> String {
        format!("{}:{}", self.chrom, self.pos)
    }
}

/// Statistics about MNP decomposition
#[derive(Debug, Clone, Default)]
pub struct MnpStats {
    /// Number of MNPs found in the original VCF
    pub mnps_found: usize,
    /// Number of SNPs generated from MNP decomposition
    pub snps_from_mnps: usize,
    /// Original MNP strings for reporting (e.g., "AT→GC at 12345")
    pub mnp_details: Vec<String>,
}

/// Collection of variants from a single VCF file
#[derive(Debug, Clone)]
pub struct VcfFile {
    /// Path to the VCF file
    pub path: String,
    /// Sample name
    pub sample: String,
    /// Pipeline name (snippy, cfsan, etc.)
    pub pipeline: String,
    /// All variant calls (including decomposed MNPs)
    pub variants: Vec<VariantCall>,
    /// Variants indexed by position (chrom:pos -> variant)
    pub by_position: HashMap<String, VariantCall>,
    /// Statistics about MNP decomposition
    pub mnp_stats: MnpStats,
}

impl VcfFile {
    /// Get only SNPs
    pub fn snps(&self) -> Vec<&VariantCall> {
        self.variants.iter().filter(|v| v.is_snp()).collect()
    }

    /// Get only PASS variants
    pub fn passed(&self) -> Vec<&VariantCall> {
        self.variants.iter().filter(|v| v.is_pass()).collect()
    }

    /// Get only PASS SNPs
    pub fn passed_snps(&self) -> Vec<&VariantCall> {
        self.variants.iter().filter(|v| v.is_snp() && v.is_pass()).collect()
    }

    /// Check if a position has a variant
    pub fn has_position(&self, chrom: &str, pos: usize) -> bool {
        let key = format!("{}:{}", chrom, pos);
        self.by_position.contains_key(&key)
    }

    /// Get variant at position
    pub fn get_position(&self, chrom: &str, pos: usize) -> Option<&VariantCall> {
        let key = format!("{}:{}", chrom, pos);
        self.by_position.get(&key)
    }

    /// Check if any MNPs were decomposed
    pub fn has_mnp_decomposition(&self) -> bool {
        self.mnp_stats.mnps_found > 0
    }

    /// Get a warning message about MNP decomposition, if any
    pub fn mnp_warning(&self) -> Option<String> {
        if self.mnp_stats.mnps_found > 0 {
            Some(format!(
                "{} MNP(s) were decomposed into {} individual SNPs",
                self.mnp_stats.mnps_found,
                self.mnp_stats.snps_from_mnps
            ))
        } else {
            None
        }
    }
}

/// Parse a VCF file, automatically decomposing MNPs into individual SNPs
pub fn parse_vcf(path: &Path, pipeline: &str) -> Result<VcfFile> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("Failed to read VCF file: {}", path.display()))?;

    // Extract sample name from filename
    let sample = extract_sample_name(path);

    let mut variants = Vec::new();
    let mut by_position = HashMap::new();
    let mut mnp_stats = MnpStats::default();

    for line in content.lines() {
        // Skip header lines
        if line.starts_with('#') {
            continue;
        }

        // Parse variant line
        if let Some(variant) = parse_vcf_line(line, &sample) {
            // Check if this is an MNP that needs decomposition
            if variant.is_mnp() {
                let decomposed = variant.decompose_mnp();
                if !decomposed.is_empty() {
                    // Track MNP statistics
                    mnp_stats.mnps_found += 1;
                    mnp_stats.snps_from_mnps += decomposed.len();
                    mnp_stats.mnp_details.push(format!(
                        "{}→{} at {}:{}",
                        variant.ref_allele, variant.alt_allele, variant.chrom, variant.pos
                    ));

                    log::debug!(
                        "Decomposed MNP {}→{} at {}:{} into {} SNPs",
                        variant.ref_allele, variant.alt_allele, variant.chrom, variant.pos,
                        decomposed.len()
                    );

                    // Add decomposed SNPs
                    for snp in decomposed {
                        let key = snp.position_key();
                        by_position.insert(key, snp.clone());
                        variants.push(snp);
                    }
                }
            } else {
                // Regular variant (SNP, indel, etc.)
                let key = variant.position_key();
                by_position.insert(key, variant.clone());
                variants.push(variant);
            }
        }
    }

    // Log MNP summary if any were found
    if mnp_stats.mnps_found > 0 {
        log::info!(
            "VCF {}: Found {} MNPs, decomposed into {} SNPs",
            path.display(), mnp_stats.mnps_found, mnp_stats.snps_from_mnps
        );
    }

    Ok(VcfFile {
        path: path.display().to_string(),
        sample,
        pipeline: pipeline.to_string(),
        variants,
        by_position,
        mnp_stats,
    })
}

/// Parse a single VCF line
fn parse_vcf_line(line: &str, sample: &str) -> Option<VariantCall> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        return None;
    }

    let chrom = fields[0].to_string();
    let pos: usize = fields[1].parse().ok()?;
    let ref_allele = fields[3].to_string();
    let alt_allele = fields[4].to_string();
    let qual: f64 = fields[5].parse().unwrap_or(0.0);
    let filter = fields[6].to_string();
    let info = fields[7].to_string();

    // Parse INFO field for DP and AF
    let depth = parse_info_field(&info, "DP")
        .or_else(|| parse_info_field(&info, "ADP"));
    let allele_freq = parse_info_float(&info, "AF");

    // Parse additional quality metrics from INFO (Snippy format)
    let ref_obs = parse_info_field(&info, "RO");
    let alt_obs = parse_info_field(&info, "AO");
    let alt_qual = parse_info_field(&info, "QA");

    // Parse FORMAT/SAMPLE fields if present
    let format_data = if fields.len() >= 10 {
        parse_format_fields(fields[8], fields[9])
    } else {
        FormatData::default()
    };

    // Prefer FORMAT values over INFO values
    let depth = format_data.depth.or(depth);
    let allele_freq = format_data.allele_freq.or(allele_freq);

    Some(VariantCall {
        chrom,
        pos,
        ref_allele,
        alt_allele,
        qual,
        filter,
        depth,
        allele_freq,
        genotype: format_data.genotype,
        info,
        sample: sample.to_string(),
        // Additional quality metrics
        ref_obs: format_data.ref_obs.or(ref_obs),
        alt_obs: format_data.alt_obs.or(alt_obs),
        alt_qual: format_data.alt_qual.or(alt_qual),
        genotype_qual: format_data.genotype_qual,
        avg_base_qual: format_data.avg_base_qual,
        pvalue: format_data.pvalue,
    })
}

/// Parse an integer field from INFO
fn parse_info_field(info: &str, key: &str) -> Option<u32> {
    for part in info.split(';') {
        if let Some(value) = part.strip_prefix(&format!("{}=", key)) {
            return value.split(',').next()?.parse().ok();
        }
    }
    None
}

/// Parse a float field from INFO
fn parse_info_float(info: &str, key: &str) -> Option<f64> {
    for part in info.split(';') {
        if let Some(value) = part.strip_prefix(&format!("{}=", key)) {
            return value.split(',').next()?.parse().ok();
        }
    }
    None
}

/// Parsed FORMAT field data
#[derive(Debug, Default)]
struct FormatData {
    genotype: Option<String>,
    depth: Option<u32>,
    allele_freq: Option<f64>,
    ref_obs: Option<u32>,
    alt_obs: Option<u32>,
    alt_qual: Option<u32>,
    genotype_qual: Option<u32>,
    avg_base_qual: Option<u32>,
    pvalue: Option<f64>,
}

/// Parse FORMAT and SAMPLE fields
fn parse_format_fields(format: &str, sample_data: &str) -> FormatData {
    let format_keys: Vec<&str> = format.split(':').collect();
    let sample_values: Vec<&str> = sample_data.split(':').collect();

    let mut data = FormatData::default();

    for (i, key) in format_keys.iter().enumerate() {
        if i >= sample_values.len() {
            break;
        }
        let value = sample_values[i];

        match *key {
            "GT" => data.genotype = Some(value.to_string()),
            "DP" => data.depth = value.parse().ok(),
            "AD" => {
                // Parse allele depth: ref,alt -> calculate AF
                let parts: Vec<u32> = value.split(',')
                    .filter_map(|v| v.parse().ok())
                    .collect();
                if parts.len() >= 2 {
                    let total: u32 = parts.iter().sum();
                    if total > 0 {
                        data.allele_freq = Some(parts[1] as f64 / total as f64);
                    }
                    if data.depth.is_none() {
                        data.depth = Some(total);
                    }
                    // Also extract ref/alt observations from AD
                    data.ref_obs = Some(parts[0]);
                    data.alt_obs = Some(parts[1]);
                }
            }
            "AF" | "FREQ" => {
                // FREQ is CFSAN format like "100%"
                let clean_value = value.trim_end_matches('%');
                if let Ok(v) = clean_value.parse::<f64>() {
                    data.allele_freq = Some(if value.contains('%') { v / 100.0 } else { v });
                }
            }
            // Snippy fields
            "RO" => data.ref_obs = value.parse().ok(),
            "AO" => data.alt_obs = value.parse().ok(),
            "QA" => data.alt_qual = value.parse().ok(),
            // CFSAN fields
            "GQ" => data.genotype_qual = value.parse().ok(),
            "ABQ" => data.avg_base_qual = value.parse().ok(),
            "PVAL" => data.pvalue = value.parse().ok(),
            "RD" => {
                // CFSAN reference depth
                if data.ref_obs.is_none() {
                    data.ref_obs = value.parse().ok();
                }
            }
            _ => {}
        }
    }

    data
}

/// Extract sample name from VCF file path
fn extract_sample_name(path: &Path) -> String {
    let filename = path.file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("unknown");

    // Try to extract sample name from common patterns
    // snippy: sample_snippy.vcf or sample.vcf
    // cfsan: var.flt.vcf in sample directory

    if filename == "var.flt.vcf" || filename == "consensus.vcf" {
        // CFSAN pattern: get parent directory name
        path.parent()
            .and_then(|p| p.file_name())
            .and_then(|n| n.to_str())
            .unwrap_or("unknown")
            .to_string()
    } else {
        // Snippy pattern: remove common suffixes
        filename
            .trim_end_matches(".vcf")
            .trim_end_matches(".vcf.gz")
            .trim_end_matches("_snippy")
            .trim_end_matches(".flt")
            .trim_end_matches(".raw")
            .to_string()
    }
}

/// Load multiple VCF files from a directory
pub fn load_vcfs_from_dir(dir: &Path, pipeline: &str) -> Result<Vec<VcfFile>> {
    let mut vcfs = Vec::new();

    if !dir.exists() {
        anyhow::bail!("VCF directory does not exist: {}", dir.display());
    }

    // For CFSAN, look for var.flt.vcf in sample subdirectories
    if pipeline == "cfsan" {
        let samples_dir = dir.join("samples");
        if samples_dir.exists() {
            for entry in std::fs::read_dir(&samples_dir)? {
                let entry = entry?;
                let vcf_path = entry.path().join("var.flt.vcf");
                if vcf_path.exists() {
                    match parse_vcf(&vcf_path, pipeline) {
                        Ok(vcf) => vcfs.push(vcf),
                        Err(e) => log::warn!("Failed to parse {}: {}", vcf_path.display(), e),
                    }
                }
            }
        }
    }

    // Look for VCF files directly in directory
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().map(|e| e == "vcf").unwrap_or(false) {
            match parse_vcf(&path, pipeline) {
                Ok(vcf) => vcfs.push(vcf),
                Err(e) => log::warn!("Failed to parse {}: {}", path.display(), e),
            }
        }
    }

    // Also check subdirectories for *_snippy/*.vcf pattern
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() {
            let dirname = path.file_name()
                .and_then(|n| n.to_str())
                .unwrap_or("");
            if dirname.ends_with("_snippy") {
                // Look for multiple possible VCF patterns:
                // 1. snps.vcf (standard snippy output)
                // 2. <dirname>.vcf (e.g., TE15676_snippy.vcf)
                // 3. Any *.vcf that's not .raw.vcf or .subs.vcf
                let standard_names = ["snps.vcf", "snps.filt.vcf"];
                let mut found = false;

                for vcf_name in &standard_names {
                    let vcf_path = path.join(vcf_name);
                    if vcf_path.exists() {
                        match parse_vcf(&vcf_path, pipeline) {
                            Ok(vcf) => vcfs.push(vcf),
                            Err(e) => log::warn!("Failed to parse {}: {}", vcf_path.display(), e),
                        }
                        found = true;
                        break;
                    }
                }

                // Try <dirname>.vcf pattern
                if !found {
                    let vcf_path = path.join(format!("{}.vcf", dirname));
                    if vcf_path.exists() {
                        match parse_vcf(&vcf_path, pipeline) {
                            Ok(vcf) => vcfs.push(vcf),
                            Err(e) => log::warn!("Failed to parse {}: {}", vcf_path.display(), e),
                        }
                        found = true;
                    }
                }

                // Fallback: look for any .vcf file (excluding raw, subs)
                if !found {
                    if let Ok(entries) = std::fs::read_dir(&path) {
                        for subentry in entries.flatten() {
                            let subpath = subentry.path();
                            let filename = subpath.file_name()
                                .and_then(|n| n.to_str())
                                .unwrap_or("");
                            if filename.ends_with(".vcf")
                                && !filename.ends_with(".raw.vcf")
                                && !filename.ends_with(".subs.vcf")
                                && !filename.ends_with(".filt.vcf") {
                                match parse_vcf(&subpath, pipeline) {
                                    Ok(vcf) => vcfs.push(vcf),
                                    Err(e) => log::warn!("Failed to parse {}: {}", subpath.display(), e),
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    Ok(vcfs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_vcf_line() {
        let line = "AL591824.1\t32859\t.\tG\tT\t2577.1\tPASS\tDP=87;AF=1.0\tGT:DP\t1/1:87";
        let variant = parse_vcf_line(line, "test_sample").unwrap();

        assert_eq!(variant.chrom, "AL591824.1");
        assert_eq!(variant.pos, 32859);
        assert_eq!(variant.ref_allele, "G");
        assert_eq!(variant.alt_allele, "T");
        assert!((variant.qual - 2577.1).abs() < 0.1);
        assert_eq!(variant.filter, "PASS");
        assert_eq!(variant.depth, Some(87));
        assert!(variant.is_snp());
        assert!(variant.is_pass());
    }

    #[test]
    fn test_is_snp() {
        let snp = VariantCall {
            chrom: "chr1".to_string(),
            pos: 100,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            qual: 100.0,
            filter: "PASS".to_string(),
            depth: Some(50),
            allele_freq: Some(1.0),
            genotype: Some("1/1".to_string()),
            info: String::new(),
            sample: "test".to_string(),
            ref_obs: None,
            alt_obs: None,
            alt_qual: None,
            genotype_qual: None,
            avg_base_qual: None,
            pvalue: None,
        };
        assert!(snp.is_snp());

        let indel = VariantCall {
            ref_allele: "AT".to_string(),
            alt_allele: "A".to_string(),
            ..snp.clone()
        };
        assert!(!indel.is_snp());
        assert!(indel.is_indel());
    }

    #[test]
    fn test_is_mnp() {
        let mnp = VariantCall {
            chrom: "chr1".to_string(),
            pos: 100,
            ref_allele: "AT".to_string(),
            alt_allele: "GC".to_string(),
            qual: 100.0,
            filter: "PASS".to_string(),
            depth: Some(50),
            allele_freq: Some(1.0),
            genotype: Some("1/1".to_string()),
            info: String::new(),
            sample: "test".to_string(),
            ref_obs: None,
            alt_obs: None,
            alt_qual: None,
            genotype_qual: None,
            avg_base_qual: None,
            pvalue: None,
        };
        assert!(mnp.is_mnp());
        assert!(!mnp.is_snp());
        assert!(!mnp.is_indel());
    }

    #[test]
    fn test_decompose_mnp() {
        let mnp = VariantCall {
            chrom: "chr1".to_string(),
            pos: 100,
            ref_allele: "ATG".to_string(),
            alt_allele: "GCC".to_string(),
            qual: 150.0,
            filter: "PASS".to_string(),
            depth: Some(60),
            allele_freq: Some(1.0),
            genotype: Some("1/1".to_string()),
            info: "DP=60".to_string(),
            sample: "test".to_string(),
            ref_obs: None,
            alt_obs: None,
            alt_qual: None,
            genotype_qual: None,
            avg_base_qual: None,
            pvalue: None,
        };

        let decomposed = mnp.decompose_mnp();

        // ATG→GCC: A→G at 100, T→C at 101, G→C at 102
        assert_eq!(decomposed.len(), 3);

        assert_eq!(decomposed[0].pos, 100);
        assert_eq!(decomposed[0].ref_allele, "A");
        assert_eq!(decomposed[0].alt_allele, "G");
        assert!(decomposed[0].is_snp());

        assert_eq!(decomposed[1].pos, 101);
        assert_eq!(decomposed[1].ref_allele, "T");
        assert_eq!(decomposed[1].alt_allele, "C");

        assert_eq!(decomposed[2].pos, 102);
        assert_eq!(decomposed[2].ref_allele, "G");
        assert_eq!(decomposed[2].alt_allele, "C");

        // Check that quality and depth are preserved
        assert!((decomposed[0].qual - 150.0).abs() < 0.1);
        assert_eq!(decomposed[0].depth, Some(60));
    }

    #[test]
    fn test_decompose_mnp_partial_change() {
        // MNP where only some bases differ: ATG→ATC (only G→C)
        let mnp = VariantCall {
            chrom: "chr1".to_string(),
            pos: 100,
            ref_allele: "ATG".to_string(),
            alt_allele: "ATC".to_string(),
            qual: 100.0,
            filter: "PASS".to_string(),
            depth: Some(50),
            allele_freq: None,
            genotype: None,
            info: String::new(),
            sample: "test".to_string(),
            ref_obs: None,
            alt_obs: None,
            alt_qual: None,
            genotype_qual: None,
            avg_base_qual: None,
            pvalue: None,
        };

        let decomposed = mnp.decompose_mnp();

        // Only G→C at position 102 differs
        assert_eq!(decomposed.len(), 1);
        assert_eq!(decomposed[0].pos, 102);
        assert_eq!(decomposed[0].ref_allele, "G");
        assert_eq!(decomposed[0].alt_allele, "C");
    }

    #[test]
    fn test_non_mnp_decompose_returns_empty() {
        let snp = VariantCall {
            chrom: "chr1".to_string(),
            pos: 100,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            qual: 100.0,
            filter: "PASS".to_string(),
            depth: Some(50),
            allele_freq: None,
            genotype: None,
            info: String::new(),
            sample: "test".to_string(),
            ref_obs: None,
            alt_obs: None,
            alt_qual: None,
            genotype_qual: None,
            avg_base_qual: None,
            pvalue: None,
        };

        // SNPs should return empty decomposition
        assert!(snp.decompose_mnp().is_empty());
    }
}
