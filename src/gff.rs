//! GFF/GTF annotation parsing and gene zone analysis
//!
//! Identifies which annotated features (genes, CDS) are affected by alignment gaps.

use crate::gaps::Region;
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// A genomic feature from GFF/GTF
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Feature {
    /// Chromosome/contig name
    pub seqid: String,
    /// Feature type (gene, CDS, exon, etc.)
    pub feature_type: String,
    /// Start position (1-based)
    pub start: usize,
    /// End position (1-based, inclusive)
    pub end: usize,
    /// Strand (+, -, or .)
    pub strand: char,
    /// Gene name (if available)
    pub gene_name: Option<String>,
    /// Locus tag
    pub locus_tag: Option<String>,
    /// Product description
    pub product: Option<String>,
    /// Feature ID
    pub id: Option<String>,
}

impl Feature {
    /// Feature length in bp
    pub fn length(&self) -> usize {
        self.end.saturating_sub(self.start) + 1
    }

    /// Check if this feature overlaps with a region
    pub fn overlaps(&self, region: &Region) -> bool {
        // Convert to 0-based for comparison with Region
        let feat_start = self.start.saturating_sub(1);
        let feat_end = self.end;
        feat_start < region.end && region.start < feat_end
    }

    /// Calculate overlap amount with a region (in bp)
    pub fn overlap_bp(&self, region: &Region) -> usize {
        let feat_start = self.start.saturating_sub(1);
        let feat_end = self.end;

        if feat_start >= region.end || region.start >= feat_end {
            return 0;
        }

        let overlap_start = feat_start.max(region.start);
        let overlap_end = feat_end.min(region.end);
        overlap_end.saturating_sub(overlap_start)
    }

    /// Display name for the feature
    pub fn display_name(&self) -> String {
        self.gene_name
            .clone()
            .or_else(|| self.locus_tag.clone())
            .or_else(|| self.id.clone())
            .unwrap_or_else(|| format!("{}:{}-{}", self.seqid, self.start, self.end))
    }
}

/// Parsed GFF annotation
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Annotation {
    /// All features
    pub features: Vec<Feature>,
    /// Features indexed by type
    pub by_type: HashMap<String, Vec<usize>>,
    /// Source file path
    pub source_path: String,
    /// Reference sequence ID (from GFF)
    pub sequence_id: Option<String>,
    /// Total genes
    pub total_genes: usize,
    /// Total CDS
    pub total_cds: usize,
}

impl Annotation {
    /// Parse a GFF/GTF file
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)
            .with_context(|| format!("Failed to open GFF file: {}", path.display()))?;
        let reader = BufReader::new(file);

        let mut features = Vec::new();
        let mut by_type: HashMap<String, Vec<usize>> = HashMap::new();
        let mut sequence_id = None;

        for line in reader.lines() {
            let line = line?;

            // Skip comments and empty lines
            if line.starts_with('#') || line.trim().is_empty() {
                // Extract sequence-region if present
                if line.starts_with("##sequence-region") {
                    let parts: Vec<&str> = line.split_whitespace().collect();
                    if parts.len() >= 2 {
                        sequence_id = Some(parts[1].to_string());
                    }
                }
                continue;
            }

            if let Some(feature) = parse_gff_line(&line) {
                let idx = features.len();
                by_type
                    .entry(feature.feature_type.clone())
                    .or_default()
                    .push(idx);
                features.push(feature);
            }
        }

        let total_genes = by_type.get("gene").map(|v| v.len()).unwrap_or(0);
        let total_cds = by_type.get("CDS").map(|v| v.len()).unwrap_or(0);

        Ok(Annotation {
            features,
            by_type,
            source_path: path.display().to_string(),
            sequence_id,
            total_genes,
            total_cds,
        })
    }

    /// Get all genes
    pub fn genes(&self) -> Vec<&Feature> {
        self.by_type
            .get("gene")
            .map(|indices| indices.iter().map(|&i| &self.features[i]).collect())
            .unwrap_or_default()
    }

    /// Get all CDS
    pub fn cds(&self) -> Vec<&Feature> {
        self.by_type
            .get("CDS")
            .map(|indices| indices.iter().map(|&i| &self.features[i]).collect())
            .unwrap_or_default()
    }

    /// Find features overlapping with gap regions
    pub fn features_in_gaps(&self, gaps: &[Region], feature_types: &[&str]) -> Vec<AffectedFeature> {
        let mut affected = Vec::new();

        for gap in gaps {
            for feat_type in feature_types {
                if let Some(indices) = self.by_type.get(*feat_type) {
                    for &idx in indices {
                        let feature = &self.features[idx];
                        let overlap = feature.overlap_bp(gap);
                        if overlap > 0 {
                            let pct = (overlap as f64 / feature.length() as f64) * 100.0;
                            affected.push(AffectedFeature {
                                feature: feature.clone(),
                                gap_region: gap.clone(),
                                overlap_bp: overlap,
                                overlap_pct: pct,
                            });
                        }
                    }
                }
            }
        }

        // Sort by overlap percentage (most affected first)
        affected.sort_by(|a, b| b.overlap_pct.partial_cmp(&a.overlap_pct).unwrap_or(std::cmp::Ordering::Equal));

        // Deduplicate by feature name (keep highest overlap)
        let mut seen = std::collections::HashSet::new();
        affected.retain(|af| seen.insert(af.feature.display_name()));

        affected
    }
}

/// A feature affected by a gap
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AffectedFeature {
    /// The affected feature
    pub feature: Feature,
    /// The gap region causing the effect
    pub gap_region: Region,
    /// Overlap in bp
    pub overlap_bp: usize,
    /// Overlap as percentage of feature length
    pub overlap_pct: f64,
}

impl AffectedFeature {
    /// Is this feature completely removed (100% overlap)?
    pub fn is_completely_removed(&self) -> bool {
        self.overlap_pct >= 99.9
    }

    /// Is this feature partially affected (>50% overlap)?
    pub fn is_significantly_affected(&self) -> bool {
        self.overlap_pct >= 50.0
    }
}

/// Summary of gene zone analysis for a sample
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneZoneAnalysis {
    /// Sample name
    pub sample: String,
    /// Total genes in annotation
    pub total_genes: usize,
    /// Genes affected by gaps
    pub affected_genes: usize,
    /// Genes completely removed (100% in gap)
    pub removed_genes: usize,
    /// Genes significantly affected (>50% in gap)
    pub significantly_affected: usize,
    /// List of affected features
    pub affected_features: Vec<AffectedFeature>,
    /// Total CDS in annotation
    pub total_cds: usize,
    /// CDS affected by gaps
    pub affected_cds: usize,
}

/// Parse a single GFF line
fn parse_gff_line(line: &str) -> Option<Feature> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 9 {
        return None;
    }

    let seqid = fields[0].to_string();
    let feature_type = fields[2].to_string();
    let start: usize = fields[3].parse().ok()?;
    let end: usize = fields[4].parse().ok()?;
    let strand = fields[6].chars().next().unwrap_or('.');

    // Parse attributes
    let attributes = fields[8];
    let gene_name = extract_attribute(attributes, "gene")
        .or_else(|| extract_attribute(attributes, "Name"));
    let locus_tag = extract_attribute(attributes, "locus_tag");
    let product = extract_attribute(attributes, "product");
    let id = extract_attribute(attributes, "ID");

    Some(Feature {
        seqid,
        feature_type,
        start,
        end,
        strand,
        gene_name,
        locus_tag,
        product,
        id,
    })
}

/// Extract an attribute value from GFF attributes field
fn extract_attribute(attributes: &str, key: &str) -> Option<String> {
    for attr in attributes.split(';') {
        let attr = attr.trim();
        // GFF3 format: key=value
        if let Some(value) = attr.strip_prefix(&format!("{}=", key)) {
            // URL decode common sequences
            let value = value
                .replace("%3B", ";")
                .replace("%2C", ",")
                .replace("%3D", "=")
                .replace("%26", "&")
                .replace("%25", "%");
            return Some(value);
        }
        // GTF format: key "value"
        if attr.starts_with(key) {
            let parts: Vec<&str> = attr.splitn(2, ' ').collect();
            if parts.len() == 2 {
                let value = parts[1].trim_matches('"');
                return Some(value.to_string());
            }
        }
    }
    None
}

/// Analyze gene zones for all samples
pub fn analyze_gene_zones(
    annotation: &Annotation,
    sample_gaps: &[(String, Vec<Region>)],
) -> Vec<GeneZoneAnalysis> {
    sample_gaps
        .iter()
        .map(|(sample, gaps)| {
            let affected = annotation.features_in_gaps(gaps, &["gene", "CDS"]);

            // Count stats before moving affected
            let affected_genes = affected.iter().filter(|af| af.feature.feature_type == "gene").count();
            let affected_cds = affected.iter().filter(|af| af.feature.feature_type == "CDS").count();
            let removed_genes = affected.iter()
                .filter(|af| af.feature.feature_type == "gene" && af.is_completely_removed())
                .count();
            let significantly_affected = affected.iter()
                .filter(|af| af.feature.feature_type == "gene" && af.is_significantly_affected())
                .count();

            GeneZoneAnalysis {
                sample: sample.clone(),
                total_genes: annotation.total_genes,
                affected_genes,
                removed_genes,
                significantly_affected,
                affected_features: affected,
                total_cds: annotation.total_cds,
                affected_cds,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_feature_overlap() {
        let feature = Feature {
            seqid: "chr1".to_string(),
            feature_type: "gene".to_string(),
            start: 100,  // 1-based
            end: 200,
            strand: '+',
            gene_name: Some("testgene".to_string()),
            locus_tag: None,
            product: None,
            id: None,
        };

        // Region is 0-based
        let region = Region::new(99, 150);  // 0-based: 99-150, overlaps with feature 100-200
        assert!(feature.overlaps(&region));
        assert_eq!(feature.overlap_bp(&region), 51);  // positions 99-150 overlap with 99-200

        let region_no_overlap = Region::new(0, 50);
        assert!(!feature.overlaps(&region_no_overlap));
        assert_eq!(feature.overlap_bp(&region_no_overlap), 0);
    }

    #[test]
    fn test_parse_attribute() {
        let attrs = "ID=gene-lmo0001;Name=dnaA;locus_tag=lmo0001;product=DnaA protein";
        assert_eq!(extract_attribute(attrs, "Name"), Some("dnaA".to_string()));
        assert_eq!(extract_attribute(attrs, "locus_tag"), Some("lmo0001".to_string()));
        assert_eq!(extract_attribute(attrs, "product"), Some("DnaA protein".to_string()));
    }
}
