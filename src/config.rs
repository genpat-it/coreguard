//! Configuration parsing for CoreGuard
//!
//! Parses YAML configuration files for pipeline-agnostic SNP comparison.

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

/// Main configuration structure
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    /// Reference genome
    pub reference: ReferenceConfig,

    /// Sample definitions
    #[serde(default)]
    pub samples: HashMap<String, SampleConfig>,

    /// Pipeline definitions
    #[serde(default)]
    pub pipelines: HashMap<String, PipelineConfig>,

    /// Optional settings
    #[serde(default)]
    pub options: Options,
}

/// Reference genome configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceConfig {
    /// Path to reference FASTA
    pub path: String,

    /// Display label (optional)
    #[serde(default)]
    pub label: Option<String>,
}

/// Sample configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleConfig {
    /// Display label (optional)
    #[serde(default)]
    pub label: Option<String>,

    /// Paired-end or single-end FASTQ reads (optional)
    #[serde(default)]
    pub reads: Vec<String>,
}

/// Pipeline configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineConfig {
    /// Display label (optional)
    #[serde(default)]
    pub label: Option<String>,

    /// Command line used to run this pipeline (optional, for documentation)
    #[serde(default)]
    pub command: Option<String>,

    /// Mark this pipeline as ground truth (baseline for comparison)
    #[serde(default)]
    pub ground_truth: bool,

    /// Sample files for this pipeline
    #[serde(default)]
    pub samples: HashMap<String, PipelineSampleFiles>,
}

/// Files for a sample in a pipeline
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineSampleFiles {
    /// VCF file (optional) - for SNPs
    #[serde(default)]
    pub vcf: Option<String>,

    /// BAM file (optional) - for gaps/coverage
    #[serde(default)]
    pub bam: Option<String>,
}

/// Optional settings
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Options {
    /// Minimum depth to consider position covered
    #[serde(default = "default_min_depth")]
    pub min_depth: usize,

    /// Minimum SNP quality to include
    #[serde(default = "default_min_qual")]
    pub min_qual: f64,

    /// Include insertions/deletions
    #[serde(default)]
    pub include_indels: bool,
}

fn default_min_depth() -> usize {
    1
}

fn default_min_qual() -> f64 {
    20.0
}

impl Config {
    /// Load configuration from YAML file
    pub fn from_yaml<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let content = std::fs::read_to_string(path)
            .with_context(|| format!("Failed to read config file: {}", path.display()))?;

        let config: Config = serde_yaml::from_str(&content)
            .with_context(|| format!("Failed to parse YAML config: {}", path.display()))?;

        config.validate()?;
        Ok(config)
    }

    /// Validate configuration
    pub fn validate(&self) -> Result<()> {
        // Check reference exists
        if !Path::new(&self.reference.path).exists() {
            anyhow::bail!("Reference file not found: {}", self.reference.path);
        }

        // Check at least one pipeline defined
        if self.pipelines.is_empty() {
            anyhow::bail!("At least one pipeline must be defined");
        }

        // Check pipeline files exist
        for (pipeline_id, pipeline) in &self.pipelines {
            for (sample_id, files) in &pipeline.samples {
                if let Some(vcf) = &files.vcf {
                    if !Path::new(vcf).exists() {
                        anyhow::bail!(
                            "VCF file not found for pipeline '{}', sample '{}': {}",
                            pipeline_id,
                            sample_id,
                            vcf
                        );
                    }
                }
                if let Some(bam) = &files.bam {
                    if !Path::new(bam).exists() {
                        anyhow::bail!(
                            "BAM file not found for pipeline '{}', sample '{}': {}",
                            pipeline_id,
                            sample_id,
                            bam
                        );
                    }
                }
                // At least one of vcf or bam must be provided
                if files.vcf.is_none() && files.bam.is_none() {
                    anyhow::bail!(
                        "No VCF or BAM file specified for pipeline '{}', sample '{}'",
                        pipeline_id,
                        sample_id
                    );
                }
            }
        }

        Ok(())
    }

    /// Get all unique sample IDs across all pipelines
    pub fn all_sample_ids(&self) -> Vec<String> {
        let mut ids: Vec<String> = self
            .pipelines
            .values()
            .flat_map(|p| p.samples.keys().cloned())
            .collect();
        ids.sort();
        ids.dedup();
        ids
    }

    /// Get all pipeline IDs
    pub fn all_pipeline_ids(&self) -> Vec<String> {
        let mut ids: Vec<String> = self.pipelines.keys().cloned().collect();
        ids.sort();
        ids
    }

    /// Get display label for a sample
    #[allow(dead_code)]
    pub fn sample_label(&self, sample_id: &str) -> String {
        self.samples
            .get(sample_id)
            .and_then(|s| s.label.clone())
            .unwrap_or_else(|| sample_id.to_string())
    }

    /// Get display label for a pipeline
    #[allow(dead_code)]
    pub fn pipeline_label(&self, pipeline_id: &str) -> String {
        self.pipelines
            .get(pipeline_id)
            .and_then(|p| p.label.clone())
            .unwrap_or_else(|| pipeline_id.to_string())
    }

    /// Get the ground truth pipeline ID (if any)
    pub fn ground_truth_pipeline(&self) -> Option<String> {
        self.pipelines
            .iter()
            .find(|(_, p)| p.ground_truth)
            .map(|(id, _)| id.clone())
    }

    /// Check if a pipeline is the ground truth
    pub fn is_ground_truth(&self, pipeline_id: &str) -> bool {
        self.pipelines
            .get(pipeline_id)
            .map(|p| p.ground_truth)
            .unwrap_or(false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_yaml() {
        let yaml = r#"
reference:
  path: test.fa
  label: "Test Reference"

samples:
  sample1:
    label: "Sample One"

pipelines:
  snippy:
    label: "Snippy v4.6"
    samples:
      sample1:
        vcf: snippy/s1.vcf
        bam: snippy/s1.bam
"#;
        let config: Config = serde_yaml::from_str(yaml).unwrap();
        assert_eq!(config.reference.path, "test.fa");
        assert_eq!(config.reference.label, Some("Test Reference".to_string()));
        assert!(config.pipelines.contains_key("snippy"));
    }
}
