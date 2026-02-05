//! Core SNP parser for CoreGuard TSV format.
//!
//! CoreGuard uses a standard TSV format identical to Snippy's `core.tab`:
//! `CHR\tPOS\tREF\tsample1\tsample2\t...`
//!
//! For pipelines that don't natively output this format (CFSAN, GATK, etc.),
//! use `coreguard convert <format>` to generate a compatible file.

mod coreguard_tsv;

use std::collections::HashMap;
use std::path::Path;

pub use coreguard_tsv::CoreGuardTsvParser;

/// Core SNP data extracted from a pipeline's native output file.
#[derive(Debug, Clone)]
pub struct CoreSnpData {
    /// Pipeline positions with per-sample alleles
    pub positions: Vec<CoreSnpPosition>,
    /// True if the format includes per-sample alleles (e.g., Snippy core.tab)
    pub has_alleles: bool,
    /// Number of discriminating positions (at least 2 samples differ).
    /// Exact for formats with alleles, conservative for position-only formats.
    pub discriminating_count: usize,
}

/// A single position from a core SNP file.
#[derive(Debug, Clone)]
pub struct CoreSnpPosition {
    /// 0-based genomic position
    pub pos: usize,
    /// Reference allele (if available)
    pub ref_allele: Option<String>,
    /// Sample -> allele mapping (empty if format lacks allele info)
    pub alleles: HashMap<String, String>,
    /// Samples that carry a SNP at this position
    pub samples_with_snp: Vec<String>,
}

/// Trait for parsing pipeline-specific core SNP output files.
///
/// Each implementation handles one file format. The parser is selected
/// automatically via [`parse_core_snps()`] based on [`can_parse()`].
///
/// # Example
/// ```ignore
/// pub struct MyFormatParser;
///
/// impl CoreSnpParser for MyFormatParser {
///     fn format_name(&self) -> &str { "my_format.tsv" }
///     fn can_parse(&self, path: &Path) -> bool { /* peek at header */ }
///     fn parse(&self, path: &Path) -> anyhow::Result<CoreSnpData> { /* ... */ }
/// }
/// ```
pub trait CoreSnpParser {
    /// Human-readable name of the format (e.g., "snippycore.tab")
    fn format_name(&self) -> &str;

    /// Check if this parser can handle the given file (peek at header/content).
    fn can_parse(&self, path: &Path) -> bool;

    /// Parse the file and return core SNP data.
    fn parse(&self, path: &Path) -> anyhow::Result<CoreSnpData>;
}

/// Return the CoreGuard TSV parser.
fn all_parsers() -> Vec<Box<dyn CoreSnpParser>> {
    vec![
        Box::new(CoreGuardTsvParser),
    ]
}

/// Try all registered parsers and return the first successful result.
pub fn parse_core_snps(path: &str) -> anyhow::Result<CoreSnpData> {
    let path = Path::new(path);
    for parser in all_parsers() {
        if parser.can_parse(path) {
            log::info!("Detected core SNP format: {}", parser.format_name());
            return parser.parse(path);
        }
    }
    anyhow::bail!("No parser found for core SNP file: {}", path.display())
}
