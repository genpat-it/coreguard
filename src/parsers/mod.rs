//! Core SNP parser plugin system.
//!
//! Each pipeline format (Snippy, CFSAN, etc.) has its own parser module.
//! To add a new format:
//! 1. Create `src/parsers/my_format.rs` implementing [`CoreSnpParser`]
//! 2. Register it in [`all_parsers()`] below
//!
//! The dispatcher [`parse_core_snps()`] tries each parser in order and uses
//! the first one that accepts the file.

mod snippy;
mod cfsan;

use std::collections::HashMap;
use std::path::Path;

// Re-export parser implementations for direct use if needed
pub use snippy::SnippyCoreTabParser;
pub use cfsan::CfsanSnplistParser;

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

/// Return all registered parsers, in priority order.
///
/// Parsers are tried in order; the first one whose [`CoreSnpParser::can_parse()`]
/// returns `true` wins. Put more specific parsers first.
fn all_parsers() -> Vec<Box<dyn CoreSnpParser>> {
    vec![
        Box::new(SnippyCoreTabParser),
        Box::new(CfsanSnplistParser),
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
