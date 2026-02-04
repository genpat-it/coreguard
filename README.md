# CoreGuard

Pipeline-agnostic SNP comparison tool for bacterial genomics.

## Overview

CoreGuard compares SNP calls from multiple pipelines (Snippy, CFSAN, GATK, etc.), helping identify discrepancies and artifacts in variant calling. Optionally, a reference alignment can be used to detect coverage gaps and validate SNP calls.

**Live Viewer:** [https://genpat-it.github.io/coreguard/](https://genpat-it.github.io/coreguard/)

### Demo

Click **"Load Demo"** in the viewer to explore a pre-loaded dataset:

- **Organism**: *Listeria monocytogenes* (4 samples)
- **Pipelines compared**: Snippy, CFSAN SNP Pipeline
- **Reference alignment**: minimap2 (BAM pileup without variant calling)

The demo illustrates how different pipelines can produce varying SNP calls on the same data, and how CoreGuard helps identify these discrepancies.

## Why CoreGuard?

Different SNP pipelines can produce different results on the same data. This matters because:

- **Outbreak investigation**: SNP distances determine if isolates belong to the same outbreak. A difference of 1-2 SNPs can change epidemiological conclusions.
- **Pipeline validation**: When adopting a new pipeline, you need to verify it produces consistent results with established methods.
- **Troubleshooting**: When results seem wrong, CoreGuard helps visualize exactly where and why pipelines disagree.
- **Quality control**: Identify problematic samples or genomic regions where pipelines consistently fail.

## Features

### Core Features
- **Multi-pipeline comparison**: Compare SNP calls from any number of pipelines side-by-side
- **Reference alignment support** (optional): Designate a baseline BAM alignment for gap detection and validation
- **Dashboard viewer**: KPI dashboard with pre-computed statistics, runs entirely in the browser via WebAssembly
- **Pre-computed statistics**: All KPIs computed at CLI time — the viewer loads instantly with no recomputation

### SNP Analysis
- **MNP decomposition**: Multi-nucleotide polymorphisms (e.g., `TTGGCG→CCGGCT`) are automatically decomposed into individual SNPs
- **Pre-computed distance matrices**: Display pairwise SNP distance matrices from each pipeline's native output
- **Reference discriminating SNPs vs pipelines**: Compare reference alignment discriminating positions against each pipeline's core SNP output, with per-position drill-down detail

### Per-Pipeline KPI Dashboard
For each pipeline, CoreGuard computes:

- **Gap-Intersect / Gap-Union**: Two gap-exclusion strategies
  - *Gap-Intersect*: exclude positions where ALL samples have a gap (permissive)
  - *Gap-Union*: exclude positions where ANY sample has a gap (restrictive)
- **Usable Space**: Reference length minus excluded gap positions
- **Total SNPs**: Positions in usable space where at least one sample has an alt allele
- **Consensus SNPs**: All samples agree on the same alt allele (non-discriminating)
- **Discriminating SNPs**: Samples differ — these contribute to Hamming distance
- **Missing Calls**: Positions where some samples have a call but others don't (variant calling inconsistency)

#### Discriminating SNP Breakdown
For non-reference pipelines, discriminating SNPs are classified by heuristic:
- **Gap-affected**: At least one sample was skipped due to gap — partial comparison
- **Ref-consensus**: Reference alignment shows all samples agree — likely variant calling artifact
- **Majority-rule**: All but one sample agree — likely single-sample VC miss
- **Confirmed**: Genuine disagreement between samples

#### Pairwise Statistics
- **Avg Usable Space / Discriminating SNPs**: Averaged across all sample pairs (Gap-Union)
- **Min / Median / Max**: Distribution of pairwise discriminating SNPs
- **Per-Sample Table**: For each sample, average usable space and discriminating SNPs across all pairs involving that sample

#### BAM Pileup SNP Calling (Reference Alignment)
The reference pipeline calls SNPs from BAM pileup using majority vote:
1. Count reads at each position (A, C, G, T)
2. If total depth < `min_depth` → **Gap**
3. If majority base < `min_consensus` (default 80%) → **Ambiguous** (skipped)
4. If majority base ≠ reference → **SNP**

### Viewer Panels
- **Description**: Project description in Markdown (from config)
- **Pipelines**: Pipeline metadata (labels, commands, data types)
- **Statistics**: KPI dashboard with reference metrics and reference discriminating SNPs vs pipelines
- **SNP Distance Matrix**: Pre-computed distance matrices from each pipeline

### User Interface
- **Dark/Light theme**: Toggle between themes
- **Collapsible panels**: Organize information in expandable sections
- **Info icons**: Hover/click for detailed explanations of each metric

## Installation

```bash
# Clone the repository
git clone https://github.com/genpat-it/coreguard.git
cd coreguard

# Build with Rust
cargo build --release

# Binary will be at ./target/release/coreguard
```

## Usage

### 1. Create a configuration file

```yaml
# project.yaml
reference:
  path: reference.fasta
  label: "My Reference Genome"

# Optional: project description (inline markdown or path to .md file)
description: "## My Study\nComparison of Snippy vs CFSAN on *Listeria* dataset."
# Or: description: docs/study_description.md

samples:
  sample1:
    label: "Sample 1"   # optional display label
  sample2: {}
  sample3: {}

pipelines:
  # Reference alignment (BAM only — gap detection and SNP validation baseline)
  minimap2:
    label: "Reference (minimap2)"
    command: "minimap2 -ax sr -t 8 ref.fa reads_1.fq.gz reads_2.fq.gz | samtools sort -o out.bam"
    ground_truth: true
    samples:
      sample1:
        bam: alignments/sample1.bam
      sample2:
        bam: alignments/sample2.bam
      sample3:
        bam: alignments/sample3.bam

  # SNP pipeline with core_snps + BAM (recommended)
  snippy:
    label: "Snippy v4.6"
    command: "snippy --ref reference.fa --R1 reads_1.fq.gz --R2 reads_2.fq.gz --outdir out"
    distance_matrix: snippy/core.distances.tsv     # pre-computed distance matrix (optional)
    core_snps: snippy/core.tab                      # core SNP output file (recommended)
    samples:
      sample1:
        bam: snippy/sample1/snps.bam
      sample2:
        bam: snippy/sample2/snps.bam
      sample3:
        bam: snippy/sample3/snps.bam

  # Another pipeline
  cfsan:
    label: "CFSAN SNP Pipeline"
    command: "cfsan_snp_pipeline run -m soft -o output reference.fasta"
    distance_matrix: cfsan/snp_distance_matrix.tsv  # pre-computed distance matrix (optional)
    core_snps: cfsan/snplist.txt                     # core SNP positions (recommended)
    samples:
      sample1:
        bam: cfsan/sample1/reads.sorted.bam
      sample2:
        bam: cfsan/sample2/reads.sorted.bam
      sample3:
        bam: cfsan/sample3/reads.sorted.bam

options:
  min_depth: 1
  min_qual: 0
  include_indels: false
```

### 2. Generate the comparison report

```bash
# JSON format (human-readable)
coreguard compare --config project.yaml -o report.json

# Compact JSON with gzip (smaller file)
coreguard compare --config project.yaml -o report.json.gz --gzip --compact

# Binary format with gzip compression (recommended, ~10x faster loading)
coreguard compare --config project.yaml -o report.bin.gz --binary --gzip
```

### 3. Visualize the results

Open [https://genpat-it.github.io/coreguard/](https://genpat-it.github.io/coreguard/) and drag & drop your report file.

Supported formats: `.json`, `.json.gz`, `.bin`, `.bin.gz`

## Configuration Reference

### Global Options

| Option | Description | Default |
|--------|-------------|---------|
| `min_depth` | Minimum read depth to consider position covered | 1 |
| `min_qual` | Minimum VCF QUAL score (0 = no filtering; recommended since variant callers already apply their own filters) | 0 |
| `min_consensus` | Minimum fraction of reads agreeing on a base for reference pileup (0.0-1.0) | 0.8 |
| `include_indels` | Include insertions/deletions | false |

### Pipeline Options

| Option | Description | Default |
|--------|-------------|---------|
| `label` | Display name in the viewer | pipeline ID |
| `command` | Command line used (shown in viewer) | - |
| `ground_truth` | Mark as reference alignment (BAM-only baseline) | false |
| `distance_matrix` | Path to pre-computed SNP distance matrix (TSV) | - |
| `core_snps` | Path to core SNP output (snippy `core.tab` or CFSAN `snplist.txt`) | - |

### Per-Sample Options

| Option | Description |
|--------|-------------|
| `bam` | Path to BAM file for coverage/gap detection (required) |
| `vcf` | Path to VCF file with SNP calls (optional, only needed if no `core_snps`) |
| `label` | Display name for the sample (optional) |

### Description

The `description` field accepts either inline markdown or a path to a `.md`/`.txt` file:

```yaml
# Inline
description: "## My Study\nComparing pipelines on *Listeria* outbreak data."

# File reference
description: docs/study_description.md
```

### Core SNPs

The `core_snps` field points to a pipeline's native core SNP output. When provided, CoreGuard uses it for both pipeline statistics and the "Reference Discriminating SNPs vs Pipelines" comparison, making per-sample VCF files unnecessary.

Supported formats:

- **Snippy**: `core.tab` — TSV with CHR, POS, REF, and per-sample alleles
- **CFSAN**: `snplist.txt` — one position per line (or `snpma.fasta` for alleles)

#### Adding a New Core SNP Parser

CoreGuard uses a plugin system (`CoreSnpParser` trait in `src/parsers/mod.rs`) for parsing pipeline-specific core SNP output files. To add support for a new format:

```rust
pub trait CoreSnpParser {
    fn format_name(&self) -> &str;
    fn can_parse(&self, path: &Path) -> bool;
    fn parse(&self, path: &Path) -> anyhow::Result<CoreSnpData>;
}
```

1. Implement the trait for your format
2. Register it in `parse_core_snps()` (add to the `parsers` vector)

Current implementations: `SnippyCoreTabParser` (Snippy `core.tab`), `CfsanSnplistParser` (CFSAN `snplist.txt`/`snpma.fasta`).

### Distance Matrix

The `distance_matrix` field points to a TSV file with pairwise SNP distances between samples, as produced by the pipeline itself. These are displayed in the viewer's "SNP Distance Matrix" panel.

## Output Formats

| Format | Extension | Pros | Cons |
|--------|-----------|------|------|
| JSON | `.json` | Human-readable, debuggable | Large files, slower parsing |
| JSON gzip | `.json.gz` | Smaller download | Still slower parsing |
| Binary | `.bin` | Fast parsing (~10x) | Not human-readable |
| Binary gzip | `.bin.gz` | Fast + small | Not human-readable |

**Recommendation**: Use `--binary --gzip` for production, JSON for debugging.

## Technology

- **CLI**: Rust (fast BAM/core_snps processing, pre-computes all KPIs)
- **Viewer**: Vanilla JavaScript dashboard
- **WASM**: Rust compiled to WebAssembly (for backward compatibility with v1 reports)
- **No backend required**: All processing happens client-side

## Dependencies

The CLI requires `samtools` to be installed and available in PATH for BAM depth calculation.

## Citation

*Paper in preparation*

## License

MIT License - see [LICENSE](LICENSE)

## Author

GenPat Team

Contact: genpat@izs.it
