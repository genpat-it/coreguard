# CoreGuard

Pipeline-agnostic SNP comparison tool for bacterial genomics.

## Overview

CoreGuard compares SNP calls from multiple pipelines (Snippy, CFSAN, GATK, etc.), helping identify discrepancies and artifacts in variant calling. Optionally, a ground truth alignment can be used to detect coverage gaps and validate SNP calls.

**Live Viewer:** [https://genpat-it.github.io/coreguard/](https://genpat-it.github.io/coreguard/)

### Demo

Click **"Load Demo"** in the viewer to explore a pre-loaded dataset:

- **Organism**: *Listeria monocytogenes* (4 samples)
- **Pipelines compared**: Snippy, CFSAN SNP Pipeline
- **Ground Truth**: minimap2 alignment (BAM pileup without variant calling)

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
- **Ground truth support** (optional): Designate a baseline BAM alignment for gap detection and validation
- **Interactive viewer**: Navigate through genomic positions with zoom/pan on HTML5 Canvas
- **Client-side processing**: Viewer runs entirely in the browser via WebAssembly (no server needed)

### SNP Analysis
- **Smart filtering**: Filter by consensus, discordant calls, gaps, and pipeline-specific SNPs
- **MNP decomposition**: Multi-nucleotide polymorphisms (e.g., `TTGGCG→CCGGCT`) are automatically decomposed into individual SNPs
- **Pre-computed distance matrices**: Display pairwise SNP distance matrices from each pipeline's native output
- **GT discriminating SNPs vs pipelines**: Compare ground truth discriminating positions against each pipeline's core SNP output

### Per-Pipeline KPI Dashboard
For each pipeline (GT and VCF), CoreGuard computes:

- **Gap-Intersect / Gap-Union**: Two gap-exclusion strategies
  - *Gap-Intersect*: exclude positions where ALL samples have a gap (permissive)
  - *Gap-Union*: exclude positions where ANY sample has a gap (restrictive)
- **Usable Space**: Reference length minus excluded gap positions
- **Total SNPs**: Positions in usable space where at least one sample has an alt allele
- **Consensus SNPs**: All samples agree on the same alt allele (non-discriminating)
- **Discriminating SNPs**: Samples differ — these contribute to Hamming distance
- **Missing VCF Calls**: Positions where some samples have a VCF call but others don't (variant calling inconsistency)

#### Discriminating SNP Breakdown (VCF pipelines)
For non-GT pipelines, discriminating SNPs are classified by heuristic:
- **Gap-affected**: At least one sample was skipped due to gap — partial comparison
- **GT-consensus**: GT pileup shows all samples agree — likely variant calling artifact
- **Majority-rule**: All but one sample agree — likely single-sample VC miss
- **Confirmed**: Genuine disagreement between samples

#### Pairwise Statistics
- **Avg Usable Space / Discriminating SNPs**: Averaged across all sample pairs (Gap-Union)
- **Min / Median / Max**: Distribution of pairwise discriminating SNPs
- **Per-Sample Table**: For each sample, average usable space and discriminating SNPs across all pairs involving that sample

#### BAM Pileup SNP Calling (Ground Truth)
The ground truth pipeline calls SNPs from BAM pileup using majority vote:
1. Count reads at each position (A, C, G, T)
2. If total depth < `min_depth` → **Gap**
3. If majority base < `min_consensus` (default 80%) → **Ambiguous** (skipped)
4. If majority base ≠ reference → **SNP**

### Viewer Panels
- **Description**: Project description in Markdown (from config)
- **Pipelines**: Pipeline metadata (labels, commands, data types)
- **Statistics**: KPI dashboard with GT metrics and GT discriminating SNPs vs pipelines
- **SNP Distance Matrix**: Pre-computed distance matrices from each pipeline
- **Genome Overview**: Canvas visualization of SNPs and gaps across the reference
- **View Settings**: Row visibility, nucleotide display options, legend
- **Filters & Navigation**: Filter by pipeline SNPs, consensus, discordant calls, gaps; AND/OR logic; go-to position

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
  # Ground truth alignment (BAM only, no VCF)
  minimap2:
    label: "Ground Truth (minimap2)"
    command: "minimap2 -ax sr -t 8 ref.fa reads_1.fq.gz reads_2.fq.gz | samtools sort -o out.bam"
    ground_truth: true
    samples:
      sample1:
        bam: alignments/sample1.bam
      sample2:
        bam: alignments/sample2.bam
      sample3:
        bam: alignments/sample3.bam

  # SNP pipeline with VCF + BAM + optional pipeline outputs
  snippy:
    label: "Snippy v4.6"
    command: "snippy --ref reference.fa --R1 reads_1.fq.gz --R2 reads_2.fq.gz --outdir out"
    distance_matrix: snippy/core.distances.tsv     # pre-computed distance matrix (optional)
    core_snps: snippy/core.tab                      # core SNP output file (optional)
    samples:
      sample1:
        vcf: snippy/sample1/snps.vcf
        bam: snippy/sample1/snps.bam
      sample2:
        vcf: snippy/sample2/snps.vcf
        bam: snippy/sample2/snps.bam
      sample3:
        vcf: snippy/sample3/snps.vcf
        bam: snippy/sample3/snps.bam

  # Another pipeline (VCF only - no gap info)
  cfsan:
    label: "CFSAN SNP Pipeline"
    command: "cfsan_snp_pipeline run -m soft -o output reference.fasta"
    distance_matrix: cfsan/snp_distance_matrix.tsv  # pre-computed distance matrix (optional)
    core_snps: cfsan/snplist.txt                     # core SNP positions (optional)
    samples:
      sample1:
        vcf: cfsan/sample1/var.flt.vcf
      sample2:
        vcf: cfsan/sample2/var.flt.vcf
      sample3:
        vcf: cfsan/sample3/var.flt.vcf

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
| `min_consensus` | Minimum fraction of reads agreeing on a base for GT pileup (0.0-1.0) | 0.8 |
| `include_indels` | Include insertions/deletions | false |

### Pipeline Options

| Option | Description | Default |
|--------|-------------|---------|
| `label` | Display name in the viewer | pipeline ID |
| `command` | Command line used (shown in viewer) | - |
| `ground_truth` | Mark as ground truth (BAM-only baseline) | false |
| `distance_matrix` | Path to pre-computed SNP distance matrix (TSV) | - |
| `core_snps` | Path to core SNP output (snippy `core.tab` or CFSAN `snplist.txt`) | - |

### Per-Sample Options

| Option | Description |
|--------|-------------|
| `vcf` | Path to VCF file with SNP calls |
| `bam` | Path to BAM file for coverage/gap detection |
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

The `core_snps` field points to a pipeline's native core SNP output. Supported formats:

- **Snippy**: `core.tab` — TSV with CHR, POS, REF, and per-sample alleles
- **CFSAN**: `snplist.txt` — one position per line

These are used to compute the "GT Discriminating SNPs vs Pipelines" comparison, showing how many ground truth discriminating positions each pipeline captures.

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

- **CLI**: Rust (fast VCF/BAM processing)
- **Viewer**: Vanilla JavaScript + HTML5 Canvas
- **Computation**: Rust compiled to WebAssembly (runs entirely in browser)
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
