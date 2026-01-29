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
- **Consensus SNPs**: Find high-confidence variants where all pipelines agree.

## Features

### Core Features
- **Multi-pipeline comparison**: Compare SNP calls from any number of pipelines side-by-side
- **Ground truth support** (optional): Designate a baseline BAM alignment for gap detection and validation
- **Interactive viewer**: Navigate through genomic positions with zoom/pan on HTML5 Canvas
- **Client-side processing**: Viewer runs entirely in the browser via WebAssembly (no server needed)

### SNP Analysis
- **Smart filtering**: Filter by consensus, discordant calls, gaps, and pipeline-specific SNPs
- **MNP decomposition**: Multi-nucleotide polymorphisms (e.g., `TTGGCG→CCGGCT`) are automatically decomposed into individual SNPs
- **Distance matrix**: Calculate pairwise SNP distances between samples per pipeline
- **Polymorphic sites**: Identify positions where samples differ from each other

### Statistics & KPIs
- **Per-Pipeline Statistics**: SNPs, gaps, #All (core), #Consensus counts per pipeline
- **Per-Sample Statistics**: Detailed breakdown per sample with agreement metrics
- **Cross-Pipeline Concordance**: Pairwise comparison of SNP positions between pipelines
- **Cross-Pipeline Gap Analysis**: SNPs from each pipeline falling in gap regions of other pipelines
- **Coverage Statistics**: Depth, quality, and consensus metrics per sample/pipeline

### Ground Truth Comparison
- **SNPs in GT Gaps**: Detect SNPs called in low-coverage regions of the ground truth
- **GT SNPs Filtered**: Identify ground truth SNPs missed by other pipelines
- **BAM Pileup Warning**: Visual indicator when GT SNPs come from BAM pileup (no variant calling)

### User Interface
- **Dark/Light theme**: Toggle between themes
- **Collapsible panels**: Organize information in expandable sections
- **Info icons**: Hover/click for detailed explanations of each metric
- **Export options**: Download distance matrices and statistics

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

samples:
  sample1: {}
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

  # SNP pipeline with VCF + BAM
  snippy:
    label: "Snippy v4.6"
    command: "snippy --ref reference.fa --R1 reads_1.fq.gz --R2 reads_2.fq.gz --outdir out"
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
    samples:
      sample1:
        vcf: cfsan/sample1/var.flt.vcf
      sample2:
        vcf: cfsan/sample2/var.flt.vcf
      sample3:
        vcf: cfsan/sample3/var.flt.vcf

options:
  min_depth: 1
  min_qual: 20
  include_indels: false
```

### 2. Generate the comparison report

```bash
# JSON format (human-readable)
coreguard compare --config project.yaml -o report.json

# Binary format with gzip compression (recommended for large datasets, ~10x faster loading)
coreguard compare --config project.yaml -o report.bin.gz --binary --gzip

# Compact JSON with gzip
coreguard compare --config project.yaml -o report.json.gz --gzip --compact
```

### 3. Visualize the results

Open [https://genpat-it.github.io/coreguard/](https://genpat-it.github.io/coreguard/) and drag & drop your report file.

Supported formats: `.json`, `.json.gz`, `.bin`, `.bin.gz`

## Configuration Options

| Option | Description | Default |
|--------|-------------|---------|
| `min_depth` | Minimum read depth to consider position covered | 1 |
| `min_qual` | Minimum SNP quality score to include | 20 |
| `include_indels` | Include insertions/deletions | false |

### Pipeline Options

| Option | Description | Default |
|--------|-------------|---------|
| `ground_truth` | Mark pipeline as ground truth (BAM-only baseline) | false |
| `label` | Display name in the viewer | pipeline ID |
| `command` | Command line used (shown in viewer for reproducibility) | - |

### Per-Sample Options

| Option | Description |
|--------|-------------|
| `vcf` | Path to VCF file with SNP calls |
| `bam` | Path to BAM file for coverage/gap detection |

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

Andrea De Ruvo - GenPat Team, IZSNT

Contact: genpat@izs.it
