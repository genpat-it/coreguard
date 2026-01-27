# CoreGuard

Pipeline-agnostic SNP comparison tool for bacterial genomics.

## Overview

CoreGuard compares SNP calls from multiple pipelines (Snippy, CFSAN, GATK, etc.), helping identify discrepancies and artifacts in variant calling. Optionally, a ground truth alignment can be used to detect coverage gaps.

**Live Viewer:** [https://genpat-it.github.io/coreguard/](https://genpat-it.github.io/coreguard/)

## Why CoreGuard?

Different SNP pipelines can produce different results on the same data. This matters because:

- **Outbreak investigation**: SNP distances determine if isolates belong to the same outbreak. A difference of 1-2 SNPs can change epidemiological conclusions.
- **Pipeline validation**: When adopting a new pipeline, you need to verify it produces consistent results with established methods.
- **Troubleshooting**: When results seem wrong, CoreGuard helps visualize exactly where and why pipelines disagree.
- **Quality control**: Identify problematic samples or genomic regions where pipelines consistently fail.
- **Consensus SNPs**: Find high-confidence variants where all pipelines agree.

## Features

- **Multi-pipeline comparison**: Compare SNP calls from different pipelines side-by-side
- **Ground truth support** (optional): Designate a baseline BAM alignment for gap detection
- **Interactive viewer**: Navigate through genomic positions with zoom/pan
- **Smart filtering**: Filter by consensus, discordant calls, gaps, and pipeline-specific SNPs
- **Distance matrix**: Calculate pairwise SNP distances between samples
- **Client-side processing**: Viewer runs entirely in the browser via WebAssembly

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

  # Another pipeline (VCF only)
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
coreguard compare --config project.yaml -o report.json
```

### 3. Visualize the results

Open [https://genpat-it.github.io/coreguard/](https://genpat-it.github.io/coreguard/) and load your `report.json` file.

## Configuration Options

| Option | Description | Default |
|--------|-------------|---------|
| `min_depth` | Minimum read depth to consider position covered | 1 |
| `min_qual` | Minimum SNP quality score to include | 20 |
| `include_indels` | Include insertions/deletions | false |
| `ground_truth` | Mark pipeline as ground truth (BAM-only baseline) | false |
| `command` | Command line used to run the pipeline (shown in viewer) | - |

## Technology

- **CLI**: Rust
- **Viewer**: Vanilla JavaScript + HTML5 Canvas
- **Computation**: Rust compiled to WebAssembly (runs in browser)

## Citation

*Paper in progress*

## License

MIT License

## Author

GenPat Team - genpat@izs.it
