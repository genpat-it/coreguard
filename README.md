# CoreGuard

Pipeline-agnostic SNP comparison tool for bacterial genomics.

## Overview

CoreGuard compares SNP calls from multiple pipelines (Snippy, CFSAN, GATK, etc.), helping identify discrepancies and artifacts in variant calling. Optionally, a reference alignment can be used to detect coverage gaps and validate SNP calls.

**Live Viewer:** [https://genpat-it.github.io/coreguard/](https://genpat-it.github.io/coreguard/)

### Demo Datasets

The viewer includes 5 pre-loaded demo datasets. Click the demo buttons to explore:

| Demo | Organism | Samples | Pipelines | Description |
|------|----------|---------|-----------|-------------|
| **Demo 1** | *Listeria monocytogenes* | 4 | Snippy, CFSAN, SPANDx | Small dataset, clone group + outlier (EGD-e ref) |
| **Demo 2** | *Listeria monocytogenes* | 53 | Snippy, CFSAN, SPANDx | Large outbreak dataset |
| **Demo 3** | *Brucella melitensis* | 17 | Snippy, CFSAN, SPANDx | Different organism |
| **Demo 4** | West Nile Virus | 10 | Snippy, CFSAN, SPANDx | Viral genome (small reference) |
| **Demo 5** | *Listeria monocytogenes* | 4 | Snippy, CFSAN, SPANDx | Same as Demo 1 with F2365 ref (different serovar) |

All demos use **minimap2** as the reference alignment (BAM pileup without variant calling) for ground truth comparison, plus **Snippy v4.6.0**, **CFSAN SNP Pipeline v2.2.1**, and **SPANDx v4.0.5** as variant calling pipelines.

The demos illustrate how different pipelines can produce varying SNP calls on the same data, and how CoreGuard helps identify these discrepancies. Demo 5 specifically demonstrates the critical impact of reference genome choice on SNP calling resolution.

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
    reference: true
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
# Dashboard mode (recommended) - lightweight output for web viewing
coreguard compare --config project.yaml -o report.json.gz --dashboard

# Full mode - includes raw data for advanced analysis
coreguard compare --config project.yaml -o report.json.gz

# JSON format (human-readable, for debugging)
coreguard compare --config project.yaml -o report.json
```

**Dashboard mode** (`--dashboard`) omits raw gap/SNP data and keeps only pre-computed statistics. This reduces file size by 10-100x while maintaining full dashboard functionality.

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
| `reference` | Mark as reference alignment (BAM-only baseline) | false |
| `gaps_only` | Only load gaps from BAM, skip SNP pileup (reduces output size) | false |
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

### Core SNPs Format

The `core_snps` field points to a **CoreGuard TSV** file — a simple tab-separated format with genomic positions and per-sample alleles:

```tsv
CHR	POS	REF	sample1	sample2	sample3	sample4
CP014790.1	16686	C	C	C	T	C
CP014790.1	17153	C	C	T	C	C
CP014790.1	23349	T	T	T	A	T
```

- `CHR` — Chromosome/contig name
- `POS` — Genomic position (1-based)
- `REF` — Reference allele
- `<sample>` — Sample allele (`-` for gap, `N` for no call)

**Snippy** `core.tab` output is already in this format — use it directly.

**CFSAN** requires conversion:

```bash
coreguard convert --from cfsan-snpma \
  -i snpma.fasta \
  --snplist snplist.txt \
  --reference-snp referenceSNP.fasta \
  -o cfsan_core_snps.tsv
```

**SPANDx** — convert GATK VariantsToTable output:

```bash
coreguard convert --from spandx-vcf-table \
  -i out.vcf.table \
  -o spandx_core_snps.tsv
```

Then reference the converted file in your YAML:

```yaml
pipelines:
  cfsan:
    label: "CFSAN 2.2.1"
    core_snps: cfsan_core_snps.tsv
  spandx:
    label: "SPANDx v4.0.5"
    core_snps: spandx_core_snps.tsv
```

### Distance Matrix

The `distance_matrix` field points to a TSV file with pairwise SNP distances between samples. CoreGuard can compute this automatically from various sources:

**From Nexus SNP matrix** (e.g., SPANDx `Ortho_SNP_matrix.nex`):

```bash
coreguard convert --from nexus \
  -i Ortho_SNP_matrix.nex \
  -o distances.tsv
```

**From CoreGuard TSV core_snps**:

```bash
coreguard convert --from core-snps \
  -i core_snps.tsv \
  -o distances.tsv
```

These computed distances are displayed in the viewer's "SNP Distance Matrix" panel.

## Output Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| JSON gzip | `.json.gz` | Recommended for production |
| JSON | `.json` | Human-readable, for debugging |

### Dashboard vs Full Mode

| Mode | Flag | Output Size | Use Case |
|------|------|-------------|----------|
| **Dashboard** | `--dashboard` | 10-100x smaller | Web viewing, sharing |
| **Full** | (default) | Large | Advanced analysis, raw data access |

**Recommendation**: Use `--dashboard` for web viewing. Use full mode only if you need raw gap/SNP data for custom analysis.

## Technology

- **CLI**: Rust with rayon parallelization (fast BAM/core_snps processing, pre-computes all KPIs)
- **Viewer**: Vanilla JavaScript dashboard with pre-computed statistics
- **No backend required**: All processing happens at CLI time; viewer loads instantly

## Dependencies

The CLI requires `samtools` to be installed and available in PATH for BAM depth calculation.

## Citation

*Paper in preparation*

## License

MIT License - see [LICENSE](LICENSE)

## Author

GenPat Team

Contact: genpat@izs.it
