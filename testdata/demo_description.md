# Listeria monocytogenes Outbreak Analysis

This demo compares SNP calls from two widely-used bacterial genomics pipelines on a set of *Listeria monocytogenes* isolates from an IZSAM 2020 outbreak investigation.

## Dataset

- **4 samples**: TE15676, TE15677, TE15678 (clonal group), TE8064 (outlier)
- **Reference genome**: NCBI assembly used for mapping
- **Sequencing**: Illumina paired-end reads

## Pipelines Compared

1. **Snippy v4.6** - Popular SNP calling pipeline using bwa + freebayes
2. **CFSAN SNP Pipeline** - FDA's pipeline for foodborne pathogen analysis
3. **Ground Truth (minimap2)** - BAM pileup without variant calling filters (baseline)

## Key Observations

The clonal samples (TE15676/77/78) show **0-1 SNP differences** in core genes, consistent with outbreak clustering. TE8064 is an **outlier (~300 SNPs distant)** representing a different strain.

This comparison reveals:
- Positions where pipelines **agree** (high-confidence SNPs)
- Positions where pipelines **disagree** (potential artifacts)
- **Gap regions** where coverage is insufficient for reliable calls

## Purpose

Use this demo to explore:
- How different pipelines can produce varying SNP counts
- The importance of ground truth validation
- How CoreGuard helps identify discrepancies
