# Variant Quality Control and Annotation Pipeline 2024 summer

This repository contains a comprehensive pipeline developed for performing high-quality variant calling and annotation of genetic data stored in Variant Call Format (VCF) files. The pipeline is designed to meet the variant-level and sample-level quality control (QC) standards commonly required for downstream genetic analysis, including disease association studies and population genomics.

## Project Overview

In the Dr. Gupta Laboratory, we focus on leveraging bioinformatics tools to analyze next-generation sequencing (NGS) data, particularly for identifying and annotating genetic variants. These variants can provide insights into disease mechanisms, ancestry, and population genetics. This project implements a streamlined workflow for handling VCF files, performing variant-level QC, and annotating variants using tools like `snpEff` and `SnpSift`.

## Pipeline Steps

The pipeline includes multiple steps, each designed to filter and annotate genetic variants to ensure the highest quality for downstream analysis:

### 1. Variant-Level Quality Control

This step applies various filters to the variants to ensure high-quality calls. The filters include:

- **Maximum missing genotype rate of 10%**: Variants with more than 10% missing genotypes are excluded.
- **Biallelic variants only**: Only variants with two alleles (biallelic) are retained.
- **Minimum read depth ≥ 10X**: Variants with fewer than 10 reads across samples are filtered out.
- **Genotype quality ≥ 90**: Variants with a genotype quality score lower than 90 are excluded.
- **SNP-specific filters**:
  - **QD > 2.0**: Quality by depth ratio.
  - **FS < 60**: Phred-scaled p-value using Fisher's exact test for strand bias.
  - **MQRankSum > -12.5**: Test for strand bias of mapping quality.
  - **ReadPosRankSum > -8.0**: Test for strand bias of read position.
  - **SOR ≤ 3**: Symmetry of Odds Ratio test for strand bias.
  - **Inbreeding Coefficient > -0.8**: Ensures variants have a reasonable inbreeding coefficient.
- **Indel-specific filters**:
  - **QD > 2.0**: Quality by depth ratio.
  - **FS < 200**: Fisher strand test for strand bias.
  - **ReadPosRankSum > -20.0**: Test for bias in read position at the ends of reads.
  - **SOR < 10.0**: Symmetry of Odds Ratio test for strand bias.
- **Hardy-Weinberg Equilibrium (HWE) p-value ≥ 0.001**: Ensures that variants are in Hardy-Weinberg equilibrium.

### 2. Variant Annotation

The pipeline uses two tools to annotate the filtered variants:

- **snpEff**: Annotates variants by predicting their effects on genes, transcripts, and proteins based on the GRCh38 reference genome.
- **SnpSift**: Further refines the annotations using the dbNSFP database to provide detailed information on functional predictions, population frequencies, and clinical relevance.

### 3. Sample-Level Quality Control

This step filters samples based on various criteria to ensure high-quality and representative data:

- **PCA and Admixture Analysis**: Identifies individuals with mixed ancestry and removes them from the analysis.
- **Kinship Filtering**: Removes related individuals by excluding one of each pair with a kinship coefficient > 0.0884.
- **Sex Inference**: Uses Y chromosome data (when available) to infer and verify the biological sex of samples.
- **Outlier Detection**: Removes samples with extreme values in key metrics, defined as samples falling outside mean ± 3 * standard deviation.
  - **Key metrics include**:
    - **Expected Ts/Tv ratio**: Expected transition/transversion ratios for whole-exome sequencing (WES) and whole-genome sequencing (WGS).
    - **Ancestral vs derived Ts/Tv ratio**: A separate ratio for ancestral and derived alleles.

### 4. Compression and Indexing

After the annotations are completed, the final VCF file is compressed and indexed using `bgzip` and `tabix` for efficient storage and querying.

## Example Command

To run the pipeline, use the following command:

```bash
./variant_quality_control.sh input_vcf.gz filtered_output_vcf.gz snpEff_output_vcf.gz snpSift_output_vcf.gz

```
