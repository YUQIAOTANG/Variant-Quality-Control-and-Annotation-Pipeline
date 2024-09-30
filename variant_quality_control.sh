#!/bin/bash

# Usage: ./variant_quality_control.sh input_vcf.gz filtered_output_vcf.gz snpEff_output_vcf.gz snpSift_output_vcf.gz

# 1. Variant-level quality control:
#   a. Retain variants that meet the following criteria:
#      i.   Maximum missing genotype rate of 10%
#      ii.  Biallelic variants
#      iii. Minimum read depth ≥ 10X
#      iv.  Filter status is PASS, genotype quality ≥ 90
#      v.   SNP-specific filters: QD > 2.0, FS < 60, MQRankSum > -12.5, ReadPosRankSum > -8.0, SOR ≤ 3, InbreedingCoeff > -0.8
#      vi.  Indel-specific filters: QD > 2.0, FS < 200, ReadPosRankSum > -20.0, SOR < 10.0
#      vii. Hardy-Weinberg Equilibrium p-value ≥ 0.001
#   b. Annotate variants using snpEff and SnpSift

# Binary definitions
SNPEFF_BIN=./snpEff/snpEff.jar
SNPEFF_DATABASE=GRCh38.86
SNPSIFT_BIN=./snpEff/SnpSift.jar
SNPSIFT_DATABASE=dbNSFP4.7.gz # Adjust this database if using the hg19 reference genome

# Input parameter definitions
INPUT_VCF_GZ=${1}
QC_VCF_GZ=${2}
SNPEFF_ANNOTATED_VCF=${3}
SNPSIFT_ANNOTATED_VCF=${4}

# Quality control filters using bcftools for improved efficiency and speed
bcftools view -m2 -M2 -v snps,indels -f PASS -Ou $INPUT_VCF_GZ | \
# Set genotypes with DP<10 or GQ<90 to missing
bcftools filter -e 'FMT/DP<10 || FMT/GQ<90' -S . -Ou | \
# Calculate missingness and filter variants with missing rate >10%
bcftools +fill-tags -Ob -o temp.vcf.gz -- -t all,F_MISSING | \
bcftools view -i 'F_MISSING<=0.1' -Ou temp.vcf.gz | \
# Calculate HWE and filter variants with p-value < 0.001
bcftools +fill-tags -Ob -o temp2.vcf.gz -- -t all,HWE | \
bcftools view -i 'HWE>=0.001' -Ou temp2.vcf.gz | \
# Apply SNP and Indel specific filters
bcftools view -i '((TYPE="snp" && QD>2.0 && FS<60 && MQRankSum>-12.5 && ReadPosRankSum>-8.0 && SOR<=3 && InbreedingCoeff>-0.8) || (TYPE="indel" && QD>2.0 && FS<200 && ReadPosRankSum>-20.0 && SOR<10.0))' -Oz -o $QC_VCF_GZ

# Clean up temporary files
rm temp.vcf.gz temp2.vcf.gz

# Annotation with snpEff and SnpSift
java -Xmx40g -jar $SNPEFF_BIN $SNPEFF_DATABASE $QC_VCF_GZ > $SNPEFF_ANNOTATED_VCF
java -Xmx40g -jar $SNPSIFT_BIN dbnsfp -v -db $SNPSIFT_DATABASE \
-f hg19_chr,"hg19_pos(1-based)",aaref,aaalt,aapos,genename,Ancestral_allele,SIFT_pred,Polyphen2_HVAR_pred,CADD_phred,gnomAD_exomes_NFE_AF,gnomAD_exomes_POPMAX_AF,gnomAD_genomes_NFE_AF,gnomAD_genomes_POPMAX_AF,clinvar_hgvs,Interpro_domain,GTEx_V7_gene,GTEx_V7_tissue \
-a -m $SNPEFF_ANNOTATED_VCF  > $SNPSIFT_ANNOTATED_VCF

# Compress and index the final annotated VCF
bgzip $SNPSIFT_ANNOTATED_VCF
tabix "${SNPSIFT_ANNOTATED_VCF}.gz"
