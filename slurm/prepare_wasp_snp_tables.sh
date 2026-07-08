#!/bin/bash
# ============================================================================
# Generate per-chromosome SNP tables for WASP from cohort.snps.vcf.gz
# ============================================================================
# WASP's find_intersecting_snps.py expects SNP data in per-chromosome
# text files (one file per chromosome), gzipped, with three tab-separated
# columns: position, reference allele, alternate allele. This script
# generates that format from the joint-genotyped SNP VCF.
#
# Run once after 06_joint_genotype.slurm and before 07_wasp_intersecting_snps.slurm.
# Fast; takes <5 minutes.
# ============================================================================

set -euo pipefail

module load bcftools/1.17

PROJECT_ROOT=<!-- FILL IN: /path/to/project/root -->
VCF=${PROJECT_ROOT}/analysis/variants/cohort.snps.vcf.gz
OUTDIR=${PROJECT_ROOT}/analysis/wasp/snp_tables

mkdir -p "${OUTDIR}"

# Extract per-chromosome SNP tables
# Format: position<TAB>ref<TAB>alt (one file per chromosome, gzipped)
for CHR in chr{1..22} chrX chrY chrM; do
    echo "Extracting SNPs for ${CHR}..."
    bcftools view -r "${CHR}" "${VCF}" 2>/dev/null \
        | bcftools query -f '%POS\t%REF\t%ALT\n' \
        | gzip -c \
        > "${OUTDIR}/${CHR}.snps.txt.gz"
    N=$(zcat "${OUTDIR}/${CHR}.snps.txt.gz" | wc -l)
    echo "  ${CHR}: ${N} SNPs"
done

echo "Done. Per-chromosome SNP tables written to ${OUTDIR}"
