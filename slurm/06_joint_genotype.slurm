#!/bin/bash
#SBATCH --job-name=joint_gt
#SBATCH --output=slurm-logs/joint_gt_%j.out
#SBATCH --error=slurm-logs/joint_gt_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
# ============================================================================
# Cohort joint genotyping and SNP selection
# ============================================================================
# Combines per-sample GVCFs (from 05_haplotypecaller_gvcf.slurm) into a
# single cohort VCF, performs joint genotyping, and selects biallelic SNPs
# for downstream ASE analysis.
#
# Three GATK steps:
#   1. CombineGVCFs: merges per-sample GVCFs into a single multi-sample GVCF
#   2. GenotypeGVCFs: performs joint genotyping across the cohort
#   3. SelectVariants: restricts to biallelic SNPs
#
# This is not an array job (single cohort-wide operation across all samples).
#
# References:
#   McKenna A, et al. Genome Res. 2010;20:1297-1303.
#   DePristo MA, et al. Nat Genet. 2011;43:491-498.
#   Van der Auwera GA, et al. Curr Protoc Bioinformatics. 2013;43:11.10.1.
# ============================================================================

set -euo pipefail

# ---- Module environment ----------------------------------------------------
module load gatk/4.4.0.0
module load htslib/1.17         # bgzip, tabix

# ---- Configuration ---------------------------------------------------------
THREADS=${SLURM_CPUS_PER_TASK:-4}
PROJECT_ROOT=<!-- FILL IN: /path/to/project/root -->
WORKDIR=${PROJECT_ROOT}/analysis
GVCFDIR=${WORKDIR}/gvcf
OUTDIR=${WORKDIR}/variants
REF=${PROJECT_ROOT}/reference/GRCh38.primary_assembly.genome.fa

mkdir -p "${OUTDIR}" slurm-logs

echo "[$(date)] Starting cohort joint genotyping"

# ---- Preflight: locate per-sample GVCFs -----------------------------------
GVCFS=( $(ls -1 ${GVCFDIR}/*.g.vcf.gz) )

if [[ ${#GVCFS[@]} -eq 0 ]]; then
    echo "[ERROR] No per-sample GVCFs found in ${GVCFDIR}" >&2
    echo "  Run 05_haplotypecaller_gvcf.slurm before this script." >&2
    exit 1
fi

echo "  Found ${#GVCFS[@]} GVCF files:"
for f in "${GVCFS[@]}"; do
    echo "    ${f}"
done

if [[ ! -f "${REF}" ]]; then
    echo "[ERROR] Reference FASTA not found: ${REF}" >&2
    exit 2
fi

# ---- Output paths ---------------------------------------------------------
COMBINED=${OUTDIR}/cohort.g.vcf.gz
COHORT_VCF=${OUTDIR}/cohort.vcf.gz
SNP_VCF=${OUTDIR}/cohort.snps.vcf.gz

# ---- Step 1: CombineGVCFs -------------------------------------------------
# Merges per-sample GVCFs into a single multi-sample GVCF. For small cohorts
# (<10 samples), CombineGVCFs is preferred over GenomicsDBImport for
# simplicity. For larger cohorts, switch to GenomicsDBImport.
echo "[$(date)] Step 1: CombineGVCFs"
gatk --java-options "-Xmx24g" CombineGVCFs \
     -R "${REF}" \
     $(for f in "${GVCFS[@]}"; do echo -n "--variant ${f} "; done) \
     -O "${COMBINED}"

if [[ ! -f "${COMBINED}" ]]; then
    echo "[ERROR] CombineGVCFs did not produce output: ${COMBINED}" >&2
    exit 3
fi

# ---- Step 2: GenotypeGVCFs ------------------------------------------------
# Performs joint genotyping across the combined GVCF, producing a cohort VCF
# with genotype calls for every sample at every variant site.
echo "[$(date)] Step 2: GenotypeGVCFs"
gatk --java-options "-Xmx24g" GenotypeGVCFs \
     -R "${REF}" \
     -V "${COMBINED}" \
     -O "${COHORT_VCF}"

if [[ ! -f "${COHORT_VCF}" ]]; then
    echo "[ERROR] GenotypeGVCFs did not produce output: ${COHORT_VCF}" >&2
    exit 4
fi

N_ALL_VARIANTS=$(zcat "${COHORT_VCF}" | grep -vc "^#" || true)
echo "  Cohort VCF variants (all types): ${N_ALL_VARIANTS}"

# ---- Step 3: SelectVariants (biallelic SNPs only) -------------------------
# Restricts to biallelic SNP-type variants for downstream ASE analysis.
# --select-type-to-include SNP: excludes indels, mixed sites, symbolic alleles
# --restrict-alleles-to BIALLELIC: excludes multiallelic sites (needed for
#                                  clean ref/alt counting in ASEReadCounter)
echo "[$(date)] Step 3: SelectVariants (biallelic SNPs)"
gatk --java-options "-Xmx16g" SelectVariants \
     -R "${REF}" \
     -V "${COHORT_VCF}" \
     --select-type-to-include SNP \
     --restrict-alleles-to BIALLELIC \
     -O "${SNP_VCF%.gz}"

# ---- Compress and index ---------------------------------------------------
# SelectVariants may write uncompressed VCF depending on version behavior.
# bgzip + tabix ensures compressed, indexed output for downstream tools.
if [[ -f "${SNP_VCF%.gz}" ]]; then
    bgzip -f "${SNP_VCF%.gz}"
fi
tabix -f -p vcf "${SNP_VCF}"

if [[ ! -f "${SNP_VCF}" ]] || [[ ! -f "${SNP_VCF}.tbi" ]]; then
    echo "[ERROR] SNP VCF or its tabix index not produced: ${SNP_VCF}" >&2
    exit 5
fi

# ---- Also index the cohort VCF (pre-filter, for downstream audit) --------
if [[ ! -f "${COHORT_VCF}.tbi" ]]; then
    tabix -f -p vcf "${COHORT_VCF}"
fi

# ---- Report summary -------------------------------------------------------
N_SNPS=$(zcat "${SNP_VCF}" | grep -vc "^#" || true)
N_SAMPLES=$(zcat "${SNP_VCF}" | grep "^#CHROM" | awk -F'\t' '{print NF - 9}')
echo "[$(date)] Joint genotyping complete"
echo "  Combined GVCF: ${COMBINED}"
echo "  Cohort VCF (all variants): ${COHORT_VCF}"
echo "  SNP VCF (biallelic SNPs): ${SNP_VCF}"
echo "  Total biallelic SNPs: ${N_SNPS}"
echo "  Samples in VCF: ${N_SAMPLES}"

# ============================================================================
# Downstream use:
#   - Per-chromosome SNP tables for WASP (07_wasp_intersecting_snps.slurm)
#     built from ${SNP_VCF} via a helper script.
#   - Allelic read counting with ASEReadCounter (10_asereadcounter.slurm)
#     uses ${SNP_VCF} directly.
#   - The audit workflows (11_gene_body_audit.slurm and scripts/benchmark_audit.R)
#     query ${SNP_VCF} for het SNP counts per gene.
# ============================================================================
