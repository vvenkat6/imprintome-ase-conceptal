#!/bin/bash
#SBATCH --job-name=ase_counter
#SBATCH --output=slurm-logs/ase_counter_%A_%a.out
#SBATCH --error=slurm-logs/ase_counter_%A_%a.err
#SBATCH --array=1-6
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
# ============================================================================
# Allelic read counting at heterozygous SNPs with GATK ASEReadCounter
# ============================================================================
# For each of six samples, counts reference and alternate allele reads at
# every heterozygous SNP in the cohort SNP VCF, using the WASP-filtered BAM
# (from 09_wasp_filter_merge.slurm) as input.
#
# Output is a per-sample table with columns:
#   contig, position, variantID, refAllele, altAllele, refCount, altCount,
#   totalCount, lowMAPQDepth, lowBaseQDepth, rawDepth, otherBases,
#   improperPairs
#
# These per-sample SNP-level tables are the input to gene-level aggregation
# (scripts/snp_to_gene_mapping.R and scripts/gene_level_ase_stats.R).
#
# ASEReadCounter thresholds:
#   --min-depth 10             minimum total depth to report a SNP
#   --min-mapping-quality 20   minimum MAPQ per read (excludes multi-mappers
#                              with MAPQ 0 and low-quality alignments)
#   --min-base-quality 20      minimum base quality at the SNP position
#
# References:
#   McKenna A, et al. Genome Res. 2010;20:1297-1303.
#   Castel SE, et al. Tools and best practices for data processing in
#     allelic expression analysis. Genome Biol. 2015;16:195.
#
# Sample manifest: config/samples.tsv (columns: sample_id, tissue, conceptus)
# ============================================================================

set -euo pipefail

# ---- Module environment ----------------------------------------------------
module load gatk/4.4.0.0
module load samtools/1.17

# ---- Configuration ---------------------------------------------------------
PROJECT_ROOT=<!-- FILL IN: /path/to/project/root -->
WORKDIR=${PROJECT_ROOT}/analysis
WASP_DIR=${WORKDIR}/wasp/final
VCF_DIR=${WORKDIR}/variants
OUTDIR=${WORKDIR}/ase_counts
REF=${PROJECT_ROOT}/reference/GRCh38.primary_assembly.genome.fa
COHORT_VCF=${VCF_DIR}/cohort.snps.vcf.gz
SAMPLES_TSV=${PROJECT_ROOT}/config/samples.tsv

mkdir -p "${OUTDIR}" slurm-logs

# ---- Resolve sample for this array task -----------------------------------
SAMPLE=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==i {print $1}' \
         "${SAMPLES_TSV}")

if [[ -z "${SAMPLE}" ]]; then
    echo "[ERROR] Could not resolve sample for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

echo "[$(date)] Starting ASEReadCounter for sample: ${SAMPLE}"

# ---- Input and output paths -----------------------------------------------
WASP_BAM=${WASP_DIR}/${SAMPLE}.wasp.sorted.bam
OUTTAB=${OUTDIR}/${SAMPLE}.ase_counts.tsv

# ---- Preflight checks -----------------------------------------------------
if [[ ! -f "${WASP_BAM}" ]]; then
    echo "[ERROR] WASP-filtered BAM not found: ${WASP_BAM}" >&2
    echo "  Run 09_wasp_filter_merge.slurm before this script." >&2
    exit 2
fi

if [[ ! -f "${COHORT_VCF}" ]]; then
    echo "[ERROR] Cohort SNP VCF not found: ${COHORT_VCF}" >&2
    echo "  Run 06_joint_genotype.slurm before this script." >&2
    exit 3
fi

if [[ ! -f "${COHORT_VCF}.tbi" ]]; then
    echo "[ERROR] Cohort SNP VCF tabix index not found: ${COHORT_VCF}.tbi" >&2
    exit 4
fi

if [[ ! -f "${REF}" ]]; then
    echo "[ERROR] Reference FASTA not found: ${REF}" >&2
    exit 5
fi

# ---- Ensure BAM index present ---------------------------------------------
if [[ ! -f "${WASP_BAM}.bai" ]] && [[ ! -f "${WASP_BAM%.bam}.bai" ]]; then
    echo "[INFO] BAM index not found; creating with samtools"
    samtools index "${WASP_BAM}"
fi

echo "  Input BAM: ${WASP_BAM}"
echo "  Cohort SNP VCF: ${COHORT_VCF}"
echo "  Reference: ${REF}"
echo "  Output table: ${OUTTAB}"

# ---- Run ASEReadCounter ---------------------------------------------------
# --min-depth 10:               require at least 10 total reads at the SNP
#                               position after filtering
# --min-mapping-quality 20:     exclude reads with MAPQ < 20 (removes STAR-
#                               assigned MAPQ 0 multimappers and low-quality
#                               alignments)
# --min-base-quality 20:        require Phred base quality >= 20 at the SNP
#                               position
# --verbosity ERROR:            suppress verbose per-site info messages
#
# GATK ASEReadCounter automatically:
#   - Uses only reads passing standard duplicate-read filter (MarkDuplicates
#     flag 0x400 excluded)
#   - Considers only heterozygous genotypes in the input VCF (0/1 or 1/0)
#     for this sample; homozygous sites are skipped
#   - Reports separately reads that failed each filter (lowMAPQDepth,
#     lowBaseQDepth) for QC inspection
echo "[$(date)] Running GATK ASEReadCounter"
gatk --java-options "-Xmx12g" ASEReadCounter \
     -R "${REF}" \
     -I "${WASP_BAM}" \
     -V "${COHORT_VCF}" \
     -O "${OUTTAB}" \
     --min-depth 10 \
     --min-mapping-quality 20 \
     --min-base-quality 20 \
     --verbosity ERROR

# ---- Verify output --------------------------------------------------------
if [[ ! -f "${OUTTAB}" ]]; then
    echo "[ERROR] ASEReadCounter did not produce output: ${OUTTAB}" >&2
    exit 6
fi

# ---- Report summary -------------------------------------------------------
# Table has one header line + one row per informative SNP
N_ROWS=$(wc -l < "${OUTTAB}")
N_SNPS=$((N_ROWS - 1))

echo "[$(date)] Done ${SAMPLE}"
echo "  Output table: ${OUTTAB}"
echo "  Informative SNPs (passed all filters): ${N_SNPS}"

# ============================================================================
# Next steps (R scripts in scripts/, not SLURM):
#   - snp_to_gene_mapping.R      map SNPs to gene bodies, exclude multi-gene
#   - gene_level_ase_stats.R     per-gene aggregation, binomial + Fisher + BH
#   - icr_annotation.R           annotate with tier 1/2 ICR overlap and distance
#   - prioritization_scoring.R   compute prioritization score, assign tiers
#   - benchmark_audit.R          cross-reference geneimprint catalog
# ============================================================================
