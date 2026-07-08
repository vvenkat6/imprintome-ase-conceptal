#!/bin/bash
#SBATCH --job-name=rg_dedup
#SBATCH --output=slurm-logs/rg_dedup_%A_%a.out
#SBATCH --error=slurm-logs/rg_dedup_%A_%a.err
#SBATCH --array=1-6
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=06:00:00
# ============================================================================
# Add read groups and mark PCR/optical duplicates
# ============================================================================
# For each of six samples, adds read group metadata required by GATK
# HaplotypeCaller (AddOrReplaceReadGroups), then marks PCR and optical
# duplicates with MarkDuplicates. Deduplicated BAM files feed variant
# discovery (05_haplotypecaller_gvcf.slurm) and WASP correction
# (07_wasp_intersecting_snps.slurm).
#
# References:
#   McKenna A, et al. The Genome Analysis Toolkit: a MapReduce framework for
#     analyzing next-generation DNA sequencing data. Genome Res. 2010;20:1297-1303.
#   DePristo MA, et al. A framework for variation discovery and genotyping
#     using next-generation DNA sequencing data. Nat Genet. 2011;43:491-498.
#
# Sample manifest: config/samples.tsv (columns: sample_id, tissue, conceptus)
# ============================================================================

set -euo pipefail

# ---- Module environment ----------------------------------------------------
module load gatk/4.4.0.0
module load samtools/1.17

# ---- Configuration ---------------------------------------------------------
THREADS=${SLURM_CPUS_PER_TASK:-4}
PROJECT_ROOT=<!-- FILL IN: /path/to/project/root -->
WORKDIR=${PROJECT_ROOT}/analysis
STAR_DIR=${WORKDIR}/star/secondpass
OUTDIR=${WORKDIR}/dedup
SAMPLES_TSV=${PROJECT_ROOT}/config/samples.tsv

mkdir -p "${OUTDIR}" slurm-logs

# ---- Resolve sample for this array task -----------------------------------
SAMPLE=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==i {print $1}' \
         "${SAMPLES_TSV}")

if [[ -z "${SAMPLE}" ]]; then
    echo "[ERROR] Could not resolve sample for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

echo "[$(date)] Starting read-group and dedup for sample: ${SAMPLE}"

# ---- Input and output paths -----------------------------------------------
INBAM=${STAR_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam
TMPBAM=${OUTDIR}/${SAMPLE}.rg.bam
DEDUP_BAM=${OUTDIR}/${SAMPLE}.dedup.bam
METRICS=${OUTDIR}/${SAMPLE}.dedup.metrics.txt

if [[ ! -f "${INBAM}" ]]; then
    echo "[ERROR] Input BAM not found: ${INBAM}" >&2
    echo "  Run 03_star_secondpass.slurm before this script." >&2
    exit 2
fi

echo "  Input BAM: ${INBAM}"
echo "  Output BAM: ${DEDUP_BAM}"
echo "  Metrics file: ${METRICS}"

# ---- Step 1: Add or replace read groups -----------------------------------
# GATK HaplotypeCaller requires @RG headers with matching read tags. STAR
# does not add these by default; AddOrReplaceReadGroups populates them.
#
# Read group fields:
#   RGID: read group ID (per-sample; must be unique across samples in the cohort)
#   RGLB: library ID (arbitrary; used by MarkDuplicates for library-aware dedup)
#   RGPL: platform (ILLUMINA)
#   RGPU: platform unit (flowcell.lane.barcode; arbitrary if not available)
#   RGSM: sample name (used as the column header in VCF; must match manifest)
echo "[$(date)] Step 1: AddOrReplaceReadGroups"
gatk --java-options "-Xmx20g" AddOrReplaceReadGroups \
     -I "${INBAM}" \
     -O "${TMPBAM}" \
     -RGID "${SAMPLE}" \
     -RGLB lib1 \
     -RGPL illumina \
     -RGPU unit1 \
     -RGSM "${SAMPLE}"

# ---- Step 2: Mark duplicates ----------------------------------------------
# MarkDuplicates flags PCR and optical duplicates by 5' alignment coordinate.
# Duplicates are not removed (flag bit 0x400 is set) so they can be excluded
# later with -F 260 in samtools or the DuplicateReadFilter in GATK.
#
# --VALIDATION_STRINGENCY SILENT: STAR occasionally emits alignments with
#   header inconsistencies that MarkDuplicates warns about; SILENT suppresses
#   those non-fatal warnings without affecting the result.
# --CREATE_INDEX true: writes a .bai alongside the output BAM
echo "[$(date)] Step 2: MarkDuplicates"
gatk --java-options "-Xmx20g" MarkDuplicates \
     -I "${TMPBAM}" \
     -O "${DEDUP_BAM}" \
     -M "${METRICS}" \
     --CREATE_INDEX true \
     --VALIDATION_STRINGENCY SILENT

# ---- Clean up intermediate BAM --------------------------------------------
rm -f "${TMPBAM}"

# ---- Fallback: ensure BAM is indexed --------------------------------------
if [[ ! -f "${DEDUP_BAM}.bai" ]] && [[ ! -f "${DEDUP_BAM%.bam}.bai" ]]; then
    echo "[INFO] MarkDuplicates did not create index; running samtools index"
    samtools index "${DEDUP_BAM}"
fi

# ---- Report duplicate metrics ---------------------------------------------
# MarkDuplicates metrics file format has a header block followed by a
# tab-separated table; extract the PERCENT_DUPLICATION field for this sample.
if [[ -f "${METRICS}" ]]; then
    PCT_DUP=$(grep -A 2 "^LIBRARY" "${METRICS}" | tail -1 | awk -F'\t' '{print $9}')
    echo "  Percent duplication: ${PCT_DUP}"
fi

N_READS=$(samtools view -c "${DEDUP_BAM}")
echo "[$(date)] Done ${SAMPLE}"
echo "  Deduplicated BAM: ${DEDUP_BAM}"
echo "  Total alignments (including flagged duplicates): ${N_READS}"

# ============================================================================
# Next steps:
#   - 05_haplotypecaller_gvcf.slurm (variant discovery, uses ${DEDUP_BAM})
#   - 07_wasp_intersecting_snps.slurm (WASP correction, also uses ${DEDUP_BAM})
# ============================================================================
