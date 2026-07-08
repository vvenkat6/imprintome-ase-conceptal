#!/bin/bash
#SBATCH --job-name=star_pass2
#SBATCH --output=slurm-logs/star_pass2_%A_%a.out
#SBATCH --error=slurm-logs/star_pass2_%A_%a.err
#SBATCH --array=1-6
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=12:00:00
# ============================================================================
# STAR second-pass alignment using merged splice junctions
# ============================================================================
# For each of six samples, aligns paired-end 150 bp reads to GRCh38 using
# STAR with the merged splice-junction database from all six first-pass
# alignments (see 02_merge_splice_junctions.sh). This pass produces
# coordinate-sorted BAM files that feed downstream deduplication, variant
# calling, and WASP correction.
#
# Reference:
#   Dobin A, et al. STAR: ultrafast universal RNA-seq aligner.
#   Bioinformatics. 2013;29:15-21.
#
# Sample manifest: config/samples.tsv (columns: sample_id, tissue, conceptus)
# ============================================================================

set -euo pipefail

# ---- Module environment ----------------------------------------------------
module load star/2.7.10a

# ---- Configuration ---------------------------------------------------------
THREADS=${SLURM_CPUS_PER_TASK:-8}
PROJECT_ROOT=<!-- FILL IN: /path/to/project/root -->
FASTQ_DIR=${PROJECT_ROOT}/fastq
STAR_INDEX=${PROJECT_ROOT}/reference/star_index_GRCh38_v38
WORKDIR=${PROJECT_ROOT}/analysis/star
MERGED_SJ=${WORKDIR}/all_SJ.out.tab
SAMPLES_TSV=${PROJECT_ROOT}/config/samples.tsv

mkdir -p ${WORKDIR}/secondpass slurm-logs

# ---- Sanity checks ---------------------------------------------------------
if [[ ! -f "${MERGED_SJ}" ]]; then
    echo "[ERROR] Merged splice junction file not found: ${MERGED_SJ}" >&2
    echo "  Run 02_merge_splice_junctions.sh before this script." >&2
    exit 1
fi

# ---- Resolve sample for this array task -----------------------------------
SAMPLE=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==i {print $1}' \
         "${SAMPLES_TSV}")

if [[ -z "${SAMPLE}" ]]; then
    echo "[ERROR] Could not resolve sample for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 2
fi

echo "[$(date)] Starting STAR second-pass for sample: ${SAMPLE}"

# ---- Input FASTQ files -----------------------------------------------------
R1=${FASTQ_DIR}/${SAMPLE}_L001_R1_001.fastq.gz
R2=${FASTQ_DIR}/${SAMPLE}_L001_R2_001.fastq.gz
OUTPREFIX=${WORKDIR}/secondpass/${SAMPLE}.

if [[ ! -f "${R1}" ]] || [[ ! -f "${R2}" ]]; then
    echo "[ERROR] FASTQ pair not found for ${SAMPLE}" >&2
    echo "  Expected R1: ${R1}" >&2
    echo "  Expected R2: ${R2}" >&2
    exit 3
fi

echo "  R1: ${R1}"
echo "  R2: ${R2}"
echo "  STAR index: ${STAR_INDEX}"
echo "  Merged splice junctions: ${MERGED_SJ}"
echo "  Output prefix: ${OUTPREFIX}"
echo "  Threads: ${THREADS}"

# ---- Run STAR second pass -------------------------------------------------
# --sjdbFileChrStartEnd: use merged splice junctions from all six samples'
#                        first-pass alignments as the junction database
# --outSAMtype BAM SortedByCoordinate: write coordinate-sorted BAM
# --quantMode TranscriptomeSAM GeneCounts: also output transcriptome-aligned
#                        BAM and per-gene read counts (STAR built-in
#                        quantification; not used downstream but retained
#                        for potential future analysis)
STAR --runThreadN ${THREADS} \
     --genomeDir ${STAR_INDEX} \
     --sjdbFileChrStartEnd ${MERGED_SJ} \
     --readFilesIn ${R1} ${R2} \
     --readFilesCommand zcat \
     --outFileNamePrefix ${OUTPREFIX} \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM GeneCounts

# ---- Verify output BAM ----------------------------------------------------
FINAL_BAM=${OUTPREFIX}Aligned.sortedByCoord.out.bam
if [[ ! -f "${FINAL_BAM}" ]]; then
    echo "[ERROR] Expected output BAM not created: ${FINAL_BAM}" >&2
    exit 4
fi

N_READS=$(module load samtools/1.17 && samtools view -c "${FINAL_BAM}")
echo "[$(date)] Finished STAR second-pass for ${SAMPLE}"
echo "  Output BAM: ${FINAL_BAM}"
echo "  Total alignments: ${N_READS}"

# ============================================================================
# Next step: run 04_readgroups_dedup.slurm to add read groups and mark
# duplicates on ${FINAL_BAM}.
# ============================================================================
