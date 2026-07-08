#!/bin/bash
#SBATCH --job-name=star_pass1
#SBATCH --output=slurm-logs/star_pass1_%A_%a.out
#SBATCH --error=slurm-logs/star_pass1_%A_%a.err
#SBATCH --array=1-6
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=08:00:00
# ============================================================================
# STAR first-pass alignment for splice-junction discovery
# ============================================================================
# For each of six samples, aligns paired-end 150 bp reads to GRCh38 using
# STAR in first-pass mode. Only splice junctions (SJ.out.tab) are retained;
# no BAM output is written at this stage. Junctions from all six samples are
# subsequently merged (02_merge_splice_junctions.sh) into a single
# `all_SJ.out.tab` file that feeds the STAR second pass.
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
SAMPLES_TSV=${PROJECT_ROOT}/config/samples.tsv

mkdir -p ${WORKDIR}/firstpass slurm-logs

# ---- Resolve sample for this array task -----------------------------------
# Read the SLURM_ARRAY_TASK_ID-th data row from the sample manifest.
# samples.tsv has a header row, so we add 1 to skip it.
SAMPLE=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==i {print $1}' \
         "${SAMPLES_TSV}")

if [[ -z "${SAMPLE}" ]]; then
    echo "[ERROR] Could not resolve sample for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

echo "[$(date)] Starting STAR first-pass for sample: ${SAMPLE}"

# ---- Input FASTQ files -----------------------------------------------------
R1=${FASTQ_DIR}/${SAMPLE}_L001_R1_001.fastq.gz
R2=${FASTQ_DIR}/${SAMPLE}_L001_R2_001.fastq.gz
OUTPREFIX=${WORKDIR}/firstpass/${SAMPLE}.

if [[ ! -f "${R1}" ]] || [[ ! -f "${R2}" ]]; then
    echo "[ERROR] FASTQ pair not found for ${SAMPLE}" >&2
    echo "  Expected R1: ${R1}" >&2
    echo "  Expected R2: ${R2}" >&2
    exit 2
fi

echo "  R1: ${R1}"
echo "  R2: ${R2}"
echo "  STAR index: ${STAR_INDEX}"
echo "  Output prefix: ${OUTPREFIX}"
echo "  Threads: ${THREADS}"

# ---- Run STAR first pass ---------------------------------------------------
# --outSAMtype None: no BAM output; only splice junctions are needed here
# --outSJfilterReads All: include all reads when reporting splice junctions
STAR --runThreadN ${THREADS} \
     --genomeDir ${STAR_INDEX} \
     --readFilesIn ${R1} ${R2} \
     --readFilesCommand zcat \
     --outFileNamePrefix ${OUTPREFIX} \
     --outSAMtype None \
     --outSJfilterReads All

echo "[$(date)] Finished STAR first-pass for ${SAMPLE}"
echo "  Splice junctions written to: ${OUTPREFIX}SJ.out.tab"

# ============================================================================
# Next step: run 02_merge_splice_junctions.sh to combine all six samples'
# SJ.out.tab files into a single all_SJ.out.tab for the second-pass alignment.
# ============================================================================

# ============================================================================
# STAR index construction (run once before this script)
# ============================================================================
# The STAR index at ${STAR_INDEX} was built with:
#
#   STAR --runMode genomeGenerate \
#        --runThreadN 16 \
#        --genomeDir ${STAR_INDEX} \
#        --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
#        --sjdbGTFfile gencode.v38.annotation.gtf \
#        --sjdbOverhang 149
#
# Reference sources:
#   GRCh38 primary assembly FASTA: GENCODE release 38
#   Annotation GTF: gencode.v38.annotation.gtf
#     https://www.gencodegenes.org/human/release_38.html
# ============================================================================
