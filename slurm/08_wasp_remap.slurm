#!/bin/bash
#SBATCH --job-name=wasp_remap
#SBATCH --output=slurm-logs/wasp_remap_%A_%a.out
#SBATCH --error=slurm-logs/wasp_remap_%A_%a.err
#SBATCH --array=1-6
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
# ============================================================================
# WASP step 2: remap allele-swapped reads with STAR
# ============================================================================
# For each of six samples, realigns the allele-swapped FASTQ files produced
# by 07_wasp_intersecting_snps.slurm using STAR with the same reference and
# splice-junction database used for the original second-pass alignment.
# This is the critical step of WASP's mapping-bias correction: reads with
# their SNP allele swapped are remapped, and any read whose remapped
# alignment differs from the original position is flagged as mapping-bias-
# sensitive and discarded in step 09.
#
# Using STAR (not BWA) matches the original alignment method, which is
# required for WASP's assumption that mapping differences reflect allele
# bias rather than aligner differences.
#
# Reference:
#   van de Geijn B, McVicker G, Gilad Y, Pritchard JK. WASP: allele-specific
#   software for robust molecular quantitative trait locus discovery.
#   Nat Methods. 2015;12:1061-1063.
#
# Sample manifest: config/samples.tsv (columns: sample_id, tissue, conceptus)
# ============================================================================

set -euo pipefail

# ---- Module environment ----------------------------------------------------
module load star/2.7.10a
module load samtools/1.17

# ---- Configuration ---------------------------------------------------------
THREADS=${SLURM_CPUS_PER_TASK:-8}
PROJECT_ROOT=<!-- FILL IN: /path/to/project/root -->
WORKDIR=${PROJECT_ROOT}/analysis
STAR_INDEX=${PROJECT_ROOT}/reference/star_index_GRCh38_v38
MERGED_SJ=${WORKDIR}/star/all_SJ.out.tab
INTERSECT_DIR=${WORKDIR}/wasp/intersecting
OUTDIR_BASE=${WORKDIR}/wasp/remap
SAMPLES_TSV=${PROJECT_ROOT}/config/samples.tsv

# Node-local scratch for intermediate STAR output (large; deleted at end)
NODE_TMP_ROOT=${TMPDIR:-/tmp}

mkdir -p "${OUTDIR_BASE}" slurm-logs

# ---- Resolve sample for this array task -----------------------------------
SAMPLE=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==i {print $1}' \
         "${SAMPLES_TSV}")

if [[ -z "${SAMPLE}" ]]; then
    echo "[ERROR] Could not resolve sample for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

echo "[$(date)] Starting WASP remap (STAR) for sample: ${SAMPLE}"

# ---- Input and output paths -----------------------------------------------
SAMPLE_INTERSECT_DIR=${INTERSECT_DIR}/${SAMPLE}
SAMPLE_OUTDIR=${OUTDIR_BASE}/${SAMPLE}
NODE_TMP=${NODE_TMP_ROOT}/wasp_remap_${SAMPLE}_${SLURM_JOB_ID}

mkdir -p "${SAMPLE_OUTDIR}" "${NODE_TMP}"

# ---- Preflight: locate WASP remap FASTQ files -----------------------------
# WASP find_intersecting_snps.py names remap FASTQs variably across versions.
# Support common patterns.
R1_FASTQ=$(ls ${SAMPLE_INTERSECT_DIR}/*remap.fq1.gz \
              ${SAMPLE_INTERSECT_DIR}/*remap_R1.fq.gz \
              ${SAMPLE_INTERSECT_DIR}/*_R1.fastq.gz \
              ${SAMPLE_INTERSECT_DIR}/*R1.fq.gz 2>/dev/null | head -n1 || true)
R2_FASTQ=$(ls ${SAMPLE_INTERSECT_DIR}/*remap.fq2.gz \
              ${SAMPLE_INTERSECT_DIR}/*remap_R2.fq.gz \
              ${SAMPLE_INTERSECT_DIR}/*_R2.fastq.gz \
              ${SAMPLE_INTERSECT_DIR}/*R2.fq.gz 2>/dev/null | head -n1 || true)

if [[ -z "${R1_FASTQ}" ]] || [[ ! -f "${R1_FASTQ}" ]]; then
    echo "[ERROR] Remap R1 FASTQ not found in ${SAMPLE_INTERSECT_DIR}" >&2
    echo "  Files present:" >&2
    ls -lh "${SAMPLE_INTERSECT_DIR}" >&2
    exit 2
fi

if [[ -z "${R2_FASTQ}" ]] || [[ ! -f "${R2_FASTQ}" ]]; then
    echo "[ERROR] Remap R2 FASTQ not found in ${SAMPLE_INTERSECT_DIR}" >&2
    ls -lh "${SAMPLE_INTERSECT_DIR}" >&2
    exit 3
fi

if [[ ! -f "${MERGED_SJ}" ]]; then
    echo "[ERROR] Merged splice junction file not found: ${MERGED_SJ}" >&2
    echo "  Run 02_merge_splice_junctions.sh before this script." >&2
    exit 4
fi

echo "  Remap R1: ${R1_FASTQ}"
echo "  Remap R2: ${R2_FASTQ}"
echo "  STAR index: ${STAR_INDEX}"
echo "  Merged splice junctions: ${MERGED_SJ}"
echo "  Node-local scratch: ${NODE_TMP}"

# ---- Run STAR remap -------------------------------------------------------
# STAR parameters match the original second-pass alignment
# (03_star_secondpass.slurm) so that mapping differences between the
# original and remapped alignments reflect allele-swap effects, not aligner
# parameter differences.
#
# --outSAMtype BAM Unsorted: write unsorted BAM (sort separately with
#   samtools to control memory usage; STAR's internal sort can exceed
#   node memory on large libraries)
echo "[$(date)] Running STAR remap"
STAR --runThreadN ${THREADS} \
     --genomeDir "${STAR_INDEX}" \
     --sjdbFileChrStartEnd "${MERGED_SJ}" \
     --readFilesIn "${R1_FASTQ}" "${R2_FASTQ}" \
     --readFilesCommand zcat \
     --outFileNamePrefix "${NODE_TMP}/remap." \
     --outSAMtype BAM Unsorted

# ---- Locate unsorted BAM (STAR naming varies by version) ------------------
UNSORTED_BAM=$(ls ${NODE_TMP}/remap.Aligned.out.bam \
                  ${NODE_TMP}/remap.Aligned.Unsorted.out.bam 2>/dev/null | head -n1 || true)

if [[ -z "${UNSORTED_BAM}" ]] || [[ ! -f "${UNSORTED_BAM}" ]]; then
    echo "[ERROR] STAR did not produce an unsorted BAM in ${NODE_TMP}" >&2
    ls -lh "${NODE_TMP}" >&2
    exit 5
fi
echo "  STAR unsorted BAM: ${UNSORTED_BAM}"

# ---- Sort with samtools ---------------------------------------------------
# Sorting is done separately from STAR to control memory and use node-local
# scratch for temporary sort files. Per-thread memory kept modest to avoid
# exceeding total allocation.
export TMPDIR="${NODE_TMP}"
PER_THREAD_MEM=4G

FINAL_BAM=${SAMPLE_OUTDIR}/remap.Aligned.sortedByCoord.out.bam
echo "[$(date)] Sorting with samtools (threads=${THREADS}, per-thread=${PER_THREAD_MEM})"
samtools sort -@ "${THREADS}" -m "${PER_THREAD_MEM}" \
              -o "${FINAL_BAM}" "${UNSORTED_BAM}"

# ---- Index final BAM ------------------------------------------------------
samtools index "${FINAL_BAM}"

# ---- Copy STAR logs to persistent output directory -----------------------
for f in remap.Log.final.out remap.Log.out remap.SJ.out.tab; do
    if [[ -f "${NODE_TMP}/${f}" ]]; then
        cp "${NODE_TMP}/${f}" "${SAMPLE_OUTDIR}/" || true
    fi
done

# ---- Clean up node-local scratch ------------------------------------------
rm -rf "${NODE_TMP}"

# ---- Verify final BAM and report summary ----------------------------------
if [[ ! -f "${FINAL_BAM}" ]]; then
    echo "[ERROR] Final remapped BAM not produced: ${FINAL_BAM}" >&2
    exit 6
fi

N_MAPPED=$(samtools view -c -F 4 "${FINAL_BAM}")
echo "[$(date)] Done ${SAMPLE}"
echo "  Remapped BAM: ${FINAL_BAM}"
echo "  Mapped reads in remap: ${N_MAPPED}"

# ============================================================================
# Next step: run 09_wasp_filter_merge.slurm to compare original and remapped
# alignments with filter_remapped_reads.py, then merge with the keep BAM to
# produce the final WASP-filtered BAM.
# ============================================================================
