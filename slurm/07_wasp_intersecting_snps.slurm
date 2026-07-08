#!/bin/bash
#SBATCH --job-name=wasp_intersect
#SBATCH --output=slurm-logs/wasp_intersect_%A_%a.out
#SBATCH --error=slurm-logs/wasp_intersect_%A_%a.err
#SBATCH --array=1-6
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=08:00:00
# ============================================================================
# WASP step 1: identify reads overlapping heterozygous SNPs
# ============================================================================
# For each of six samples, uses WASP's find_intersecting_snps.py to separate
# reads that overlap heterozygous SNPs (candidates for allele-swap remapping)
# from reads that do not (retained directly for the final WASP-filtered BAM).
#
# Inputs:
#   - Deduplicated BAM (from 04_readgroups_dedup.slurm)
#   - Per-chromosome SNP tables generated from cohort.snps.vcf.gz
#     (produced by prepare_wasp_snp_tables.sh; see companion script)
#
# Outputs (written to sample-specific subdirectory):
#   ${SAMPLE}.dedup.keep.bam       Reads not overlapping any het SNP
#   ${SAMPLE}.dedup.to.remap.bam   Reads overlapping het SNPs
#   ${SAMPLE}.dedup.remap.fq.gz    Allele-swapped FASTQ for remapping
#
# Reference:
#   van de Geijn B, McVicker G, Gilad Y, Pritchard JK. WASP: allele-specific
#   software for robust molecular quantitative trait locus discovery.
#   Nat Methods. 2015;12:1061-1063.
#
# WASP repository: https://github.com/bmvdgeijn/WASP
# ============================================================================

set -euo pipefail

# ---- Module environment ----------------------------------------------------
module load python/3.9.7
module load samtools/1.17

# WASP is installed from source (not module-available). See environment/modules.txt.
WASP=<!-- FILL IN: /home3/vvenka23/bin/WASP -->

# ---- Configuration ---------------------------------------------------------
THREADS=${SLURM_CPUS_PER_TASK:-4}
PROJECT_ROOT=<!-- FILL IN: /path/to/project/root -->
WORKDIR=${PROJECT_ROOT}/analysis
DEDUP_DIR=${WORKDIR}/dedup
SNP_DIR=${WORKDIR}/wasp/snp_tables
OUTDIR_BASE=${WORKDIR}/wasp/intersecting
SAMPLES_TSV=${PROJECT_ROOT}/config/samples.tsv

mkdir -p "${OUTDIR_BASE}" slurm-logs

# ---- Resolve sample for this array task -----------------------------------
SAMPLE=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==i {print $1}' \
         "${SAMPLES_TSV}")

if [[ -z "${SAMPLE}" ]]; then
    echo "[ERROR] Could not resolve sample for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

echo "[$(date)] Starting WASP find_intersecting_snps for sample: ${SAMPLE}"

# ---- Input and output paths -----------------------------------------------
INBAM=${DEDUP_DIR}/${SAMPLE}.dedup.bam
OUTDIR=${OUTDIR_BASE}/${SAMPLE}

mkdir -p "${OUTDIR}"

# ---- Preflight checks -----------------------------------------------------
if [[ ! -f "${INBAM}" ]]; then
    echo "[ERROR] Input deduplicated BAM not found: ${INBAM}" >&2
    echo "  Run 04_readgroups_dedup.slurm before this script." >&2
    exit 2
fi

if [[ ! -d "${SNP_DIR}" ]]; then
    echo "[ERROR] WASP SNP table directory not found: ${SNP_DIR}" >&2
    echo "  Run prepare_wasp_snp_tables.sh (see below) before this script." >&2
    exit 3
fi

N_SNP_FILES=$(ls -1 ${SNP_DIR}/*.snps.txt.gz 2>/dev/null | wc -l)
if [[ ${N_SNP_FILES} -eq 0 ]]; then
    echo "[ERROR] No per-chromosome SNP tables found in ${SNP_DIR}" >&2
    exit 4
fi
echo "  Found ${N_SNP_FILES} per-chromosome SNP table(s)"

if [[ ! -x "${WASP}/mapping/find_intersecting_snps.py" ]] && \
   [[ ! -f "${WASP}/mapping/find_intersecting_snps.py" ]]; then
    echo "[ERROR] WASP script not found: ${WASP}/mapping/find_intersecting_snps.py" >&2
    exit 5
fi

echo "  Input BAM: ${INBAM}"
echo "  SNP tables: ${SNP_DIR}"
echo "  Output directory: ${OUTDIR}"

# ---- Ensure BAM index present ---------------------------------------------
if [[ ! -f "${INBAM}.bai" ]] && [[ ! -f "${INBAM%.bam}.bai" ]]; then
    echo "[INFO] BAM index not found; creating with samtools"
    samtools index "${INBAM}"
fi

# ---- Run find_intersecting_snps.py ----------------------------------------
# --is_paired_end: handle read pairs (both mates evaluated together)
# --is_sorted: BAM is coordinate-sorted (from STAR + MarkDuplicates)
# --output_dir: where to write .keep.bam, .to.remap.bam, .remap.fq.gz
# --snp_dir: directory containing per-chromosome SNP tables
echo "[$(date)] Running find_intersecting_snps.py"
python3 "${WASP}/mapping/find_intersecting_snps.py" \
        --is_paired_end \
        --is_sorted \
        --output_dir "${OUTDIR}" \
        --snp_dir "${SNP_DIR}" \
        "${INBAM}"

# ---- Verify expected outputs ----------------------------------------------
KEEP_BAM=$(ls ${OUTDIR}/*keep.bam 2>/dev/null | head -n1 || true)
TO_REMAP_BAM=$(ls ${OUTDIR}/*to.remap.bam 2>/dev/null | head -n1 || true)

if [[ -z "${KEEP_BAM}" ]] || [[ ! -f "${KEEP_BAM}" ]]; then
    echo "[ERROR] keep BAM not produced in ${OUTDIR}" >&2
    ls -lh "${OUTDIR}" >&2
    exit 6
fi

if [[ -z "${TO_REMAP_BAM}" ]] || [[ ! -f "${TO_REMAP_BAM}" ]]; then
    echo "[ERROR] to.remap BAM not produced in ${OUTDIR}" >&2
    ls -lh "${OUTDIR}" >&2
    exit 7
fi

# ---- Report summary -------------------------------------------------------
N_KEEP=$(samtools view -c "${KEEP_BAM}")
N_TO_REMAP=$(samtools view -c "${TO_REMAP_BAM}")
N_TOTAL=$((N_KEEP + N_TO_REMAP))

echo "[$(date)] Done ${SAMPLE}"
echo "  Keep BAM (no het SNP overlap): ${KEEP_BAM}"
echo "    Reads: ${N_KEEP}"
echo "  To-remap BAM (het SNP overlap): ${TO_REMAP_BAM}"
echo "    Reads: ${N_TO_REMAP}"
echo "  Total reads processed: ${N_TOTAL}"
if [[ ${N_TOTAL} -gt 0 ]]; then
    PCT=$(awk -v r=${N_TO_REMAP} -v t=${N_TOTAL} 'BEGIN {printf "%.2f", 100*r/t}')
    echo "  Fraction to remap: ${PCT}%"
fi

# ============================================================================
# Next step: run 08_wasp_remap.slurm to realign the to-remap reads with STAR
# after allele swap.
# ============================================================================
