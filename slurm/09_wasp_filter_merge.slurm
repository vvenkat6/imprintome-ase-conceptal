#!/bin/bash
#SBATCH --job-name=wasp_filter
#SBATCH --output=slurm-logs/wasp_filter_%A_%a.out
#SBATCH --error=slurm-logs/wasp_filter_%A_%a.err
#SBATCH --array=1-6
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
# ============================================================================
# WASP step 3: filter remapped reads and merge to produce final WASP BAM
# ============================================================================
# For each of six samples, compares the original (to-remap) and remapped
# alignments produced by 07 and 08 using WASP's filter_remapped_reads.py.
# Reads whose remapped alignment does not return to the same genomic
# coordinates as the original are discarded (mapping-bias-sensitive).
# Reads that remap to identical coordinates are retained.
#
# The final WASP-filtered BAM is produced by merging:
#   1. Keep BAM: reads that never overlapped any het SNP (from step 07)
#   2. Retained remapped BAM: reads that survived the coordinate-return
#      filter (from filter_remapped_reads.py)
#
# The merged BAM is coordinate-sorted, indexed, and used as input to
# ASEReadCounter (step 10) and the gene-body audit (step 11).
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
module load python/3.9.7
module load samtools/1.17

WASP=<!-- FILL IN: /home3/vvenka23/bin/WASP -->

# ---- Configuration ---------------------------------------------------------
THREADS=${SLURM_CPUS_PER_TASK:-4}
PROJECT_ROOT=<!-- FILL IN: /path/to/project/root -->
WORKDIR=${PROJECT_ROOT}/analysis
INTERSECT_DIR=${WORKDIR}/wasp/intersecting
REMAP_DIR=${WORKDIR}/wasp/remap
OUTDIR=${WORKDIR}/wasp/final
SAMPLES_TSV=${PROJECT_ROOT}/config/samples.tsv

mkdir -p "${OUTDIR}" slurm-logs

# ---- Resolve sample for this array task -----------------------------------
SAMPLE=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==i {print $1}' \
         "${SAMPLES_TSV}")

if [[ -z "${SAMPLE}" ]]; then
    echo "[ERROR] Could not resolve sample for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

echo "[$(date)] Starting WASP filter + merge for sample: ${SAMPLE}"

# ---- Input paths ----------------------------------------------------------
SAMPLE_INTERSECT_DIR=${INTERSECT_DIR}/${SAMPLE}
SAMPLE_REMAP_DIR=${REMAP_DIR}/${SAMPLE}

# WASP output file names vary; use robust ls detection
KEEP_BAM=$(ls ${SAMPLE_INTERSECT_DIR}/*keep.bam 2>/dev/null | head -n1 || true)
TO_REMAP_BAM=$(ls ${SAMPLE_INTERSECT_DIR}/*to.remap.bam 2>/dev/null | head -n1 || true)
REMAPPED_BAM=${SAMPLE_REMAP_DIR}/remap.Aligned.sortedByCoord.out.bam

# ---- Preflight checks -----------------------------------------------------
if [[ -z "${KEEP_BAM}" ]] || [[ ! -f "${KEEP_BAM}" ]]; then
    echo "[ERROR] Keep BAM not found in ${SAMPLE_INTERSECT_DIR}" >&2
    ls -lh "${SAMPLE_INTERSECT_DIR}" >&2
    exit 2
fi

if [[ -z "${TO_REMAP_BAM}" ]] || [[ ! -f "${TO_REMAP_BAM}" ]]; then
    echo "[ERROR] To-remap BAM not found in ${SAMPLE_INTERSECT_DIR}" >&2
    ls -lh "${SAMPLE_INTERSECT_DIR}" >&2
    exit 3
fi

if [[ ! -f "${REMAPPED_BAM}" ]]; then
    echo "[ERROR] Remapped BAM not found: ${REMAPPED_BAM}" >&2
    echo "  Run 08_wasp_remap.slurm before this script." >&2
    exit 4
fi

if [[ ! -f "${WASP}/mapping/filter_remapped_reads.py" ]]; then
    echo "[ERROR] WASP script not found: ${WASP}/mapping/filter_remapped_reads.py" >&2
    exit 5
fi

echo "  Keep BAM: ${KEEP_BAM}"
echo "  To-remap BAM: ${TO_REMAP_BAM}"
echo "  Remapped BAM: ${REMAPPED_BAM}"

# ---- Output paths ---------------------------------------------------------
REMAPPED_KEEP_BAM=${SAMPLE_REMAP_DIR}/${SAMPLE}.remap.keep.bam
MERGED_UNSORTED=${OUTDIR}/${SAMPLE}.wasp.unsorted.bam
FINAL_SORTED=${OUTDIR}/${SAMPLE}.wasp.sorted.bam

# ---- Step 1: filter_remapped_reads.py -------------------------------------
# Compares each read in the to-remap BAM against its remapped alignment.
# Reads whose remapped coordinates match the original are retained; reads
# that remap to a different position (mapping-bias-sensitive) are discarded.
#
# Positional arguments: to_remap_bam remap_bam keep_bam
echo "[$(date)] Step 1: filter_remapped_reads.py"
python3 "${WASP}/mapping/filter_remapped_reads.py" \
        "${TO_REMAP_BAM}" \
        "${REMAPPED_BAM}" \
        "${REMAPPED_KEEP_BAM}"

if [[ ! -f "${REMAPPED_KEEP_BAM}" ]]; then
    echo "[ERROR] filter_remapped_reads.py did not produce ${REMAPPED_KEEP_BAM}" >&2
    exit 6
fi

# ---- Step 2: merge keep BAM + retained remapped BAM -----------------------
# The two BAMs contain non-overlapping sets of reads:
#   Keep BAM: reads that never overlapped a het SNP (from step 07)
#   Remapped keep BAM: reads that overlapped a het SNP AND survived the
#                      coordinate-return filter
echo "[$(date)] Step 2: samtools merge"
samtools merge -f "${MERGED_UNSORTED}" "${KEEP_BAM}" "${REMAPPED_KEEP_BAM}"

# ---- Step 3: coordinate-sort and index the merged BAM ---------------------
echo "[$(date)] Step 3: samtools sort + index"
samtools sort -@ ${THREADS} -o "${FINAL_SORTED}" "${MERGED_UNSORTED}"
samtools index "${FINAL_SORTED}"

# ---- Clean up intermediate unsorted BAM -----------------------------------
rm -f "${MERGED_UNSORTED}"

# ---- Verify final BAM -----------------------------------------------------
if [[ ! -f "${FINAL_SORTED}" ]]; then
    echo "[ERROR] Final WASP-filtered BAM not produced: ${FINAL_SORTED}" >&2
    exit 7
fi

# ---- Report summary -------------------------------------------------------
N_KEEP=$(samtools view -c "${KEEP_BAM}")
N_TO_REMAP=$(samtools view -c "${TO_REMAP_BAM}")
N_REMAP_KEPT=$(samtools view -c "${REMAPPED_KEEP_BAM}")
N_FINAL=$(samtools view -c "${FINAL_SORTED}")

if [[ ${N_TO_REMAP} -gt 0 ]]; then
    PCT_RETAINED=$(awk -v k=${N_REMAP_KEPT} -v t=${N_TO_REMAP} \
                       'BEGIN {printf "%.2f", 100*k/t}')
else
    PCT_RETAINED="N/A"
fi

echo "[$(date)] Done ${SAMPLE}"
echo "  Reads never overlapping het SNPs (retained directly): ${N_KEEP}"
echo "  Reads overlapping het SNPs (submitted for remap): ${N_TO_REMAP}"
echo "  Reads surviving coordinate-return filter: ${N_REMAP_KEPT}"
echo "  Fraction of remapped reads retained: ${PCT_RETAINED}%"
echo "  Final WASP-filtered BAM: ${FINAL_SORTED}"
echo "  Total reads in final BAM: ${N_FINAL}"

# ============================================================================
# Next steps:
#   - 10_asereadcounter.slurm (allelic read counting from ${FINAL_SORTED})
#   - 11_gene_body_audit.slurm (gene-body coverage audit)
# ============================================================================
