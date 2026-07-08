#!/bin/bash
# ============================================================================
# Merge per-sample splice junctions from STAR first pass into a single file
# ============================================================================
# STAR's first-pass alignment (01_star_firstpass.slurm) produces one
# SJ.out.tab per sample. This script concatenates them, filters out unreliable
# junctions (novel junctions supported by fewer than 3 uniquely mapped reads,
# junctions on chrM, and non-canonical strand annotations), and writes a
# single merged file used by the STAR second pass (03_star_secondpass.slurm).
#
# STAR SJ.out.tab format (tab-separated, no header):
#   1: chromosome
#   2: first intronic base (1-based)
#   3: last intronic base (1-based)
#   4: strand (0=undefined, 1=+, 2=-)
#   5: intron motif (0=non-canonical, 1=GT/AG, 2=CT/AC, 3=GC/AG,
#                    4=CT/GC, 5=AT/AC, 6=GT/AT)
#   6: annotated (0=novel, 1=annotated)
#   7: number of uniquely mapping reads crossing the junction
#   8: number of multi-mapping reads crossing the junction
#   9: maximum spliced alignment overhang
#
# This is not a SLURM job; run on a login node or as part of an interactive
# session. Takes <1 minute for six samples.
# ============================================================================

set -euo pipefail

# ---- Configuration ---------------------------------------------------------
PROJECT_ROOT=<!-- FILL IN: /path/to/project/root -->
WORKDIR=${PROJECT_ROOT}/analysis/star
FIRSTPASS_DIR=${WORKDIR}/firstpass
OUTPUT=${WORKDIR}/all_SJ.out.tab

# ---- Sanity checks ---------------------------------------------------------
if [[ ! -d "${FIRSTPASS_DIR}" ]]; then
    echo "[ERROR] First-pass output directory not found: ${FIRSTPASS_DIR}" >&2
    exit 1
fi

SJ_FILES=(${FIRSTPASS_DIR}/*.SJ.out.tab)

if [[ ${#SJ_FILES[@]} -eq 0 ]] || [[ ! -f "${SJ_FILES[0]}" ]]; then
    echo "[ERROR] No SJ.out.tab files found in ${FIRSTPASS_DIR}" >&2
    exit 2
fi

echo "[$(date)] Merging splice junctions from ${#SJ_FILES[@]} samples"
for f in "${SJ_FILES[@]}"; do
    n_junctions=$(wc -l < "${f}")
    echo "  ${f}: ${n_junctions} junctions"
done

# ---- Merge and filter ------------------------------------------------------
# Filtering rules:
#   - Retain all annotated junctions (column 6 == 1) regardless of read support
#   - Retain novel junctions (column 6 == 0) with >=3 uniquely mapping reads
#     (column 7 >= 3)
#   - Exclude junctions on chrM (mitochondrial genome)
#   - Exclude junctions with undefined strand (column 4 == 0) unless annotated
#   - Deduplicate identical junction coordinates across samples
#
# These filters follow STAR two-pass best practice (Dobin lab recommendations)
# and reduce the risk of introducing spurious novel junctions into the
# second-pass alignment.
cat "${SJ_FILES[@]}" \
    | awk 'BEGIN{FS=OFS="\t"} 
           $1 != "chrM" && \
           ( $6 == 1 || ($6 == 0 && $7 >= 3) ) && \
           ( $4 != 0 || $6 == 1 )' \
    | sort -k1,1 -k2,2n -k3,3n -u \
    > "${OUTPUT}"

N_MERGED=$(wc -l < "${OUTPUT}")
echo "[$(date)] Merged file written: ${OUTPUT}"
echo "  Total unique junctions after filtering: ${N_MERGED}"

# ---- Summary --------------------------------------------------------------
N_ANNOTATED=$(awk '$6 == 1' "${OUTPUT}" | wc -l)
N_NOVEL=$(awk '$6 == 0' "${OUTPUT}" | wc -l)
echo "  Annotated junctions: ${N_ANNOTATED}"
echo "  Novel junctions (>=3 uniquely mapping reads): ${N_NOVEL}"

# ============================================================================
# Next step: run 03_star_secondpass.slurm using ${OUTPUT} as the
# --sjdbFileChrStartEnd input.
# ============================================================================
