#!/bin/bash
#SBATCH --job-name=gene_audit
#SBATCH --output=slurm-logs/gene_audit_%j.out
#SBATCH --error=slurm-logs/gene_audit_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
# ============================================================================
# Gene-body read count and heterozygous SNP audit
# ============================================================================
# Quantifies gene-body read coverage and heterozygous SNP counts across all
# six WASP-filtered BAMs for a curated set of established imprinted genes.
# Two outputs support the manuscript:
#
#   1. Supplementary Table S2: per-sample read counts for the 13 benchmark
#      imprinted genes that yielded ASE-informative reads in only one tissue
#      (drives the depth-limited vs SNP-coverage-limited classification in
#      Figure 2 and the Discussion).
#
#   2. Supplementary Table S3: extended audit of established human imprinted
#      genes from geneimprint.com not recovered in the ASE table, classifying
#      each into: not expressed / no het SNPs / nested (multi-gene excluded)
#      / below aggregation thresholds / not annotated.
#
# This is not an array job; the audit operates across all six samples for
# each gene in the input list.
#
# Input: config/benchmark_gene_list.txt (one gene symbol per line)
# ============================================================================

set -euo pipefail

# ---- Module environment ----------------------------------------------------
module load samtools/1.17
module load bcftools/1.17

# ---- Configuration ---------------------------------------------------------
PROJECT_ROOT=<!-- FILL IN: /path/to/project/root -->
WORKDIR=${PROJECT_ROOT}/analysis
WASP_DIR=${WORKDIR}/wasp/final
VCF=${WORKDIR}/variants/cohort.snps.vcf.gz
GTF=${PROJECT_ROOT}/reference/gencode.v38.annotation.gtf.gz
GENE_LIST=${PROJECT_ROOT}/config/benchmark_gene_list.txt
OUTDIR=${WORKDIR}/audit

mkdir -p "${OUTDIR}" slurm-logs

# ---- Sample manifest: WASP-filtered BAMs in the exact order used in the
# ---- manuscript column ordering ------------------------------------------
# Order preserved for reproducibility with Supplementary Table S2 as shown
# in the paper (sample columns: 21423L, 21424B, 21424L, 21562L, 22288B, 22316B)
BAMS=(${WASP_DIR}/21423L_S1.wasp.sorted.bam
      ${WASP_DIR}/21424B_S2.wasp.sorted.bam
      ${WASP_DIR}/21424L_S3.wasp.sorted.bam
      ${WASP_DIR}/21562L_S4.wasp.sorted.bam
      ${WASP_DIR}/22288B_S5.wasp.sorted.bam
      ${WASP_DIR}/22316B_S6.wasp.sorted.bam)

# ---- Preflight checks -----------------------------------------------------
if [[ ! -f "${GTF}" ]]; then
    echo "[ERROR] GENCODE GTF not found: ${GTF}" >&2
    exit 1
fi

if [[ ! -f "${GENE_LIST}" ]]; then
    echo "[ERROR] Gene list not found: ${GENE_LIST}" >&2
    exit 2
fi

if [[ ! -f "${VCF}" ]] || [[ ! -f "${VCF}.tbi" ]]; then
    echo "[ERROR] Cohort SNP VCF or its tabix index not found: ${VCF}" >&2
    exit 3
fi

for bam in "${BAMS[@]}"; do
    if [[ ! -f "${bam}" ]]; then
        echo "[ERROR] WASP-filtered BAM not found: ${bam}" >&2
        exit 4
    fi
done

echo "[$(date)] Starting gene-body audit"
echo "  Gene list: ${GENE_LIST}"
echo "  BAMs: ${#BAMS[@]} samples"
echo "  Output directory: ${OUTDIR}"

# ---- Step 1: extract gene coordinates from GENCODE GTF --------------------
# Match gene symbols against the GTF's gene_name field. Requires gawk for
# the 3-argument match() function; falls back to a sed-based extractor if
# gawk is unavailable.
GENE_COORDS=${OUTDIR}/gene_coords.tsv

echo "[$(date)] Step 1: extract gene coordinates from GTF"
if command -v gawk >/dev/null 2>&1; then
    zcat "${GTF}" \
      | awk '$3 == "gene"' \
      | gawk -F'\t' -v OFS='\t' '{
          match($9, /gene_name "([^"]+)"/, arr);
          print arr[1], $1 ":" $4 "-" $5
        }' \
      | grep -F -w -f "${GENE_LIST}" \
      | sort -u -k1,1 \
      > "${GENE_COORDS}"
else
    zcat "${GTF}" \
      | awk '$3 == "gene"' \
      | awk -F'\t' -v OFS='\t' '{ print $9, $1 ":" $4 "-" $5 }' \
      | sed -E 's/.*gene_name "([^"]+)".*\t/\1\t/' \
      | grep -F -w -f "${GENE_LIST}" \
      | sort -u -k1,1 \
      > "${GENE_COORDS}"
fi

N_FOUND=$(wc -l < "${GENE_COORDS}")
N_REQUESTED=$(sort -u "${GENE_LIST}" | wc -l)
echo "  Requested genes: ${N_REQUESTED}"
echo "  Found in GTF: ${N_FOUND}"

# Report genes not matched
NOT_FOUND=${OUTDIR}/genes_not_in_gtf.txt
comm -23 <(sort -u "${GENE_LIST}") <(cut -f1 "${GENE_COORDS}" | sort -u) > "${NOT_FOUND}"
N_NOT_FOUND=$(wc -l < "${NOT_FOUND}")
if [[ ${N_NOT_FOUND} -gt 0 ]]; then
    echo "  Genes not found in GTF: ${N_NOT_FOUND} (see ${NOT_FOUND})"
fi

# ---- Step 2: per-gene per-sample read counts and het SNP counts ----------
AUDIT_TSV=${OUTDIR}/gene_body_audit.tsv

echo "[$(date)] Step 2: gene-body read counts and het SNP counts"

# Build header row
{
    printf "gene\tregion"
    for bam in "${BAMS[@]}"; do
        printf "\t%s" "$(basename ${bam} .wasp.sorted.bam)"
    done
    printf "\tn_het_snps\n"

    # Iterate over each gene
    while IFS=$'\t' read -r GENE REGION; do
        printf "%s\t%s" "${GENE}" "${REGION}"

        # Read counts per sample
        for bam in "${BAMS[@]}"; do
            count=$(samtools view -c -F 260 "${bam}" "${REGION}")
            printf "\t%s" "${count}"
        done

        # Heterozygous SNP count from cohort VCF for this region
        # GT[*]="het" matches any sample being heterozygous at that site
        n_het=$(bcftools view -r "${REGION}" "${VCF}" 2>/dev/null \
                | bcftools view -i 'GT[*]="het"' 2>/dev/null \
                | grep -vc '^#' || echo 0)
        printf "\t%s\n" "${n_het:-0}"

    done < "${GENE_COORDS}"
} > "${AUDIT_TSV}"

# ---- Verify output --------------------------------------------------------
if [[ ! -s "${AUDIT_TSV}" ]]; then
    echo "[ERROR] Audit table is empty: ${AUDIT_TSV}" >&2
    exit 5
fi

N_AUDITED=$(($(wc -l < "${AUDIT_TSV}") - 1))
echo "[$(date)] Done"
echo "  Audit table: ${AUDIT_TSV}"
echo "  Genes audited: ${N_AUDITED}"

# ============================================================================
# Next step: run scripts/benchmark_audit.R to classify each gene into:
#   - Recovered in Table 1 (present in ASE gene table)
#   - Nested (lost to multi-gene SNP exclusion)
#   - Expressed with het SNPs but not in ASE table (needs pipeline investigation)
#   - Expressed with no het SNPs in cohort
#   - Expressed with 1-2 het SNPs (below aggregation thresholds)
#   - Not expressed in these tissues
#
# Downstream outputs:
#   - Supplementary Table S2: tissue-restricted benchmark subset
#   - Supplementary Table S3: extended imprinted-gene audit
#   - Figure 2 (scripts/figure2_readcount_plot.R)
# ============================================================================
