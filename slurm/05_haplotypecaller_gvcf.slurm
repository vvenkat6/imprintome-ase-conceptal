#!/bin/bash
#SBATCH --job-name=hc_gvcf
#SBATCH --output=slurm-logs/hc_gvcf_%A_%a.out
#SBATCH --error=slurm-logs/hc_gvcf_%A_%a.err
#SBATCH --array=1-6
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
# ============================================================================
# Per-sample variant calling: SplitNCigarReads + HaplotypeCaller (GVCF mode)
# ============================================================================
# For each of six samples, prepares RNA-seq alignments for variant calling
# and produces a per-sample GVCF. Two GATK steps:
#
#   1. SplitNCigarReads: RNA-seq alignments contain N in the CIGAR string
#      where reads span splice junctions. HaplotypeCaller cannot process
#      these directly. SplitNCigarReads splits each spliced read into
#      per-exon segments, reassigns mapping quality (255 -> 60), and hard-
#      clips overhangs into intronic regions.
#
#   2. HaplotypeCaller (in GVCF mode with -ERC GVCF): performs local
#      de novo assembly and emits a genomic VCF that includes reference
#      confidence blocks. GVCFs from all samples are subsequently combined
#      and jointly genotyped (06_joint_genotype.slurm).
#
# References:
#   McKenna A, et al. Genome Res. 2010;20:1297-1303.
#   DePristo MA, et al. Nat Genet. 2011;43:491-498.
#   Van der Auwera GA, et al. Curr Protoc Bioinformatics. 2013;43:11.10.1.
#
# Sample manifest: config/samples.tsv (columns: sample_id, tissue, conceptus)
# ============================================================================

set -euo pipefail

# ---- Module environment ----------------------------------------------------
module load gatk/4.4.0.0
module load samtools/1.17

# ---- Configuration ---------------------------------------------------------
THREADS=${SLURM_CPUS_PER_TASK:-8}
PROJECT_ROOT=<!-- FILL IN: /path/to/project/root -->
WORKDIR=${PROJECT_ROOT}/analysis
DEDUP_DIR=${WORKDIR}/dedup
OUTDIR=${WORKDIR}/gvcf
REF=${PROJECT_ROOT}/reference/GRCh38.primary_assembly.genome.fa
SAMPLES_TSV=${PROJECT_ROOT}/config/samples.tsv

# Scratch directory for intermediate split BAMs (large; delete after use)
TMPDIR_ROOT=${WORKDIR}/tmp
mkdir -p "${OUTDIR}" "${TMPDIR_ROOT}" slurm-logs

# ---- Resolve sample for this array task -----------------------------------
SAMPLE=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==i {print $1}' \
         "${SAMPLES_TSV}")

if [[ -z "${SAMPLE}" ]]; then
    echo "[ERROR] Could not resolve sample for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

echo "[$(date)] Starting per-sample variant calling for sample: ${SAMPLE}"

# ---- Input and output paths -----------------------------------------------
INBAM=${DEDUP_DIR}/${SAMPLE}.dedup.bam
TMP_SPLIT=${TMPDIR_ROOT}/${SAMPLE}.split.bam
OUTGVCF=${OUTDIR}/${SAMPLE}.g.vcf.gz

if [[ ! -f "${INBAM}" ]]; then
    echo "[ERROR] Input deduplicated BAM not found: ${INBAM}" >&2
    echo "  Run 04_readgroups_dedup.slurm before this script." >&2
    exit 2
fi

if [[ ! -f "${REF}" ]]; then
    echo "[ERROR] Reference FASTA not found: ${REF}" >&2
    exit 3
fi

# ---- Ensure BAM index is present ------------------------------------------
if [[ ! -f "${INBAM}.bai" ]] && [[ ! -f "${INBAM%.bam}.bai" ]]; then
    echo "[INFO] BAM index not found; creating with samtools"
    samtools index "${INBAM}"
fi

# ---- Ensure reference FASTA has .fai and .dict indices --------------------
# GATK requires both a samtools .fai index and a Picard-style .dict file.
if [[ ! -f "${REF}.fai" ]]; then
    echo "[INFO] Reference .fai index not found; creating"
    samtools faidx "${REF}"
fi
if [[ ! -f "${REF%.fa}.dict" ]] && [[ ! -f "${REF%.fasta}.dict" ]]; then
    echo "[INFO] Reference .dict not found; creating with GATK"
    gatk CreateSequenceDictionary -R "${REF}"
fi

echo "  Input BAM: ${INBAM}"
echo "  Reference: ${REF}"
echo "  Output GVCF: ${OUTGVCF}"
echo "  Threads: ${THREADS}"

# ---- Step 1: SplitNCigarReads --------------------------------------------
# Splits reads with N in CIGAR into separate per-exon reads. Required for
# RNA-seq variant calling because HaplotypeCaller does not process
# N-containing CIGAR strings correctly.
echo "[$(date)] Step 1: SplitNCigarReads"
gatk --java-options "-Xmx40g" SplitNCigarReads \
     -R "${REF}" \
     -I "${INBAM}" \
     -O "${TMP_SPLIT}"

if [[ ! -f "${TMP_SPLIT}" ]]; then
    echo "[ERROR] SplitNCigarReads did not produce output BAM" >&2
    exit 4
fi

# ---- Step 2: HaplotypeCaller (GVCF mode) ----------------------------------
# -ERC GVCF: emit reference confidence blocks (required for joint genotyping
#            in step 06)
# --native-pair-hmm-threads: parallelize the pair-HMM step of local assembly
#
# No hard filters (QD, FS, MQ, ReadPosRankSum, VQSR) are applied here. All
# GATK default variant-quality behavior is retained. If additional filtering
# becomes necessary for a revision, apply it downstream on the joint VCF.
echo "[$(date)] Step 2: HaplotypeCaller (GVCF)"
gatk --java-options "-Xmx40g -XX:ParallelGCThreads=${THREADS}" HaplotypeCaller \
     -R "${REF}" \
     -I "${TMP_SPLIT}" \
     -O "${OUTGVCF}" \
     -ERC GVCF \
     --native-pair-hmm-threads ${THREADS}

# ---- Verify output GVCF ---------------------------------------------------
if [[ ! -f "${OUTGVCF}" ]]; then
    echo "[ERROR] HaplotypeCaller did not produce output GVCF: ${OUTGVCF}" >&2
    exit 5
fi

# ---- Clean up intermediate split BAM --------------------------------------
# The split BAM is large and not needed once the GVCF has been produced.
rm -f "${TMP_SPLIT}"
rm -f "${TMP_SPLIT}.bai" "${TMP_SPLIT%.bam}.bai"

# ---- Report GVCF summary --------------------------------------------------
N_VARIANTS=$(zcat "${OUTGVCF}" | grep -vc "^#" || true)
echo "[$(date)] Done ${SAMPLE}"
echo "  Output GVCF: ${OUTGVCF}"
echo "  Total GVCF records: ${N_VARIANTS}"

# ============================================================================
# Next step: run 06_joint_genotype.slurm to combine all six GVCFs and
# perform cohort joint genotyping.
# ============================================================================
