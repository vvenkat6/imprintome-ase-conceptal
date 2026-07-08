# Tables

Manuscript tables and supplementary tables for [_A conservative allele-specific expression framework prioritizes candidate imprinted loci in human conceptal brain and liver using the human imprintome_](../README.md).

<!-- FILL IN: journal, year, DOI once accepted -->

## Contents

| File | Description | Manuscript reference |
|---|---|---|
| `table1_benchmarks.tsv` | Established human imprinted genes recovered with ASE-informative reads (n=30) | Table 1 |
| `table2_prioritized.tsv` | Prioritized novel autosomal candidates from ICR-anchored ranking (n=26) | Table 2 |
| `table3_conflicting_unknown.tsv` | Genes with ASE evidence and imprinting-relevant genomic context but unresolved human imprinting status (n=8) | Table 3 |
| `supp_table1_scoring_rubric.tsv` | Prioritization scoring rubric and priority tier assignment thresholds | Supplementary Table S1 |
| `supp_table2_tissue_restricted_readcounts.tsv` | Per-sample gene-body read counts for the 13 tissue-restricted benchmark imprinted genes | Supplementary Table S2 |
| `supp_table3_extended_audit.tsv` | Full audit of established human imprinted genes queried from the geneimprint.com catalog, with per-gene classification | Supplementary Table S3 |

## Column definitions

### `table1_benchmarks.tsv`

Established imprinted genes recovered in the ASE data, grouped by sample support tier (present in all 6 / 5 / 4 samples; brain-only; liver-only) and sorted by mean imbalance within tier.

| Column | Definition |
|---|---|
| `gene` | HGNC gene symbol |
| `chrom` | Chromosome (without `chr` prefix) |
| `n_samples` | Number of ASE-informative samples (out of 6) |
| `n_brain` | Number of ASE-informative brain samples (out of 3) |
| `n_liver` | Number of ASE-informative liver samples (out of 3) |
| `tissue_class` | `both_tissues`, `brain_only`, or `liver_only` |
| `mean_imbalance` | Mean gene-level allelic imbalance across informative samples, defined as `abs(ref - alt) / (ref + alt)` |
| `min_imbalance` | Minimum gene-level allelic imbalance across informative samples |
| `direction_consistency` | `max(n_ref_dominant, n_alt_dominant) / n_informative_direction` |
| `padj` | Benjamini-Hochberg-adjusted Fisher combined *P*-value |
| `tier1_overlap` | Direct gene-body overlap with a tier 1 (n=332) ICR from Jima et al. 2022 |
| `tier1_dist_kb` | Distance in kilobases to the nearest tier 1 ICR (0 = overlap) |
| `tier2_dist_kb` | Distance in kilobases to the nearest tier 2 ICR (0 = overlap) |
| `icrs_1mb` | Combined tier 1 + tier 2 ICR count within 1 Mb of the gene body |
| `paired_21424` | Counted-allele concordance between brain and liver of sample 21424: `same`, `different`, or `missing` |

### `table2_prioritized.tsv`

Novel autosomal candidate genes from the ICR-anchored prioritization framework, grouped into four tiers: high-priority novel candidates (n=6), moderate-priority follow-up (n=8), lower-priority follow-up (n=9), and non-autosomal confounders retained for illustration (n=3).

| Column | Definition |
|---|---|
| `rank` | Overall priority rank across all 26 candidates |
| `gene` | HGNC gene symbol |
| `chrom` | Chromosome |
| `priority_tier` | `high_priority_novel_candidate`, `moderate_priority_followup`, `lower_priority_followup`, or `non_autosomal_confounder` |
| `priority_score` | Total priority score from the additive rubric (Supplementary Table S1) |
| `n_samples`, `n_brain`, `n_liver`, `mean_imbalance`, `min_imbalance`, `direction_consistency`, `padj` | Same definitions as Table 1 |
| `tier1_dist_kb`, `tier2_dist_kb`, `icrs_1mb` | Same definitions as Table 1 |
| `icr_context` | Local ICR density class: `ICR-rich` (≥5), `Moderate cluster` (3–4), `Sparse` (2), `None` (<2) |
| `paired_21424` | Same definition as Table 1 |

### `table3_conflicting_unknown.tsv`

Genes present in the ASE data with ASE evidence and imprinting-relevant genomic context but classified as "Conflicting Data" or "Unknown" in geneimprint.com's Homo sapiens catalog rather than "Imprinted." Sorted by counted-allele direction consistency to group variable-direction candidates (more compatible with parent-of-origin effects) separately from perfect-direction candidates (more compatible with stable cis-regulatory ASE).

| Column | Definition |
|---|---|
| `rank` | Priority rank |
| `gene` | HGNC gene symbol |
| `chrom` | Chromosome |
| `geneimprint_status` | `Conflicting Data` or `Unknown` (from geneimprint.com) |
| Standard ASE columns | Same definitions as Tables 1 and 2 |
| `framework_interpretation` | Text interpretation applying the paper's framework (variable direction → parent-of-origin-compatible; perfect direction → cis-regulatory-compatible) |

### `supp_table1_scoring_rubric.tsv`

Point values for each feature in the prioritization score, plus tier-assignment thresholds. Used by the prioritization workflow to compute `priority_score` for every gene.

| Column | Definition |
|---|---|
| `feature` | Feature category (e.g., ICR overlap, sample support, imbalance) |
| `threshold` | Value or range that triggers the associated point award |
| `points` | Points contributed to the total priority score |

### `supp_table2_tissue_restricted_readcounts.tsv`

Per-sample gene-body read counts for the 13 tissue-restricted benchmark imprinted genes, obtained via `samtools view -c -F 260` on WASP-filtered BAMs. Classification (depth-limited vs SNP-coverage-limited) is included as the final column.

| Column | Definition |
|---|---|
| `gene` | HGNC gene symbol |
| `locus` | GRCh38 gene coordinates in `chr:start-end` format |
| `ase_tissue` | Tissue in which the gene was ASE-informative (`Brain` or `Liver`) |
| Per-sample columns | Read counts (primary alignments) overlapping the gene body in each of the six WASP-filtered BAMs |
| `median_brain` | Median read count across brain samples |
| `median_liver` | Median read count across liver samples |
| `ratio` | Median count in ASE-informative tissue / median count in uninformative tissue |
| `classification` | `Depth-limited` (>5× tissue difference) or `SNP-coverage-limited` (≤5×) |

### `supp_table3_extended_audit.tsv`

Full audit output for established human imprinted genes queried from the geneimprint.com Homo sapiens catalog. Includes recovered genes, unrecovered genes with reason for absence, and per-sample read counts.

| Column | Definition |
|---|---|
| `gene` | HGNC gene symbol (or alias as matched in GENCODE v38) |
| `region` | GRCh38 gene coordinates |
| `chromosome_cluster` | Imprinted cluster or genomic region (e.g., `11p15`, `15q11`, `14q32`) |
| Per-sample read count columns | From `samtools view -c -F 260` |
| `n_het_snps` | Number of heterozygous SNP sites in the gene body across the cohort |
| `recovered_in_ase_table` | `TRUE` if the gene has an entry in the master ASE gene table; `FALSE` otherwise |
| `classification` | One of: `recovered`, `nested_multi_gene_excluded`, `expressed_with_het_snps_not_in_ase_table`, `expressed_no_het_snps_in_cohort`, `expressed_below_aggregation_thresholds`, `not_expressed`, `not_annotated_in_gencode` |
| `nested_in` | For genes with `classification = nested_multi_gene_excluded`, the gene whose body it lies within |
| `notes` | Free-text notes (e.g., tissue-specificity, alias resolution, curator observations) |

## Table generation

All tables in this directory were produced by downstream R analysis of the SLURM pipeline outputs. The SLURM pipeline (see [`../slurm/`](../slurm/)) produces per-sample ASE counts (`analysis/ase_counts/`), WASP-filtered BAMs (`analysis/wasp/final/`), and the audit TSV (`analysis/audit/gene_body_audit.tsv`). The R analysis that transforms these into the manuscript tables is described in the manuscript Methods and Supplementary Methods 2, with the exact R environment documented in [`../environment/R_sessionInfo.txt`](../environment/R_sessionInfo.txt). R analysis code is available on request from the repository maintainer.

## Citation

If you use tables from this directory, please cite:

> Venkat V, Skaar D, Planchart A, Jirtle RL, Murphy SK, Tzeng J-Y, Hoyo C. *A conservative allele-specific expression framework prioritizes candidate imprinted loci in human conceptal brain and liver using the human imprintome.* <!-- FILL IN: journal, year, volume, pages, DOI -->

Also cite the underlying human imprintome resource:

> Jima DD, Skaar DA, Planchart A, et al. (2022) Genomic map of candidate human imprint control regions: the imprintome. *Epigenetics* 17(13):1920-1943. doi:10.1080/15592294.2022.2091815
