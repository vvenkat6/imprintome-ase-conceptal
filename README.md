# imprintome-ase-conceptal

Allele-specific expression (ASE) analysis of human conceptal brain and liver, anchored on the human imprintome map of candidate imprinting control regions (ICRs) from Jima et al. (2022).

This repository contains the full pipeline, analysis scripts, and output tables supporting the manuscript:

> Venkat V, Skaar D, Planchart A, Jirtle RL, Murphy SK, Tzeng J-Y, Hoyo C. *A conservative allele-specific expression framework prioritizes candidate imprinted loci in human conceptal brain and liver using the human imprintome.* (Manuscript in preparation.)

<!-- FILL IN: journal, year, DOI once accepted -->

## Overview

The pipeline processes six RNA-seq libraries (three brain, three liver) from five conceptuses previously used to define the human imprintome via whole-genome bisulfite sequencing. Reads are aligned to GRCh38 with STAR in two-pass mode, processed for variant discovery with GATK, corrected for reference mapping bias with WASP, and quantified for allelic imbalance with GATK ASEReadCounter. Gene-level ASE is combined across samples with Fisher's method, annotated against tier 1 (n=332) and tier 2 (n=1,488) ICRs from the human imprintome, and prioritized with a prespecified additive scoring rubric.

See Figure 1 of the manuscript for the workflow diagram; the Mermaid source is in `scripts/figure1_workflow.mmd`.

## Repository structure

```
imprintome-ase-conceptal/
├── README.md                    This file
├── LICENSE                      MIT license
├── .gitignore
├── environment/                 Software version pinning
│   ├── modules.txt              HPC module load commands
│   └── R_sessionInfo.txt        R package versions
├── slurm/                       SLURM scripts, run in numeric order
│   ├── 01_star_firstpass.slurm
│   ├── 02_merge_splice_junctions.sh
│   ├── 03_star_secondpass.slurm
│   ├── 04_readgroups_dedup.slurm
│   ├── 05_haplotypecaller_gvcf.slurm
│   ├── 06_joint_genotype.slurm
│   ├── 07_wasp_intersecting_snps.slurm
│   ├── 08_wasp_remap.slurm
│   ├── 09_wasp_filter_merge.slurm
│   ├── 10_asereadcounter.slurm
│   └── 11_gene_body_audit.slurm
├── config/                      Configuration and metadata
│   ├── samples.tsv              Sample manifest
│   └── benchmark_gene_list.txt  Established human imprinted genes queried
└── tables/                      Final output tables (Table 1, 2, 3 + supplemental)
    └── README.md
```

## Pipeline order

The scripts in `slurm/` are numbered and should be run in order. Each script processes all six samples as an array job (`--array=1-6`), reading sample IDs from `config/samples.tsv`.

| Step | Script | Purpose | Input | Output |
|---|---|---|---|---|
| 1 | `01_star_firstpass.slurm` | STAR first-pass alignment for splice-junction discovery | FASTQ | Per-sample `SJ.out.tab` |
| 2 | `02_merge_splice_junctions.sh` | Merge per-sample splice junctions | `SJ.out.tab` × 6 | `all_SJ.out.tab` |
| 3 | `03_star_secondpass.slurm` | STAR second-pass alignment with merged junctions | FASTQ + `all_SJ.out.tab` | Coordinate-sorted BAM |
| 4 | `04_readgroups_dedup.slurm` | Add read groups and mark duplicates | Sorted BAM | Deduplicated BAM |
| 5 | `05_haplotypecaller_gvcf.slurm` | Per-sample variant calling in GVCF mode | Deduplicated BAM | Per-sample GVCF |
| 6 | `06_joint_genotype.slurm` | Cohort joint genotyping, SNP selection | Six GVCFs | `cohort.snps.vcf.gz` |
| 7 | `07_wasp_intersecting_snps.slurm` | Identify SNP-overlapping reads for remapping | Deduplicated BAM + SNP tables | Keep BAM + to-remap BAM |
| 8 | `08_wasp_remap.slurm` | Remap allele-swapped reads with STAR | To-remap BAM | Remapped BAM |
| 9 | `09_wasp_filter_merge.slurm` | Filter and merge remapped reads | Keep + remapped BAMs | WASP-filtered BAM |
| 10 | `10_asereadcounter.slurm` | Count reference and alternate alleles at het SNPs | WASP BAM + VCF | Per-sample ASE table |
| 11 | `11_gene_body_audit.slurm` | Gene-body read counts for audit (Figure 2, Supp Table S2, S3) | WASP BAMs + gene coordinates | Per-gene per-sample read count table |

After the SLURM pipeline, the R scripts in `scripts/` produce the analysis tables and figures:

| Script | Purpose | Output |
|---|---|---|
| `snp_to_gene_mapping.R` | Assign SNPs to gene bodies, exclude multi-gene SNPs | SNP-to-gene table |
| `gene_level_ase_stats.R` | Aggregate to gene level, binomial + Fisher + BH FDR | Master gene ASE table |
| `icr_annotation.R` | Annotate genes with tier 1/2 ICR overlap, distance, density | ICR-annotated gene table |
| `prioritization_scoring.R` | Additive prioritization scoring (see Supp Table S1) | Table 2 |
| `benchmark_audit.R` | Cross-check against geneimprint.com human catalog | Table 1, Table 3, Supp Table S3 |
| `figure2_readcount_plot.R` | Dot plot of tissue-restricted benchmark read counts | Figure 2 |

## Data

Raw and processed sequencing data are deposited at NCBI Gene Expression Omnibus under accession [GSE288468](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288468).

ICR interval files (tier 1 and tier 2) are from Jima et al. (2022), *Epigenetics*, 17(13):1920-1943. See Supplementary Materials of that publication for BED files.

Reference genome and annotation:
- GRCh38 primary assembly FASTA (GENCODE)
- GENCODE v38 comprehensive annotation GTF

## Software versions

<!-- FILL IN: run `<tool> --version` on your HPC and paste the exact output here -->

| Tool | Version | Purpose |
|---|---|---|
| STAR | 2.7.10a | RNA-seq alignment |
| GATK | 4.4.0.0 | Variant calling, ASEReadCounter, MarkDuplicates |
| samtools | 1.17 | BAM manipulation |
| bcftools | 1.17 | VCF manipulation |
| WASP | <!-- FILL IN: commit hash --> | Reference mapping-bias correction |
| Python | 3.9+ | WASP scripts |
| R | 4.2+ | Statistical analysis and figures |

Full HPC module load commands are in `environment/modules.txt`. R package versions are in `environment/R_sessionInfo.txt`.

## Reproducing the analysis

The full pipeline was run on the North Carolina State University Hazel HPC cluster with SLURM. The scripts in `slurm/` are written for that environment and will need adaptation for other HPC systems (path variables, module load commands, memory/thread allocations).

To reproduce the analysis from raw FASTQ:

1. Download raw FASTQ from GEO GSE288468 into a directory specified as `FASTQ_DIR` in the SLURM scripts.
2. Build the STAR genome index from GRCh38 + GENCODE v38 GTF (see `slurm/01_star_firstpass.slurm` header for parameters).
3. Edit paths in each SLURM script to match your environment.
4. Submit scripts 01 through 11 in order.
5. Run the R scripts in `scripts/` from the R console after the SLURM pipeline completes.

To reproduce only the analysis (starting from pre-computed ASE counts and ICR intervals), skip to the R scripts in `scripts/`. The final output tables in `tables/` are provided for reference and comparison.

## Citation

If you use code or output from this repository, please cite:

> Venkat V, Skaar D, Planchart A, Jirtle RL, Murphy SK, Tzeng J-Y, Hoyo C. *A conservative allele-specific expression framework prioritizes candidate imprinted loci in human conceptal brain and liver using the human imprintome.* <!-- FILL IN: journal, year, volume, pages, DOI -->

Also cite the underlying human imprintome resource:

> Jima DD, Skaar DA, Planchart A, Motsinger-Reif A, Cevik SE, Park SS, Cowley M, Wright F, House J, Liu A, Jirtle RL, Hoyo C. (2022) Genomic map of candidate human imprint control regions: the imprintome. *Epigenetics* 17(13):1920-1943. doi:10.1080/15592294.2022.2091815

An archival snapshot of this repository at the time of publication is available at Zenodo: <!-- FILL IN: DOI once generated at submission -->

## Contact

- **Scientific correspondence:** David Skaar, Department of Biological Sciences, North Carolina State University 
- **Repository and code questions:** Vaishnavi Venkat, Bioinformatics Research Center, North Carolina State University 

For questions about the human imprintome map itself, contact the corresponding authors of Jima et al. (2022).

## License

This code is released under the MIT License (see `LICENSE`). The analysis outputs and derived tables are released under CC-BY-4.0. Please cite the manuscript above when reusing.

## Acknowledgments

This work was supported by NIH grants R01HD098857, R01ES093351, R01MD017696, and RF1AG074328/R01AG074328. Sequencing was performed at the North Carolina State University Genome Services Laboratory. HPC resources were provided by the NCSU Hazel cluster and BRC cluster.
