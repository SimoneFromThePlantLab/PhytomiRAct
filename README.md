# PhytomiRAct

**PhytomiRAct** is a standalone R framework to infer tissue-level plant microRNA (miRNA) activity from annotated single-cell RNA-seq datasets.

Most plant single-cell RNA-seq datasets profile polyadenylated mRNAs rather than mature small RNAs. PhytomiRAct therefore estimates miRNA activity indirectly by combining transcriptome-wide prediction of miRNA target sites with tissue-level expression deviations of predicted target genes.

## Overview

PhytomiRAct takes three main inputs:

1. an annotated **Seurat** object containing raw counts and tissue or cell-type metadata;
2. a mature miRNA sequence;
3. a species-matched reference cDNA transcriptome FASTA file used to scan for putative miRNA target sites.

For each miRNA, PhytomiRAct:

1. scans the reference cDNA transcriptome for candidate miRNA target sites;
2. collapses predicted target sites to gene-level target scores;
3. computes tissue-level expression deviations from single-cell mRNA profiles;
4. derives a miRNA activity score for each annotated tissue or cell type;
5. generates UMAP activity maps and ranked tissue-level summaries.

## Repository structure

```text
PhytomiRAct/
├── PhytomiRAct.R                  # Main standalone R script
├── examples/
│   └── example_run.R              # Minimal example for running the framework
├── data-raw/
│   ├── dataset_accessions.tsv     # Dataset accession template
│   ├── miRNA_sequences.tsv        # miRNA sequence table template
│   └── reference_transcriptomes.md# Notes on reference cDNA FASTA files
├── manuscript_runs/
│   └── README.md                  # How to organize manuscript-specific runs
├── results/
│   └── README.md                  # Placeholder for output files
├── CITATION.cff
├── LICENSE
├── .gitignore
└── README.md
```

## Installation

PhytomiRAct is currently distributed as a standalone R script. Clone or download this repository, then source the main script from R:

```r
source("PhytomiRAct.R")
```

Required R packages:

```r
install.packages(c(
  "Seurat",
  "Matrix",
  "dplyr",
  "tidyr",
  "tibble",
  "matrixStats",
  "readr",
  "stringr",
  "ggplot2"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("Biostrings")
```

## Quick start

```r
library(Seurat)
source("PhytomiRAct.R")

seurat_obj <- readRDS("your_annotated_seurat_object.rds")

res <- runPhytomiRAct(
  seurat_obj = seurat_obj,
  mirna_sequence = "MATURE_MIRNA_SEQUENCE",
  mirna_name = "miRNA_name",
  tissue_col = "tissue",
  cdna_fasta = "species_reference_cdna.fa",
  output_dir = "PhytomiRAct_results",
  assay = NULL,
  reduction = "umap",
  min_counts_per_cell = 1000,
  target_stringency = "permissive",
  use_precomputed_sites = TRUE,
  max_cells_per_tissue = NULL,
  downsample_seed = 123
)
```

## Input requirements

### 1. Seurat object

The input object must be a Seurat object and must contain:

- raw counts in the selected assay;
- a metadata column containing tissue or cell-type annotations;
- a dimensional reduction, such as `umap`, if UMAP activity visualization is required.

By default, PhytomiRAct expects the annotation column to be named `tissue`, but this can be changed using `tissue_col`.

### 2. Mature miRNA sequence

The mature miRNA sequence should be provided as a character string. Both RNA-style sequences with `U` and DNA-style sequences with `T` are accepted, because `U` is internally converted to `T`.

### 3. Reference cDNA FASTA

The reference FASTA should correspond to the cDNA transcriptome of the species under analysis. FASTA headers should preferably follow this pipe-separated format:

```text
>gene_id|transcript_id|gene_symbol|strand
```

Example:

```text
>AT1G01010|AT1G01010.1|NAC001|1
```

The first field, `gene_id`, must match the gene identifiers used as row names in the Seurat object.

## Output files

For each run, PhytomiRAct writes the following files to `output_dir`:

```text
<miRNA>_<stringency>_target_sites.tsv
<miRNA>_<stringency>_gene_targets.csv
delta_by_tissue.csv
activity_by_tissue.csv
UMAP_activity_tissue_level.tiff
dotplot_activity_tissue_level.tiff
seurat_with_activity.rds
```

The function also returns an R list containing the modified Seurat object, target-site table, gene-level target scores, tissue-level activity table, plots, output directory and run parameters.

## Target-detection stringency

PhytomiRAct currently supports two settings:

- `permissive`: allows broader target discovery;
- `strict`: applies more restrictive mismatch and wobble constraints.

The selected setting is controlled with:

```r
target_stringency = "permissive"
# or
 target_stringency = "strict"
```

## Data availability

Large single-cell datasets, Seurat objects, FASTA files and intermediate output files are not stored in this repository. Dataset accessions and reference transcriptome sources should be reported in `data-raw/dataset_accessions.tsv` and `data-raw/reference_transcriptomes.md`.

## Citation

If you use PhytomiRAct, please cite the associated manuscript and this repository. Citation information is provided in `CITATION.cff`.

## License

This project is released under the MIT License. See `LICENSE` for details.
