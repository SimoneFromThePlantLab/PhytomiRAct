# PhytomiRAct

**Target-guided inference of tissue-level miRNA activity from plant single-cell transcriptomes**

PhytomiRAct is an R-based workflow for inferring plant microRNA (miRNA) activity from annotated single-cell RNA-seq datasets. The tool was designed for plant scRNA-seq atlases in which mature miRNAs are usually not directly measured, but their regulatory activity can be inferred from the expression behaviour of predicted target transcripts.

PhytomiRAct takes as input an annotated Seurat object, a mature miRNA sequence, and a species-matched reference cDNA transcriptome FASTA file. It scans reference transcripts for candidate plant-compatible miRNA target sites, collapses site-level predictions into gene-level target scores, and integrates these scores with tissue-resolved expression deviations to generate tissue-level miRNA activity profiles.

## Key features

- Infers tissue- or cell-type-level miRNA activity from plant single-cell mRNA data.
- Works directly with annotated Seurat objects.
- Uses mature miRNA sequences and species-specific cDNA transcriptomes for target-site scanning.
- Implements permissive and strict plant-oriented target-detection settings.
- Produces ranked tissue-level activity summaries and UMAP activity maps.
- Exports predicted target sites, gene-level target scores, expression deviation tables, activity scores, and Seurat objects containing inferred activity metadata.
- Designed for reuse across plant species with available scRNA-seq data and reference transcriptomes.

## Repository structure

```text
PhytomiRAct/
├── PhytomiRAct.R
├── README.md
├── LICENSE
├── CITATION.cff
├── NEWS.md
├── examples/
│   └── example_run.R
├── data-raw/
│   ├── miRNA_sequences.tsv
│   ├── dataset_accessions.tsv
│   └── reference_transcriptomes.md
├── manuscript_runs/
│   └── README.md
└── results/
    └── README.md
```

Large Seurat objects, reference FASTA files, and generated result files are not included in the repository. Dataset accession numbers, miRNA sequences, and reference transcriptome information are provided in `data-raw/`.

## Installation

PhytomiRAct is distributed as a standalone R script. Clone or download this repository and source the script in R.

```r
git clone https://github.com/SimoneFromThePlantLab/PhytomiRAct.git
```

In R:

```r
source("PhytomiRAct.R")
```

## Required R packages

PhytomiRAct requires the following R packages:

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

## Required inputs

PhytomiRAct requires three main inputs:

1. **Annotated Seurat object**
   - Raw counts must be available in the selected assay.
   - The object must contain dimensional reduction coordinates, typically UMAP.
   - The metadata must include a tissue or cell-type annotation column, for example `tissue`.

2. **Mature miRNA sequence**
   - RNA sequence in 5′ to 3′ orientation.
   - Uracil bases are automatically converted to thymine before transcriptome scanning.

3. **Species-specific reference cDNA FASTA file**
   - The reference FASTA should correspond to the annotated cDNA transcriptome of the species under analysis.
   - This is not a transcriptome reconstructed from the single-cell dataset.
   - Gene identifiers in the FASTA headers must match the gene identifiers used as row names in the Seurat object.

Recommended FASTA header format:

```text
>gene_id|transcript_id|gene_symbol|strand
```

Example:

```text
>AT1G01010|AT1G01010.1|NAC001|1
```

See `data-raw/reference_transcriptomes.md` for reference transcriptomes used in the manuscript analyses.

## Quick start

```r
library(Seurat)

source("PhytomiRAct.R")

seurat_obj <- readRDS("your_annotated_seurat_object.rds")

res <- runPhytomiRAct(
  seurat_obj = seurat_obj,
  mirna_sequence = "UGCCAAAGGAGAUUUGCCCAG",
  mirna_name = "osa-miR399e-3p",
  tissue_col = "tissue",
  cdna_fasta = "reference_cdna_transcriptome.fa",
  output_dir = "results/osa_miR399e_3p",
  assay = "RNA",
  reduction = "umap",
  min_counts_per_cell = 1000,
  target_stringency = "strict",
  use_precomputed_sites = TRUE
)
```

## Main function

```r
runPhytomiRAct(
  seurat_obj,
  mirna_sequence,
  mirna_name,
  tissue_col = "tissue",
  cdna_fasta,
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

### Main arguments

| Argument | Description |
|---|---|
| `seurat_obj` | Annotated Seurat object containing raw counts and tissue or cell-type metadata. |
| `mirna_sequence` | Mature miRNA sequence in RNA alphabet. |
| `mirna_name` | Name used for output files and plot titles. |
| `tissue_col` | Metadata column containing biological tissue or cell-type annotations. |
| `cdna_fasta` | Species-specific reference cDNA FASTA file. |
| `output_dir` | Directory where results will be written. |
| `assay` | Seurat assay to use. If `NULL`, the default assay is used. |
| `reduction` | Dimensional reduction used for plotting, usually `umap`. |
| `min_counts_per_cell` | Minimum total counts required for a cell to be retained. |
| `target_stringency` | Target detection mode: `permissive` or `strict`. |
| `use_precomputed_sites` | If `TRUE`, previously generated target-site tables are reused when present. |
| `max_cells_per_tissue` | Optional downsampling limit per tissue or cell type. |
| `downsample_seed` | Random seed used for optional downsampling. |

## Target-detection stringency

PhytomiRAct implements two target-detection modes:

- `strict`: more conservative target-site detection.
- `permissive`: broader target-site detection allowing more pairing events.

Both modes use a plant-oriented complementarity model that considers extended miRNA-target pairing, mismatch position, G:U wobble pairs, and pairing around the canonical cleavage region.

## Output files

Each run generates the following files in `output_dir`:

| File | Description |
|---|---|
| `<miRNA>_<stringency>_target_sites.tsv` | Predicted transcript-level target sites. |
| `<miRNA>_<stringency>_gene_targets.csv` | Gene-level target scores collapsed from site-level predictions. |
| `delta_by_tissue.csv` | Gene-level expression deviation values by tissue or cell type. |
| `activity_by_tissue.csv` | Ranked tissue-level miRNA activity scores. |
| `UMAP_activity_tissue_level.tiff` | UMAP projection coloured by inferred tissue-level miRNA activity. |
| `dotplot_activity_tissue_level.tiff` | Ranked tissue-level activity dot plot. |
| `seurat_with_activity.rds` | Seurat object with inferred activity added to metadata. |

The dot plot reports activity rankings only; it does not include significance asterisks.

## Manuscript analyses

The manuscript analyses used PhytomiRAct on Arabidopsis, rice, and maize single-cell datasets.

### Arabidopsis thaliana

- `ath-miR399` in shoot and root datasets.
- `ath-miR165/166` in root dataset.

### Oryza sativa

- `osa-miR393a`
- `osa-miR399e-3p`
- `osa-miR166a-3p`

### Zea mays

- `zma-miR390a-5p`
- `zma-miR399a-3p`
- `zma-miR166j-k-n-3p`

The mature miRNA sequences used in the analyses are listed in `data-raw/miRNA_sequences.tsv`.
Dataset accessions are listed in `data-raw/dataset_accessions.tsv`.
Reference transcriptome details are listed in `data-raw/reference_transcriptomes.md`.

## Interpreting PhytomiRAct scores

PhytomiRAct estimates regulatory activity from the collective behaviour of predicted target transcripts. It does not directly measure mature miRNA abundance.

Activity values should therefore be interpreted as relative tissue-level rankings within each miRNA-dataset combination. They are not intended as absolute activity values for direct comparison across unrelated miRNAs, species, or datasets.

## Data availability

Large input files are not distributed through this repository. Users should download the original single-cell datasets and reference transcriptomes from the sources listed in `data-raw/dataset_accessions.tsv` and `data-raw/reference_transcriptomes.md`.

## License

This project is released under the MIT License. See `LICENSE` for details.

## Contact

For questions, suggestions, or issues, please open a GitHub issue or contact the corresponding author listed in the manuscript.
