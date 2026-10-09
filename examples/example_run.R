############################################################
# Minimal PhytomiRAct example
############################################################

library(Seurat)

# Source the standalone framework
source("PhytomiRAct.R")

# Load an annotated Seurat object
# The object should contain raw counts, a UMAP reduction, and a tissue/cell-type metadata column.
seurat_obj <- readRDS("your_annotated_seurat_object.rds")

# Run PhytomiRAct
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

# Inspect tissue-level activity
head(res$activity_by_tissue)
