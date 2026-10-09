# Manuscript-specific runs

This folder can be used to store scripts that reproduce the analyses reported in the manuscript.

Suggested organization:

```text
01_run_Arabidopsis_miR399.R
02_run_Arabidopsis_miR165_166.R
03_run_rice_benchmark_miRNAs.R
04_run_rice_conserved_miRNAs.R
05_run_maize_benchmark_miRNAs.R
06_run_maize_conserved_miRNAs.R
```

Large input files such as Seurat objects, FASTA files and output folders should not be committed to GitHub. Instead, report public accessions and download instructions in `data-raw/` and in the manuscript Data Availability statement.
