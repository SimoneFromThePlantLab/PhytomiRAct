# Reference cDNA transcriptomes

PhytomiRAct requires a species-matched reference cDNA transcriptome FASTA file for target-site scanning.

The reference FASTA is not a transcriptome reconstructed from the single-cell dataset. It should correspond to the annotated reference cDNA transcriptome of the species under analysis.

Recommended FASTA header format:

```text
>gene_id|transcript_id|gene_symbol|strand
```

Example:

```text
>AT1G01010|AT1G01010.1|NAC001|1
```

The `gene_id` field must match the gene identifiers used as row names in the Seurat object.

## Transcriptomes used in the manuscript

| Species | Reference transcriptome | Source/version | FASTA used | Notes |
|---|---|---|---|---|
| *Arabidopsis thaliana* | TAIR/Araport cDNA reference | TAIR10 genome; Araport11 annotation | `cDNA_arabidopsis_TAIR10_canonical.fa` | Canonical Arabidopsis cDNA reference used to scan transcripts for ath-miR399 and ath-miR165/166 target sites. Headers were formatted as `gene_id|transcript_id|gene_symbol|strand`, with `ATxGxxxxx` identifiers matching the Seurat object row names. |
| *Oryza sativa* | Ensembl Plants cDNA reference | IRGSP-1.0 | `cDNA_rice_IRGSP1.0_EnsemblPlants_cdna_all_PhytomiRAct.fa` | Species-matched rice cDNA transcriptome used to scan reference transcripts for osa-miR393a, osa-miR399e-3p and osa-miR166a-3p target sites. Raw Ensembl Plants headers were reformatted for PhytomiRAct compatibility. |
| *Zea mays* | B73 RefGen v4 cDNA reference | Zm-B73-REFERENCE-GRAMENE-4.0 / B73 RefGen v4 / Zm00001d gene model series | `cDNA_maize_B73_RefGen_v4_PhytomiRAct.fa` or `cDNA_maize_B73_RefGen_v4_canonical_like_PhytomiRAct.fa` | Species-matched maize cDNA transcriptome used to scan reference transcripts for zma-miR390a-5p, zma-miR399a-3p and zma-miR166j-k-n-3p target sites. Headers were reformatted for PhytomiRAct compatibility; the canonical-like version keeps the longest transcript per gene to reduce redundancy. |

## Header conversion notes

For each species, FASTA headers should be converted before running PhytomiRAct so that the first pipe-delimited field corresponds to the same gene identifiers used in the Seurat expression matrix. This matching step is required to join predicted target sites with expression-deviation values.

Example compatible headers:

```text
>AT1G01010|AT1G01010.1|NAC001|1
>Os01g0208000|Os01t0208000-01|OsUMAMIT4|1
>Zm00001d020622|Zm00001d020622_T004|NA|NA
```

Large reference FASTA files are not included in this repository. Users should download species-specific cDNA FASTA files from the corresponding reference resources and reformat headers as described above.
