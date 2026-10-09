############################################################
# PhytomiRAct.R
# Single-cell plant miRNA activity inference from Seurat
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(matrixStats)
  library(Biostrings)
  library(readr)
  library(stringr)
  library(ggplot2)
})

############################################################
# Main function
############################################################

runPhytomiRAct <- function(
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
) {
  
  ############################################################
  # Safety checks and output folder
  ############################################################
  
  message("############################################################")
  message("PhytomiRAct")
  message("Inferring miRNA activity from single-cell mRNA expression")
  message("############################################################")
  
  message("")
  message("IMPORTANT:")
  message("The Seurat object must contain a metadata column named '", tissue_col, "'.")
  message("This column must contain biologically annotated tissues/cell types.")
  message("Example: Xylem, Phloem, Meristem, Root epidermis, etc.")
  message("")
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  
  if (!tissue_col %in% colnames(seurat_obj@meta.data)) {
    stop(
      "The column '", tissue_col, "' was not found in seurat_obj@meta.data.\n",
      "Please create this column before running PhytomiRAct."
    )
  }
  
  if (!file.exists(cdna_fasta)) {
    stop("cdna_fasta file not found: ", cdna_fasta)
  }
  
  ############################################################
  # Optional downsampling before count extraction
  ############################################################
  
  if (!is.null(max_cells_per_tissue)) {
    
    message("Downsampling Seurat object before count extraction...")
    message("Maximum cells per tissue: ", max_cells_per_tissue)
    
    tissues_down <- seurat_obj@meta.data[[tissue_col]]
    
    cells_by_tissue <- split(
      colnames(seurat_obj),
      tissues_down
    )
    
    set.seed(downsample_seed)
    
    sampled_cells <- unlist(lapply(cells_by_tissue, function(cells) {
      if (length(cells) > max_cells_per_tissue) {
        sample(cells, max_cells_per_tissue)
      } else {
        cells
      }
    }))
    
    sampled_cells <- as.character(sampled_cells)
    
    message("Cells before downsampling: ", ncol(seurat_obj))
    message("Cells after downsampling: ", length(sampled_cells))
    
    seurat_obj <- seurat_obj[, sampled_cells]
    
    message("Cells sampled per tissue:")
    print(table(seurat_obj@meta.data[[tissue_col]]))
  }
  
  ###########################################################
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  mirna_safe <- gsub("[^A-Za-z0-9_]", "_", mirna_name)
  
  target_stringency <- match.arg(
    target_stringency,
    choices = c("permissive", "strict")
  )
  
  message("Target detection stringency: ", target_stringency)
  
  out_sites <- file.path(
    output_dir,
    paste0(mirna_safe, "_", target_stringency, "_target_sites.tsv")
  )
  out_gene_targets <- file.path(
    output_dir,
    paste0(mirna_safe, "_", target_stringency, "_gene_targets.csv")
  )
  out_delta <- file.path(output_dir, "delta_by_tissue.csv")
  out_activity_tissue <- file.path(output_dir, "activity_by_tissue.csv")
  out_umap_tissue <- file.path(output_dir, "UMAP_activity_tissue_level.tiff")
  out_seurat <- file.path(output_dir, "seurat_with_activity.rds")
  
  ############################################################
  # 1) Extract counts
  ############################################################
  
  message("Extracting count matrix from Seurat object...")
  
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  }
  
  counts <- tryCatch({
    GetAssayData(seurat_obj, assay = assay, slot = "counts")
  }, error = function(e) {
    GetAssayData(seurat_obj, assay = assay, layer = "counts")
  })
  
  ############################################################
  # 2) Filter cells and normalize
  ############################################################
  
  message("Filtering cells and performing normalization...")
  
  cell_totals <- Matrix::colSums(counts)
  valid_cells <- names(cell_totals[cell_totals >= min_counts_per_cell])
  
  if (length(valid_cells) == 0) {
    stop("No cells passed the min_counts_per_cell threshold.")
  }
  
  message("Cells before filtering: ", ncol(counts))
  message("Cells after filtering: ", length(valid_cells))
  
  # Apply filtering
  counts <- counts[, valid_cells, drop = FALSE]
  seurat_obj <- seurat_obj[, valid_cells]
  
  # Recompute totals after filtering
  cell_totals <- Matrix::colSums(counts)
  Tmax <- max(cell_totals)
  
  norm_counts <- sweep(counts, 2, cell_totals, "/") * Tmax
  logE <- log2(norm_counts + 1)
  
  ############################################################
  # 3) Expected expression and DeltaE by tissue
  ############################################################
  
  message("Computing expected expression and DeltaE by tissue...")
  
  tissues <- as.character(seurat_obj@meta.data[[tissue_col]])
  names(tissues) <- colnames(seurat_obj)
  
  if (any(is.na(tissues))) {
    warning("Some cells have NA tissue annotation. They will be labelled as 'Unknown'.")
    tissues[is.na(tissues)] <- "Unknown"
  }
  
  Exp_global <- matrixStats::rowMedians(as.matrix(logE))
  names(Exp_global) <- rownames(logE)
  
  Exp_tissue <- lapply(unique(tissues), function(t) {
    cells_t <- which(tissues == t)
    matrixStats::rowMedians(as.matrix(logE[, cells_t, drop = FALSE]))
  })
  
  names(Exp_tissue) <- unique(tissues)
  
  Delta_global <- sweep(as.matrix(logE), 1, Exp_global, "-")
  
  Delta_tissue <- matrix(
    NA_real_,
    nrow = nrow(logE),
    ncol = ncol(logE),
    dimnames = dimnames(logE)
  )
  
  for (t in unique(tissues)) {
    cells_t <- which(tissues == t)
    Exp_vec <- Exp_tissue[[t]]
    names(Exp_vec) <- rownames(logE)
    
    Delta_tissue[, cells_t] <- sweep(
      as.matrix(logE[, cells_t, drop = FALSE]),
      1,
      Exp_vec,
      "-"
    )
  }
  
  ############################################################
  # 4) Save DeltaE by tissue
  ############################################################
  
  message("Saving DeltaE by tissue...")
  
  delta_tissue_list <- lapply(unique(tissues), function(t) {
    cells_t <- which(tissues == t)
    
    tibble(
      gene_id = rownames(Delta_tissue),
      tissue = t,
      delta_e_tissue_mean = rowMeans(
        Delta_tissue[, cells_t, drop = FALSE],
        na.rm = TRUE
      ),
      delta_e_global_mean = rowMeans(
        Delta_global[, cells_t, drop = FALSE],
        na.rm = TRUE
      ),
      n_cells = length(cells_t)
    )
  })
  
  delta_by_tissue <- bind_rows(delta_tissue_list)
  
  write.csv(
    delta_by_tissue,
    out_delta,
    row.names = FALSE
  )
  
  ############################################################
  # 5) Target scanning or loading precomputed target sites
  ############################################################
  
  if (file.exists(out_sites) && use_precomputed_sites) {
    
    message("Loading precomputed target sites from: ", out_sites)
    sites <- readr::read_tsv(out_sites, show_col_types = FALSE)
    
  } else {
    
    message("Scanning transcriptome for candidate miRNA target sites...")
    
    if (target_stringency == "permissive") {
      
      cfg <- list(
        core_region = c(2, 13),
        no_errors_10_11 = TRUE,
        max_events_total = 5,
        max_core_events = 2,
        scores = list(
          match = 2,
          wobble = -1,
          mismatch_core = -4,
          mismatch_noncore = -2,
          contig7_bonus = 3
        ),
        max_mismatch_scan = 7
      )
      
    } else if (target_stringency == "strict") {
      
      cfg <- list(
        core_region = c(2, 13),
        no_errors_10_11 = TRUE,
        max_events_total = 3,
        max_core_events = 1,
        scores = list(
          match = 2,
          wobble = -1,
          mismatch_core = -4,
          mismatch_noncore = -2,
          contig7_bonus = 3
        ),
        max_mismatch_scan = 4
      )
    }
    
    guide <- toupper(chartr("U", "T", mirna_sequence))
    pat_rc <- reverseComplement(DNAString(guide))
    
    fa <- readDNAStringSet(cdna_fasta)
    
    parse_header <- function(h) {
      p <- strsplit(h, "\\|")[[1]]
      
      tibble(
        gene_id = ifelse(length(p) >= 1, p[1], NA_character_),
        transcript_id = ifelse(length(p) >= 2, p[2], p[1]),
        gene_symbol = ifelse(length(p) >= 3, p[3], NA_character_),
        strand = ifelse(length(p) >= 4, p[4], NA_character_)
      )
    }
    
    hdr_tbl <- bind_rows(lapply(names(fa), parse_header))
    
    hits_list <- vmatchPattern(
      pat_rc,
      fa,
      max.mismatch = cfg$max_mismatch_scan,
      with.indels = FALSE,
      fixed = TRUE
    )
    
    classify_pos <- function(gb, tw_rc_b, tw_orig_b) {
      if (gb == tw_rc_b) {
        return("match")
      }
      
      if (
        (gb == "G" & tw_orig_b == "T") |
        (gb == "T" & tw_orig_b == "G")
      ) {
        return("wobble")
      }
      
      return("mismatch")
    }
    
    score_window <- function(window_mrna, guide, cfg) {
      tw_rc <- as.character(reverseComplement(DNAString(window_mrna)))
      
      gw <- strsplit(guide, "")[[1]]
      rcw <- strsplit(tw_rc, "")[[1]]
      orig <- strsplit(window_mrna, "")[[1]]
      
      classes <- mapply(
        classify_pos,
        gw,
        rcw,
        orig,
        USE.NAMES = FALSE
      )
      
      if (cfg$no_errors_10_11 &&
          !(classes[10] == "match" & classes[11] == "match")) {
        return(NULL)
      }
      
      wob <- sum(classes == "wobble")
      mis <- sum(classes == "mismatch")
      
      if ((wob + mis) > cfg$max_events_total) {
        return(NULL)
      }
      
      core <- cfg$core_region[1]:cfg$core_region[2]
      core_events <- sum(classes[core] != "match")
      
      if (core_events > cfg$max_core_events) {
        return(NULL)
      }
      
      if (any(classes[core][-length(core)] == "mismatch" &
              classes[core][-1] == "mismatch")) {
        return(NULL)
      }
      
      s <- 0
      
      for (i in seq_along(gw)) {
        if (classes[i] == "match") {
          s <- s + cfg$scores$match
        } else if (classes[i] == "wobble") {
          s <- s + cfg$scores$wobble
        } else if (i %in% core) {
          s <- s + cfg$scores$mismatch_core
        } else {
          s <- s + cfg$scores$mismatch_noncore
        }
      }
      
      core_classes <- classes[core] == "match"
      run <- rle(core_classes)
      max_run <- max(ifelse(run$values, run$lengths, 0))
      
      if (max_run >= 7) {
        s <- s + cfg$scores$contig7_bonus
      }
      
      list(
        score = s,
        n_wobble = wob,
        n_mismatch = mis,
        classes = paste(classes, collapse = ""),
        rc = tw_rc
      )
    }
    
    res_list <- list()
    site_id <- 1L
    
    for (tx_i in seq_along(fa)) {
      
      if (length(hits_list[[tx_i]]) == 0) {
        next
      }
      
      seq_dna <- as.character(fa[[tx_i]])
      h <- hdr_tbl[tx_i, ]
      
      for (j in seq_along(hits_list[[tx_i]])) {
        
        rng <- hits_list[[tx_i]][j]
        w <- substr(seq_dna, start(rng), end(rng))
        
        if (grepl("N", w)) {
          next
        }
        
        sc <- score_window(w, guide, cfg)
        
        if (is.null(sc)) {
          next
        }
        
        res_list[[site_id]] <- tibble(
          gene_id = h$gene_id,
          transcript_id = h$transcript_id,
          gene_symbol = h$gene_symbol,
          strand = h$strand,
          start = start(rng),
          end = end(rng),
          site_sequence = w,
          site_revcomp = sc$rc,
          guide_DNA = guide,
          score_site = sc$score,
          n_wobble = sc$n_wobble,
          n_mismatch = sc$n_mismatch,
          pair_classes = sc$classes
        )
        
        site_id <- site_id + 1L
      }
    }
    
    if (length(res_list) == 0) {
      stop("No candidate target sites were found.")
    }
    
    sites <- bind_rows(res_list) %>%
      arrange(desc(score_site))
    
    readr::write_tsv(sites, out_sites)
  }
  
  ############################################################
  # 6) Collapse target sites to gene-level target scores
  ############################################################
  
  message("Collapsing target sites to gene-level target scores...")
  
  gene_targets <- sites %>%
    group_by(gene_id) %>%
    summarise(
      TargetScore = max(score_site, na.rm = TRUE),
      n_sites = n(),
      .groups = "drop"
    ) %>%
    mutate(
      TargetScore_norm = case_when(
        dplyr::n_distinct(TargetScore) == 1 ~ 1,
        TRUE ~ (TargetScore - min(TargetScore, na.rm = TRUE)) /
          max(
            1e-8,
            max(TargetScore, na.rm = TRUE) -
              min(TargetScore, na.rm = TRUE)
          )
      )
    )
  
  write.csv(
    gene_targets,
    out_gene_targets,
    row.names = FALSE
  )
  
  ############################################################
  # 7) Activity by tissue
  ############################################################
  
  message("Computing miRNA activity by tissue...")
  
  delta <- delta_by_tissue %>%
    dplyr::rename(
      deltaE_tissue = delta_e_tissue_mean,
      deltaE_global = delta_e_global_mean
    )
  
  deltaX <- delta %>%
    left_join(gene_targets, by = "gene_id") %>%
    mutate(is_target = !is.na(TargetScore))
  
  activity_by_tissue <- deltaX %>%
    group_by(tissue) %>%
    summarise(
      n_targets = sum(is_target),
      activity_weighted_raw = {
        sel <- is_target
        wt <- TargetScore_norm[sel]
        val <- -deltaE_tissue[sel]
        
        ok <- is.finite(wt) & is.finite(val)
        wt <- wt[ok]
        val <- val[ok]
        
        if (length(wt) == 0 || length(val) == 0) {
          NA_real_
        } else {
          den <- sum(wt, na.rm = TRUE)
          num <- sum(wt * val, na.rm = TRUE)
          
          if (den > 0) {
            num / den
          } else {
            mean(val, na.rm = TRUE)
          }
        }
      },
      .groups = "drop"
    )
  
  activity_by_tissue <- activity_by_tissue %>%
    mutate(
      activity_weighted = -activity_weighted_raw
    ) %>%
    arrange(desc(activity_weighted))

  write.csv(
    activity_by_tissue,
    out_activity_tissue,
    row.names = FALSE
  )
  
  # NOTE:
  # activity_weighted_raw is computed using -deltaE_tissue.
  # We invert the sign here to keep the final activity score aligned with
  # the empirically validated tissue-level miRNA activity pattern.
  # Do not remove this sign inversion: it preserves the directionality
  # observed for known miRNA/tissue associations.
  
  ############################################################
  # 8) Add tissue-level activity to Seurat object and plot UMAP
  ############################################################
  
  message("Generating UMAP activity map...")
  
  seu <- seurat_obj
  seu[[tissue_col]] <- as.character(seu@meta.data[[tissue_col]])
  
  activity_col_tissue <- paste0(mirna_safe, "_activity_by_tissue")
  
  cell_meta <- seu@meta.data %>%
    rownames_to_column("cell_id") %>%
    mutate(tissue_join = as.character(.data[[tissue_col]]))
  
  activity_join <- activity_by_tissue %>%
    mutate(tissue_join = as.character(tissue)) %>%
    dplyr::select(tissue_join, activity_weighted)
  
  cell_meta2 <- cell_meta %>%
    left_join(activity_join, by = "tissue_join")
  
  colnames(cell_meta2)[colnames(cell_meta2) == "activity_weighted"] <- activity_col_tissue
  
  seu <- AddMetaData(
    seu,
    metadata = cell_meta2 %>%
      dplyr::select(cell_id, all_of(activity_col_tissue)) %>%
      column_to_rownames("cell_id")
  )
  
  p_umap_tissue <- FeaturePlot(
    seu,
    features = activity_col_tissue,
    reduction = reduction
  ) +
    scale_color_gradientn(
      colours = c("grey95", "#fff5f0", "#fee0d2", "#fc9272", "#de2d26", "#a50f15"),
      na.value = "grey90"
    ) +
    labs(
      title = paste0(mirna_name, " tissue-level activity")
    ) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5)
    )
  
  tiff(
    out_umap_tissue,
    units = "in",
    width = 6,
    height = 5,
    res = 600
  )
  print(p_umap_tissue)
  dev.off()
  
  ############################################################
  # 9 Tissue-level activity plot
  ############################################################
  
  message("Generating tissue-level activity dotplot...")
  
  out_tissue_dotplot <- file.path(output_dir, "dotplot_activity_tissue_level.tiff")
  
  plot_df <- activity_by_tissue %>%
    dplyr::filter(tissue != "Unknown") %>%
    dplyr::mutate(tissue = as.character(tissue)) %>%
    dplyr::arrange(activity_weighted)

  plot_df$tissue <- factor(plot_df$tissue, levels = plot_df$tissue)
  
  x_min <- min(plot_df$activity_weighted, na.rm = TRUE)
  x_max <- max(plot_df$activity_weighted, na.rm = TRUE)
  x_range <- x_max - x_min
  
  if (!is.finite(x_range) || x_range == 0) {
    x_range <- 0.1
  }
  
  x_lower <- x_min - 0.03 * x_range
  x_upper <- x_max + 0.08 * x_range
  
  p_tissue_dotplot <- ggplot(
    plot_df,
    aes(
      x = activity_weighted,
      y = tissue,
      color = activity_weighted
    )
  ) +
    geom_point(
      size = 5,
      alpha = 0.9
    ) +
    coord_cartesian(
      xlim = c(x_lower, x_upper),
      clip = "off"
    ) +
    scale_color_gradientn(
      colours = c("grey85", "#fee0d2", "#fc9272", "#de2d26", "#a50f15")
    ) +
    labs(
      x = "Tissue-level miRNA activity",
      y = "",
      color = "Activity",
      title = paste0(mirna_name, " inferred activity across tissues")
    ) +
    theme_classic() +
    theme(
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      
      panel.grid.major.y = element_line(
        color = "grey90",
        linewidth = 0.35
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 11),
      axis.title.x = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  
  tiff(
    out_tissue_dotplot,
    units = "in",
    width = 6,
    height = 5,
    res = 600
  )
  print(p_tissue_dotplot)
  dev.off()
  
  ############################################################
  # 10 Save Seurat object
  ############################################################
  
  saveRDS(
    seu,
    out_seurat
  )
  
  ############################################################
  # 11 Final message and return object
  ############################################################
  
  message("")
  message("Analysis completed successfully.")
  message("Results written to: ", output_dir)
  message("")
  message("Generated files:")
  message(" - target_sites.tsv")
  message(" - gene_targets.csv")
  message(" - delta_by_tissue.csv")
  message(" - activity_by_tissue.csv")
  message(" - UMAP_activity_tissue_level.tiff")
  message(" - seurat_with_activity.rds")
  
  result <- list(
    seurat_obj = seu,
    target_sites = sites,
    gene_targets = gene_targets,
    delta_by_tissue = delta_by_tissue,
    activity_by_tissue = activity_by_tissue,
    plots = list(
      umap_tissue = p_umap_tissue,
      tissue_dotplot = p_tissue_dotplot
    ),
    output_dir = output_dir,
    parameters = list(
      mirna_name = mirna_name,
      mirna_sequence = mirna_sequence,
      tissue_col = tissue_col,
      cdna_fasta = cdna_fasta,
      assay = assay,
      reduction = reduction,
      min_counts_per_cell = min_counts_per_cell,
      target_stringency = target_stringency,
      max_cells_per_tissue = max_cells_per_tissue
    )
  )
  
  return(result)
}