
## ============================================================================
## NO3b v3: DESeq2-only from saved pseudobulk matrix
##
## Purpose:
##   Recalculate DEG for syn21788402 EC/SFG using saved pseudobulk matrix.
##   Avoid Seurat/dplyr/ggplot2/S7 conflicts.
##
## Important fix:
##   Some cell types work with DESeq2 default/parametric.
##   Some cell types fail under parametric but work with mean.
##   Therefore, this version automatically tries:
##       fitType = "parametric" -> "mean" -> "local"
##
## Output:
##   DEG_*_AD_vs_Control_allgenes.csv
##   INTERMEDIATE_plotdata_*_UBL3_*.csv
##   NO3b_v3_DESeq2_only_summary_EC_SFG.csv
## ============================================================================

rm(list = ls())
gc()

Sys.setenv(LANG = "en")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(DESeq2)
})

dataset_tag <- "syn21788402"

regions <- list(
  EC  = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/NO3_UBL3_Boxplots_DESeq2_byDonor_EC",
  SFG = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/NO3_UBL3_Boxplots_DESeq2_byDonor_SFG"
)

summary_fp <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/NO3b_v3_DESeq2_only_summary_EC_SFG.csv"

safe_n <- function(tab, nm) {
  if (nm %in% names(tab)) as.integer(tab[[nm]]) else 0L
}

run_deseq_with_fallback <- function(y, cd) {
  
  fit_types <- c("parametric", "mean", "local")
  
  last_error <- NA_character_
  
  for (ft in fit_types) {
    
    ans <- tryCatch({
      
      dds <- DESeqDataSetFromMatrix(
        countData = y,
        colData   = cd,
        design    = ~ group
      )
      
      dds <- DESeq(
        dds,
        quiet   = TRUE,
        sfType  = "poscounts",
        fitType = ft
      )
      
      res <- results(
        dds,
        contrast             = c("group", "AD", "Control"),
        independentFiltering = FALSE,
        cooksCutoff          = FALSE
      )
      
      list(
        ok = TRUE,
        fitType_used = ft,
        dds = dds,
        res = res,
        message = NA_character_
      )
      
    }, error = function(e) {
      
      list(
        ok = FALSE,
        fitType_used = ft,
        dds = NULL,
        res = NULL,
        message = conditionMessage(e)
      )
      
    })
    
    if (ans$ok) {
      return(ans)
    } else {
      last_error <- paste0("fitType=", ft, ": ", ans$message)
    }
  }
  
  list(
    ok = FALSE,
    fitType_used = NA_character_,
    dds = NULL,
    res = NULL,
    message = last_error
  )
}

all_summary <- data.frame()

for (region_tag in names(regions)) {
  
  out_dir <- regions[[region_tag]]
  
  pb_fp <- file.path(out_dir, "INTERMEDIATE_pseudobulk_matrix_celltype6_donor.rds")
  cd_fp <- file.path(out_dir, "INTERMEDIATE_pseudobulk_coldata_celltype6_donor_group.csv")
  
  cat("\n############################################################\n")
  cat("Region:", region_tag, "\n")
  
  if (!file.exists(pb_fp)) {
    cat("  !! Missing pseudobulk matrix:\n  ", pb_fp, "\n")
    
    all_summary <- rbind(
      all_summary,
      data.frame(
        region = region_tag,
        celltype6 = NA_character_,
        Ctrl = NA_integer_,
        AD = NA_integer_,
        status = "MISSING_pseudobulk_matrix",
        fitType_used = NA_character_,
        UBL3_log2FC = NA_real_,
        UBL3_padj = NA_real_,
        DEG_file = NA_character_,
        plotdata_file = NA_character_,
        error = pb_fp
      )
    )
    
    next
  }
  
  if (!file.exists(cd_fp)) {
    cat("  !! Missing pseudobulk coldata:\n  ", cd_fp, "\n")
    
    all_summary <- rbind(
      all_summary,
      data.frame(
        region = region_tag,
        celltype6 = NA_character_,
        Ctrl = NA_integer_,
        AD = NA_integer_,
        status = "MISSING_pseudobulk_coldata",
        fitType_used = NA_character_,
        UBL3_log2FC = NA_real_,
        UBL3_padj = NA_real_,
        DEG_file = NA_character_,
        plotdata_file = NA_character_,
        error = cd_fp
      )
    )
    
    next
  }
  
  pb  <- readRDS(pb_fp)
  pbm <- read.csv(cd_fp, stringsAsFactors = FALSE)
  
  required_cols <- c("key", "donor", "group", "celltype6")
  missing_cols <- setdiff(required_cols, colnames(pbm))
  
  if (length(missing_cols) > 0) {
    stop("Missing columns in coldata: ", paste(missing_cols, collapse = ", "))
  }
  
  pbm <- pbm[pbm$key %in% colnames(pb), , drop = FALSE]
  
  gene_row <- if ("UBL3" %in% rownames(pb)) {
    "UBL3"
  } else if ("ENSG00000122042" %in% rownames(pb)) {
    "ENSG00000122042"
  } else {
    NA_character_
  }
  
  if (is.na(gene_row)) {
    cat("  !! UBL3 not found in pseudobulk rownames. Skip.\n")
    
    all_summary <- rbind(
      all_summary,
      data.frame(
        region = region_tag,
        celltype6 = NA_character_,
        Ctrl = NA_integer_,
        AD = NA_integer_,
        status = "UBL3_not_found",
        fitType_used = NA_character_,
        UBL3_log2FC = NA_real_,
        UBL3_padj = NA_real_,
        DEG_file = NA_character_,
        plotdata_file = NA_character_,
        error = "UBL3 not found in rownames(pb)"
      )
    )
    
    next
  }
  
  dm <- unique(pbm[, c("donor", "group")])
  
  cat("  Donor group table:\n")
  print(table(dm$group))
  
  for (ct in sort(unique(pbm$celltype6))) {
    
    pbm_ct <- pbm[pbm$celltype6 == ct, , drop = FALSE]
    
    cols <- pbm_ct$key
    grp_ct <- factor(pbm_ct$group, levels = c("Control", "AD"))
    tab <- table(grp_ct)
    
    n_ctrl <- safe_n(tab, "Control")
    n_ad   <- safe_n(tab, "AD")
    
    deg_fp <- file.path(
      out_dir,
      paste0(
        "DEG_",
        dataset_tag, "_",
        region_tag, "_",
        ct,
        "_AD_vs_Control_allgenes.csv"
      )
    )
    
    plotdata_fp <- file.path(
      out_dir,
      paste0(
        "INTERMEDIATE_plotdata_",
        dataset_tag, "_",
        region_tag, "_UBL3_",
        ct,
        ".csv"
      )
    )
    

    if (file.exists(deg_fp)) file.remove(deg_fp)
    if (file.exists(plotdata_fp)) file.remove(plotdata_fp)
    
    if (length(cols) < 4 || n_ctrl < 2 || n_ad < 2) {
      
      cat(sprintf(
        "  Skip %-20s (Ctrl=%d, AD=%d)\n",
        ct, n_ctrl, n_ad
      ))
      
      all_summary <- rbind(
        all_summary,
        data.frame(
          region = region_tag,
          celltype6 = ct,
          Ctrl = n_ctrl,
          AD = n_ad,
          status = "SKIP_low_n",
          fitType_used = NA_character_,
          UBL3_log2FC = NA_real_,
          UBL3_padj = NA_real_,
          DEG_file = deg_fp,
          plotdata_file = plotdata_fp,
          error = NA_character_
        )
      )
      
      next
    }
    
    ymat <- as.matrix(pb[, cols, drop = FALSE])
    
    y <- matrix(
      as.integer(round(as.numeric(ymat))),
      nrow = nrow(ymat),
      ncol = ncol(ymat),
      dimnames = dimnames(ymat)
    )
    
    keep <- rowSums(y >= 10) >= min(n_ctrl, n_ad)
    keep[rownames(y) == gene_row] <- TRUE
    y <- y[keep, , drop = FALSE]
    
    cd <- data.frame(
      group = grp_ct,
      row.names = cols
    )
    
    ans <- run_deseq_with_fallback(y, cd)
    
    if (!ans$ok) {
      
      cat(sprintf(
        "  X DESeq failed for %-20s (Ctrl=%d, AD=%d): %s\n",
        ct, n_ctrl, n_ad, ans$message
      ))
      
      all_summary <- rbind(
        all_summary,
        data.frame(
          region = region_tag,
          celltype6 = ct,
          Ctrl = n_ctrl,
          AD = n_ad,
          status = "FAILED",
          fitType_used = NA_character_,
          UBL3_log2FC = NA_real_,
          UBL3_padj = NA_real_,
          DEG_file = deg_fp,
          plotdata_file = plotdata_fp,
          error = ans$message
        )
      )
      
      next
    }
    
    dds <- ans$dds
    res <- ans$res
    
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    
    res_df <- res_df[, c(
      "gene",
      "log2FoldChange",
      "padj",
      "pvalue",
      "baseMean"
    )]
    
    write.csv(res_df, deg_fp, row.names = FALSE)
    
    norm <- counts(dds, normalized = TRUE)
    
    dfp <- data.frame(
      donor = pbm_ct$donor,
      group = factor(pbm_ct$group, levels = c("AD", "Control")),
      value = as.numeric(norm[gene_row, cols]),
      stringsAsFactors = FALSE
    )
    
    write.csv(dfp, plotdata_fp, row.names = FALSE)
    
    u_fc <- as.numeric(res[gene_row, "log2FoldChange"])
    u_p  <- as.numeric(res[gene_row, "padj"])
    
    cat(sprintf(
      "  OK %-20s (Ctrl=%d, AD=%d)  fitType=%s  UBL3 log2FC=%.3f padj=%s\n",
      ct,
      n_ctrl,
      n_ad,
      ans$fitType_used,
      u_fc,
      ifelse(is.na(u_p), "NA", format(u_p, digits = 3, scientific = TRUE))
    ))
    
    all_summary <- rbind(
      all_summary,
      data.frame(
        region = region_tag,
        celltype6 = ct,
        Ctrl = n_ctrl,
        AD = n_ad,
        status = "OK",
        fitType_used = ans$fitType_used,
        UBL3_log2FC = u_fc,
        UBL3_padj = u_p,
        DEG_file = deg_fp,
        plotdata_file = plotdata_fp,
        error = NA_character_
      )
    )
  }
}

write.csv(all_summary, summary_fp, row.names = FALSE)

cat("\n============================================================\n")
cat("NO3b v3 finished.\n")
cat("Summary saved to:\n")
cat(summary_fp, "\n")
cat("============================================================\n")

cat("\nStatus table:\n")
print(table(all_summary$region, all_summary$status, useNA = "ifany"))

cat("\nFitType used:\n")
print(all_summary[, c("region", "celltype6", "status", "fitType_used", "UBL3_log2FC", "UBL3_padj")])


expected <- expand.grid(
  region = names(regions),
  celltype6 = c(
    "Astrocytes",
    "Endothelial",
    "Excitatory neurons",
    "Inhibitory neurons",
    "Microglia",
    "Oligodendrocytes"
  ),
  stringsAsFactors = FALSE
)

expected$file <- mapply(
  function(region, ct) {
    file.path(
      regions[[region]],
      paste0(
        "DEG_",
        dataset_tag, "_",
        region, "_",
        ct,
        "_AD_vs_Control_allgenes.csv"
      )
    )
  },
  expected$region,
  expected$celltype6
)

expected$exists <- file.exists(expected$file)

cat("\nDEG file existence check:\n")
print(expected[, c("region", "celltype6", "exists")])

missing_expected <- expected[!expected$exists, , drop = FALSE]

if (nrow(missing_expected) == 0) {

} else {

  print(missing_expected[, c("region", "celltype6", "file")])
}




############################################################

# build_SuppData2_DEG_NBDV1senseirevised.R
#






#

#   GSE157827_AD_vs_Control
#   GSE174367_AD_vs_Control
#   syn52082747_AD_vs_Control
#   syn21788402_EC_AD_vs_Control
#   syn21788402_SFG_AD_vs_Control
#   syn52082747_FTD_vs_Control
#   syn52082747_PSP_vs_Control
############################################################

###############################

###############################
rm(list = ls())
options(scipen = 999)
options(stringsAsFactors = FALSE)

###############################

###############################
need_cran_pkgs <- c("openxlsx")
need_bioc_pkgs <- c("AnnotationDbi", "org.Hs.eg.db")

for (pkg in need_cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in need_bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(openxlsx)
library(AnnotationDbi)
library(org.Hs.eg.db)

###############################

###############################
output_dir <- "D:/RNA/supptable/supply4/Results/NBDV1senseirevised"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

per_sheet_csv_dir <- file.path(output_dir, "per_sheet_csv")
if (!dir.exists(per_sheet_csv_dir)) dir.create(per_sheet_csv_dir, recursive = TRUE)

###############################

###############################
dataset_cfg <- data.frame(
  sheet_prefix = c(
    "GSE157827",
    "GSE174367",
    "syn21788402_EC",
    "syn21788402_SFG",
    "syn52082747"
  ),
  Dataset = c(
    "GSE157827",
    "GSE174367",
    "syn21788402",
    "syn21788402",
    "syn52082747"
  ),
  Region = c(
    "Middle frontal gyrus",
    "Prefrontal cortex",
    "Entorhinal cortex",
    "Superior frontal gyrus",
    "Primary visual cortex (V1)"
  ),
  result_dir = c(
    "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results/NOxx_GSE157827_UBL3_Boxplots_DESeq2_byDonor_ADvsControl",
    "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/NO3_GSE174367_UBL3_Boxplots_DESeq2_byDonor_ADvsControl",
    "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/NO3_UBL3_Boxplots_DESeq2_byDonor_EC",
    "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/NO3_UBL3_Boxplots_DESeq2_byDonor_SFG",
    "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO4/NO4_05_boxplots_deseq2_byDonor_6panel"
  ),
  qc_file = c(
    "QC_donor_counts_by_celltype6_by_group.csv",
    "QC_donor_counts_by_celltype6_by_group.csv",
    "QC_donor_counts_by_celltype6_by_group.csv",
    "QC_donor_counts_by_celltype6_by_group.csv",
    "QC_donor_counts_by_celltype6_by_group4.csv"
  ),
  stringsAsFactors = FALSE
)

###############################

###############################

clean_names <- function(x) {
  x <- trimws(x)
  x <- gsub("\ufeff", "", x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x <- tolower(x)
  return(x)
}

safe_read_csv <- function(file) {
  try_list <- list(
    list(fileEncoding = "UTF-8-BOM"),
    list(fileEncoding = "UTF-8"),
    list(fileEncoding = "GB18030"),
    list()
  )
  for (opt in try_list) {
    tmp <- try(
      do.call(read.csv, c(list(file = file, stringsAsFactors = FALSE,
                               check.names = FALSE), opt)),
      silent = TRUE
    )
    if (!inherits(tmp, "try-error")) return(tmp)
  }

}

standardize_celltype <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  x0[is.na(x0)] <- ""
  x0 <- gsub("_", " ", x0)
  x0 <- gsub("-", " ", x0)
  x0 <- gsub("\\s+", " ", x0)
  out <- rep(NA_character_, length(x0))
  out[grepl("astro", x0)] <- "Astrocytes"
  out[grepl("endo|endothelial", x0)] <- "Endothelial"
  out[grepl("excit", x0)] <- "Excitatory neurons"
  out[grepl("inhib", x0)] <- "Inhibitory neurons"
  out[grepl("microgl", x0)] <- "Microglia"
  out[grepl("oligo", x0)] <- "Oligodendrocytes"
  return(out)
}

standardize_group <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y <- gsub("[[:space:]]+", "", y)
  y[y %in% c("CTRL", "HC")] <- "CONTROL"

  return(y)
}

extract_comparison_from_filename <- function(file) {
  nm <- basename(file)
  possible_comparisons <- c("AD_vs_Control", "FTD_vs_Control",
                            "PSP_vs_Control", "CTE_vs_Control")
  hits <- possible_comparisons[
    sapply(possible_comparisons, function(cc) grepl(cc, nm, ignore.case = TRUE))
  ]
  if (length(hits) == 1) return(hits[1])

  if (length(hits) > 1) {

    return(hits[1])
  }
}

extract_celltype_from_filename <- function(file) {
  nm <- basename(file)
  nm2 <- sub("\\.csv$", "", nm, ignore.case = TRUE)
  ct <- standardize_celltype(nm2)
  if (length(unique(ct[!is.na(ct)])) == 0) {

  }
  ct_unique <- unique(ct[!is.na(ct)])
  return(ct_unique[1])
}

strip_ensembl_version <- function(x) {
  x <- as.character(x); x <- sub("\\..*$", "", x); return(x)
}

annotate_gene_columns <- function(gene_vec) {
  raw_gene <- trimws(as.character(gene_vec))
  raw_gene[raw_gene == ""] <- NA_character_

  ensembl_id  <- rep(NA_character_, length(raw_gene))
  gene_symbol <- rep(NA_character_, length(raw_gene))

  raw_gene_upper <- toupper(raw_gene)
  is_ensembl <- !is.na(raw_gene_upper) &
    grepl("^ENSG[0-9]+(\\.[0-9]+)?$", raw_gene_upper)


  if (any(is_ensembl)) {
    ens_clean <- strip_ensembl_version(raw_gene_upper[is_ensembl])
    ensembl_id[is_ensembl] <- ens_clean
    ens_keys <- unique(ens_clean[!is.na(ens_clean)])
    if (length(ens_keys) > 0) {
      symbol_map <- AnnotationDbi::mapIds(
        x = org.Hs.eg.db, keys = ens_keys,
        column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"
      )
      gene_symbol[is_ensembl] <- unname(symbol_map[ensembl_id[is_ensembl]])
    }
  }


  idx_non_ensembl <- which(!is_ensembl & !is.na(raw_gene))
  if (length(idx_non_ensembl) > 0) {
    gene_symbol[idx_non_ensembl] <- raw_gene[idx_non_ensembl]
    sym_keys <- unique(raw_gene[idx_non_ensembl])
    symbol_to_ens <- try(
      AnnotationDbi::mapIds(x = org.Hs.eg.db, keys = sym_keys,
                             column = "ENSEMBL", keytype = "SYMBOL",
                             multiVals = "first"),
      silent = TRUE
    )
    if (!inherits(symbol_to_ens, "try-error")) {
      ensembl_id[idx_non_ensembl] <- unname(symbol_to_ens[raw_gene[idx_non_ensembl]])
    }
  }

  data.frame(
    Ensembl_Gene_ID = ensembl_id,
    Gene_symbol     = gene_symbol,
    stringsAsFactors = FALSE
  )
}

read_deg_table <- function(file) {
  deg <- safe_read_csv(file)


  original_colnames <- names(deg)
  names(deg) <- clean_names(names(deg))

  gene_col <- which(names(deg) %in% c(
    "gene", "genes", "geneid", "gene_id",
    "symbol", "ensembl", "ensembl_id", "ensembl_gene_id"))
  if (length(gene_col) == 0) {
    first_vals <- as.character(deg[[1]][1:min(50, nrow(deg))])
    if (any(grepl("^ENSG|^[A-Za-z]", first_vals))) gene_col <- 1


  }
  gene_col <- gene_col[1]

  lfc_col <- which(grepl("^log2foldchange$|^log2foldc$|^log2fc$|log2fold", names(deg)))

  lfc_col <- lfc_col[1]

  pvalue_col <- which(grepl("^pvalue$|^p_value$|^pval$|^p_val$", names(deg)))

  pvalue_col <- pvalue_col[1]

  padj_col <- which(grepl("^padj$|^adjp$|^fdr$|adjusted", names(deg)))

  padj_col <- padj_col[1]

  basemean_col <- which(grepl("^basemean$|^base_mean$", names(deg)))

  basemean_col <- basemean_col[1]

  data.frame(
    Raw_Gene = as.character(deg[[gene_col]]),
    baseMean = suppressWarnings(as.numeric(as.character(deg[[basemean_col]]))),
    log2FoldChange = suppressWarnings(as.numeric(as.character(deg[[lfc_col]]))),
    pvalue = suppressWarnings(as.numeric(as.character(deg[[pvalue_col]]))),
    padj   = suppressWarnings(as.numeric(as.character(deg[[padj_col]]))),
    stringsAsFactors = FALSE
  )
}

read_qc_donor_table <- function(qc_file) {
  qc <- safe_read_csv(qc_file)

  names(qc) <- clean_names(names(qc))


  drop_cols <- c()
  for (j in seq_along(qc)) {
    nm <- names(qc)[j]
    if (nm %in% c("x", "x1", "row_names", "row_names_", "row.names",
                  "unnamed_0", "unnamed_1", "1")) {
      tmp_num <- suppressWarnings(as.numeric(as.character(qc[[j]])))
      if (sum(!is.na(tmp_num)) == nrow(qc)) drop_cols <- c(drop_cols, j)
    }
  }
  if (length(drop_cols) > 0) qc <- qc[, -drop_cols, drop = FALSE]


  char_cols <- which(sapply(qc, function(v) is.character(v) || is.factor(v)))
  celltype_col <- NA
  for (j in char_cols) {
    vals <- tolower(as.character(qc[[j]]))
    if (any(grepl("astro|endo|excit|inhib|microgl|oligo", vals))) {
      celltype_col <- j; break
    }
  }



  wide_group_cols <- names(qc)[
    grepl("^(ad|control|cte|ftd|psp|pid|ctrl|hc)$", names(qc))
  ]

  if (length(wide_group_cols) > 0) {
    out_list <- lapply(wide_group_cols, function(grp_nm) {
      data.frame(
        Cell_type = standardize_celltype(as.character(qc[[celltype_col]])),
        Group_raw = grp_nm,
        n_donors = suppressWarnings(as.integer(as.character(qc[[grp_nm]]))),
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, out_list)
  } else {

    group_col <- NA
    for (j in char_cols) {
      if (j == celltype_col) next
      vals <- tolower(as.character(qc[[j]]))
      if (any(grepl("control|ad|cte|ftd|psp|pid|ctrl|hc", vals))) {
        group_col <- j; break
      }
    }


    score <- rep(0, ncol(qc))
    for (j in seq_along(qc)) {
      if (j %in% c(celltype_col, group_col)) next
      tmp_num <- suppressWarnings(as.numeric(as.character(qc[[j]])))
      if (sum(!is.na(tmp_num)) == nrow(qc)) score[j] <- score[j] + 1
      if (grepl("donor|count|n_donor|n_donors|n$", names(qc)[j])) score[j] <- score[j] + 2
    }

    count_col <- which.max(score)

    out <- data.frame(
      Cell_type = standardize_celltype(as.character(qc[[celltype_col]])),
      Group_raw = as.character(qc[[group_col]]),
      n_donors  = suppressWarnings(as.integer(as.character(qc[[count_col]]))),
      stringsAsFactors = FALSE
    )
  }

  out$Group_std <- standardize_group(out$Group_raw)
  out <- out[!is.na(out$Cell_type) & !is.na(out$n_donors), , drop = FALSE]
  out
}

pick_one_count <- function(v, where_text = "") {
  v <- unique(v[!is.na(v)])
  if (length(v) == 0) return(NA_integer_)
  if (length(v) > 1) {


  }
  as.integer(v[1])
}

get_donor_counts <- function(qc_df, comparison, cell_type, qc_file = "") {
  parts <- strsplit(comparison, "_vs_", fixed = TRUE)[[1]]
  disease_label <- standardize_group(parts[1])
  control_label <- standardize_group(parts[2])
  sub_df <- qc_df[qc_df$Cell_type == cell_type, , drop = FALSE]

  n_disease <- pick_one_count(
    sub_df$n_donors[sub_df$Group_std == disease_label],
    where_text = paste(qc_file, "|", cell_type, "|", disease_label)
  )
  n_control <- pick_one_count(
    sub_df$n_donors[sub_df$Group_std == control_label],
    where_text = paste(qc_file, "|", cell_type, "|", control_label)
  )
  if (is.na(n_disease) || is.na(n_control)) {

            "\ncomparison = ", comparison, "\ncell_type = ", cell_type,
            "\nn_disease = ", n_disease, " n_control = ", n_control)
  }
  c(disease = n_disease, control = n_control)
}

make_sheet_name <- function(prefix, comparison) {
  x <- paste(prefix, comparison, sep = "_")
  x <- gsub("[:\\\\/?*\\[\\]]", "_", x)
  substr(x, 1, 31)
}

###############################

###############################
celltype_order <- c(
  "Astrocytes", "Endothelial", "Excitatory neurons",
  "Inhibitory neurons", "Microglia", "Oligodendrocytes"
)
comparison_order <- c(
  "AD_vs_Control", "FTD_vs_Control",
  "PSP_vs_Control", "CTE_vs_Control"
)

###############################

###############################
sheet_data_list <- list()
summary_list <- list()

for (i in seq_len(nrow(dataset_cfg))) {
  cfg <- dataset_cfg[i, ]
  cat("\n====================================================\n")


  cat("====================================================\n")



  deg_files <- list.files(path = cfg$result_dir, pattern = "^DEG_.*\\.csv$",
                          full.names = TRUE)



  qc_path <- file.path(cfg$result_dir, cfg$qc_file)

  qc_df <- read_qc_donor_table(qc_path)


  all_comparisons_found <- unique(
    vapply(deg_files, extract_comparison_from_filename, character(1))
  )
  all_comparisons_found <- comparison_order[
    comparison_order %in% all_comparisons_found
  ]


  for (comp in all_comparisons_found) {
    cat("\n---- comparison:", comp, "----\n")

    files_this_comp <- deg_files[
      vapply(deg_files, extract_comparison_from_filename, character(1)) == comp
    ]
    celltypes_this_comp <- unique(
      vapply(files_this_comp, extract_celltype_from_filename, character(1))
    )
    missing_ct <- setdiff(celltype_order, celltypes_this_comp)
    if (length(missing_ct) > 0) {

              "\nsheet: ", make_sheet_name(cfg$sheet_prefix, comp))
    }

    comp_table_list <- list()

    for (f in files_this_comp) {

      cell_type_now <- extract_celltype_from_filename(f)
      deg_std <- read_deg_table(f)

      gene_anno <- annotate_gene_columns(deg_std$Raw_Gene)
      deg_std$Ensembl_Gene_ID <- gene_anno$Ensembl_Gene_ID
      deg_std$Gene_symbol     <- gene_anno$Gene_symbol

      donor_counts <- get_donor_counts(qc_df, comp, cell_type_now, qc_path)

      deg_std$Dataset <- cfg$Dataset
      deg_std$Region <- cfg$Region
      deg_std$Comparison <- comp
      deg_std$Cell_type <- cell_type_now
      deg_std$n_disease_donors <- donor_counts["disease"]
      deg_std$n_control_donors <- donor_counts["control"]
      deg_std$method <- "donor-level pseudobulk DESeq2"

      deg_std <- deg_std[, c(
        "Dataset", "Region", "Comparison", "Cell_type",
        "Ensembl_Gene_ID", "Gene_symbol",
        "baseMean", "log2FoldChange", "pvalue", "padj",
        "n_disease_donors", "n_control_donors", "method"
      )]
      comp_table_list[[basename(f)]] <- deg_std

      mapped_symbol_n <- sum(!is.na(deg_std$Gene_symbol) & deg_std$Gene_symbol != "")
      mapped_symbol_pct <- round(100 * mapped_symbol_n / nrow(deg_std), 2)

      summary_list[[length(summary_list) + 1]] <- data.frame(
        sheet_name = make_sheet_name(cfg$sheet_prefix, comp),
        Dataset = cfg$Dataset, Region = cfg$Region,
        Comparison = comp, Cell_type = cell_type_now,
        n_genes = nrow(deg_std),
        n_symbol_mapped = mapped_symbol_n,
        pct_symbol_mapped = mapped_symbol_pct,
        n_disease_donors = donor_counts["disease"],
        n_control_donors = donor_counts["control"],
        deg_file = f, qc_file = qc_path,
        stringsAsFactors = FALSE
      )
    }

    comp_df <- do.call(rbind, comp_table_list)
    comp_df$Cell_type <- factor(comp_df$Cell_type, levels = celltype_order)
    comp_df <- comp_df[
      order(comp_df$Cell_type, comp_df$Ensembl_Gene_ID, comp_df$Gene_symbol),
      , drop = FALSE
    ]
    comp_df$Cell_type <- as.character(comp_df$Cell_type)

    sheet_name <- make_sheet_name(cfg$sheet_prefix, comp)
    sheet_data_list[[sheet_name]] <- comp_df

  }
}

###############################

###############################
all_df <- do.call(rbind, sheet_data_list)
summary_df <- do.call(rbind, summary_list)

all_csv_path <- file.path(output_dir,
  "Supplementary_Data_2_DESeq2_ALL_NBDV1senseirevised.csv")
summary_csv_path <- file.path(output_dir,
  "Supplementary_Data_2_build_summary_NBDV1senseirevised.csv")

write.csv(all_df, file = all_csv_path, row.names = FALSE, na = "")
write.csv(summary_df, file = summary_csv_path, row.names = FALSE, na = "")




###############################

###############################
for (sheet_nm in names(sheet_data_list)) {
  out_csv <- file.path(per_sheet_csv_dir, paste0(sheet_nm, ".csv"))
  write.csv(sheet_data_list[[sheet_nm]], file = out_csv,
            row.names = FALSE, na = "")
}


###############################
## 9. Excel workbook
###############################
xlsx_path <- file.path(output_dir,
  "Supplementary_Data_2_complete_omics_pseudobulk_DESeq2_NBDV1senseirevised.xlsx")

wb <- createWorkbook()
header_style <- createStyle(textDecoration = "bold",
                            halign = "center", valign = "center")

for (sheet_nm in names(sheet_data_list)) {
  addWorksheet(wb, sheet_nm)
  writeData(wb, sheet = sheet_nm, x = sheet_data_list[[sheet_nm]])
  addStyle(wb, sheet = sheet_nm, style = header_style,
           rows = 1, cols = 1:ncol(sheet_data_list[[sheet_nm]]),
           gridExpand = TRUE)
  freezePane(wb, sheet = sheet_nm, firstRow = TRUE)
  setColWidths(wb, sheet = sheet_nm,
               cols = 1:ncol(sheet_data_list[[sheet_nm]]),
               widths = c(14, 28, 18, 22, 20, 16, 12, 16, 12, 12, 16, 16, 28))
}

saveWorkbook(wb, file = xlsx_path, overwrite = TRUE)


###############################

###############################
cat("\n====================================================\n")
cat("🎉 Supplementary Data 2 (NBDV1senseirevised) 完成\n")
cat("====================================================\n")



cat("  1) ", xlsx_path, "\n", sep = "")
cat("  2) ", all_csv_path, "\n", sep = "")
cat("  3) ", summary_csv_path, "\n", sep = "")

cat("====================================================\n")
