############################################################
# 代码名称:
# build_SuppData2_DEG_NBDV1senseirevised.R
#
# 改动 vs 原 build_SuppData2_from_existing_DEG_tables.R:
#   - dataset_cfg 删 GSE155114 和 GSE261807 两行 (CTE 数据集排除)
#   - output_dir → D:/RNA/supptable/supply4/Results/NBDV1senseirevised
#   - 输出文件名加 _NBDV1senseirevised 后缀
#   - 其余 helper / DESeq2 表读取 / ENSG->symbol 注释 / QC donor 解析
#     / Excel 写入逻辑全部完整保留
#
# 输出 sheet 结构 (7 cohort, 与 SF3 / SuppTable SUMO 一致):
#   GSE157827_AD_vs_Control
#   GSE174367_AD_vs_Control
#   syn52082747_AD_vs_Control
#   syn21788402_EC_AD_vs_Control
#   syn21788402_SFG_AD_vs_Control
#   syn52082747_FTD_vs_Control
#   syn52082747_PSP_vs_Control
############################################################

###############################
## 0. 基础设置
###############################
rm(list = ls())
options(scipen = 999)
options(stringsAsFactors = FALSE)

###############################
## 1. 包加载
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
## 2. 输出目录 (NBDV1senseirevised)
###############################
output_dir <- "D:/RNA/supptable/supply4/Results/NBDV1senseirevised"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

per_sheet_csv_dir <- file.path(output_dir, "per_sheet_csv")
if (!dir.exists(per_sheet_csv_dir)) dir.create(per_sheet_csv_dir, recursive = TRUE)

###############################
## 3. ★ 数据集 cfg —— 5 行 (CTE 两条已删, 实际产 7 sheets)
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
## 4. 辅助函数 (与原版完全一致)
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
  stop("读取 csv 失败: ", file)
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
  y[y %in% c("PID")] <- "FTD"   # syn52082747 PiD 按 FTD 处理
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
  if (length(hits) == 0) stop("无法从文件名识别 comparison: ", nm)
  if (length(hits) > 1) {
    warning("识别到多个 comparison, 取第一个: ", nm)
    return(hits[1])
  }
}

extract_celltype_from_filename <- function(file) {
  nm <- basename(file)
  nm2 <- sub("\\.csv$", "", nm, ignore.case = TRUE)
  ct <- standardize_celltype(nm2)
  if (length(unique(ct[!is.na(ct)])) == 0) {
    stop("无法从文件名识别细胞类型: ", nm)
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

  ## 情况 1: ENSG → 反查 symbol
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

  ## 情况 2: 不是 ENSG → 当成 symbol, 反查 Ensembl
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
  if (nrow(deg) == 0 || ncol(deg) == 0) stop("DEG 文件为空: ", file)

  original_colnames <- names(deg)
  names(deg) <- clean_names(names(deg))

  gene_col <- which(names(deg) %in% c(
    "gene", "genes", "geneid", "gene_id",
    "symbol", "ensembl", "ensembl_id", "ensembl_gene_id"))
  if (length(gene_col) == 0) {
    first_vals <- as.character(deg[[1]][1:min(50, nrow(deg))])
    if (any(grepl("^ENSG|^[A-Za-z]", first_vals))) gene_col <- 1
    else stop("无法识别 gene 列: ", file,
              "\n原始列名: ", paste(original_colnames, collapse = ", "))
  }
  gene_col <- gene_col[1]

  lfc_col <- which(grepl("^log2foldchange$|^log2foldc$|^log2fc$|log2fold", names(deg)))
  if (length(lfc_col) == 0) stop("无法识别 log2FoldChange 列: ", file)
  lfc_col <- lfc_col[1]

  pvalue_col <- which(grepl("^pvalue$|^p_value$|^pval$|^p_val$", names(deg)))
  if (length(pvalue_col) == 0) stop("无法识别 pvalue 列: ", file)
  pvalue_col <- pvalue_col[1]

  padj_col <- which(grepl("^padj$|^adjp$|^fdr$|adjusted", names(deg)))
  if (length(padj_col) == 0) stop("无法识别 padj 列: ", file)
  padj_col <- padj_col[1]

  basemean_col <- which(grepl("^basemean$|^base_mean$", names(deg)))
  if (length(basemean_col) == 0) stop("无法识别 baseMean 列: ", file)
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
  if (nrow(qc) == 0 || ncol(qc) == 0) stop("QC donor 文件为空: ", qc_file)
  names(qc) <- clean_names(names(qc))

  ## 去掉行号列
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

  ## 识别 cell type 列
  char_cols <- which(sapply(qc, function(v) is.character(v) || is.factor(v)))
  celltype_col <- NA
  for (j in char_cols) {
    vals <- tolower(as.character(qc[[j]]))
    if (any(grepl("astro|endo|excit|inhib|microgl|oligo", vals))) {
      celltype_col <- j; break
    }
  }
  if (is.na(celltype_col)) stop("无法识别 QC donor 表 cell type 列: ", qc_file)

  ## 宽表识别
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
    ## 长表
    group_col <- NA
    for (j in char_cols) {
      if (j == celltype_col) next
      vals <- tolower(as.character(qc[[j]]))
      if (any(grepl("control|ad|cte|ftd|psp|pid|ctrl|hc", vals))) {
        group_col <- j; break
      }
    }
    if (is.na(group_col)) stop("无法识别 QC donor group 列: ", qc_file)

    score <- rep(0, ncol(qc))
    for (j in seq_along(qc)) {
      if (j %in% c(celltype_col, group_col)) next
      tmp_num <- suppressWarnings(as.numeric(as.character(qc[[j]])))
      if (sum(!is.na(tmp_num)) == nrow(qc)) score[j] <- score[j] + 1
      if (grepl("donor|count|n_donor|n_donors|n$", names(qc)[j])) score[j] <- score[j] + 2
    }
    if (max(score) == 0) stop("无法识别 QC donor 计数列: ", qc_file)
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
    warning("识别到多个 donor 计数, 取第一个。位置: ", where_text,
            "\n候选: ", paste(v, collapse = ", "))
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
    warning("QC donor 表中未完整找到 donor 数:\nQC = ", qc_file,
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
## 5. 固定顺序
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
## 6. 主循环
###############################
sheet_data_list <- list()
summary_list <- list()

for (i in seq_len(nrow(dataset_cfg))) {
  cfg <- dataset_cfg[i, ]
  cat("\n====================================================\n")
  cat("处理数据集:", cfg$sheet_prefix, "\n")
  cat("结果目录:", cfg$result_dir, "\n")
  cat("====================================================\n")

  if (!dir.exists(cfg$result_dir)) stop("结果目录不存在: ", cfg$result_dir)

  deg_files <- list.files(path = cfg$result_dir, pattern = "^DEG_.*\\.csv$",
                          full.names = TRUE)
  if (length(deg_files) == 0) stop("没有找到 DEG_*.csv: ", cfg$result_dir)
  cat("DEG 文件数量:", length(deg_files), "\n")

  qc_path <- file.path(cfg$result_dir, cfg$qc_file)
  if (!file.exists(qc_path)) stop("没有找到 QC donor 文件: ", qc_path)
  qc_df <- read_qc_donor_table(qc_path)
  cat("QC donor 读取成功:", qc_path, "\n")

  all_comparisons_found <- unique(
    vapply(deg_files, extract_comparison_from_filename, character(1))
  )
  all_comparisons_found <- comparison_order[
    comparison_order %in% all_comparisons_found
  ]
  cat("识别 comparison:", paste(all_comparisons_found, collapse = ", "), "\n")

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
      warning("缺少细胞类型: ", paste(missing_ct, collapse = ", "),
              "\nsheet: ", make_sheet_name(cfg$sheet_prefix, comp))
    }

    comp_table_list <- list()

    for (f in files_this_comp) {
      cat("读取:", basename(f), "\n")
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
    cat("✓ sheet:", sheet_name, "| rows:", nrow(comp_df), "\n")
  }
}

###############################
## 7. 输出总长表 + 摘要
###############################
all_df <- do.call(rbind, sheet_data_list)
summary_df <- do.call(rbind, summary_list)

all_csv_path <- file.path(output_dir,
  "Supplementary_Data_2_DESeq2_ALL_NBDV1senseirevised.csv")
summary_csv_path <- file.path(output_dir,
  "Supplementary_Data_2_build_summary_NBDV1senseirevised.csv")

write.csv(all_df, file = all_csv_path, row.names = FALSE, na = "")
write.csv(summary_df, file = summary_csv_path, row.names = FALSE, na = "")

cat("\n✓ 总长表:", all_csv_path, "\n")
cat("✓ 构建摘要:", summary_csv_path, "\n")

###############################
## 8. 每个 sheet 单独 csv 备份
###############################
for (sheet_nm in names(sheet_data_list)) {
  out_csv <- file.path(per_sheet_csv_dir, paste0(sheet_nm, ".csv"))
  write.csv(sheet_data_list[[sheet_nm]], file = out_csv,
            row.names = FALSE, na = "")
}
cat("✓ 单 sheet csv:", per_sheet_csv_dir, "\n")

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
cat("\n✓ Excel workbook:", xlsx_path, "\n")

###############################
## 10. 结束汇总
###############################
cat("\n====================================================\n")
cat("🎉 Supplementary Data 2 (NBDV1senseirevised) 完成\n")
cat("====================================================\n")
cat("Sheet 数量:", length(sheet_data_list), "(7 cohort, CTE 已排除)\n")
cat("总基因行数:", nrow(all_df), "\n\n")
cat("主要输出:\n")
cat("  1) ", xlsx_path, "\n", sep = "")
cat("  2) ", all_csv_path, "\n", sep = "")
cat("  3) ", summary_csv_path, "\n", sep = "")
cat("  4) ", per_sheet_csv_dir, "  (每 sheet csv 备份)\n", sep = "")
cat("====================================================\n")
