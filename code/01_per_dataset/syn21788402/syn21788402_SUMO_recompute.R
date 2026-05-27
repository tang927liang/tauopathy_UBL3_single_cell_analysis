










一、各个SUMO的值
## ============================================================================
## NO2 SUMO overlap histograms (FIXED 完整版)
## 改动（仅 2 处，其余逻辑不变）：
##  1) SFG 的 braak 分组：原来是 bs>0 -> AD（把 Braak2 错算进 AD，变成 3 vs 7）。
##     按你的设计改为：Braak0 -> Control(3人)，Braak6 -> AD(3人)，Braak2 -> 排除(NA)。
##     => SFG 变回正确的 3 vs 3。EC 走 meta，不受影响。
##  2) 画图包进 tryCatch：ggplot2 4.0/S7 在 SFG 偶发崩溃；stat_input 在画图前已存盘，
##     画图失败也不影响 SUMO 表数据。另加 suppressWarnings 读 .rds（伪警告）。
## ============================================================================
rm(list = ls()); gc()

# 关键：先指定 Seurat 所在的干净 library
clean_lib <- "D:/Rlibs/R45_seurat_clean"
.libPaths(c(clean_lib, .libPaths()))

Sys.setenv(LANG = "en")
options(stringsAsFactors = FALSE)
SEED <- 20251023
set.seed(SEED)

# 如果缺少绘图包，自动安装到 clean_lib
pkgs <- c("Seurat", "SeuratObject", "Matrix", "dplyr", "ggplot2", "ragg")

missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, lib = clean_lib, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ragg)
})

cat("Current .libPaths():\n")
print(.libPaths())

cat("\nPackage versions:\n")
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("SeuratObject:", as.character(packageVersion("SeuratObject")), "\n")
cat("Matrix:", as.character(packageVersion("Matrix")), "\n")
cat("dplyr:", as.character(packageVersion("dplyr")), "\n")
cat("ggplot2:", as.character(packageVersion("ggplot2")), "\n")
cat("ragg:", as.character(packageVersion("ragg")), "\n")

dataset_tag <- "syn21788402"
gene_list   <- c("SUMO1", "SUMO2", "SUMO3")
binwidth    <- 0.2
disease     <- "AD"

dataset_info <- list(
  EC = list(
    root_dir = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify",
    obj_fp   = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds",
    meta_fp  = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepP_syn21788402_matched_cells_meta.csv",
    group_mode = "meta"
  ),
  SFG = list(
    root_dir = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify",
    obj_fp   = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepH_syn21788402_SFG_obj_celltype6.rds",
    meta_fp  = NA,
    group_mode = "braak"
  )
)

get_counts_matrix_allcells <- function(obj, assay = "RNA", log_fp = NULL, out_dir = NULL) {
  a <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(a), error = function(e) character(0))
  counts_layers <- layers[grepl("^counts", layers)]
  if (length(counts_layers) == 0) {
    m2 <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "counts"), error = function(e) NULL)
    if (!is.null(m2)) return(m2)
    stop("ERROR: cannot get counts.")
  }
  mats <- list()
  for (ly in counts_layers) {
    m <- tryCatch(SeuratObject::LayerData(a, layer = ly), error = function(e) NULL)
    if (!is.null(m)) mats[[ly]] <- m
  }
  if (length(mats) == 0) stop("ERROR: LayerData failed.")
  ref_genes <- rownames(mats[[1]])
  for (k in names(mats)) {
    if (!identical(rownames(mats[[k]]), ref_genes)) {
      m0 <- mats[[k]]
      m_aligned <- Matrix::Matrix(0, nrow = length(ref_genes), ncol = ncol(m0), sparse = TRUE)
      rownames(m_aligned) <- ref_genes; colnames(m_aligned) <- colnames(m0)
      common <- intersect(ref_genes, rownames(m0))
      if (length(common) > 0) m_aligned[common, ] <- m0[common, , drop = FALSE]
      mats[[k]] <- m_aligned
    }
  }
  mat_all <- if (length(mats) == 1) mats[[1]] else Reduce(Matrix::cbind2, mats)
  if (!is.null(colnames(mat_all))) {
    dup <- duplicated(colnames(mat_all))
    if (any(dup)) mat_all <- mat_all[, !dup, drop = FALSE]
  }
  all_cells <- colnames(obj)
  miss_cells <- setdiff(all_cells, colnames(mat_all))
  if (length(miss_cells) > 0) {
    m_fill <- Matrix::Matrix(0, nrow = nrow(mat_all), ncol = length(miss_cells), sparse = TRUE)
    rownames(m_fill) <- rownames(mat_all); colnames(m_fill) <- miss_cells
    mat_all <- Matrix::cbind2(mat_all, m_fill)
  }
  j <- match(all_cells, colnames(mat_all))
  if (anyNA(j)) stop("ERROR: cells not aligned.")
  mat_all <- mat_all[, j, drop = FALSE]
  return(mat_all)
}

locate_gene_row <- function(counts_mat, gene_symbol) {
  rn <- rownames(counts_mat)
  if (gene_symbol %in% rn) return(gene_symbol)
  idx_ci <- which(toupper(rn) == toupper(gene_symbol))
  if (length(idx_ci) == 1) return(rn[idx_ci[1]])
  stop("ERROR: gene not found: ", gene_symbol)
}

assign_group_from_meta <- function(md0, meta_fp, out_dir, region_tag, gene_symbol) {
  if (!file.exists(meta_fp)) stop("ERROR: meta file missing: ", meta_fp)
  meta_raw <- read.csv(meta_fp, stringsAsFactors = FALSE)
  if (!all(c("sample", "group") %in% colnames(meta_raw))) stop("ERROR: meta missing sample/group.")
  ctrl_alias <- c("NC","Control","CTRL","Ctr","CTR","Normal","N","control","ctrl","ctr","nc","normal")
  grp_raw <- trimws(as.character(meta_raw$group))
  meta_raw$group4 <- ifelse(grp_raw %in% ctrl_alias, "Control", grp_raw)
  map_by_sample <- meta_raw %>%
    transmute(sample_join = trimws(as.character(sample)), group4 = group4) %>%
    filter(!is.na(sample_join), sample_join != "", group4 %in% c("AD", "Control")) %>%
    distinct()
  join_key <- c("SampleID", "sample", "orig.ident")[c("SampleID", "sample", "orig.ident") %in% colnames(md0)][1]
  if (is.na(join_key)) stop("ERROR: object meta missing SampleID/sample/orig.ident.")
  md0$sample_join_in_obj <- trimws(as.character(md0[[join_key]]))
  md1 <- dplyr::left_join(md0, map_by_sample, by = c("sample_join_in_obj" = "sample_join"))
  match_rate <- mean(!is.na(md1$group4))
  md1 <- md1[!is.na(md1$group4) & md1$group4 %in% c("AD", "Control"), , drop = FALSE]
  rownames(md1) <- md1$cell_id
  list(md = md1, join_key = join_key, match_rate = match_rate, mode = "meta")
}

## ★★★ 改动 1：braak 分组 —— Braak0->Control, Braak6->AD, 其它(含Braak2)->排除 ★★★
assign_group_from_braak <- function(md0, out_dir, region_tag, gene_symbol) {
  if (!("BraakStage" %in% colnames(md0))) stop("ERROR: no BraakStage column.")
  bs <- suppressWarnings(as.numeric(as.character(md0$BraakStage)))
  ## 你的设计：只比较两端，排除中间阶段（如 Braak2）
  group4 <- ifelse(is.na(bs), NA,
                   ifelse(bs == 0, "Control",
                          ifelse(bs == 6, "AD", NA)))
  md0$group4 <- group4
  md1 <- md0[!is.na(md0$group4) & md0$group4 %in% c("AD", "Control"), , drop = FALSE]
  rownames(md1) <- md1$cell_id
  write.csv(data.frame(BraakStage = bs, assigned_group = group4),
            file.path(out_dir, paste0("CHECK_", gene_symbol, "_braak_assignment_", region_tag, ".csv")),
            row.names = FALSE)
  list(md = md1, join_key = "BraakStage(0=Ctrl,6=AD)", match_rate = mean(!is.na(group4)), mode = "braak")
}

run_one_gene_one_region <- function(obj_fp, meta_fp, region_tag, out_root, gene_symbol, group_mode) {

  out_dir <- file.path(out_root, paste0("NO2_overlap_hist_density_SUMO_GSE157827style_", region_tag))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  log_fp <- file.path(out_dir, paste0("NO2_log_", region_tag, "_", gene_symbol, ".txt"))
  sink(log_fp); on.exit(sink(), add = TRUE)

  cat("==== START ====\n")
  cat("Region:", region_tag, "\nGene:", gene_symbol, "\ngroup_mode:", group_mode, "\n")
  cat("Time:", as.character(Sys.time()), "\nobj_fp :", obj_fp, "\nout_dir:", out_dir, "\n\n")

  if (!file.exists(obj_fp)) stop("ERROR: object file missing: ", obj_fp)
  obj <- suppressWarnings(readRDS(obj_fp))      ## 忽略伪警告
  DefaultAssay(obj) <- "RNA"

  md0 <- obj@meta.data
  md0$cell_id <- colnames(obj); rownames(md0) <- md0$cell_id
  if (!("celltype6" %in% colnames(md0))) stop("ERROR: no celltype6.")
  md0$celltype6 <- trimws(as.character(md0$celltype6))

  if ("PatientID" %in% colnames(md0)) { md0$donor <- trimws(as.character(md0$PatientID)); donor_source <- "PatientID"
  } else if ("sample" %in% colnames(md0)) { md0$donor <- trimws(as.character(md0$sample)); donor_source <- "sample"
  } else if ("SampleID" %in% colnames(md0)) { md0$donor <- trimws(as.character(md0$SampleID)); donor_source <- "SampleID"
  } else stop("ERROR: no donor column.")
  cat("donor_source =", donor_source, "\n\n")

  if (group_mode == "meta") {
    group_res <- assign_group_from_meta(md0, meta_fp, out_dir, region_tag, gene_symbol)
  } else if (group_mode == "braak") {
    group_res <- assign_group_from_braak(md0, out_dir, region_tag, gene_symbol)
  } else stop("ERROR: unknown group_mode.")

  md1 <- group_res$md; rownames(md1) <- md1$cell_id
  cat("group_assign_mode =", group_res$mode, "\njoin_key/source =", group_res$join_key,
      "\nmatch_rate =", sprintf("%.4f", group_res$match_rate), "\n\n")

  if (!all(c("AD","Control") %in% unique(md1$group4))) {
    cat("Group counts:\n"); print(table(md1$group4)); cat("\n")
    stop("ERROR: not both AD and Control present.")
  }
  cells_keep <- intersect(rownames(md1), colnames(obj))
  if (length(cells_keep) == 0) stop("ERROR: cells_keep=0.")
  md <- md1[cells_keep, , drop = FALSE]

  rna_counts_all <- get_counts_matrix_allcells(obj, "RNA", log_fp = log_fp, out_dir = out_dir)
  rna_counts <- rna_counts_all[, cells_keep, drop = FALSE]
  stopifnot(identical(colnames(rna_counts), rownames(md)))

  gene_row <- locate_gene_row(rna_counts, gene_symbol)
  lib_size <- Matrix::colSums(rna_counts)
  expr <- log1p((as.numeric(rna_counts[gene_row, , drop = TRUE]) / pmax(lib_size, 1)) * 1e4)

  df_all <- data.frame(expr = expr, donor = md$donor, group4 = md$group4,
                       celltype6 = md$celltype6, stringsAsFactors = FALSE)
  df0 <- df_all[df_all$expr > 0, , drop = FALSE]
  saveRDS(df0, file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_df0_exprGT0_", region_tag, ".rds")))
  write.csv(df0, file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_df0_exprGT0_", region_tag, ".csv")), row.names = FALSE)

  don_all <- unique(df0[, c("donor", "group4")])
  qc_don <- as.data.frame(table(don_all$group4), stringsAsFactors = FALSE)
  colnames(qc_don) <- c("group4", "n_donors")
  write.csv(qc_don, file.path(out_dir, paste0("QC_", gene_symbol, "_donors_by_group4_", region_tag, ".csv")), row.names = FALSE)

  plot_one_unit <- function(unit = c("donor", "cell")) {
    unit <- match.arg(unit)
    df2 <- df0 %>% filter(group4 %in% c(disease, "Control")) %>%
      mutate(group = ifelse(group4 == disease, disease, "Control"))
    don_pair <- unique(df2[, c("donor", "group")])
    n_dis <- sum(don_pair$group == disease); n_ctl <- sum(don_pair$group == "Control")
    lab_dis <- sprintf("%s\n(n = %d)", disease, n_dis)
    lab_ctl <- sprintf("Control\n(n = %d)", n_ctl)
    df2$group_lab <- factor(ifelse(df2$group == disease, lab_dis, lab_ctl), levels = c(lab_dis, lab_ctl))

    if (unit == "donor") {
      stat_input <- df2 %>% group_by(celltype6, donor, group_lab) %>%
        summarise(val = median(expr), .groups = "drop")
    } else {
      stat_input <- df2 %>% transmute(celltype6 = celltype6, group_lab = group_lab, val = expr)
    }
    ## stat_input 在画图前就存盘 —— SUMO 表只依赖它
    saveRDS(stat_input, file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_stat_input_AD_vs_Control_", unit, "_", region_tag, ".rds")))
    write.csv(stat_input, file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_stat_input_AD_vs_Control_", unit, "_", region_tag, ".csv")), row.names = FALSE)

    n_by_celltype <- stat_input %>% group_by(celltype6, group_lab) %>% summarise(n = n(), .groups = "drop")
    write.csv(n_by_celltype, file.path(out_dir, paste0("CHECK_", gene_symbol, "_n_by_celltype_AD_vs_Control_", unit, "_", region_tag, ".csv")), row.names = FALSE)

    stats <- stat_input %>% group_by(celltype6) %>%
      summarise(p_raw = tryCatch(wilcox.test(val ~ group_lab, exact = FALSE)$p.value,
                                 error = function(e) NA_real_), .groups = "drop")
    stats$padj <- p.adjust(stats$p_raw, method = "BH")
    stats$label <- sprintf("Mann-Whitney U\nPadj=%.2e", stats$padj)
    stats$x <- Inf; stats$y <- Inf
    write.csv(stats, file.path(out_dir, paste0("STATS_", gene_symbol, "_AD_vs_Control_", unit, "_", region_tag, ".csv")), row.names = FALSE)

    ## ★★★ 改动 2：画图包进 tryCatch（崩了也不影响上面已存的 stat_input / 表） ★★★
    tryCatch({
      fill_vals <- c("red", "blue"); names(fill_vals) <- c(lab_dis, lab_ctl)
      p <- ggplot(df2, aes(x = expr, y = after_stat(density), fill = group_lab)) +
        geom_histogram(binwidth = binwidth, alpha = 0.7, position = "identity", colour = NA) +
        facet_wrap(~celltype6, scales = "free_y") +
        scale_fill_manual(values = fill_vals, drop = FALSE) +
        geom_label(data = stats, inherit.aes = FALSE, aes(x = x, y = y, label = label),
                   hjust = 1.02, vjust = 1.02, size = 2.3, linewidth = 0, fill = "white", alpha = 0.7) +
        labs(title = paste0(gene_symbol, " expression per cell type (only expressed cells): ",
                            disease, " vs Control (", unit, "-level)"),
             x = paste0(gene_symbol, " log1p(CP10k)"), y = "Density", fill = "group_lab") +
        theme_bw() + theme(plot.margin = margin(10, 25, 10, 10)) + coord_cartesian(clip = "off")
      out_png <- file.path(out_dir, paste0(dataset_tag, "_", region_tag, "_", gene_symbol,
                           "_OverlapHistDensity_AD_vs_Control_", unit, ".png"))
      ragg::agg_png(out_png, width = 10, height = 6, units = "in", res = 1000, background = "white")
      print(p); dev.off()
      out_pdf <- file.path(out_dir, paste0(dataset_tag, "_", region_tag, "_", gene_symbol,
                           "_OverlapHistDensity_AD_vs_Control_", unit, ".pdf"))
      ggsave(out_pdf, p, width = 10, height = 6, device = cairo_pdf)
      cat("saved:", out_png, "\nsaved:", out_pdf, "\n")
    }, error = function(e) {
      try(while (!is.null(dev.list())) dev.off(), silent = TRUE)
      cat("PLOT failed for", unit, "(stat_input already saved, table NOT affected):",
          conditionMessage(e), "\n")
    })
  }

  plot_one_unit("cell")
  plot_one_unit("donor")

  cat("\n==== END ====\n")
  message("DONE region=", region_tag, " | gene=", gene_symbol)
}

for (region_tag in names(dataset_info)) {
  gc()
  for (gene_symbol in gene_list) {
    run_one_gene_one_region(
      obj_fp     = dataset_info[[region_tag]]$obj_fp,
      meta_fp    = dataset_info[[region_tag]]$meta_fp,
      region_tag = region_tag,
      out_root   = dataset_info[[region_tag]]$root_dir,
      gene_symbol= gene_symbol,
      group_mode = dataset_info[[region_tag]]$group_mode
    )
  }
}



































二、合并纳入的4个数据的SUMO123值





# ========================================================================
# 脚本名称: make_SuppTable_SUMO_NBDV1senseirevised.R  (FIXED)
#
# 目的:
#   汇总 7 个 disease/control 比较 (CTE 两条已删) 的 SUMO1/2/3 donor-level
#   comparator results, 输出 NBDV1senseirevised 版 Supplementary Table。
#
# 【本次唯一修改】read_donor_median_tbl() 里取"每 donor 中位数"那一列的方式：
#   原 guess_value_col(prefer=c("expr","median","value","log1p","cp10k")) 匹配不到
#   名为 `val` 的中位数列（"value"≠"val"），于是退到第一个数值列；而 syn21788402
#   的 `donor` 列是纯数字编号("1","2","8","9","10")会被 readr 当成数值列，结果
#   误取 donor 编号、算出 donor 编号的中位数（EC AD={8,9,10}→9，Control={1,2,3}→2）。
#   修复：排除 donor/id/计数类数值列后，优先取 val/median/expr/value，取不到再退第一列。
#   （powered 数据集不受影响：它们 donor 列非数字，本就不会被误选。）
# ========================================================================

suppressPackageStartupMessages({
  library(fs); library(readr); library(dplyr); library(purrr)
  library(stringr); library(stringi); library(tidyr); library(tibble)
  library(openxlsx); library(glue)
})

# -----------------------------
# 一、用户参数区
# -----------------------------
out_dir  <- "D:/RNA/supptable/Supp5SUMO/Results/NBDV1senseirevised"
out_xlsx <- path(out_dir,
  "Supplementary_Table_SUMO_donor_level_comparator_results_NBDV1senseirevised.xlsx")
dir_create(out_dir, recurse = TRUE)

## ★ 7 个比较 (CTE 的 H/I 两行已删)
comparison_cfg <- tribble(
  ~panel, ~dataset,       ~disease, ~region,                       ~comparison,       ~case_label, ~control_label, ~n_case_total, ~n_control_total, ~source_dir,
  "A",   "GSE157827",    "AD",     "Middle frontal gyrus",        "AD_vs_Control",   "AD",       "Control",      12,            9,                "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results",
  "B",   "GSE174367",    "AD",     "Prefrontal cortex",           "AD_vs_Control",   "AD",       "Control",      11,            7,                "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results",
  "C",   "syn52082747",  "AD",     "Primary visual cortex (V1)",  "AD_vs_Control",   "AD",       "Control",      10,            10,               "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results",
  "D",   "syn21788402",  "AD",     "Entorhinal cortex",           "AD_vs_Control",   "AD",       "Control",      3,             3,                "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/NO2_overlap_hist_density_SUMO_GSE157827style_EC",
  "E",   "syn21788402",  "AD",     "Superior frontal gyrus",      "AD_vs_Control",   "AD",       "Control",      3,             3,                "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/NO2_overlap_hist_density_SUMO_GSE157827style_SFG",
  "F",   "syn52082747",  "FTD",    "Primary visual cortex (V1)",  "FTD_vs_Control",  "FTD",      "Control",      9,             10,               "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results",
  "G",   "syn52082747",  "PSP",    "Primary visual cortex (V1)",  "PSP_vs_Control",  "PSP",      "Control",      11,            10,               "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results"
)

genes <- c("SUMO1", "SUMO2", "SUMO3")
celltype_order <- c(
  "Astrocytes", "Endothelial", "Excitatory neurons",
  "Inhibitory neurons", "Microglia", "Oligodendrocytes"
)

# -----------------------------
# 二、基础辅助函数 (与原版一致)
# -----------------------------
to_utf8_safe <- function(x) {
  x <- as.character(x)
  x <- iconv(x, from = "", to = "UTF-8", sub = "")
  x <- stringi::stri_enc_toutf8(x, is_unknown_8bit = TRUE, validate = TRUE)
  x <- gsub("[[:cntrl:]]", " ", x)
  x <- gsub("\\s+", " ", x)
  x[is.na(x)] <- ""
  x
}

sanitize_df_for_excel <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.factor(col)) col <- as.character(col)
    if (inherits(col, c("POSIXct", "POSIXt"))) col <- as.character(col)
    if (is.character(col)) col <- to_utf8_safe(col)
    col
  })
  names(df) <- to_utf8_safe(names(df))
  df
}

check_bad_encoding_column <- function(df) {
  for (nm in names(df)) {
    col <- df[[nm]]
    if (is.factor(col)) col <- as.character(col)
    if (is.character(col)) {
      test <- try(stringi::stri_length(col), silent = TRUE)
      if (inherits(test, "try-error")) {
        cat("BAD COLUMN:", nm, "\n"); print(utils::head(col, 20))
        stop(paste("Encoding problem in column:", nm))
      }
    }
  }
  invisible(TRUE)
}

resolve_existing_dir <- function(root) {
  cand <- unique(c(
    root,
    str_replace(root, fixed("sn_scRNA"), "sn_RNA"),
    str_replace(root, fixed("sn_RNA"), "sn_scRNA")
  ))
  hit <- cand[dir_exists(cand)]
  if (length(hit) >= 1) {
    if (normalizePath(hit[1], winslash = "/", mustWork = TRUE) !=
        normalizePath(root, winslash = "/", mustWork = FALSE)) {
      message("路径自动修正: ", root, "  -->  ", hit[1])
    }
    return(hit[1])
  }
  stop("目录不存在 (已尝试 sn_scRNA / sn_RNA 替换): ", root)
}

find_one_file <- function(root, patterns, recurse = TRUE) {
  root <- resolve_existing_dir(root)
  files <- dir_ls(root, recurse = recurse, type = "file")
  for (pat in patterns) {
    hit <- files[str_detect(path_file(files), regex(pat, ignore_case = TRUE))]
    if (length(hit) == 1) return(hit)
    if (length(hit) > 1) {
      hit2 <- hit[str_detect(hit, regex("SUMO", ignore_case = TRUE))]
      if (length(hit2) == 1) return(hit2)
      hit3 <- hit[str_detect(hit, regex("donorMWU|donorMedian|cellHist", ignore_case = TRUE))]
      if (length(hit3) == 1) return(hit3)
      stop(glue("目录 {root} 中模式 {pat} 命中多个文件:\n{paste(hit, collapse = '\n')}"))
    }
  }
  stop(glue("目录 {root} 中找不到目标。尝试的模式:\n{paste(patterns, collapse = '\n')}"))
}

clean_group_label <- function(x) {
  x %>% as.character() %>% str_replace_all("[\r\n]+", " ") %>%
    str_replace_all("\\s+", " ") %>% str_trim()
}

extract_group_key <- function(x) {
  x <- clean_group_label(x)
  x <- str_replace(x, "\\s*\\(n\\s*=.*$", "")
  str_trim(x)
}

extract_total_n_from_group_label <- function(x) {
  x <- clean_group_label(x)
  out <- str_match(x, "\\(n\\s*=\\s*([0-9]+)\\)")[, 2]
  suppressWarnings(as.integer(out))
}

guess_group_col <- function(df) {
  cand <- names(df)[str_detect(names(df),
    regex("group|group_lab|group4|condition|disease", ignore_case = TRUE))]
  if (length(cand) == 0) stop("无法识别 group 列。列名: ", paste(names(df), collapse = ", "))
  cand[1]
}

guess_celltype_col <- function(df) {
  cand <- names(df)[str_detect(names(df),
    regex("celltype|cell_type|celltype6", ignore_case = TRUE))]
  if (length(cand) == 0) stop("无法识别 cell type 列。列名: ", paste(names(df), collapse = ", "))
  cand[1]
}

guess_value_col <- function(df, prefer = c("expr","median","value","cp10k","log1p","n_","count")) {
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  if (length(num_cols) == 0) stop("数据中无数值列。")
  num_cols2 <- setdiff(num_cols, c("x", "y"))
  if (length(num_cols2) == 0) num_cols2 <- num_cols
  for (p in prefer) {
    hit <- num_cols2[str_detect(num_cols2, regex(p, ignore_case = TRUE))]
    if (length(hit) >= 1) return(hit[1])
  }
  num_cols2[1]
}

# ★FIX 专用：稳健地选出"每 donor 中位数"那一列（排除会被当成数值的 donor/id/计数列）
pick_donor_median_col <- function(df) {
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  num_cols <- setdiff(num_cols, c("x", "y"))
  # 排除 donor 编号 / id / 计数 / 索引类数值列
  bad <- grepl("donor|^id$|_id$|index|^n$|count|cells|rank|^row",
               num_cols, ignore.case = TRUE)
  cand <- num_cols[!bad]
  if (length(cand) == 0) cand <- num_cols   # 兜底
  # 在剩余列里优先取像"中位数 / 表达值"的列
  hit <- cand[grepl("val|median|expr|value|cp10k|log1p", cand, ignore.case = TRUE)]
  vc <- if (length(hit) >= 1) hit[1] else cand[1]
  if (is.na(vc) || length(vc) == 0)
    stop("找不到 donor 中位数值列 (期望 val / median_expr 之类的数值列)。列名: ",
         paste(names(df), collapse = ", "))
  vc
}

read_csv_safely <- function(file) readr::read_csv(file, show_col_types = FALSE, progress = FALSE)

map_case_control <- function(group_values, case_label, control_label) {
  g <- extract_group_key(group_values)
  case_when(
    str_to_upper(g) == str_to_upper(case_label)    ~ "case",
    str_to_upper(g) == str_to_upper(control_label) ~ "control",
    TRUE ~ NA_character_
  )
}

# -----------------------------
# 三、各类输入文件的读取器
# -----------------------------
read_stats_tbl <- function(source_dir, gene, comparison) {
  f <- find_one_file(source_dir, c(
    glue("^STATS_{gene}_{comparison}_MWU_donorMedian_BH.*\\.csv$"),
    glue("^STATS_{gene}_{comparison}_donor.*\\.csv$"),
    glue("^STATS_{gene}_MWU_donorMedian_BH.*\\.csv$")
  ), recurse = TRUE)
  df <- read_csv_safely(f)
  cell_col  <- guess_celltype_col(df)
  p_col     <- names(df)[str_detect(names(df), regex("p_raw|pvalue|p_value|^p$", ignore_case = TRUE))][1]
  padj_col  <- names(df)[str_detect(names(df), regex("padj|fdr|adj", ignore_case = TRUE))][1]
  label_col <- names(df)[str_detect(names(df), regex("label", ignore_case = TRUE))][1]
  if (is.na(p_col) || is.na(padj_col)) stop("统计文件无法识别 p_raw / padj 列: ", f)
  df %>% transmute(
    cell_type = .data[[cell_col]],
    p_raw     = .data[[p_col]],
    padj      = .data[[padj_col]],
    stat_label = if (!is.na(label_col)) .data[[label_col]] else NA_character_
  )
}

read_n_donors_expr_tbl <- function(source_dir, gene, comparison, case_label, control_label) {
  f <- find_one_file(source_dir, c(
    glue("^CHECK_{gene}_n_donors_by_celltype_{comparison}.*\\.csv$"),
    glue("^CHECK_{gene}_n_by_celltype_{comparison}_donor.*\\.csv$"),
    glue("^CHECK_{gene}_n_donors_by_celltype_.*\\.csv$"),
    glue("^CHECK_{gene}_n_by_celltype_.*donor.*\\.csv$")
  ), recurse = TRUE)
  df <- read_csv_safely(f)
  cell_col  <- guess_celltype_col(df)
  group_col <- guess_group_col(df)
  value_col <- guess_value_col(df, prefer = c("n_donors", "exprGT0", "count", "n_"))
  df %>% transmute(
    cell_type = .data[[cell_col]],
    group_raw = .data[[group_col]],
    group_key = map_case_control(.data[[group_col]], case_label, control_label),
    n_donors_exprGT0 = .data[[value_col]],
    total_n_from_label = extract_total_n_from_group_label(.data[[group_col]])
  ) %>% filter(!is.na(group_key))
}

read_n_cells_expr_tbl <- function(source_dir, gene, comparison, case_label, control_label) {
  source_dir <- resolve_existing_dir(source_dir)
  files <- dir_ls(source_dir, recurse = TRUE, type = "file")
  possible_patterns <- c(
    glue("^CHECK_{gene}_n_cells_by_celltype_{comparison}.*\\.csv$"),
    glue("^CHECK_{gene}_n_by_celltype_{comparison}_cell.*\\.csv$"),
    glue("^CHECK_{gene}_n_cells_by_celltype_.*\\.csv$"),
    glue("^CHECK_{gene}_n_by_celltype_.*cell.*\\.csv$")
  )
  for (pat in possible_patterns) {
    hit <- files[str_detect(path_file(files), regex(pat, ignore_case = TRUE))]
    if (length(hit) == 1) {
      df <- read_csv_safely(hit)
      cell_col  <- guess_celltype_col(df)
      group_col <- guess_group_col(df)
      value_col <- guess_value_col(df, prefer = c("n_cells", "exprGT0", "count", "n_"))
      return(df %>% transmute(
        cell_type = .data[[cell_col]],
        group_key = map_case_control(.data[[group_col]], case_label, control_label),
        n_cells_exprGT0 = .data[[value_col]]
      ) %>% filter(!is.na(group_key)))
    }
  }
  f <- find_one_file(source_dir, c(
    glue("^INTERMEDIATE_{gene}_df2_{comparison}_exprGT0.*\\.csv$"),
    glue("^INTERMEDIATE_{gene}_df0_exprGT0.*\\.csv$"),
    glue("^INTERMEDIATE_{gene}_df2_.*exprGT0.*\\.csv$")
  ), recurse = TRUE)
  df <- read_csv_safely(f)
  cell_col  <- guess_celltype_col(df)
  group_col <- guess_group_col(df)
  df %>% transmute(
    cell_type = .data[[cell_col]],
    group_key = map_case_control(.data[[group_col]], case_label, control_label)
  ) %>% filter(!is.na(group_key)) %>%
    count(cell_type, group_key, name = "n_cells_exprGT0")
}

read_donor_median_tbl <- function(source_dir, gene, comparison, case_label, control_label) {
  f <- find_one_file(source_dir, c(
    glue("^INTERMEDIATE_{gene}_stat_input_{comparison}_donor.*\\.csv$"),
    glue("^INTERMEDIATE_{gene}_stat_input_.*donorMedian.*\\.csv$"),
    glue("^INTERMEDIATE_{gene}_stat_input_.*donor.*\\.csv$")
  ), recurse = TRUE)
  df <- read_csv_safely(f)
  cell_col  <- guess_celltype_col(df)
  group_col <- guess_group_col(df)
  # ★FIX: 用稳健的列选择，避免把数字编号的 donor 列当成中位数值列
  value_col <- pick_donor_median_col(df)
  df %>% transmute(
    cell_type = .data[[cell_col]],
    group_key = map_case_control(.data[[group_col]], case_label, control_label),
    donor_median_value = .data[[value_col]]
  ) %>% filter(!is.na(group_key)) %>%
    group_by(cell_type, group_key) %>%
    summarise(donor_median_group = median(donor_median_value, na.rm = TRUE), .groups = "drop")
}

# -----------------------------
# 四、单个 comparison × gene 汇总 (与原版一致)
# -----------------------------
summarise_one_comparison_gene <- function(cfg_row, gene) {
  src             <- cfg_row$source_dir[[1]]
  disease         <- cfg_row$disease[[1]]
  dataset         <- cfg_row$dataset[[1]]
  region          <- cfg_row$region[[1]]
  comparison      <- cfg_row$comparison[[1]]
  case_label      <- cfg_row$case_label[[1]]
  control_label   <- cfg_row$control_label[[1]]
  panel           <- cfg_row$panel[[1]]
  n_case_total    <- cfg_row$n_case_total[[1]]
  n_control_total <- cfg_row$n_control_total[[1]]

  stats_tbl  <- read_stats_tbl(src, gene, comparison)
  ndonor_tbl <- read_n_donors_expr_tbl(src, gene, comparison, case_label, control_label)
  ncell_tbl  <- read_n_cells_expr_tbl(src, gene, comparison, case_label, control_label)
  dmed_tbl   <- read_donor_median_tbl(src, gene, comparison, case_label, control_label)

  ndonor_wide <- ndonor_tbl %>%
    select(cell_type, group_key, n_donors_exprGT0) %>%
    pivot_wider(names_from = group_key, values_from = n_donors_exprGT0, names_prefix = "n_") %>%
    rename(n_case_donors_exprGT0 = n_case, n_control_donors_exprGT0 = n_control)

  ncell_wide <- ncell_tbl %>%
    select(cell_type, group_key, n_cells_exprGT0) %>%
    pivot_wider(names_from = group_key, values_from = n_cells_exprGT0, names_prefix = "n_") %>%
    rename(n_case_cells_exprGT0 = n_case, n_control_cells_exprGT0 = n_control)

  dmed_wide <- dmed_tbl %>%
    pivot_wider(names_from = group_key, values_from = donor_median_group, names_prefix = "median_") %>%
    rename(donor_median_case = median_case, donor_median_control = median_control)

  full_join(stats_tbl, ndonor_wide, by = "cell_type") %>%
    full_join(ncell_wide, by = "cell_type") %>%
    full_join(dmed_wide, by = "cell_type") %>%
    mutate(
      panel = panel, disease = disease, dataset = dataset, region = region,
      comparison = comparison, gene = gene,
      case_label = case_label, control_label = control_label,
      n_case_total_donors = n_case_total, n_control_total_donors = n_control_total,
      delta_case_minus_control = donor_median_case - donor_median_control,
      endpoint            = "conditional expression among gene-positive cells (UMI > 0)",
      donor_summary_stat  = "median expr per donor within cell type",
      statistical_test    = "Mann-Whitney U / Wilcoxon rank-sum test",
      multiple_testing    = "Benjamini-Hochberg across six cell types within each cohort"
    ) %>%
    select(
      panel, disease, dataset, region, comparison, gene, cell_type,
      case_label, control_label,
      n_case_total_donors, n_control_total_donors,
      n_case_donors_exprGT0, n_control_donors_exprGT0,
      n_case_cells_exprGT0, n_control_cells_exprGT0,
      donor_median_case, donor_median_control, delta_case_minus_control,
      p_raw, padj, stat_label,
      endpoint, donor_summary_stat, statistical_test, multiple_testing
    )
}

# -----------------------------
# 五、汇总所有 (7) comparison × 3 gene = 21 组合
# -----------------------------
message("正在处理 7 个 comparison × 3 SUMO 基因 = 21 组合...")
all_tbl <- map_dfr(seq_len(nrow(comparison_cfg)), function(i) {
  cfg_row <- comparison_cfg[i, , drop = FALSE]
  map_dfr(genes, ~ summarise_one_comparison_gene(cfg_row, .x))
}) %>%
  mutate(
    cell_type = factor(cell_type, levels = celltype_order),
    gene      = factor(gene, levels = genes),
    panel     = factor(panel, levels = comparison_cfg$panel)
  ) %>%
  arrange(panel, gene, cell_type) %>%
  mutate(
    cell_type = as.character(cell_type),
    gene      = as.character(gene),
    panel     = as.character(panel)
  )

all_tbl <- sanitize_df_for_excel(all_tbl)
if ("stat_label" %in% names(all_tbl)) {
  all_tbl$stat_label <- iconv(all_tbl$stat_label, from = "", to = "ASCII//TRANSLIT", sub = "")
  all_tbl$stat_label[is.na(all_tbl$stat_label)] <- ""
}
all_tbl <- sanitize_df_for_excel(all_tbl)
check_bad_encoding_column(all_tbl)

# ★FIX 的健全性自检：donor 中位数应在 log1p(CP10k) 量级（约 0–5），若仍出现 >5 的
#   "整数常数"，多半还是取错了列，提前报警。
.bad_med <- with(all_tbl, donor_median_case[is.finite(donor_median_case) & donor_median_case > 5])
if (length(.bad_med) > 0) {
  warning("⚠ 检测到 donor_median_case > 5 的异常值（疑似又取到 donor 编号列），请核对：\n  ",
          paste(round(unique(.bad_med), 2), collapse = ", "))
}

message(sprintf("✅ 总表生成: %d 行 × %d 列", nrow(all_tbl), ncol(all_tbl)))

# -----------------------------
# 六、写 Excel 工作簿
# -----------------------------
wb <- createWorkbook()

## README ("nine" → "seven")
addWorksheet(wb, "README")
readme_lines <- tibble(
  Item = c(
    "Workbook purpose",
    "Table title",
    "Number of comparisons",
    "Cohort exclusions",
    "Primary endpoint",
    "Unit of inference",
    "Statistical test",
    "Multiple-testing correction",
    "n_case_total_donors / n_control_total_donors",
    "n_case_donors_exprGT0 / n_control_donors_exprGT0",
    "n_case_cells_exprGT0 / n_control_cells_exprGT0",
    "donor_median_case / donor_median_control",
    "delta_case_minus_control",
    "Source of this table"
  ),
  Description = c(
    "Supplementary Table for SUMO1/2/3 donor-level comparator results (NBDV1 senseirevised, CTE excluded).",
    "SUMO1/2/3 donor-level comparator results across seven disease-control comparisons.",
    "Seven (7) disease-vs-control comparisons spanning AD, FTD, and PSP cohorts (CTE cohorts excluded in this revision).",
    "GSE155114 (CTE) and GSE261807 (CTE) cohorts were excluded from this revised analysis to focus the manuscript on AD and other tauopathy cohorts.",
    "Conditional expression among gene-positive cells (UMI > 0).",
    "Donor.",
    "Mann-Whitney U / Wilcoxon rank-sum test on donor-wise median expression within cell type.",
    "Benjamini-Hochberg correction across six major cell types within each cohort.",
    "Total donor counts used for the comparison, recorded according to the current manuscript / finalized cohort definition.",
    "Number of donors contributing expr > 0 cells for the specified gene and cell type.",
    "Number of expr > 0 cells included in the conditional-expression analysis for the specified gene and cell type.",
    "Median of donor-wise median expression values for case / control groups.",
    "donor_median_case - donor_median_control.",
    "This workbook is generated from intermediate/stat/check CSV files in the original dataset result directories."
  )
)
readme_lines <- sanitize_df_for_excel(readme_lines)
writeData(wb, "README", readme_lines)
setColWidths(wb, "README", cols = 1:2, widths = c(28, 110))
freezePane(wb, "README", firstActiveRow = 2)

## ALL_SUMO_long
addWorksheet(wb, "ALL_SUMO_long")
writeData(wb, "ALL_SUMO_long", all_tbl)
freezePane(wb, "ALL_SUMO_long", firstActiveRow = 2)
setColWidths(wb, "ALL_SUMO_long", cols = 1:ncol(all_tbl), widths = "auto")

## 每个 comparison 一个 sheet (7 sheets)
for (i in seq_len(nrow(comparison_cfg))) {
  rowi <- comparison_cfg[i, ]
  sheet_name <- glue("{rowi$panel}_{rowi$dataset}_{rowi$disease}")
  sheet_name <- to_utf8_safe(sheet_name)
  sheet_name <- str_sub(sheet_name, 1, 31)

  sub_tbl <- all_tbl %>%
    filter(panel == rowi$panel) %>%
    arrange(gene, factor(cell_type, levels = celltype_order))
  sub_tbl <- sanitize_df_for_excel(sub_tbl)
  check_bad_encoding_column(sub_tbl)

  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, sub_tbl)
  freezePane(wb, sheet_name, firstActiveRow = 2)
  setColWidths(wb, sheet_name, cols = 1:ncol(sub_tbl), widths = "auto")
}

## 表头样式
headerStyle <- createStyle(textDecoration = "bold", halign = "center", valign = "center")
sheet_names_written <- c(
  "README", "ALL_SUMO_long",
  vapply(seq_len(nrow(comparison_cfg)), function(i) {
    nm <- glue("{comparison_cfg$panel[i]}_{comparison_cfg$dataset[i]}_{comparison_cfg$disease[i]}")
    nm <- to_utf8_safe(nm); str_sub(nm, 1, 31)
  }, FUN.VALUE = character(1))
)
for (s in unique(sheet_names_written)) {
  try(addStyle(wb, s, headerStyle, rows = 1, cols = 1:50,
               gridExpand = TRUE, stack = TRUE), silent = TRUE)
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

# -----------------------------
# 七、汇总
# -----------------------------
message("\n========================================================")
message("✅ Supplementary Table (NBDV1senseirevised) 已生成")
message("========================================================")
message("文件: ", out_xlsx)
message("总比较数: ", nrow(comparison_cfg), " (CTE 已排除)")
message("总 sheet 数: ", length(unique(sheet_names_written)),
        " (README + ALL_SUMO_long + 7 个分 comparison)")
message("总数据行: ", nrow(all_tbl),
        " (= 7 comparison × 3 gene × 6 celltype)")
