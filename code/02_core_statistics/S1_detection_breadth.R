###############################################################################
## NBD revision — donor-level UBL3 analytical framework
##
## Version: v6 — CTE cohorts (GSE155114, GSE261807) removed from the analytical
##          set on principled grounds (Methods §2.1). The within-unit BH 
##          correction is independent across units, so removing CTE units does
##          NOT alter the 3 syn52082747-V1 positive results — see
##          comparison verification in /mnt/user-data/outputs/.
##
## Diff vs FULL9 version (v5):
##   - dataset_configs: removed the trailing two list() entries for 
##     GSE155114 and GSE261807.

##   - Everything else (helpers, statistics, BH logic) is identical.
###############################################################################

rm(list = ls()); gc()
SEED <- 20251023; set.seed(SEED)
Sys.setenv(LANG = "en")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(dplyr); library(tidyr); library(qs)
})

# ============================================================
# Chapter 1: paths
# ============================================================
out_root <- "D:/RNA/supptable/SupTable_S5_donor_level_combined"
res_dir  <- file.path(out_root, "Results_7units_CTE_removed")
for (d in c(res_dir,
            file.path(res_dir, "01_per_donor"),
            file.path(res_dir, "02_per_unit_test"),
            file.path(res_dir, "Logs"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_fp <- file.path(res_dir, "Logs", "S5_master_log_7units.txt")
sink(log_fp, split = TRUE)
cat("==== START (7-unit, CTE removed) ====\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("SEED:", SEED, "\n\n")
sink()

# ============================================================
# Chapter 2: per-cohort raw object paths
# ============================================================
ec_path    <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds"
sfg_path   <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepH_syn21788402_SFG_obj_celltype6.rds"
syn52_path <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3/stepH_slim_uncompressed.rds"

# ============================================================
# Chapter 3: helper functions (UNCHANGED from v5)
# ============================================================
load_obj <- function(path, format) if (format == "qs") qs::qread(path) else readRDS(path)

do_step <- function(label, fn) {
  cat("    [step]", label, "...")
  out <- tryCatch(fn(),
                  error = function(e) {
                    cat(" FAIL\n")
                    cat("      error: ", conditionMessage(e), "\n")
                    cat("      step:  ", label, "\n")
                    stop(paste0("STEP_FAILED: ", label, " - ", conditionMessage(e)),
                         call. = FALSE)
                  })
  cat(" ok\n")
  out
}

get_counts_robust <- function(obj, assay = "RNA", verbose = TRUE) {
  a <- obj[[assay]]
  ncells <- ncol(obj)
  if (verbose) {
    cat("      assay class:", paste(class(a), collapse = "/"), "\n")
    layers_now <- tryCatch(SeuratObject::Layers(a), error = function(e) character(0))
    cat("      layers (n=", length(layers_now), ")\n", sep = "")
    cat("      target ncells:", ncells, "\n")
  }
  layers_all <- tryCatch(SeuratObject::Layers(a), error = function(e) character(0))
  count_layers <- layers_all[grepl("^counts", layers_all)]
  if (length(count_layers) > 1) {
    if (verbose) cat("      Path-multi: concatenating", length(count_layers), "counts layers\n")
    mats <- list()
    for (ly in count_layers) {
      mm <- tryCatch(SeuratObject::LayerData(a, layer = ly), error = function(e) NULL)
      if (!is.null(mm) && ncol(mm) > 0) mats[[ly]] <- mm
    }
    if (length(mats) >= 1) {
      all_genes <- unique(unlist(lapply(mats, rownames)))
      mats_aligned <- list()
      for (nm in names(mats)) {
        m0 <- mats[[nm]]
        if (identical(rownames(m0), all_genes)) {
          mats_aligned[[nm]] <- m0
        } else {
          m_new <- Matrix::Matrix(0, nrow = length(all_genes), ncol = ncol(m0), sparse = TRUE)
          rownames(m_new) <- all_genes; colnames(m_new) <- colnames(m0)
          common_g <- intersect(all_genes, rownames(m0))
          m_new[common_g, ] <- m0[common_g, , drop = FALSE]
          mats_aligned[[nm]] <- m_new
        }
      }
      mat_all <- Reduce(Matrix::cbind2, mats_aligned)
      dup <- duplicated(colnames(mat_all))
      if (any(dup)) mat_all <- mat_all[, !dup, drop = FALSE]
      miss <- setdiff(colnames(obj), colnames(mat_all))
      if (length(miss) > 0) {
        m_fill <- Matrix::Matrix(0, nrow = nrow(mat_all), ncol = length(miss), sparse = TRUE)
        rownames(m_fill) <- rownames(mat_all); colnames(m_fill) <- miss
        mat_all <- Matrix::cbind2(mat_all, m_fill)
      }
      mat_all <- mat_all[, colnames(obj), drop = FALSE]
      if (verbose) cat("      => Path-multi ok, dim=", nrow(mat_all), "x", ncol(mat_all), "\n")
      return(mat_all)
    }
  }
  m <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "counts"),
                error = function(e) NULL)
  if (!is.null(m) && nrow(m) > 0 && ncol(m) == ncells) {
    if (verbose) cat("      => Path1 (GetAssayData) ok, dim=", nrow(m), "x", ncol(m), "\n")
    return(m)
  }
  m <- tryCatch(SeuratObject::LayerData(a, layer = "counts"),
                error = function(e) NULL)
  if (!is.null(m) && nrow(m) > 0 && ncol(m) == ncells) {
    if (verbose) cat("      => Path2 (LayerData) ok, dim=", nrow(m), "x", ncol(m), "\n")
    return(m)
  }
  m <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "data"),
                error = function(e) NULL)
  if (!is.null(m) && ncol(m) == ncells) {
    if (verbose) cat("      => Path-data slot ok\n")
    return(m)
  }
  stop("get_counts_robust: all paths failed")
}

locate_UBL3 <- function(counts_mat) {
  rn <- as.character(rownames(counts_mat))
  if ("UBL3" %in% rn) return("UBL3")
  if ("ENSG00000122042" %in% rn) return("ENSG00000122042")
  rs <- sub("\\.\\d+$", "", rn)
  h <- which(rs == "ENSG00000122042")
  if (length(h) >= 1) return(rn[h[1]])
  NA_character_
}

cliffs_delta <- function(x_d, x_c) {
  if (length(x_d) < 2 || length(x_c) < 2) return(NA_real_)
  dm <- outer(x_d, x_c, FUN = "-")
  (sum(dm > 0) - sum(dm < 0)) / (length(x_d) * length(x_c))
}

es_wilcox <- function(x_d, x_c) {
  res <- list(p_raw = NA_real_, hl = NA_real_, lo = NA_real_, hi = NA_real_,
              ci_method = "wilcox_HL")
  if (length(x_d) < 2 || length(x_c) < 2) return(res)
  wt <- tryCatch(suppressWarnings(wilcox.test(x_d, x_c, exact = FALSE,
                                              conf.int = TRUE, conf.level = 0.95)),
                 error = function(e) NULL)
  if (!is.null(wt)) {
    res$p_raw <- wt$p.value; res$hl <- as.numeric(wt$estimate)
    res$lo <- wt$conf.int[1]; res$hi <- wt$conf.int[2]
  }
  res
}

es_bootstrap <- function(x_d, x_c, B = 2000) {
  res <- list(p_raw = NA_real_, hl = NA_real_, lo = NA_real_, hi = NA_real_,
              ci_method = "bootstrap_HL")
  if (length(x_d) < 2 || length(x_c) < 2) return(res)
  res$hl <- median(outer(x_d, x_c, FUN = "-"))
  hl_b <- replicate(B, {
    xd <- sample(x_d, replace = TRUE); xc <- sample(x_c, replace = TRUE)
    median(outer(xd, xc, FUN = "-"))
  })
  res$lo <- as.numeric(quantile(hl_b, 0.025, na.rm = TRUE))
  res$hi <- as.numeric(quantile(hl_b, 0.975, na.rm = TRUE))
  res
}

# ============================================================
# Chapter 4: group resolvers (UNCHANGED from v5)
# ============================================================
resolver_braakstage <- function(md, disease_label = "AD") {
  bs <- trimws(as.character(md$BraakStage))
  out <- rep(NA_character_, length(bs))
  out[bs %in% c("0","Braak 0","BraakStage_0","Braak0")] <- "Control"
  out[bs %in% c("6","Braak 6","BraakStage_6","Braak6")] <- disease_label
  out
}
resolver_group_simple <- function(md, group_col, disease_aliases, disease_label) {
  g <- trimws(as.character(md[[group_col]]))
  ctrl_alias <- c("Control","CTRL","Ctr","CTR","NC","Normal",
                  "control","ctrl","ctr","nc","normal")
  out <- rep(NA_character_, length(g))
  out[g %in% ctrl_alias] <- "Control"
  out[g %in% disease_aliases] <- disease_label
  out
}

# ============================================================
# Chapter 5: 7-unit configuration (CTE cohorts removed)
# ============================================================
# Rationale for CTE removal (Methods 2.1): the two CTE cohorts (GSE155114,
# GSE261807) differ from the syn52082747 V1 cohort along three simultaneous
# and confounded dimensions - source dataset, anatomical region (deep frontal
# white matter and dorsolateral frontal cortex, respectively, vs primary
# visual cortex) and sequencing platform - so any negative finding in those
# cohorts cannot be attributed unambiguously to CTE biology rather than to
# inter-dataset, inter-regional or platform differences.
# ------------------------------------------------------------
dataset_configs <- list(
  list(unit_name = "GSE157827", file_format = "rds",
       stepH_path = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results/stepH_obj_celltype6_named.rds",
       disease_label = "AD", region = "Middle frontal gyrus",
       group_resolver = function(md) resolver_group_simple(md, "group", c("AD"), "AD"),
       celltype_col = "celltype6", donor_col = "sample"),
  list(unit_name = "GSE174367", file_format = "rds",
       stepH_path = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/stepH_obj_celltype6_named.rds",
       disease_label = "AD", region = "Prefrontal cortex",
       group_resolver = function(md) resolver_group_simple(md, "Diagnosis", c("AD"), "AD"),
       celltype_col = "celltype6", donor_col = "SampleID"),
  list(unit_name = "syn21788402_EC", file_format = "rds", stepH_path = ec_path,
       disease_label = "AD", region = "Entorhinal cortex",
       group_resolver = function(md) resolver_braakstage(md, "AD"),
       celltype_col = "celltype6", donor_col = "SampleID"),
  list(unit_name = "syn21788402_SFG", file_format = "rds", stepH_path = sfg_path,
       disease_label = "AD", region = "Superior frontal gyrus",
       group_resolver = function(md) resolver_braakstage(md, "AD"),
       celltype_col = "celltype6", donor_col = "SampleID"),
  list(unit_name = "syn52082747_AD", file_format = "rds", stepH_path = syn52_path,
       disease_label = "AD", region = "Primary visual cortex (V1)",
       group_resolver = function(md) resolver_group_simple(md, "group4", c("AD"), "AD"),
       celltype_col = "celltype6", donor_col = "autopsy_id"),
  list(unit_name = "syn52082747_FTD", file_format = "rds", stepH_path = syn52_path,
       disease_label = "FTD", region = "Primary visual cortex (V1)",
       group_resolver = function(md) resolver_group_simple(md, "group4", c("FTD"), "FTD"),
       celltype_col = "celltype6", donor_col = "autopsy_id"),
  list(unit_name = "syn52082747_PSP", file_format = "rds", stepH_path = syn52_path,
       disease_label = "PSP", region = "Primary visual cortex (V1)",
       group_resolver = function(md) resolver_group_simple(md, "group4", c("PSP"), "PSP"),
       celltype_col = "celltype6", donor_col = "autopsy_id")
  # ---- REMOVED on principled grounds (see Methods 2.1) ----
  # list(unit_name = "GSE155114", file_format = "qs", ...  CTE, DLF white matter)
  # list(unit_name = "GSE261807", file_format = "qs", ...  CTE, DL frontal cortex)
)

# ============================================================
# Chapter 6: process_unit() (UNCHANGED from v5)
# ============================================================
process_unit <- function(cfg) {
  cat("\n========== ", cfg$unit_name, " (", cfg$disease_label, ") ==========\n", sep = "")
  if (!file.exists(cfg$stepH_path)) {
    cat("File does not exist: ", cfg$stepH_path, "\n"); return(NULL)
  }
  obj <- do_step("A. load object", function() load_obj(cfg$stepH_path, cfg$file_format))
  obj <- do_step("A2. set default assay", function() { DefaultAssay(obj) <- "RNA"; obj })
  cat("    cells =", ncol(obj), "\n")
  md_full <- do_step("B. extract meta.data", function() obj@meta.data)
  md_full <- do_step("B1. resolve group_std", function() {
    md_full$group_std <- cfg$group_resolver(md_full); md_full
  })
  if (!cfg$donor_col %in% colnames(md_full)) {
    cat("donor_col '", cfg$donor_col, "' not found\n", sep = ""); return(NULL)
  }
  md_full <- do_step("B2-3. cast celltype/donor as character", function() {
    md_full$celltype6_std <- as.character(md_full[[cfg$celltype_col]])
    md_full$donor_std     <- trimws(as.character(md_full[[cfg$donor_col]]))
    md_full
  })
  md <- do_step("C. filter keep mask", function() {
    keep <- !is.na(md_full$group_std) & !is.na(md_full$donor_std) & md_full$donor_std != "" &
      !is.na(md_full$celltype6_std) & md_full$celltype6_std != ""
    md_full[keep, , drop = FALSE]
  })
  if (sum(md$group_std == cfg$disease_label) == 0) {
    cat("No disease cells\n"); return(NULL)
  }
  cells_keep <- as.character(rownames(md))
  cat("    after filter: cells =", length(cells_keep),
      "| disease donors =",
      length(unique(md$donor_std[md$group_std == cfg$disease_label])),
      "| control donors =",
      length(unique(md$donor_std[md$group_std == "Control"])), "\n")
  rna_counts_full <- do_step("E. extract full obj counts",
                             function() get_counts_robust(obj, "RNA", verbose = TRUE))
  rna_counts <- do_step("E2. column-index counts by cells_keep", function() {
    cells_in_counts <- intersect(cells_keep, colnames(rna_counts_full))
    if (length(cells_in_counts) == 0) stop("no cells_keep matching counts columns")
    rna_counts_full[, cells_in_counts, drop = FALSE]
  })
  cat("    counts subset dim:", nrow(rna_counts), "x", ncol(rna_counts), "\n")
  md <- md[colnames(rna_counts), , drop = FALSE]
  ubl3_row <- do_step("F. locate UBL3", function() locate_UBL3(rna_counts))
  if (is.na(ubl3_row)) { cat("UBL3 not found\n"); return(NULL) }
  cat("    UBL3 row:", ubl3_row, "\n")
  lib_size    <- do_step("G1. compute lib_size", function() Matrix::colSums(rna_counts))
  ubl3_counts <- do_step("G2. extract UBL3 row",
                         function() as.numeric(rna_counts[ubl3_row, , drop = TRUE]))
  ubl3_log1p  <- do_step("G3. log-norm",
                         function() log1p((ubl3_counts / pmax(lib_size, 1)) * 1e4))
  df <- do_step("H. build df", function() {
    cell_ids <- as.character(colnames(rna_counts))
    md_rn    <- as.character(rownames(md))
    idx      <- match(cell_ids, md_rn)
    out <- data.frame(
      cell = cell_ids, ubl3_count = ubl3_counts, ubl3_lognorm = ubl3_log1p,
      donor = as.character(md$donor_std[idx]),
      group = as.character(md$group_std[idx]),
      celltype = as.character(md$celltype6_std[idx]),
      stringsAsFactors = FALSE)
    out[!is.na(out$donor) & !is.na(out$group) & !is.na(out$celltype), , drop = FALSE]
  })
  cat("    df rows after match cleanup:", nrow(df), "\n")
  donor_summary <- do_step("I. donor x celltype summarise", function() {
    df %>% group_by(celltype, donor, group) %>%
      summarise(n_cells_total = n(), n_cells_ubl3pos = sum(ubl3_count > 0),
                prop_ubl3pos = sum(ubl3_count > 0) / n(),
                median_lognorm_ubl3pos = ifelse(sum(ubl3_count > 0) > 0,
                                                median(ubl3_lognorm[ubl3_count > 0]),
                                                NA_real_),
                .groups = "drop") %>%
      mutate(unit = cfg$unit_name, disease = cfg$disease_label, region = cfg$region) %>%
      select(unit, disease, region, celltype, donor, group,
             n_cells_total, n_cells_ubl3pos, prop_ubl3pos, median_lognorm_ubl3pos) %>%
      as.data.frame()
  })
  test_df <- do_step("J. per cell-type test", function() {
    cts <- sort(unique(donor_summary$celltype))
    rows <- list()
    for (ct in cts) {
      sub <- donor_summary[donor_summary$celltype == ct, , drop = FALSE]
      xp_d <- sub$prop_ubl3pos[sub$group == cfg$disease_label]
      xp_c <- sub$prop_ubl3pos[sub$group == "Control"]
      sub_e <- sub[!is.na(sub$median_lognorm_ubl3pos), , drop = FALSE]
      xe_d <- sub_e$median_lognorm_ubl3pos[sub_e$group == cfg$disease_label]
      xe_c <- sub_e$median_lognorm_ubl3pos[sub_e$group == "Control"]
      n_pd <- length(xp_d); n_pc <- length(xp_c)
      n_ed <- length(xe_d); n_ec <- length(xe_c)
      status_p <- if (n_pd >= 4 && n_pc >= 4) "YES_wilcox" else "DESCRIPTIVE_n<4"
      status_e <- if (n_ed >= 4 && n_ec >= 4) "YES_wilcox" else "DESCRIPTIVE_n<4"
      es_p <- if (status_p == "YES_wilcox") es_wilcox(xp_d, xp_c) else es_bootstrap(xp_d, xp_c)
      es_e <- if (status_e == "YES_wilcox") es_wilcox(xe_d, xe_c) else es_bootstrap(xe_d, xe_c)
      rows[[length(rows)+1]] <- data.frame(
        unit = cfg$unit_name, disease = cfg$disease_label,
        region = cfg$region, celltype = ct,
        prop_n_donors_disease = n_pd, prop_n_donors_control = n_pc,
        prop_median_disease = if(n_pd>0) median(xp_d) else NA_real_,
        prop_median_control = if(n_pc>0) median(xp_c) else NA_real_,
        prop_HL_estimate = es_p$hl, prop_HL_ci95_low = es_p$lo, prop_HL_ci95_high = es_p$hi,
        prop_cliffs_delta = cliffs_delta(xp_d, xp_c),
        prop_p_raw = es_p$p_raw, prop_padj_BH = NA_real_,
        prop_test_status = status_p, prop_ci_method = es_p$ci_method,
        expr_n_donors_disease = n_ed, expr_n_donors_control = n_ec,
        expr_median_disease = if(n_ed>0) median(xe_d) else NA_real_,
        expr_median_control = if(n_ec>0) median(xe_c) else NA_real_,
        expr_HL_estimate = es_e$hl, expr_HL_ci95_low = es_e$lo, expr_HL_ci95_high = es_e$hi,
        expr_cliffs_delta = cliffs_delta(xe_d, xe_c),
        expr_p_raw = es_e$p_raw, expr_padj_BH = NA_real_,
        expr_test_status = status_e, expr_ci_method = es_e$ci_method,
        stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
  })
  test_df <- do_step("K. BH correction (within-unit, across 6 cell types)", function() {
    td <- test_df
    idx_p <- which(td$prop_test_status == "YES_wilcox" & !is.na(td$prop_p_raw))
    if (length(idx_p) > 0)
      td$prop_padj_BH[idx_p] <- p.adjust(td$prop_p_raw[idx_p], method = "BH")
    idx_e <- which(td$expr_test_status == "YES_wilcox" & !is.na(td$expr_p_raw))
    if (length(idx_e) > 0)
      td$expr_padj_BH[idx_e] <- p.adjust(td$expr_p_raw[idx_e], method = "BH")
    td
  })
  list(donor_summary = donor_summary, test = test_df)
}

# ============================================================
# Chapter 7: run 7 units
# ============================================================
all_donor <- list(); all_test <- list()
failed_units <- character(0)

for (cfg in dataset_configs) {
  res <- tryCatch(process_unit(cfg),
                  error = function(e) {
                    cat("Unit ", cfg$unit_name, " failed: ", 
                        conditionMessage(e), "\n", sep = ""); NULL
                  })
  if (is.null(res)) { failed_units <- c(failed_units, cfg$unit_name); next }
  all_donor[[cfg$unit_name]] <- res$donor_summary
  all_test [[cfg$unit_name]] <- res$test
  write.csv(res$donor_summary,
            file.path(res_dir, "01_per_donor",
                      paste0("per_donor_summary_", cfg$unit_name, ".csv")),
            row.names = FALSE)
  write.csv(res$test,
            file.path(res_dir, "02_per_unit_test",
                      paste0("test_", cfg$unit_name, ".csv")),
            row.names = FALSE)
}

# ============================================================
# Chapter 8: combined master outputs
# ============================================================
combined_donor <- do.call(rbind, all_donor)
combined_test  <- do.call(rbind, all_test)

# Outputs renamed FULL9 -> FINAL_7units
write.csv(combined_donor,
          file.path(res_dir, "SupTable_S5A_per_donor_data_FINAL_7units.csv"),
          row.names = FALSE)
write.csv(combined_test,
          file.path(res_dir, "SupTable_S5B_donor_level_test_FINAL_7units.csv"),
          row.names = FALSE)

# ============================================================
# Chapter 9: final summary
# ============================================================
sink(log_fp, append = TRUE)
cat("\n\n================ FINAL SUMMARY (7-unit version) ================\n")
cat("Units processed: ", length(all_donor), "/ 7\n")
if (length(failed_units) > 0)
  cat("Units FAILED:    ", paste(failed_units, collapse = ", "), "\n")
cat("Combined rows (Panel A):", nrow(combined_donor), "\n")
cat("Combined rows (Panel B):", nrow(combined_test), "\n")
cat("\n--- Detection-breadth (prop) endpoint ---\n")
cat("Wilcoxon tests:    ", sum(combined_test$prop_test_status == "YES_wilcox"), "\n")
cat("Descriptive (n<4): ", sum(combined_test$prop_test_status == "DESCRIPTIVE_n<4"), "\n")
cat("padj_BH < 0.05:    ", sum(combined_test$prop_padj_BH < 0.05, na.rm = TRUE), "\n")
cat("Min padj_BH:       ", min(combined_test$prop_padj_BH, na.rm = TRUE), "\n\n")
cat("--- Conditional-expression endpoint ---\n")
cat("Wilcoxon tests:    ", sum(combined_test$expr_test_status == "YES_wilcox"), "\n")
cat("Descriptive (n<4): ", sum(combined_test$expr_test_status == "DESCRIPTIVE_n<4"), "\n")
cat("padj_BH < 0.05:    ", sum(combined_test$expr_padj_BH < 0.05, na.rm = TRUE), "\n")
cat("Min padj_BH:       ", min(combined_test$expr_padj_BH, na.rm = TRUE), "\n\n")
cat("==== sessionInfo ====\n"); print(sessionInfo())
cat("\n==== END ====\n")
sink()

cat("\nDONE\n")
cat("Main outputs (7-unit, CTE removed):\n")
cat("  -", file.path(res_dir, "SupTable_S5A_per_donor_data_FINAL_7units.csv"), "\n")
cat("  -", file.path(res_dir, "SupTable_S5B_donor_level_test_FINAL_7units.csv"), "\n")
