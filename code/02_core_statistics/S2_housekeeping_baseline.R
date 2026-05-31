###############################################################################
## NBD revision #13 — Extended housekeeping baseline (11 genes)
##
## Version: 7-unit (CTE cohorts removed)
##
## Diff vs FULL9 version:
##   - dataset_configs: removed the trailing two list() entries for
##     GSE155114 and GSE261807 (CTE cohorts).
##   - Output subdirectory: Results_extended11 -> Results_extended11_7units
##   - Output filenames: FULL9 -> FINAL_7units
##   - Everything else (helpers, gene list, statistics, BH logic, spotlight)
##     is identical.
##
## Note: BH correction is within-unit and within-gene (across 6 cell types).
##       Removing CTE units therefore does NOT alter any other unit's BH
##       results, nor the spotlight (syn52082747 AD/PSP neurons), which is
##       unchanged by definition.
##
##   原 6 (传统 housekeeping):
##     GAPDH, ACTB, B2M, HPRT1, PPIA, TBP
##   新加 5 (Eisenberg & Levanon 2013 transcriptome-stability genes):
##     C1ORF43, VCP, PSMB4, SNRPD3, REEP5
##
##   特别意义：VCP 和 PSMB4 是 UPS (ubiquitin-proteasome system) 成员，
##             与 UBL3 同通路 — 它们的稳定证明 UBL3 不是 UPS 整体上调
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
# 第 1 章：路径
# ============================================================
out_root <- "D:/RNA/supptable/SupTable_S6_housekeeping_baseline"
res_dir  <- file.path(out_root, "Results_extended11_7units")  # 新子目录 (7-unit, CTE 移除)
for (d in c(res_dir,
            file.path(res_dir, "01_per_donor_per_gene"),
            file.path(res_dir, "02_per_unit_per_gene_test"),
            file.path(res_dir, "Logs"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_fp <- file.path(res_dir, "Logs", "S6_extended11_7units_log.txt")
sink(log_fp, split = TRUE)
cat("==== START (7-unit, CTE removed) ====\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("SEED:", SEED, "\n\n"); sink()

# ============================================================
# 第 2 章：11 个管家基因（symbol + ENSG）
# ============================================================
HK_GENES <- list(
  # ---- 传统管家基因（前 6 个，与原版一致）----
  GAPDH   = list(symbol = "GAPDH",   ensg = "ENSG00000111640"),
  ACTB    = list(symbol = "ACTB",    ensg = "ENSG00000075624"),
  B2M     = list(symbol = "B2M",     ensg = "ENSG00000166710"),
  HPRT1   = list(symbol = "HPRT1",   ensg = "ENSG00000165704"),
  PPIA    = list(symbol = "PPIA",    ensg = "ENSG00000196262"),
  TBP     = list(symbol = "TBP",     ensg = "ENSG00000112592"),
  # ---- Eisenberg & Levanon 2013 跨组织最稳定基因（新加 5 个）----
  C1ORF43 = list(symbol = "C1orf43", ensg = "ENSG00000143612"),
  VCP     = list(symbol = "VCP",     ensg = "ENSG00000165280"),
  PSMB4   = list(symbol = "PSMB4",   ensg = "ENSG00000159377"),
  SNRPD3  = list(symbol = "SNRPD3",  ensg = "ENSG00000100028"),
  REEP5   = list(symbol = "REEP5",   ensg = "ENSG00000129625")
)

# ============================================================
# 第 3 章：数据集路径
# ============================================================
ec_path    <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds"
sfg_path   <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepH_syn21788402_SFG_obj_celltype6.rds"
syn52_path <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3/stepH_slim_uncompressed.rds"

# ============================================================
# 第 4 章：辅助函数（与之前完全一致）
# ============================================================
load_obj <- function(path, format) if (format == "qs") qs::qread(path) else readRDS(path)

do_step <- function(label, fn) {
  cat("    [step]", label, "...")
  out <- tryCatch(fn(),
                  error = function(e) {
                    cat(" FAIL\n      错误信息:", conditionMessage(e), "\n      步骤:", label, "\n")
                    stop(paste0("STEP_FAILED: ", label, " — ", conditionMessage(e)), call. = FALSE)
                  })
  cat(" ok\n"); out
}

get_counts_robust <- function(obj, assay = "RNA", verbose = TRUE) {
  a <- obj[[assay]]; ncells <- ncol(obj)
  if (verbose) {
    cat("      assay class:", paste(class(a), collapse = "/"), "\n")
    layers_now <- tryCatch(SeuratObject::Layers(a), error = function(e) character(0))
    cat("      layers (n=", length(layers_now), ")\n", sep = "")
  }
  layers_all <- tryCatch(SeuratObject::Layers(a), error = function(e) character(0))
  count_layers <- layers_all[grepl("^counts", layers_all)]
  if (length(count_layers) > 1) {
    if (verbose) cat("      Path-multi: 拼接", length(count_layers), "个 counts layer\n")
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
        if (identical(rownames(m0), all_genes)) mats_aligned[[nm]] <- m0
        else {
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
      if (verbose) cat("      => Path-multi ok\n")
      return(mat_all)
    }
  }
  m <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "counts"), error = function(e) NULL)
  if (!is.null(m) && nrow(m) > 0 && ncol(m) == ncells) {
    if (verbose) cat("      => Path1 ok\n"); return(m)
  }
  m <- tryCatch(SeuratObject::LayerData(a, layer = "counts"), error = function(e) NULL)
  if (!is.null(m) && nrow(m) > 0 && ncol(m) == ncells) {
    if (verbose) cat("      => Path2 ok\n"); return(m)
  }
  m <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "data"), error = function(e) NULL)
  if (!is.null(m) && ncol(m) == ncells) {
    if (verbose) cat("      => Path-data ok\n"); return(m)
  }
  stop("get_counts_robust: 全部路径失败")
}

# 找指定基因（symbol + ensg + 大小写宽容 + 版本号去除）
locate_gene <- function(counts_mat, symbol, ensg) {
  rn <- as.character(rownames(counts_mat))
  if (symbol %in% rn) return(symbol)
  hit <- which(toupper(rn) == toupper(symbol))
  if (length(hit) >= 1) return(rn[hit[1]])
  if (ensg %in% rn) return(ensg)
  rs <- sub("\\.\\d+$", "", rn)
  h <- which(rs == ensg)
  if (length(h) >= 1) return(rn[h[1]])
  NA_character_
}

cliffs_delta <- function(x_d, x_c) {
  if (length(x_d) < 2 || length(x_c) < 2) return(NA_real_)
  dm <- outer(x_d, x_c, FUN = "-")
  (sum(dm > 0) - sum(dm < 0)) / (length(x_d) * length(x_c))
}

es_wilcox <- function(x_d, x_c) {
  res <- list(p_raw = NA_real_, hl = NA_real_, lo = NA_real_, hi = NA_real_, ci_method = "wilcox_HL")
  if (length(x_d) < 2 || length(x_c) < 2) return(res)
  wt <- tryCatch(suppressWarnings(wilcox.test(x_d, x_c, exact = FALSE, conf.int = TRUE, conf.level = 0.95)),
                 error = function(e) NULL)
  if (!is.null(wt)) {
    res$p_raw <- wt$p.value; res$hl <- as.numeric(wt$estimate)
    res$lo <- wt$conf.int[1]; res$hi <- wt$conf.int[2]
  }
  res
}

es_bootstrap <- function(x_d, x_c, B = 2000) {
  res <- list(p_raw = NA_real_, hl = NA_real_, lo = NA_real_, hi = NA_real_, ci_method = "bootstrap_HL")
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
# 第 5 章：Group resolvers
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
  ctrl_alias <- c("Control","CTRL","Ctr","CTR","NC","Normal","control","ctrl","ctr","nc","normal")
  out <- rep(NA_character_, length(g))
  out[g %in% ctrl_alias] <- "Control"
  out[g %in% disease_aliases] <- disease_label
  out
}

# ============================================================
# 第 6 章：7-unit 配置（CTE cohorts removed）
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
  # list(unit_name = "GSE155114", file_format = "qs", ... CTE, DLF white matter)
  # list(unit_name = "GSE261807", file_format = "qs", ... CTE, DL frontal cortex)
)

# ============================================================
# 第 7 章：处理一个 unit（11 基因都算） — 完全沿用原逻辑
# ============================================================
process_unit <- function(cfg) {
  cat("\n========== ", cfg$unit_name, " (", cfg$disease_label, ") ==========\n", sep = "")
  if (!file.exists(cfg$stepH_path)) { cat("File does not exist\n"); return(NULL) }
  
  obj <- do_step("A. load object", function() load_obj(cfg$stepH_path, cfg$file_format))
  obj <- do_step("A2. set default assay", function() { DefaultAssay(obj) <- "RNA"; obj })
  
  md_full <- do_step("B. meta processing", function() {
    md0 <- obj@meta.data
    md0$group_std    <- cfg$group_resolver(md0)
    md0$celltype6_std<- as.character(md0[[cfg$celltype_col]])
    md0$donor_std    <- trimws(as.character(md0[[cfg$donor_col]]))
    md0
  })
  
  md <- do_step("C. filter", function() {
    keep <- !is.na(md_full$group_std) & !is.na(md_full$donor_std) & md_full$donor_std != "" &
      !is.na(md_full$celltype6_std) & md_full$celltype6_std != ""
    md_full[keep, , drop = FALSE]
  })
  if (sum(md$group_std == cfg$disease_label) == 0) { cat("No disease cells\n"); return(NULL) }
  
  cells_keep <- as.character(rownames(md))
  cat("    cells_keep =", length(cells_keep), "\n")
  
  rna_counts_full <- do_step("E. extract counts", function() get_counts_robust(obj, "RNA", verbose = TRUE))
  rna_counts <- do_step("E2. column-index by cells_keep", function() {
    cells_in <- intersect(cells_keep, colnames(rna_counts_full))
    rna_counts_full[, cells_in, drop = FALSE]
  })
  md <- md[colnames(rna_counts), , drop = FALSE]
  
  gene_results <- list()
  
  for (gene_label in names(HK_GENES)) {
    g <- HK_GENES[[gene_label]]
    gene_row <- locate_gene(rna_counts, g$symbol, g$ensg)
    if (is.na(gene_row)) {
      cat("    [WARN]", gene_label, "not found, skipped\n"); next
    }
    cat("    [gene]", gene_label, "row=", gene_row, "\n")
    
    gene_counts <- as.numeric(rna_counts[gene_row, , drop = TRUE])
    cell_ids <- as.character(colnames(rna_counts))
    md_rn <- as.character(rownames(md))
    idx <- match(cell_ids, md_rn)
    
    df <- data.frame(
      cell = cell_ids, gene_count = gene_counts,
      donor = as.character(md$donor_std[idx]),
      group = as.character(md$group_std[idx]),
      celltype = as.character(md$celltype6_std[idx]),
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$donor) & !is.na(df$group) & !is.na(df$celltype), , drop = FALSE]
    
    donor_summary <- df %>% group_by(celltype, donor, group) %>%
      summarise(n_cells_total = n(), n_cells_pos = sum(gene_count > 0),
                prop_pos = sum(gene_count > 0) / n(), .groups = "drop") %>%
      mutate(unit = cfg$unit_name, disease = cfg$disease_label,
             region = cfg$region, gene = gene_label) %>%
      select(unit, disease, region, gene, celltype, donor, group,
             n_cells_total, n_cells_pos, prop_pos) %>%
      as.data.frame()
    
    cts <- sort(unique(donor_summary$celltype))
    rows <- list()
    for (ct in cts) {
      sub <- donor_summary[donor_summary$celltype == ct, , drop = FALSE]
      x_d <- sub$prop_pos[sub$group == cfg$disease_label]
      x_c <- sub$prop_pos[sub$group == "Control"]
      n_d <- length(x_d); n_c <- length(x_c)
      status <- if (n_d >= 4 && n_c >= 4) "YES_wilcox" else "DESCRIPTIVE_n<4"
      es <- if (status == "YES_wilcox") es_wilcox(x_d, x_c) else es_bootstrap(x_d, x_c)
      rows[[length(rows)+1]] <- data.frame(
        unit = cfg$unit_name, disease = cfg$disease_label, region = cfg$region,
        gene = gene_label, celltype = ct,
        n_donors_disease = n_d, n_donors_control = n_c,
        median_disease = if(n_d>0) median(x_d) else NA_real_,
        median_control = if(n_c>0) median(x_c) else NA_real_,
        HL_estimate = es$hl, HL_ci95_low = es$lo, HL_ci95_high = es$hi,
        cliffs_delta = cliffs_delta(x_d, x_c),
        p_raw = es$p_raw, padj_BH = NA_real_,
        test_status = status, ci_method = es$ci_method,
        stringsAsFactors = FALSE
      )
    }
    test_df <- do.call(rbind, rows)
    # within-(unit x gene) BH correction across 6 cell types
    idx_y <- which(test_df$test_status == "YES_wilcox" & !is.na(test_df$p_raw))
    if (length(idx_y) > 0) test_df$padj_BH[idx_y] <- p.adjust(test_df$p_raw[idx_y], method = "BH")
    
    gene_results[[gene_label]] <- list(donor_summary = donor_summary, test = test_df)
  }
  list(gene_results = gene_results)
}

# ============================================================
# 第 8 章：跑全部 7 个 unit
# ============================================================
all_donor <- list(); all_test <- list(); failed_units <- character(0)

for (cfg in dataset_configs) {
  res <- tryCatch(process_unit(cfg),
                  error = function(e) {
                    cat("Unit ", cfg$unit_name, " failed: ", conditionMessage(e), "\n", sep = ""); NULL
                  })
  if (is.null(res)) { failed_units <- c(failed_units, cfg$unit_name); next }
  for (gene_label in names(res$gene_results)) {
    gr <- res$gene_results[[gene_label]]
    key <- paste0(cfg$unit_name, "__", gene_label)
    all_donor[[key]] <- gr$donor_summary
    all_test [[key]] <- gr$test
    write.csv(gr$donor_summary,
              file.path(res_dir, "01_per_donor_per_gene",
                        paste0("per_donor_", cfg$unit_name, "_", gene_label, ".csv")),
              row.names = FALSE)
    write.csv(gr$test,
              file.path(res_dir, "02_per_unit_per_gene_test",
                        paste0("test_", cfg$unit_name, "_", gene_label, ".csv")),
              row.names = FALSE)
  }
}

combined_donor <- do.call(rbind, all_donor)
combined_test  <- do.call(rbind, all_test)

# Output filenames renamed FULL9 -> FINAL_7units
write.csv(combined_donor,
          file.path(res_dir, "SupTable_S6A_per_donor_housekeeping_FINAL_7units_11genes.csv"),
          row.names = FALSE)
write.csv(combined_test,
          file.path(res_dir, "SupTable_S6B_donor_level_test_housekeeping_FINAL_7units_11genes.csv"),
          row.names = FALSE)

# ============================================================
# 第 9 章：摘要 + spotlight (UNCHANGED — spotlight only uses syn52082747)
# ============================================================
sink(log_fp, append = TRUE)
cat("\n\n================ FINAL SUMMARY (7-unit, 11 genes) ================\n")
cat("Units processed:", length(dataset_configs) - length(failed_units), "/", length(dataset_configs), "\n\n")

cat("--- All tests ---\n")
cat("Wilcoxon tests:    ", sum(combined_test$test_status == "YES_wilcox"), "\n")
cat("Descriptive (n<4): ", sum(combined_test$test_status == "DESCRIPTIVE_n<4"), "\n")
cat("padj_BH < 0.05:    ", sum(combined_test$padj_BH < 0.05, na.rm = TRUE), "\n\n")

cat("--- Per-gene summary ---\n")
gene_summary <- combined_test %>% group_by(gene) %>%
  summarise(n_test = sum(test_status == "YES_wilcox"),
            n_sig  = sum(padj_BH < 0.05, na.rm = TRUE),
            min_padj = min(padj_BH, na.rm = TRUE), .groups = "drop") %>%
  arrange(min_padj)
print(as.data.frame(gene_summary))

cat("\n\n*** SPOTLIGHT - syn52082747 neurons, 11 housekeeping genes ***\n")
cat("(Spotlight is unchanged by CTE removal: only uses syn52082747 AD/PSP neurons)\n\n")
spotlight <- combined_test %>%
  filter(unit %in% c("syn52082747_AD", "syn52082747_PSP")) %>%
  filter(celltype %in% c("Excitatory neurons", "Inhibitory neurons")) %>%
  arrange(unit, celltype, gene) %>%
  select(unit, disease, celltype, gene, n_donors_disease, n_donors_control,
         HL_estimate, cliffs_delta, p_raw, padj_BH, test_status)
print(as.data.frame(spotlight))

write.csv(spotlight,
          file.path(res_dir, "SPOTLIGHT_syn52082747_neurons_housekeeping_11genes_7units.csv"),
          row.names = FALSE)

cat("\n\n*** Critical question: across 11 genes x 3 critical comparisons (= 33 tests),\n")
cat("how many reach padj < 0.05? Lower is better. ***\n")
n_sig_critical <- sum(spotlight$padj_BH < 0.05, na.rm = TRUE)
cat("ANSWER:", n_sig_critical, "/ 33\n")
cat("\nNote: this number is IDENTICAL to the FULL9 version because the spotlight\n")
cat("uses only syn52082747 AD + PSP neurons (unaffected by CTE removal).\n")
sink()

cat("\nDONE\n")
cat("Main outputs (7-unit, CTE removed):\n")
cat("  -", file.path(res_dir, "SupTable_S6B_donor_level_test_housekeeping_FINAL_7units_11genes.csv"), "\n")
cat("  -", file.path(res_dir, "SPOTLIGHT_syn52082747_neurons_housekeeping_11genes_7units.csv"), "\n")
