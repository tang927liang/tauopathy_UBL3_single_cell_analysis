###############################################################################
## S4 / 04 - Manual recomputation from cell layer (Panel D)
## Supplementary Table S4 (robustness/sensitivity) | MN submission
##
## INPUT: syn52082747 Seurat object (.rds) + S1A + S1B. OUTPUT: PanelD_*.csv.
## Run order: 01 -> 02 -> 03 -> 04 -> 05. Set REPO; outputs go to S4_robustness/.
###############################################################################
rm(list = ls()); gc()
SEED <- 20251023
set.seed(SEED)
options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
  library(Seurat)
  library(dplyr)
})

# ============================================================
# Path configuration (CONFIRMED to match NBD revision source script)
# ============================================================
## === 改这里：代码库根目录 ===
REPO <- "D:/RNA/Code/UBL3_tauopathy"

INPUT_SYN_OBJ <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3/stepH_slim_uncompressed.rds"
FILE_FORMAT   <- "rds"     # confirmed .rds (NOT .qs)

INPUT_S5A     <- file.path(REPO, "output", "stats_tables", "S1_detection_breadth", "S1A_per_donor_data.csv")
INPUT_S5B     <- file.path(REPO, "output", "stats_tables", "S1_detection_breadth", "S1B_donor_level_test.csv")
OUTPUT_DIR    <- file.path(REPO, "output", "stats_tables", "S4_robustness")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

TOL_PROP <- 1e-9   # per-donor prop tolerance vs S5A
TOL_STAT <- 1e-6   # test-statistic tolerance vs S5B

# ============================================================
# Helper: control alias list (verbatim from NBD revision, line 196-197)
# ============================================================
CTRL_ALIAS <- c("Control","CTRL","Ctr","CTR","NC","Normal",
                "control","ctrl","ctr","nc","normal")

# Three focal recomputations (each is a (unit, focal_celltype) pair)
focal <- data.frame(
  unit_label     = c("syn52082747_AD",
                     "syn52082747_PSP",
                     "syn52082747_PSP"),
  disease_label  = c("AD", "PSP", "PSP"),
  focal_celltype = c("Excitatory neurons",
                     "Excitatory neurons",
                     "Inhibitory neurons"),
  stringsAsFactors = FALSE
)

# ============================================================
# Step 1. Load deposited cell-level object
# ============================================================
if (!file.exists(INPUT_SYN_OBJ)) stop("syn52082747 object not found: ", INPUT_SYN_OBJ)
cat("Loading syn52082747 cell-level object:\n  ", INPUT_SYN_OBJ, "\n")

if (FILE_FORMAT == "qs") {
  if (!requireNamespace("qs", quietly = TRUE)) stop("Install 'qs' package.")
  obj <- qs::qread(INPUT_SYN_OBJ)
} else {
  obj <- readRDS(INPUT_SYN_OBJ)
}
DefaultAssay(obj) <- "RNA"
cat("  Total cells:", ncol(obj), "\n")

# ============================================================
# Step 2. Combine multi-layer counts (Seurat v5) — mirrors
#         get_counts_robust() in NBD revision script
# ============================================================
get_counts_matrix_allcells <- function(obj, assay = "RNA") {
  a <- obj[[assay]]
  layers_all   <- tryCatch(SeuratObject::Layers(a), error = function(e) character(0))
  count_layers <- layers_all[grepl("^counts", layers_all)]
  cat("  counts layers found: ", length(count_layers),
      " (", paste(count_layers, collapse=", "), ")\n", sep = "")

  if (length(count_layers) > 1) {
    mats <- list()
    for (ly in count_layers) {
      mm <- tryCatch(SeuratObject::LayerData(a, layer = ly),
                     error = function(e) NULL)
      if (!is.null(mm) && ncol(mm) > 0) mats[[ly]] <- mm
    }
    all_genes <- unique(unlist(lapply(mats, rownames)))
    mats_aligned <- list()
    for (nm in names(mats)) {
      m0 <- mats[[nm]]
      if (identical(rownames(m0), all_genes)) {
        mats_aligned[[nm]] <- m0
      } else {
        m_new <- Matrix::Matrix(0, nrow = length(all_genes),
                                ncol = ncol(m0), sparse = TRUE)
        rownames(m_new) <- all_genes
        colnames(m_new) <- colnames(m0)
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
      m_fill <- Matrix::Matrix(0, nrow = nrow(mat_all),
                               ncol = length(miss), sparse = TRUE)
      rownames(m_fill) <- rownames(mat_all)
      colnames(m_fill) <- miss
      mat_all <- Matrix::cbind2(mat_all, m_fill)
    }
    return(mat_all[, colnames(obj), drop = FALSE])
  }
  m <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "counts"),
                error = function(e) NULL)
  if (!is.null(m) && ncol(m) == ncol(obj)) return(m)
  m <- tryCatch(SeuratObject::LayerData(a, layer = "counts"),
                error = function(e) NULL)
  if (!is.null(m) && ncol(m) == ncol(obj)) return(m)
  stop("Cannot extract counts matrix.")
}

# ============================================================
# Step 3. Locate UBL3 row — mirrors locate_UBL3() in NBD revision
# ============================================================
locate_UBL3 <- function(counts_mat) {
  rn <- as.character(rownames(counts_mat))
  if ("UBL3" %in% rn) return("UBL3")
  if ("ENSG00000122042" %in% rn) return("ENSG00000122042")
  rs <- sub("\\.\\d+$", "", rn)
  h <- which(rs == "ENSG00000122042")
  if (length(h) >= 1) return(rn[h[1]])
  NA_character_
}

cat("\nExtracting counts matrix and locating UBL3...\n")
rna_counts <- get_counts_matrix_allcells(obj, "RNA")
cat("  counts dim:", nrow(rna_counts), "x", ncol(rna_counts), "\n")

ubl3_row <- locate_UBL3(rna_counts)
if (is.na(ubl3_row)) stop("UBL3 not found in counts rownames.")
cat("  UBL3 row:", ubl3_row, "\n")

ubl3_counts <- as.numeric(rna_counts[ubl3_row, , drop = TRUE])
names(ubl3_counts) <- colnames(rna_counts)
ubl3_pos <- ubl3_counts > 0

cat("  Total cells with UBL3 count > 0:", sum(ubl3_pos),
    "(", round(100 * mean(ubl3_pos), 2), "%)\n", sep = "")

# ============================================================
# Step 4. Build cell metadata (mirrors NBD revision lines 262-281)
# ============================================================
md_full <- obj@meta.data
need_cols <- c("celltype6", "autopsy_id", "group4")
miss <- setdiff(need_cols, colnames(md_full))
if (length(miss) > 0) stop("obj meta missing columns: ",
                            paste(miss, collapse = ", "))

cell_md <- data.frame(
  cell      = rownames(md_full),
  celltype6 = as.character(md_full$celltype6),
  donor_std = trimws(as.character(md_full$autopsy_id)),
  group_raw = trimws(as.character(md_full$group4)),
  stringsAsFactors = FALSE
)

# Group resolver — mirrors resolver_group_simple() VERBATIM
# (Control alias list per NBD revision line 196-197)
cell_md$group_std <- NA_character_
cell_md$group_std[cell_md$group_raw %in% CTRL_ALIAS] <- "Control"
cell_md$group_std[cell_md$group_raw == "AD"]        <- "AD"
cell_md$group_std[cell_md$group_raw == "FTD"]       <- "FTD"
cell_md$group_std[cell_md$group_raw == "PSP"]       <- "PSP"

# Keep mask — mirrors NBD revision lines 274-278
keep_mask <- !is.na(cell_md$group_std) &
             !is.na(cell_md$donor_std) & cell_md$donor_std != "" &
             !is.na(cell_md$celltype6) & cell_md$celltype6 != ""

cat("\nMeta filtering (mirrors NBD revision keep mask):\n")
cat("  total cells :", nrow(cell_md), "\n")
cat("  kept cells  :", sum(keep_mask), "\n")
cat("  dropped     :", sum(!keep_mask),
    "(NA group/donor/celltype)\n", sep = "")

cell_md       <- cell_md[keep_mask, , drop = FALSE]
cell_md$ubl3_pos <- as.logical(ubl3_pos[cell_md$cell])

cat("  Group_std distribution (kept cells):\n")
print(table(cell_md$group_std))

# ============================================================
# Step 5. Three focal recomputations
# ============================================================
s5a <- read.csv(INPUT_S5A, stringsAsFactors = FALSE)
s5b <- read.csv(INPUT_S5B, stringsAsFactors = FALSE)

cliffs_delta <- function(x, y) {
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  d <- outer(x, y, "-")
  (sum(d > 0) - sum(d < 0)) / (length(x) * length(y))
}
hodges_lehmann <- function(x, y) {
  # v4 fix: use wilcox.test(..., conf.int=TRUE)$estimate (R's pseudomedian /
  # uniroot algorithm) — identical to the function used to generate S5B
  # (NBD revision script line 161). The previous median-of-pairwise-differences
  # formula is theoretically equivalent but differs slightly with ties.
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  wt <- tryCatch(
    suppressWarnings(wilcox.test(x, y, exact = FALSE,
                                  conf.int = TRUE, conf.level = 0.95)),
    error = function(e) NULL
  )
  if (is.null(wt)) return(NA_real_)
  as.numeric(wt$estimate)
}
wilcox_p <- function(x, y) {
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(x, y, exact = FALSE))$p.value
}

celltypes_all <- c("Astrocytes", "Endothelial",
                   "Excitatory neurons", "Inhibitory neurons",
                   "Microglia", "Oligodendrocytes")

per_donor_rows   <- list()
test_summary_rows <- list()

for (k in seq_len(nrow(focal))) {
  unit_label   <- focal$unit_label[k]
  disease_lab  <- focal$disease_label[k]
  focal_ct     <- focal$focal_celltype[k]

  cat("\n--------------------------------------------------\n")
  cat(sprintf("Focal #%d: %s / %s (disease=%s)\n",
              k, unit_label, focal_ct, disease_lab))
  cat("--------------------------------------------------\n")

  # Restrict to disease + Control cells (mirrors per-unit resolver result)
  cells_unit <- cell_md[cell_md$group_std %in% c(disease_lab, "Control"), ,
                        drop = FALSE]

  cat("  Cells in this unit (disease + Control): ", nrow(cells_unit), "\n", sep="")

  # Per-donor x celltype summary (mirrors NBD revision lines 318-330)
  donor_summary_re <- cells_unit %>%
    dplyr::group_by(celltype6, donor_std, group_std) %>%
    dplyr::summarise(
      n_cells_total      = dplyr::n(),
      n_cells_ubl3pos    = sum(ubl3_pos),
      prop_ubl3pos       = sum(ubl3_pos) / dplyr::n(),
      .groups            = "drop"
    ) %>%
    as.data.frame()

  # Now: compare focal cell-type per-donor rows against S5A
  expected <- s5a[s5a$unit == unit_label & s5a$celltype == focal_ct,
                  c("donor", "group", "n_cells_total",
                    "n_cells_ubl3pos", "prop_ubl3pos")]

  for (i in seq_len(nrow(expected))) {
    don <- expected$donor[i]
    grp <- expected$group[i]

    row_re <- donor_summary_re[donor_summary_re$donor_std == don &
                               donor_summary_re$celltype6 == focal_ct, ]
    if (nrow(row_re) == 0) {
      # Donor has no cells of this celltype in recomputation — record as FAIL
      per_donor_rows[[length(per_donor_rows) + 1]] <- data.frame(
        focal_id              = paste0("focal", k),
        unit                  = unit_label,
        celltype              = focal_ct,
        donor                 = don,
        group                 = grp,
        n_cells_total_S5A     = expected$n_cells_total[i],
        n_cells_total_recomp  = 0,
        n_cells_ubl3pos_S5A   = expected$n_cells_ubl3pos[i],
        n_cells_ubl3pos_recomp = 0,
        prop_ubl3pos_S5A      = expected$prop_ubl3pos[i],
        prop_ubl3pos_recomp   = NA_real_,
        abs_diff_prop         = NA_real_,
        PASS                  = FALSE,
        stringsAsFactors = FALSE
      )
      cat(sprintf("  FAIL: donor %s — no recomputed row for %s\n", don, focal_ct))
      next
    }

    d_n_total <- row_re$n_cells_total  - expected$n_cells_total[i]
    d_n_pos   <- row_re$n_cells_ubl3pos - expected$n_cells_ubl3pos[i]
    d_prop    <- row_re$prop_ubl3pos    - expected$prop_ubl3pos[i]
    pass <- (d_n_total == 0) && (d_n_pos == 0) && (abs(d_prop) < TOL_PROP)

    per_donor_rows[[length(per_donor_rows) + 1]] <- data.frame(
      focal_id              = paste0("focal", k),
      unit                  = unit_label,
      celltype              = focal_ct,
      donor                 = don,
      group                 = grp,
      n_cells_total_S5A     = expected$n_cells_total[i],
      n_cells_total_recomp  = row_re$n_cells_total,
      n_cells_ubl3pos_S5A   = expected$n_cells_ubl3pos[i],
      n_cells_ubl3pos_recomp = row_re$n_cells_ubl3pos,
      prop_ubl3pos_S5A      = expected$prop_ubl3pos[i],
      prop_ubl3pos_recomp   = row_re$prop_ubl3pos,
      abs_diff_prop         = abs(d_prop),
      PASS                  = pass,
      stringsAsFactors = FALSE
    )

    if (!pass) {
      cat(sprintf("  FAIL: donor %s (group %s): n_total %d vs %d  n_pos %d vs %d  prop %.6f vs %.6f\n",
                  don, grp, row_re$n_cells_total, expected$n_cells_total[i],
                  row_re$n_cells_ubl3pos, expected$n_cells_ubl3pos[i],
                  row_re$prop_ubl3pos, expected$prop_ubl3pos[i]))
    }
  }

  # Re-run Wilcoxon for focal celltype using recomputed donor props
  xp_d <- donor_summary_re$prop_ubl3pos[donor_summary_re$celltype6 == focal_ct &
                                         donor_summary_re$group_std == disease_lab]
  xp_c <- donor_summary_re$prop_ubl3pos[donor_summary_re$celltype6 == focal_ct &
                                         donor_summary_re$group_std == "Control"]
  cat(sprintf("  Recomputed N: %s=%d, Control=%d\n",
              disease_lab, length(xp_d), length(xp_c)))

  # Recompute Wilcoxon + Cliff's delta + HL for focal ct
  re_p     <- wilcox_p(xp_d, xp_c)
  re_delta <- cliffs_delta(xp_d, xp_c)
  re_HL    <- hodges_lehmann(xp_d, xp_c)

  # Compute within-unit BH (across 6 cell types of this unit) for focal ct
  ct_p_vec <- numeric(length(celltypes_all))
  for (ci in seq_along(celltypes_all)) {
    ct <- celltypes_all[ci]
    xd <- donor_summary_re$prop_ubl3pos[donor_summary_re$celltype6 == ct &
                                        donor_summary_re$group_std == disease_lab]
    xc <- donor_summary_re$prop_ubl3pos[donor_summary_re$celltype6 == ct &
                                        donor_summary_re$group_std == "Control"]
    ct_p_vec[ci] <- wilcox_p(xd, xc)
  }
  ct_padj <- p.adjust(ct_p_vec, method = "BH")
  i_focal <- which(celltypes_all == focal_ct)
  re_padj <- ct_padj[i_focal]

  # Match against S5B
  ref <- s5b[s5b$unit == unit_label & s5b$celltype == focal_ct, ]
  stopifnot(nrow(ref) == 1L)

  d_p     <- abs(re_p     - ref$prop_p_raw)
  d_padj  <- abs(re_padj  - ref$prop_padj_BH)
  d_delta <- abs(re_delta - ref$prop_cliffs_delta)
  d_HL    <- abs(re_HL    - ref$prop_HL_estimate)

  pass_p     <- d_p     < TOL_STAT
  pass_padj  <- d_padj  < TOL_STAT
  pass_delta <- d_delta < TOL_STAT
  pass_HL    <- d_HL    < TOL_STAT

  cat(sprintf("\n  Wilcoxon recomputation match against S5B:\n"))
  cat(sprintf("    raw p   : %.6f vs %.6f   diff=%.3e   %s\n",
              re_p, ref$prop_p_raw, d_p, ifelse(pass_p, "PASS", "FAIL")))
  cat(sprintf("    padj_BH : %.6f vs %.6f   diff=%.3e   %s\n",
              re_padj, ref$prop_padj_BH, d_padj, ifelse(pass_padj, "PASS", "FAIL")))
  cat(sprintf("    delta   : %+.6f vs %+.6f   diff=%.3e   %s\n",
              re_delta, ref$prop_cliffs_delta, d_delta,
              ifelse(pass_delta, "PASS", "FAIL")))
  cat(sprintf("    HL      : %+.6f vs %+.6f   diff=%.3e   %s\n",
              re_HL, ref$prop_HL_estimate, d_HL, ifelse(pass_HL, "PASS", "FAIL")))

  test_summary_rows[[length(test_summary_rows) + 1]] <- data.frame(
    focal_id           = paste0("focal", k),
    unit               = unit_label,
    celltype           = focal_ct,
    n_disease_recomp   = length(xp_d),
    n_control_recomp   = length(xp_c),
    recomp_raw_p       = re_p,
    S5B_raw_p          = ref$prop_p_raw,
    abs_diff_raw_p     = d_p,
    PASS_raw_p         = pass_p,
    recomp_padj_BH     = re_padj,
    S5B_padj_BH        = ref$prop_padj_BH,
    abs_diff_padj_BH   = d_padj,
    PASS_padj_BH       = pass_padj,
    recomp_cliffs_delta = re_delta,
    S5B_cliffs_delta   = ref$prop_cliffs_delta,
    abs_diff_delta     = d_delta,
    PASS_delta         = pass_delta,
    recomp_HL_estimate = re_HL,
    S5B_HL_estimate    = ref$prop_HL_estimate,
    abs_diff_HL        = d_HL,
    PASS_HL            = pass_HL,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Step 6. Write Panel D outputs
# ============================================================
per_donor_df    <- do.call(rbind, per_donor_rows)
test_summary_df <- do.call(rbind, test_summary_rows)

per_out  <- file.path(OUTPUT_DIR,
                      "SupTable_S4_PanelD_manual_recomputation_per_donor.csv")
test_out <- file.path(OUTPUT_DIR,
                      "SupTable_S4_PanelD_manual_recomputation_test_summary.csv")
write.csv(per_donor_df, per_out, row.names = FALSE)
write.csv(test_summary_df, test_out, row.names = FALSE)
cat("\n[OK] Wrote per-donor recomputation:", per_out, "\n")
cat("[OK] Wrote test-level summary:",       test_out, "\n")

# ============================================================
# Step 7. Final PASS / FAIL verdict
# ============================================================
total_donor_rows <- nrow(per_donor_df)
passed_donor    <- sum(per_donor_df$PASS)
all_test_pass   <- all(test_summary_df$PASS_raw_p,
                       test_summary_df$PASS_padj_BH,
                       test_summary_df$PASS_delta,
                       test_summary_df$PASS_HL)

cat("\n", strrep("=", 75), "\n", sep = "")
cat("FINAL VERDICT\n")
cat(strrep("=", 75), "\n", sep = "")
cat(sprintf("Per-donor recomputation:   %d / %d donor-celltype rows PASS\n",
            passed_donor, total_donor_rows))
cat(sprintf("Test-level recomputation:  %s\n",
            ifelse(all_test_pass, "ALL 3 FOCAL TESTS PASS",
                   "SOME FOCAL TESTS FAIL — see test summary csv")))

if (passed_donor == total_donor_rows && all_test_pass) {
  cat("\n>>> Manual recomputation from cell-level counts is fully consistent\n")
  cat("    with S5A and S5B. AI-code verification PASSED.\n")
} else {
  cat("\n>>> Discrepancies detected. Investigate before publication.\n")
}

cat("\n[DONE] Script 04 v3 finished.\n")
