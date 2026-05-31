###############################################################################
## S4 / 03 - Leave-one-donor-out stability (Panel C)
## Supplementary Table S4 (robustness/sensitivity) | MN submission
##
## INPUT: S1A + S1B. OUTPUT: PanelC_*.csv.
## Run order: 01 -> 02 -> 03 -> 04 -> 05. Set REPO; outputs go to S4_robustness/.
###############################################################################
rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)
SEED <- 20251023
set.seed(SEED)

# ============================================================
# Path configuration
# ============================================================
## === 改这里：代码库根目录 ===
REPO <- "D:/RNA/Code/UBL3_tauopathy"

INPUT_S5A  <- file.path(REPO, "output", "stats_tables", "S1_detection_breadth", "S1A_per_donor_data.csv")
INPUT_S5B  <- file.path(REPO, "output", "stats_tables", "S1_detection_breadth", "S1B_donor_level_test.csv")
OUTPUT_DIR <- file.path(REPO, "output", "stats_tables", "S4_robustness")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Step 1. Load data
# ============================================================
if (!file.exists(INPUT_S5A)) stop("S5A not found: ", INPUT_S5A)
if (!file.exists(INPUT_S5B)) stop("S5B not found: ", INPUT_S5B)
s5a <- read.csv(INPUT_S5A, stringsAsFactors = FALSE)
s5b <- read.csv(INPUT_S5B, stringsAsFactors = FALSE)
cat("Loaded S5A:", nrow(s5a), "rows\n")
cat("Loaded S5B:", nrow(s5b), "rows\n")

# ============================================================
# Step 2. Helpers
# ============================================================
cliffs_delta <- function(x, y) {
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  d <- outer(x, y, "-")
  (sum(d > 0) - sum(d < 0)) / (length(x) * length(y))
}
hodges_lehmann <- function(x, y) {
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  median(as.vector(outer(x, y, "-")))
}
wilcox_p <- function(x, y) {
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(x, y, exact = FALSE))$p.value
}

celltypes <- c("Astrocytes", "Endothelial",
               "Excitatory neurons", "Inhibitory neurons",
               "Microglia", "Oligodendrocytes")

# ============================================================
# Step 3. Three focal findings
# ============================================================
focal <- data.frame(
  unit            = c("syn52082747_AD", "syn52082747_PSP", "syn52082747_PSP"),
  disease         = c("AD", "PSP", "PSP"),
  focal_celltype  = c("Excitatory neurons", "Excitatory neurons",
                      "Inhibitory neurons"),
  stringsAsFactors = FALSE
)

# ============================================================
# Step 4. LOO loop
# ============================================================
per_iter_rows <- list()
summary_rows  <- list()

for (k in seq_len(nrow(focal))) {
  unit_k    <- focal$unit[k]
  disease_k <- focal$disease[k]
  focal_k   <- focal$focal_celltype[k]

  cat("\n--------------------------------------------------\n")
  cat(sprintf("Focal #%d: %s / %s\n", k, unit_k, focal_k))
  cat("--------------------------------------------------\n")

  unit_data <- s5a[s5a$unit == unit_k, ]
  donors_in_unit <- unique(unit_data[, c("donor", "group")])
  donors_in_unit <- donors_in_unit[order(donors_in_unit$group,
                                         donors_in_unit$donor), ]
  cat("  Donors in unit:", nrow(donors_in_unit),
      "  (groups:", paste(names(table(donors_in_unit$group)),
                         table(donors_in_unit$group), sep = "=",
                         collapse = ", "), ")\n")

  i_focal <- which(celltypes == focal_k)

  # ---- Baseline (no LOO) ----
  baseline_p     <- numeric(length(celltypes))
  baseline_delta <- numeric(length(celltypes))
  baseline_HL    <- numeric(length(celltypes))
  for (ci in seq_along(celltypes)) {
    ct  <- celltypes[ci]
    sub <- unit_data[unit_data$celltype == ct, ]
    x   <- sub$prop_ubl3pos[sub$group == disease_k]
    y   <- sub$prop_ubl3pos[sub$group == "Control"]
    baseline_p[ci]     <- wilcox_p(x, y)
    baseline_delta[ci] <- cliffs_delta(x, y)
    baseline_HL[ci]    <- hodges_lehmann(x, y)
  }
  baseline_padj <- p.adjust(baseline_p, method = "BH")

  # Baseline check vs S5B (HARD STOP on mismatch)
  ref <- s5b[s5b$unit == unit_k & s5b$celltype == focal_k, ]
  if (nrow(ref) != 1L) {
    stop("Cannot uniquely find focal row in S5B: ", unit_k, " / ", focal_k)
  }
  if (abs(baseline_p[i_focal] - ref$prop_p_raw) > 1e-6) {
    stop("Baseline raw p does not match S5B prop_p_raw for ",
         unit_k, " / ", focal_k,
         " (recomputed=", baseline_p[i_focal],
         "; S5B=", ref$prop_p_raw, ")")
  }
  if (abs(baseline_padj[i_focal] - ref$prop_padj_BH) > 1e-6) {
    stop("Baseline within-unit BH does not match S5B for ",
         unit_k, " / ", focal_k)
  }
  if (abs(baseline_delta[i_focal] - ref$prop_cliffs_delta) > 1e-6) {
    stop("Baseline Cliff's delta does not match S5B for ",
         unit_k, " / ", focal_k)
  }
  cat(sprintf(
    "  Baseline (no LOO):  raw p=%.4f  within-BH=%.4f  delta=%+.3f  HL=%+.4f  [match S5B]\n",
    baseline_p[i_focal], baseline_padj[i_focal],
    baseline_delta[i_focal], baseline_HL[i_focal]))

  # ---- LOO iterations ----
  for (d_idx in seq_len(nrow(donors_in_unit))) {
    drop_donor <- donors_in_unit$donor[d_idx]
    drop_group <- donors_in_unit$group[d_idx]
    loo_data <- unit_data[unit_data$donor != drop_donor, ]

    p_vec     <- numeric(length(celltypes))
    delta_vec <- numeric(length(celltypes))
    HL_vec    <- numeric(length(celltypes))
    n_dis_vec <- integer(length(celltypes))
    n_ctl_vec <- integer(length(celltypes))
    for (ci in seq_along(celltypes)) {
      ct  <- celltypes[ci]
      sub <- loo_data[loo_data$celltype == ct, ]
      x   <- sub$prop_ubl3pos[sub$group == disease_k]
      y   <- sub$prop_ubl3pos[sub$group == "Control"]
      p_vec[ci]     <- wilcox_p(x, y)
      delta_vec[ci] <- cliffs_delta(x, y)
      HL_vec[ci]    <- hodges_lehmann(x, y)
      n_dis_vec[ci] <- length(x)
      n_ctl_vec[ci] <- length(y)
    }
    padj_vec <- p.adjust(p_vec, method = "BH")

    per_iter_rows[[length(per_iter_rows) + 1]] <- data.frame(
      positive_id              = paste0("focal", k),
      unit                     = unit_k,
      focal_celltype           = focal_k,
      dropped_donor            = drop_donor,
      dropped_donor_group      = drop_group,
      n_disease                = n_dis_vec[i_focal],
      n_control                = n_ctl_vec[i_focal],
      raw_p                    = p_vec[i_focal],
      padj_within_unit_BH      = padj_vec[i_focal],
      cliffs_delta             = delta_vec[i_focal],
      HL_estimate              = HL_vec[i_focal],
      direction_positive_HL    = HL_vec[i_focal]    > 0,
      direction_positive_delta = delta_vec[i_focal] > 0,
      stringsAsFactors = FALSE
    )
  }

  # ---- Summary for this focal ----
  this_iters <- do.call(rbind, per_iter_rows)
  this_iters <- this_iters[this_iters$positive_id == paste0("focal", k), ]

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    positive_id                  = paste0("focal", k),
    unit                         = unit_k,
    focal_celltype               = focal_k,
    n_LOO_iterations             = nrow(this_iters),

    # Direction stability — report BOTH metrics
    pct_direction_positive_HL    = 100 *
      mean(this_iters$direction_positive_HL,    na.rm = TRUE),
    pct_direction_positive_delta = 100 *
      mean(this_iters$direction_positive_delta, na.rm = TRUE),

    # Effect size stability — PRIMARY decision criterion
    cliffs_delta_min             = min(this_iters$cliffs_delta, na.rm = TRUE),
    cliffs_delta_median          = median(this_iters$cliffs_delta, na.rm = TRUE),
    cliffs_delta_max             = max(this_iters$cliffs_delta, na.rm = TRUE),
    HL_estimate_min              = min(this_iters$HL_estimate, na.rm = TRUE),
    HL_estimate_median           = median(this_iters$HL_estimate, na.rm = TRUE),
    HL_estimate_max              = max(this_iters$HL_estimate, na.rm = TRUE),

    # p / padj fluctuation — reported for completeness, NOT decision criterion
    raw_p_min                    = min(this_iters$raw_p, na.rm = TRUE),
    raw_p_median                 = median(this_iters$raw_p, na.rm = TRUE),
    raw_p_max                    = max(this_iters$raw_p, na.rm = TRUE),
    padj_BH_min                  = min(this_iters$padj_within_unit_BH, na.rm = TRUE),
    padj_BH_median               = median(this_iters$padj_within_unit_BH, na.rm = TRUE),
    padj_BH_max                  = max(this_iters$padj_within_unit_BH, na.rm = TRUE),
    n_padj_lt_0_05               = sum(this_iters$padj_within_unit_BH < 0.05,
                                        na.rm = TRUE),
    pct_padj_lt_0_05             = 100 *
      mean(this_iters$padj_within_unit_BH < 0.05, na.rm = TRUE),

    stringsAsFactors = FALSE
  )

  cat(sprintf("  LOO direction positive (HL>0):    %d / %d  (%.1f%%)\n",
              sum(this_iters$direction_positive_HL, na.rm = TRUE),
              nrow(this_iters),
              100 * mean(this_iters$direction_positive_HL, na.rm = TRUE)))
  cat(sprintf("  LOO direction positive (delta>0): %d / %d  (%.1f%%)\n",
              sum(this_iters$direction_positive_delta, na.rm = TRUE),
              nrow(this_iters),
              100 * mean(this_iters$direction_positive_delta, na.rm = TRUE)))
  cat(sprintf("  Cliff's delta range: [%+.3f, %+.3f]   median %+.3f\n",
              min(this_iters$cliffs_delta, na.rm = TRUE),
              max(this_iters$cliffs_delta, na.rm = TRUE),
              median(this_iters$cliffs_delta, na.rm = TRUE)))
  cat(sprintf("  within-unit BH padj range: [%.4f, %.4f]   median %.4f   (%d/%d iters padj<0.05)\n",
              min(this_iters$padj_within_unit_BH, na.rm = TRUE),
              max(this_iters$padj_within_unit_BH, na.rm = TRUE),
              median(this_iters$padj_within_unit_BH, na.rm = TRUE),
              sum(this_iters$padj_within_unit_BH < 0.05, na.rm = TRUE),
              nrow(this_iters)))
}

# ============================================================
# Step 5. Write outputs
# ============================================================
per_iter_df <- do.call(rbind, per_iter_rows)
summary_df  <- do.call(rbind, summary_rows)

iter_out <- file.path(OUTPUT_DIR,
                      "SupTable_S4_PanelC_LOO_per_iteration.csv")
sum_out  <- file.path(OUTPUT_DIR,
                      "SupTable_S4_PanelC_LOO_summary.csv")
write.csv(per_iter_df, iter_out, row.names = FALSE)
write.csv(summary_df,  sum_out,  row.names = FALSE)
cat("\n[OK] Wrote:", iter_out, "\n")
cat("[OK] Wrote:", sum_out, "\n")

cat("\nDecision criterion (per design):\n",
    "  - PRIMARY: stability of direction (HL>0 / delta>0) and effect size\n",
    "    magnitude (Cliff's delta) across LOO iterations.\n",
    "  - Iteration-level padj fluctuation is reported for completeness\n",
    "    only; padj will increase when N decreases, so iteration-level\n",
    "    padj is NOT the LOO decision criterion.\n",
    sep = "")

cat("\n[DONE] Script 03 finished.\n")
