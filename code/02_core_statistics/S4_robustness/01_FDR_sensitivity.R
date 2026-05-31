###############################################################################
## S4 / 01 - Multiple-testing sensitivity (Panel A)
## Supplementary Table S4 (robustness/sensitivity) | MN submission
##
## INPUT: S1B detection-breadth table. OUTPUT: PanelA_*.csv.
## Run order: 01 -> 02 -> 03 -> 04 -> 05. Set REPO; outputs go to S4_robustness/.
###############################################################################
rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ============================================================
# Path configuration (EDIT if needed)
# ============================================================
## === 改这里：代码库根目录 ===
REPO <- "D:/RNA/Code/UBL3_tauopathy"

INPUT_S5B  <- file.path(REPO, "output", "stats_tables", "S1_detection_breadth", "S1B_donor_level_test.csv")
OUTPUT_DIR <- file.path(REPO, "output", "stats_tables", "S4_robustness")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Step 1. Load and validate S5B
# ============================================================
if (!file.exists(INPUT_S5B)) stop("S5B not found at: ", INPUT_S5B)
s5b <- read.csv(INPUT_S5B, stringsAsFactors = FALSE)

required_cols <- c("unit", "disease", "region", "celltype",
                   "prop_p_raw", "prop_padj_BH", "prop_test_status",
                   "expr_p_raw", "expr_padj_BH", "expr_test_status",
                   "prop_cliffs_delta",
                   "prop_HL_estimate", "prop_HL_ci95_low", "prop_HL_ci95_high")
missing <- setdiff(required_cols, names(s5b))
if (length(missing) > 0) {
  stop("S5B is missing required columns: ", paste(missing, collapse = ", "))
}

cat("Loaded S5B:", nrow(s5b), "rows (expected 42)\n")
stopifnot(nrow(s5b) == 42L)

# ============================================================
# Step 2. Define inferential set and validate
# ============================================================
unpowered_units <- c("syn21788402_EC", "syn21788402_SFG")
syn_v1_units    <- c("syn52082747_AD", "syn52082747_FTD", "syn52082747_PSP")
s5b$in_inferential_set <- !(s5b$unit %in% unpowered_units)

idx_prop_powered <- which(s5b$in_inferential_set &
                          s5b$prop_test_status == "YES_wilcox")
idx_expr_powered <- which(s5b$in_inferential_set &
                          s5b$expr_test_status == "YES_wilcox")
idx_syn18        <- which(s5b$unit %in% syn_v1_units &
                          s5b$prop_test_status == "YES_wilcox")

cat("Powered detection-breadth tests:",  length(idx_prop_powered), "(expect 30)\n")
cat("Powered conditional-expression tests:", length(idx_expr_powered), "(expect 30)\n")
cat("syn52082747 V1 detection-breadth tests:", length(idx_syn18), "(expect 18)\n")

stopifnot(length(idx_prop_powered) == 30L)
stopifnot(length(idx_expr_powered) == 30L)
stopifnot(length(idx_syn18)        == 18L)

stopifnot(!anyNA(s5b$prop_p_raw[idx_prop_powered]))
stopifnot(!anyNA(s5b$expr_p_raw[idx_expr_powered]))

# ============================================================
# Step 3. Build sensitivity table
# ============================================================
out <- data.frame(
  unit                                  = s5b$unit,
  disease                               = s5b$disease,
  region                                = s5b$region,
  celltype                              = s5b$celltype,
  in_inferential_set                    = ifelse(s5b$in_inferential_set,
                                                 "YES", "NO_descriptive_only"),

  # Effect size (carried from S5B for one-stop view)
  prop_cliffs_delta                     = s5b$prop_cliffs_delta,
  prop_HL_estimate                      = s5b$prop_HL_estimate,
  prop_HL_ci95_low                      = s5b$prop_HL_ci95_low,
  prop_HL_ci95_high                     = s5b$prop_HL_ci95_high,

  # Detection-breadth (prop)
  prop_p_raw                            = s5b$prop_p_raw,
  prop_padj_within_unit_BH_primary      = s5b$prop_padj_BH,
  prop_padj_syn52082747_V1_BH18         = NA_real_,
  prop_padj_syn52082747_V1_bonferroni18 = NA_real_,
  prop_padj_global_BH_across30          = NA_real_,
  prop_padj_bonferroni_across30         = NA_real_,

  # Conditional-expression (expr)
  expr_p_raw                            = s5b$expr_p_raw,
  expr_padj_within_unit_BH_primary      = s5b$expr_padj_BH,
  expr_padj_global_BH_across30          = NA_real_,
  expr_padj_bonferroni_across30         = NA_real_,

  # Combined 60 (worst-case stress test)
  combined_prop_padj_global_BH_60       = NA_real_,
  combined_prop_padj_bonferroni_60      = NA_real_,
  combined_expr_padj_global_BH_60       = NA_real_,
  combined_expr_padj_bonferroni_60      = NA_real_,

  stringsAsFactors = FALSE
)

# (i) syn52082747 V1 dataset-level BH18 (intermediate sensitivity)
out$prop_padj_syn52082747_V1_BH18[idx_syn18] <-
  p.adjust(s5b$prop_p_raw[idx_syn18], method = "BH")
out$prop_padj_syn52082747_V1_bonferroni18[idx_syn18] <-
  p.adjust(s5b$prop_p_raw[idx_syn18], method = "bonferroni")

# (ii)-(iii) global BH/Bonferroni across 30 detection-breadth
out$prop_padj_global_BH_across30[idx_prop_powered] <-
  p.adjust(s5b$prop_p_raw[idx_prop_powered], method = "BH")
out$prop_padj_bonferroni_across30[idx_prop_powered] <-
  p.adjust(s5b$prop_p_raw[idx_prop_powered], method = "bonferroni")

# conditional-expression
out$expr_padj_global_BH_across30[idx_expr_powered] <-
  p.adjust(s5b$expr_p_raw[idx_expr_powered], method = "BH")
out$expr_padj_bonferroni_across30[idx_expr_powered] <-
  p.adjust(s5b$expr_p_raw[idx_expr_powered], method = "bonferroni")

# (iv)-(v) Combined 60-test family
p_combined   <- c(s5b$prop_p_raw[idx_prop_powered],
                  s5b$expr_p_raw[idx_expr_powered])
padj_BH_60   <- p.adjust(p_combined, method = "BH")
padj_bonf_60 <- p.adjust(p_combined, method = "bonferroni")

n_prop <- length(idx_prop_powered)
out$combined_prop_padj_global_BH_60[idx_prop_powered]   <- padj_BH_60[1:n_prop]
out$combined_prop_padj_bonferroni_60[idx_prop_powered]  <- padj_bonf_60[1:n_prop]
out$combined_expr_padj_global_BH_60[idx_expr_powered]   <- padj_BH_60[(n_prop+1):length(p_combined)]
out$combined_expr_padj_bonferroni_60[idx_expr_powered]  <- padj_bonf_60[(n_prop+1):length(p_combined)]

# ============================================================
# Step 4. Write Panel A
# ============================================================
out_path <- file.path(OUTPUT_DIR,
                      "SupTable_S4_PanelA_alternative_FDR_corrections.csv")
write.csv(out, out_path, row.names = FALSE)
cat("\n[OK] Wrote:", out_path, "\n")

# ============================================================
# Step 5. Focal-finding summary
# (Three focal within-unit FDR-positive findings, NOT prespecified-positive)
# ============================================================
focal <- list(
  c("syn52082747_AD",  "Excitatory neurons"),
  c("syn52082747_PSP", "Excitatory neurons"),
  c("syn52082747_PSP", "Inhibitory neurons")
)

summary_rows <- list()
cat("\n", strrep("=", 90), "\n", sep = "")
cat("Detection-breadth padj for the three focal findings under all schemes:\n")
cat(strrep("=", 90), "\n", sep = "")

for (pos in focal) {
  row <- out[out$unit == pos[1] & out$celltype == pos[2], ]
  cat(sprintf("\n%s / %s   (Cliff delta=%+.3f, HL=%+.4f [%+.4f, %+.4f])\n",
              pos[1], pos[2],
              row$prop_cliffs_delta,
              row$prop_HL_estimate,
              row$prop_HL_ci95_low, row$prop_HL_ci95_high))
  cat(sprintf("  raw p:                              %.4f\n", row$prop_p_raw))
  cat(sprintf("  within-unit BH6 (PRIMARY):          %.4f%s\n",
              row$prop_padj_within_unit_BH_primary,
              ifelse(row$prop_padj_within_unit_BH_primary < 0.05, "  *", "")))
  cat(sprintf("  syn52082747 V1 BH18 (intermediate): %.4f%s\n",
              row$prop_padj_syn52082747_V1_BH18,
              ifelse(row$prop_padj_syn52082747_V1_BH18 < 0.05, "  *", "")))
  cat(sprintf("  global BH30 (conservative):         %.4f%s\n",
              row$prop_padj_global_BH_across30,
              ifelse(row$prop_padj_global_BH_across30 < 0.05, "  *", "")))
  cat(sprintf("  Bonferroni x30:                     %.4f%s\n",
              row$prop_padj_bonferroni_across30,
              ifelse(row$prop_padj_bonferroni_across30 < 0.05, "  *", "")))
  cat(sprintf("  combined BH60 (worst-case):         %.4f%s\n",
              row$combined_prop_padj_global_BH_60,
              ifelse(row$combined_prop_padj_global_BH_60 < 0.05, "  *", "")))
  cat(sprintf("  combined Bonferroni x60:            %.4f%s\n",
              row$combined_prop_padj_bonferroni_60,
              ifelse(row$combined_prop_padj_bonferroni_60 < 0.05, "  *", "")))

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    unit                            = pos[1],
    celltype                        = pos[2],
    prop_cliffs_delta               = row$prop_cliffs_delta,
    prop_HL_estimate                = row$prop_HL_estimate,
    prop_HL_ci95_low                = row$prop_HL_ci95_low,
    prop_HL_ci95_high               = row$prop_HL_ci95_high,
    raw_p                           = row$prop_p_raw,
    within_unit_BH6_primary         = row$prop_padj_within_unit_BH_primary,
    syn52082747_V1_BH18             = row$prop_padj_syn52082747_V1_BH18,
    syn52082747_V1_bonferroni18     = row$prop_padj_syn52082747_V1_bonferroni18,
    global_BH30                     = row$prop_padj_global_BH_across30,
    Bonferroni30                    = row$prop_padj_bonferroni_across30,
    combined_BH60                   = row$combined_prop_padj_global_BH_60,
    combined_Bonferroni60           = row$combined_prop_padj_bonferroni_60,
    stringsAsFactors = FALSE
  )
}

summary_df  <- do.call(rbind, summary_rows)
summary_out <- file.path(OUTPUT_DIR,
                         "SupTable_S4_PanelA_summary_focal_findings.csv")
write.csv(summary_df, summary_out, row.names = FALSE)
cat("\n[OK] Wrote summary:", summary_out, "\n")

# Survival counts under each scheme
cat("\nSurvival count (padj < 0.05) of detection-breadth tests:\n")
for (col in c("prop_padj_within_unit_BH_primary",
              "prop_padj_syn52082747_V1_BH18",
              "prop_padj_global_BH_across30",
              "prop_padj_bonferroni_across30",
              "combined_prop_padj_global_BH_60",
              "combined_prop_padj_bonferroni_60")) {
  n_sig <- sum(out[[col]] < 0.05, na.rm = TRUE)
  total <- sum(!is.na(out[[col]]))
  cat(sprintf("  %-45s: %d / %d tests\n", col, n_sig, total))
}

cat("\n[DONE] Script 01 finished.\n")
