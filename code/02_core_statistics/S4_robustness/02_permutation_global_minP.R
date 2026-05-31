###############################################################################
## S4 / 02 - Family-wise permutation null, global min-P (Panel B)
## Supplementary Table S4 (robustness/sensitivity) | MN submission
##
## INPUT: S1A + S1B. OUTPUT: PanelB_*.csv (incl. 10,000-perm null).
## Run order: 01 -> 02 -> 03 -> 04 -> 05. Set REPO; outputs go to S4_robustness/.
###############################################################################
rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

SEED   <- 20251023
N_PERM <- 10000L
set.seed(SEED)

# ============================================================
# Path configuration (EDIT if needed)
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
# Step 2. Test ordering (5 powered units x 6 celltypes = 30)
# ============================================================
powered_units <- c("GSE157827", "GSE174367",
                   "syn52082747_AD", "syn52082747_FTD", "syn52082747_PSP")
celltypes <- c("Astrocytes", "Endothelial",
               "Excitatory neurons", "Inhibitory neurons",
               "Microglia", "Oligodendrocytes")
unit_disease <- c(GSE157827       = "AD",
                  GSE174367       = "AD",
                  syn52082747_AD  = "AD",
                  syn52082747_FTD = "FTD",
                  syn52082747_PSP = "PSP")
syn_units <- c("syn52082747_AD", "syn52082747_FTD", "syn52082747_PSP")

test_grid <- expand.grid(unit = powered_units, celltype = celltypes,
                         stringsAsFactors = FALSE)
test_grid <- test_grid[order(match(test_grid$unit, powered_units),
                             match(test_grid$celltype, celltypes)), ]
rownames(test_grid) <- NULL
stopifnot(nrow(test_grid) == 30L)

# ============================================================
# Step 3. Build donor blocks from ALL cell types with consistency check
# ============================================================
make_donor_block <- function(dat, expected_n = NULL, block_name = "") {
  # 1) Aggregate to one row per (donor, group) — checks no donor has >1 group
  block <- unique(dat[, c("donor", "group")])
  check <- aggregate(group ~ donor, data = block,
                     FUN = function(z) length(unique(z)))
  if (any(check$group != 1)) {
    bad <- check$donor[check$group != 1]
    stop("[", block_name, "] inconsistent group labels for donors: ",
         paste(bad, collapse = ", "))
  }
  block <- block[!duplicated(block$donor), ]
  block <- block[order(block$group, block$donor), ]
  rownames(block) <- NULL
  if (!is.null(expected_n)) stopifnot(nrow(block) == expected_n)
  block
}

syn_donor_block <- make_donor_block(
  s5a[s5a$unit %in% syn_units, ],
  expected_n = 40L, block_name = "syn52082747"
)
gse157827_donor_block <- make_donor_block(
  s5a[s5a$unit == "GSE157827", ],
  expected_n = 21L, block_name = "GSE157827"
)
gse174367_donor_block <- make_donor_block(
  s5a[s5a$unit == "GSE174367", ],
  expected_n = 18L, block_name = "GSE174367"
)

cat("\nBlock 1 (syn52082747 V1):\n");  print(table(syn_donor_block$group))
cat("Block 2 (GSE157827):\n");         print(table(gse157827_donor_block$group))
cat("Block 3 (GSE174367):\n");         print(table(gse174367_donor_block$group))

# ============================================================
# Step 4. Verify syn52082747 shared donor-celltype prop values
#         are identical across the 3 V1 units
# ============================================================
dup_check <- aggregate(
  prop_ubl3pos ~ donor + celltype,
  data = s5a[s5a$unit %in% syn_units, ],
  FUN  = function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
)
max_diff <- max(dup_check$prop_ubl3pos, na.rm = TRUE)
cat(sprintf("\n[Check] Max within-(donor,celltype) prop range across syn V1 units: %.3e\n",
            max_diff))
if (max_diff > 1e-12) {
  stop("Shared syn52082747 donor-celltype prop_ubl3pos values differ across units. ",
       "Cannot safely aggregate.")
}

# ============================================================
# Step 5. Build wide per-donor prop_ubl3pos tables
# ============================================================
syn_wide <- aggregate(prop_ubl3pos ~ donor + celltype,
                      data = s5a[s5a$unit %in% syn_units, ],
                      FUN  = function(x) x[1])
syn_prop <- reshape(syn_wide, timevar = "celltype", idvar = "donor",
                    direction = "wide")
names(syn_prop) <- sub("^prop_ubl3pos\\.", "", names(syn_prop))
stopifnot(nrow(syn_prop) == 40L)

gse157827_prop <- reshape(s5a[s5a$unit == "GSE157827",
                              c("donor", "celltype", "prop_ubl3pos")],
                          timevar = "celltype", idvar = "donor",
                          direction = "wide")
names(gse157827_prop) <- sub("^prop_ubl3pos\\.", "", names(gse157827_prop))
stopifnot(nrow(gse157827_prop) == 21L)

gse174367_prop <- reshape(s5a[s5a$unit == "GSE174367",
                              c("donor", "celltype", "prop_ubl3pos")],
                          timevar = "celltype", idvar = "donor",
                          direction = "wide")
names(gse174367_prop) <- sub("^prop_ubl3pos\\.", "", names(gse174367_prop))
stopifnot(nrow(gse174367_prop) == 18L)

# ============================================================
# Step 6. Wilcoxon helper and 30-test compute function
# ============================================================
wilcox_p <- function(x, y) {
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(x, y, exact = FALSE))$p.value
}

compute_30_p <- function(syn_lab, gse157_lab, gse174_lab) {
  p_vec <- numeric(nrow(test_grid))
  for (i in seq_len(nrow(test_grid))) {
    unit    <- test_grid$unit[i]
    ct      <- test_grid$celltype[i]
    disease <- unit_disease[unit]
    if (unit %in% syn_units) {
      labels <- syn_lab;       props <- syn_prop
    } else if (unit == "GSE157827") {
      labels <- gse157_lab;    props <- gse157827_prop
    } else if (unit == "GSE174367") {
      labels <- gse174_lab;    props <- gse174367_prop
    } else stop("Unknown unit: ", unit)

    dis_ids <- labels$donor[labels$group == disease]
    ctl_ids <- labels$donor[labels$group == "Control"]
    x <- props[match(dis_ids, props$donor), ct]; x <- x[!is.na(x)]
    y <- props[match(ctl_ids, props$donor), ct]; y <- y[!is.na(y)]
    p_vec[i] <- wilcox_p(x, y)
  }
  p_vec
}

# ============================================================
# Step 7. Observed P (sanity check vs S5B)
# ============================================================
observed_p <- compute_30_p(syn_donor_block,
                           gse157827_donor_block,
                           gse174367_donor_block)

s5b_observed_p <- numeric(nrow(test_grid))
for (i in seq_len(nrow(test_grid))) {
  row <- s5b[s5b$unit == test_grid$unit[i] &
             s5b$celltype == test_grid$celltype[i], ]
  s5b_observed_p[i] <- row$prop_p_raw
}
max_pdiff <- max(abs(observed_p - s5b_observed_p), na.rm = TRUE)
cat(sprintf("\n[Sanity] Max |recomputed obs p - S5B prop_p_raw|: %.3e\n", max_pdiff))
if (max_pdiff > 1e-6) {
  stop("Recomputed observed p does NOT match S5B prop_p_raw. ",
       "Investigate before trusting permutation results.")
}

# ============================================================
# Step 8. Permutation loop
# ============================================================
cat("\nRunning", N_PERM, "permutations...\n")
t0 <- Sys.time()

perm_p <- matrix(NA_real_, nrow = N_PERM, ncol = nrow(test_grid))

for (iter in seq_len(N_PERM)) {
  syn_perm <- syn_donor_block
  syn_perm$group <- sample(syn_perm$group)

  g157_perm <- gse157827_donor_block
  g157_perm$group <- sample(g157_perm$group)

  g174_perm <- gse174367_donor_block
  g174_perm$group <- sample(g174_perm$group)

  perm_p[iter, ] <- compute_30_p(syn_perm, g157_perm, g174_perm)

  if (iter %% 500 == 0) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    eta <- el / iter * (N_PERM - iter)
    cat(sprintf("  iter %d / %d   elapsed %.0fs   eta %.0fs\n",
                iter, N_PERM, el, eta))
  }
}
cat(sprintf("\nPermutation done. Total: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ============================================================
# Step 9. Empirical P values with (count + 1) / (N + 1) correction
# ============================================================
emp_per <- numeric(nrow(test_grid))
for (j in seq_along(emp_per)) {
  ok <- !is.na(perm_p[, j])
  emp_per[j] <- (sum(perm_p[ok, j] <= observed_p[j]) + 1) / (sum(ok) + 1)
}

null_minP <- apply(perm_p, 1, function(row) min(row, na.rm = TRUE))
ok_min    <- !is.na(null_minP)
emp_fw    <- numeric(nrow(test_grid))
for (j in seq_along(emp_fw)) {
  emp_fw[j] <- (sum(null_minP[ok_min] <= observed_p[j]) + 1) /
               (sum(ok_min) + 1)
}

# ============================================================
# Step 10. Assemble Panel B
# ============================================================
prop_padj_within <- numeric(nrow(test_grid))
prop_delta       <- numeric(nrow(test_grid))
prop_HL          <- numeric(nrow(test_grid))
for (i in seq_len(nrow(test_grid))) {
  row <- s5b[s5b$unit == test_grid$unit[i] &
             s5b$celltype == test_grid$celltype[i], ]
  prop_padj_within[i] <- row$prop_padj_BH
  prop_delta[i]       <- row$prop_cliffs_delta
  prop_HL[i]          <- row$prop_HL_estimate
}

panelB <- data.frame(
  unit                                  = test_grid$unit,
  disease                               = unname(unit_disease[test_grid$unit]),
  celltype                              = test_grid$celltype,
  observed_prop_p_raw                   = observed_p,
  observed_prop_padj_within_unit_BH     = prop_padj_within,
  observed_prop_cliffs_delta            = prop_delta,
  observed_prop_HL_estimate             = prop_HL,
  empirical_p_per_comparison            = emp_per,
  empirical_p_familyWise_globalMinP     = emp_fw,
  n_permutations                        = N_PERM,
  is_focal_within_unit_FDR_positive     = with(test_grid,
    (unit == "syn52082747_AD"  & celltype == "Excitatory neurons") |
    (unit == "syn52082747_PSP" & celltype == "Excitatory neurons") |
    (unit == "syn52082747_PSP" & celltype == "Inhibitory neurons")),
  stringsAsFactors = FALSE
)

out_path <- file.path(OUTPUT_DIR,
                      "SupTable_S4_PanelB_permutation_results.csv")
write.csv(panelB, out_path, row.names = FALSE)
cat("\n[OK] Wrote:", out_path, "\n")

null_df  <- data.frame(iteration = seq_len(N_PERM), null_min_p = null_minP)
null_out <- file.path(OUTPUT_DIR,
                      "SupTable_S4_PanelB_null_distribution_global_minP.csv")
write.csv(null_df, null_out, row.names = FALSE)
cat("[OK] Wrote null distribution:", null_out, "\n")

# ============================================================
# Step 11. Console summary
# ============================================================
cat("\n", strrep("=", 90), "\n", sep = "")
cat("Permutation results for three focal within-unit FDR-positive findings ",
    "(N_PERM = ", N_PERM, ")\n", sep = "")
cat(strrep("=", 90), "\n", sep = "")
pos_rows <- panelB[panelB$is_focal_within_unit_FDR_positive, ]
print(pos_rows[, c("unit", "celltype",
                   "observed_prop_p_raw",
                   "observed_prop_padj_within_unit_BH",
                   "observed_prop_cliffs_delta",
                   "empirical_p_per_comparison",
                   "empirical_p_familyWise_globalMinP")],
      row.names = FALSE, digits = 4)

cat("\nInterpretation:\n",
    "  - empirical_p_per_comparison supports the donor-label sensitivity\n",
    "    of THIS individual test (does NOT correct for multiple testing).\n",
    "  - empirical_p_familyWise_globalMinP is the permutation-based\n",
    "    family-wise (study-wide) error rate across 30 powered detection-\n",
    "    breadth tests; it does NOT prove or disprove the individual test.\n",
    sep = "")

cat("\n[DONE] Script 02 finished.\n")
