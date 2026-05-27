# =============================================================================
# Build Supplementary Table S4 workbook for Molecular Neurodegeneration
#
# Output:
#   D:/RNA/supptable/SupTable_S4_robustness_sensitivity/result/MN/
#   Supplementary_Table_S4_robustness_sensitivity.xlsx
#
# Workbook structure:
#   1. README
#   2. PanelA_full
#   3. PanelA_focal_summary
#   4. PanelB_perm
#   5. PanelB_null
#   6. PanelC_summary
#   7. PanelC_iter
#   8. PanelD_donor
#   9. PanelD_test
# =============================================================================

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
Sys.setenv(LANG = "en")

# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------
required_pkgs <- c("readr", "dplyr", "openxlsx", "tibble", "stringr")

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(openxlsx)
  library(tibble)
  library(stringr)
})

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------
SRC_DIR <- "D:/RNA/supptable/SupTable_S4_robustness_sensitivity/result"

OUT_DIR <- "D:/RNA/supptable/SupTable_S4_robustness_sensitivity/result/MN"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_XLSX <- file.path(
  OUT_DIR,
  "Supplementary_Table_S4_robustness_sensitivity.xlsx"
)

# -----------------------------------------------------------------------------
# 2. Whether to include extra old/duplicate PanelA 3-positive summary
# -----------------------------------------------------------------------------
# Keep FALSE for the final nine-sheet submission workbook.
# If set TRUE, the workbook becomes 10 sheets and the legend should be changed
# from "Nine-sheet workbook" to "Ten-sheet workbook".
INCLUDE_EXTRA_PANELA_3POSITIVES <- FALSE

# -----------------------------------------------------------------------------
# 3. Sheet map
# -----------------------------------------------------------------------------
sheet_map <- tibble::tribble(
  ~sheet_name,             ~file_name,                                                   ~description,
  "PanelA_full",           "SupTable_S4_PanelA_alternative_FDR_corrections.csv",          "Panel A full multiple-testing sensitivity table for all 30 powered detection-breadth comparisons.",
  "PanelA_focal_summary",  "SupTable_S4_PanelA_summary_focal_findings.csv",              "Panel A focal-subset summary for the three syn52082747 V1 cortical-neuron findings.",
  "PanelB_perm",           "SupTable_S4_PanelB_permutation_results.csv",                 "Panel B donor-label permutation results, including per-comparison and family-wise empirical P values.",
  "PanelB_null",           "SupTable_S4_PanelB_null_distribution_global_minP.csv",       "Panel B global minimum P null distribution from donor-label permutations.",
  "PanelC_summary",        "SupTable_S4_PanelC_LOO_summary.csv",                         "Panel C leave-one-donor-out summary for the three focal findings.",
  "PanelC_iter",           "SupTable_S4_PanelC_LOO_per_iteration.csv",                   "Panel C leave-one-donor-out per-iteration donor-influence results.",
  "PanelD_donor",          "SupTable_S4_PanelD_manual_recomputation_per_donor.csv",      "Panel D independent manual recomputation per-donor table.",
  "PanelD_test",           "SupTable_S4_PanelD_manual_recomputation_test_summary.csv",   "Panel D independent manual recomputation per-test summary with PASS flags."
)

if (INCLUDE_EXTRA_PANELA_3POSITIVES) {
  sheet_map <- bind_rows(
    sheet_map,
    tibble::tibble(
      sheet_name  = "PanelA_3positives",
      file_name   = "SupTable_S4_PanelA_summary_3positives.csv",
      description = "Extra Panel A three-positive-finding summary. Include only if changing the workbook legend to ten sheets."
    )
  )
}

# -----------------------------------------------------------------------------
# 4. Check input files
# -----------------------------------------------------------------------------
sheet_map <- sheet_map %>%
  mutate(file_path = file.path(SRC_DIR, file_name))

missing_files <- sheet_map %>%
  filter(!file.exists(file_path))

if (nrow(missing_files) > 0) {
  message("Missing input files:")
  print(missing_files[, c("sheet_name", "file_name", "file_path")])
  stop("Some required CSV files are missing. Please check SRC_DIR.")
}

# -----------------------------------------------------------------------------
# 5. Read CSV files
# -----------------------------------------------------------------------------
read_csv_safe <- function(fp) {
  readr::read_csv(
    fp,
    show_col_types = FALSE,
    locale = readr::locale(encoding = "UTF-8"),
    guess_max = 100000,
    progress = FALSE
  )
}

data_list <- lapply(sheet_map$file_path, read_csv_safe)
names(data_list) <- sheet_map$sheet_name

dim_check <- tibble::tibble(
  sheet_name = names(data_list),
  n_rows = vapply(data_list, nrow, integer(1)),
  n_cols = vapply(data_list, ncol, integer(1)),
  source_file = sheet_map$file_name
)

message("\nInput table dimensions:")
print(dim_check)

# -----------------------------------------------------------------------------
# 6. Light sanity checks
# -----------------------------------------------------------------------------
check_expected_rows <- function(sheet, expected_n) {
  if (sheet %in% names(data_list)) {
    observed_n <- nrow(data_list[[sheet]])
    if (!identical(observed_n, expected_n)) {
      warning(sprintf(
        "%s has %d rows, but expected %d rows based on the planned legend.",
        sheet, observed_n, expected_n
      ))
    }
  }
}

check_expected_rows("PanelA_full", 30)
check_expected_rows("PanelD_donor", 62)
check_expected_rows("PanelD_test", 3)

# -----------------------------------------------------------------------------
# 7. Create README sheet content
# -----------------------------------------------------------------------------
readme_main <- tibble::tribble(
  ~Item, ~Description,
  "Supplementary table title",
  "Supplementary Table S4. Statistical robustness and sensitivity analyses for the three focal V1 cortical-neuron findings.",
  "Workbook purpose",
  "Nine-sheet workbook providing four orthogonal robustness checks for the three syn52082747 V1 cortical-neuron focal findings: AD excitatory neurons, PSP excitatory neurons, and PSP inhibitory neurons.",
  "Companion manuscript section",
  "Statistical robustness and sensitivity analyses.",
  "Companion figure",
  "Figure 4.",
  "Panel A",
  "Multiple-testing sensitivity. Five Benjamini-Hochberg / Bonferroni schemes are reported: within-unit BH6 primary, syn52082747 V1 dataset-level BH18, study-wide BH30, study-wide Bonferroni30, and combined BH60.",
  "Inferential set",
  "The 30 powered detection-breadth tests comprise five inferentially powered analytical units: GSE157827 AD MFG, GSE174367 AD PFC, syn52082747 AD V1, syn52082747 PSP V1, and syn52082747 FTD/PiD V1. Each contributes six harmonized cell types.",
  "Descriptive set",
  "The 12 syn21788402 EC and SFG comparisons with n = 3 vs 3 are descriptive bootstrap-CI mode and are not included in the powered inferential set.",
  "Panel B",
  "Donor-label permutation with 10,000 iterations. Two-level reporting includes per-comparison empirical P values, family-wise global min-P empirical P values, and syn52082747 four-group joint permutation preserving the shared-control structure.",
  "Panel C",
  "Leave-one-donor-out donor-influence analysis. Per-iteration disease-Control direction, Cliff's delta range, Hodges-Lehmann range, and BH-adjusted P-value range are reported for the three focal findings.",
  "Panel D",
  "Manual recomputation from the analysis-ready cell-level count layer generated in this study. Independent recomputation is provided per donor and per test, with PASS flags for raw P, FDR-adjusted P, Cliff's delta, and Hodges-Lehmann estimates against Supplementary Table S1/S5B within numerical tolerance.",
  "Data lineage scripts",
  "01_FDR_robustness; 02_permutation_global_minP; 03_LOO_donor_influence; 04_manual_recomputation_from_celllayer.",
  "Random seed",
  "SEED = 42 where applicable.",
  "Submission file",
  "Supplementary_Table_S4_robustness_sensitivity.xlsx.",
  "Created",
  as.character(Sys.time())
)

readme_sheet_map <- sheet_map %>%
  select(sheet_name, file_name, description)

# -----------------------------------------------------------------------------
# 8. Helper functions for formatting
# -----------------------------------------------------------------------------
calc_widths <- function(df, min_width = 8, max_width = 38, sample_n = 200) {
  widths <- vapply(seq_along(df), function(j) {
    txt <- c(names(df)[j], as.character(utils::head(df[[j]], sample_n)))
    txt[is.na(txt)] <- ""
    w <- max(nchar(txt, type = "width"), na.rm = TRUE) + 2
    w <- max(min_width, min(max_width, w))
    if (is.numeric(df[[j]])) {
      w <- min(w, 16)
    }
    w
  }, numeric(1))
  widths
}

safe_table_name <- function(sheet_name) {
  x <- gsub("[^A-Za-z0-9_]", "_", sheet_name)
  x <- paste0("tbl_", x)
  substr(x, 1, 255)
}

# -----------------------------------------------------------------------------
# 9. Build workbook
# -----------------------------------------------------------------------------
wb <- openxlsx::createWorkbook(creator = "Setou Lab")
openxlsx::modifyBaseFont(wb, fontName = "Arial", fontSize = 10)

style_title <- openxlsx::createStyle(
  fontName = "Arial",
  fontSize = 13,
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  halign = "left",
  valign = "center"
)

style_header <- openxlsx::createStyle(
  fontName = "Arial",
  fontSize = 10,
  fontColour = "#FFFFFF",
  textDecoration = "bold",
  fgFill = "#1F4E79",
  halign = "center",
  valign = "center",
  border = "Bottom",
  borderColour = "#808080"
)

style_body <- openxlsx::createStyle(
  fontName = "Arial",
  fontSize = 9,
  valign = "top"
)

style_wrap <- openxlsx::createStyle(
  fontName = "Arial",
  fontSize = 9,
  wrapText = TRUE,
  valign = "top"
)

# -----------------------------------------------------------------------------
# 10. README sheet
# -----------------------------------------------------------------------------
openxlsx::addWorksheet(wb, "README", gridLines = FALSE)

openxlsx::writeData(
  wb,
  "README",
  "Supplementary Table S4. Statistical robustness and sensitivity analyses",
  startRow = 1,
  startCol = 1
)

openxlsx::addStyle(
  wb,
  "README",
  style_title,
  rows = 1,
  cols = 1,
  gridExpand = TRUE
)

openxlsx::writeDataTable(
  wb,
  "README",
  readme_main,
  startRow = 3,
  startCol = 1,
  tableStyle = "TableStyleMedium2",
  tableName = "tbl_README_main",
  withFilter = FALSE
)

start_row_map <- nrow(readme_main) + 6

openxlsx::writeData(
  wb,
  "README",
  "Workbook sheet map",
  startRow = start_row_map,
  startCol = 1
)

openxlsx::addStyle(
  wb,
  "README",
  style_title,
  rows = start_row_map,
  cols = 1,
  gridExpand = TRUE
)

openxlsx::writeDataTable(
  wb,
  "README",
  readme_sheet_map,
  startRow = start_row_map + 2,
  startCol = 1,
  tableStyle = "TableStyleMedium9",
  tableName = "tbl_README_sheet_map",
  withFilter = TRUE
)

openxlsx::writeData(
  wb,
  "README",
  "Input table dimensions",
  startRow = start_row_map + nrow(readme_sheet_map) + 5,
  startCol = 1
)

openxlsx::addStyle(
  wb,
  "README",
  style_title,
  rows = start_row_map + nrow(readme_sheet_map) + 5,
  cols = 1,
  gridExpand = TRUE
)

openxlsx::writeDataTable(
  wb,
  "README",
  dim_check,
  startRow = start_row_map + nrow(readme_sheet_map) + 7,
  startCol = 1,
  tableStyle = "TableStyleMedium4",
  tableName = "tbl_README_dimensions",
  withFilter = TRUE
)

openxlsx::setColWidths(wb, "README", cols = 1, widths = 28)
openxlsx::setColWidths(wb, "README", cols = 2, widths = 110)
openxlsx::setColWidths(wb, "README", cols = 3:5, widths = 32)
openxlsx::addStyle(
  wb,
  "README",
  style_wrap,
  rows = 1:200,
  cols = 1:5,
  gridExpand = TRUE,
  stack = TRUE
)
openxlsx::freezePane(wb, "README", firstActiveRow = 3)

# -----------------------------------------------------------------------------
# 11. Data sheets
# -----------------------------------------------------------------------------
for (sh in names(data_list)) {
  df <- data_list[[sh]]
  
  openxlsx::addWorksheet(wb, sh, gridLines = FALSE)
  
  openxlsx::writeDataTable(
    wb,
    sh,
    df,
    startRow = 1,
    startCol = 1,
    tableStyle = "TableStyleMedium2",
    tableName = safe_table_name(sh),
    withFilter = TRUE
  )
  
  if (ncol(df) > 0) {
    openxlsx::addStyle(
      wb,
      sh,
      style_header,
      rows = 1,
      cols = 1:ncol(df),
      gridExpand = TRUE,
      stack = TRUE
    )
    
    openxlsx::addStyle(
      wb,
      sh,
      style_body,
      rows = 2:(nrow(df) + 1),
      cols = 1:ncol(df),
      gridExpand = TRUE,
      stack = TRUE
    )
    
    widths <- calc_widths(df)
    openxlsx::setColWidths(wb, sh, cols = 1:ncol(df), widths = widths)
  }
  
  openxlsx::freezePane(wb, sh, firstActiveRow = 2)
}

# -----------------------------------------------------------------------------
# 12. Save workbook
# -----------------------------------------------------------------------------
if (file.exists(OUT_XLSX)) {
  file.remove(OUT_XLSX)
}

openxlsx::saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)

message("\nDONE.")
message("Workbook saved to:")
message(OUT_XLSX)

message("\nWorkbook sheets:")
print(names(wb))

message("\nFinal sheet count: ", length(names(wb)))
if (!INCLUDE_EXTRA_PANELA_3POSITIVES && length(names(wb)) == 9) {
  message("PASS: Workbook has 9 sheets, consistent with the planned legend.")
} else if (INCLUDE_EXTRA_PANELA_3POSITIVES && length(names(wb)) == 10) {
  message("NOTE: Workbook has 10 sheets. Please change the legend from 'Nine-sheet workbook' to 'Ten-sheet workbook'.")
} else {
  warning("Unexpected sheet count. Please check workbook structure.")
}