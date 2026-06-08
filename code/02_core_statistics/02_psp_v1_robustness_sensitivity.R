###############################################################################
# Supplementary Table S2 Route C robustness/sensitivity workbook
# v2 aligned to final Figure 4 v26 panels:
#   A multiple-testing sensitivity
#   B detection-count quasibinomial sensitivity
#   C leave-one-donor-out stability
#   D neuronal subcluster localization
###############################################################################

rm(list = ls()); gc()
Sys.setenv(LANG = "en")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(openxlsx)
})

base_dir <- "D:/RNA/2026063Molecular Neurodegeneration"
fig4_dir <- file.path(base_dir, "Figure4", "results", "RouteC_20260606_v26_2x2_axis_title_spacing")
out_dir  <- file.path(base_dir, "Supplementary_Table_S2", "results", "RouteC_20260606_robustness_sensitivity_aligned_to_Fig4_v26")
script_dir <- file.path(base_dir, "Supplementary_Table_S2", "R")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(name) {
  path <- file.path(fig4_dir, name)
  if (!file.exists(path)) stop("Missing Figure 4 source file: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

panelA <- read_csv("Figure4_panelA_multiple_testing_source.csv")
panelB <- read_csv("Figure4_panelB_detection_count_source.csv")
panelC_iter <- read_csv("Figure4_panelC_leave_one_out_source.csv")
panelC_sum  <- read_csv("Figure4_panelC_leave_one_out_summary.csv")
panelD <- read_csv("Figure4_panelD_subcluster_source.csv")

readme <- data.frame(
  Item = c(
    "Workbook title",
    "Companion figure",
    "Purpose",
    "Statistical unit",
    "Panel A",
    "Panel B",
    "Panel C",
    "Panel D",
    "Multiplicity note",
    "Source files",
    "Created"
  ),
  Description = c(
    "Supplementary Table S2. Robustness and sensitivity analyses for the PSP V1 cortical-neuron UBL3 candidate signal.",
    "Figure 4, Route C final v26.",
    "Source tables for robustness and sensitivity checks of PSP V1 excitatory- and inhibitory-neuron UBL3 detection-breadth findings.",
    "Donor is the statistical unit for Wilcoxon/Hodges-Lehmann analyses; detection-count models use donor-level binomial counts with quasibinomial dispersion.",
    "Multiple-testing sensitivity for the two PSP V1 cortical-neuron findings across increasingly broad correction families.",
    "Detection-count quasibinomial sensitivity using UBL3-positive nuclei counts and total nuclei per donor; reports odds ratios, 95% confidence intervals, raw P values and q values.",
    "Leave-one-donor-out donor-level robustness for the two PSP V1 cortical-neuron findings.",
    "Exploratory source-label neuronal subcluster localization within PSP V1 excitatory and inhibitory neurons; this does not replace the six-class primary framework.",
    "Within-unit q values refer to BH correction across the six major cell classes in the relevant disease-region unit unless otherwise stated; broader correction families are shown as sensitivity checks.",
    fig4_dir,
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  ),
  stringsAsFactors = FALSE
)

workbook_map <- data.frame(
  sheet_name = c("README", "Workbook_map", "PanelA_multiple_testing", "PanelB_detection_count", "PanelC_LOO_summary", "PanelC_LOO_iterations", "PanelD_subclusters", "Source_manifest"),
  figure_panel = c("", "", "Figure 4A", "Figure 4B", "Figure 4C", "Figure 4C", "Figure 4D", ""),
  source_file = c("", "", "Figure4_panelA_multiple_testing_source.csv", "Figure4_panelB_detection_count_source.csv", "Figure4_panelC_leave_one_out_summary.csv", "Figure4_panelC_leave_one_out_source.csv", "Figure4_panelD_subcluster_source.csv", ""),
  description = c(
    "Workbook description and definitions.",
    "Sheet map.",
    "Multiple-testing sensitivity results.",
    "Detection-count quasibinomial model results.",
    "Leave-one-donor-out summary.",
    "Leave-one-donor-out iteration-level results.",
    "Exploratory neuronal subcluster localization results.",
    "Input file paths and sizes."
  ),
  stringsAsFactors = FALSE
)

source_files <- c(
  "Figure4_panelA_multiple_testing_source.csv",
  "Figure4_panelB_detection_count_source.csv",
  "Figure4_panelC_leave_one_out_summary.csv",
  "Figure4_panelC_leave_one_out_source.csv",
  "Figure4_panelD_subcluster_source.csv",
  "Figure4_RouteC_count_subclusters_legend_draft.txt",
  "Figure4_RouteC_count_subclusters_file_sizes.csv"
)
source_manifest <- data.frame(
  file_name = source_files,
  path = file.path(fig4_dir, source_files),
  exists = file.exists(file.path(fig4_dir, source_files)),
  size_bytes = ifelse(file.exists(file.path(fig4_dir, source_files)), file.info(file.path(fig4_dir, source_files))$size, NA),
  stringsAsFactors = FALSE
)

wb <- createWorkbook()
add_sheet <- function(name, df, freeze = TRUE) {
  addWorksheet(wb, name)
  writeData(wb, name, df)
  if (freeze && nrow(df) > 0) freezePane(wb, name, firstActiveRow = 2)
  setColWidths(wb, name, cols = 1:max(1, ncol(df)), widths = "auto")
}

add_sheet("README", readme)
add_sheet("Workbook_map", workbook_map)
add_sheet("PanelA_multiple_testing", panelA)
add_sheet("PanelB_detection_count", panelB)
add_sheet("PanelC_LOO_summary", panelC_sum)
add_sheet("PanelC_LOO_iterations", panelC_iter)
add_sheet("PanelD_subclusters", panelD)
add_sheet("Source_manifest", source_manifest)

out_xlsx <- file.path(out_dir, "Supplementary_Table_S2_robustness_sensitivity_RouteC_aligned_to_Figure4_v26.xlsx")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)

writeLines(capture.output(sessionInfo()), file.path(script_dir, "sessionInfo_Supplementary_Table_S2_aligned_20260606_v2.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_Supplementary_Table_S2_aligned_20260606_v2.txt"))

cat("WROTE", out_xlsx, "\n")
cat("SHEETS", paste(names(wb), collapse=", "), "\n")
cat("SIZE_MB", sprintf("%.3f", file.info(out_xlsx)$size/1024^2), "\n")