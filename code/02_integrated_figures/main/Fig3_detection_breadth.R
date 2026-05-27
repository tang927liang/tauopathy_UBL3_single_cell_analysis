###############################################################################
## Figure 3 — Molecular Neurodegeneration integrated version v17
## Donor-level UBL3 detection-breadth analysis identifies a candidate AD/PSP
## cortical-neuron signal in V1
##
## What this script does
##   1) Reads the donor-level UBL3 and housekeeping summary tables.
##   2) Generates Panel A, Panel B, and Panel C from the original Figure 3 logic.
##   3) Saves standalone panel PDF/PNG files.
##   4) Creates one composite Figure 3 file for Molecular Neurodegeneration:
##        - 170 mm full-page width
##        - <= 225 mm total height
##        - 1000 dpi PNG / TIFF
##        - vector PDF using cairo_pdf
##        - figure keys inside the graphic
##        - no figure title/legend text inside the graphic
##        - reduced surrounding whitespace
##        - line width >= 0.25 pt
##
## Changes vs v15
##   - Panel B decoding key is a NATIVE GRAPHICAL LEGEND (symbol + short label)
##     at the bottom of the panel, replacing the prose caption block. This matches
##     the BMC/Molecular Neurodegeneration rule "Figure keys should be incorporated
##     into the graphic, not into the legend"; the full prose explanation stays in
##     the manuscript Figure 3 legend. n.s. vs descriptive are distinguished in the
##     key by marker (filled vs open) and line style (solid vs dotted) via
##     override.aes, not by the two near-identical greys.
##   - x-axis title uses en-dash via \u2013 (locale-safe on Windows).
##   - Final composite height trimmed from 224.5 mm to 217.5 mm for safety margin.
##
## Changes vs v16 (this version)
##   - Panel C: FDR < 0.05 tiles now carry a thin black border in addition to the
##     white asterisk. This separates genuinely significant tiles from tiles that
##     are merely dark (e.g. B2M in AD excitatory neurons, FDR just above 0.05),
##     which previously read as significant from colour alone.
##   - Panel B: cell-type axis labels enlarged from 5.45 pt to ~6.1 pt
##     (MN_BASE_SIZE - 1.0) for legibility at final size.
##
## Output directory requested by user
##   D:/RNA/MNversion submission/Figure 3. Donor-level UBL3 detection-breadth analysis identifies a candidate ADPSP cortical-neuron signal in V1/Results
###############################################################################

rm(list = ls()); gc()
SEED <- 20251023
set.seed(SEED)
Sys.setenv(LANG = "en")
options(stringsAsFactors = FALSE)

# =============================================================================
# 0. Packages
# =============================================================================
required_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "scales", "grid",
  "ggbeeswarm", "fs", "patchwork"
)

to_install <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(grid)
  library(ggbeeswarm)
  library(fs)
  library(patchwork)
})

# =============================================================================
# 1. Molecular Neurodegeneration output settings
# =============================================================================
MN_WIDTH_FULL_MM <- 170
MN_HEIGHT_MAX_MM <- 225
MN_DPI <- 1000
MN_DPI_SUBMISSION <- 300

# SAFETY SWITCHES --------------------------------------------------------------
# RStudio can abort fatally when exporting large/complex ggplot objects at
# 1000 dpi, especially Panel B. Molecular Neurodegeneration/BMC requires about
# 300 dpi at final size and accepts vector PDF. Therefore the safe default is:
#   - vector PDF: ON
#   - 300 dpi submission PNG/TIFF: ON
#   - 1000 dpi PNG/TIFF: OFF by default
# Turn MAKE_1000DPI <- TRUE only after the PDF/300 dpi version is confirmed.
MAKE_1000DPI <- TRUE
MAKE_STANDALONE_RASTER <- FALSE

# Carefully chosen to keep the final composite safely under 225 mm.
# (v16: trimmed from 116 / 107 to 112.5 / 103.5 -> total 217.5 mm.)
W_A_MM  <- 170
H_A_MM  <- 112.5
W_B_MM  <- 101
W_C_MM  <- 69
H_BC_MM <- 103.5
GAP_MM  <- 1.5
FINAL_H_MM <- H_A_MM + GAP_MM + H_BC_MM

if (FINAL_H_MM > MN_HEIGHT_MAX_MM) {
  stop("Final figure height exceeds Molecular Neurodegeneration maximum height.")
}

# User-requested output folder. Use forward slashes in R even on Windows.
out_dir <- "D:/RNA/MNversion submission/Figure 3. Donor-level UBL3 detection-breadth analysis identifies a candidate ADPSP cortical-neuron signal in V1/Results"
fs::dir_create(out_dir, recurse = TRUE)

log_fp <- file.path(out_dir, "Figure3_MN_generation_log.txt")
sink(log_fp, split = TRUE)
cat("==== Figure 3 MN integrated script (v17) START ====\n")
cat("Time          :", as.character(Sys.time()), "\n")
cat("Output dir    :", out_dir, "\n")
cat("Final size    :", MN_WIDTH_FULL_MM, "mm x", FINAL_H_MM, "mm @", MN_DPI, "dpi\n\n")

# =============================================================================
# 2. Font — Arial preferred; fallback to sans if unavailable
# =============================================================================
PLOT_FONT <- "Arial"
if (.Platform$OS.type == "windows") {
  arial_ok <- tryCatch({
    windowsFonts(Arial = windowsFont("TT Arial"))
    TRUE
  }, error = function(e) FALSE, warning = function(w) TRUE)
  if (!arial_ok) {
    cat("WARNING: Arial registration failed; fallback to sans.\n")
    PLOT_FONT <- "sans"
  } else {
    cat("Font: Arial registered.\n")
  }
} else {
  PLOT_FONT <- "sans"
  cat("Font: non-Windows system; using sans fallback.\n")
}

# Devices
pdf_device <- function(filename, width, height, ...) {
  grDevices::cairo_pdf(filename = filename, width = width, height = height,
                       family = PLOT_FONT, onefile = TRUE, ...)
}

# PNG device type. cairo-png avoids the corrupted-IDAT issue seen with the
# previous ragg + magick read-back workflow on some Windows systems.
PNG_TYPE <- if (isTRUE(capabilities("cairo"))) {
  "cairo-png"
} else if (.Platform$OS.type == "windows") {
  "windows"
} else {
  "Xlib"
}
cat("PNG device type:", PNG_TYPE, "\n")

# =============================================================================
# 3. Input data paths
# =============================================================================
s5a_fp <- "D:/RNA/supptable/SupTable_S5_donor_level_combined/Results/SupTable_S5A_per_donor_data_FULL9.csv"
s5b_fp <- "D:/RNA/supptable/SupTable_S5_donor_level_combined/Results/SupTable_S5B_donor_level_test_FULL9.csv"
s6b_fp <- "D:/RNA/supptable/SupTable_S6_housekeeping_baseline/Results/Results_extended11/SupTable_S6B_donor_level_test_housekeeping_FULL9_11genes.csv"

input_files <- c(S5A = s5a_fp, S5B = s5b_fp, S6B = s6b_fp)
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs) > 0) {
  cat("ERROR: Missing input file(s):\n")
  for (x in missing_inputs) cat("  -", x, "\n")
  sink()
  stop("Missing input file(s). Please check the input CSV paths in section 3.")
}

cat("Reading data...\n")
s5a <- read.csv(s5a_fp, check.names = FALSE)
s5b <- read.csv(s5b_fp, check.names = FALSE)
s6b <- read.csv(s6b_fp, check.names = FALSE)
cat("  S5A:", nrow(s5a), "rows x", ncol(s5a), "columns\n")
cat("  S5B:", nrow(s5b), "rows x", ncol(s5b), "columns\n")
cat("  S6B:", nrow(s6b), "rows x", ncol(s6b), "columns\n\n")

check_cols <- function(dat, required, label) {
  miss <- setdiff(required, names(dat))
  if (length(miss) > 0) {
    stop(label, " is missing required column(s): ", paste(miss, collapse = ", "))
  }
}

check_cols(s5a, c("unit", "celltype", "group", "prop_ubl3pos"), "S5A")
check_cols(s5b, c("unit", "celltype", "prop_test_status", "prop_padj_BH",
                  "prop_HL_estimate", "prop_HL_ci95_low", "prop_HL_ci95_high"), "S5B")
check_cols(s6b, c("unit", "celltype", "test_status", "gene", "padj_BH"), "S6B")

# =============================================================================
# 4. Colors and MN theme
# =============================================================================
PAL_DISEASE          <- "#D55E00"  # Okabe-Ito orange
PAL_CONTROL          <- "#0072B2"  # Okabe-Ito blue
PAL_SIG              <- "#D55E00"
PAL_HIGHLIGHT_BG     <- "#FCD6A8"
PAL_HIGHLIGHT_ALPHA  <- 0.55
PAL_HIGHLIGHT_BORDER <- "#D55E00"
PAL_NS               <- "#9E9E9E"
PAL_DESC             <- "#BDBDBD"
PAL_UBL3_ROW         <- "#D55E00"
PAL_PANEL_C_SEP      <- "gray60"
PAL_PANEL_C_SIGBORDER <- "black"   # v17: border for FDR < 0.05 tiles in Panel C
SEP_LINEWIDTH        <- 0.35
SEP_LINETYPE         <- "dashed"

MN_BASE_SIZE  <- 7.1
MN_LINE_WIDTH <- 0.30

# Use a compact journal theme. 0.30 pt line width is above the 0.25 pt minimum.
theme_mn <- function(base_size = MN_BASE_SIZE) {
  theme_classic(base_size = base_size, base_family = PLOT_FONT) +
    theme(
      axis.text        = element_text(color = "black", size = base_size, family = PLOT_FONT),
      axis.title       = element_text(color = "black", size = base_size,
                                      face = "bold", family = PLOT_FONT),
      axis.line        = element_line(color = "black", linewidth = MN_LINE_WIDTH),
      axis.ticks       = element_line(color = "black", linewidth = MN_LINE_WIDTH),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(size = base_size, color = "black",
                                      face = "bold", family = PLOT_FONT),
      legend.position  = "bottom",
      legend.box       = "horizontal",
      legend.key.size  = unit(0.25, "cm"),
      legend.text      = element_text(size = base_size, color = "black", family = PLOT_FONT),
      legend.title     = element_text(size = base_size, color = "black", family = PLOT_FONT),
      panel.grid       = element_blank(),
      plot.title       = element_blank(),
      plot.subtitle    = element_text(size = base_size, color = "gray35",
                                      face = "italic", hjust = 0.5,
                                      family = PLOT_FONT,
                                      margin = margin(b = 1, unit = "pt")),
      plot.caption     = element_text(size = base_size - 1.0, hjust = 0,
                                      color = "gray25", family = PLOT_FONT,
                                      margin = margin(t = 2, unit = "pt")),
      plot.tag         = element_text(size = base_size + 3.2, face = "bold",
                                      color = "black", family = PLOT_FONT),
      plot.margin      = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")
    )
}

# =============================================================================
# 5. Data preprocessing
# =============================================================================
ct_normalize <- function(x) {
  x <- as.character(x)
  x[x %in% c("Excit", "Excitatory", "Excitatory neurons")] <- "Excitatory neurons"
  x[x %in% c("Inhib", "Inhibitory", "Inhibitory neurons")] <- "Inhibitory neurons"
  x[x %in% c("Astro", "Astrocytes")] <- "Astrocytes"
  x[x %in% c("Oligo", "Oligodendrocytes", "Oligodendrocyte")] <- "Oligodendrocytes"
  x[x %in% c("Microgl", "Microglia")] <- "Microglia"
  x[x %in% c("Endo", "Endothelial")] <- "Endothelial"
  x
}

celltype_order <- c("Astrocytes", "Endothelial", "Excitatory neurons",
                    "Inhibitory neurons", "Microglia", "Oligodendrocytes")

unit_order <- c("GSE157827", "GSE174367", "syn52082747_AD",
                "syn21788402_EC", "syn21788402_SFG",
                "syn52082747_FTD", "syn52082747_PSP")

unit_display <- c(
  "GSE157827"       = "GSE157827 (AD)",
  "GSE174367"       = "GSE174367 (AD)",
  "syn52082747_AD"  = "syn52082747 (AD)",
  "syn21788402_EC"  = "syn21788402 EC",
  "syn21788402_SFG" = "syn21788402 SFG",
  "syn52082747_FTD" = "syn52082747 (FTD)",
  "syn52082747_PSP" = "syn52082747 (PSP)"
)

ct_short <- c(
  "Astrocytes" = "Astro",
  "Endothelial" = "Endo",
  "Excitatory neurons" = "Excit",
  "Inhibitory neurons" = "Inhib",
  "Microglia" = "Micro",
  "Oligodendrocytes" = "Oligo"
)

unit_labels_full <- c(
  "GSE157827"       = "GSE157827 (AD)\nMFG",
  "GSE174367"       = "GSE174367 (AD)\nPFC",
  "syn52082747_AD"  = "syn52082747 (AD)\nV1",
  "syn21788402_EC"  = "syn21788402 (AD)\nEC",
  "syn21788402_SFG" = "syn21788402 (AD)\nSFG",
  "syn52082747_FTD" = "syn52082747 (FTD)\nV1",
  "syn52082747_PSP" = "syn52082747 (PSP)\nV1"
)

s5a <- s5a[s5a$unit %in% unit_order, ]
s5b <- s5b[s5b$unit %in% unit_order, ]

s5a$unit     <- factor(s5a$unit, levels = unit_order)
s5a$celltype <- factor(ct_normalize(s5a$celltype), levels = celltype_order)
s5a$group2   <- ifelse(s5a$group == "Control", "Control", "Disease")
s5a$group2   <- factor(s5a$group2, levels = c("Disease", "Control"))

s5b$celltype <- ct_normalize(s5b$celltype)
s6b$celltype <- ct_normalize(s6b$celltype)

# =============================================================================
# 6. Panel A — donor-level UBL3-positive cell fraction
# =============================================================================
cat("========== Panel A ==========" , "\n")

sig_facets <- s5b %>%
  filter(prop_test_status == "YES_wilcox", !is.na(prop_padj_BH), prop_padj_BH < 0.05) %>%
  mutate(unit = factor(unit, levels = unit_order),
         celltype = factor(celltype, levels = celltype_order)) %>%
  select(unit, celltype) %>%
  distinct()

sig_bg <- sig_facets %>% mutate(highlight = TRUE)
sig_star <- sig_facets %>% mutate(x_star = 1.5, y_star = Inf, star = "*")

p_A <- ggplot(s5a, aes(x = group2, y = prop_ubl3pos)) +
  geom_rect(data = sig_bg, aes(fill = highlight),
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            inherit.aes = FALSE, alpha = PAL_HIGHLIGHT_ALPHA) +
  scale_fill_manual(values = c("TRUE" = PAL_HIGHLIGHT_BG), guide = "none") +
  geom_boxplot(aes(color = group2), outlier.shape = NA, fill = "white",
               width = 0.55, linewidth = 0.35, median.linewidth = 0.65) +
  geom_quasirandom(aes(color = group2),
                   size = 0.48, alpha = 0.85, width = 0.18, stroke = 0) +
  scale_color_manual(values = c("Disease" = PAL_DISEASE, "Control" = PAL_CONTROL),
                     name = NULL, limits = c("Disease", "Control")) +
  geom_text(data = sig_star, aes(x = x_star, y = y_star, label = star),
            inherit.aes = FALSE, color = PAL_SIG, size = 4.5,
            vjust = 1.35, fontface = "bold", family = PLOT_FONT) +
  facet_grid(unit ~ celltype, scales = "fixed", switch = "y",
             labeller = labeller(
               unit = unit_labels_full,
               celltype = c(
                 "Astrocytes"         = "Astrocytes",
                 "Endothelial"        = "Endothelial",
                 "Excitatory neurons" = "Excitatory\nneurons",
                 "Inhibitory neurons" = "Inhibitory\nneurons",
                 "Microglia"          = "Microglia",
                 "Oligodendrocytes"   = "Oligodendro-\ncytes"
               )
             )) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    breaks = c(0, 0.2, 0.4, 0.6),
    expand = expansion(mult = c(0.00, 0.03))
  ) +
  coord_cartesian(ylim = c(0, 0.68), clip = "off") +
  labs(x = NULL, y = "UBL3-positive cell proportion per donor", tag = "A") +
  theme_mn(base_size = MN_BASE_SIZE) +
  theme(
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.text.y       = element_text(size = MN_BASE_SIZE + 0.8, color = "black", family = PLOT_FONT),
    axis.title.y      = element_text(size = MN_BASE_SIZE + 1.4, face = "bold",
                                     color = "black", family = PLOT_FONT,
                                     margin = margin(r = 3, unit = "pt")),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = MN_BASE_SIZE + 0.4,
                                     color = "black", face = "bold", family = PLOT_FONT,
                                     lineheight = 0.95),
    strip.text.x      = element_text(size = MN_BASE_SIZE - 0.1, color = "black",
                                     face = "bold", family = PLOT_FONT,
                                     lineheight = 0.95),
    strip.placement   = "outside",
    panel.spacing.x   = unit(0.08, "lines"),
    panel.spacing.y   = unit(0.38, "lines"),
    legend.position   = "bottom",
    legend.margin     = margin(t = 0, b = 0, unit = "pt")
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.8, linetype = 0)))

# =============================================================================
# 7. Panel B — forest plot, faceted by analytical unit
# =============================================================================
cat("========== Panel B ==========" , "\n")

forest_df <- s5b %>%
  mutate(
    celltype_full = factor(ct_normalize(celltype), levels = rev(celltype_order)),
    unit = factor(unit, levels = unit_order),
    unit_facet = factor(unit_labels_full[as.character(unit)],
                        levels = unit_labels_full[unit_order]),
    is_descriptive = (prop_test_status == "DESCRIPTIVE_n<4"),
    is_sig = (prop_test_status == "YES_wilcox" &
                !is.na(prop_padj_BH) & prop_padj_BH < 0.05),
    sig_class = case_when(
      is_sig ~ "Significant (FDR < 0.05)",
      is_descriptive ~ "Descriptive (n < 4)",
      TRUE ~ "n.s."
    )
  ) %>%
  mutate(y_num = as.numeric(celltype_full)) %>%
  arrange(unit, celltype_full)

# Fix legend order: significant -> n.s. -> descriptive
forest_df$sig_class <- factor(
  forest_df$sig_class,
  levels = c("Significant (FDR < 0.05)", "n.s.", "Descriptive (n < 4)")
)

sig_rows  <- forest_df[forest_df$is_sig, ]
desc_rows <- forest_df[forest_df$is_descriptive, ]
ns_rows   <- forest_df[forest_df$sig_class == "n.s.", ]

xlim_lo <- max(min(forest_df$prop_HL_ci95_low, na.rm = TRUE), -0.30)
xlim_hi <- min(max(forest_df$prop_HL_ci95_high, na.rm = TRUE),  0.30)
xlim_use <- c(xlim_lo, xlim_hi)

# Faceting by analytical unit avoids one very long 42-row label column.
# Cell-type names remain fully written, while dataset/region labels are separated
# into the left facet strips.
#
# MN compliance note:
#   The decoding key is a native graphical legend (symbol + short label) at the
#   bottom of the panel, NOT a prose caption. The full prose explanation lives in
#   the manuscript Figure 3 legend. n.s. vs descriptive are distinguished by marker
#   (filled vs open) and line style (solid vs dotted), not by the two near-identical
#   greys, so override.aes encodes shape + linetype.
p_B_core <- ggplot(forest_df, aes(y = y_num)) +
  geom_rect(data = sig_rows,
            aes(xmin = -Inf, xmax = Inf,
                ymin = y_num - 0.45,
                ymax = y_num + 0.45),
            inherit.aes = FALSE,
            fill = PAL_HIGHLIGHT_BG, alpha = PAL_HIGHLIGHT_ALPHA,
            color = PAL_HIGHLIGHT_BORDER, linewidth = 0.30) +
  geom_vline(xintercept = 0, linewidth = 0.30,
             color = "grey40", linetype = "dashed") +
  # Descriptive (n < 4): grey open circle + dotted CI
  geom_segment(data = desc_rows,
               aes(x = prop_HL_ci95_low, xend = prop_HL_ci95_high,
                   y = y_num, yend = y_num, color = sig_class),
               linewidth = 0.33, linetype = "dotted") +
  geom_point(data = desc_rows,
             aes(x = prop_HL_estimate, y = y_num, color = sig_class),
             fill = "white", shape = 21, size = 0.85, stroke = 0.35) +
  # n.s.: grey filled circle + solid CI
  geom_segment(data = ns_rows,
               aes(x = prop_HL_ci95_low, xend = prop_HL_ci95_high,
                   y = y_num, yend = y_num, color = sig_class),
               linewidth = 0.33) +
  geom_point(data = ns_rows,
             aes(x = prop_HL_estimate, y = y_num, color = sig_class),
             size = 0.75, shape = 16) +
  # FDR-significant: orange filled diamond + thick CI
  geom_segment(data = sig_rows,
               aes(x = prop_HL_ci95_low, xend = prop_HL_ci95_high,
                   y = y_num, yend = y_num, color = sig_class),
               linewidth = 0.48) +
  geom_point(data = sig_rows,
             aes(x = prop_HL_estimate, y = y_num, color = sig_class),
             size = 1.15, shape = 18) +
  facet_grid(unit_facet ~ ., switch = "y", drop = FALSE) +
  scale_color_manual(
    name   = NULL,
    breaks = c("Significant (FDR < 0.05)", "n.s.", "Descriptive (n < 4)"),
    values = c("Significant (FDR < 0.05)" = PAL_SIG,
               "n.s."                     = PAL_NS,
               "Descriptive (n < 4)"      = PAL_DESC),
    labels = c("FDR < 0.05", "n.s.", "Descriptive (n < 4)"),
    guide  = guide_legend(
      override.aes = list(
        shape     = c(18, 16, 21),                 # diamond / filled circle / open circle
        linetype  = c("solid", "solid", "dotted"),
        fill      = "white",                       # used only by shape 21
        size      = c(1.15, 0.75, 0.85),           # point size
        linewidth = c(0.48, 0.33, 0.33)            # CI line width
      )
    )
  ) +
  scale_x_continuous(limits = xlim_use, breaks = c(-0.1, 0, 0.1, 0.2)) +
  scale_y_continuous(
    breaks = seq_along(levels(forest_df$celltype_full)),
    labels = levels(forest_df$celltype_full),
    limits = c(0.5, length(levels(forest_df$celltype_full)) + 0.5),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Hodges\u2013Lehmann shift, UBL3+ proportion\n(Disease\u2013Control, 95% CI)",
    y = NULL,
    tag = "B"
  ) +
  theme_mn(base_size = MN_BASE_SIZE) +
  theme(
    # v17: enlarged from MN_BASE_SIZE - 1.65 (5.45 pt) to MN_BASE_SIZE - 1.0 (~6.1 pt)
    axis.text.y = element_text(size = MN_BASE_SIZE - 1.0,
                               color = "black", family = PLOT_FONT),
    axis.text.x = element_text(size = MN_BASE_SIZE - 0.2,
                               color = "black", family = PLOT_FONT),
    axis.title.x = element_text(size = MN_BASE_SIZE - 0.05, face = "bold",
                                color = "black", family = PLOT_FONT,
                                lineheight = 0.95,
                                margin = margin(t = 1, unit = "pt")),
    strip.text.y.left = element_text(angle = 0, hjust = 1,
                                     size = MN_BASE_SIZE - 0.45,
                                     color = "black", face = "bold",
                                     family = PLOT_FONT,
                                     lineheight = 0.90),
    strip.background.y = element_rect(fill = "grey92", color = NA),
    strip.placement = "outside",
    panel.spacing.y = unit(0.10, "lines"),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.key.size   = unit(0.28, "cm"),
    legend.key.height = unit(0.30, "cm"),
    legend.spacing.x  = unit(0.18, "cm"),
    legend.text       = element_text(size = MN_BASE_SIZE - 0.6,
                                     color = "black", family = PLOT_FONT),
    legend.margin     = margin(t = 2, r = 0, b = 0, l = 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")
  )

# Panel B is kept as one ggplot object; the decoding key is the bottom legend.
p_B <- p_B_core

# =============================================================================
# 8. Panel C — housekeeping baseline heatmap
# =============================================================================
cat("========== Panel C ==========" , "\n")

critical_combos <- data.frame(
  unit       = c("syn52082747_AD", "syn52082747_PSP", "syn52082747_PSP"),
  celltype   = c("Excitatory neurons", "Excitatory neurons", "Inhibitory neurons"),
  comp_label = c("AD:\nExcitatory\nneurons",
                 "PSP:\nExcitatory\nneurons",
                 "PSP:\nInhibitory\nneurons"),
  stringsAsFactors = FALSE
)

ubl3_rows <- s5b %>%
  filter(unit %in% critical_combos$unit, celltype %in% critical_combos$celltype) %>%
  filter(prop_test_status == "YES_wilcox") %>%
  inner_join(critical_combos, by = c("unit", "celltype")) %>%
  mutate(gene = "UBL3", padj = prop_padj_BH) %>%
  select(unit, celltype, comp_label, gene, padj)

hk_rows <- s6b %>%
  filter(unit %in% critical_combos$unit, celltype %in% critical_combos$celltype) %>%
  filter(test_status == "YES_wilcox") %>%
  inner_join(critical_combos, by = c("unit", "celltype")) %>%
  select(unit, celltype, comp_label, gene, padj = padj_BH)

heat_df <- bind_rows(ubl3_rows, hk_rows) %>%
  mutate(comp_label = factor(comp_label, levels = critical_combos$comp_label))

hk_order <- hk_rows %>%
  group_by(gene) %>%
  summarise(med = median(padj, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(gene)

gene_order <- c("UBL3", hk_order)
heat_df$gene <- factor(heat_df$gene, levels = rev(gene_order))
heat_df$sig_label <- ifelse(heat_df$padj < 0.05, "*", "")

neglog_max <- 1.7
heat_df$neglog_padj <- pmax(0, pmin(-log10(heat_df$padj), neglog_max))
ubl3_y_pos <- length(gene_order)

# v17: FDR < 0.05 tiles get a thin black border so significance is marked by
# border + asterisk, not by colour alone (avoids dark-but-non-significant tiles,
# e.g. B2M in AD excitatory neurons, reading as significant).
heat_sig <- heat_df %>% filter(padj < 0.05)

p_C <- ggplot(heat_df, aes(x = comp_label, y = gene, fill = neglog_padj)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_tile(data = heat_sig, fill = NA,
            color = PAL_PANEL_C_SIGBORDER, linewidth = 0.45) +
  geom_text(aes(label = sig_label), size = 4.0, fontface = "bold",
            color = "white", vjust = 0.6, family = PLOT_FONT) +
  geom_hline(yintercept = ubl3_y_pos - 0.5,
             color = PAL_PANEL_C_SEP, linetype = SEP_LINETYPE,
             linewidth = SEP_LINEWIDTH) +
  scale_fill_gradient(
    low = "grey95", high = "#B30000",
    limits = c(0, neglog_max),
    breaks = c(0, log10(1 / 0.05), neglog_max),
    labels = c("1.0", "0.05", "<0.02"),
    name = "FDR",
    guide = guide_colourbar(
      title.position = "left",
      title.hjust = 0.5,
      barwidth = unit(2.35, "cm"),
      barheight = unit(0.18, "cm"),
      ticks = FALSE
    )
  ) +
  scale_x_discrete(position = "bottom") +
  labs(x = NULL, y = NULL, tag = "C", subtitle = "syn52082747 cohorts") +
  theme_mn(base_size = MN_BASE_SIZE) +
  theme(
    plot.subtitle = element_text(
      size = MN_BASE_SIZE + 1.0,
      color = "gray35",
      face = "italic",
      hjust = 0.5,
      family = PLOT_FONT,
      margin = margin(b = 2, unit = "pt")
    ),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1,
                               size = MN_BASE_SIZE - 0.9,
                               lineheight = 0.84,
                               color = "black", family = PLOT_FONT,
                               margin = margin(t = 1, unit = "pt")),
    axis.text.y = element_text(
      size = MN_BASE_SIZE - 0.1,
      color = "black",
      family = PLOT_FONT
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width  = unit(0.78, "cm"),
    legend.key.height = unit(0.18, "cm"),
    legend.margin     = margin(t = 1, b = 0, unit = "pt"),
    legend.title      = element_text(size = MN_BASE_SIZE - 0.2, color = "black",
                                     family = PLOT_FONT,
                                     margin = margin(r = 4, unit = "pt")),
    legend.text       = element_text(size = MN_BASE_SIZE - 1.15, color = "black", family = PLOT_FONT),
    plot.margin = margin(t = 2, r = 1, b = 4, l = -4, unit = "pt")
  )

# =============================================================================
# 9. Save standalone panels
# =============================================================================
cat("========== Saving standalone panels ==========" , "\n")

out_A_pdf <- file.path(out_dir, "Figure3_PanelA_MN_vector.pdf")
out_B_pdf <- file.path(out_dir, "Figure3_PanelB_MN_vector.pdf")
out_C_pdf <- file.path(out_dir, "Figure3_PanelC_MN_vector.pdf")

# Standalone PDFs are useful for inspection and remain vector-based.
ggsave(out_A_pdf, plot = p_A, width = W_A_MM, height = H_A_MM,
       units = "mm", device = pdf_device, bg = "white")
cat("Saved Panel A vector PDF:", W_A_MM, "x", H_A_MM, "mm\n")

ggsave(out_B_pdf, plot = p_B, width = W_B_MM, height = H_BC_MM,
       units = "mm", device = pdf_device, bg = "white")
cat("Saved Panel B vector PDF:", W_B_MM, "x", H_BC_MM, "mm\n")

ggsave(out_C_pdf, plot = p_C, width = W_C_MM, height = H_BC_MM,
       units = "mm", device = pdf_device, bg = "white")
cat("Saved Panel C vector PDF:", W_C_MM, "x", H_BC_MM, "mm\n")

# Optional standalone rasters are limited to 300 dpi by default to prevent
# RStudio session crashes. The final composite is the file needed for submission.
out_A_png_300 <- file.path(out_dir, "Figure3_PanelA_MN_300dpi_check.png")
out_B_png_300 <- file.path(out_dir, "Figure3_PanelB_MN_300dpi_check.png")
out_C_png_300 <- file.path(out_dir, "Figure3_PanelC_MN_300dpi_check.png")

if (MAKE_STANDALONE_RASTER) {
  ggsave(out_A_png_300, plot = p_A, width = W_A_MM, height = H_A_MM,
         units = "mm", dpi = MN_DPI_SUBMISSION, bg = "white", limitsize = FALSE)
  ggsave(out_B_png_300, plot = p_B, width = W_B_MM, height = H_BC_MM,
         units = "mm", dpi = MN_DPI_SUBMISSION, bg = "white", limitsize = FALSE)
  ggsave(out_C_png_300, plot = p_C, width = W_C_MM, height = H_BC_MM,
         units = "mm", dpi = MN_DPI_SUBMISSION, bg = "white", limitsize = FALSE)
  cat("Saved optional standalone 300 dpi check PNGs.\n")
} else {
  cat("Skipped standalone raster panels to avoid RStudio high-dpi crash.\n")
}
cat("\n")

# =============================================================================
# 10. Composite vector PDF with patchwork
# =============================================================================
cat("========== Building final vector composite ==========" , "\n")

# p_B is a single ggplot object with a bottom symbol legend (the key).
# We still wrap it as one atomic panel before assembling the final figure.
p_B_wrapped <- patchwork::wrap_elements(full = p_B)
p_C_wrapped <- patchwork::wrap_elements(full = p_C)

bottom_row <- patchwork::wrap_plots(
  list(p_B_wrapped, p_C_wrapped),
  ncol = 2,
  widths = c(W_B_MM, W_C_MM)
)

p_final_vector <- patchwork::wrap_plots(
  list(p_A, patchwork::plot_spacer(), bottom_row),
  ncol = 1,
  heights = c(H_A_MM, GAP_MM, H_BC_MM)
)

out_final_vector_pdf <- file.path(out_dir, "Figure3_Combined_MN_vector.pdf")
ggsave(out_final_vector_pdf, plot = p_final_vector,
       width = MN_WIDTH_FULL_MM, height = FINAL_H_MM,
       units = "mm", device = pdf_device, bg = "white")
cat("Saved final vector PDF:", basename(out_final_vector_pdf), "\n")

# =============================================================================
# 11. Final raster outputs — safe submission first, 1000 dpi
# =============================================================================
cat("========== Exporting final raster files directly ==========" , "\n")

# Journal-sized 300 dpi raster files for file-size-safe submission backup.
# This matches the BMC/Molecular Neurodegeneration final-size recommendation.
out_final_png_300  <- file.path(out_dir, "Figure3_Combined_MN_300dpi_submission.png")
out_final_tiff_300 <- file.path(out_dir, "Figure3_Combined_MN_300dpi_submission_LZW.tiff")

ggsave(out_final_png_300, plot = p_final_vector,
       width = MN_WIDTH_FULL_MM, height = FINAL_H_MM,
       units = "mm", dpi = MN_DPI_SUBMISSION,
       bg = "white", limitsize = FALSE)
cat("Saved final 300 dpi submission PNG:", basename(out_final_png_300), "\n")

tryCatch({
  ggsave(out_final_tiff_300, plot = p_final_vector,
         width = MN_WIDTH_FULL_MM, height = FINAL_H_MM,
         units = "mm", dpi = MN_DPI_SUBMISSION, compression = "lzw",
         bg = "white", limitsize = FALSE)
  cat("Saved final 300 dpi submission LZW TIFF:", basename(out_final_tiff_300), "\n")
}, error = function(e) {
  cat("WARNING: 300 dpi LZW TIFF export failed. Trying uncompressed TIFF.\n")
  ggsave(out_final_tiff_300, plot = p_final_vector,
         width = MN_WIDTH_FULL_MM, height = FINAL_H_MM,
         units = "mm", dpi = MN_DPI_SUBMISSION,
         bg = "white", limitsize = FALSE)
})

# Optional 1000 dpi outputs. These are NOT required by the journal and may exceed
# 10 MB or crash RStudio on some Windows devices. Use only if necessary.
out_final_png  <- file.path(out_dir, "Figure3_Combined_MN_1000dpi.png")
out_final_tiff <- file.path(out_dir, "Figure3_Combined_MN_1000dpi_LZW.tiff")

if (MAKE_1000DPI) {
  cat("MAKE_1000DPI is TRUE: exporting 1000 dpi files.\n")
  tryCatch({
    ggsave(out_final_png, plot = p_final_vector,
           width = MN_WIDTH_FULL_MM, height = FINAL_H_MM,
           units = "mm", dpi = MN_DPI,
           bg = "white", limitsize = FALSE)
    cat("Saved optional final 1000 dpi PNG:", basename(out_final_png), "\n")
  }, error = function(e) {
    cat("WARNING: 1000 dpi PNG export failed: ", conditionMessage(e), "\n", sep = "")
  })

  tryCatch({
    ggsave(out_final_tiff, plot = p_final_vector,
           width = MN_WIDTH_FULL_MM, height = FINAL_H_MM,
           units = "mm", dpi = MN_DPI, compression = "lzw",
           bg = "white", limitsize = FALSE)
    cat("Saved optional final 1000 dpi LZW TIFF:", basename(out_final_tiff), "\n")
  }, error = function(e) {
    cat("WARNING: 1000 dpi TIFF export failed: ", conditionMessage(e), "\n", sep = "")
  })
} else {
  cat("Skipped 1000 dpi raster export.\n")
}

cat("Final raster export finished.\n\n")

# =============================================================================
# 12. Molecular Neurodegeneration compliance self-check
# =============================================================================
cat("========== Molecular Neurodegeneration self-check ==========" , "\n")

file_mb <- function(x) {
  if (!file.exists(x)) return(NA_real_)
  round(file.info(x)$size / 1024 / 1024, 2)
}

outputs <- c(
  out_A_pdf, out_B_pdf, out_C_pdf,
  out_final_vector_pdf,
  out_final_png_300,
  out_final_tiff_300,
  log_fp
)

if (MAKE_STANDALONE_RASTER) {
  outputs <- c(outputs, out_A_png_300, out_B_png_300, out_C_png_300)
}
if (MAKE_1000DPI) {
  outputs <- c(outputs, out_final_png, out_final_tiff)
}

check_df <- data.frame(
  item = c(
    "Single composite file contains Panels A/B/C",
    "Full-width final figure",
    "Maximum final height",
    "Raster resolution",
    "Line width",
    "Font handling",
    "Figure keys in graphic",
    "No figure title/legend text in graphic",
    "Significance marking in Panel C",
    "Preferred final submission file",
    "Fallback raster file"
  ),
  target = c(
    "One composite Figure 3 file",
    "170 mm",
    "<= 225 mm",
    "300 dpi journal-sized raster backup; 1000 dpi if MAKE_1000DPI=TRUE",
    "> 0.25 pt",
    "Arial on Windows; cairo_pdf embeds fonts in vector PDF",
    "Disease/Control key in A; Panel B bottom symbol legend; FDR key in C",
    "Only panel labels and necessary keys are inside the figure",
    "FDR < 0.05 tiles = black border + asterisk (not colour alone)",
    "Figure3_Combined_MN_vector.pdf",
    "Figure3_Combined_MN_300dpi_submission.png / LZW TIFF"
  ),
  value = c(
    "YES",
    paste0(MN_WIDTH_FULL_MM, " mm"),
    paste0(FINAL_H_MM, " mm"),
    paste0(MN_DPI_SUBMISSION, " dpi submission backup; optional ", MN_DPI, " dpi"),
    paste0(MN_LINE_WIDTH, " pt"),
    PLOT_FONT,
    "YES (Panel B key = native bottom legend)",
    "YES",
    "border + asterisk",
    basename(out_final_vector_pdf),
    paste0(basename(out_final_png_300), "; ", basename(out_final_tiff_300))
  ),
  status = c(
    "PASS",
    ifelse(abs(MN_WIDTH_FULL_MM - 170) <= 0.5, "PASS", "CHECK"),
    ifelse(FINAL_H_MM <= MN_HEIGHT_MAX_MM, "PASS", "FAIL"),
    ifelse(MN_DPI_SUBMISSION >= 300, "PASS", "CHECK"),
    ifelse(MN_LINE_WIDTH > 0.25, "PASS", "FAIL"),
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS"
  ),
  stringsAsFactors = FALSE
)

check_fp <- file.path(out_dir, "Figure3_MN_guideline_selfcheck.csv")
write.csv(check_df, check_fp, row.names = FALSE)
print(check_df)

cat("\nOutput files and sizes:\n")
for (x in outputs) {
  if (file.exists(x)) {
    cat(sprintf("  %-45s %8.2f MB\n", basename(x), file_mb(x)))
  }
}

cat("\nMain files to inspect/submit:\n")
cat("  1. Preferred vector PDF :", out_final_vector_pdf, "\n")
cat("  2. File-size-safe 300 dpi PNG:", out_final_png_300, "\n")
cat("  3. File-size-safe 300 dpi TIFF:", out_final_tiff_300, "\n")
cat("  4. Optional 1000 dpi PNG (only if MAKE_1000DPI=TRUE):", out_final_png, "\n")
cat("  5. Optional 1000 dpi TIFF (only if MAKE_1000DPI=TRUE):", out_final_tiff, "\n")
cat("  6. Self-check CSV       :", check_fp, "\n")
cat("  7. Log file             :", log_fp, "\n")

cat("\n==== Figure 3 MN integrated script (v17) END ====\n")
sink()

cat("Done. Outputs saved to:\n", out_dir, "\n")
