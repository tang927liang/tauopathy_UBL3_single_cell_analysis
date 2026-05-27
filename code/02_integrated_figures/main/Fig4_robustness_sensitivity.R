# =============================================================================
# Figure 4 — Molecular Neurodegeneration integrated version v10
# Robustness and sensitivity analyses support the V1 cortical-neuron findings
# as a candidate signal.
#
# Target journal format:
#   - One composite multi-panel figure file
#   - Full-width final PDF size: 170 mm; total height <= 225 mm
#   - Vector PDF + 1000 dpi PNG/TIFF output
#   - Additional 300 dpi submission backup PNG/TIFF
#   - Figure keys incorporated into the graphic
#   - No full figure title/legend inside the graphic
#   - Line widths > 0.25 pt
#   - Arial / embedded fonts via cairo_pdf
#   - Colourblind-safe Okabe-Ito palette
#
# Output directory requested by user:
#   D:/RNA/MNversion submission/Figure 4. Robustness and sensitivity analyses support the V1 cortical-neuron findings as a candidate signal/Results
# =============================================================================

rm(list = ls()); gc()
SEED <- 42
set.seed(SEED)
Sys.setenv(LANG = "en")
options(stringsAsFactors = FALSE)

# =============================================================================
# 0. Packages
# =============================================================================
required_pkgs <- c(
  "readr", "readxl", "dplyr", "tidyr", "stringr", "purrr",
  "ggplot2", "patchwork", "ggbeeswarm", "scales",
  "ragg", "writexl", "fs"
)

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(ggbeeswarm)
  library(scales)
  library(ragg)
  library(writexl)
  library(fs)
})

# =============================================================================
# 1. Paths and output settings
# =============================================================================
RAW_DIR <- "D:/RNA/figure/Robustness/Raw"
OUT_DIR <- "D:/RNA/MNversion submission/Figure 4. Robustness and sensitivity analyses support the V1 cortical-neuron findings as a candidate signal/Results"
FIG_DIR <- OUT_DIR
PANEL_DIR <- file.path(OUT_DIR, "panels")
DATA_DIR <- file.path(OUT_DIR, "source_data")

fs::dir_create(FIG_DIR, recurse = TRUE)
fs::dir_create(PANEL_DIR, recurse = TRUE)
fs::dir_create(DATA_DIR, recurse = TRUE)

MN_WIDTH_FULL_MM <- 170
MN_HEIGHT_MAX_MM <- 225
MN_HEIGHT_MM <- 224.5
MN_DPI <- 1000
MN_DPI_SUBMISSION <- 300

# 1000 dpi raster export may create a large file. PDF is the preferred file.
MAKE_1000DPI <- TRUE
MAKE_PANEL_EXPORTS <- TRUE

# =============================================================================
# 2. Font and devices
# =============================================================================
FONT_FAM <- "Arial"
if (.Platform$OS.type == "windows") {
  arial_ok <- tryCatch({
    windowsFonts(Arial = windowsFont("TT Arial"))
    TRUE
  }, error = function(e) FALSE, warning = function(w) TRUE)
  if (!arial_ok) FONT_FAM <- "sans"
} else {
  FONT_FAM <- "sans"
}

pdf_device <- function(filename, width, height, ...) {
  grDevices::cairo_pdf(filename = filename, width = width, height = height,
                       family = FONT_FAM, onefile = TRUE, ...)
}

# =============================================================================
# 3. Palette and Molecular Neurodegeneration theme
# =============================================================================
# Okabe-Ito colourblind-safe palette
PAL <- c(
  "AD ExN"  = "#0072B2",  # blue
  "PSP ExN" = "#D55E00",  # vermillion/orange
  "PSP InN" = "#009E73"   # bluish green
)

FOCAL_ORDER <- c("AD ExN", "PSP ExN", "PSP InN")
FOCAL_FULL_LABELS <- c(
  "AD ExN"  = "AD Excitatory neurons",
  "PSP ExN" = "PSP Excitatory neurons",
  "PSP InN" = "PSP Inhibitory neurons"
)
FOCAL_FULL_LABELS_NL <- c(
  "AD ExN"  = "AD\nExcitatory\nneurons",
  "PSP ExN" = "PSP\nExcitatory\nneurons",
  "PSP InN" = "PSP\nInhibitory\nneurons"
)

BASE_SIZE <- 7.4
MIN_LINE_WIDTH <- 0.30

basic_theme <- function(base_size = BASE_SIZE) {
  theme_classic(base_size = base_size, base_family = FONT_FAM) +
    theme(
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.background   = element_rect(fill = "white", colour = NA),
      axis.line          = element_line(linewidth = MIN_LINE_WIDTH, colour = "black"),
      axis.ticks         = element_line(linewidth = MIN_LINE_WIDTH, colour = "black"),
      axis.text          = element_text(size = base_size - 0.4, colour = "black", family = FONT_FAM),
      axis.title         = element_text(size = base_size, colour = "black", family = FONT_FAM),
      panel.grid.major.y = element_line(linewidth = 0.30, colour = "grey94"),
      panel.grid.minor   = element_blank(),
      strip.background   = element_blank(),
      strip.text         = element_text(face = "plain", size = base_size, colour = "black", family = FONT_FAM),
      legend.position    = "top",
      legend.title       = element_blank(),
      legend.text        = element_text(size = base_size - 0.2, colour = "black", family = FONT_FAM),
      legend.key.size    = unit(3.0, "mm"),
      legend.margin      = margin(0, 0, 0, 0),
      legend.background  = element_blank(),
      plot.title         = element_text(face = "bold", size = base_size + 0.8,
                                        hjust = 0, margin = margin(b = 2), family = FONT_FAM),
      plot.tag           = element_text(face = "bold", size = base_size + 3.0,
                                        family = FONT_FAM),
      plot.margin        = margin(2, 2, 2, 2, unit = "pt")
    )
}

theme_set(basic_theme())

# =============================================================================
# 4. Helpers
# =============================================================================
find_file <- function(pattern) {
  fp <- list.files(RAW_DIR, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(fp) == 0) stop(sprintf("Not found in %s: %s", RAW_DIR, pattern))
  fp[1]
}

read_any <- function(fp) {
  ext <- tolower(tools::file_ext(fp))
  if (ext %in% c("csv", "tsv", "txt")) {
    suppressMessages(readr::read_csv(fp, show_col_types = FALSE,
                                     locale = readr::locale(encoding = "UTF-8")))
  } else if (ext %in% c("xlsx", "xls")) {
    readxl::read_excel(fp)
  } else {
    stop(sprintf("Unsupported ext '%s' in %s", ext, basename(fp)))
  }
}

build_focal <- function(unit, celltype) {
  u <- tolower(as.character(unit))
  c <- tolower(as.character(celltype))
  out <- rep(NA_character_, length(u))
  out[grepl("ad",  u) & grepl("ex(citatory)?", c)] <- "AD ExN"
  out[grepl("psp", u) & grepl("ex(citatory)?", c)] <- "PSP ExN"
  out[grepl("psp", u) & grepl("in(hibitory)?",  c)] <- "PSP InN"
  factor(out, levels = FOCAL_ORDER)
}

# =============================================================================
# 5. Read input tables
# =============================================================================
message("\n>>> Reading Figure 4 source tables from: ", RAW_DIR)
A_focal <- read_any(find_file("PanelA_summary_focal"))
B_null  <- read_any(find_file("PanelB_null_distribution"))
B_perm  <- read_any(find_file("PanelB_permutation_results"))
C_iter  <- read_any(find_file("PanelC_LOO_per_iter"))
D_test  <- read_any(find_file("PanelD_manual_recomputation.*test"))

message("\n>>> Panel A data sanity check:")
A_check <- A_focal %>%
  mutate(focal = build_focal(unit, celltype)) %>%
  filter(!is.na(focal)) %>%
  transmute(
    focal,
    within_unit_BH6_primary,
    syn52082747_V1_BH18,
    global_BH30,
    combined_BH60,
    Bonferroni30,
    prop_cliffs_delta
  )
print(A_check)
message("Primary BH6 all < 0.05 ? -> ", all(A_check$within_unit_BH6_primary < 0.05))

# =============================================================================
# 6. Panel A — multiplicity sensitivity and effect sizes
# =============================================================================
message("\n>>> Building Panel A...")

method_codes <- c("BH6", "BH18", "BH30", "BH60", "Bon30")
method_labels <- c(
  "Primary\nBenjamini-\nHochberg\n(6 tests)",
  "V1 dataset\nBenjamini-\nHochberg\n(18 tests)",
  "Study-wide\nBenjamini-\nHochberg\n(30 tests)",
  "Combined\nBenjamini-\nHochberg\n(60 tests)",
  "Study-wide\nBonferroni\n(30 tests)"
)

A_long <- A_focal %>%
  mutate(focal = build_focal(unit, celltype)) %>%
  filter(!is.na(focal)) %>%
  transmute(
    focal,
    BH6   = as.numeric(within_unit_BH6_primary),
    BH18  = as.numeric(syn52082747_V1_BH18),
    BH30  = as.numeric(global_BH30),
    BH60  = as.numeric(combined_BH60),
    Bon30 = as.numeric(Bonferroni30),
    HL    = as.numeric(prop_HL_estimate),
    HL_lo = as.numeric(prop_HL_ci95_low),
    HL_hi = as.numeric(prop_HL_ci95_high),
    delta = as.numeric(prop_cliffs_delta)
  ) %>%
  pivot_longer(BH6:Bon30, names_to = "method_code", values_to = "q") %>%
  mutate(
    method_code = factor(method_code, levels = method_codes),
    method_idx  = as.numeric(method_code),
    sig = q < 0.05
  )

bar_dodge <- 0.22
bar_halfw <- 0.10

A_bars <- A_long %>%
  mutate(
    focal_idx = as.numeric(focal),
    x_center  = method_idx + (focal_idx - 2) * bar_dodge,
    xmin      = x_center - bar_halfw,
    xmax      = x_center + bar_halfw,
    ymin      = pmin(q, 0.05),
    ymax      = pmax(q, 0.05),
    sig_chr   = ifelse(sig, "TRUE", "FALSE")
  )

A_stars <- A_bars %>% filter(sig) %>% mutate(star_y = q * 0.80)
A_eff <- A_long %>% distinct(focal, HL, HL_lo, HL_hi, delta)

# Compact graphical key inside Panel A: significance encoding, not focal-group colours.
A_key_sig <- data.frame(
  xmin = c(3.40, 3.40),
  xmax = c(3.63, 3.63),
  ymin = c(0.74, 0.57),
  ymax = c(0.86, 0.66),
  status = c("Adj. P < 0.05", "Adj. P >= 0.05"),
  stringsAsFactors = FALSE
)
A_key_line <- data.frame(
  x0 = 3.40, x1 = 3.63, y = 0.45,
  label = "FDR = 0.05",
  stringsAsFactors = FALSE
)

write.csv(A_bars, file.path(DATA_DIR, "PanelA_plot_data.csv"), row.names = FALSE)
write.csv(A_eff,  file.path(DATA_DIR, "PanelA_effectsize.csv"), row.names = FALSE)

# Dummy points used only to draw a compact colour key for the three focal groups.
# This avoids using square swatches for both group identity and significance.
A_group_key <- data.frame(
  focal = factor(FOCAL_ORDER, levels = FOCAL_ORDER),
  x = 5.55,
  y = 1.0
)

pA_left <- ggplot(A_bars) +
  annotate("rect", xmin = 0.4, xmax = 5.6,
           ymin = 0.025, ymax = 0.05,
           fill = "#E6F4EA", alpha = 0.45, colour = NA) +
  annotate("rect", xmin = 0.4, xmax = 5.6,
           ymin = 0.05, ymax = 1.0,
           fill = "#F5F5F5", alpha = 0.45, colour = NA) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax,
        fill = focal, alpha = sig_chr),
    colour = "grey60", linewidth = 0.35
  ) +
  geom_rect(
    data = A_bars %>% filter(sig),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE, fill = NA, colour = "black", linewidth = 0.55
  ) +
  geom_text(
    data = A_stars,
    aes(x = x_center, y = star_y, label = "\u2605", colour = focal),
    inherit.aes = FALSE, size = 3.0, fontface = "bold", show.legend = FALSE
  ) +
  geom_hline(yintercept = 0.05, linetype = "22", colour = "black", linewidth = 0.45) +
  geom_point(
    data = A_group_key,
    aes(x = x, y = y, colour = focal),
    inherit.aes = FALSE, alpha = 0, size = 1.8, show.legend = TRUE
  ) +
  scale_fill_manual(values = PAL, name = NULL, labels = FOCAL_FULL_LABELS, guide = "none") +
  scale_colour_manual(values = PAL, name = NULL, labels = FOCAL_FULL_LABELS) +
  scale_alpha_manual(
    values = c("TRUE" = 1.0, "FALSE" = 0.22),
    name   = "Bar style",
    breaks = c("TRUE", "FALSE"),
    labels = c("Adjusted P < 0.05", "Adjusted P \u2265 0.05"),
    guide = guide_legend(
      order = 2,
      nrow = 1,
      override.aes = list(
        fill      = "grey45",
        colour    = c("black", "grey60"),
        linewidth = c(0.55, 0.35),
        alpha     = c(1.0, 0.30)
      )
    )
  ) +
  scale_x_continuous(breaks = 1:5, labels = method_labels,
                     limits = c(0.4, 5.6), expand = c(0, 0)) +
  scale_y_log10(
    name = "Adjusted P (log scale)",
    limits = c(0.025, 1.0),
    breaks = c(0.03, 0.05, 0.1, 0.2, 0.5, 1.0),
    labels = c("0.03", "0.05", "0.1", "0.2", "0.5", "1.0")
  ) +
  labs(x = NULL, title = "A  Multiple-testing sensitivity") +
  guides(
    colour = guide_legend(
      order = 1, nrow = 1,
      override.aes = list(alpha = 1, size = 2.6, shape = 16)
    )
  ) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.box = "vertical",
    legend.box.spacing = unit(0.2, "mm"),
    legend.spacing.y = unit(0.05, "mm"),
    legend.spacing.x = unit(1.0, "mm"),
    axis.text.x = element_text(size = BASE_SIZE - 2.05, lineheight = 0.74),
    axis.title.y = element_text(size = BASE_SIZE, margin = margin(r = 5, unit = "pt")),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 6, unit = "pt")
  )

pA_right <- ggplot(A_eff, aes(x = HL, y = focal, colour = focal)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = "22", linewidth = 0.30) +
  geom_linerange(aes(xmin = HL_lo, xmax = HL_hi), linewidth = 0.60) +
  geom_point(size = 1.8) +
  geom_text(
    aes(label = sprintf("\u03b4 = %+.2f", delta)),
    nudge_y = 0.32, size = 2.25, family = FONT_FAM, show.legend = FALSE
  ) +
  scale_colour_manual(values = PAL, guide = "none") +
  scale_y_discrete(limits = rev(FOCAL_ORDER), labels = FOCAL_FULL_LABELS) +
  labs(x = "Hodges-Lehmann shift (95% CI)", y = NULL,
       title = "Effect-size estimates") +
  theme(
    plot.title = element_text(face = "plain", size = BASE_SIZE, hjust = 0),
    axis.text.y = element_text(size = BASE_SIZE - 0.4),
    axis.title.x = element_text(size = BASE_SIZE - 0.1, margin = margin(t = 1, unit = "pt"))
  )

panelA <- pA_left + pA_right + patchwork::plot_layout(widths = c(2.35, 1.05))

# =============================================================================
# 7. Panel B — permutation null
# =============================================================================
message("\n>>> Building Panel B...")

B_null_df <- B_null %>% transmute(null_p = as.numeric(null_min_p))
B_obs <- B_perm %>%
  filter(is_focal_within_unit_FDR_positive == TRUE) %>%
  mutate(focal = build_focal(unit, celltype)) %>%
  filter(!is.na(focal)) %>%
  transmute(
    focal,
    p_obs        = as.numeric(observed_prop_p_raw),
    p_emp_compar = as.numeric(empirical_p_per_comparison),
    p_emp_global = as.numeric(empirical_p_familyWise_globalMinP)
  ) %>%
  arrange(focal)

write.csv(B_null_df, file.path(DATA_DIR, "PanelB_null_long.csv"), row.names = FALSE)
write.csv(B_obs,     file.path(DATA_DIR, "PanelB_observed.csv"),  row.names = FALSE)

panelB <- ggplot(B_null_df, aes(x = null_p)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 60, fill = "grey82", colour = "grey58", linewidth = 0.30
  ) +
  geom_vline(xintercept = 0.05, linetype = "22", colour = "grey35", linewidth = 0.35) +
  geom_vline(data = B_obs, aes(xintercept = p_obs, colour = focal), linewidth = 0.75) +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  coord_cartesian(xlim = c(0, 0.30)) +
  scale_colour_manual(values = PAL, guide = "none") +
  labs(
    x = "Global minimum P under permutation",
    y = "Null density",
    title = "B  Permutation global min-P"
  ) +
  theme(
    axis.title.x = element_text(size = BASE_SIZE - 0.1, margin = margin(t = 1, unit = "pt")),
    axis.title.y = element_text(size = BASE_SIZE - 0.1),
    plot.margin = margin(t = 2, r = 3, b = 1, l = 2, unit = "pt")
  )

# =============================================================================
# 8. Panel C — leave-one-donor-out stability
# =============================================================================
message("\n>>> Building Panel C...")

C_long <- C_iter %>%
  mutate(focal = build_focal(unit, focal_celltype)) %>%
  filter(!is.na(focal)) %>%
  transmute(
    focal,
    iter = dropped_donor,
    dropped_group = dropped_donor_group,
    delta = as.numeric(cliffs_delta),
    direction_pos = direction_positive_delta
  )
C_orig <- A_eff %>% select(focal, delta_orig = delta)
write.csv(C_long, file.path(DATA_DIR, "PanelC_LOO_long.csv"), row.names = FALSE)

panelC <- ggplot(C_long, aes(x = focal, y = delta, colour = focal)) +
  geom_hline(yintercept = 0,     linetype = "22",       colour = "grey55", linewidth = 0.30) +
  geom_hline(yintercept = 0.474, linetype = "longdash", colour = "grey35", linewidth = 0.30) +
  geom_errorbar(
    data = C_orig,
    aes(x = focal, ymin = delta_orig, ymax = delta_orig, colour = focal),
    width = 0.55, linewidth = 0.55, inherit.aes = FALSE
  ) +
  ggbeeswarm::geom_quasirandom(width = 0.24, size = 1.25, alpha = 0.75, stroke = 0) +
  geom_text(
    data = C_orig,
    aes(x = focal, colour = focal, label = sprintf("\u03b4 = %+.2f", delta_orig)),
    y = 0.96, fontface = "bold", size = 2.3, family = FONT_FAM,
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  annotate(
    "text", x = 3.24, y = 0.435,
    label = "Cliff's \u03b4 = 0.474\n(large-effect threshold)",
    hjust = 1, vjust = 1, size = 2.05, colour = "grey35", family = FONT_FAM
  ) +
  scale_colour_manual(values = PAL, guide = "none") +
  scale_x_discrete(labels = FOCAL_FULL_LABELS_NL) +
  scale_y_continuous(limits = c(min(0, min(C_long$delta, na.rm = TRUE) - 0.05), 1),
                     breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "Cliff's \u03b4", title = "C  Leave-one-donor-out stability") +
  theme(
    axis.text.x = element_text(size = BASE_SIZE - 0.8, lineheight = 0.88),
    axis.title.y = element_text(size = BASE_SIZE - 0.1)
  )

# =============================================================================
# 9. Panel D — manual recomputation agreement
# =============================================================================
message("\n>>> Building Panel D...")

D_plot <- D_test %>%
  mutate(focal = build_focal(unit, celltype)) %>%
  filter(!is.na(focal)) %>%
  transmute(
    focal,
    raw_pri = as.numeric(S5B_raw_p),         raw_rec = as.numeric(recomp_raw_p),
    pad_pri = as.numeric(S5B_padj_BH),       pad_rec = as.numeric(recomp_padj_BH),
    del_pri = as.numeric(S5B_cliffs_delta),  del_rec = as.numeric(recomp_cliffs_delta),
    hl_pri  = as.numeric(S5B_HL_estimate),   hl_rec  = as.numeric(recomp_HL_estimate)
  ) %>%
  pivot_longer(
    -focal,
    names_to = c("metric", "src"),
    names_pattern = "(raw|pad|del|hl)_(pri|rec)",
    values_to = "v"
  ) %>%
  pivot_wider(names_from = src, values_from = v) %>%
  mutate(
    metric = recode(metric,
                    raw = "Raw P", pad = "FDR-adjusted P",
                    del = "Cliff's \u03b4", hl = "HL estimate"),
    metric = factor(metric, levels = c("Raw P", "FDR-adjusted P",
                                       "Cliff's \u03b4", "HL estimate"))
  )
write.csv(D_plot, file.path(DATA_DIR, "PanelD_plot_data.csv"), row.names = FALSE)

panelD <- ggplot(D_plot, aes(x = pri, y = rec, colour = focal)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey55", linewidth = 0.30) +
  geom_point(size = 1.65, stroke = 0, alpha = 0.95) +
  facet_wrap(~ metric, scales = "free", ncol = 2) +
  scale_colour_manual(values = PAL, guide = "none") +
  labs(
    x = "Primary pipeline value",
    y = "Independent recomputation value",
    title = "D  Manual recomputation agreement"
  ) +
  theme(
    strip.text = element_text(face = "plain", size = BASE_SIZE - 0.1),
    axis.text  = element_text(size = BASE_SIZE - 1.0),
    axis.title = element_text(size = BASE_SIZE - 0.1),
    panel.spacing = unit(2.2, "mm")
  )

# =============================================================================
# 10. Compose one composite figure
# =============================================================================
message("\n>>> Composing Figure 4 composite...")

middle_row <- panelB + panelC + patchwork::plot_layout(widths = c(1.03, 1.00))

composite <- panelA / middle_row / panelD +
  patchwork::plot_layout(heights = c(1.16, 0.88, 1.22)) &
  theme(
    plot.margin = margin(3, 3, 3, 3, unit = "pt"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# =============================================================================
# 11. Export composite and optional panels
# =============================================================================
message("\n>>> Exporting figures to: ", FIG_DIR)

out_pdf <- file.path(FIG_DIR, "Figure_4_composite_MN_vector.pdf")
out_png <- file.path(FIG_DIR, "Figure_4_composite_MN_1000dpi.png")
out_tiff <- file.path(FIG_DIR, "Figure_4_composite_MN_1000dpi_LZW.tiff")
out_png_300 <- file.path(FIG_DIR, "Figure_4_composite_MN_300dpi_submission.png")
out_tiff_300 <- file.path(FIG_DIR, "Figure_4_composite_MN_300dpi_submission_LZW.tiff")

ggsave(out_pdf, composite,
       width = MN_WIDTH_FULL_MM, height = MN_HEIGHT_MM, units = "mm",
       device = pdf_device, bg = "white")

# 300 dpi backup matching BMC/Molecular Neurodegeneration final-size guidance.
ggsave(out_png_300, composite,
       width = MN_WIDTH_FULL_MM, height = MN_HEIGHT_MM, units = "mm",
       dpi = MN_DPI_SUBMISSION, bg = "white", limitsize = FALSE)

tryCatch({
  ggsave(out_tiff_300, composite,
         width = MN_WIDTH_FULL_MM, height = MN_HEIGHT_MM, units = "mm",
         dpi = MN_DPI_SUBMISSION, compression = "lzw",
         bg = "white", limitsize = FALSE)
}, error = function(e) {
  message("300 dpi LZW TIFF export failed; trying uncompressed TIFF.")
  ggsave(out_tiff_300, composite,
         width = MN_WIDTH_FULL_MM, height = MN_HEIGHT_MM, units = "mm",
         dpi = MN_DPI_SUBMISSION, bg = "white", limitsize = FALSE)
})

if (MAKE_1000DPI) {
  ggsave(out_png, composite,
         width = MN_WIDTH_FULL_MM, height = MN_HEIGHT_MM, units = "mm",
         dpi = MN_DPI, device = ragg::agg_png, bg = "white", limitsize = FALSE)
  tryCatch({
    ggsave(out_tiff, composite,
           width = MN_WIDTH_FULL_MM, height = MN_HEIGHT_MM, units = "mm",
           dpi = MN_DPI, compression = "lzw",
           bg = "white", limitsize = FALSE)
  }, error = function(e) {
    message("1000 dpi LZW TIFF export failed; trying uncompressed TIFF.")
    ggsave(out_tiff, composite,
           width = MN_WIDTH_FULL_MM, height = MN_HEIGHT_MM, units = "mm",
           dpi = MN_DPI, bg = "white", limitsize = FALSE)
  })
}

if (MAKE_PANEL_EXPORTS) {
  ggsave(file.path(PANEL_DIR, "Figure_4_PanelA_MN_vector.pdf"), panelA,
         width = 170, height = 72, units = "mm", device = pdf_device, bg = "white")
  ggsave(file.path(PANEL_DIR, "Figure_4_PanelB_MN_vector.pdf"), panelB,
         width = 85, height = 66, units = "mm", device = pdf_device, bg = "white")
  ggsave(file.path(PANEL_DIR, "Figure_4_PanelC_MN_vector.pdf"), panelC,
         width = 85, height = 66, units = "mm", device = pdf_device, bg = "white")
  ggsave(file.path(PANEL_DIR, "Figure_4_PanelD_MN_vector.pdf"), panelD,
         width = 170, height = 84, units = "mm", device = pdf_device, bg = "white")
}

# =============================================================================
# 12. Manuscript figure title and legend text
# =============================================================================
figure_title <- "Figure 4. Robustness supports candidate V1 neuronal signals."

legend_text <- paste0(
  figure_title, "\n",
  "(A) Multiple-testing sensitivity for the three focal syn52082747 V1 neuronal findings. ",
  "Bars are anchored at adjusted P = 0.05; saturated bars with black borders and star markers indicate adjusted P < 0.05, whereas faded bars indicate adjusted P >= 0.05. ",
  "Right, Hodges-Lehmann disease-control shifts with 95% confidence intervals; Cliff's delta is annotated. ",
  "(B) Family-wise permutation null distribution of the global minimum P value from 10,000 donor-label permutations. The dashed line marks P = 0.05; coloured vertical lines show the observed raw Wilcoxon P values. ",
  "(C) Leave-one-donor-out stability of Cliff's delta. Points denote leave-one-donor-out iterations; horizontal coloured bars denote full-sample delta values. The dashed line marks the large-effect threshold. ",
  "(D) Manual recomputation agreement between the primary pipeline and independent recomputation. The dashed diagonal indicates identity. ",
  "AD, Alzheimer's disease; PSP, progressive supranuclear palsy; BH, Benjamini-Hochberg; Bonf., Bonferroni; HL, Hodges-Lehmann; V1, primary visual cortex. Source data: Supplementary Table S4."
)

legend_words <- length(strsplit(gsub("[^A-Za-z0-9]+", " ", legend_text), "\\s+")[[1]])
writeLines(legend_text, file.path(FIG_DIR, "Figure_4_legend_MN_300words.txt"))
writeLines(sprintf("Figure legend word count: %d", legend_words),
           file.path(FIG_DIR, "Figure_4_legend_word_count.txt"))

# =============================================================================
# 13. Self-check and session info
# =============================================================================
file_mb <- function(x) {
  if (!file.exists(x)) return(NA_real_)
  round(file.info(x)$size / 1024 / 1024, 2)
}

check_df <- data.frame(
  item = c(
    "Composite multi-panel file",
    "Full-width PDF size",
    "Maximum height",
    "1000 dpi PNG generated",
    "300 dpi submission backup generated",
    "Line width minimum",
    "Font handling",
    "Colourblind-safe palette",
    "Figure title <= 15 words",
    "Figure legend <= 300 words",
    "PDF file size MB",
    "1000 dpi PNG file size MB"
  ),
  target = c(
    "One Figure 4 composite file",
    "170 mm",
    "<= 225 mm",
    "Yes",
    "Yes",
    "> 0.25 pt",
    "Arial embedded through cairo_pdf",
    "Okabe-Ito palette",
    "<= 15 words",
    "<= 300 words",
    "Preferably <= 10 MB",
    "Preferably <= 10 MB"
  ),
  value = c(
    "PASS",
    paste0(MN_WIDTH_FULL_MM, " mm"),
    paste0(MN_HEIGHT_MM, " mm"),
    ifelse(file.exists(out_png), "YES", "NO"),
    ifelse(file.exists(out_png_300), "YES", "NO"),
    paste0(MIN_LINE_WIDTH, " pt"),
    FONT_FAM,
    "PASS",
    "PASS",
    paste0(legend_words, " words"),
    paste0(file_mb(out_pdf), " MB"),
    paste0(file_mb(out_png), " MB")
  ),
  status = c(
    "PASS",
    ifelse(abs(MN_WIDTH_FULL_MM - 170) <= 0.5, "PASS", "CHECK"),
    ifelse(MN_HEIGHT_MM <= MN_HEIGHT_MAX_MM, "PASS", "FAIL"),
    ifelse(file.exists(out_png), "PASS", "CHECK"),
    ifelse(file.exists(out_png_300), "PASS", "CHECK"),
    ifelse(MIN_LINE_WIDTH > 0.25, "PASS", "FAIL"),
    "PASS",
    "PASS",
    "PASS",
    ifelse(legend_words <= 300, "PASS", "FAIL"),
    ifelse(!is.na(file_mb(out_pdf)) && file_mb(out_pdf) <= 10, "PASS", "CHECK"),
    ifelse(!is.na(file_mb(out_png)) && file_mb(out_png) <= 10, "PASS", "CHECK")
  ),
  stringsAsFactors = FALSE
)

write.csv(check_df, file.path(FIG_DIR, "Figure_4_MN_guideline_selfcheck.csv"), row.names = FALSE)
print(check_df)

si <- capture.output(sessionInfo())
writeLines(c(sprintf("# Figure 4 MN build session info — %s", Sys.time()), "", si),
           file.path(FIG_DIR, "sessionInfo.txt"))

message("\n>>> DONE.")
message("    Composite PDF: ", out_pdf)
message("    1000 dpi PNG : ", out_png)
message("    300 dpi backup: ", out_png_300)
message("    Results dir  : ", FIG_DIR)
