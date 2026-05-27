###############################################################################
## Figure 1 - Study workflow and dataset selection
## UBL3 in tauopathy | Molecular Neurodegeneration submission
##
## PRISMA-style flowchart (6 phases). Self-contained: NO external data input.
## Output: 170 x 215 mm; vector PDF (cairo_pdf, fonts embedded) + 1000 dpi PNG.
##
## HOW TO RUN:
##   1) Edit the REPO line below so it points to your repository folder.
##   2) Source this whole file (RStudio: Ctrl+Shift+S, or the "Source" button).
##   Output is written automatically to  <REPO>/output/figures/Fig1/
###############################################################################
rm(list = ls()); gc()

## === 改这一行就够了:指向你的代码库根目录 =====================================
REPO <- "D:/RNA/Code/UBL3_tauopathy"
## ============================================================================

SEED <- 20251023; set.seed(SEED)
Sys.setenv(LANG = "en"); options(stringsAsFactors = FALSE)
for (pkg in c("ggplot2", "ragg")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
suppressPackageStartupMessages({
  library(ggplot2); library(ragg); library(grid)
})
if (.Platform$OS.type == "windows") {
  tryCatch(windowsFonts(Arial = windowsFont("Arial")),
           error = function(e) message("Arial \u5B57\u4F53\u6CE8\u518C\u8DF3\u8FC7"))
}
###############################################################################
# 0. \u8F93\u51FA\u8DEF\u5F84
###############################################################################
out_dir <- file.path(REPO, "output", "figures", "Fig1")   # 由 REPO 拼出绝对路径
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
###############################################################################
# 1. Canvas + \u8272\u76F2\u53CB\u597D\u914D\u8272 (Okabe-Ito \u542F\u53D1)
###############################################################################
fig_w_mm  <- 170
fig_h_mm  <- 215
base_fontfamily <- "Arial"
png_dpi   <- 1000
# 6 \u4E2A phase \u989C\u8272 \u2014 \u5DF2\u9A8C\u8BC1\u5BF9 protan/deutan/tritan \u53EF\u533A\u5206
col_id    <- "#034748"   # dark teal     IDENTIFICATION
col_scr   <- "#117A8A"   # medium teal   SCREENING
col_elig  <- "#C7901C"   # dark amber    ELIGIBILITY
col_comp  <- "#B0533A"   # vermillion    COMPARABILITY  (replaces pure red)
col_incl  <- "#117A4F"   # bluish green  INCLUDED       (replaces pure green)
col_anal  <- "#1F3A6B"   # navy          ANALYSIS
# \u5185\u5BB9\u533A\u989C\u8272
col_excl_fill   <- "#FCEAE3"; col_excl_border <- "#B0533A"
col_para_fill   <- "#D6E5EE"; col_para_border <- "#117A8A"
col_diam_fill   <- "#FCEFCC"; col_diam_border <- "#C7901C"
col_box_fill    <- "#FFFFFF"; col_box_border  <- "#2A2A2A"
col_anal_fill   <- "#E5E9F2"
###############################################################################
# 2. Phase \u8272\u6761
###############################################################################
phase_x_left  <- 3
phase_x_right <- 20
phases <- data.frame(
  label = c("IDENTIFICATION","SCREENING","ELIGIBILITY",
            "COMPARABILITY","INCLUDED","ANALYSIS"),
  # v5: COMPARABILITY 24mm 容纳 13 字竖排 (v4 仅 19mm → "Y" 被剪)
  y_top = c(213, 158, 120,  87, 62, 45),
  y_bot = c(159, 121,  88,  63, 46,  1),
  fill  = c(col_id, col_scr, col_elig, col_comp, col_incl, col_anal),
  stringsAsFactors = FALSE
)
phases$y_mid <- (phases$y_top + phases$y_bot) / 2
# 略小字号，避免 "COMPARABILITY" / "INCLUDED" 竖排超出短 band
phase_font_size <- 2.55   # ≈ 7.25 pt
###############################################################################
# 3. \u4E3B\u6D41\u7A0B\u6846\uFF08\u4E2D\u5217 cx = 73\uFF09
###############################################################################
main_cx <- 73
main_boxes <- data.frame(
  id     = 1:9,
  cx     = rep(main_cx, 9),
  cy     = c(204, 184, 166, 148, 130, 112,  94,  74,  55),
  w      = c( 92,  92,  74,  92,  30,  92,  30,  86,  92),
  h      = c( 15,  15,  10,  15,  12,  15,  12,  15,  15),
  shape  = c("rect","para","rect","rect","diamond",
             "rect","diamond","rect","rect"),
  border = c(col_id, col_para_border, col_box_border, col_box_border,
             col_diam_border, col_box_border, col_diam_border,
             col_box_border, col_incl),
  fill   = c(col_box_fill, col_para_fill, col_box_fill, col_box_fill,
             col_diam_fill, col_box_fill, col_diam_fill,
             col_box_fill, col_box_fill),
  title  = c(
    "Pre-specified scope & inclusion criteria",
    "Records identified from public repositories",
    "Records de-duplicated across repositories",
    "Records screened on title and metadata",
    "Screening\npassed?",
    "Full eligibility assessment of candidate datasets",
    "Eligible for\nanalysis?",
    "Cross-cohort comparability assessment",
    "Included in donor-level UBL3 re-analysis"
  ),
  detail = c(
    "Human cortical tauopathies \u2014 AD, PSP, PiD (FTD-tau), and CTE;\nsc/snRNA-seq; matched non-diseased controls; \u22653 donors per group.",
    "GEO \u00B7 ArrayExpress (BioStudies) \u00B7 AMP-AD Synapse \u00B7 CZ CELLxGENE\nSearch completed 23 June 2025; refreshed 4 August 2025",
    "",
    "Retain: human; brain/CNS; sc/snRNA-seq;\nprimary tau-related neuropathology",
    "",
    "Counts + sample-level metadata accessible; matched non-diseased\ncontrols; \u22653 biological replicates per group. Six cohorts retained.",
    "",
    "Cohorts differing in dataset source, anatomical region or platform\nvs the V1 anchor (syn52082747) were not retained for primary inference.",
    "4 cohorts \u2192 7 donor \u00D7 disease\u2013region analytical units\nAD, PSP, PiD (FTD-tau);   MFG, PFC, EC, SFG, V1"
  ),
  stringsAsFactors = FALSE
)
###############################################################################
# 4. \u53F3\u4FA7 Exclusion \u6846\uFF08\u4E0E\u5BF9\u5E94 diamond / box \u5728\u540C\u4E00 y \u4E0A\u5BF9\u9F50\uFF09
###############################################################################
excl_cx <- 148
excl_boxes <- data.frame(
  id    = 1:3,
  cx    = c(excl_cx, excl_cx, excl_cx),
  cy    = c(130, 94, 74),         # 与 Diamond 5 / Diamond 7 / Box 8 同 y
  w     = c( 38, 38, 38),
  h     = c( 18, 20, 18),         # v5: Excl 2 16→20, Excl 3 15→18 防止重叠
  title = c(
    "Excluded at screening",
    "Excluded after\neligibility",
    "Not retained for\nprimary inference"
  ),
  detail = c(
    "\u2022 Non-human\n\u2022 Non-brain / non-CNS\n\u2022 Not sc/snRNA-seq\n\u2022 Non-tauopathy pathology",
    "\u2022 No matched controls\n\u2022 <3 biological reps\n\u2022 Missing metadata\n(none excluded here)",
    "Dataset/region/platform:\n\u2022 GSE155114 (CTE)\n\u2022 GSE261807 (CTE)"
  ),
  title_lines  = c(1, 2, 2),
  detail_lines = c(4, 4, 3),
  stringsAsFactors = FALSE
)
###############################################################################
# 5. ANALYSIS \u5BB9\u5668 + 2 \u5B50\u6846 (Step1 / Step2)
#    \u5BB9\u5668 cx=85, w=120 -> \u5DE6\u7F18 x=25 (\u8DDD phase band right=20 \u67095mm \u95F4\u9699)
#    Step 1 cx=56, w=50  -> \u5DE6\u7F18 x=31 (\u5BB9\u5668\u5185\u90E8\u67096mm padding)
#    Step 2 cx=114, w=50 -> \u53F3\u7F18 x=139 (\u5BB9\u5668\u5185\u90E8\u67096mm padding)
###############################################################################
anal_box <- data.frame(cx = 85, cy = 22, w = 120, h = 42)
anal_subs <- data.frame(
  id     = 1:2,
  cx     = c( 56, 114),
  cy     = c( 18,  18),
  w      = c( 50,  50),
  h      = c( 24,  24),
  title  = c("Step 1\nDetection breadth",
             "Step 2\nConditional UBL3 expression"),
  detail = c("fraction of UBL3\u207A cells per\ndonor \u00D7 disease\u2013region \u00D7 cell-type",
             "among UBL3\u207A cells"),
  stringsAsFactors = FALSE
)
###############################################################################
# 6. \u8F85\u52A9\uFF1A\u5E73\u884C\u56DB\u8FB9\u5F62 / \u83F1\u5F62 polygon
###############################################################################
build_para_poly <- function(rows) {
  do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
    cx <- rows$cx[i]; cy <- rows$cy[i]
    w  <- rows$w[i];  h  <- rows$h[i]
    skew <- 4
    data.frame(
      x = c(cx - w/2 + skew, cx + w/2 + skew, cx + w/2 - skew, cx - w/2 - skew),
      y = c(cy + h/2, cy + h/2, cy - h/2, cy - h/2),
      group_id   = i,
      fill_col   = rows$fill[i],
      border_col = rows$border[i]
    )
  }))
}
build_diam_poly <- function(rows) {
  do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
    cx <- rows$cx[i]; cy <- rows$cy[i]
    w  <- rows$w[i];  h  <- rows$h[i]
    data.frame(
      x = c(cx, cx + w/2, cx, cx - w/2),
      y = c(cy + h/2, cy, cy - h/2, cy),
      group_id   = i,
      fill_col   = rows$fill[i],
      border_col = rows$border[i]
    )
  }))
}
rect_boxes <- main_boxes[main_boxes$shape == "rect", ]
para_boxes <- main_boxes[main_boxes$shape == "para", ]
diam_boxes <- main_boxes[main_boxes$shape == "diamond", ]
para_poly  <- build_para_poly(para_boxes)
diam_poly  <- build_diam_poly(diam_boxes)
###############################################################################
# 7. \u7BAD\u5934\uFF1A\u4E25\u683C\u5BF9\u9F50
###############################################################################
# 7a) \u4E3B\u5782\u76F4\u7BAD\u5934\uFF08\u4E0A\u6846\u5E95 \u2192 \u4E0B\u6846\u9876\uFF09
arrows_v <- data.frame(
  from_id = 1:8,
  to_id   = 2:9
)
arrows_v$x_start <- main_boxes$cx[arrows_v$from_id]
arrows_v$y_start <- main_boxes$cy[arrows_v$from_id] - main_boxes$h[arrows_v$from_id]/2
arrows_v$x_end   <- main_boxes$cx[arrows_v$to_id]
arrows_v$y_end   <- main_boxes$cy[arrows_v$to_id]   + main_boxes$h[arrows_v$to_id]/2 + 0.4
# 7b) Box 9 \u2192 ANALYSIS \u5BB9\u5668\u9876\u90E8\uFF08\u4E2D\u5FC3\u5BF9\u9F50\u5230\u5BB9\u5668\u8FB9\uFF09
arrow_to_anal <- data.frame(
  x_start = main_cx, y_start = main_boxes$cy[9] - main_boxes$h[9]/2,
  x_end   = main_cx, y_end   = anal_box$cy + anal_box$h/2 + 0.4
)
# 7c) \u6C34\u5E73\u7BAD\u5934\uFF08\u4E25\u683C y \u5BF9\u9F50\uFF09\uFF1A
#     \u83F1\u5F62 5 / \u83F1\u5F62 7 / Box 8 -> Excl 1 / Excl 2 / Excl 3
arrows_h <- data.frame(
  x_start = c(main_cx + 30/2,            # Diamond 5 \u53F3\u8FB9\u7F18 (w=30)
              main_cx + 30/2,            # Diamond 7
              main_cx + 86/2),           # Box 8 (w=86)
  y_start = c(130, 94, 74),              # \u4E0E excl_boxes$cy \u5B8C\u5168\u4E00\u81F4
  x_end   = c(excl_cx - 38/2 - 0.3,      # Excl box \u5DE6\u8FB9\u7F18\u63A5\u8FD1
              excl_cx - 38/2 - 0.3,
              excl_cx - 38/2 - 0.3),
  y_end   = c(130, 94, 74)               # \u4E25\u683C\u6C34\u5E73 (y_start == y_end)
)
# 7d) Step1 \u2192 Step2
arrow_step <- data.frame(
  x_start = anal_subs$cx[1] + anal_subs$w[1]/2,
  y_start = anal_subs$cy[1],
  x_end   = anal_subs$cx[2] - anal_subs$w[2]/2 - 0.3,
  y_end   = anal_subs$cy[2]
)
###############################################################################
# 8. Yes / No \u6807\u7B7E\u4F4D\u7F6E
###############################################################################
arrow_labels <- data.frame(
  x     = c(main_cx + 2.5,            # \u83F1\u5F62 5 "Yes" \u5728\u5782\u7EBF\u53F3\u8FB9
            main_cx + 30/2 + 8,        # \u83F1\u5F62 5 "No" \u5728\u6C34\u5E73\u7BAD\u4E0A\u65B9
            main_cx + 2.5,            # \u83F1\u5F62 7 "Yes"
            main_cx + 30/2 + 8),       # \u83F1\u5F62 7 "No"
  y     = c(122, 132,
             86,  96),
  label = c("Yes", "No",
            "Yes", "No")
)
###############################################################################
# 9. \u7ED8\u56FE
###############################################################################
# v5: 从框边推算 title / detail y 位置，确保几何健壮
# title 从框顶向下: 1.8mm padding + title_block_height/2
# detail 从框底向上: 1.8mm padding + detail_block_height/2
# 每行约 2.5mm (size 2.55 + lineheight ≈ 1.0)
line_h_mm <- 2.5
top_pad   <- 1.8
bot_pad   <- 1.8
excl_boxes$title_y  <- excl_boxes$cy + excl_boxes$h/2 - top_pad -
                       (excl_boxes$title_lines * line_h_mm)/2
excl_boxes$detail_y <- excl_boxes$cy - excl_boxes$h/2 + bot_pad +
                       (excl_boxes$detail_lines * line_h_mm)/2
# \u4E3B\u6D41\u6846 title / detail y\uFF1A\u6709\u6210\u5BF9 vs \u53EA\u6709\u6807\u9898
main_boxes$has_detail <- nchar(main_boxes$detail) > 0
main_boxes$title_y  <- ifelse(main_boxes$has_detail,
                              main_boxes$cy + main_boxes$h/4 + 0.7,
                              main_boxes$cy)
main_boxes$detail_y <- main_boxes$cy - main_boxes$h/4 - 0.3
p <- ggplot() +
  # ---- Phase \u8272\u6761 ----
  geom_rect(data = phases,
            aes(xmin = phase_x_left, xmax = phase_x_right,
                ymin = y_bot, ymax = y_top, fill = fill),
            color = NA) +
  geom_text(data = phases,
            aes(x = (phase_x_left + phase_x_right)/2, y = y_mid, label = label),
            angle = 90, color = "white",
            size = phase_font_size, fontface = "bold",
            family = base_fontfamily) +
  # ---- ANALYSIS \u5BB9\u5668\uFF08\u80CC\u666F\u5E95\u5C42\uFF09 ----
  geom_rect(data = anal_box,
            aes(xmin = cx - w/2, xmax = cx + w/2,
                ymin = cy - h/2, ymax = cy + h/2),
            fill = col_anal_fill, color = col_anal, linewidth = 0.55) +
  annotate("text", x = anal_box$cx, y = anal_box$cy + anal_box$h/2 - 3.2,
           label = "Primary donor-level inference framework",
           size = 3.1, fontface = "bold", color = col_anal,
           family = base_fontfamily) +
  # ---- ANALYSIS \u5B50\u6846 ----
  geom_rect(data = anal_subs,
            aes(xmin = cx - w/2, xmax = cx + w/2,
                ymin = cy - h/2, ymax = cy + h/2),
            fill = col_box_fill, color = col_anal, linewidth = 0.5) +
  geom_text(data = anal_subs,
            aes(x = cx, y = cy + 5, label = title),
            size = 2.7, fontface = "bold", color = col_anal,
            family = base_fontfamily, lineheight = 1.0) +
  geom_text(data = anal_subs,
            aes(x = cx, y = cy - 5, label = detail),
            size = 2.35, color = "grey20",
            family = base_fontfamily, lineheight = 1.05) +
  # ---- \u4E3B\u6D41 \u77E9\u5F62\u6846 ----
  geom_rect(data = rect_boxes,
            aes(xmin = cx - w/2, xmax = cx + w/2,
                ymin = cy - h/2, ymax = cy + h/2,
                fill = fill, color = border),
            linewidth = 0.5) +
  # ---- \u5E73\u884C\u56DB\u8FB9\u5F62 ----
  geom_polygon(data = para_poly,
               aes(x = x, y = y, group = group_id,
                   fill = fill_col, color = border_col),
               linewidth = 0.5) +
  # ---- \u83F1\u5F62 ----
  geom_polygon(data = diam_poly,
               aes(x = x, y = y, group = group_id,
                   fill = fill_col, color = border_col),
               linewidth = 0.5) +
  # ---- \u4E3B\u6D41 \u6807\u9898 ----
  geom_text(data = main_boxes,
            aes(x = cx, y = title_y, label = title),
            size = 2.85, fontface = "bold", color = "grey15",
            family = base_fontfamily, lineheight = 1.0) +
  # ---- \u4E3B\u6D41 detail ----
  geom_text(data = main_boxes[main_boxes$has_detail, ],
            aes(x = cx, y = detail_y, label = detail),
            size = 2.3, color = "grey25",
            family = base_fontfamily, lineheight = 1.05) +
  # ---- Exclusion \u6846 ----
  geom_rect(data = excl_boxes,
            aes(xmin = cx - w/2, xmax = cx + w/2,
                ymin = cy - h/2, ymax = cy + h/2),
            fill = col_excl_fill, color = col_excl_border, linewidth = 0.5) +
  geom_text(data = excl_boxes,
            aes(x = cx, y = title_y, label = title),
            size = 2.55, fontface = "bold", color = col_excl_border,
            family = base_fontfamily, lineheight = 1.0) +
  geom_text(data = excl_boxes,
            aes(x = cx, y = detail_y, label = detail),
            size = 2.2, color = "grey20",
            family = base_fontfamily, lineheight = 1.05) +
  # ---- \u4E3B\u5782\u76F4\u7BAD\u5934 ----
  geom_segment(data = arrows_v,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               linewidth = 0.5, color = "grey20") +
  geom_segment(data = arrow_to_anal,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               linewidth = 0.5, color = "grey20") +
  # ---- \u6C34\u5E73\u7BAD\u5934\u5230 excl ----
  geom_segment(data = arrows_h,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               linewidth = 0.5, color = "grey20") +
  # ---- Step1 \u2192 Step2 ----
  geom_segment(data = arrow_step,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               linewidth = 0.5, color = "grey20") +
  # ---- Yes/No \u6807\u7B7E ----
  geom_text(data = arrow_labels,
            aes(x = x, y = y, label = label),
            size = 2.4, fontface = "italic", color = "grey25",
            family = base_fontfamily) +
  scale_fill_identity() +
  scale_color_identity() +
  coord_fixed(xlim = c(0, fig_w_mm), ylim = c(0, fig_h_mm), expand = FALSE) +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    plot.background = element_rect(fill = "white", color = NA)
  )
###############################################################################
# 10. \u4FDD\u5B58 (PNG 1000 dpi + cairo_pdf \u77E2\u91CF + \u5B57\u4F53\u5D4C\u5165)
###############################################################################
fig_w_in <- fig_w_mm / 25.4
fig_h_in <- fig_h_mm / 25.4
out_png <- file.path(out_dir, "Figure1_workflow.png")
out_pdf <- file.path(out_dir, "Figure1_workflow.pdf")
# PNG : ragg \u9AD8\u6E05
ragg::agg_png(out_png, width = fig_w_in, height = fig_h_in,
              units = "in", res = png_dpi, background = "white")
print(p); dev.off()
# PDF : cairo_pdf \u77E2\u91CF + \u5B57\u4F53\u5D4C\u5165
if (capabilities("cairo")) {
  ggsave(out_pdf, p, width = fig_w_in, height = fig_h_in,
         units = "in", device = cairo_pdf, bg = "white")
} else {
  ggsave(out_pdf, p, width = fig_w_in, height = fig_h_in,
         units = "in", device = "pdf", bg = "white")
  # \u540E\u7EED\u5B57\u4F53\u5D4C\u5165\u53EF\u7528 embedFonts(out_pdf) \u624B\u52A8\u5904\u7406
}
png_size_mb <- file.size(out_png) / 1024 / 1024
pdf_size_kb <- file.size(out_pdf) / 1024
message("\n\u2705 Figure 1 \u5DF2\u751F\u6210\uFF1A")
message("   PNG : ", out_png, sprintf("  (%.2f MB)", png_size_mb))
message("   PDF : ", out_pdf, sprintf("  (%.1f KB)", pdf_size_kb))
message("   \u5C3A\u5BF8\uFF1A", fig_w_mm, " mm \u00D7 ", fig_h_mm, " mm")
if (png_size_mb > 10) {
  message("   \u26A0 PNG \u8D85\u8FC7 MN \u4E3B\u56FE 10 MB \u4E0A\u9650\uFF1B\u63D0\u4EA4\u7528 PDF (\u77E2\u91CF\u3001\u6781\u5C0F)")
}
