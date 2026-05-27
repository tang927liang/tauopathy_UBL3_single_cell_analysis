###############################################################################
# Figure 5
# Donor-level conditional UBL3 expression is preserved in the three focal
# tauopathy compartments.
#
#   Panel A: cell-level UBL3 density histograms (disease vs control), 3 rows
#            (AD-Excitatory, PSP-Excitatory, PSP-Inhibitory; shared y-axis).
#   Panel B: donor-level median UBL3 strip dotplots for the same 3 comparisons
#            (shared log1p(CP10k) y-axis; within-unit BH-FDR Mann-Whitney U).
#
# Input   : syn52082747 Seurat object (.rds) at the data-project path in
#           section 0 below (confirm it exists before running).
# Outputs : <REPO>/output/figures/Fig5/ (1000 dpi PNG + vector cairo_pdf).
# MN specs: 170 x 180 mm, Arial, line widths >0.25 pt, colourblind-safe palette.
###############################################################################

rm(list = ls()); gc()
SEED <- 20251023; set.seed(SEED)

# --- Paths -------------------------------------------------------------------
# REPO : repository root; figure outputs go under here (same as Fig2-Fig4 / S3).
# Input Seurat objects and metadata CSVs are read from their data-project
# locations below; those are source data and are left unchanged.
REPO <- "D:/RNA/Code/UBL3_tauopathy"
Sys.setenv(LANG = "en"); options(stringsAsFactors = FALSE)

for (pkg in c("lemon", "ggtext", "cowplot")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(dplyr);  library(ggplot2);      library(ragg)
  library(grid);   library(lemon);        library(ggtext); library(cowplot)
  library(org.Hs.eg.db); library(AnnotationDbi)
})
if (.Platform$OS.type == "windows") {
  tryCatch(windowsFonts(Arial = windowsFont("Arial")),
           error = function(e) message("Arial 字体注册跳过"))
}

###############################################################################
# 0. 路径 + 参数
###############################################################################
res_dir <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3"
obj_fp  <- file.path(res_dir, "stepH_slim_uncompressed.rds")
stopifnot(file.exists(obj_fp))

out_dir <- file.path(REPO, "output", "figures", "Fig5")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

positive_set <- list(
  list(disease = "AD",  celltype = "Excitatory neurons"),
  list(disease = "PSP", celltype = "Excitatory neurons"),
  list(disease = "PSP", celltype = "Inhibitory neurons")
)

gene_symbol  <- "UBL3"
compare_ctrl <- "Control"
binwidth_use <- 0.15
png_dpi      <- 1000
base_fontfamily <- "Arial"
palette_groups <- c("AD"="#E69F00","Control"="#56B4E9","FTD"="#009E73","PSP"="#CC79A7")

###############################################################################
# 1. 工具函数
###############################################################################
normalize_celltype6 <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("Astro","Astrocyte","Astrocytes")]             <- "Astrocytes"
  x[x %in% c("Endo","Endothelial","Endothelial cells")]     <- "Endothelial"
  x[x %in% c("Excit","Excitatory","Excitatory neurons")]    <- "Excitatory neurons"
  x[x %in% c("Inhib","Inhibitory","Inhibitory neurons")]    <- "Inhibitory neurons"
  x[x %in% c("Microgl","Micro","Microglia")]                <- "Microglia"
  x[x %in% c("Oligo","Oligodendrocyte","Oligodendrocytes")] <- "Oligodendrocytes"
  factor(x, levels = c("Astrocytes","Endothelial","Excitatory neurons",
                       "Inhibitory neurons","Microglia","Oligodendrocytes"))
}

get_counts_matrix_allcells <- function(obj, assay = "RNA") {
  a <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(a), error = function(e) character(0))
  count_layers <- layers[grepl("^counts", layers)]
  if (length(count_layers) > 0) {
    mats <- list()
    for (ly in count_layers) {
      m <- tryCatch(SeuratObject::LayerData(a, layer = ly), error = function(e) NULL)
      if (!is.null(m) && !is.null(dim(m)) && length(dim(m)) == 2 && ncol(m) > 0) mats[[ly]] <- m
    }
    if (length(mats) == 0) stop("counts layers 无法读取。")
    ref_genes <- rownames(mats[[1]])
    for (nm in names(mats)) {
      if (!identical(rownames(mats[[nm]]), ref_genes)) {
        m0 <- mats[[nm]]
        m_align <- Matrix::Matrix(0, nrow = length(ref_genes), ncol = ncol(m0), sparse = TRUE)
        rownames(m_align) <- ref_genes; colnames(m_align) <- colnames(m0)
        common <- intersect(ref_genes, rownames(m0))
        m_align[common, ] <- m0[common, , drop = FALSE]
        mats[[nm]] <- m_align
      }
    }
    mat_all <- if (length(mats) == 1) mats[[1]] else Reduce(Matrix::cbind2, mats)
    dup <- duplicated(colnames(mat_all))
    if (any(dup)) mat_all <- mat_all[, !dup, drop = FALSE]
    miss <- setdiff(colnames(obj), colnames(mat_all))
    if (length(miss) > 0) {
      m_fill <- Matrix::Matrix(0, nrow = nrow(mat_all), ncol = length(miss), sparse = TRUE)
      rownames(m_fill) <- rownames(mat_all); colnames(m_fill) <- miss
      mat_all <- Matrix::cbind2(mat_all, m_fill)
    }
    return(mat_all[, colnames(obj), drop = FALSE])
  }
  m4 <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "counts"),
                 error = function(e) NULL)
  if (!is.null(m4)) return(m4)
  stop("无法获取 counts 矩阵。")
}

locate_gene_row <- function(counts_mat, gene_symbol) {
  rn <- rownames(counts_mat)
  if (gene_symbol %in% rn) return(gene_symbol)
  hit_ci <- which(toupper(rn) == toupper(gene_symbol))
  if (length(hit_ci) >= 1) return(rn[hit_ci[1]])
  ens_tbl <- tryCatch(
    AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbol,
                          keytype = "SYMBOL", columns = "ENSEMBL"),
    error = function(e) NULL)
  if (!is.null(ens_tbl)) {
    ens_ids <- unique(ens_tbl$ENSEMBL[!is.na(ens_tbl$ENSEMBL)])
    h1 <- intersect(ens_ids, rn); if (length(h1) >= 1) return(h1[1])
    rn_strip <- sub("\\.\\d+$", "", rn); h2 <- which(rn_strip %in% ens_ids)
    if (length(h2) >= 1) return(rn[h2[1]])
  }
  stop(paste0("counts 行名中找不到基因：", gene_symbol))
}

###############################################################################
# 2. 读 syn52082747 + 准备数据
###############################################################################
obj <- readRDS(obj_fp); DefaultAssay(obj) <- "RNA"
md  <- obj@meta.data
ct_col  <- intersect(c("celltype6","celltype","cell_type"), colnames(md))[1]
grp_col <- intersect(c("group4","group2","group","Group","Diagnosis","diagnosis"), colnames(md))[1]
don_col <- intersect(c("autopsy_id","donor","SampleID","sample","orig.ident"), colnames(md))[1]
stopifnot(!is.na(ct_col), !is.na(grp_col), !is.na(don_col))

ctrl_alias <- c("Control","CTRL","Ctr","CTR","NC","Normal","N",
                "control","ctrl","ctr","nc","normal")
md$celltype6 <- normalize_celltype6(md[[ct_col]])
md$group_std <- ifelse(trimws(as.character(md[[grp_col]])) %in% ctrl_alias,
                       compare_ctrl, trimws(as.character(md[[grp_col]])))
md$donor <- trimws(as.character(md[[don_col]]))
obj@meta.data <- md

counts_mat <- get_counts_matrix_allcells(obj, assay = "RNA")
gene_row   <- locate_gene_row(counts_mat, gene_symbol)
lib_size   <- Matrix::colSums(counts_mat)
expr_vec   <- log1p((as.numeric(counts_mat[gene_row, , drop = TRUE]) /
                      pmax(lib_size, 1)) * 1e4)

df_all <- data.frame(
  expr = expr_vec, donor = md$donor, group = md$group_std,
  celltype = md$celltype6, stringsAsFactors = FALSE)

###############################################################################
# 3. Within-unit BH FDR 计算
###############################################################################
compute_fdr_table <- function(df_all, disease, control = compare_ctrl) {
  df <- df_all %>% filter(group %in% c(disease, control), expr > 0)
  donor_med <- df %>%
    group_by(celltype, donor, group) %>%
    summarise(median_expr = median(expr), .groups = "drop")
  stats <- donor_med %>%
    group_by(celltype) %>%
    summarise(
      p_raw = {
        g <- group; v <- median_expr
        if (length(unique(g)) < 2) NA_real_
        else tryCatch(wilcox.test(v ~ g, exact = FALSE)$p.value,
                      error = function(e) NA_real_)
      }, .groups = "drop")
  stats$padj <- p.adjust(stats$p_raw, method = "BH")
  stats
}

fdr_AD  <- compute_fdr_table(df_all, "AD")
fdr_PSP <- compute_fdr_table(df_all, "PSP")
fdr_lookup <- list(AD = fdr_AD, PSP = fdr_PSP)

###############################################################################
# 4. ★ v3 关键：预计算 3 个 hist + 3 个 dot 的全局 Y 范围
###############################################################################
# 收集 3 个 (disease, celltype) 对应的所有 cell-level expr (>0)
collect_data_for_panels <- function(df_all, positive_set) {
  hist_data_list <- list()
  dot_data_list  <- list()
  for (cfg in positive_set) {
    key <- paste0(cfg$disease, "_", cfg$celltype)
    
    # Cell-level data for histograms
    df_pair <- df_all %>%
      filter(group %in% c(cfg$disease, compare_ctrl),
             celltype == cfg$celltype, expr > 0)
    hist_data_list[[key]] <- df_pair
    
    # Donor-level data for dotplots
    donor_median <- df_pair %>%
      group_by(donor, group) %>%
      summarise(median_expr = median(expr), .groups = "drop")
    dot_data_list[[key]] <- donor_median
  }
  list(hist = hist_data_list, dot = dot_data_list)
}

panel_data <- collect_data_for_panels(df_all, positive_set)

# 计算 3 个 histogram 的 density Y max
compute_hist_y_max <- function(hist_data_list, binwidth) {
  y_max_all <- 0
  for (df_pair in hist_data_list) {
    for (grp in unique(df_pair$group)) {
      x_sub <- df_pair$expr[df_pair$group == grp]
      if (length(x_sub) < 2) next
      h <- hist(x_sub, breaks = seq(0, ceiling(max(x_sub)/binwidth)*binwidth + binwidth, by = binwidth),
                plot = FALSE)
      y_max_all <- max(y_max_all, h$density, na.rm = TRUE)
    }
  }
  ceiling(y_max_all * 10) / 10 * 1.05   # round up + 5% headroom
}

hist_y_max <- compute_hist_y_max(panel_data$hist, binwidth_use)
message(sprintf("📐 Histograms shared Y max: %.3f", hist_y_max))

# 计算 3 个 dotplot 的 UBL3 log1p(CP10k) global Y range
all_dot_vals <- unlist(lapply(panel_data$dot, function(d) d$median_expr))
dot_y_min <- min(all_dot_vals, na.rm = TRUE)
dot_y_max <- max(all_dot_vals, na.rm = TRUE)
dot_range <- dot_y_max - dot_y_min
dot_y_lim <- c(dot_y_min - dot_range * 0.08, dot_y_max + dot_range * 0.12)
message(sprintf("📐 Dotplots shared Y range: [%.3f, %.3f]", dot_y_lim[1], dot_y_lim[2]))

###############################################################################
# 5. 主题
###############################################################################
hist_theme <- theme_bw(base_size = 7.5, base_family = base_fontfamily) +
  theme(
    plot.title       = element_text(face = "bold", size = 7.8, hjust = 0,
                                    margin = margin(b = 2)),
    axis.title.x     = element_text(size = 7.2, margin = margin(t = 2)),
    axis.title.y     = element_text(size = 7.2, margin = margin(r = 2)),
    axis.text        = element_text(size = 6.6, colour = "black"),
    axis.ticks       = element_line(colour = "black", linewidth = 0.30),
    panel.border     = element_rect(colour = "grey35", fill = NA, linewidth = 0.35),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.22),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom", legend.direction = "horizontal",
    legend.title     = element_blank(), legend.text = element_text(size = 6.8),
    legend.key.width = unit(0.7,"lines"), legend.key.height = unit(0.7,"lines"),
    legend.margin    = margin(t = -2, b = 0),
    plot.margin      = margin(t = 4, r = 5, b = 2, l = 4)
  )

dot_theme <- theme_bw(base_size = 7.5, base_family = base_fontfamily) +
  theme(
    plot.title       = element_text(face = "bold", size = 7.8, hjust = 0,
                                    margin = margin(b = 2)),
    plot.subtitle    = element_text(face = "italic", size = 7.0, colour = "grey25",
                                    hjust = 0, margin = margin(b = 3)),
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(size = 7.2, margin = margin(r = 3)),
    axis.text.x      = element_text(size = 7.0, colour = "black"),
    axis.text.y      = element_text(size = 6.8, colour = "black"),
    axis.ticks       = element_line(colour = "black", linewidth = 0.30),
    panel.border     = element_rect(colour = "grey35", fill = NA, linewidth = 0.35),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.22),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom", legend.direction = "horizontal",
    legend.title     = element_blank(), legend.text = element_text(size = 6.8),
    legend.key.width = unit(0.7,"lines"), legend.key.height = unit(0.7,"lines"),
    legend.margin    = margin(t = -2, b = 0),
    plot.margin      = margin(t = 4, r = 5, b = 2, l = 4)
  )

###############################################################################
# 6. Sub-plot 生成函数（注入全局 Y limits）
###############################################################################
make_hist <- function(df_pair, disease, celltype, y_max_shared) {
  donor_pair <- unique(df_pair[, c("donor","group")])
  n_dis <- sum(donor_pair$group == disease)
  n_ctl <- sum(donor_pair$group == compare_ctrl)
  lab_dis <- sprintf("%s (n=%d)", disease, n_dis)
  lab_ctl <- sprintf("%s (n=%d)", compare_ctrl, n_ctl)
  df_ctrl <- df_pair[df_pair$group == compare_ctrl,, drop = FALSE]
  df_dis  <- df_pair[df_pair$group == disease,     , drop = FALSE]
  
  ggplot() +
    geom_histogram(data = df_ctrl,
      aes(x = expr, y = after_stat(density), fill = lab_ctl, color = lab_ctl),
      binwidth = binwidth_use, position = "identity",
      alpha = 0.92, linewidth = 0.20, boundary = 0, closed = "left") +
    geom_histogram(data = df_dis,
      aes(x = expr, y = after_stat(density), fill = lab_dis, color = lab_dis),
      binwidth = binwidth_use, position = "identity",
      alpha = 0.78, linewidth = 0.20, boundary = 0, closed = "left") +
    scale_fill_manual(values = setNames(c("white","grey65"), c(lab_dis, lab_ctl)),
                      breaks = c(lab_dis, lab_ctl)) +
    scale_color_manual(values = setNames(c("black","grey40"), c(lab_dis, lab_ctl)),
                       breaks = c(lab_dis, lab_ctl)) +
    # ★ 共享 Y 上限
    scale_y_continuous(limits = c(0, y_max_shared),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = paste0(gene_symbol, " log1p(CP10k)"), y = "Density",
         title = sprintf("%s — %s", disease, celltype)) +
    guides(color = "none",
           fill = guide_legend(nrow = 1, byrow = TRUE,
                               override.aes = list(color = c("black","grey40"),
                                                   fill  = c("white","grey65")))) +
    hist_theme
}

make_dot <- function(donor_median, disease, celltype, fdr_value, y_lim_shared) {
  n_dis <- sum(donor_median$group == disease)
  n_ctl <- sum(donor_median$group == compare_ctrl)
  lab_dis <- sprintf("%s (n=%d)", disease, n_dis)
  lab_ctl <- sprintf("%s (n=%d)", compare_ctrl, n_ctl)
  fdr_str <- if (is.finite(fdr_value)) sprintf("MWU FDR = %.2e", fdr_value) else "MWU FDR = NA"
  
  summary_bars <- donor_median %>%
    group_by(group) %>%
    summarise(med = median(median_expr),
              q1 = quantile(median_expr, 0.25),
              q3 = quantile(median_expr, 0.75), .groups = "drop")
  donor_median$group <- factor(donor_median$group, levels = c(disease, compare_ctrl))
  summary_bars$group <- factor(summary_bars$group, levels = c(disease, compare_ctrl))
  
  ggplot() +
    geom_linerange(data = summary_bars,
                   aes(x = group, ymin = q1, ymax = q3),
                   colour = "grey45", linewidth = 0.38) +
    geom_errorbar(data = summary_bars,
                  aes(x = group, ymin = med, ymax = med),
                  width = 0.45, linewidth = 0.75, colour = "black") +
    geom_point(data = donor_median,
               aes(x = group, y = median_expr, fill = group),
               position = position_jitter(width = 0.14, height = 0, seed = SEED),
               size = 1.8, shape = 21, stroke = 0.30,
               colour = "black", alpha = 0.92) +
    scale_fill_manual(values = palette_groups,
                      breaks = c(disease, compare_ctrl),
                      labels = c(lab_dis, lab_ctl)) +
    scale_x_discrete(labels = c(disease, compare_ctrl)) +
    # ★ 共享 Y 范围
    scale_y_continuous(limits = y_lim_shared, expand = c(0, 0)) +
    labs(y = paste0(gene_symbol, " log1p(CP10k)"),
         title = sprintf("%s — %s", disease, celltype),
         subtitle = fdr_str) +
    coord_cartesian(clip = "off") +
    dot_theme
}

###############################################################################
# 7. 生成 6 个 sub-plot
###############################################################################
hist_plots <- list()
dot_plots  <- list()
for (i in seq_along(positive_set)) {
  cfg <- positive_set[[i]]
  key <- paste0(cfg$disease, "_", cfg$celltype)
  
  fdr_tbl <- fdr_lookup[[cfg$disease]]
  this_fdr <- fdr_tbl$padj[as.character(fdr_tbl$celltype) == cfg$celltype]
  if (length(this_fdr) == 0) this_fdr <- NA_real_
  message(sprintf("📌 row %d: %s — %s | BH FDR = %.3e",
                  i, cfg$disease, cfg$celltype, this_fdr))
  
  hist_plots[[i]] <- make_hist(panel_data$hist[[key]], cfg$disease, cfg$celltype, hist_y_max)
  dot_plots[[i]]  <- make_dot (panel_data$dot[[key]],  cfg$disease, cfg$celltype, this_fdr, dot_y_lim)
}

###############################################################################
# 8. 组合 — A/B 标签变小、定位到 corner 不压图
###############################################################################
panel_A <- plot_grid(hist_plots[[1]], hist_plots[[2]], hist_plots[[3]],
                     ncol = 1, align = "v")
panel_B <- plot_grid(dot_plots[[1]], dot_plots[[2]], dot_plots[[3]],
                     ncol = 1, align = "v")

combined <- plot_grid(
  panel_A, panel_B,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 10,                     # ★ 缩小: 14 -> 10
  label_fontfamily = base_fontfamily,
  label_fontface = "bold",
  label_x = 0.01,                      # ★ 略左移
  label_y = 0.995,                     # ★ 顶部留白
  hjust = 0, vjust = 1,                # 锚定左上角
  rel_widths = c(1, 1)
)

###############################################################################
# 9. 保存（PNG 1000dpi + 矢量 PDF）
###############################################################################
fig_width_mm  <- 170
fig_height_mm <- 180
fig_width_in  <- fig_width_mm  / 25.4
fig_height_in <- fig_height_mm / 25.4

out_png <- file.path(out_dir, "Figure5_MN_v3_1000dpi.png")
out_pdf <- file.path(out_dir, "Figure5_MN_v3_vector.pdf")

ragg::agg_png(out_png, width = fig_width_in, height = fig_height_in,
              units = "in", res = png_dpi, background = "white")
print(combined); dev.off()

if (capabilities("cairo")) {
  ggsave(out_pdf, combined, width = fig_width_in, height = fig_height_in,
         units = "in", device = cairo_pdf, bg = "white")
} else {
  ggsave(out_pdf, combined, width = fig_width_in, height = fig_height_in,
         units = "in", device = "pdf", bg = "white")
}

message("\n🎉 Figure 5 (MN v3) 已生成：")
message("   PNG  : ", out_png)
message("   PDF  : ", out_pdf)
message("   尺寸 : ", fig_width_mm, " mm × ", fig_height_mm, " mm")
message("   Y 轴统一: histogram [0, ", round(hist_y_max, 3),
        "], dotplot [", round(dot_y_lim[1], 3), ", ", round(dot_y_lim[2], 3), "]")

