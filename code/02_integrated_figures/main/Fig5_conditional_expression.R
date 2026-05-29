###############################################################################
# 统一脚本：Figure 5 + Supplementary Figure S4（MN 投稿 unified v2 版）
#
# 【v2 相对 v1 的改动】（修复 Section B 崩溃）
#   报错: ! `.homonyms` must be a string or character vector.
#   原因: patchwork 的 `wrap_plots() + plot_annotation() & theme()` 这套写法
#         在同一 R session 里第 2 次调用会触发 vctrs/patchwork 主题合并 bug
#         （与 ragg 'res' 锁死同类的状态污染问题）。v1 在 Stage A 用 cowplot
#         跑通 Figure 5 + Stage B Section A 用 patchwork 也跑通，但 Section B
#         作为 patchwork 的第二次调用就翻车。
#   修法: Stage B 的 section 组装全部改用 cowplot::plot_grid()，删掉所有
#         patchwork::plot_annotation() 和 `&` 操作。section 标题用
#         cowplot::ggdraw() + draw_label() 渲染；A/B 标签直接用 plot_grid 的
#         labels 参数。这样 Stage A 和 Stage B 都只用 cowplot 组装，状态干净。
#   其他逻辑保持 v1 不变（共享主题 / 函数 / 调色板 / descriptive 处理 /
#   稳健 PNG 导出 / 缓存）。
#
# 【v1 设计目标】（保持不变）
#   - 同一份 hist_theme / dot_theme（base_size = 7.5）
#   - 同一份 palette_groups（AD=#E69F00 / Control=#56B4E9 / FTD=#009E73 / PSP=#CC79A7）
#   - 同一份 make_hist() / make_dot()（descriptive: n<4 自动显示 "Descriptive only"）
#   - 同一份 binwidth（0.15）/ point size（1.8）/ shape（圆点带描边）/ Arial 字体
#   - syn52082747 两阶段共享缓存
#   - 稳健 PNG 导出（grDevices::png type="cairo" + tryCatch + safe_close_all_devices）
#
# 【输出】
#   Figure 5 → D:/RNA/MNversion submission/Figure 5.../Results/
#     - Figure5_MN_unified_1000dpi.png
#     - Figure5_MN_unified_vector.pdf
#   Supp Fig S4 → D:/RNA/MNversion submission/Supplementary Figure S4.../Results2/
#     - SuppS4_SectionA-G_{cohort}_{disease}_1000dpi.png  (7 个)
#     - SuppS4_SectionA-G_{cohort}_{disease}_vector.pdf   (7 个)
#     - Supplementary_Figure_S4_MN_unified_ALL_vector.pdf  ← 投稿主用
#     - Supplementary_Figure_S4_MN_unified_ALL_1000dpi.png
###############################################################################

rm(list = ls()); gc()
SEED <- 20251023; set.seed(SEED)
Sys.setenv(LANG = "en"); options(stringsAsFactors = FALSE)

for (pkg in c("lemon", "ggtext", "cowplot", "magick", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(dplyr);  library(ggplot2);      library(ragg)
  library(grid);   library(lemon);        library(ggtext); library(cowplot)
  library(patchwork)   # 仍然加载（合并多页 PDF 时可能用），但 section 组装不再用
  library(org.Hs.eg.db); library(AnnotationDbi)
})
if (.Platform$OS.type == "windows") {
  tryCatch(windowsFonts(Arial = windowsFont("Arial")),
           error = function(e) message("Arial 字体注册跳过"))
}

###############################################################################
# 0. 输出目录
###############################################################################
FIG5_DIR <- "D:/RNA/MNversion submission/Figure 5. Donor-level conditional UBL3 expression preserves baseline magnitude across tauopathies/Results"
SUPS4_DIR <- "D:/RNA/MNversion submission/Supplementary Figure S4. negative overlap and dotplot/Results2"
dir.create(FIG5_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(SUPS4_DIR, showWarnings = FALSE, recursive = TRUE)

###############################################################################
# 1. 共享参数
###############################################################################
gene_symbol     <- "UBL3"
compare_ctrl    <- "Control"
binwidth_use    <- 0.15
png_dpi         <- 1000
base_fontfamily <- "Arial"

MIN_DONORS_PER_GROUP <- 1
DOT_POINT_SIZE       <- 1.8

palette_groups  <- c("AD"="#E69F00", "Control"="#56B4E9",
                     "FTD"="#009E73", "PSP"="#CC79A7")

all_celltype_levels <- c("Astrocytes","Endothelial","Excitatory neurons",
                         "Inhibitory neurons","Microglia","Oligodendrocytes")

ctrl_alias <- c("Control","CTRL","Ctr","CTR","NC","Normal","N",
                "control","ctrl","ctr","nc","normal")

###############################################################################
# 2. cohort 配置
###############################################################################
USER_CONFIG <- list(
  syn52082747 = list(
    rds_path     = "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3/stepH_slim_uncompressed.rds",
    meta_csv     = NA_character_,
    use_csv      = FALSE,
    ct_col_pref  = c("celltype6","celltype","cell_type"),
    grp_col_pref = c("group4","group2","group","Group","Diagnosis","diagnosis"),
    don_col_pref = c("autopsy_id","donor","SampleID","sample","orig.ident","subject")
  ),
  GSE157827 = list(
    rds_path     = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results/stepH_obj_celltype6_named.rds",
    meta_csv     = NA_character_,
    use_csv      = FALSE,
    ct_col_pref  = c("celltype6","celltype","cell_type","celltype6_named","celltype6_manual"),
    grp_col_pref = c("group4","group2","group","Group","diagnosis","Dx","clinical_diagnosis"),
    don_col_pref = c("autopsy_id","donor","Donor","sample","Sample","orig.ident","patient","subject")
  ),
  GSE174367 = list(
    rds_path     = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/stepH_obj_celltype6_named.rds",
    meta_csv     = NA_character_,
    use_csv      = FALSE,
    ct_col_pref  = c("celltype6","celltype","cell_type","celltype6_named","celltype6_manual"),
    grp_col_pref = c("group4","group2","group","Group","diagnosis","Dx","clinical_diagnosis"),
    don_col_pref = c("autopsy_id","donor","Donor","sample","Sample","orig.ident","patient","subject")
  ),
  syn21788402_EC = list(
    rds_path     = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds",
    meta_csv     = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepP_syn21788402_matched_cells_meta.csv",
    use_csv      = TRUE,
    ct_col_pref  = c("celltype6","celltype","cell_type","CellType","broad_celltype"),
    don_col_pref = c("PatientID","autopsy_id","donor","donor_id","SampleID","sample","SubID","subject","orig.ident"),
    join_key_pref = c("SampleID","sample","orig.ident")
  ),
  syn21788402_SFG = list(
    rds_path     = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepH_syn21788402_SFG_obj_celltype6.rds",
    meta_csv     = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepP_syn21788402_SFG_matched_cells_meta.csv",
    use_csv      = TRUE,
    ct_col_pref  = c("celltype6","celltype","cell_type","CellType","broad_celltype"),
    don_col_pref = c("PatientID","autopsy_id","donor","donor_id","SampleID","sample","SubID","subject","orig.ident"),
    join_key_pref = c("SampleID","sample","orig.ident")
  )
)

###############################################################################
# 3. 共享工具函数
###############################################################################
normalize_celltype6 <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("Astro","Astrocyte","Astrocytes")]             <- "Astrocytes"
  x[x %in% c("Endo","Endothelial","Endothelial cells")]     <- "Endothelial"
  x[x %in% c("Excit","Excitatory","Excitatory neurons")]    <- "Excitatory neurons"
  x[x %in% c("Inhib","Inhibitory","Inhibitory neurons")]    <- "Inhibitory neurons"
  x[x %in% c("Microgl","Micro","Microglia")]                <- "Microglia"
  x[x %in% c("Oligo","Oligodendrocyte","Oligodendrocytes")] <- "Oligodendrocytes"
  factor(x, levels = all_celltype_levels)
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

safe_close_all_devices <- function() {
  while (!is.null(dev.list())) {
    try(dev.off(), silent = TRUE)
  }
}

###############################################################################
# 4. 共享 cohort loader（带缓存）
###############################################################################
cohort_cache <- list()

load_cohort_df <- function(cohort_key) {
  if (!is.null(cohort_cache[[cohort_key]])) {
    message(sprintf("  ▶ 命中缓存: %s", cohort_key))
    return(cohort_cache[[cohort_key]])
  }
  
  cfg <- USER_CONFIG[[cohort_key]]
  stopifnot(file.exists(cfg$rds_path))
  
  message(sprintf("  ▼ 读对象: %s", cohort_key))
  obj <- readRDS(cfg$rds_path); DefaultAssay(obj) <- "RNA"
  md  <- obj@meta.data
  md$cell_id <- colnames(obj)
  rownames(md) <- md$cell_id
  
  ct_col  <- intersect(cfg$ct_col_pref,  colnames(md))[1]
  don_col <- intersect(cfg$don_col_pref, colnames(md))[1]
  stopifnot(!is.na(ct_col), !is.na(don_col))
  
  if (cfg$use_csv) {
    join_key <- intersect(cfg$join_key_pref, colnames(md))[1]
    stopifnot(!is.na(join_key), file.exists(cfg$meta_csv))
    
    meta_raw <- read.csv(cfg$meta_csv, stringsAsFactors = FALSE)
    stopifnot(all(c("sample","group") %in% colnames(meta_raw)))
    
    grp_raw <- trimws(as.character(meta_raw$group))
    meta_raw$group_std <- ifelse(grp_raw %in% ctrl_alias, compare_ctrl, grp_raw)
    
    map_by_sample <- meta_raw %>%
      transmute(sample_join = trimws(as.character(sample)),
                group_std   = group_std) %>%
      filter(!is.na(sample_join), sample_join != "") %>%
      distinct()
    
    md$sample_join_in_obj <- trimws(as.character(md[[join_key]]))
    md1 <- dplyr::left_join(md, map_by_sample,
                            by = c("sample_join_in_obj" = "sample_join"))
    rownames(md1) <- md1$cell_id
    
    md1$celltype6 <- normalize_celltype6(md1[[ct_col]])
    md1$donor     <- trimws(as.character(md1[[don_col]]))
    
    keep <- !is.na(md1$group_std) & md1$group_std != "" &
      !is.na(md1$celltype6) & !is.na(md1$donor) & md1$donor != ""
    md1 <- md1[keep, , drop = FALSE]
    cells_keep <- md1$cell_id
    obj_meta <- md1
  } else {
    grp_col <- intersect(cfg$grp_col_pref, colnames(md))[1]
    stopifnot(!is.na(grp_col))
    
    md$celltype6 <- normalize_celltype6(md[[ct_col]])
    md$group_std <- ifelse(trimws(as.character(md[[grp_col]])) %in% ctrl_alias,
                           compare_ctrl, trimws(as.character(md[[grp_col]])))
    md$donor <- trimws(as.character(md[[don_col]]))
    
    keep <- !is.na(md$group_std) & md$group_std != "" &
      !is.na(md$celltype6) & !is.na(md$donor) & md$donor != ""
    md1 <- md[keep, , drop = FALSE]
    cells_keep <- md1$cell_id
    obj_meta <- md1
  }
  
  counts_mat <- get_counts_matrix_allcells(obj, assay = "RNA")
  gene_row   <- locate_gene_row(counts_mat, gene_symbol)
  j_keep     <- match(cells_keep, colnames(counts_mat))
  cm_sub     <- counts_mat[, j_keep, drop = FALSE]
  lib_size   <- Matrix::colSums(cm_sub)
  expr_vec   <- log1p((as.numeric(cm_sub[gene_row, , drop = TRUE]) /
                         pmax(lib_size, 1)) * 1e4)
  
  df <- data.frame(
    expr = expr_vec, donor = obj_meta$donor, group = obj_meta$group_std,
    celltype = obj_meta$celltype6, stringsAsFactors = FALSE
  )
  message(sprintf("    Cells loaded: %d | celltype col=%s | donor col=%s",
                  nrow(df), ct_col, don_col))
  rm(obj, counts_mat, cm_sub, md); gc()
  cohort_cache[[cohort_key]] <<- df
  df
}

compute_fdr_table <- function(df_all, disease) {
  df <- df_all %>% filter(group %in% c(disease, compare_ctrl), expr > 0)
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

###############################################################################
# 5. 共享 themes
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
# 6. 共享 sub-plot 生成函数
###############################################################################
make_hist <- function(df_pair, disease, celltype, y_max_shared) {
  if (nrow(df_pair) == 0) return(NULL)
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
  if (nrow(donor_median) == 0) return(NULL)
  n_dis <- sum(donor_median$group == disease)
  n_ctl <- sum(donor_median$group == compare_ctrl)
  lab_dis <- sprintf("%s (n=%d)", disease, n_dis)
  lab_ctl <- sprintf("%s (n=%d)", compare_ctrl, n_ctl)
  
  is_descriptive <- (n_dis < 4 || n_ctl < 4)
  if (is_descriptive) {
    fdr_str <- sprintf("Descriptive only (n = %d vs %d)", n_dis, n_ctl)
  } else {
    fdr_str <- if (is.finite(fdr_value)) sprintf("MWU FDR = %.2e", fdr_value) else "MWU FDR = NA"
  }
  
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
               size = DOT_POINT_SIZE, shape = 21, stroke = 0.30,
               colour = "black", alpha = 0.92) +
    scale_fill_manual(values = palette_groups,
                      breaks = c(disease, compare_ctrl),
                      labels = c(lab_dis, lab_ctl)) +
    scale_x_discrete(labels = c(disease, compare_ctrl)) +
    scale_y_continuous(limits = y_lim_shared, expand = c(0, 0)) +
    labs(y = paste0(gene_symbol, " log1p(CP10k)"),
         title = sprintf("%s — %s", disease, celltype),
         subtitle = fdr_str) +
    coord_cartesian(clip = "off") +
    dot_theme
}

compute_hist_y_max <- function(hist_data_list, binwidth) {
  y_max_all <- 0
  for (df_pair in hist_data_list) {
    if (is.null(df_pair) || nrow(df_pair) == 0) next
    for (grp in unique(df_pair$group)) {
      x_sub <- df_pair$expr[df_pair$group == grp]
      if (length(x_sub) < 2) next
      h <- hist(x_sub,
                breaks = seq(0, ceiling(max(x_sub)/binwidth)*binwidth + binwidth, by = binwidth),
                plot = FALSE)
      y_max_all <- max(y_max_all, h$density, na.rm = TRUE)
    }
  }
  if (y_max_all <= 0) y_max_all <- 1
  ceiling(y_max_all * 10) / 10 * 1.05
}

compute_dot_y_lim <- function(dot_data_list) {
  all_vals <- unlist(lapply(dot_data_list,
                            function(d) if (is.null(d)) NULL else d$median_expr))
  if (length(all_vals) == 0) return(c(0, 1))
  v_min <- min(all_vals, na.rm = TRUE)
  v_max <- max(all_vals, na.rm = TRUE)
  v_range <- v_max - v_min
  if (v_range <= 0) v_range <- max(abs(v_max), 1)
  c(v_min - v_range * 0.08, v_max + v_range * 0.12)
}

save_png_pdf <- function(plot_obj, png_path, pdf_path, w_in, h_in, dpi) {
  safe_close_all_devices()
  png_ok <- tryCatch({
    grDevices::png(filename = png_path, width = w_in, height = h_in,
                   units = "in", res = dpi, type = "cairo", bg = "white")
    print(plot_obj); dev.off()
    TRUE
  }, error = function(e) {
    message(sprintf("  ⚠ PNG 导出失败 (%s)，ggsave 兜底", conditionMessage(e)))
    safe_close_all_devices()
    tryCatch({
      ggsave(png_path, plot_obj, width = w_in, height = h_in,
             units = "in", dpi = dpi, bg = "white")
      TRUE
    }, error = function(e2) { safe_close_all_devices(); FALSE })
  })
  safe_close_all_devices()
  pdf_ok <- tryCatch({
    if (capabilities("cairo")) {
      ggsave(pdf_path, plot_obj, width = w_in, height = h_in,
             units = "in", device = cairo_pdf, bg = "white")
    } else {
      ggsave(pdf_path, plot_obj, width = w_in, height = h_in,
             units = "in", device = "pdf", bg = "white")
    }
    TRUE
  }, error = function(e) {
    message(sprintf("  ⚠ PDF 导出失败: %s", conditionMessage(e)))
    safe_close_all_devices(); FALSE
  })
  safe_close_all_devices(); gc()
  c(png = png_ok, pdf = pdf_ok)
}

###############################################################################
# v2 新增：section page 组装函数（纯 cowplot，无 patchwork &）
# 把 hist_list + dot_list + section 标题组装成一页完整 section page
###############################################################################
build_section_page_cowplot <- function(hist_list, dot_list, hdr_text,
                                       font_family = base_fontfamily) {
  # 左列：垂直堆叠所有 histogram
  left_col  <- cowplot::plot_grid(plotlist = unname(hist_list),
                                  ncol = 1, align = "v")
  # 右列：垂直堆叠所有 dotplot
  right_col <- cowplot::plot_grid(plotlist = unname(dot_list),
                                  ncol = 1, align = "v")
  
  # 2 列并排 + A/B 标签（沿用 Figure 5 的 cowplot::plot_grid 风格）
  inner_body <- cowplot::plot_grid(
    left_col, right_col,
    ncol = 2,
    labels = c("A", "B"),
    label_size = 10,
    label_fontfamily = font_family,
    label_fontface = "bold",
    label_x = 0.01, label_y = 0.998,
    hjust = 0, vjust = 1,
    rel_widths = c(1, 1)
  )
  
  # 顶部 section 标题（cowplot::ggdraw + draw_label，避开 patchwork plot_annotation）
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(
      hdr_text,
      fontfamily = font_family, fontface = "bold", size = 8.5,
      x = 0.02, hjust = 0, y = 0.5
    )
  
  # 标题在顶 + 内容在下
  cowplot::plot_grid(
    title_grob, inner_body,
    ncol = 1,
    rel_heights = c(0.04, 0.96)
  )
}


###############################################################################
# =============================================================================
# STAGE A：Figure 5
# =============================================================================
###############################################################################
message("\n\n#############################################################")
message("# STAGE A: Figure 5 — Donor-level conditional UBL3 expression")
message("#############################################################")

positive_set <- list(
  list(disease = "AD",  celltype = "Excitatory neurons"),
  list(disease = "PSP", celltype = "Excitatory neurons"),
  list(disease = "PSP", celltype = "Inhibitory neurons")
)

df_v1 <- load_cohort_df("syn52082747")
fdr_AD  <- compute_fdr_table(df_v1, "AD")
fdr_PSP <- compute_fdr_table(df_v1, "PSP")
fdr_lookup <- list(AD = fdr_AD, PSP = fdr_PSP)

hist_data_f5 <- list()
dot_data_f5  <- list()
for (cfg in positive_set) {
  key <- paste0(cfg$disease, "_", cfg$celltype)
  df_pair <- df_v1 %>%
    filter(group %in% c(cfg$disease, compare_ctrl),
           celltype == cfg$celltype, expr > 0)
  hist_data_f5[[key]] <- df_pair
  
  donor_med <- df_pair %>%
    group_by(donor, group) %>%
    summarise(median_expr = median(expr), .groups = "drop")
  dot_data_f5[[key]] <- donor_med
}

hist_y_max_f5 <- compute_hist_y_max(hist_data_f5, binwidth_use)
dot_y_lim_f5  <- compute_dot_y_lim(dot_data_f5)
message(sprintf("📐 Figure 5 Y 范围: hist [0, %.3f], dot [%.3f, %.3f]",
                hist_y_max_f5, dot_y_lim_f5[1], dot_y_lim_f5[2]))

hist_plots_f5 <- list(); dot_plots_f5 <- list()
for (i in seq_along(positive_set)) {
  cfg <- positive_set[[i]]
  key <- paste0(cfg$disease, "_", cfg$celltype)
  fdr_tbl <- fdr_lookup[[cfg$disease]]
  this_fdr <- fdr_tbl$padj[as.character(fdr_tbl$celltype) == cfg$celltype]
  if (length(this_fdr) == 0) this_fdr <- NA_real_
  message(sprintf("  📌 row %d: %s — %s | BH FDR = %.3e",
                  i, cfg$disease, cfg$celltype, this_fdr))
  hist_plots_f5[[i]] <- make_hist(hist_data_f5[[key]], cfg$disease, cfg$celltype, hist_y_max_f5)
  dot_plots_f5[[i]]  <- make_dot (dot_data_f5[[key]],  cfg$disease, cfg$celltype, this_fdr, dot_y_lim_f5)
}

panel_A_f5 <- cowplot::plot_grid(hist_plots_f5[[1]], hist_plots_f5[[2]], hist_plots_f5[[3]],
                                 ncol = 1, align = "v")
panel_B_f5 <- cowplot::plot_grid(dot_plots_f5[[1]], dot_plots_f5[[2]], dot_plots_f5[[3]],
                                 ncol = 1, align = "v")
fig5 <- cowplot::plot_grid(
  panel_A_f5, panel_B_f5,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 10,
  label_fontfamily = base_fontfamily,
  label_fontface = "bold",
  label_x = 0.01, label_y = 0.995,
  hjust = 0, vjust = 1,
  rel_widths = c(1, 1)
)

fig5_w_mm <- 170; fig5_h_mm <- 180
fig5_w_in <- fig5_w_mm / 25.4; fig5_h_in <- fig5_h_mm / 25.4

fig5_png <- file.path(FIG5_DIR, "Figure5_MN_unified_1000dpi.png")
fig5_pdf <- file.path(FIG5_DIR, "Figure5_MN_unified_vector.pdf")

status_f5 <- save_png_pdf(fig5, fig5_png, fig5_pdf, fig5_w_in, fig5_h_in, png_dpi)
if (all(status_f5)) {
  message(sprintf("  ✅ Figure 5 saved (%dmm × %dmm) → %s",
                  fig5_w_mm, fig5_h_mm, FIG5_DIR))
} else {
  message("  ⚠ Figure 5 部分导出失败")
}


###############################################################################
# =============================================================================
# STAGE B：Supplementary Figure S4（7 sections，纯 cowplot 组装）
# =============================================================================
###############################################################################
message("\n\n#############################################################")
message("# STAGE B: Supplementary Figure S4 — comprehensive coverage")
message("#############################################################")

sections_config <- list(
  A = list(cohort_key="GSE157827",       disease="AD",  region="Middle frontal gyrus",
           excluded_celltypes = character(0)),
  B = list(cohort_key="GSE174367",       disease="AD",  region="Prefrontal cortex",
           excluded_celltypes = character(0)),
  C = list(cohort_key="syn52082747",     disease="AD",  region="Primary visual cortex (V1)",
           excluded_celltypes = c("Excitatory neurons")),
  D = list(cohort_key="syn21788402_EC",  disease="AD",  region="Entorhinal cortex",
           excluded_celltypes = character(0)),
  E = list(cohort_key="syn21788402_SFG", disease="AD",  region="Superior frontal gyrus",
           excluded_celltypes = character(0)),
  F = list(cohort_key="syn52082747",     disease="FTD", region="Primary visual cortex (V1)",
           excluded_celltypes = character(0)),
  G = list(cohort_key="syn52082747",     disease="PSP", region="Primary visual cortex (V1)",
           excluded_celltypes = c("Excitatory neurons","Inhibitory neurons"))
)

section_pages <- list()

for (section_letter in names(sections_config)) {
  scfg <- sections_config[[section_letter]]
  message(sprintf("\n========== Section %s: %s | %s | %s ==========",
                  section_letter, scfg$cohort_key, scfg$disease, scfg$region))
  
  df_all <- load_cohort_df(scfg$cohort_key)
  ct_to_plot <- setdiff(all_celltype_levels, scfg$excluded_celltypes)
  ct_to_plot <- intersect(ct_to_plot, levels(df_all$celltype))
  fdr_tbl <- compute_fdr_table(df_all, scfg$disease)
  
  hist_data_sec <- list(); dot_data_sec <- list()
  for (ct in ct_to_plot) {
    df_pair <- df_all %>%
      filter(group %in% c(scfg$disease, compare_ctrl),
             celltype == ct, expr > 0)
    if (nrow(df_pair) < 2) next
    hist_data_sec[[ct]] <- df_pair
    dot_data_sec[[ct]]  <- df_pair %>%
      group_by(donor, group) %>%
      summarise(median_expr = median(expr), .groups = "drop")
  }
  hist_y_max_sec <- compute_hist_y_max(hist_data_sec, binwidth_use)
  dot_y_lim_sec  <- compute_dot_y_lim(dot_data_sec)
  message(sprintf("  📐 Section Y 范围: hist [0, %.3f], dot [%.3f, %.3f]",
                  hist_y_max_sec, dot_y_lim_sec[1], dot_y_lim_sec[2]))
  
  hist_list <- list(); dot_list <- list()
  for (ct in ct_to_plot) {
    donor_med <- dot_data_sec[[ct]]
    if (is.null(donor_med)) next
    n_by_grp <- table(donor_med$group)
    has_both <- all(c(scfg$disease, compare_ctrl) %in% names(n_by_grp))
    if (!has_both ||
        any(n_by_grp[c(scfg$disease, compare_ctrl)] < MIN_DONORS_PER_GROUP)) {
      message("  ⤷ skip ", ct, " (某组 donor 数不足)"); next
    }
    fdr_val <- fdr_tbl$padj[as.character(fdr_tbl$celltype) == ct]
    if (length(fdr_val) == 0) fdr_val <- NA_real_
    h <- make_hist(hist_data_sec[[ct]], scfg$disease, ct, hist_y_max_sec)
    d <- make_dot (donor_med, scfg$disease, ct, fdr_val, dot_y_lim_sec)
    if (!is.null(h) && !is.null(d)) {
      hist_list[[ct]] <- h
      dot_list[[ct]]  <- d
    }
  }
  if (length(hist_list) == 0) {
    message("  ⚠ Section ", section_letter, " 没有数据可画"); next
  }
  
  # cohort_key 末尾的 _EC/_SFG 后缀不进标题展示
  cohort_display <- sub("_(EC|SFG)$", "", scfg$cohort_key)
  hdr_text <- sprintf("Section %s. %s — %s (%s)",
                      section_letter, cohort_display, scfg$disease, scfg$region)
  
  # ====== v2 关键：用 cowplot 组装 section page（不再用 patchwork & ）======
  section_page <- build_section_page_cowplot(hist_list, dot_list, hdr_text,
                                             font_family = base_fontfamily)
  section_pages[[section_letter]] <- section_page
  
  n_rows    <- length(hist_list)
  page_w_mm <- 170
  page_h_mm <- min(225, 22 + n_rows * 33)
  page_w_in <- page_w_mm / 25.4
  page_h_in <- page_h_mm / 25.4
  
  png_path <- file.path(SUPS4_DIR,
                        sprintf("SuppS4_Section%s_%s_%s_1000dpi.png",
                                section_letter, scfg$cohort_key, scfg$disease))
  pdf_path <- file.path(SUPS4_DIR,
                        sprintf("SuppS4_Section%s_%s_%s_vector.pdf",
                                section_letter, scfg$cohort_key, scfg$disease))
  
  status_sec <- save_png_pdf(section_page, png_path, pdf_path,
                             page_w_in, page_h_in, png_dpi)
  if (all(status_sec)) {
    message(sprintf("  ✅ Saved: %s (%dmm × %.0fmm)",
                    basename(pdf_path), page_w_mm, page_h_mm))
  } else {
    message("  ⚠ Section ", section_letter, " 部分导出失败")
  }
}

# ============ 合并多页 PDF + 1000dpi PNG ============
all_pdf <- file.path(SUPS4_DIR, "Supplementary_Figure_S4_MN_unified_ALL_vector.pdf")
safe_close_all_devices()
if (length(section_pages) > 0 && capabilities("cairo")) {
  tryCatch({
    cairo_pdf(all_pdf, width = 170/25.4, height = 225/25.4,
              onefile = TRUE, family = base_fontfamily)
    for (sl in names(section_pages)) print(section_pages[[sl]])
    dev.off()
    message("\n🎉 Supp S4 合并多页 PDF: ", all_pdf)
  }, error = function(e) {
    message("\n⚠ 合并 PDF 失败: ", conditionMessage(e))
    safe_close_all_devices()
  })
}

all_png <- file.path(SUPS4_DIR, "Supplementary_Figure_S4_MN_unified_ALL_1000dpi.png")
if (requireNamespace("magick", quietly = TRUE) && length(section_pages) > 0) {
  message("🔧 用 magick 拼接合并 1000dpi PNG ...")
  png_files_all <- list.files(SUPS4_DIR,
                              pattern = "^SuppS4_Section[A-G]_.*_1000dpi\\.png$", full.names = TRUE)
  letters_found <- sub("^SuppS4_Section([A-G]).*$", "\\1", basename(png_files_all))
  png_files_all <- png_files_all[order(letters_found)]
  if (length(png_files_all) > 0) {
    imgs <- magick::image_read(png_files_all)
    combined_img <- magick::image_append(imgs, stack = TRUE)
    info <- magick::image_info(combined_img)
    magick::image_write(combined_img, path = all_png,
                        format = "png", density = "1000x1000", quality = 100)
    size_mb <- file.size(all_png) / 1024 / 1024
    message(sprintf("🎉 Supp S4 合并总图 PNG: %s (%.2f MB)", all_png, size_mb))
    if (size_mb > 20) {
      message("   ⚠ 超过 MN 20 MB 上限，投稿用矢量 PDF")
    } else {
      message("   ✅ 符合 MN 20 MB 上限")
    }
  }
}

###############################################################################
# 最终总结
###############################################################################
message("\n\n=============================================================")
message("✅ 全部完成（v2，纯 cowplot 组装）。两张图保存在：")
message("  Figure 5    : ", FIG5_DIR)
message("    - Figure5_MN_unified_1000dpi.png")
message("    - Figure5_MN_unified_vector.pdf")
message("  Supp Fig S4 : ", SUPS4_DIR)
message("    - SuppS4_SectionA-G_{cohort}_{disease}_1000dpi.png / vector.pdf")
message("    - Supplementary_Figure_S4_MN_unified_ALL_vector.pdf  ← 投稿主用")
message("    - Supplementary_Figure_S4_MN_unified_ALL_1000dpi.png")
message("=============================================================\n")
