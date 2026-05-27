###############################################################################
# Supplementary Figure S4
# Donor-level conditional UBL3 expression for all comparisons other than the
# three focal compartments in Figure 5.
#
#   Per analytical unit (7 sections, A-G): cell-level UBL3 density histograms
#   (Panel A) and donor-level median strip dotplots (Panel B) for every cell
#   type, with within-unit BH-FDR Mann-Whitney U. The three focal compartments
#   shown in Figure 5 (syn52082747 AD/PSP V1 neurons) are excluded here.
#
# NOTE: this figure was numbered 'S3' in earlier drafts; it is Supplementary
#       Figure S4 in the current manuscript. Output filenames use S4.
#
# Inputs  : 5 Seurat .rds objects (+ 2 syn21788402 metadata CSVs) at the
#           data-project paths in USER_CONFIG below (confirm before running).
# Outputs : <REPO>/output/figures/S4/ (per-section vector PDF + 1000 dpi PNG,
#           a merged multi-page vector PDF, and a stitched 1000 dpi PNG).
# MN specs: 170 mm wide, <=225 mm tall per page, vector cairo_pdf (fonts
#           embedded) + 1000 dpi PNG, line widths >0.25 pt, colourblind-safe,
#           supplementary file-size limit 20 MB.
###############################################################################

rm(list = ls()); gc()
SEED <- 20251023; set.seed(SEED)

# --- Paths -------------------------------------------------------------------
# REPO : repository root; figure outputs go under here (same as Fig2-Fig4 / S3).
# Input Seurat objects and metadata CSVs are read from their data-project
# locations below; those are source data and are left unchanged.
REPO <- "D:/RNA/Code/UBL3_tauopathy"
Sys.setenv(LANG = "en"); options(stringsAsFactors = FALSE)

# ★ v3 新增：magick 加入自动安装清单
for (pkg in c("lemon", "ggtext", "cowplot", "magick")) {
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
# =====================  USER CONFIG (已根据你脚本填好)  ========================
# 如有路径需要改，只改这块
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

# 7-section 定义
sections_config <- list(
  A = list(cohort_key="GSE157827",       disease="AD",  region="Middle frontal gyrus",
           excluded_celltypes = character(0)),
  B = list(cohort_key="GSE174367",       disease="AD",  region="Prefrontal cortex",
           excluded_celltypes = character(0)),
  C = list(cohort_key="syn52082747",     disease="AD",  region="Primary visual cortex (V1)",
           excluded_celltypes = c("Excitatory neurons")),  # 阳性→Fig 5
  D = list(cohort_key="syn21788402_EC",  disease="AD",  region="Entorhinal cortex",
           excluded_celltypes = character(0)),
  E = list(cohort_key="syn21788402_SFG", disease="AD",  region="Superior frontal gyrus",
           excluded_celltypes = character(0)),
  F = list(cohort_key="syn52082747",     disease="FTD", region="Primary visual cortex (V1)",
           excluded_celltypes = character(0)),
  G = list(cohort_key="syn52082747",     disease="PSP", region="Primary visual cortex (V1)",
           excluded_celltypes = c("Excitatory neurons","Inhibitory neurons"))  # 阳性→Fig 5
)

out_dir <- file.path(REPO, "output", "figures", "S4")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

###############################################################################
# 通用参数
###############################################################################
gene_symbol     <- "UBL3"
compare_ctrl    <- "Control"
binwidth_use    <- 0.15
png_dpi         <- 1000
base_fontfamily <- "Arial"
palette_groups  <- c("AD"="#E69F00","Control"="#56B4E9","FTD"="#009E73","PSP"="#CC79A7")
all_celltype_levels <- c("Astrocytes","Endothelial","Excitatory neurons",
                         "Inhibitory neurons","Microglia","Oligodendrocytes")
ctrl_alias <- c("Control","CTRL","Ctr","CTR","NC","Normal","N",
                "control","ctrl","ctr","nc","normal")

###############################################################################
# 工具函数（与 Fig 5 一致）
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

###############################################################################
# Cohort loader：处理 2 种数据来源模式
###############################################################################
load_cohort_df <- function(cohort_key) {
  cfg <- USER_CONFIG[[cohort_key]]
  stopifnot(file.exists(cfg$rds_path))
  
  message(sprintf("  ▼ 读对象: %s", cohort_key))
  obj <- readRDS(cfg$rds_path); DefaultAssay(obj) <- "RNA"
  md  <- obj@meta.data
  md$cell_id <- colnames(obj)
  rownames(md) <- md$cell_id
  
  # 列名自动识别
  ct_col  <- intersect(cfg$ct_col_pref,  colnames(md))[1]
  don_col <- intersect(cfg$don_col_pref, colnames(md))[1]
  stopifnot(!is.na(ct_col), !is.na(don_col))
  
  if (cfg$use_csv) {
    # ====== Pattern B: 外部 CSV 映射 sample→group ======
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
    # ====== Pattern A: obj@meta.data 直接有 group ======
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
  
  # 计算 UBL3 表达
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
# Themes (较 Fig 5 略紧凑因 supp 多面板)
###############################################################################
hist_theme <- theme_bw(base_size = 7.0, base_family = base_fontfamily) +
  theme(
    plot.title       = element_text(face = "bold", size = 7.2, hjust = 0,
                                    margin = margin(b = 2)),
    axis.title.x     = element_text(size = 6.6, margin = margin(t = 2)),
    axis.title.y     = element_text(size = 6.6, margin = margin(r = 2)),
    axis.text        = element_text(size = 6.0, colour = "black"),
    axis.ticks       = element_line(colour = "black", linewidth = 0.30),
    panel.border     = element_rect(colour = "grey35", fill = NA, linewidth = 0.30),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.20),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom", legend.direction = "horizontal",
    legend.title     = element_blank(), legend.text = element_text(size = 6.2),
    legend.key.width = unit(0.6,"lines"), legend.key.height = unit(0.6,"lines"),
    legend.margin    = margin(t = -3, b = 0),
    plot.margin      = margin(t = 3, r = 4, b = 1, l = 3)
  )

dot_theme <- theme_bw(base_size = 7.0, base_family = base_fontfamily) +
  theme(
    plot.title       = element_text(face = "bold", size = 7.2, hjust = 0,
                                    margin = margin(b = 2)),
    plot.subtitle    = element_text(face = "italic", size = 6.4, colour = "grey25",
                                    hjust = 0, margin = margin(b = 2)),
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(size = 6.6, margin = margin(r = 3)),
    axis.text.x      = element_text(size = 6.4, colour = "black"),
    axis.text.y      = element_text(size = 6.2, colour = "black"),
    axis.ticks       = element_line(colour = "black", linewidth = 0.30),
    panel.border     = element_rect(colour = "grey35", fill = NA, linewidth = 0.30),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.20),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom", legend.direction = "horizontal",
    legend.title     = element_blank(), legend.text = element_text(size = 6.2),
    legend.key.width = unit(0.6,"lines"), legend.key.height = unit(0.6,"lines"),
    legend.margin    = margin(t = -3, b = 0),
    plot.margin      = margin(t = 3, r = 4, b = 1, l = 3)
  )

###############################################################################
# Sub-plot 生成函数（注入 section 内统一 Y 范围）
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
               size = 1.6, shape = 21, stroke = 0.30,
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

###############################################################################
# 计算 section 内的统一 Y 范围
###############################################################################
compute_section_y_limits <- function(df_all, disease, ct_list, binwidth) {
  y_max_hist <- 0
  dot_vals <- numeric(0)
  for (ct in ct_list) {
    df_pair <- df_all %>%
      filter(group %in% c(disease, compare_ctrl),
             celltype == ct, expr > 0)
    if (nrow(df_pair) < 2) next
    for (grp in unique(df_pair$group)) {
      x_sub <- df_pair$expr[df_pair$group == grp]
      if (length(x_sub) < 2) next
      h <- hist(x_sub,
                breaks = seq(0, ceiling(max(x_sub)/binwidth)*binwidth + binwidth, by = binwidth),
                plot = FALSE)
      y_max_hist <- max(y_max_hist, h$density, na.rm = TRUE)
    }
    donor_med <- df_pair %>%
      group_by(donor, group) %>%
      summarise(median_expr = median(expr), .groups = "drop")
    dot_vals <- c(dot_vals, donor_med$median_expr)
  }
  
  if (y_max_hist <= 0) y_max_hist <- 1
  y_max_hist <- ceiling(y_max_hist * 10) / 10 * 1.05
  
  if (length(dot_vals) > 0) {
    dot_min <- min(dot_vals, na.rm = TRUE)
    dot_max <- max(dot_vals, na.rm = TRUE)
    dot_range <- dot_max - dot_min
    dot_lim <- c(dot_min - dot_range * 0.08, dot_max + dot_range * 0.12)
  } else {
    dot_lim <- c(0, 1)
  }
  
  list(hist_y_max = y_max_hist, dot_y_lim = dot_lim)
}

###############################################################################
# 主循环：逐 section 处理 + 保存
###############################################################################
section_pages <- list()
cohort_cache  <- list()  # 缓存 df_all 避免重复读 RDS

for (section_letter in names(sections_config)) {
  scfg <- sections_config[[section_letter]]
  
  message(sprintf("\n========== Section %s: %s | %s | %s ==========",
                  section_letter, scfg$cohort_key, scfg$disease, scfg$region))
  
  if (is.null(cohort_cache[[scfg$cohort_key]])) {
    cohort_cache[[scfg$cohort_key]] <- load_cohort_df(scfg$cohort_key)
  }
  df_all <- cohort_cache[[scfg$cohort_key]]
  
  ct_to_plot <- setdiff(all_celltype_levels, scfg$excluded_celltypes)
  ct_to_plot <- intersect(ct_to_plot, levels(df_all$celltype))
  
  fdr_tbl <- compute_fdr_table(df_all, scfg$disease)
  
  ylim <- compute_section_y_limits(df_all, scfg$disease, ct_to_plot, binwidth_use)
  message(sprintf("  📐 Section Y range: hist [0, %.3f], dot [%.3f, %.3f]",
                  ylim$hist_y_max, ylim$dot_y_lim[1], ylim$dot_y_lim[2]))
  
  hist_list <- list(); dot_list  <- list()
  for (ct in ct_to_plot) {
    df_pair <- df_all %>%
      filter(group %in% c(scfg$disease, compare_ctrl),
             celltype == ct, expr > 0)
    if (nrow(df_pair) == 0) next
    
    donor_med <- df_pair %>%
      group_by(donor, group) %>%
      summarise(median_expr = median(expr), .groups = "drop")
    
    fdr_val <- fdr_tbl$padj[as.character(fdr_tbl$celltype) == ct]
    if (length(fdr_val) == 0) fdr_val <- NA_real_
    
    h <- make_hist(df_pair, scfg$disease, ct, ylim$hist_y_max)
    d <- make_dot (donor_med, scfg$disease, ct, fdr_val, ylim$dot_y_lim)
    if (!is.null(h)) hist_list[[ct]] <- h
    if (!is.null(d)) dot_list[[ct]]  <- d
  }
  if (length(hist_list) == 0) {
    message("  ⚠ Section ", section_letter, " 没有数据可画")
    next
  }
  
  panel_A <- plot_grid(plotlist = hist_list, ncol = 1, align = "v")
  panel_B <- plot_grid(plotlist = dot_list,  ncol = 1, align = "v")
  
  hdr_text <- sprintf("Section %s. %s — %s (%s)",
                      section_letter, scfg$cohort_key, scfg$disease, scfg$region)
  header <- ggdraw() +
    draw_label(hdr_text, x = 0.02, hjust = 0,
               fontfamily = base_fontfamily, fontface = "bold", size = 8.5)
  
  section_body <- plot_grid(
    panel_A, panel_B, ncol = 2,
    labels = c("A","B"),
    label_size = 10,
    label_fontfamily = base_fontfamily,
    label_fontface = "bold",
    label_x = 0.01, label_y = 0.995,
    hjust = 0, vjust = 1,
    rel_widths = c(1, 1)
  )
  
  section_page <- plot_grid(header, section_body, ncol = 1,
                            rel_heights = c(0.035, 1))
  section_pages[[section_letter]] <- section_page
  
  n_rows <- length(hist_list)
  page_w_mm <- 170
  page_h_mm <- min(225, 22 + n_rows * 32)
  page_w_in <- page_w_mm / 25.4
  page_h_in <- page_h_mm / 25.4
  
  png_path <- file.path(out_dir,
    sprintf("SuppS4_Section%s_%s_%s_1000dpi.png",
            section_letter, scfg$cohort_key, scfg$disease))
  pdf_path <- file.path(out_dir,
    sprintf("SuppS4_Section%s_%s_%s_vector.pdf",
            section_letter, scfg$cohort_key, scfg$disease))
  
  ragg::agg_png(png_path, width = page_w_in, height = page_h_in,
                units = "in", res = png_dpi, background = "white")
  print(section_page); dev.off()
  
  if (capabilities("cairo")) {
    ggsave(pdf_path, section_page, width = page_w_in, height = page_h_in,
           units = "in", device = cairo_pdf, bg = "white")
  } else {
    ggsave(pdf_path, section_page, width = page_w_in, height = page_h_in,
           units = "in", device = "pdf", bg = "white")
  }
  message("  ✅ Saved: ", basename(pdf_path),
          sprintf(" (%dmm × %.0fmm)", page_w_mm, page_h_mm))
}

###############################################################################
# 合并多页 PDF（所有 sections 一文件，矢量）
###############################################################################
all_pdf <- file.path(out_dir, "SupplementaryFigureS4_MN_v3_ALL_vector.pdf")
if (length(section_pages) > 0 && capabilities("cairo")) {
  cairo_pdf(all_pdf, width = 170/25.4, height = 225/25.4,
            onefile = TRUE, family = base_fontfamily)
  for (sl in names(section_pages)) print(section_pages[[sl]])
  dev.off()
  message("\n🎉 多页合并 PDF: ", all_pdf)
}

###############################################################################
# ★ v3 新增：合并 1000dpi PNG（垂直拼接已生成的 section PNG）
###############################################################################
all_png <- file.path(out_dir, "SupplementaryFigureS4_MN_v3_ALL_1000dpi.png")

if (requireNamespace("magick", quietly = TRUE) && length(section_pages) > 0) {
  message("\n🔧 用 magick 拼接合并 1000dpi PNG ...")
  
  # 按 A→G 顺序取出已生成的 section PNG
  png_files_all <- list.files(out_dir,
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
    message(sprintf("🎉 合并总图 PNG: %s", all_png))
    message(sprintf("   尺寸: %d × %d px (@1000dpi → %.0f × %.0f mm)",
                    info$width, info$height,
                    info$width * 25.4 / 1000, info$height * 25.4 / 1000))
    message(sprintf("   文件大小: %.2f MB", size_mb))
    
    if (size_mb > 20) {
      message("   ⚠ 超过 MN supp 单文件 20 MB 上限")
      message("   → 投稿时建议改用矢量 PDF (极小, 完美无损)")
      message("   → 本 PNG 留作本地存档/快速浏览用")
    } else {
      message("   ✅ 文件大小符合 MN supp 20 MB 上限")
    }
  } else {
    message("⚠ 未找到 section PNG, 跳过合并 PNG 生成")
  }
} else {
  message("\n⚠ magick 包不可用, 跳过合并 PNG 生成")
}

###############################################################################
# 最终总结
###############################################################################
message("\n✅ 完成。所有 Supp S4 文件保存在：\n   ", out_dir)
message("   - SuppS4_SectionX_<cohort>_<disease>_1000dpi.png  (各 section, 1000dpi)")
message("   - SuppS4_SectionX_<cohort>_<disease>_vector.pdf   (各 section, 矢量 PDF)")
message("   - SupplementaryFigureS4_MN_v3_ALL_vector.pdf      (合并多页 PDF, 矢量)")
message("   - SupplementaryFigureS4_MN_v3_ALL_1000dpi.png     (合并总图 PNG, 1000dpi)")

