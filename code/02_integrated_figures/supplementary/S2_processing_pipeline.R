###############################################################################
# Supplementary Figure S2
# Per-cohort processing-pipeline evidence for the two cohorts processed from
# raw count matrices in this study (GSE157827, GSE174367).
#
#   Panel A: per-nucleus QC violins before/after filtering (nFeature, nCount,
#            percent.mt) with applied thresholds and retained fraction.
#   Panel B: PCA elbow plots (PCs 1-50; first 20 retained), shared y-axis.
#   Panel C: unsupervised-cluster UMAPs at the chosen resolution.
#
# Inputs  : per-cohort QC / integration / clustering objects at the data-project
#           paths in section 1 below (confirm they exist before running).
# Outputs : <REPO>/output/figures/S2/ (composite vector PDF + 1000 dpi PNG,
#           plus a 600 dpi fallback if the PNG exceeds 20 MB; legend sidecar).
# MN specs: 170 mm wide, <=225 mm tall, vector cairo_pdf + 1000 dpi PNG,
#           line widths >0.25 pt, colourblind-safe palettes.
###############################################################################

SEED <- 20251023
set.seed(SEED)

# --- Paths -------------------------------------------------------------------
# REPO : repository root; figure outputs go under here (same as Fig2-Fig4 / S3).
# Input Seurat objects and metadata CSVs are read from their data-project
# locations below; those are source data and are left unchanged.
REPO <- "D:/RNA/Code/UBL3_tauopathy"

###############################################################################
# 0. 加载包
###############################################################################
need_pkgs <- c(
  "Seurat", "SeuratObject", "Matrix", "qs", "ggplot2", "patchwork",
  "scales", "ragg", "scattermore", "shadowtext", "png", "grid"
)

to_install <- need_pkgs[
  !vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(to_install)) {
  install.packages(to_install, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(qs)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ragg)
  library(scattermore)
  library(shadowtext)
  library(grid)
})

if (.Platform$OS.type == "windows") {
  tryCatch(
    windowsFonts(Arial = windowsFont("Arial")),
    error = function(e) NULL
  )
}

PLOT_FONT <- "Arial"

###############################################################################
# 1. 路径与 cohort 信息
###############################################################################
# A: QC 用 merged raw 对象
paths_qc <- list(
  GSE157827 = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results/GSE157827_merged_with_group_CLEAN.qs",
  GSE174367 = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/GSE174367_merged_raw_with_meta_v2_CLEAN.qs"
)

# B: Elbow 用 stepE 整合对象
paths_pca <- list(
  GSE157827 = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results/stepE_integrated_ALL_cca_v2.rds",
  GSE174367 = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/stepE_integrated_ALL_cca_GSE174367.rds"
)

# C: Cluster UMAP 用 stepF 聚类对象
paths_umap <- list(
  GSE157827 = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results/stepF_afterPCA_graph_umap_res1_v3style.rds",
  GSE174367 = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/stepF_afterPCA_graph_umap_res0.8_dims1_20_GSE174367.rds"
)

COHORTS <- c("GSE157827", "GSE174367")

meta <- list(
  GSE157827 = list(
    line1   = "GSE157827 | AD",
    line2   = "Middle frontal gyrus",
    n_pc    = 20,
    res     = "1",
    mt_mode = "ensg"
  ),
  GSE174367 = list(
    line1   = "GSE174367 | AD",
    line2   = "Prefrontal cortex",
    n_pc    = 20,
    res     = "0.8",
    mt_mode = "symbol"
  )
)

OUT_DIR <- file.path(REPO, "output", "figures", "S2")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FIG_TAG <- "SuppFigureS2_MN_v4_final"

###############################################################################
# 2. 配色与尺寸
###############################################################################
# 色盲安全配色
VIO_FILL   <- "#4E79A7"
VIO_STROKE <- "grey25"
THR_COLOR  <- "#E15759"
COL_USED   <- "#4E79A7"
COL_UNUSED <- "grey72"

# MN 尺寸
# 总高度 = 92 + 48 + 72 + 2*3 = 218 mm, 小于 225 mm
FIG_W_MM <- 170
H_A_MM   <- 92
H_B_MM   <- 48
H_C_MM   <- 72
GAP_MM   <- 3
FIG_H_MM <- H_A_MM + H_B_MM + H_C_MM + 2 * GAP_MM

PANELC_DPI   <- 600
png_limit_mb <- 20

# QC 阈值
THR_NFEATURE <- 200
THR_NCOUNT   <- 20000
THR_PCTMT    <- 20

# GSE157827 线粒体基因，ENSEMBL rownames
mt_genes_ensg <- c(
  "ENSG00000198888", "ENSG00000198727", "ENSG00000198804",
  "ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
  "ENSG00000198695", "ENSG00000198712", "ENSG00000198899",
  "ENSG00000198938", "ENSG00000198840", "ENSG00000198763",
  "ENSG00000210107", "ENSG00000210112", "ENSG00000210117",
  "ENSG00000210127", "ENSG00000210133", "ENSG00000210140",
  "ENSG00000210144", "ENSG00000210151", "ENSG00000210156",
  "ENSG00000210160", "ENSG00000210164", "ENSG00000210169",
  "ENSG00000210174", "ENSG00000210179", "ENSG00000210184",
  "ENSG00000210189", "ENSG00000210194", "ENSG00000210199",
  "ENSG00000210204", "ENSG00000210209", "ENSG00000210214",
  "ENSG00000210219", "ENSG00000228253", "ENSG00000228630",
  "ENSG00000210130"
)

###############################################################################
# 3. 通用 helpers
###############################################################################
load_obj <- function(fp) {
  if (!file.exists(fp)) {
    stop("❌ 缺失文件: ", fp)
  }
  
  if (grepl("\\.qs$", fp, ignore.case = TRUE)) {
    qs::qread(fp)
  } else {
    readRDS(fp)
  }
}

fmt_int <- function(x) {
  formatC(x, format = "d", big.mark = ",")
}

get_counts_v5 <- function(obj, assay = "RNA") {
  DefaultAssay(obj) <- assay
  
  m <- try(
    SeuratObject::GetAssayData(obj, assay = assay, layer = "counts"),
    silent = TRUE
  )
  
  if (!inherits(m, "try-error") && !is.null(m) && ncol(m) > 0) {
    return(m)
  }
  
  stop("无法取 counts layer")
}

rebuild_singlelayer_counts_if_needed <- function(obj) {
  a <- obj[["RNA"]]
  
  lyr <- try(SeuratObject::Layers(a), silent = TRUE)
  
  if (inherits(lyr, "try-error") || length(lyr) <= 1) {
    return(obj)
  }
  
  m1 <- try(SeuratObject::LayerData(a, layer = lyr[1]), silent = TRUE)
  
  if (inherits(m1, "try-error") || is.null(m1) || ncol(m1) == 0) {
    return(obj)
  }
  
  if (length(intersect(colnames(m1), colnames(obj))) > 0) {
    return(obj)
  }
  
  message("  ⚠ multi-layer counts → 重建单层对象")
  
  mats <- list()
  
  for (L in lyr) {
    mm <- try(SeuratObject::LayerData(a, layer = L), silent = TRUE)
    
    if (!inherits(mm, "try-error") && !is.null(mm) && ncol(mm) > 0) {
      mats[[L]] <- mm
    }
  }
  
  stopifnot(length(mats) > 0)
  
  g1 <- rownames(mats[[1]])
  stopifnot(all(vapply(mats, function(x) identical(rownames(x), g1), logical(1))))
  
  obj_new <- CreateSeuratObject(
    counts       = do.call(cbind, mats),
    project      = "QC_rebuilt",
    min.cells    = 0,
    min.features = 0
  )
  
  DefaultAssay(obj_new) <- "RNA"
  
  obj_new
}

###############################################################################
# 4. 共享主题
###############################################################################
theme_base_mn <- function() {
  theme_classic(base_size = 8, base_family = PLOT_FONT) +
    theme(
      plot.title = element_text(
        face   = "bold",
        hjust  = 0,
        size   = 8.5,
        margin = margin(b = 0)
      ),
      plot.subtitle = element_text(
        hjust  = 0,
        size   = 7.2,
        color  = "grey30",
        margin = margin(b = 2, t = 0)
      ),
      axis.title = element_text(
        size  = 7.2,
        color = "black"
      ),
      axis.text = element_text(
        color = "black",
        size  = 6.8
      ),
      axis.line = element_line(
        linewidth = 0.3,
        color     = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.3,
        color     = "black"
      ),
      axis.ticks.length = unit(1.8, "pt"),
      plot.margin = margin(3, 5, 1, 3)
    )
}

###############################################################################
# 5. Panel A — QC violins
###############################################################################
make_violin <- function(values, ylab, title_main, title_sub, thr = NULL, ylims = NULL) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  
  if (length(values) < 10) {
    return(
      ggplot() +
        theme_void() +
        labs(title = title_main, subtitle = title_sub)
    )
  }
  
  v <- values
  
  if (!is.null(ylims)) {
    v <- pmin(pmax(v, ylims[1]), ylims[2])
  }
  
  p <- ggplot(data.frame(v = v), aes(x = "", y = v)) +
    geom_violin(
      fill      = VIO_FILL,
      color     = VIO_STROKE,
      linewidth = 0.35,
      trim      = TRUE,
      alpha     = 0.85
    ) +
    geom_boxplot(
      width         = 0.10,
      fill          = "white",
      color         = VIO_STROKE,
      linewidth     = 0.3,
      outlier.shape = NA
    ) +
    stat_summary(
      fun    = median,
      geom   = "point",
      shape  = 21,
      size   = 1.2,
      color  = VIO_STROKE,
      fill   = "white",
      stroke = 0.3
    ) +
    labs(
      title    = title_main,
      subtitle = title_sub,
      x        = NULL,
      y        = ylab
    ) +
    scale_y_continuous(
      labels = scales::label_comma(accuracy = 1),
      limits = ylims,
      expand = expansion(mult = c(0.02, 0.08)),
      oob    = scales::squish
    ) +
    theme_base_mn() +
    coord_cartesian(clip = "on")
  
  if (!is.null(thr)) {
    p <- p +
      geom_hline(
        yintercept = thr,
        linetype   = "dashed",
        linewidth  = 0.35,
        color      = THR_COLOR
      )
  }
  
  p
}

# 修正版 QC strip:
# 三行显示，避免 cohort 标题、脑区和 retention 信息重叠
qc_strip <- function(line1, region, retained) {
  ggplot() +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    annotate(
      "text",
      x        = 0.02,
      y        = 0.92,
      label    = line1,
      hjust    = 0,
      vjust    = 1,
      fontface = "bold",
      size     = 3.0,
      family   = PLOT_FONT
    ) +
    annotate(
      "text",
      x      = 0.02,
      y      = 0.57,
      label  = region,
      hjust  = 0,
      vjust  = 1,
      size   = 2.55,
      color  = "grey30",
      family = PLOT_FONT
    ) +
    annotate(
      "text",
      x      = 0.02,
      y      = 0.25,
      label  = retained,
      hjust  = 0,
      vjust  = 1,
      size   = 2.45,
      color  = "grey30",
      family = PLOT_FONT
    ) +
    coord_cartesian(clip = "off") +
    theme(
      plot.margin = margin(0, 4, 4, 0)
    )
}

build_qc_block <- function(cohort) {
  m <- meta[[cohort]]
  
  cat("══ QC:", cohort, "══\n")
  
  obj <- load_obj(paths_qc[[cohort]])
  DefaultAssay(obj) <- "RNA"
  
  if (m$mt_mode == "ensg") {
    obj <- rebuild_singlelayer_counts_if_needed(obj)
  }
  
  mat <- get_counts_v5(obj)
  
  nFeature <- Matrix::colSums(mat > 0)
  nCount   <- Matrix::colSums(mat)
  
  is_mt <- if (m$mt_mode == "ensg") {
    rownames(mat) %in% mt_genes_ensg
  } else {
    grepl("^MT-", rownames(mat))
  }
  
  tot <- Matrix::colSums(mat)
  
  pct_mt <- ifelse(
    tot > 0,
    100 * Matrix::colSums(mat[is_mt, , drop = FALSE]) / tot,
    0
  )
  
  keep <- (nFeature > THR_NFEATURE) &
    (nCount < THR_NCOUNT) &
    (pct_mt < THR_PCTMT)
  
  n_b      <- length(nFeature)
  n_a      <- sum(keep)
  pct_kept <- 100 * n_a / n_b
  
  cat(sprintf(
    "  cells: %s → %s (%.1f%% kept)\n",
    fmt_int(n_b), fmt_int(n_a), pct_kept
  ))
  
  cap99 <- function(x, d) {
    x <- x[is.finite(x)]
    
    if (!length(x)) {
      return(d)
    }
    
    max(quantile(x, 0.995, na.rm = TRUE) * 1.15, d)
  }
  
  # 每个 cohort 内 before / after 使用同一 y-axis
  # 不同 cohort 之间允许按各自数据范围缩放，避免图形被压扁
  lim_nf <- c(0, cap99(nFeature, 1000))
  lim_nc <- c(0, cap99(nCount, 5000))
  lim_pm <- c(0, 25)
  
  bf <- function(metric, ylab, thr, ylims, phase, is_col1) {
    make_violin(
      values     = metric,
      ylab       = ylab,
      title_main = "",
      title_sub  = if (is_col1) phase else "",
      thr        = thr,
      ylims      = ylims
    )
  }
  
  row_before <- (
    bf(nFeature, "nFeature_RNA", THR_NFEATURE, lim_nf, "Before QC", TRUE) |
      bf(nCount, "nCount_RNA", THR_NCOUNT, lim_nc, "Before QC", FALSE) |
      bf(pct_mt, "percent.mt", THR_PCTMT, lim_pm, "Before QC", FALSE)
  )
  
  row_after <- (
    bf(nFeature[keep], "nFeature_RNA", THR_NFEATURE, lim_nf, "After QC", TRUE) |
      bf(nCount[keep], "nCount_RNA", THR_NCOUNT, lim_nc, "After QC", FALSE) |
      bf(pct_mt[keep], "percent.mt", THR_PCTMT, lim_pm, "After QC", FALSE)
  )
  
  strip <- qc_strip(
    line1    = m$line1,
    region   = m$line2,
    retained = sprintf(
      "retained %.1f%% (%s/%s nuclei)",
      pct_kept,
      fmt_int(n_a),
      fmt_int(n_b)
    )
  )
  
  rm(obj, mat, nFeature, nCount, pct_mt, tot, keep)
  gc()
  
  strip / row_before / row_after +
    plot_layout(heights = c(0.44, 1, 1.08))
}

###############################################################################
# 6. Panel B — PCA elbow plots
###############################################################################
ensure_pca50 <- function(obj) {
  stopifnot("integrated" %in% Assays(obj))
  
  if (!("pca" %in% Reductions(obj)) || length(obj[["pca"]]@stdev) < 50) {
    set.seed(SEED)
    obj <- ScaleData(
      object  = obj,
      assay   = "integrated",
      verbose = FALSE
    )
    
    set.seed(SEED)
    obj <- RunPCA(
      object  = obj,
      assay   = "integrated",
      npcs    = 50,
      verbose = FALSE
    )
  }
  
  obj
}

build_elbow <- function(cohort) {
  m <- meta[[cohort]]
  
  cat("══ Elbow:", cohort, "══\n")
  
  obj <- ensure_pca50(load_obj(paths_pca[[cohort]]))
  
  sdev <- obj[["pca"]]@stdev
  nd   <- min(50, length(sdev))
  
  df <- data.frame(
    PC    = 1:nd,
    stdev = sdev[1:nd]
  )
  
  df$status <- factor(
    ifelse(df$PC <= m$n_pc, "used", "unused"),
    levels = c("used", "unused")
  )
  
  rm(obj)
  gc()
  
  ggplot(df, aes(x = PC, y = stdev)) +
    geom_line(
      linewidth = 0.3,
      color     = "grey75",
      alpha     = 0.8
    ) +
    geom_point(
      aes(color = status, size = status, alpha = status),
      shape = 16
    ) +
    geom_vline(
      xintercept = m$n_pc + 0.5,
      linetype   = "dashed",
      linewidth  = 0.35,
      color      = THR_COLOR
    ) +
    scale_color_manual(
      values = c(used = COL_USED, unused = COL_UNUSED),
      guide  = "none"
    ) +
    scale_size_manual(
      values = c(used = 1.2, unused = 0.7),
      guide  = "none"
    ) +
    scale_alpha_manual(
      values = c(used = 1.0, unused = 0.65),
      guide  = "none"
    ) +
    scale_x_continuous(
      breaks = c(1, 10, 20, 30, 40, 50),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, 16),
      breaks = c(0, 5, 10, 15),
      labels = c("0", "5", "10", "15"),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    labs(
      title    = m$line1,
      subtitle = paste0(m$line2, " · first ", m$n_pc, " PCs used"),
      x        = "PC",
      y        = "Standard deviation"
    ) +
    theme_base_mn() +
    theme(
      panel.grid.major.y = element_line(
        linewidth = 0.2,
        color     = "grey92"
      ),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

###############################################################################
# 7. Panel C — Cluster UMAP
###############################################################################
build_palette_50 <- function() {
  core24 <- c(
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728",
    "#9467BD", "#8C564B", "#E377C2", "#7F7F7F",
    "#BCBD22", "#17BECF", "#AEC7E8", "#FFBB78",
    "#98DF8A", "#FF9896", "#C5B0D5", "#C49C94",
    "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5",
    "#393B79", "#E7CB94", "#843C39", "#8CA252"
  )
  
  c(
    core24,
    grDevices::hcl(
      h = seq(15, 345, length.out = 26),
      c = 80,
      l = 55
    )
  )
}

compute_centroids <- function(emb, clusters, min_dist = 0.90, iter = 1500, step = 0.09) {
  df <- data.frame(
    x   = emb[, 1],
    y   = emb[, 2],
    clu = as.character(clusters),
    stringsAsFactors = FALSE
  )
  
  cen <- aggregate(cbind(x, y) ~ clu, data = df, FUN = median)
  cen <- cen[order(suppressWarnings(as.numeric(cen$clu))), ]
  
  if (nrow(cen) <= 1) {
    return(cen)
  }
  
  xy <- as.matrix(cen[, c("x", "y")])
  n  <- nrow(xy)
  
  set.seed(SEED)
  xy <- xy + matrix(rnorm(2 * n, 0, 0.015), ncol = 2)
  
  k <- 0L
  
  while (k < iter) {
    k <- k + 1L
    
    dx <- outer(xy[, 1], xy[, 1], "-")
    dy <- outer(xy[, 2], xy[, 2], "-")
    
    d <- sqrt(dx^2 + dy^2)
    
    cm <- d < min_dist
    diag(cm) <- FALSE
    
    if (!any(cm)) {
      break
    }
    
    ds   <- pmax(d, 1e-3)
    push <- pmax(min_dist - d, 0) * step
    push[!cm] <- 0
    
    xy[, 1] <- xy[, 1] + rowSums((dx / ds) * push)
    xy[, 2] <- xy[, 2] + rowSums((dy / ds) * push)
  }
  
  cen$x <- xy[, 1]
  cen$y <- xy[, 2]
  
  cen
}

make_cluster_umap <- function(obj, title_main, title_sub, color_map) {
  emb <- Seurat::Embeddings(obj, "umap")
  clu <- as.character(Idents(obj))
  
  df <- data.frame(
    UMAP_1 = emb[, 1],
    UMAP_2 = emb[, 2],
    clu    = clu
  )
  
  df$color <- color_map[df$clu]
  
  cen <- compute_centroids(
    emb       = emb,
    clusters  = clu,
    min_dist  = 0.90,
    iter      = 1500,
    step      = 0.09
  )
  
  ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
    scattermore::geom_scattermore(
      aes(color = color),
      pointsize = 3.5,
      pixels    = c(2000, 2000),
      alpha     = 0.85
    ) +
    scale_color_identity() +
    shadowtext::geom_shadowtext(
      data        = cen,
      aes(x = x, y = y, label = clu),
      size        = 2.4,
      color       = "black",
      bg.color    = "white",
      bg.r        = 0.22,
      fontface    = "bold",
      family      = PLOT_FONT,
      inherit.aes = FALSE
    ) +
    labs(
      title    = title_main,
      subtitle = title_sub,
      x        = "UMAP_1",
      y        = "UMAP_2"
    ) +
    theme_base_mn() +
    theme(
      panel.grid      = element_blank(),
      legend.position = "none"
    )
}

###############################################################################
# 8. 构建 Panel A
###############################################################################
cat("\n########## 构建 Panel A: QC violins ##########\n")

panelA <- (
  build_qc_block("GSE157827") |
    build_qc_block("GSE174367")
) +
  plot_layout(widths = c(1, 1)) &
  theme(
    plot.margin = margin(5, 7, 4, 7)
  )

###############################################################################
# 9. 构建 Panel B
###############################################################################
cat("\n########## 构建 Panel B: PCA elbow plots ##########\n")

panelB <- (
  build_elbow("GSE157827") |
    build_elbow("GSE174367")
) +
  plot_layout(widths = c(1, 1)) &
  theme(
    plot.margin = margin(4, 7, 4, 7)
  )

###############################################################################
# 10. 构建 Panel C
###############################################################################
cat("\n########## 构建 Panel C: Cluster UMAPs ##########\n")

objsC <- lapply(COHORTS, function(c) {
  o <- load_obj(paths_umap[[c]])
  
  if (!"seurat_clusters" %in% colnames(o@meta.data)) {
    stop(c, " 无 seurat_clusters")
  }
  
  Idents(o) <- "seurat_clusters"
  
  o
})

names(objsC) <- COHORTS

all_levels <- unique(
  unlist(lapply(objsC, function(o) levels(Idents(o))))
)

all_levels <- all_levels[
  order(suppressWarnings(as.numeric(all_levels)))
]

stopifnot(length(all_levels) <= 50)

color_map <- setNames(
  build_palette_50()[seq_along(all_levels)],
  all_levels
)

cat("  跨 2 数据集 unique clusters:", length(all_levels), "\n")

umapC <- lapply(COHORTS, function(c) {
  o <- objsC[[c]]
  m <- meta[[c]]
  
  nclu <- length(levels(Idents(o)))
  
  make_cluster_umap(
    obj        = o,
    title_main = m$line1,
    title_sub  = sprintf(
      "%s (res = %s; %d clusters)",
      m$line2,
      m$res,
      nclu
    ),
    color_map  = color_map
  )
})

panelC_patch <- (
  umapC[[1]] |
    umapC[[2]]
) +
  plot_layout(widths = c(1, 1)) &
  theme(
    plot.margin = margin(4, 7, 4, 7)
  )

rm(objsC)
gc()

# Panel C 单独渲染成 PNG，再作为 rasterGrob 嵌入。
# 这样可以避免 scattermore 在 cowplot/patchwork 组合时偶发的 embedded nul bug。
c_png <- file.path(
  OUT_DIR,
  paste0(FIG_TAG, "_panelC_raster_source.png")
)

ragg::agg_png(
  filename   = c_png,
  width      = FIG_W_MM / 25.4,
  height     = H_C_MM / 25.4,
  units      = "in",
  res        = PANELC_DPI,
  background = "white"
)

print(panelC_patch)
dev.off()

gC <- grid::rasterGrob(
  png::readPNG(c_png),
  interpolate = TRUE
)

###############################################################################
# 11. 合成 A/B/C composite figure — grid 稳定导出版
###############################################################################
# 重要：
# 这里不再使用 wrap_elements(full = panelA) / wrap_elements(full = panelB) 的方式导出 PDF。
# 原因是复杂 patchwork + rasterGrob 在部分环境下可能导致 PDF 只显示局部 panel。
# 这里将 Panel A/B 转为 grob，然后在固定 170 × 218 mm 画布上手动排版。
# Panel letter A/B/C 也手动绘制，避免 patchwork 自动给内部小图加过多编号。

cat("\n########## 准备 composite grobs ##########\n")

gA <- patchwork::patchworkGrob(panelA)
gB <- patchwork::patchworkGrob(panelB)

draw_SF2_grid <- function() {
  grid::grid.newpage()
  
  main_layout <- grid::grid.layout(
    nrow = 5,
    ncol = 1,
    heights = grid::unit(
      c(H_A_MM, GAP_MM, H_B_MM, GAP_MM, H_C_MM),
      "mm"
    )
  )
  
  grid::pushViewport(
    grid::viewport(
      layout = main_layout,
      width  = grid::unit(FIG_W_MM, "mm"),
      height = grid::unit(FIG_H_MM, "mm")
    )
  )
  
  draw_panel_row <- function(grob_obj, row_id, tag_label) {
    grid::pushViewport(
      grid::viewport(
        layout.pos.row = row_id,
        layout.pos.col = 1,
        clip = "off"
      )
    )
    
    grid::grid.draw(grob_obj)
    
    grid::grid.text(
      label = tag_label,
      x = grid::unit(2.0, "mm"),
      y = grid::unit(1, "npc") - grid::unit(1.5, "mm"),
      just = c("left", "top"),
      gp = grid::gpar(
        fontfamily = PLOT_FONT,
        fontface   = "bold",
        fontsize   = 11,
        col        = "black"
      )
    )
    
    grid::upViewport()
  }
  
  draw_panel_row(gA, 1, "A")
  draw_panel_row(gB, 3, "B")
  draw_panel_row(gC, 5, "C")
  
  grid::upViewport()
}

###############################################################################
# 12. 导出 PDF 和 PNG — 稳定版
###############################################################################
W_in <- FIG_W_MM / 25.4
H_in <- FIG_H_MM / 25.4

out_pdf <- file.path(
  OUT_DIR,
  paste0(FIG_TAG, "_vector.pdf")
)

out_png_1000 <- file.path(
  OUT_DIR,
  paste0(FIG_TAG, "_1000dpi.png")
)

out_png_600 <- file.path(
  OUT_DIR,
  paste0(FIG_TAG, "_600dpi_under20MB.png")
)

cat("\n💾 Exporting vector PDF using grid device...\n")

pdf_ok <- tryCatch(
  {
    grDevices::cairo_pdf(
      filename = out_pdf,
      width    = W_in,
      height   = H_in,
      family   = "Arial",
      bg       = "white"
    )
    
    draw_SF2_grid()
    dev.off()
    
    TRUE
  },
  error = function(e) {
    try(dev.off(), silent = TRUE)
    message("PDF export failed: ", conditionMessage(e))
    FALSE
  }
)

cat("💾 Exporting 1000 dpi PNG using grid device...\n")

png_ok <- tryCatch(
  {
    ragg::agg_png(
      filename   = out_png_1000,
      width      = W_in,
      height     = H_in,
      units      = "in",
      res        = 1000,
      background = "white"
    )
    
    draw_SF2_grid()
    dev.off()
    
    TRUE
  },
  error = function(e) {
    try(dev.off(), silent = TRUE)
    message("1000 dpi PNG export failed: ", conditionMessage(e))
    FALSE
  }
)

png_1000_mb <- if (file.exists(out_png_1000)) {
  file.info(out_png_1000)$size / 1024^2
} else {
  NA_real_
}

out_png_final <- out_png_1000
png_final_mb  <- png_1000_mb

if (isTRUE(png_ok) && !is.na(png_1000_mb) && png_1000_mb > png_limit_mb) {
  cat(sprintf(
    "⚠ 1000 dpi PNG %.1f MB > %d MB. Exporting additional 600 dpi PNG for size control...\n",
    png_1000_mb,
    png_limit_mb
  ))
  
  ragg::agg_png(
    filename   = out_png_600,
    width      = W_in,
    height     = H_in,
    units      = "in",
    res        = 600,
    background = "white"
  )
  
  draw_SF2_grid()
  dev.off()
  
  out_png_final <- out_png_600
  png_final_mb  <- file.info(out_png_600)$size / 1024^2
}

###############################################################################
# 13. Figure title and legend draft
###############################################################################
title_txt <- "Supplementary Figure S2. Per-cohort processing-pipeline evidence from raw count matrices."

legend_txt <- paste(
  "Processing-pipeline diagnostics for the two cohorts reprocessed from raw count matrices in this study (GSE157827, middle frontal gyrus; GSE174367, prefrontal cortex), shown one cohort per column.",
  "(A) Per-nucleus quality control before and after filtering, shown as violin plots of the number of detected genes (nFeature_RNA), total UMI counts (nCount_RNA), and the percentage of mitochondrial reads (percent.mt). Dashed red lines mark the applied thresholds (nFeature > 200, nCount < 20,000, percent.mt < 20%), and the retained fraction is annotated. QC y-axes are scaled by cohort and metric.",
  "(B) PCA elbow plots showing the standard deviation explained by principal components 1-50. The first 20 components retained for downstream neighbour-graph construction are highlighted in blue, and the cut-off is marked by the dashed red line. Y-axes are shown on the same scale for the two cohorts.",
  "(C) UMAP embeddings coloured by unsupervised cluster at the resolution used for annotation, with cluster identities labelled in-plot.",
  "Colours are colour-blind-safe, and figure keys are shown within the panels.",
  sep = "\n\n"
)

title_n_words <- lengths(gregexpr("\\S+", title_txt))
legend_n_words <- lengths(gregexpr("\\S+", legend_txt))

legend_file <- file.path(
  OUT_DIR,
  paste0(FIG_TAG, "_legend.txt")
)

writeLines(
  c(
    title_txt,
    "",
    legend_txt,
    "",
    sprintf(
      "[words: title=%d, legend=%d]",
      title_n_words,
      legend_n_words
    )
  ),
  con = legend_file,
  useBytes = TRUE
)

###############################################################################
# 14. 完成提示
###############################################################################
cat("\n🎉 Supplementary Figure S2 final version completed\n")
cat("   PDF: ", out_pdf, "\n")
cat("   PNG 1000 dpi: ", out_png_1000, sprintf(" (%.2f MB)\n", png_1000_mb))

if (!identical(out_png_final, out_png_1000)) {
  cat("   PNG size-controlled final: ", out_png_final, sprintf(" (%.2f MB)\n", png_final_mb))
}

cat(sprintf(
  "   Size: %d × %d mm, Panel A=%d mm, Panel B=%d mm, Panel C=%d mm\n",
  FIG_W_MM,
  FIG_H_MM,
  H_A_MM,
  H_B_MM,
  H_C_MM
))

cat("   Legend: ", legend_file, "\n")

cat("\n检查重点:\n")
cat("   1) 打开 PDF，确认 A/B/C 三个 panel 都完整显示。\n")
cat("   2) Panel B 左右 Y 轴应统一为 0, 5, 10, 15。\n")
cat("   3) Panel A 顶部 cohort 信息不应重叠。\n")
cat("   4) Panel C UMAP 标签应保持可读。\n")

