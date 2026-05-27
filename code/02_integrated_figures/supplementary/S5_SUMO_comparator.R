###############################################################################
# Supplementary Figure S5
# SUMO-family comparator histograms across the seven analytical units
# (cell-type-resolved).
#
#   Multi-page composite: facet_nested(unit + cell type ~ gene); columns =
#   SUMO1/SUMO2/SUMO3, rows = 7 analytical units x 6 cell types. Each panel is
#   the disease-vs-control overlapping density histogram of non-zero nuclei
#   (no statistical test; a distributional reference for the UBL3 analyses).
#
# Inputs  : 5 Seurat .rds objects (+ 2 syn21788402 metadata CSVs) at the
#           data-project paths in UNITS below (confirm they exist before running).
# Outputs : <REPO>/output/figures/S5/ (per-page vector PDF + 1000 dpi PNG,
#           a merged composite vector PDF, and a stitched composite PNG;
#           legend sidecar).
# MN specs: 170 mm wide, <=225 mm tall per page, vector cairo_pdf (fonts
#           embedded) + 1000 dpi PNG, line widths >0.25 pt, colourblind-safe,
#           supplementary file-size limit 20 MB.
###############################################################################

SEED <- 20251023; set.seed(SEED)

# --- Paths -------------------------------------------------------------------
# REPO : repository root; figure outputs go under here (same as Fig2-Fig4 / S3).
# Input Seurat objects and metadata CSVs are read from their data-project
# locations below; those are source data and are left unchanged.
REPO <- "D:/RNA/Code/UBL3_tauopathy"
Sys.setenv(LANG = "en"); options(stringsAsFactors = FALSE)

need_pkgs <- c("Seurat","SeuratObject","Matrix","dplyr","ggplot2","ragg",
               "grid","org.Hs.eg.db","AnnotationDbi","scales","ggh4x","pdftools")
to_install <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(dplyr)
  library(ggplot2); library(ragg); library(grid); library(scales)
  library(org.Hs.eg.db); library(AnnotationDbi); library(ggh4x)
})
if (.Platform$OS.type == "windows")
  tryCatch(windowsFonts(Arial = windowsFont("Arial")), error = function(e) NULL)
PLOT_FONT <- "Arial"

###############################################################################
# 0. 配置
###############################################################################
OUT_DIR <- file.path(REPO, "output", "figures", "S5")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIG_TAG <- "SuppFigS5_SUMO_byCelltype_MN_v10"

GENES        <- c("SUMO1", "SUMO2", "SUMO3")
BINWIDTH     <- 0.15
CTRL_NAME    <- "Control"
CTRL_ALIAS   <- c("NC","Control","CTRL","Ctr","CTR","Normal","N",
                  "control","ctrl","ctr","nc","normal")

# ── 排版常量 (可调) ─────────────────────────────────────────────────────────
FIG_W_MM        <- 170     # 全页宽 (MN 全栏)
UNITS_PER_PAGE  <- 2       # 每页放几个单元: 2 → 4 页(2+2+2+1, 行更高更不挤); 3 → 3 页(较挤)
ROW_MM          <- 15.0    # 每个细胞类型行的目标高度 (mm); 调大更松, 调小更紧
OVERHEAD_MM     <- 24      # 顶部 gene strip + 底部 x 轴/图例 的固定开销
PAGE_H_CAP_MM   <- 224     # 每页高度上限 (<225)
png_limit_mb    <- 20

# 每个单元配置 (与 v9 完全一致)
UNITS <- list(
  list(unit="GSE157827", disease="AD", region="MFG",
       obj="D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results/stepH_obj_celltype6_named.rds",
       group_mode="col",
       group_cands=c("group4","group2","group","Group","diagnosis","Dx","clinical_diagnosis"),
       donor_cands=c("autopsy_id","donor","Donor","sample","Sample","orig.ident","patient","subject")),
  list(unit="GSE174367", disease="AD", region="PFC",
       obj="D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/stepH_obj_celltype6_named.rds",
       group_mode="col",
       group_cands=c("group4","group2","group","Diagnosis","diagnosis","Group","Dx","clinical_diagnosis"),
       donor_cands=c("autopsy_id","donor","Donor","sample","Sample","orig.ident","patient","subject")),
  list(unit="syn21788402_EC", disease="AD", region="EC",
       obj="D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds",
       group_mode="csv",
       meta_fp="D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepP_syn21788402_matched_cells_meta.csv",
       donor_cands=c("PatientID","autopsy_id","donor","sample","SampleID","orig.ident")),
  list(unit="syn21788402_SFG", disease="AD", region="SFG",
       obj="D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepH_syn21788402_SFG_obj_celltype6.rds",
       group_mode="csv",
       meta_fp="D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepP_syn21788402_SFG_matched_cells_meta.csv",
       donor_cands=c("PatientID","autopsy_id","donor","sample","SampleID","orig.ident")),
  list(unit="syn52082747_AD",  disease="AD",  region="V1",
       obj="D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3/stepH_slim_uncompressed.rds",
       group_mode="col", group_cands=c("group4"), donor_cands=c("autopsy_id")),
  list(unit="syn52082747_PSP", disease="PSP", region="V1",
       obj="D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3/stepH_slim_uncompressed.rds",
       group_mode="col", group_cands=c("group4"), donor_cands=c("autopsy_id")),
  list(unit="syn52082747_FTD", disease="FTD", region="V1",   # 图内标 FTD (正文 PiD)
       obj="D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3/stepH_slim_uncompressed.rds",
       group_mode="col", group_cands=c("group4"), donor_cands=c("autopsy_id"))
)
REGION_FULL <- c(MFG="Middle frontal gyrus", PFC="Prefrontal cortex",
                 EC="Entorhinal cortex", SFG="Superior frontal gyrus",
                 V1="Primary visual cortex")

# ── celltype6 归一化 + 显示标签 (移植自原始 per-dataset 脚本) ────────────────
CT_LEVELS <- c("Astrocytes","Endothelial","Excitatory neurons",
               "Inhibitory neurons","Microglia","Oligodendrocytes")
normalize_celltype6 <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("Astro","Astrocyte","Astrocytes")]                 <- "Astrocytes"
  x[x %in% c("Endo","Endothelial","Endothelial cells")]         <- "Endothelial"
  x[x %in% c("Excit","Excitatory","Excitatory neurons")]        <- "Excitatory neurons"
  x[x %in% c("Inhib","Inhibitory","Inhibitory neurons")]        <- "Inhibitory neurons"
  x[x %in% c("Microgl","Micro","Microglia")]                    <- "Microglia"
  x[x %in% c("Oligo","Oligodendrocyte","Oligodendrocytes")]     <- "Oligodendrocytes"
  factor(x, levels = CT_LEVELS)
}
CT_LABEL <- c("Astrocytes"="Astrocytes","Endothelial"="Endothelial",
              "Excitatory neurons"="Excitatory\nneurons",
              "Inhibitory neurons"="Inhibitory\nneurons",
              "Microglia"="Microglia","Oligodendrocytes"="Oligodendrocytes")
CT_cands <- c("celltype6","celltype","cell_type","celltype6_named","celltype6_manual")

###############################################################################
# 1. helpers (与 v9 一致)
###############################################################################
get_counts_matrix_allcells <- function(obj, assay = "RNA") {
  DefaultAssay(obj) <- assay
  a <- obj[[assay]]; n_obj <- ncol(obj)
  m <- tryCatch({ a2 <- SeuratObject::JoinLayers(a)
                  SeuratObject::LayerData(a2, layer = "counts") }, error = function(e) NULL)
  if (is.null(m) || ncol(m) < n_obj * 0.9) {       # Seurat v5 分层未合并 → 手动 cbind
    lyrs  <- tryCatch(SeuratObject::Layers(a), error = function(e) character(0))
    clyrs <- grep("^counts", lyrs, value = TRUE)
    if (length(clyrs) >= 1) {
      mats <- lapply(clyrs, function(L)
        tryCatch(SeuratObject::LayerData(a, layer = L), error = function(e) NULL))
      mats <- mats[!vapply(mats, is.null, logical(1))]
      if (length(mats)) {
        m2 <- do.call(cbind, mats); m2 <- m2[, !duplicated(colnames(m2)), drop = FALSE]
        if (is.null(m) || ncol(m2) > ncol(m)) m <- m2
      }
    }
  }
  if (is.null(m)) m <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "counts"),
                                error = function(e) NULL)
  if (is.null(m) || length(dim(m)) != 2) stop("无法获取 counts 矩阵")
  cells <- intersect(colnames(obj), colnames(m)); m[, cells, drop = FALSE]
}
ENS_FALLBACK <- c(SUMO1="ENSG00000116030", SUMO2="ENSG00000188612",
                  SUMO3="ENSG00000184900", UBL3="ENSG00000122042")
locate_gene_row <- function(rn, gene_symbol) {
  if (gene_symbol %in% rn) return(gene_symbol)
  ci <- which(toupper(rn) == toupper(gene_symbol)); if (length(ci)) return(rn[ci[1]])
  rn_strip <- sub("\\.\\d+$", "", rn)
  ens <- tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys=gene_symbol,
            keytype="SYMBOL", columns="ENSEMBL")$ENSEMBL, error=function(e) NULL)
  ens <- unique(c(ens[!is.na(ens)], unname(ENS_FALLBACK[gene_symbol]))); ens <- ens[!is.na(ens)]
  if (length(ens)) {
    h1 <- intersect(ens, rn); if (length(h1)) return(h1[1])
    h2 <- which(rn_strip %in% ens); if (length(h2)) return(rn[h2[1]])
  }
  NA_character_
}
find_col <- function(md, cands) { hit <- intersect(cands, colnames(md)); if (length(hit)) hit[1] else NA_character_ }

# 返回每个细胞的 group2 / donor / celltype6
resolve_group_donor_ct <- function(obj, u) {
  md <- obj@meta.data
  dcol <- find_col(md, u$donor_cands); if (is.na(dcol)) stop(u$unit, ": 找不到 donor 列")
  ctcol <- find_col(md, CT_cands);     if (is.na(ctcol)) stop(u$unit, ": 找不到 celltype 列")
  donor <- trimws(as.character(md[[dcol]]))
  ct    <- normalize_celltype6(md[[ctcol]])
  if (u$group_mode == "col") {
    gcol <- find_col(md, u$group_cands); if (is.na(gcol)) stop(u$unit, ": 找不到 group 列")
    g <- trimws(as.character(md[[gcol]])); g[g %in% CTRL_ALIAS] <- CTRL_NAME
    grp <- ifelse(g == u$disease, "Disease", ifelse(g == CTRL_NAME, "Control", NA_character_))
  } else {
    if (!file.exists(u$meta_fp)) stop(u$unit, ": 缺 meta CSV ", u$meta_fp)
    mr <- read.csv(u$meta_fp, stringsAsFactors = FALSE)
    if (!all(c("sample","group") %in% colnames(mr))) stop(u$unit, ": CSV 缺 sample/group 列")
    jk <- find_col(md, c("SampleID","sample","orig.ident"))
    if (is.na(jk)) stop(u$unit, ": 对象缺 join key")
    obj_sample <- trimws(as.character(md[[jk]]))
    gr <- trimws(as.character(mr$group)); gr[gr %in% CTRL_ALIAS] <- CTRL_NAME
    map <- setNames(gr, trimws(as.character(mr$sample)))
    g <- unname(map[obj_sample])
    grp <- ifelse(g == u$disease, "Disease", ifelse(g == CTRL_NAME, "Control", NA_character_))
    if (all(is.na(donor))) donor <- obj_sample
  }
  data.frame(donor=donor, group2=grp, celltype6=ct,
             stringsAsFactors = FALSE, row.names = colnames(obj))
}

###############################################################################
# 2. 逐单元 × SUMO 计算 expr>0 长表 (含 celltype6); 无任何统计
###############################################################################
big_df <- list(); unit_lab <- character(0)   # unit_lab[u$unit] = 带 donor 数的左侧标题

for (u in UNITS) {
  cat("══", u$unit, "══\n")
  tryCatch({
    obj <- readRDS(u$obj); DefaultAssay(obj) <- "RNA"
    gd  <- resolve_group_donor_ct(obj, u)
    keep <- !is.na(gd$group2) & !is.na(gd$donor) & gd$donor != "" & !is.na(gd$celltype6)
    nD <- length(unique(gd$donor[keep & gd$group2=="Disease"]))
    nC <- length(unique(gd$donor[keep & gd$group2=="Control"]))
    cat(sprintf("   donors → %s(n=%d)  Control(n=%d)\n", u$disease, nD, nC))
    if (!any(keep)) stop("0 细胞匹配")

    # 左侧 unit 标题: 3 行 (datasetID | disease / 脑区(缩写) / donor 数)
    dataset <- sub("_.*$", "", u$unit)
    unit_lab[u$unit] <- sprintf("%s | %s\n%s (%s)\n%s n=%d \u00B7 Ctrl n=%d",
                                dataset, u$disease, REGION_FULL[[u$region]], u$region,
                                u$disease, nD, nC)

    mat <- get_counts_matrix_allcells(obj)
    if (ncol(mat) < ncol(obj) * 0.9)
      cat(sprintf("   \u26A0 矩阵仅 %d 细胞 < 对象 %d (分层可能未合并!)\n", ncol(mat), ncol(obj)))
    gd  <- gd[colnames(mat), , drop = FALSE]
    keep <- !is.na(gd$group2) & !is.na(gd$donor) & gd$donor != "" & !is.na(gd$celltype6)
    mat <- mat[, keep, drop = FALSE]; gd <- gd[keep, , drop = FALSE]
    lib <- Matrix::colSums(mat); rn <- rownames(mat)
    cat(sprintf("   matrix: %d genes × %d cells\n", nrow(mat), ncol(mat)))

    for (g in GENES) {
      gr <- locate_gene_row(rn, g)
      if (is.na(gr)) { cat("   \u26A0 找不到", g, "\n"); next }
      expr <- log1p((as.numeric(mat[gr, ]) / pmax(lib, 1)) * 1e4)
      d <- data.frame(expr=expr, group2=gd$group2, celltype6=gd$celltype6,
                      unit=u$unit, gene=factor(g, levels=GENES), stringsAsFactors=FALSE)
      d <- d[is.finite(d$expr) & d$expr > 0, , drop = FALSE]
      big_df[[paste(u$unit, g)]] <- d
    }
    rm(obj, mat, lib); gc()
  }, error = function(e) cat("   \u274C", u$unit, "失败:", conditionMessage(e), "\n"))
}
big_df <- do.call(rbind, big_df)

# unit 转 factor (按 UNITS 顺序, label 带 donor 数)
unit_order <- vapply(UNITS, function(u) u$unit, character(1))
unit_order <- unit_order[unit_order %in% names(unit_lab)]
lvl <- unname(unit_lab[unit_order])
big_df$unit      <- factor(unname(unit_lab[big_df$unit]), levels = lvl)
big_df$celltype6 <- factor(as.character(big_df$celltype6), levels = CT_LEVELS)

# 补齐空面板骨架: 保证每个 unit×celltype6×gene 组合都存在 (NA 行不画柱, 但保留面板,
#   否则某 unit 若缺某细胞类型, facet_nested 会少一行导致整列错位)
skel <- expand.grid(unit = levels(big_df$unit), celltype6 = CT_LEVELS,
                    gene = GENES, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
skel$expr <- NA_real_; skel$group2 <- "Disease"
skel$unit      <- factor(skel$unit,      levels = levels(big_df$unit))
skel$celltype6 <- factor(skel$celltype6, levels = CT_LEVELS)
skel$gene      <- factor(skel$gene,      levels = GENES)
big_df <- rbind(big_df, skel[, names(big_df)])

###############################################################################
# 3. 单页绘图函数 (一页 = 若干个 unit)
###############################################################################
# 安全的 y 刻度: 空面板(NA/Inf 范围)返回空, 正常面板给 ~3 个 pretty 断点
.bp3 <- scales::breaks_pretty(n = 3)
y_breaks_safe <- function(lims) {
  if (length(lims) != 2 || any(!is.finite(lims)) || lims[1] == lims[2]) return(numeric(0))
  .bp3(lims)
}
base_theme <- theme_bw(base_size = 7.2, base_family = PLOT_FONT) +
  theme(
    axis.title        = element_text(size = 8),
    axis.text         = element_text(size = 6, color = "black"),
    axis.text.y       = element_text(size = 5.6, color = "black"),
    axis.line         = element_line(linewidth = 0.3, color = "black"),
    axis.ticks        = element_line(linewidth = 0.3, color = "black"),
    strip.text.x      = element_text(face = "bold", size = 8.5),
    strip.text.y.left  = element_text(size = 6.4, angle = 0, lineheight = 1.0, hjust = 0),
    strip.background  = element_rect(fill = "grey92", colour = "grey35", linewidth = 0.3),
    strip.placement   = "outside",
    panel.grid.major  = element_line(colour = "grey90", linewidth = 0.22),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(colour = "grey35", fill = NA, linewidth = 0.3),
    panel.spacing.x   = unit(0.35, "lines"),
    panel.spacing.y   = unit(0.45, "lines"),
    legend.position   = "bottom", legend.title = element_blank(),
    legend.text       = element_text(size = 7.5),
    legend.key.size   = unit(3.4, "mm"),
    plot.margin       = margin(3, 5, 3, 3)
  )

make_page <- function(unit_labels) {
  dsub <- big_df[big_df$unit %in% unit_labels, , drop = FALSE]
  dsub$unit <- droplevels(factor(dsub$unit, levels = unit_labels))
  dctl <- dsub[dsub$group2 == "Control", ]; ddis <- dsub[dsub$group2 == "Disease", ]
  ggplot() +
    geom_histogram(data = dctl, aes(x=expr, y=after_stat(density), fill="Control", color="Control"),
                   binwidth=BINWIDTH, position="identity", alpha=0.92, linewidth=0.3, boundary=0, closed="left", na.rm=TRUE) +
    geom_histogram(data = ddis, aes(x=expr, y=after_stat(density), fill="Disease", color="Disease"),
                   binwidth=BINWIDTH, position="identity", alpha=0.78, linewidth=0.3, boundary=0, closed="left", na.rm=TRUE) +
    ggh4x::facet_nested(unit + celltype6 ~ gene, scales = "free_y", switch = "y",
                        labeller = labeller(celltype6 = CT_LABEL),
                        nest_line = element_line(linewidth = 0.3, colour = "grey35")) +
    scale_fill_manual(values = c(Disease="white", Control="grey65"),
                      breaks = c("Disease","Control"), labels = c("Disease","Control")) +
    scale_color_manual(values = c(Disease="black", Control="grey40"), guide = "none") +
    scale_y_continuous(breaks = y_breaks_safe,
                       expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = "SUMO  log1p(CP10k)", y = "Density") +
    guides(fill = guide_legend(nrow = 1, override.aes = list(
             alpha = c(1,1), color = c("black","grey40"),
             fill = c("white","grey65"), linewidth = 0.35))) +
    base_theme
}

###############################################################################
# 4. 分页 → 每页单独 PDF + 1000 dpi PNG → 合并成一个 composite PDF
###############################################################################
pages <- split(lvl, ceiling(seq_along(lvl) / UNITS_PER_PAGE))
W_in  <- FIG_W_MM / 25.4
page_pdfs <- character(0)

for (pi in seq_along(pages)) {
  ul <- pages[[pi]]
  n_rows  <- length(ul) * length(CT_LEVELS)                 # 该页细胞类型行数
  page_h  <- min(PAGE_H_CAP_MM, OVERHEAD_MM + n_rows * ROW_MM)
  H_in    <- page_h / 25.4
  p <- make_page(ul)

  pdf_p <- file.path(OUT_DIR, sprintf("%s_p%d.pdf", FIG_TAG, pi))
  png_p <- file.path(OUT_DIR, sprintf("%s_p%d_1000dpi.png", FIG_TAG, pi))
  cat(sprintf("\n💾 page %d: %d units, %d rows, %.0f mm 高\n", pi, length(ul), n_rows, page_h))

  grDevices::cairo_pdf(pdf_p, width = W_in, height = H_in, family = "Arial", bg = "white")
  print(p); dev.off()
  page_pdfs <- c(page_pdfs, pdf_p)

  ok <- tryCatch({ ragg::agg_png(png_p, width=W_in, height=H_in, units="in", res=1000, background="white")
                   print(p); dev.off(); TRUE }, error=function(e){ try(dev.off(),silent=TRUE); FALSE })
  mb <- if (file.exists(png_p)) file.info(png_p)$size/1024^2 else NA
  if (isTRUE(ok) && !is.na(mb) && mb > png_limit_mb) {
    cat(sprintf("   ⚠ PNG %.1f MB > %d, 回退 600 dpi\n", mb, png_limit_mb))
    ragg::agg_png(png_p, width=W_in, height=H_in, units="in", res=600, background="white")
    print(p); dev.off()
  }
}

# 合并所有页 PDF 为一个 composite (页高可不同, 完全合法)
out_pdf <- file.path(OUT_DIR, paste0(FIG_TAG, "_composite_vector.pdf"))
merged <- tryCatch({ pdftools::pdf_combine(input = page_pdfs, output = out_pdf); TRUE },
                   error = function(e) { cat("   ⚠ pdf_combine 失败:", conditionMessage(e), "\n"); FALSE })

###############################################################################
# 5. 标题 + legend 草稿 (无统计版本)
###############################################################################
title_txt <- "Cell-type-resolved SUMO1/2/3 expression distributions across seven analytical units"
legend_txt <- paste(
  "Cell-level expression distributions of SUMO1, SUMO2 and SUMO3 (columns) shown as overlapping density histograms of disease versus control nuclei, separated by major cell type (rows) within each of the seven analytical units.",
  "Only nuclei with non-zero expression are plotted; expression is log1p of counts-per-10,000. Disease nuclei are white with black outlines and control nuclei grey; the number of disease and control donors is given beside each unit. No statistical comparison is shown: these panels are a cell-type-resolved distributional reference for the donor-level UBL3 analyses in the main figures.",
  "Colours are colour-blind-safe and the key is shown within the figure. The syn52082747 frontotemporal-lobar group (Pick's disease in the main text) is labelled 'FTD' here.",
  sep = "\n\n")
writeLines(c(title_txt, "", legend_txt, "",
             sprintf("[words: title=%d, legend=%d]",
                     lengths(gregexpr("\\S+", title_txt)), lengths(gregexpr("\\S+", legend_txt)))),
           file.path(OUT_DIR, paste0(FIG_TAG, "_legend.txt")), useBytes = TRUE)

cat("\n🎉 完成\n   composite PDF:", out_pdf, if (isTRUE(merged)) "(已合并)" else "(合并失败, 用分页 PDF)",
    "\n   分页 PNG:", length(pages), "张  目录:", OUT_DIR, "\n")

## 把已生成的分页 1000dpi PNG 竖向拼成一张 composite PNG (无需重跑流程)
if (!requireNamespace("magick", quietly = TRUE)) install.packages("magick")
library(magick)

OUT_DIR <- file.path(REPO, "output", "figures", "S5")
FIG_TAG <- "SuppFigS5_SUMO_byCelltype_MN_v10"   # 你当前文件名前缀就是这个

# 按页号 p1,p2,... 顺序收集分页 PNG
pngs <- list.files(OUT_DIR, pattern = sprintf("^%s_p\\d+_1000dpi\\.png$", FIG_TAG), full.names = TRUE)
pngs <- pngs[order(as.integer(sub(".*_p(\\d+)_1000dpi\\.png$", "\\1", basename(pngs))))]
cat("拼接顺序:\n"); print(basename(pngs))

comp <- image_append(image_read(pngs), stack = TRUE)   # stack=TRUE = 竖向拼
out  <- file.path(OUT_DIR, paste0(FIG_TAG, "_composite_1000dpi.png"))
image_write(comp, out, format = "png")

mb <- file.info(out)$size / 1024^2
cat(sprintf("\u2705 写出 %s (%.1f MB)\n", out, mb))
if (mb > 20) cat("\u26A0 >20 MB, 补充图上限 20MB; 投稿就用 composite PDF, 或把上面 image_read 换成 image_read(pngs) |> image_resize(\"60%\")\n")

