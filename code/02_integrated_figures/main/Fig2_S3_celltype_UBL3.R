###############################################################################
# Figure 2  +  Supplementary Figure S3  (MN submission, v6)
#
# Cell-type harmonization and UBL3 detectability
#   (a) DotPlot   (b) UMAP  — SAME a/b structure for BOTH figures
#
# 【这一版的核心改动】
#   - Panel (b) UMAP 不再读预渲染 PDF, 改为在脚本内从 Seurat 对象现算
#     (与 (a) DotPlot 共用同一次对象加载), 全部为 ggplot (scattermore+shadowtext).
#   - (b) 用 *一个* 共享 Cell type 图例 + *一个* 共享 UBL3 色条 (放在 panel b 底部),
#     与 (a) 右侧的共享 color+size 图例对称 —— 不再每个 cohort 各带一套图例.
#   - 每个 cohort 内不再有 A/B 子标签; 全图只有 2 个面板号: A(上, DotPlot 行) / B(下, UMAP 行).
#   - UMAP 绘图核心 (scattermore pointsize / shadowtext / Okabe-Ito 配色 /
#     YlOrRd UBL3 梯度 / UBL3_SHARED_MAX=4 / coord_cartesian 锁轴) 与
#     <dataset>UBL3umap.R (v15) 完全一致, 仅去掉各自图例并抽出共享图例.
#
# 【产出两张图 (各一个 composite)】
#   Figure 2 (2 个代表性阳性单元):
#       syn52082747_AD (AD, V1) , syn52082747_PSP (PSP, V1)
#       -> D:/RNA/MNversion submission/Figure 2. Cell-type harmonization and
#          UBL3 detectability across seven analytical units/Results
#   Supplementary Figure S3 (其余 5 个单元):
#       GSE157827(AD,MFG), GSE174367(AD,PFC), syn21788402_EC(AD,EC),
#       syn21788402_SFG(AD,SFG), syn52082747_FTD(FTD,V1)
#       -> D:/RNA/MNversion submission/Supplementary Figure S3. negative UMAP/Results
#
# 【MN 合规】170 mm 宽; <=225 mm 高; 矢量 cairo_pdf (字体内嵌) + 1000 dpi PNG;
#   线宽 >=0.3 pt; 色盲安全 (Okabe-Ito + YlOrRd + 蓝白橙发散); 图内含 key;
#   主图 <=10 MB, 补充图 <=20 MB; figure title <=15 words; legend <=300 words(写入 sidecar).
#
# 【关键可调旋钮 (出图后微调用)】见各处标注 "<<< TUNE"
###############################################################################

rm(list = ls()); gc()
SEED <- 20251023; set.seed(SEED)
Sys.setenv(LANG = "en"); options(stringsAsFactors = FALSE)

# =============================================================================
# 0. 包 + 字体
# =============================================================================
need_pkgs <- c("Seurat","SeuratObject","Matrix","dplyr","ggplot2","patchwork",
               "scales","cowplot","ragg","fs","scattermore","shadowtext","magick","png")
to_install <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(dplyr)
  library(ggplot2); library(patchwork); library(scales); library(grid)
  library(cowplot); library(ragg); library(fs); library(magick)
  library(scattermore); library(shadowtext)
})

PLOT_FONT <- "Arial"
if (.Platform$OS.type == "windows") {
  ok <- tryCatch({ windowsFonts(Arial = windowsFont("TT Arial")); TRUE },
                 error = function(e) FALSE, warning = function(w) TRUE)
  if (!ok) { PLOT_FONT <- "sans"; cat("\u26A0 Arial \u6CE8\u518C\u5931\u8D25 -> sans\n") } else cat("\u2705 Arial OK\n")
} else PLOT_FONT <- "sans"

# =============================================================================
# 1. 常量: 字号 / 线宽 / 配色 / UMAP 绘图参数
# =============================================================================
NBD_BASE_SIZE  <- 7      # 主字号 (title/y-axis/legend)
NBD_XAXIS_SIZE <- 6      # dotplot x 轴 marker 名
NBD_LINE_WIDTH <- 0.3    # 线宽 (>0.25 pt)
OUTER_TAG_SIZE <- 12     # A / B 面板号
COHORT_TITLE_SIZE <- 7   # cohort 小标题
SUBPLOT_TITLE_SIZE <- 6.5 # "Cell type" / "UBL3" 子标题

# --- UMAP scattermore 参数 (与 <dataset>UBL3umap.R v15 一致) ---
PIXELS_GRID  <- c(2400, 2400)
PT_CELL_MAIN <- 5.2      # <<< TUNE: cohort 变小后若点太大可降到 3.5~4
PT_UBL3_MAIN <- 3.6
PT_UBL3_BG   <- 3.8
BG_GREY      <- "grey92"   # 与 v15 单图脚本一致 (背景细胞云浅灰)

# --- Okabe-Ito cell-type 配色 (UMAP, 与 v15 一致) ---
celltype6_levels_std <- c("Astrocytes","Excitatory neurons","Inhibitory neurons",
                          "Microglia","Endothelial","Oligodendrocytes")
celltype6_palette_NBD <- c(
  "Astrocytes" = "#E69F00", "Excitatory neurons" = "#009E73",
  "Inhibitory neurons" = "#0072B2", "Microglia" = "#56B4E9",
  "Endothelial" = "#D55E00", "Oligodendrocytes" = "#CC79A7"
)

# --- UBL3 共享色条 (与 v15 一致) ---
UBL3_SHARED_MAX <- 4.0   # <<< 若任一 cohort UBL3 max 超 4, 改成 5 并重跑
ubl3_breaks  <- seq(0, UBL3_SHARED_MAX, by = 1)
ubl3_labels  <- as.character(ubl3_breaks)
ubl3_palette <- c("#FFFFCC","#FFEDA0","#FED976","#FEB24C","#FD8D3C",
                  "#FC4E2A","#E31A1C","#BD0026","#800026")

# --- dotplot y 轴 cell-type 缩写 ---
celltype_order  <- c("Astrocytes","Endothelial","Excitatory neurons",
                     "Inhibitory neurons","Microglia","Oligodendrocytes")
celltype_abbrev <- c("Astrocytes"="Astro","Endothelial"="Endo",
                     "Excitatory neurons"="Excit","Inhibitory neurons"="Inhib",
                     "Microglia"="Micro","Oligodendrocytes"="Oligo")

# UBL3 canonical Ensembl gene ID (GRCh38, chr13) = ENSG00000122042
#   之前误写成 ENSG00000146521 (别的基因!), 导致 ENSEMBL-keyed 的 GSE157827 取错行 → UBL3 全灰。
#   优先现场查 org.Hs.eg.db (与 GSE155114 / v15 模板一致); 不可用则用下方正确常量。
UBL3_ENSG <- "ENSG00000122042"
if (requireNamespace("org.Hs.eg.db", quietly=TRUE) && requireNamespace("AnnotationDbi", quietly=TRUE)) {
  .e <- tryCatch(
    AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys="UBL3",
                          keytype="SYMBOL", columns="ENSEMBL")$ENSEMBL[1],
    error = function(e) NA_character_)
  if (!is.na(.e) && nzchar(.e)) UBL3_ENSG <- .e
}
message("✅ UBL3 Ensembl ID 使用: ", UBL3_ENSG)

# =============================================================================
# 2. Markers (18) + dataset configs (7 units)
# =============================================================================
MARKERS <- list(
  list(symbol="SLC17A7", ensg="ENSG00000104888"), list(symbol="RBFOX3", ensg="ENSG00000167281"),
  list(symbol="SATB2",   ensg="ENSG00000119042"), list(symbol="GAD1",   ensg="ENSG00000128683"),
  list(symbol="GAD2",    ensg="ENSG00000136750"), list(symbol="SLC32A1",ensg="ENSG00000101438"),
  list(symbol="AQP4",    ensg="ENSG00000171885"), list(symbol="GFAP",   ensg="ENSG00000131095"),
  list(symbol="SLC1A2",  ensg="ENSG00000110436"), list(symbol="MBP",    ensg="ENSG00000197971"),
  list(symbol="MOG",     ensg="ENSG00000204655"), list(symbol="PLP1",   ensg="ENSG00000123560"),
  list(symbol="P2RY12",  ensg="ENSG00000169313"), list(symbol="CX3CR1", ensg="ENSG00000168329"),
  list(symbol="C1QA",    ensg="ENSG00000173372"), list(symbol="CLDN5",  ensg="ENSG00000184113"),
  list(symbol="FLT1",    ensg="ENSG00000102755"), list(symbol="VWF",    ensg="ENSG00000110799")
)

ec_path    <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds"
sfg_path   <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepH_syn21788402_SFG_obj_celltype6.rds"
syn52_path <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3/stepH_slim_uncompressed.rds"

resolver_braak <- function(md, lab="AD") {
  bs <- trimws(as.character(md$BraakStage)); out <- rep(NA_character_, length(bs))
  out[bs %in% c("0","Braak 0","BraakStage_0","Braak0")] <- "Control"
  out[bs %in% c("6","Braak 6","BraakStage_6","Braak6")] <- lab; out
}
resolver_simple <- function(md, col, dz, lab) {
  g <- trimws(as.character(md[[col]]))
  ctrl <- c("Control","CTRL","Ctr","CTR","NC","Normal","control","ctrl","ctr","nc","normal")
  out <- rep(NA_character_, length(g)); out[g %in% ctrl] <- "Control"; out[g %in% dz] <- lab; out
}

# umap_group4: UMAP 取细胞方式 (照搬各自 v15 单图脚本)
#   syn52082747 三个 -> 按 group4 == 疾病 subset; 其余四个 -> NA (用全部细胞, 不 subset)
CFG_ALL <- list(
  GSE157827      = list(unit="GSE157827", path="D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results/stepH_obj_celltype6_named.rds",
                        disease="AD", region="Middle frontal gyrus",
                        resolver=function(md) resolver_simple(md,"group",c("AD"),"AD"), ct_col="celltype6", umap_group4=NA),
  GSE174367      = list(unit="GSE174367", path="D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/stepH_obj_celltype6_named.rds",
                        disease="AD", region="Prefrontal cortex",
                        resolver=function(md) resolver_simple(md,"Diagnosis",c("AD"),"AD"), ct_col="celltype6", umap_group4=NA),
  syn52082747_AD = list(unit="syn52082747_AD", path=syn52_path, disease="AD", region="Primary visual cortex (V1)",
                        resolver=function(md) resolver_simple(md,"group4",c("AD"),"AD"), ct_col="celltype6", umap_group4="AD"),
  syn21788402_EC = list(unit="syn21788402_EC", path=ec_path, disease="AD", region="Entorhinal cortex",
                        resolver=function(md) resolver_braak(md,"AD"), ct_col="celltype6", umap_group4=NA),
  syn21788402_SFG= list(unit="syn21788402_SFG", path=sfg_path, disease="AD", region="Superior frontal gyrus",
                        resolver=function(md) resolver_braak(md,"AD"), ct_col="celltype6", umap_group4=NA),
  syn52082747_FTD= list(unit="syn52082747_FTD", path=syn52_path, disease="FTD", region="Primary visual cortex (V1)",
                        resolver=function(md) resolver_simple(md,"group4",c("FTD"),"FTD"), ct_col="celltype6", umap_group4="FTD"),
  syn52082747_PSP= list(unit="syn52082747_PSP", path=syn52_path, disease="PSP", region="Primary visual cortex (V1)",
                        resolver=function(md) resolver_simple(md,"group4",c("PSP"),"PSP"), ct_col="celltype6", umap_group4="PSP")
)

UNITS_FIG2 <- c("syn52082747_AD","syn52082747_PSP")
UNITS_SF3  <- c("GSE157827","GSE174367","syn21788402_EC","syn21788402_SFG","syn52082747_FTD")

# =============================================================================
# 3. Helpers (load / cell-type 标准化 / marker 定位 / counts)
# =============================================================================
ct_normalize <- function(x) {
  x <- as.character(x)
  x[x %in% c("Excit","Excitatory","Excitatory neurons")] <- "Excitatory neurons"
  x[x %in% c("Inhib","Inhibitory","Inhibitory neurons")] <- "Inhibitory neurons"
  x[x %in% c("Astro","Astrocytes")] <- "Astrocytes"
  x[x %in% c("Oligo","Oligodendrocytes","Oligodendrocyte")] <- "Oligodendrocytes"
  x[x %in% c("Microgl","Microglia")] <- "Microglia"
  x[x %in% c("Endo","Endothelial")] <- "Endothelial"; x
}
standardize_celltype6 <- function(x) {
  x <- as.character(x)
  x[x %in% c("Astrocytes","Astro","Astrocyte")] <- "Astrocytes"
  x[x %in% c("Excitatory neurons","Excitatory","Ex_neuron","Excit")] <- "Excitatory neurons"
  x[x %in% c("Inhibitory neurons","Inhibitory","Inh_neuron","Inhib")] <- "Inhibitory neurons"
  x[x %in% c("Microglia","Micro","Microgl")] <- "Microglia"
  x[x %in% c("Endothelial","Endothelial cells","Endo")] <- "Endothelial"
  x[x %in% c("Oligodendrocytes","Oligo","Oligodendrocyte","Oligodendro")] <- "Oligodendrocytes"
  factor(x, levels = celltype6_levels_std)
}
locate_marker <- function(rn, symbol, ensg) {
  if (symbol %in% rn) return(symbol)
  hit <- which(toupper(rn) == toupper(symbol)); if (length(hit)) return(rn[hit[1]])
  if (ensg %in% rn) return(ensg)
  h <- which(sub("\\.\\d+$","",rn) == ensg); if (length(h)) return(rn[h[1]])
  NA_character_
}
get_counts_robust <- function(obj, assay = "RNA") {
  a <- obj[[assay]]; ncells <- ncol(obj)
  layers_all <- tryCatch(SeuratObject::Layers(a), error=function(e) character(0))
  cl <- layers_all[grepl("^counts", layers_all)]
  if (length(cl) > 1) {
    mats <- list()
    for (ly in cl) { mm <- tryCatch(SeuratObject::LayerData(a, layer=ly), error=function(e) NULL)
                     if (!is.null(mm) && ncol(mm) > 0) mats[[ly]] <- mm }
    if (length(mats)) {
      allg <- unique(unlist(lapply(mats, rownames))); ali <- list()
      for (nm in names(mats)) { m0 <- mats[[nm]]
        if (identical(rownames(m0), allg)) ali[[nm]] <- m0 else {
          mn <- Matrix::Matrix(0, nrow=length(allg), ncol=ncol(m0), sparse=TRUE)
          rownames(mn) <- allg; colnames(mn) <- colnames(m0)
          cg <- intersect(allg, rownames(m0)); mn[cg,] <- m0[cg,,drop=FALSE]; ali[[nm]] <- mn } }
      mat <- Reduce(Matrix::cbind2, ali); mat <- mat[, !duplicated(colnames(mat)), drop=FALSE]
      miss <- setdiff(colnames(obj), colnames(mat))
      if (length(miss)) { mf <- Matrix::Matrix(0, nrow=nrow(mat), ncol=length(miss), sparse=TRUE)
        rownames(mf) <- rownames(mat); colnames(mf) <- miss; mat <- Matrix::cbind2(mat, mf) }
      return(mat[, colnames(obj), drop=FALSE])
    }
  }
  m <- tryCatch(Seurat::GetAssayData(obj, assay=assay, slot="counts"), error=function(e) NULL)
  if (!is.null(m) && nrow(m) && ncol(m)==ncells) return(m)
  m <- tryCatch(SeuratObject::LayerData(a, layer="counts"), error=function(e) NULL)
  if (!is.null(m) && nrow(m) && ncol(m)==ncells) return(m)
  m <- tryCatch(Seurat::GetAssayData(obj, assay=assay, slot="data"), error=function(e) NULL)
  if (!is.null(m) && ncol(m)==ncells) return(m)
  stop("get_counts_robust failed")
}
display_id <- function(u) sub("_(EC|SFG|AD|FTD|PSP)$","",u)

# =============================================================================
# 4. process_unit: 一次加载 -> 同时产出 dotplot 数据 + UMAP 数据
# =============================================================================
process_unit <- function(cfg) {
  cat("\n========== ", cfg$unit, " (", cfg$disease, ") ==========\n", sep="")
  stopifnot(file.exists(cfg$path))
  obj <- readRDS(cfg$path); DefaultAssay(obj) <- "RNA"
  cat("  total cells =", ncol(obj), "\n")
  md <- obj@meta.data
  group_std <- cfg$resolver(md)
  ct_raw    <- as.character(md[[cfg$ct_col]])

  # ---- counts (once) ----
  counts <- get_counts_robust(obj, "RNA")
  rn  <- rownames(counts)
  lib <- Matrix::colSums(counts); lib_safe <- pmax(lib, 1)

  # ---- (a) dotplot df: disease + control 细胞 ----
  keep_dot <- !is.na(group_std) & !is.na(ct_raw) & ct_raw != ""
  mk_rows <- vapply(MARKERS, function(m) locate_marker(rn, m$symbol, m$ensg), character(1))
  mk_syms <- vapply(MARKERS, function(m) m$symbol, character(1))
  found <- !is.na(mk_rows); mk_rows <- mk_rows[found]; mk_syms <- mk_syms[found]
  cat("  markers found:", sum(found), "/", length(MARKERS), "\n")
  dm <- log1p(Matrix::t(Matrix::t(counts[mk_rows,,drop=FALSE]) / lib_safe) * 1e4)
  dm <- as.matrix(dm[, keep_dot, drop=FALSE]); rownames(dm) <- mk_syms
  ctf <- factor(ct_normalize(ct_raw[keep_dot]), levels = celltype_order)
  res <- list()
  for (lvl in levels(ctf)) {
    idx <- which(ctf == lvl & !is.na(ctf)); if (!length(idx)) next
    sm <- dm[, idx, drop=FALSE]
    for (g in mk_syms) { e <- as.numeric(sm[g,])
      res[[length(res)+1]] <- data.frame(gene=g, celltype=as.character(lvl),
        avg_exp=mean(e,na.rm=TRUE), pct_exp=mean(e>0,na.rm=TRUE)*100, stringsAsFactors=FALSE) }
  }
  dot_df <- do.call(rbind, res)
  dot_df$avg_exp_scaled <- NA_real_
  for (g in unique(dot_df$gene)) { ix <- which(dot_df$gene==g); v <- dot_df$avg_exp[ix]
    dot_df$avg_exp_scaled[ix] <- if (length(v)>=2 && sd(v,na.rm=TRUE)>0) (v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE) else 0 }
  dot_df$avg_exp_scaled <- pmin(pmax(dot_df$avg_exp_scaled,-2.5),2.5)
  dot_df$gene     <- factor(dot_df$gene, levels = mk_syms)
  dot_df$celltype <- factor(dot_df$celltype, levels = rev(celltype_order))

  # ---- (b) UMAP df: disease-only 细胞 ----
  #   UBL3 完全照搬 v15 单图脚本: 逐 counts layer 求 lib_size 与 UBL3 行, gene_row 先 ENSEMBL 后 symbol
  rna <- obj[["RNA"]]
  ubl3_layers <- grep("^counts", SeuratObject::Layers(rna), value=TRUE)
  if (length(ubl3_layers)==0) ubl3_layers <- "counts"
  m0 <- SeuratObject::LayerData(rna, layer=ubl3_layers[1])
  gene_row <- if (UBL3_ENSG %in% rownames(m0)) UBL3_ENSG else
              if ("UBL3" %in% rownames(m0)) "UBL3" else locate_marker(rownames(m0), "UBL3", UBL3_ENSG)
  if (is.na(gene_row)) stop("UBL3 not found for ", cfg$unit)
  cells_all <- colnames(obj)
  lib_v15  <- setNames(rep(0, length(cells_all)), cells_all)
  ubl3_cnt <- setNames(rep(0, length(cells_all)), cells_all)
  for (ly in ubl3_layers) {
    m <- SeuratObject::LayerData(rna, layer=ly); cn <- colnames(m)
    if (length(cn)==0) next
    lib_v15[cn]  <- Matrix::colSums(m)
    ubl3_cnt[cn] <- as.numeric(m[gene_row, ])
  }
  ubl3_log1p <- log1p((ubl3_cnt / pmax(lib_v15, 1)) * 1e4)
  emb <- Embeddings(obj, "umap")
  # cell type 列: v15 用 celltype6 / celltype7 / celltype 依次匹配
  ct_col_use <- intersect(c("celltype6","celltype7","celltype"), colnames(obj@meta.data))[1]
  if (is.na(ct_col_use)) ct_col_use <- cfg$ct_col
  ct6 <- standardize_celltype6(obj@meta.data[[ct_col_use]])
  # UMAP 取细胞: syn52082747 按 group4 subset; 其余四个用全部细胞 (各自 v15 一致)
  if (!is.null(cfg$umap_group4) && !is.na(cfg$umap_group4)) {
    g4  <- as.character(obj$group4)
    sel <- !is.na(g4) & g4 == cfg$umap_group4 & !is.na(ct6)
  } else {
    sel <- !is.na(ct6)
  }
  cells_dz <- colnames(obj)[sel]
  df_all <- data.frame(UMAP_1=emb[cells_dz,1], UMAP_2=emb[cells_dz,2],
                       celltype6=ct6[sel], UBL3=ubl3_log1p[cells_dz])
  df_all <- df_all[!is.na(df_all$celltype6), ]
  cat("  UMAP (disease) cells:", nrow(df_all), " | UBL3 max =", round(max(df_all$UBL3,na.rm=TRUE),2), "\n")
  df_centroid <- df_all %>% group_by(celltype6) %>%
    summarise(UMAP_1=median(UMAP_1,na.rm=TRUE), UMAP_2=median(UMAP_2,na.rm=TRUE), .groups="drop")
  xr <- range(df_all$UMAP_1,na.rm=TRUE); yr <- range(df_all$UMAP_2,na.rm=TRUE)
  xpad <- diff(xr)*0.04; ypad <- diff(yr)*0.04

  rm(obj, counts, dm); gc()
  list(unit=cfg$unit, line1=paste0(display_id(cfg$unit)," | ",cfg$disease), line2=cfg$region,
       dot_df=dot_df, df_all=df_all, df_centroid=df_centroid,
       xlim=c(xr[1]-xpad, xr[2]+xpad), ylim=c(yr[1]-ypad, yr[2]+ypad))
}

# =============================================================================
# 5. 绘图函数
# =============================================================================
theme_dotplot <- function() {
  theme_classic(base_size = NBD_BASE_SIZE, base_family = PLOT_FONT) +
    theme(
      axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=NBD_XAXIS_SIZE,
                                 color="black", family=PLOT_FONT, face="italic"),
      axis.text.y = element_text(size=NBD_BASE_SIZE, color="black", family=PLOT_FONT),
      axis.title  = element_blank(),
      axis.line   = element_line(color="black", linewidth=NBD_LINE_WIDTH),
      axis.ticks  = element_line(color="black", linewidth=NBD_LINE_WIDTH),
      axis.ticks.length = unit(0.06,"cm"),
      plot.title  = element_text(size=NBD_BASE_SIZE, face="bold", hjust=0, color="black",
                                 family=PLOT_FONT, lineheight=1.10, margin=margin(b=1.5,unit="pt")),
      plot.title.position = "panel",
      panel.grid  = element_blank(), panel.background = element_blank(),
      plot.margin = margin(t=2.5,r=1,b=3.5,l=1,unit="mm"),  # 上下留白: 防竖排基因名被下排标题截断
      legend.key.size = unit(0.30,"cm"),
      legend.text  = element_text(size=NBD_BASE_SIZE, family=PLOT_FONT),
      legend.title = element_text(size=NBD_BASE_SIZE, family=PLOT_FONT, lineheight=1.05)
    )
}
draw_dotplot <- function(df, title, show_legend = FALSE) {
  p <- ggplot(df, aes(x=gene, y=celltype)) +
    geom_point(aes(size=pct_exp, color=avg_exp_scaled), shape=16) +
    scale_color_gradient2(low="#0072B2", mid="#FFFFFF", high="#D55E00", midpoint=0,
                          limits=c(-2.5,2.5), breaks=c(-2,-1,0,1,2), name="Avg expr\n(z-score)") +
    scale_size_continuous(range=c(0.3,2.4), limits=c(0,100), breaks=c(25,50,75,100),
                          name="% cells\nexpressing") +
    scale_x_discrete(expand=expansion(mult=c(0.04,0.04))) +
    scale_y_discrete(expand=expansion(mult=c(0.07,0.07))) +   # Y 轴全称(不缩写)
    labs(title=title) + theme_dotplot()
  if (show_legend) {
    p + guides(color=guide_colorbar(barwidth=unit(0.32,"cm"), barheight=unit(1.6,"cm"), title.position="top"),
               size =guide_legend(title.position="top")) +
      theme(legend.position="right", legend.box="horizontal", legend.direction="vertical",
            legend.key.size=unit(0.45,"cm"), legend.spacing.x=unit(5,"mm"),
            legend.title=element_text(face="plain"))
  } else p + theme(legend.position="none")
}

umap_theme <- theme_classic(base_size=8, base_family=PLOT_FONT) +
  theme(
    plot.title = element_text(face="bold", hjust=0, size=SUBPLOT_TITLE_SIZE),  # 左对齐到 y 轴
    plot.title.position = "panel",
    axis.title = element_text(size=6.5, color="black"),
    axis.text  = element_text(size=6, color="black"),
    axis.line  = element_line(linewidth=0.35, color="black"),
    axis.ticks = element_line(linewidth=0.30, color="black"),
    legend.title = element_text(size=NBD_BASE_SIZE, face="bold"),
    legend.text  = element_text(size=NBD_BASE_SIZE),
    legend.key.size = unit(3,"mm"),
    plot.margin = margin(1,1,1,1,"mm"), panel.grid = element_blank()
  )

draw_umap_ct <- function(df_all, df_centroid, xlim, ylim, show_legend=FALSE) {
  set.seed(SEED); df_c <- df_all[sample.int(nrow(df_all)), ]
  p <- ggplot(df_c, aes(UMAP_1, UMAP_2, color=celltype6)) +
    scattermore::geom_scattermore(pointsize=PT_CELL_MAIN, pixels=PIXELS_GRID, alpha=1) +
    # centroid 文字标签已移除(改由 panel b 底部共享 Cell type 图例对应), 避免压字重叠
    scale_color_manual(values=celltype6_palette_NBD, breaks=celltype6_levels_std,
                       limits=celltype6_levels_std, drop=FALSE, name="Cell type",
                       guide=guide_legend(nrow=2, byrow=TRUE,
                                          override.aes=list(size=3, alpha=1, shape=16))) +
    ggtitle("Cell type") + labs(x="UMAP_1", y="UMAP_2") +
    coord_cartesian(xlim=xlim, ylim=ylim, expand=FALSE) + umap_theme
  if (show_legend) p + theme(legend.position="bottom", legend.justification="center")
  else p + theme(legend.position="none")
}
draw_umap_ubl3 <- function(df_all, xlim, ylim, show_legend=FALSE) {
  d <- df_all[order(df_all$UBL3, na.last=FALSE), ]
  df_bg <- subset(d, UBL3<=0 | is.na(UBL3)); df_pos <- subset(d, UBL3>0)
  p <- ggplot() +
    scattermore::geom_scattermore(data=df_bg, aes(UMAP_1,UMAP_2), color=BG_GREY,
                                  pointsize=PT_UBL3_BG, pixels=PIXELS_GRID, alpha=1) +
    scattermore::geom_scattermore(data=df_pos, aes(UMAP_1,UMAP_2,color=UBL3),
                                  pointsize=PT_UBL3_MAIN, pixels=PIXELS_GRID, alpha=1) +
    scale_color_gradientn(colours=ubl3_palette, limits=c(0,UBL3_SHARED_MAX),
                          breaks=ubl3_breaks, labels=ubl3_labels, oob=scales::squish,
                          name="UBL3 log1p(CP10k)",
                          guide=guide_colorbar(barwidth=unit(22,"mm"), barheight=unit(1.8,"mm"),
                                               title.position="top", title.hjust=0.5,
                                               ticks.linewidth=0.3, frame.linewidth=0.3,
                                               frame.colour="grey30")) +
    ggtitle("UBL3 (>0 highlighted)") + labs(x="UMAP_1", y="UMAP_2") +
    coord_cartesian(xlim=xlim, ylim=ylim, expand=FALSE) + umap_theme
  if (show_legend) p + theme(legend.position="bottom", legend.justification="center")
  else p + theme(legend.position="none")
}

# ---- 共享图例 (无 scattermore 的 dummy 图 -> 抽 legend, 纯矢量, 不触发 scattermore bug) ----
legend_celltype <- function() {
  d <- data.frame(x=seq_along(celltype6_levels_std), y=1,
                  celltype6=factor(celltype6_levels_std, levels=celltype6_levels_std))
  p <- ggplot(d, aes(x, y, color=celltype6)) + geom_point(size=2) +
    scale_color_manual(values=celltype6_palette_NBD, breaks=celltype6_levels_std,
                       limits=celltype6_levels_std, drop=FALSE, name="Cell type",
                       guide=guide_legend(nrow=2, byrow=TRUE, override.aes=list(size=3, shape=16))) +
    theme_void(base_family=PLOT_FONT) +
    theme(legend.position="bottom", legend.title=element_text(size=NBD_BASE_SIZE, face="bold"),
          legend.text=element_text(size=NBD_BASE_SIZE), legend.key.size=unit(3,"mm"))
  cowplot::get_legend(p)
}
legend_ubl3 <- function() {
  d <- data.frame(x=1:2, y=1, UBL3=c(0, UBL3_SHARED_MAX))
  p <- ggplot(d, aes(x, y, color=UBL3)) + geom_point(size=2) +
    scale_color_gradientn(colours=ubl3_palette, limits=c(0, UBL3_SHARED_MAX),
                          breaks=ubl3_breaks, labels=ubl3_labels, oob=scales::squish,
                          name="UBL3 log1p(CP10k)",
                          guide=guide_colorbar(barwidth=unit(22,"mm"), barheight=unit(1.8,"mm"),
                                               title.position="top", title.hjust=0.5,
                                               ticks.linewidth=0.3, frame.linewidth=0.3, frame.colour="grey30")) +
    theme_void(base_family=PLOT_FONT) +
    theme(legend.position="bottom", legend.title=element_text(size=NBD_BASE_SIZE, face="bold"),
          legend.text=element_text(size=NBD_BASE_SIZE))
  cowplot::get_legend(p)
}

make_tag <- function(letter, y = 0.95) {        # <<< TUNE y: 让字母底部与小标题底部齐平
  cowplot::ggdraw() +
    cowplot::draw_label(letter, x=0.5, y=y, hjust=0.5, vjust=1,
                        fontfamily=PLOT_FONT, fontface="bold", color="black", size=OUTER_TAG_SIZE) +
    theme(plot.background=element_rect(fill="white",color=NA),
          panel.background=element_rect(fill="white",color=NA))
}
white_spacer <- function() cowplot::ggdraw() +
  theme(plot.background=element_rect(fill="white",color=NA),
        panel.background=element_rect(fill="white",color=NA), plot.margin=margin(0,0,0,0))

# =============================================================================
# 6. Builder: 给定单元列表 -> 组装 (a) DotPlot + (b) UMAP 复合图 (MN 输出)
# =============================================================================
build_fig <- function(units, out_dir, fig_tag, figure_title, legend_text,
                      dot_ncol, umap_ncol, fig_w, fig_h, a_h, b_h, gap = 4,
                      png_limit_mb = 10, a_tag_y = 0.99, b_tag_y = 0.99,  # <<< TUNE 字母高度
                      umap_frac = 0.80, dot_rel_widths = NULL, img_h = 0.80) {
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  cat("\n##################  BUILD:", fig_tag, " (", length(units), "units)  ##################\n")

  proc <- lapply(units, function(u) process_unit(CFG_ALL[[u]]))
  names(proc) <- units

  # ---- (a) DotPlot: 各单元点图(无图例) + 一个共享 color+size 图例 ----
  dot_plots <- lapply(proc, function(p) draw_dotplot(p$dot_df, paste0(p$line1,"\n",p$line2), show_legend=FALSE))
  dot_legend <- cowplot::get_legend(draw_dotplot(proc[[1]]$dot_df, "x", show_legend=TRUE))
  panel_a <- if (is.null(dot_rel_widths))
    cowplot::plot_grid(plotlist=c(dot_plots, list(dot_legend)), ncol=dot_ncol, align="hv", axis="tblr")
  else
    cowplot::plot_grid(plotlist=c(dot_plots, list(dot_legend)), ncol=dot_ncol, align="hv", axis="tblr",
                       rel_widths=dot_rel_widths)

  # ---- (b) UMAP: 每个单元的 [Cell type | UBL3] 先用 agg_png *单独* 渲染成 PNG
  #      (scattermore 在独立设备里渲染才稳定; 在复合 grob 内部会触发 rgb nul bug),
  #      再把 PNG 摆进网格. 共享图例用无 scattermore 的 dummy 图抽取 (纯矢量). ----
  tmp_dir <- file.path(out_dir, "_tmp_umap_png"); dir.create(tmp_dir, recursive=TRUE, showWarnings=FALSE)
  PAIR_W_MM <- 120; PAIR_H_MM <- 56; PAIR_DPI <- 600   # <<< TUNE: pair 长宽比(~2.14) 与分辨率
  blocks <- lapply(seq_along(proc), function(i) {
    p_ct   <- draw_umap_ct(proc[[i]]$df_all, proc[[i]]$df_centroid, proc[[i]]$xlim, proc[[i]]$ylim, FALSE)
    p_ubl3 <- draw_umap_ubl3(proc[[i]]$df_all, proc[[i]]$xlim, proc[[i]]$ylim, FALSE)
    pair   <- p_ct + p_ubl3 + patchwork::plot_layout(widths = c(1, 1))
    png_i  <- file.path(tmp_dir, paste0(names(proc)[i], "_pair.png"))
    ragg::agg_png(png_i, width=PAIR_W_MM/25.4, height=PAIR_H_MM/25.4, units="in", res=PAIR_DPI, background="white")
    print(pair); dev.off()
    # PNG 用 png::readPNG + rasterGrob 摆放 (轻量), 不走 magick -> 避免多图在 cairo_pdf 里重采样卡死
    base <- if (requireNamespace("png", quietly=TRUE)) {
      .raw <- png::readPNG(png_i)
      cowplot::ggdraw() + cowplot::draw_grob(grid::rasterGrob(
        .raw, x=grid::unit(0,"npc"), y=grid::unit(img_h,"npc"),
        width=grid::unit(1,"npc"), height=grid::unit(img_h,"npc"),
        hjust=0, vjust=1, interpolate=TRUE))   # 顶+左对齐: 图贴着标题
    } else {
      cowplot::ggdraw() + cowplot::draw_image(png_i, x=0, y=0, width=1, height=img_h,
                                              hjust=0, vjust=0, halign=0, valign=1)
    }
    base +
      cowplot::draw_label(paste0(proc[[i]]$line1, "\n", proc[[i]]$line2),   # 2 行, 与 (a) 点图标题一致
                          x=0.02, y=0.97, hjust=0, vjust=1, lineheight=1.05,  # <<< TUNE y: cohort 标题高度
                          fontfamily=PLOT_FONT, fontface="bold", size=COHORT_TITLE_SIZE, color="black")
  })
  umap_grid <- cowplot::plot_grid(plotlist=blocks, ncol=umap_ncol)

  ct_legend   <- legend_celltype()
  ubl3_legend <- legend_ubl3()
  legend_row  <- cowplot::plot_grid(ct_legend, ubl3_legend, nrow=1, rel_widths=c(1.5, 1))
  panel_b <- cowplot::plot_grid(umap_grid, legend_row, ncol=1, rel_heights=c(umap_frac, 1-umap_frac))

  # ---- 外层 A / B + 垂直拼接 ----
  pA <- cowplot::plot_grid(make_tag("A", a_tag_y), panel_a, ncol=2, rel_widths=c(0.032, 0.968))
  pB <- cowplot::plot_grid(make_tag("B", b_tag_y), panel_b, ncol=2, rel_widths=c(0.032, 0.968))
  fig <- cowplot::plot_grid(pA, white_spacer(), pB, ncol=1, rel_heights=c(a_h, gap, b_h))

  # ---- 导出: 矢量 PDF + 1000 dpi PNG (+300 dpi 备份) ----
  W <- fig_w/25.4; H <- fig_h/25.4
  out_pdf <- file.path(out_dir, paste0(fig_tag, "_MN_v15_vector.pdf"))
  out_png <- file.path(out_dir, paste0(fig_tag, "_MN_v15_1000dpi.png"))
  out_png300 <- file.path(out_dir, paste0(fig_tag, "_MN_v15_300dpi.png"))
  gc()
  cat("\U0001F4BE PDF...\n")
  tryCatch({ grDevices::cairo_pdf(out_pdf, width=W, height=H, family="Arial", bg="white")
             print(fig); dev.off() }, error=function(e) cat("\u274C PDF:", e$message,"\n"))
  save_png <- function(path, res) {
    ok <- tryCatch({ ragg::agg_png(path, width=W, height=H, units="in", res=res, background="white")
                     print(fig); dev.off(); TRUE },
                   error=function(e){ try(grDevices::dev.off(), silent=TRUE)
                                      cat("\u274C PNG", res, "dpi:", e$message, "\n"); FALSE })
    isTRUE(ok) && file.exists(path) && file.size(path) > 0
  }
  cat("\U0001F4BE PNG 1000dpi...\n")
  if (!save_png(out_png, 1000)) {
    cat("\u21BB 1000 dpi \u672A\u751F\u6210, \u56DE\u9000 600 dpi ->", basename(out_png), "\n")
    save_png(out_png, 600)
  }
  save_png(out_png300, 300)

  # ---- source data + legend sidecar ----
  combined <- do.call(rbind, lapply(names(proc), function(u){ d <- proc[[u]]$dot_df; d$unit <- u; d }))
  write.csv(combined, file.path(out_dir, paste0(fig_tag, "_PanelA_dotplot_data.csv")), row.names=FALSE)
  full_legend <- paste0(figure_title, "\n", legend_text)
  nwords <- length(strsplit(trimws(gsub("[^A-Za-z0-9]+"," ", full_legend))," ")[[1]])
  writeLines(full_legend, file.path(out_dir, paste0(fig_tag, "_legend_MN.txt")))

  png_mb <- if (file.exists(out_png)) file.size(out_png)/1024/1024 else NA
  pdf_kb <- if (file.exists(out_pdf)) file.size(out_pdf)/1024 else NA
  cat("\n\u2705", fig_tag, "DONE\n")
  cat("  PDF:", out_pdf, sprintf("(%.1f KB)\n", pdf_kb))
  cat("  PNG:", out_png, sprintf("(%.2f MB; MN limit %d MB)\n", png_mb, png_limit_mb))
  if (!is.na(png_mb) && png_mb > png_limit_mb)
    cat("  \u26A0 PNG >", png_limit_mb, "MB -> \u4F18\u5148\u63D0\u4EA4\u77E2\u91CF PDF, \u6216\u964D\u5230 600 dpi\n")
  cat("  Legend words:", nwords, "/ 300\n")
  invisible(list(fig=fig, png_mb=png_mb))
}

# =============================================================================
# 7. RUN — Figure 2 (2 units)
# =============================================================================
fig2_dir <- "D:/RNA/MNversion submission/Figure 2. Cell-type harmonization and UBL3 detectability across seven analytical units/Results"

fig2_title  <- "Figure 2. Cell-type harmonization and UBL3 detectability in representative analytical units."
fig2_legend <- paste0(
  "Two representative syn52082747 primary visual cortex (V1) units are shown: ",
  "Alzheimer's disease (AD) and progressive supranuclear palsy (PSP); the remaining five units appear in Supplementary Figure S3. ",
  "(a) Dotplots of 18 canonical lineage markers (columns) across six harmonized cell types ",
  "(rows; Astro, astrocytes; Endo, endothelial; Excit, excitatory neurons; Inhib, inhibitory neurons; Micro, microglia; Oligo, oligodendrocytes). ",
  "Dot size, percentage of cells expressing the marker; dot colour, mean log-normalized expression z-scored across cell types. ",
  "(b) UMAP embeddings coloured by harmonized cell type (left) and by UBL3 expression (right; cells with UBL3 > 0 highlighted on a shared log1p(CP10k) colour scale, grey = non-detecting). ",
  "A single shared cell-type legend and UBL3 colour bar apply to all UMAP panels. UBL3-detecting cells are present across neuronal and glial lineages. ",
  "V1, primary visual cortex. Source data, Supplementary Table S3."
)

build_fig(units=UNITS_FIG2, out_dir=fig2_dir, fig_tag="Figure2",
          figure_title=fig2_title, legend_text=fig2_legend,
          dot_ncol=3, umap_ncol=2,                       # a: 2 点图+图例 一行; b: 2 个 cohort 块 一行
          fig_w=170, fig_h=130, a_h=62, b_h=64, gap=4,   # 消除 panel B 大留白; a_h 略加防截断
          umap_frac=0.75, img_h=0.79, dot_rel_widths=c(1, 1, 0.5), png_limit_mb=10)

# =============================================================================
# 8. RUN — Supplementary Figure S3 (5 units)
# =============================================================================
sf3_dir <- "D:/RNA/MNversion submission/Supplementary Figure S3. negative UMAP/Results"

sf3_title  <- "Supplementary Figure S3. Cell-type harmonization and UBL3 detectability in the remaining five analytical units."
sf3_legend <- paste0(
  "The five non-representative units: GSE157827 (AD, middle frontal gyrus), GSE174367 (AD, prefrontal cortex), ",
  "syn21788402 (AD, entorhinal cortex), syn21788402 (AD, superior frontal gyrus) and syn52082747 (FTD, primary visual cortex). ",
  "(a) Dotplots of 18 canonical lineage markers (columns) across six harmonized cell types ",
  "(rows; Astro, astrocytes; Endo, endothelial; Excit, excitatory neurons; Inhib, inhibitory neurons; Micro, microglia; Oligo, oligodendrocytes). ",
  "Dot size, percentage of cells expressing the marker; dot colour, mean log-normalized expression z-scored across cell types. ",
  "(b) UMAP embeddings coloured by harmonized cell type (left of each pair) and by UBL3 expression (right; cells with UBL3 > 0 highlighted on a shared log1p(CP10k) colour scale, grey = non-detecting). ",
  "A single shared cell-type legend and UBL3 colour bar apply to all UMAP panels. FTD, frontotemporal dementia. Source data, Supplementary Table S3."
)

build_fig(units=UNITS_SF3, out_dir=sf3_dir, fig_tag="SuppFigureS3",
          figure_title=sf3_title, legend_text=sf3_legend,
          dot_ncol=2, umap_ncol=3,                       # a: 5 点图+图例 (2 列, 更宽以容全称, 3 行); b: 5 块 (3 列, 2 行)
          fig_w=170, fig_h=216, a_h=132, b_h=80, gap=4,  # 全称更宽 -> 点图 2 列; 总高 <225
          umap_frac=0.81, img_h=0.78, png_limit_mb=20)

cat("\n\U0001F389 ALL DONE: Figure 2 + Supplementary Figure S3\n")
