###############################################################################
# 01_per_dataset / syn21788402  (Alzheimer's disease; entorhinal cortex EC +
# superior frontal gyrus SFG)  -- PER-DATASET DOWNSTREAM / UMAP SCRIPT
#
# IMPORTANT: this script CONSUMES the harmonized 6-cell-type objects; it does
#            NOT create them. It rebuilds a UMAP from the authors' aligned object
#            and imports the celltype6 labels FROM the pre-existing stepH objects
#            (attach_celltype6_from_stepH -> readRDS) to render the per-dataset
#            cell-type / UBL3 UMAPs.
#
# Reads (does NOT write):
#   .../syn21788402/resultsmodify/stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds
#   .../syn21788402/resultsmodify/stepH_syn21788402_SFG_obj_celltype6.rds
# These two objects are the inputs consumed by code/02_integrated_figures/.
#
# Source     : Synapse accession syn21788402 (authors' QC'd / aligned / clustered
#              objects, e.g. sce.EC.scAlign.assigned.rds, + matched-cells CSVs).
# Paths      : all paths point to the data-project location on disk; source data
#              and intermediate objects are not stored in this repository.
# Environment: R + Seurat (v5) + SingleCellExperiment; see sessionInfo logs.
#
# NOTE: the upstream script that actually CREATES the two stepH objects above
#       (author object -> marker scoring -> 7 classes -> celltype6 ->
#       saveRDS stepH) is a SEPARATE file and should be added to this folder for
#       full reproducibility. Kept faithful to the script as run.
###############################################################################

#直接用原文提供的质控后的EC区的：sce.EC.scAlign.assigned.rds：进行 UMAP + 6 类 cell type + UBL3 对照 vs AD
############################################################
## syn21788402 EC_allCells：对应 GSE157827 的第 6 + 7 章
##
## 起点：rawdata/sce.EC.scAlign.assigned.rds
##   - 作者已完成 QC + CCA/scAlign 对齐 + clustering
##   - reducedDims: CCA / CCA.ALIGNED / TSNE
##   - colData: SampleID, BraakStage, clusterAssignment, clusterCellType …
##
## 目标：
##   1) 使用 SCE.counts + CCA.ALIGNED 建一个 Seurat UMAP 对象
##   2) 用 GSE157827 模版的 marker 打分 → 7 类 → 合并为 6 大类：
##        Astrocytes / Excitatory neurons / Microglia /
##        Endothelial / Inhibitory neurons / Oligodendrocytes
##   3) 从 counts 计算 UBL3 的 log1p(CPM)
##      （公式完全与 GSE157827 第 7 章一致）
##   4) 画 “Cell type UMAP + UBL3 UMAP(>0)” 两联图
############################################################

###############################################################################
# syn21788402：SFG 与 EC 两个部位
# 第 7 章（2）UMAP：6 个细胞类型 + UBL3（>0 高亮）
# 要求：严格模仿 GSE157827 模版（主题、轴标签、配色、legend两行、布局、导出PNG+PDF）
#
# 输出目录（按你要求固定）：
#   - SFG：D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/results_SFG
#   - EC ：D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/results
#
# 关键处理：
#   - SFG：使用你现成的 stepF（UMAP）+ stepH（celltype6 + counts）
#   - EC ：由于旧 stepF_EC 可能 R 版本不兼容（unknown type 212），本脚本会
#          直接从 rawdata/sce.EC.scAlign.assigned.rds 重建一个新的 stepF_EC（UMAP），
#          然后从你指定 stepH_EC 搬入 celltype6，再按模版出图。
###############################################################################

## ============================================================
## 0. 基本设置：随机数种子 + 加载 R 包（尽量贴近你模版）
## ============================================================
SEED <- 20251023
set.seed(SEED)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(data.table)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(SingleCellExperiment)  # ★EC 重建 stepF 需要
})

## ============================================================
## 0.1 全局统一的 celltype6 名称 + 颜色（模版核心）
## ============================================================
celltype6_levels_std <- c(
  "Astrocytes",
  "Excitatory neurons",
  "Microglia",
  "Endothelial",
  "Inhibitory neurons",
  "Oligodendrocytes"
)

celltype6_palette_std <- c(
  "Astrocytes"          = "#FF8E8E",
  "Excitatory neurons"  = "#09BB3C",
  "Microglia"           = "#36B0E1",
  "Endothelial"         = "#B8A109",
  "Inhibitory neurons"  = "#00BFC4",
  "Oligodendrocytes"    = "#E16AFC"
)

standardize_celltype6 <- function(x) {
  x <- as.character(x)
  x[x %in% c("Astrocytes", "Astro", "Astrocyte")] <- "Astrocytes"
  x[x %in% c("Excitatory neurons", "Excitatory", "Ex_neuron", "Excit")] <- "Excitatory neurons"
  x[x %in% c("Inhibitory neurons", "Inhibitory", "Inh_neuron", "Inhib")] <- "Inhibitory neurons"
  x[x %in% c("Microglia", "Micro", "Microgl")] <- "Microglia"
  x[x %in% c("Endothelial", "Endothelial cells", "Endo")] <- "Endothelial"
  x[x %in% c("Oligodendrocytes", "Oligo", "Oligodendrocyte", "Oligodendro")] <- "Oligodendrocytes"
  factor(x, levels = celltype6_levels_std)
}

## ============================================================
## 0.2 统一 UMAP 视觉主题（严格模版）
## ============================================================
umap_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    legend.title = element_text(face = "bold"),
    panel.border = element_blank(),
    panel.grid   = element_blank()
  )

## ============================================================
## 0.3 兼容函数：Seurat v4/v5 读取 counts（避免 Assay5 差异）
## ============================================================
get_counts_compat <- function(obj, assay = "RNA") {
  DefaultAssay(obj) <- assay
  tryCatch(
    GetAssayData(obj, assay = assay, slot = "counts"),
    error = function(e) {
      message("GetAssayData(slot='counts') 出错，改用 LayerData(layer='counts')")
      LayerData(obj[[assay]], layer = "counts")
    }
  )
}

###############################################################################
# 1) 工具函数：把 celltype6 从 stepH 搬到 stepF（按细胞名对齐）
###############################################################################
attach_celltype6_from_stepH <- function(objF, stepH_rds) {
  objH <- readRDS(stepH_rds)
  if (!"celltype6" %in% colnames(objH@meta.data)) stop("stepH 中没有 celltype6：", stepH_rds)
  
  common <- intersect(colnames(objF), colnames(objH))
  if (length(common) == 0) stop("stepF 与 stepH 细胞名无法对齐（common=0）")
  
  objF$celltype6 <- objH@meta.data[colnames(objF), "celltype6", drop = TRUE]
  objF
}

###############################################################################
# 2) 工具函数：从 SCE 重建 stepF（解决 EC stepF 的 unknown type 212）
###############################################################################
rebuild_stepF_from_sce <- function(sce_rds, region_tag = "EC", seed = 20251023) {
  sce <- readRDS(sce_rds)
  
  if (!"counts" %in% assayNames(sce)) stop("SCE 中没有 counts assay：", sce_rds)
  counts <- assay(sce, "counts")
  
  # 优先 CCA.ALIGNED，其次 CCA
  use_red <- NULL
  if ("CCA.ALIGNED" %in% reducedDimNames(sce)) {
    use_red <- "CCA.ALIGNED"
  } else if ("CCA" %in% reducedDimNames(sce)) {
    use_red <- "CCA"
  } else {
    stop("SCE 中找不到 CCA.ALIGNED 或 CCA：", sce_rds)
  }
  
  emb <- reducedDim(sce, use_red)
  emb <- emb[colnames(sce), , drop = FALSE]
  colnames(emb) <- paste0("CCAALN_", seq_len(ncol(emb)))
  
  # 建 Seurat
  obj <- CreateSeuratObject(
    counts       = counts,
    project      = paste0("syn21788402_", region_tag),
    min.cells    = 0,
    min.features = 0
  )
  
  # 搬 meta（可选，但建议）
  meta <- as.data.frame(colData(sce))
  meta <- meta[colnames(obj), , drop = FALSE]
  obj  <- AddMetaData(obj, metadata = meta)
  
  # 挂 embedding
  obj[[use_red]] <- CreateDimReducObject(
    embeddings = emb[colnames(obj), , drop = FALSE],
    key       = "CCAALN_",
    assay     = "RNA"
  )
  
  # RunUMAP
  set.seed(seed)
  obj <- RunUMAP(
    obj,
    reduction = use_red,
    dims      = 1:ncol(emb),
    verbose   = FALSE
  )
  
  obj
}

###############################################################################
# 3) 核心函数：严格按 GSE157827 模版画两联 UMAP 并输出 PNG+PDF
###############################################################################
plot_umap_2panel_likeGSE157827 <- function(region = c("SFG","EC"),
                                           res_dir,
                                           obj_umap,        # ★已经含 umap reduction + celltype6 的 Seurat 对象
                                           counts_obj_rds) { # ★stepH 对象（提供 counts，用于重算 UBL3_log1p）
  region <- match.arg(region)
  if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)
  
  if (!"umap" %in% Reductions(obj_umap)) stop("对象里没有 reduction='umap'")
  if (!"celltype6" %in% colnames(obj_umap@meta.data)) stop("对象 meta.data 中没有 celltype6")
  
  # 读取 counts 对象（严格按模版重算 UBL3_log1p）
  obj_expr <- readRDS(counts_obj_rds)
  DefaultAssay(obj_expr) <- "RNA"
  
  rna_counts <- get_counts_compat(obj_expr, "RNA")
  cat(">>>", region, "counts 维度：", nrow(rna_counts), "基因 ×", ncol(rna_counts), "细胞\n")
  
  # 找 UBL3 ENSEMBL
  map_ubl3 <- AnnotationDbi::select(
    org.Hs.eg.db, keys = "UBL3", keytype = "SYMBOL", columns = "ENSEMBL"
  )
  ubl3_id <- map_ubl3$ENSEMBL[1]
  cat(">>>", region, "UBL3 ENSG ID:", ubl3_id, "\n")
  
  # 定位 UBL3 行：优先 ENSG，否则 SYMBOL
  if (!is.na(ubl3_id) && ubl3_id %in% rownames(rna_counts)) {
    gene_row <- ubl3_id
  } else if ("UBL3" %in% rownames(rna_counts)) {
    gene_row <- "UBL3"
  } else {
    stop(region, "：counts 中找不到 UBL3（既无 ENSEMBL 也无 SYMBOL）")
  }
  cat(">>>", region, "实际使用的 UBL3 行名：", gene_row, "\n")
  
  # 对齐细胞
  common_cells <- intersect(colnames(rna_counts), colnames(obj_umap))
  cat(">>>", region, "共有细胞数：", length(common_cells), "\n")
  if (length(common_cells) == 0) stop(region, "：common_cells=0，counts 与 UMAP 对象无法对齐")
  
  # 模版公式：log1p((raw/libsize)*1e4)
  raw_vec  <- as.numeric(rna_counts[gene_row, common_cells, drop = FALSE][1, ])
  lib_size <- Matrix::colSums(rna_counts[, common_cells, drop = FALSE])
  UBL3_norm <- log1p((raw_vec / lib_size) * 1e4)
  
  # 写入 UMAP 对象（按 UMAP 细胞顺序）
  ubl3_for_umap <- rep(NA_real_, ncol(obj_umap))
  names(ubl3_for_umap) <- colnames(obj_umap)
  ubl3_for_umap[common_cells] <- UBL3_norm
  obj_umap$UBL3_log1p <- ubl3_for_umap
  
  cat(">>>", region, "UBL3_log1p 范围：",
      min(obj_umap$UBL3_log1p, na.rm = TRUE), "~",
      max(obj_umap$UBL3_log1p, na.rm = TRUE), "\n")
  
  # celltype6 标准化（模版）
  obj_umap$celltype6 <- standardize_celltype6(obj_umap$celltype6)
  
  cat(">>>", region, "标准化后的 celltype6：\n")
  print(table(obj_umap$celltype6, useNA = "ifany"))
  
  ## -------- Panel A：Cell type UMAP（严格模版）--------
  if ("integrated" %in% Assays(obj_umap)) {
    DefaultAssay(obj_umap) <- "integrated"
  }
  
  p_celltype <- DimPlot(
    obj_umap,
    reduction  = "umap",
    group.by   = "celltype6",
    label      = TRUE,
    label.size = 4.5,
    repel      = TRUE,
    raster     = FALSE,
    pt.size    = 0.1
  ) +
    scale_color_manual(
      values = celltype6_palette_std,
      breaks = celltype6_levels_std,
      limits = celltype6_levels_std,
      drop   = TRUE
    ) +
    ggtitle("Cell Type UMAP") +
    labs(color = "Cell type") +
    xlab("UMAP_1") + ylab("UMAP_2") +
    umap_theme +
    theme(
      legend.position = "bottom",
      legend.box      = "horizontal"
    ) +
    guides(
      color = guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(size = 5)
      )
    )
  
  ## -------- Panel B：UBL3 expression UMAP（>0 高亮，严格模版）--------
  DefaultAssay(obj_umap) <- "RNA"
  
  umap_coords <- Embeddings(obj_umap, "umap")
  df_umap <- data.frame(
    umap_1 = umap_coords[, 1],
    umap_2 = umap_coords[, 2],
    UBL3   = obj_umap$UBL3_log1p
  )
  df_bg  <- subset(df_umap, is.na(UBL3) | UBL3 <= 0)
  df_pos <- subset(df_umap, UBL3 >  0)
  
  break_vals   <- c(0, 1, 2, 3, 4)
  break_labels <- c("0", "1", "2", "3", "4+")
  
  p_ubl3 <- ggplot() +
    geom_point(
      data  = df_bg,
      aes(x = umap_1, y = umap_2),
      color = "grey95",
      size  = 0.15,
      alpha = 0.4
    ) +
    geom_point(
      data  = df_pos,
      aes(x = umap_1, y = umap_2, color = UBL3),
      size  = 0.25,
      alpha = 0.9
    ) +
    scale_color_gradientn(
      colours = c("#FEE0D2", "#FC9272", "#CB181D"),
      limits  = c(0, max(df_umap$UBL3, na.rm = TRUE)),
      breaks  = break_vals,
      labels  = break_labels,
      name    = "UBL3"
    ) +
    ggtitle("UBL3 expression (>0 highlighted)") +
    xlab("UMAP_1") + ylab("UMAP_2") +
    umap_theme +
    theme(
      legend.position = "bottom",
      legend.box      = "horizontal"
    )
  
  ## -------- 拼图并输出（严格模版：12×6，PNG 300dpi + PDF）--------
  p_AB <- p_celltype + p_ubl3 +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(tag_levels = "A")
  
  out_png <- file.path(res_dir, paste0("Fig_UMAP_Celltype_UBL3_syn21788402_", region, ".png"))
  out_pdf <- file.path(res_dir, paste0("Fig_UMAP_Celltype_UBL3_syn21788402_", region, ".pdf"))
  
  ggsave(out_png, p_AB, width = 12, height = 6, dpi = 300)
  ggsave(out_pdf, p_AB, width = 12, height = 6)
  
  cat("🎉", region, "已输出 PNG：", out_png, "\n")
  cat("🎉", region, "已输出 PDF：", out_pdf, "\n")
  
  invisible(p_AB)
}

###############################################################################
# 4) 主流程：一次性同时出 SFG + EC，并分别保存到两个文件夹
###############################################################################

## ---------- 固定工程路径 ----------
proj_dir  <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402"
raw_dir   <- file.path(proj_dir, "rawdata")
res_SFG   <- file.path(proj_dir, "results_SFG")  # 你指定
res_EC    <- file.path(proj_dir, "results")      # 你指定

if (!dir.exists(res_SFG)) dir.create(res_SFG, recursive = TRUE)
if (!dir.exists(res_EC))  dir.create(res_EC,  recursive = TRUE)

## =========================
## 4.1 SFG：直接用你现成的 stepF + stepH
## =========================
stepF_SFG <- file.path(res_SFG, "stepF_syn21788402_SFG_counts_umap.rds")
stepH_SFG <- file.path(res_SFG, "stepH_syn21788402_SFG_obj_celltype6.rds")

stopifnot(file.exists(stepF_SFG))
stopifnot(file.exists(stepH_SFG))

# 读 stepF（UMAP）
objSFG <- readRDS(stepF_SFG)

# 如果 stepF 里没有 celltype6，就从 stepH 搬入（你之前就是这个问题）
if (!"celltype6" %in% colnames(objSFG@meta.data)) {
  cat(">>> SFG：stepF 中缺少 celltype6，正在从 stepH 搬入...\n")
  objSFG <- attach_celltype6_from_stepH(objSFG, stepH_SFG)
  
  # 保存一个修复版（可选，但推荐，后面就不怕再丢）
  fp_fix <- file.path(res_SFG, "stepF_syn21788402_SFG_counts_umap_with_celltype6.rds")
  saveRDS(objSFG, fp_fix, compress = "xz")
  cat("✅ SFG：已保存修复版 stepF：", fp_fix, "\n")
}

# 出图（严格模版）
plot_umap_2panel_likeGSE157827(
  region        = "SFG",
  res_dir       = res_SFG,
  obj_umap      = objSFG,
  counts_obj_rds= stepH_SFG
)

## =========================
## 4.2 EC：从 SCE 重建 stepF（避免 unknown type 212）
## =========================
sce_EC_rds <- file.path(raw_dir, "sce.EC.scAlign.assigned.rds")  # 你截图中确认存在
stepH_EC   <- file.path(res_EC, "stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds")  # 你指定

stopifnot(file.exists(sce_EC_rds))
stopifnot(file.exists(stepH_EC))

cat("\n>>> EC：开始从 SCE 重建 stepF（避免旧 stepF 的 unknown type 212）...\n")
objEC <- rebuild_stepF_from_sce(
  sce_rds    = sce_EC_rds,
  region_tag = "EC",
  seed       = SEED
)

# 从 stepH 搬入 celltype6（保证与既有定义一致）
objEC <- attach_celltype6_from_stepH(objEC, stepH_EC)

# 保存重建版 stepF（可选但推荐）
stepF_EC_rebuilt <- file.path(res_EC, "stepF_syn21788402_EC_counts_umap_rebuilt_with_celltype6.rds")
saveRDS(objEC, stepF_EC_rebuilt, compress = "xz")
cat("✅ EC：已保存重建版 stepF：", stepF_EC_rebuilt, "\n")

# 出图（严格模版）
plot_umap_2panel_likeGSE157827(
  region        = "EC",
  res_dir       = res_EC,
  obj_umap      = objEC,
  counts_obj_rds= stepH_EC
)

cat("\n✅ syn21788402：SFG 与 EC 两个部位的两联 UMAP 已同时完成并分别输出到指定文件夹。\n")











#第 2章（0）：6 个细胞类型 × 每个样本的直方图（Y = cell 数）
## ===========================
## 2.0 EC：6 个细胞类型 × 每个样本的 UBL3>0 直方图（Y = Cell count）
## ===========================
SEED <- 20251023; set.seed(SEED)

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(dplyr)
})

res_dir_EC <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/results"

## ------ 1) 读回 EC 的 counts 和 meta（只含 Braak0/6 = AD / NC）------
cnt_EC <- readRDS(file.path(res_dir_EC,
                            "stepP_syn21788402_counts_raw_matched.rds"))
md_EC  <- read.csv(file.path(res_dir_EC,
                             "stepP_syn21788402_matched_cells_meta.csv"),
                   row.names = 1, stringsAsFactors = FALSE)

## 把 NC 改成 Control，再设定因子顺序（确保 AD 在前）
md_EC$group[md_EC$group == "NC"] <- "Control"
md_EC$group <- factor(md_EC$group, levels = c("AD","Control"))

## 统一 celltype6 名称
md_EC$celltype6 <- standardize_celltype6(md_EC$celltype6)

cat(">>> [EC 10.1] 每组细胞数：\n")
print(table(md_EC$group, useNA = "ifany"))

## ------ 2) 计算 UBL3 的 log1p(CPM) ------
if (!"UBL3" %in% rownames(cnt_EC)) {
  stop("EC counts 中找不到 UBL3 行，请检查。")
}

raw_vec  <- as.numeric(cnt_EC["UBL3", ])
lib_size <- Matrix::colSums(cnt_EC)

UBL3_log1p_EC <- log1p((raw_vec / lib_size) * 1e4)

## ------ 3) 构建每细胞数据框 df_all_EC ------
df_all_EC <- data.frame(
  cell      = colnames(cnt_EC),
  UBL3      = UBL3_log1p_EC,
  raw       = raw_vec,
  sample    = md_EC$sample,
  celltype6 = md_EC$celltype6,
  group     = md_EC$group,
  stringsAsFactors = FALSE
)

## 只保留 UBL3 > 0 的细胞
df_pos_EC <- subset(df_all_EC, UBL3 > 0 & !is.na(celltype6) & !is.na(group))
cat(">>> [EC 10.1] UBL3>0 的细胞数：", nrow(df_pos_EC), "\n")

## ------ 4) 画图：每个 celltype 一张，按样本分面 ------
## 颜色：AD 橙色，Control 蓝色（与 GSE157827 一致）
group_colors_EC <- c(AD = "#D55E00", Control = "#0072B2")

df_pos_EC$sample_group <- with(df_pos_EC, paste0(sample, "_", group))
df_pos_EC$sample_group <- factor(df_pos_EC$sample_group,
                                 levels = sort(unique(df_pos_EC$sample_group)))

celltypes_EC <- sort(unique(df_pos_EC$celltype6))

out_dir1_EC <- file.path(res_dir_EC,
                         "Fig10_1_UBL3_hist_per_celltype_noZero_syn21788402_EC")
dir.create(out_dir1_EC, showWarnings = FALSE, recursive = TRUE)

sanitize_name <- function(x) {
  x <- gsub(" ", "_", x)
  x <- gsub("/", "_", x)
  x
}

for (ct in celltypes_EC) {
  df_ct <- df_pos_EC[df_pos_EC$celltype6 == ct, ]
  if (nrow(df_ct) == 0) next
  
  cat(">>> [EC 10.1] 绘制细胞类型：", ct, "；细胞数 =", nrow(df_ct), "\n")
  
  p_ct <- ggplot(df_ct, aes(x = UBL3, fill = group)) +
    geom_histogram(
      bins     = 40,
      position = "identity",
      alpha    = 0.7,
      color    = "grey30"
    ) +
    facet_wrap(~ sample_group, ncol = 3) +
    scale_fill_manual(values = group_colors_EC, name = "Group") +
    labs(
      title = paste0("UBL3>0 distribution in ", ct, " (syn21788402 EC)"),
      x     = "UBL3 expression (log1p CPM)",
      y     = "Cell count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 16),
      strip.text      = element_text(face = "bold", size = 10),
      legend.position = "top",
      legend.title    = element_text(face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = file.path(out_dir1_EC,
                         paste0("Fig10_1_UBL3_noZero_", sanitize_name(ct), "_syn21788402_EC.png")),
    plot   = p_ct,
    width  = 10,
    height = 8,
    dpi    = 300
  )
}

cat("✅ [EC 10.1] EC 6 个 celltype × 单样本直方图已生成：", out_dir1_EC, "\n")





#SFG的单样本直方图
# 第 2 章（0）：SFG – 6 个细胞类型 × 每个样本的直方图（Y = Cell count）
## ===========================
## 10.1 SFG：6 个细胞类型 × 每个样本的 UBL3>0 直方图（Y = Cell count）
## ===========================
SEED <- 20251023; set.seed(SEED)

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(dplyr)
})

res_dir_SFG <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/results_SFG"

## ------ 1) 读回 SFG counts + meta（只含 Braak0/6 = AD / NC）------
cnt_SFG <- readRDS(file.path(res_dir_SFG,
                             "stepP_syn21788402_SFG_counts_raw_matched.rds"))
md_SFG  <- read.csv(file.path(res_dir_SFG,
                              "stepP_syn21788402_SFG_matched_cells_meta.csv"),
                    row.names = 1, stringsAsFactors = FALSE)

## 把 NC 改成 Control，再设定因子顺序
md_SFG$group[md_SFG$group == "NC"] <- "Control"
md_SFG$group <- factor(md_SFG$group, levels = c("AD","Control"))

## 统一 celltype6 名称
md_SFG$celltype6 <- standardize_celltype6(md_SFG$celltype6)

cat(">>> [SFG 10.1] 每组细胞数：\n")
print(table(md_SFG$group, useNA = "ifany"))

## ------ 2) 计算 UBL3 的 log1p(CPM) ------
if (!"UBL3" %in% rownames(cnt_SFG)) {
  stop("SFG counts 中找不到 UBL3 行，请检查。")
}

raw_vec  <- as.numeric(cnt_SFG["UBL3", ])
lib_size <- Matrix::colSums(cnt_SFG)

UBL3_log1p_SFG <- log1p((raw_vec / lib_size) * 1e4)

## ------ 3) 构建每细胞数据框 df_all_SFG ------
df_all_SFG <- data.frame(
  cell      = colnames(cnt_SFG),
  UBL3      = UBL3_log1p_SFG,
  raw       = raw_vec,
  sample    = md_SFG$sample,
  celltype6 = md_SFG$celltype6,
  group     = md_SFG$group,
  stringsAsFactors = FALSE
)

## 只保留 UBL3>0 的细胞
df_pos_SFG <- subset(df_all_SFG, UBL3 > 0 & !is.na(celltype6) & !is.na(group))
cat(">>> [SFG 10.1] UBL3>0 的细胞数：", nrow(df_pos_SFG), "\n")

## ------ 4) 画图：每个 celltype 一张，按样本分面 ------
group_colors_SFG <- c(AD = "#D55E00", Control = "#0072B2")

df_pos_SFG$sample_group <- with(df_pos_SFG, paste0(sample, "_", group))
df_pos_SFG$sample_group <- factor(df_pos_SFG$sample_group,
                                  levels = sort(unique(df_pos_SFG$sample_group)))

celltypes_SFG <- sort(unique(df_pos_SFG$celltype6))

out_dir1_SFG <- file.path(res_dir_SFG,
                          "Fig10_1_UBL3_hist_per_celltype_noZero_syn21788402_SFG")
dir.create(out_dir1_SFG, showWarnings = FALSE, recursive = TRUE)

sanitize_name <- function(x) {
  x <- gsub(" ", "_", x)
  x <- gsub("/", "_", x)
  x
}

for (ct in celltypes_SFG) {
  df_ct <- df_pos_SFG[df_pos_SFG$celltype6 == ct, ]
  if (nrow(df_ct) == 0) next
  
  cat(">>> [SFG 10.1] 绘制细胞类型：", ct, "；细胞数 =", nrow(df_ct), "\n")
  
  p_ct <- ggplot(df_ct, aes(x = UBL3, fill = group)) +
    geom_histogram(
      bins     = 40,
      position = "identity",
      alpha    = 0.7,
      color    = "grey30"
    ) +
    facet_wrap(~ sample_group, ncol = 3) +
    scale_fill_manual(values = group_colors_SFG, name = "Group") +
    labs(
      title = paste0("UBL3>0 distribution in ", ct, " (syn21788402 SFG)"),
      x     = "UBL3 expression (log1p CPM)",
      y     = "Cell count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 16),
      strip.text      = element_text(face = "bold", size = 10),
      legend.position = "top",
      legend.title    = element_text(face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = file.path(out_dir1_SFG,
                         paste0("Fig10_1_UBL3_noZero_", sanitize_name(ct), "_syn21788402_SFG.png")),
    plot   = p_ct,
    width  = 10,
    height = 8,
    dpi    = 300
  )
}

cat("✅ [SFG 10.1] SFG 6 个 celltype × 单样本直方图已生成：", out_dir1_SFG, "\n")






#第 2 章（1）：重叠直方图（Y = cell count，含检验方法 + FDR）
###############################################################################
## NO1_syn21788402_EC_SFG_overlap_hist_count.R
###############################################################################

###############################################################################
## NO1_syn21788402_EC_SFG_OverlapHistCount_ADvsControl_cell_and_donor_level.R
##
## 【严格模仿 GSE157827 模版】
##  - overlap hist（Y=Cell count）
##  - facet_wrap(~celltype6, scales=free_y)
##  - 右上角统计标签（Mann–Whitney U + BH校正Padj）
##  - 两张图：cell-level + donor-level
##  - UBL3 表达：log1p((raw/lib_size)*1e4) = log1p(CP10k)
##
## 【syn21788402 特殊点】
##  - 两个脑区：EC 与 SFG，相当于两套分析
##  - 组别在 stepP 的 meta csv 里是 AD / NC
##  - 但你要求：最终所有图里对照组名称固定显示为 Control
##    => NC/CTRL/... 都统一映射为 Control（只改变显示与统计标签，不改变数据含义）
##
## 【最关键的稳定性修复（不偏离模版，仅改变实现方式）】
##  - Seurat v5 的 Assay5 多 layers 在 subset(obj) 后，counts layers 可能不同步
##    导致 .subscript.2ary(x,i,j) out of bounds
##  - 解决：不要先 subset Seurat 对象再取 counts
##          而是：
##          (1) 从“原始 obj（未subset）”合并得到全量 rna_counts_all
##          (2) 用 cells_keep 在矩阵层面取子集 rna_counts <- rna_counts_all[, cells_keep]
##          (3) 后续流程完全照 GSE157827 模版（CP10k/log1p/画图/统计）
###############################################################################

rm(list=ls()); gc()
Sys.setenv(LANG="en")
SEED <- 20251023; set.seed(SEED)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(ragg)
})

## =========================
## 0) 全局参数（与你模版一致）
## =========================
dataset_tag <- "syn21788402"
gene_symbol <- "UBL3"
binwidth    <- 0.2
disease     <- "AD"

## =========================
## 1) 两个脑区的根目录（严格不弄乱）
## =========================
root_EC  <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/results"
root_SFG <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/results_SFG"

## RDS（对象）
obj_fp_EC  <- file.path(root_EC,  "stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds")
obj_fp_SFG <- file.path(root_SFG, "stepH_syn21788402_SFG_obj_celltype6.rds")
stopifnot(file.exists(obj_fp_EC))
stopifnot(file.exists(obj_fp_SFG))

## stepP meta（你截图里：有 sample 列 + group(AD/NC)）
meta_fp_EC  <- file.path(root_EC,  "stepP_syn21788402_matched_cells_meta.csv")
meta_fp_SFG <- file.path(root_SFG, "stepP_syn21788402_SFG_matched_cells_meta.csv")
stopifnot(file.exists(meta_fp_EC))
stopifnot(file.exists(meta_fp_SFG))

## =========================
## 2) GSE157827 模版同款：合并 counts layers（Seurat v5 多 layers 兼容）
##    （这里保留模版逻辑：合并所有 counts.* → gene×cell 矩阵）
## =========================
get_counts_matrix_allcells <- function(obj, assay="RNA", log_fp=NULL, out_dir=NULL) {
  
  a <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(a), error=function(e) character(0))
  counts_layers <- layers[grepl("^counts", layers)]
  
  ## 如果没有 layers，就退回传统 counts slot
  if (length(counts_layers) == 0) {
    m2 <- tryCatch(Seurat::GetAssayData(obj, assay=assay, slot="counts"), error=function(e) NULL)
    if (!is.null(m2)) return(m2)
    stop("❌ 无法获取 counts：没有 counts.* layers，GetAssayData(counts) 也失败。")
  }
  
  if (!is.null(log_fp)) {
    sink(log_fp, append=TRUE)
    cat("Detected counts layers:\n"); print(counts_layers)
    sink()
  }
  
  mats <- list()
  for (ly in counts_layers) {
    m <- tryCatch(SeuratObject::LayerData(a, layer=ly), error=function(e) NULL)
    if (is.null(m)) next
    if (!is.null(dim(m)) && length(dim(m))==2) {
      mats[[ly]] <- m
      if (!is.null(log_fp)) {
        sink(log_fp, append=TRUE)
        cat("layer:", ly, " dim=", paste(dim(m), collapse="x"), "\n")
        sink()
      }
    }
  }
  if (length(mats)==0) stop("❌ counts layers 存在，但读取 LayerData 全部失败。")
  
  ## 基因对齐：以第一个 layer 的 gene 顺序为标准
  ref_genes <- rownames(mats[[1]])
  if (is.null(ref_genes) || length(ref_genes)==0) stop("❌ counts layer 的 rownames 为空。")
  
  for (k in names(mats)) {
    if (!identical(rownames(mats[[k]]), ref_genes)) {
      m0 <- mats[[k]]
      m_aligned <- Matrix::Matrix(0, nrow=length(ref_genes), ncol=ncol(m0), sparse=TRUE)
      rownames(m_aligned) <- ref_genes
      colnames(m_aligned) <- colnames(m0)
      common <- intersect(ref_genes, rownames(m0))
      if (length(common) > 0) m_aligned[common, ] <- m0[common, , drop=FALSE]
      mats[[k]] <- m_aligned
    }
  }
  
  ## 合并所有 layer 的列（cells）
  mat_all <- if (length(mats)==1) mats[[1]] else Reduce(Matrix::cbind2, mats)
  
  ## 去除重复 cell
  if (!is.null(colnames(mat_all))) {
    dup <- duplicated(colnames(mat_all))
    if (any(dup)) mat_all <- mat_all[, !dup, drop=FALSE]
  }
  
  ## 对齐到 obj 的 cell 顺序（用 match 数值索引，避免 Matrix 字符索引越界）
  all_cells <- colnames(obj)
  if (length(all_cells) == 0) stop("❌ obj 没有细胞（ncol(obj)=0）。")
  if (is.null(colnames(mat_all))) stop("❌ mat_all 没有 colnames，无法对齐。")
  
  miss_cells <- setdiff(all_cells, colnames(mat_all))
  if (length(miss_cells) > 0) {
    if (!is.null(out_dir)) {
      write.csv(data.frame(missing_cells=miss_cells),
                file.path(out_dir, "CHECK_missing_cells_after_merge.csv"),
                row.names=FALSE)
    }
    ## 补 0 列（确保所有 obj cells 都在 mat_all 里）
    m_fill <- Matrix::Matrix(0, nrow=nrow(mat_all), ncol=length(miss_cells), sparse=TRUE)
    rownames(m_fill) <- rownames(mat_all)
    colnames(m_fill) <- miss_cells
    mat_all <- Matrix::cbind2(mat_all, m_fill)
  }
  
  j <- match(all_cells, colnames(mat_all))
  if (anyNA(j)) {
    bad <- all_cells[is.na(j)]
    if (!is.null(out_dir)) {
      write.csv(data.frame(cells_not_found=bad),
                file.path(out_dir, "CHECK_cells_not_found_in_mat_all.csv"),
                row.names=FALSE)
    }
    stop("❌ counts 合并后仍有 cells 对不上（见 CHECK_cells_not_found_in_mat_all.csv）。")
  }
  
  mat_all <- mat_all[, j, drop=FALSE]
  return(mat_all)
}

## =========================
## 3) 单脑区分析函数（EC 或 SFG）
##    输出严格落在 out_root 的子目录，不会混乱
## =========================
run_overlap_hist_count_one_region <- function(obj_fp, meta_fp, region_tag, out_root) {
  
  ## 输出目录
  out_dir <- file.path(out_root, paste0("NO1_overlap_hist_count_GSE157827style_", region_tag))
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  
  log_fp <- file.path(out_dir, paste0("NO1_log_", region_tag, ".txt"))
  sink(log_fp); on.exit(sink(), add=TRUE)
  
  cat("==== START ====\n")
  cat("Region:", region_tag, "\n")
  cat("Time:", as.character(Sys.time()), "\n")
  cat("obj_fp :", obj_fp, "\n")
  cat("meta_fp:", meta_fp, "\n")
  cat("out_dir:", out_dir, "\n")
  cat("SEED:", SEED, "\n\n")
  
  ## ------------------------------------------------------------
  ## 3.1 读对象（先不 subset，避免 Seurat v5 layers 同步 bug）
  ## ------------------------------------------------------------
  obj <- readRDS(obj_fp)
  DefaultAssay(obj) <- "RNA"
  
  ## meta：注意，obj@meta.data 的 rownames 理论上是 cell barcode
  md0 <- obj@meta.data
  if (is.null(rownames(md0)) || length(rownames(md0)) != ncol(obj)) {
    ## 双保险：强制用 colnames(obj) 当 cell_id
    md0$cell_id <- colnames(obj)
    rownames(md0) <- md0$cell_id
  } else {
    md0$cell_id <- rownames(md0)
  }
  
  cat("Object loaded. cells=", ncol(obj), " genes=", nrow(obj), "\n")
  cat("meta columns:\n"); print(colnames(md0)); cat("\n")
  
  ## ------------------------------------------------------------
  ## 3.2 从 stepP meta.csv 构建 sample->group4 映射（AD/NC->Control）
  ## ------------------------------------------------------------
  meta_raw <- read.csv(meta_fp, stringsAsFactors=FALSE)
  if (!all(c("sample","group") %in% colnames(meta_raw))) {
    stop("❌ meta_fp 缺少 sample/group 列：", meta_fp)
  }
  
  ## group=AD/NC，但图里对照固定叫 Control => NC->Control
  ctrl_alias <- c("NC","Control","CTRL","Ctr","CTR","Normal","N","control","ctrl","ctr","nc","normal")
  grp_raw <- trimws(as.character(meta_raw$group))
  meta_raw$group4 <- ifelse(grp_raw %in% ctrl_alias, "Control", grp_raw)
  
  map_by_sample <- meta_raw %>%
    transmute(sample_join = trimws(as.character(sample)),
              group4 = group4) %>%
    filter(!is.na(sample_join), sample_join != "", group4 %in% c("AD","Control")) %>%
    distinct()
  
  ## 一个 sample 只能对应一个组别
  if (any(table(map_by_sample$sample_join) > 1)) {
    bad <- names(which(table(map_by_sample$sample_join) > 1))
    write.csv(data.frame(sample=bad),
              file.path(out_dir, "CHECK_sample_maps_to_multiple_group.csv"),
              row.names=FALSE)
    stop("❌ sample 对应多个 group（见 CHECK_sample_maps_to_multiple_group.csv）。")
  }
  
  cat("Loaded sample->group map. n_samples=", nrow(map_by_sample), "\n")
  cat("Map group counts:\n"); print(table(map_by_sample$group4)); cat("\n")
  
  ## ------------------------------------------------------------
  ## 3.3 在对象 meta 中构建 join 键：优先 SampleID（EC 的 SampleID 是 EC2/EC9 等样本号）
  ## ------------------------------------------------------------
  join_key <- c("SampleID","sample","orig.ident")[c("SampleID","sample","orig.ident") %in% colnames(md0)][1]
  if (is.na(join_key)) {
    stop("❌ 对象 meta 中缺少 SampleID/sample/orig.ident，无法与 meta_fp 的 sample 对齐。")
  }
  
  md0$sample_join_in_obj <- trimws(as.character(md0[[join_key]]))
  
  ## 用 dplyr join 会丢 rownames，所以先确保 cell_id 在列里（已做）
  md1 <- dplyr::left_join(md0, map_by_sample, by=c("sample_join_in_obj"="sample_join"))
  
  match_rate <- mean(!is.na(md1$group4))
  cat("Join used: object.", join_key, " -> meta.sample\n", sep="")
  cat("Match rate =", sprintf("%.4f", match_rate), "\n\n")
  
  ## 匹配不到 group 的细胞：写检查表，但不 stop（因为 stepP 是 matched 子集）
  if (any(is.na(md1$group4))) {
    bad_vals <- sort(unique(md1$sample_join_in_obj[is.na(md1$group4)]))
    write.csv(data.frame(join_key=join_key, join_value=bad_vals),
              file.path(out_dir, "CHECK_cells_without_group_after_join.csv"),
              row.names=FALSE)
    cat("⚠ unmapped sample values saved. unmapped cells will be dropped.\n\n")
  }
  
  ## 只保留 AD/Control 的细胞
  md1 <- md1[!is.na(md1$group4) & md1$group4 %in% c("AD","Control"), , drop=FALSE]
  
  ## 必须仍包含两组
  if (!all(c("AD","Control") %in% unique(md1$group4))) {
    cat("Group counts after drop:\n"); print(table(md1$group4)); cat("\n")
    stop("❌ 丢弃未映射细胞后不再同时包含 AD 与 Control。")
  }
  
  ## ------------------------------------------------------------
  ## 3.4 donor + celltype6（按模版必需字段）
  ## ------------------------------------------------------------
  if (!("celltype6" %in% colnames(md1))) stop("❌ meta.data 中没有 celltype6。")
  md1$celltype6 <- trimws(as.character(md1$celltype6))
  
  ## donor：优先 PatientID，否则退回 sample（在 syn21788402 里 sample 通常就是 donor）
  if ("PatientID" %in% colnames(md1)) {
    md1$donor <- trimws(as.character(md1$PatientID))
  } else if ("sample" %in% colnames(md1)) {
    md1$donor <- trimws(as.character(md1$sample))
  } else if ("SampleID" %in% colnames(md1)) {
    md1$donor <- trimws(as.character(md1$SampleID))
  } else {
    stop("❌ 找不到 donor 列（PatientID/sample/SampleID）。")
  }
  
  ## 以 cell_id 恢复 rownames（用于后面 cells_keep）
  rownames(md1) <- md1$cell_id
  
  ## 得到最终要分析的细胞集合（cell barcode）
  cells_keep <- intersect(rownames(md1), colnames(obj))
  if (length(cells_keep) == 0) stop("❌ cells_keep=0。请检查 join_key 与 CHECK_cells_without_group_after_join.csv")
  
  ## 只保留这些细胞的 meta（注意顺序要与 cells_keep 一致）
  md <- md1[cells_keep, , drop=FALSE]
  
  ## ------------------------------------------------------------
  ## ★关键修复：counts 不对 obj subset，而是矩阵层面 subset
  ## ------------------------------------------------------------
  ## 1) 先合并“全量 counts”（未 subset 的 obj）
  rna_counts_all <- get_counts_matrix_allcells(obj, "RNA", log_fp=log_fp, out_dir=out_dir)
  
  ## 2) 矩阵层面取子集（完全等价于 subset 后取 counts，但更稳）
  miss_in_counts <- setdiff(cells_keep, colnames(rna_counts_all))
  if (length(miss_in_counts) > 0) {
    write.csv(data.frame(cells_missing_in_counts=miss_in_counts),
              file.path(out_dir, "CHECK_cells_missing_in_counts_after_join.csv"),
              row.names=FALSE)
    stop("❌ 有 cells_keep 不在 counts 列名中（见 CHECK_cells_missing_in_counts_after_join.csv）。")
  }
  
  rna_counts <- rna_counts_all[, cells_keep, drop=FALSE]
  
  ## counts 与 md 再次一致性检查
  stopifnot(identical(colnames(rna_counts), rownames(md)))
  
  cat("Counts(all) dim:", paste(dim(rna_counts_all), collapse=" x "), "\n")
  cat("Counts(subset) dim:", paste(dim(rna_counts), collapse=" x "), "\n\n")
  
  ## ------------------------------------------------------------
  ## 3.5 UBL3 log1p(CP10k)（严格模版）
  ## ------------------------------------------------------------
  gene_row <- if ("UBL3" %in% rownames(rna_counts)) "UBL3" else
    if ("ENSG00000122042" %in% rownames(rna_counts)) "ENSG00000122042" else NA_character_
  if (is.na(gene_row)) {
    write.csv(data.frame(head_rownames=head(rownames(rna_counts), 200)),
              file.path(out_dir, "CHECK_counts_rownames_head200.csv"),
              row.names=FALSE)
    stop("❌ counts 行名中找不到 UBL3（见 CHECK_counts_rownames_head200.csv）。")
  }
  
  lib_size <- Matrix::colSums(rna_counts)
  expr <- log1p((as.numeric(rna_counts[gene_row, , drop=TRUE]) / pmax(lib_size, 1)) * 1e4)
  
  df_all <- data.frame(
    expr      = expr,
    donor     = md$donor,
    group4    = md$group4,      ## 已经是 AD/Control（NC 已映射为 Control）
    celltype6 = md$celltype6,
    stringsAsFactors = FALSE
  )
  
  ## 只保留表达细胞
  df0 <- df_all[df_all$expr > 0, , drop=FALSE]
  saveRDS(df0, file.path(out_dir, "INTERMEDIATE_df0_exprGT0.rds"))
  write.csv(df0, file.path(out_dir, "INTERMEDIATE_df0_exprGT0.csv"), row.names=FALSE)
  
  cat("Expressed cells (expr>0):", nrow(df0), "\n\n")
  
  ## donor QC（补充材料）
  don_all <- unique(df0[, c("donor","group4")])
  qc_don <- as.data.frame(table(don_all$group4), stringsAsFactors=FALSE)
  colnames(qc_don) <- c("group4","n_donors")
  write.csv(qc_don, file.path(out_dir, "QC_donors_by_group4.csv"), row.names=FALSE)
  
  ## ------------------------------------------------------------
  ## 3.6 绘图函数：cell-level / donor-level（严格模版）
  ## ------------------------------------------------------------
  plot_one_unit <- function(unit=c("donor","cell")) {
    
    unit <- match.arg(unit)
    
    df2 <- df0 %>%
      filter(group4 %in% c(disease, "Control")) %>%
      mutate(group = ifelse(group4 == disease, disease, "Control"))
    
    ## legend donor n（expr>0 口径）
    don_pair <- unique(df2[, c("donor","group")])
    n_dis <- sum(don_pair$group == disease)
    n_ctl <- sum(don_pair$group == "Control")
    
    lab_dis <- sprintf("%s\n(n = %d)", disease,  n_dis)
    lab_ctl <- sprintf("Control\n(n = %d)",     n_ctl)
    
    df2$group_lab <- factor(ifelse(df2$group==disease, lab_dis, lab_ctl),
                            levels=c(lab_dis, lab_ctl))
    
    ## 统计输入：donor-level 用 median；cell-level 用每细胞 expr
    if (unit == "donor") {
      stat_input <- df2 %>%
        group_by(celltype6, donor, group_lab) %>%
        summarise(val = median(expr), .groups="drop")
    } else {
      stat_input <- df2 %>%
        transmute(celltype6 = celltype6, group_lab = group_lab, val = expr)
    }
    
    saveRDS(stat_input, file.path(out_dir, paste0("INTERMEDIATE_stat_input_AD_vs_Control_", unit, ".rds")))
    write.csv(stat_input, file.path(out_dir, paste0("INTERMEDIATE_stat_input_AD_vs_Control_", unit, ".csv")),
              row.names=FALSE)
    
    ## MWU + BH（每个 celltype 一次）
    stats <- stat_input %>%
      group_by(celltype6) %>%
      summarise(
        p_raw = tryCatch(wilcox.test(val ~ group_lab, exact=FALSE)$p.value,
                         error=function(e) NA_real_),
        .groups="drop"
      )
    
    stats$padj  <- p.adjust(stats$p_raw, method="BH")
    stats$label <- sprintf("Mann–Whitney U\nPadj=%.2e", stats$padj)
    stats$x <- Inf; stats$y <- Inf
    write.csv(stats, file.path(out_dir, paste0("STATS_AD_vs_Control_", unit, ".csv")), row.names=FALSE)
    
    ## 颜色：红/蓝（严格模版，不变灰）
    fill_vals <- c("red","blue")
    names(fill_vals) <- c(lab_dis, lab_ctl)
    
    ## 画图：Y=Cell count
    p <- ggplot(df2, aes(x=expr, fill=group_lab)) +
      geom_histogram(binwidth=binwidth, alpha=0.7, position="identity", colour=NA) +
      facet_wrap(~celltype6, scales="free_y") +
      scale_fill_manual(values=fill_vals, drop=FALSE) +
      geom_label(
        data=stats, inherit.aes=FALSE,
        aes(x=x, y=y, label=label),
        hjust=1.02, vjust=1.02, size=2.3,
        label.size=0, fill="white", alpha=0.7
      ) +
      labs(
        title = paste0(gene_symbol,
                       " expression per cell type (only expressed cells): ",
                       disease, " vs Control (", unit, "-level)"),
        x = paste0(gene_symbol, " log1p(CP10k)"),
        y = "Cell count",
        fill = "group_lab"
      ) +
      theme_bw() +
      theme(plot.margin = margin(10,25,10,10)) +
      coord_cartesian(clip="off")
    
    out_png <- file.path(out_dir,
                         paste0(dataset_tag, "_", region_tag, "_", gene_symbol,
                                "_OverlapHistCount_AD_vs_Control_", unit, ".png"))
    
    ragg::agg_png(out_png, width=10, height=6, units="in", res=300, background="white")
    print(p)
    dev.off()
    
    cat("✅ saved:", out_png, "\n")
  }
  
  ## 输出两张：cell-level + donor-level
  plot_one_unit("cell")
  plot_one_unit("donor")
  
  cat("\n==== sessionInfo ====\n")
  print(sessionInfo())
  cat("==== END ====\n")
  
  message("🎉 DONE region=", region_tag, " | out_dir=", out_dir)
}

## =========================
## 4) 开跑：先 EC 再 SFG（输出严格落在各自根目录）
## =========================
run_overlap_hist_count_one_region(obj_fp_EC,  meta_fp_EC,  "EC",  root_EC)
run_overlap_hist_count_one_region(obj_fp_SFG, meta_fp_SFG, "SFG", root_SFG)










#🔹第 2 章（2）：重叠直方图（Y = density）
###############################################################################
## NO2_syn21788402_EC_SFG_OverlapHistDensity_ADvsControl_cell_and_donor_level.R
##
## 【第2章目标】（严格模仿 GSE157827 模版：Overlap Hist (Density)）
##  - 两个部位：EC 与 SFG 各自独立跑（输出目录严格分开，避免混）
##  - 每个部位输出两张图（共4张）：
##      (1) cell-level：每个 cell 一个观测 → MWU → BH(FDR)
##      (2) donor-level：每 donor 先 median → MWU → BH(FDR)
##  - 只用表达细胞 expr>0
##  - UBL3：log1p((raw/lib_size)*1e4) (CP10K) ——与 GSE157827 模版一致
##  - 直方图：Y轴 = Density（y=after_stat(density)）
##  - 图中对照组名称固定显示 Control（数据里为 NC，也会映射为 Control）
##
## 【关键稳定性修复：完全继承第1章跑通方案，不改变分析逻辑】
##  - 不对 Seurat 对象 subset 后再读 layers
##  - 先合并全量 counts：rna_counts_all
##  - 再用 cells_keep 在矩阵层面取子集：rna_counts <- rna_counts_all[, cells_keep]
##  - 避免 Seurat v5 Assay5 多 layers subset 不同步导致 out-of-bounds
###############################################################################

rm(list=ls()); gc()
Sys.setenv(LANG="en")
SEED <- 20251023; set.seed(SEED)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(ragg)
})

## =========================
## 0) 全局参数（与你模版一致）
## =========================
dataset_tag <- "syn21788402"
gene_symbol <- "UBL3"
binwidth    <- 0.2
disease     <- "AD"

## =========================
## 1) 两个部位路径（严格不弄乱）
## =========================
root_EC  <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/results"
root_SFG <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/results_SFG"

obj_fp_EC  <- file.path(root_EC,  "stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds")
obj_fp_SFG <- file.path(root_SFG, "stepH_syn21788402_SFG_obj_celltype6.rds")
stopifnot(file.exists(obj_fp_EC))
stopifnot(file.exists(obj_fp_SFG))

meta_fp_EC  <- file.path(root_EC,  "stepP_syn21788402_matched_cells_meta.csv")
meta_fp_SFG <- file.path(root_SFG, "stepP_syn21788402_SFG_matched_cells_meta.csv")
stopifnot(file.exists(meta_fp_EC))
stopifnot(file.exists(meta_fp_SFG))

## =========================
## 2) 合并 counts layers（Seurat v5 多 layers 兼容）
##    说明：这段与第1章一致，只负责把 counts.* layers 合并成 gene×cell
## =========================
get_counts_matrix_allcells <- function(obj, assay="RNA", log_fp=NULL, out_dir=NULL) {
  
  a <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(a), error=function(e) character(0))
  counts_layers <- layers[grepl("^counts", layers)]
  
  ## 无 layers：退回传统 counts slot
  if (length(counts_layers) == 0) {
    m2 <- tryCatch(Seurat::GetAssayData(obj, assay=assay, slot="counts"), error=function(e) NULL)
    if (!is.null(m2)) return(m2)
    stop("❌ 无法获取 counts：没有 counts.* layers，GetAssayData(counts) 也失败。")
  }
  
  if (!is.null(log_fp)) {
    sink(log_fp, append=TRUE)
    cat("Detected counts layers:\n"); print(counts_layers)
    sink()
  }
  
  mats <- list()
  for (ly in counts_layers) {
    m <- tryCatch(SeuratObject::LayerData(a, layer=ly), error=function(e) NULL)
    if (is.null(m)) next
    if (!is.null(dim(m)) && length(dim(m))==2) mats[[ly]] <- m
  }
  if (length(mats)==0) stop("❌ counts layers 存在，但读取 LayerData 全部失败。")
  
  ## 基因对齐
  ref_genes <- rownames(mats[[1]])
  for (k in names(mats)) {
    if (!identical(rownames(mats[[k]]), ref_genes)) {
      m0 <- mats[[k]]
      m_aligned <- Matrix::Matrix(0, nrow=length(ref_genes), ncol=ncol(m0), sparse=TRUE)
      rownames(m_aligned) <- ref_genes
      colnames(m_aligned) <- colnames(m0)
      common <- intersect(ref_genes, rownames(m0))
      if (length(common) > 0) m_aligned[common, ] <- m0[common, , drop=FALSE]
      mats[[k]] <- m_aligned
    }
  }
  
  mat_all <- if (length(mats)==1) mats[[1]] else Reduce(Matrix::cbind2, mats)
  
  ## 去除重复 cell
  if (!is.null(colnames(mat_all))) {
    dup <- duplicated(colnames(mat_all))
    if (any(dup)) mat_all <- mat_all[, !dup, drop=FALSE]
  }
  
  ## 对齐到 obj 当前细胞顺序
  all_cells <- colnames(obj)
  if (is.null(colnames(mat_all))) stop("❌ mat_all 没有 colnames，无法对齐。")
  
  miss_cells <- setdiff(all_cells, colnames(mat_all))
  if (length(miss_cells) > 0) {
    if (!is.null(out_dir)) {
      write.csv(data.frame(missing_cells=miss_cells),
                file.path(out_dir, "CHECK_missing_cells_after_merge.csv"),
                row.names=FALSE)
    }
    m_fill <- Matrix::Matrix(0, nrow=nrow(mat_all), ncol=length(miss_cells), sparse=TRUE)
    rownames(m_fill) <- rownames(mat_all)
    colnames(m_fill) <- miss_cells
    mat_all <- Matrix::cbind2(mat_all, m_fill)
  }
  
  j <- match(all_cells, colnames(mat_all))
  if (anyNA(j)) {
    bad <- all_cells[is.na(j)]
    if (!is.null(out_dir)) {
      write.csv(data.frame(cells_not_found=bad),
                file.path(out_dir, "CHECK_cells_not_found_in_mat_all.csv"),
                row.names=FALSE)
    }
    stop("❌ counts 合并后仍有 cells 对不上（见 CHECK_cells_not_found_in_mat_all.csv）。")
  }
  
  mat_all <- mat_all[, j, drop=FALSE]
  return(mat_all)
}

## =========================
## 3) 单部位函数：输出 Density overlap hist（cell+donor）
## =========================
run_overlap_hist_density_one_region <- function(obj_fp, meta_fp, region_tag, out_root) {
  
  ## 输出目录（第2章）
  out_dir <- file.path(out_root, paste0("NO2_overlap_hist_density_GSE157827style_", region_tag))
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  
  log_fp <- file.path(out_dir, paste0("NO2_log_", region_tag, ".txt"))
  sink(log_fp); on.exit(sink(), add=TRUE)
  
  cat("==== START ====\n")
  cat("Region:", region_tag, "\n")
  cat("Time:", as.character(Sys.time()), "\n")
  cat("obj_fp :", obj_fp, "\n")
  cat("meta_fp:", meta_fp, "\n")
  cat("out_dir:", out_dir, "\n")
  cat("SEED:", SEED, "\n\n")
  
  ## ------------------------------------------------------------
  ## 3.1 读对象（不 subset）
  ## ------------------------------------------------------------
  obj <- readRDS(obj_fp)
  DefaultAssay(obj) <- "RNA"
  
  md0 <- obj@meta.data
  ## 保证 cell_id 与 cell barcode 一一对应
  md0$cell_id <- colnames(obj)
  rownames(md0) <- md0$cell_id
  
  if (!("celltype6" %in% colnames(md0))) stop("❌ meta.data 中没有 celltype6。")
  md0$celltype6 <- trimws(as.character(md0$celltype6))
  
  ## donor：优先 PatientID，否则退回 sample/SampleID
  if ("PatientID" %in% colnames(md0)) {
    md0$donor <- trimws(as.character(md0$PatientID))
    donor_source <- "PatientID"
  } else if ("sample" %in% colnames(md0)) {
    md0$donor <- trimws(as.character(md0$sample))
    donor_source <- "sample"
  } else if ("SampleID" %in% colnames(md0)) {
    md0$donor <- trimws(as.character(md0$SampleID))
    donor_source <- "SampleID"
  } else {
    stop("❌ 找不到 donor 列（PatientID/sample/SampleID）。")
  }
  cat("donor_source =", donor_source, "\n\n")
  
  ## ------------------------------------------------------------
  ## 3.2 从 stepP meta.csv 构建 sample->group4（AD/NC->Control）
  ## ------------------------------------------------------------
  meta_raw <- read.csv(meta_fp, stringsAsFactors=FALSE)
  if (!all(c("sample","group") %in% colnames(meta_raw))) stop("❌ meta_fp 缺少 sample/group。")
  
  ctrl_alias <- c("NC","Control","CTRL","Ctr","CTR","Normal","N","control","ctrl","ctr","nc","normal")
  grp_raw <- trimws(as.character(meta_raw$group))
  meta_raw$group4 <- ifelse(grp_raw %in% ctrl_alias, "Control", grp_raw)
  
  map_by_sample <- meta_raw %>%
    transmute(sample_join = trimws(as.character(sample)),
              group4 = group4) %>%
    filter(!is.na(sample_join), sample_join != "", group4 %in% c("AD","Control")) %>%
    distinct()
  
  if (any(table(map_by_sample$sample_join) > 1)) {
    bad <- names(which(table(map_by_sample$sample_join) > 1))
    write.csv(data.frame(sample=bad),
              file.path(out_dir, "CHECK_sample_maps_to_multiple_group.csv"),
              row.names=FALSE)
    stop("❌ sample 对应多个 group（见 CHECK_sample_maps_to_multiple_group.csv）。")
  }
  
  ## ------------------------------------------------------------
  ## 3.3 对象侧 join 键：优先 SampleID（样本号）对 meta 的 sample
  ## ------------------------------------------------------------
  join_key <- c("SampleID","sample","orig.ident")[c("SampleID","sample","orig.ident") %in% colnames(md0)][1]
  if (is.na(join_key)) stop("❌ 对象 meta 缺少 SampleID/sample/orig.ident，无法对齐 meta.sample。")
  
  md0$sample_join_in_obj <- trimws(as.character(md0[[join_key]]))
  
  ## dplyr join 会丢 rownames，但我们有 cell_id 列可恢复
  md1 <- dplyr::left_join(md0, map_by_sample, by=c("sample_join_in_obj"="sample_join"))
  
  match_rate <- mean(!is.na(md1$group4))
  cat("Join key =", join_key, " | match_rate =", sprintf("%.4f", match_rate), "\n\n")
  
  ## 记录 join 不上的 sample 值（不 stop，丢弃继续）
  if (any(is.na(md1$group4))) {
    bad_vals <- sort(unique(md1$sample_join_in_obj[is.na(md1$group4)]))
    write.csv(data.frame(join_key=join_key, join_value=bad_vals),
              file.path(out_dir, "CHECK_cells_without_group_after_join.csv"),
              row.names=FALSE)
    cat("⚠ unmapped sample values saved; unmapped cells will be dropped.\n\n")
  }
  
  ## 保留 AD/Control
  md1 <- md1[!is.na(md1$group4) & md1$group4 %in% c("AD","Control"), , drop=FALSE]
  rownames(md1) <- md1$cell_id
  
  if (!all(c("AD","Control") %in% unique(md1$group4))) {
    cat("Group counts after drop:\n"); print(table(md1$group4)); cat("\n")
    stop("❌ 丢弃未映射细胞后不再同时包含 AD 与 Control。")
  }
  
  ## 最终 cells_keep
  cells_keep <- intersect(rownames(md1), colnames(obj))
  if (length(cells_keep)==0) stop("❌ cells_keep=0。请检查 join_key 与 CHECK_cells_without_group_after_join.csv")
  
  md <- md1[cells_keep, , drop=FALSE]
  
  ## ------------------------------------------------------------
  ## ★关键稳定方案：先合并全量 counts，再矩阵层面 subset
  ## ------------------------------------------------------------
  rna_counts_all <- get_counts_matrix_allcells(obj, "RNA", log_fp=log_fp, out_dir=out_dir)
  
  miss_in_counts <- setdiff(cells_keep, colnames(rna_counts_all))
  if (length(miss_in_counts) > 0) {
    write.csv(data.frame(cells_missing_in_counts=miss_in_counts),
              file.path(out_dir, "CHECK_cells_missing_in_counts_after_join.csv"),
              row.names=FALSE)
    stop("❌ 有 cells_keep 不在 counts 列名中（见 CHECK_cells_missing_in_counts_after_join.csv）。")
  }
  
  rna_counts <- rna_counts_all[, cells_keep, drop=FALSE]
  stopifnot(identical(colnames(rna_counts), rownames(md)))
  
  ## ------------------------------------------------------------
  ## 3.4 UBL3 log1p(CP10k)（严格模版）
  ## ------------------------------------------------------------
  gene_row <- if ("UBL3" %in% rownames(rna_counts)) "UBL3" else
    if ("ENSG00000122042" %in% rownames(rna_counts)) "ENSG00000122042" else NA_character_
  if (is.na(gene_row)) stop("❌ counts 行名中找不到 UBL3（UBL3 或 ENSG00000122042）。")
  
  lib_size <- Matrix::colSums(rna_counts)
  expr <- log1p((as.numeric(rna_counts[gene_row, , drop=TRUE]) / pmax(lib_size,1)) * 1e4)
  
  df_all <- data.frame(
    expr      = expr,
    donor     = md$donor,
    group4    = md$group4,      ## AD/Control（NC已映射）
    celltype6 = md$celltype6,
    stringsAsFactors = FALSE
  )
  
  df0 <- df_all[df_all$expr > 0, , drop=FALSE]  ## only expressed cells
  saveRDS(df0, file.path(out_dir, "INTERMEDIATE_df0_exprGT0.rds"))
  write.csv(df0, file.path(out_dir, "INTERMEDIATE_df0_exprGT0.csv"), row.names=FALSE)
  
  ## donor QC（补充材料）
  don_all <- unique(df0[, c("donor","group4")])
  qc_don <- as.data.frame(table(don_all$group4), stringsAsFactors=FALSE)
  colnames(qc_don) <- c("group4","n_donors")
  write.csv(qc_don, file.path(out_dir, "QC_donors_by_group4.csv"), row.names=FALSE)
  
  ## ------------------------------------------------------------
  ## 3.5 画图函数：Density overlap hist（严格复刻 syn520 NO4_03 风格）
  ## ------------------------------------------------------------
  plot_one_unit <- function(unit=c("donor","cell")) {
    
    unit <- match.arg(unit)
    
    df2 <- df0 %>%
      filter(group4 %in% c(disease,"Control")) %>%
      mutate(group = ifelse(group4==disease, disease, "Control"))
    
    ## legend：donor n（expr>0口径）
    don_pair <- unique(df2[, c("donor","group")])
    n_dis <- sum(don_pair$group == disease)
    n_ctl <- sum(don_pair$group == "Control")
    
    lab_dis <- sprintf("%s\n(n = %d)", disease,  n_dis)
    lab_ctl <- sprintf("Control\n(n = %d)",     n_ctl)
    
    ## ★关键：levels 必须和 fill_vals 的 names 一致，否则颜色会变灰
    df2$group_lab <- factor(ifelse(df2$group==disease, lab_dis, lab_ctl),
                            levels=c(lab_dis, lab_ctl))
    
    ## 统计输入
    if (unit=="donor") {
      stat_input <- df2 %>%
        group_by(celltype6, donor, group_lab) %>%
        summarise(val = median(expr), .groups="drop")
    } else {
      stat_input <- df2 %>%
        transmute(celltype6=celltype6, group_lab=group_lab, val=expr)
    }
    
    ## 保存补充材料
    saveRDS(stat_input, file.path(out_dir, paste0("INTERMEDIATE_stat_input_AD_vs_Control_", unit, ".rds")))
    write.csv(stat_input, file.path(out_dir, paste0("INTERMEDIATE_stat_input_AD_vs_Control_", unit, ".csv")),
              row.names=FALSE)
    
    ## 每面板样本量检查（补充材料/核查用）
    n_by_celltype <- stat_input %>%
      group_by(celltype6, group_lab) %>%
      summarise(n=n(), .groups="drop")
    write.csv(n_by_celltype, file.path(out_dir, paste0("CHECK_n_by_celltype_AD_vs_Control_", unit, ".csv")),
              row.names=FALSE)
    
    ## MWU + BH
    stats <- stat_input %>%
      group_by(celltype6) %>%
      summarise(
        p_raw = tryCatch(wilcox.test(val ~ group_lab, exact=FALSE)$p.value,
                         error=function(e) NA_real_),
        .groups="drop"
      )
    stats$padj <- p.adjust(stats$p_raw, method="BH")
    stats$label <- sprintf("Mann–Whitney U\nPadj=%.2e", stats$padj)
    stats$x <- Inf; stats$y <- Inf
    write.csv(stats, file.path(out_dir, paste0("STATS_AD_vs_Control_", unit, ".csv")), row.names=FALSE)
    
    ## 颜色：红/蓝
    fill_vals <- c("red","blue")
    names(fill_vals) <- c(lab_dis, lab_ctl)
    
    ## ★Y轴=Density：y=after_stat(density)
    p <- ggplot(df2, aes(x=expr, y=after_stat(density), fill=group_lab)) +
      geom_histogram(binwidth=binwidth, alpha=0.7, position="identity", colour=NA) +
      facet_wrap(~celltype6, scales="free_y") +
      scale_fill_manual(values=fill_vals, drop=FALSE) +
      geom_label(
        data=stats, inherit.aes=FALSE,
        aes(x=x, y=y, label=label),
        hjust=1.02, vjust=1.02, size=2.3,
        label.size=0, fill="white", alpha=0.7
      ) +
      labs(
        title = paste0(gene_symbol,
                       " expression per cell type (only expressed cells): ",
                       disease, " vs Control (", unit, "-level)"),
        x = paste0(gene_symbol, " log1p(CP10k)"),
        y = "Density",
        fill = "group_lab"
      ) +
      theme_bw() +
      theme(plot.margin = margin(10,25,10,10)) +
      coord_cartesian(clip="off")
    
    ## 输出文件（带 region_tag，避免 EC/SFG 覆盖）
    out_png <- file.path(out_dir,
                         paste0(dataset_tag,"_",region_tag,"_",gene_symbol,
                                "_OverlapHistDensity_AD_vs_Control_", unit, ".png"))
    
    ragg::agg_png(out_png, width=10, height=6, units="in", res=300, background="white")
    print(p)
    dev.off()
    
    cat("✅ saved:", out_png, "\n")
    sink(log_fp, append=TRUE)
    cat("Saved figure:", out_png, "\n")
    sink()
  }
  
  ## 输出两张：cell-level + donor-level
  plot_one_unit("cell")
  plot_one_unit("donor")
  
  cat("\n==== sessionInfo ====\n")
  print(sessionInfo())
  cat("==== END ====\n")
  
  message("🎉 DONE region=", region_tag, " | out_dir=", out_dir)
}

## =========================
## 4) 开跑：EC 与 SFG（分部位分别保存）
## =========================
run_overlap_hist_density_one_region(obj_fp_EC,  meta_fp_EC,  "EC",  root_EC)
run_overlap_hist_density_one_region(obj_fp_SFG, meta_fp_SFG, "SFG", root_SFG)





#箱线图
#6个细胞类型的
###############################################################################
## NO3_syn21788402_EC_SFG_UBL3_Boxplots_DESeq2_byDonor.R
##
## 【第3章目标】（严格模仿 GSE157827 pseudo-bulk DESeq2 by donor 模版）
##  - 对 EC 与 SFG 两个脑区分别做：
##    1) donor 为统计单位：先 pseudo-bulk 到 celltype6 × donor（counts 累加）
##    2) 每个 celltype6 跑 DESeq2（AD vs Control）
##    3) 输出：
##       a) 每个 celltype 一张 UBL3 donor-level 箱线图（四分位数）
##       b) 每个 celltype 一份全基因 DEG 表（AD vs Control）
##       c) 可选：6-panel 总图（2×3）
##
## 【分组规则】原始 meta 里是 AD/NC
##  - 内部统一：NC/CTRL/... -> Control
##  - 图上对照固定显示 Control（满足你的要求）
##
## 【两脑区输出目录严格分开，不弄乱】
##  - EC : .../results/NO3_UBL3_Boxplots_DESeq2_byDonor_EC/
##  - SFG: .../results_SFG/NO3_UBL3_Boxplots_DESeq2_byDonor_SFG/
##
## 【必须核对的中间结果（脚本会自动输出）】
##  - CHECK_sample_maps_to_multiple_group.csv（若 sample 对应多个组会 stop）
##  - CHECK_cells_without_group_after_join.csv（哪些样本号没法匹配，会被丢弃）
##  - QC_donor_counts_by_group.csv
##  - QC_donor_counts_by_celltype6_by_group.csv
##  - QC_cells_by_celltype6_by_group.csv
##  - INTERMEDIATE_pseudobulk_coldata_celltype6_donor_group.csv
##  - INTERMEDIATE_pseudobulk_matrix_celltype6_donor.rds
###############################################################################

rm(list = ls()); gc()
Sys.setenv(LANG = "en")
SEED <- 20251023; set.seed(SEED)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ragg)
})

## =========================
## 0) 全局参数
## =========================
dataset_tag <- "syn21788402"
gene_symbol <- "UBL3"

## 颜色：与 GSE157827 模版一致
pal2 <- c(AD="#D24B40", Control="#2C7FB8")

## 是否输出每脑区的 6-panel 总图（你如果只要12张单图，可改 FALSE）
MAKE_6PANEL <- TRUE

## =========================
## 1) 两脑区路径（严格不弄乱）
## =========================
root_EC  <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/results"
root_SFG <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/results_SFG"

obj_fp_EC  <- file.path(root_EC,  "stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds")
obj_fp_SFG <- file.path(root_SFG, "stepH_syn21788402_SFG_obj_celltype6.rds")
stopifnot(file.exists(obj_fp_EC), file.exists(obj_fp_SFG))

meta_fp_EC  <- file.path(root_EC,  "stepP_syn21788402_matched_cells_meta.csv")
meta_fp_SFG <- file.path(root_SFG, "stepP_syn21788402_SFG_matched_cells_meta.csv")
stopifnot(file.exists(meta_fp_EC), file.exists(meta_fp_SFG))

## =========================
## 2) counts 多 layers 合并函数（与第1章/第2章同款稳定方案）
## =========================
get_counts_matrix_allcells <- function(obj, assay="RNA", log_fp=NULL, out_dir=NULL) {
  
  a <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(a), error=function(e) character(0))
  counts_layers <- layers[grepl("^counts", layers)]
  
  if (length(counts_layers) == 0) {
    m2 <- tryCatch(Seurat::GetAssayData(obj, assay=assay, slot="counts"), error=function(e) NULL)
    if (!is.null(m2)) return(m2)
    stop("❌ 无法获取 counts：没有 counts.* layers，GetAssayData(counts) 也失败。")
  }
  
  if (!is.null(log_fp)) {
    sink(log_fp, append=TRUE)
    cat("Detected counts layers:\n"); print(counts_layers)
    sink()
  }
  
  mats <- list()
  for (ly in counts_layers) {
    m <- tryCatch(SeuratObject::LayerData(a, layer=ly), error=function(e) NULL)
    if (is.null(m)) next
    if (!is.null(dim(m)) && length(dim(m))==2) mats[[ly]] <- m
  }
  if (length(mats)==0) stop("❌ counts layers 存在，但读取 LayerData 全部失败。")
  
  ref_genes <- rownames(mats[[1]])
  for (k in names(mats)) {
    if (!identical(rownames(mats[[k]]), ref_genes)) {
      m0 <- mats[[k]]
      m_aligned <- Matrix::Matrix(0, nrow=length(ref_genes), ncol=ncol(m0), sparse=TRUE)
      rownames(m_aligned) <- ref_genes
      colnames(m_aligned) <- colnames(m0)
      common <- intersect(ref_genes, rownames(m0))
      if (length(common) > 0) m_aligned[common, ] <- m0[common, , drop=FALSE]
      mats[[k]] <- m_aligned
    }
  }
  
  mat_all <- if (length(mats)==1) mats[[1]] else Reduce(Matrix::cbind2, mats)
  
  ## 去重 cell
  if (!is.null(colnames(mat_all))) {
    dup <- duplicated(colnames(mat_all))
    if (any(dup)) mat_all <- mat_all[, !dup, drop=FALSE]
  }
  
  ## 对齐 obj cell 顺序（match 数值索引）
  all_cells <- colnames(obj)
  miss_cells <- setdiff(all_cells, colnames(mat_all))
  if (length(miss_cells) > 0) {
    if (!is.null(out_dir)) {
      write.csv(data.frame(missing_cells=miss_cells),
                file.path(out_dir, "CHECK_missing_cells_after_merge.csv"),
                row.names=FALSE)
    }
    m_fill <- Matrix::Matrix(0, nrow=nrow(mat_all), ncol=length(miss_cells), sparse=TRUE)
    rownames(m_fill) <- rownames(mat_all)
    colnames(m_fill) <- miss_cells
    mat_all <- Matrix::cbind2(mat_all, m_fill)
  }
  
  j <- match(all_cells, colnames(mat_all))
  if (anyNA(j)) {
    bad <- all_cells[is.na(j)]
    if (!is.null(out_dir)) {
      write.csv(data.frame(cells_not_found=bad),
                file.path(out_dir, "CHECK_cells_not_found_in_mat_all.csv"),
                row.names=FALSE)
    }
    stop("❌ counts 合并后仍有 cells 对不上（见 CHECK_cells_not_found_in_mat_all.csv）。")
  }
  
  mat_all <- mat_all[, j, drop=FALSE]
  return(mat_all)
}

## =========================
## 3) 单脑区主函数：pseudo-bulk + DESeq2 + UBL3 箱线图
## =========================
run_boxplots_deseq2_one_region <- function(obj_fp, meta_fp, region_tag, out_root) {
  
  out_dir <- file.path(out_root, paste0("NO3_UBL3_Boxplots_DESeq2_byDonor_", region_tag))
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  
  log_fp <- file.path(out_dir, paste0("NO3_log_", region_tag, ".txt"))
  sink(log_fp); on.exit(sink(), add=TRUE)
  
  cat("==== START ====\n")
  cat("Region:", region_tag, "\n")
  cat("Time:", as.character(Sys.time()), "\n")
  cat("obj_fp :", obj_fp, "\n")
  cat("meta_fp:", meta_fp, "\n")
  cat("out_dir:", out_dir, "\n")
  cat("SEED:", SEED, "\n\n")
  
  ## ------------------------------------------------------------
  ## 3.1 读对象（不 subset，避免 v5 layers bug）
  ## ------------------------------------------------------------
  obj <- readRDS(obj_fp)
  DefaultAssay(obj) <- "RNA"
  
  md0 <- obj@meta.data
  md0$cell_id <- colnames(obj)      ## 统一使用 cell barcode
  rownames(md0) <- md0$cell_id
  
  if (!("celltype6" %in% colnames(md0))) stop("❌ meta.data 中没有 celltype6。")
  md0$celltype6 <- trimws(as.character(md0$celltype6))
  
  ## donor：优先 PatientID，否则退回 sample/SampleID
  if ("PatientID" %in% colnames(md0)) {
    md0$donor <- trimws(as.character(md0$PatientID))
    donor_source <- "PatientID"
  } else if ("sample" %in% colnames(md0)) {
    md0$donor <- trimws(as.character(md0$sample))
    donor_source <- "sample"
  } else if ("SampleID" %in% colnames(md0)) {
    md0$donor <- trimws(as.character(md0$SampleID))
    donor_source <- "SampleID"
  } else {
    stop("❌ 找不到 donor 列（PatientID/sample/SampleID）。")
  }
  cat("donor_source =", donor_source, "\n\n")
  
  ## ------------------------------------------------------------
  ## 3.2 从 stepP meta.csv 构建 sample->group4（AD/NC->Control）
  ## ------------------------------------------------------------
  meta_raw <- read.csv(meta_fp, stringsAsFactors=FALSE)
  if (!all(c("sample","group") %in% colnames(meta_raw))) stop("❌ meta_fp 缺少 sample/group。")
  
  ctrl_alias <- c("NC","Control","CTRL","Ctr","CTR","Normal","N","control","ctrl","ctr","nc","normal")
  grp_raw <- trimws(as.character(meta_raw$group))
  meta_raw$group4 <- ifelse(grp_raw %in% ctrl_alias, "Control", grp_raw)
  
  map_by_sample <- meta_raw %>%
    transmute(sample_join = trimws(as.character(sample)),
              group4 = group4) %>%
    filter(!is.na(sample_join), sample_join != "", group4 %in% c("AD","Control")) %>%
    distinct()
  
  ## sample 一对一校验（否则 stop）
  if (any(table(map_by_sample$sample_join) > 1)) {
    bad <- names(which(table(map_by_sample$sample_join) > 1))
    write.csv(data.frame(sample=bad),
              file.path(out_dir, "CHECK_sample_maps_to_multiple_group.csv"),
              row.names=FALSE)
    stop("❌ sample 对应多个 group（见 CHECK_sample_maps_to_multiple_group.csv）。")
  }
  
  ## ------------------------------------------------------------
  ## 3.3 join：对象侧用 SampleID（样本号）优先对 meta.sample
  ## ------------------------------------------------------------
  join_key <- c("SampleID","sample","orig.ident")[c("SampleID","sample","orig.ident") %in% colnames(md0)][1]
  if (is.na(join_key)) stop("❌ 对象 meta 缺少 SampleID/sample/orig.ident，无法对齐 meta.sample。")
  
  md0$sample_join_in_obj <- trimws(as.character(md0[[join_key]]))
  md1 <- dplyr::left_join(md0, map_by_sample, by=c("sample_join_in_obj"="sample_join"))
  
  ## join 不上的记录下来，但不 stop（matched 子集）
  if (any(is.na(md1$group4))) {
    bad_vals <- sort(unique(md1$sample_join_in_obj[is.na(md1$group4)]))
    write.csv(data.frame(join_key=join_key, join_value=bad_vals),
              file.path(out_dir, "CHECK_cells_without_group_after_join.csv"),
              row.names=FALSE)
    cat("⚠ unmapped sample values saved; unmapped cells will be dropped.\n\n")
  }
  
  ## 只保留 AD/Control
  md1 <- md1[!is.na(md1$group4) & md1$group4 %in% c("AD","Control"), , drop=FALSE]
  rownames(md1) <- md1$cell_id
  
  if (!all(c("AD","Control") %in% unique(md1$group4))) {
    cat("Group counts after drop:\n"); print(table(md1$group4)); cat("\n")
    stop("❌ 丢弃未映射细胞后不再同时包含 AD 与 Control。")
  }
  
  ## cells_keep（必须是 cell barcode）
  cells_keep <- intersect(rownames(md1), colnames(obj))
  if (length(cells_keep)==0) stop("❌ cells_keep=0。请检查 join_key 与 CHECK_cells_without_group_after_join.csv")
  md <- md1[cells_keep, , drop=FALSE]
  
  ## ------------------------------------------------------------
  ## 3.4 QC：donor->group 一对一校验（必须）
  ## ------------------------------------------------------------
  donor_map <- unique(md[, c("donor","group4")])
  if (any(table(donor_map$donor) > 1)) {
    bad <- donor_map[donor_map$donor %in% names(which(table(donor_map$donor) > 1)), ]
    write.csv(bad, file.path(out_dir, "CHECK_donor_maps_to_multiple_group.csv"), row.names=FALSE)
    stop("❌ donor 对应多个 group，已输出 CHECK_donor_maps_to_multiple_group.csv")
  }
  
  write.csv(as.data.frame(table(donor_map$group4)),
            file.path(out_dir, "QC_donor_counts_by_group.csv"),
            row.names=FALSE)
  
  tmp_donor_ct <- unique(md[, c("donor","celltype6","group4")])
  write.csv(as.data.frame(table(tmp_donor_ct$celltype6, tmp_donor_ct$group4)),
            file.path(out_dir, "QC_donor_counts_by_celltype6_by_group.csv"),
            row.names=FALSE)
  
  write.csv(as.data.frame(table(md$celltype6, md$group4)),
            file.path(out_dir, "QC_cells_by_celltype6_by_group.csv"),
            row.names=FALSE)
  
  ## ------------------------------------------------------------
  ## 3.5 counts：先全量合并，再矩阵层面 subset（稳定方案）
  ## ------------------------------------------------------------
  cnt_all <- get_counts_matrix_allcells(obj, "RNA", log_fp=log_fp, out_dir=out_dir)
  
  ## 子集到 cells_keep，并保证顺序一致
  cnt <- cnt_all[, cells_keep, drop=FALSE]
  stopifnot(identical(colnames(cnt), rownames(md)))
  
  ## ------------------------------------------------------------
  ## 3.6 伪 bulk：celltype6 × donor（稀疏聚合）
  ## ------------------------------------------------------------
  pb_key <- paste(md$celltype6, md$donor, sep="__")
  grp <- factor(pb_key, levels=unique(pb_key))
  
  ## M：cell × (celltype6__donor) 的指示矩阵
  M <- Matrix::sparseMatrix(
    i = seq_along(grp),
    j = as.integer(grp),
    x = 1,
    dims = c(length(grp), length(levels(grp))),
    dimnames = list(rownames(md), levels(grp))
  )
  
  ## gene×cell 乘 cell×pb -> gene×pb
  pb <- cnt %*% M
  pb <- as(pb, "dgCMatrix")
  
  pb_meta <- data.frame(
    key       = colnames(pb),
    celltype6 = sub("__.*","", colnames(pb)),
    donor     = sub(".*__","", colnames(pb)),
    group     = donor_map$group4[match(sub(".*__","", colnames(pb)), donor_map$donor)],
    stringsAsFactors = FALSE
  )
  
  saveRDS(pb, file.path(out_dir, "INTERMEDIATE_pseudobulk_matrix_celltype6_donor.rds"))
  write.csv(pb_meta, file.path(out_dir, "INTERMEDIATE_pseudobulk_coldata_celltype6_donor_group.csv"),
            row.names=FALSE)
  
  ## ------------------------------------------------------------
  ## 3.7 找 UBL3 行名（symbol 或 ENSG）
  ## ------------------------------------------------------------
  gene_row <- if ("UBL3" %in% rownames(pb)) "UBL3" else
    if ("ENSG00000122042" %in% rownames(pb)) "ENSG00000122042" else NA_character_
  if (is.na(gene_row)) {
    write.csv(data.frame(head_rownames=head(rownames(pb), 200)),
              file.path(out_dir, "CHECK_pseudobulk_rownames_head200.csv"),
              row.names=FALSE)
    stop("❌ pseudo-bulk 矩阵行名里找不到 UBL3（见 CHECK_pseudobulk_rownames_head200.csv）。")
  }
  write.csv(data.frame(gene="UBL3", gene_row=gene_row),
            file.path(out_dir, "CHECK_geneSymbol_to_rowname.csv"),
            row.names=FALSE)
  
  ## ------------------------------------------------------------
  ## 3.8 每个 celltype：DESeq2 + UBL3 箱线图 + 全基因 DEG
  ## ------------------------------------------------------------
  plots <- list()
  
  for (ct in sort(unique(pb_meta$celltype6))) {
    
    cols <- pb_meta$key[pb_meta$celltype6 == ct]
    if (length(cols) < 4) {
      cat("⚠ Skip celltype (too few pseudo-bulk columns):", ct, " n=", length(cols), "\n")
      next
    }
    
    ## DESeq2 colData：Control 为 reference
    coldata <- data.frame(
      group = factor(pb_meta$group[pb_meta$celltype6 == ct], levels=c("Control","AD")),
      row.names = cols
    )
    
    ## counts：整数
    y <- round(as.matrix(pb[, cols, drop=FALSE]))
    storage.mode(y) <- "integer"
    
    dds <- DESeqDataSetFromMatrix(y, coldata, design=~group)
    dds <- DESeq(dds, quiet=TRUE)
    
    res <- results(dds, contrast=c("group","AD","Control"))
    
    ## 输出全基因 DEG（补充材料）
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df <- res_df[, c("gene","log2FoldChange","padj","pvalue","baseMean")]
    write.csv(res_df,
              file.path(out_dir, paste0("DEG_", dataset_tag, "_", region_tag, "_", ct, "_AD_vs_Control_allgenes.csv")),
              row.names=FALSE)
    
    ## UBL3 的 log2FC/padj（用于 subtitle）
    u_log2fc <- as.numeric(res[gene_row, "log2FoldChange"])
    u_padj   <- as.numeric(res[gene_row, "padj"])
    subtitle <- sprintf("AD vs Control : log2FC=%.3f, padj=%s",
                        u_log2fc,
                        ifelse(is.na(u_padj), "NA", format(u_padj, digits=3, scientific=TRUE)))
    
    ## donor-level normalized counts
    norm <- counts(dds, normalized=TRUE)
    
    dfp <- data.frame(
      donor = pb_meta$donor[pb_meta$celltype6==ct],
      group = factor(pb_meta$group[pb_meta$celltype6==ct], levels=c("AD","Control")),
      value = as.numeric(norm[gene_row, cols]),
      stringsAsFactors = FALSE
    )
    write.csv(dfp,
              file.path(out_dir, paste0("INTERMEDIATE_plotdata_", dataset_tag, "_", region_tag, "_UBL3_", ct, ".csv")),
              row.names=FALSE)
    
    ## 画箱线图（四分位数）+ donor 点（模板风格）
    p <- ggplot(dfp, aes(x=group, y=value, fill=group)) +
      geom_boxplot(width=0.55, outlier.shape=NA, alpha=0.95, colour="grey15") +
      geom_point(position=position_jitter(width=0.10),
                 size=2.6, alpha=0.9, shape=21, stroke=0.4, colour="grey10") +
      scale_fill_manual(values=pal2, drop=FALSE) +
      labs(
        title = paste0("UBL3 in ", ct, " (pseudo-bulk per donor)"),
        subtitle = subtitle,
        y = "Normalized counts",
        x = NULL
      ) +
      theme_bw(base_size=13) +
      theme(
        plot.title.position = "plot",
        plot.title    = element_text(face="bold", size=22, margin=margin(b=4)),
        plot.subtitle = element_text(size=12, colour="grey15", margin=margin(b=8)),
        legend.position = "none"
      )
    
    out_png <- file.path(out_dir, paste0("Plot_UBL3_", dataset_tag, "_", region_tag, "_", ct, "_AD_vs_Control_DESeq2_pretty.png"))
    ragg::agg_png(out_png, width=7.6, height=5.6, units="in", res=300, background="white")
    print(p)
    dev.off()
    
    plots[[ct]] <- p
  }
  
  ## ------------------------------------------------------------
  ## 3.9 6-panel 总图（可选）
  ## ------------------------------------------------------------
  if (MAKE_6PANEL) {
    ## 用当前 celltype6 的排序（如果你想强制固定顺序，可以在这里写死 order）
    panel_order <- sort(intersect(names(plots), unique(pb_meta$celltype6)))
    if (length(panel_order) == 6) {
      p_all <- patchwork::wrap_plots(plots[panel_order], ncol=3)
      out_panel <- file.path(out_dir, paste0("Plot_UBL3_", dataset_tag, "_", region_tag, "_6panel.png"))
      ragg::agg_png(out_panel, width=16, height=9, units="in", res=300, background="white")
      print(p_all)
      dev.off()
    } else {
      cat("⚠ 6-panel 未生成：实际 plots 数量=", length(panel_order), "（可能某些 celltype donor 太少被跳过）\n")
    }
  }
  
  cat("\n==== sessionInfo ====\n")
  print(sessionInfo())
  cat("==== END ====\n")
  
  message("🎉 DONE region=", region_tag, " | out_dir=", out_dir)
}

## =========================
## 4) 开跑：EC 与 SFG（分别保存，不混）
## =========================
run_boxplots_deseq2_one_region(obj_fp_EC,  meta_fp_EC,  "EC",  root_EC)
run_boxplots_deseq2_one_region(obj_fp_SFG, meta_fp_SFG, "SFG", root_SFG)










#整体水平箱线图
###############################################################################
## NO4_syn21788402_EC_SFG_WholeCell_ADvsControl_perDonor_UBL3_boxplot_Assay5Layers.R
##
## 【严格模仿 GSE157827 NO12 模版】
##  - counts -> CP10k -> log1p -> donor mean
##  - donor 为统计单位：Wilcoxon + log2FC
##  - 图风格/颜色/标题与模式图一致
##
## 【syn21788402 特殊点】
##  - 两个脑区：EC 与 SFG 分别跑、分别保存
##  - meta 的 group 是 AD/NC，但图中对照必须显示 Control（NC->Control）
##  - 分组信息来自 stepP matched_cells_meta.csv（sample + group）
##
## 【输出】
##  - 每脑区 1 张图（共2张）：
##      Fig_NO4_UBL3_wholecell_perDonor_log1pCP10k_boxplot_ADvsControl.png
##  - 中间表：
##      UBL3_wholecell_perDonor_CP10k_mean_ADvsControl.csv
##      UBL3_wholecell_stats_ADvsControl.csv
##      QC_donor_counts_by_group.csv
###############################################################################

###############################################################################
## NO4_syn21788402_EC_SFG_WholeCell_ADvsControl_perDonor_UBL3_boxplot_SAFEcountsLayers.R
##
## 【严格模仿 GSE157827 NO12 模版分析流程】
##  - counts -> CP10k -> log1p -> donor mean
##  - donor 为统计单位：Wilcoxon + log2FC
##  - 图风格/颜色/标题与模式图一致
##
## 【syn21788402 特殊点】
##  - 两个脑区：EC 与 SFG 分别跑、分别保存（不混）
##  - stepP meta 的 group 是 AD/NC，但图中对照必须显示 Control（NC->Control）
##  - 分组信息来自 stepP matched_cells_meta.csv（sample + group）
##
## 【关键修复：避免 subscript out of bounds】
##  - Seurat v5 counts layers 里，可能存在某些 layer 没有 UBL3 行
##  - 仍按 NO12 思路逐 layer 遍历，但取 UBL3 counts 时：
##      若该 layer 有 target 行 -> 正常取
##      若没有 -> 该 layer 的 UBL3 counts 视为 0，并记录缺失 layer
##  - 不改变统计逻辑，只让脚本稳定跑完
###############################################################################

###############################################################################
## NO4_syn21788402_EC_SFG_WholeCell_ADvsControl_perDonor_UBL3_boxplot_FINALSAFE.R
##
## 【严格模仿 GSE157827 NO12】counts -> CP10k -> log1p -> donor mean
## 【统计】Wilcoxon + log2FC
## 【输出】EC 1张 + SFG 1张（分别保存，不混）
##
## 【关键修复】Seurat v5 多 counts layers 中取 UBL3 行时：
##   - 不使用 m["UBL3", ] 这种字符索引（容易 subscript out of bounds）
##   - 改用 idx <- match(target, rownames(m)) 再 m[idx, ] 数值索引（最稳）
##   - 任何 layer 出错会写入 DEBUG_counts_layers_trace.txt
###############################################################################

###############################################################################
## NO4_syn21788402_EC_SFG_WholeCell_ADvsControl_perDonor_UBL3_boxplot_COUNTSMERGE.R
##
## 【分析逻辑严格不变（等价于 NO12）】
##  counts -> CP10k -> log1p -> donor mean -> Wilcoxon + log2FC -> boxplot
##
## 【与 NO12 的唯一区别（为了解决你现在的 out-of-bounds）】
##  - NO12：逐 counts.* layer 遍历取 lib_size 与 UBL3 counts（在 syn21788402 EC 上会触发 Matrix 越界）
##  - 本脚本：先把 counts.* layers 合并成一个稳定的 gene×cell counts 矩阵，再在矩阵上计算
##    => 仍然是 counts 流程，数值完全等价，但更稳
##
## 【syn21788402 特殊点】
##  - 两个脑区：EC 与 SFG 分别跑、分别保存（不混）
##  - stepP meta 的 group 是 AD/NC，但图中对照必须显示 Control（NC->Control）
##  - 分组来源：stepP matched_cells_meta.csv（sample + group）
##
## 【输出】每脑区 1 张图（共2张）+ 中间表 + stats
###############################################################################

###############################################################################
## NO4_syn21788402_EC_SFG_WholeCell_ADvsControl_perDonor_UBL3_boxplot_FINAL.R
##
## 逻辑：counts -> CP10k -> log1p -> donor mean -> Wilcoxon + log2FC -> boxplot
## 组别：AD/NC -> AD/Control（图中固定显示 Control）
## 输出：EC 1张 + SFG 1张，分别保存在各自脑区目录内
##
## 修复点：避免 dplyr::count() 报 “group2 is list”
##  - 强制把 group2 / autopsy_id 等列 unlist 成普通向量
###############################################################################

rm(list=ls()); gc()
Sys.setenv(LANG="en")
SEED <- 20251023; set.seed(SEED)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(ragg)
})

dataset_tag <- "syn21788402"
gene_symbol <- "UBL3"

root_EC  <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/results"
root_SFG <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/results_SFG"

obj_fp_EC  <- file.path(root_EC,  "stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds")
obj_fp_SFG <- file.path(root_SFG, "stepH_syn21788402_SFG_obj_celltype6.rds")
stopifnot(file.exists(obj_fp_EC), file.exists(obj_fp_SFG))

meta_fp_EC  <- file.path(root_EC,  "stepP_syn21788402_matched_cells_meta.csv")
meta_fp_SFG <- file.path(root_SFG, "stepP_syn21788402_SFG_matched_cells_meta.csv")
stopifnot(file.exists(meta_fp_EC), file.exists(meta_fp_SFG))

## ===== counts 合并（你前面章节已验证稳定）=====
get_counts_matrix_allcells <- function(obj, assay="RNA") {
  a <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(a), error=function(e) character(0))
  counts_layers <- layers[grepl("^counts", layers)]
  if (length(counts_layers) == 0) {
    m2 <- tryCatch(Seurat::GetAssayData(obj, assay=assay, slot="counts"), error=function(e) NULL)
    if (!is.null(m2)) return(m2)
    stop("❌ 无法获取 counts。")
  }
  mats <- list()
  for (ly in counts_layers) {
    m <- tryCatch(SeuratObject::LayerData(a, layer=ly), error=function(e) NULL)
    if (!is.null(m)) mats[[ly]] <- m
  }
  ref <- rownames(mats[[1]])
  mats <- lapply(mats, function(m){
    if (identical(rownames(m), ref)) return(m)
    m2 <- Matrix::Matrix(0, nrow=length(ref), ncol=ncol(m), sparse=TRUE)
    rownames(m2) <- ref; colnames(m2) <- colnames(m)
    common <- intersect(ref, rownames(m))
    if (length(common) > 0) m2[common, ] <- m[common, , drop=FALSE]
    m2
  })
  mat_all <- if (length(mats)==1) mats[[1]] else Reduce(Matrix::cbind2, mats)
  dup <- duplicated(colnames(mat_all))
  if (any(dup)) mat_all <- mat_all[, !dup, drop=FALSE]
  mat_all <- mat_all[, colnames(obj), drop=FALSE]
  mat_all
}

run_wholecell_box_one_region <- function(obj_fp, meta_fp, region_tag, out_root) {
  
  out_dir <- file.path(out_root, paste0("NO4_WholeCell_ADvsControl_perDonor_", region_tag))
  dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
  
  ## 1) 读对象 + meta
  obj <- readRDS(obj_fp)
  DefaultAssay(obj) <- "RNA"
  md <- obj@meta.data
  md$cell_id <- colnames(obj)
  rownames(md) <- md$cell_id
  
  ## 2) 从 stepP meta 构建 sample->group2（AD/NC->Control）
  meta_raw <- read.csv(meta_fp, stringsAsFactors=FALSE)
  stopifnot(all(c("sample","group") %in% colnames(meta_raw)))
  
  ctrl_alias <- c("NC","Control","CTRL","Ctr","CTR","Normal","N","control","ctrl","ctr","nc","normal")
  meta_raw$group2 <- ifelse(trimws(as.character(meta_raw$group)) %in% ctrl_alias,
                            "Control", trimws(as.character(meta_raw$group)))
  
  map_by_sample <- meta_raw %>%
    transmute(sample_join = trimws(as.character(sample)),
              group2 = as.character(group2)) %>%
    filter(!is.na(sample_join), sample_join != "", group2 %in% c("AD","Control")) %>%
    distinct()
  
  ## 3) join_key：优先 SampleID
  join_key <- c("SampleID","sample","orig.ident")[c("SampleID","sample","orig.ident") %in% colnames(md)][1]
  if (is.na(join_key)) stop("❌ 对象 meta 缺少 SampleID/sample/orig.ident。")
  md$sample_join_in_obj <- trimws(as.character(md[[join_key]]))
  
  md2 <- dplyr::left_join(md, map_by_sample, by=c("sample_join_in_obj"="sample_join"))
  
  ## 丢弃 join 不到的 cells
  md2 <- md2[!is.na(md2$group2) & md2$group2 %in% c("AD","Control"), , drop=FALSE]
  rownames(md2) <- md2$cell_id
  
  if (!all(c("AD","Control") %in% unique(as.character(md2$group2)))) {
    stop("❌ 丢弃未映射细胞后不再同时包含 AD 与 Control。")
  }
  
  ## donor：优先 PatientID，否则 sample_join_in_obj
  if ("PatientID" %in% colnames(md2)) {
    md2$autopsy_id <- trimws(as.character(md2$PatientID))
  } else {
    md2$autopsy_id <- trimws(as.character(md2$sample_join_in_obj))
  }
  
  cells_keep <- intersect(rownames(md2), colnames(obj))
  md2 <- md2[cells_keep, , drop=FALSE]
  
  ## 4) counts 矩阵层面取子集，计算 CP10k/log1p
  cnt_all <- get_counts_matrix_allcells(obj, "RNA")
  cnt <- cnt_all[, cells_keep, drop=FALSE]
  
  gene_row <- if ("UBL3" %in% rownames(cnt)) "UBL3" else
    if ("ENSG00000122042" %in% rownames(cnt)) "ENSG00000122042" else NA_character_
  if (is.na(gene_row)) stop("❌ counts 里找不到 UBL3。")
  
  lib_size <- Matrix::colSums(cnt)
  ubl3_counts <- as.numeric(cnt[gene_row, , drop=TRUE])
  cp10k <- (ubl3_counts / pmax(lib_size, 1)) * 10000
  ubl3_log1p_cp10k <- log1p(cp10k)
  
  ## 5) donor 汇总（每 donor=1点）
  df_cells <- data.frame(
    autopsy_id = as.character(md2$autopsy_id),
    group2     = as.character(md2$group2),
    ubl3_cp10k = as.numeric(cp10k),
    ubl3_log   = as.numeric(ubl3_log1p_cp10k),
    stringsAsFactors = FALSE
  )
  
  ## ★关键修复：防止 group2/autopsy_id 变成 list 列
  df_cells$autopsy_id <- as.character(unlist(df_cells$autopsy_id, use.names=FALSE))
  df_cells$group2     <- as.character(unlist(df_cells$group2, use.names=FALSE))
  
  df_cells <- df_cells %>% filter(!is.na(autopsy_id), autopsy_id != "", !is.na(group2))
  
  df_donor <- df_cells %>%
    group_by(autopsy_id, group2) %>%
    summarise(
      n_cells = dplyr::n(),
      mean_cp10k = mean(ubl3_cp10k, na.rm=TRUE),
      mean_log1p_cp10k = mean(ubl3_log, na.rm=TRUE),
      .groups="drop"
    )
  
  ## 再次防止 list 列
  df_donor$autopsy_id <- as.character(unlist(df_donor$autopsy_id, use.names=FALSE))
  df_donor$group2     <- as.character(unlist(df_donor$group2, use.names=FALSE))
  
  ## donor QC（这里用 base::table，彻底避开 dplyr::count 的 list-column 坑）
  donor_tab <- as.data.frame(table(unique(df_donor[, c("autopsy_id","group2")])$group2),
                             stringsAsFactors = FALSE)
  colnames(donor_tab) <- c("group2","n_donors")
  write.csv(donor_tab, file.path(out_dir, "QC_donor_counts_by_group.csv"), row.names=FALSE)
  
  write.csv(df_donor,
            file.path(out_dir, "UBL3_wholecell_perDonor_CP10k_mean_ADvsControl.csv"),
            row.names=FALSE)
  
  ## 6) 统计：Wilcoxon + log2FC
  d1 <- df_donor %>% filter(group2=="AD")
  d0 <- df_donor %>% filter(group2=="Control")
  
  mu1 <- mean(d1$mean_cp10k, na.rm=TRUE)
  mu0 <- mean(d0$mean_cp10k, na.rm=TRUE)
  log2FC <- log2((mu1 + 1e-8) / (mu0 + 1e-8))
  
  pval <- if (nrow(d1)>=2 && nrow(d0)>=2) {
    wilcox.test(d1$mean_log1p_cp10k, d0$mean_log1p_cp10k, exact=FALSE)$p.value
  } else NA_real_
  
  stats_1 <- data.frame(
    contrast="AD_vs_Control",
    n_AD=nrow(d1),
    n_Control=nrow(d0),
    mean_cp10k_AD=mu1,
    mean_cp10k_Control=mu0,
    log2FC=log2FC,
    p_wilcox=pval,
    expr_source="merged_counts_matrix",
    used_object=basename(obj_fp),
    join_key_used=join_key,
    gene_row=gene_row,
    stringsAsFactors=FALSE
  )
  write.csv(stats_1, file.path(out_dir, "UBL3_wholecell_stats_ADvsControl.csv"), row.names=FALSE)
  
  ## 7) 画图（严格对齐你的模式图）
  group_order_plot <- c("AD","Control")
  df_donor$group_plot <- factor(df_donor$group2, levels=group_order_plot)
  
  fmt_p <- function(p) if (is.na(p)) "NA" else formatC(p, format="e", digits=2)
  sub1 <- sprintf("AD vs Control:  log2FC=%s, p=%s",
                  if (is.na(log2FC)) "NA" else sprintf("%.3f", log2FC),
                  fmt_p(pval))
  
  pal2_plot <- c(AD="#D24B40", Control="#2C7FB8")
  
  theme_sci <- theme_bw(base_size = 16) +
    theme(
      plot.title.position = "plot",
      plot.title    = element_text(face="bold", size=20, margin=margin(b=4)),
      plot.subtitle = element_text(size=12, colour="grey25", margin=margin(b=8)),
      panel.border  = element_rect(colour="grey25", fill=NA, linewidth=0.8),
      panel.grid.major.y = element_line(linewidth=0.28, linetype="dashed", colour="grey88"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.y = element_text(margin = margin(r=10)),
      legend.position = "none",
      plot.margin = margin(t=10, r=16, b=8, l=10)
    )
  
  p1 <- ggplot(df_donor, aes(x=group_plot, y=mean_log1p_cp10k, fill=group_plot)) +
    geom_boxplot(width=0.55, outlier.shape=NA, linewidth=1.0,
                 alpha=0.96, colour="grey15", median.linewidth=1.6) +
    geom_point(position=position_jitter(width=0.10, height=0),
               size=2.6, alpha=0.9, shape=21, stroke=0.5, colour="grey10") +
    scale_fill_manual(values=pal2_plot) +
    labs(
      title = "UBL3 expression per donor (Whole cells)",
      subtitle = sub1,
      x = NULL,
      y = "Mean UBL3 log1p(CP10k)"
    ) +
    theme_sci +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.10)))
  
  fig_fp <- file.path(out_dir, "Fig_NO4_UBL3_wholecell_perDonor_log1pCP10k_boxplot_ADvsControl.png")
  ragg::agg_png(fig_fp, width=7.8, height=5.6, units="in", res=450, background="white")
  print(p1)
  dev.off()
  
  cat("✅ saved: ", fig_fp, "\n")
  invisible(fig_fp)
}

## =========================
## 3) 开跑：EC 与 SFG（分别保存，不混）
## =========================
run_wholecell_box_one_region(obj_fp_EC,  meta_fp_EC,  "EC",  root_EC)
run_wholecell_box_one_region(obj_fp_SFG, meta_fp_SFG, "SFG", root_SFG)
