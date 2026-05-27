###############################################################################
# 01_per_dataset / syn52082747  (AD / PSP / FTD vs Control, primary visual
# cortex V1)
# Author-processed snRNA-seq harmonization: relabels the author's annotated
# object to the 6 harmonized cell types and standardizes diagnostic groups.
#
# Source     : Synapse accession syn52082747. Starts from the authors' processed
#              object (NO2_step7_obj_original_backup.rds) plus the authors'
#              subcluster metadata CSV (liger_subcluster_metadata_v2.csv); no
#              re-QC/clustering.
# Pipeline   : align cell barcodes to the author metadata -> derive group4
#              (AD / PSP / FTD / Control) -> marker-score clusters to the 6
#              harmonized cell types -> compute UBL3 log1p(CP10k) -> write the
#              supplementary count tables -> save a slimmed object.
# Produces   : stepH_slim_uncompressed.rds  (the object read by the integrated
#              figure/stat scripts for all syn52082747 units).
# Paths      : all paths point to the data-project location on disk. Source
#              objects are obtained from syn52082747; intermediate .rds objects
#              live with the data, not in this repository. Confirm before run.
# Environment: R + Seurat (v5) + SeuratObject + Matrix + data.table + dplyr;
#              see the sessionInfo logs this script writes for exact versions.
#
# NOTE: kept faithful to the script as run. Sections after the celltype6
#       hand-off overlap the integrated figures/tables and are retained only
#       for provenance.
###############################################################################

###############################################################################
## NO3A_make_stepH_slim_syn52082747.R
##
## 作用（只做这一件事）：
##  1) 读取起点对象 NO2_step7_obj_original_backup.rds
##  2) 读取作者 CSV（liger_subcluster_metadata_v2.csv）
##  3) 生成 group4（AD/PSP/FTD/Control）
##  4) 用 marker→cluster 打分得到 celltype6（模仿 GSE157827 第6章精神）
##  5) 计算 UBL3_log1pCP10K（严格：log1p((raw/lib_size)*1e4)）
##  6) 生成补充材料用的统计表（group/celltype 的细胞数+donor数+sample数）
##  7) 生成“瘦身版 StepH”并保存（稳定、不再坏）
##
## 输出：
##  - stepH_slim.rds
##  - Tables_*.csv
###############################################################################

rm(list = ls()); gc()  # 清空环境+释放内存
Sys.setenv(LANG="en")  # 避免中文乱码
SEED <- 20251023; set.seed(SEED)  # 固定随机种子（保证结果可复现）


suppressPackageStartupMessages({
  library(ggrepel)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(dplyr)
})

## ===== 路径（你只改这里）=====
raw_dir <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/rawdata"     # 原始数据路径
res_dir <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3"  # 结果输出路径
dir.create(res_dir, recursive=TRUE, showWarnings=FALSE)

obj_fp  <- file.path(raw_dir, "NO2_step7_obj_original_backup.rds") # 原始Seurat对象
meta_fp <- file.path(raw_dir, "liger_subcluster_metadata_v2.csv")   # 作者提供的元数据CSV
stopifnot(file.exists(obj_fp), file.exists(meta_fp))


## ============ 0) 环境记录（Table S3 用） ============
sink(file.path(res_dir, "sessionInfo_NO1_step1.txt"))
cat("=== Time ===\n"); print(Sys.time())
cat("\nSeurat: "); print(packageVersion("Seurat"))
cat("SeuratObject: "); print(packageVersion("SeuratObject"))
cat("\n=== sessionInfo ===\n"); print(sessionInfo())
sink()  


## ===== 1) 读对象（只读一次）=====
cat("▶ Read obj...\n")
obj <- readRDS(obj_fp)
stopifnot("RNA" %in% Assays(obj), "umap" %in% Reductions(obj))    # 检查关键assay/降维结果

## ===== 2) 读 CSV + 对齐 UMI =====
cat("▶ Read CSV...\n")
meta_csv <- data.table::fread(meta_fp, data.table=FALSE)
meta_csv$UMI <- as.character(meta_csv$UMI)

strip_barcode <- function(x){  # 标准化barcode格式（移除前缀/后缀）
  x <- as.character(x)
  x <- gsub("^.+_", "", x)   # 移除"样本名_"前缀
  x <- gsub("-1$", "", x)     # 移除末尾"-1"
  x
}

idx <- match(strip_barcode(colnames(obj)), strip_barcode(meta_csv$UMI))    # 匹配细胞barcode
meta_aligned <- meta_csv[idx, , drop=FALSE]     # 按Seurat对象的细胞顺序对齐元数据
rownames(meta_aligned) <- colnames(obj)      # 元数据行名=细胞barcode（关键对齐）

## ===== 3) group4 =====
normalize_group4 <- function(x) {   
  x0 <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(x0))
  out[grepl("ALZ|ALZH|ALZHEIMER|\\bAD\\b", x0)] <- "AD"
  out[grepl("\\bPSP\\b", x0)] <- "PSP"
  out[grepl("FTD|BVFTD|PICK|PICK'S|\\bPID\\b|FTLD", x0)] <- "FTD"
  out[grepl("CTRL|CONTROL|\\bNC\\b|NORMAL|HEALTH|NO DEMENTIA", x0)] <- "Control"
  out
}
dx_col <- if ("npdx1" %in% colnames(meta_aligned)) "npdx1" else "clinical_dx"

obj$group4     <- normalize_group4(meta_aligned[[dx_col]])
obj$sample     <- as.character(meta_aligned$sample)
obj$autopsy_id <- as.character(meta_aligned$autopsy_id)
obj$region     <- as.character(meta_aligned$region)

## ===== 4) celltype6（marker→cluster 打分）=====
if ("seurat_clusters" %in% colnames(obj@meta.data)) Idents(obj) <- "seurat_clusters"

obj$cluster_id <- paste0("c", as.character(Idents(obj)))
clu_levels_c   <- paste0("c", levels(Idents(obj)))

celltype6_levels_std <- c(
  "Astrocytes","Excitatory neurons","Microglia",
  "Endothelial","Inhibitory neurons","Oligodendrocytes"
)

markers_ref <- list(
  Astro   = c("AQP4","GFAP","ALDH1L1","SLC1A3","GPC5","RYR3"),
  Endo    = c("CLDN5","KDR","FLT1","PECAM1","ABCB1"),
  Excit   = c("SLC17A7","CAMK2A","TBR1","CBLN2","LDB2"),
  Inhib   = c("GAD1","GAD2","SLC6A1","PVALB","SST"),
  Microgl = c("CX3CR1","P2RY12","AIF1","C3","LRMDA"),
  Oligo   = c("MBP","PLP1","MOG","MOBP","ST18"),
  Peri    = c("PDGFRB","RGS5","MCAM","ACTA2")
)

DefaultAssay(obj) <- "RNA"
if (nrow(LayerData(obj[["RNA"]], layer="data")) == 0) {
  obj <- NormalizeData(obj, normalization.method="LogNormalize", scale.factor=1e4, verbose=FALSE)
}

score_one_group <- function(genes) {
  g <- intersect(genes, rownames(obj))   # 只保留数据集中存在的marker
  if (length(g) == 0) {
    out <- rep(NA_real_, length(clu_levels_c)); names(out) <- clu_levels_c; return(out)
  }
  av <- AverageExpression(obj, features=g, assays="RNA", layer="data", group.by="cluster_id", verbose=FALSE)$RNA    # 聚类平均表达
  sc <- colMeans(av, na.rm=TRUE)    # 每个聚类的marker平均得分
  sc2 <- sc[clu_levels_c]; names(sc2) <- clu_levels_c   # 按聚类顺序返回
  sc2
}

avg_by_cluster <- sapply(markers_ref, score_one_group)   # 所有细胞类型的聚类得分矩阵
tmp <- avg_by_cluster; tmp[is.na(tmp)] <- -Inf
lab_7 <- colnames(tmp)[max.col(tmp, ties.method="first")]   # 7类→合并周细胞到内皮
names(lab_7) <- rownames(tmp)
lab_7[lab_7=="Peri"] <- "Endo"

short2nice <- c(
  Astro="Astrocytes",
  Endo="Endothelial",
  Excit="Excitatory neurons",
  Inhib="Inhibitory neurons",
  Microgl="Microglia",
  Oligo="Oligodendrocytes"
)
lab_6 <- lab_7
lab_6[lab_6 %in% names(short2nice)] <- short2nice[lab_6[lab_6 %in% names(short2nice)]]
lab_6 <- factor(lab_6, levels=celltype6_levels_std)

map_cluster_to_celltype6 <- setNames(as.character(lab_6), names(lab_6))
celltype6_vec <- map_cluster_to_celltype6[obj$cluster_id]
names(celltype6_vec) <- colnames(obj)
obj$celltype6 <- factor(celltype6_vec, levels=celltype6_levels_std)

## ===== 5) 计算 UBL3_log1pCP10K（严格公式）=====
rna_counts <- LayerData(obj[["RNA"]], layer="counts") # 提取原始count矩阵
stopifnot("UBL3" %in% rownames(rna_counts))  # UBL3原始表达量
raw_ubl3 <- as.numeric(rna_counts["UBL3", ])  
lib_size <- Matrix::colSums(rna_counts)    # 每个细胞的总文库大小
obj$UBL3_log1pCP10K <- log1p((raw_ubl3 / lib_size) * 1e4)   # 严格按公式计算

write.csv(data.frame(UBL3_min=min(obj$UBL3_log1pCP10K,na.rm=TRUE),
                     UBL3_max=max(obj$UBL3_log1pCP10K,na.rm=TRUE)),
          file.path(res_dir, "Tables_UBL3_range.csv"), row.names=FALSE)

## ===== 6) 输出补充材料统计（只统计非NA的 group4/celltype6）=====
#核心输出：补充材料需要的 “分组 - 细胞类型” 维度的样本量统计，包含 3 个关键指标：n_cells：细胞数量；n_donors：独立捐赠者数量（去重 autopsy_id）；n_samples：独立样本数量（去重 sample）。
obj_stat <- subset(obj, subset = !is.na(group4) & !is.na(celltype6) & !is.na(sample) & !is.na(autopsy_id))
df <- data.frame(
  group4=obj_stat$group4,
  celltype6=obj_stat$celltype6,
  sample=obj_stat$sample,
  autopsy_id=obj_stat$autopsy_id,
  stringsAsFactors=FALSE
)

tabA <- df %>% group_by(group4) %>% summarise(
  n_cells=n(), n_donors=n_distinct(autopsy_id), n_samples=n_distinct(sample), .groups="drop"
)
write.csv(tabA, file.path(res_dir, "Tables_group4_nCells_nDonors_nSamples.csv"), row.names=FALSE)

tabB <- df %>% group_by(group4, celltype6) %>% summarise(
  n_cells=n(), n_donors=n_distinct(autopsy_id), n_samples=n_distinct(sample), .groups="drop"
)
write.csv(tabB, file.path(res_dir, "Tables_group4_by_celltype6_nCells_nDonors_nSamples.csv"), row.names=FALSE)

## ===== 7) DietSeurat 瘦身（关键：以后都用这个 slim 对象）=====
## 说明：
## - counts=TRUE：保证你后面还能按同公式算 UBL3 / 做箱线图
## - data=FALSE：更省；如果你后面要 FindAllMarkers，再改成 TRUE
obj_slim <- DietSeurat(
  obj,
  assays="RNA",   # 仅保留RNA assay
  dimreducs=c("pca","umap","tsne"),   # 保留降维结果
  graphs=NULL,    # 移除聚类图（节省空间）
  counts=TRUE,    # 保留原始count（后续可重新计算UBL3）
  data=FALSE,    # 移除归一化后的数据（节省空间）
  scale.data=FALSE     # 移除缩放数据（节省空间）
)

## =========================
## 用最稳方式保存 slim：先写 C:\Temp（不压缩）→ 读回 → copy 到 D盘 → md5校验
## =========================

tmp_dir <- "C:/Temp"
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

tmp_fp  <- file.path(tmp_dir, "stepH_slim_uncompressed_tmp.rds")
final_fp <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3/stepH_slim_uncompressed.rds"

## 1) 生成 slim（Seurat v5 推荐 layers 参数）
##    - 你后续算 UBL3 只需要 counts，所以保留 counts layer
obj_slim <- DietSeurat(
  obj,
  assays     = "RNA",
  dimreducs  = c("pca","umap","tsne"),
  graphs     = NULL,
  layers     = c("counts"),   # ✅ v5写法：只保留 counts
  scale.data = FALSE
)

gc()

## 2) 先写到 C:\Temp，且不压缩（最稳）
cat("▶ 保存到 C:/Temp（不压缩，最稳）...\n")
saveRDS(obj_slim, tmp_fp, compress = FALSE, version = 2)

## 3) 立刻读回验证（关键）
cat("▶ 读回验证...\n")
test <- readRDS(tmp_fp)
rm(test); gc()

## 4) copy 到 D盘结果目录
cat("▶ copy 到最终目录...\n")
ok <- file.copy(tmp_fp, final_fp, overwrite = TRUE)
stopifnot(ok)

## 5) md5 校验 copy 前后是否一致
m1 <- tools::md5sum(tmp_fp)
m2 <- tools::md5sum(final_fp)
stopifnot(identical(unname(m1), unname(m2)))

cat("✅ slim 保存成功且可读：", final_fp, "\n")
cat("   size(GB)=", round(file.info(final_fp)$size/1024^3, 2), "\n")














#第7章 UMAP
###############################################################################
## NO7_UMAP_syn52082747_GSE157827style_fromSlim.R
##
## 目标（严格模仿 GSE157827 模版）：
##  - 使用 stepH（含 umap + celltype6 + group4）
##  - 计算 UBL3：log1p((raw/lib_size)*1e4)  # log1p(CP10k)
##  - 画两联图：Panel A = celltype6 UMAP；Panel B = UBL3>0 高亮 UMAP
##  - 输出 4 张：Total + AD + PSP + FTD
##  - 固定 celltype6 顺序与配色（与你模版一致）
##  - UBL3 色条范围固定为“全局最大值”（可比性）
###############################################################################

rm(list = ls()); gc()
Sys.setenv(LANG="en")
SEED <- 20251023; set.seed(SEED)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(ggrepel)
  
  library(patchwork)
})

## =========================
## 1) 路径
## =========================
res_dir <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3"
obj_fp  <- file.path(res_dir, "stepH_slim_uncompressed.rds")
stopifnot(file.exists(obj_fp))

out_dir <- file.path(res_dir, "NO7_UMAP_fromSlim")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## =========================
## 2) 读取对象 + 核查
## =========================
cat("▶ read:", obj_fp, "\n")
obj <- readRDS(obj_fp)

stopifnot("umap" %in% Reductions(obj))
stopifnot(all(c("group4","celltype6") %in% colnames(obj@meta.data)))
stopifnot("counts" %in% Layers(obj[["RNA"]]))

## =========================
## 3) 固定 celltype6 顺序 + 配色（模版核心）
## =========================
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

## 只保留非NA（避免图里出现 NA）
obj <- subset(obj, subset = !is.na(group4) & !is.na(celltype6))
obj$celltype6 <- factor(as.character(obj$celltype6), levels = celltype6_levels_std)

## 组顺序（Control 在最后）
group_order <- c("AD","FTD","PSP","Control")
obj$group4 <- factor(as.character(obj$group4), levels = group_order)

## =========================
## 4) 计算 UBL3（严格模版公式）
## =========================
DefaultAssay(obj) <- "RNA"
rna_counts <- LayerData(obj[["RNA"]], layer = "counts")

if (!("UBL3" %in% rownames(rna_counts))) stop("❌ counts 行名中找不到 UBL3")

raw_ubl3 <- as.numeric(rna_counts["UBL3", ])
lib_size <- Matrix::colSums(rna_counts)
obj$UBL3_log1pCP10K <- log1p((raw_ubl3 / lib_size) * 1e4)

## 全局范围（所有图固定色条上限）
UBL3_max_global <- max(obj$UBL3_log1pCP10K, na.rm = TRUE)
write.csv(data.frame(UBL3_min=min(obj$UBL3_log1pCP10K, na.rm=TRUE),
                     UBL3_max=UBL3_max_global),
          file.path(out_dir, "Check_UBL3_range_global.csv"),
          row.names = FALSE)

## ===== UBL3 色条刻度（严格模版）=====
break_vals   <- c(0, 1, 2, 3, 4)
break_labels <- c("0", "1", "2", "3", "4+")

## =========================
## 5) 更稳的绘图函数：完全用 ggplot 画 Panel A/B
##     - 避免 DimPlot / FeaturePlot 在循环里触发 locked binding
## =========================

## 如果全局环境里有人把 name 锁了，先清掉（安全）
if (exists("name", envir = .GlobalEnv, inherits = FALSE)) {
  if (bindingIsLocked("name", .GlobalEnv)) unlockBinding("name", .GlobalEnv)
  rm(name, envir = .GlobalEnv)
}

## 提前把“全量数据框”准备好：一次提取，后面只过滤
um <- Embeddings(obj, "umap")
df_umap_all <- data.frame(
  umap_1   = um[,1],
  umap_2   = um[,2],
  celltype6 = obj$celltype6,
  group4    = obj$group4,
  UBL3      = obj$UBL3_log1pCP10K,
  stringsAsFactors = FALSE
)

## 主题（与你原来一致）
## =========================
## 主题（与你原来一致）
## =========================
umap_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face="bold", hjust=0.5, size=18),
    axis.title   = element_text(size=14),
    axis.text    = element_text(size=12),
    legend.title = element_text(face="bold"),
    panel.grid   = element_blank()
  )

## ============================================================
## ✅ 修改点：在函数里新增参数 show_title
##   - show_title=TRUE  ：显示左上角 “AD (6 cell types)” 这种总标题
##   - show_title=FALSE ：不显示（你要的）
## ============================================================
plot_umap_pair_df <- function(df_sub, tag, out_prefix, show_title = TRUE) {
  
  ## ---------- 计算标签位置（用 median，位置更稳定） ----------
  label_df <- df_sub %>%
    dplyr::filter(!is.na(celltype6)) %>%
    dplyr::group_by(celltype6) %>%
    dplyr::summarise(
      x = median(umap_1, na.rm = TRUE),
      y = median(umap_2, na.rm = TRUE),
      .groups = "drop"
    )
  
  ## ---------- Panel A：Cell type（用 ggplot 画点 + 标签） ----------
  pA <- ggplot(df_sub, aes(x = umap_1, y = umap_2, color = celltype6)) +
    geom_point(size = 0.10, alpha = 0.90) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(x = x, y = y, label = celltype6),
      inherit.aes = FALSE,
      size = 4.5,            # ✅ 与模版 label.size=4.5 对齐
      fontface = "plain",    # ✅ 关键：不要粗体（像 GSE157827）
      color = "black",
      box.padding = 0.30,
      point.padding = 0.15,
      max.overlaps = Inf,
      seed = SEED
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
      legend.box      = "horizontal",
      text = element_text(face = "plain")  # ✅ 防止全局变粗
    ) +
    guides(
      color = guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(size = 5)
      )
    )
  
  ## ---------- Panel B：UBL3>0 高亮 ----------
  df_bg  <- df_sub[df_sub$UBL3 <= 0 | is.na(df_sub$UBL3), ]
  df_pos <- df_sub[df_sub$UBL3 >  0 & !is.na(df_sub$UBL3), ]
  
  pB <- ggplot() +
    geom_point(
      data = df_bg,
      aes(x = umap_1, y = umap_2),
      color = "grey95",
      size  = 0.15,
      alpha = 0.40
    ) +
    geom_point(
      data = df_pos,
      aes(x = umap_1, y = umap_2, color = UBL3),
      size  = 0.25,
      alpha = 0.90
    ) +
    scale_color_gradientn(
      colours = c("#FEE0D2", "#FC9272", "#CB181D"),
      limits  = c(0, UBL3_max_global),
      breaks  = break_vals,
      labels  = break_labels,
      name    = "UBL3"
    ) +
    ggtitle("UBL3 expression (>0 highlighted)") +
    xlab("UMAP_1") + ylab("UMAP_2") +
    umap_theme +
    theme(legend.position = "bottom")
  
  ## ---------- 拼图并输出 ----------
  ## ✅ 修改点：疾病图不显示左上角 title（AD/FTD/PSP 那句）
  if (show_title) {
    pAB <- pA + pB +
      plot_layout(widths = c(1, 1)) +
      plot_annotation(title = tag, tag_levels = "A")
  } else {
    pAB <- pA + pB +
      plot_layout(widths = c(1, 1)) +
      plot_annotation(tag_levels = "A") +
      plot_annotation(theme = theme(
        plot.title = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
      ))
  }
  
  ggsave(file.path(out_dir, paste0(out_prefix, ".png")), pAB, width = 12, height = 6, dpi = 300)
  ggsave(file.path(out_dir, paste0(out_prefix, ".pdf")), pAB, width = 12, height = 6)
  cat("✅ 输出：", out_prefix, "\n")
}

## ============================================================
## ✅ 调用部分修改：
##   - Total：show_title=TRUE（保留 Total (6 cell types)）
##   - AD/PSP/FTD：show_title=FALSE（去掉图2那句）
## ============================================================

plot_umap_pair_df(df_umap_all,
                  "Total (6 cell types)",
                  "Fig3style_UMAP_Total_syn52082747",
                  show_title = TRUE)

for (g in c("AD","PSP","FTD")) {
  df_g <- df_umap_all[df_umap_all$group4 == g, ]
  plot_umap_pair_df(df_g,
                    paste0(g, " (6 cell types)"),
                    paste0("Fig3style_UMAP_", g, "_syn52082747"),
                    show_title = FALSE)
}








#第10章（1）直方图，这快代码老是有问题做不出6张完整图，所以另外还有1块代码保存在这个R文件夹内：NO3单个样本直方图2
###############################################################################
## NO10_Hist_syn52082747_GSE157827style_fromSlim.R
##
###############################################################################
## NO10_Hist_syn52082747_GSE157827style_AUTO.R
##
## 目标：严格按 GSE157827 模版的数据处理流程做“样本分面直方图”
## - UBL3：log1p((raw/lib_size)*1e4)  # log1p(CP10k)  ✅与模板一致
## - UBL3>0 过滤 ✅与模板一致
## - 每个 celltype6 输出 1 张图
## - facet_wrap(~ sample_group, ncol=5) ✅模板风格
##
## 重要：你的电脑 ggplot2/SeuratObject 会触发 locked binding
## - 脚本会优先尝试 ggplot2 输出（完全照抄模板）
## - 如果 ggplot2 保存失败，则自动改用 lattice 输出（布局一致、稳定）
###############################################################################

SEED <- 20251023; set.seed(SEED)
Sys.setenv(LANG="en")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

## ===========================
## 0) 路径（你的 syn52082747）
## ===========================
res_dir <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3"
obj_fp  <- file.path(res_dir, "stepH_slim_uncompressed.rds")
stopifnot(file.exists(obj_fp))

out_dir <- file.path(res_dir, "Fig_UBL3_hist_per_celltype_noZero2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_fp <- file.path(out_dir, "log_NO10_env.txt")

## ===========================
## 1) 读入对象（模板对应 obj_meta + obj_expr）
##    你这里 slim 对象已经同时包含表达与meta，因此都从一个对象取
## ===========================
obj_meta <- readRDS(obj_fp)
obj_expr <- obj_meta
DefaultAssay(obj_expr) <- "RNA"

## 必要字段检查
stopifnot(all(c("sample","group4","celltype6") %in% colnames(obj_meta@meta.data)))
stopifnot("counts" %in% Layers(obj_expr[["RNA"]]))

## 去 NA（避免 facet 出现 NA 面板）
obj_meta <- subset(obj_meta, subset = !is.na(sample) & !is.na(group4) & !is.na(celltype6))
obj_expr <- obj_meta  # 同步

## ===========================
## 2) 从 counts 里算 UBL3 的 log1p(CP10k)（逐行对齐模板）
## ===========================
rna_counts <- tryCatch(
  GetAssayData(obj_expr, assay = "RNA", slot = "counts"),
  error = function(e) LayerData(obj_expr[["RNA"]], layer = "counts")
)

## --- 找 UBL3 的 ENSEMBL ID（模板做法）---
ubl3_id <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = "UBL3",
  keytype = "SYMBOL",
  columns = "ENSEMBL"
)$ENSEMBL[1]

## --- 兼容 ENSEMBL / SYMBOL 行名（模板做法）---
gene_row <- if (!is.na(ubl3_id) && ubl3_id %in% rownames(rna_counts)) {
  ubl3_id
} else if ("UBL3" %in% rownames(rna_counts)) {
  "UBL3"
} else {
  stop("在 counts 中找不到 UBL3（无 ENSEMBL 也无 SYMBOL）")
}

## --- 对齐细胞（模板做法）---
common_cells <- intersect(colnames(obj_expr), colnames(obj_meta))
if (length(common_cells) == 0) stop("两个对象没有共同细胞名")

raw_vec  <- as.numeric(rna_counts[gene_row, common_cells])
lib_size <- Matrix::colSums(rna_counts[, common_cells, drop = FALSE])

## ✅ 这是你最关心的核心：与模板完全一致
UBL3_norm <- log1p((raw_vec / lib_size) * 1e4)  # log1p (CP10k)

meta_sub <- obj_meta@meta.data[common_cells, c("sample", "celltype6", "group4")]

df_all <- data.frame(
  cell      = common_cells,
  UBL3      = UBL3_norm,
  raw       = raw_vec,
  sample    = meta_sub$sample,
  celltype6 = meta_sub$celltype6,
  group     = meta_sub$group4,   # 模板叫 group，这里用 group4 赋值给 group
  stringsAsFactors = FALSE
)

## 只保留 UBL3 > 0（模板一致）
df_pos <- subset(df_all, UBL3 > 0)


## 3) 画图（模板风格）- 已改好版本
##    目标：让 legend（图例）里 Control 永远排在最后
##         并且在直方图重叠时，Control 最后绘制（视觉上在最上层）
## ===========================

## ---------------------------
## (A) 统一规定组别顺序（最关键！）
## 这会同时影响：
## 1) 图例（legend）的顺序
## 2)（可选）直方图重叠时的绘制顺序
## ---------------------------
group_levels <- c("AD", "FTD", "PSP", "Control")  # Control 放最后

## 强制把 df_pos$group 变成 factor，并按 group_levels 排序
## 如果 df_pos$group 还是字符型，ggplot 会按字母顺序/出现顺序，导致 Control 不在最后
df_pos$group <- factor(as.character(df_pos$group), levels = group_levels)

## ---------------------------
## (B) 颜色表：必须与 group_levels 对应
## 注意：names(group_colors) 的顺序也会影响 lattice 的 auto.key 显示顺序
## ---------------------------
group_colors <- c(
  AD      = "#D55E00",
  FTD     = "#009E73",
  PSP     = "#CC79A7",
  Control = "#0072B2"
)

## ---------------------------
## (C) facet 标签：sample_group = sample + "_" + group
## 说明：每个面板是一个 sample_group
## ---------------------------
df_pos$sample_group <- with(df_pos, paste0(sample, "_", group))

## ---------------------------
## (D) 面板顺序：按组 → 按样本数字 → 按样本名
## 说明：你原先想要“更贴近你之前要求”的布局，这里保留
## ---------------------------
df_pos$sample_num <- suppressWarnings(as.numeric(gsub("\\D+", "", df_pos$sample)))
tmp <- unique(df_pos[, c("sample","sample_num","group","sample_group")])

## 关键：tmp$group 也强制为 factor，并按 group_levels（Control 最后）
tmp$group <- factor(as.character(tmp$group), levels = group_levels)

## 排序：先组别（AD/FTD/PSP/Control），再样本数字，再样本名
tmp <- tmp[order(tmp$group, tmp$sample_num, tmp$sample), ]
ord_levels <- as.character(tmp$sample_group)

## 让 facet 面板严格按 ord_levels 顺序出现
df_pos$sample_group <- factor(as.character(df_pos$sample_group), levels = ord_levels)

## ---------------------------
## (E) 细胞类型列表
## ---------------------------
celltypes <- sort(unique(df_pos$celltype6))

## ---------------------------
## (F) 文件名安全化（防止空格/斜杠导致存图失败）
## ---------------------------
sanitize_name <- function(x) {
  x <- gsub(" ", "_", x)
  x <- gsub("/", "_", x)
  x
}

## ---------------------------
## (G) 记录环境日志（方便追查 locked binding 等问题）
## ---------------------------
sink(log_fp)
cat("SEED:", SEED, "\n")
cat("n_df_pos:", nrow(df_pos), "\n")
cat("groups (levels):", paste(levels(df_pos$group), collapse=", "), "\n\n")
print(sessionInfo())
sink()

## ===========================
## 4) 输出函数：优先 ggplot2，失败则 lattice
## ===========================




#开始的画图版本，先不删
plot_with_ggplot <- function(df_ct, ct) {
  
  ## -------------------------
  ## (1) 在每个 celltype 子集里再次强制 group 顺序
  ##     防止 df_ct 被筛选后 drop 掉 levels 或变字符
  ## -------------------------
  df_ct$group <- factor(as.character(df_ct$group), levels = group_levels)
  
  ## -------------------------
  ## (2) 控制“重叠时的绘制顺序”
  ## 你用 position="identity" 会重叠：
  ## - ggplot 的绘制顺序跟数据行顺序有关
  ## - 我们把 df_ct 按 group_levels 排序
  ## - 因为 Control 在 levels 最后，所以它会最后画，视觉上更“在上面”
  ## -------------------------
  df_ct <- df_ct[order(df_ct$group), ]
  
  ## -------------------------
  ## (3) 画图
  ## -------------------------
  p_ct <- ggplot(df_ct, aes(x = UBL3, fill = group)) +
    geom_histogram(
      bins     = 40,
      position = "identity",  # 重叠叠加
      alpha    = 0.7,
      color    = "grey30"
    ) +
    facet_wrap(~ sample_group, ncol = 5) +
    
    ## -----------------------
  ## (4) legend 顺序：用 breaks 显式指定
  ## breaks = group_levels → 图例显示顺序固定为 AD→FTD→PSP→Control
  ## drop = FALSE → 即使某个 celltype 子集里缺某组，也保留图例顺序（可按需关掉）
  ## -----------------------
  scale_fill_manual(
    values = group_colors,
    breaks = group_levels,
    drop   = FALSE,
    name   = "Group"
  ) +
    
    labs(
      title = paste0("UBL3>0 distribution in ", ct),
      x     = "UBL3 expression (log1p (CP10k))",
      y     = "Cell count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 16),
      strip.text       = element_text(face = "bold", size = 10),
      legend.position  = "top",
      legend.title     = element_text(face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ## 输出 png
  out_png <- file.path(out_dir, paste0("UBL3_noZero_", sanitize_name(ct), ".png"))
  ggsave(filename = out_png, plot = p_ct, width = 12, height = 8, dpi = 300)
  
  return(out_png)
}

plot_with_lattice <- function(df_ct, ct) {
  
  suppressPackageStartupMessages({
    library(lattice)
    library(grDevices)
  })
  
  ## -------------------------
  ## (1) 依旧强制 group 顺序（lattice 也需要）
  ## -------------------------
  df_ct$group <- factor(as.character(df_ct$group), levels = group_levels)
  
  ## -------------------------
  ## (2) bins = 40 → breaks 要 41 个点
  ## -------------------------
  rng <- range(df_ct$UBL3, finite = TRUE)
  breaks <- seq(from = rng[1], to = rng[2], length.out = 41)
  
  ## -------------------------
  ## (3) 只保留出现的面板（sample_group）
  ##     面板顺序仍然尽量沿用 ord_levels
  ## -------------------------
  lev_ct <- intersect(ord_levels, unique(as.character(df_ct$sample_group)))
  df_ct$sample_group <- factor(as.character(df_ct$sample_group), levels = lev_ct)
  
  ## -------------------------
  ## (4) lattice 的图例顺序：
  ##     auto.key$text = group_levels 或 names(group_colors)
  ##     这里用 group_levels 更直观（确保 Control 最后）
  ## -------------------------
  p <- lattice::histogram(
    ~ UBL3 | sample_group,
    data   = df_ct,
    layout = c(5, ceiling(nlevels(df_ct$sample_group) / 5)),
    breaks = breaks,
    type   = "count",
    xlab   = "UBL3 expression (log1p (CP10k))",
    ylab   = "Cell count",
    main   = paste0("UBL3>0 distribution in ", ct),
    strip  = strip.custom(bg = "grey90", par.strip.text = list(cex = 0.75, font = 2)),
    par.settings = list(
      superpose.polygon = list(
        col    = adjustcolor(group_colors, alpha.f = 0.70),
        border = "grey30"
      )
    ),
    auto.key = list(
      space       = "top",
      columns     = 4,
      rectangles  = TRUE,
      points      = FALSE,
      lines       = FALSE,
      text        = group_levels  # ✅ 保证图例顺序 Control 最后
    ),
    panel = function(x, subscripts, ...) {
      
      ## panel 内每个面板的数据是单一 group（因为你 sample_group = sample + "_" + group）
      ## 所以这里取该面板对应的 group，决定填充颜色
      g <- unique(df_ct$group[subscripts])
      g <- as.character(g[1])
      
      col_fill <- adjustcolor(group_colors[g], alpha.f = 0.70)
      panel.histogram(x, breaks = breaks, type = "count", col = col_fill, border = "grey30")
    }
  )
  
  out_png <- file.path(out_dir, paste0("UBL3_noZero_", sanitize_name(ct), ".png"))
  out_pdf <- file.path(out_dir, paste0("UBL3_noZero_", sanitize_name(ct), ".pdf"))
  
  png(out_png, width = 12, height = 8, units = "in", res = 300); print(p); dev.off()
  pdf(out_pdf, width = 12, height = 8); print(p); dev.off()
  
  return(out_png)
}

## ===========================
## 5) 循环输出每个 celltype 的图
## ===========================
for (ct in celltypes) {
  
  df_ct <- df_pos[df_pos$celltype6 == ct, ]
  if (nrow(df_ct) == 0) next
  
  ok <- TRUE
  outp <- NA_character_
  
  ## 优先 ggplot（模板）
  tryCatch({
    outp <- plot_with_ggplot(df_ct, ct)
  }, error = function(e) {
    ok <<- FALSE
    message("⚠ ggplot2 输出失败（将自动改用 lattice）: ", conditionMessage(e))
  })
  
  ## 若 ggplot 失败 → lattice
  if (!ok) {
    outp <- plot_with_lattice(df_ct, ct)
  }
  
  cat("✅ done:", ct, " -> ", outp, "\n")
}

cat("✅ NO10 完成：直方图输出到 ", out_dir, "\n")
cat("📌 环境日志：", log_fp, "\n")















#10.2 重叠直方图（3组在1个图上）
##
###############################################################################
## NO10_2_OverlapHist_syn52082747_SAMPLELEVEL_MAINFIG_CLEAN.R
##
## 改动点（按你反馈）：
##  1) panel 右上角 label 太长 → 改为极简：只显示 “KW  FDR=xx”
##  2) 各组样本数不再写在每个panel里 → 写到图例里：AD(n=10)...
##  3) 柱子不要发白/透明 → 提高 alpha（0.6），边框更细
##
## 不变点：
##  - UBL3 计算公式：log1p((raw/lib_size)*1e4)  # log1p(CP10k)
##  - 统计：样本级 median 主分析 + Kruskal–Wallis + BH-FDR
##  - 图形结构：overlap hist + facet(celltype6, free_y)
###############################################################################

SEED <- 20251023; set.seed(SEED)
Sys.setenv(LANG="en")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

## ============================================================================
## 0) 路径与参数
## ============================================================================
dataset_tag <- "syn52082747"
res_dir     <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3"
obj_file    <- file.path(res_dir, "stepH_slim_uncompressed.rds")
stopifnot(file.exists(obj_file))

out_dir <- file.path(res_dir, "NO10_2_overlap_hist_sampleLevel_MAINFIG_CLEAN")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gene_symbol <- "UBL3"
binwidth    <- 0.2

group_order <- c("AD","FTD","PSP","Control")

## 高级但清晰的四色（你之前那套保留）
group_colors <- c(
  AD      = "#E69F00",
  FTD     = "#009E73",
  PSP     = "#CC79A7",
  Control = "#4C72B0"
)

out_fig_png <- file.path(out_dir, paste0(dataset_tag, "_", gene_symbol, "_overlapHist_4groups_SAMPLELEVEL_CLEAN.png"))
out_fig_pdf <- file.path(out_dir, paste0(dataset_tag, "_", gene_symbol, "_overlapHist_4groups_SAMPLELEVEL_CLEAN.pdf"))

out_stats_median <- file.path(out_dir, paste0(dataset_tag, "_UBL3_condExpr_median_sampleLevel_stats.csv"))
out_summ_expr    <- file.path(out_dir, paste0(dataset_tag, "_UBL3_condExpr_sampleSummaries.csv"))
out_log          <- file.path(out_dir, "NO10_2_QC_log.txt")

## ============================================================================
## 1) 读对象 + 清理 NA
## ============================================================================
cat("▶ read:", obj_file, "\n")
obj <- readRDS(obj_file)
DefaultAssay(obj) <- "RNA"

stopifnot(all(c("sample","group4","celltype6") %in% colnames(obj@meta.data)))
stopifnot("counts" %in% Layers(obj[["RNA"]]))

obj <- subset(obj, subset = !is.na(sample) & !is.na(group4) & !is.na(celltype6))

## ============================================================================
## 2) 找 UBL3 行名（ENSG优先，其次SYMBOL）
## ============================================================================
gene_ensg <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = gene_symbol,
  keytype = "SYMBOL",
  columns = "ENSEMBL"
)$ENSEMBL[1]

rna_counts <- tryCatch(
  GetAssayData(obj, assay = "RNA", slot = "counts"),
  error = function(e) LayerData(obj[["RNA"]], layer = "counts")
)

features <- rownames(rna_counts)

gene_row <- if (!is.na(gene_ensg) && gene_ensg %in% features) {
  gene_ensg
} else if (gene_symbol %in% features) {
  gene_symbol
} else {
  stop("❌ counts 中找不到 UBL3（ENSG 或 SYMBOL 都不在行名）")
}
cat("✅ 使用基因行为:", gene_row, "\n")

## ============================================================================
## 3) 计算 UBL3 = log1p(CP10k)（严格模板公式）
## ============================================================================
cells <- colnames(obj)

raw_vec  <- as.numeric(rna_counts[gene_row, cells])
lib_size <- Matrix::colSums(rna_counts[, cells, drop = FALSE])

if (any(lib_size <= 0, na.rm = TRUE)) stop("❌ lib_size<=0，counts异常。")

## ✅ 模板核心公式（不动）
expr_cp10k_log1p <- log1p((raw_vec / lib_size) * 1e4)

if (any(!is.finite(expr_cp10k_log1p))) stop("❌ expr 出现 NaN/Inf。")

meta <- obj@meta.data[cells, c("sample","group4","celltype6")]
df_all <- data.frame(
  cell      = cells,
  expr      = as.numeric(expr_cp10k_log1p),
  raw       = raw_vec,
  sample    = as.character(meta$sample),
  celltype6 = as.character(meta$celltype6),
  group     = as.character(meta$group4),
  stringsAsFactors = FALSE
)

df_all$group <- factor(df_all$group, levels = group_order)

## 条件表达：只保留 UBL3>0 细胞
df_pos <- df_all %>% filter(expr > 0)
if (any(df_pos$expr <= 0, na.rm = TRUE)) stop("❌ df_pos 过滤异常（仍有<=0）。")

celltypes <- sort(unique(df_pos$celltype6))

## ============================================================================
## 4) 样本级汇总：median(expr)（主分析）
## ============================================================================
summ_expr <- df_pos %>%
  group_by(celltype6, sample, group) %>%
  summarise(
    n_cells_pos = n(),
    median_expr = median(expr, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summ_expr, out_summ_expr, row.names = FALSE)

## ============================================================================
## 5) 样本级检验：Kruskal–Wallis + BH-FDR（跨celltype）
## ============================================================================
stats_list <- lapply(sort(unique(summ_expr$celltype6)), function(ct) {
  
  d <- summ_expr[summ_expr$celltype6 == ct, , drop = FALSE]
  d$group <- factor(as.character(d$group), levels = group_order)
  
  n_by_group <- table(d$group)
  
  ## 默认多组KW
  p_raw <- tryCatch(
    kruskal.test(median_expr ~ group, data = d)$p.value,
    error = function(e) NA_real_
  )
  
  data.frame(
    celltype6 = ct,
    method    = "Kruskal–Wallis (sample-level, median)",
    p_raw     = p_raw,
    n_AD      = ifelse("AD"      %in% names(n_by_group), as.integer(n_by_group[["AD"]]),      0L),
    n_FTD     = ifelse("FTD"     %in% names(n_by_group), as.integer(n_by_group[["FTD"]]),     0L),
    n_PSP     = ifelse("PSP"     %in% names(n_by_group), as.integer(n_by_group[["PSP"]]),     0L),
    n_Control = ifelse("Control" %in% names(n_by_group), as.integer(n_by_group[["Control"]]), 0L),
    stringsAsFactors = FALSE
  )
})

stats_df <- do.call(rbind, stats_list)
stats_df$p_FDR <- p.adjust(stats_df$p_raw, method = "BH")

write.csv(stats_df, out_stats_median, row.names = FALSE)

## ============================================================================
## 6) 图例里写清楚各组样本数（避免每个panel重复写n）
##    注意：不同celltype里样本数可能略不同（例如microglia少一个FTD）
##    为了图例统一且不混乱，这里用“总体（全celltype合并）样本数”
## ============================================================================
overall_n <- summ_expr %>%
  distinct(sample, group) %>%
  count(group) %>%
  right_join(data.frame(group = factor(group_order, levels=group_order)), by="group") %>%
  mutate(n = ifelse(is.na(n), 0L, n))

legend_labels <- setNames(
  paste0(as.character(overall_n$group), " (n=", overall_n$n, ")"),
  as.character(overall_n$group)
)

## ============================================================================
## 7) panel右上角 label：只保留短信息（不挡柱子）
##    - 只显示：FDR（以及KW缩写）
## ============================================================================
stats_df$label_short <- sprintf("KW,FDR=%.2e", stats_df$p_FDR)

label_df <- stats_df %>%
  transmute(
    celltype6 = celltype6,
    label     = label_short,
    x_pos     = Inf,
    y_pos     = Inf
  )

## ============================================================================
## 8) 作图（主图：4组重叠直方图）
##    - alpha 提高到 0.6，避免你说的“发白透明”
##    - label 放右上角空白区（Inf定位 + margin）
## ============================================================================
p_main <- ggplot(df_pos, aes(x = expr, fill = group)) +
  geom_histogram(
    binwidth  = binwidth,
    alpha     = 0.60,          # ✅ 不要太透明
    position  = "identity",
    color     = "grey20",
    linewidth = 0.20
  ) +
  facet_wrap(~ celltype6, scales = "free_y") +
  scale_fill_manual(values = group_colors, labels = legend_labels, name = "Group") +
  geom_label(
    data        = label_df,
    inherit.aes = FALSE,
    aes(x = x_pos, y = y_pos, label = label),
    hjust      = 1.02,
    vjust      = 1.10,
    size       = 3.0,
    label.size = 0,
    fill       = "white",
    alpha      = 0.80
  ) +
  labs(
    title = paste0(gene_symbol, " expression per cell type (only expressed cells)"),
    x     = paste0(gene_symbol, " (log1p(CP10k))"),
    y     = "Cell count"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill="grey92", color="grey70"),
    strip.text       = element_text(face="bold"),
    plot.margin      = margin(10, 30, 10, 10),
    legend.position  = "right",
    legend.title     = element_text(face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

## ============================================================================
## 9) 保存（ggsave失败则设备兜底）
## ============================================================================
save_plot_safe <- function(p, out_png, out_pdf, w=10, h=6, dpi=300) {
  ok <- TRUE
  tryCatch({
    ggsave(out_png, p, width=w, height=h, dpi=dpi)
    ggsave(out_pdf, p, width=w, height=h)
  }, error = function(e) {
    ok <<- FALSE
    message("⚠ ggsave失败，改用设备兜底：", conditionMessage(e))
  })
  if (!ok) {
    tryCatch({
      png(out_png, width=w, height=h, units="in", res=dpi); print(p); dev.off()
      pdf(out_pdf, width=w, height=h); print(p); dev.off()
    }, error = function(e) {
      message("❌ 兜底保存也失败：", conditionMessage(e))
    })
  }
}
save_plot_safe(p_main, out_fig_png, out_fig_pdf, w=10, h=6, dpi=300)

## ============================================================================
## 10) 输出核查日志（保证你能核对正确性）
## ============================================================================
sink(out_log)
cat("SEED:", SEED, "\n")
cat("gene_symbol:", gene_symbol, "\n")
cat("gene_ensg:", gene_ensg, "\n")
cat("gene_row used:", gene_row, "\n\n")
cat("Total cells:", nrow(df_all), "\n")
cat("Expressed cells (expr>0):", nrow(df_pos), "\n\n")
cat("Cells(expr>0) by group:\n"); print(table(df_pos$group))
cat("\nOverall sample counts used in legend:\n"); print(overall_n)
cat("\nStats table (median, sample-level):\n"); print(stats_df)
cat("\nNOTE: 若多个celltype的FDR非常接近/相同，是BH校正+原始p值偏大造成的正常现象。\n")
cat("\nSessionInfo:\n"); print(sessionInfo())
sink()

cat("✅ 主图已生成：", out_fig_png, "\n")
cat("✅ 样本级统计表：", out_stats_median, "\n")
cat("✅ 样本级汇总表：", out_summ_expr, "\n")
cat("✅ 核查日志：", out_log, "\n")

























#重叠直方图，3组图各自分开，（Y轴=cell count），检验方法按照：sample统计

###############################################################################
## NO10_2_OverlapHist_syn52082747_3Pairs_vsControl_SAMPLELEVEL_GSE157827STRICT.R
##
###############################################################################
## NO10_2_OverlapHist_syn52082747_3Pairs_vsControl_SAMPLELEVEL_GSE157827_FINAL2.R
## - 3张图：AD/FTD/PSP vs Control
## - 样本级 median + Mann–Whitney U + BH-FDR
## - 图形严格走 GSE157827 风格（hist alpha=0.7 facet theme_bw）
## - 关键修复：给每个panel加右/上留白，把label放到留白区（不遮挡柱子）
###############################################################################

###############################################################################
## NO10_2_OverlapHist_syn52082747_3Pairs_vsControl_SAMPLELEVEL_GSE157827_CORNERLABEL.R
## - 3张图：AD/FTD/PSP vs Control
## - 统计：sample-level median + Mann–Whitney U + BH-FDR（不变）
## - 图形：保持 GSE157827 风格
## - 关键：label 固定在每个 panel 的“最右上角”(x=Inf,y=Inf)，不再用 hist 算坐标
###############################################################################

SEED <- 20251023; set.seed(SEED)
Sys.setenv(LANG="en")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

dataset_tag <- "syn52082747"
res_dir     <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3"
obj_file    <- file.path(res_dir, "stepH_slim_uncompressed.rds")
stopifnot(file.exists(obj_file))

out_dir <- file.path(res_dir, "NO10_2_overlap_hist_sampleLevel_3Pairs_vsControl_GSE157827_CORNERLABEL")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gene_symbol <- "UBL3"
binwidth    <- 0.2
pairs       <- c("AD","FTD","PSP")

## 读对象
obj <- readRDS(obj_file)
DefaultAssay(obj) <- "RNA"
stopifnot(all(c("sample","group4","celltype6") %in% colnames(obj@meta.data)))
stopifnot("counts" %in% Layers(obj[["RNA"]]))
obj <- subset(obj, subset = !is.na(sample) & !is.na(group4) & !is.na(celltype6))

## 找 UBL3 行名
gene_ensg <- AnnotationDbi::select(
  org.Hs.eg.db, keys=gene_symbol, keytype="SYMBOL", columns="ENSEMBL"
)$ENSEMBL[1]

rna_counts <- tryCatch(
  GetAssayData(obj, assay="RNA", slot="counts"),
  error=function(e) LayerData(obj[["RNA"]], layer="counts")
)
features <- rownames(rna_counts)
gene_row <- if (!is.na(gene_ensg) && gene_ensg %in% features) gene_ensg else if (gene_symbol %in% features) gene_symbol else stop("❌ counts中找不到UBL3")

## 计算 log1p(CP10k)
cells <- colnames(obj)
raw_vec  <- as.numeric(rna_counts[gene_row, cells])
lib_size <- Matrix::colSums(rna_counts[, cells, drop=FALSE])
if (any(lib_size <= 0, na.rm=TRUE)) stop("❌ lib_size<=0")
expr <- log1p((raw_vec / lib_size) * 1e4)

meta <- obj@meta.data[cells, c("sample","group4","celltype6")]
df_all <- data.frame(
  expr      = as.numeric(expr),
  sample    = as.character(meta$sample),
  celltype6 = as.character(meta$celltype6),
  group4    = trimws(as.character(meta$group4)),
  stringsAsFactors = FALSE
)

## 只用表达细胞
df_pos_all <- df_all[df_all$expr > 0, , drop=FALSE]
if (nrow(df_pos_all) == 0) stop("❌ expr>0 后无细胞")

run_pair <- function(disease){
  
  ## -------------------------
  ## 1) 取 pair 数据（细胞层）
  ## -------------------------
  expr_df2 <- df_pos_all %>%
    filter(group4 %in% c(disease,"Control")) %>%
    mutate(
      group = ifelse(group4==disease, disease, "Control"),
      group = factor(group, levels=c(disease,"Control"))
    )
  
  celltypes <- sort(unique(expr_df2$celltype6))
  
  ## -------------------------
  ## 2) 样本级汇总（median）
  ## -------------------------
  summ_expr <- expr_df2 %>%
    group_by(celltype6, sample, group) %>%
    summarise(median_expr = median(expr, na.rm=TRUE), .groups="drop")
  
  write.csv(
    summ_expr,
    file.path(out_dir, paste0(dataset_tag,"_",gene_symbol,"_sampleSummaries_exprGT0_",disease,"_vs_Control.csv")),
    row.names = FALSE
  )
  
  ## -------------------------
  ## 3) 每个 celltype：Mann–Whitney U + BH-FDR
  ## -------------------------
  stats_list <- lapply(celltypes, function(ct){
    d <- summ_expr[summ_expr$celltype6==ct, , drop=FALSE]
    d$group <- factor(as.character(d$group), levels=c(disease,"Control"))
    
    p_raw <- tryCatch(
      wilcox.test(median_expr ~ group, data=d, exact=FALSE)$p.value,
      error=function(e) NA_real_
    )
    
    data.frame(celltype6=ct, p_raw=p_raw, stringsAsFactors=FALSE)
  })
  stats_df <- do.call(rbind, stats_list)
  stats_df$p_FDR <- NA_real_
  ok <- which(is.finite(stats_df$p_raw))
  if (length(ok)>0) stats_df$p_FDR[ok] <- p.adjust(stats_df$p_raw[ok], method="BH")
  
  stats_df$label <- ifelse(
    is.finite(stats_df$p_FDR),
    sprintf("Mann–Whitney U\nFDR=%.2e", stats_df$p_FDR),
    "Mann–Whitney U\nFDR=NA"
  )
  
  write.csv(
    stats_df,
    file.path(out_dir, paste0(dataset_tag,"_",gene_symbol,"_stats_sampleMedian_MWU_",disease,"_vs_Control.csv")),
    row.names = FALSE
  )
  
  ## -------------------------
  ## 4) label_df：固定右上角坐标 (Inf,Inf)
  ##    注意：每个 panel 一条 label
  ## -------------------------
  label_df <- stats_df %>%
    transmute(
      celltype6 = celltype6,
      label     = label,
      x_pos     = Inf,
      y_pos     = Inf
    )
  
  ## -------------------------
  ## 5) 作图（GSE157827风格 + 右上角固定label）
  ## -------------------------
  p_hist <- ggplot(expr_df2, aes(x=expr, fill=group)) +
    geom_histogram(binwidth=binwidth, alpha=0.7, position="identity") +
    facet_wrap(~celltype6, scales="free_y") +
    scale_fill_manual(values=setNames(c("red","blue"), c(disease,"Control"))) +
    geom_label(
      data        = label_df,
      inherit.aes = FALSE,
      aes(x=x_pos, y=y_pos, label=label),
      hjust      = 1.02,     # 让label略微“贴外侧”，更像右上角
      vjust      = 1.02,
      size       = 2.5,
      label.size = 0,
      fill       = "white",
      alpha      = 0.7
    ) +
    labs(
      title = paste0(gene_symbol," expression per cell type (only expressed cells): ",disease," vs Control"),
      x     = paste0(gene_symbol," (log1p(CP10k))"),
      y     = "Cell count"
    ) +
    theme_bw() +
    theme(
      plot.margin = margin(10, 25, 10, 10)  # 给右侧留点空间，避免被裁切
    ) +
    coord_cartesian(clip="off")             # 允许label略微超出panel边界也不被裁切
  
  out_png <- file.path(out_dir, paste0(dataset_tag,"_",gene_symbol,"_hist_nonparam_",disease,"_vs_Control.png"))
  out_pdf <- file.path(out_dir, paste0(dataset_tag,"_",gene_symbol,"_hist_nonparam_",disease,"_vs_Control.pdf"))
  ggsave(out_png, p_hist, width=10, height=6, dpi=300)
  ggsave(out_pdf, p_hist, width=10, height=6)
  
  cat("✅ Saved:", out_png, "\n")
}

for (disease in pairs) run_pair(disease)

cat("\n🎉 完成：3张图已生成（CORNERLABEL），label固定每个panel最右上角，不会遮挡。\n")










#重叠直方图，Y轴是cell
##############################################################################
## NO4_01_OverlapHist_3Pairs_compare_donor_vs_cell_FINAL.R
##
## donor-level（人，autopsy_id） + cell-level（细胞）重叠直方图
## n 的定义与原文一致：donor = 41 (AD10 / Control10 / FTD9 / PSP11)
###############################################################################

rm(list=ls()); gc()
Sys.setenv(LANG="en")
SEED <- 20251023; set.seed(SEED)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

## =========================
## 0) 路径
## =========================
base_dir <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747"
obj_fp   <- file.path(base_dir, "results/NO3/stepH_slim_uncompressed.rds")
stopifnot(file.exists(obj_fp))

out_dir <- file.path(base_dir, "results/NO4/NO4_01_overlap_hist_compare")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

dataset_tag <- "syn52082747"
gene_symbol <- "UBL3"
binwidth    <- 0.2
pairs       <- c("AD","FTD","PSP")

log_fp <- file.path(out_dir, "NO4_01_log.txt")
sink(log_fp)
cat("==== NO4_01 START ====\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("obj_fp:", obj_fp, "\n")
cat("out_dir:", out_dir, "\n")
cat("SEED:", SEED, "\n\n")
sink()

## =========================
## 1) 读对象 + 基本核查
## =========================
obj <- readRDS(obj_fp)
DefaultAssay(obj) <- "RNA"

md <- obj@meta.data
stopifnot(all(c("sample","group4","celltype6","autopsy_id") %in% colnames(md)))

md$sample     <- trimws(as.character(md$sample))
md$donor      <- trimws(as.character(md$autopsy_id))  # ★ donor 定义
md$group4     <- trimws(as.character(md$group4))
md$celltype6  <- as.character(md$celltype6)

## donor 总数（原文口径）
don_all <- unique(md[, c("donor","group4")])
write.csv(
  as.data.frame(table(don_all$group4), stringsAsFactors=FALSE),
  file.path(out_dir, "QC_donors_by_group4.csv"),
  row.names=FALSE
)

sink(log_fp, append=TRUE)
cat("Total donors (all cells):", nrow(don_all), "\n")
cat("Donors by group4:\n"); print(table(don_all$group4))
cat("\n")
sink()

## =========================
## 2) 计算 UBL3 log1p(CP10k)
## =========================
rna_counts <- LayerData(obj[["RNA"]], layer="counts")

ubl3_ensg <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys="UBL3", keytype="SYMBOL", columns="ENSEMBL"
)$ENSEMBL[1]

gene_row <- if (!is.na(ubl3_ensg) && ubl3_ensg %in% rownames(rna_counts)) {
  ubl3_ensg
} else if ("UBL3" %in% rownames(rna_counts)) {
  "UBL3"
} else stop("❌ counts 行名中找不到 UBL3")

lib_size <- Matrix::colSums(rna_counts)
expr <- log1p((as.numeric(rna_counts[gene_row, ]) / pmax(lib_size,1)) * 1e4)

df0 <- data.frame(
  expr      = expr,
  sample    = md$sample,
  donor     = md$donor,
  group4    = md$group4,
  celltype6 = md$celltype6,
  stringsAsFactors = FALSE
)

## 只用表达细胞画分布
df0 <- df0[df0$expr > 0, ]

sink(log_fp, append=TRUE)
cat("Expressed cells (expr>0):", nrow(df0), "\n")
cat("Unique donors (expr>0):", length(unique(df0$donor)), "\n\n")
sink()

## =========================
## 3) 绘图函数
## =========================
plot_one_pair <- function(disease, unit=c("donor","cell")) {
  
  unit <- match.arg(unit)
  
  df2 <- df0 %>%
    filter(group4 %in% c(disease,"Control")) %>%
    mutate(group = ifelse(group4==disease, disease, "Control"))
  
  ## donor n（原文口径）
  don_pair <- unique(df2[, c("donor","group")])
  n_dis <- sum(don_pair$group == disease)
  n_ctl <- sum(don_pair$group == "Control")
  
  lab_dis <- sprintf("%s\n(n = %d)", disease, n_dis)
  lab_ctl <- sprintf("Control\n(n = %d)", n_ctl)
  
  df2$group_lab <- factor(
    ifelse(df2$group==disease, lab_dis, lab_ctl),
    levels=c(lab_dis, lab_ctl)
  )
  
  ## 统计输入
  if (unit=="donor") {
    stat_input <- df2 %>%
      group_by(celltype6, donor, group_lab) %>%
      summarise(val=median(expr), .groups="drop")
  } else {
    stat_input <- df2 %>%
      transmute(celltype6=celltype6, group_lab=group_lab, val=expr)
  }
  
  ## MWU + BH
  stats <- stat_input %>%
    group_by(celltype6) %>%
    summarise(
      p_raw = tryCatch(
        wilcox.test(val ~ group_lab, exact=FALSE)$p.value,
        error=function(e) NA_real_
      ),
      .groups="drop"
    )
  stats$padj <- p.adjust(stats$p_raw, method="BH")
  stats$label <- sprintf("Mann–Whitney U\nPadj=%.2e", stats$padj)
  stats$x <- Inf; stats$y <- Inf
  
  write.csv(
    stats,
    file.path(out_dir,
              paste0(dataset_tag,"_",gene_symbol,"_Stats_",
                     disease,"_vs_Control_",unit,".csv")),
    row.names=FALSE
  )
  
  fill_vals <- c("red","blue")
  names(fill_vals) <- c(lab_dis, lab_ctl)
  
  p <- ggplot(df2, aes(x=expr, fill=group_lab)) +
    geom_histogram(binwidth=binwidth, alpha=0.7, position="identity") +
    facet_wrap(~celltype6, scales="free_y") +
    scale_fill_manual(values=fill_vals, drop=FALSE) +
    geom_label(
      data=stats, inherit.aes=FALSE,
      aes(x=x, y=y, label=label),
      hjust=1.02, vjust=1.02, size=2.3,
      label.size=0, fill="white", alpha=0.7
    ) +
    labs(
      title = paste0(
        gene_symbol,
        " expression per cell type (only expressed cells): ",
        disease, " vs Control (", unit, "-level)"
      ),
      x = paste0(gene_symbol, " log1p(CP10k)"),
      y = "Cell count"
    ) +
    theme_bw() +
    theme(plot.margin = margin(10,25,10,10)) +
    coord_cartesian(clip="off")
  
  out_png <- file.path(
    out_dir,
    paste0(dataset_tag,"_",gene_symbol,
           "_OverlapHist_",disease,"_vs_Control_",unit,".png")
  )
  ggsave(out_png, p, width=10, height=6, dpi=300)
  cat("✅ saved:", out_png, "\n")
}

for (d in pairs) {
  plot_one_pair(d, "donor")  # ★原文对齐的人
  plot_one_pair(d, "cell")   # 分布对照
}

sink(log_fp, append=TRUE)
cat("\n==== sessionInfo ====\n")
print(sessionInfo())
cat("==== NO4_01 END ====\n")
sink()

cat("\n🎉 NO4_01 完成：donor-level + cell-level 共 6 张图\n输出目录：\n",
    out_dir, "\n", sep="")











#NO4_03：density 重叠图（6 张）+ 中间结果全保存
###############################################################################
###############################################################################
## NO4_03_OverlapHist_Density_3Pairs_compare_donor_vs_cell_FINAL.R
##
## 仍然是“柱状重叠直方图”（红/蓝，和 NO4_01 一样的格式）
## 但 Y 轴改为 Density（after_stat(density)）
## 输出：AD/FTD/PSP vs Control × (donor-level + cell-level) = 6 张
## 并保存中间结果与统计表（补充材料用）
###############################################################################

rm(list=ls()); gc()
Sys.setenv(LANG="en")
SEED <- 20251023; set.seed(SEED)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

## =========================
## 0) 路径
## =========================
base_dir <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747"
obj_fp   <- file.path(base_dir, "results/NO3/stepH_slim_uncompressed.rds")
stopifnot(file.exists(obj_fp))

out_dir <- file.path(base_dir, "results/NO4/NO4_03_overlap_hist_density_compare")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

dataset_tag <- "syn52082747"
gene_symbol <- "UBL3"
binwidth    <- 0.2
pairs       <- c("AD","FTD","PSP")

log_fp <- file.path(out_dir, "NO4_03_log.txt")
sink(log_fp)
cat("==== NO4_03 START ====\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("obj_fp:", obj_fp, "\n")
cat("out_dir:", out_dir, "\n")
cat("SEED:", SEED, "\n\n")
sink()

## =========================
## 1) 读对象 + donor 口径核查
## =========================
obj <- readRDS(obj_fp)
DefaultAssay(obj) <- "RNA"

md <- obj@meta.data
stopifnot(all(c("sample","group4","celltype6","autopsy_id") %in% colnames(md)))
stopifnot("counts" %in% Layers(obj[["RNA"]]))

md$sample     <- trimws(as.character(md$sample))
md$donor      <- trimws(as.character(md$autopsy_id))  # ★ donor（原文口径）
md$group4     <- trimws(as.character(md$group4))
md$celltype6  <- as.character(md$celltype6)

don_all <- unique(md[, c("donor","group4")])
qc_donors_by_group4 <- as.data.frame(table(don_all$group4), stringsAsFactors=FALSE)
colnames(qc_donors_by_group4) <- c("group4","n_donors")
write.csv(qc_donors_by_group4, file.path(out_dir, "QC_donors_by_group4.csv"), row.names=FALSE)

sink(log_fp, append=TRUE)
cat("Total donors (all cells):", nrow(don_all), "\n")
cat("Donors by group4:\n"); print(table(don_all$group4))
cat("\n")
sink()

## =========================
## 2) 计算 UBL3 log1p(CP10k) 并生成 df0（expr>0）
## =========================
rna_counts <- LayerData(obj[["RNA"]], layer="counts")

ubl3_ensg <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys="UBL3", keytype="SYMBOL", columns="ENSEMBL"
)$ENSEMBL[1]

gene_row <- if (!is.na(ubl3_ensg) && ubl3_ensg %in% rownames(rna_counts)) {
  ubl3_ensg
} else if ("UBL3" %in% rownames(rna_counts)) {
  "UBL3"
} else stop("❌ counts 行名中找不到 UBL3")

lib_size <- Matrix::colSums(rna_counts)
expr <- log1p((as.numeric(rna_counts[gene_row, ]) / pmax(lib_size,1)) * 1e4)

df_all <- data.frame(
  expr      = expr,
  sample    = md$sample,
  donor     = md$donor,
  group4    = md$group4,
  celltype6 = md$celltype6,
  stringsAsFactors = FALSE
)

df0 <- df_all[df_all$expr > 0, ]

## 保存中间结果（补充材料用）
saveRDS(df0, file.path(out_dir, "INTERMEDIATE_df0_exprGT0.rds"))
write.csv(df0, file.path(out_dir, "INTERMEDIATE_df0_exprGT0.csv"), row.names=FALSE)

sink(log_fp, append=TRUE)
cat("Expressed cells (expr>0):", nrow(df0), "\n")
cat("Unique donors in expr>0:", length(unique(df0$donor)), "\n")
cat("\n")
sink()

## =========================
## 3) 单个 pair：柱状重叠直方图（Y=density）
## =========================
plot_one_pair_hist_density <- function(disease, unit=c("donor","cell")) {
  
  unit <- match.arg(unit)
  
  df2 <- df0 %>%
    filter(group4 %in% c(disease,"Control")) %>%
    mutate(group = ifelse(group4==disease, disease, "Control"))
  
  ## 保存 pair 内细胞数据（中间结果）
  saveRDS(df2, file.path(out_dir, paste0("INTERMEDIATE_df2_", disease, "_vs_Control_exprGT0.rds")))
  
  ## donor n（图例用原文口径）
  don_pair <- unique(df2[, c("donor","group")])
  n_dis <- sum(don_pair$group == disease)
  n_ctl <- sum(don_pair$group == "Control")
  
  lab_dis <- sprintf("%s\n(n = %d)", disease, n_dis)
  lab_ctl <- sprintf("Control\n(n = %d)", n_ctl)
  
  df2$group_lab <- factor(
    ifelse(df2$group==disease, lab_dis, lab_ctl),
    levels=c(lab_dis, lab_ctl)
  )
  
  ## 统计输入（用于 MWU + BH）
  if (unit=="donor") {
    ## 每 donor 每 celltype 一条：median(expr)
    stat_input <- df2 %>%
      group_by(celltype6, donor, group_lab) %>%
      summarise(val=median(expr), .groups="drop")
  } else {
    stat_input <- df2 %>%
      transmute(celltype6=celltype6, group_lab=group_lab, val=expr)
  }
  
  write.csv(
    stat_input,
    file.path(out_dir, paste0("INTERMEDIATE_stat_input_", disease, "_vs_Control_", unit, ".csv")),
    row.names=FALSE
  )
  
  ## 每 celltype 每组 n（补充材料）
  n_by_celltype <- stat_input %>%
    group_by(celltype6, group_lab) %>%
    summarise(n = n(), .groups="drop")
  write.csv(
    n_by_celltype,
    file.path(out_dir, paste0("CHECK_n_by_celltype_", disease, "_vs_Control_", unit, ".csv")),
    row.names=FALSE
  )
  
  ## MWU + BH（按 celltype6）
  stats <- stat_input %>%
    group_by(celltype6) %>%
    summarise(
      p_raw = tryCatch(
        wilcox.test(val ~ group_lab, exact=FALSE)$p.value,
        error=function(e) NA_real_
      ),
      .groups="drop"
    )
  stats$padj <- p.adjust(stats$p_raw, method="BH")
  stats$label <- ifelse(
    is.finite(stats$padj),
    sprintf("Mann–Whitney U\nPadj=%.2e", stats$padj),
    "Mann–Whitney U\nPadj=NA"
  )
  stats$x <- Inf; stats$y <- Inf
  
  write.csv(
    stats,
    file.path(out_dir, paste0("STATS_", disease, "_vs_Control_", unit, ".csv")),
    row.names=FALSE
  )
  
  ## 颜色与外观：保持你之前“直方图版”的红/蓝格式
  fill_vals <- c("red","blue")
  names(fill_vals) <- c(lab_dis, lab_ctl)
  
  ## ★关键：y = after_stat(density)
  p <- ggplot(df2, aes(x = expr, y = after_stat(density), fill = group_lab)) +   # ★加 y=after_stat(density)
    geom_histogram(
      binwidth = binwidth,
      alpha = 0.7,
      position = "identity",
      colour = NA   # ★明确不要柱子外侧边框线
    ) +
    facet_wrap(~celltype6, scales="free_y") +
    scale_fill_manual(values=fill_vals, drop=FALSE) +
    geom_label(
      data=stats, inherit.aes=FALSE,
      aes(x=x, y=y, label=label),
      hjust=1.02, vjust=1.02, size=2.3,
      label.size=0, fill="white", alpha=0.7
    ) +
    labs(
      title = paste0(
        gene_symbol,
        " expression per cell type (only expressed cells): ",
        disease, " vs Control (", unit, "-level)"
      ),
      x = paste0(gene_symbol, " log1p(CP10k)"),
      y = "Density"   # ★改 y 轴标题
    ) +
    theme_bw() +
    theme(plot.margin = margin(10,25,10,10)) +
    coord_cartesian(clip="off")
  
  
  out_png <- file.path(
    out_dir,
    paste0(dataset_tag,"_",gene_symbol,
           "_OverlapHistDensity_", disease, "_vs_Control_", unit, ".png")
  )
  ggsave(out_png, p, width=10, height=6, dpi=300)
  
  sink(log_fp, append=TRUE)
  cat("---- Pair:", disease, "vs Control | unit=", unit, "\n", sep="")
  cat("Legend donor n: ", disease, "=", n_dis, " | Control=", n_ctl, "\n", sep="")
  cat("Saved figure:", out_png, "\n\n")
  sink()
  
  cat("✅ saved:", out_png, "\n")
}

for (d in pairs) {
  plot_one_pair_hist_density(d, "donor")
  plot_one_pair_hist_density(d, "cell")
}

sink(log_fp, append=TRUE)
cat("\n==== sessionInfo ====\n")
print(sessionInfo())
cat("==== NO4_03 END ====\n")
sink()

cat("\n🎉 NO4_03 完成：柱状重叠直方图（Y=density）6 张 + 中间结果已保存\n输出目录：\n",
    out_dir, "\n", sep="")
cat("📌 日志：", log_fp, "\n")








###############################################################################
## NO4_XX_syn52082747_SUMO_OverlapHistDensity_3Pairs_cellHist_donorMWU.R
##
## 【目标】syn52082747：3种病（AD/FTD/PSP）分别 vs Control
##  - 基因：SUMO1 / SUMO2 / SUMO3（3个基因都跑）
##  - 输出图：3个基因 × 3个pair = 9张
##
## 【重要口径（避免混淆）】
## 1) 直方图（图形形状）：cell-level（每个细胞一行），只用 expr>0 的细胞
##    - Y轴：Density（after_stat(density)）
## 2) 统计检验：donor-level（每 donor×celltype 汇总 median(expr)）
##    - 按 celltype6 分面分别做 Mann–Whitney U
##    - BH 校正得到 Padj，并标在每个分面右上角
## 3) 不输出 donor-level 直方图（只出主图 9张）
##
## 【表达计算】严格按 CP10k：
##    expr = log1p( (raw_counts / lib_size) * 1e4 )
##
## 【输入对象】固定：
##    D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3/stepH_slim_uncompressed.rds
##
## 【输出目录】建议新建一个独立目录，避免覆盖你 NO4_03：
##    D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO4/NO4_XX_overlap_hist_density_SUMOref_cellHist_donorMWU
##
## 【中间结果/核查输出】（补充材料与自查用）
##  - QC_donors_by_group4.csv：全局 donor 数
##  - 每个基因：INTERMEDIATE_{GENE}_df0_exprGT0.rds/csv（expr>0 全部细胞）
##  - 每个基因×pair：
##      * INTERMEDIATE_{GENE}_df2_{DIS}_vs_Control_exprGT0.rds/csv（pair内细胞数据）
##      * INTERMEDIATE_{GENE}_stat_input_{DIS}_vs_Control_donorMedian.rds/csv（donor median）
##      * CHECK_{GENE}_n_donors_by_celltype_{DIS}_vs_Control_exprGT0.csv（每celltype donor数）
##      * STATS_{GENE}_{DIS}_vs_Control_MWU_donorMedian_BH.csv（p_raw/padj）
##      * 图：{GENE}_{DIS}_vs_Control_cellHist_MWUdonorMedian.png
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
  ## 如果 rownames 是 ENSEMBL，下面两个包可用于 SYMBOL→ENSEMBL 映射
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

## =========================
## 0) 路径与参数
## =========================
base_dir <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747"
obj_fp   <- file.path(base_dir, "results/NO3/stepH_slim_uncompressed.rds")
stopifnot(file.exists(obj_fp))

out_dir <- file.path(base_dir, "results/NO4/NO4_XX_overlap_hist_density_SUMOref_cellHist_donorMWU")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

dataset_tag <- "syn52082747"
pairs       <- c("AD","FTD","PSP")              # 3种病分别 vs Control
gene_list   <- c("SUMO1","SUMO2","SUMO3")       # 3个候选内参基因
binwidth    <- 0.2

log_fp <- file.path(out_dir, "NO4_XX_log_SUMO.txt")
sink(log_fp)
cat("==== START ====\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("obj_fp:", obj_fp, "\n")
cat("out_dir:", out_dir, "\n")
cat("SEED:", SEED, "\n")
cat("pairs:", paste(pairs, collapse=", "), "\n")
cat("genes:", paste(gene_list, collapse=", "), "\n")
cat("binwidth:", binwidth, "\n\n")
sink()

## =========================
## 1) 读对象 + meta 列检查（保持你原 syn520 口径）
## =========================
obj <- readRDS(obj_fp)
DefaultAssay(obj) <- "RNA"
md <- obj@meta.data

## 这些列是你 syn520 既有流程中一直使用的口径
need_cols <- c("sample","group4","celltype6","autopsy_id")
if (!all(need_cols %in% colnames(md))) {
  sink(log_fp, append=TRUE)
  cat("❌ meta.data 缺少必要列。当前列名：\n")
  print(colnames(md))
  sink()
  stop(paste0("❌ meta.data 需要列：", paste(need_cols, collapse=", ")))
}

md$sample    <- trimws(as.character(md$sample))
md$donor     <- trimws(as.character(md$autopsy_id))  # ★ donor 口径：autopsy_id
md$group4    <- trimws(as.character(md$group4))
md$celltype6 <- as.character(md$celltype6)

## 统一 Control 的写法（防止出现 CTRL/NC 等）
ctrl_alias <- c("Control","CTRL","Ctr","CTR","NC","Normal","N","control","ctrl","ctr","nc","normal")
md$group4 <- ifelse(md$group4 %in% ctrl_alias, "Control", md$group4)

## 全局 donor QC（不依赖 expr>0，按全细胞）
don_all <- unique(md[, c("donor","group4")])
qc_don <- as.data.frame(table(don_all$group4), stringsAsFactors=FALSE)
colnames(qc_don) <- c("group4","n_donors")
write.csv(qc_don, file.path(out_dir, "QC_donors_by_group4.csv"), row.names=FALSE)

sink(log_fp, append=TRUE)
cat("Total donors (all cells):", nrow(don_all), "\n")
cat("Donors by group4:\n"); print(table(don_all$group4))
cat("\n")
sink()

## =========================
## 2) 获取 counts（Seurat v5 多 layers 兼容）
## =========================
get_counts_matrix_allcells <- function(obj, assay="RNA") {
  a <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(a), error=function(e) character(0))
  counts_layers <- layers[grepl("^counts", layers)]
  
  if (length(counts_layers) > 0) {
    sink(log_fp, append=TRUE)
    cat("Detected counts layers:\n"); print(counts_layers)
    sink()
    
    mats <- list()
    for (ly in counts_layers) {
      m <- tryCatch(SeuratObject::LayerData(a, layer=ly), error=function(e) NULL)
      if (is.null(m)) next
      if (!is.null(dim(m)) && length(dim(m))==2) mats[[ly]] <- m
    }
    if (length(mats)==0) stop("❌ counts layers 存在，但读取 LayerData 失败。")
    
    ## 对齐基因顺序
    ref_genes <- rownames(mats[[1]])
    for (k in names(mats)) {
      if (!identical(rownames(mats[[k]]), ref_genes)) {
        m0 <- mats[[k]]
        m_aligned <- Matrix::Matrix(0, nrow=length(ref_genes), ncol=ncol(m0), sparse=TRUE)
        rownames(m_aligned) <- ref_genes
        colnames(m_aligned) <- colnames(m0)
        common <- intersect(ref_genes, rownames(m0))
        m_aligned[common, ] <- m0[common, , drop=FALSE]
        mats[[k]] <- m_aligned
      }
    }
    
    mat_all <- if (length(mats)==1) mats[[1]] else Reduce(Matrix::cbind2, mats)
    
    ## 去重 cell
    dup <- duplicated(colnames(mat_all))
    if (any(dup)) {
      sink(log_fp, append=TRUE)
      cat("⚠ duplicated cells across layers:", sum(dup), " -> keep first\n")
      sink()
      mat_all <- mat_all[, !dup, drop=FALSE]
    }
    
    ## 补齐缺失 cell（以 obj 的细胞为准）
    all_cells <- colnames(obj)
    miss_cells <- setdiff(all_cells, colnames(mat_all))
    if (length(miss_cells) > 0) {
      sink(log_fp, append=TRUE)
      cat("⚠ missing cells in merged counts:", length(miss_cells), " -> fill zeros\n")
      sink()
      m_fill <- Matrix::Matrix(0, nrow=nrow(mat_all), ncol=length(miss_cells), sparse=TRUE)
      rownames(m_fill) <- rownames(mat_all)
      colnames(m_fill) <- miss_cells
      mat_all <- Matrix::cbind2(mat_all, m_fill)
    }
    
    mat_all <- mat_all[, all_cells, drop=FALSE]
    return(mat_all)
  }
  
  ## v4 兜底
  m2 <- tryCatch(Seurat::GetAssayData(obj, assay=assay, slot="counts"), error=function(e) NULL)
  if (!is.null(m2)) return(m2)
  
  stop("❌ 无法获取 counts。")
}

rna_counts <- get_counts_matrix_allcells(obj, "RNA")
stopifnot(ncol(rna_counts) == ncol(obj))

lib_size <- Matrix::colSums(rna_counts)

sink(log_fp, append=TRUE)
cat("Counts dim:", paste(dim(rna_counts), collapse=" x "), "\n\n")
sink()

## =========================
## 3) 基因行名定位（优先用 org.Hs.eg.db 做 SYMBOL→ENSEMBL；不“瞎编”ID）
## =========================
locate_gene_row <- function(counts_mat, gene_symbol) {
  rn <- rownames(counts_mat)
  
  ## 1) 先尝试 SYMBOL 直接存在
  if (gene_symbol %in% rn) return(gene_symbol)
  
  ## 2) 大小写不敏感匹配
  idx_ci <- which(toupper(rn) == toupper(gene_symbol))
  if (length(idx_ci) == 1) return(rn[idx_ci[1]])
  
  ## 3) 尝试 SYMBOL→ENSEMBL（如果 rownames 是 ENSG/ENSEMBL）
  ens_tbl <- tryCatch(
    AnnotationDbi::select(org.Hs.eg.db,
                          keys=gene_symbol, keytype="SYMBOL", columns=c("ENSEMBL")),
    error=function(e) NULL
  )
  if (!is.null(ens_tbl) && "ENSEMBL" %in% colnames(ens_tbl)) {
    ens_ids <- unique(ens_tbl$ENSEMBL[!is.na(ens_tbl$ENSEMBL)])
    m1 <- intersect(ens_ids, rn)
    if (length(m1) >= 1) return(m1[1])
    
    ## 兼容 ENSG.版本号（如果 counts 行名带 .数字）
    rn_strip <- sub("\\.\\d+$", "", rn)
    hit <- which(rn_strip %in% ens_ids)
    if (length(hit) >= 1) return(rn[hit[1]])
  }
  
  ## 4) 兜底：给出候选，便于你手动检查 rownames 格式
  cand <- rn[grep(gene_symbol, rn, ignore.case=TRUE)]
  stop(paste0(
    "❌ counts 行名中找不到基因：", gene_symbol, "\n",
    "请检查 rownames(rna_counts) 是否为 SYMBOL 或 ENSEMBL。\n",
    "grep 候选（前 20 个）：", paste(head(cand, 20), collapse=", ")
  ))
}

## =========================
## 4) 画图函数：单个基因 × 单个疾病pair（cell-hist + donor-MWU）
## =========================
plot_one_pair_for_gene <- function(df0, gene_symbol, disease) {
  
  ## 4.1 取出 pair 的表达细胞数据（expr>0 已经在 df0 保证）
  df2 <- df0 %>%
    filter(group4 %in% c(disease, "Control")) %>%
    mutate(group = ifelse(group4 == disease, disease, "Control"))
  
  ## 保存 pair 内细胞数据（中间结果，补充材料追溯）
  saveRDS(df2, file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_df2_", disease, "_vs_Control_exprGT0.rds")))
  write.csv(df2, file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_df2_", disease, "_vs_Control_exprGT0.csv")),
            row.names=FALSE)
  
  ## 4.2 legend 的 donor n（按该 pair 且 expr>0 的 donor 计数，口径与你前面一致）
  don_pair <- unique(df2[, c("donor","group")])
  n_dis <- sum(don_pair$group == disease)
  n_ctl <- sum(don_pair$group == "Control")
  
  lab_dis <- sprintf("%s\n(n = %d)", disease, n_dis)
  lab_ctl <- sprintf("Control\n(n = %d)", n_ctl)
  
  df2$group_lab <- factor(ifelse(df2$group == disease, lab_dis, lab_ctl),
                          levels=c(lab_dis, lab_ctl))
  
  ## 4.3 donor-level 统计输入：每 donor×celltype 的 median(expr)
  stat_input <- df2 %>%
    group_by(celltype6, donor, group_lab) %>%
    summarise(val = median(expr), .groups="drop")
  
  saveRDS(stat_input,
          file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_stat_input_", disease, "_vs_Control_donorMedian.rds")))
  write.csv(stat_input,
            file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_stat_input_", disease, "_vs_Control_donorMedian.csv")),
            row.names=FALSE)
  
  ## donor 数核查（每 celltype×组）
  n_don <- stat_input %>%
    group_by(celltype6, group_lab) %>%
    summarise(n_donors_with_exprGT0 = n(), .groups="drop")
  write.csv(n_don,
            file.path(out_dir, paste0("CHECK_", gene_symbol, "_n_donors_by_celltype_", disease, "_vs_Control_exprGT0.csv")),
            row.names=FALSE)
  
  ## 4.4 MWU + BH（按 celltype6）
  stats <- stat_input %>%
    group_by(celltype6) %>%
    summarise(
      p_raw = {
        g <- group_lab; v <- val
        ## 若某 celltype 只有一组 donor 有值，则无法检验，记 NA
        if (length(unique(g)) < 2) NA_real_
        else tryCatch(wilcox.test(v ~ g, exact=FALSE)$p.value, error=function(e) NA_real_)
      },
      .groups="drop"
    )
  
  stats$padj <- p.adjust(stats$p_raw, method="BH")
  stats$label <- ifelse(
    is.finite(stats$padj),
    sprintf("Mann–Whitney U\nPadj=%.2e", stats$padj),
    "Mann–Whitney U\nPadj=NA"
  )
  stats$x <- Inf; stats$y <- Inf
  
  write.csv(stats,
            file.path(out_dir, paste0("STATS_", gene_symbol, "_", disease, "_vs_Control_MWU_donorMedian_BH.csv")),
            row.names=FALSE)
  
  ## 4.5 颜色：红/蓝（与之前保持一致）
  fill_vals <- c("red","blue")
  names(fill_vals) <- c(lab_dis, lab_ctl)
  
  ## 4.6 绘图：cell-level 直方图（Y = Density），但标签来自 donor-level MWU
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
      ## 标题要写清楚：直方图是 cell-level，检验是 donor-level（防止你自己/审稿人误解）
      title = paste0(gene_symbol,
                     " expression per cell type (only expressed cells): ",
                     disease, " vs Control (donor-level)"),
      x = paste0(gene_symbol, " log1p(CP10k)"),
      y = "Density",
      fill = "group_lab"
    ) +
    theme_bw() +
    theme(plot.margin = margin(10,25,10,10)) +
    coord_cartesian(clip="off")
  
  out_png <- file.path(out_dir,
                       paste0(dataset_tag, "_", gene_symbol,
                              "_OverlapHistDensity_", disease, "_vs_Control_cellHist_MWUdonorMedian.png"))
  
  ragg::agg_png(out_png, width=10, height=6, units="in", res=300, background="white")
  print(p)
  dev.off()
  
  sink(log_fp, append=TRUE)
  cat("Saved figure:", out_png, "\n")
  cat("Pair:", disease, "vs Control | donor n (expr>0): ", disease, "=", n_dis, " | Control=", n_ctl, "\n\n", sep="")
  sink()
  
  cat("✅ saved: ", out_png, "\n", sep="")
  invisible(out_png)
}

## =========================
## 5) 主循环：3个基因 → 3个pair（共9张）
## =========================
summary_rows <- list()

for (gene_symbol in gene_list) {
  
  sink(log_fp, append=TRUE)
  cat("==== Gene START:", gene_symbol, "====\n")
  sink()
  
  ## 5.1 计算该基因 expr（CP10k）并构建 df0（expr>0）
  gene_row <- locate_gene_row(rna_counts, gene_symbol)
  
  raw_vec <- as.numeric(rna_counts[gene_row, , drop=TRUE])
  expr <- log1p((raw_vec / pmax(lib_size, 1)) * 1e4)
  
  df_all <- data.frame(
    expr      = expr,
    sample    = md$sample,
    donor     = md$donor,
    group4    = md$group4,
    celltype6 = md$celltype6,
    stringsAsFactors = FALSE
  )
  
  ## 只保留表达阳性细胞
  df0 <- df_all[df_all$expr > 0, , drop=FALSE]
  
  saveRDS(df0, file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_df0_exprGT0.rds")))
  write.csv(df0, file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_df0_exprGT0.csv")),
            row.names=FALSE)
  
  sink(log_fp, append=TRUE)
  cat("gene_row used:", gene_row, "\n")
  cat("Expressed cells (expr>0):", nrow(df0), "\n\n")
  sink()
  
  ## 5.2 对 3 个疾病分别画图（每个 disease vs Control）
  for (d in pairs) {
    plot_one_pair_for_gene(df0, gene_symbol, d)
    
    ## 汇总表：读回 stats 文件并合并（便于你最后选“最像内参”的基因）
    st_fp <- file.path(out_dir, paste0("STATS_", gene_symbol, "_", d, "_vs_Control_MWU_donorMedian_BH.csv"))
    st <- read.csv(st_fp, stringsAsFactors=FALSE)
    st$gene <- gene_symbol
    st$pair <- paste0(d, "_vs_Control")
    summary_rows[[paste(gene_symbol, d, sep="__")]] <- st[, c("gene","pair","celltype6","p_raw","padj","label")]
  }
  
  sink(log_fp, append=TRUE)
  cat("==== Gene END:", gene_symbol, "====\n\n")
  sink()
}

summary_df <- bind_rows(summary_rows)
write.csv(summary_df, file.path(out_dir, "SUMMARY_SUMO_genes_MWU_BH_by_celltype_and_pair.csv"),
          row.names=FALSE)

## =========================
## 6) sessionInfo 写入日志
## =========================
sink(log_fp, append=TRUE)
cat("\n==== sessionInfo ====\n")
print(sessionInfo())
cat("==== END ====\n")
sink()

cat("\n🎉 DONE.\nOutput dir:\n", out_dir, "\n", sep="")
cat("Log:\n", log_fp, "\n", sep="")
cat("图数量应为：3 genes × 3 pairs = 9 张（cell-hist + donor-MWU 标签）。\n")


























###############################################################################
###############################################################################
## NO4_04_Boxplots_DESeq2_byDonor_6panel_FINAL_CN.R
###############################################################################
###############################################################################
## NO4_04_Boxplots_DESeq2_byDonor_6panel_FINAL_CN_sparsePB.R
##
## 【你要的最终目标】
##  - donor 粒度（autopsy_id）pseudo-bulk + DESeq2（每个 celltype6 单独跑）
##  - 每个 celltype6 输出：
##      1) 1 张 donor-level 箱线图（配色/布局/标题/副标题格式固定）
##      2) 3 张 DEG 表：AD/FTD/PSP vs Control（每个 celltype 一套）
##  - 额外输出：
##      - 6-panel 总图（2×3）
##      - QC 表（donor 数、celltype×group 的 donor 数、celltype×group 的细胞数）
##      - pseudo-bulk 中间结果（pb 矩阵 + pb_meta）
##      - 每张箱线图的绘图数据表（INTERMEDIATE_plotdata_*.csv）
##      - 运行日志（NO4_04_log.txt）
##
## 【为什么要用 sparse pseudo-bulk】
##  - rowsum 需要 dense matrix，会提示分配上百 GB 内存并失败
##  - 因此这里用稀疏 one-hot 矩阵聚合：pb = cnt %*% M
##
## 【输入】
##  - stepH_slim_uncompressed.rds（必须包含 sample/group4/celltype6/autopsy_id）
##
## 【输出目录】
##  D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO4/NO4_05_boxplots_deseq2_byDonor_6panel
###############################################################################
###############################################################################

rm(list=ls()); gc()
Sys.setenv(LANG="en")
SEED <- 20251023; set.seed(SEED)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

## =========================
## 0) 路径与全局参数
## =========================
base_dir <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747"
obj_fp   <- file.path(base_dir, "results/NO3/stepH_slim_uncompressed.rds")
stopifnot(file.exists(obj_fp))

out_dir <- file.path(base_dir, "results/NO4/NO4_05_boxplots_deseq2_byDonor_6panel")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

dataset_tag <- "syn52082747"
gene_symbol <- "UBL3"

## 配色：保持你前面那套（红/绿/紫/蓝）
pal4 <- c(
  AD      = "#D24B40",
  FTD     = "#009E73",
  PSP     = "#CC79A7",
  Control = "#2C7FB8"
)

## 画图 X 轴顺序（你要求固定）
group_order_plot  <- c("AD","FTD","PSP","Control")

## DESeq2 的 reference（Control）
group_order_deseq <- c("Control","AD","FTD","PSP")

## 三个对照（写到 subtitle 里）
contrasts <- list(
  c("group4","AD","Control"),
  c("group4","FTD","Control"),
  c("group4","PSP","Control")
)

## 日志
log_fp <- file.path(out_dir, "NO4_04_log.txt")
sink(log_fp)
cat("==== NO4_04 START ====\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("obj_fp:", obj_fp, "\n")
cat("out_dir:", out_dir, "\n")
cat("SEED:", SEED, "\n\n")
sink()

## 为了避免你找不到输出：在控制台也打印一次 out_dir
cat("✅ out_dir = ", out_dir, "\n", sep="")

## =========================
## 1) 读对象 + 基本核查
## =========================
cat("▶ read:", obj_fp, "\n")
obj <- readRDS(obj_fp)
DefaultAssay(obj) <- "RNA"

md <- obj@meta.data
stopifnot(all(c("sample","group4","celltype6","autopsy_id") %in% colnames(md)))
stopifnot("counts" %in% Layers(obj[["RNA"]]))

## 统一成字符（避免 factor/list-column/空格）
md$sample    <- trimws(as.character(md$sample))
md$group4    <- trimws(as.character(md$group4))
md$celltype6 <- trimws(as.character(md$celltype6))

## donor 口径：autopsy_id（原文 41 人就是这个）
md$donor     <- trimws(as.character(md$autopsy_id))

## 只保留关键列不缺失的细胞（pseudo-bulk 必须干净）
keep_cells <- rownames(md)[
  !is.na(md$donor) & md$donor != "" &
    !is.na(md$group4) & md$group4 != "" &
    !is.na(md$celltype6) & md$celltype6 != ""
]
md2 <- md[keep_cells, , drop=FALSE]

sink(log_fp, append=TRUE)
cat("Total cells in object:", ncol(obj), "\n")
cat("Cells after filtering donor/group4/celltype6 not NA:", nrow(md2), "\n")
cat("Unique donors:", length(unique(md2$donor)), "\n")
cat("Donors by group4:\n"); print(table(unique(md2[,c("donor","group4")])$group4))
cat("Celltype6 counts (cell-level):\n"); print(table(md2$celltype6))
cat("\n")
sink()

## =========================
## 2) donor -> group4 一对一校验（顶刊必做）
##    不满足就停止并输出冲突表
## =========================
donor2grp <- unique(md2[, c("donor","group4"), drop=FALSE])

tab_donor_groups <- table(donor2grp$donor)         # 每 donor 出现多少行（理想=1）
bad_donors <- names(tab_donor_groups[tab_donor_groups > 1])

if (length(bad_donors) > 0) {
  bad_tbl <- donor2grp[donor2grp$donor %in% bad_donors, , drop=FALSE]
  write.csv(bad_tbl, file.path(out_dir, "CHECK_donor_maps_to_multiple_group4.csv"), row.names=FALSE)
  stop("❌ donor 对应多个 group4（标签混乱）。已输出：CHECK_donor_maps_to_multiple_group4.csv")
}

## QC：donor 数（总表）
qc_donor_counts <- as.data.frame(table(donor2grp$group4), stringsAsFactors = FALSE)
colnames(qc_donor_counts) <- c("group4","n_donors")
qc_donor_counts <- qc_donor_counts[order(qc_donor_counts$group4), ]
write.csv(qc_donor_counts, file.path(out_dir, "QC_donor_counts_by_group4.csv"), row.names=FALSE)

## =========================
## 3) pseudo-bulk：celltype6 × donor（counts 求和）
##    ★用稀疏聚合，不转 dense，避免 100GB+ 内存爆炸
## =========================
cnt <- LayerData(obj[["RNA"]], layer="counts")[, rownames(md2), drop=FALSE]

## UBL3 行名兼容（用于后续画图/汇总）
ubl3_ensg <- AnnotationDbi::select(
  org.Hs.eg.db, keys="UBL3", keytype="SYMBOL", columns="ENSEMBL"
)$ENSEMBL[1]

gene_row <- if (!is.na(ubl3_ensg) && ubl3_ensg %in% rownames(cnt)) {
  ubl3_ensg
} else if (gene_symbol %in% rownames(cnt)) {
  gene_symbol
} else if ("ENSG00000122042" %in% rownames(cnt)) {
  "ENSG00000122042"
} else {
  stop("❌ counts 行名中找不到 UBL3 或 ENSG00000122042")
}

## 为每个细胞建立分组键：celltype6__donor
pb_key <- paste(md2$celltype6, md2$donor, sep="__")

## 固定 group 顺序（保证可复现）
grp_levels <- sort(unique(pb_key))
grp <- factor(pb_key, levels = grp_levels)

## 构造稀疏 one-hot：cells × groups
## 行名必须与 cnt 的列名一致（我们前面用了 rownames(md2) 取列，因此严格对应）
M <- sparseMatrix(
  i = seq_along(grp),
  j = as.integer(grp),
  x = 1,
  dims = c(length(grp), length(grp_levels)),
  dimnames = list(rownames(md2), grp_levels)
)

## 关键：genes × groups  = (genes × cells) %*% (cells × groups)
pb <- cnt %*% M
pb <- as(pb, "dgCMatrix")

cat("✅ pseudo-bulk built (sparse): genes=", nrow(pb), " cols=", ncol(pb), "\n")

## pseudo-bulk 列注释表：每列属于哪个 celltype/donor/group
pb_meta <- data.frame(
  key       = colnames(pb),
  celltype6 = sub("__.*","", colnames(pb)),
  donor     = sub(".*__","", colnames(pb)),
  stringsAsFactors = FALSE
)

map_grp <- setNames(donor2grp$group4, donor2grp$donor)
pb_meta$group4 <- unname(map_grp[pb_meta$donor])

## 保存中间结果（补充材料/复现）
saveRDS(pb, file.path(out_dir, "INTERMEDIATE_pseudobulk_matrix_genes_by_celltype6__donor_dgCMatrix.rds"))
write.csv(pb_meta, file.path(out_dir, "INTERMEDIATE_pseudobulk_coldata_celltype6_donor_group4.csv"), row.names=FALSE)

## QC：每 celltype × group 的 donor 数
qc_donor_by_ct <- as.data.frame(table(pb_meta$celltype6, pb_meta$group4), stringsAsFactors = FALSE)
colnames(qc_donor_by_ct) <- c("celltype6","group4","n_donors")
write.csv(qc_donor_by_ct, file.path(out_dir, "QC_donor_counts_by_celltype6_by_group4.csv"), row.names=FALSE)

## QC：每 celltype × group 的细胞数（cell-level）
qc_cells_by_ct <- as.data.frame(table(md2$celltype6, md2$group4), stringsAsFactors = FALSE)
colnames(qc_cells_by_ct) <- c("celltype6","group4","cells_n")
write.csv(qc_cells_by_ct, file.path(out_dir, "QC_cells_by_celltype6_by_group4.csv"), row.names=FALSE)

## =========================
## 4) 工具函数：从 DESeq2::results 取 UBL3 指标（避免 S4 转换坑）
## =========================
get_ubl3_stats <- function(dds, con, gene_row) {
  res <- DESeq2::results(dds, contrast = con)
  list(
    log2FC = as.numeric(res[gene_row, "log2FoldChange"]),
    padj   = as.numeric(res[gene_row, "padj"]),
    pvalue = as.numeric(res[gene_row, "pvalue"])
  )
}

fmt_line <- function(dds, con, gene_row) {
  st <- get_ubl3_stats(dds, con, gene_row)
  l2s <- ifelse(is.na(st$log2FC), "NA", sprintf("%.3f", st$log2FC))
  pjs <- ifelse(is.na(st$padj),  "NA", format(st$padj, digits = 3, scientific = TRUE))
  sprintf("%s vs %s : log2FC=%s, padj=%s", con[2], con[3], l2s, pjs)
}

## =========================
## 5) 每个 celltype6：DESeq2 + 箱线图 + DEG 输出
## =========================
celltypes_all <- sort(unique(pb_meta$celltype6))

plots_single <- list()
res_ubl3_list <- list()

for (ct in celltypes_all) {
  
  cat("\n===========================\n")
  cat("▶ Celltype:", ct, " | time=", as.character(Sys.time()), "\n", sep="")
  
  cols <- pb_meta$key[pb_meta$celltype6 == ct]
  
  ## colData（group4），Control 做 reference
  cd_group <- pb_meta$group4[pb_meta$celltype6 == ct]
  cd_group <- factor(as.character(cd_group), levels = group_order_deseq)
  coldata <- data.frame(group4 = cd_group, row.names = cols)
  
  ## 取出该 celltype 的 pseudo-bulk counts
  ## 注意：cols 数量=donor 数（约10个），转 dense 很小，安全
  y <- pb[, cols, drop=FALSE]
  y_int <- round(as.matrix(y))
  storage.mode(y_int) <- "integer"
  
  ## DESeq2
  dds <- DESeqDataSetFromMatrix(countData = y_int, colData = coldata, design = ~ group4)
  dds <- DESeq(dds, quiet = TRUE)
  
  ## ------------------------------------------------------------
  ## 【DEG 输出】每个 celltype 输出 3 张 DEG 表（AD/FTD/PSP vs Control）
  ## ------------------------------------------------------------
  for (con in contrasts) {
    
    grp1 <- con[2]
    grp0 <- con[3]
    
    res <- DESeq2::results(dds, contrast = con)
    
    ## 转成 data.frame 并加 gene 列（此时是 genes×1 的 S4 结果转表，允许）
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    
    ## 常用列（补充材料最友好）
    keep_cols <- c("gene","log2FoldChange","padj","pvalue","baseMean")
    res_df <- res_df[, keep_cols]
    
    ## 按 padj 排序（显著的在前）
    res_df <- res_df[order(res_df$padj), ]
    
    out_deg <- file.path(out_dir, paste0("DEG_", gsub("[ /]","_",ct), "_", grp1, "_vs_", grp0, ".csv"))
    write.csv(res_df, out_deg, row.names = FALSE)
    
    ## 强制验收：写完立刻确认文件存在（避免“写了但没落盘”的错觉）
    ok <- file.exists(out_deg)
    cat(">>> DEG saved:", out_deg, " | exists=", ok, "\n")
    if (!ok) stop("❌ write.csv 后文件不存在：", out_deg)
  }
  
  ## ------------------------------------------------------------
  ## subtitle（三行对照：UBL3 的 log2FC + padj）
  ## ------------------------------------------------------------
  subtitle3 <- paste(
    fmt_line(dds, contrasts[[1]], gene_row),
    fmt_line(dds, contrasts[[2]], gene_row),
    fmt_line(dds, contrasts[[3]], gene_row),
    sep = "\n"
  )
  
  ## normalized counts（用于 donor-level 箱线图）
  norm <- counts(dds, normalized = TRUE)
  
  ## donor-level 作图数据（这里只画 UBL3）
  dfp <- data.frame(
    donor  = pb_meta$donor[pb_meta$celltype6 == ct],
    group4 = factor(pb_meta$group4[pb_meta$celltype6 == ct], levels = group_order_plot),
    value  = as.numeric(norm[gene_row, cols]),
    stringsAsFactors = FALSE
  )
  
  ## 保存绘图数据（补充材料/复现）
  write.csv(
    dfp,
    file.path(out_dir, paste0("INTERMEDIATE_plotdata_", gene_symbol, "_", gsub("[ /]","_",ct), "_byDonor.csv")),
    row.names = FALSE
  )
  
  ## UBL3 单基因汇总（补充材料汇总表）
  one <- do.call(rbind, lapply(contrasts, function(con0){
    st <- get_ubl3_stats(dds, con0, gene_row)
    data.frame(
      celltype6 = ct,
      contrast  = paste0(con0[2], "_vs_", con0[3]),
      log2FC    = st$log2FC,
      padj      = st$padj,
      pvalue    = st$pvalue,
      stringsAsFactors = FALSE
    )
  }))
  res_ubl3_list[[ct]] <- one
  
  ## ------------------------------------------------------------
  ## 画箱线图：颜色/布局/标题/副标题与你示例一致
  ## ------------------------------------------------------------
  p <- ggplot(dfp, aes(x = group4, y = value, fill = group4)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.95, colour = "grey15") +
    geom_point(
      position = position_jitter(width = 0.10),
      size = 2.6, alpha = 0.90,
      shape = 21, stroke = 0.4, colour = "grey10"
    ) +
    scale_fill_manual(values = pal4, drop = FALSE) +
    labs(
      title    = paste0(gene_symbol, " in ", ct, " (pseudo-bulk per donor)"),
      subtitle = subtitle3,
      x        = NULL,
      y        = "Normalized counts"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title.position = "plot",
      plot.title    = element_text(face = "bold", size = 16, margin = margin(b = 4)),
      plot.subtitle = element_text(size = 8.8, colour = "grey15", margin = margin(b = 6)),
      legend.position = "none"
    )
  
  ## 单图保存
  out_png <- file.path(out_dir, paste0(dataset_tag, "_", gene_symbol, "_Box_DESeq2_byDonor_", gsub("[ /]","_",ct), ".png"))
  ggsave(out_png, p, width = 6.8, height = 4.4, dpi = 300)
  
  plots_single[[ct]] <- p
  cat("✅ saved:", out_png, "\n")
}

## =========================
## 6) 输出：UBL3 汇总表（补充材料）
## =========================
res_ubl3 <- do.call(rbind, res_ubl3_list)
write.csv(res_ubl3, file.path(out_dir, paste0("DESeq2_", gene_symbol, "_summary_byCelltype6_byDonor.csv")),
          row.names = FALSE)

## =========================
## 7) 6-panel 总图（2×3）
## =========================
panel_order <- c("Astrocytes","Endothelial","Excitatory neurons",
                 "Inhibitory neurons","Microglia","Oligodendrocytes")
panel_order <- panel_order[panel_order %in% names(plots_single)]

p_all <- wrap_plots(plots_single[panel_order], ncol = 3)

out_all <- file.path(out_dir, paste0(dataset_tag, "_", gene_symbol, "_Box_DESeq2_byDonor_6panel.png"))
ggsave(out_all, p_all, width = 14, height = 7.5, dpi = 300)
cat("✅ saved 6-panel:", out_all, "\n")

## =========================
## 8) 日志 + sessionInfo + 输出文件清单（让你一眼看到 DEG 是否生成）
## =========================
sink(log_fp, append=TRUE)
cat("\n==== Files in out_dir (DEG_*.csv) ====\n")
print(list.files(out_dir, pattern="^DEG_.*\\.csv$", full.names=FALSE))
cat("\n==== sessionInfo ====\n")
print(sessionInfo())
cat("\n==== NO4_04 END ====\n")
cat("Time:", as.character(Sys.time()), "\n")
sink()

cat("\n🎉 NO4_04 完成：箱线图 + 6-panel + DEG 表（每 celltype 3 张）\n")
cat("📁 输出目录：", out_dir, "\n", sep="")
cat("🔎 DEG 文件数：", length(list.files(out_dir, pattern="^DEG_.*\\.csv$")), "\n", sep="")
###############################################################################














#SUMO当内参基因，做6个细胞类型的箱线图
###############################################################################
###############################################################################
## NO11_SUMO_ref_select_and_boxplots_FIXED.R
##
## 【目的】
## 1) donor-level pseudo-bulk（celltype6__autopsy_id）基础上，
##    用 DESeq2（design = ~ group4，Control 作为 reference）
##    计算 SUMO1 / SUMO2 / SUMO3 在每个 celltype 中的差异：
##      - AD vs Control
##      - FTD vs Control
##      - PSP vs Control
##
## 2) 自动选择“最稳定”的 SUMO 内参基因：
##    - 理想内参：在任何 celltype、任何 contrast 下都“不显著”(padj 大)、且效应小(|log2FC|小)
##    - 综合评分：
##        (A) padj 越大越好
##        (B) |log2FC| 越小越好
##
## 3) 用和 UBL3 完全同样的风格输出箱线图：
##    - 标题：GENE in Celltype (pseudo-bulk per donor)
##    - subtitle 三行：AD/FTD/PSP vs Control 的 log2FC/padj
##    - y 轴：DESeq2 normalized counts
##    - 四组顺序：AD, FTD, PSP, Control（Control 最右）
##
## 【输入文件（你已生成）】
## - NO11_Boxplots_4groups_DESeq2_byDonor/pbulk_counts_celltype6__autopsy_id.rds
## - NO11_Boxplots_4groups_DESeq2_byDonor/pbulk_coldata_celltype6__autopsy_id.csv
##
## 【输出目录】
## - NO11_Boxplots_4groups_DESeq2_byDonor/
##     ├─ SUMO_REF_selected/                 （最佳 SUMO 的箱线图）
##     ├─ SUMO_REF_all3_optional/            （可选：SUMO1/2/3 全部箱线图）
##     ├─ SUMO_stats_all.csv                 （每个 SUMO × celltype × contrast 的统计）
##     └─ SUMO_stability_ranking.csv         （稳定性排名 + 选择结果）
###############################################################################

## ========= 0) 基础设置 =========
Sys.setenv(LANG = "en")
options(stringsAsFactors = FALSE)

## ========= 0.1) 固定库路径（避免子进程找不到包）=========
lib <- "C:/Users/setou/AppData/Local/R/win-library/4.4"
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(ragg)     # ✅ 直接加载
})


SEED <- 20251023
set.seed(SEED)

## ========= 0.2) 【关键修复】locked binding / match.arg 异常防护 =========
## 你遇到的错误：
##   Error in match.arg(method) : cannot change value of locked binding for 'choices'
## 常见原因：
##   - 全局环境里出现了一个名字叫 choices 的对象，而且被 lockBinding() 锁住
##   - 或 match.arg 被某些包/脚本覆盖（极少见，但一旦发生会很玄学）
## 解决策略：
##   1) 强制使用 base::match.arg
##   2) 如果全局环境存在 choices 且被锁：unlockBinding + rm
fix_locked_choices <- function(verbose = TRUE) {
  
  ## (1) 强制 match.arg 是 base 版本
  ##     说明：如果 match.arg 被覆盖，DESeq2 内部调用 match.arg 会走错函数
  if (!identical(get("match.arg", mode = "function"), base::match.arg)) {
    if (verbose) message("⚠ match.arg 被覆盖/遮蔽：强制恢复为 base::match.arg")
    match.arg <- base::match.arg
  }
  
  ## (2) 清理/解锁全局 choices（最常见的罪魁祸首）
  if (exists("choices", envir = .GlobalEnv, inherits = FALSE)) {
    locked <- try(bindingIsLocked("choices", .GlobalEnv), silent = TRUE)
    
    ## 若被锁则解锁
    if (!inherits(locked, "try-error") && isTRUE(locked)) {
      if (verbose) message("⚠ 发现全局 choices 被 lockBinding：正在 unlockBinding() ...")
      unlockBinding("choices", .GlobalEnv)
    }
    
    ## 删除全局 choices
    if (verbose) message("🧹 删除全局对象 choices（避免 match.arg 路径被污染）")
    rm("choices", envir = .GlobalEnv)
  }
  
  ## (3) 轻量提示：如果你全局环境里有这些对象，最好别用同名（不强制删除）
  for (nm in c("method", "arg", "layer")) {
    if (exists(nm, envir = .GlobalEnv, inherits = FALSE)) {
      if (verbose) message("ℹ 注意：全局环境存在对象名：", nm, "（一般不建议与函数参数同名）")
    }
  }
  
  invisible(TRUE)
}

## 脚本开始先修复一次
fix_locked_choices(verbose = TRUE)

## ========= 1) 路径（只改这一处即可） =========
res_dir <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3"
out_dir <- file.path(res_dir, "NO11_Boxplots_4groups_DESeq2_byDonor")

pbulk_counts_fp  <- file.path(out_dir, "pbulk_counts_celltype6__autopsy_id.rds")
pbulk_coldata_fp <- file.path(out_dir, "pbulk_coldata_celltype6__autopsy_id.csv")

stopifnot(file.exists(pbulk_counts_fp), file.exists(pbulk_coldata_fp))

## 输出目录
dir_sel  <- file.path(out_dir, "SUMO_REF_selected")
dir_all3 <- file.path(out_dir, "SUMO_REF_all3_optional")
dir.create(dir_sel,  showWarnings = FALSE, recursive = TRUE)
dir.create(dir_all3, showWarnings = FALSE, recursive = TRUE)

## ========= 2) 读入 pseudo-bulk =========
counts_pbulk <- readRDS(pbulk_counts_fp)  # genes × (celltype__donor)
pb_annot     <- read.csv(pbulk_coldata_fp, row.names = 1, stringsAsFactors = FALSE)

celltypes <- sort(unique(pb_annot$celltype6))
cat("✅ celltypes =", paste(celltypes, collapse = " | "), "\n")

## ========= 3) 分组与配色（与你 UBL3 一致） =========
group_order_plot  <- c("AD", "FTD", "PSP", "Control")     # 图上 Control 最右
group_order_deseq <- c("Control", "AD", "FTD", "PSP")     # DESeq2 中 Control 做 reference

pal4 <- c(
  AD      = "#D24B40",
  Control = "#2C7FB8",
  FTD     = "#009E73",
  PSP     = "#CC79A7"
)

## ========= 4) SUMO 候选 =========
sumo_genes <- c("SUMO1", "SUMO2", "SUMO3")

missing_sumo <- sumo_genes[!sumo_genes %in% rownames(counts_pbulk)]
if (length(missing_sumo) > 0) {
  stop("❌ counts_pbulk 行名中找不到这些 SUMO：", paste(missing_sumo, collapse = ", "),
       "\n请确认 pseudo-bulk 的行名是否为 SYMBOL。")
}

## contrasts：固定 3 个
contrasts <- list(
  c("group4", "AD",  "Control"),
  c("group4", "FTD", "Control"),
  c("group4", "PSP", "Control")
)

## ========= 5) 稳定写 PNG：原子写入 + 校验 + 重试（避免黑图/缺图） =========
## ========= 5) 稳定写 PNG（强烈推荐版） =========
## 核心思想：
##  - 先写到 tempdir()（避免被资源管理器缩略图/杀毒/云同步抢占）
##  - 成功后再复制到目标路径（目标路径被占用也不影响绘图）
##  - 文件大小阈值提高，避免“半张图”也被当成功
save_png_atomic <- function(filename, plot, width=8.4, height=5.6, dpi=320, retries=5) {
  
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  
  ## ✅ 临时输出：写到系统 temp（最不容易被占用/锁定）
  tmpfile <- file.path(
    tempdir(),
    paste0("tmp_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", sample.int(1e9, 1), ".png")
  )
  
  for (k in seq_len(retries)) {
    
    if (file.exists(tmpfile)) file.remove(tmpfile)
    
    ok <- FALSE
    try({
      ## ✅ 强制 ragg：稳定且跨平台
      ragg::agg_png(tmpfile, width=width, height=height, units="in", res=dpi, background="white")
      print(plot)
      dev.off()
      
      ## ✅ 更严格的完整性判定：至少 80KB（你这种大图通常远大于此）
      if (file.exists(tmpfile) && is.finite(file.info(tmpfile)$size) && file.info(tmpfile)$size > 80*1024) {
        ok <- TRUE
      }
    }, silent = TRUE)
    
    if (ok) {
      ## ✅ 复制到目标路径（先删再写，避免旧文件残留）
      if (file.exists(filename)) file.remove(filename)
      ok2 <- file.copy(tmpfile, filename, overwrite = TRUE)
      
      ## ✅ 再做一次校验：目标文件也要足够大
      if (isTRUE(ok2) && file.exists(filename) && file.info(filename)$size > 80*1024) {
        file.remove(tmpfile)
        gc(verbose = FALSE)  # ✅ 每张图后清理一下内存/设备残留
        return(TRUE)
      }
    }
    
    Sys.sleep(0.4)
  }
  
  if (file.exists(tmpfile)) file.remove(tmpfile)
  gc(verbose = FALSE)
  return(FALSE)
}



## ========= 6) 生成 subtitle（三行 log2FC/padj） =========
## 注意：results(dds) 行名是基因名（与 counts 一致）
fmt_one <- function(dds, con, gene_row) {
  case <- con[2]; ctrl <- con[3]
  
  ## results 也可能报错（比如某组没有样本/对比无意义），所以做 tryCatch
  res <- tryCatch(
    as.data.frame(results(dds, contrast = con)),
    error = function(e) NULL
  )
  
  if (is.null(res) || !(gene_row %in% rownames(res))) {
    return(sprintf("%s vs %s: log2FC=NA, padj=NA", case, ctrl))
  }
  
  l2 <- res[gene_row, "log2FoldChange"]
  pj <- res[gene_row, "padj"]
  
  pj_str <- if (is.na(pj)) "NA" else format(pj, digits = 3, scientific = TRUE)
  l2_str <- if (is.na(l2)) "NA" else sprintf("%.3f", l2)
  
  sprintf("%s vs %s: log2FC=%s, padj=%s", case, ctrl, l2_str, pj_str)
}

## ========= 7) 画图函数（与 UBL3 同风格） =========
pretty_box <- function(dfp, title, subtitle) {
  dfp$group4 <- factor(as.character(dfp$group4), levels = group_order_plot)
  
  ggplot(dfp, aes(x = group4, y = value, fill = group4)) +
    geom_boxplot(
      width = 0.55,
      outlier.shape = NA,
      linewidth = 1.15,
      median.linewidth = 1.9,
      alpha = 0.96,
      colour = "grey15"
    ) +
    geom_point(
      position = position_jitter(width = 0.10, height = 0),
      size = 2.6,
      alpha = 0.90,
      shape = 21,
      stroke = 0.5,
      colour = "grey10"
    ) +
    scale_fill_manual(values = pal4) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Normalized counts") +
    theme_bw(base_size = 16) +
    theme(
      plot.title.position = "plot",
      plot.title    = element_text(face = "bold", size = 20, margin = margin(b = 4)),
      plot.subtitle = element_text(size = 12, colour = "grey20", margin = margin(b = 8)),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

## ========= 8) 主计算：收集 SUMO 的 log2FC / padj =========
stats_list <- list()

for (ct in celltypes) {
  
  ## ★ 每个 celltype 开始前都再修复一次（防止某一步污染全局环境）
  fix_locked_choices(verbose = FALSE)
  
  cat("\n====================================================\n")
  cat(">>> DESeq2 | celltype =", ct, "\n")
  
  ## --- 8.1 取出该 celltype 的 donor 列 ---
  keep_cols <- rownames(pb_annot)[pb_annot$celltype6 == ct & !is.na(pb_annot$group4)]
  
  ## donor 太少就跳过（避免 DESeq2 不稳定）
  if (length(keep_cols) < 8) {
    cat("⚠ donor 太少，跳过：", ct, "\n")
    next
  }
  
  y  <- counts_pbulk[, keep_cols, drop = FALSE]
  cd <- pb_annot[keep_cols, , drop = FALSE]
  cd$group4 <- factor(as.character(cd$group4), levels = group_order_deseq)
  
  ## --- 8.2 过滤低表达：但强制保留 SUMO1/2/3 ---
  keep_genes <- rowSums(y >= 10) >= 3
  for (g in sumo_genes) keep_genes <- keep_genes | (rownames(y) == g)
  y <- y[keep_genes, , drop = FALSE]
  
  ## --- 8.3 counts 必须是整数矩阵（DESeq2 要求）---
  ## 注意：这里转 as.matrix 可能占内存，但每个 celltype 的列数不大（donor 级别），通常可接受
  y_int <- round(as.matrix(y))
  storage.mode(y_int) <- "integer"
  dimnames(y_int) <- dimnames(y)
  
  cd_use <- data.frame(group4 = cd$group4, row.names = rownames(cd), check.names = FALSE)
  
  ## --- 8.4 构建 DESeq2 对象 ---
  dds <- DESeqDataSetFromMatrix(countData = y_int, colData = cd_use, design = ~ group4)
  
  ## --- 8.5 运行 DESeq：用 tryCatch 防止某个 celltype 直接把全脚本炸停 ---
  dds <- tryCatch(
    DESeq(dds, fitType = "parametric", quiet = TRUE),
    error = function(e) {
      message("❌ DESeq2 失败 | celltype=", ct, " | ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(dds)) next
  
  ## --- 8.6 对每个 SUMO × contrast 记录 stats ---
  for (g in sumo_genes) {
    for (con in contrasts) {
      case <- con[2]; ctrl <- con[3]
      
      res <- tryCatch(
        as.data.frame(results(dds, contrast = con)),
        error = function(e) NULL
      )
      
      stats_list[[length(stats_list) + 1]] <- data.frame(
        gene      = g,
        celltype6 = ct,
        contrast  = paste0(case, "_vs_", ctrl),
        log2FC    = if (is.null(res)) NA_real_ else res[g, "log2FoldChange"],
        padj      = if (is.null(res)) NA_real_ else res[g, "padj"],
        stringsAsFactors = FALSE
      )
    }
  }
  
  ## --- 8.7 取 normalized counts，准备画箱线图 ---
  ## DESeq2 的 normalized counts：counts(dds, normalized=TRUE)
  norm_mat <- tryCatch(
    counts(dds, normalized = TRUE),
    error = function(e) NULL
  )
  if (is.null(norm_mat)) {
    message("⚠ 取 normalized counts 失败，跳过画图 | celltype=", ct)
    next
  }
  
  ## --- 8.8 对 SUMO1/2/3 画图（可选输出全部）---
  for (g in sumo_genes) {
    
    ## 组装绘图数据：每个 donor 一行
    dfp <- data.frame(
      sample = colnames(norm_mat),
      value  = as.numeric(norm_mat[g, colnames(norm_mat)]),
      stringsAsFactors = FALSE
    )
    dfp$group4    <- pb_annot[dfp$sample, "group4"]
    dfp$celltype6 <- ct
    
    ## subtitle 三行（AD/FTD/PSP vs Control）
    sub3 <- paste(
      fmt_one(dds, contrasts[[1]], g),
      fmt_one(dds, contrasts[[2]], g),
      fmt_one(dds, contrasts[[3]], g),
      sep = "\n"
    )
    
    p <- pretty_box(
      dfp,
      title    = sprintf("%s in %s (pseudo-bulk per donor)", g, ct),
      subtitle = sub3
    )
    
    ## 输出文件名（保证可读）
    fn <- sprintf("BOX_%s__%s__4groups_byDonor.png", g, gsub("[ /]", "_", ct))
    fp <- file.path(dir_all3, fn)
    
    ok <- save_png_atomic(fp, p, width = 8.4, height = 5.6, dpi = 320, retries = 5)
    
    if (!ok) message("⚠ 写图失败（已重试）: ", fp)
  }
}

## ========= 9) 汇总统计表 =========
stats_df <- bind_rows(stats_list)

## 保存每个 SUMO × celltype × contrast 的统计
stats_fp <- file.path(out_dir, "SUMO_stats_all.csv")
write.csv(stats_df, stats_fp, row.names = FALSE)
cat("\n✅ 已保存：", stats_fp, "\n")

## ========= 10) 自动选择最稳定 SUMO（综合评分） =========
## 评分思路：
## - 对每个 gene：把所有 celltype×contrast 的结果合并
## - padj 越大越稳定：我们用 -log10(padj) 的反方向不太直观
##   这里直接做一个简单稳定性分：
##     score = median(padj, na.rm=TRUE)  -  median(|log2FC|, na.rm=TRUE)
##   -> padj 大加分，|log2FC| 大扣分
## 你后续如果想更严格也可以加“最差情况”惩罚（min padj / max |log2FC|）
rank_df <- stats_df %>%
  mutate(absL2 = abs(log2FC)) %>%
  group_by(gene) %>%
  summarise(
    n_total = n(),
    n_NA    = sum(is.na(padj) | is.na(log2FC)),
    padj_median = suppressWarnings(median(padj, na.rm = TRUE)),
    absL2_median = suppressWarnings(median(absL2, na.rm = TRUE)),
    padj_min = suppressWarnings(min(padj, na.rm = TRUE)),
    absL2_max = suppressWarnings(max(absL2, na.rm = TRUE)),
    score = padj_median - absL2_median,
    .groups = "drop"
  ) %>%
  arrange(desc(score), desc(padj_median), absL2_median)

rank_fp <- file.path(out_dir, "SUMO_stability_ranking.csv")
write.csv(rank_df, rank_fp, row.names = FALSE)
cat("✅ 已保存：", rank_fp, "\n")

best_sumo <- rank_df$gene[1]
cat("🏆 最稳定 SUMO（综合评分最高）=", best_sumo, "\n")

## ========= 11) 仅复制“最佳 SUMO”的箱线图到 SUMO_REF_selected =========
## 说明：上面我们把 SUMO1/2/3 全部画到 dir_all3
##      这里把 best_sumo 对应的 6 个 celltype 图复制到 dir_sel
all_png <- list.files(dir_all3, pattern = "\\.png$", full.names = TRUE)

best_png <- all_png[grepl(paste0("BOX_", best_sumo, "__"), basename(all_png))]
if (length(best_png) == 0) {
  message("⚠ 没找到最佳 SUMO 的图，可能前面绘图失败：", best_sumo)
} else {
  file.copy(best_png, file.path(dir_sel, basename(best_png)), overwrite = TRUE)
  cat("✅ 已复制最佳 SUMO 图到：", dir_sel, "\n")
}

cat("\n🎉 NO11 完成。\n输出目录：", out_dir, "\n")
###############################################################################

































#整体水平的箱线图
###############################################################################
## NO12_WholeCell_UBL3_4groups_perDonor_CP10k_Wilcoxon.R
##
## 【本章目标】
## A) 全细胞总体水平（不分 celltype）评估 UBL3 表达趋势（4组：AD/FTD/PSP/Control）
##    - 表达度量：counts-based 的 UBL3_log1p(CP10k)（与 UMAP 一致）
##    - 每个点 = 1 个 donor（autopsy_id）
##    - donor 层面均值：先计算每个细胞的 CP10k，再按 donor 求 mean
##    - 统计检验：Wilcoxon 秩和检验（疾病 vs Control）
##    - 效应量：log2FC 基于 donor 层面的线性 CP10k 均值
##
## B) 同章继续：SUMO1 作为内参候选（4组），按 celltype6 分 6 张 donor 箱线图
##    - 同样用 counts-based 的 log1p(CP10k)
##
## 【输入】
## - Seurat 对象（你 NO3/NO2 里生成的那个 slim 对象或 celltype6 对象）
##   必须包含：
##     meta.data: group4, autopsy_id, celltype6
##     RNA counts layer（counts）
##
## 【输出】
## - NO12_WholeCell_4groups_perDonor/
##     ├─ UBL3_wholecell_perDonor_CP10k_mean.csv
##     ├─ UBL3_wholecell_stats_3contrasts.csv
##     ├─ Fig_NO12_UBL3_wholecell_perDonor_CP10k_boxplot_4groups.png
##     └─ SUMO1_celltype6_perDonor/
##          ├─ BOX_SUMO1__Astrocytes__4groups_byDonor.png
##          └─ ...（6张）
###############################################################################

## ========= 0) 基础设置 =========
SEED <- 20251023; set.seed(SEED)
Sys.setenv(LANG = "en")
options(stringsAsFactors = FALSE)

## ✅ 不要在脚本里 install.packages()（会干扰图形设备 & 文件写入）
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(ragg)     # 用于稳定输出 PNG
})

## ========= 1) 路径（只改这里） =========
res_dir <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3"

## 你常用的对象路径（按你现有文件改一个能读到的）
## ✅ 推荐用你“轻量 slim”对象（不会爆内存）
cand_obj <- c(
  file.path(res_dir, "stepH_slim_uncompressed.rds"),
  file.path(res_dir, "NO2_step7_obj_celltype6_UBL3_named.rds"),
  file.path(res_dir, "stepH_obj_celltype6_named.rds")
)
obj_fp <- cand_obj[file.exists(cand_obj)][1]
stopifnot(length(obj_fp) == 1)

## 输出目录
out_dir <- file.path(res_dir, "NO12_WholeCell_4groups_perDonor")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ========= 2) 读入对象 + 最小自检 =========
obj <- readRDS(obj_fp)
DefaultAssay(obj) <- "RNA"
message("Loaded: ", basename(obj_fp), " | cells=", ncol(obj))

## 必备列检查（缺一不可）
need_cols <- c("group4","autopsy_id","celltype6")
miss_cols <- need_cols[!need_cols %in% colnames(obj@meta.data)]
if (length(miss_cols) > 0) {
  stop("❌ meta.data 缺少列：", paste(miss_cols, collapse=", "),
       "\n请确认你读入的是带 group4/autopsy_id/celltype6 的对象。")
}

## 组别顺序固定（Control 做 reference）
group_order_plot  <- c("AD","FTD","PSP","Control")
group_order_stats <- c("Control","AD","FTD","PSP")

obj$group4 <- factor(as.character(obj$group4), levels = group_order_stats)

## 快速核对每组 donor 数
donor_tab <- obj@meta.data %>%
  distinct(autopsy_id, group4) %>%
  count(group4) %>%
  tidyr::complete(group4 = factor(group_order_stats, levels = group_order_stats), fill = list(n=0))
print(donor_tab)

## ========= 3) 取 counts-based UBL3_log1p(CP10k)（与 UMAP 一致） =========
## 关键点：
## - 必须用 RNA@counts（原始计数）
## - CP10k：每个细胞 counts / library_size * 10000
## - 再 log1p 得到 UBL3_log1p(CP10k)

rna <- obj[["RNA"]]

## 确保 counts layer 存在
if (!"counts" %in% SeuratObject::Layers(rna)) {
  ## Seurat v5 有时 layers 未 join，做一次 JoinLayers
  rna <- SeuratObject::JoinLayers(rna)
  obj[["RNA"]] <- rna
}

mat_counts <- SeuratObject::LayerData(obj[["RNA"]], layer = "counts")  # dgCMatrix

## 找 UBL3 行（SYMBOL 或 ENSG）
target <- if ("UBL3" %in% rownames(mat_counts)) "UBL3" else
  if ("ENSG00000122042" %in% rownames(mat_counts)) "ENSG00000122042" else NULL
stopifnot(!is.null(target))

## 每个细胞总 counts（library size）
lib_size <- Matrix::colSums(mat_counts)
stopifnot(length(lib_size) == ncol(obj))

## UBL3 counts（每细胞）
ubl3_counts <- as.numeric(mat_counts[target, ])
stopifnot(length(ubl3_counts) == ncol(obj))

## 线性 CP10k（避免除0）
cp10k <- (ubl3_counts / pmax(lib_size, 1)) * 10000

## log1p(CP10k)（论文定义/与你 UMAP 一致）
ubl3_log1p_cp10k <- log1p(cp10k)

## ========= 4) donor 层面汇总：每个 donor=1点 =========
df_cells <- data.frame(
  autopsy_id = as.character(obj$autopsy_id),
  group4     = as.character(obj$group4),
  ubl3_cp10k = cp10k,                # 线性（用于 log2FC）
  ubl3_log   = ubl3_log1p_cp10k,      # log1p（用于“趋势/描述”，但统计按论文是 donor mean 后检验）
  stringsAsFactors = FALSE
)

## 去掉缺失 donor / group 的细胞
df_cells <- df_cells %>% filter(!is.na(autopsy_id), !is.na(group4))

## 每个 donor 求均值：注意这里同时保存线性CP10k均值 & log1p均值
df_donor <- df_cells %>%
  group_by(autopsy_id, group4) %>%
  summarise(
    n_cells = dplyr::n(),
    mean_cp10k = mean(ubl3_cp10k, na.rm = TRUE),   # ✅ log2FC 用它
    mean_log1p_cp10k = mean(ubl3_log, na.rm = TRUE),
    .groups = "drop"
  )

## 画图用的组顺序（Control 最右）
df_donor$group4_plot <- factor(df_donor$group4, levels = group_order_plot)

## 保存 donor 表（留档）
write.csv(df_donor,
          file.path(out_dir, "UBL3_wholecell_perDonor_CP10k_mean.csv"),
          row.names = FALSE)

message("✅ donor points = ", nrow(df_donor))
print(table(df_donor$group4))

## ========= 5) 统计：3个疾病 vs Control 的 Wilcoxon + log2FC =========
## 说明：
## - Wilcoxon：用 donor 的 mean_log1p_cp10k（按论文“同一表达度量”）
## - log2FC：用 donor 的 mean_cp10k（线性均值）计算 log2(meanDisease/meanControl)

calc_one <- function(disease) {
  d1 <- df_donor %>% filter(group4 == disease)
  d0 <- df_donor %>% filter(group4 == "Control")
  
  ## 样本太少就返回 NA（避免 wilcox 报错）
  if (nrow(d1) < 2 || nrow(d0) < 2) {
    return(data.frame(
      contrast = paste0(disease, "_vs_Control"),
      n_disease = nrow(d1),
      n_control = nrow(d0),
      mean_cp10k_disease = mean(d1$mean_cp10k, na.rm=TRUE),
      mean_cp10k_control = mean(d0$mean_cp10k, na.rm=TRUE),
      log2FC = NA_real_,
      p_wilcox = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  mu1 <- mean(d1$mean_cp10k, na.rm=TRUE)
  mu0 <- mean(d0$mean_cp10k, na.rm=TRUE)
  log2FC <- log2((mu1 + 1e-8) / (mu0 + 1e-8))
  
  pval <- wilcox.test(d1$mean_log1p_cp10k, d0$mean_log1p_cp10k, exact = FALSE)$p.value
  
  data.frame(
    contrast = paste0(disease, "_vs_Control"),
    n_disease = nrow(d1),
    n_control = nrow(d0),
    mean_cp10k_disease = mu1,
    mean_cp10k_control = mu0,
    log2FC = log2FC,
    p_wilcox = pval,
    stringsAsFactors = FALSE
  )
}

stats_3 <- bind_rows(lapply(c("AD","FTD","PSP"), calc_one))
write.csv(stats_3, file.path(out_dir, "UBL3_wholecell_stats_3contrasts.csv"), row.names = FALSE)
print(stats_3)

## 为图的 subtitle 组装三行文本
fmt_p <- function(p) if (is.na(p)) "NA" else formatC(p, format="e", digits=2)
sub3 <- paste(
  sprintf("AD vs Control:  log2FC=%s, p=%s",  sprintf("%.3f", stats_3$log2FC[stats_3$contrast=="AD_vs_Control"]),  fmt_p(stats_3$p_wilcox[stats_3$contrast=="AD_vs_Control"])),
  sprintf("FTD vs Control: log2FC=%s, p=%s",  sprintf("%.3f", stats_3$log2FC[stats_3$contrast=="FTD_vs_Control"]), fmt_p(stats_3$p_wilcox[stats_3$contrast=="FTD_vs_Control"])),
  sprintf("PSP vs Control: log2FC=%s, p=%s",  sprintf("%.3f", stats_3$log2FC[stats_3$contrast=="PSP_vs_Control"]), fmt_p(stats_3$p_wilcox[stats_3$contrast=="PSP_vs_Control"])),
  sep = "\n"
)

## ========= 6) 画图：4组 donor 箱线图（不画白色菱形均值点） =========
pal4_plot <- c(
  AD      = "#D24B40",
  FTD     = "#009E73",
  PSP     = "#CC79A7",
  Control = "#2C7FB8"
)

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

p1 <- ggplot(df_donor, aes(x = group4_plot, y = mean_log1p_cp10k, fill = group4_plot)) +
  geom_boxplot(width=0.55, outlier.shape=NA, linewidth=1.0,
               alpha=0.96, colour="grey15", median.linewidth=1.6) +
  geom_point(position = position_jitter(width=0.10, height=0),
             size=2.6, alpha=0.9, shape=21, stroke=0.5, colour="grey10") +
  scale_fill_manual(values = pal4_plot) +
  labs(
    title = "UBL3 expression per donor (Whole cells)",
    subtitle = sub3,
    x = NULL,
    y = "Mean UBL3 log1p(CP10k)"
  ) +
  theme_sci +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10)))

## ========= 7) 稳定保存 PNG（避免黑图：写到 temp 再 copy） =========
save_png_atomic <- function(filename, plot, width=7.8, height=5.6, dpi=450, retries=5) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  
  tmpfile <- file.path(
    tempdir(),
    paste0("tmp_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", sample.int(1e9, 1), ".png")
  )
  
  for (k in seq_len(retries)) {
    if (file.exists(tmpfile)) file.remove(tmpfile)
    
    ok <- FALSE
    try({
      ragg::agg_png(tmpfile, width=width, height=height, units="in", res=dpi, background="white")
      print(plot)
      dev.off()
      
      if (file.exists(tmpfile) && is.finite(file.info(tmpfile)$size) && file.info(tmpfile)$size > 80*1024) ok <- TRUE
    }, silent = TRUE)
    
    if (ok) {
      if (file.exists(filename)) file.remove(filename)
      ok2 <- file.copy(tmpfile, filename, overwrite = TRUE)
      if (isTRUE(ok2) && file.exists(filename) && file.info(filename)$size > 80*1024) {
        file.remove(tmpfile)
        gc(FALSE)
        return(TRUE)
      }
    }
    Sys.sleep(0.3)
  }
  
  if (file.exists(tmpfile)) file.remove(tmpfile)
  gc(FALSE)
  return(FALSE)
}

fig1_fp <- file.path(out_dir, "Fig_NO12_UBL3_wholecell_perDonor_log1pCP10k_boxplot_4groups.png")
ok <- save_png_atomic(fig1_fp, p1, width=7.8, height=5.6, dpi=450, retries=5)
if (!ok) message("⚠ 写图失败：", fig1_fp) else message("✅ saved: ", fig1_fp)

###############################################################################
## Part B：同章继续做 SUMO1（4组）按 celltype6 的 donor 箱线图（counts-based log1p(CP10k)）
###############################################################################

###############################################################################
## NO12B_WholeCell_SUMO1_4groups_perDonor_CP10k_Wilcoxon.R
##
## 【本脚本目的】
##  - 在“全细胞总体水平（不分 celltype）”比较 SUMO1 的表达稳定性
##  - 采用与 UMAP 完全一致的 counts-based 指标：SUMO1_log1p(CP10k)
##  - 每个点 = 1 个 donor（autopsy_id）
##  - 统计检验：Wilcoxon 秩和检验（疾病组 vs Control）
##  - 效应量 log2FC：基于 donor 层面的线性 CP10k 均值计算：
##        log2FC = log2(mean_CP10k_disease / mean_CP10k_control)
##
## 【为什么这么做？】
##  - SUMO1 作为内参候选，应该在不同疾病组之间“尽可能不变”
##  - 你论文方法写的是 counts-based log1p(CP10k)，所以图的 Y 轴也必须一致
##  - log2FC 用线性 CP10k（避免在 log 尺度上直接算 FC 的歧义）
##
## 【输入要求】
##  - Seurat 对象中必须有：
##      meta.data: group4, autopsy_id
##      RNA assay: counts layer（原始计数）
##
## 【输出】
##  - 目录：<res_dir>/NO12_WholeCell_4groups_perDonor/
##      SUMO1_wholecell_perDonor_CP10k_mean.csv
##      SUMO1_wholecell_stats_3contrasts.csv
##      Fig_NO12_SUMO1_wholecell_perDonor_log1pCP10k_boxplot_4groups.png
###############################################################################

## ========= 0) 基础设置 =========
SEED <- 20251023; set.seed(SEED)
Sys.setenv(LANG = "en")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(ragg)    # ✅ 稳定写 PNG（避免偶发黑图）
})

## ========= 1) 路径（只改这里） =========
res_dir <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3"

## 你的对象文件候选（选第一个存在的）
cand_obj <- c(
  file.path(res_dir, "stepH_slim_uncompressed.rds"),
  file.path(res_dir, "NO2_step7_obj_celltype6_UBL3_named.rds"),
  file.path(res_dir, "stepH_obj_celltype6_named.rds")
)
obj_fp <- cand_obj[file.exists(cand_obj)][1]
stopifnot(length(obj_fp) == 1)

## 输出目录：与 NO12 同一个总目录（保证章节一致）
out_dir <- file.path(res_dir, "NO12_WholeCell_4groups_perDonor")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ========= 2) 读入对象 + 最小自检 =========
obj <- readRDS(obj_fp)
DefaultAssay(obj) <- "RNA"
message("Loaded: ", basename(obj_fp), " | cells=", ncol(obj))

## 必备列检查
need_cols <- c("group4","autopsy_id")
miss_cols <- need_cols[!need_cols %in% colnames(obj@meta.data)]
if (length(miss_cols) > 0) {
  stop("❌ meta.data 缺少列：", paste(miss_cols, collapse=", "),
       "\n请确认你读入的是带 group4/autopsy_id 的对象。")
}

## 组别顺序固定：Control 做 reference（用于统计）
group_order_stats <- c("Control","AD","FTD","PSP")
## 图上顺序固定：Control 最右（与你前面的 UBL3/NO11 一致）
group_order_plot  <- c("AD","FTD","PSP","Control")

obj$group4 <- factor(as.character(obj$group4), levels = group_order_stats)

## donor 数快速核对（每个 donor 只计一次）
donor_tab <- obj@meta.data %>%
  distinct(autopsy_id, group4) %>%
  count(group4)
message("Donor counts per group:")
print(donor_tab)

## ========= 3) 取 counts-based SUMO1_log1p(CP10k) =========
## 核心逻辑：
## - 使用 RNA@counts（原始计数）
## - 每个细胞总 counts = library size
## - SUMO1_CP10k = SUMO1_counts / library_size * 10000
## - SUMO1_log1p(CP10k) = log1p(SUMO1_CP10k)

rna <- obj[["RNA"]]

## 确保 counts layer 存在（Seurat v5 有时需要 JoinLayers）
if (!"counts" %in% SeuratObject::Layers(rna)) {
  rna <- SeuratObject::JoinLayers(rna)
  obj[["RNA"]] <- rna
}

mat_counts <- SeuratObject::LayerData(obj[["RNA"]], layer = "counts")  # dgCMatrix

## 找 SUMO1 行（SYMBOL）
target <- if ("SUMO1" %in% rownames(mat_counts)) "SUMO1" else NULL
stopifnot(!is.null(target))

## 每个细胞总 counts（library size）
lib_size <- Matrix::colSums(mat_counts)
stopifnot(length(lib_size) == ncol(obj))

## SUMO1 counts（每细胞）
sumo1_counts <- as.numeric(mat_counts[target, ])
stopifnot(length(sumo1_counts) == ncol(obj))

## 线性 CP10k（避免除 0）
sumo1_cp10k <- (sumo1_counts / pmax(lib_size, 1)) * 10000

## log1p(CP10k)（用于图 & Wilcoxon）
sumo1_log1p_cp10k <- log1p(sumo1_cp10k)

## ========= 4) donor 层面汇总：每个 donor=1点 =========
df_cells <- data.frame(
  autopsy_id = as.character(obj$autopsy_id),
  group4     = as.character(obj$group4),
  sumo1_cp10k = sumo1_cp10k,               # 线性：用于 log2FC
  sumo1_log   = sumo1_log1p_cp10k,         # log：用于图 & Wilcoxon
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(autopsy_id), !is.na(group4))

## donor 层面均值：同时保存线性均值 & log均值
df_donor <- df_cells %>%
  group_by(autopsy_id, group4) %>%
  summarise(
    n_cells = dplyr::n(),
    mean_cp10k = mean(sumo1_cp10k, na.rm = TRUE),
    mean_log1p_cp10k = mean(sumo1_log, na.rm = TRUE),
    .groups = "drop"
  )

## 图上组顺序（Control 最右）
df_donor$group4_plot <- factor(df_donor$group4, levels = group_order_plot)

## 保存 donor 表（留档可复现）
write.csv(df_donor,
          file.path(out_dir, "SUMO1_wholecell_perDonor_CP10k_mean.csv"),
          row.names = FALSE)

message("✅ donor points = ", nrow(df_donor))
print(table(df_donor$group4))

## ========= 5) 统计：3个疾病 vs Control（Wilcoxon + log2FC） =========
## - Wilcoxon：比较 donor 的 mean_log1p_cp10k（非参数）
## - log2FC：用 donor 的 mean_cp10k（线性均值）算 log2(meanDisease/meanControl)

calc_one <- function(disease) {
  d1 <- df_donor %>% filter(group4 == disease)
  d0 <- df_donor %>% filter(group4 == "Control")
  
  ## donor 太少时返回 NA（避免 wilcox 报错）
  if (nrow(d1) < 2 || nrow(d0) < 2) {
    return(data.frame(
      contrast = paste0(disease, "_vs_Control"),
      n_disease = nrow(d1),
      n_control = nrow(d0),
      mean_cp10k_disease = mean(d1$mean_cp10k, na.rm=TRUE),
      mean_cp10k_control = mean(d0$mean_cp10k, na.rm=TRUE),
      log2FC = NA_real_,
      p_wilcox = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  mu1 <- mean(d1$mean_cp10k, na.rm=TRUE)
  mu0 <- mean(d0$mean_cp10k, na.rm=TRUE)
  log2FC <- log2((mu1 + 1e-8) / (mu0 + 1e-8))
  
  pval <- wilcox.test(d1$mean_log1p_cp10k, d0$mean_log1p_cp10k, exact = FALSE)$p.value
  
  data.frame(
    contrast = paste0(disease, "_vs_Control"),
    n_disease = nrow(d1),
    n_control = nrow(d0),
    mean_cp10k_disease = mu1,
    mean_cp10k_control = mu0,
    log2FC = log2FC,
    p_wilcox = pval,
    stringsAsFactors = FALSE
  )
}

stats_3 <- bind_rows(lapply(c("AD","FTD","PSP"), calc_one))
write.csv(stats_3,
          file.path(out_dir, "SUMO1_wholecell_stats_3contrasts.csv"),
          row.names = FALSE)
print(stats_3)

## subtitle 三行文本（用于图上展示）
fmt_p <- function(p) if (is.na(p)) "NA" else formatC(p, format="e", digits=2)
get_row <- function(con) stats_3[stats_3$contrast == con, , drop=FALSE]

r1 <- get_row("AD_vs_Control")
r2 <- get_row("FTD_vs_Control")
r3 <- get_row("PSP_vs_Control")

sub3 <- paste(
  sprintf("AD vs Control:  log2FC=%s, p=%s",  sprintf("%.3f", r1$log2FC), fmt_p(r1$p_wilcox)),
  sprintf("FTD vs Control: log2FC=%s, p=%s",  sprintf("%.3f", r2$log2FC), fmt_p(r2$p_wilcox)),
  sprintf("PSP vs Control: log2FC=%s, p=%s",  sprintf("%.3f", r3$log2FC), fmt_p(r3$p_wilcox)),
  sep = "\n"
)

## ========= 6) 画图：4组 donor 箱线图（不画白色菱形均值点） =========
pal4_plot <- c(
  AD      = "#D24B40",
  FTD     = "#009E73",
  PSP     = "#CC79A7",
  Control = "#2C7FB8"
)

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

p <- ggplot(df_donor, aes(x = group4_plot, y = mean_log1p_cp10k, fill = group4_plot)) +
  geom_boxplot(width=0.55, outlier.shape=NA, linewidth=1.0,
               alpha=0.96, colour="grey15", median.linewidth=1.6) +
  geom_point(position = position_jitter(width=0.10, height=0),
             size=2.6, alpha=0.9, shape=21, stroke=0.5, colour="grey10") +
  scale_fill_manual(values = pal4_plot) +
  labs(
    title = "SUMO1 expression per donor (Whole cells)",
    subtitle = sub3,
    x = NULL,
    y = "Mean SUMO1 log1p(CP10k) "
  ) +
  theme_sci +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10)))

## ========= 7) 稳定保存 PNG（避免偶发黑图：写到 temp 再 copy） =========
save_png_atomic <- function(filename, plot, width=7.8, height=5.6, dpi=450, retries=5) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  
  tmpfile <- file.path(
    tempdir(),
    paste0("tmp_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", sample.int(1e9, 1), ".png")
  )
  
  for (k in seq_len(retries)) {
    if (file.exists(tmpfile)) file.remove(tmpfile)
    
    ok <- FALSE
    try({
      ragg::agg_png(tmpfile, width=width, height=height, units="in", res=dpi, background="white")
      print(plot)
      dev.off()
      
      ## 80KB 阈值：防止半张图/黑图误判为成功
      if (file.exists(tmpfile) && is.finite(file.info(tmpfile)$size) && file.info(tmpfile)$size > 80*1024) ok <- TRUE
    }, silent = TRUE)
    
    if (ok) {
      if (file.exists(filename)) file.remove(filename)
      ok2 <- file.copy(tmpfile, filename, overwrite = TRUE)
      if (isTRUE(ok2) && file.exists(filename) && file.info(filename)$size > 80*1024) {
        file.remove(tmpfile)
        gc(FALSE)
        return(TRUE)
      }
    }
    Sys.sleep(0.3)
  }
  
  if (file.exists(tmpfile)) file.remove(tmpfile)
  gc(FALSE)
  return(FALSE)
}

fig_fp <- file.path(out_dir, "Fig_NO12_SUMO1_wholecell_perDonor_log1pCP10k_boxplot_4groups.png")
ok <- save_png_atomic(fig_fp, p, width=7.8, height=5.6, dpi=450, retries=5)
if (!ok) message("⚠ 写图失败：", fig_fp) else message("✅ saved: ", fig_fp)

message("\n🎉 NO12B（SUMO1 whole-cell 4groups）完成。输出：", out_dir)
###############################################################################















#
#
## ===========================
## B. 生成投稿核对表（QC后/当前对象）
mv_fp <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3/QC_tables_for_submission/Sample_Group4_MajorityVote.csv"
mv <- read.csv(mv_fp, stringsAsFactors = FALSE)

cat("Total unique samples in MajorityVote:", length(unique(mv$sample)), "\n\n")

cat("Correct sample counts by diagnosis (majority vote):\n")
print(table(mv$group4_majority))

cat("\nThese are the 4 mixed-label samples (flag_multi_group=TRUE):\n")
print(mv[mv$flag_multi_group == TRUE, ])


rc_fp <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3/QC_tables_for_submission/Sample_Region_Counts.csv"
rc <- read.csv(rc_fp, stringsAsFactors = FALSE)

cat("Region -> number of samples with cells (cells_n > 0):\n")
region_sample_n <- aggregate(sample ~ region, data = rc[rc$cells_n > 0, ],
                             FUN = function(x) length(unique(x)))
colnames(region_sample_n)[2] <- "n_samples_with_cells"
print(region_sample_n)

cat("\nRegion -> total cells:\n")
region_cells <- aggregate(cells_n ~ region, data = rc, sum)
print(region_cells)

cat("\nHow many regions each sample covers (cells_n > 0):\n")
sample_region_n <- aggregate(region ~ sample, data = rc[rc$cells_n > 0, ],
                             FUN = function(x) length(unique(x)))
colnames(sample_region_n)[2] <- "n_regions_with_cells"
print(table(sample_region_n$n_regions_with_cells))
