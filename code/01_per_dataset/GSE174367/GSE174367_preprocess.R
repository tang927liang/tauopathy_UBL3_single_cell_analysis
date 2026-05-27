###############################################################################
# 01_per_dataset / GSE174367  (Alzheimer's disease, prefrontal cortex)
# From-raw snRNA-seq preprocessing: builds the harmonized 6-cell-type Seurat
# object consumed by code/02_integrated_figures/.
#
# Source     : GEO accession GSE174367 (Cell Ranger aggr filtered_feature_bc
#              matrix .h5 + snRNA-seq cell metadata CSV).
# Pipeline   : build object from h5 + meta -> QC (nFeature_RNA > 200,
#              nCount_RNA < 20000, percent.mt < 20) -> CCA integration
#              (dims = 1:20, nfeatures = 2000) -> PCA -> clustering + UMAP
#              (dims = 1:20) -> cluster-marker annotation collapsed to 6 cell
#              types (same template as GSE157827).
# Produces   : stepH_obj_celltype6_named.rds  (the object read by the integrated
#              figure/stat scripts; also underpins Supplementary Figure S2).
# Paths      : all paths point to the data-project location on disk. Raw inputs
#              are public at GSE174367; intermediate .rds objects live with the
#              data, not in this repository. Confirm paths before running.
# Environment: R 4.4.3 + Seurat 4.3.0 + SeuratObject 5.2.0 + Matrix 1.7-4 +
#              data.table 1.17.8.
#
# NOTE: kept faithful to the script as run. Sections after the celltype6
#       hand-off (per-dataset cell-type / UBL3 UMAPs) overlap Figure 2 /
#       Supplementary Figure S3 and are retained only for provenance.
###############################################################################

#GSE174367
#先读取文件，看格式

# ## 设置路径为你存放 GSE174367 的文件夹
setwd("D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367")

## 读取 meta 数据（R 会自动识别 .gz）
meta <- read.csv("GSE174367_snRNA-seq_cell_meta.csv.gz")

## 查看数据基本情况
dim(meta)
head(meta)

## 统计每个样本的诊断类型（AD / Control）
table(meta$Diagnosis, meta$SampleID)

## 统计每个 SampleID 的诊断（确认 AD 11 vs Control 7）
sample_diag <- tapply(meta$Diagnosis, meta$SampleID, unique)
table(sample_diag)


















#一，合并原始数据和meta信息得到含所有样本和分组信息的rds大文件
###############################################################################
# GSE174367：从头复现 GSE157827 风格的 snRNA 分析流水线
# 环境：R 4.4.3 + Seurat 4.3.0 + SeuratObject 5.2.0 + Matrix 1.7-4 + data.table 1.17.8
# 目录：D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367
###############################################################################

## ========== 第 0 章：基础设置（固定随机种子 + 路径 + 环境记录） ==========

# 统一随机种子（整套分析都用这个）
SEED <- 20251023
set.seed(SEED)

# 路径（请确认和你的文件夹一致）
base_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367"
raw_dir  <- base_dir
res_dir  <- file.path(base_dir, "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# 加载必要 R 包
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

# 小工具：记录当前环境信息（R / Seurat / Matrix 等），方便以后复现
log_env <- function(res_dir, step_name = "step_env") {
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
  log_file <- file.path(res_dir, paste0(step_name, "_sessionInfo.txt"))
  SEED <- 20251023
  set.seed(SEED)
  sink(log_file)
  cat("SEED =", SEED, "\n")
  print(sessionInfo())
  sink()
  message("✅ 环境信息已保存到: ", log_file)
}

# 记录一次基础环境
log_env(res_dir, step_name = "step0_setup")


## ========== 第 1 章：读入 h5 + meta → 合并成大 Seurat 对象 ==========

# 1.1 读入 cellranger aggr 的 filtered_feature_bc_matrix.h5 -----------------
h5_file <- file.path(raw_dir, "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5")
mat <- Read10X_h5(h5_file)
cat("表达矩阵维度（genes x cells）：", nrow(mat), "x", ncol(mat), "\n")

# 1.2 创建最原始的 Seurat 对象（不额外过滤）------------------------------
obj <- CreateSeuratObject(
  counts       = mat,
  project      = "GSE174367",
  min.cells    = 0,
  min.features = 0
)
cat("初始 Seurat 对象：", nrow(obj), "genes x", ncol(obj), "cells\n")

# 1.3 读入 cell-level metadata ---------------------------------------------
meta_file <- file.path(raw_dir, "GSE174367_snRNA-seq_cell_meta.csv.gz")
meta <- fread(meta_file)
cat("meta 行列数：", nrow(meta), "x", ncol(meta), "\n")
print(head(meta))

# 1.4 修正：h5 比 meta 多 298 个细胞 → 只保留交集 ------------------------
# 确认 meta 中所有 Barcode 都在 h5 中
stopifnot(all(meta$Barcode %in% colnames(obj)))

extra_cells  <- setdiff(colnames(obj), meta$Barcode)
cat("h5 中多出的无 meta 细胞数：", length(extra_cells), "\n")
common_cells <- intersect(colnames(obj), meta$Barcode)
cat("交集细胞数：", length(common_cells), "\n")

# 只保留有 meta 的细胞
obj <- subset(obj, cells = common_cells)
cat("subset 后对象维度：", nrow(obj), "genes x", ncol(obj), "cells\n")

# 对齐 meta 行顺序
meta2 <- meta[match(colnames(obj), meta$Barcode), ]
stopifnot(all(meta2$Barcode == colnames(obj)))

# 1.5 统一细胞命名：SampleID_barcode（与 GSE157827 完全同风格）---------
new_cell_names <- paste0(meta2$SampleID, "_", colnames(obj))
stopifnot(sum(duplicated(new_cell_names)) == 0)

colnames(obj) <- new_cell_names

# meta2 去掉 Barcode，并用新细胞名作为行名
meta2$Barcode <- NULL
rownames(meta2) <- colnames(obj)
stopifnot(nrow(meta2) == ncol(obj))

# 1.6 将全部 metadata 加入 Seurat 对象 -------------------------------
obj <- AddMetaData(obj, metadata = as.data.frame(meta2))

# 创建统一分组变量 group（Control / AD）
obj$group <- factor(obj$Diagnosis, levels = c("Control", "AD"))

cat("按 group 统计细胞数：\n")
print(table(obj$group))
cat("按 SampleID × group 统计：\n")
print(table(obj$SampleID, obj$group))

# 检查重复名
cat("重复细胞名数：", sum(duplicated(colnames(obj))), "\n")
cat("重复基因名数：", sum(duplicated(rownames(obj))), "\n")

# 1.7 保存合并后的大对象（正常版 + 兼容旧 R 的 v2 版）----------------
saveRDS(obj, file.path(res_dir, "GSE174367_merged_raw_with_meta.rds"))
# 如需兼容旧 R，可再保存一份 version=2
saveRDS(
  obj,
  file.path(res_dir, "GSE174367_merged_raw_with_meta_v2.rds"),
  version = 2
)

# 输出每个样本 × 组别的细胞数
cell_summary <- as.data.frame(table(obj$SampleID, obj$group))
colnames(cell_summary) <- c("SampleID", "group", "cell_n")
write.csv(cell_summary,
          file.path(res_dir, "GSE174367_sample_cell_counts_by_group.csv"),
          row.names = FALSE)

cat("✅ 第 1 章完成：已保存 GSE174367_merged_raw_with_meta.rds\n")
log_env(res_dir, step_name = "step1_merge")














### ========== 第 2 章：QC（两步）—— 单层 counts + 过滤前/后小提琴图 ==========

# 如果是新会话，可以用上一章的 rds 重新读入：
# obj <- readRDS(file.path(res_dir, "GSE174367_merged_raw_with_meta.rds"))
DefaultAssay(obj) <- "RNA"

# 2.1 从 counts 层取出矩阵（GSE174367 本身就只有单层）------------------
mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
cat("QC 使用的 counts 维度：", nrow(mat), "genes x", ncol(mat), "cells\n")

# 2.2 计算三大 QC 指标：nFeature_RNA / nCount_RNA / percent.mt ----------
obj$nFeature_RNA <- Matrix::colSums(mat > 0)
obj$nCount_RNA   <- Matrix::colSums(mat)

# 线粒体基因（GSE174367 行名为 symbol，使用 MT- 前缀）
is_mt <- grepl("^MT-", rownames(mat))
obj$percent.mt <- Matrix::colSums(mat[is_mt, , drop = FALSE]) /
  Matrix::colSums(mat) * 100

summary(obj$percent.mt)

# 2.3 过滤前 QC 小提琴图 ----------------------------------------------
png(file.path(res_dir, "stepB_QC_violin_before_filter.png"),
    width = 1600, height = 600)
print(
  VlnPlot(obj,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          ncol = 3, pt.size = 0, raster = TRUE)
)
dev.off()

writeLines(
  sprintf("Before filter: genes=%d, cells=%d", nrow(obj), ncol(obj)),
  file.path(res_dir, "stepB_QC_sizes.txt")
)

# 2.4 按 GSE157827 / 原文风格的阈值进行二次 QC -----------------------
#   nFeature_RNA > 200
#   nCount_RNA   < 20000
#   percent.mt   < 20
keep <- (obj$nFeature_RNA > 200) &
  (obj$nCount_RNA   < 20000) &
  (obj$percent.mt   < 20)

cat("过滤后保留细胞数：", sum(keep), "/", ncol(obj), "\n")
obj_flt <- subset(obj, cells = colnames(obj)[keep])

# 2.5 过滤后 QC 小提琴图 ----------------------------------------------
png(file.path(res_dir, "stepC_QC_violin_after_filter.png"),
    width = 1600, height = 600)
print(
  VlnPlot(obj_flt,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          ncol = 3, pt.size = 0, raster = TRUE)
)
dev.off()

# 保存过滤后对象（对应 stepC_filtered_obj）
saveRDS(obj_flt, file.path(res_dir, "stepC_filtered_obj.rds"))

writeLines(
  sprintf("After filter: genes=%d, cells=%d",
          nrow(obj_flt), ncol(obj_flt)),
  file.path(res_dir, "stepC_QC_sizes.txt")
)

# 导出过滤前后样本/组别细胞数量
write.csv(as.data.frame(table(obj$SampleID)),
          file.path(res_dir, "stepC_cells_by_sample_before.csv"),
          row.names = FALSE)
write.csv(as.data.frame(table(obj_flt$SampleID)),
          file.path(res_dir, "stepC_cells_by_sample_after.csv"),
          row.names = FALSE)
write.csv(as.data.frame(table(obj$group)),
          file.path(res_dir, "stepC_cells_by_group_before.csv"),
          row.names = FALSE)
write.csv(as.data.frame(table(obj_flt$group)),
          file.path(res_dir, "stepC_cells_by_group_after.csv"),
          row.names = FALSE)

cat("✅ 第 2 章完成：QC 过滤后对象已保存为 stepC_filtered_obj.rds\n")
log_env(res_dir, step_name = "step2_QC")















###############################################################################
# GSE174367: 对齐 GSE157827 第 3–5 章的整合流程
#   1) 按 SampleID 拆分样本 → 每样本 Normalize + HVG(1000)
#   2) 把 11 个 AD 样本分成 AD1/AD2 两组 + 所有 Control 样本为一组
#   3) 三组内部各做一次 CCA 整合 (dims=1:20, nfeatures=2000)
#   4) 三组再一起 CCA 整合 → 得到总 integrated assay
#   5) 用 integrated assay 做 Scale → PCA(50) → JackStraw → 
#      Neighbors/Clusters(res=1) → UMAP(dims=1:20, v3 风格)
###############################################################################

## ========= 通用设置：随机种子 + 路径 + 环境日志函数 =========
SEED <- 20251023; set.seed(SEED)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

base_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367"
res_dir  <- file.path(base_dir, "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# 环境记录函数（每个大步骤开头调用一次，保存 sessionInfo）
log_env <- function(res_dir, step_name = "step_env") {
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
  log_file <- file.path(res_dir, paste0(step_name, "_sessionInfo.txt"))
  SEED <- 20251023
  set.seed(SEED)
  sink(log_file)
  cat("SEED =", SEED, "\n")
  print(sessionInfo())
  sink()
  message("✅ 环境信息已保存到: ", log_file)
}

log_env(res_dir, "stepD_split_normalize_HVG")  # 记录当前环境


###############################################################################
# 第 D 章：按 SampleID 拆分样本 → Normalize(LogNormalize) + HVG(1000)
###############################################################################

# 从“过滤后的大对象”开始（和 GSE157827 一样）
obj_flt <- readRDS(file.path(res_dir, "stepC_filtered_obj.rds"))
DefaultAssay(obj_flt) <- "RNA"

# 按 SampleID 拆分成 18 个样本（11 AD + 7 Control）
obj_list <- SplitObject(obj_flt, split.by = "SampleID")

# 对每个样本单独进行 Normalize + FindVariableFeatures(vst, 1000)
for (nm in names(obj_list)) {
  message("Processing sample: ", nm)
  obj_list[[nm]] <- NormalizeData(
    obj_list[[nm]],
    normalization.method = "LogNormalize",
    scale.factor = 1e4,
    verbose = FALSE
  )
  obj_list[[nm]] <- FindVariableFeatures(
    obj_list[[nm]],
    selection.method = "vst",
    nfeatures = 1000,
    verbose = FALSE
  )
}

# 保存拆分 + 标准化 + HVG 列表，便于后续断点续跑
saveRDS(
  obj_list,
  file.path(res_dir, "stepD_split_normalized_vst1000_GSE174367.rds")
)
cat("✅ 已保存: stepD_split_normalized_vst1000_GSE174367.rds\n")


###############################################################################
# 第 E1 章：三组内部 CCA 整合（AD1 / AD2 / Control）
###############################################################################
log_env(res_dir, "stepE_groupwise_cca")

# 重新读入刚才的列表（防止中途重启）
obj_list <- readRDS(file.path(res_dir, "stepD_split_normalized_vst1000_GSE174367.rds"))

# 查看有哪些 SampleID，确认无误
print(sort(names(obj_list)))
# 预期：
# "Sample-100" "Sample-17" "Sample-19" "Sample-22" "Sample-27" "Sample-33"
# "Sample-37" "Sample-43" "Sample-45" "Sample-46" "Sample-47" "Sample-50"
# "Sample-52" "Sample-58" "Sample-66" "Sample-82" "Sample-90" "Sample-96"

# 按 AD/Control 把 18 个样本分成 3 组（类似 GSE157827 的 AD1/AD2/NC）
g_ad1 <- c("Sample-17", "Sample-19", "Sample-22", "Sample-27", "Sample-33")
g_ad2 <- c("Sample-37", "Sample-43", "Sample-45", "Sample-46", "Sample-47", "Sample-50")
g_ctrl <- c("Sample-52", "Sample-58", "Sample-66", "Sample-82", "Sample-90", "Sample-96", "Sample-100")

# 小工具：对一组样本做 CCA 整合（参数完全照 GSE157827）
integrate_by_cca <- function(xlist, outfile, seed = SEED) {
  set.seed(seed)
  feats <- SelectIntegrationFeatures(object.list = xlist, nfeatures = 2000)
  
  set.seed(seed)
  anchors <- FindIntegrationAnchors(
    object.list     = xlist,
    anchor.features = feats,
    dims            = 1:20,
    reduction       = "cca",
    verbose         = TRUE
  )
  
  set.seed(seed)
  obj_int <- IntegrateData(anchorset = anchors, dims = 1:20)  # 得到 integrated assay
  
  saveRDS(obj_int, outfile)
  invisible(obj_int)
}

# 对三组分别做 CCA 整合并保存
obj_ad1_int <- integrate_by_cca(
  obj_list[g_ad1],
  file.path(res_dir, "stepE_integrated_AD_part1_GSE174367.rds")
)
obj_ad2_int <- integrate_by_cca(
  obj_list[g_ad2],
  file.path(res_dir, "stepE_integrated_AD_part2_GSE174367.rds")
)
obj_ctrl_int <- integrate_by_cca(
  obj_list[g_ctrl],
  file.path(res_dir, "stepE_integrated_CTRL_GSE174367.rds")
)

cat("✅ 三组内部整合完成并保存：AD_part1 / AD_part2 / CTRL\n")


###############################################################################
# 第 E2 章：把三段整合对象再整合为一个总对象（仍用 CCA，dims=1:20）
###############################################################################
###############################################################################
# 重新运行 E2：三组整合 → ALL 整合 → 保存为 version=2 → 当场测试
###############################################################################
# 重新运行 E2：三组整合 → ALL 整合 → 保存为 version=2 → 当场测试
###############################################################################

SEED <- 20251023; set.seed(SEED)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results"

# 读取 E1 的三个部分
ad1  <- readRDS(file.path(res_dir, "stepE_integrated_AD_part1_GSE174367.rds"))
ad2  <- readRDS(file.path(res_dir, "stepE_integrated_AD_part2_GSE174367.rds"))
ctrl <- readRDS(file.path(res_dir, "stepE_integrated_CTRL_GSE174367.rds"))

cat("三个部分维度：\n")
print(dim(ad1))
print(dim(ad2))
print(dim(ctrl))

# ----------- E2 开始：总整合 -----------
set.seed(SEED)
feats2 <- SelectIntegrationFeatures(object.list = list(ad1, ad2, ctrl), nfeatures = 2000)

set.seed(SEED)
anchors2 <- FindIntegrationAnchors(
  object.list     = list(ad1, ad2, ctrl),
  anchor.features = feats2,
  dims            = 1:20,
  reduction       = "cca",
  verbose         = TRUE
)

set.seed(SEED)
obj_int_all <- IntegrateData(anchorset = anchors2, dims = 1:20)

# ----------- 保存为 version=2（更稳） -----------
saveRDS(
  obj_int_all,
  file.path(res_dir, "stepE_integrated_ALL_cca_GSE174367.rds"),
  version = 2,
  compress = "xz"
)

cat("✅ 已保存：stepE_integrated_ALL_cca_GSE174367.rds (version=2)\n")

# ----------- 当场测试能否读回（最关键） -----------
obj_test <- readRDS(file.path(res_dir, "stepE_integrated_ALL_cca_GSE174367.rds"))
cat("读回整合对象的维度：\n")
print(dim(obj_test))



# 这个对象中：
#   assay "integrated" = 批次校正后的表达（用于 PCA/UMAP/聚类）
#   assay "RNA"        = 原始表达（用于找 marker / UBL3 表达）


















###############################################################################
#最开始是按照GSE157827UMAP参数设置（dims=1:20、resolution=1、n.neighbors=30、min.dist=0.3）进行PCA和聚类及UMAP
# GSE174367 - 第 F 章: PCA(50) + JackStraw + 聚类 + UMAP（Seurat v3 风格）
###############################################################################
#根据后面跑出来的stepF_JackStraw_PC_pvalues_GSE174367和stepF_ElbowPlot_50PCs_GSE174367.png），不合适再重新调上面4个参数
SEED <- 20251023; set.seed(SEED)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# 环境记录（方便以后复现）
log_env <- function(res_dir, step_name = "step_env") {
  log_file <- file.path(res_dir, paste0(step_name, "_sessionInfo.txt"))
  SEED <- 20251023
  set.seed(SEED)
  sink(log_file)
  cat("SEED =", SEED, "\n")
  print(sessionInfo())
  sink()
  message("✅ 环境信息已保存到: ", log_file)
}

log_env(res_dir, "stepF_PCA_JackStraw_UMAP")


###############################################################################
# 1. 读取 ALL 整合对象，用 integrated 做下游分析
###############################################################################

obj <- readRDS(file.path(res_dir, "stepE_integrated_ALL_cca_GSE174367.rds"))

# 确认 assay
cat("Assays in obj: ", paste(Assays(obj), collapse = ", "), "\n")
stopifnot("integrated" %in% Assays(obj), "RNA" %in% Assays(obj))

# 之后所有降维/聚类都在 integrated 上做
DefaultAssay(obj) <- "integrated"


###############################################################################
# 2. ScaleData + PCA(50)
###############################################################################

obj <- ScaleData(obj, verbose = FALSE)

set.seed(SEED)
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

# 看一下前 5 个 PC 的 top genes，确认一下有没有明显的技术噪声
print(obj[["pca"]], dims = 1:5, nfeatures = 5)




###############################################################################
# 3. JackStraw(50 PCs) + ElbowPlot
###############################################################################

library(future)
plan(sequential)  # 单进程
options(
  future.globals.maxSize = Inf,          # 不再限制 global 大小
  future.rng.onMisuse    = "ignore"      # 关闭那个随机数 warning
)

set.seed(SEED)
obj <- JackStraw(
  obj,
  reduction      = "pca",
  dims           = 50,
  num.replicate  = 50,   # 比 100 少一点，减少负担
  verbose        = FALSE
)
obj <- ScoreJackStraw(obj, dims = 1:50)

# 保存一下
saveRDS(
  obj,
  file.path(res_dir, "stepF_afterPCA_JackStraw_GSE174367.rds"),
  version = 2,
  compress = "xz"
)


# JackStraw p 值图：用于决定保留多少 PC
obj <- readRDS(file.path(res_dir, "stepF_afterPCA_JackStraw_GSE174367.rds"))

png(file.path(res_dir, "stepF_JackStraw_PC_pvalues_GSE174367.png"),
    width = 1600, height = 900)
print(JackStrawPlot(obj, dims = 1:50))
dev.off()

png(file.path(res_dir, "stepF_ElbowPlot_50PCs_GSE174367.png"),
    width = 1200, height = 800)
print(ElbowPlot(obj, ndims = 50))
dev.off()


# 保存 PCA + JackStraw 后的检查点（可选）
saveRDS(
  obj,
  file.path(res_dir, "stepF_afterPCA_JackStraw_GSE174367.rds"),
  version = 2,
  compress = "xz"
)
cat("✅ PCA(50) + JackStraw 已完成并保存检查点\n")
# JackStraw p 值图
png(file.path(res_dir, "stepF_JackStraw_PC_pvalues_GSE174367.png"),
    width = 1600, height = 900)
print(JackStrawPlot(obj, dims = 1:50))
dev.off()

# ElbowPlot
png(file.path(res_dir, "stepF_ElbowPlot_50PCs_GSE174367.png"),
    width = 1200, height = 800)
print(ElbowPlot(obj, ndims = 50))
dev.off()






###############################################################################
# 4. 选择 PC 数：dims_use <- 1:20 （先照 GSE157827，一会儿看图再决定要不要改）
###############################################################################

dims_use <- 1:20   # **关键参数 1：PC 数量**


###############################################################################
# 5. 构图（Neighbors）、聚类（resolution=1）、UMAP（v3 风格参数）
###############################################################################

# 5.1 Neighbors 图
set.seed(SEED)
obj <- FindNeighbors(
  obj,
  dims    = dims_use,
  k.param = 20,          # 邻居数量（图结构）——通常不需要改
  verbose = FALSE
)

# 5.2 图聚类（Louvain；resolution 控制簇数量）
set.seed(SEED)
obj <- FindClusters(
  obj,
  resolution = 1,        # **关键参数 2：cluster 粒度**
  algorithm = 1,         # Louvain
  n.start   = 10,
  n.iter    = 10,
  verbose   = FALSE
)

cat("seurat_clusters 统计：\n")
print(table(obj$seurat_clusters))

# 5.3 UMAP（v3 风格参数，主要是 n.neighbors / min.dist 影响布局）
set.seed(SEED)
obj <- RunUMAP(
  obj,
  reduction    = "pca",
  dims         = dims_use,
  umap.method  = "uwot",
  metric       = "cosine",
  n.neighbors  = 30,     # **关键参数 3：UMAP 邻居数**
  min.dist     = 0.3,    # **关键参数 4：UMAP 簇的紧密程度**
  spread       = 1,
  init         = "spectral",
  n.components = 2,
  verbose      = FALSE
)

# 5.4 画 UMAP 聚类分布
png(file.path(res_dir, "stepF_UMAP_clusters_res1_GSE174367.png"),
    width = 2000, height = 1200)
print(
  DimPlot(obj, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.2) +
    NoLegend()
)
dev.off()

# 5.5 保存最终 StepF 对象
saveRDS(
  obj,
  file.path(res_dir, "stepF_afterPCA_graph_umap_res1_GSE174367.rds"),
  version = 2,
  compress = "xz"
)

cat("✅ StepF 完成：stepF_afterPCA_graph_umap_res1_GSE174367.rds 已保存\n")
log_env(res_dir, "stepF_done")








#验证
#--------------验证按照GSE157827PCA后得到聚类的正确性UMAP参数设置（dims=1:20、resolution=1、n.neighbors=30、min.dist=0.3）是否合适
#如何“专业地验证上面PCA后得到聚类的正确性”？方法 1：ARI（Adjusted Rand Index）评估 cluster 与作者 Cell.Type 的一致性
# 先安装mclust包（只需安装一次）
install.packages("mclust")

# 再加载包

library(mclust)
ref_labels <- as.integer(factor(obj$Cell.Type))
clu_labels <- as.integer(obj$seurat_clusters)
ARI <- adjustedRandIndex(ref_labels, clu_labels)
ARI


#✔ ✔ 方法 3：Marker DotPlot（生物学验证）

DefaultAssay(obj) <- "RNA"

markers_ref <- list(
  Astro   = c("AQP4","GFAP","ALDH1L1","SLC1A3","ADGRV1","GPC5","RYR3"), 
  Endo    = c("CLDN5","KDR","FLT1","PECAM1","ABCB1","EBF1"),
  Excit   = c("CAMK2A","SLC17A7","TBR1","CBLN2","LDB2"),
  Inhib   = c("GAD1","GAD2","SLC6A1","LHFPL3","PCDH15"),
  Microgl = c("C3","CX3CR1","P2RY12","AIF1","DOCK8","LRMDA"),
  Oligo   = c("MBP","MOG","PLP1","MOBP","ST18")
)

features_all <- unique(unlist(markers_ref))

DotPlot(obj, features = features_all,
        group.by = "seurat_clusters") +
  RotatedAxis()




#看按照GSE157827模版代码的UMAP参数设置（dims=1:20、resolution=1、n.neighbors=30、min.dist=0.3）进行调整对6个细胞类型UMAP的影响

#做一个小网格搜索：dims 和 resolution 组合里，谁的 ARI 更高？

#现在 PCA 已经算好了（50 PCs），我们只需要在同一对象上反复跑 Neighbors/Clusters，不用重复 PCA，很快的。

#下面这段代码帮你在 两种 dims × 三种 resolution 上算 ARI：
library(mclust)

res_dir  <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results"
obj_base <- readRDS(file.path(res_dir, "stepF_afterPCA_JackStraw_GSE174367.rds"))
DefaultAssay(obj_base) <- "integrated"

# 我们只关心 “前20个PC”和“前25个PC” 两种情况
param_grid <- expand.grid(
  max_dim    = c(20, 25),          # 20 vs 25
  resolution = c(0.8, 1.0, 1.2)    # 三种分辨率
)

results <- list()

for (i in seq_len(nrow(param_grid))) {
  max_dim <- param_grid$max_dim[i]
  res_use <- param_grid$resolution[i]
  dims_use <- 1:max_dim
  
  cat(">>> Testing dims 1:", max_dim, "resolution =", res_use, "\n")
  
  obj <- obj_base
  
  # Neighbors + clustering（不重复 PCA）
  set.seed(SEED)
  obj <- FindNeighbors(obj, dims = dims_use, k.param = 20, verbose = FALSE)
  set.seed(SEED)
  obj <- FindClusters(obj,
                      resolution = res_use,
                      algorithm  = 1,
                      n.start    = 10,
                      n.iter     = 10,
                      verbose    = FALSE)
  
  # ARI：和 Cell.Type 的一致性
  ref_labels <- as.integer(factor(obj$Cell.Type))
  clu_labels <- as.integer(obj$seurat_clusters)
  ARI <- adjustedRandIndex(ref_labels, clu_labels)
  
  # silhouette：在 PCA 空间评估簇分离度
  emb <- Embeddings(obj, "pca")[, dims_use]
  sil <- cluster::silhouette(as.integer(obj$seurat_clusters), dist(emb))
  sil_avg <- summary(sil)$avg.width
  
  results[[i]] <- data.frame(
    max_dim    = max_dim,
    resolution = res_use,
    ARI        = ARI,
    sil_avg    = sil_avg
  )
}

res_df <- do.call(rbind, results)
res_df

#验证结束，最终是只把resolution = 0.8    # 从1变为0.8，这一个变化
#最终结果建议：最终推荐参数（强烈建议使用）：
#✔ dims_use = 1:20  （保持不动）
#✔ resolution = 0.8    # 从1变为0.8
#✔ neighbors = 30 （保持不动）
#✔ min.dist = 0.3（保持不动）











#经过验证GSE157827的UMAP参数不太合适，找到真真的参数后，重新进行聚类+UMAP。JackStraw 前的PCA不用重新跑

###############################################################################
# GSE174367 - StepF 最终版
# 参数：dims_use = 1:20, resolution = 0.8, n.neighbors = 30, min.dist = 0.3
# 从 JackStraw 后的对象继续：stepF_afterPCA_JackStraw_GSE174367.rds
###############################################################################

SEED <- 20251023
set.seed(SEED)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(mclust)     # 计算 ARI
})

res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# 小工具：记录环境信息，方便以后复现
log_env <- function(res_dir, step_name = "step_env") {
  log_file <- file.path(res_dir, paste0(step_name, "_sessionInfo.txt"))
  SEED <- 20251023
  set.seed(SEED)
  sink(log_file)
  cat("SEED =", SEED, "\n")
  print(sessionInfo())
  sink()
  message("✅ 环境信息已保存到: ", log_file)
}

log_env(res_dir, "stepF_final_dims1_20_res0.8")


###############################################################################
# 1. 读取 JackStraw 后的对象（已经完成 ScaleData + PCA(50) + JackStraw）
###############################################################################

obj <- readRDS(file.path(res_dir, "stepF_afterPCA_JackStraw_GSE174367.rds"))
DefaultAssay(obj) <- "integrated"

# 我们使用前 20 个 PC
dims_use <- 1:20
cat("使用的 PCA 维度: 1-", max(dims_use), "\n")


###############################################################################
# 2. Neighbors + 聚类（resolution = 0.8）
###############################################################################

set.seed(SEED)
obj <- FindNeighbors(obj, dims = dims_use, k.param = 20, verbose = FALSE)

set.seed(SEED)
obj <- FindClusters(
  obj,
  resolution = 0.8,   # ★ 最终推荐的分辨率
  algorithm  = 1,     # Louvain
  n.start    = 10,
  n.iter     = 10,
  verbose    = FALSE
)

cat("seurat_clusters 统计（dims=1:20, res=0.8）：\n")
print(table(obj$seurat_clusters))


###############################################################################
# 3. UMAP（v3 风格参数）
###############################################################################

set.seed(SEED)
obj <- RunUMAP(
  obj,
  reduction    = "pca",
  dims         = dims_use,
  umap.method  = "uwot",
  metric       = "cosine",
  n.neighbors  = 30,   # 保持与 GSE157827 一致
  min.dist     = 0.3,  # 保持与 GSE157827 一致
  spread       = 1,
  init         = "spectral",
  n.components = 2,
  verbose      = FALSE
)

# 保存 UMAP 聚类图
png(file.path(res_dir, "stepF_UMAP_clusters_res0.8_dims1_20_GSE174367.png"),
    width = 2000, height = 1200)
print(
  DimPlot(obj, reduction = "umap",
          label = TRUE, label.size = 8, pt.size = 0.2) + NoLegend()
)
dev.off()
cat("✅ 已保存 UMAP 图: stepF_UMAP_clusters_res0.8_dims1_20_GSE174367.png\n")


###############################################################################
# 4. 计算 ARI（与作者 Cell.Type 的一致性）+ silhouette（几何分离度）
###############################################################################

# ARI
ref_labels <- as.integer(factor(obj$Cell.Type))
clu_labels <- as.integer(obj$seurat_clusters)
ARI <- adjustedRandIndex(ref_labels, clu_labels)

# silhouette（在 PCA 空间，使用 cluster::silhouette）
emb <- Embeddings(obj, "pca")[, dims_use]
sil <- cluster::silhouette(as.integer(obj$seurat_clusters), dist(emb))
sil_avg <- summary(sil)$avg.width

cat(sprintf("🔹 ARI (Cell.Type vs clusters) = %.3f\n", ARI))
cat(sprintf("🔹 平均 silhouette (dims=1:20) = %.3f\n", sil_avg))

# 把指标写入一个 txt，方便以后查
metrics_file <- file.path(res_dir, "stepF_metrics_dims1_20_res0.8_GSE174367.txt")
writeLines(
  c(
    paste0("dims_use = 1:20"),
    paste0("resolution = 0.8"),
    paste0("ARI = ", ARI),
    paste0("silhouette_avg = ", sil_avg)
  ),
  con = metrics_file
)
cat("✅ 已保存聚类评估指标到: ", metrics_file, "\n")


###############################################################################
# 5. 生物学验证：使用 GSE157827 的 markers_ref 画 DotPlot
###############################################################################

DefaultAssay(obj) <- "RNA"

markers_ref <- list(
  Astro   = c("AQP4","GFAP","ALDH1L1","SLC1A3","ADGRV1","GPC5","RYR3"), 
  Endo    = c("CLDN5","KDR","FLT1","PECAM1","ABCB1","EBF1"),
  Excit   = c("CAMK2A","SLC17A7","TBR1","CBLN2","LDB2"),
  Inhib   = c("GAD1","GAD2","SLC6A1","LHFPL3","PCDH15"),
  Microgl = c("C3","CX3CR1","P2RY12","AIF1","DOCK8","LRMDA"),
  Oligo   = c("MBP","MOG","PLP1","MOBP","ST18")
)

features_all <- unique(unlist(markers_ref))

png(file.path(res_dir, "stepF_DotPlot_markers_res0.8_dims1_20_GSE174367.png"),
    width = 2200, height = 1200)
print(
  DotPlot(obj, features = features_all,
          group.by = "seurat_clusters") + RotatedAxis()
)
dev.off()
cat("✅ 已保存 marker DotPlot: stepF_DotPlot_markers_res0.8_dims1_20_GSE174367.png\n")


###############################################################################
# 6. 保存最终 StepF 对象（用于之后的 UBL3 / AD vs Control 分析）
###############################################################################

saveRDS(
  obj,
  file.path(res_dir, "stepF_afterPCA_graph_umap_res0.8_dims1_20_GSE174367.rds"),
  version = 2,
  compress = "xz"
)

cat("✅ 最终对象已保存: stepF_afterPCA_graph_umap_res0.8_dims1_20_GSE174367.rds\n")
log_env(res_dir, "stepF_final_done")
###############################################################################














#这个代码是改成7个细胞类型（和原文一样）
###############################################################################
# 第 6 章：按簇找 marker + 利用 canonical markers 自动标注 6大细胞类型
# 数据集：GSE174367   对象起点：stepF_afterPCA_graph_umap_res0.8_dims1_20_GSE174367.rds
# 结果对象：stepH_obj_celltype6_named.rds
###############################################################################

SEED <- 20251023
set.seed(SEED)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(mclust)  # 之前用过，可留着
  library(plyr)
  library(ggplot2)
})

res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# 环境记录函数（保证可复现性）
log_env <- function(res_dir, step_name = "step_env") {
  log_file <- file.path(res_dir, paste0(step_name, "_sessionInfo.txt"))
  SEED <- 20251023
  set.seed(SEED)
  sink(log_file)
  cat("SEED =", SEED, "\n")
  print(sessionInfo())
  sink()
  message("✅ 环境信息已保存到: ", log_file)
}
log_env(res_dir, "stepG_markers_and_celltype6")



## 6.1 读取上一步的整合对象 -----------------------------------------------

obj <- readRDS(file.path(
  res_dir,
  "stepF_afterPCA_graph_umap_res0.8_dims1_20_GSE174367.rds"
))

# 用 integrated assay（和 GSE157827 一样）
DefaultAssay(obj) <- "integrated"

# 当前身份设为 seurat_clusters
Idents(obj) <- "seurat_clusters"

cat("当前簇数量：", nlevels(Idents(obj)), "\n")
print(table(Idents(obj)))


###############################################################################
# 6.2 每簇找 marker（与 GSE157827/原文一致：Wilcoxon + logFC≥0.25 + adj.P<0.1）
###############################################################################

## 关掉 future 并行，防止 future.globals.maxSize 报错
plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)  # 预防以后如果又启用并行

markers_all <- FindAllMarkers(
  obj,                     # Seurat 对象
  only.pos        = FALSE, # 同时找上调和下调基因
  test.use        = "wilcox",
  logfc.threshold = 0.25,  # 最小 log2(倍数变化)
  min.pct         = 0.1,   # 至少在10%的细胞里表达
  verbose         = TRUE
)

cat("FindAllMarkers 返回行数：", nrow(markers_all), "\n")

if (nrow(markers_all) == 0) {
  stop("❌ FindAllMarkers 没有找到任何 DE 基因，多半是之前并行/future 的问题，请检查。")
}

# p 值校正后是否 <0.1（和 GSE157827 模板一致，也符合原文）
markers_all$pass_0.1 <- markers_all$p_val_adj < 0.1

# 兼容不同 Seurat 版本：有的叫 avg_log2FC，有的叫 avg_logFC
cn <- colnames(markers_all)
if (!"avg_log2FC" %in% cn && "avg_logFC" %in% cn) {
  markers_all$avg_log2FC <- markers_all$avg_logFC
}
if (!"avg_log2FC" %in% colnames(markers_all)) {
  stop("❌ markers_all 中既没有 avg_log2FC 也没有 avg_logFC，请用 colnames(markers_all) 检查。")
}

# 保存完整结果（含上下调基因），CSV 文件可在 Excel 打开
write.csv(
  markers_all,
  file.path(res_dir, "stepG_FindAllMarkers_wilcox_logfc0.25_GSE174367.csv"),
  row.names = FALSE
)
cat("✅ 已保存所有 marker：stepG_FindAllMarkers_wilcox_logfc0.25_GSE174367.csv\n")


###############################################################################
# 6.3 每簇取前 20 个正向 marker，画 DotPlot / Heatmap 辅助命名
# —— 写法和 GSE157827 完全一致，只是文件名加上 GSE174367
###############################################################################

# 只保留上调的基因（avg_log2FC > 0 且通过显著性）
top20 <- subset(markers_all, avg_log2FC > 0 & pass_0.1)

# 按簇号和 logFC 从大到小排序
top20 <- top20[order(top20$cluster, -top20$avg_log2FC), ]

# 每个簇取前 20 个上调基因（视觉展示用）
top20 <- do.call(rbind, by(top20, top20$cluster, head, n = 20))

feat20 <- unique(top20$gene)
cat("用于 DotPlot/Heatmap 的 marker 基因数：", length(feat20), "\n")

# DotPlot（每个簇前 20 个 marker）
png(file.path(res_dir, "stepG_DotPlot_top20_per_cluster_GSE174367.png"),
    width = 2200, height = 1400, res = 180)
print(DotPlot(obj, features = feat20) + RotatedAxis())
dev.off()

# Heatmap（同样这些 marker）
png(file.path(res_dir, "stepG_Heatmap_top20_per_cluster_GSE174367.png"),
    width = 2200, height = 1400, res = 180)
print(DoHeatmap(obj, features = feat20, raster = TRUE))
dev.off()

cat("✅ 已输出每簇 top20 marker 的 DotPlot 和 Heatmap（GSE174367 版本）\n")
###############################################################################

###############################################################################


# 6.4 基于 canonical markers 计算“簇 × 细胞大类”的平均分数
# 注意：GSE174367 的行名已经是基因 symbol，不需要再做 ENSG→SYMBOL 映射
###############################################################################

# 使用 RNA assay 的 log-normalized 数据
DefaultAssay(obj) <- "RNA"
mat_sym <- GetAssayData(obj, slot = "data")  # 行=SYMBOL, 列=cell
sym_all <- rownames(mat_sym)

clu <- Idents(obj)
clu_levels <- levels(clu)

# GSE157827 使用的 6+1 类 marker 列表（7类：ODC、INH、EX、OPC、ASC、MG、PER/END）
markers_ref <- list(
  "Astrocytes",
  "Excitatory neurons",
  "Microglia",
  "Endothelial",
  "Inhibitory neurons",
  "Oligodendrocytes"
)


score_one_group <- function(genes) {
  g <- intersect(genes, sym_all)
  if (length(g) == 0) {
    return(setNames(rep(NA_real_, length(clu_levels)), clu_levels))
  }
  per_cell <- Matrix::colMeans(mat_sym[g, , drop = FALSE])
  tapply(per_cell, INDEX = clu, FUN = mean, na.rm = TRUE)[clu_levels]
}

avg_by_cluster <- sapply(markers_ref, score_one_group)
stopifnot(identical(rownames(avg_by_cluster), clu_levels))

# 若存在 NA，改为 -Inf，便于 max.col 选择最大类
tmp <- avg_by_cluster
tmp[is.na(tmp)] <- -Inf

# 7 类标签（含 Peri），用于审计
lab_6 <- colnames(tmp)[max.col(tmp, ties.method = "first")]
names(lab_6) <- rownames(tmp)



# 保存簇得分 + 标签表，便于检查/以后回看
# 保存簇得分 + 7 类短标签
df_scores <- data.frame(
  cluster = rownames(avg_by_cluster),
  label_6 = lab_6,
  avg_by_cluster,
  check.names = FALSE
)

write.csv(
  df_scores,
  file.path(res_dir, "stepH_cluster_scores_label6_GSE174367.csv"),
  row.names = FALSE
)
cat("✅ 已保存 cluster×7大类 分数表：stepH_cluster_scores_label7E174367.csv\n")


###############################################################################

## 6.5 把簇标签写回 Seurat 对象：celltype7（短标签）
obj$celltype6_short <- plyr::mapvalues(
  Idents(obj),
  from = names(lab_6),
  to   = unname(lab_6)
)

## 6.6 变成人类可读的长名字
nice_levels  <- c("Astrocytes",
                  "Excitatory neurons",
                  "Microglia",
                  "Inhibitory neurons",
                  "Oligodendrocytes",
                  "Endothelial")

short2nice   <- c(
  Astro   = "Astrocytes",
  Excit   = "Excitatory neurons",
  Microgl = "Microglia",
  Endo    = "Endothelial",
  Inhib   = "Inhibitory neurons"，
  Oligo   = "Oligodendrocytes",
)

val <- as.character(obj$celltype7_short)
val[val %in% names(short2nice)] <- short2nice[val[val %in% names(short2nice)]]
obj$celltype7 <- factor(val, levels = nice_levels)

cat("7 大类分布：\n")
print(table(obj$celltype7))


# 若没有 sample 列，就用 SampleID 补一个
if (!"sample" %in% colnames(obj@meta.data) && "SampleID" %in% colnames(obj@meta.data)) {
  obj$sample <- obj$SampleID
}

# 6 类的细胞数 & 占比
counts6<- table(obj$celltype7)
props6  <- round(100 * prop.table(counts7), 1)

cat("6 大类 cell 数:\n")
print(counts6)
cat("6 大类 cell 占比(%):\n")
print(props6)

write.csv(
  data.frame(celltype = names(counts6),
             n_cells = as.integer(counts6),
             percent = as.numeric(props6)),
  file.path(res_dir, "stepH_celltype7_counts_percent_GSE174367.csv"),
  row.names = FALSE
)

# 画 Cell type UMAP（长名字版本）
pal7 <- c(
  "Astrocytes"             = "#FF8E8E",  # 粉红（原样）
  "Excitatory neurons"     = "#09BB3C",  # 亮绿（原样）
  "Microglia"              = "#36B0E1",  # 蓝（原样）
  "Endothelial"            = "#B8A109",  # 橄榄黄（原样）
  "Inhibitory neurons"     = "#00BFC4",  # 青蓝（原样）
  "Oligodendrocytes"       = "#E16AFC",  # 粉紫（原样）
)

png(file.path(res_dir, "stepH_umap_celltype7_named_GSE174367.png"),
    width = 2000, height = 1400, res = 180)
print(
  DimPlot(obj, reduction = "umap", group.by = "celltype7",
          label = TRUE, label.size = 5, repel = TRUE) +
    scale_color_manual(values = pal7, drop = FALSE) +
    labs(color = "Cell type")
)
dev.off()
cat("✅ 已保存 7 大类命名版 UMAP：stepH_umap_celltype7_named_GSE174367.png\n")

# 若存在 sample 信息，再画每个样本的 7 类比例
if ("sample" %in% colnames(obj@meta.data)) {
  df_bar <- as.data.frame(table(sample = obj$sample, celltype = obj$celltype7))
  df_bar <- within(df_bar, {
    total_by_sample <- ave(Freq, sample, FUN = sum)
    percent <- 100 * Freq / total_by_sample
  })
  png(file.path(res_dir, "stepH_bar_sample_celltype6_GSE174367.png"),
      width = 1800, height = 1200, res = 180)
  print(
    ggplot(df_bar, aes(x = sample, y = percent, fill = celltype)) +
      geom_bar(stat = "identity", width = 0.9) +
      scale_fill_manual(values = pal7, drop = FALSE) +
      coord_flip() +
      labs(x = NULL, y = "Proportion (%)", fill = "Cell type") +
      theme_bw(base_size = 12)
  )
  dev.off()
  cat("✅ 已保存样本 × 6比例条形图：stepH_bar_sample_celltype6_GSE174367.png\n")
}


###############################################################################
# 6.7 保存最终对象：stepH_obj_celltype6_named.rds
###############################################################################

saveRDS(
  obj,
  file.path(res_dir, "stepH_obj_celltype6_named.rds"),
  version = 2,
  compress = "xz"
)


cat("🎉 最终对象已保存: stepH_obj_celltype6_named.rds\n")
cat("   之后的 Cell type UMAP + UBL3 UMAP + pseudo-bulk DESeq2 分析均以此对象为起点。\n")

log_env(res_dir, "stepH_celltype7_done")
###############################################################################









###############################################################################
# 第 7 章（2）UMAP：6 个细胞类型 + UBL3（>0 高亮）
# 数据集：GSE174367 版本（按 GSE157827 模版改写）
# 起点对象：
#   - stepH_obj_celltype6_named.rds（带 UMAP + celltype6）
#   - stepC_filtered_obj.rds（QC 后 RNA counts）
###############################################################################

## ============================================================
## 0. 基本设置：随机数种子 + 加载 R 包
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
})

## ============================================================
## 0.1 全局统一的 celltype6 名称 + 颜色（模版核心）
## ============================================================

## 6 个“标准名字”，以后所有 snRNA 数据都朝这个名字对齐
celltype6_levels_std <- c(
  "Astrocytes",
  "Excitatory neurons",
  "Microglia",
  "Endothelial",
  "Inhibitory neurons",
  "Oligodendrocytes"
)

## 和 GSE157827 / 174367 一样的统一配色
celltype6_palette_std <- c(
  "Astrocytes"          = "#FF8E8E",  # 粉红
  "Excitatory neurons"  = "#09BB3C",  # 亮绿
  "Microglia"           = "#36B0E1",  # 蓝
  "Endothelial"         = "#B8A109",  # 橄榄黄
  "Inhibitory neurons"  = "#00BFC4",  # 青蓝
  "Oligodendrocytes"    = "#E16AFC"   # 粉紫
)



## 把不同写法统一成“标准名字”的函数
## 如果 GSE174367 里有别的写法，可以在下面补充映射
standardize_celltype6 <- function(x) {
  x <- as.character(x)
  x[x %in% c("Astrocytes","Astro","Astrocyte")] <- "Astrocytes"
  x[x %in% c("Excitatory neurons","Excit")]      <- "Excitatory neurons"
  x[x %in% c("Inhibitory neurons","Inhib")]      <- "Inhibitory neurons"
  x[x %in% c("Microglia","Microgl")]            <- "Microglia"
  x[x %in% c("Oligodendrocytes","Oligo")]       <- "Oligodendrocytes"
  x[x %in% c("OPCs","OPC")]                     <- "OPCs"
  x[x %in% c("Pericytes/Endothelial","PerEnd","PER/END")] <- "Pericytes/Endothelial"
  factor(x, levels = celltype6_levels_std)
}


## ============================================================
## 1. 路径与对象：读入 UMAP 对象 & QC 后 counts 对象
## ============================================================

## ⚠ 这里改成 GSE174367 的 results 目录
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results"

## 1.1 stepH 对象：已经整合 + 聚类 + celltype6 + UMAP
obj <- readRDS(file.path(res_dir, "stepH_obj_celltype6_named.rds"))

## 1.2 stepC 对象：QC 后的 RNA counts（单层），不含 UMAP
obj_expr <- readRDS(file.path(res_dir, "stepC_filtered_obj.rds"))
DefaultAssay(obj_expr) <- "RNA"

## ============================================================
## 2. 从 QC 后 counts 计算 UBL3 的 log1p(CPM)，构建 df_all
## ============================================================

## 2.1 counts 矩阵：行 = 基因 (ENSG/SYMBOL)，列 = 细胞
rna_counts <- tryCatch(
  GetAssayData(obj_expr, assay = "RNA", slot = "counts"),
  error = function(e) {
    message("GetAssayData(slot='counts') 出错，改用 LayerData(slot='counts')")
    LayerData(obj_expr[["RNA"]], layer = "counts")
  }
)
cat("counts 矩阵维度：", nrow(rna_counts), "基因 ×", ncol(rna_counts), "细胞\n")

## 2.2 找到 UBL3 的 ENSEMBL ID
map_ubl3 <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = "UBL3",
  keytype  = "SYMBOL",
  columns  = "ENSEMBL"
)
ubl3_id <- map_ubl3$ENSEMBL[1]
cat("UBL3 ENSG ID:", ubl3_id, "\n")

## 2.3 在 counts 中定位 UBL3 一行（优先 ENSG，不行用 SYMBOL）
if (!is.null(ubl3_id) && ubl3_id %in% rownames(rna_counts)) {
  gene_row <- ubl3_id
} else if ("UBL3" %in% rownames(rna_counts)) {
  gene_row <- "UBL3"
} else {
  stop("在 counts 中找不到 UBL3（既无 ENSEMBL 也无 SYMBOL）")
}
cat("实际使用的 UBL3 行名为：", gene_row, "\n")

## 2.4 对齐细胞名：只保留同时出现在 obj_expr 和 obj 里的细胞
cells_expr   <- colnames(obj_expr)
cells_umap   <- colnames(obj)
common_cells <- intersect(cells_expr, cells_umap)
cat("两个对象共有细胞数：", length(common_cells), "\n")

## 2.5 提取 UBL3 counts & library size
rna_sub  <- rna_counts[gene_row, common_cells, drop = FALSE]  # 1 × N
raw_vec  <- as.numeric(rna_sub[1, ])                          # UBL3 原始 counts
lib_size <- Matrix::colSums(rna_counts[, common_cells, drop = FALSE])

## 2.6 按 Seurat::NormalizeData(LogNormalize) 计算 log1p(CPM)
UBL3_norm <- log1p((raw_vec / lib_size) * 1e4)

## 2.7 从 stepH 对象拿 meta 信息（sample / celltype6 / group）
meta_sub <- obj@meta.data[common_cells, c("sample", "celltype6", "group")]

## 2.8 构建 df_all：直方图 & KS 检验用的整合数据框
df_all <- data.frame(
  cell      = common_cells,
  UBL3      = UBL3_norm,
  raw       = raw_vec,
  sample    = meta_sub$sample,
  celltype6 = meta_sub$celltype6,
  group     = meta_sub$group,
  stringsAsFactors = FALSE
)

df_all <- df_all[!is.na(df_all$sample) &
                   !is.na(df_all$celltype6) &
                   !is.na(df_all$group), ]

cat("df_all 总细胞数：", nrow(df_all), "\n")
cat("细胞类型（原始命名）：", paste(sort(unique(df_all$celltype6)), collapse = ", "), "\n")

## ============================================================
## 3. 数据一致性检查：将 UBL3 加入 UMAP 对象
## ============================================================

## 3.1 按 UMAP 对象的细胞顺序排列 UBL3
ubl3_for_umap <- df_all$UBL3[match(colnames(obj), df_all$cell)]

## 3.2 加入 Seurat 对象 meta.data
obj$UBL3_log1p <- ubl3_for_umap

cat("UBL3_log1p 在 obj 中的范围：",
    min(obj$UBL3_log1p, na.rm = TRUE), " ~ ",
    max(obj$UBL3_log1p, na.rm = TRUE), "\n")

## ============================================================
## 4. 统一 UMAP 视觉主题 + celltype6 标准化（模版通用部分）
## ============================================================

# 统一绘图主题
umap_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    legend.title = element_text(face = "bold"),
    panel.border = element_blank(),
    panel.grid   = element_blank()
  )

# ★★ 关键一步：把当前数据集的 celltype6 统一成“标准名字” + 固定顺序
obj$celltype6 <- standardize_celltype6(obj$celltype6)

cat("标准化后的 celltype6：\n")
print(table(obj$celltype6, useNA = "ifany"))
# 如果看到有 <NA>，就说明有没映射的名称，需要回到 standardize_celltype6() 里补充。

###############################################################################
# 4. Panel A：Cell type UMAP（使用统一模版颜色）
###############################################################################

DefaultAssay(obj) <- "integrated"

p_celltype <- DimPlot(
  obj,
  reduction  = "umap",
  group.by   = "celltype6",
  label      = TRUE,
  label.size = 4.5,
  repel      = TRUE,
  raster     = FALSE,   # 强制用矢量点
  pt.size    = 0.1      # 可按需要微调
) +
  scale_color_manual(
    values = celltype6_palette_std,
    breaks = celltype6_levels_std,
    limits = celltype6_levels_std,
    drop   = TRUE
  ) +
  ggtitle("Cell Type UMAP") +
  labs(color = "Cell type") +
  umap_theme +
  theme(
    legend.position = "bottom",
    legend.box      = "horizontal"
  ) +
  guides(
    color = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(size = 5)  # 图例里的点大小
    )
  )

###############################################################################
# 5. Panel B：UBL3 expression UMAP（>0 高亮）
###############################################################################

DefaultAssay(obj) <- "RNA"

umap_coords <- Embeddings(obj, "umap")
df_umap <- data.frame(
  umap_1 = umap_coords[, 1],
  umap_2 = umap_coords[, 2],
  UBL3   = obj$UBL3_log1p
)

df_bg  <- subset(df_umap, UBL3 <= 0)
df_pos <- subset(df_umap, UBL3 >  0)

break_vals   <- c(0, 1, 2, 3, 4)
break_labels <- c("0", "1", "2", "3", "4+")

p_ubl3_highlight <- ggplot() +
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
    limits  = c(0, max(df_all$UBL3)),
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

###############################################################################
# 6. 拼图并导出：Panel A + Panel B
###############################################################################

p_AB_high <- p_celltype + p_ubl3_highlight +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A")

out_file <- file.path(res_dir, "Fig_UBL3_Celltype_and_UMAP_highlight_GSE174367.png")
ggsave(out_file, p_AB_high, width = 12, height = 6, dpi = 300)
cat("🎉 已输出两联图（PNG）：", out_file, "\n")

out_pdf <- file.path(res_dir, "Fig_UBL3_Celltype_and_UMAP_highlight_GSE174367.pdf")
ggsave(out_pdf, p_AB_high, width = 12, height = 6)
cat("🎉 已输出两联图（PDF）：", out_pdf, "\n")













## ============================================================
## 第 2 章（0）：6 张【每个细胞类型 × 每个样本】直方图（y=Cell count）
## ===2 章（0）：UBL3>0 每个细胞类型的直方图（按样本分面）
## ===========================
SEED <- 20251023; set.seed(SEED)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

## 按需要改成 GSE174367 的结果路径
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results"

## 1. 读入对象 -------------------------------------------------
obj_meta <- readRDS(file.path(res_dir, "stepH_obj_celltype6_named.rds"))
obj_expr <- readRDS(file.path(res_dir, "stepC_filtered_obj.rds"))
DefaultAssay(obj_expr) <- "RNA"

## 2. 从 counts 里算 UBL3 的 log1p(CPM) -----------------------
rna_counts <- tryCatch(
  GetAssayData(obj_expr, assay = "RNA", slot = "counts"),
  error = function(e) LayerData(obj_expr[["RNA"]], layer = "counts")
)

# 找 UBL3 的 ENSEMBL ID
ubl3_id <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = "UBL3",
  keytype = "SYMBOL",
  columns = "ENSEMBL"
)$ENSEMBL[1]

# 兼容行名用 SYMBOL 的情况
gene_row <- if (ubl3_id %in% rownames(rna_counts)) {
  ubl3_id
} else if ("UBL3" %in% rownames(rna_counts)) {
  "UBL3"
} else {
  stop("在 counts 中找不到 UBL3（无 ENSEMBL 也无 SYMBOL）")
}

# 对齐细胞
common_cells <- intersect(colnames(obj_expr), colnames(obj_meta))
if (length(common_cells) == 0) stop("两个对象没有共同细胞名")

raw_vec  <- as.numeric(rna_counts[gene_row, common_cells])
lib_size <- Matrix::colSums(rna_counts[, common_cells, drop = FALSE])
UBL3_norm <- log1p((raw_vec / lib_size) * 1e4)  # log1p(CPM)

meta_sub <- obj_meta@meta.data[common_cells, c("sample", "celltype6", "group")]

df_all <- data.frame(
  cell      = common_cells,
  UBL3      = UBL3_norm,
  raw       = raw_vec,
  sample    = meta_sub$sample,
  celltype6 = meta_sub$celltype6,
  group     = meta_sub$group,
  stringsAsFactors = FALSE
)

## 只保留 UBL3 > 0 的细胞
df_pos <- subset(df_all, UBL3 > 0)

## 3. 画图 -----------------------------------------------------
# 颜色：AD 橙、Control 蓝
group_colors <- c(AD = "#D55E00", Control = "#0072B2")

# facet 用的标签：比如 AD1_AD / NC14_Control
df_pos$sample_group <- with(df_pos, paste0(sample, "_", group))
df_pos$sample_group <- factor(df_pos$sample_group,
                              levels = sort(unique(df_pos$sample_group)))

celltypes <- sort(unique(df_pos$celltype6))

# 输出目录
out_dir <- file.path(res_dir, "Fig_UBL3_hist_per_celltype_noZero_GSE174367")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sanitize_name <- function(x) {
  x <- gsub(" ", "_", x)
  x <- gsub("/", "_", x)
  x
}

for (ct in celltypes) {
  df_ct <- df_pos[df_pos$celltype6 == ct, ]
  if (nrow(df_ct) == 0) next
  
  p_ct <- ggplot(df_ct, aes(x = UBL3, fill = group)) +
    geom_histogram(
      bins     = 40,
      position = "identity",
      alpha    = 0.7,
      color    = "grey30"
    ) +
    facet_wrap(~ sample_group, ncol = 5) +
    scale_fill_manual(values = group_colors, name = "Group") +
    labs(
      title = paste0("UBL3>0 distribution in ", ct),
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
    filename = file.path(out_dir,
                         paste0("UBL3_noZero_", sanitize_name(ct), "_GSE174367.png")),
    plot   = p_ct,
    width  = 12,
    height = 8,
    dpi    = 300
  )
}

cat("✅ 已为 GSE174367 生成 UBL3>0 直方图（每个细胞类型一张）：", out_dir, "\n")








#第 2 章（1）：1 张 6 个细胞类型的 overlap 直方图（y=Cell count，右上角标检验）
## ============================================================
## 第 2章（1）：overlap 直方图 + 固定 Mann–Whitney U 检验
## y 轴 = Cell count，右上角显示检验方法和 FDR
## 适用于 GSE174367
## ============================================================
###############################################################################
## NO1_GSE174367_OverlapHistCount_AD_vs_Control_cell_and_donor_level.R
##
## 【第1章目标】（严格模仿 GSE157827 模版）
##  - 使用 6 个细胞类型（celltype6）做 overlap 直方图
##  - Y轴 = Cell count（不是 density）
##  - 只使用“表达细胞”（expr > 0）
##  - 同时输出两张图：
##      (1) cell-level：每个 cell 作为观测 → Mann–Whitney U → BH(FDR)
##      (2) donor-level：每个 donor 先汇总(默认 median) → MWU → BH(FDR)
##
## 【UBL3表达计算】严格按模版：
##    expr = log1p( (raw_counts / lib_size) * 1e4 )   # CP10K
##    X轴显示：UBL3 log1p(CP10k)
##
## 【输入对象】强制只用：
##    D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/stepH_obj_celltype6_named.rds
##
## 【输出目录】全部保存在：
##    D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/NO1_overlap_hist_count_GSE157827style
##
## 【补充材料/中间结果】自动保存（若已存在则跳过）：
##  - 每组 donor 数：QC_donors_by_group4.csv
##  - 表达细胞 df0：INTERMEDIATE_df0_exprGT0.csv/.rds
##  - 统计输入表：INTERMEDIATE_stat_input_*.csv/.rds
##  - 每 celltype 的 p/padj 表：STATS_*.csv
##  - 日志与 sessionInfo：NO1_log.txt
###############################################################################

## =========================
## 0) 清理环境 + 基本设置
## =========================
rm(list = ls()); gc()

## 强制英文环境，避免某些系统输出中文导致日志/解析不一致
Sys.setenv(LANG = "en")

## 固定随机种子，保证可复现（即使本章本身随机性不强，也建议保留）
SEED <- 20251023
set.seed(SEED)

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
## 1) 路径参数（★强制只用 celltype6 对象★）
## =========================
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results"

## ★强制只用这个对象：不存在就立刻报错，不允许自动 fallback 到 celltype7
obj_fp <- file.path(res_dir, "stepH_obj_celltype6_named.rds")
stopifnot(file.exists(obj_fp))

## 输出目录（本章统一放这里）
out_dir <- file.path(res_dir, "NO1_overlap_hist_count_GSE157827style")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## 数据集/基因/绘图参数
dataset_tag <- "GSE174367"
gene_symbol <- "UBL3"
binwidth    <- 0.2       # 直方图 bin 宽度（与模版一致）
disease     <- "AD"      # 本章固定 AD vs Control

## 日志文件：把所有关键检查、统计、sessionInfo 都写进去，方便投稿补充材料
log_fp <- file.path(out_dir, "NO1_log.txt")
sink(log_fp)
cat("==== START ====\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("SEED:", SEED, "\n")
cat("res_dir:", res_dir, "\n")
cat("obj_fp:", obj_fp, "\n")
cat("out_dir:", out_dir, "\n\n")
sink()

## =========================
## 2) 读入对象 + 强制检查 celltype6 必须存在且只有 6 类
## =========================
obj <- readRDS(obj_fp)
DefaultAssay(obj) <- "RNA"

md <- obj@meta.data

## 2.1 强制 celltype6 列存在
if (!("celltype6" %in% colnames(md))) {
  sink(log_fp, append = TRUE)
  cat("❌ meta.data 列名如下（找不到 celltype6）：\n")
  print(colnames(md))
  sink()
  stop("❌ 该对象 meta.data 中缺少 celltype6 列。请回到上游命名步骤修正。")
}

## 2.2 强制 celltype6 正好 6 类（否则 facet 会不等于6）
ct_levels <- sort(unique(as.character(md$celltype6)))
n_ct <- length(ct_levels)

sink(log_fp, append = TRUE)
cat("celltype6 unique =", n_ct, "\n")
print(ct_levels)
cat("\n")
sink()

if (n_ct != 6) {
  stop("❌ celltype6 不是 6 类（请检查上游 celltype6 的命名/过滤是否正确）。")
}

## 2.3 自动识别 group 列（优先 group4）
grp_col_candidates <- c("group4","group2","group","Group","diagnosis","Dx","clinical_diagnosis")
grp_col <- grp_col_candidates[grp_col_candidates %in% colnames(md)][1]
if (is.na(grp_col)) {
  sink(log_fp, append = TRUE)
  cat("❌ meta.data 列名如下（找不到 group 列候选）：\n")
  print(colnames(md))
  sink()
  stop("❌ meta.data 找不到 group 列（group4/group/...）。")
}

## 2.4 自动识别 donor 列（优先 autopsy_id，其次 donor/sample）
donor_col_candidates <- c("autopsy_id","donor","Donor","sample","Sample","orig.ident","patient","subject")
don_col <- donor_col_candidates[donor_col_candidates %in% colnames(md)][1]
if (is.na(don_col)) {
  sink(log_fp, append = TRUE)
  cat("❌ meta.data 列名如下（找不到 donor 列候选）：\n")
  print(colnames(md))
  sink()
  stop("❌ meta.data 找不到 donor 列（autopsy_id/sample/...）。")
}

## 2.5 标准化 group：把各种 Control/NC/CTRL 等统一为 "Control"
grp_raw <- trimws(as.character(md[[grp_col]]))
ctrl_alias <- c("Control","CTRL","Ctr","CTR","NC","Normal","N","control","ctrl","ctr","nc","normal")
grp_std <- ifelse(grp_raw %in% ctrl_alias, "Control", grp_raw)

## 2.6 写入统一字段（后续所有分析都用这三个标准列名）
md$celltype6_std <- as.character(md$celltype6)                 # 明确：来自 celltype6
md$group4_std    <- grp_std                                    # 明确：统一 Control 命名
md$donor_std     <- trimws(as.character(md[[don_col]]))         # 明确：donor id

## 2.7 检查是否同时存在 AD 和 Control
if (!all(c(disease, "Control") %in% unique(md$group4_std))) {
  sink(log_fp, append = TRUE)
  cat("⚠ group unique after standardization:\n")
  print(sort(unique(md$group4_std)))
  sink()
  stop("❌ 分组里没有同时包含 AD 和 Control（或 Control 同义词未被识别）。")
}

## 2.8 统计 donor 数（以 donor×group 唯一对为准，作为 legend 的 n 口径）
don_all <- unique(md[, c("donor_std","group4_std")])
qc_don <- as.data.frame(table(don_all$group4_std), stringsAsFactors = FALSE)
colnames(qc_don) <- c("group4","n_donors")

qc_don_fp <- file.path(out_dir, "QC_donors_by_group4.csv")
if (!file.exists(qc_don_fp)) write.csv(qc_don, qc_don_fp, row.names = FALSE)

sink(log_fp, append = TRUE)
cat("Loaded cells:", ncol(obj), " genes:", nrow(obj), "\n")
cat("group col:", grp_col, " | donor col:", don_col, "\n")
cat("Donors by group4_std:\n"); print(table(don_all$group4_std))
cat("\n")
sink()

## =========================
## 3) 获取 counts（兼容 Seurat v5 多 layers）
## =========================
## 说明：
##  - Seurat v5 可能把 counts 分散在多个 counts.* layers
##  - 我们需要把所有 counts 合并成一个 gene×cell 矩阵，并严格按 obj 的 cell 顺序对齐
get_counts_matrix_allcells <- function(obj, assay = "RNA") {
  a <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(a), error = function(e) character(0))
  counts_layers <- layers[grepl("^counts", layers)]
  
  ## 情况A：存在 counts layers
  if (length(counts_layers) > 0) {
    sink(log_fp, append = TRUE)
    cat("Detected counts layers:\n"); print(counts_layers)
    sink()
    
    mats <- list()
    for (ly in counts_layers) {
      m <- tryCatch(SeuratObject::LayerData(a, layer = ly), error = function(e) NULL)
      if (is.null(m)) next
      if (!is.null(dim(m)) && length(dim(m)) == 2) mats[[ly]] <- m
    }
    if (length(mats) == 0) stop("❌ counts layers 存在，但读取 LayerData 失败。")
    
    ## 基因对齐：以第一个 layer 的 genes 作为参考
    ref_genes <- rownames(mats[[1]])
    for (k in names(mats)) {
      if (!identical(rownames(mats[[k]]), ref_genes)) {
        m0 <- mats[[k]]
        m_aligned <- Matrix::Matrix(0, nrow = length(ref_genes), ncol = ncol(m0), sparse = TRUE)
        rownames(m_aligned) <- ref_genes
        colnames(m_aligned) <- colnames(m0)
        common <- intersect(ref_genes, rownames(m0))
        m_aligned[common, ] <- m0[common, , drop = FALSE]
        mats[[k]] <- m_aligned
      }
    }
    
    ## 合并列（cells）
    mat_all <- if (length(mats) == 1) mats[[1]] else Reduce(Matrix::cbind2, mats)
    
    ## 去除重复 cell（防止多层重复）
    dup <- duplicated(colnames(mat_all))
    if (any(dup)) {
      sink(log_fp, append = TRUE)
      cat("⚠ duplicated cells across layers:", sum(dup), " -> keep first\n")
      sink()
      mat_all <- mat_all[, !dup, drop = FALSE]
    }
    
    ## 补齐缺失 cell（理论上不应发生；发生则补零保证对齐）
    all_cells <- colnames(obj)
    miss_cells <- setdiff(all_cells, colnames(mat_all))
    if (length(miss_cells) > 0) {
      sink(log_fp, append = TRUE)
      cat("⚠ missing cells in merged counts:", length(miss_cells), " -> fill zeros\n")
      sink()
      m_fill <- Matrix::Matrix(0, nrow = nrow(mat_all), ncol = length(miss_cells), sparse = TRUE)
      rownames(m_fill) <- rownames(mat_all)
      colnames(m_fill) <- miss_cells
      mat_all <- Matrix::cbind2(mat_all, m_fill)
    }
    
    ## 强制列顺序与 obj 一致（非常关键，否则 meta 对不上表达）
    mat_all <- mat_all[, colnames(obj), drop = FALSE]
    return(mat_all)
  }
  
  ## 情况B：没有 layers，尝试传统 counts slot
  m2 <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "counts"), error = function(e) NULL)
  if (!is.null(m2)) return(m2)
  
  stop("❌ 无法获取 counts（既无 counts layers，也无法 GetAssayData counts）。")
}

rna_counts <- get_counts_matrix_allcells(obj, "RNA")

## 强制检查：counts 的 cell 数必须与对象一致
stopifnot(ncol(rna_counts) == ncol(obj))
stopifnot(identical(colnames(rna_counts), colnames(obj)))

sink(log_fp, append = TRUE)
cat("Merged counts dim:", paste(dim(rna_counts), collapse = " x "), "\n\n")
sink()

## =========================
## 4) 计算 UBL3 log1p(CP10k) 并构建表达细胞 df0(expr>0)
## =========================
## 4.1 找到 UBL3 对应行名：
##     - 优先 UBL3 gene symbol
##     - 若是 Ensembl 行名，则尝试 ENSG00000122042
gene_row <- if ("UBL3" %in% rownames(rna_counts)) {
  "UBL3"
} else if ("ENSG00000122042" %in% rownames(rna_counts)) {
  "ENSG00000122042"
} else {
  NA_character_
}
if (is.na(gene_row)) stop("❌ counts 行名中找不到 UBL3（UBL3 或 ENSG00000122042）。")

## 4.2 计算每个 cell 的 library size（总 counts）
lib_size <- Matrix::colSums(rna_counts)

## 4.3 取出 UBL3 raw counts（每个 cell 一个数）
raw_vec <- as.numeric(rna_counts[gene_row, , drop = TRUE])

## 4.4 严格按模版计算 CP10K 并 log1p
expr <- log1p((raw_vec / pmax(lib_size, 1)) * 1e4)

## 4.5 合并 meta 信息（注意：必须和 counts 的列顺序一致，我们已强制检查过）
df_all <- data.frame(
  expr      = expr,
  donor     = md$donor_std,
  group4    = md$group4_std,
  celltype6 = md$celltype6_std,
  stringsAsFactors = FALSE
)

## 4.6 只保留“表达细胞”（expr > 0）
df0 <- df_all[df_all$expr > 0, , drop = FALSE]

## 4.7 保存 df0（补充材料可复现；若文件已存在则跳过）
df0_rds <- file.path(out_dir, "INTERMEDIATE_df0_exprGT0.rds")
df0_csv <- file.path(out_dir, "INTERMEDIATE_df0_exprGT0.csv")

if (!file.exists(df0_rds)) saveRDS(df0, df0_rds)
if (!file.exists(df0_csv)) write.csv(df0, df0_csv, row.names = FALSE)

sink(log_fp, append = TRUE)
cat("Expressed cells (expr>0):", nrow(df0), "\n")
cat("Donors in expressed cells (unique donor×group pairs):\n")
print(table(unique(df0[, c("donor","group4")])$group4))
cat("\n")
sink()

## =========================
## 5) 绘图函数：overlap hist（Y=Cell count）+ MWU(BH) + 右上角统计标签
## =========================
## unit="cell"：每个 cell 一个观测
## unit="donor"：每个 donor 在每个 celltype 先聚合成一个值（默认 median）
plot_one_unit <- function(unit = c("donor", "cell"),
                          donor_summary_fun = c("median","mean")) {
  
  unit <- match.arg(unit)
  donor_summary_fun <- match.arg(donor_summary_fun)
  
  ## 5.1 只保留 AD vs Control（和模版一致）
  df2 <- df0 %>%
    filter(group4 %in% c(disease, "Control")) %>%
    mutate(group = ifelse(group4 == disease, disease, "Control"))
  
  ## 5.2 legend 的 n：按“表达细胞 df0”口径统计 donor 数（而不是全体细胞）
  don_pair <- unique(df2[, c("donor","group")])
  n_dis <- sum(don_pair$group == disease)
  n_ctl <- sum(don_pair$group == "Control")
  
  ## legend 显示两行：组名 + (n = donor数)
  lab_dis <- sprintf("%s\n(n = %d)", disease,  n_dis)
  lab_ctl <- sprintf("Control\n(n = %d)",     n_ctl)
  
  ## 设置为 factor，保证 legend 顺序固定：AD 在上/左，Control 在下/右
  df2$group_lab <- factor(ifelse(df2$group == disease, lab_dis, lab_ctl),
                          levels = c(lab_dis, lab_ctl))
  
  ## 5.3 构建统计输入表（stat_input）
  ##     - donor-level：每个 celltype×donor×group → 一个 val（median/mean）
  ##     - cell-level：每个 cell 直接用 expr 作为 val
  if (unit == "donor") {
    if (donor_summary_fun == "median") {
      stat_input <- df2 %>%
        group_by(celltype6, donor, group_lab) %>%
        summarise(val = median(expr), .groups = "drop")
    } else {
      stat_input <- df2 %>%
        group_by(celltype6, donor, group_lab) %>%
        summarise(val = mean(expr), .groups = "drop")
    }
  } else {
    stat_input <- df2 %>%
      transmute(celltype6 = celltype6, group_lab = group_lab, val = expr)
  }
  
  ## 5.4 保存统计输入表（补充材料可以直接复算 MWU）
  stat_rds <- file.path(out_dir, paste0("INTERMEDIATE_stat_input_AD_vs_Control_", unit, ".rds"))
  stat_csv <- file.path(out_dir, paste0("INTERMEDIATE_stat_input_AD_vs_Control_", unit, ".csv"))
  
  if (!file.exists(stat_rds)) saveRDS(stat_input, stat_rds)
  if (!file.exists(stat_csv)) write.csv(stat_input, stat_csv, row.names = FALSE)
  
  ## 5.5 对每个 celltype 做 MWU，并进行 BH(FDR) 校正
  stats <- stat_input %>%
    group_by(celltype6) %>%
    summarise(
      p_raw = tryCatch(wilcox.test(val ~ group_lab, exact = FALSE)$p.value,
                       error = function(e) NA_real_),
      .groups = "drop"
    )
  
  stats$padj  <- p.adjust(stats$p_raw, method = "BH")
  stats$label <- sprintf("Mann–Whitney U\nPadj=%.2e", stats$padj)
  
  ## 为了把标签放在每个 facet 的右上角：
  ## - x=Inf, y=Inf
  ## - coord_cartesian(clip="off")
  ## - geom_label hjust/vjust > 1 稍微往外推
  stats$x <- Inf
  stats$y <- Inf
  
  ## 保存统计结果表（补充材料）
  stats_fp <- file.path(out_dir, paste0("STATS_AD_vs_Control_", unit, ".csv"))
  if (!file.exists(stats_fp)) write.csv(stats, stats_fp, row.names = FALSE)
  
  ## 5.6 颜色：严格红/蓝（与模版一致）
  fill_vals <- c("red","blue")
  names(fill_vals) <- c(lab_dis, lab_ctl)
  
  ## 5.7 构建 overlap 直方图（Y=Cell count）
  ## 注意：不写 y=after_stat(density)，默认就是 count
  p <- ggplot(df2, aes(x = expr, fill = group_lab)) +
    geom_histogram(binwidth = binwidth, alpha = 0.7, position = "identity", colour = NA) +
    facet_wrap(~ celltype6, scales = "free_y") +
    scale_fill_manual(values = fill_vals, drop = FALSE) +
    geom_label(
      data = stats, inherit.aes = FALSE,
      aes(x = x, y = y, label = label),
      hjust = 1.02, vjust = 1.02, size = 2.3,
      label.size = 0, fill = "white", alpha = 0.7
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
    theme(plot.margin = margin(10, 25, 10, 10)) +
    coord_cartesian(clip = "off")
  
  ## 5.8 输出文件名（与模版一致）
  out_png <- file.path(
    out_dir,
    paste0(dataset_tag, "_", gene_symbol,
           "_OverlapHistCount_AD_vs_Control_", unit, ".png")
  )
  
  ## 若已存在则跳过（你要求：已经生成的就跳过）
  if (file.exists(out_png)) {
    cat("↪ skip (already exists): ", out_png, "\n", sep = "")
    return(invisible(out_png))
  }
  
  ## 使用 ragg 输出高质量 png（论文/补充材料更稳定）
  ragg::agg_png(out_png, width = 10, height = 6, units = "in", res = 300, background = "white")
  print(p)
  dev.off()
  
  cat("✅ saved: ", out_png, "\n", sep = "")
  return(invisible(out_png))
}

## =========================
## 6) 输出两张图：cell-level + donor-level
## =========================
## 先 cell-level（每个 cell 一个观测）
plot_one_unit("cell")

## 再 donor-level（每个 donor 聚合为一个值：默认 median）
plot_one_unit("donor")  # donor_summary_fun 默认 median（与你要求一致）

## =========================
## 7) 把 sessionInfo 写入日志（补充材料常要求）
## =========================
sink(log_fp, append = TRUE)
cat("\n==== sessionInfo ====\n")
print(sessionInfo())
cat("==== END ====\n")
sink()

cat("\n🎉 DONE. Output dir:\n", out_dir, "\n", sep = "")







#第2章 2 ，Y轴 = Density，
###############################################################################
## NO2_GSE174367_OverlapHistDensity_AD_vs_Control_cell_and_donor_level.R
##
## 【第2章目标】（严格模仿 GSE157827 模版 / 你刚刚那两张图的风格）
##  - 使用 6 个细胞类型（celltype6）做 overlap 直方图
##  - Y轴 = Density（使用 after_stat(density)）
##  - 只使用“表达细胞”（expr > 0）
##  - 同时输出两张图：
##      (1) cell-level：每个 cell 作为观测 → Mann–Whitney U → BH(FDR)
##      (2) donor-level：每个 donor 先汇总(默认 median) → MWU → BH(FDR)
##
## 【UBL3表达计算】严格按模版：
##    expr = log1p( (raw_counts / lib_size) * 1e4 )   # CP10K
##    X轴显示：UBL3 log1p(CP10k)
##
## 【输入对象】强制只用：
##    D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/stepH_obj_celltype6_named.rds
##
## 【输出目录】全部保存在：
##    D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/NO2_overlap_hist_density_GSE157827style
##
## 【补充材料/中间结果】自动保存（若已存在则跳过）：
##  - 每组 donor 数：QC_donors_by_group4.csv
##  - 表达细胞 df0：INTERMEDIATE_df0_exprGT0.csv/.rds
##  - df2（AD vs Control 且 expr>0）：INTERMEDIATE_df2_AD_vs_Control_exprGT0.csv/.rds
##  - 统计输入表：INTERMEDIATE_stat_input_*.csv/.rds
##  - 每 celltype 的 p/padj 表：STATS_*.csv
##  - 每 celltype 的样本量检查表：CHECK_n_by_celltype_*.csv
##  - 日志与 sessionInfo：NO2_log.txt
###############################################################################

## =========================
## 0) 清理环境 + 基本设置
## =========================
rm(list = ls()); gc()
Sys.setenv(LANG = "en")

SEED <- 20251023
set.seed(SEED)
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
## 1) 路径参数（★强制只用 celltype6 对象★）
## =========================
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results"

obj_fp <- file.path(res_dir, "stepH_obj_celltype6_named.rds")
stopifnot(file.exists(obj_fp))

out_dir <- file.path(res_dir, "NO2_overlap_hist_density_GSE157827style")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dataset_tag <- "GSE174367"
gene_symbol <- "UBL3"
binwidth    <- 0.2
disease     <- "AD"

log_fp <- file.path(out_dir, "NO2_log.txt")
sink(log_fp)
cat("==== START ====\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("SEED:", SEED, "\n")
cat("obj_fp:", obj_fp, "\n")
cat("out_dir:", out_dir, "\n\n")
sink()

## =========================
## 2) 读对象 + 强制检查 celltype6 必须存在且只有 6 类
## =========================
obj <- readRDS(obj_fp)
DefaultAssay(obj) <- "RNA"
md <- obj@meta.data

if (!("celltype6" %in% colnames(md))) {
  sink(log_fp, append = TRUE)
  cat("❌ meta.data 列名如下（找不到 celltype6）：\n")
  print(colnames(md))
  sink()
  stop("❌ 该对象 meta.data 中缺少 celltype6 列。")
}

ct_levels <- sort(unique(as.character(md$celltype6)))
n_ct <- length(ct_levels)

sink(log_fp, append = TRUE)
cat("celltype6 unique =", n_ct, "\n")
print(ct_levels)
cat("\n")
sink()

if (n_ct != 6) stop("❌ celltype6 不是 6 类（请检查上游 celltype6 命名/过滤）。")

## 自动识别 group / donor 列（与你模版一致）
grp_col_candidates <- c("group4","group2","group","Group","diagnosis","Dx","clinical_diagnosis")
grp_col <- grp_col_candidates[grp_col_candidates %in% colnames(md)][1]
if (is.na(grp_col)) stop("❌ meta.data 找不到 group 列（group4/group/...）。")

donor_col_candidates <- c("autopsy_id","donor","Donor","sample","Sample","orig.ident","patient","subject")
don_col <- donor_col_candidates[donor_col_candidates %in% colnames(md)][1]
if (is.na(don_col)) stop("❌ meta.data 找不到 donor 列（autopsy_id/sample/...）。")

## 统一 Control/NC 名称
grp_raw <- trimws(as.character(md[[grp_col]]))
ctrl_alias <- c("Control","CTRL","Ctr","CTR","NC","Normal","N","control","ctrl","ctr","nc","normal")
grp_std <- ifelse(grp_raw %in% ctrl_alias, "Control", grp_raw)

## 写入标准列名（后续统一使用）
md$celltype6_std <- as.character(md$celltype6)
md$group4_std    <- grp_std
md$donor_std     <- trimws(as.character(md[[don_col]]))

if (!all(c(disease, "Control") %in% unique(md$group4_std))) {
  sink(log_fp, append = TRUE)
  cat("⚠ group unique after standardization:\n")
  print(sort(unique(md$group4_std)))
  sink()
  stop("❌ 分组里没有同时包含 AD 和 Control。")
}

## donor QC（补充材料）
don_all <- unique(md[, c("donor_std","group4_std")])
qc_don <- as.data.frame(table(don_all$group4_std), stringsAsFactors = FALSE)
colnames(qc_don) <- c("group4","n_donors")

qc_don_fp <- file.path(out_dir, "QC_donors_by_group4.csv")
if (!file.exists(qc_don_fp)) write.csv(qc_don, qc_don_fp, row.names = FALSE)

sink(log_fp, append = TRUE)
cat("Loaded cells:", ncol(obj), " genes:", nrow(obj), "\n")
cat("celltype6 fixed at 6 panels.\n")
cat("group col:", grp_col, " | donor col:", don_col, "\n\n")
sink()

## =========================
## 3) 获取 counts（Seurat v5 多 layers 兼容）
## =========================
get_counts_matrix_allcells <- function(obj, assay = "RNA") {
  a <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(a), error = function(e) character(0))
  counts_layers <- layers[grepl("^counts", layers)]
  
  if (length(counts_layers) > 0) {
    sink(log_fp, append = TRUE)
    cat("Detected counts layers:\n"); print(counts_layers)
    sink()
    
    mats <- list()
    for (ly in counts_layers) {
      m <- tryCatch(SeuratObject::LayerData(a, layer = ly), error = function(e) NULL)
      if (is.null(m)) next
      if (!is.null(dim(m)) && length(dim(m)) == 2) mats[[ly]] <- m
    }
    if (length(mats) == 0) stop("❌ counts layers 存在，但读取 LayerData 失败。")
    
    ref_genes <- rownames(mats[[1]])
    for (k in names(mats)) {
      if (!identical(rownames(mats[[k]]), ref_genes)) {
        m0 <- mats[[k]]
        m_aligned <- Matrix::Matrix(0, nrow = length(ref_genes), ncol = ncol(m0), sparse = TRUE)
        rownames(m_aligned) <- ref_genes
        colnames(m_aligned) <- colnames(m0)
        common <- intersect(ref_genes, rownames(m0))
        m_aligned[common, ] <- m0[common, , drop = FALSE]
        mats[[k]] <- m_aligned
      }
    }
    
    mat_all <- if (length(mats) == 1) mats[[1]] else Reduce(Matrix::cbind2, mats)
    
    dup <- duplicated(colnames(mat_all))
    if (any(dup)) {
      sink(log_fp, append = TRUE)
      cat("⚠ duplicated cells across layers:", sum(dup), " -> keep first\n")
      sink()
      mat_all <- mat_all[, !dup, drop = FALSE]
    }
    
    all_cells <- colnames(obj)
    miss_cells <- setdiff(all_cells, colnames(mat_all))
    if (length(miss_cells) > 0) {
      sink(log_fp, append = TRUE)
      cat("⚠ missing cells in merged counts:", length(miss_cells), " -> fill zeros\n")
      sink()
      m_fill <- Matrix::Matrix(0, nrow = nrow(mat_all), ncol = length(miss_cells), sparse = TRUE)
      rownames(m_fill) <- rownames(mat_all)
      colnames(m_fill) <- miss_cells
      mat_all <- Matrix::cbind2(mat_all, m_fill)
    }
    
    mat_all <- mat_all[, all_cells, drop = FALSE]
    return(mat_all)
  }
  
  m2 <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "counts"), error = function(e) NULL)
  if (!is.null(m2)) return(m2)
  
  stop("❌ 无法获取 counts。")
}

rna_counts <- get_counts_matrix_allcells(obj, "RNA")
stopifnot(ncol(rna_counts) == ncol(obj))

sink(log_fp, append = TRUE)
cat("Merged counts dim:", paste(dim(rna_counts), collapse = " x "), "\n\n")
sink()

## =========================
## 4) 计算 UBL3 log1p(CP10k) + df0(expr>0)
## =========================
gene_row <- if ("UBL3" %in% rownames(rna_counts)) "UBL3" else
  if ("ENSG00000122042" %in% rownames(rna_counts)) "ENSG00000122042" else NA_character_
if (is.na(gene_row)) stop("❌ counts 行名中找不到 UBL3（UBL3 或 ENSG00000122042）。")

lib_size <- Matrix::colSums(rna_counts)
raw_vec  <- as.numeric(rna_counts[gene_row, , drop = TRUE])
expr <- log1p((raw_vec / pmax(lib_size, 1)) * 1e4)

df_all <- data.frame(
  expr      = expr,
  donor     = md$donor_std,
  group4    = md$group4_std,
  celltype6 = md$celltype6_std,
  stringsAsFactors = FALSE
)

df0 <- df_all[df_all$expr > 0, , drop = FALSE]   # only expressed cells

df0_rds <- file.path(out_dir, "INTERMEDIATE_df0_exprGT0.rds")
df0_csv <- file.path(out_dir, "INTERMEDIATE_df0_exprGT0.csv")
if (!file.exists(df0_rds)) saveRDS(df0, df0_rds)
if (!file.exists(df0_csv)) write.csv(df0, df0_csv, row.names = FALSE)

sink(log_fp, append = TRUE)
cat("Expressed cells (expr>0):", nrow(df0), "\n\n")
sink()

## =========================
## 5) 绘图函数：Overlap Hist (Density) + MWU(BH) + 右上角统计标签
## =========================
plot_one_unit <- function(unit = c("donor","cell")) {
  
  unit <- match.arg(unit)
  
  ## 5.1 只做 AD vs Control
  df2 <- df0 %>%
    filter(group4 %in% c(disease, "Control")) %>%
    mutate(group = ifelse(group4 == disease, disease, "Control"))
  
  ## 保存 df2（便于补充材料追溯）
  df2_rds <- file.path(out_dir, "INTERMEDIATE_df2_AD_vs_Control_exprGT0.rds")
  df2_csv <- file.path(out_dir, "INTERMEDIATE_df2_AD_vs_Control_exprGT0.csv")
  if (!file.exists(df2_rds)) saveRDS(df2, df2_rds)
  if (!file.exists(df2_csv)) write.csv(df2, df2_csv, row.names = FALSE)
  
  ## 5.2 legend 的 n：用“表达细胞df0口径”统计 donor 数（与你上一章一致）
  don_pair <- unique(df2[, c("donor","group")])
  n_dis <- sum(don_pair$group == disease)
  n_ctl <- sum(don_pair$group == "Control")
  
  lab_dis <- sprintf("%s\n(n = %d)", disease,  n_dis)
  lab_ctl <- sprintf("Control\n(n = %d)",     n_ctl)
  
  ## ★关键：levels 必须与 fill_vals 名称一致，避免颜色退回灰色
  df2$group_lab <- factor(ifelse(df2$group == disease, lab_dis, lab_ctl),
                          levels = c(lab_dis, lab_ctl))
  
  ## 5.3 统计输入：donor-level vs cell-level（检验方法不变：MWU）
  if (unit == "donor") {
    stat_input <- df2 %>%
      group_by(celltype6, donor, group_lab) %>%
      summarise(val = median(expr), .groups = "drop")
  } else {
    stat_input <- df2 %>%
      transmute(celltype6 = celltype6, group_lab = group_lab, val = expr)
  }
  
  ## 保存统计输入（补充材料）
  stat_rds <- file.path(out_dir, paste0("INTERMEDIATE_stat_input_AD_vs_Control_", unit, ".rds"))
  stat_csv <- file.path(out_dir, paste0("INTERMEDIATE_stat_input_AD_vs_Control_", unit, ".csv"))
  if (!file.exists(stat_rds)) saveRDS(stat_input, stat_rds)
  if (!file.exists(stat_csv)) write.csv(stat_input, stat_csv, row.names = FALSE)
  
  ## 样本量检查表（每个 celltype×组）
  n_by_celltype <- stat_input %>%
    group_by(celltype6, group_lab) %>%
    summarise(n = n(), .groups = "drop")
  check_fp <- file.path(out_dir, paste0("CHECK_n_by_celltype_AD_vs_Control_", unit, ".csv"))
  if (!file.exists(check_fp)) write.csv(n_by_celltype, check_fp, row.names = FALSE)
  
  ## 5.4 MWU + BH(FDR)：每个 celltype 一次检验
  stats <- stat_input %>%
    group_by(celltype6) %>%
    summarise(
      p_raw = tryCatch(wilcox.test(val ~ group_lab, exact = FALSE)$p.value,
                       error = function(e) NA_real_),
      .groups = "drop"
    )
  stats$padj  <- p.adjust(stats$p_raw, method = "BH")
  stats$label <- sprintf("Mann–Whitney U\nPadj=%.2e", stats$padj)
  stats$x <- Inf; stats$y <- Inf
  
  stats_fp <- file.path(out_dir, paste0("STATS_AD_vs_Control_", unit, ".csv"))
  if (!file.exists(stats_fp)) write.csv(stats, stats_fp, row.names = FALSE)
  
  ## 5.5 颜色（严格红/蓝）
  fill_vals <- c("red", "blue")
  names(fill_vals) <- c(lab_dis, lab_ctl)
  
  ## 5.6 绘图：Y轴=Density（关键差异：y=after_stat(density)）
  p <- ggplot(df2, aes(x = expr, y = after_stat(density), fill = group_lab)) +
    geom_histogram(binwidth = binwidth, alpha = 0.7, position = "identity", colour = NA) +
    facet_wrap(~ celltype6, scales = "free_y") +
    scale_fill_manual(values = fill_vals, drop = FALSE) +
    geom_label(
      data = stats, inherit.aes = FALSE,
      aes(x = x, y = y, label = label),
      hjust = 1.02, vjust = 1.02, size = 2.3,
      label.size = 0, fill = "white", alpha = 0.7
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
    theme(plot.margin = margin(10, 25, 10, 10)) +
    coord_cartesian(clip = "off")
  
  ## 5.7 输出文件（已存在则跳过）
  out_png <- file.path(
    out_dir,
    paste0(dataset_tag, "_", gene_symbol,
           "_OverlapHistDensity_AD_vs_Control_", unit, ".png")
  )
  
  if (file.exists(out_png)) {
    cat("↪ skip (already exists): ", out_png, "\n", sep = "")
    return(invisible(out_png))
  }
  
  ragg::agg_png(out_png, width = 10, height = 6, units = "in", res = 300, background = "white")
  print(p)
  dev.off()
  
  cat("✅ saved: ", out_png, "\n", sep = "")
  
  sink(log_fp, append = TRUE)
  cat("---- AD vs Control | unit=", unit, "\n", sep = "")
  cat("Legend donor n: AD=", n_dis, " | Control=", n_ctl, "\n", sep = "")
  cat("Saved figure:", out_png, "\n\n")
  sink()
  
  invisible(out_png)
}

## =========================
## 6) 输出两张图：cell-level + donor-level
## =========================
plot_one_unit("cell")
plot_one_unit("donor")

## =========================
## 7) sessionInfo 写入日志（补充材料常要求）
## =========================
sink(log_fp, append = TRUE)
cat("\n==== sessionInfo ====\n")
print(sessionInfo())
cat("==== END ====\n")
sink()

cat("\n🎉 DONE. Output dir:\n", out_dir, "\n", sep = "")
cat("Log:\n", log_fp, "\n", sep = "")














###############################################################################
## NO2_GSE174367_SUMO_OverlapHistDensity_cellHist_donorMWU.R
##
## 【目标】（严格对齐你 GSE157827 的最新风格）
##  - 6 个细胞类型（celltype6）分面 overlap 直方图
##  - Y轴 = Density（after_stat(density)）
##  - 只使用“表达细胞”（expr > 0）
##  - 只输出 1 种图：
##      ✅ cell-level 直方图（用每个 cell 的 expr 画分布形态）
##      ✅ donor-level 统计（每 donor×celltype 先汇总 median → MWU → BH）
##    ※不再输出 donor-level 直方图（会稀疏且不利于形态对比）
##
## 【表达计算】严格按你模版：
##    expr = log1p( (raw_counts / lib_size) * 1e4 )   # CP10K
##    X轴显示：SUMO1 log1p(CP10k) / SUMO2 ... / SUMO3 ...
##
## 【输入对象】强制只用：
##    D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/stepH_obj_celltype6_named.rds
##
## 【输出目录】全部保存在：
##    D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/NO2_overlap_hist_density_GSE157827style_SUMOref
##
## 【中间结果/核查】每个基因都会输出：
##  - QC_donors_by_group4.csv（全局 donor 数）
##  - INTERMEDIATE_{GENE}_df0_exprGT0.csv/.rds（expr>0 全部细胞）
##  - INTERMEDIATE_{GENE}_df2_AD_vs_Control_exprGT0.csv/.rds（AD vs Control 且 expr>0）
##  - CHECK_{GENE}_n_cells_by_celltype_AD_vs_Control_exprGT0.csv（每 celltype 的表达细胞数）
##  - INTERMEDIATE_{GENE}_stat_input_donorMedian.csv/.rds（donor-level 统计输入）
##  - CHECK_{GENE}_n_donors_by_celltype_AD_vs_Control_exprGT0.csv（每 celltype 的 donor 数）
##  - STATS_{GENE}_MWU_donorMedian_BH_by_celltype.csv（p/padj）
##  - 图：{GENE}_OverlapHistDensity_AD_vs_Control_cellHist_MWUdonorMedian.png
##  - 日志与 sessionInfo：NO2_log_SUMO.txt
###############################################################################

## =========================
## 0) 清理环境 + 基本设置
## =========================
rm(list = ls()); gc()
Sys.setenv(LANG = "en")

SEED <- 20251023
set.seed(SEED)
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
## 1) 路径参数（★强制只用 celltype6 对象★）
## =========================
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results"
obj_fp  <- file.path(res_dir, "stepH_obj_celltype6_named.rds")
stopifnot(file.exists(obj_fp))

## 为避免覆盖你已有的 UBL3 NO2 输出，这里单独建 SUMOref 输出目录
out_dir <- file.path(res_dir, "NO2_overlap_hist_density_GSE157827style_SUMOref")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dataset_tag <- "GSE174367"
binwidth    <- 0.2
disease     <- "AD"

## 默认跑 SUMO1/2/3；如果你只想先跑 SUMO1，把 gene_list 改成 c("SUMO1")
gene_list <- c("SUMO1", "SUMO2", "SUMO3")

log_fp <- file.path(out_dir, "NO2_log_SUMO.txt")
sink(log_fp)
cat("==== START ====\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("SEED:", SEED, "\n")
cat("obj_fp:", obj_fp, "\n")
cat("out_dir:", out_dir, "\n")
cat("binwidth:", binwidth, "\n")
cat("genes:", paste(gene_list, collapse = ", "), "\n\n")
sink()

## =========================
## 2) 读对象 + 强制检查 celltype6 必须存在且只有 6 类
## =========================
obj <- readRDS(obj_fp)
DefaultAssay(obj) <- "RNA"
md <- obj@meta.data

if (!("celltype6" %in% colnames(md))) {
  sink(log_fp, append = TRUE)
  cat("❌ meta.data 列名如下（找不到 celltype6）：\n")
  print(colnames(md))
  sink()
  stop("❌ 该对象 meta.data 中缺少 celltype6 列。")
}

ct_levels <- sort(unique(as.character(md$celltype6)))
n_ct <- length(ct_levels)

sink(log_fp, append = TRUE)
cat("celltype6 unique =", n_ct, "\n")
print(ct_levels)
cat("\n")
sink()

if (n_ct != 6) stop("❌ celltype6 不是 6 类（请检查上游 celltype6 命名/过滤）。")

## 自动识别 group / donor 列（与你模版一致）
grp_col_candidates <- c("group4","group2","group","Group","diagnosis","Dx","clinical_diagnosis")
grp_col <- grp_col_candidates[grp_col_candidates %in% colnames(md)][1]
if (is.na(grp_col)) stop("❌ meta.data 找不到 group 列（group4/group2/...）。")

donor_col_candidates <- c("autopsy_id","donor","Donor","sample","Sample","orig.ident","patient","subject")
don_col <- donor_col_candidates[donor_col_candidates %in% colnames(md)][1]
if (is.na(don_col)) stop("❌ meta.data 找不到 donor 列（autopsy_id/sample/...）。")

## 统一 Control/NC 名称（保持口径一致）
grp_raw <- trimws(as.character(md[[grp_col]]))
ctrl_alias <- c("Control","CTRL","Ctr","CTR","NC","Normal","N","control","ctrl","ctr","nc","normal")
grp_std <- ifelse(grp_raw %in% ctrl_alias, "Control", grp_raw)

## 标准列名（后续统一使用）
md$celltype6_std <- as.character(md$celltype6)
md$group4_std    <- grp_std
md$donor_std     <- trimws(as.character(md[[don_col]]))

if (!all(c(disease, "Control") %in% unique(md$group4_std))) {
  sink(log_fp, append = TRUE)
  cat("⚠ group unique after standardization:\n")
  print(sort(unique(md$group4_std)))
  sink()
  stop("❌ 分组里没有同时包含 AD 和 Control。")
}

## donor QC（补充材料：全局 donor 数）
don_all <- unique(md[, c("donor_std","group4_std")])
qc_don <- as.data.frame(table(don_all$group4_std), stringsAsFactors = FALSE)
colnames(qc_don) <- c("group4","n_donors")
write.csv(qc_don, file.path(out_dir, "QC_donors_by_group4.csv"), row.names = FALSE)

sink(log_fp, append = TRUE)
cat("Loaded cells:", ncol(obj), " genes:", nrow(obj), "\n")
cat("celltype6 fixed at 6 panels.\n")
cat("group col:", grp_col, " | donor col:", don_col, "\n\n")
sink()

## =========================
## 3) 获取 counts（Seurat v5 多 layers 兼容）
## =========================
get_counts_matrix_allcells <- function(obj, assay = "RNA") {
  a <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(a), error = function(e) character(0))
  counts_layers <- layers[grepl("^counts", layers)]
  
  if (length(counts_layers) > 0) {
    sink(log_fp, append = TRUE)
    cat("Detected counts layers:\n"); print(counts_layers)
    sink()
    
    mats <- list()
    for (ly in counts_layers) {
      m <- tryCatch(SeuratObject::LayerData(a, layer = ly), error = function(e) NULL)
      if (is.null(m)) next
      if (!is.null(dim(m)) && length(dim(m)) == 2) mats[[ly]] <- m
    }
    if (length(mats) == 0) stop("❌ counts layers 存在，但读取 LayerData 失败。")
    
    ## 对齐 gene 顺序
    ref_genes <- rownames(mats[[1]])
    for (k in names(mats)) {
      if (!identical(rownames(mats[[k]]), ref_genes)) {
        m0 <- mats[[k]]
        m_aligned <- Matrix::Matrix(0, nrow = length(ref_genes), ncol = ncol(m0), sparse = TRUE)
        rownames(m_aligned) <- ref_genes
        colnames(m_aligned) <- colnames(m0)
        common <- intersect(ref_genes, rownames(m0))
        m_aligned[common, ] <- m0[common, , drop = FALSE]
        mats[[k]] <- m_aligned
      }
    }
    
    mat_all <- if (length(mats) == 1) mats[[1]] else Reduce(Matrix::cbind2, mats)
    
    ## 去掉重复 cell
    dup <- duplicated(colnames(mat_all))
    if (any(dup)) {
      sink(log_fp, append = TRUE)
      cat("⚠ duplicated cells across layers:", sum(dup), " -> keep first\n")
      sink()
      mat_all <- mat_all[, !dup, drop = FALSE]
    }
    
    ## 补齐缺失 cell（用 0 填充）
    all_cells <- colnames(obj)
    miss_cells <- setdiff(all_cells, colnames(mat_all))
    if (length(miss_cells) > 0) {
      sink(log_fp, append = TRUE)
      cat("⚠ missing cells in merged counts:", length(miss_cells), " -> fill zeros\n")
      sink()
      m_fill <- Matrix::Matrix(0, nrow = nrow(mat_all), ncol = length(miss_cells), sparse = TRUE)
      rownames(m_fill) <- rownames(mat_all)
      colnames(m_fill) <- miss_cells
      mat_all <- Matrix::cbind2(mat_all, m_fill)
    }
    
    mat_all <- mat_all[, all_cells, drop = FALSE]
    return(mat_all)
  }
  
  ## v4 兜底
  m2 <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "counts"), error = function(e) NULL)
  if (!is.null(m2)) return(m2)
  
  stop("❌ 无法获取 counts。")
}

rna_counts <- get_counts_matrix_allcells(obj, "RNA")
stopifnot(ncol(rna_counts) == ncol(obj))

sink(log_fp, append = TRUE)
cat("Merged counts dim:", paste(dim(rna_counts), collapse = " x "), "\n\n")
sink()

lib_size <- Matrix::colSums(rna_counts)

## =========================
## 4) 基因行名定位（不“瞎编”ENSG；仅在可用时尝试 SYMBOL→ENSEMBL）
## =========================
locate_gene_row <- function(counts_mat, gene_symbol) {
  rn <- rownames(counts_mat)
  
  ## 1) 直接匹配（rownames 是 SYMBOL）
  if (gene_symbol %in% rn) return(gene_symbol)
  
  ## 2) 大小写不敏感匹配
  idx_ci <- which(toupper(rn) == toupper(gene_symbol))
  if (length(idx_ci) == 1) return(rn[idx_ci[1]])
  
  ## 3) 若 rownames 是 ENSG：且你安装了 org.Hs.eg.db，则尝试 SYMBOL→ENSEMBL
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    ens_tbl <- tryCatch(
      AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                            keys = gene_symbol, keytype = "SYMBOL", columns = c("ENSEMBL")),
      error = function(e) NULL
    )
    if (!is.null(ens_tbl) && "ENSEMBL" %in% colnames(ens_tbl)) {
      ens_ids <- unique(ens_tbl$ENSEMBL[!is.na(ens_tbl$ENSEMBL)])
      
      ## 精确 ENSG
      m1 <- intersect(ens_ids, rn)
      if (length(m1) >= 1) return(m1[1])
      
      ## ENSG.版本号
      rn_strip <- sub("\\.\\d+$", "", rn)
      hit <- which(rn_strip %in% ens_ids)
      if (length(hit) >= 1) return(rn[hit[1]])
    }
  }
  
  cand <- rn[grep(gene_symbol, rn, ignore.case = TRUE)]
  stop(paste0(
    "❌ counts 行名中找不到基因：", gene_symbol, "\n",
    "请检查 rownames 是否为 SYMBOL 或 ENSG。\n",
    "grep 候选（前 20 个）：", paste(head(cand, 20), collapse = ", ")
  ))
}

## =========================
## 5) 单基因流程：cell-level 直方图 + donor-level MWU(BH)
## =========================
run_one_gene <- function(gene_symbol) {
  
  sink(log_fp, append = TRUE)
  cat("==== Gene:", gene_symbol, "====\n")
  sink()
  
  gene_row <- locate_gene_row(rna_counts, gene_symbol)
  sink(log_fp, append = TRUE); cat("gene_row used:", gene_row, "\n"); sink()
  
  ## 5.1 计算 expr = log1p(CP10k)
  raw_vec <- as.numeric(rna_counts[gene_row, , drop = TRUE])
  expr <- log1p((raw_vec / pmax(lib_size, 1)) * 1e4)
  
  df_all <- data.frame(
    expr      = expr,
    donor     = md$donor_std,
    group4    = md$group4_std,
    celltype6 = md$celltype6_std,
    stringsAsFactors = FALSE
  )
  
  ## 只保留表达阳性细胞（expr>0）
  df0 <- df_all[df_all$expr > 0, , drop = FALSE]
  saveRDS(df0, file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_df0_exprGT0.rds")))
  write.csv(df0, file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_df0_exprGT0.csv")), row.names = FALSE)
  
  ## AD vs Control
  df2 <- df0 %>%
    filter(group4 %in% c(disease, "Control")) %>%
    mutate(group = ifelse(group4 == disease, disease, "Control"))
  
  saveRDS(df2, file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_df2_AD_vs_Control_exprGT0.rds")))
  write.csv(df2, file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_df2_AD_vs_Control_exprGT0.csv")), row.names = FALSE)
  
  ## 5.2 legend donor n（按“有表达细胞”的 donor 计数，口径与 GSE157827 一致）
  don_pair <- unique(df2[, c("donor", "group")])
  n_dis <- sum(don_pair$group == disease)
  n_ctl <- sum(don_pair$group == "Control")
  
  lab_dis <- sprintf("%s\n(n = %d)", disease,  n_dis)
  lab_ctl <- sprintf("Control\n(n = %d)",     n_ctl)
  
  df2$group_lab <- factor(ifelse(df2$group == disease, lab_dis, lab_ctl),
                          levels = c(lab_dis, lab_ctl))
  
  ## 5.3 细胞数核查（每 celltype×组：表达细胞数量）
  n_cells <- df2 %>%
    group_by(celltype6, group_lab) %>%
    summarise(n_cells_exprGT0 = n(), .groups = "drop")
  write.csv(n_cells,
            file.path(out_dir, paste0("CHECK_", gene_symbol, "_n_cells_by_celltype_AD_vs_Control_exprGT0.csv")),
            row.names = FALSE)
  
  ## 5.4 donor-level 统计输入：每 donor×celltype 汇总 median(expr)
  stat_input <- df2 %>%
    group_by(celltype6, donor, group_lab) %>%
    summarise(val = median(expr), .groups = "drop")
  
  saveRDS(stat_input,
          file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_stat_input_donorMedian.rds")))
  write.csv(stat_input,
            file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_stat_input_donorMedian.csv")),
            row.names = FALSE)
  
  ## donor 数核查（每 celltype×组：进入统计的 donor 数）
  n_don <- stat_input %>%
    group_by(celltype6, group_lab) %>%
    summarise(n_donors_with_exprGT0 = n(), .groups = "drop")
  write.csv(n_don,
            file.path(out_dir, paste0("CHECK_", gene_symbol, "_n_donors_by_celltype_AD_vs_Control_exprGT0.csv")),
            row.names = FALSE)
  
  ## 5.5 MWU + BH（严格 donor-level median）
  stats <- stat_input %>%
    group_by(celltype6) %>%
    summarise(
      p_raw = {
        g <- group_lab; v <- val
        if (length(unique(g)) < 2) NA_real_
        else tryCatch(wilcox.test(v ~ g, exact = FALSE)$p.value, error = function(e) NA_real_)
      },
      .groups = "drop"
    )
  
  stats$padj  <- p.adjust(stats$p_raw, method = "BH")
  stats$label <- sprintf("Mann–Whitney U\nPadj=%.2e", stats$padj)
  stats$x <- Inf; stats$y <- Inf
  
  write.csv(stats,
            file.path(out_dir, paste0("STATS_", gene_symbol, "_MWU_donorMedian_BH_by_celltype.csv")),
            row.names = FALSE)
  
  ## 5.6 绘图：cell-level overlap hist（Y=Density）+ donor-level 标签
  fill_vals <- c("red", "blue")
  names(fill_vals) <- c(lab_dis, lab_ctl)
  
  p <- ggplot(df2, aes(x = expr, y = after_stat(density), fill = group_lab)) +
    geom_histogram(binwidth = binwidth, alpha = 0.7, position = "identity", colour = NA) +
    facet_wrap(~celltype6, scales = "free_y") +
    scale_fill_manual(values = fill_vals, drop = FALSE) +
    geom_label(
      data = stats, inherit.aes = FALSE,
      aes(x = x, y = y, label = label),
      hjust = 1.02, vjust = 1.02, size = 2.3,
      label.size = 0, fill = "white", alpha = 0.7
    ) +
    labs(
      title = paste0(gene_symbol,
                     " expression per cell type (only expressed cells): ",
                     disease, " vs Control (donor-level)"),
      x = paste0(gene_symbol, " log1p(CP10k)"),
      y = "Density",
      fill = "group_lab"
    ) +
    theme_bw() +
    theme(plot.margin = margin(10, 25, 10, 10)) +
    coord_cartesian(clip = "off")
  
  out_png <- file.path(out_dir,
                       paste0(dataset_tag, "_", gene_symbol,
                              "_OverlapHistDensity_AD_vs_Control_cellHist_MWUdonorMedian.png"))
  
  ragg::agg_png(out_png, width = 10, height = 6, units = "in", res = 300, background = "white")
  print(p)
  dev.off()
  
  sink(log_fp, append = TRUE)
  cat("Saved figure:", out_png, "\n")
  cat("Legend donor n (expr>0 donors): AD=", n_dis, " | Control=", n_ctl, "\n", sep = "")
  cat("Expressed cells (expr>0):", nrow(df0), "\n\n")
  sink()
  
  invisible(out_png)
}

## =========================
## 6) 批量运行 SUMO1/2/3 + 汇总统计表
## =========================
stats_long <- list()

for (g in gene_list) {
  run_one_gene(g)
  st_fp <- file.path(out_dir, paste0("STATS_", g, "_MWU_donorMedian_BH_by_celltype.csv"))
  st <- read.csv(st_fp, stringsAsFactors = FALSE)
  st$gene <- g
  stats_long[[g]] <- st
}

stats_long_df <- bind_rows(stats_long) %>%
  select(gene, celltype6, p_raw, padj, label)

write.csv(stats_long_df,
          file.path(out_dir, "SUMMARY_SUMO_genes_MWU_BH_by_celltype.csv"),
          row.names = FALSE)

## =========================
## 7) sessionInfo 写入日志
## =========================
sink(log_fp, append = TRUE)
cat("\n==== sessionInfo ====\n")
print(sessionInfo())
cat("==== END ====\n")
sink()

cat("\n🎉 DONE. Output dir:\n", out_dir, "\n", sep = "")
cat("Log:\n", log_fp, "\n", sep = "")






#第3章，箱线图
###############################################################################
## NO3_GSE174367_UBL3_Boxplots_DESeq2_byDonor_ADvsControl_FINAL.R
##
## 【数据集】GSE174367：AD vs Control（Control/NC/CTRL 等会统一映射为 Control）
## 【统计单位】donor（先 pseudo-bulk 到 celltype6×donor，再 DESeq2）
##
## 【本章做什么】（严格模仿 GSE157827 模版风格：颜色/标题/点+箱线/字幕log2FC+padj）
##  1) 读取 stepH_obj_celltype6_named.rds（★强制使用 celltype6，且必须正好 6 类）
##  2) 自动识别 group 列 + donor 列；统一 group 为 AD/Control；仅保留这两组
##  3) donor -> group 一对一校验（不满足直接 stop，并输出检查表）
##  4) pseudo-bulk：按 celltype6×donor 聚合 counts（稀疏矩阵乘法，速度快、内存省）
##  5) 对每个 celltype6：
##     - DESeq2：AD vs Control（输出全基因 DEG 表）
##     - 画 UBL3 donor-level 箱线图（Normalized counts；四分位数箱体；点为 donor）
##  6) 输出 6-panel 总图（2×3），风格对齐 GSE157827
##  7) 保存所有关键中间结果，便于核查与补充材料投稿
##
## 【输出目录】全部保存在：
##    D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/NO3_GSE174367_UBL3_Boxplots_DESeq2_byDonor_ADvsControl
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
## 0) 路径（★强制只用 celltype6 对象★）
## =========================
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results"
obj_fp  <- file.path(res_dir, "stepH_obj_celltype6_named.rds")
stopifnot(file.exists(obj_fp))

out_dir <- file.path(res_dir, "NO3_GSE174367_UBL3_Boxplots_DESeq2_byDonor_ADvsControl")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_fp <- file.path(out_dir, "NO3_log.txt")
sink(log_fp)
cat("==== START ====\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("SEED:", SEED, "\n")
cat("obj :", obj_fp, "\n")
cat("out :", out_dir, "\n\n")
sink()

dataset_tag <- "GSE174367"
gene_symbol <- "UBL3"
disease     <- "AD"

## 颜色：与模版一致（红=AD，蓝=Control）
pal2 <- c(AD = "#D24B40", Control = "#2C7FB8")

## =========================
## 1) 读对象 + meta
## =========================
obj <- readRDS(obj_fp)
DefaultAssay(obj) <- "RNA"
md0 <- obj@meta.data

## =========================
## 2) 强制 celltype6 存在且正好 6 类（确保 6-panel）
## =========================
if (!("celltype6" %in% colnames(md0))) {
  sink(log_fp, append = TRUE)
  cat("❌ meta.data 列名如下（找不到 celltype6）：\n")
  print(colnames(md0))
  sink()
  stop("❌ 该对象缺少 celltype6 列。")
}

ct_levels <- sort(unique(trimws(as.character(md0$celltype6))))
n_ct <- length(ct_levels)

sink(log_fp, append = TRUE)
cat("celltype6 unique =", n_ct, "\n")
print(ct_levels)
cat("\n")
sink()

if (n_ct != 6) stop("❌ celltype6 不是 6 类，请检查上游 celltype6 命名/过滤。")

md0$celltype6_std <- trimws(as.character(md0$celltype6))

## =========================
## 3) 自动识别 group 列，并统一成 AD/Control（只保留这两组）
## =========================
grp_candidates <- c("group4","group2","group","Group","diagnosis","Dx","clinical_diagnosis")
grp_col <- grp_candidates[grp_candidates %in% colnames(md0)][1]
if (is.na(grp_col)) {
  stop("❌ 找不到分组列（group4/group2/group/diagnosis/...）。当前列名：\n",
       paste(colnames(md0), collapse = ", "))
}

grp_raw <- trimws(as.character(md0[[grp_col]]))
ctrl_alias <- c("Control","CTRL","Ctr","CTR","NC","Normal","N","control","ctrl","ctr","nc","normal")
md0$group_std <- ifelse(grp_raw %in% ctrl_alias, "Control", grp_raw)

## 只保留 AD/Control（其它组直接剔除，避免污染）
keep <- md0$group_std %in% c("AD","Control")
md  <- md0[keep, , drop = FALSE]

## =========================
## 4) 自动识别 donor 列（优先 autopsy_id）
## =========================
don_candidates <- c("autopsy_id","donor","Donor","sample","Sample","orig.ident","patient","subject")
don_col <- don_candidates[don_candidates %in% colnames(md)][1]
if (is.na(don_col)) {
  stop("❌ 找不到 donor 列（autopsy_id/donor/sample/...）。当前列名：\n",
       paste(colnames(md), collapse = ", "))
}

md$donor_std <- trimws(as.character(md[[don_col]]))

## 去掉空 donor / 空 celltype
md <- md[md$donor_std != "" & !is.na(md$donor_std) &
           md$celltype6_std != "" & !is.na(md$celltype6_std), , drop = FALSE]

## 记录：本次到底用了哪几列
sink(log_fp, append = TRUE)
cat("Used columns:\n")
cat("  group column   =", grp_col, "\n")
cat("  donor column   =", don_col, "\n")
cat("  celltype column= celltype6\n\n")
cat("Group counts (cells):\n"); print(table(md$group_std))
cat("Celltype counts (cells):\n"); print(table(md$celltype6_std))
cat("Unique donors:", length(unique(md$donor_std)), "\n\n")
sink()

## =========================
## 5) donor -> group 一对一校验（必须！）
## =========================
donor_map <- unique(md[, c("donor_std","group_std")])
if (any(table(donor_map$donor_std) > 1)) {
  bad <- donor_map[donor_map$donor_std %in% names(which(table(donor_map$donor_std) > 1)), ]
  write.csv(bad, file.path(out_dir, "CHECK_donor_maps_to_multiple_group.csv"), row.names = FALSE)
  stop("❌ donor 对应多个 group，已输出 CHECK_donor_maps_to_multiple_group.csv")
}

## donor 数（补充材料）
qc_donor_counts_by_group <- as.data.frame(table(donor_map$group_std), stringsAsFactors = FALSE)
colnames(qc_donor_counts_by_group) <- c("group", "n_donors")
write.csv(qc_donor_counts_by_group,
          file.path(out_dir, "QC_donor_counts_by_group.csv"),
          row.names = FALSE)

## 每个 celltype6 中，每组有多少 donor（补充材料）
tmp_donor_ct <- unique(md[, c("donor_std","celltype6_std","group_std")])
qc_donor_counts_by_celltype <- as.data.frame(table(tmp_donor_ct$celltype6_std, tmp_donor_ct$group_std),
                                             stringsAsFactors = FALSE)
colnames(qc_donor_counts_by_celltype) <- c("celltype6","group","n_donors")
write.csv(qc_donor_counts_by_celltype,
          file.path(out_dir, "QC_donor_counts_by_celltype6_by_group.csv"),
          row.names = FALSE)

## 每个 celltype6 中，每组有多少 cells（补充材料）
qc_cells_by_celltype <- as.data.frame(table(md$celltype6_std, md$group_std), stringsAsFactors = FALSE)
colnames(qc_cells_by_celltype) <- c("celltype6","group","n_cells")
write.csv(qc_cells_by_celltype,
          file.path(out_dir, "QC_cells_by_celltype6_by_group.csv"),
          row.names = FALSE)

## =========================
## 6) 获取 counts（兼容 Seurat v5 多 layers）
## =========================
get_counts <- function(obj, assay = "RNA") {
  a <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(a), error = function(e) character(0))
  cl <- layers[grepl("^counts", layers)]
  
  if (length(cl) > 0) {
    mats <- list()
    for (x in cl) {
      m <- tryCatch(SeuratObject::LayerData(a, layer = x), error = function(e) NULL)
      if (!is.null(m)) mats[[x]] <- m
    }
    if (length(mats) == 0) stop("❌ counts layers 存在，但读取 LayerData 失败。")
    
    ref <- rownames(mats[[1]])
    mats <- lapply(mats, function(m){
      if (identical(rownames(m), ref)) return(m)
      m2 <- Matrix::Matrix(0, nrow = length(ref), ncol = ncol(m), sparse = TRUE)
      rownames(m2) <- ref; colnames(m2) <- colnames(m)
      common <- intersect(ref, rownames(m))
      m2[common, ] <- m[common, , drop = FALSE]
      m2
    })
    
    m <- if (length(mats) == 1) mats[[1]] else Reduce(Matrix::cbind2, mats)
    
    ## 去重（防止 layer 合并后重复 cell）
    dup <- duplicated(colnames(m))
    if (any(dup)) m <- m[, !dup, drop = FALSE]
    
    ## 强制对齐 obj cell 顺序
    m[, colnames(obj), drop = FALSE]
  } else {
    Seurat::GetAssayData(obj, assay = assay, slot = "counts")
  }
}

cnt_all <- get_counts(obj, "RNA")

## 只取本章保留的 cells（AD/Control 且 meta 已过滤）
## 注意：md 的行名应是 cell barcode；确保与 counts colnames 匹配
stopifnot(all(rownames(md) %in% colnames(cnt_all)))
cnt <- cnt_all[, rownames(md), drop = FALSE]

sink(log_fp, append = TRUE)
cat("Counts dim (genes x kept_cells):", paste(dim(cnt), collapse = " x "), "\n\n")
sink()

## =========================
## 7) 确认 UBL3 的行名（symbol 或 ENSG）
## =========================
gene_row <- if ("UBL3" %in% rownames(cnt)) {
  "UBL3"
} else if ("ENSG00000122042" %in% rownames(cnt)) {
  "ENSG00000122042"
} else {
  NA_character_
}
if (is.na(gene_row)) stop("❌ counts 里找不到 UBL3（UBL3 或 ENSG00000122042）。")

write.csv(data.frame(gene = "UBL3", gene_row = gene_row),
          file.path(out_dir, "CHECK_geneSymbol_to_rowname.csv"),
          row.names = FALSE)

## =========================
## 8) pseudo-bulk：celltype6 × donor 聚合 counts（稀疏聚合）
## =========================
## 8.1 为每个 cell 定义 pseudo-bulk 组键：celltype6__donor
pb_key <- paste(md$celltype6_std, md$donor_std, sep = "__")

## 8.2 构建“cell -> pseudo-bulk列”的稀疏指示矩阵 M
##     M 的行=细胞，列=每个 celltype6×donor 的组合
grp <- factor(pb_key, levels = unique(pb_key))

M <- Matrix::sparseMatrix(
  i = seq_along(grp),
  j = as.integer(grp),
  x = 1,
  dims = c(length(grp), length(levels(grp))),
  dimnames = list(rownames(md), levels(grp))
)

## 8.3 聚合：gene×cell 乘 cell×pb -> gene×pb
pb <- cnt %*% M
pb <- as(pb, "dgCMatrix")

## 8.4 构建 pseudo-bulk 的列注释（coldata）
pb_meta <- data.frame(
  key       = colnames(pb),
  celltype6 = sub("__.*", "", colnames(pb)),
  donor     = sub(".*__", "", colnames(pb)),
  group     = donor_map$group_std[match(sub(".*__", "", colnames(pb)), donor_map$donor_std)],
  stringsAsFactors = FALSE
)

## 中间结果保存（补充材料/复现）
saveRDS(pb, file.path(out_dir, "INTERMEDIATE_pseudobulk_matrix_celltype6_byDonor.rds"))
write.csv(pb_meta, file.path(out_dir, "INTERMEDIATE_pseudobulk_coldata_celltype6_byDonor.csv"), row.names = FALSE)

## =========================
## 9) 每个 celltype6：DESeq2 + UBL3 箱线图 + 全基因 DEG 表
## =========================
plots <- list()

## 为了让 6-panel 排版稳定：使用 celltype6 当前的水平顺序
celltype_order <- ct_levels

for (ct in celltype_order) {
  
  ## 9.1 取出该 celltype 的 pseudo-bulk 列
  cols <- pb_meta$key[pb_meta$celltype6 == ct]
  if (length(cols) < 4) {
    ## donor 太少时 DESeq2 可能不稳，这里仍然记录并跳过
    sink(log_fp, append = TRUE)
    cat("⚠ Skip celltype (too few pseudo-bulk columns):", ct, " n=", length(cols), "\n")
    sink()
    next
  }
  
  ## 9.2 DESeq2 的 colData：group（Control 为 reference）
  coldata <- data.frame(
    group = factor(pb_meta$group[pb_meta$celltype6 == ct], levels = c("Control", "AD")),
    row.names = cols
  )
  
  ## 9.3 counts：必须是整数矩阵
  y <- round(as.matrix(pb[, cols, drop = FALSE]))
  storage.mode(y) <- "integer"
  
  ## 9.4 跑 DESeq2
  dds <- DESeqDataSetFromMatrix(y, coldata, design = ~ group)
  dds <- DESeq(dds, quiet = TRUE)
  
  ## 9.5 输出全基因 DEG（补充材料投稿常用）
  res_all <- results(dds, contrast = c("group", "AD", "Control"))
  res_df <- as.data.frame(res_all)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[, c("gene","log2FoldChange","padj","pvalue","baseMean")]
  
  deg_fp <- file.path(out_dir, paste0("DEG_", dataset_tag, "_", ct, "_AD_vs_Control.csv"))
  write.csv(res_df, deg_fp, row.names = FALSE)
  
  ## 9.6 取 UBL3 的 log2FC / padj（用于图 subtitle）
  u_log2fc <- as.numeric(res_all[gene_row, "log2FoldChange"])
  u_padj   <- as.numeric(res_all[gene_row, "padj"])
  
  subtitle <- sprintf("AD vs Control : log2FC=%.3f, padj=%s",
                      u_log2fc,
                      ifelse(is.na(u_padj), "NA", format(u_padj, digits = 3, scientific = TRUE)))
  
  ## 9.7 作图数据：使用 DESeq2 的 normalized counts（donor-level）
  norm <- DESeq2::counts(dds, normalized = TRUE)
  
  dfp <- data.frame(
    donor = pb_meta$donor[pb_meta$celltype6 == ct],
    group = factor(pb_meta$group[pb_meta$celltype6 == ct], levels = c("AD","Control")),
    value = as.numeric(norm[gene_row, cols]),
    stringsAsFactors = FALSE
  )
  
  ## 保存每张图的作图数据（补充材料/复现）
  plotdat_fp <- file.path(out_dir, paste0("INTERMEDIATE_plotdata_", dataset_tag, "_UBL3_", ct, ".csv"))
  write.csv(dfp, plotdat_fp, row.names = FALSE)
  
  ## 9.8 画箱线图（四分位箱体）+ donor 点（与模版一致）
  ##     - box：geom_boxplot（默认四分位）
  ##     - 点：geom_point + jitter（每个点 = donor）
  ##     - 颜色：AD 红、Control 蓝
  ##     - 标题：UBL3 in {ct} (pseudo-bulk per donor)
  p <- ggplot(dfp, aes(x = group, y = value, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.95, colour = "grey15") +
    geom_point(
      position = position_jitter(width = 0.10),
      size = 2.6, alpha = 0.9, shape = 21, stroke = 0.4, colour = "grey10"
    ) +
    scale_fill_manual(values = pal2, drop = FALSE) +
    labs(
      title = paste0("UBL3 in ", ct, " (pseudo-bulk per donor)"),
      subtitle = subtitle,
      y = "Normalized counts",
      x = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title.position = "plot",
      plot.title    = element_text(face = "bold", size = 22, margin = margin(b = 4)),
      plot.subtitle = element_text(size = 12, colour = "grey15", margin = margin(b = 8)),
      legend.position = "none"
    )
  
  ## 9.9 输出单张图（已存在则跳过，避免重复生成）
  out_png <- file.path(out_dir, paste0("UBL3_Box_byDonor_", dataset_tag, "_", ct, ".png"))
  if (!file.exists(out_png)) {
    ragg::agg_png(out_png, width = 10, height = 7, units = "in", res = 300, background = "white")
    print(p)
    dev.off()
  }
  
  plots[[ct]] <- p
  
  sink(log_fp, append = TRUE)
  cat("Saved:", out_png, "\n")
  cat("DEG  :", deg_fp, "\n")
  cat("Plot data:", plotdat_fp, "\n\n")
  sink()
}

## =========================
## 10) 6-panel 总图（2×3）——风格对齐 GSE157827
## =========================
panel_order <- celltype_order[celltype_order %in% names(plots)]
if (length(panel_order) != 6) {
  sink(log_fp, append = TRUE)
  cat("⚠ panel_order length != 6. Available plots:\n")
  print(names(plots))
  cat("panel_order used:\n")
  print(panel_order)
  sink()
  stop("❌ 6-panel 总图无法生成：某些 celltype 没有成功生成图（常见原因：该 celltype donor 太少）。")
}

p_all <- patchwork::wrap_plots(plots[panel_order], ncol = 3)

out_panel <- file.path(out_dir, paste0("UBL3_Box_byDonor_", dataset_tag, "_6panel.png"))
if (!file.exists(out_panel)) {
  ragg::agg_png(out_panel, width = 16, height = 9, units = "in", res = 300, background = "white")
  print(p_all)
  dev.off()
}

sink(log_fp, append = TRUE)
cat("Saved 6-panel:", out_panel, "\n\n")
cat("==== sessionInfo ====\n")
print(sessionInfo())
cat("==== END ====\n")
sink()

cat("\n🎉 DONE\nOutput dir:\n", out_dir, "\n", sep = "")
cat("Log:\n", log_fp, "\n", sep = "")







#
###############################################################################
## NO4_GSE174367_WholeCell_ADvsControl_perDonor_UBL3_boxplot_Assay5Layers.R
##
## 【本章目标】
##  - 统计单位：donor（每个 donor = 1 个点）
##  - 表达计算：counts -> CP10k -> log1p（✅严格沿用你给的模版）
##      CP10k = (UBL3_counts / lib_size) * 1e4
##      log1pCP10k = log1p(CP10k)
##  - donor 汇总：对每个 donor 计算 mean_cp10k 与 mean_log1pCP10k
##  - 检验：Wilcoxon（对 mean_log1pCP10k） + log2FC（对 mean_cp10k 的均值）
##  - 画图：箱线图（四分位数）+ donor 散点；颜色/标题/布局与模版一致
##  - Seurat v5 / Assay5 多 layers：直接遍历 counts.* layers 拼回全细胞向量（不 JoinLayers）
##
## 【输入对象】（你这次固定用 celltype6 对象也可以做 whole-cell）
##   D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/stepH_obj_celltype6_named.rds
##
## 【输出目录】：
##   D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/NO4_WholeCell_ADvsControl_perDonor
##
## 【补充材料/中间结果】会输出：
##   - UBL3_wholecell_perDonor_CP10k_mean_ADvsControl.csv（每 donor 1 行）
##   - UBL3_wholecell_stats_ADvsControl.csv（统计结果）
##   - CHECK_counts_layers.csv（counts layers 列表）
##   - CHECK_geneSymbol_to_rowname.csv（UBL3 行名）
##   - NO4_log.txt（日志+sessionInfo）
###############################################################################

rm(list = ls()); gc()

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
  library(ragg)
  library(tidyr)
})

## ========= 1) 路径（只改这里） =========
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results"

## ★本章：GSE174367 固定用这个对象（你说要继续用 stepH_obj_celltype6_named.rds）
obj_fp <- file.path(res_dir, "stepH_obj_celltype6_named.rds")
stopifnot(file.exists(obj_fp))

out_dir <- file.path(res_dir, "NO4_WholeCell_ADvsControl_perDonor")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_fp <- file.path(out_dir, "NO4_log.txt")
sink(log_fp)
cat("==== START ====\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("SEED:", SEED, "\n")
cat("obj_fp:", obj_fp, "\n")
cat("out_dir:", out_dir, "\n\n")
sink()

## ========= 2) 读入对象 + 最小自检 =========
obj <- readRDS(obj_fp)
DefaultAssay(obj) <- "RNA"
message("Loaded: ", basename(obj_fp), " | cells=", ncol(obj))

md <- obj@meta.data

## 自动识别：组别列 + donor列（与模版一致）
group_candidates <- c("group", "group2", "group4", "diagnosis", "Dx", "clinical_diagnosis")
donor_candidates <- c("autopsy_id", "sample", "donor", "individual", "subject", "subj", "case", "patient_id", "orig.ident")

group_col <- group_candidates[group_candidates %in% colnames(md)][1]
donor_col <- donor_candidates[donor_candidates %in% colnames(md)][1]
if (is.na(group_col)) stop("❌ meta.data 找不到组别列。候选：", paste(group_candidates, collapse=", "))
if (is.na(donor_col)) stop("❌ meta.data 找不到 donor 列。候选：", paste(donor_candidates, collapse=", "))

sink(log_fp, append = TRUE)
cat("✅ group_col = ", group_col, "\n", sep="")
cat("✅ donor_col = ", donor_col, "\n\n", sep="")
sink()

## 统一字段名（贴近你模版）
obj$group2 <- as.character(md[[group_col]])
obj$autopsy_id <- as.character(md[[donor_col]])

## 去掉 NA
keep <- !is.na(obj$group2) & !is.na(obj$autopsy_id) & obj$autopsy_id != ""
obj <- subset(obj, cells = colnames(obj)[keep])

## 统一 Control/NC/CTRL 等名称为 Control
ctrl_alias <- c("Control","CTRL","Ctr","CTR","NC","Normal","N","control","ctrl","ctr","nc","normal")
obj$group2 <- ifelse(obj$group2 %in% ctrl_alias, "Control", as.character(obj$group2))

## 只保留 AD / Control
obj$group2 <- ifelse(obj$group2 %in% c("AD","Control"), obj$group2, NA_character_)
obj <- subset(obj, subset = !is.na(group2))

## 顺序：
## - 统计：Control 作为 reference（和模版一致）
## - 作图：AD 左、Control 右（和你展示的图一致）
group_order_plot  <- c("AD","Control")
group_order_stats <- c("Control","AD")
obj$group2 <- factor(as.character(obj$group2), levels = group_order_stats)

## donor 数核对（补充材料）
## donor 数核对（修复：group2/autopsy_id 可能是 list 列）
md2 <- obj@meta.data
obj$group2     <- as.character(unlist(obj$group2, use.names = FALSE))
obj$autopsy_id <- as.character(unlist(obj$autopsy_id, use.names = FALSE))


donor_tab <- md2 %>%
  dplyr::select(autopsy_id, group2) %>%
  dplyr::filter(!is.na(autopsy_id), autopsy_id != "", !is.na(group2), group2 != "") %>%
  dplyr::distinct() %>%
  dplyr::count(group2, name = "n") %>%
  tidyr::complete(group2 = factor(group_order_stats, levels = group_order_stats),
                  fill = list(n = 0))

print(donor_tab)


## ========= 3) Seurat v5/Assay5：从多 layers 的 counts.* 拼回每细胞 CP10k / log1p =========
rna <- obj[["RNA"]]

## 注意：不同 Seurat 版本可能不是 Assay5；但 counts layers 逻辑仍可用
all_layers <- tryCatch(SeuratObject::Layers(rna), error = function(e) character(0))
count_layers <- grep("^counts(\\.|$)", all_layers, value = TRUE)   # 兼容 counts 和 counts.*

if (length(count_layers) == 0) {
  stop("❌ 没找到任何 counts 或 counts.* layers：无法按 counts->CP10k 流程计算。")
}

write.csv(data.frame(count_layers = count_layers),
          file.path(out_dir, "CHECK_counts_layers.csv"),
          row.names = FALSE)

## 3.1 确定 UBL3 的行名（SYMBOL 或 ENSG）
m0 <- SeuratObject::LayerData(rna, layer = count_layers[1])

target <- if ("UBL3" %in% rownames(m0)) "UBL3" else
  if ("ENSG00000122042" %in% rownames(m0)) "ENSG00000122042" else NULL
stopifnot(!is.null(target))

write.csv(data.frame(gene = "UBL3", gene_row = target),
          file.path(out_dir, "CHECK_geneSymbol_to_rowname.csv"),
          row.names = FALSE)

## 3.2 初始化“全长向量”（以 obj 当前细胞顺序为准）
cells_all <- colnames(obj)
lib_size <- setNames(rep(NA_real_, length(cells_all)), cells_all)
ubl3_counts <- setNames(rep(NA_real_, length(cells_all)), cells_all)

## 3.3 遍历每个 counts layer，填充对应细胞的 lib_size 和 ubl3_counts
for (ly in count_layers) {
  m <- SeuratObject::LayerData(rna, layer = ly)   # dgCMatrix（genes x cells-of-layer）
  cn <- colnames(m)
  if (length(cn) == 0) next
  
  ## 每个细胞总 counts
  lib <- Matrix::colSums(m)
  ## 每个细胞 UBL3 counts
  uct <- as.numeric(m[target, , drop = TRUE])
  
  ## 按 cell name 对齐写入
  lib_size[cn] <- as.numeric(lib)
  ubl3_counts[cn] <- uct
}

## 3.4 兜底：没覆盖到的细胞填 0（通常不该发生，但写上防止 NA 传播）
lib_size[is.na(lib_size)] <- 0
ubl3_counts[is.na(ubl3_counts)] <- 0

## 3.5 计算 CP10k 与 log1p(CP10k)
cp10k <- (ubl3_counts / pmax(lib_size, 1)) * 10000
ubl3_log1p_cp10k <- log1p(cp10k)

expr_source <- "counts_layers"  # 记录表达来源（补充材料）

## ========= 4) donor 层面汇总：每个 donor = 1 点 =========
df_cells <- data.frame(
  autopsy_id = as.character(obj$autopsy_id),
  group2     = as.character(obj$group2),
  ubl3_cp10k = as.numeric(cp10k),
  ubl3_log   = as.numeric(ubl3_log1p_cp10k),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(autopsy_id), autopsy_id != "", !is.na(group2))

## donor 汇总（严格按模版：mean_cp10k 用于 log2FC；mean_log 用于检验+作图）
df_donor <- df_cells %>%
  group_by(autopsy_id, group2) %>%
  summarise(
    n_cells = dplyr::n(),
    mean_cp10k = mean(ubl3_cp10k, na.rm = TRUE),
    mean_log1p_cp10k = mean(ubl3_log, na.rm = TRUE),
    .groups = "drop"
  )

## 作图顺序：AD 左、Control 右
df_donor$group_plot <- factor(df_donor$group2, levels = group_order_plot)

write.csv(df_donor,
          file.path(out_dir, "UBL3_wholecell_perDonor_CP10k_mean_ADvsControl.csv"),
          row.names = FALSE)

message("✅ donor points = ", nrow(df_donor))
print(table(df_donor$group2))

sink(log_fp, append = TRUE)
cat("donor points:", nrow(df_donor), "\n")
cat("donor group table:\n"); print(table(df_donor$group2))
cat("\n")
sink()

## ========= 5) 统计：Wilcoxon + log2FC =========
d1 <- df_donor %>% filter(group2 == "AD")
d0 <- df_donor %>% filter(group2 == "Control")

mu1 <- mean(d1$mean_cp10k, na.rm = TRUE)
mu0 <- mean(d0$mean_cp10k, na.rm = TRUE)

## log2FC（用 mean_cp10k 的均值，避免 log1p 值做 fold-change 的不直观）
log2FC <- log2((mu1 + 1e-8) / (mu0 + 1e-8))

## Wilcoxon：对 donor 的 mean_log1p_cp10k 做检验（与你模版一致）
pval <- if (nrow(d1) >= 2 && nrow(d0) >= 2) {
  wilcox.test(d1$mean_log1p_cp10k, d0$mean_log1p_cp10k, exact = FALSE)$p.value
} else NA_real_

stats_1 <- data.frame(
  contrast = "AD_vs_Control",
  n_AD = nrow(d1),
  n_Control = nrow(d0),
  mean_cp10k_AD = mu1,
  mean_cp10k_Control = mu0,
  log2FC = log2FC,
  p_wilcox = pval,
  expr_source = expr_source,
  used_object = basename(obj_fp),
  stringsAsFactors = FALSE
)

write.csv(stats_1,
          file.path(out_dir, "UBL3_wholecell_stats_ADvsControl.csv"),
          row.names = FALSE)

print(stats_1)

fmt_p <- function(p) if (is.na(p)) "NA" else formatC(p, format = "e", digits = 2)
sub1 <- sprintf("AD vs Control:  log2FC=%s, p=%s",
                if (is.na(log2FC)) "NA" else sprintf("%.3f", log2FC),
                fmt_p(pval))

## ========= 6) 画图：颜色/标题/布局与你模版一致 =========
pal2_plot <- c(AD = "#D24B40", Control = "#2C7FB8")

theme_sci <- theme_bw(base_size = 16) +
  theme(
    plot.title.position = "plot",
    plot.title    = element_text(face = "bold", size = 20, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 12, colour = "grey25", margin = margin(b = 8)),
    panel.border  = element_rect(colour = "grey25", fill = NA, linewidth = 0.8),
    panel.grid.major.y = element_line(linewidth = 0.28, linetype = "dashed", colour = "grey88"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 16, b = 8, l = 10)
  )

p1 <- ggplot(df_donor, aes(x = group_plot, y = mean_log1p_cp10k, fill = group_plot)) +
  geom_boxplot(
    width = 0.55, outlier.shape = NA, linewidth = 1.0,
    alpha = 0.96, colour = "grey15", median.linewidth = 1.6
  ) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0),
    size = 2.6, alpha = 0.9, shape = 21, stroke = 0.5, colour = "grey10"
  ) +
  scale_fill_manual(values = pal2_plot, drop = FALSE) +
  labs(
    title = "UBL3 expression per donor (Whole cells)",
    subtitle = sub1,
    x = NULL,
    y = "Mean UBL3 log1p(CP10k)"
  ) +
  theme_sci +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10)))

## ========= 7) 保存 PNG =========
fig1_fp <- file.path(out_dir, "Fig_NO4_UBL3_wholecell_perDonor_log1pCP10k_boxplot_ADvsControl.png")

ragg::agg_png(fig1_fp, width = 7.8, height = 5.6, units = "in", res = 450, background = "white")
print(p1)
dev.off()

message("✅ saved: ", fig1_fp)

sink(log_fp, append = TRUE)
cat("\n==== sessionInfo ====\n")
print(sessionInfo())
cat("==== END ====\n")
sink()

message("DONE. out_dir = ", out_dir)
