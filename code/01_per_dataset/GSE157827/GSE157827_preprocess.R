###############################################################################
# 01_per_dataset / GSE157827  (Alzheimer's disease, middle frontal gyrus)
# From-raw snRNA-seq preprocessing: builds the harmonized 6-cell-type Seurat
# object consumed by code/02_integrated_figures/.
#
# Source     : GEO accession GSE157827 (raw 10x triplets: matrix.mtx /
#              features.tsv / barcodes.tsv per GSM sample).
# Pipeline   : merge all samples -> QC (nFeature_RNA > 200, nCount_RNA < 20000,
#              percent.mt < 20) -> per-sample Normalize + HVG(1000) -> CCA
#              integration (FindIntegrationAnchors dims = 1:20, IntegrateData
#              dims = 1:20, nfeatures = 2000) -> PCA -> clustering (resolution 1)
#              + UMAP (dims = 1:20) -> cluster markers (Wilcoxon, logFC >= 0.25,
#              adj.P < 0.1) -> 7-class annotation collapsed to 6 cell types.
# Produces   : stepH_obj_celltype6_named.rds  (the object read by the integrated
#              figure/stat scripts; also underpins Supplementary Figure S2).
# Paths      : all paths point to the data-project location on disk. Raw inputs
#              are public at GSE157827; intermediate .rds objects live with the
#              data, not in this repository. Confirm paths before running.
# Environment: R 4.5.1 + Seurat (v5) + SeuratObject + Matrix + data.table.
#
# NOTE: kept faithful to the script as run. Sections after the celltype6
#       hand-off (per-dataset cell-type / UBL3 UMAPs) overlap Figure 2 /
#       Supplementary Figure S3 and are retained only for provenance.
###############################################################################

#GSE157827第 0 章：基础准备
#第 1 章：读入所有样本的 10x 三件套 → 做成一个大 Seurat 对象
#--先手动（readMM）合并1个小样本GSM4775561_AD1的3件套试下（[64-bit] C:\Program Files\R\R-4.5.1）----

# 设置工作路径
setwd("D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/Raw_data")

# 加载包
library(Matrix)
library(Seurat)

# 读取矩阵文件
mat <- readMM("GSM4775561_AD1_matrix.mtx.gz")

# 读取特征和细胞文件
features <- read.delim("GSM4775561_AD1_features.tsv.gz", header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.delim("GSM4775561_AD1_barcodes.tsv.gz", header = FALSE, stringsAsFactors = FALSE)

#检查维度是否对齐
dim(mat)
nrow(features)
nrow(barcodes)

# 添加行列名（这里用 ENSG ID）
rownames(mat) <- features$V1  # 第一列是 ENSG ID
colnames(mat) <- barcodes$V1

# 检查重复或缺失
sum(duplicated(features$V1))   # 重复的 ENSG 数量（应=0）
sum(duplicated(features$V2))   # 重复的 symbol 数量（会>0，正常）
sum(is.na(features$V2) | features$V2 == "")   # symbol 是否有空值

#随手看几个值确认
head(rownames(mat))
head(features$V2)
head(colnames(mat))

# 保存为 R 对象
saveRDS(mat, file = "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results/GSM4775561_AD1_counts_ENSG.rds")







# ---再把GEO157827中所有的样本的3件套合并在一起成一个rds文件----


# 第1步设定输入/输出目录 + 日志
# install.packages(c("Seurat","Matrix","data.table","stringr"))
# 固定随机种子（UMAP/聚类等用得到；贯穿全流程保证复现）
SEED <- 20251023; set.seed(SEED)

# 输入/输出路径
raw_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/Raw_data"
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# 必要包：Matrix 读稀疏矩阵；data.table 快速读文件；Seurat 用于合并
library(Matrix); library(data.table); library(Seurat)

# 记录一次会话信息，保证完全复现（版本差异一眼可查）
sink(file.path(res_dir, "merge_sessionInfo.txt"))
cat("SEED =", SEED, "\n"); print(sessionInfo()); sink()



#第2步找到所有样本前缀（自动）
mtx_files <- list.files(raw_dir, pattern = "_matrix\\.mtx(\\.gz)?$", full.names = TRUE)
prefixes  <- sort(sub("_matrix\\.mtx(\\.gz)?$", "", basename(mtx_files)))
prefixes




#第3步定义一个超短函数：读一个样本（行名=ENSG）
read_one <- function(pfx){
  # 路径
  f_mtx <- file.path(raw_dir, paste0(pfx, "_matrix.mtx.gz"))
  f_fea <- file.path(raw_dir, paste0(pfx, "_features.tsv.gz"))
  f_bar <- file.path(raw_dir, paste0(pfx, "_barcodes.tsv.gz"))
  
  # 读三件套
  mat <- as(readMM(f_mtx), "dgCMatrix")                          # 稀疏计数矩阵
  fea <- fread(f_fea, header = FALSE)                             # V1=ENSG, V2=symbol
  bar <- fread(f_bar, header = FALSE)                             # 条形码
  
  # 对齐维度 + 命名（行=ENSG；列=barcode）
  stopifnot(nrow(fea) == nrow(mat), nrow(bar) == ncol(mat))
  rownames(mat) <- fea$V1                                         # 行名=ENSG（唯一稳定）
  colnames(mat) <- bar$V1
  
  # 建最原始 Seurat 对象（不做任何过滤/归一化，忠实保留原始计数）
  so <- CreateSeuratObject(counts = mat, project = "GSE157827",
                           min.cells = 0, min.features = 0)
  so$gsm    <- sub("_.*$", "", pfx)                               # 如 GSM4775561
  so$sample <- sub("^.*?_", "", pfx)                              # 如 AD1
  colnames(so) <- paste0(so$sample, "_", colnames(so))            # 细胞名唯一：sample_barcode
  so
}





#第4步读所有样本 → 合并 → 自检 → 保存（全在这几行）
# 读所有样本（每读一个会打印维度，便于你看进度）
# 读所有样本（过程会打印每个样本名，方便你看进度）
objs <- lapply(prefixes, function(pfx){ message("Reading ", pfx); read_one(pfx) })
names(objs) <- prefixes

# 合并（只有一个样本就直接取它）
merged <- if (length(objs) == 1) objs[[1]] else Reduce(function(a, b) merge(a, y = b), objs)

# 自检：总基因数×细胞数 + 各样本细胞数（检查是否合理）
cat("Genes x Cells:", nrow(merged), "x", ncol(merged), "\n")
print(table(merged$sample))

# 保存合并后的“原始计数大对象”（后续分析都基于它）
saveRDS(merged, file.path(res_dir, "GSE157827_merged_ENSG_raw_seurat.rds"))



# 各样本细胞数汇总
cell_summary <- as.data.frame(table(merged$sample))
write.csv(cell_summary, file.path(res_dir, "sample_cell_counts.csv"), row.names = FALSE)

# 检查是否存在重复的细胞名（应该为 0）
sum(duplicated(colnames(merged)))

# 检查是否存在重复的基因名（我们用 ENSG，所以也应该为 0）
sum(duplicated(rownames(merged)))

#结束






#第 2 章：合并 GEO metadata（gsm / sample / group）
#-----再把所有样本的rds文件和其medata文件合在一起----
# ---- 路径设置 ----
res_dir  <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"
meta_csv <- file.path(res_dir, "GSE157827_sample_metadata_from_series_matrix.csv")
merged_rds <- file.path(res_dir, "GSE157827_merged_ENSG_raw_seurat.rds")

# ---- 加载包 ----
library(Seurat)
library(data.table)

# ---- 读取对象与meta ----
obj  <- readRDS(merged_rds)
meta <- fread(meta_csv)

# ---- 保留关键列并统一命名 ----
meta_use <- meta[, .(gsm = GSM, sample = sample_suggest, group = disease_group)]

# 检查一下
print(meta_use)

# ---- 合并到 Seurat 对象的 meta.data ----
obj@meta.data <- merge(obj@meta.data, meta_use, by = "gsm", all.x = TRUE)

# ---- 验证合并是否成功 ----
cat("基因数 × 细胞数：", nrow(obj), "×", ncol(obj), "\n")
cat("每组细胞数量：\n")
print(table(obj$group))

# ---- 保存增强版对象 ----
saveRDS(obj, file.path(res_dir, "GSE157827_merged_with_group.rds"))
cat("✅ 已保存：GSE157827_merged_with_group.rds\n")


#这会已经成功得到一个含所有样本和分组信息的rds大文件





#第 3 章：重新整理 counts（单层）+ QC（两步）
#----------------按照原文复现单细胞分析流程（合并好的大矩阵 → Seuratv5→ Cell Ranger 初筛，所以从这步开始做2次QC→每样本 Normalize + HVG(1000) → 样本整合（Anchors dims=1:20 → IntegrateData dims=1:20） → 
→聚类（resolution=1）+ UMAP(dims=1:20) → 找每个簇的 marker（logfc≥0.25；Wilcoxon；adjP<0.1）------------



#第1步 固定随机种子 + 加载环境（可复现的基础设置）
# —— 固定随机种子（整篇分析都用这个）——
SEED <- 20251023; set.seed(SEED)

res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"
suppressPackageStartupMessages({library(Seurat); library(data.table)})

# 读入你合并好的大对象
obj <- readRDS(file.path(res_dir, "GSE157827_merged_with_group.rds"))
DefaultAssay(obj) <- "RNA"

# 1) 规范出唯一的样本列：sample
if (!"sample" %in% colnames(obj@meta.data)) {
  obj$sample <- if ("sample.x" %in% colnames(obj@meta.data)) obj$sample.x else obj$sample.y
}
# 2) 分组列：group（已经有就保留）
stopifnot("group" %in% colnames(obj@meta.data))

# 3) 清理不再需要的列（避免 .x/.y 继续干扰）
for (cl in c("sample.x","sample.y","orig.ident")) {
  if (cl %in% colnames(obj@meta.data)) obj@meta.data[[cl]] <- NULL
}

# 4) 快速自检 + 导出统计（保证“固定”成功）
cat("Genes x Cells:", nrow(obj), "x", ncol(obj), "\n")
print(head(obj@meta.data[, c("gsm","sample","group")], 3))
fwrite(as.data.table(table(obj$sample)), file.path(res_dir, "step0_cells_by_sample_beforeQC.csv"))
fwrite(as.data.table(table(obj$group)),  file.path(res_dir, "step0_cells_by_group_beforeQC.csv"))

# 5) 保存“已固定meta”的对象
saveRDS(obj, file.path(res_dir, "step0_meta_fixed.rds"))
cat("✅ meta 已规范并保存：step0_meta_fixed.rds\n")



#第2步｜准备/补齐 QC 指标并画图（不过滤，仅查看）
# 固定随机种子，保证可复现
SEED <- 20251023; set.seed(SEED)

library(Matrix)
library(Seurat)
library(SeuratObject)

res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"

# 读你现有的大对象（就是“179392 细胞”的那个）
obj_old <- readRDS(file.path(res_dir, "GSE157827_merged_with_group.rds"))
DefaultAssay(obj_old) <- "RNA"
a <- obj_old[["RNA"]]

# 1) 取出每个 layer 的矩阵与条形码
lyr <- Layers(a)
message("共发现 ", length(lyr), " 个 layer.")
mats <- list()
barcodes_all <- character(0)

for (L in lyr) {
  m <- LayerData(a, layer = L)         # 稀疏矩阵（行=基因 ENSG，列=该样本细胞）
  if (ncol(m) == 0) next
  mats[[L]] <- m
  barcodes_all <- c(barcodes_all, colnames(m))
  cat(sprintf("  layer: %-20s | genes=%d cells=%d | 示例条形码=%s\n",
              substr(L,1,20), nrow(m), ncol(m), colnames(m)[1]))
}

# 2) 自检：所有 layer 的基因顺序必须完全一致（这在 cellranger 产物上通常成立）
genes1 <- rownames(mats[[1]])
stopifnot(all(vapply(mats, function(mm) identical(rownames(mm), genes1), logical(1))))

# 3) 合并为一个大 counts（按列拼接）
counts_combined <- do.call(cbind, mats)
stopifnot(ncol(counts_combined) == length(barcodes_all))
stopifnot(identical(colnames(counts_combined), barcodes_all))

cat("合并完成：", nrow(counts_combined), "genes ×", ncol(counts_combined), "cells\n")

# 4) 用合并好的 counts 重建一个“干净”的 Seurat 对象（单 assay / 单 layer）
obj <- CreateSeuratObject(
  counts = counts_combined,
  project = "GSE157827",
  min.cells = 0, min.features = 0
)
# 现在 obj 的细胞名就是 "AD1_AAAC..." / "NC15_TTT..." 等真实条形码
head(colnames(obj))




# 从条形码拆出 sample 与 group
# 条形码形如 "AD1_AAACCCAAG..."；样本名在 "_" 之前
samples <- sub("_.*$", "", colnames(obj))           # "AD1","AD2","NC14", ...
groups  <- ifelse(grepl("^AD", samples), "AD", "Control")

obj$sample <- samples
obj$group  <- groups

# 快速核对
table(obj$sample)[1:10]
table(obj$group)

# 保存一步，便于回滚
saveRDS(obj, file.path(res_dir, "stepA_rebuilt_single_counts.rds"))



# 人类线粒体基因（Ensembl ID；和你的行名 ENSG 对上）
mt_genes_ensg <- c(
  "ENSG00000198888","ENSG00000198727","ENSG00000198804","ENSG00000198886",
  "ENSG00000212907","ENSG00000198786","ENSG00000198695","ENSG00000198712",
  "ENSG00000198899","ENSG00000198938","ENSG00000198840","ENSG00000198763",
  "ENSG00000210107","ENSG00000210112","ENSG00000210117","ENSG00000210127",
  "ENSG00000210133","ENSG00000210140","ENSG00000210144","ENSG00000210151",
  "ENSG00000210156","ENSG00000210160","ENSG00000210164","ENSG00000210169",
  "ENSG00000210174","ENSG00000210179","ENSG00000210184","ENSG00000210189",
  "ENSG00000210194","ENSG00000210199","ENSG00000210204","ENSG00000210209",
  "ENSG00000210214","ENSG00000210219","ENSG00000228253","ENSG00000228630",
  "ENSG00000210130"
)

# 1) 直接基于 obj@assays$RNA@counts 计算（现在就是单层 counts）
mat <- GetAssayData(obj, assay = "RNA", layer = "counts")  # v5 写法
# 若你的 SeuratObject 版本提示 layer 参数不支持，可用：
# mat <- obj@assays$RNA@layers[["counts"]]

# 2) 基础 QC 指标
obj[["nFeature_RNA"]] <- Matrix::colSums(mat > 0)
obj[["nCount_RNA"]]   <- Matrix::colSums(mat)

is_mt <- rownames(mat) %in% mt_genes_ensg
mt_counts    <- Matrix::colSums(mat[is_mt, , drop = FALSE])
total_counts <- Matrix::colSums(mat)
obj[["percent.mt"]] <- ifelse(total_counts > 0, 100 * mt_counts / total_counts, 0)

summary(obj$percent.mt)

# 3) 小提琴图（过滤前）
png(file.path(res_dir, "stepB_QC_violin_before_filter.png"), width=1600, height=600)
print(VlnPlot(obj, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
              ncol = 3, pt.size = 0, raster = TRUE))
dev.off()

# 记录基线规模
writeLines(sprintf("Before filter: genes=%d, cells=%d",
                   nrow(obj), ncol(obj)),
           file.path(res_dir, "stepB_QC_sizes.txt"))




#-------用原文阈值过滤（完全按论文2次QC），并保存“过滤后”对象------
# 论文阈值（复现，不改动）
keep <- (obj$nFeature_RNA > 200) &
  (obj$nCount_RNA   < 20000) &
  (obj$percent.mt   < 20)

cat("过滤后保留细胞数：", sum(keep), " / ", ncol(obj), "\n")

obj_flt <- subset(obj, cells = colnames(obj)[keep])

# 小提琴图（过滤后）
png(file.path(res_dir, "stepC_QC_violin_after_filter.png"), width=1600, height=600)
print(VlnPlot(obj_flt, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
              ncol = 3, pt.size = 0, raster = TRUE))
dev.off()

# 保存过滤后对象
saveRDS(obj_flt, file.path(res_dir, "stepC_filtered_obj.rds"))

# 导出按样本与分组的细胞数，方便对照文章
write.csv(as.data.frame(table(obj$sample, useNA="ifany")),
          file.path(res_dir, "stepC_cells_by_sample_before.csv"), row.names = FALSE)
write.csv(as.data.frame(table(obj_flt$sample, useNA="ifany")),
          file.path(res_dir, "stepC_cells_by_sample_after.csv"),  row.names = FALSE)

write.csv(as.data.frame(table(obj$group, useNA="ifany")),
          file.path(res_dir, "stepC_cells_by_group_before.csv"), row.names = FALSE)
write.csv(as.data.frame(table(obj_flt$group, useNA="ifany")),
          file.path(res_dir, "stepC_cells_by_group_after.csv"),  row.names = FALSE)

# 基因数对照（过滤后，通常会因“全零行”被自动去除而下降）
writeLines(sprintf("After filter: genes=%d, cells=%d",
                   nrow(obj_flt), ncol(obj_flt)),
           file.path(res_dir, "stepC_QC_sizes.txt"))











#上面步骤的质控小提琴图打不开了，所以：
#重新做标准的质控前后的小提琴图2张，按照GSE174367的模版图片格式
###############################################################################
###############################################################################
# GSE157827：按 GSE174367 模版风格输出 QC 前/后两张小提琴图（可直接运行）
# 输出目录：D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results
###############################################################################
# GSE157827：按 GSE174367 模版风格输出 QC 前/后两张小提琴图（稳定版，可直接运行）
# 输入：step0_meta_fixed.rds（优先）或 GSE157827_merged_with_group.rds
# 输出：stepB_QC_violin_before_filter.png / stepC_QC_violin_after_filter.png
###############################################################################
# GSE157827：按 GSE174367 模版风格输出 QC 前/后两张小提琴图（可直接运行）
# 输出目录：D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results
# 特点：
#   - step0_meta_fixed.rds 若有损坏/读取异常，会自动回退到 merged_with_group.rds
#   - 从多 layer 合并重建单层 counts → 画 QC before/after 两张图（模版风格）
###############################################################################

SEED <- 20251023
set.seed(SEED)

res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

# ---------- 安全读 RDS：失败返回 error ----------
safe_read <- function(p){
  tryCatch(readRDS(p), error = function(e) e)
}

# ---------- 0) 选择一个能读的输入对象 ----------
rds_step0  <- file.path(res_dir, "step0_meta_fixed.rds")
rds_merged <- file.path(res_dir, "GSE157827_merged_with_group.rds")

obj_old <- NULL

if (file.exists(rds_step0)) {
  tmp <- safe_read(rds_step0)
  if (!inherits(tmp, "error")) {
    obj_old <- tmp
    message("✅ Loaded: step0_meta_fixed.rds")
  } else {
    message("⚠️ step0_meta_fixed.rds 读入失败：", conditionMessage(tmp))
  }
}

if (is.null(obj_old)) {
  if (!file.exists(rds_merged)) stop("step0 失败且 merged_with_group.rds 不存在。")
  tmp <- safe_read(rds_merged)
  if (inherits(tmp, "error")) stop("merged_with_group.rds 也读入失败：", conditionMessage(tmp))
  obj_old <- tmp
  message("✅ Loaded: GSE157827_merged_with_group.rds")
}

DefaultAssay(obj_old) <- "RNA"
a <- obj_old[["RNA"]]

# ---------- 1) 多 layer → 合并成单层 counts ----------
lyr <- Layers(a)
if (length(lyr) <= 0) stop("未发现 RNA layers，无法重建 counts。")

mats <- lapply(lyr, function(L) LayerData(a, layer = L))
mats <- mats[sapply(mats, function(m) ncol(m) > 0)]
if (length(mats) == 0) stop("layers 全部为空（0 cells），无法重建 counts。")

genes1 <- rownames(mats[[1]])
ok <- all(vapply(mats, function(mm) identical(rownames(mm), genes1), logical(1)))
if (!ok) stop("不同 layer 的基因顺序不一致，不能直接 cbind。")

counts_combined <- do.call(cbind, mats)
cat("Rebuilt counts:", nrow(counts_combined), "genes x", ncol(counts_combined), "cells\n")

obj <- CreateSeuratObject(
  counts = counts_combined,
  project = "GSE157827",
  min.cells = 0, min.features = 0
)

# 保存一份兼容版（避免下次再遇到读不进）
saveRDS(obj, file.path(res_dir, "stepA_rebuilt_single_counts_v2.rds"), version = 2)
message("✅ Saved: stepA_rebuilt_single_counts_v2.rds")

# ---------- 2) 统一 Identity（X轴只显示一个“GSE157827”，和模版一致） ----------
obj$orig.ident <- "GSE157827"
Idents(obj) <- "orig.ident"
DefaultAssay(obj) <- "RNA"

# ---------- 3) 计算 QC 指标 ----------
mat <- GetAssayData(obj, assay = "RNA", layer = "counts")

mt_genes_ensg <- c(
  "ENSG00000198888","ENSG00000198727","ENSG00000198804","ENSG00000198886",
  "ENSG00000212907","ENSG00000198786","ENSG00000198695","ENSG00000198712",
  "ENSG00000198899","ENSG00000198938","ENSG00000198840","ENSG00000198763",
  "ENSG00000210107","ENSG00000210112","ENSG00000210117","ENSG00000210127",
  "ENSG00000210133","ENSG00000210140","ENSG00000210144","ENSG00000210151",
  "ENSG00000210156","ENSG00000210160","ENSG00000210164","ENSG00000210169",
  "ENSG00000210174","ENSG00000210179","ENSG00000210184","ENSG00000210189",
  "ENSG00000210194","ENSG00000210199","ENSG00000210204","ENSG00000210209",
  "ENSG00000210214","ENSG00000210219","ENSG00000228253","ENSG00000228630",
  "ENSG00000210130"
)

obj$nFeature_RNA <- Matrix::colSums(mat > 0)
obj$nCount_RNA   <- Matrix::colSums(mat)

is_mt <- rownames(mat) %in% mt_genes_ensg
obj$percent.mt <- Matrix::colSums(mat[is_mt, , drop = FALSE]) /
  Matrix::colSums(mat) * 100

# ---------- 4) 过滤前 QC 小提琴图（stepB） ----------
png(file.path(res_dir, "stepB_QC_violin_before_filter.png"),
    width = 1600, height = 600)
print(
  VlnPlot(obj,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          ncol = 3, pt.size = 0, raster = TRUE)
)
dev.off()

writeLines(sprintf("Before filter: genes=%d, cells=%d", nrow(obj), ncol(obj)),
           file.path(res_dir, "stepB_QC_sizes.txt"))
message("✅ Saved: stepB_QC_violin_before_filter.png")

# ---------- 5) 按论文阈值过滤 + 过滤后 QC 小提琴图（stepC） ----------
keep <- (obj$nFeature_RNA > 200) &
  (obj$nCount_RNA   < 20000) &
  (obj$percent.mt   < 20)

cat("Cells kept:", sum(keep), "/", ncol(obj), "\n")

obj_flt <- subset(obj, cells = colnames(obj)[keep])
Idents(obj_flt) <- "orig.ident"

png(file.path(res_dir, "stepC_QC_violin_after_filter.png"),
    width = 1600, height = 600)
print(
  VlnPlot(obj_flt,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          ncol = 3, pt.size = 0, raster = TRUE)
)
dev.off()

saveRDS(obj_flt, file.path(res_dir, "stepC_filtered_obj_from_v2.rds"), version = 2)

writeLines(sprintf("After filter: genes=%d, cells=%d", nrow(obj_flt), ncol(obj_flt)),
           file.path(res_dir, "stepC_QC_sizes.txt"))
message("✅ Saved: stepC_QC_violin_after_filter.png")

message("✅ Done. All outputs in: ", res_dir)
























#第 4 章：按样本 Normalize + HVG → 三组 CCA 整合 → 全部整合
#✔️ 1. 每个样本 Normalize（LogNormalize）
#✔️ 2. 每个样本找高变基因（1000）
#✔️ 3. 样本分三组（AD1、AD2、NC）
#✔️ 4. 每组内部做一次 CCA integration
#✔️ 5. 三组再一起做一次 CCA integration
#✔️ 6. 得到最终 integrated assay 用于 UMAP / 聚类
#-------每样本 Normalize + HVG(1000) → 样本整合（Anchors dims=1:20 → IntegrateData dims=1:20）→聚类（resolution=1）+ UMAP(dims=1:20) → 找每个簇的 marker-----
## ========= 统一设置 =========
#将所有样本分3组再整合+Normalize +HVF → CCA 锚点 → 组内整合 → 保存
SEED <- 20251023; set.seed(SEED)
rm(list = ls())   # 清除环境中所有变量
gc()              # 清理内存垃圾
# ——清除可能被锁的符号——
if (exists("method", envir = .GlobalEnv, inherits = FALSE)) rm(method, envir = .GlobalEnv)
if (exists("FUN",    envir = .GlobalEnv, inherits = FALSE)) rm(FUN,    envir = .GlobalEnv)
if (exists("i",      envir = .GlobalEnv, inherits = FALSE)) rm(i,      envir = .GlobalEnv)

# ==== 只加载必要包，减少冲突 ====
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

# ==== 路径 ====
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"

# ==== 从“过滤后的大对象”新鲜开始 ====
obj_flt <- readRDS(file.path(res_dir, "stepC_filtered_obj.rds"))
DefaultAssay(obj_flt) <- "RNA"
# === 按 sample 拆分成 21 个样本 ===
obj_list <- SplitObject(obj_flt, split.by = "sample")





# === 对每个样本单独进行标准化 & 识别高变基因 ===
for (nm in names(obj_list)) {
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
saveRDS(obj_list, file.path(res_dir, "stepD_split_normalized_vst1000_clean.rds"))
cat("✅ 已保存: stepD_split_normalized_vst1000_clean.rds\n")



#分组整合（CCA，dims=1:20，nfeatures=2000）
 #读入拆分后的对象列表
SEED <- 20251023; set.seed(SEED) # 固定随机种子（保证你我多次运行一致）
suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(Matrix)})

res_dir  <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"
## ========= 读入上一步“拆分并标准化”的对象列表 =========
## 这个文件是你刚刚保存好的：每个样本已经 LogNormalize + HVF(vst, n=1000)
obj_list <- readRDS(file.path(res_dir, "stepD_split_normalized_vst1000_clean.rds"))

# 三组划分（便于稳健 CCA）
g_ad1 <- c("AD1","AD2","AD4","AD5","AD6","AD8","AD9")  # AD 第一组
g_ad2 <- c("AD10","AD13","AD19","AD20","AD21")          # AD 第二组
g_nc  <- c("NC3","NC7","NC11","NC12","NC14","NC15","NC16","NC17","NC18")  # 对照组



#写一个小函数（按照原文参数做 CCA 整合 + 保存）
integrate_by_cca <- function(xlist, outfile, seed=SEED) {
  set.seed(seed) # 固定随机性（影响锚点匹配/UMAP初始等）
  # 1) 选择整合特征（原文：nfeatures=2000）
  feats <- SelectIntegrationFeatures(object.list = xlist, nfeatures = 2000)
  set.seed(seed)  # 再固定一次
  # 2) 找锚点（原文：CCA；dims=1:20）
  anchors <- FindIntegrationAnchors(
    object.list     = xlist,
    anchor.features = feats,
    dims            = 1:20,
    reduction       = "cca",
    verbose         = TRUE
  )
  set.seed(seed) # 再固定一次
  # 3) 整合（把批次效应消掉，生成 integrated assay）
  obj_int <- IntegrateData(anchorset = anchors, dims = 1:20)  # integrated assay
  # 4) 落盘，便于断点续跑/复现
  saveRDS(obj_int, outfile)
  invisible(obj_int)
}


#分三段整合并保存
obj_ad1_int <- integrate_by_cca(obj_list[g_ad1],
                                file.path(res_dir, "stepE_integrated_AD_part1.rds"))
obj_ad2_int <- integrate_by_cca(obj_list[g_ad2],
                                file.path(res_dir, "stepE_integrated_AD_part2.rds"))
obj_nc_int  <- integrate_by_cca(obj_list[g_nc],
                                file.path(res_dir, "stepE_integrated_NC.rds"))
cat("✅ 三段整合完成并保存\n")









#-------分组整合（原文：Seurat v3 的 CCA，原文参数：dims=1:20，nfeatures=2000）把三段整合对象合并为一个总对象（仍用 CCA、dims=1:20）-------
## 固定随机性 + 只加载必需包
SEED <- 20251023; set.seed(SEED)
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(Matrix) })

res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"

## 读取你刚刚保存的三段整合对象
ad1 <- readRDS(file.path(res_dir, "stepE_integrated_AD_part1.rds"))
ad2 <- readRDS(file.path(res_dir, "stepE_integrated_AD_part2.rds"))
nc  <- readRDS(file.path(res_dir, "stepE_integrated_NC.rds"))

## 按原文参数：先选整合特征（nfeatures=2000）
set.seed(SEED)
feats2 <- SelectIntegrationFeatures(object.list = list(ad1, ad2, nc), nfeatures = 2000)

## 按原文参数：CCA 找锚点（dims=1:20）
set.seed(SEED)
anchors2 <- FindIntegrationAnchors(
  object.list     = list(ad1, ad2, nc),
  anchor.features = feats2,
  dims            = 1:20,
  reduction       = "cca",
  verbose         = TRUE
)

## 按原文参数：整合（生成 integrated assay）
set.seed(SEED)
obj <- IntegrateData(anchorset = anchors2, dims = 1:20)

## 保存总整合对象
saveRDS(obj, file.path(res_dir, "stepE_integrated_ALL_cca.rds"))
cat("✅ 已保存：stepE_integrated_ALL_cca.rds\n")
#上面文件包含有 21 个样本整合后的表达矩阵，assay "integrated"：批次校正后的表达（比较重要），assay "RNA"：原始表达









#第 5 章：降维（PCA）+ JackStraw + 聚类 + UMAP
# 第 5 章：（1）当前情况 ：R: 4.5.1 (ucrt) Seurat: 5.3.0 SeuratObject: 5.1.99.9000（开发版） Matrix: 1.7-4

#
#-------缩放 → PCA(npca=50) → JackStraw 选前20PC → KNN/聚类(res=1) → UMAP(1:20)-------
#Seurat v5 风格（第一套）：
#读取整合对象 + 基本体检
## ——固定随机种子（可复现）——
SEED <- 20251023
set.seed(SEED)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(qs)     # 你已安装
})

res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"

## 优先读 RDS(v2)，失败就读 qs
obj <- tryCatch(
  readRDS(file.path(res_dir, "stepE_integrated_ALL_cca_v2.rds")),
  error = function(e) {
    message("RDS 读取失败，改读 QS：", conditionMessage(e))
    qs::qread(file.path(res_dir, "stepE_integrated_ALL_cca.qs"))
  }
)

## 体检：看一下有哪些 assay
cat("Assays in obj: ", paste(Assays(obj), collapse = ", "), "\n")
## 按 Seurat/论文惯例：integrated 用于降维/聚类/UMAP；RNA 用于找 marker/差异
stopifnot("integrated" %in% Assays(obj), "RNA" %in% Assays(obj))



#PCA（50 PCs）+ JackStraw 评估 + 选前 20 PCs
## ——用 integrated 做降维——
DefaultAssay(obj) <- "integrated"

## 标准化（集成矩阵），PCA=50
set.seed(SEED)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

## JackStraw（可能较慢，但与原文一致）
set.seed(SEED)
obj <- JackStraw(obj, num.replicate = 100, dims = 50)
obj <- ScoreJackStraw(obj, dims = 1:50)

## 导出评估图（保存到文件夹）
png(file.path(res_dir, "stepF_JackStraw_PC1_50.png"), width = 1600, height = 900)
print(JackStrawPlot(obj, dims = 1:50))
dev.off()

png(file.path(res_dir, "stepF_PC_Elbow.png"), width = 1000, height = 800)
print(ElbowPlot(obj, ndims = 50))
dev.off()

## 保存检查点（防止重算）
saveRDS(obj, file.path(res_dir, "stepF_after_PCA_JackStraw_v2.rds"), version = 2, compress = "xz")
qs::qsave(obj, file.path(res_dir, "stepF_after_PCA_JackStraw.qs"))
#两个文件保存的是同一个时间点的 obj（PCA+JackStraw 之后的状态），内容完全一致，只是存储格式不同
cat("✅ PCA(50) + JackStraw 完成并已保存\n")





#图谱构建（Neighbors/Clusters/UMAP）
## 以 1:20 PCs 做邻居图 + 图聚类 + UMAP
set.seed(SEED)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)

## 论文：resolution = 1
obj <- FindClusters(obj, resolution = 1)

## 论文：UMAP 用 dims = 1:20
set.seed(SEED)
obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)

## 快速可视化保存
png(file.path(res_dir, "stepF_UMAP_clusters.png"), width = 1600, height = 900)
print(DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()

## 保存检查点
saveRDS(obj, file.path(res_dir, "stepF_cluster_umap_v2.rds"), version = 2, compress = "xz")
#（Seurat v5 style）
qs::qsave(obj, file.path(res_dir, "stepF_cluster_umap.qs"))
cat("✅ 图谱构建（Neighbors/Clusters/UMAP）完成并已保存\n")









#第 5 章：（2）为了逼近原文 Seurat v3，又跑了一次：
## =========因为上面Seurat V5跑的UMAP和原文不一样，所以用原文的Seurat V3 （只是默认参数设置和V3一样，其实还是V5的版本）=========

## ========= 固定随机性 & 只加载必须包 =========
SEED <- 20251023; set.seed(SEED)
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(Matrix) })

res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"

## 读取整合对象（优先 v2 RDS；没有就读 qs）
obj <- readRDS(file.path(res_dir, "stepE_integrated_ALL_cca_v2.rds"))
# 如果上面这一行报错，可改用：
# library(qs)
# obj <- qread(file.path(res_dir, "stepE_integrated_ALL_cca.qs"))

## 1) 用 integrated assay 做下游分析（论文就是这么做的）
DefaultAssay(obj) <- "integrated"

## 2) 尺度化（v3 风格：center=TRUE, scale=TRUE）
obj <- ScaleData(obj, verbose = FALSE)

## 3) PCA：保留 50 PCs（论文 npcs=50）
set.seed(SEED)
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

## 4) JackStraw：可选但论文画过；这里跑 100 次以省时间（更慢的话你已跑过一次就可以跳过）
set.seed(SEED)
obj <- JackStraw(obj, reduction = "pca", dims = 50, num.replicate = 100, verbose = FALSE)
obj <- ScoreJackStraw(obj, dims = 1:50)
# 保存一张图，方便查 20 PCs 是否合适
png(file.path(res_dir, "stepF_JackStraw_PC_pvalues.png"), width=1200, height=800)
print(JackStrawPlot(obj, dims = 1:50))
dev.off()

## 5) 按论文使用前 20 PCs
## 你已经完成到 RunPCA 了，直接继续：
dims_use <- 1:20                   # 👉 要改就改这里的dims_use数字（比如改成1:30）

## （可选）用肘部图简单确认一下 20 PCs 是合理的
png(file.path(res_dir, "stepF_ElbowPlot_50PCs.png"), width=1200, height=800)
print(ElbowPlot(obj, ndims = 50))
dev.off()

## 构图 → 聚类 → UMAP（已按我们之前的 v3 风格固定好参数）
set.seed(SEED)
obj <- FindNeighbors(obj, dims = dims_use, k.param = 20, verbose = FALSE)

set.seed(SEED)
obj <- FindClusters(obj, resolution = 1, algorithm = 1, n.start = 10, n.iter = 10, verbose = FALSE)
# 👉 要改就改这里的resolution 数字（（比如改成0.6或1.5））
set.seed(SEED)
obj <- RunUMAP(
  obj, reduction = "pca", dims = dims_use,
  umap.method = "uwot", metric = "cosine",
  n.neighbors = 30, min.dist = 0.3, spread = 1,    # 👉 要改就改这里的neighbors数字（比如改成20或40）
  init = "spectral", n.components = 2, verbose = FALSE     # 👉 要改就改这里的min.dist数字（比如改成0.2或0.5）
)

## 保存图 & 对象
png(file.path(res_dir, "stepF_umap_clusters_res1_v3style.png"), width=2000, height=1200)
print(DimPlot(obj, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.2) + NoLegend())
dev.off()

saveRDS(obj, file.path(res_dir, "stepF_afterPCA_graph_umap_res1_v3style.rds"), version = 2, compress = "xz")
#（更接近论文，最终使用做 downstream（第 6、7、8、9 章）




  


#读取对象 + 设置 integrated assay
# 固定随机种子，保证每次结果一致
SEED <- 20251023; set.seed(SEED)

# 加载必须的包（Seurat 是主角）
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

# 设置你的结果文件夹路径
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"

# 读取你上一步保存的整合对象
# （如果你保存的是 43 簇版本就换那行）
obj <- readRDS(file.path(res_dir, "stepF_afterPCA_graph_umap_res1_v3style.rds"))

# 指定后续分析都用“integrated”数据
# （就是整合后去掉批次效应的矩阵）
DefaultAssay(obj) <- "integrated"

# 设置当前身份为“聚类编号”，即每个簇的 ID
Idents(obj) <- "seurat_clusters"










#第 6 章：找 marker + 标注 cell type（6 大类）
#------每簇找 marker（与原文一致：Wilcoxon + logFC≥0.25 + adj.P<0.1）------
# 用 Wilcoxon 秩和检验找每个簇与其它所有簇的差异基因
# 论文参数：logFC 阈值 0.25，检验方法 wilcox
markers_all <- FindAllMarkers(
  obj,                     # 你的 Seurat 对象
  only.pos = FALSE,         # 同时找上调和下调基因
  test.use = "wilcox",      # 检验方法
  logfc.threshold = 0.25,   # 最小 log2(倍数变化)
  min.pct = 0.1,            # 至少在10%的细胞里表达
  verbose = FALSE           # 不输出太多过程信息
)

# 加一列：p 值校正后是否 <0.1（论文阈值）
markers_all$pass_0.1 <- markers_all$p_val_adj < 0.1

# 保存完整结果（含上下调基因），CSV 文件可在 Excel 打开
write.csv(markers_all,
          file.path(res_dir, "stepG_FindAllMarkers_wilcox_logfc0.25.csv"),
          row.names = FALSE)


## 取每簇前20个正向marker，方便画图和命名
# 只保留上调的基因（avg_log2FC > 0 且通过显著性）
top20 <- subset(markers_all, avg_log2FC > 0 & pass_0.1)

# 按簇号和 logFC 从大到小排序
top20 <- top20[order(top20$cluster, -top20$avg_log2FC), ]

# 每个簇取前 20 个上调基因（视觉展示用）
top20 <- do.call(rbind, by(top20, top20$cluster, head, n = 20))

# 画点图（每个簇前20个 marker）
png(file.path(res_dir, "stepG_DotPlot_top20_per_cluster.png"),
    width=2200, height=1400, res=180)
print(DotPlot(obj, features = unique(top20$gene)) + RotatedAxis())
dev.off()

# 画热图（同样这些 marker）
png(file.path(res_dir, "stepG_Heatmap_top20_per_cluster.png"),
    width=2200, height=1400, res=180)
print(DoHeatmap(obj, features = unique(top20$gene), raster = TRUE))
dev.off()









#查看共是否是和原文一样43个簇
length(levels(Idents(obj)))  # 应该返回 43

length(levels(Idents(obj)))
nlevels(Idents(obj))










#------根据 marker 表达手动命名细胞类型-------
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
})

# 1) 取 integrated 的数据层（v5 推荐 LayerData）
mat_en <- LayerData(obj[["integrated"]], layer = "data")  # 行=ENSG，列=cell
ensg <- rownames(mat_en)

# 2) ENSG → SYMBOL 映射
map_df <- AnnotationDbi::select(org.Hs.eg.db, keys = ensg,
                                keytype = "ENSEMBL", columns = "SYMBOL")
sym_by_row <- map_df$SYMBOL[match(ensg, map_df$ENSEMBL)]
keep <- !is.na(sym_by_row)
mat_en <- mat_en[keep, , drop = FALSE]
sym_by_row <- sym_by_row[keep]

# 3) 聚合到“符号行”
sym_levels <- sort(unique(sym_by_row))
row_index  <- match(sym_by_row, sym_levels)
G <- sparseMatrix(i = row_index, j = seq_along(row_index),
                  x = 1, dims = c(length(sym_levels), length(row_index)))
mat_sym <- G %*% mat_en
rownames(mat_sym) <- sym_levels

# 4) 经典 marker（论文里点名的加上）
markers_ref <- list(
  Astro   = c("AQP4","GFAP","ALDH1L1","SLC1A3","ADGRV1","GPC5","RYR3"),
  Endo    = c("CLDN5","KDR","FLT1","PECAM1","ABCB1","EBF1"),
  Excit   = c("CAMK2A","SLC17A7","TBR1","CBLN2","LDB2"),
  Inhib   = c("GAD1","GAD2","SLC6A1","LHFPL3","PCDH15"),
  Microgl = c("C3","CX3CR1","P2RY12","AIF1","DOCK8","LRMDA"),
  Oligo   = c("MBP","MOG","PLP1","MOBP","ST18"),
  # 这组是我之前单列的“周细胞”，论文未单列：等会儿会合并进 Endo
  Peri    = c("PDGFRB","RGS5","MCAM","ACTA2")
)

# 5) 计算“每个簇 × 每个大类”的平均 marker 分数
clu <- Idents(obj); clu_levels <- levels(clu)

score_one_group <- function(genes) {
  g <- intersect(genes, rownames(mat_sym))
  if (length(g) == 0) return(setNames(rep(NA_real_, length(clu_levels)), clu_levels))
  per_cell <- Matrix::colMeans(mat_sym[g, , drop = FALSE])
  tapply(per_cell, INDEX = clu, FUN = mean, na.rm = TRUE)[clu_levels]
}

avg_by_cluster <- sapply(markers_ref, score_one_group)  # 行=簇；列=细胞大类(含Peri)
stopifnot(identical(rownames(avg_by_cluster), clu_levels))

# 6) 先在“7 类”框架下选择每簇的最佳标签（便于审计）
tmp <- avg_by_cluster; tmp[is.na(tmp)] <- -Inf
lab_7 <- colnames(tmp)[max.col(tmp, ties.method = "first")]
names(lab_7) <- rownames(tmp)

# 7) 为了严格对齐原文：把 Peri 合并到 Endo → 得到“6 大类”标签
lab_6 <- lab_7
lab_6[lab_6 == "Peri"] <- "Endo"  # 合并
lab_6 <- factor(lab_6, levels = c("Astro","Endo","Excit","Inhib","Microgl","Oligo"))

# 8) 写回对象并保存两套标签（便于回看）
obj$celltype7 <- plyr::mapvalues(Idents(obj), from = names(lab_7), to = unname(lab_7))
obj$celltype6 <- plyr::mapvalues(Idents(obj), from = names(lab_6), to = unname(lab_6))

saveRDS(obj, file.path(res_dir, "stepH_obj_labeled_celltype7_celltype6.rds"),
        version = 2, compress = "xz")
#UBL3 UMAP，histogram，KS density curves，pseudo-bulk DESeq2都是基于这个对象。

# 9) 导出“每簇的分数和两套标签”CSV（你要的“中间步骤表”）
df_scores <- data.frame(cluster = rownames(avg_by_cluster),
                        label_7 = lab_7,
                        label_6 = lab_6,
                        avg_by_cluster,
                        check.names = FALSE)
#统计表
write.csv(df_scores,
          file.path(res_dir, "stepH_cluster_scores_label7_label6.csv"),
          row.names = FALSE)

# 10) 快速看分布（应为 6 大类）
# 看 6 大类的分布
table(obj$celltype6)

# 看 “簇 × 6 大类” 的交叉频数（每个簇主要对应哪一大类一目了然）
table(Idents(obj), obj$celltype6)

# 只要对象在，随时能复画 6 大类 UMAP
png(file.path(res_dir, "stepH_umap_celltype6.png"), width=2000, height=1400, res=180)
print(DimPlot(obj, reduction = "umap", group.by = "celltype6",
              label = TRUE, label.size = 5))
dev.off()






#------之前的UMAP生成的是只带数字的6种细胞簇，这会把它转成各自名字的UMAP------
## 固定随机数 & 必要包
SEED <- 20251023; set.seed(SEED)
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(plyr)
})

## 路径 & 读取你刚保存过的对象（含 UMAP/聚类）
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"
obj <- readRDS(file.path(res_dir, "stepH_obj_labeled_celltype7_celltype6.rds"))

## ── 1) 把数字代码 → 变成文字标签（并设定想要的顺序）──
## 有些时候 mapvalues 写入 factor 会变成整数编码；这里统一转成字符
## 我们定义 6 大类的人类可读标签
nice_levels  <- c("Astrocytes","Endothelial","Excitatory neurons",
                  "Inhibitory neurons","Microglia","Oligodendrocytes")
short2nice   <- c(Astro="Astrocytes", Endo="Endothelial", Excit="Excitatory neurons",
                  Inhib="Inhibitory neurons", Microgl="Microglia", Oligo="Oligodendrocytes")

## 如果 celltype6 是数字/杂乱，先尽力识别；最稳妥的做法：重新从簇→标签的映射文件重建
csv_map <- file.path(res_dir, "stepH_cluster_scores_label7_label6.csv")
if (file.exists(csv_map)) {
  map_df <- read.csv(csv_map, stringsAsFactors = FALSE)
  ## map_df$label_6 目前是短标签（Astro/Endo/...）；转成人类可读
  lab6_map <- setNames(unname(short2nice[ map_df$label_6 ]), map_df$cluster)  # 名字=簇号
  ## 把“簇身份”映射成 6 大类（字符）
  obj$celltype6 <- plyr::mapvalues(Idents(obj), from = names(lab6_map),
                                   to   = as.character(lab6_map))
} else {
  ## 如果映射表不存在，就把现有的 meta 列尝试转成字符，再套用中英文对照
  val <- as.character(obj$celltype6)
  ## 若里面是短码（Astro/Endo/...），替换成全名；若是全名就会保持不变
  val[val %in% names(short2nice)] <- short2nice[val[val %in% names(short2nice)]]
  obj$celltype6 <- val
}

## 设定有序因子，保证图例顺序与论文一致
obj$celltype6 <- factor(obj$celltype6, levels = nice_levels)

## ── 2) 快速检查 & 导出计数/百分比表 ──
## 每类细胞的细胞数
counts6 <- table(obj$celltype6)
print(counts6)

## 百分比（四舍五入到 1 位小数）
props6 <- round(100 * prop.table(counts6), 1)
print(props6)

## 保存 CSV：细胞数 & 占比
write.csv(data.frame(celltype = names(counts6),
                     n_cells = as.integer(counts6),
                     percent = as.numeric(props6)),
          file.path(res_dir, "stepH_celltype6_counts_percent.csv"),
          row.names = FALSE)

## ── 3) 重新画 6 大类 UMAP，显示文字图例与标签 ──
## 自定义一组稳定且易读的颜色（按 Astro/Endo/Excit/Inhib/Microgl/Oligo 的顺序）
pal6 <- c(Astrocytes="#E76F51", Endothelial="#2A9D8F", `Excitatory neurons`="#457B9D",
          `Inhibitory neurons`="#F4A261", Microglia="#8D99AE", Oligodendrocytes="#6D597A")

png(file.path(res_dir, "stepH_umap_celltype6_named.png"), width=2000, height=1400, res=180)
print(
  DimPlot(obj, reduction = "umap", group.by = "celltype6",
          label = TRUE, label.size = 5, repel = TRUE) +
    scale_color_manual(values = pal6, drop = FALSE) +
    labs(color = "Cell type")
)
dev.off()

## ── 4) 再给你一个“每个样本的 6 类比例”堆叠条形图（对齐原文 Fig S1D）──
if ("sample" %in% colnames(obj@meta.data)) {
  df_bar <- as.data.frame(table(sample = obj$sample, celltype = obj$celltype6))
  df_bar <- within(df_bar, {
    total_by_sample <- ave(Freq, sample, FUN = sum)
    percent <- 100 * Freq / total_by_sample
  })
  png(file.path(res_dir, "stepH_bar_sample_celltype6.png"), width=1800, height=1200, res=180)
  print(
    ggplot(df_bar, aes(x = sample, y = percent, fill = celltype)) +
      geom_bar(stat = "identity", width = 0.9) +
      scale_fill_manual(values = pal6, drop = FALSE) +
      coord_flip() +
      labs(x = NULL, y = "Proportion (%)", fill = "Cell type") +
      theme_bw(base_size = 12)
  )
  dev.off()
}
#UMAP图中点的大小与以下有关：没有显式设 pt.size（→ 使用 Seurat 默认，大约 0.2）；没有设 raster（→ 对十几万细胞会自动 raster，但颜色仍然比较正常）
#没有加 theme_bw，用的是 Seurat 的默认 theme（其实是 ggplot2 的 theme_grey 做了一点修改）
## 最后保存对象（带改好的标签）
saveRDS(obj, file.path(res_dir, "stepH_obj_celltype6_named.rds"), version = 2, compress = "xz")
cat("✅ 已把 6 大类改成文字标签、导出计数/比例，并重画 UMAP/（可选）样本堆叠图。\n")







#第 1 章UMAP
#直接从最终对象"stepH_obj_celltype6_named.rds"画“Cell type UMAP + UBL3 UMAP”（两联图），
#完全不需要从第一行开始重跑。你已经有了“整合 + 聚类 + 细胞类型注释 + UMAP”的最终对象，只要从这个 RDS 继续往下画图就行。

#第一步，先获取 UBL3 的全局真实最大值（用于 colorbar）（画UMAP图右侧的颜色条（Color Bar 或 Color Scale），也称为色阶条。）

## ============================================================
## 0. 基本设置：随机数种子 + 加载 R 包   （保证结果可复现）
## ============================================================







#第7章（2）UMAP ，6个细胞类型和UBL3的，只画UBL3＞0的
#还有验证UMAP的UBL3数据是不是和后面直方图的UBL3表达数据是否一致
###############################################################################
# 第 7 章（2）UMAP：6 个细胞类型 + UBL3（>0 高亮）
# 数据集：GSE157827
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
## 以后如果别的 GSE 里有新的写法，就在这里继续添加映射就行
standardize_celltype6 <- function(x) {
  x <- as.character(x)
  
  # 下面这些右边的字符串，请根据你 stepH 对象里真实的 celltype6 名称再调整 / 添加
  # 可以先运行 sort(unique(obj$celltype6)) 看看有哪些名字
  
  # Astrocytes
  x[x %in% c("Astrocytes", "Astro", "Astrocyte")] <- "Astrocytes"
  
  # Excitatory neurons
  x[x %in% c("Excitatory neurons", "Excitatory", "Ex_neuron", "Excit")] <- "Excitatory neurons"
  
  # Inhibitory neurons
  x[x %in% c("Inhibitory neurons", "Inhibitory", "Inh_neuron", "Inhib")] <- "Inhibitory neurons"
  
  # Microglia
  x[x %in% c("Microglia", "Micro"， "Microgl")] <- "Microglia"
  
  # Endothelial
  x[x %in% c("Endothelial", "Endothelial cells", "Endo")] <- "Endothelial"
  
  # Oligodendrocytes
  x[x %in% c("Oligodendrocytes", "Oligo", "Oligodendrocyte", "Oligodendro")] <- "Oligodendrocytes"
  
  factor(x, levels = celltype6_levels_std)
}

## ============================================================
## 1. 路径与对象：读入 UMAP 对象 & QC 后 counts 对象
## ============================================================
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"

## 1.1 stepH 对象：已经整合 + 聚类 + celltype6 + UMAP
obj <- readRDS(file.path(res_dir, "stepH_obj_celltype6_named.rds"))

## 1.2 stepC 对象：QC 后的 RNA counts（单层），不含 UMAP
obj_expr <- readRDS(file.path(res_dir, "stepC_filtered_obj.rds"))
DefaultAssay(obj_expr) <- "RNA"

## ============================================================
## 2. 从 QC 后 counts 计算 UBL3 的 log1p(CPM)，构建 df_all
## ============================================================

## 2.1 counts 矩阵：行 = 基因 (ENSG)，列 = 细胞
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
if (ubl3_id %in% rownames(rna_counts)) {
  gene_row <- ubl3_id
} else if ("UBL3" %in% rownames(rna_counts)) {
  gene_row <- "UBL3"
} else {
  stop("在 counts 中找不到 UBL3（既无 ENSEMBL 也无 SYMBOL）")
}
cat("实际使用的 UBL3 行名为：", gene_row, "\n")

## 2.4 对齐细胞名：只保留同时出现在 obj_expr 和 obj 里的细胞
cells_expr <- colnames(obj_expr)
cells_umap <- colnames(obj)
common_cells <- intersect(cells_expr, cells_umap)
cat("两个对象共有细胞数：", length(common_cells), "\n")

## 2.5 提取 UBL3 counts & library size
rna_sub  <- rna_counts[gene_row, common_cells, drop = FALSE]  # 1 × N
raw_vec  <- as.numeric(rna_sub[1, ])                          # UBL3 原始 counts
lib_size <- Matrix::colSums(rna_counts[, common_cells, drop = FALSE])

## 2.6 按 Seurat::NormalizeData(LogNormalize) 计算 
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
#如果看到有 <NA> 或者有一些特别的名字没有被映射，就回到最前面的 standardize_celltype6() 里，把这些名字加到相应的一行即可。
###############################################################################



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
  raster     = FALSE,   # ★ 关闭自动栅格化，强制用普通点来画
  pt.size    = 0.1      # ★ 根据需要调小一点点大小（可选）
) +
  scale_color_manual(
    values = celltype6_palette_std,
    breaks = celltype6_levels_std,
    limits = celltype6_levels_std,
    drop   = TRUE
  ) +
  ggtitle("Cell Type UMAP") +
  labs(color = "Cell type") +
  xlab("UMAP_1") + ylab("UMAP_2") +  # 新增：强制轴标签为大写
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
  xlab("UMAP_1") + ylab("UMAP_2") +  # 修改：小写umap_1 → 大写UMAP_1
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

out_file <- file.path(res_dir, "Fig_UBL3_Celltype_and_UMAP_highlight_GSE157827.png")
ggsave(out_file, p_AB_high, width = 12, height = 6, dpi = 300)
cat("🎉 已输出两联图（PNG）：", out_file, "\n")

out_pdf <- file.path(res_dir, "Fig_UBL3_Celltype_and_UMAP_highlight_GSE157827.pdf")
ggsave(out_pdf, p_AB_high, width = 12, height = 6)
cat("🎉 已输出两联图（PDF）：", out_pdf, "\n")














#第2章 ，直方图图

#第 2 章（0）：UBL3 的直方图，单个样本的直方图
#直方图（一般用UMAP确定好细胞类型后的那个rds文件做）
#

## ===========================
## Part 1: UBL3>0 每个细胞类型的直方图（按样本分面）
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

## 按需要改成自己的路径（C:/ 或 D:/ 都可以）
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"

## 1. 读入对象 -------------------------------------------------
obj_meta <- readRDS(file.path(res_dir, "stepH_obj_celltype6_named.rds"))
obj_expr <- readRDS(file.path(res_dir, "stepC_filtered_obj.rds"))
DefaultAssay(obj_expr) <- "RNA"

## 2. 从 counts 里算 UBL3 的 log1p (CP10k) -----------------------
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
UBL3_norm <- log1p((raw_vec / lib_size) * 1e4)  # log1p (CP10k)

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

# facet 用的标签：AD1_AD / NC14_Control 等
df_pos$sample_group <- with(df_pos, paste0(sample, "_", group))
df_pos$sample_group <- factor(df_pos$sample_group,
                              levels = sort(unique(df_pos$sample_group)))

celltypes <- sort(unique(df_pos$celltype6))

# 输出目录
out_dir <- file.path(res_dir, "Fig_UBL3_hist_per_celltype_noZero")
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
      x     = "UBL3 expression (log1p (CP10k))",
      y     = "Cell count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
      strip.text    = element_text(face = "bold", size = 10),
      legend.position = "top",
      legend.title    = element_text(face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = file.path(out_dir,
                         paste0("UBL3_noZero_", sanitize_name(ct), ".png")),
    plot   = p_ct,
    width  = 12,
    height = 8,
    dpi    = 300
  )
}

cat("✅ 已生成 6 张 UBL3>0 直方图：", out_dir, "\n")















## 第2章（1），Y轴是cell，的overlap 直方图。UBL3 计算：仍然是 log1p((raw/lib_size)*1e4)（CP10K），X轴应为
## ============================================================
## 你要做的**第10章（2）**是：
#图的格式：继续严格模仿 GSE157827 模版那张“overlap 直方图（Y轴=cell count）+ facet_wrap(~celltype6, scales=free_y) + 右上角统计标签”的样式
###############################################################################
## GSE157827：Overlap Hist (Cell count) | AD vs Control | donor-level + cell-level
## 风格对齐你刚刚的 density 图（红/蓝、legend n、label 右上角、6分面）
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
## 0) 路径
## =========================
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"
obj_fp  <- file.path(res_dir, "stepH_obj_celltype6_named.rds")
stopifnot(file.exists(obj_fp))

out_dir <- file.path(res_dir, "NOxx_overlap_hist_count_compare_syn520style")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

dataset_tag <- "GSE157827"
gene_symbol <- "UBL3"
binwidth    <- 0.2
disease     <- "AD"

log_fp <- file.path(out_dir, "NOxx_log.txt")
sink(log_fp)
cat("==== START ====\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("obj_fp:", obj_fp, "\n")
cat("out_dir:", out_dir, "\n")
cat("SEED:", SEED, "\n\n")
sink()

## =========================
## 1) 读对象 + 自动识别列
## =========================
obj <- readRDS(obj_fp)
DefaultAssay(obj) <- "RNA"
md <- obj@meta.data

ct_col_candidates <- c("celltype6","celltype","cell_type","celltype6_named","celltype6_manual")
ct_col <- ct_col_candidates[ct_col_candidates %in% colnames(md)][1]
if (is.na(ct_col)) stop("❌ meta.data 找不到 celltype6 列。")

grp_col_candidates <- c("group4","group2","group","Group","diagnosis","Dx","clinical_diagnosis")
grp_col <- grp_col_candidates[grp_col_candidates %in% colnames(md)][1]
if (is.na(grp_col)) stop("❌ meta.data 找不到 group 列。")

donor_col_candidates <- c("autopsy_id","donor","Donor","sample","Sample","orig.ident","patient","subject")
don_col <- donor_col_candidates[donor_col_candidates %in% colnames(md)][1]
if (is.na(don_col)) stop("❌ meta.data 找不到 donor 列。")

grp_raw <- trimws(as.character(md[[grp_col]]))
ctrl_alias <- c("Control","CTRL","Ctr","CTR","NC","Normal","N","control","ctrl","ctr","nc","normal")
grp_std <- ifelse(grp_raw %in% ctrl_alias, "Control", grp_raw)

md$celltype6 <- as.character(md[[ct_col]])
md$group4    <- grp_std
md$donor     <- trimws(as.character(md[[don_col]]))

if (!all(c("AD","Control") %in% unique(md$group4))) {
  sink(log_fp, append=TRUE)
  cat("⚠ group unique:\n"); print(sort(unique(md$group4)))
  sink()
  stop("❌ 分组里没有同时包含 AD 和 Control（或 NC→Control）。")
}

don_all <- unique(md[, c("donor","group4")])
qc_don <- as.data.frame(table(don_all$group4), stringsAsFactors=FALSE)
colnames(qc_don) <- c("group4","n_donors")
write.csv(qc_don, file.path(out_dir, "QC_donors_by_group4.csv"), row.names=FALSE)

sink(log_fp, append=TRUE)
cat("Loaded cells:", ncol(obj), " genes:", nrow(obj), "\n")
cat("celltype col:", ct_col, " | group col:", grp_col, " | donor col:", don_col, "\n")
cat("Donors by group4:\n"); print(table(don_all$group4))
cat("\n")
sink()

## =========================
## 2) 合并 counts layers（Seurat v5 多 layers 兼容）
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
    if (length(mats)==0) stop("❌ counts layers 存在，但读取 layerData 失败。")
    
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
    
    dup <- duplicated(colnames(mat_all))
    if (any(dup)) {
      sink(log_fp, append=TRUE)
      cat("⚠ duplicated cells across layers:", sum(dup), " -> keep first\n")
      sink()
      mat_all <- mat_all[, !dup, drop=FALSE]
    }
    
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
  
  m2 <- tryCatch(Seurat::GetAssayData(obj, assay=assay, slot="counts"), error=function(e) NULL)
  if (!is.null(m2)) return(m2)
  
  stop("❌ 无法获取 counts。")
}

rna_counts <- get_counts_matrix_allcells(obj, "RNA")
stopifnot(ncol(rna_counts) == ncol(obj))

sink(log_fp, append=TRUE)
cat("Merged counts dim:", paste(dim(rna_counts), collapse=" x "), "\n\n")
sink()

## =========================
## 3) 计算 UBL3 log1p(CP10k) + df0(expr>0)
## =========================
gene_row <- if ("UBL3" %in% rownames(rna_counts)) "UBL3" else
  if ("ENSG00000122042" %in% rownames(rna_counts)) "ENSG00000122042" else NA_character_
if (is.na(gene_row)) stop("❌ counts 行名中找不到 UBL3（UBL3 或 ENSG00000122042）")

lib_size <- Matrix::colSums(rna_counts)
expr <- log1p((as.numeric(rna_counts[gene_row, , drop=TRUE]) / pmax(lib_size,1)) * 1e4)

df_all <- data.frame(
  expr      = expr,
  donor     = md$donor,
  group4    = md$group4,
  celltype6 = md$celltype6,
  stringsAsFactors = FALSE
)

df0 <- df_all[df_all$expr > 0, ]  # only expressed cells
saveRDS(df0, file.path(out_dir, "INTERMEDIATE_df0_exprGT0.rds"))
write.csv(df0, file.path(out_dir, "INTERMEDIATE_df0_exprGT0.csv"), row.names=FALSE)

sink(log_fp, append=TRUE)
cat("Expressed cells (expr>0):", nrow(df0), "\n\n")
sink()

## =========================
## 4) 绘图（Y轴=Cell count）+ MWU(BH)
## =========================
plot_one_unit <- function(unit=c("donor","cell")) {
  
  unit <- match.arg(unit)
  
  df2 <- df0 %>%
    filter(group4 %in% c(disease,"Control")) %>%
    mutate(group = ifelse(group4==disease, disease, "Control"))
  
  ## legend：donor n（按 expr>0 的数据口径）
  don_pair <- unique(df2[, c("donor","group")])
  n_dis <- sum(don_pair$group == disease)
  n_ctl <- sum(don_pair$group == "Control")
  
  lab_dis <- sprintf("%s\n(n = %d)", disease,  n_dis)
  lab_ctl <- sprintf("Control\n(n = %d)",     n_ctl)
  
  df2$group_lab <- factor(ifelse(df2$group==disease, lab_dis, lab_ctl),
                          levels=c(lab_dis, lab_ctl))
  
  ## 统计输入
  if (unit=="donor") {
    stat_input <- df2 %>%
      group_by(celltype6, donor, group_lab) %>%
      summarise(val=median(expr), .groups="drop")
  } else {
    stat_input <- df2 %>%
      transmute(celltype6=celltype6, group_lab=group_lab, val=expr)
  }
  
  ## 保存中间结果（补充材料）
  saveRDS(stat_input, file.path(out_dir, paste0("INTERMEDIATE_stat_input_AD_vs_Control_", unit, ".rds")))
  write.csv(stat_input, file.path(out_dir, paste0("INTERMEDIATE_stat_input_AD_vs_Control_", unit, ".csv")),
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
  
  ## 颜色（确保不灰）
  fill_vals <- c("red","blue")
  names(fill_vals) <- c(lab_dis, lab_ctl)
  
  ## ★Y轴=Cell count（不指定 y=after_stat(density)）
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
                       paste0(dataset_tag,"_",gene_symbol,
                              "_OverlapHistCount_AD_vs_Control_", unit, ".png"))
  
  ragg::agg_png(out_png, width=10, height=6, units="in", res=300, background="white")
  print(p)
  dev.off()
  
  cat("✅ saved:", out_png, "\n")
}

## 输出两张：cell-level + donor-level
plot_one_unit("cell")
plot_one_unit("donor")

sink(log_fp, append=TRUE)
cat("\n==== sessionInfo ====\n")
print(sessionInfo())
cat("==== END ====\n")
sink()

cat("\n🎉 DONE. Output dir:\n", out_dir, "\n", sep="")









#第2章，2 重叠直方图
#donor-level  ，cell-level  2种方法的各自的 重叠直方图  ，Y轴是  Density  
###############################################################################
###############################################################################
###############################################################################
## GSE157827：Overlap Hist (Density) | AD vs Control | donor-level + cell-level
## 风格严格对齐 syn52082747 NO4_03（图3/4）
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
## 0) 路径
## =========================
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"
obj_fp  <- file.path(res_dir, "stepH_obj_celltype6_named.rds")
stopifnot(file.exists(obj_fp))

out_dir <- file.path(res_dir, "NOxx_overlap_hist_density_compare_syn520style")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

dataset_tag <- "GSE157827"
gene_symbol <- "UBL3"
binwidth    <- 0.2
disease     <- "AD"

log_fp <- file.path(out_dir, "NOxx_log.txt")
sink(log_fp)
cat("==== START ====\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("obj_fp:", obj_fp, "\n")
cat("out_dir:", out_dir, "\n")
cat("SEED:", SEED, "\n\n")
sink()

## =========================
## 1) 读对象 + 自动识别列
## =========================
obj <- readRDS(obj_fp)
DefaultAssay(obj) <- "RNA"
md <- obj@meta.data

ct_col_candidates <- c("celltype6","celltype","cell_type","celltype6_named","celltype6_manual")
ct_col <- ct_col_candidates[ct_col_candidates %in% colnames(md)][1]
if (is.na(ct_col)) stop("❌ meta.data 找不到 celltype6 列。")

grp_col_candidates <- c("group4","group2","group","Group","diagnosis","Dx","clinical_diagnosis")
grp_col <- grp_col_candidates[grp_col_candidates %in% colnames(md)][1]
if (is.na(grp_col)) stop("❌ meta.data 找不到 group 列。")

donor_col_candidates <- c("autopsy_id","donor","Donor","sample","Sample","orig.ident","patient","subject")
don_col <- donor_col_candidates[donor_col_candidates %in% colnames(md)][1]
if (is.na(don_col)) stop("❌ meta.data 找不到 donor 列。")

## 统一 Control/NC 名称
grp_raw <- trimws(as.character(md[[grp_col]]))
ctrl_alias <- c("Control","CTRL","Ctr","CTR","NC","Normal","N","control","ctrl","ctr","nc","normal")
grp_std <- ifelse(grp_raw %in% ctrl_alias, "Control", grp_raw)

md$celltype6 <- as.character(md[[ct_col]])
md$group4    <- grp_std
md$donor     <- trimws(as.character(md[[don_col]]))

if (!all(c("AD","Control") %in% unique(md$group4))) {
  sink(log_fp, append=TRUE)
  cat("⚠ group unique:\n"); print(sort(unique(md$group4)))
  sink()
  stop("❌ 分组里没有同时包含 AD 和 Control（或 NC→Control）。")
}

## donor QC（补充材料）
don_all <- unique(md[, c("donor","group4")])
qc_don <- as.data.frame(table(don_all$group4), stringsAsFactors=FALSE)
colnames(qc_don) <- c("group4","n_donors")
write.csv(qc_don, file.path(out_dir, "QC_donors_by_group4.csv"), row.names=FALSE)

sink(log_fp, append=TRUE)
cat("Loaded cells:", ncol(obj), " genes:", nrow(obj), "\n")
cat("celltype col:", ct_col, " | group col:", grp_col, " | donor col:", don_col, "\n")
cat("Donors by group4:\n"); print(table(don_all$group4))
cat("\n")
sink()

## =========================
## 2) 合并 counts layers（Seurat v5 多 layers 兼容）
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
    if (length(mats)==0) stop("❌ counts layers 存在，但读取 layerData 失败。")
    
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
    
    dup <- duplicated(colnames(mat_all))
    if (any(dup)) {
      sink(log_fp, append=TRUE)
      cat("⚠ duplicated cells across layers:", sum(dup), " -> keep first\n")
      sink()
      mat_all <- mat_all[, !dup, drop=FALSE]
    }
    
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

sink(log_fp, append=TRUE)
cat("Merged counts dim:", paste(dim(rna_counts), collapse=" x "), "\n\n")
sink()

## =========================
## 3) 计算 UBL3 log1p(CP10k) + df0(expr>0)
## =========================
gene_row <- if ("UBL3" %in% rownames(rna_counts)) "UBL3" else
  if ("ENSG00000122042" %in% rownames(rna_counts)) "ENSG00000122042" else NA_character_
if (is.na(gene_row)) stop("❌ counts 行名中找不到 UBL3（UBL3 或 ENSG00000122042）")

lib_size <- Matrix::colSums(rna_counts)
expr <- log1p((as.numeric(rna_counts[gene_row, , drop=TRUE]) / pmax(lib_size,1)) * 1e4)

df_all <- data.frame(
  expr      = expr,
  donor     = md$donor,
  group4    = md$group4,
  celltype6 = md$celltype6,
  stringsAsFactors = FALSE
)

df0 <- df_all[df_all$expr > 0, ]   # only expressed cells
saveRDS(df0, file.path(out_dir, "INTERMEDIATE_df0_exprGT0.rds"))
write.csv(df0, file.path(out_dir, "INTERMEDIATE_df0_exprGT0.csv"), row.names=FALSE)

sink(log_fp, append=TRUE)
cat("Expressed cells (expr>0):", nrow(df0), "\n")
cat("Unique donors in expr>0:", length(unique(df0$donor)), "\n\n")
sink()

## =========================
## 4) 画图函数：严格复刻 syn520 NO4_03（关键点：颜色/标题/legend/n）
## =========================
plot_one_unit <- function(unit=c("donor","cell")) {
  
  unit <- match.arg(unit)
  
  ## AD vs Control
  df2 <- df0 %>%
    filter(group4 %in% c(disease,"Control")) %>%
    mutate(group = ifelse(group4==disease, disease, "Control"))
  
  ## pair 中间结果保存
  saveRDS(df2, file.path(out_dir, "INTERMEDIATE_df2_AD_vs_Control_exprGT0.rds"))
  write.csv(df2, file.path(out_dir, "INTERMEDIATE_df2_AD_vs_Control_exprGT0.csv"), row.names=FALSE)
  
  ## donor n（legend 必须写出来）
  don_pair <- unique(df2[, c("donor","group")])
  n_dis <- sum(don_pair$group == disease)
  n_ctl <- sum(don_pair$group == "Control")
  
  lab_dis <- sprintf("%s\n(n = %d)", disease,  n_dis)
  lab_ctl <- sprintf("Control\n(n = %d)",     n_ctl)
  
  ## ★关键：group_lab 的 levels 必须与 fill_vals 的 names 精确一致，否则会变灰
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
  
  ## ★关键：颜色名字必须=levels(df2$group_lab)，否则 ggplot 会退回默认灰
  fill_vals <- c("red","blue")
  names(fill_vals) <- c(lab_dis, lab_ctl)
  
  ## ★标题只保留一行（不写 subtitle）
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
      fill = "group_lab"   # syn520 的 legend 标题就是 group_lab
    ) +
    theme_bw() +
    theme(plot.margin = margin(10,25,10,10)) +
    coord_cartesian(clip="off")
  
  out_png <- file.path(out_dir,
                       paste0(dataset_tag,"_",gene_symbol,
                              "_OverlapHistDensity_AD_vs_Control_", unit, ".png"))
  
  ragg::agg_png(out_png, width=10, height=6, units="in", res=300, background="white")
  print(p)
  dev.off()
  
  cat("✅ saved:", out_png, "\n")
  
  sink(log_fp, append=TRUE)
  cat("---- AD vs Control | unit=", unit, "\n", sep="")
  cat("Legend donor n: AD=", n_dis, " | Control=", n_ctl, "\n", sep="")
  cat("Saved figure:", out_png, "\n\n")
  sink()
}

plot_one_unit("cell")
plot_one_unit("donor")

sink(log_fp, append=TRUE)
cat("\n==== sessionInfo ====\n")
print(sessionInfo())
cat("==== END ====\n")
sink()

cat("\n🎉 DONE. Output dir:\n", out_dir, "\n", sep="")
cat("Log:\n", log_fp, "\n", sep="")










#第2章，内参基因SUMO的 重叠直方图
#donor-level  ，重叠直方图  ，Y轴是Density  ，
###############################################################################
## GSE157827：SUMO(1/2/3) Overlap Hist (Density) | AD vs Control | donor-level
## 风格严格对齐你给的 UBL3 脚本与图：facet、颜色(red/blue)、legend(n=donor)、
## X/Y 轴命名、标题、右上角 MWU Padj 标签、free_y、ragg 输出等。
##
## 重要：本脚本仅做 donor-level（用每个 donor 在每个 celltype 的 median(expr)）
##      的重叠直方图；不做 cell-level（避免伪重复）。
###############################################################################

###############################################################################
## GSE157827：SUMO1/2/3 Overlap Hist (Density) | AD vs Control
##
## 重要说明（避免混淆）：
## - Y轴是 Density：只是直方图的归一化方式（面积≈1），与 cell-level / donor-level 无关
## - 本脚本“图形形状”使用 cell-level（每个细胞一行）的 expr 分布（仅 expr>0）
## - 右上角统计（MWU）严格使用 donor-level：每 donor×celltype 先汇总成一个值（median）
##   → 这样既能看“分布形态”，又能避免伪重复做推断
##
## 输出内容（仅这些）：
## 1) SUMO1/2/3：cell-level overlap histogram（Y=Density） + donor-level MWU(BH) 标签
## 2) 中间结果：df0/df2/stat_input/stats/n_by_celltype/summary + log
## 不输出：donor-median 的箱线/散点图
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
## 0) 路径与参数
## =========================
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"
obj_fp  <- file.path(res_dir, "stepH_obj_celltype6_named.rds")
stopifnot(file.exists(obj_fp))

## 输出目录（仍在原 results 下）
out_dir <- file.path(res_dir, "NOxx_overlap_hist_density_SUMO_cellHist_donorMWU_syn520style")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

dataset_tag <- "GSE157827"
disease     <- "AD"
binwidth    <- 0.2
gene_list   <- c("SUMO1","SUMO2","SUMO3")

log_fp <- file.path(out_dir, "NOxx_log_SUMO_cellHist_donorMWU.txt")
sink(log_fp)
cat("==== START ====\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("obj_fp:", obj_fp, "\n")
cat("out_dir:", out_dir, "\n")
cat("SEED:", SEED, "\n")
cat("binwidth:", binwidth, "\n")
cat("genes:", paste(gene_list, collapse=", "), "\n\n")
sink()

## =========================
## 1) 读对象 + 自动识别 meta 列
## =========================
obj <- readRDS(obj_fp)
DefaultAssay(obj) <- "RNA"
md <- obj@meta.data

ct_col_candidates <- c("celltype6","celltype","cell_type","celltype6_named","celltype6_manual")
ct_col <- ct_col_candidates[ct_col_candidates %in% colnames(md)][1]
if (is.na(ct_col)) stop("❌ meta.data 找不到 celltype 列。")

grp_col_candidates <- c("group4","group2","group","Group","diagnosis","Dx","clinical_diagnosis")
grp_col <- grp_col_candidates[grp_col_candidates %in% colnames(md)][1]
if (is.na(grp_col)) stop("❌ meta.data 找不到 group 列。")

donor_col_candidates <- c("autopsy_id","donor","Donor","sample","Sample","orig.ident","patient","subject")
don_col <- donor_col_candidates[donor_col_candidates %in% colnames(md)][1]
if (is.na(don_col)) stop("❌ meta.data 找不到 donor 列。")

## 统一 Control/NC 名称（与 UBL3 脚本一致）
grp_raw <- trimws(as.character(md[[grp_col]]))
ctrl_alias <- c("Control","CTRL","Ctr","CTR","NC","Normal","N","control","ctrl","ctr","nc","normal")
grp_std <- ifelse(grp_raw %in% ctrl_alias, "Control", grp_raw)

md$celltype6 <- as.character(md[[ct_col]])
md$group4    <- grp_std
md$donor     <- trimws(as.character(md[[don_col]]))

if (!all(c(disease, "Control") %in% unique(md$group4))) {
  sink(log_fp, append=TRUE)
  cat("⚠ group unique:\n"); print(sort(unique(md$group4)))
  sink()
  stop("❌ 分组里没有同时包含 AD 和 Control（或 NC→Control）。")
}

## donor QC（全局 donor 数）
don_all <- unique(md[, c("donor","group4")])
qc_don <- as.data.frame(table(don_all$group4), stringsAsFactors=FALSE)
colnames(qc_don) <- c("group4","n_donors")
write.csv(qc_don, file.path(out_dir, "QC_donors_by_group4.csv"), row.names=FALSE)

sink(log_fp, append=TRUE)
cat("Loaded cells:", ncol(obj), " genes:", nrow(obj), "\n")
cat("celltype col:", ct_col, " | group col:", grp_col, " | donor col:", don_col, "\n")
cat("Donors by group4:\n"); print(table(don_all$group4))
cat("\n")
sink()

## =========================
## 2) 合并 counts layers（Seurat v5 多 layers 兼容；与 UBL3 脚本一致）
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
    if (length(mats)==0) stop("❌ counts layers 存在，但读取 layerData 失败。")
    
    ## 对齐 gene 顺序（防止不同 layer 行顺序不一致）
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
    
    ## 去掉重复 cell（多 layer 合并后可能出现重复列名）
    dup <- duplicated(colnames(mat_all))
    if (any(dup)) {
      sink(log_fp, append=TRUE)
      cat("⚠ duplicated cells across layers:", sum(dup), " -> keep first\n")
      sink()
      mat_all <- mat_all[, !dup, drop=FALSE]
    }
    
    ## 补齐缺失 cell（如果 obj 中有 cell 不在 mat_all 里，用 0 填充）
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
    
    ## 按 obj 的 cell 顺序对齐
    mat_all <- mat_all[, colnames(obj), drop=FALSE]
    return(mat_all)
  }
  
  ## Seurat v4 兜底
  m2 <- tryCatch(Seurat::GetAssayData(obj, assay=assay, slot="counts"), error=function(e) NULL)
  if (!is.null(m2)) return(m2)
  
  stop("❌ 无法获取 counts。")
}

rna_counts <- get_counts_matrix_allcells(obj, "RNA")
stopifnot(ncol(rna_counts) == ncol(obj))

sink(log_fp, append=TRUE)
cat("Merged counts dim:", paste(dim(rna_counts), collapse=" x "), "\n\n")
sink()

## 每个 cell 的 library size（总 counts），用于 CP10k 归一化
lib_size <- Matrix::colSums(rna_counts)

## =========================
## 3) 基因行名定位函数（不“瞎编”ENSG；仅在可用时尝试 SYMBOL→ENSEMBL）
## =========================
locate_gene_row <- function(counts_mat, gene_symbol) {
  rn <- rownames(counts_mat)
  
  ## 1) 直接匹配（最常见：rownames 是 SYMBOL）
  if (gene_symbol %in% rn) return(gene_symbol)
  
  ## 2) 大小写不敏感匹配
  idx_ci <- which(toupper(rn) == toupper(gene_symbol))
  if (length(idx_ci) == 1) return(rn[idx_ci[1]])
  
  ## 3) 如果 rownames 是 ENSG：且你安装了 org.Hs.eg.db，则尝试 SYMBOL→ENSEMBL
  if (requireNamespace("org.Hs.eg.db", quietly=TRUE) &&
      requireNamespace("AnnotationDbi", quietly=TRUE)) {
    
    ens_tbl <- tryCatch(
      AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                            keys=gene_symbol, keytype="SYMBOL", columns=c("ENSEMBL")),
      error=function(e) NULL
    )
    
    if (!is.null(ens_tbl) && "ENSEMBL" %in% colnames(ens_tbl)) {
      ens_ids <- unique(ens_tbl$ENSEMBL[!is.na(ens_tbl$ENSEMBL)])
      
      ## 3a) 精确 ENSG 匹配
      m1 <- intersect(ens_ids, rn)
      if (length(m1) >= 1) return(m1[1])
      
      ## 3b) ENSG.版本号（去掉末尾 .数字 再匹配）
      rn_strip <- sub("\\.\\d+$", "", rn)
      hit <- which(rn_strip %in% ens_ids)
      if (length(hit) >= 1) return(rn[hit[1]])
    }
  }
  
  ## 4) 兜底：给出 grep 候选，方便你检查 rownames 格式
  cand <- rn[grep(gene_symbol, rn, ignore.case=TRUE)]
  stop(paste0(
    "❌ counts 行名中找不到基因：", gene_symbol, "\n",
    "请检查 rownames(rna_counts) 是否为 SYMBOL 或 ENSG。\n",
    "grep 候选（前 20 个）：", paste(head(cand, 20), collapse=", ")
  ))
}

## =========================
## 4) 单基因：cell-level 直方图 + donor-level MWU（median）+ 输出核查/统计
## =========================
run_one_gene <- function(gene_symbol) {
  
  sink(log_fp, append=TRUE)
  cat("==== Gene:", gene_symbol, "====\n")
  sink()
  
  ## 4.1 找到该基因在 counts 里的行名
  gene_row <- locate_gene_row(rna_counts, gene_symbol)
  sink(log_fp, append=TRUE); cat("gene_row used:", gene_row, "\n"); sink()
  
  ## 4.2 计算 log1p(CP10k)
  ## expr = log1p( (counts_gene / library_size) * 1e4 )
  expr <- log1p((as.numeric(rna_counts[gene_row, , drop=TRUE]) / pmax(lib_size, 1)) * 1e4)
  
  df_all <- data.frame(
    expr      = expr,
    donor     = md$donor,
    group4    = md$group4,
    celltype6 = md$celltype6,
    stringsAsFactors = FALSE
  )
  
  ## 仅保留表达阳性细胞（expr>0）
  df0 <- df_all[df_all$expr > 0, ]
  
  ## 中间结果：df0（所有组）
  saveRDS(df0,  file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_df0_exprGT0.rds")))
  write.csv(df0, file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_df0_exprGT0.csv")),
            row.names=FALSE)
  
  ## 只取 AD vs Control
  df2 <- df0 %>%
    filter(group4 %in% c(disease, "Control")) %>%
    mutate(group = ifelse(group4 == disease, disease, "Control"))
  
  ## 中间结果：df2（AD vs Control）
  saveRDS(df2,  file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_df2_AD_vs_Control_exprGT0.rds")))
  write.csv(df2, file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_df2_AD_vs_Control_exprGT0.csv")),
            row.names=FALSE)
  
  ## legend 显示 donor n（按 df2 中“有 expr>0 细胞”的 donor 计数；与 UBL3 逻辑一致）
  don_pair <- unique(df2[, c("donor","group")])
  n_dis <- sum(don_pair$group == disease)
  n_ctl <- sum(don_pair$group == "Control")
  
  lab_dis <- sprintf("%s\n(n = %d)", disease,  n_dis)
  lab_ctl <- sprintf("Control\n(n = %d)",     n_ctl)
  
  ## ★关键：levels 与颜色 names 必须一致，否则 ggplot 会变灰
  df2$group_lab <- factor(ifelse(df2$group == disease, lab_dis, lab_ctl),
                          levels=c(lab_dis, lab_ctl))
  
  ## 4.3 donor-level 统计输入：每 donor×celltype 对 expr 取 median
  stat_input <- df2 %>%
    group_by(celltype6, donor, group_lab) %>%
    summarise(val = median(expr), .groups="drop")
  
  saveRDS(stat_input,  file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_stat_input_donorMedian.rds")))
  write.csv(stat_input, file.path(out_dir, paste0("INTERMEDIATE_", gene_symbol, "_stat_input_donorMedian.csv")),
            row.names=FALSE)
  
  ## 核查：每个 celltype 每组进入统计的 donor 数（用于确认没有异常大量丢失）
  n_by_celltype <- stat_input %>%
    group_by(celltype6, group_lab) %>%
    summarise(n_donors_with_exprGT0 = n(), .groups="drop")
  write.csv(n_by_celltype,
            file.path(out_dir, paste0("CHECK_", gene_symbol, "_n_donors_by_celltype_AD_vs_Control_donor.csv")),
            row.names=FALSE)
  
  ## 4.4 MWU（严格 donor-level median）+ BH
  stats <- stat_input %>%
    group_by(celltype6) %>%
    summarise(
      p_raw = {
        g <- group_lab; v <- val
        ## 某 celltype 若只有一组 donor 有值，则无法做 MWU，记为 NA
        if (length(unique(g)) < 2) NA_real_
        else tryCatch(wilcox.test(v ~ g, exact=FALSE)$p.value, error=function(e) NA_real_)
      },
      .groups="drop"
    )
  
  stats$padj  <- p.adjust(stats$p_raw, method="BH")
  stats$label <- sprintf("Mann–Whitney U\nPadj=%.2e", stats$padj)
  stats$x <- Inf; stats$y <- Inf
  
  write.csv(stats,
            file.path(out_dir, paste0("STATS_", gene_symbol, "_MWU_donorMedian_BH_by_celltype.csv")),
            row.names=FALSE)
  
  ## 4.5 作图：cell-level overlap histogram（Y=Density）+ 右上角 donor-level MWU 标签
  fill_vals <- c("red","blue")
  names(fill_vals) <- c(lab_dis, lab_ctl)
  
  ## 标题建议写清楚“图=cell-level，检验=donor-level”，防止你自己/审稿人误解
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
                              "_OverlapHistDensity_AD_vs_Control_cellHist_MWUdonorMedian.png"))
  
  ragg::agg_png(out_png, width=10, height=6, units="in", res=300, background="white")
  print(p)
  dev.off()
  
  sink(log_fp, append=TRUE)
  cat("Saved figure:", out_png, "\n")
  cat("Legend donor n (expr>0 donors): AD=", n_dis, " | Control=", n_ctl, "\n", sep="")
  cat("Expressed cells (expr>0):", nrow(df0), "\n\n")
  sink()
  
  invisible(list(fig=out_png, stats=stats))
}

## =========================
## 5) 批量运行 SUMO1/2/3 + 汇总统计表
## =========================
stats_long <- list()

for (g in gene_list) {
  run_one_gene(g)
  st_fp <- file.path(out_dir, paste0("STATS_", g, "_MWU_donorMedian_BH_by_celltype.csv"))
  st <- read.csv(st_fp, stringsAsFactors=FALSE)
  st$gene <- g
  stats_long[[g]] <- st
}

stats_long_df <- bind_rows(stats_long) %>%
  select(gene, celltype6, p_raw, padj, label)

write.csv(stats_long_df,
          file.path(out_dir, "SUMMARY_SUMO_genes_MWU_BH_by_celltype.csv"),
          row.names=FALSE)

sink(log_fp, append=TRUE)
cat("\n==== sessionInfo ====\n")
print(sessionInfo())
cat("==== END ====\n")
sink()

cat("\n🎉 DONE.\nOutput dir:\n", out_dir, "\n", sep="")
cat("Log:\n", log_fp, "\n", sep="")
cat("仅输出：cell-level 重叠直方图（expr>0，Y=Density）+ donor-level MWU(BH) 标签。\n")









#第3章，箱线图
#第3章1，6个细胞类型的箱线图
###############################################################################
## NOxx_GSE157827_pseudobulk_DESeq2_byDonor_UBL3_SUMO_ADvsControl.R
##
## 【数据集】GSE157827：只有 AD vs Control（NC 已映射为 Control）
## 【统计单位】donor（先 pseudo-bulk 到 celltype6×donor，再 DESeq2）
##
## A) UBL3：
##   - 6 个 celltype6：每个 celltype 输出 1 张 donor-level 箱线图（模板风格）
##   - 每个 celltype 输出 1 张 DEG 表：AD vs Control（全基因）
##   - 输出 6-panel 总图（2×3）
##
## B) SUMO 内参筛选：
##   - 候选：SUMO1/2/3
##   - 规则：在每个 celltype6 中，AD vs Control 的 padj > 0.05
##   - 通过者：输出同款 6 张箱线图 + 6-panel 总图
##
## 【必须核查的中间结果】
##   - donor->group 一对一映射（不满足直接 stop）
##   - QC_donor_counts_by_group.csv
##   - QC_donor_counts_by_celltype6_by_group.csv
##   - QC_cells_by_celltype6_by_group.csv
##   - INTERMEDIATE_pseudobulk_coldata_celltype6_donor_group.csv
###############################################################################
## NOxx_GSE157827_UBL3_Boxplots_DESeq2_byDonor_ADvsControl_FINAL.R
##
## 【适用数据】GSE157827（只有 AD vs Control）
## 【统计单位】donor（先 pseudo-bulk，再 DESeq2）
##
## 【本脚本做什么】
##  1) 使用你已经核查通过的映射：celltype6 → celltype7（Astro/Endo/…）
##  2) 对每个 celltype（6 个）：
##     - 进行 donor-level pseudo-bulk
##     - 跑 DESeq2（AD vs Control）
##     - 输出：
##         a) 1 张 UBL3 箱线图
##         b) 1 张 DEG 表（全基因）
##  3) 输出 6-panel 总图
##  4) 保存所有关键中间结果，便于核查与补充材料
###############################################################################

###############################################################################
## NOxx_GSE157827_UBL3_Boxplots_DESeq2_byDonor_ADvsControl_FIXED.R
##
## 修复点：
##  1) 自动识别分组列（避免 md$group4 不存在导致 logical(0)）
##  2) celltype 使用 celltype7（你已核查）
##  3) 仅保留 AD vs Control
###############################################################################

rm(list = ls()); gc()
Sys.setenv(LANG = "en")
set.seed(20251023)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

## =========================
## 0) 路径
## =========================
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"
obj_fp  <- file.path(res_dir, "stepH_obj_labeled_celltype7_celltype6.rds")
stopifnot(file.exists(obj_fp))

out_dir <- file.path(res_dir, "NOxx_GSE157827_UBL3_Boxplots_DESeq2_byDonor_ADvsControl")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_fp <- file.path(out_dir, "NOxx_log.txt")
sink(log_fp)
cat("==== START ====\n")
cat("Time:", Sys.time(), "\n")
cat("obj :", obj_fp, "\n")
cat("out :", out_dir, "\n\n")
sink()

## =========================
## 1) 读对象 + 取 meta
## =========================
obj <- readRDS(obj_fp)
DefaultAssay(obj) <- "RNA"
md <- obj@meta.data

## =========================
## 2) celltype：直接用你已核查的 celltype7（名称）
## =========================
if (!("celltype7" %in% colnames(md))) {
  stop("❌ meta.data 中没有 celltype7 列。请检查对象列名。")
}
md$celltype6 <- trimws(as.character(md$celltype7))

## =========================
## 3) 自动识别 group 列，并统一成 AD/Control
## =========================
## 你对象里可能不是 group4，所以必须自动找
grp_candidates <- c("group4","group2","group","Group","diagnosis","Dx","clinical_diagnosis")
grp_col <- grp_candidates[grp_candidates %in% colnames(md)][1]
if (is.na(grp_col)) {
  stop("❌ 找不到分组列（group4/group2/group/diagnosis/...）。当前列名：\n",
       paste(colnames(md), collapse=", "))
}

grp_raw <- trimws(as.character(md[[grp_col]]))

## 把 NC 等映射成 Control
ctrl_alias <- c("Control","CTRL","Ctr","CTR","NC","Normal","N","control","ctrl","ctr","nc","normal")
md$group <- ifelse(grp_raw %in% ctrl_alias, "Control", grp_raw)

## 只保留 AD/Control（如果还有其它组一律扔掉）
md <- md[md$group %in% c("AD","Control"), , drop = FALSE]

## =========================
## 4) donor 列（优先 autopsy_id；如果没有则自动找）
## =========================
don_candidates <- c("autopsy_id","donor","sample","orig.ident","patient","subject")
don_col <- don_candidates[don_candidates %in% colnames(md)][1]
if (is.na(don_col)) {
  stop("❌ 找不到 donor 列（autopsy_id/donor/sample/...）。当前列名：\n",
       paste(colnames(md), collapse=", "))
}
md$donor <- trimws(as.character(md[[don_col]]))

## 去掉空 donor / 空 celltype
md <- md[md$donor != "" & !is.na(md$donor) & md$celltype6 != "" & !is.na(md$celltype6), , drop=FALSE]

## 记录：本次到底用了哪几列
sink(log_fp, append=TRUE)
cat("Used columns:\n")
cat("  group column   =", grp_col, "\n")
cat("  donor column   =", don_col, "\n")
cat("  celltype column= celltype7\n\n")
cat("Group counts:\n"); print(table(md$group))
cat("Celltype counts:\n"); print(table(md$celltype6))
cat("Unique donors:", length(unique(md$donor)), "\n\n")
sink()

## =========================
## 5) donor -> group 一对一校验（必须）
## =========================
donor_map <- unique(md[, c("donor","group")])
if (any(table(donor_map$donor) > 1)) {
  bad <- donor_map[donor_map$donor %in% names(which(table(donor_map$donor) > 1)), ]
  write.csv(bad, file.path(out_dir, "CHECK_donor_maps_to_multiple_group.csv"), row.names=FALSE)
  stop("❌ donor 对应多个 group，已输出 CHECK_donor_maps_to_multiple_group.csv")
}

write.csv(as.data.frame(table(donor_map$group)),
          file.path(out_dir, "QC_donor_counts_by_group.csv"),
          row.names=FALSE)

tmp_donor_ct <- unique(md[, c("donor","celltype6","group")])
write.csv(as.data.frame(table(tmp_donor_ct$celltype6, tmp_donor_ct$group)),
          file.path(out_dir, "QC_donor_counts_by_celltype6_by_group.csv"),
          row.names=FALSE)

write.csv(as.data.frame(table(md$celltype6, md$group)),
          file.path(out_dir, "QC_cells_by_celltype6_by_group.csv"),
          row.names=FALSE)

## =========================
## 6) counts 获取（兼容 v5 多 layers）
## =========================
get_counts <- function(obj) {
  a <- obj[["RNA"]]
  layers <- tryCatch(Layers(a), error=function(e) character(0))
  cl <- layers[grepl("^counts", layers)]
  
  if (length(cl) > 0) {
    mats <- lapply(cl, function(x) LayerData(a, layer=x))
    ref <- rownames(mats[[1]])
    mats <- lapply(mats, function(m){
      m2 <- Matrix(0, nrow=length(ref), ncol=ncol(m), sparse=TRUE)
      rownames(m2) <- ref; colnames(m2) <- colnames(m)
      common <- intersect(ref, rownames(m))
      m2[common, ] <- m[common, , drop=FALSE]
      m2
    })
    m <- if (length(mats)==1) mats[[1]] else Reduce(Matrix::cbind2, mats)
    m[, colnames(obj), drop=FALSE]
  } else {
    GetAssayData(obj, slot="counts")
  }
}

cnt_all <- get_counts(obj)
cnt <- cnt_all[, rownames(md), drop=FALSE]

## =========================
## 7) 找 UBL3 行名
## =========================
ensg <- AnnotationDbi::select(org.Hs.eg.db, keys="UBL3", keytype="SYMBOL", columns="ENSEMBL")$ENSEMBL[1]
gene_row <- if (!is.na(ensg) && ensg %in% rownames(cnt)) ensg else "UBL3"
if (!(gene_row %in% rownames(cnt))) stop("❌ counts 里找不到 UBL3（symbol/ensg 都不匹配）")

write.csv(data.frame(gene="UBL3", gene_row=gene_row),
          file.path(out_dir, "CHECK_geneSymbol_to_rowname.csv"),
          row.names=FALSE)

## =========================
## 8) pseudo-bulk：celltype6 × donor（稀疏聚合）
## =========================
pb_key <- paste(md$celltype6, md$donor, sep="__")
grp <- factor(pb_key, levels=unique(pb_key))

M <- sparseMatrix(
  i = seq_along(grp),
  j = as.integer(grp),
  x = 1,
  dims = c(length(grp), length(levels(grp))),
  dimnames = list(rownames(md), levels(grp))
)

pb <- cnt %*% M
pb <- as(pb, "dgCMatrix")

pb_meta <- data.frame(
  key       = colnames(pb),
  celltype6 = sub("__.*","", colnames(pb)),
  donor     = sub(".*__","", colnames(pb)),
  group     = donor_map$group[match(sub(".*__","", colnames(pb)), donor_map$donor)],
  stringsAsFactors = FALSE
)

saveRDS(pb, file.path(out_dir, "INTERMEDIATE_pseudobulk_matrix.rds"))
write.csv(pb_meta, file.path(out_dir, "INTERMEDIATE_pseudobulk_coldata.csv"), row.names=FALSE)

## =========================
## 9) 每个 celltype：DESeq2 + UBL3 箱线图 + DEG
## =========================
pal2 <- c(AD="#D24B40", Control="#2C7FB8")
plots <- list()

for (ct in sort(unique(pb_meta$celltype6))) {
  
  cols <- pb_meta$key[pb_meta$celltype6 == ct]
  
  coldata <- data.frame(
    group = factor(pb_meta$group[pb_meta$celltype6 == ct], levels=c("Control","AD")),
    row.names = cols
  )
  
  y <- round(as.matrix(pb[, cols, drop=FALSE]))
  storage.mode(y) <- "integer"
  
  dds <- DESeqDataSetFromMatrix(y, coldata, design=~group)
  dds <- DESeq(dds, quiet=TRUE)
  
  res <- results(dds, contrast=c("group","AD","Control"))
  
  ## 输出 DEG（全基因）
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[,c("gene","log2FoldChange","padj","pvalue","baseMean")]
  write.csv(res_df,
            file.path(out_dir, paste0("DEG_UBL3_", ct, "_AD_vs_Control.csv")),
            row.names=FALSE)
  
  ## 作图数据（donor-level）
  norm <- counts(dds, normalized=TRUE)
  dfp <- data.frame(
    donor = pb_meta$donor[pb_meta$celltype6==ct],
    group = factor(pb_meta$group[pb_meta$celltype6==ct], levels=c("AD","Control")),
    value = as.numeric(norm[gene_row, cols]),
    stringsAsFactors = FALSE
  )
  write.csv(dfp,
            file.path(out_dir, paste0("INTERMEDIATE_plotdata_UBL3_", ct, ".csv")),
            row.names=FALSE)
  
  subtitle <- sprintf("AD vs Control : log2FC=%.3f, padj=%s",
                      as.numeric(res[gene_row,"log2FoldChange"]),
                      ifelse(is.na(res[gene_row,"padj"]), "NA", format(res[gene_row,"padj"], digits=3, scientific=TRUE)))
  
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
      plot.title    = element_text(face = "bold", size = 16, margin = margin(b = 4)),
      plot.subtitle = element_text(size = 8.8, colour = "grey15", margin = margin(b = 6)),
      legend.position = "none"
    )
  
  
  out_png <- file.path(out_dir, paste0("UBL3_Box_byDonor_", ct, ".png"))
  ggsave(out_png, p, width=6.5, height=4.2, dpi=300)
  
  plots[[ct]] <- p
}

## 6-panel（按你 celltype7 的顺序）
panel_order <- c("Astro","Endo","Excit","Inhib","Microgl","Oligo")
panel_order <- panel_order[panel_order %in% names(plots)]
if (length(panel_order)==0) {
  stop("❌ 6-panel 为空：plots 的名字是：", paste(names(plots), collapse=", "))
}

p_all <- wrap_plots(plots[panel_order], ncol=3)
ggsave(file.path(out_dir, "UBL3_Box_byDonor_6panel.png"),
       p_all, width=14, height=7.5, dpi=300)

sink(log_fp, append=TRUE)
cat("\n==== END ====\n")
print(sessionInfo())
sink()

cat("\n🎉 完成：已自动识别 group 列，使用 celltype7 名称，按 donor 做 pseudo-bulk + DESeq2 + 箱线图 + DEG。\n")
cat("输出目录：", out_dir, "\n", sep="")










#SUMO 内参筛选 + 箱线图 + 导出数值（AD vs Control）
###############################################################################
## NOxx2_GSE157827_SUMO_housekeeping_screen_and_boxplots_byDonor_ADvsControl.R
##
## 【输入】优先复用上一章已生成的 pseudo-bulk：
##   - INTERMEDIATE_pseudobulk_matrix.rds
##   - INTERMEDIATE_pseudobulk_coldata.csv
## 若不存在，则自动从对象重新构建（与 UBL3 章一致）
##
## 【输出】
##  1) CHECK_SUMO_stability_ADvsControl_byCelltype.csv   （每 celltype 的 log2FC/pvalue/padj）
##  2) CHECK_SUMO_passed_genes.txt                      （通过的 SUMO 列表）
##  3) 对每个通过的 SUMO：
##      - 6 张箱线图：SUMO*_Box_byDonor_<celltype>.png
##      - 6-panel：SUMO*_Box_byDonor_6panel.png
##      - 每个 celltype 的 donor 值表：VALUES_SUMO*_byDonor_<celltype>.csv
###############################################################################

###############################################################################
## NO11_GSE157827_SUMO_ref_select_and_boxplots_ADvsControl_FINAL.R
##
## 【GSE157827 专用：只有 AD vs Control】
## 【统计单位：donor】pseudo-bulk（celltype6__donor）+ DESeq2（每个 celltype 单独跑）
##
## 1) 计算 SUMO1/2/3 在每个 celltype 的 AD vs Control 差异（log2FC, pvalue, padj）
## 2) 自动筛选最稳定 SUMO 内参：
##    - 硬条件：6 个 celltype 全部 padj > 0.05
##    - 且效应足够小：6 个 celltype 全部 |log2FC| < 0.10   （可改阈值）
##    - 再按“稳定性评分”排序，选第一名
## 3) 用与 UBL3 完全一致的风格画箱线图（donor-level）
##    - title：GENE in Celltype (pseudo-bulk per donor)
##    - subtitle：AD vs Control : log2FC=..., padj=...
##    - y 轴：DESeq2 normalized counts
## 4) 输出被选中 SUMO 的 donor 值表（每 celltype 一张 CSV）
##
## 【输入】
##  优先复用 UBL3 章输出（你已经跑通）：
##   - NOxx_GSE157827_UBL3_Boxplots_DESeq2_byDonor_ADvsControl/INTERMEDIATE_pseudobulk_matrix.rds
##   - NOxx_GSE157827_UBL3_Boxplots_DESeq2_byDonor_ADvsControl/INTERMEDIATE_pseudobulk_coldata.csv
##  若不存在，会从对象 stepH_obj_labeled_celltype7_celltype6.rds 自动重建 pseudo-bulk
###############################################################################

rm(list=ls()); gc()
Sys.setenv(LANG="en")
set.seed(20251023)
options(stringsAsFactors = FALSE)

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
  library(ragg)
})

## =========================
## 0) 路径与参数（只改这里也行）
## =========================
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"
obj_fp  <- file.path(res_dir, "stepH_obj_labeled_celltype7_celltype6.rds")
stopifnot(file.exists(obj_fp))

## 你 UBL3 那章的输出目录（用来复用 pseudo-bulk）
ubl3_dir <- file.path(res_dir, "NOxx_GSE157827_UBL3_Boxplots_DESeq2_byDonor_ADvsControl")
stopifnot(dir.exists(ubl3_dir))

## 本章输出目录
out_dir <- file.path(res_dir, "NO11_SUMO_REF_select_ADvsControl")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

dir_sel  <- file.path(out_dir, "SUMO_REF_selected")   # 最佳 SUMO 图
dir_all3 <- file.path(out_dir, "SUMO_REF_all3")       # 可选：SUMO1/2/3 全部图（本脚本默认都画，便于核查）
dir.create(dir_sel,  recursive=TRUE, showWarnings=FALSE)
dir.create(dir_all3, recursive=TRUE, showWarnings=FALSE)

## SUMO 候选
sumo_genes <- c("SUMO1","SUMO2","SUMO3")

## 筛选阈值（你要更宽松就把 0.10 改成 0.20）
TH_padj <- 0.05
TH_absL2 <- 0.10

## celltype 固定顺序（与你前面一致）
panel_order <- c("Astro","Endo","Excit","Inhib","Microgl","Oligo")

## 配色（两组）
pal2 <- c(
  AD      = "#D24B40",
  Control = "#2C7FB8"
)

## 日志
log_fp <- file.path(out_dir, "NO11_log.txt")
sink(log_fp)
cat("==== NO11 START ====\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("obj_fp :", obj_fp, "\n")
cat("ubl3_dir:", ubl3_dir, "\n")
cat("out_dir :", out_dir, "\n")
cat("TH_padj :", TH_padj, " | TH_absL2 :", TH_absL2, "\n\n")
sink()

cat("✅ out_dir = ", out_dir, "\n", sep="")

## =========================
## 1) 读入/重建 pseudo-bulk（优先复用）
## =========================
pb_fp     <- file.path(ubl3_dir, "INTERMEDIATE_pseudobulk_matrix.rds")
pbmeta_fp <- file.path(ubl3_dir, "INTERMEDIATE_pseudobulk_coldata.csv")

if (file.exists(pb_fp) && file.exists(pbmeta_fp)) {
  pb <- readRDS(pb_fp)
  pb_annot <- read.csv(pbmeta_fp, stringsAsFactors=FALSE)
  sink(log_fp, append=TRUE)
  cat("✅ Reuse pseudo-bulk from UBL3 step.\n")
  cat("pb dim:", paste(dim(pb), collapse=" x "), "\n\n")
  sink()
} else {
  ## 如果你不小心删了中间文件，则自动从对象重建（与 UBL3 章一致）
  sink(log_fp, append=TRUE)
  cat("⚠ pseudo-bulk not found in UBL3 dir -> rebuild from object.\n\n")
  sink()
  
  obj <- readRDS(obj_fp); DefaultAssay(obj) <- "RNA"
  md <- obj@meta.data
  
  ## celltype 用 celltype7（你已核查映射）
  md$celltype6 <- trimws(as.character(md$celltype7))
  
  ## 自动找 group 列
  grp_candidates <- c("group4","group2","group","Group","diagnosis","Dx","clinical_diagnosis")
  grp_col <- grp_candidates[grp_candidates %in% colnames(md)][1]
  if (is.na(grp_col)) stop("❌ 找不到分组列。")
  
  grp_raw <- trimws(as.character(md[[grp_col]]))
  ctrl_alias <- c("Control","CTRL","Ctr","CTR","NC","Normal","N","control","ctrl","ctr","nc","normal")
  md$group <- ifelse(grp_raw %in% ctrl_alias, "Control", grp_raw)
  
  ## donor 列
  don_candidates <- c("autopsy_id","donor","sample","orig.ident","patient","subject")
  don_col <- don_candidates[don_candidates %in% colnames(md)][1]
  if (is.na(don_col)) stop("❌ 找不到 donor 列。")
  md$donor <- trimws(as.character(md[[don_col]]))
  
  ## 只留 AD/Control
  md <- md[md$group %in% c("AD","Control"), ]
  md <- md[md$donor != "" & md$celltype6 != "", ]
  
  ## donor->group 一对一
  donor_map <- unique(md[,c("donor","group")])
  if (any(table(donor_map$donor) > 1)) stop("❌ donor 对应多个 group。")
  
  ## counts 获取（兼容 v5 多 layers）
  get_counts <- function(obj) {
    a <- obj[["RNA"]]
    layers <- tryCatch(Layers(a), error=function(e) character(0))
    cl <- layers[grepl("^counts", layers)]
    if (length(cl) > 0) {
      mats <- lapply(cl, function(x) LayerData(a, layer=x))
      ref <- rownames(mats[[1]])
      mats <- lapply(mats, function(m){
        m2 <- Matrix(0, nrow=length(ref), ncol=ncol(m), sparse=TRUE)
        rownames(m2) <- ref; colnames(m2) <- colnames(m)
        common <- intersect(ref, rownames(m))
        m2[common, ] <- m[common, , drop=FALSE]
        m2
      })
      m <- if (length(mats)==1) mats[[1]] else Reduce(Matrix::cbind2, mats)
      m[, colnames(obj), drop=FALSE]
    } else {
      GetAssayData(obj, slot="counts")
    }
  }
  cnt_all <- get_counts(obj)
  cnt <- cnt_all[, rownames(md), drop=FALSE]
  
  ## pseudo-bulk（celltype6__donor）
  pb_key <- paste(md$celltype6, md$donor, sep="__")
  grp <- factor(pb_key, levels=unique(pb_key))
  M <- sparseMatrix(
    i=seq_along(grp),
    j=as.integer(grp),
    x=1,
    dims=c(length(grp), length(levels(grp))),
    dimnames=list(rownames(md), levels(grp))
  )
  pb <- as(cnt %*% M, "dgCMatrix")
  
  pb_annot <- data.frame(
    key       = colnames(pb),
    celltype6 = sub("__.*","", colnames(pb)),
    donor     = sub(".*__","", colnames(pb)),
    group     = donor_map$group[match(sub(".*__","", colnames(pb)), donor_map$donor)],
    stringsAsFactors=FALSE
  )
  
  saveRDS(pb, file.path(out_dir, "INTERMEDIATE_pseudobulk_matrix.rds"))
  write.csv(pb_annot, file.path(out_dir, "INTERMEDIATE_pseudobulk_coldata.csv"), row.names=FALSE)
}

## =========================
## 2) SUMO 行名定位（SYMBOL/ENSG 兼容）
## =========================
symbol_to_row <- function(symbol, mat) {
  ensg <- tryCatch({
    AnnotationDbi::select(org.Hs.eg.db, keys=symbol, keytype="SYMBOL", columns="ENSEMBL")$ENSEMBL[1]
  }, error=function(e) NA_character_)
  if (!is.na(ensg) && ensg %in% rownames(mat)) return(ensg)
  if (symbol %in% rownames(mat)) return(symbol)
  return(NA_character_)
}

sumo_rows <- sapply(sumo_genes, symbol_to_row, mat=pb)
sumo_rows <- sumo_rows[!is.na(sumo_rows)]

write.csv(
  data.frame(gene=names(sumo_rows), gene_row=unname(sumo_rows)),
  file.path(out_dir, "CHECK_SUMO_geneSymbol_to_rowname.csv"),
  row.names=FALSE
)

if (length(sumo_rows)==0) stop("❌ SUMO1/2/3 都无法在 pb 行名中匹配。")

## =========================
## 3) 统计：每个 celltype 跑 DESeq2，提取 SUMO 的 log2FC/pvalue/padj
## =========================
celltypes <- sort(unique(pb_annot$celltype6))
stats_list <- list()

for (ct in celltypes) {
  
  cols <- pb_annot$key[pb_annot$celltype6 == ct]
  gvec <- pb_annot$group[pb_annot$celltype6 == ct]
  
  ## 必须两组都存在
  if (!all(c("AD","Control") %in% unique(gvec))) next
  
  y <- round(as.matrix(pb[, cols, drop=FALSE]))
  storage.mode(y) <- "integer"
  
  dds <- DESeqDataSetFromMatrix(
    countData = y,
    colData = data.frame(group=factor(gvec, levels=c("Control","AD"))),
    design = ~ group
  )
  dds <- DESeq(dds, quiet=TRUE)
  res <- results(dds, contrast=c("group","AD","Control"))
  
  for (g in names(sumo_rows)) {
    gr <- sumo_rows[[g]]
    stats_list[[length(stats_list)+1]] <- data.frame(
      gene      = g,
      gene_row  = gr,
      celltype6 = ct,
      log2FC    = as.numeric(res[gr, "log2FoldChange"]),
      pvalue    = as.numeric(res[gr, "pvalue"]),
      padj      = as.numeric(res[gr, "padj"]),
      stringsAsFactors=FALSE
    )
  }
}

stats_df <- bind_rows(stats_list)

stats_fp <- file.path(out_dir, "SUMO_stats_all_ADvsControl.csv")
write.csv(stats_df, stats_fp, row.names=FALSE)

## =========================
## 4) 筛选 + 排名（严格：6 个 celltype 都满足 padj>0.05 且 |log2FC|<阈值）
## =========================
rank_df <- stats_df %>%
  mutate(absL2 = abs(log2FC)) %>%
  group_by(gene) %>%
  summarise(
    n_celltypes = n(),
    pass_all_celltypes_padj = all(is.finite(padj) & padj > TH_padj),
    pass_all_celltypes_absL2 = all(is.finite(absL2) & absL2 < TH_absL2),
    padj_median = median(padj, na.rm=TRUE),
    absL2_median = median(absL2, na.rm=TRUE),
    padj_min = min(padj, na.rm=TRUE),
    absL2_max = max(absL2, na.rm=TRUE),
    ## 稳定性评分：padj 越大越好，|log2FC| 越小越好
    score = padj_median - absL2_median,
    .groups="drop"
  ) %>%
  arrange(desc(pass_all_celltypes_padj), desc(pass_all_celltypes_absL2), desc(score))

rank_fp <- file.path(out_dir, "SUMO_stability_ranking_ADvsControl.csv")
write.csv(rank_df, rank_fp, row.names=FALSE)

## 选“同时通过两条硬条件”的第一名；若一个都没通过，则降级：只要求 padj>0.05
best_gene <- NA_character_

cand_strict <- rank_df %>% filter(pass_all_celltypes_padj & pass_all_celltypes_absL2)
if (nrow(cand_strict) > 0) {
  best_gene <- cand_strict$gene[1]
  mode_sel <- "STRICT(padj>0.05 & |log2FC|<TH)"
} else {
  cand_padj <- rank_df %>% filter(pass_all_celltypes_padj)
  if (nrow(cand_padj) > 0) {
    best_gene <- cand_padj$gene[1]
    mode_sel <- "RELAX(padj>0.05 only)"
  } else {
    ## 最后兜底：谁 score 最高就选谁（但会在日志里明确提示）
    best_gene <- rank_df$gene[1]
    mode_sel <- "FALLBACK(score only)"
  }
}

writeLines(best_gene, con=file.path(out_dir, "SUMO_REF_selected_gene.txt"))

sink(log_fp, append=TRUE)
cat("Selected SUMO =", best_gene, " | mode =", mode_sel, "\n\n")
cat("Ranking table:\n"); print(rank_df)
cat("\n")
sink()

cat("🏆 Selected SUMO =", best_gene, " | mode =", mode_sel, "\n")

## =========================
## 5) 画图：最佳 SUMO（6 张 + 6-panel）并导出 donor 值（VALUES_*.csv）
## =========================
best_row <- sumo_rows[[best_gene]]

## 稳定保存 PNG：避免黑图
## =========================
## 稳定写 PNG：先写 temp，再复制到目标路径
## - temp 文件名避免使用 Sys.time()（会带冒号/空格，Windows 易出问题）
## - 校验文件大小，避免“空文件/半文件”也当成功
## - 重试，避免资源管理器/杀毒/同步软件短暂占用导致失败
## =========================
save_png_atomic <- function(filename, plot, width=6.5, height=4.2, dpi=300, retries=5) {
  
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
      
      ## 文件完整性：至少 80KB（你的箱线图通常远大于此）
      if (file.exists(tmpfile) && is.finite(file.info(tmpfile)$size) && file.info(tmpfile)$size > 80*1024) {
        ok <- TRUE
      }
    }, silent = TRUE)
    
    if (ok) {
      ## 复制到目标路径（先删旧文件，避免被占用导致 copy 失败）
      if (file.exists(filename)) file.remove(filename)
      ok2 <- file.copy(tmpfile, filename, overwrite = TRUE)
      
      if (isTRUE(ok2) && file.exists(filename) && file.info(filename)$size > 80*1024) {
        file.remove(tmpfile)
        gc(FALSE)
        return(TRUE)
      }
    }
    
    Sys.sleep(0.4)  # 给系统一点时间释放文件句柄
  }
  
  if (file.exists(tmpfile)) file.remove(tmpfile)
  gc(FALSE)
  return(FALSE)
}


plots <- list()

for (ct in celltypes) {
  
  cols <- pb_annot$key[pb_annot$celltype6 == ct]
  gvec <- pb_annot$group[pb_annot$celltype6 == ct]
  if (!all(c("AD","Control") %in% unique(gvec))) next
  
  y <- round(as.matrix(pb[, cols, drop=FALSE]))
  storage.mode(y) <- "integer"
  
  dds <- DESeqDataSetFromMatrix(
    countData = y,
    colData = data.frame(group=factor(gvec, levels=c("Control","AD"))),
    design = ~ group
  )
  dds <- DESeq(dds, quiet=TRUE)
  res <- results(dds, contrast=c("group","AD","Control"))
  
  ## donor-level normalized counts（输出你要的“值”）
  norm <- counts(dds, normalized=TRUE)
  dfp <- data.frame(
    donor = pb_annot$donor[pb_annot$celltype6==ct],
    group = factor(pb_annot$group[pb_annot$celltype6==ct], levels=c("AD","Control")),
    value = as.numeric(norm[best_row, cols]),
    stringsAsFactors=FALSE
  )
  
  write.csv(dfp, file.path(out_dir, paste0("VALUES_", best_gene, "_byDonor_", ct, ".csv")),
            row.names=FALSE)
  
  subtitle <- sprintf("AD vs Control : log2FC=%.3f, padj=%s",
                      as.numeric(res[best_row,"log2FoldChange"]),
                      ifelse(is.na(res[best_row,"padj"]), "NA", format(res[best_row,"padj"], digits=3, scientific=TRUE)))
  
  ## 图风格：与 UBL3 同款（标题更重点）
  p <- ggplot(dfp, aes(x=group, y=value, fill=group)) +
    geom_boxplot(width=0.55, outlier.shape=NA, alpha=0.95, colour="grey15") +
    geom_point(position=position_jitter(width=0.10),
               size=2.6, alpha=0.90, shape=21, stroke=0.4, colour="grey10") +
    scale_fill_manual(values=pal2, drop=FALSE) +
    labs(
      title    = paste0(best_gene, " in ", ct, " (pseudo-bulk per donor)"),
      subtitle = subtitle,
      x = NULL, y = "Normalized counts"
    ) +
    theme_bw(base_size=13) +
    theme(
      plot.title.position = "plot",
      plot.title    = element_text(face="bold", size=16, margin=margin(b=4)),
      plot.subtitle = element_text(size=8.8, colour="grey15", margin=margin(b=6)),
      legend.position="none"
    )
  
  fp <- file.path(dir_sel, paste0(best_gene, "_Box_byDonor_", ct, ".png"))
  ok <- save_png_atomic(fp, p, width=6.5, height=4.2, dpi=300, retries=5)
  if (!ok) message("⚠ 写图失败（重试后仍失败）：", fp)
  
  
  plots[[ct]] <- p
}

## 6-panel 总图
ord <- panel_order[panel_order %in% names(plots)]
if (length(ord) > 0) {
  p_all <- wrap_plots(plots[ord], ncol=3)
  ok <- save_png_atomic(file.path(dir_sel, paste0(best_gene, "_Box_byDonor_6panel.png")),
                        p_all, width=14, height=7.5, dpi=300, retries=5)
  if (!ok) message("⚠ 写6-panel失败：", file.path(dir_sel, paste0(best_gene, "_Box_byDonor_6panel.png")))
}

sink(log_fp, append=TRUE)
cat("\n==== NO11 END ====\n")
print(sessionInfo())
sink()

cat("\n🎉 NO11 完成：SUMO 内参筛选 + 核实 +（最佳 SUMO）箱线图 + donor 值输出\n")
cat("📁 输出目录：", out_dir, "\n", sep="")













#整体水平的箱线图
###############################################################################
################################################################################
## NO12_GSE157827_WholeCell_ADvsControl_perDonor_UBL3_boxplot_Assay5Layers.R
##
## 适用场景：
##  - Seurat v5, RNA assay 为 Assay5
##  - counts/data 被拆成多个 layers（counts.1.1, data.1.1 ...）
##  - GetAssayData(slot="counts"/"data") 返回 NULL（正常）
##
## 目标：
##  - 严格按 syn52082747 流程：counts -> CP10k -> log1p -> donor mean
##  - donor 为统计单位：Wilcoxon + log2FC
##  - 出图风格/颜色/标题与模式图一致
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
  library(ragg)
  library(tidyr)
})

## ========= 1) 路径（只改这里） =========
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"

cand_obj <- c(file.path(res_dir, "stepH_obj_labeled_celltype7_celltype6.rds"))
obj_fp <- cand_obj[file.exists(cand_obj)][1]
stopifnot(length(obj_fp) == 1)

out_dir <- file.path(res_dir, "NO12_WholeCell_ADvsControl_perDonor")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ========= 2) 读入对象 + 最小自检 =========
obj <- readRDS(obj_fp)
DefaultAssay(obj) <- "RNA"
message("Loaded: ", basename(obj_fp), " | cells=", ncol(obj))

md <- obj@meta.data

## 自动识别：组别列 + donor列（你这里是 group / sample）
group_candidates <- c("group", "group2", "group4", "diagnosis", "clinical_diagnosis")
donor_candidates <- c("sample", "autopsy_id", "donor", "individual", "subject", "subj", "case", "patient_id")

group_col <- group_candidates[group_candidates %in% colnames(md)][1]
donor_col <- donor_candidates[donor_candidates %in% colnames(md)][1]
if (is.na(group_col)) stop("❌ meta.data 找不到组别列（AD/Control）。候选：", paste(group_candidates, collapse=", "))
if (is.na(donor_col)) stop("❌ meta.data 找不到 donor 列。候选：", paste(donor_candidates, collapse=", "))

message("✅ group_col = ", group_col)
message("✅ donor_col = ", donor_col)

## 统一字段名（贴近 syn 模版）
obj$group2 <- as.character(md[[group_col]])
obj$autopsy_id <- as.character(md[[donor_col]])

## 去掉 NA
keep <- !is.na(obj$group2) & !is.na(obj$autopsy_id)
obj <- subset(obj, cells = colnames(obj)[keep])

## 只保留 AD / Control
obj$group2 <- ifelse(obj$group2 %in% c("AD","Control"), obj$group2, NA_character_)
obj <- subset(obj, subset = !is.na(group2))

## 顺序：统计 Control reference；作图 AD 左、Control 右
group_order_plot  <- c("AD","Control")
group_order_stats <- c("Control","AD")
obj$group2 <- factor(as.character(obj$group2), levels = group_order_stats)

## donor 数核对
donor_tab <- obj@meta.data %>%
  distinct(autopsy_id, group2) %>%
  count(group2) %>%
  tidyr::complete(group2 = factor(group_order_stats, levels = group_order_stats), fill = list(n=0))
print(donor_tab)

## ========= 3) Seurat v5/Assay5：从多 layers 的 counts.* 拼回每细胞 CP10k / log1p =========
## 核心思想：
##  - RNA 是 Assay5，counts/data 不在 slot，而在 layers
##  - 每个 counts.* layer 覆盖一部分 cells（列名就是 cell barcodes）
##  - 我们遍历所有 counts.* layers：
##      对每个 layer：
##        - 取列总和作为该 layer 细胞的 library size
##        - 取 UBL3 这一行作为该 layer 细胞的 UBL3 counts
##      然后“按 cell name 对齐”拼成全长向量（长度 = ncol(obj)）
##  - 完全不需要 JoinLayers（避开你之前 FUN locked 报错）

rna <- obj[["RNA"]]
stopifnot(inherits(rna, "Assay5"))

all_layers <- SeuratObject::Layers(rna)

## 3.1 找所有 counts.* layers（严格走 counts 流程）
count_layers <- grep("^counts\\.", all_layers, value = TRUE)

if (length(count_layers) == 0) {
  stop("❌ 没找到任何 counts.* layers：无法按 counts->CP10k 流程计算。")
}
message("✅ counts layers = ", length(count_layers))

## 3.2 确定 UBL3 的行名（SYMBOL 或 ENSG）
## 注意：不同 layer 的 features 应该一致，我们用第一个 layer 检查 rownames
m0 <- SeuratObject::LayerData(rna, layer = count_layers[1])
target <- if ("UBL3" %in% rownames(m0)) "UBL3" else
  if ("ENSG00000122042" %in% rownames(m0)) "ENSG00000122042" else NULL
stopifnot(!is.null(target))
message("✅ UBL3 row = ", target)

## 3.3 初始化“全长向量”（以 obj 当前细胞顺序为准）
cells_all <- colnames(obj)
lib_size <- setNames(rep(NA_real_, length(cells_all)), cells_all)
ubl3_counts <- setNames(rep(NA_real_, length(cells_all)), cells_all)

## 3.4 遍历每个 counts layer，填充对应细胞的 lib_size 和 ubl3_counts
for (ly in count_layers) {
  m <- SeuratObject::LayerData(rna, layer = ly)   # dgCMatrix（genes x cells-of-layer）
  cn <- colnames(m)
  if (length(cn) == 0) next
  
  ## 该 layer 每个细胞的总 counts
  lib <- Matrix::colSums(m)
  ## 该 layer 每个细胞的 UBL3 counts
  uct <- as.numeric(m[target, ])
  
  ## 写入全局向量（按 cell name 对齐）
  lib_size[cn] <- as.numeric(lib)
  ubl3_counts[cn] <- uct
}

## 3.5 兜底：如果某些细胞没被任何 counts layer 覆盖，填 0（一般不该发生）
lib_size[is.na(lib_size)] <- 0
ubl3_counts[is.na(ubl3_counts)] <- 0

## 3.6 计算 CP10k 与 log1p(CP10k)
cp10k <- (ubl3_counts / pmax(lib_size, 1)) * 10000
ubl3_log1p_cp10k <- log1p(cp10k)

## 记录来源：严格 counts multi-layer
expr_source <- "counts_layers"

## ========= 4) donor 层面汇总：每个 donor=1点 =========
df_cells <- data.frame(
  autopsy_id = as.character(obj$autopsy_id),
  group2     = as.character(obj$group2),
  ubl3_cp10k = as.numeric(cp10k),
  ubl3_log   = as.numeric(ubl3_log1p_cp10k),
  stringsAsFactors = FALSE
) %>% filter(!is.na(autopsy_id), !is.na(group2))

df_donor <- df_cells %>%
  group_by(autopsy_id, group2) %>%
  summarise(
    n_cells = dplyr::n(),
    mean_cp10k = mean(ubl3_cp10k, na.rm = TRUE),          # ✅ log2FC 用
    mean_log1p_cp10k = mean(ubl3_log, na.rm = TRUE),      # ✅ Wilcoxon + 图用
    .groups = "drop"
  )

## 作图顺序：AD 左，Control 右
df_donor$group_plot <- factor(df_donor$group2, levels = group_order_plot)

write.csv(df_donor,
          file.path(out_dir, "UBL3_wholecell_perDonor_CP10k_mean_ADvsControl.csv"),
          row.names = FALSE)

message("✅ donor points = ", nrow(df_donor))
print(table(df_donor$group2))

## ========= 5) 统计：AD vs Control 的 Wilcoxon + log2FC =========
d1 <- df_donor %>% filter(group2 == "AD")
d0 <- df_donor %>% filter(group2 == "Control")

mu1 <- mean(d1$mean_cp10k, na.rm=TRUE)
mu0 <- mean(d0$mean_cp10k, na.rm=TRUE)
log2FC <- log2((mu1 + 1e-8) / (mu0 + 1e-8))
pval <- if (nrow(d1) >= 2 && nrow(d0) >= 2) {
  wilcox.test(d1$mean_log1p_cp10k, d0$mean_log1p_cp10k, exact = FALSE)$p.value
} else NA_real_

stats_1 <- data.frame(
  contrast="AD_vs_Control",
  n_AD=nrow(d1),
  n_Control=nrow(d0),
  mean_cp10k_AD=mu1,
  mean_cp10k_Control=mu0,
  log2FC=log2FC,
  p_wilcox=pval,
  expr_source=expr_source,
  used_object=basename(obj_fp),
  stringsAsFactors = FALSE
)
write.csv(stats_1, file.path(out_dir, "UBL3_wholecell_stats_ADvsControl.csv"), row.names = FALSE)
print(stats_1)

fmt_p <- function(p) if (is.na(p)) "NA" else formatC(p, format="e", digits=2)
sub1 <- sprintf("AD vs Control:  log2FC=%s, p=%s",
                if (is.na(log2FC)) "NA" else sprintf("%.3f", log2FC),
                fmt_p(pval))

## ========= 6) 画图：颜色/标题/布局与你模式图一致 =========
pal2_plot <- c(
  AD      = "#D24B40",
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

p1 <- ggplot(df_donor, aes(x = group_plot, y = mean_log1p_cp10k, fill = group_plot)) +
  geom_boxplot(width=0.55, outlier.shape=NA, linewidth=1.0,
               alpha=0.96, colour="grey15", median.linewidth=1.6) +
  geom_point(position = position_jitter(width=0.10, height=0),
             size=2.6, alpha=0.9, shape=21, stroke=0.5, colour="grey10") +
  scale_fill_manual(values = pal2_plot) +
  labs(
    title = "UBL3 expression per donor (Whole cells)",
    subtitle = sub1,
    x = NULL,
    y = "Mean UBL3 log1p(CP10k)"
  ) +
  theme_sci +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10)))

## ========= 7) 保存 PNG（稳定写图） =========
fig1_fp <- file.path(out_dir, "Fig_NO12_UBL3_wholecell_perDonor_log1pCP10k_boxplot_ADvsControl.png")

ragg::agg_png(fig1_fp, width=7.8, height=5.6, units="in", res=450, background="white")
print(p1)
dev.off()

message("✅ saved: ", fig1_fp)
message("DONE. out_dir = ", out_dir)
###############################################################################
