###############################################################################
# 01_per_dataset / syn21788402  (Leng et al. 2021; entorhinal cortex EC +
# superior frontal gyrus SFG)
# Canonical creation script: builds the harmonized 6-cell-type Seurat objects
# end-to-end from the authors' public SingleCellExperiment objects. Consumed by
# code/02_integrated_figures/.
#
# Source  : Synapse syn21788402 (Leng et al. 2021) --
#           sce.EC/SFG.scAlign.assigned.rds (counts + scran logcounts +
#           colData[clusterCellType, clusterAssignment, SampleID, BraakStage] +
#           reducedDims[CCA, CCA.ALIGNED]); QC / alignment / clustering /
#           per-cell annotation were all performed by the authors.
# Method  : build a Seurat object from counts (authors' scran logcounts used as
#           the data layer, not re-normalized); attach CCA.ALIGNED and RunUMAP
#           (seed.use = 42, matching the authors' @commands); map the authors'
#           clusterCellType (7 classes) to the 6 harmonized classes, with
#           OPC -> Oligodendrocytes (oligodendrocyte lineage). clusterCellType is
#           retained for traceability.
# Outputs : resultsmodify/stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds
#           resultsmodify/stepH_syn21788402_SFG_obj_celltype6.rds
#           (these resultsmodify/ objects are the canonical inputs read by the
#            integrated figure/stat scripts.)
#
# METHOD NOTE vs the other datasets: GSE157827 / GSE174367 / syn52082747 assign
#   the 6 cell types by de-novo cluster-marker scoring; syn21788402 instead reuses
#   the authors' published annotation (clusterCellType) mapped to the same 6
#   classes, because this dataset is distributed already annotated.
# ENV NOTE: the `.libPaths()/clean_lib` line below points to a local clean Seurat
#   5.5.0 library; adjust or remove it for your environment. Built/tested with
#   Seurat 5.5.0 / SeuratObject 5.4.0.
#
# Paths point to the data-project location on disk; source data and the resulting
# objects are not stored in this repository. Confirm paths before running.
###############################################################################
# syn21788402 (Leng et al. 2021) — stepH 生成流水线
# 修正版：适配 Seurat 5.5.0 / SeuratObject 5.4.0
###############################################################################

# =========================
# 0. 使用干净 Seurat library
# =========================
clean_lib <- "D:/Rlibs/R45_seurat_clean"
.libPaths(c(clean_lib, .libPaths()))

SEED <- 42
set.seed(SEED)

# 尽量保持经典 v3 Assay 风格
options(Seurat.object.assay.version = "v3")

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
  library(Seurat)
  library(SingleCellExperiment)
  library(plyr)
})

cat("当前 R library 路径:\n")
print(.libPaths())

cat("\n当前关键包版本:\n")
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("SeuratObject:", as.character(packageVersion("SeuratObject")), "\n")
cat("Matrix:", as.character(packageVersion("Matrix")), "\n")
cat("SingleCellExperiment:", as.character(packageVersion("SingleCellExperiment")), "\n")

# =========================
# 1. 路径设置
# =========================
in_dir  <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/rawdata"
out_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ec_file  <- file.path(in_dir, "sce.EC.scAlign.assigned.rds")
sfg_file <- file.path(in_dir, "sce.SFG.scAlign.assigned.rds")

cat("\n检查输入文件:\n")
cat("EC file exists:", file.exists(ec_file), "\n")
cat("SFG file exists:", file.exists(sfg_file), "\n")

if (!file.exists(ec_file)) {
  stop("找不到 EC 输入文件: ", ec_file)
}
if (!file.exists(sfg_file)) {
  stop("找不到 SFG 输入文件: ", sfg_file)
}

# =========================
# 2. 统一 6 大类命名 + 映射
# =========================
celltype6_levels_std <- c(
  "Astrocytes",
  "Excitatory neurons",
  "Microglia",
  "Endothelial",
  "Inhibitory neurons",
  "Oligodendrocytes"
)

short_levels_std <- c(
  "Astro",
  "Excit",
  "Microgl",
  "Endo",
  "Inhib",
  "Oligo"
)

# clusterCellType（作者 7 类）→ celltype6_short
# OPC 归入 Oligodendrocytes
clusterCellType_to_short <- c(
  Astro = "Astro",
  Exc   = "Excit",
  Micro = "Microgl",
  Endo  = "Endo",
  Inh   = "Inhib",
  OPC   = "Oligo",
  Oligo = "Oligo"
)

short2nice <- c(
  Astro   = "Astrocytes",
  Excit   = "Excitatory neurons",
  Microgl = "Microglia",
  Endo    = "Endothelial",
  Inhib   = "Inhibitory neurons",
  Oligo   = "Oligodendrocytes"
)

###############################################################################
# 3. 主函数：从 SCE 重建 Seurat stepH 对象
###############################################################################
process_region <- function(region, sce_file, out_file) {

  message("\n============================================================")
  message("Processing region: ", region)
  message("Input file: ", sce_file)
  message("Output file: ", out_file)
  message("============================================================")

  stopifnot(file.exists(sce_file))

  # ---- 读取 SCE ----
  sce <- readRDS(sce_file)

  # ---- 检查 SCE 内容 ----
  if (!"counts" %in% assayNames(sce)) {
    stop(region, ": SCE 中没有 assay 'counts'")
  }
  if (!"logcounts" %in% assayNames(sce)) {
    stop(region, ": SCE 中没有 assay 'logcounts'")
  }
  if (!"CCA.ALIGNED" %in% reducedDimNames(sce)) {
    stop(region, ": SCE 中没有 reducedDim 'CCA.ALIGNED'")
  }
  if (!"clusterCellType" %in% colnames(colData(sce))) {
    stop(region, ": colData 中没有 clusterCellType")
  }
  if (!"clusterAssignment" %in% colnames(colData(sce))) {
    stop(region, ": colData 中没有 clusterAssignment")
  }
  if (!"SampleID" %in% colnames(colData(sce))) {
    stop(region, ": colData 中没有 SampleID")
  }

  counts_mat    <- assay(sce, "counts")
  logcounts_mat <- assay(sce, "logcounts")
  meta_df       <- as.data.frame(colData(sce))
  emb_aligned   <- as.matrix(reducedDim(sce, "CCA.ALIGNED"))

  # ---- 基础一致性检查 ----
  if (!identical(colnames(counts_mat), rownames(meta_df))) {
    stop(region, ": counts_mat 的 colnames 与 meta_df 的 rownames 不一致")
  }

  if (!identical(colnames(counts_mat), colnames(logcounts_mat))) {
    stop(region, ": counts_mat 与 logcounts_mat 的细胞名不一致")
  }

  if (!identical(rownames(counts_mat), rownames(logcounts_mat))) {
    stop(region, ": counts_mat 与 logcounts_mat 的基因名不一致")
  }

  message("genes x cells: ", nrow(counts_mat), " x ", ncol(counts_mat))
  message("CCA.ALIGNED dims: ", ncol(emb_aligned))

  # =========================
  # 3.1 处理 gene names
  # =========================
  # Seurat 不允许 feature names 中含有 "_"
  # 如果不提前处理，CreateSeuratObject 会自动把 "_" 改成 "-"
  # 但 logcounts_mat 仍然保留原名，后续 SetAssayData 会不一致
  gene_names_clean <- gsub("_", "-", rownames(counts_mat))

  if (anyDuplicated(gene_names_clean)) {
    message(region, ": '_' 替换为 '-' 后产生重复 gene names，使用 make.unique() 处理")
    gene_names_clean <- make.unique(gene_names_clean, sep = ".dup")
  }

  rownames(counts_mat)    <- gene_names_clean
  rownames(logcounts_mat) <- gene_names_clean

  # 再次确认
  if (!identical(rownames(counts_mat), rownames(logcounts_mat))) {
    stop(region, ": gene names 清理后 counts 与 logcounts 仍然不一致")
  }

  # =========================
  # 3.2 创建 Seurat 对象
  # =========================
  obj <- CreateSeuratObject(
    counts = counts_mat,
    project = paste0("syn21788402_", region),
    min.cells = 0,
    min.features = 0
  )

  DefaultAssay(obj) <- "RNA"

  # data 层直接使用 SCE 的 scran logcounts
  # SeuratObject >= 5.0 必须使用 layer = "data"，不能再用 slot = "data"
  obj <- SetAssayData(
    object = obj,
    assay = "RNA",
    layer = "data",
    new.data = logcounts_mat
  )

  # =========================
  # 3.3 添加 metadata
  # =========================
  obj <- AddMetaData(obj, metadata = meta_df)
  obj$orig.ident <- factor(as.character(obj$SampleID))

  # =========================
  # 3.4 添加 CCA.ALIGNED reduction 并重跑 UMAP
  # =========================
  rownames(emb_aligned) <- colnames(obj)
  colnames(emb_aligned) <- paste0("CCAALN_", seq_len(ncol(emb_aligned)))

  obj[["CCA.ALIGNED"]] <- CreateDimReducObject(
    embeddings = emb_aligned,
    key = "CCAALN_",
    assay = "RNA"
  )

  obj <- RunUMAP(
    object = obj,
    reduction = "CCA.ALIGNED",
    dims = 1:ncol(emb_aligned),
    seed.use = SEED,
    verbose = FALSE
  )

  # =========================
  # 3.5 cluster_use + Idents
  # =========================
  obj$cluster_use <- obj$clusterAssignment
  Idents(obj) <- obj$clusterAssignment

  # =========================
  # 3.6 映射 6 大类
  # =========================
  cct <- as.character(obj$clusterCellType)
  short <- unname(clusterCellType_to_short[cct])

  if (any(is.na(short))) {
    bad <- unique(cct[is.na(short)])
    stop(region, ": 以下 clusterCellType 没有映射关系: ", paste(bad, collapse = ", "))
  }

  obj$celltype6_short <- factor(short, levels = short_levels_std)
  obj$celltype6 <- factor(
    unname(short2nice[as.character(obj$celltype6_short)]),
    levels = celltype6_levels_std
  )

  # =========================
  # 3.7 核查
  # =========================
  cat("\nclusterCellType x celltype6:\n")
  print(table(obj$clusterCellType, obj$celltype6))

  cat("\ncelltype6 counts:\n")
  print(table(obj$celltype6))

  stopifnot(sum(table(obj$celltype6)) == ncol(obj))

  # OPC 应归入 Oligodendrocytes
  stopifnot(
    sum(obj$clusterCellType == "OPC") ==
      sum(obj$celltype6 == "Oligodendrocytes") -
      sum(obj$clusterCellType == "Oligo")
  )

  # =========================
  # 3.8 保存
  # =========================
  saveRDS(obj, out_file)

  message("\n✅ Saved: ", out_file)
  message("Object cells: ", ncol(obj))
  message("Object genes: ", nrow(obj))

  # 保存后立即测试是否可以重新读入
  test_obj <- readRDS(out_file)
  message("✅ Reload test passed: ", out_file)
  message("Reloaded cells: ", ncol(test_obj))
  message("Reloaded genes: ", nrow(test_obj))
  rm(test_obj)

  invisible(obj)
}

###############################################################################
# 4. 分别处理 EC 和 SFG
###############################################################################

obj_EC <- process_region(
  region = "EC",
  sce_file = ec_file,
  out_file = file.path(
    out_dir,
    "stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds"
  )
)

obj_SFG <- process_region(
  region = "SFG",
  sce_file = sfg_file,
  out_file = file.path(
    out_dir,
    "stepH_syn21788402_SFG_obj_celltype6.rds"
  )
)

message("\n============================================================")
message("完成。两个 stepH 文件已写入:")
message(out_dir)
message("============================================================")

# =========================
# 5. 保存 sessionInfo
# =========================
sink(file.path(out_dir, "sessionInfo_stepH_reconstruction.txt"))
print(sessionInfo())
sink()

message("✅ sessionInfo saved: ", file.path(out_dir, "sessionInfo_stepH_reconstruction.txt"))

# =========================
# 6. 最终检查
# =========================
cat("\n最终对象检查:\n")
cat("EC cells:", ncol(obj_EC), " genes:", nrow(obj_EC), "\n")
cat("SFG cells:", ncol(obj_SFG), " genes:", nrow(obj_SFG), "\n")

cat("\nEC celltype6:\n")
print(table(obj_EC$celltype6))

cat("\nSFG celltype6:\n")
print(table(obj_SFG$celltype6))

cat("\n全部完成。\n")
