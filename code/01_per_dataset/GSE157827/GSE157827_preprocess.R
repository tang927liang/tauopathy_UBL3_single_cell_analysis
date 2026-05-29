###############################################################################
# GSE157827_preprocess.R
#
# Dataset: GSE157827, Alzheimer’s disease, middle frontal gyrus
# Purpose: raw-to-StepH pipeline with optional verified StepF checkpoint mode
#
# Main output:
#   stepH_obj_celltype6_named.rds
#
# Original environment used for the submitted result:
#   R 4.5.1
#   Seurat 5.3.0
#   SeuratObject 5.1.99.9000
#   Matrix 1.7-4
#   future 1.67.0
#
# Important:
#   Exact clustering is version-sensitive. For exact reproduction of the
#   submitted StepH object, set USE_VERIFIED_STEPF_CHECKPOINT <- TRUE and place
#   the verified 43-cluster StepF object in the results directory.
#   If the checkpoint is absent, the script reruns the full raw-to-StepH pipeline.
###############################################################################

rm(list = ls())
gc()
Sys.setenv(LANG = "en")
options(stringsAsFactors = FALSE)

SEED <- 20251023
set.seed(SEED)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(plyr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

raw_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/Raw_data"
res_dir <- "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Reproducibility mode
# ---------------------------------------------------------------------------
# TRUE  = if the verified original StepF checkpoint exists, start from it and
#         regenerate StepH annotation/output. This is recommended for exact
#         reproduction of the submitted object.
# FALSE = always rerun raw-to-StepF processing before StepH annotation.
USE_VERIFIED_STEPF_CHECKPOINT <- TRUE

verified_stepF <- file.path(res_dir, "stepF_cluster_umap_v2.rds")

sink(file.path(res_dir, "step0_sessionInfo.txt"))
cat("SEED =", SEED, "\n")
print(sessionInfo())
sink()

get_counts_safe <- function(seu, assay = "RNA") {
  out <- tryCatch(
    GetAssayData(seu, assay = assay, layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(out)) {
    out <- GetAssayData(seu, assay = assay, slot = "counts")
  }
  return(out)
}

if (USE_VERIFIED_STEPF_CHECKPOINT && file.exists(verified_stepF)) {
  message("Using verified StepF checkpoint: ", verified_stepF)
  obj_int <- readRDS(verified_stepF)
} else {

###############################################################################
# 1. Read raw 10x matrices and merge all samples
###############################################################################

read_one_sample <- function(prefix) {
  f_mtx <- file.path(raw_dir, paste0(prefix, "_matrix.mtx.gz"))
  f_fea <- file.path(raw_dir, paste0(prefix, "_features.tsv.gz"))
  f_bar <- file.path(raw_dir, paste0(prefix, "_barcodes.tsv.gz"))
  
  mat <- as(readMM(f_mtx), "dgCMatrix")
  fea <- fread(f_fea, header = FALSE)
  bar <- fread(f_bar, header = FALSE)
  
  stopifnot(nrow(fea) == nrow(mat))
  stopifnot(nrow(bar) == ncol(mat))
  
  rownames(mat) <- fea$V1
  colnames(mat) <- bar$V1
  
  sample_id <- sub("^.*?_", "", prefix)
  
  so <- CreateSeuratObject(
    counts = mat,
    project = "GSE157827",
    min.cells = 0,
    min.features = 0
  )
  
  so$gsm <- sub("_.*$", "", prefix)
  so$sample <- sample_id
  so$group <- ifelse(grepl("^AD", sample_id), "AD", "Control")
  colnames(so) <- paste0(sample_id, "_", colnames(so))
  
  so
}

mtx_files <- list.files(raw_dir, pattern = "_matrix\\.mtx(\\.gz)?$", full.names = TRUE)
prefixes <- sort(sub("_matrix\\.mtx(\\.gz)?$", "", basename(mtx_files)))

objs <- lapply(prefixes, function(x) {
  message("Reading: ", x)
  read_one_sample(x)
})
names(objs) <- prefixes

obj <- Reduce(function(x, y) merge(x, y = y), objs)

saveRDS(obj, file.path(res_dir, "GSE157827_merged_ENSG_raw_seurat.rds"))

write.csv(as.data.frame(table(obj$sample)),
          file.path(res_dir, "sample_cell_counts.csv"),
          row.names = FALSE)

###############################################################################
# 2. Rebuild single-layer counts and perform QC
###############################################################################

DefaultAssay(obj) <- "RNA"

mat <- get_counts_safe(obj, assay = "RNA")

obj$nFeature_RNA <- Matrix::colSums(mat > 0)
obj$nCount_RNA <- Matrix::colSums(mat)

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

obj$percent.mt <- 100 * Matrix::colSums(mat[rownames(mat) %in% mt_genes_ensg, , drop = FALSE]) /
  Matrix::colSums(mat)

writeLines(
  sprintf("Before filter: genes=%d, cells=%d", nrow(obj), ncol(obj)),
  file.path(res_dir, "stepB_QC_sizes.txt")
)

png(file.path(res_dir, "stepB_QC_violin_before_filter.png"),
    width = 1600, height = 600)
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0, raster = TRUE))
dev.off()

keep <- obj$nFeature_RNA > 200 &
  obj$nCount_RNA < 20000 &
  obj$percent.mt < 20

obj_flt <- subset(obj, cells = colnames(obj)[keep])

writeLines(
  sprintf("After filter: genes=%d, cells=%d", nrow(obj_flt), ncol(obj_flt)),
  file.path(res_dir, "stepC_QC_sizes.txt")
)

png(file.path(res_dir, "stepC_QC_violin_after_filter.png"),
    width = 1600, height = 600)
print(VlnPlot(obj_flt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3, pt.size = 0, raster = TRUE))
dev.off()

saveRDS(obj_flt, file.path(res_dir, "stepC_filtered_obj.rds"))

write.csv(as.data.frame(table(obj_flt$sample)),
          file.path(res_dir, "stepC_cells_by_sample_after.csv"),
          row.names = FALSE)
write.csv(as.data.frame(table(obj_flt$group)),
          file.path(res_dir, "stepC_cells_by_group_after.csv"),
          row.names = FALSE)

###############################################################################
# 3. Per-sample normalization and CCA integration
###############################################################################

obj_list <- SplitObject(obj_flt, split.by = "sample")

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

g_ad1 <- c("AD1","AD2","AD4","AD5","AD6","AD8","AD9")
g_ad2 <- c("AD10","AD13","AD19","AD20","AD21")
g_nc  <- c("NC3","NC7","NC11","NC12","NC14","NC15","NC16","NC17","NC18")

integrate_by_cca <- function(xlist, outfile) {
  set.seed(SEED)
  feats <- SelectIntegrationFeatures(object.list = xlist, nfeatures = 2000)
  
  set.seed(SEED)
  anchors <- FindIntegrationAnchors(
    object.list = xlist,
    anchor.features = feats,
    dims = 1:20,
    reduction = "cca",
    verbose = TRUE
  )
  
  set.seed(SEED)
  obj_int <- IntegrateData(anchorset = anchors, dims = 1:20)
  
  saveRDS(obj_int, outfile, version = 2, compress = "xz")
  obj_int
}

ad1 <- integrate_by_cca(obj_list[g_ad1],
                        file.path(res_dir, "stepE_integrated_AD_part1.rds"))
ad2 <- integrate_by_cca(obj_list[g_ad2],
                        file.path(res_dir, "stepE_integrated_AD_part2.rds"))
nc  <- integrate_by_cca(obj_list[g_nc],
                        file.path(res_dir, "stepE_integrated_NC.rds"))

rm(obj_list)
gc()

set.seed(SEED)
feats2 <- SelectIntegrationFeatures(object.list = list(ad1, ad2, nc), nfeatures = 2000)

set.seed(SEED)
anchors2 <- FindIntegrationAnchors(
  object.list = list(ad1, ad2, nc),
  anchor.features = feats2,
  dims = 1:20,
  reduction = "cca",
  verbose = TRUE
)

set.seed(SEED)
obj_int <- IntegrateData(anchorset = anchors2, dims = 1:20)

saveRDS(obj_int, file.path(res_dir, "stepE_integrated_ALL_cca.rds"),
        version = 2, compress = "xz")

###############################################################################
# 4. PCA, clustering, and UMAP
###############################################################################

DefaultAssay(obj_int) <- "integrated"

set.seed(SEED)
obj_int <- ScaleData(obj_int, verbose = FALSE)

set.seed(SEED)
obj_int <- RunPCA(obj_int, npcs = 50, verbose = FALSE)

set.seed(SEED)
obj_int <- JackStraw(obj_int, reduction = "pca", dims = 50,
                     num.replicate = 100, verbose = FALSE)
obj_int <- ScoreJackStraw(obj_int, dims = 1:50)

png(file.path(res_dir, "stepF_JackStraw_PC1_50.png"),
    width = 1600, height = 900)
print(JackStrawPlot(obj_int, dims = 1:50))
dev.off()

png(file.path(res_dir, "stepF_ElbowPlot_50PCs.png"),
    width = 1200, height = 800)
print(ElbowPlot(obj_int, ndims = 50))
dev.off()

set.seed(SEED)
obj_int <- FindNeighbors(obj_int, dims = 1:20, k.param = 20, verbose = FALSE)

set.seed(SEED)
obj_int <- FindClusters(obj_int, resolution = 1, algorithm = 1,
                        n.start = 10, n.iter = 10, verbose = FALSE)

set.seed(SEED)
obj_int <- RunUMAP(
  obj_int,
  reduction = "pca",
  dims = 1:20,
  umap.method = "uwot",
  metric = "cosine",
  n.neighbors = 30,
  min.dist = 0.3,
  spread = 1,
  init = "spectral",
  n.components = 2,
  verbose = FALSE
)

png(file.path(res_dir, "stepF_umap_clusters_res1_v3style.png"),
    width = 2000, height = 1200)
print(DimPlot(obj_int, reduction = "umap", label = TRUE,
              label.size = 4, pt.size = 0.2) + NoLegend())
dev.off()

saveRDS(obj_int, file.path(res_dir, "stepF_cluster_umap_v2.rds"),
        version = 2, compress = "xz")
}

###############################################################################
# 5. Marker-based cell-type annotation
###############################################################################

obj <- obj_int
DefaultAssay(obj) <- "integrated"
Idents(obj) <- "seurat_clusters"

markers_all <- FindAllMarkers(
  obj,
  only.pos = FALSE,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  verbose = FALSE
)

markers_all$pass_0.1 <- markers_all$p_val_adj < 0.1

write.csv(markers_all,
          file.path(res_dir, "stepG_FindAllMarkers_wilcox_logfc0.25.csv"),
          row.names = FALSE)

mat_en <- LayerData(obj[["integrated"]], layer = "data")
ensg <- rownames(mat_en)

map_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensg,
  keytype = "ENSEMBL",
  columns = "SYMBOL"
)

sym_by_row <- map_df$SYMBOL[match(ensg, map_df$ENSEMBL)]
keep_gene <- !is.na(sym_by_row)

mat_en <- mat_en[keep_gene, , drop = FALSE]
sym_by_row <- sym_by_row[keep_gene]

sym_levels <- sort(unique(sym_by_row))
row_index <- match(sym_by_row, sym_levels)

G <- sparseMatrix(
  i = row_index,
  j = seq_along(row_index),
  x = 1,
  dims = c(length(sym_levels), length(row_index))
)

mat_sym <- G %*% mat_en
rownames(mat_sym) <- sym_levels

markers_ref <- list(
  Astro   = c("AQP4","GFAP","ALDH1L1","SLC1A3","ADGRV1","GPC5","RYR3"),
  Endo    = c("CLDN5","KDR","FLT1","PECAM1","ABCB1","EBF1"),
  Excit   = c("CAMK2A","SLC17A7","TBR1","CBLN2","LDB2"),
  Inhib   = c("GAD1","GAD2","SLC6A1","LHFPL3","PCDH15"),
  Microgl = c("C3","CX3CR1","P2RY12","AIF1","DOCK8","LRMDA"),
  Oligo   = c("MBP","MOG","PLP1","MOBP","ST18"),
  Peri    = c("PDGFRB","RGS5","MCAM","ACTA2")
)

clu <- Idents(obj)
clu_levels <- levels(clu)

score_one_group <- function(genes) {
  g <- intersect(genes, rownames(mat_sym))
  if (length(g) == 0) {
    return(setNames(rep(NA_real_, length(clu_levels)), clu_levels))
  }
  per_cell <- Matrix::colMeans(mat_sym[g, , drop = FALSE])
  tapply(per_cell, INDEX = clu, FUN = mean, na.rm = TRUE)[clu_levels]
}

avg_by_cluster <- sapply(markers_ref, score_one_group)
tmp <- avg_by_cluster
tmp[is.na(tmp)] <- -Inf

lab_7 <- colnames(tmp)[max.col(tmp, ties.method = "first")]
names(lab_7) <- rownames(tmp)

lab_6 <- lab_7
lab_6[lab_6 == "Peri"] <- "Endo"

obj$celltype7 <- plyr::mapvalues(Idents(obj), from = names(lab_7), to = unname(lab_7))
obj$celltype6 <- plyr::mapvalues(Idents(obj), from = names(lab_6), to = unname(lab_6))

df_scores <- data.frame(
  cluster = rownames(avg_by_cluster),
  label_7 = lab_7,
  label_6 = lab_6,
  avg_by_cluster,
  check.names = FALSE
)

write.csv(df_scores,
          file.path(res_dir, "stepH_cluster_scores_label7_label6.csv"),
          row.names = FALSE)

saveRDS(obj, file.path(res_dir, "stepH_obj_labeled_celltype7_celltype6.rds"),
        version = 2, compress = "xz")

###############################################################################
# 6. Convert abbreviated labels to harmonized six-cell-type names
###############################################################################

nice_levels <- c(
  "Astrocytes",
  "Endothelial",
  "Excitatory neurons",
  "Inhibitory neurons",
  "Microglia",
  "Oligodendrocytes"
)

short2nice <- c(
  Astro   = "Astrocytes",
  Endo    = "Endothelial",
  Excit   = "Excitatory neurons",
  Inhib   = "Inhibitory neurons",
  Microgl = "Microglia",
  Oligo   = "Oligodendrocytes"
)

map_df <- read.csv(file.path(res_dir, "stepH_cluster_scores_label7_label6.csv"),
                   stringsAsFactors = FALSE)

lab6_map <- setNames(unname(short2nice[map_df$label_6]), map_df$cluster)

obj$celltype6 <- plyr::mapvalues(
  Idents(obj),
  from = names(lab6_map),
  to = as.character(lab6_map)
)

obj$celltype6 <- factor(obj$celltype6, levels = nice_levels)

write.csv(
  data.frame(
    celltype = names(table(obj$celltype6)),
    n_cells = as.integer(table(obj$celltype6)),
    percent = round(100 * prop.table(table(obj$celltype6)), 2)
  ),
  file.path(res_dir, "stepH_celltype6_counts_percent.csv"),
  row.names = FALSE
)

pal6 <- c(
  Astrocytes = "#E76F51",
  Endothelial = "#2A9D8F",
  `Excitatory neurons` = "#457B9D",
  `Inhibitory neurons` = "#F4A261",
  Microglia = "#8D99AE",
  Oligodendrocytes = "#6D597A"
)

png(file.path(res_dir, "stepH_umap_celltype6_named.png"),
    width = 2000, height = 1400, res = 180)
print(
  DimPlot(obj, reduction = "umap", group.by = "celltype6",
          label = TRUE, label.size = 5, repel = TRUE) +
    scale_color_manual(values = pal6, drop = FALSE) +
    labs(color = "Cell type")
)
dev.off()

saveRDS(obj, file.path(res_dir, "stepH_obj_celltype6_named.rds"),
        version = 2, compress = "xz")

cat("DONE: saved stepH_obj_celltype6_named.rds\n")
cat("Cells:", ncol(obj), "\n")
cat("Clusters:", length(unique(obj$seurat_clusters)), "\n")
print(table(obj$celltype6))
