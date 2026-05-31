###############################################################################
# Reproduce GSE174367 upstream object:
#   stepH_obj_celltype6_named.rds
#
# This one-file script is the compact GitHub version of the original long
# "#GSE174367 全程代码.docx" workflow.
#
# Required raw files from GEO:
#   GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5
#   GSE174367_snRNA-seq_cell_meta.csv.gz
#
# Default input folder for reviewers:
#   ./data/GSE174367
# Default output folder:
#   ./results/GSE174367_stepH
#
# You can override paths without editing the script:
#   Sys.setenv(GSE174367_RAW_DIR = "D:/path/to/GSE174367")
#   Sys.setenv(GSE174367_OUT_DIR = "D:/path/to/output")
#
# For a fast local rerun from an existing StepF object:
#   Sys.setenv(GSE174367_STEPF_RDS = "D:/codex/174367/stepF_afterPCA_graph_umap_res0.8_dims1_20_GSE174367.rds")
#
# Historical environment:
#   R 4.4.3
#   Seurat 4.3.0
#   SeuratObject 5.2.0
#   Matrix 1.7-4
#   data.table 1.17.8
###############################################################################

SEED <- 20251023
set.seed(SEED)

RAW_DIR <- Sys.getenv("GSE174367_RAW_DIR", unset = file.path(getwd(), "data", "GSE174367"))
OUT_DIR <- Sys.getenv("GSE174367_OUT_DIR", unset = file.path(getwd(), "results", "GSE174367_stepH"))
STEPF_RDS <- Sys.getenv("GSE174367_STEPF_RDS", unset = "")
RUN_FIND_MARKERS <- tolower(Sys.getenv("GSE174367_RUN_FIND_MARKERS", unset = "true")) %in%
  c("1", "true", "yes", "y")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)

load_required_packages <- function(extra = character()) {
  pkgs <- unique(c(
    "Seurat",
    "SeuratObject",
    "Matrix",
    "data.table",
    "future",
    "plyr",
    "ggplot2",
    extra
  ))
  suppressPackageStartupMessages({
    for (pkg in pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Required R package is missing: ", pkg, call. = FALSE)
      }
      library(pkg, character.only = TRUE)
    }
  })
}

load_required_packages(extra = c("cluster"))

log_env <- function(step_name) {
  log_file <- file.path(OUT_DIR, paste0(step_name, "_sessionInfo.txt"))
  set.seed(SEED)
  sink(log_file)
  cat("SEED =", SEED, "\n")
  cat("RAW_DIR =", RAW_DIR, "\n")
  cat("OUT_DIR =", OUT_DIR, "\n")
  print(sessionInfo())
  sink()
  message("Saved sessionInfo: ", log_file)
}

get_assay_data_compat <- function(object, assay = NULL, slot = "data", layer = slot) {
  if (is.null(assay)) assay <- DefaultAssay(object)
  out <- tryCatch(
    GetAssayData(object, assay = assay, slot = slot),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)
  GetAssayData(object, assay = assay, layer = layer)
}

save_rds_v2_xz <- function(object, filename) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  saveRDS(object, filename, version = 2, compress = "xz")
  message("Saved: ", filename)
}

celltype6_levels <- c(
  "Astrocytes",
  "Endothelial",
  "Excitatory neurons",
  "Inhibitory neurons",
  "Microglia",
  "Oligodendrocytes"
)

celltype6_palette <- c(
  "Astrocytes" = "#FF8E8E",
  "Endothelial" = "#B8A109",
  "Excitatory neurons" = "#09BB3C",
  "Inhibitory neurons" = "#00BFC4",
  "Microglia" = "#36B0E1",
  "Oligodendrocytes" = "#E16AFC"
)

short_to_nice <- c(
  Astro = "Astrocytes",
  Endo = "Endothelial",
  Excit = "Excitatory neurons",
  Inhib = "Inhibitory neurons",
  Microgl = "Microglia",
  Oligo = "Oligodendrocytes"
)

# Reviewed cluster-level labels from the original manuscript workflow.
# celltype7 is retained only as an audit trail. The final manuscript class is
# celltype6, where Peri and Endo are merged into Endothelial.
cluster_label_map <- data.frame(
  cluster = as.character(0:25),
  label_7 = c(
    "Oligo", "Oligo", "Oligo", "Oligo", "Microgl", "Astro",
    "Inhib", "Astro", "Oligo", "Excit", "Oligo", "Inhib",
    "Inhib", "Microgl", "Oligo", "Excit", "Inhib", "Inhib",
    "Oligo", "Peri", "Inhib", "Inhib", "Endo", "Oligo",
    "Excit", "Inhib"
  ),
  label_6 = c(
    "Oligo", "Oligo", "Oligo", "Oligo", "Microgl", "Astro",
    "Inhib", "Astro", "Oligo", "Excit", "Oligo", "Inhib",
    "Inhib", "Microgl", "Oligo", "Excit", "Inhib", "Inhib",
    "Oligo", "Endo", "Inhib", "Inhib", "Endo", "Oligo",
    "Excit", "Inhib"
  ),
  stringsAsFactors = FALSE
)

###############################################################################
# Raw data to StepF
###############################################################################

run_raw_to_stepf <- function() {
  log_env("step0_setup")

  h5_file <- file.path(RAW_DIR, "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5")
  meta_file <- file.path(RAW_DIR, "GSE174367_snRNA-seq_cell_meta.csv.gz")
  if (!file.exists(h5_file)) stop("Missing h5 file: ", h5_file, call. = FALSE)
  if (!file.exists(meta_file)) stop("Missing metadata file: ", meta_file, call. = FALSE)

  message("Step 1: read raw h5 and metadata")
  mat <- Read10X_h5(h5_file)
  obj <- CreateSeuratObject(
    counts = mat,
    project = "GSE174367",
    min.cells = 0,
    min.features = 0
  )

  meta <- data.table::fread(meta_file)
  stopifnot(all(meta$Barcode %in% colnames(obj)))
  common_cells <- intersect(colnames(obj), meta$Barcode)
  obj <- subset(obj, cells = common_cells)

  meta2 <- meta[match(colnames(obj), meta$Barcode), ]
  stopifnot(all(meta2$Barcode == colnames(obj)))

  new_cell_names <- paste0(meta2$SampleID, "_", colnames(obj))
  stopifnot(!anyDuplicated(new_cell_names))
  colnames(obj) <- new_cell_names

  meta2$Barcode <- NULL
  rownames(meta2) <- colnames(obj)
  obj <- AddMetaData(obj, metadata = as.data.frame(meta2))
  obj$group <- factor(obj$Diagnosis, levels = c("Control", "AD"))

  saveRDS(obj, file.path(OUT_DIR, "GSE174367_merged_raw_with_meta.rds"))
  saveRDS(obj, file.path(OUT_DIR, "GSE174367_merged_raw_with_meta_v2.rds"), version = 2)
  log_env("step1_merge")

  message("Step 2: QC filtering")
  DefaultAssay(obj) <- "RNA"
  counts <- get_assay_data_compat(obj, assay = "RNA", slot = "counts")
  obj$nFeature_RNA <- Matrix::colSums(counts > 0)
  obj$nCount_RNA <- Matrix::colSums(counts)
  is_mt <- grepl("^MT-", rownames(counts))
  obj$percent.mt <- Matrix::colSums(counts[is_mt, , drop = FALSE]) /
    Matrix::colSums(counts) * 100

  png(file.path(OUT_DIR, "figures", "stepB_QC_violin_before_filter.png"),
      width = 1600, height = 600)
  print(VlnPlot(
    obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0,
    raster = TRUE
  ))
  dev.off()

  keep <- (obj$nFeature_RNA > 200) &
    (obj$nCount_RNA < 20000) &
    (obj$percent.mt < 20)
  obj_flt <- subset(obj, cells = colnames(obj)[keep])

  png(file.path(OUT_DIR, "figures", "stepC_QC_violin_after_filter.png"),
      width = 1600, height = 600)
  print(VlnPlot(
    obj_flt,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0,
    raster = TRUE
  ))
  dev.off()

  saveRDS(obj_flt, file.path(OUT_DIR, "stepC_filtered_obj.rds"))
  writeLines(
    sprintf("Before filter: genes=%d, cells=%d", nrow(obj), ncol(obj)),
    file.path(OUT_DIR, "stepB_QC_sizes.txt")
  )
  writeLines(
    sprintf("After filter: genes=%d, cells=%d", nrow(obj_flt), ncol(obj_flt)),
    file.path(OUT_DIR, "stepC_QC_sizes.txt")
  )
  write.csv(as.data.frame(table(obj$SampleID)),
            file.path(OUT_DIR, "tables", "stepC_cells_by_sample_before.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(table(obj_flt$SampleID)),
            file.path(OUT_DIR, "tables", "stepC_cells_by_sample_after.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(table(obj$group)),
            file.path(OUT_DIR, "tables", "stepC_cells_by_group_before.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(table(obj_flt$group)),
            file.path(OUT_DIR, "tables", "stepC_cells_by_group_after.csv"),
            row.names = FALSE)

  rm(obj, mat, counts)
  gc()
  log_env("step2_QC")

  message("Step D: split by SampleID, LogNormalize, HVG(vst, 1000)")
  DefaultAssay(obj_flt) <- "RNA"
  obj_list <- SplitObject(obj_flt, split.by = "SampleID")
  for (nm in names(obj_list)) {
    message("  sample: ", nm)
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
  saveRDS(obj_list, file.path(OUT_DIR, "stepD_split_normalized_vst1000_GSE174367.rds"))
  rm(obj_flt)
  gc()
  log_env("stepD_split_normalize_HVG")

  message("Step E: hierarchical CCA integration")
  g_ad1 <- c("Sample-17", "Sample-19", "Sample-22", "Sample-27", "Sample-33")
  g_ad2 <- c("Sample-37", "Sample-43", "Sample-45", "Sample-46", "Sample-47", "Sample-50")
  g_ctrl <- c("Sample-52", "Sample-58", "Sample-66", "Sample-82", "Sample-90", "Sample-96", "Sample-100")

  integrate_by_cca <- function(xlist, outfile, seed = SEED) {
    set.seed(seed)
    feats <- SelectIntegrationFeatures(object.list = xlist, nfeatures = 2000)
    set.seed(seed)
    anchors <- FindIntegrationAnchors(
      object.list = xlist,
      anchor.features = feats,
      dims = 1:20,
      reduction = "cca",
      verbose = TRUE
    )
    set.seed(seed)
    obj_int <- IntegrateData(anchorset = anchors, dims = 1:20)
    saveRDS(obj_int, outfile)
    obj_int
  }

  obj_ad1_int <- integrate_by_cca(
    obj_list[g_ad1],
    file.path(OUT_DIR, "stepE_integrated_AD_part1_GSE174367.rds")
  )
  obj_ad2_int <- integrate_by_cca(
    obj_list[g_ad2],
    file.path(OUT_DIR, "stepE_integrated_AD_part2_GSE174367.rds")
  )
  obj_ctrl_int <- integrate_by_cca(
    obj_list[g_ctrl],
    file.path(OUT_DIR, "stepE_integrated_CTRL_GSE174367.rds")
  )
  rm(obj_list)
  gc()
  log_env("stepE_groupwise_cca")

  set.seed(SEED)
  feats2 <- SelectIntegrationFeatures(
    object.list = list(obj_ad1_int, obj_ad2_int, obj_ctrl_int),
    nfeatures = 2000
  )
  set.seed(SEED)
  anchors2 <- FindIntegrationAnchors(
    object.list = list(obj_ad1_int, obj_ad2_int, obj_ctrl_int),
    anchor.features = feats2,
    dims = 1:20,
    reduction = "cca",
    verbose = TRUE
  )
  set.seed(SEED)
  obj <- IntegrateData(anchorset = anchors2, dims = 1:20)
  save_rds_v2_xz(obj, file.path(OUT_DIR, "stepE_integrated_ALL_cca_GSE174367.rds"))
  rm(obj_ad1_int, obj_ad2_int, obj_ctrl_int, anchors2)
  gc()
  log_env("stepE_all_cca")

  message("Step F: PCA, JackStraw, Louvain clustering and UMAP")
  DefaultAssay(obj) <- "integrated"
  obj <- ScaleData(obj, verbose = FALSE)
  set.seed(SEED)
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

  future::plan(sequential)
  options(
    future.globals.maxSize = Inf,
    future.rng.onMisuse = "ignore"
  )
  set.seed(SEED)
  obj <- JackStraw(
    obj,
    reduction = "pca",
    dims = 50,
    num.replicate = 50,
    verbose = FALSE
  )
  obj <- ScoreJackStraw(obj, dims = 1:50)
  save_rds_v2_xz(obj, file.path(OUT_DIR, "stepF_afterPCA_JackStraw_GSE174367.rds"))

  png(file.path(OUT_DIR, "figures", "stepF_JackStraw_PC_pvalues_GSE174367.png"),
      width = 1600, height = 900)
  print(JackStrawPlot(obj, dims = 1:50))
  dev.off()

  png(file.path(OUT_DIR, "figures", "stepF_ElbowPlot_50PCs_GSE174367.png"),
      width = 1200, height = 800)
  print(ElbowPlot(obj, ndims = 50))
  dev.off()

  dims_use <- 1:20
  set.seed(SEED)
  obj <- FindNeighbors(obj, dims = dims_use, k.param = 20, verbose = FALSE)
  set.seed(SEED)
  obj <- FindClusters(
    obj,
    resolution = 0.8,
    algorithm = 1,
    n.start = 10,
    n.iter = 10,
    verbose = FALSE
  )
  set.seed(SEED)
  obj <- RunUMAP(
    obj,
    reduction = "pca",
    dims = dims_use,
    umap.method = "uwot",
    metric = "cosine",
    n.neighbors = 30,
    min.dist = 0.3,
    spread = 1,
    init = "spectral",
    n.components = 2,
    verbose = FALSE
  )

  png(file.path(OUT_DIR, "figures", "stepF_UMAP_clusters_res0.8_dims1_20_GSE174367.png"),
      width = 2000, height = 1200)
  print(DimPlot(
    obj,
    reduction = "umap",
    label = TRUE,
    label.size = 8,
    pt.size = 0.2
  ) + NoLegend())
  dev.off()

  ari <- NA_real_
  if (requireNamespace("mclust", quietly = TRUE)) {
    ref_labels <- as.integer(factor(obj$Cell.Type))
    clu_labels <- as.integer(obj$seurat_clusters)
    ari <- mclust::adjustedRandIndex(ref_labels, clu_labels)
  } else {
    message("Package mclust is not installed; ARI metric will be recorded as NA.")
  }
  emb <- Embeddings(obj, "pca")[, dims_use]
  sil <- cluster::silhouette(as.integer(obj$seurat_clusters), dist(emb))
  sil_avg <- summary(sil)$avg.width
  writeLines(
    c(
      "dims_use = 1:20",
      "resolution = 0.8",
      paste0("ARI = ", ari),
      paste0("silhouette_avg = ", sil_avg)
    ),
    con = file.path(OUT_DIR, "tables", "stepF_metrics_dims1_20_res0.8_GSE174367.txt")
  )

  save_rds_v2_xz(obj, file.path(OUT_DIR, "stepF_afterPCA_graph_umap_res0.8_dims1_20_GSE174367.rds"))
  log_env("stepF_final_done")
  obj
}

###############################################################################
# StepG/StepH annotation
###############################################################################

annotate_stepH <- function(obj) {
  DefaultAssay(obj) <- "integrated"
  Idents(obj) <- "seurat_clusters"

  message("Cluster counts:")
  print(table(Idents(obj)))

  if (RUN_FIND_MARKERS) {
    message("Step G: FindAllMarkers")
    future::plan("sequential")
    options(future.globals.maxSize = 8 * 1024^3)
    markers_all <- FindAllMarkers(
      obj,
      only.pos = FALSE,
      test.use = "wilcox",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      verbose = TRUE
    )

    if (nrow(markers_all) == 0) stop("FindAllMarkers returned zero rows.", call. = FALSE)
    if (!"avg_log2FC" %in% colnames(markers_all) && "avg_logFC" %in% colnames(markers_all)) {
      markers_all$avg_log2FC <- markers_all$avg_logFC
    }
    markers_all$pass_0.1 <- markers_all$p_val_adj < 0.1
    write.csv(
      markers_all,
      file.path(OUT_DIR, "tables", "stepG_FindAllMarkers_wilcox_logfc0.25_GSE174367.csv"),
      row.names = FALSE
    )

    top20 <- subset(markers_all, avg_log2FC > 0 & pass_0.1)
    top20 <- top20[order(top20$cluster, -top20$avg_log2FC), ]
    top20 <- do.call(rbind, by(top20, top20$cluster, head, n = 20))
    feat20 <- unique(top20$gene)

    png(file.path(OUT_DIR, "figures", "stepG_DotPlot_top20_per_cluster_GSE174367.png"),
        width = 2200, height = 1400, res = 180)
    print(DotPlot(obj, features = feat20) + RotatedAxis())
    dev.off()

    png(file.path(OUT_DIR, "figures", "stepG_Heatmap_top20_per_cluster_GSE174367.png"),
        width = 2200, height = 1400, res = 180)
    print(DoHeatmap(obj, features = feat20, raster = TRUE))
    dev.off()
  }

  observed_clusters <- sort(unique(as.character(obj@meta.data$seurat_clusters)))
  missing_clusters <- setdiff(observed_clusters, cluster_label_map$cluster)
  if (length(missing_clusters) > 0) {
    stop("No reviewed cell-type label for cluster(s): ",
         paste(missing_clusters, collapse = ", "), call. = FALSE)
  }

  label7_by_cluster <- setNames(cluster_label_map$label_7, cluster_label_map$cluster)
  label6_by_cluster <- setNames(cluster_label_map$label_6, cluster_label_map$cluster)
  cluster_vec <- as.character(obj@meta.data$seurat_clusters)

  label_7 <- unname(label7_by_cluster[cluster_vec])
  label_6_short <- unname(label6_by_cluster[cluster_vec])
  label_6_nice <- unname(short_to_nice[label_6_short])

  obj$celltype7 <- factor(
    label_7,
    levels = c("Oligo", "Microgl", "Astro", "Inhib", "Excit", "Peri", "Endo")
  )
  obj$celltype6 <- factor(label_6_nice, levels = celltype6_levels)
  if (!"sample" %in% colnames(obj@meta.data) && "SampleID" %in% colnames(obj@meta.data)) {
    obj$sample <- obj$SampleID
  }

  write.csv(
    cluster_label_map,
    file.path(OUT_DIR, "tables", "stepH_cluster_label_map_GSE174367.csv"),
    row.names = FALSE
  )

  counts6 <- table(obj$celltype6)
  props6 <- round(100 * prop.table(counts6), 1)
  counts6_df <- data.frame(
    celltype = names(counts6),
    n_cells = as.integer(counts6),
    percent = as.numeric(props6)
  )
  write.csv(
    counts6_df,
    file.path(OUT_DIR, "tables", "stepH_celltype6_counts_percent_GSE174367.csv"),
    row.names = FALSE
  )

  counts7 <- table(obj$celltype7)
  props7 <- round(100 * prop.table(counts7), 1)
  write.csv(
    data.frame(
      celltype = names(counts7),
      n_cells = as.integer(counts7),
      percent = as.numeric(props7)
    ),
    file.path(OUT_DIR, "tables", "stepH_celltype7_counts_percent_GSE174367.csv"),
    row.names = FALSE
  )

  png(file.path(OUT_DIR, "figures", "stepH_umap_celltype6_named_GSE174367.png"),
      width = 2000, height = 1400, res = 180)
  print(
    DimPlot(
      obj,
      reduction = "umap",
      group.by = "celltype6",
      label = TRUE,
      label.size = 5,
      repel = TRUE
    ) +
      scale_color_manual(values = celltype6_palette, drop = FALSE) +
      labs(color = "Cell type")
  )
  dev.off()

  if ("sample" %in% colnames(obj@meta.data)) {
    df_bar <- as.data.frame(table(sample = obj$sample, celltype = obj$celltype6))
    df_bar <- within(df_bar, {
      total_by_sample <- ave(Freq, sample, FUN = sum)
      percent <- 100 * Freq / total_by_sample
    })

    png(file.path(OUT_DIR, "figures", "stepH_bar_sample_celltype6_GSE174367.png"),
        width = 1800, height = 1200, res = 180)
    print(
      ggplot(df_bar, aes(x = sample, y = percent, fill = celltype)) +
        geom_bar(stat = "identity", width = 0.9) +
        scale_fill_manual(values = celltype6_palette, drop = FALSE) +
        coord_flip() +
        labs(x = NULL, y = "Proportion (%)", fill = "Cell type") +
        theme_bw(base_size = 12)
    )
    dev.off()
  }

  message("Final celltype6 counts:")
  print(counts6_df)

  DefaultAssay(obj) <- "RNA"
  save_rds_v2_xz(obj, file.path(OUT_DIR, "stepH_obj_celltype6_named.rds"))
  log_env("stepH_celltype6_done")
  obj
}

###############################################################################
# Main
###############################################################################

if (nzchar(STEPF_RDS)) {
  if (!file.exists(STEPF_RDS)) stop("GSE174367_STEPF_RDS does not exist: ", STEPF_RDS, call. = FALSE)
  message("Reading existing StepF object: ", STEPF_RDS)
  obj_stepf <- readRDS(STEPF_RDS)
} else {
  obj_stepf <- run_raw_to_stepf()
}

obj_stepH <- annotate_stepH(obj_stepf)
message("Done. Final object: ", file.path(OUT_DIR, "stepH_obj_celltype6_named.rds"))
