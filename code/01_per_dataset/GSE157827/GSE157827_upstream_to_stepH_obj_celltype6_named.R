# GSE157827 upstream pipeline to stepH_obj_celltype6_named.rds.
#
# This is a single-file, GitHub-ready version of the upstream workflow used to
# regenerate the manuscript object:
#
#   stepH_obj_celltype6_named.rds
#
# It combines the cleaned Step A-H code and the verified cluster-to-cell-type
# labels. Parameters match the manuscript rerun used for downstream validation.
#
# Typical rerun from the merged Seurat object:
#   Rscript GSE157827_upstream_to_stepH_obj_celltype6_named.R
#
# Verified shortcut from an existing StepF object:
#   set START_FROM_STEPF=true and STEPF_RDS=/path/to/stepF_cluster_umap_v2.rds
#
# Path configuration can be edited below or supplied as environment variables.

rm(list = ls())
gc()
Sys.setenv(LANG = "en")
options(stringsAsFactors = FALSE)

###############################################################################
# 0. User configuration
###############################################################################

SEED <- as.integer(Sys.getenv("SEED", unset = "20251023"))
set.seed(SEED)

project_root <- Sys.getenv("PROJECT_ROOT", unset = "D:/codex/157827")
output_dir <- Sys.getenv("OUTPUT_DIR", unset = file.path(project_root, "redo"))
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

merged_rds <- Sys.getenv(
  "MERGED_RDS",
  unset = file.path(project_root, "GSE157827_merged_with_group.rds")
)
stepf_rds <- Sys.getenv(
  "STEPF_RDS",
  unset = file.path(output_dir, "stepF_cluster_umap_v2.rds")
)
if (!file.exists(stepf_rds)) {
  stepf_rds <- file.path(project_root, "stepF_cluster_umap_v2.rds")
}

final_rds <- Sys.getenv(
  "FINAL_RDS",
  unset = file.path(output_dir, "stepH_obj_celltype6_named.rds")
)

env_flag <- function(name, default = FALSE) {
  val <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
  val %in% c("1", "true", "yes", "y")
}

# START_FROM_STEPF=true reproduces StepH directly from a StepF checkpoint. This
# is the exact mode used for the final validation against the historical object.
START_FROM_STEPF <- env_flag("START_FROM_STEPF", default = FALSE)

# Marker discovery/figures do not change the final object. They are optional
# because FindAllMarkers and heatmaps are slow for this object.
RUN_FIND_MARKERS <- env_flag("RUN_FIND_MARKERS", default = FALSE)
MAKE_FIGURES <- env_flag("MAKE_FIGURES", default = TRUE)
SAVE_INTERMEDIATE_RDS <- env_flag("SAVE_INTERMEDIATE_RDS", default = TRUE)
SAVE_QS <- env_flag("SAVE_QS", default = FALSE)
DRY_RUN_NO_SAVE <- env_flag("DRY_RUN_NO_SAVE", default = FALSE)
STOP_ON_EXPECTED_COUNT_MISMATCH <- env_flag(
  "STOP_ON_EXPECTED_COUNT_MISMATCH",
  default = TRUE
)

###############################################################################
# 1. Packages and compatibility helpers
###############################################################################

required_pkgs <- c(
  "Matrix",
  "SeuratObject",
  "Seurat",
  "future",
  "ggplot2"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
  library(Seurat)
  library(future)
  library(ggplot2)
})

if (SAVE_QS && !requireNamespace("qs", quietly = TRUE)) {
  stop("SAVE_QS=true but the qs package is not installed.", call. = FALSE)
}

log_env <- function(out_dir, step_name = "session") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  log_file <- file.path(out_dir, paste0(step_name, "_sessionInfo.txt"))
  sink(log_file)
  cat("SEED =", SEED, "\n")
  print(sessionInfo())
  sink()
  message("Saved sessionInfo: ", log_file)
}

get_assay_data_compat <- function(object, assay = NULL, slot = "data", layer = slot) {
  if (is.null(assay)) {
    assay <- DefaultAssay(object)
  }
  out <- tryCatch(
    GetAssayData(object, assay = assay, layer = layer),
    error = function(e) NULL
  )
  if (!is.null(out)) {
    return(out)
  }
  GetAssayData(object, assay = assay, slot = slot)
}

get_layer_names_compat <- function(assay_obj) {
  out <- tryCatch(Layers(assay_obj), error = function(e) character())
  if (length(out) == 0) {
    out <- "counts"
  }
  out
}

get_layer_data_compat <- function(assay_obj, layer = "counts") {
  out <- tryCatch(LayerData(assay_obj, layer = layer), error = function(e) NULL)
  if (!is.null(out)) {
    return(out)
  }
  GetAssayData(assay_obj, slot = layer)
}

save_rds_v2 <- function(object, filename) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  compress_method <- Sys.getenv("RDS_COMPRESS", unset = "none")
  if (tolower(compress_method) %in% c("false", "none", "no", "0")) {
    compress_method <- FALSE
  }
  tmp_file <- paste0(filename, ".tmp")
  if (file.exists(tmp_file)) {
    unlink(tmp_file)
  }
  saveRDS(object, tmp_file, version = 2, compress = compress_method)
  if (file.exists(filename)) {
    unlink(filename)
  }
  ok <- file.rename(tmp_file, filename)
  if (!ok) {
    stop("Could not move temporary RDS into place: ", filename, call. = FALSE)
  }
  message("Saved: ", filename)
}

save_optional_qs <- function(object, filename) {
  if (!SAVE_QS) {
    return(invisible(FALSE))
  }
  qs::qsave(object, filename)
  message("Saved: ", filename)
  invisible(TRUE)
}

###############################################################################
# 2. Constants from the manuscript pipeline
###############################################################################

mt_genes_ensg <- c(
  "ENSG00000198888", "ENSG00000198727", "ENSG00000198804", "ENSG00000198886",
  "ENSG00000212907", "ENSG00000198786", "ENSG00000198695", "ENSG00000198712",
  "ENSG00000198899", "ENSG00000198938", "ENSG00000198840", "ENSG00000198763",
  "ENSG00000210107", "ENSG00000210112", "ENSG00000210117", "ENSG00000210127",
  "ENSG00000210133", "ENSG00000210140", "ENSG00000210144", "ENSG00000210151",
  "ENSG00000210156", "ENSG00000210160", "ENSG00000210164", "ENSG00000210169",
  "ENSG00000210174", "ENSG00000210179", "ENSG00000210184", "ENSG00000210189",
  "ENSG00000210194", "ENSG00000210199", "ENSG00000210204", "ENSG00000210209",
  "ENSG00000210214", "ENSG00000210219", "ENSG00000228253", "ENSG00000228630",
  "ENSG00000210130"
)

celltype6_levels <- c(
  "Astrocytes",
  "Endothelial",
  "Excitatory neurons",
  "Inhibitory neurons",
  "Microglia",
  "Oligodendrocytes"
)

celltype6_palette <- c(
  "Astrocytes" = "#E76F51",
  "Endothelial" = "#2A9D8F",
  "Excitatory neurons" = "#457B9D",
  "Inhibitory neurons" = "#F4A261",
  "Microglia" = "#8D99AE",
  "Oligodendrocytes" = "#6D597A"
)

short_to_nice <- c(
  Astro = "Astrocytes",
  Endo = "Endothelial",
  Excit = "Excitatory neurons",
  Inhib = "Inhibitory neurons",
  Microgl = "Microglia",
  Oligo = "Oligodendrocytes"
)

cluster_label_map <- data.frame(
  cluster = as.character(0:42),
  celltype = c(
    "Oligo", "Excit", "Oligo", "Astro", "OPC", "Excit", "Excit", "Oligo",
    "Microgl", "Inhib", "Inhib", "Excit", "Astro", "Inhib", "Excit",
    "Excit", "Inhib", "Inhib", "Excit", "Excit", "Excit", "Endo",
    "Excit", "Oligo", "Excit", "Excit", "Inhib", "Excit", "Excit",
    "Oligo", "Excit", "Oligo", "OPC", "Excit", "Oligo", "Astro",
    "Oligo", "Excit", "Astro", "OPC", "OPC", "Inhib", "OPC"
  ),
  celltype7 = c(
    "Oligo", "Excit", "Oligo", "Astro", "Inhib", "Excit", "Excit",
    "Oligo", "Microgl", "Inhib", "Inhib", "Excit", "Astro", "Inhib",
    "Excit", "Excit", "Inhib", "Inhib", "Excit", "Excit", "Excit",
    "Endo", "Excit", "Oligo", "Excit", "Excit", "Inhib", "Excit",
    "Astro", "Oligo", "Excit", "Oligo", "Excit", "Microgl", "Oligo",
    "Astro", "Oligo", "Excit", "Astro", "Inhib", "Astro", "Inhib",
    "Inhib"
  ),
  celltype6_short = c(
    "Oligo", "Excit", "Oligo", "Astro", "Inhib", "Excit", "Excit",
    "Oligo", "Microgl", "Inhib", "Inhib", "Excit", "Astro", "Inhib",
    "Excit", "Excit", "Inhib", "Inhib", "Excit", "Excit", "Excit",
    "Endo", "Excit", "Oligo", "Excit", "Excit", "Inhib", "Excit",
    "Astro", "Oligo", "Excit", "Oligo", "Excit", "Microgl", "Oligo",
    "Astro", "Oligo", "Excit", "Astro", "Inhib", "Astro", "Inhib",
    "Inhib"
  ),
  stringsAsFactors = FALSE
)
cluster_label_map$celltype6 <- unname(short_to_nice[cluster_label_map$celltype6_short])

expected_celltype6_counts <- c(
  "Astrocytes" = 19157L,
  "Endothelial" = 2420L,
  "Excitatory neurons" = 60095L,
  "Inhibitory neurons" = 34821L,
  "Microglia" = 7562L,
  "Oligodendrocytes" = 45451L
)

###############################################################################
# 3. Step A-F: preprocessing, integration, PCA, clustering and UMAP
###############################################################################

if (!START_FROM_STEPF) {
  if (!file.exists(merged_rds)) {
    stop(
      "Missing merged input object: ", merged_rds,
      "\nSet MERGED_RDS or use START_FROM_STEPF=true with STEPF_RDS.",
      call. = FALSE
    )
  }

  log_env(output_dir, "step0_setup")

  message("Reading merged object: ", merged_rds)
  obj_old <- readRDS(merged_rds)
  DefaultAssay(obj_old) <- "RNA"

  # Step A. Rebuild single-layer RNA counts and restore sample/group metadata.
  rna_assay <- obj_old[["RNA"]]
  layers <- get_layer_names_compat(rna_assay)
  mats <- list()
  barcodes_all <- character()

  for (lyr in layers) {
    mat <- get_layer_data_compat(rna_assay, layer = lyr)
    if (ncol(mat) == 0) {
      next
    }
    mats[[lyr]] <- mat
    barcodes_all <- c(barcodes_all, colnames(mat))
    message(sprintf("layer: %s | genes=%d cells=%d", lyr, nrow(mat), ncol(mat)))
  }

  if (length(mats) == 0) {
    stop("No non-empty RNA count layers were found.", call. = FALSE)
  }

  genes1 <- rownames(mats[[1]])
  same_gene_order <- vapply(mats, function(x) identical(rownames(x), genes1), logical(1))
  if (!all(same_gene_order)) {
    stop("RNA layers have different gene order; cannot cbind safely.", call. = FALSE)
  }

  counts_combined <- do.call(cbind, mats)
  stopifnot(ncol(counts_combined) == length(barcodes_all))
  stopifnot(identical(colnames(counts_combined), barcodes_all))

  obj <- CreateSeuratObject(
    counts = counts_combined,
    project = "GSE157827",
    min.cells = 0,
    min.features = 0
  )

  samples <- sub("_.*$", "", colnames(obj))
  groups <- ifelse(grepl("^AD", samples), "AD", "Control")
  obj$sample <- samples
  obj$group <- factor(groups, levels = c("Control", "AD"))

  write.csv(as.data.frame(table(obj$sample)),
            file.path(output_dir, "step0_cells_by_sample_beforeQC.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(table(obj$group)),
            file.path(output_dir, "step0_cells_by_group_beforeQC.csv"),
            row.names = FALSE)

  if (SAVE_INTERMEDIATE_RDS) {
    save_rds_v2(obj, file.path(output_dir, "stepA_rebuilt_single_counts.rds"))
  }
  rm(obj_old, rna_assay, mats, counts_combined)
  gc()

  # Step B/C. QC metrics and manuscript filtering thresholds.
  DefaultAssay(obj) <- "RNA"
  mat <- get_assay_data_compat(obj, assay = "RNA", slot = "counts", layer = "counts")

  obj$nFeature_RNA <- Matrix::colSums(mat > 0)
  obj$nCount_RNA <- Matrix::colSums(mat)

  is_mt <- rownames(mat) %in% mt_genes_ensg
  mt_counts <- Matrix::colSums(mat[is_mt, , drop = FALSE])
  total_counts <- Matrix::colSums(mat)
  obj$percent.mt <- ifelse(total_counts > 0, 100 * mt_counts / total_counts, 0)

  if (MAKE_FIGURES) {
    png(file.path(figure_dir, "stepB_QC_violin_before_filter.png"),
        width = 1600, height = 600)
    print(VlnPlot(
      obj,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      ncol = 3,
      pt.size = 0,
      raster = TRUE
    ))
    dev.off()
  }

  writeLines(
    sprintf("Before filter: genes=%d, cells=%d", nrow(obj), ncol(obj)),
    file.path(output_dir, "stepB_QC_sizes.txt")
  )

  keep <- (obj$nFeature_RNA > 200) &
    (obj$nCount_RNA < 20000) &
    (obj$percent.mt < 20)

  obj_flt <- subset(obj, cells = colnames(obj)[keep])

  if (MAKE_FIGURES) {
    png(file.path(figure_dir, "stepC_QC_violin_after_filter.png"),
        width = 1600, height = 600)
    print(VlnPlot(
      obj_flt,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      ncol = 3,
      pt.size = 0,
      raster = TRUE
    ))
    dev.off()
  }

  if (SAVE_INTERMEDIATE_RDS) {
    save_rds_v2(obj_flt, file.path(output_dir, "stepC_filtered_obj.rds"))
  }

  write.csv(as.data.frame(table(obj$sample, useNA = "ifany")),
            file.path(output_dir, "stepC_cells_by_sample_before.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(table(obj_flt$sample, useNA = "ifany")),
            file.path(output_dir, "stepC_cells_by_sample_after.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(table(obj$group, useNA = "ifany")),
            file.path(output_dir, "stepC_cells_by_group_before.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(table(obj_flt$group, useNA = "ifany")),
            file.path(output_dir, "stepC_cells_by_group_after.csv"),
            row.names = FALSE)
  writeLines(
    sprintf("After filter: genes=%d, cells=%d", nrow(obj_flt), ncol(obj_flt)),
    file.path(output_dir, "stepC_QC_sizes.txt")
  )

  rm(obj, mat)
  gc()
  log_env(output_dir, "stepC_QC")

  # Step D. Split by sample, LogNormalize, HVG(vst, 1000).
  DefaultAssay(obj_flt) <- "RNA"
  obj_list <- SplitObject(obj_flt, split.by = "sample")

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

  if (SAVE_INTERMEDIATE_RDS) {
    save_rds_v2(obj_list, file.path(output_dir, "stepD_split_normalized_vst1000_clean.rds"))
  }
  rm(obj_flt)
  gc()
  log_env(output_dir, "stepD_split_normalize_HVG")

  # Step E. Hierarchical CCA integration.
  g_ad1 <- c("AD1", "AD2", "AD4", "AD5", "AD6", "AD8", "AD9")
  g_ad2 <- c("AD10", "AD13", "AD19", "AD20", "AD21")
  g_nc <- c("NC3", "NC7", "NC11", "NC12", "NC14", "NC15", "NC16", "NC17", "NC18")

  expected_samples <- c(g_ad1, g_ad2, g_nc)
  missing_samples <- setdiff(expected_samples, names(obj_list))
  if (length(missing_samples) > 0) {
    stop("Missing sample(s) from split object: ",
         paste(missing_samples, collapse = ", "), call. = FALSE)
  }

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
    if (SAVE_INTERMEDIATE_RDS) {
      save_rds_v2(obj_int, outfile)
    }
    obj_int
  }

  obj_ad1_int <- integrate_by_cca(
    obj_list[g_ad1],
    file.path(output_dir, "stepE_integrated_AD_part1.rds")
  )
  obj_ad2_int <- integrate_by_cca(
    obj_list[g_ad2],
    file.path(output_dir, "stepE_integrated_AD_part2.rds")
  )
  obj_nc_int <- integrate_by_cca(
    obj_list[g_nc],
    file.path(output_dir, "stepE_integrated_NC.rds")
  )

  rm(obj_list)
  gc()

  set.seed(SEED)
  feats2 <- SelectIntegrationFeatures(
    object.list = list(obj_ad1_int, obj_ad2_int, obj_nc_int),
    nfeatures = 2000
  )
  set.seed(SEED)
  anchors2 <- FindIntegrationAnchors(
    object.list = list(obj_ad1_int, obj_ad2_int, obj_nc_int),
    anchor.features = feats2,
    dims = 1:20,
    reduction = "cca",
    verbose = TRUE
  )
  set.seed(SEED)
  obj <- IntegrateData(anchorset = anchors2, dims = 1:20)

  if (SAVE_INTERMEDIATE_RDS) {
    save_rds_v2(obj, file.path(output_dir, "stepE_integrated_ALL_cca_v2.rds"))
    save_optional_qs(obj, file.path(output_dir, "stepE_integrated_ALL_cca.qs"))
  }

  rm(obj_ad1_int, obj_ad2_int, obj_nc_int, anchors2)
  gc()
  log_env(output_dir, "stepE_all_cca")

  # Step F. ScaleData, PCA(50), JackStraw, Louvain clustering and UMAP.
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
    num.replicate = 100,
    verbose = FALSE
  )
  obj <- ScoreJackStraw(obj, dims = 1:50)

  if (SAVE_INTERMEDIATE_RDS) {
    save_rds_v2(obj, file.path(output_dir, "stepF_after_PCA_JackStraw_v2.rds"))
    save_optional_qs(obj, file.path(output_dir, "stepF_after_PCA_JackStraw.qs"))
  }

  if (MAKE_FIGURES) {
    png(file.path(figure_dir, "stepF_JackStraw_PC1_50.png"),
        width = 1600, height = 900)
    print(JackStrawPlot(obj, dims = 1:50))
    dev.off()

    png(file.path(figure_dir, "stepF_ElbowPlot_50PCs.png"),
        width = 1200, height = 800)
    print(ElbowPlot(obj, ndims = 50))
    dev.off()
  }

  dims_use <- 1:20

  set.seed(SEED)
  obj <- FindNeighbors(obj, dims = dims_use, k.param = 20, verbose = FALSE)

  set.seed(SEED)
  obj <- FindClusters(
    obj,
    resolution = 1,
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

  if (MAKE_FIGURES) {
    png(file.path(figure_dir, "stepF_umap_clusters_res1_v3style.png"),
        width = 2000, height = 1200)
    print(DimPlot(
      obj,
      reduction = "umap",
      label = TRUE,
      label.size = 4,
      pt.size = 0.2
    ) + NoLegend())
    dev.off()
  }

  save_rds_v2(obj, file.path(output_dir, "stepF_cluster_umap_v2.rds"))
  save_optional_qs(obj, file.path(output_dir, "stepF_cluster_umap.qs"))
  stepf_rds <- file.path(output_dir, "stepF_cluster_umap_v2.rds")
  log_env(output_dir, "stepF_final_done")
} else {
  if (!file.exists(stepf_rds)) {
    stop("START_FROM_STEPF=true but StepF input is missing: ", stepf_rds, call. = FALSE)
  }
  message("START_FROM_STEPF=true; using StepF object: ", stepf_rds)
}

###############################################################################
# 4. Step G/H: marker review and final six-cell-type annotation
###############################################################################

if (!exists("obj")) {
  message("Reading StepF object: ", stepf_rds)
  obj <- readRDS(stepf_rds)
}
stopifnot("seurat_clusters" %in% colnames(obj@meta.data))

DefaultAssay(obj) <- "integrated"
Idents(obj) <- "seurat_clusters"

message("Cluster counts:")
print(table(Idents(obj)))

if (RUN_FIND_MARKERS) {
  future::plan("sequential")
  options(future.globals.maxSize = Inf)

  markers_all <- FindAllMarkers(
    obj,
    only.pos = FALSE,
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    verbose = TRUE
  )

  if (nrow(markers_all) == 0) {
    stop("FindAllMarkers returned zero rows.", call. = FALSE)
  }

  if (!"avg_log2FC" %in% colnames(markers_all) && "avg_logFC" %in% colnames(markers_all)) {
    markers_all$avg_log2FC <- markers_all$avg_logFC
  }
  if (!"avg_log2FC" %in% colnames(markers_all)) {
    stop("Marker table lacks avg_log2FC/avg_logFC.", call. = FALSE)
  }

  markers_all$pass_0.1 <- markers_all$p_val_adj < 0.1
  marker_csv <- file.path(table_dir, "stepG_FindAllMarkers_wilcox_logfc0.25.csv")
  write.csv(markers_all, marker_csv, row.names = FALSE)
  message("Saved marker table: ", marker_csv)

  if (MAKE_FIGURES) {
    top20 <- subset(markers_all, avg_log2FC > 0 & pass_0.1)
    top20 <- top20[order(top20$cluster, -top20$avg_log2FC), ]
    top20 <- do.call(rbind, by(top20, top20$cluster, head, n = 20))
    feat20 <- unique(top20$gene)

    png(file.path(figure_dir, "stepG_DotPlot_top20_per_cluster.png"),
        width = 2200, height = 1400, res = 180)
    print(DotPlot(obj, features = feat20) + RotatedAxis())
    dev.off()

    png(file.path(figure_dir, "stepG_Heatmap_top20_per_cluster.png"),
        width = 2200, height = 1400, res = 180)
    print(DoHeatmap(obj, features = feat20, raster = TRUE))
    dev.off()
  }
} else {
  message("RUN_FIND_MARKERS=false; marker review is skipped because it does not change StepH.")
}

observed_clusters <- sort(unique(as.character(obj@meta.data$seurat_clusters)))
missing_clusters <- setdiff(observed_clusters, cluster_label_map$cluster)
if (length(missing_clusters) > 0) {
  stop("No cell-type label for cluster(s): ",
       paste(missing_clusters, collapse = ", "), call. = FALSE)
}

celltype_by_cluster <- setNames(cluster_label_map$celltype, cluster_label_map$cluster)
celltype7_by_cluster <- setNames(cluster_label_map$celltype7, cluster_label_map$cluster)
celltype6_by_cluster <- setNames(cluster_label_map$celltype6, cluster_label_map$cluster)

cluster_vec <- as.character(obj@meta.data$seurat_clusters)
obj$celltype <- unname(celltype_by_cluster[cluster_vec])
obj$celltype7 <- factor(
  unname(celltype7_by_cluster[cluster_vec]),
  levels = c("Oligo", "Excit", "Astro", "Inhib", "Microgl", "Endo")
)
obj$celltype6 <- factor(
  unname(celltype6_by_cluster[cluster_vec]),
  levels = celltype6_levels
)

write.csv(
  cluster_label_map,
  file.path(table_dir, "stepH_cluster_label_map_GSE157827.csv"),
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
  file.path(table_dir, "stepH_celltype6_counts_percent.csv"),
  row.names = FALSE
)

if (MAKE_FIGURES) {
  png(file.path(figure_dir, "stepH_umap_celltype6_named.png"),
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
}

if (MAKE_FIGURES && "sample" %in% colnames(obj@meta.data)) {
  df_bar <- as.data.frame(table(sample = obj$sample, celltype = obj$celltype6))
  df_bar <- within(df_bar, {
    total_by_sample <- ave(Freq, sample, FUN = sum)
    percent <- 100 * Freq / total_by_sample
  })

  png(file.path(figure_dir, "stepH_bar_sample_celltype6.png"),
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

observed_counts <- setNames(as.integer(counts6), names(counts6))
expected_same <- identical(observed_counts[names(expected_celltype6_counts)], expected_celltype6_counts)
if (!isTRUE(expected_same)) {
  msg <- paste0(
    "Final celltype6 counts differ from the validated manuscript object.\n",
    "Observed:\n",
    paste(capture.output(print(observed_counts)), collapse = "\n"),
    "\nExpected:\n",
    paste(capture.output(print(expected_celltype6_counts)), collapse = "\n")
  )
  if (STOP_ON_EXPECTED_COUNT_MISMATCH) {
    stop(msg, call. = FALSE)
  } else {
    warning(msg, call. = FALSE)
  }
}

if (DRY_RUN_NO_SAVE) {
  message("DRY_RUN_NO_SAVE=true; final RDS was not written.")
} else {
  save_rds_v2(obj, final_rds)
}
log_env(output_dir, "stepH_celltype6_done")

cat("\nCompleted GSE157827 upstream pipeline.\n")
cat("Final object:", final_rds, "\n")
