###############################################################################
# Rebuild syn21788402 SFG stepH Seurat object from raw SCE
#
# Output:
#   stepH_syn21788402_SFG_obj_celltype6.rds
#
# This is a standalone upstream reconstruction script prepared for code
# deposition. It follows the original syn21788402 full-code parameters:
#   - SEED = 42
#   - options(Seurat.object.assay.version = "v3")
#   - CreateSeuratObject(min.cells = 0, min.features = 0)
#   - RNA data layer copied from SCE logcounts
#   - CCA.ALIGNED reduction copied from the SCE object
#   - UMAP rerun on all CCA.ALIGNED dimensions
#   - clusterCellType OPC and Oligo are merged into Oligodendrocytes
###############################################################################

rm(list = ls())
gc()

SEED <- 42
set.seed(SEED)
options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v3")

# ---- Paths -------------------------------------------------------------------
# Override these with environment variables when running outside the original
# workstation layout.
raw_dir <- Sys.getenv(
  "SYN21788402_RAW_DIR",
  "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/rawdata"
)
out_dir <- Sys.getenv(
  "SYN21788402_SFG_OUT_DIR",
  "D:/codex/21788402/redo/SFG"
)
sce_file <- file.path(raw_dir, "sce.SFG.scAlign.assigned.rds")
out_file <- file.path(
  out_dir,
  "stepH_syn21788402_SFG_obj_celltype6.rds"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- R libraries and packages -------------------------------------------------
candidate_libs <- unique(c(
  Sys.getenv("SYN21788402_R_LIB", unset = NA_character_),
  "D:/Rlibs/R45_seurat_clean",
  "D:/R_lib_clean",
  "C:/Users/setou/Documents/R/libs-v4-clean"
))
candidate_libs <- candidate_libs[!is.na(candidate_libs)]
candidate_libs <- candidate_libs[dir.exists(candidate_libs)]
if (length(candidate_libs) > 0) {
  .libPaths(unique(c(candidate_libs, .libPaths())))
}

required_packages <- c("Matrix", "SeuratObject", "Seurat")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nCurrent .libPaths():\n",
    paste(.libPaths(), collapse = "\n")
  )
}

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
  library(Seurat)
})

cat("Current R library paths:\n")
print(.libPaths())
cat("\nKey package versions:\n")
for (pkg in required_packages) {
  cat(pkg, ":", as.character(packageVersion(pkg)), "\n")
}

# ---- Cell-type harmonization --------------------------------------------------
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

clusterCellType_to_short <- c(
  Astro = "Astro",
  Exc = "Excit",
  Micro = "Microgl",
  Endo = "Endo",
  Inh = "Inhib",
  OPC = "Oligo",
  Oligo = "Oligo"
)

short2nice <- c(
  Astro = "Astrocytes",
  Excit = "Excitatory neurons",
  Microgl = "Microglia",
  Endo = "Endothelial",
  Inhib = "Inhibitory neurons",
  Oligo = "Oligodendrocytes"
)

# ---- Helpers -----------------------------------------------------------------
extract_sce_with_api <- function(sce) {
  counts_mat <- SummarizedExperiment::assay(sce, "counts")
  logcounts_mat <- SummarizedExperiment::assay(sce, "logcounts")
  meta_df <- as.data.frame(SummarizedExperiment::colData(sce))
  emb_aligned <- as.matrix(
    SingleCellExperiment::reducedDim(sce, "CCA.ALIGNED")
  )

  list(
    counts = counts_mat,
    logcounts = logcounts_mat,
    meta = meta_df,
    cca_aligned = emb_aligned,
    reduced_dim_names = SingleCellExperiment::reducedDimNames(sce),
    assay_names = SummarizedExperiment::assayNames(sce)
  )
}

extract_sce_from_raw_slots <- function(sce) {
  assays_env <- attr(attr(sce, "assays"), ".xData")
  assay_list <- attr(get("data", assays_env), "listData")
  assay_names <- names(assay_list)
  if (!all(c("counts", "logcounts") %in% assay_names)) {
    stop(
      "SCE raw slots do not contain counts/logcounts assays. Found: ",
      paste(assay_names, collapse = ", ")
    )
  }

  counts_mat <- assay_list[["counts"]]
  logcounts_mat <- assay_list[["logcounts"]]
  if (is.null(rownames(counts_mat)) || is.null(colnames(counts_mat))) {
    stop("counts assay is missing rownames or colnames")
  }
  rownames(logcounts_mat) <- rownames(counts_mat)
  colnames(logcounts_mat) <- colnames(counts_mat)

  col_data <- attr(sce, "colData")
  meta_list <- attr(col_data, "listData")
  meta_df <- as.data.frame(meta_list, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(meta_df) <- attr(col_data, "rownames")

  reduced_list <- attr(attr(sce, "reducedDims"), "listData")
  reduced_dim_names <- names(reduced_list)
  if (!"CCA.ALIGNED" %in% reduced_dim_names) {
    stop(
      "SCE raw slots do not contain reducedDim CCA.ALIGNED. Found: ",
      paste(reduced_dim_names, collapse = ", ")
    )
  }
  emb_aligned <- as.matrix(reduced_list[["CCA.ALIGNED"]])
  rownames(emb_aligned) <- colnames(counts_mat)

  list(
    counts = counts_mat,
    logcounts = logcounts_mat,
    meta = meta_df,
    cca_aligned = emb_aligned,
    reduced_dim_names = reduced_dim_names,
    assay_names = assay_names
  )
}

extract_sce_data <- function(sce) {
  can_use_sce_api <- all(vapply(
    c("SingleCellExperiment", "SummarizedExperiment"),
    requireNamespace,
    logical(1),
    quietly = TRUE
  ))

  if (can_use_sce_api) {
    suppressPackageStartupMessages({
      library(SingleCellExperiment)
      library(SummarizedExperiment)
    })
    return(extract_sce_with_api(sce))
  }

  cat("SingleCellExperiment/SummarizedExperiment unavailable; using raw S4 slot extractor.\n")
  extract_sce_from_raw_slots(sce)
}

safe_save_rds <- function(obj, path, compress = TRUE, version = NULL) {
  tmp <- paste0(path, ".tmp")
  if (file.exists(tmp)) file.remove(tmp)

  if (is.null(version)) {
    saveRDS(obj, tmp, compress = compress)
  } else {
    saveRDS(obj, tmp, compress = compress, version = version)
  }
  test_obj <- readRDS(tmp)
  rm(test_obj)
  gc()

  ok <- file.copy(tmp, path, overwrite = TRUE)
  if (!ok) stop("Failed to move temporary RDS into place: ", path)
  file.remove(tmp)
  invisible(path)
}

write_celltype_qc <- function(obj, region, out_dir) {
  counts <- as.data.frame(table(obj$celltype6), stringsAsFactors = FALSE)
  colnames(counts) <- c("celltype6", "n_cells")
  counts$percent <- counts$n_cells / sum(counts$n_cells) * 100
  write.csv(
    counts,
    file.path(out_dir, paste0("stepH_syn21788402_", region, "_celltype6_counts_percent.csv")),
    row.names = FALSE
  )

  ctab <- as.data.frame(
    table(obj$clusterCellType, obj$celltype6),
    stringsAsFactors = FALSE
  )
  colnames(ctab) <- c("clusterCellType", "celltype6", "n_cells")
  write.csv(
    ctab,
    file.path(out_dir, paste0("stepH_syn21788402_", region, "_clusterCellType_to_celltype6.csv")),
    row.names = FALSE
  )
}

# ---- Build SFG object ---------------------------------------------------------
cat("\n============================================================\n")
cat("Processing syn21788402 SFG\n")
cat("Input file : ", sce_file, "\n", sep = "")
cat("Output file: ", out_file, "\n", sep = "")
cat("============================================================\n")
if (!file.exists(sce_file)) stop("Cannot find SFG input file: ", sce_file)

sce <- readRDS(sce_file)
dat <- extract_sce_data(sce)
rm(sce)
gc()

if (!"counts" %in% dat$assay_names) stop("SFG: SCE has no assay 'counts'")
if (!"logcounts" %in% dat$assay_names) stop("SFG: SCE has no assay 'logcounts'")
if (!"CCA.ALIGNED" %in% dat$reduced_dim_names) {
  stop("SFG: SCE has no reducedDim 'CCA.ALIGNED'")
}

required_meta <- c("clusterCellType", "clusterAssignment", "SampleID")
missing_meta <- setdiff(required_meta, colnames(dat$meta))
if (length(missing_meta) > 0) {
  stop("SFG: missing colData columns: ", paste(missing_meta, collapse = ", "))
}

counts_mat <- dat$counts
logcounts_mat <- dat$logcounts
meta_df <- dat$meta
emb_aligned <- dat$cca_aligned
rm(dat)
gc()

if (!identical(colnames(counts_mat), rownames(meta_df))) {
  stop("SFG: counts colnames and metadata rownames differ")
}
if (!identical(colnames(counts_mat), colnames(logcounts_mat))) {
  stop("SFG: counts and logcounts cell names differ")
}
if (!identical(rownames(counts_mat), rownames(logcounts_mat))) {
  stop("SFG: counts and logcounts gene names differ")
}

cat("genes x cells: ", nrow(counts_mat), " x ", ncol(counts_mat), "\n", sep = "")
cat("CCA.ALIGNED dims: ", ncol(emb_aligned), "\n", sep = "")

gene_names_clean <- gsub("_", "-", rownames(counts_mat))
if (anyDuplicated(gene_names_clean)) {
  cat("SFG: duplicate gene names after '_' -> '-', applying make.unique()\n")
  gene_names_clean <- make.unique(gene_names_clean, sep = ".dup")
}
rownames(counts_mat) <- gene_names_clean
rownames(logcounts_mat) <- gene_names_clean

obj <- CreateSeuratObject(
  counts = counts_mat,
  project = "syn21788402_SFG",
  min.cells = 0,
  min.features = 0
)
rm(counts_mat)
gc()

DefaultAssay(obj) <- "RNA"
obj <- SetAssayData(
  object = obj,
  assay = "RNA",
  layer = "data",
  new.data = logcounts_mat
)
rm(logcounts_mat)
gc()

obj <- AddMetaData(obj, metadata = meta_df)
obj$orig.ident <- factor(as.character(obj$SampleID))
rm(meta_df)
gc()

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
rm(emb_aligned)
gc()

obj$cluster_use <- obj$clusterAssignment
Idents(obj) <- obj$clusterAssignment

cct <- as.character(obj$clusterCellType)
short <- unname(clusterCellType_to_short[cct])
if (any(is.na(short))) {
  bad <- unique(cct[is.na(short)])
  stop("SFG: unmapped clusterCellType values: ", paste(bad, collapse = ", "))
}

obj$celltype6_short <- factor(short, levels = short_levels_std)
obj$celltype6 <- factor(
  unname(short2nice[as.character(obj$celltype6_short)]),
  levels = celltype6_levels_std
)

cat("\nclusterCellType x celltype6 for SFG:\n")
print(table(obj$clusterCellType, obj$celltype6))
cat("\ncelltype6 counts for SFG:\n")
print(table(obj$celltype6))

expected_counts <- c(
  "Astrocytes" = 8025,
  "Excitatory neurons" = 20301,
  "Microglia" = 4093,
  "Endothelial" = 1397,
  "Inhibitory neurons" = 7964,
  "Oligodendrocytes" = 21828
)
got_counts <- setNames(as.integer(table(obj$celltype6)), names(table(obj$celltype6)))
if (!identical(as.integer(got_counts[names(expected_counts)]), as.integer(expected_counts))) {
  stop("SFG celltype6 counts do not match expected corrected mapping")
}

stopifnot(nrow(obj) == 33694)
stopifnot(ncol(obj) == 63608)
stopifnot(ncol(Embeddings(obj, "CCA.ALIGNED")) == 64)
stopifnot(ncol(Embeddings(obj, "umap")) == 2)
stopifnot(sum(table(obj$celltype6)) == ncol(obj))

safe_save_rds(obj, out_file)
write_celltype_qc(obj, "SFG", out_dir)

session_file <- file.path(out_dir, "sessionInfo_stepH_syn21788402_SFG_rebuild.txt")
sink(session_file)
print(sessionInfo())
sink()

test_obj <- readRDS(out_file)
cat("\nReload test passed: ", out_file, "\n", sep = "")
cat("Reloaded genes x cells: ", nrow(test_obj), " x ", ncol(test_obj), "\n", sep = "")
rm(test_obj)
gc()

cat("\nFinished SFG rebuild.\n")
cat("RDS written: ", out_file, "\n", sep = "")
cat("Session info: ", session_file, "\n", sep = "")

quit(save = "no", status = 0, runLast = FALSE)
