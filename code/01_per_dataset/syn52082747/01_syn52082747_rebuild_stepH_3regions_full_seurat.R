###############################################################################
# 01_syn52082747_rebuild_stepH_3regions_full_seurat.R
#
# syn52082747 upstream reconstruction for the Route C manuscript.
#
# Purpose
#   Build the final region-resolved syn52082747 Seurat object used by the
#   manuscript directly from the public author-processed object and source
#   metadata, without requiring the local historical intermediate
#   results/NO3/stepH_slim_uncompressed.rds.
#
# Public/source inputs expected on disk
#   rawdata/NO2_step7_obj_original_backup.rds
#   rawdata/liger_subcluster_metadata_v2.csv
#
# Main output
#   results/3regions/syn52082747_3regions_stepH_slim_uncompressed_full_seurat.rds
#
# The output retains RNA counts and UMAP/PCA/tSNE reductions, adds harmonized
# six-class cell type labels, normalizes source region labels to V1/insula/preCG,
# and computes UBL3 count/log1p(CP10K)/positive metadata for downstream figures.
###############################################################################

rm(list = ls())
gc()
Sys.setenv(LANG = "en")
options(stringsAsFactors = FALSE, scipen = 999)

SEED <- 20251023
set.seed(SEED)

if (dir.exists("D:/R_lib_clean")) {
  .libPaths(c("D:/R_lib_clean", .libPaths()))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(dplyr)
})

project_dir <- Sys.getenv("SYN520_PROJECT_DIR")
if (!nzchar(project_dir)) {
  project_dir <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747"
}

raw_dir <- file.path(project_dir, "rawdata")
out_dir <- file.path(project_dir, "results", "3regions")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw_obj_fp <- file.path(raw_dir, "NO2_step7_obj_original_backup.rds")
source_meta_fp <- file.path(raw_dir, "liger_subcluster_metadata_v2.csv")
full_rds <- file.path(out_dir, "syn52082747_3regions_stepH_slim_uncompressed_full_seurat.rds")
allow_overwrite <- identical(Sys.getenv("ALLOW_OVERWRITE"), "1")

stopifnot(file.exists(raw_obj_fp), file.exists(source_meta_fp))
if (file.exists(full_rds) && !allow_overwrite) {
  stop(
    "Output already exists; refusing to overwrite. ",
    "Set ALLOW_OVERWRITE=1 only after confirming that replacement is intended: ",
    full_rds
  )
}

strip_barcode <- function(x) {
  x <- as.character(x)
  x <- gsub("^.+_", "", x)
  x <- gsub("-1$", "", x)
  x
}

normalize_group4 <- function(x) {
  x0 <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(x0))
  out[grepl("ALZ|ALZH|ALZHEIMER|\\bAD\\b", x0)] <- "AD"
  out[grepl("\\bPSP\\b", x0)] <- "PSP"
  out[grepl("FTD|BVFTD|PICK|PICK'S|\\bPID\\b|FTLD", x0)] <- "FTD"
  out[grepl("CTRL|CONTROL|\\bNC\\b|NORMAL|HEALTH|NO DEMENTIA", x0)] <- "Control"
  out
}

normalize_region3 <- function(x) {
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(z))
  out[grepl("calcarine|visual|\\bv1\\b", z)] <- "V1"
  out[grepl("insula", z)] <- "insula"
  out[grepl("precg|precentral|motor", z)] <- "preCG"
  out
}

standardize_celltype6 <- function(x) {
  z <- trimws(as.character(x))
  zl <- tolower(z)
  out <- rep(NA_character_, length(z))
  out[grepl("astro", zl)] <- "Astrocytes"
  out[grepl("endo", zl)] <- "Endothelial"
  out[grepl("excit", zl)] <- "Excitatory neurons"
  out[grepl("inhib", zl)] <- "Inhibitory neurons"
  out[grepl("micro", zl)] <- "Microglia"
  out[grepl("oligo|opc", zl)] <- "Oligodendrocytes"
  out
}

locate_gene <- function(rn, symbol, ensg = NA_character_) {
  rn0 <- as.character(rn)
  idx <- which(rn0 == symbol)
  if (length(idx) > 0) return(rn0[idx[1]])
  if (!is.na(ensg)) {
    idx <- which(rn0 == ensg)
    if (length(idx) > 0) return(rn0[idx[1]])
    idx <- which(sub("\\.\\d+$", "", rn0) == ensg)
    if (length(idx) > 0) return(rn0[idx[1]])
  }
  idx <- grep(paste0("^", symbol, "$"), rn0, ignore.case = TRUE)
  if (length(idx) > 0) return(rn0[idx[1]])
  NA_character_
}

get_layer <- function(obj, layer) {
  tryCatch(
    SeuratObject::LayerData(obj[["RNA"]], layer = layer),
    error = function(e) NULL
  )
}

message("Reading public author-processed Seurat object: ", raw_obj_fp)
obj <- readRDS(raw_obj_fp)
stopifnot(inherits(obj, "Seurat"))
stopifnot("RNA" %in% Assays(obj))
stopifnot("umap" %in% Reductions(obj))
DefaultAssay(obj) <- "RNA"

message("Reading source metadata: ", source_meta_fp)
meta_csv <- data.table::fread(source_meta_fp, data.table = FALSE)
required_meta_cols <- c("UMI", "sample", "autopsy_id", "region")
missing_meta <- setdiff(required_meta_cols, colnames(meta_csv))
if (length(missing_meta) > 0) {
  stop("Missing required metadata columns: ", paste(missing_meta, collapse = ", "))
}

dx_col <- if ("npdx1" %in% colnames(meta_csv)) {
  "npdx1"
} else if ("clinical_dx" %in% colnames(meta_csv)) {
  "clinical_dx"
} else {
  stop("Could not find npdx1 or clinical_dx in source metadata.")
}

meta_csv$UMI <- as.character(meta_csv$UMI)
idx <- match(strip_barcode(colnames(obj)), strip_barcode(meta_csv$UMI))
if (any(is.na(idx))) {
  stop("Metadata alignment failed for ", sum(is.na(idx)), " cells.")
}
meta_aligned <- meta_csv[idx, , drop = FALSE]
rownames(meta_aligned) <- colnames(obj)

obj$group4 <- normalize_group4(meta_aligned[[dx_col]])
obj$sample <- as.character(meta_aligned$sample)
obj$autopsy_id <- as.character(meta_aligned$autopsy_id)
obj$donor <- as.character(meta_aligned$autopsy_id)
obj$region_raw_syn52082747 <- as.character(meta_aligned$region)
obj$region <- as.character(meta_aligned$region)

celltype6_levels_std <- c(
  "Astrocytes", "Endothelial", "Excitatory neurons",
  "Inhibitory neurons", "Microglia", "Oligodendrocytes"
)

markers_ref <- list(
  Astro = c("AQP4", "GFAP", "ALDH1L1", "SLC1A3", "GPC5", "RYR3"),
  Endo = c("CLDN5", "KDR", "FLT1", "PECAM1", "ABCB1"),
  Excit = c("SLC17A7", "CAMK2A", "TBR1", "CBLN2", "LDB2"),
  Inhib = c("GAD1", "GAD2", "SLC6A1", "PVALB", "SST"),
  Microgl = c("CX3CR1", "P2RY12", "AIF1", "C3", "LRMDA"),
  Oligo = c("MBP", "PLP1", "MOG", "MOBP", "ST18"),
  Peri = c("PDGFRB", "RGS5", "MCAM", "ACTA2")
)

if ("seurat_clusters" %in% colnames(obj@meta.data)) {
  Idents(obj) <- "seurat_clusters"
}
obj$cluster_id <- paste0("c", as.character(Idents(obj)))
clu_levels_c <- paste0("c", levels(Idents(obj)))

data_layer <- get_layer(obj, "data")
if (is.null(data_layer) || nrow(data_layer) == 0 || ncol(data_layer) == 0) {
  message("RNA data layer absent; running LogNormalize for marker scoring.")
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
}

score_one_group <- function(genes) {
  g <- intersect(genes, rownames(obj))
  if (length(g) == 0) {
    out <- rep(NA_real_, length(clu_levels_c))
    names(out) <- clu_levels_c
    return(out)
  }
  av <- AverageExpression(
    obj,
    features = g,
    assays = "RNA",
    layer = "data",
    group.by = "cluster_id",
    verbose = FALSE
  )$RNA
  sc <- colMeans(av, na.rm = TRUE)
  sc2 <- sc[clu_levels_c]
  names(sc2) <- clu_levels_c
  sc2
}

message("Assigning six major cell classes from marker scores.")
avg_by_cluster <- sapply(markers_ref, score_one_group)
tmp <- avg_by_cluster
tmp[is.na(tmp)] <- -Inf
lab_7 <- colnames(tmp)[max.col(tmp, ties.method = "first")]
names(lab_7) <- rownames(tmp)
lab_7[lab_7 == "Peri"] <- "Endo"

short2nice <- c(
  Astro = "Astrocytes",
  Endo = "Endothelial",
  Excit = "Excitatory neurons",
  Inhib = "Inhibitory neurons",
  Microgl = "Microglia",
  Oligo = "Oligodendrocytes"
)
lab_6 <- lab_7
lab_6[lab_6 %in% names(short2nice)] <- short2nice[lab_6[lab_6 %in% names(short2nice)]]
lab_6 <- factor(lab_6, levels = celltype6_levels_std)
map_cluster_to_celltype6 <- setNames(as.character(lab_6), names(lab_6))
obj$celltype6 <- factor(map_cluster_to_celltype6[obj$cluster_id], levels = celltype6_levels_std)

rna_counts <- get_layer(obj, "counts")
if (is.null(rna_counts) || nrow(rna_counts) == 0) {
  stop("RNA counts layer not found.")
}
ubl3_row <- locate_gene(rownames(rna_counts), "UBL3", "ENSG00000122042")
if (is.na(ubl3_row)) stop("UBL3 not found in RNA counts.")
lib_size <- Matrix::colSums(rna_counts)
ubl3_count <- as.numeric(rna_counts[ubl3_row, ])
obj$UBL3_count_routeC <- ubl3_count
obj$UBL3_log1p_CP10k_routeC <- log1p((ubl3_count / pmax(lib_size, 1)) * 1e4)
obj$UBL3_positive_routeC <- ubl3_count > 0

md <- obj@meta.data
routeC_region <- normalize_region3(md$region)
routeC_celltype6 <- standardize_celltype6(md$celltype6)
keep <- !is.na(md$autopsy_id) &
  !is.na(md$group4) & md$group4 %in% c("AD", "Control", "FTD", "PSP") &
  !is.na(routeC_region) & routeC_region %in% c("V1", "insula", "preCG") &
  !is.na(routeC_celltype6)

message("Source cells: ", ncol(obj))
message("Route C retained cells: ", sum(keep))
message("Excluded cells: ", sum(!keep))

obj3 <- subset(obj, cells = colnames(obj)[keep])
md3 <- obj3@meta.data
obj3$routeC_region_raw <- as.character(md3$region)
obj3$routeC_region <- normalize_region3(md3$region)
obj3$region <- obj3$routeC_region
obj3$celltype6 <- factor(standardize_celltype6(md3$celltype6), levels = celltype6_levels_std)
obj3$donor <- as.character(md3$autopsy_id)
obj3$autopsy_id <- as.character(md3$autopsy_id)
obj3$group4 <- factor(as.character(md3$group4), levels = c("AD", "Control", "FTD", "PSP"))
obj3$routeC_unit <- paste0(
  "syn52082747_", as.character(obj3$group4), "_",
  ifelse(obj3$routeC_region == "V1", "V1", obj3$routeC_region)
)

DefaultAssay(obj3) <- "RNA"
rna_counts3 <- get_layer(obj3, "counts")
lib_size3 <- Matrix::colSums(rna_counts3)
ubl3_count3 <- as.numeric(rna_counts3[ubl3_row, ])
obj3$nCount_RNA_routeC <- as.numeric(lib_size3)
obj3$UBL3_count_routeC <- as.numeric(ubl3_count3)
obj3$UBL3_log1p_CP10k_routeC <- log1p((ubl3_count3 / pmax(lib_size3, 1)) * 1e4)
obj3$UBL3_positive_routeC <- obj3$UBL3_count_routeC > 0

message("Reducing object to RNA counts plus standard reductions for GitHub-reproducible downstream use.")
obj3 <- tryCatch(
  DietSeurat(
    obj3,
    assays = "RNA",
    dimreducs = intersect(c("pca", "tsne", "umap"), Reductions(obj3)),
    graphs = NULL,
    layers = c("counts"),
    scale.data = FALSE
  ),
  error = function(e) {
    message("DietSeurat failed; preserving object before slimming. Error: ", conditionMessage(e))
    obj3
  }
)

meta3 <- obj3@meta.data
cell_counts <- as.data.frame(
  table(meta3$group4, meta3$routeC_region, meta3$celltype6, useNA = "ifany"),
  stringsAsFactors = FALSE
)
names(cell_counts) <- c("disease", "region", "celltype6", "n_cells")
cell_counts <- cell_counts[cell_counts$n_cells > 0, ]
cell_counts <- cell_counts[order(cell_counts$disease, cell_counts$region, cell_counts$celltype6), ]

donor_region <- unique(meta3[, c("group4", "routeC_region", "donor")])
donor_counts <- as.data.frame(
  table(donor_region$group4, donor_region$routeC_region, useNA = "ifany"),
  stringsAsFactors = FALSE
)
names(donor_counts) <- c("disease", "region", "n_donors")
donor_counts <- donor_counts[donor_counts$n_donors > 0, ]
donor_counts <- donor_counts[order(donor_counts$disease, donor_counts$region), ]

region_label_summary <- as.data.frame(
  table(
    source_region = as.character(md$region[keep]),
    routeC_region = normalize_region3(md$region[keep]),
    disease = as.character(md$group4[keep]),
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)
names(region_label_summary)[4] <- "n_cells"
region_label_summary <- region_label_summary[region_label_summary$n_cells > 0, ]
region_label_summary <- region_label_summary[
  order(region_label_summary$routeC_region, region_label_summary$source_region, region_label_summary$disease),
]

object_summary <- data.frame(
  metric = c(
    "source_author_object", "source_metadata", "output_rds", "n_genes", "n_cells",
    "assays", "reductions", "default_assay", "has_RNA_counts", "has_umap", "ubl3_row"
  ),
  value = c(
    raw_obj_fp, source_meta_fp, full_rds, nrow(obj3), ncol(obj3),
    paste(Assays(obj3), collapse = ";"),
    paste(Reductions(obj3), collapse = ";"),
    DefaultAssay(obj3),
    "TRUE",
    as.character("umap" %in% Reductions(obj3)),
    ubl3_row
  ),
  stringsAsFactors = FALSE
)

data.table::fwrite(
  cell_counts,
  file.path(out_dir, "syn52082747_3regions_full_seurat_cell_counts_by_disease_region_celltype.csv")
)
data.table::fwrite(
  donor_counts,
  file.path(out_dir, "syn52082747_3regions_full_seurat_donor_counts_by_disease_region.csv")
)
data.table::fwrite(
  region_label_summary,
  file.path(out_dir, "syn52082747_3regions_full_seurat_region_label_summary.csv")
)
data.table::fwrite(
  object_summary,
  file.path(out_dir, "syn52082747_3regions_full_seurat_object_summary.csv")
)

readme <- c(
  "syn52082747 three-region full Seurat object for Route C downstream figures",
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Source author object: ", raw_obj_fp),
  paste0("Source metadata: ", source_meta_fp),
  paste0("Output RDS: ", full_rds),
  "",
  "Purpose:",
  "  This file is the final syn52082747 input for downstream Route C Figure 2 and Supplementary Figure S3.",
  "  It is generated directly from the public author-processed object and source metadata, without using the local historical NO3 StepH object.",
  "",
  "Region convention:",
  "  V1 corresponds to the calcarine cortex label in syn52082747.",
  "  Figure titles should use V1, insula and preCG.",
  "",
  "Required retained fields:",
  "  RNA assay with counts; UMAP reduction; celltype6; donor/autopsy_id; group4; routeC_region; UBL3 count/log1p(CP10K)/positive metadata."
)
writeLines(readme, file.path(out_dir, "README_syn52082747_3regions_full_seurat.txt"), useBytes = TRUE)

session_file <- file.path(out_dir, "sessionInfo_syn52082747_3regions_full_seurat.txt")
sink(session_file)
cat("=== Time ===\n")
print(Sys.time())
cat("\n=== Source/output ===\n")
print(object_summary)
cat("\n=== sessionInfo ===\n")
print(sessionInfo())
sink()

if (file.exists(full_rds) && allow_overwrite) {
  backup_rds <- paste0(full_rds, ".backup_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  message("Backing up existing output RDS to: ", backup_rds)
  ok_backup <- file.rename(full_rds, backup_rds)
  if (!ok_backup) stop("Could not back up existing output RDS: ", full_rds)
}

tmp_rds <- file.path(
  out_dir,
  paste0("syn52082747_3regions_stepH_slim_uncompressed_full_seurat.tmp_", Sys.getpid(), ".rds")
)
message("Saving temporary full Seurat RDS: ", tmp_rds)
saveRDS(obj3, tmp_rds, compress = FALSE, version = 2)

message("Read-back verification")
test <- readRDS(tmp_rds)
stopifnot(inherits(test, "Seurat"))
stopifnot("RNA" %in% Assays(test))
stopifnot("umap" %in% Reductions(test))
stopifnot(all(c(
  "celltype6", "donor", "autopsy_id", "group4", "routeC_region",
  "UBL3_log1p_CP10k_routeC", "UBL3_positive_routeC"
) %in% colnames(test@meta.data)))
rm(test)
gc()

ok <- file.rename(tmp_rds, full_rds)
if (!ok) stop("Could not move temporary RDS to final path.")

message("DONE")
message("full_rds: ", full_rds)
message("n_cells: ", ncol(obj3))
message("size_GB: ", round(file.info(full_rds)$size / 1024^3, 3))
print(object_summary)
print(cell_counts)
print(donor_counts)
