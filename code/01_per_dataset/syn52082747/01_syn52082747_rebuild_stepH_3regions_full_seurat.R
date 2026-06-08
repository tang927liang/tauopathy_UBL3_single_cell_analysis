rm(list = ls())
gc()
Sys.setenv(LANG = "en")
options(stringsAsFactors = FALSE, scipen = 999)

SEED <- 20251023
set.seed(SEED)

.libPaths(c("D:/R_lib_clean", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

project_dir <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747"
script_path <- file.path(project_dir, "R",
                         "01_syn52082747_rebuild_stepH_3regions_full_seurat.R")
source_stepH <- file.path(project_dir, "results", "NO3", "stepH_slim_uncompressed.rds")
out_dir <- file.path(project_dir, "results", "3regions")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

full_rds <- file.path(out_dir, "syn52082747_3regions_stepH_slim_uncompressed_full_seurat.rds")
if (file.exists(full_rds)) {
  stop("Output already exists; refusing to overwrite: ", full_rds)
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
  }
  idx <- grep(paste0("^", symbol, "$"), rn0, ignore.case = TRUE)
  if (length(idx) > 0) return(rn0[idx[1]])
  NA_character_
}

message("Reading source Seurat object: ", source_stepH)
stopifnot(file.exists(source_stepH))
obj <- readRDS(source_stepH)

stopifnot(inherits(obj, "Seurat"))
stopifnot("RNA" %in% Assays(obj))
stopifnot("umap" %in% Reductions(obj))

md <- obj@meta.data
required <- c("autopsy_id", "group4", "region", "celltype6")
missing <- setdiff(required, colnames(md))
if (length(missing) > 0) {
  stop("Missing required metadata columns in source object: ", paste(missing, collapse = ", "))
}

routeC_region <- normalize_region3(md$region)
routeC_celltype6 <- standardize_celltype6(md$celltype6)
keep <- !is.na(md$autopsy_id) &
  !is.na(md$group4) & md$group4 %in% c("AD", "Control", "FTD", "PSP") &
  !is.na(routeC_region) & routeC_region %in% c("V1", "insula", "preCG") &
  !is.na(routeC_celltype6)

message("Source cells: ", ncol(obj))
message("Route C retained cells: ", sum(keep))
message("Excluded cells: ", sum(!keep))

cells_keep <- colnames(obj)[keep]
obj3 <- subset(obj, cells = cells_keep)
idx <- match(colnames(obj3), rownames(md))
if (any(is.na(idx))) stop("Internal metadata alignment failed")

obj3$routeC_region_raw <- as.character(md$region[idx])
obj3$routeC_region <- normalize_region3(md$region[idx])
obj3$region <- obj3$routeC_region
obj3$celltype6 <- factor(standardize_celltype6(md$celltype6[idx]), levels = c(
  "Astrocytes",
  "Endothelial",
  "Excitatory neurons",
  "Inhibitory neurons",
  "Microglia",
  "Oligodendrocytes"
))
obj3$donor <- as.character(md$autopsy_id[idx])
obj3$autopsy_id <- as.character(md$autopsy_id[idx])
obj3$group4 <- factor(as.character(md$group4[idx]), levels = c("AD", "Control", "FTD", "PSP"))
obj3$routeC_unit <- paste0("syn52082747_", as.character(obj3$group4), "_",
                           ifelse(obj3$routeC_region == "V1", "calcarine", obj3$routeC_region))

DefaultAssay(obj3) <- "RNA"
rna_counts <- GetAssayData(obj3, assay = "RNA", layer = "counts")
ubl3_row <- locate_gene(rownames(rna_counts), "UBL3", "ENSG00000122042")
if (is.na(ubl3_row)) stop("UBL3 not found in RNA counts")
lib_size <- Matrix::colSums(rna_counts)
ubl3_count <- as.numeric(rna_counts[ubl3_row, ])
obj3$nCount_RNA_routeC <- as.numeric(lib_size)
obj3$UBL3_count_routeC <- as.numeric(ubl3_count)
obj3$UBL3_log1p_CP10k_routeC <- log1p((ubl3_count / pmax(lib_size, 1)) * 1e4)
obj3$UBL3_positive_routeC <- obj3$UBL3_count_routeC > 0

meta3 <- obj3@meta.data
cell_counts <- as.data.frame(
  table(meta3$group4, meta3$routeC_region, meta3$celltype6, useNA = "ifany"),
  stringsAsFactors = FALSE
)
names(cell_counts) <- c("disease", "region", "celltype6", "n_cells")
cell_counts <- cell_counts[cell_counts$n_cells > 0, ]
cell_counts <- cell_counts[order(cell_counts$disease, cell_counts$region, cell_counts$celltype6), ]

donor_counts <- unique(meta3[, c("group4", "routeC_region", "donor")])
donor_counts <- as.data.frame(
  table(donor_counts$group4, donor_counts$routeC_region, useNA = "ifany"),
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
region_label_summary <- region_label_summary[order(region_label_summary$routeC_region,
                                                   region_label_summary$source_region,
                                                   region_label_summary$disease), ]

object_summary <- data.frame(
  metric = c("source_stepH", "output_rds", "n_genes", "n_cells", "assays", "reductions",
             "default_assay", "has_RNA_counts", "has_umap", "ubl3_row"),
  value = c(source_stepH, full_rds, nrow(obj3), ncol(obj3),
            paste(Assays(obj3), collapse = ";"),
            paste(Reductions(obj3), collapse = ";"),
            DefaultAssay(obj3),
            "TRUE",
            as.character("umap" %in% Reductions(obj3)),
            ubl3_row),
  stringsAsFactors = FALSE
)

fwrite(cell_counts,
       file.path(out_dir, "syn52082747_3regions_full_seurat_cell_counts_by_disease_region_celltype.csv"))
fwrite(donor_counts,
       file.path(out_dir, "syn52082747_3regions_full_seurat_donor_counts_by_disease_region.csv"))
fwrite(region_label_summary,
       file.path(out_dir, "syn52082747_3regions_full_seurat_region_label_summary.csv"))
fwrite(object_summary,
       file.path(out_dir, "syn52082747_3regions_full_seurat_object_summary.csv"))

readme <- c(
  "syn52082747 three-region full Seurat object for Route C downstream figures",
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Source stepH object: ", source_stepH),
  paste0("Output RDS: ", full_rds),
  "",
  "Purpose:",
  "  This file is the final syn52082747 input for downstream Route C Figure 2 and Supplementary Figure S3.",
  "  The historical NO3 stepH object is read only once here to preserve RNA counts/data and UMAP, then region labels are normalized to V1, insula, and preCG.",
  "",
  "Region convention:",
  "  V1 corresponds to the calcarine cortex label in syn52082747.",
  "  Figure titles should use V1, not calcarine cortex/V1.",
  "",
  "Required retained fields:",
  "  RNA assay with counts; UMAP reduction; celltype6; donor/autopsy_id; group4; routeC_region; UBL3 count/log1p CP10K/positive metadata.",
  "",
  "Do not use the old NO3 stepH object as a final downstream Figure2/S3 input."
)
writeLines(readme,
           file.path(out_dir, "README_syn52082747_3regions_full_seurat.txt"),
           useBytes = TRUE)

session_file <- file.path(out_dir, "sessionInfo_syn52082747_3regions_full_seurat.txt")
sink(session_file)
cat("=== Time ===\n")
print(Sys.time())
cat("\n=== Script ===\n")
cat(normalizePath(script_path, winslash = "/", mustWork = FALSE), "\n")
cat("\n=== Source/output ===\n")
print(object_summary)
cat("\n=== sessionInfo ===\n")
print(sessionInfo())
sink()

tmp_rds <- file.path(out_dir, paste0("syn52082747_3regions_stepH_slim_uncompressed_full_seurat.tmp_",
                                     Sys.getpid(), ".rds"))
message("Saving temporary full Seurat RDS: ", tmp_rds)
saveRDS(obj3, tmp_rds, compress = FALSE, version = 2)

message("Read-back verification")
test <- readRDS(tmp_rds)
stopifnot(inherits(test, "Seurat"))
stopifnot("RNA" %in% Assays(test))
stopifnot("umap" %in% Reductions(test))
stopifnot(all(c("celltype6", "donor", "autopsy_id", "group4", "routeC_region",
                "UBL3_log1p_CP10k_routeC", "UBL3_positive_routeC") %in% colnames(test@meta.data)))
rm(test)
gc()

ok <- file.rename(tmp_rds, full_rds)
if (!ok) stop("Could not move temporary RDS to final path")

message("DONE")
message("full_rds: ", full_rds)
message("n_cells: ", ncol(obj3))
message("size_GB: ", round(file.info(full_rds)$size / 1024^3, 3))
print(object_summary)
print(cell_counts)
print(donor_counts)

