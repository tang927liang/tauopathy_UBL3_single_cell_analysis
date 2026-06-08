# make_syn52082747_stepH_slim_uncompressed.R
#
# Purpose
#   Reproduce the upstream Seurat object used by the manuscript downstream
#   analyses:
#     D:/codex/52082747/redo/stepH_slim_uncompressed.rds
#
# Provenance
#   Cleaned single-file version of the original NO3A_make_stepH_slim_syn52082747.R
#   code block extracted from 52082747.docx. The manuscript parameters are kept
#   unchanged for diagnosis mapping, marker scoring, UBL3 expression, celltype6
#   assignment, DietSeurat slimming, and uncompressed RDS version-2 saving.
#
# Inputs expected on disk
#   D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/rawdata/
#     NO2_step7_obj_original_backup.rds
#     liger_subcluster_metadata_v2.csv
#
# Main output
#   D:/codex/52082747/redo/stepH_slim_uncompressed.rds
#
# Optional historical verification target
#   D:/codex/52082747/stepH_slim_uncompressed.rds
#
# Usage
#   Rscript make_syn52082747_stepH_slim_uncompressed.R
#
# Environment switches
#   RUN_MAKE_STEPH=true|false
#   RUN_VERIFY_STEPH=true|false
#   VERIFY_HASH_RDS=true|false
#   VERIFY_LOAD_HISTORICAL_OBJECT=true|false

rm(list = ls())
gc()
Sys.setenv(LANG = "en")
options(stringsAsFactors = FALSE)

SEED <- 20251023
set.seed(SEED)

.libPaths(c("D:/R_lib_clean", .libPaths()))

project_dir <- "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747"
raw_dir <- file.path(project_dir, "rawdata")
history_dir <- "D:/codex/52082747"
redo_dir <- "D:/codex/52082747/redo"
tmp_dir <- "C:/Temp"

obj_fp <- file.path(raw_dir, "NO2_step7_obj_original_backup.rds")
meta_fp <- file.path(raw_dir, "liger_subcluster_metadata_v2.csv")

dir.create(redo_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

load_required_packages <- function(extra = character()) {
  pkgs <- unique(c(
    "Seurat",
    "SeuratObject",
    "Matrix",
    "data.table",
    "dplyr",
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

log_env <- function(out_dir, step_name) {
  log_file <- file.path(
    out_dir,
    if (grepl("\\.txt$", step_name)) step_name else paste0(step_name, "_sessionInfo.txt")
  )
  sink(log_file)
  cat("=== Time ===\n")
  print(Sys.time())
  cat("\nSEED: ", SEED, "\n", sep = "")
  cat("\nSeurat: ")
  print(packageVersion("Seurat"))
  cat("SeuratObject: ")
  print(packageVersion("SeuratObject"))
  cat("\n=== sessionInfo ===\n")
  print(sessionInfo())
  sink()
  message("Saved sessionInfo: ", log_file)
}

get_layer_data_compat <- function(object, assay = "RNA", layer = "counts") {
  assay_obj <- object[[assay]]
  out <- tryCatch(
    LayerData(assay_obj, layer = layer),
    error = function(e) NULL
  )
  if (!is.null(out)) {
    return(out)
  }
  slot <- if (identical(layer, "counts")) "counts" else "data"
  GetAssayData(object, assay = assay, slot = slot)
}

has_layer_rows <- function(object, assay = "RNA", layer = "data") {
  out <- tryCatch(
    nrow(LayerData(object[[assay]], layer = layer)),
    error = function(e) 0L
  )
  isTRUE(out > 0L)
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

celltype6_levels_std <- c(
  "Astrocytes",
  "Excitatory neurons",
  "Microglia",
  "Endothelial",
  "Inhibitory neurons",
  "Oligodendrocytes"
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

short2nice <- c(
  Astro = "Astrocytes",
  Endo = "Endothelial",
  Excit = "Excitatory neurons",
  Inhib = "Inhibitory neurons",
  Microgl = "Microglia",
  Oligo = "Oligodendrocytes"
)

write_count_percent <- function(x, out_file, label_col) {
  counts <- table(x, useNA = "ifany")
  percent <- round(100 * prop.table(counts), 2)
  df <- data.frame(
    value = names(counts),
    n_cells = as.integer(counts),
    percent = as.numeric(percent),
    stringsAsFactors = FALSE
  )
  names(df)[1] <- label_col
  df[[label_col]][df[[label_col]] == "<NA>"] <- NA_character_
  write.csv(df, out_file, row.names = FALSE)
  invisible(df)
}

make_stepH_slim <- function() {
  load_required_packages(extra = c("ggrepel"))
  stopifnot(file.exists(obj_fp), file.exists(meta_fp))
  log_env(redo_dir, "sessionInfo_NO1_step1.txt")

  message("Read object: ", obj_fp)
  obj <- readRDS(obj_fp)
  stopifnot("RNA" %in% Assays(obj), "umap" %in% Reductions(obj))

  message("Read metadata: ", meta_fp)
  meta_csv <- data.table::fread(meta_fp, data.table = FALSE)
  meta_csv$UMI <- as.character(meta_csv$UMI)

  idx <- match(strip_barcode(colnames(obj)), strip_barcode(meta_csv$UMI))
  meta_aligned <- meta_csv[idx, , drop = FALSE]
  rownames(meta_aligned) <- colnames(obj)

  dx_col <- if ("npdx1" %in% colnames(meta_aligned)) "npdx1" else "clinical_dx"
  obj$group4 <- normalize_group4(meta_aligned[[dx_col]])
  obj$sample <- as.character(meta_aligned$sample)
  obj$autopsy_id <- as.character(meta_aligned$autopsy_id)
  obj$region <- as.character(meta_aligned$region)

  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    Idents(obj) <- "seurat_clusters"
  }
  obj$cluster_id <- paste0("c", as.character(Idents(obj)))
  clu_levels_c <- paste0("c", levels(Idents(obj)))

  DefaultAssay(obj) <- "RNA"
  if (!has_layer_rows(obj, assay = "RNA", layer = "data")) {
    obj <- NormalizeData(
      obj,
      normalization.method = "LogNormalize",
      scale.factor = 1e4,
      verbose = FALSE
    )
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

  avg_by_cluster <- sapply(markers_ref, score_one_group)
  write.csv(avg_by_cluster, file.path(redo_dir, "Check_avgMarkerScores_matrix.csv"))

  tmp <- avg_by_cluster
  tmp[is.na(tmp)] <- -Inf
  lab_7_raw <- colnames(tmp)[max.col(tmp, ties.method = "first")]
  names(lab_7_raw) <- rownames(tmp)

  lab_7_merged <- lab_7_raw
  lab_7_merged[lab_7_merged == "Peri"] <- "Endo"

  lab_6 <- lab_7_merged
  lab_6[lab_6 %in% names(short2nice)] <- short2nice[lab_6[lab_6 %in% names(short2nice)]]
  lab_6 <- factor(lab_6, levels = celltype6_levels_std)

  marker_scores <- data.frame(
    cluster = rownames(avg_by_cluster),
    label_7 = unname(lab_7_raw),
    label_6 = as.character(lab_6),
    avg_by_cluster,
    row.names = NULL,
    check.names = FALSE
  )
  write.csv(
    marker_scores,
    file.path(redo_dir, "Check_markerScores_perCluster.csv"),
    row.names = FALSE
  )

  map_cluster_to_celltype6 <- setNames(as.character(lab_6), names(lab_6))
  celltype6_vec <- map_cluster_to_celltype6[obj$cluster_id]
  names(celltype6_vec) <- colnames(obj)
  obj$celltype6 <- factor(celltype6_vec, levels = celltype6_levels_std)

  write_count_percent(
    obj$group4,
    file.path(redo_dir, "Table_group4_counts_percent.csv"),
    "group4"
  )
  write_count_percent(
    obj$celltype6,
    file.path(redo_dir, "Table_celltype6_counts_percent.csv"),
    "celltype6"
  )

  group4_celltype_counts <- as.data.frame(table(obj$group4, obj$celltype6, useNA = "ifany"))
  names(group4_celltype_counts) <- c("group4", "celltype6", "Freq")
  write.csv(
    group4_celltype_counts,
    file.path(redo_dir, "Table_group4_by_celltype6_counts.csv"),
    row.names = FALSE
  )

  rna_counts <- get_layer_data_compat(obj, assay = "RNA", layer = "counts")
  stopifnot("UBL3" %in% rownames(rna_counts))
  raw_ubl3 <- as.numeric(rna_counts["UBL3", ])
  lib_size <- Matrix::colSums(rna_counts)
  obj$UBL3_log1pCP10K <- log1p((raw_ubl3 / lib_size) * 1e4)
  write.csv(
    data.frame(
      UBL3_min = min(obj$UBL3_log1pCP10K, na.rm = TRUE),
      UBL3_max = max(obj$UBL3_log1pCP10K, na.rm = TRUE)
    ),
    file.path(redo_dir, "Tables_UBL3_range.csv"),
    row.names = FALSE
  )

  obj_stat <- subset(
    obj,
    subset = !is.na(group4) & !is.na(celltype6) & !is.na(sample) & !is.na(autopsy_id)
  )
  df <- data.frame(
    group4 = obj_stat$group4,
    celltype6 = obj_stat$celltype6,
    sample = obj_stat$sample,
    autopsy_id = obj_stat$autopsy_id,
    stringsAsFactors = FALSE
  )
  tabA <- df %>%
    group_by(group4) %>%
    summarise(
      n_cells = n(),
      n_donors = n_distinct(autopsy_id),
      n_samples = n_distinct(sample),
      .groups = "drop"
    )
  write.csv(
    tabA,
    file.path(redo_dir, "Tables_group4_nCells_nDonors_nSamples.csv"),
    row.names = FALSE
  )
  tabB <- df %>%
    group_by(group4, celltype6) %>%
    summarise(
      n_cells = n(),
      n_donors = n_distinct(autopsy_id),
      n_samples = n_distinct(sample),
      .groups = "drop"
    )
  write.csv(
    tabB,
    file.path(redo_dir, "Tables_group4_by_celltype6_nCells_nDonors_nSamples.csv"),
    row.names = FALSE
  )

  rm(obj_stat, df, meta_csv, meta_aligned, rna_counts)
  gc()

  DefaultAssay(obj) <- "RNA"
  obj_slim <- DietSeurat(
    obj,
    assays = "RNA",
    dimreducs = c("pca", "umap", "tsne"),
    graphs = NULL,
    layers = c("counts"),
    scale.data = FALSE
  )
  gc()

  tmp_fp <- file.path(tmp_dir, "stepH_slim_uncompressed_tmp.rds")
  final_fp <- file.path(redo_dir, "stepH_slim_uncompressed.rds")

  message("Save temporary uncompressed RDS: ", tmp_fp)
  saveRDS(obj_slim, tmp_fp, compress = FALSE, version = 2)

  message("Read-back verification: ", tmp_fp)
  test <- readRDS(tmp_fp)
  rm(test)
  gc()

  message("Copy to final output: ", final_fp)
  ok <- file.copy(tmp_fp, final_fp, overwrite = TRUE)
  stopifnot(ok)

  m1 <- tools::md5sum(tmp_fp)
  m2 <- tools::md5sum(final_fp)
  stopifnot(identical(unname(m1), unname(m2)))

  summary_lines <- c(
    paste0("final_fp = ", final_fp),
    paste0("tmp_fp = ", tmp_fp),
    paste0("md5 = ", unname(m2)),
    paste0("size_GB = ", round(file.info(final_fp)$size / 1024^3, 2)),
    paste0("cells = ", ncol(obj_slim)),
    paste0("genes = ", nrow(obj_slim)),
    paste0("assays = ", paste(Assays(obj_slim), collapse = ", ")),
    paste0("reductions = ", paste(Reductions(obj_slim), collapse = ", ")),
    paste0("default_assay = ", DefaultAssay(obj_slim))
  )
  writeLines(summary_lines, file.path(redo_dir, "stepH_slim_uncompressed_summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")

  log_env(redo_dir, "stepH_slim_done")
  invisible(final_fp)
}

compare_csv <- function(filename) {
  redo_file <- file.path(redo_dir, filename)
  hist_file <- file.path(history_dir, filename)
  if (!file.exists(redo_file) || !file.exists(hist_file)) {
    return(data.frame(file = filename, exists = FALSE, identical = NA))
  }
  redo_df <- read.csv(redo_file, check.names = FALSE)
  hist_df <- read.csv(hist_file, check.names = FALSE)
  data.frame(
    file = filename,
    exists = TRUE,
    identical = isTRUE(all.equal(redo_df, hist_df, check.attributes = FALSE)),
    stringsAsFactors = FALSE
  )
}

summarize_obj <- function(fp) {
  message("Reading object for summary: ", fp)
  obj <- readRDS(fp)
  md <- obj@meta.data
  out <- list(
    file = fp,
    size = file.info(fp)$size,
    dim = dim(obj),
    default_assay = DefaultAssay(obj),
    assays = Assays(obj),
    reductions = Reductions(obj),
    group4 = table(md$group4, useNA = "ifany"),
    celltype6 = table(md$celltype6, useNA = "ifany"),
    group4_celltype6 = table(md$group4, md$celltype6, useNA = "ifany")
  )
  rm(obj, md)
  gc()
  out
}

verify_outputs <- function(load_historical_object = FALSE, hash_rds = TRUE) {
  load_required_packages()

  redo_stepH <- file.path(redo_dir, "stepH_slim_uncompressed.rds")
  historical_stepH <- file.path(history_dir, "stepH_slim_uncompressed.rds")
  if (!file.exists(redo_stepH)) stop("Missing regenerated object: ", redo_stepH)
  if (!file.exists(historical_stepH)) stop("Missing historical object: ", historical_stepH)

  csv_files <- c(
    "Check_avgMarkerScores_matrix.csv",
    "Check_markerScores_perCluster.csv",
    "Table_group4_counts_percent.csv",
    "Table_celltype6_counts_percent.csv",
    "Table_group4_by_celltype6_counts.csv",
    "Tables_UBL3_range.csv",
    "Tables_group4_nCells_nDonors_nSamples.csv",
    "Tables_group4_by_celltype6_nCells_nDonors_nSamples.csv"
  )
  csv_cmp <- do.call(rbind, lapply(csv_files, compare_csv))

  redo_summary <- summarize_obj(redo_stepH)
  historical_summary <- NULL
  if (load_historical_object) {
    historical_summary <- summarize_obj(historical_stepH)
  }

  rds_md5 <- data.frame(
    file = c(historical_stepH, redo_stepH),
    md5 = NA_character_,
    stringsAsFactors = FALSE
  )
  if (hash_rds) {
    rds_md5$md5 <- unname(tools::md5sum(rds_md5$file))
  }

  report <- c(
    paste0("redo_stepH = ", redo_stepH),
    paste0("historical_stepH = ", historical_stepH),
    paste0("redo_size = ", file.info(redo_stepH)$size),
    paste0("historical_size = ", file.info(historical_stepH)$size),
    paste0("same_file_size = ", identical(file.info(redo_stepH)$size, file.info(historical_stepH)$size)),
    "",
    "RDS MD5:",
    capture.output(print(rds_md5, row.names = FALSE)),
    paste0("same_rds_md5 = ", length(unique(rds_md5$md5)) == 1L),
    "",
    "CSV comparisons:",
    capture.output(print(csv_cmp, row.names = FALSE)),
    "",
    "Regenerated object summary:",
    paste0("dim = ", paste(redo_summary$dim, collapse = " x ")),
    paste0("default_assay = ", redo_summary$default_assay),
    paste0("assays = ", paste(redo_summary$assays, collapse = ", ")),
    paste0("reductions = ", paste(redo_summary$reductions, collapse = ", ")),
    "",
    "Regenerated group4 counts:",
    capture.output(print(redo_summary$group4)),
    "",
    "Regenerated celltype6 counts:",
    capture.output(print(redo_summary$celltype6))
  )

  if (!is.null(historical_summary)) {
    report <- c(
      report,
      "",
      "Historical object summary:",
      paste0("dim = ", paste(historical_summary$dim, collapse = " x ")),
      paste0("default_assay = ", historical_summary$default_assay),
      paste0("assays = ", paste(historical_summary$assays, collapse = ", ")),
      paste0("reductions = ", paste(historical_summary$reductions, collapse = ", ")),
      paste0("same_dim = ", identical(redo_summary$dim, historical_summary$dim)),
      paste0("same_default_assay = ", identical(redo_summary$default_assay, historical_summary$default_assay)),
      paste0("same_assays = ", identical(redo_summary$assays, historical_summary$assays)),
      paste0("same_reductions = ", identical(redo_summary$reductions, historical_summary$reductions)),
      paste0("same_group4_counts = ", identical(redo_summary$group4, historical_summary$group4)),
      paste0("same_celltype6_counts = ", identical(redo_summary$celltype6, historical_summary$celltype6)),
      paste0("same_group4_celltype6_counts = ", identical(redo_summary$group4_celltype6, historical_summary$group4_celltype6))
    )
  }

  writeLines(report, file.path(redo_dir, "verification_stepH_vs_historical.txt"))
  write.csv(csv_cmp, file.path(redo_dir, "verification_csv_comparisons.csv"), row.names = FALSE)
  write.csv(rds_md5, file.path(redo_dir, "verification_rds_md5.csv"), row.names = FALSE)
  cat(paste(report, collapse = "\n"))
  cat("\n")
  invisible(list(csv = csv_cmp, md5 = rds_md5, summary = redo_summary))
}

run_make <- tolower(Sys.getenv("RUN_MAKE_STEPH", "true")) %in% c("1", "true", "yes", "y")
run_verify <- tolower(Sys.getenv("RUN_VERIFY_STEPH", "true")) %in% c("1", "true", "yes", "y")
load_hist <- tolower(Sys.getenv("VERIFY_LOAD_HISTORICAL_OBJECT", "false")) %in% c("1", "true", "yes", "y")
hash_rds <- tolower(Sys.getenv("VERIFY_HASH_RDS", "true")) %in% c("1", "true", "yes", "y")

if (run_make) {
  make_stepH_slim()
}

if (run_verify) {
  verify_outputs(load_historical_object = load_hist, hash_rds = hash_rds)
}

###############################################################################
# Route C add-on: split syn52082747 into three cortical regions and write
# donor-level UBL3 detection-breadth source files.
#
# This section does not change the historical stepH reconstruction above.
# It reads the regenerated stepH object when available, otherwise the manuscript
# project copy in results/NO3, and writes all region-resolved outputs to:
#   D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/3regions
#
# Environment switches
#   RUN_MAKE_3REGIONS=true|false
#   SAVE_3REGION_SOURCE_RDS=true|false
#   SAVE_3REGION_SEURAT_RDS=true|false
#   SAVE_SPLIT_REGION_RDS=true|false
#   WRITE_CELL_METADATA=true|false
###############################################################################

find_first_col <- function(md, candidates) {
  hit <- intersect(candidates, colnames(md))
  if (length(hit)) hit[1] else NA_character_
}

infer_region_col <- function(md) {
  candidates <- c(
    "region", "Region", "brain_region", "BrainRegion", "tissue", "Tissue",
    "roi", "ROI", "structure", "Structure", "dissected_region"
  )
  hit <- find_first_col(md, candidates)
  if (!is.na(hit)) return(hit)
  target <- c("calcarine", "insula", "precg", "precentral", "v1")
  for (cn in colnames(md)) {
    x <- tolower(trimws(as.character(md[[cn]])))
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) next
    if (any(grepl(paste(target, collapse = "|"), unique(x)))) return(cn)
  }
  NA_character_
}

normalize_region3 <- function(x) {
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(z))
  out[grepl("calcarine|\\bv1\\b|visual", z)] <- "V1"
  out[grepl("insula", z)] <- "insula"
  out[grepl("precg|precentral|pre-central|motor", z)] <- "preCG"
  out
}

standardize_celltype6_routeC <- function(x) {
  z <- trimws(as.character(x))
  z[z %in% c("Astro", "Astrocyte", "Astrocytes")] <- "Astrocytes"
  z[z %in% c("Endo", "Endothelial", "Endothelial cells")] <- "Endothelial"
  z[z %in% c("Excit", "Excitatory", "Excitatory neurons")] <- "Excitatory neurons"
  z[z %in% c("Inhib", "Inhibitory", "Inhibitory neurons")] <- "Inhibitory neurons"
  z[z %in% c("Micro", "Microgl", "Microglia")] <- "Microglia"
  z[z %in% c("Oligo", "Oligodendrocyte", "Oligodendrocytes")] <- "Oligodendrocytes"
  z
}

locate_ubl3_row <- function(rownames_vec) {
  ubl3_ensg <- "ENSG00000122042"
  rn <- rownames_vec
  if ("UBL3" %in% rn) return("UBL3")
  ci <- which(toupper(rn) == "UBL3")
  if (length(ci)) return(rn[ci[1]])
  if (ubl3_ensg %in% rn) return(ubl3_ensg)
  ci <- which(sub("\\.\\d+$", "", rn) == ubl3_ensg)
  if (length(ci)) return(rn[ci[1]])
  NA_character_
}

extract_ubl3_log1p_cp10k <- function(obj, assay = "RNA") {
  DefaultAssay(obj) <- assay
  rna <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(rna), error = function(e) character(0))
  count_layers <- grep("^counts", layers, value = TRUE)
  if (!length(count_layers)) count_layers <- "counts"

  cells_all <- colnames(obj)
  lib_size <- setNames(rep(NA_real_, length(cells_all)), cells_all)
  ubl3_count <- setNames(rep(NA_real_, length(cells_all)), cells_all)
  rows_used <- character(0)

  for (ly in count_layers) {
    m <- tryCatch(SeuratObject::LayerData(rna, layer = ly), error = function(e) NULL)
    if (is.null(m) || !ncol(m)) next
    gene_row <- locate_ubl3_row(rownames(m))
    if (is.na(gene_row)) next
    cn <- intersect(colnames(m), cells_all)
    if (!length(cn)) next
    lib_size[cn] <- Matrix::colSums(m[, cn, drop = FALSE])
    ubl3_count[cn] <- as.numeric(m[gene_row, cn, drop = TRUE])
    rows_used <- c(rows_used, paste0(ly, ":", gene_row))
  }

  if (!length(rows_used)) stop("UBL3 row was not found in any counts layer.", call. = FALSE)
  expr <- log1p((ubl3_count / pmax(lib_size, 1)) * 1e4)
  list(expr = expr, lib_size = lib_size, ubl3_count = ubl3_count, rows_used = unique(rows_used))
}

hl_shift <- function(x, y) {
  if (!length(x) || !length(y)) return(NA_real_)
  median(as.vector(outer(x, y, "-")), na.rm = TRUE)
}

cliffs_delta <- function(x, y) {
  if (!length(x) || !length(y)) return(NA_real_)
  dif <- as.vector(outer(x, y, "-"))
  (sum(dif > 0, na.rm = TRUE) - sum(dif < 0, na.rm = TRUE)) / length(dif)
}

safe_wilcox <- function(x, y) {
  if (length(x) < 1L || length(y) < 1L) {
    return(list(p = NA_real_, conf.low = NA_real_, conf.high = NA_real_))
  }
  wt <- tryCatch(
    stats::wilcox.test(
      x = x, y = y, alternative = "two.sided", exact = FALSE,
      conf.int = TRUE, conf.level = 0.95
    ),
    error = function(e) NULL
  )
  if (is.null(wt)) return(list(p = NA_real_, conf.low = NA_real_, conf.high = NA_real_))
  ci <- suppressWarnings(as.numeric(wt$conf.int))
  if (length(ci) < 2L) ci <- c(NA_real_, NA_real_)
  list(p = as.numeric(wt$p.value), conf.low = ci[1], conf.high = ci[2])
}

make_syn52082747_3regions <- function() {
  load_required_packages(c("data.table", "dplyr"))
  out_dir <- file.path(project_dir, "results", "3regions")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  redo_stepH <- file.path(redo_dir, "stepH_slim_uncompressed.rds")
  project_stepH <- file.path(project_dir, "results", "NO3", "stepH_slim_uncompressed.rds")
  stepH_fp <- if (isTRUE(run_make) && file.exists(redo_stepH)) {
    redo_stepH
  } else if (file.exists(project_stepH)) {
    project_stepH
  } else {
    redo_stepH
  }
  if (!file.exists(stepH_fp)) stop("Cannot find stepH object: ", stepH_fp, call. = FALSE)

  cat("\n========== syn52082747 Route C three-region outputs ==========\n")
  cat("  stepH object:", stepH_fp, "\n")
  cat("  output dir:", out_dir, "\n")

  obj <- readRDS(stepH_fp)
  DefaultAssay(obj) <- "RNA"
  md <- obj@meta.data

  donor_col <- find_first_col(md, c("autopsy_id", "donor", "Donor", "subject", "patient"))
  group_col <- find_first_col(md, c("group4", "group", "Group", "diagnosis", "Diagnosis"))
  ct_col <- find_first_col(md, c("celltype6", "celltype", "cell_type", "celltype6_named"))
  region_col <- infer_region_col(md)
  needed <- c(donor = donor_col, group = group_col, celltype = ct_col, region = region_col)
  if (any(is.na(needed))) {
    stop("Missing required metadata columns: ",
         paste(names(needed)[is.na(needed)], collapse = ", "), call. = FALSE)
  }
  cat("  columns: ", paste(names(needed), needed, sep = "=", collapse = "; "), "\n")

  ubl3 <- extract_ubl3_log1p_cp10k(obj, "RNA")
  cat("  UBL3 rows used:", paste(ubl3$rows_used, collapse = "; "), "\n")

  meta <- data.frame(
    cell = colnames(obj),
    donor = trimws(as.character(md[[donor_col]])),
    group4 = trimws(as.character(md[[group_col]])),
    region_raw = trimws(as.character(md[[region_col]])),
    region = normalize_region3(md[[region_col]]),
    celltype6 = standardize_celltype6_routeC(md[[ct_col]]),
    nCount_RNA_routeC = as.numeric(ubl3$lib_size[colnames(obj)]),
    UBL3_count = as.numeric(ubl3$ubl3_count[colnames(obj)]),
    UBL3_log1p_CP10k = as.numeric(ubl3$expr[colnames(obj)]),
    stringsAsFactors = FALSE
  )
  meta$UBL3_positive <- is.finite(meta$UBL3_log1p_CP10k) & meta$UBL3_log1p_CP10k > 0
  meta <- meta[!is.na(meta$donor) & nzchar(meta$donor) &
                 !is.na(meta$group4) & nzchar(meta$group4) &
                 !is.na(meta$region) &
                 !is.na(meta$celltype6), , drop = FALSE]

  region_levels <- c("V1", "insula", "preCG")
  ct_levels <- c("Astrocytes", "Endothelial", "Excitatory neurons",
                 "Inhibitory neurons", "Microglia", "Oligodendrocytes")
  meta$region <- factor(meta$region, levels = region_levels)
  meta$celltype6 <- factor(meta$celltype6, levels = ct_levels)

  overview <- meta |>
    dplyr::count(group4, region, celltype6, name = "n_cells") |>
    dplyr::arrange(group4, region, celltype6)
  donor_region <- meta |>
    dplyr::count(group4, donor, region, name = "n_cells") |>
    dplyr::arrange(group4, donor, region)
  donor_ct <- meta |>
    dplyr::group_by(group4, donor, region, celltype6) |>
    dplyr::summarise(
      n_cells_total = dplyr::n(),
      n_cells_ubl3pos = sum(UBL3_positive, na.rm = TRUE),
      prop_ubl3pos = n_cells_ubl3pos / n_cells_total,
      median_ubl3_log1p_cp10k_pos = ifelse(
        n_cells_ubl3pos > 0,
        median(UBL3_log1p_CP10k[UBL3_positive], na.rm = TRUE),
        NA_real_
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(group4, donor, region, celltype6)

  diseases <- c("AD", "FTD", "PSP")
  tests <- list()
  for (dz in diseases) {
    for (rg in region_levels) {
      for (ct in ct_levels) {
        d <- donor_ct[donor_ct$region == rg & donor_ct$celltype6 == ct &
                        donor_ct$group4 %in% c(dz, "Control"), , drop = FALSE]
        x <- d$prop_ubl3pos[d$group4 == dz]
        y <- d$prop_ubl3pos[d$group4 == "Control"]
        wt <- safe_wilcox(x, y)
        tests[[length(tests) + 1L]] <- data.frame(
          analytical_unit = paste("syn52082747", dz, rg, sep = "_"),
          disease = dz,
          control = "Control",
          region = rg,
          celltype6 = ct,
          n_disease = length(x),
          n_control = length(y),
          mode = ifelse(length(x) >= 4L && length(y) >= 4L, "powered", "descriptive"),
          disease_median = ifelse(length(x), median(x, na.rm = TRUE), NA_real_),
          control_median = ifelse(length(y), median(y, na.rm = TRUE), NA_real_),
          hl_shift = hl_shift(x, y),
          ci_low = wt$conf.low,
          ci_high = wt$conf.high,
          cliffs_delta = cliffs_delta(x, y),
          p_value = wt$p,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  tests <- do.call(rbind, tests)
  tests$fdr_BH6 <- NA_real_
  for (unit in unique(tests$analytical_unit)) {
    ix <- which(tests$analytical_unit == unit & tests$mode == "powered")
    tests$fdr_BH6[ix] <- p.adjust(tests$p_value[ix], method = "BH")
  }
  tests$is_positive_BH6 <- with(tests, mode == "powered" & is.finite(fdr_BH6) &
                                  fdr_BH6 < 0.05 & hl_shift > 0)
  positives <- tests[tests$is_positive_BH6, , drop = FALSE]

  data.table::fwrite(overview, file.path(out_dir, "syn52082747_3regions_cell_counts_by_group_region_celltype.csv"))
  data.table::fwrite(donor_region, file.path(out_dir, "syn52082747_3regions_donor_region_cell_counts.csv"))
  data.table::fwrite(donor_ct, file.path(out_dir, "syn52082747_3regions_donor_celltype_UBL3_detection.csv"))
  data.table::fwrite(tests, file.path(out_dir, "syn52082747_3regions_donor_level_tests_BH6.csv"))
  data.table::fwrite(positives, file.path(out_dir, "syn52082747_3regions_positive_within_unit_findings.csv"))

  write_cell_metadata <- tolower(Sys.getenv("WRITE_CELL_METADATA", "false")) %in% c("1", "true", "yes", "y")
  if (write_cell_metadata) {
    data.table::fwrite(meta, file.path(out_dir, "syn52082747_3regions_cell_metadata_minimal.csv.gz"))
  }

  save_3region_source_rds <- tolower(Sys.getenv("SAVE_3REGION_SOURCE_RDS", "true")) %in% c("1", "true", "yes", "y")
  if (save_3region_source_rds) {
    source_payload <- list(
      created_at = as.character(Sys.time()),
      source_stepH = stepH_fp,
      region_mapping = c(
        "calcarine / V1 / visual" = "V1",
        "insula" = "insula",
        "preCG / precentral / motor" = "preCG"
      ),
      cell_metadata_minimal = meta,
      cell_counts_by_group_region_celltype = overview,
      donor_region_cell_counts = donor_region,
      donor_celltype_UBL3_detection = donor_ct,
      donor_level_tests_BH6 = tests,
      positive_within_unit_findings = positives,
      sessionInfo = capture.output(sessionInfo())
    )
    saveRDS(
      source_payload,
      file = file.path(out_dir, "syn52082747_3regions_routeC_source.rds"),
      compress = "xz",
      version = 2
    )
    rm(source_payload)
    gc()
  }

  save_3region_seurat_rds <- tolower(Sys.getenv("SAVE_3REGION_SEURAT_RDS", "false")) %in% c("1", "true", "yes", "y")
  if (save_3region_seurat_rds) {
    cells_3region <- intersect(meta$cell, colnames(obj))
    obj_3region <- subset(obj, cells = cells_3region)
    idx3 <- match(colnames(obj_3region), meta$cell)
    obj_3region$region3 <- as.character(meta$region[idx3])
    obj_3region$region_raw_routeC <- as.character(meta$region_raw[idx3])
    obj_3region$UBL3_log1p_CP10k_routeC <- meta$UBL3_log1p_CP10k[idx3]
    obj_3region$UBL3_positive_routeC <- meta$UBL3_positive[idx3]
    obj_3region$nCount_RNA_routeC <- meta$nCount_RNA_routeC[idx3]
    saveRDS(
      obj_3region,
      file = file.path(out_dir, "stepH_slim_uncompressed_syn52082747_3regions_Seurat_optional.rds"),
      compress = "xz",
      version = 2
    )
    rm(obj_3region)
    gc()
  }

  save_split_region_rds <- tolower(Sys.getenv("SAVE_SPLIT_REGION_RDS", "false")) %in% c("1", "true", "yes", "y")
  if (save_split_region_rds) {
    rds_dir <- file.path(out_dir, "split_region_rds_optional")
    dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
    meta_cells <- split(meta$cell, meta$region)
    for (rg in names(meta_cells)) {
      cells_rg <- intersect(meta_cells[[rg]], colnames(obj))
      if (!length(cells_rg)) next
      obj_rg <- subset(obj, cells = cells_rg)
      saveRDS(
        obj_rg,
        file = file.path(rds_dir, paste0("stepH_slim_uncompressed_syn52082747_", rg, ".rds")),
        compress = "xz",
        version = 2
      )
      rm(obj_rg)
      gc()
    }
  }

  readme <- c(
    "syn52082747 three-region Route C outputs",
    "",
    paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Source stepH object: ", stepH_fp),
    paste0("Output directory: ", out_dir),
    "",
    "Region mapping:",
    "  calcarine / V1 / visual -> V1",
    "  insula -> insula",
    "  preCG / precentral / motor -> preCG",
    "",
    "Files:",
    "  syn52082747_3regions_cell_counts_by_group_region_celltype.csv",
    "  syn52082747_3regions_donor_region_cell_counts.csv",
    "  syn52082747_3regions_donor_celltype_UBL3_detection.csv",
    "  syn52082747_3regions_donor_level_tests_BH6.csv",
    "  syn52082747_3regions_positive_within_unit_findings.csv",
    "  syn52082747_3regions_routeC_source.rds if SAVE_3REGION_SOURCE_RDS=true",
    "  stepH_slim_uncompressed_syn52082747_3regions_Seurat_optional.rds only if SAVE_3REGION_SEURAT_RDS=true",
    "  split_region_rds_optional/*.rds only if SAVE_SPLIT_REGION_RDS=true",
    "",
    "Statistical note:",
    "  Tests compare each disease group with Control within each region and cell type.",
    "  BH FDR is applied within each disease-region analytical unit over six cell types.",
    "  Rows with fewer than 4 donors per group are marked descriptive."
  )
  writeLines(readme, file.path(out_dir, "README_syn52082747_3regions.txt"), useBytes = TRUE)
  log_env(out_dir, "sessionInfo_syn52082747_3regions.txt")

  cat("\nThree-region Route C output summary:\n")
  cat("  cells used:", nrow(meta), "\n")
  cat("  donor-celltype rows:", nrow(donor_ct), "\n")
  cat("  test rows:", nrow(tests), "\n")
  cat("  positive BH6 rows:", nrow(positives), "\n")
  if (nrow(positives)) print(positives[, c("analytical_unit", "celltype6", "hl_shift", "fdr_BH6")])
  invisible(list(meta = meta, donor_celltype = donor_ct, tests = tests, positives = positives))
}

run_make_3regions <- tolower(Sys.getenv("RUN_MAKE_3REGIONS", "true")) %in% c("1", "true", "yes", "y")
if (run_make_3regions) {
  make_syn52082747_3regions()
}
