.libPaths(c("D:/R_lib_clean", "C:/Program Files/R/R-4.5.3/library"))

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(openxlsx)
  library(DESeq2)
})

run_started <- Sys.time()

base_results_dir <- "D:/RNA/2026063Molecular Neurodegeneration/Supplementary_Data_1/results/RouteC_20260605_pseudobulk_DESeq2"

make_safe_out_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    return(path)
  }
  existing <- list.files(path, all.files = TRUE, no.. = TRUE)
  if (length(existing) == 0) return(path)
  for (i in 2:99) {
    candidate <- paste0(path, "_v", i)
    if (!dir.exists(candidate)) {
      dir.create(candidate, recursive = TRUE, showWarnings = FALSE)
      return(candidate)
    }
    if (length(list.files(candidate, all.files = TRUE, no.. = TRUE)) == 0) return(candidate)
  }
  stop("Could not create a non-overwriting output directory.")
}

out_dir <- make_safe_out_dir(base_results_dir)

source_paths <- data.table(
  source_key = c("GSE157827", "GSE174367", "syn21788402_EC", "syn21788402_SFG", "syn52082747"),
  accession = c("GSE157827", "GSE174367", "syn21788402", "syn21788402", "syn52082747"),
  rds_path = c(
    "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results/stepH_obj_celltype6_named.rds",
    "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/stepH_obj_celltype6_named.rds",
    "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds",
    "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepH_syn21788402_SFG_obj_celltype6.rds",
    "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/3regions/syn52082747_3regions_stepH_slim_uncompressed_full_seurat.rds"
  )
)

unit_order <- data.table(
  analytical_unit_order = 1:13,
  source_key = c(
    "GSE157827", "GSE174367",
    rep("syn52082747", 3),
    "syn21788402_EC", "syn21788402_SFG",
    rep("syn52082747", 6)
  ),
  accession = c(
    "GSE157827", "GSE174367",
    rep("syn52082747", 3),
    "syn21788402", "syn21788402",
    rep("syn52082747", 6)
  ),
  disease = c("AD", "AD", "AD", "AD", "AD", "AD", "AD", "FTD", "FTD", "FTD", "PSP", "PSP", "PSP"),
  region = c("PFC", "PFC", "V1", "insula", "preCG", "EC", "SFG", "V1", "insula", "preCG", "V1", "insula", "preCG")
)
unit_order[, analytical_unit := sprintf("%s (%s), %s", accession, disease, region)]

celltype_order <- c(
  "Astrocytes", "Endothelial", "Excitatory neurons",
  "Inhibitory neurons", "Microglia", "Oligodendrocytes"
)

pick_first_col <- function(meta, candidates, required = TRUE) {
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) > 0) return(hit[1])
  if (required) stop("Missing required metadata column among: ", paste(candidates, collapse = ", "))
  NA_character_
}

clean_text <- function(x) {
  y <- as.character(x)
  trimws(y)
}

standardize_celltype <- function(x) {
  y <- clean_text(x)
  yl <- tolower(y)
  out <- rep(NA_character_, length(y))
  out[grepl("astro", yl)] <- "Astrocytes"
  out[grepl("endothel", yl)] <- "Endothelial"
  out[grepl("excit", yl)] <- "Excitatory neurons"
  out[grepl("inhibit", yl)] <- "Inhibitory neurons"
  out[grepl("micro", yl)] <- "Microglia"
  out[grepl("oligo", yl)] <- "Oligodendrocytes"
  exact <- y %in% celltype_order
  out[is.na(out) & exact] <- y[is.na(out) & exact]
  out
}

standardize_disease <- function(x) {
  y <- clean_text(x)
  yl <- toupper(y)
  out <- rep(NA_character_, length(y))
  num <- suppressWarnings(as.numeric(y))
  out[!is.na(num) & num >= 5] <- "AD"
  out[!is.na(num) & num == 0] <- "Control"
  out[grepl("CONTROL|CTRL|CON\\b|NORMAL", yl)] <- "Control"
  out[grepl("\\bAD\\b|ALZHEIMER", yl)] <- "AD"
  out[grepl("\\bFTD\\b|FRONTOTEMPORAL", yl)] <- "FTD"
  out[grepl("\\bPSP\\b|PROGRESSIVE", yl)] <- "PSP"
  out
}

standardize_region <- function(x, fixed_region = NULL) {
  if (!is.null(fixed_region)) return(rep(fixed_region, length(x)))
  y <- clean_text(x)
  yl <- tolower(y)
  out <- rep(NA_character_, length(y))
  out[grepl("v1|calcarine|primary visual|visual cortex", yl)] <- "V1"
  out[grepl("insula", yl)] <- "insula"
  out[grepl("precg|precentral|pre-central", yl)] <- "preCG"
  out[is.na(out) & y %in% c("V1", "insula", "preCG")] <- y[is.na(out) & y %in% c("V1", "insula", "preCG")]
  out
}

prepare_metadata <- function(obj, source_key) {
  meta <- as.data.table(obj@meta.data, keep.rownames = "cell")
  if (source_key == "GSE157827") {
    donor_col <- pick_first_col(meta, c("sample", "SampleID", "donor", "donor_id"))
    group_col <- pick_first_col(meta, c("group", "Diagnosis", "diagnosis", "disease"))
    region <- rep("PFC", nrow(meta))
  } else if (source_key == "GSE174367") {
    donor_col <- pick_first_col(meta, c("sample", "SampleID", "donor", "donor_id"))
    group_col <- pick_first_col(meta, c("group", "Diagnosis", "diagnosis", "disease"))
    region <- rep("PFC", nrow(meta))
  } else if (source_key == "syn21788402_EC") {
    donor_col <- pick_first_col(meta, c("PatientID", "patient_id", "SampleID", "sample", "donor", "donor_id"))
    group_col <- pick_first_col(meta, c("BraakStage", "braak_stage", "group", "Diagnosis", "diagnosis", "disease"))
    region <- rep("EC", nrow(meta))
  } else if (source_key == "syn21788402_SFG") {
    donor_col <- pick_first_col(meta, c("PatientID", "patient_id", "SampleID", "sample", "donor", "donor_id"))
    group_col <- pick_first_col(meta, c("BraakStage", "braak_stage", "group", "Diagnosis", "diagnosis", "disease"))
    region <- rep("SFG", nrow(meta))
  } else if (source_key == "syn52082747") {
    donor_col <- pick_first_col(meta, c("autopsy_id", "donor", "sample", "SampleID", "donor_id"))
    group_col <- pick_first_col(meta, c("group4", "group", "Diagnosis", "diagnosis", "disease"))
    region_col <- pick_first_col(meta, c("routeC_region", "region", "Region"))
    region <- standardize_region(meta[[region_col]])
  } else {
    stop("Unhandled source key: ", source_key)
  }
  cell_col <- pick_first_col(meta, c("celltype6", "cell_type", "celltype", "cell_type6"))
  data.table(
    cell = meta$cell,
    donor = clean_text(meta[[donor_col]]),
    disease_group = standardize_disease(meta[[group_col]]),
    region = region,
    celltype = standardize_celltype(meta[[cell_col]])
  )
}

get_count_layers <- function(obj, assay = "RNA") {
  assay_obj <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(assay_obj), error = function(e) character())
  count_layers <- grep("^counts", layers, value = TRUE)
  if (length(count_layers) == 0) count_layers <- "counts"
  count_layers
}

get_layer_matrix <- function(obj, layer, assay = "RNA") {
  if (layer == "counts") {
    return(tryCatch(
      SeuratObject::LayerData(obj, assay = assay, layer = "counts"),
      error = function(e) SeuratObject::GetAssayData(obj, assay = assay, layer = "counts")
    ))
  }
  SeuratObject::LayerData(obj, assay = assay, layer = layer)
}

get_feature_template <- function(obj, assay = "RNA") {
  count_layers <- get_count_layers(obj, assay = assay)
  mat <- get_layer_matrix(obj, count_layers[1], assay = assay)
  rownames(mat)
}

make_gene_annotation <- function(features) {
  stripped <- sub("\\.\\d+$", "", features)
  is_ensg <- grepl("^ENSG[0-9]+", stripped)
  data.table(
    raw_feature = features,
    gene_id = fifelse(is_ensg, stripped, features),
    gene_symbol = fifelse(is_ensg, NA_character_, features)
  )
}

aggregate_pseudobulk_counts <- function(obj, selected_meta, features, assay = "RNA") {
  donors <- sort(unique(selected_meta$donor))
  counts_out <- matrix(
    0,
    nrow = length(features),
    ncol = length(donors),
    dimnames = list(features, donors)
  )
  donor_by_cell <- setNames(selected_meta$donor, selected_meta$cell)
  selected_cells <- selected_meta$cell
  count_layers <- get_count_layers(obj, assay = assay)
  layer_cells_used <- 0L
  for (layer in count_layers) {
    mat <- get_layer_matrix(obj, layer, assay = assay)
    if (!identical(rownames(mat), features)) {
      stop("Feature order mismatch in layer ", layer)
    }
    cn <- intersect(colnames(mat), selected_cells)
    if (length(cn) == 0) next
    layer_cells_used <- layer_cells_used + length(cn)
    layer_donors <- donor_by_cell[cn]
    for (donor_id in unique(layer_donors)) {
      donor_cells <- cn[layer_donors == donor_id]
      if (length(donor_cells) == 1) {
        counts_out[, donor_id] <- counts_out[, donor_id] + as.numeric(mat[, donor_cells, drop = TRUE])
      } else {
        counts_out[, donor_id] <- counts_out[, donor_id] + as.numeric(Matrix::rowSums(mat[, donor_cells, drop = FALSE]))
      }
    }
  }
  list(counts = counts_out, n_cells_used = layer_cells_used)
}

run_deseq_unit_celltype <- function(count_mat, coldata) {
  keep <- rowSums(count_mat) > 0
  count_mat <- count_mat[keep, , drop = FALSE]
  if (nrow(count_mat) == 0) stop("No expressed genes after pseudobulk aggregation.")
  count_mat <- round(count_mat)
  coldata <- as.data.frame(coldata)
  rownames(coldata) <- coldata$donor
  coldata <- coldata[colnames(count_mat), , drop = FALSE]
  coldata$condition <- factor(coldata$condition, levels = c("Control", "Disease"))
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = coldata,
    design = ~condition
  )
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "Disease", "Control"), independentFiltering = FALSE)
  list(result = as.data.table(res, keep.rownames = "raw_feature"), fit = "DESeq2_default")
}

sanitize_dt <- function(dt) {
  for (nm in names(dt)) {
    if (is.character(dt[[nm]])) {
      set(dt, j = nm, value = clean_text(dt[[nm]]))
    }
  }
  dt
}

result_rows <- list()
input_checks <- list()
pseudobulk_manifest_rows <- list()
source_manifest_rows <- list()

for (src_i in seq_len(nrow(source_paths))) {
  src <- source_paths[src_i]
  message("Reading ", src$source_key, ": ", src$rds_path)
  obj <- readRDS(src$rds_path)
  meta <- prepare_metadata(obj, src$source_key)
  features <- get_feature_template(obj, assay = "RNA")
  gene_annot <- make_gene_annotation(features)
  src_units <- unit_order[source_key == src$source_key]

  source_manifest_rows[[length(source_manifest_rows) + 1L]] <- data.table(
    source_key = src$source_key,
    accession = src$accession,
    rds_path = src$rds_path,
    n_features_RNA = length(features),
    n_cells = ncol(obj),
    count_layers = paste(get_count_layers(obj, assay = "RNA"), collapse = ";"),
    metadata_columns = paste(colnames(obj@meta.data), collapse = ";")
  )

  for (u_i in seq_len(nrow(src_units))) {
    unit <- src_units[u_i]
    for (ct in celltype_order) {
      message("DESeq2: ", unit$analytical_unit, " | ", ct)
      selected <- meta[
        region == unit$region &
          disease_group %in% c(unit$disease, "Control") &
          celltype == ct &
          !is.na(donor) & donor != "" &
          !is.na(disease_group)
      ]
      donor_table <- selected[, .(
        n_cells = .N
      ), by = .(donor, disease_group)]
      n_disease <- donor_table[disease_group == unit$disease, uniqueN(donor)]
      n_control <- donor_table[disease_group == "Control", uniqueN(donor)]
      analytical_mode <- if (unit$source_key %in% c("syn21788402_EC", "syn21788402_SFG")) {
        "descriptive_only_n3v3"
      } else if (n_disease >= 4 && n_control >= 4) {
        "inferential_n_at_least_4"
      } else {
        "descriptive_low_n"
      }

      input_checks[[length(input_checks) + 1L]] <- data.table(
        analytical_unit_order = unit$analytical_unit_order,
        accession = unit$accession,
        analytical_unit = unit$analytical_unit,
        disease = unit$disease,
        region = unit$region,
        celltype = ct,
        n_disease_donors = n_disease,
        n_control_donors = n_control,
        n_cells_selected = nrow(selected),
        analytical_mode = analytical_mode
      )

      if (n_disease < 2 || n_control < 2 || nrow(selected) == 0) {
        warning("Skipping low-donor unit/cell type: ", unit$analytical_unit, " | ", ct)
        next
      }

      agg <- aggregate_pseudobulk_counts(obj, selected, features, assay = "RNA")
      count_mat <- agg$counts
      sample_info <- donor_table[, .(
        donor,
        disease_group,
        n_cells
      )]
      sample_info[, condition := fifelse(disease_group == unit$disease, "Disease", "Control")]
      sample_info <- sample_info[donor %in% colnames(count_mat)]
      setorder(sample_info, condition, donor)
      count_mat <- count_mat[, sample_info$donor, drop = FALSE]

      pseudobulk_manifest_rows[[length(pseudobulk_manifest_rows) + 1L]] <- data.table(
        analytical_unit_order = unit$analytical_unit_order,
        accession = unit$accession,
        analytical_unit = unit$analytical_unit,
        disease = unit$disease,
        region = unit$region,
        celltype = ct,
        donor = sample_info$donor,
        group = sample_info$disease_group,
        condition = sample_info$condition,
        n_cells = sample_info$n_cells,
        total_pseudobulk_counts = colSums(count_mat),
        analytical_mode = analytical_mode
      )

      fit <- tryCatch(
        run_deseq_unit_celltype(count_mat, sample_info),
        error = function(e) {
          warning("DESeq2 failed for ", unit$analytical_unit, " | ", ct, ": ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(fit)) next
      res <- fit$result
      res <- merge(res, gene_annot, by = "raw_feature", all.x = TRUE, sort = FALSE)
      res[, `:=`(
        accession = unit$accession,
        analytical_unit_order = unit$analytical_unit_order,
        analytical_unit = unit$analytical_unit,
        disease = unit$disease,
        region = unit$region,
        celltype = ct,
        n_disease_donors = n_disease,
        n_control_donors = n_control,
        analytical_mode = analytical_mode,
        provenance_note = paste0(
          "Route C donor-level pseudobulk DESeq2; ",
          if (unit$accession == "syn52082747") "syn52082747 full three-region Seurat RDS; " else "",
          "disease_vs_control"
        ),
        deseq2_fit = fit$fit
      )]
      setcolorder(res, c(
        "accession", "analytical_unit_order", "analytical_unit", "disease", "region", "celltype",
        "gene_id", "gene_symbol", "raw_feature",
        "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
        "n_disease_donors", "n_control_donors", "analytical_mode", "provenance_note", "deseq2_fit"
      ))
      result_rows[[length(result_rows) + 1L]] <- sanitize_dt(res)
      rm(agg, count_mat, fit, res)
      invisible(gc())
    }
  }
  rm(obj, meta)
  invisible(gc())
}

all_results <- rbindlist(result_rows, fill = TRUE)
input_check <- rbindlist(input_checks, fill = TRUE)
pseudobulk_manifest <- rbindlist(pseudobulk_manifest_rows, fill = TRUE)
source_manifest <- rbindlist(source_manifest_rows, fill = TRUE)

all_results[, celltype_order_index := match(celltype, celltype_order)]
setorder(all_results, analytical_unit_order, celltype_order_index, gene_id)
all_results[, celltype_order_index := NULL]

input_check[, celltype_order_index := match(celltype, celltype_order)]
setorder(input_check, analytical_unit_order, celltype_order_index)
input_check[, celltype_order_index := NULL]

pseudobulk_manifest[, celltype_order_index := match(celltype, celltype_order)]
setorder(pseudobulk_manifest, analytical_unit_order, celltype_order_index, condition, donor)
pseudobulk_manifest[, celltype_order_index := NULL]

row_count_by_unit <- all_results[, .(n_rows = .N), by = .(analytical_unit_order, analytical_unit)]
setorder(row_count_by_unit, analytical_unit_order)
row_count_by_celltype <- all_results[, .(n_rows = .N), by = .(celltype)]
row_count_by_celltype[, celltype_order_index := match(celltype, celltype_order)]
setorder(row_count_by_celltype, celltype_order_index)
row_count_by_celltype[, celltype_order_index := NULL]

readme <- data.table(
  field = c(
    "table_title",
    "analysis_type",
    "routeC_region_resolution",
    "syn52082747_source",
    "syn52082747_V1_note",
    "GSE157827_region_note",
    "syn21788402_note",
    "statistical_unit",
    "deseq2_replicate",
    "not_main_endpoint_table",
    "zenodo_note",
    "formal_storyline_note",
    "old_data1_audit",
    "total_rows",
    "run_started",
    "run_completed"
  ),
  value = c(
    "Supplementary Data 1. Donor-level pseudobulk DESeq2 differential-expression output across Route C analytical units",
    "This is complete gene-level donor-level pseudobulk DESeq2 output.",
    "Route C region-resolved structure is used; syn52082747 is split into V1, insula, and preCG.",
    "syn52082747 was read from the full three-region Seurat RDS, not the old pooled NO3 object.",
    "V1 corresponds to the calcarine cortex metadata label in syn52082747.",
    "GSE157827 is reported as PFC.",
    "syn21788402 EC/SFG are descriptive_only_n3v3.",
    "All differential-expression models use donor-level pseudobulk samples.",
    "Cells/nuclei are aggregated by donor before DESeq2 and are not treated as independent DESeq2 replicates.",
    "This table is not the UBL3 main endpoint table; UBL3 donor-level detection-breadth and conditional-expression statistics are in Supplementary Table S1.",
    "If the workbook exceeds the Molecular Neurodegeneration additional-file limit, deposit this data resource in Zenodo.",
    "The focal formal manuscript claim remains a PSP V1 cortical-neuron UBL3 detection-breadth candidate signal, without positive-cell conditional-expression upregulation.",
    "The existing root-level Supplementary_Data_1.xlsx was audited as old seven-unit/pooled syn52082747 output and was not reused as final Route C input.",
    as.character(nrow(all_results)),
    as.character(run_started),
    as.character(Sys.time())
  )
)

session_lines <- capture.output(sessionInfo())
session_dt <- data.table(line = seq_along(session_lines), text = session_lines)

csv_all <- file.path(out_dir, "Supplementary_Data_1_RouteC_DESeq2_all.csv")
csv_unit <- file.path(out_dir, "row_count_by_analytical_unit.csv")
csv_celltype <- file.path(out_dir, "row_count_by_celltype.csv")
csv_input <- file.path(out_dir, "input_check.csv")
csv_pseudobulk <- file.path(out_dir, "pseudobulk_sample_manifest.csv")
csv_source <- file.path(out_dir, "source_data_manifest.csv")

fwrite(all_results, csv_all)
fwrite(row_count_by_unit, csv_unit)
fwrite(row_count_by_celltype, csv_celltype)
fwrite(input_check, csv_input)
fwrite(pseudobulk_manifest, csv_pseudobulk)
fwrite(source_manifest, csv_source)
writeLines(readme[, paste(field, value, sep = ": ")], file.path(out_dir, "README_Supplementary_Data_1.txt"))
writeLines(session_lines, file.path(out_dir, "sessionInfo.txt"))

wb <- createWorkbook()
header_style <- createStyle(textDecoration = "bold")
add_sheet <- function(sheet, data) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet, data)
  if (ncol(data) > 0 && nrow(data) > 0) {
    addStyle(wb, sheet, header_style, rows = 1, cols = seq_len(ncol(data)), gridExpand = TRUE)
    freezePane(wb, sheet, firstRow = TRUE)
  }
}

add_sheet("README", readme)
add_sheet("row_count_by_unit", row_count_by_unit)
add_sheet("row_count_by_celltype", row_count_by_celltype)
add_sheet("input_check", input_check)
add_sheet("pseudobulk_sample_manifest", pseudobulk_manifest)
add_sheet("source_data_manifest", source_manifest)
add_sheet("session_info", session_dt)

max_rows_per_sheet <- 900000L
n_parts <- ceiling(nrow(all_results) / max_rows_per_sheet)
for (part in seq_len(n_parts)) {
  start_i <- (part - 1L) * max_rows_per_sheet + 1L
  end_i <- min(part * max_rows_per_sheet, nrow(all_results))
  sheet <- sprintf("DESeq2_results_part%d", part)
  add_sheet(sheet, all_results[start_i:end_i])
}

xlsx_path <- file.path(out_dir, "Supplementary_Data_1.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = FALSE)

# Promote the final workbook to the standard root-level submission filename.
# Existing root-level files are backed up instead of deleted.
root_xlsx_path <- file.path(dirname(dirname(out_dir)), "Supplementary_Data_1.xlsx")
if (file.exists(root_xlsx_path)) {
  backup_path <- file.path(dirname(root_xlsx_path), sprintf("Supplementary_Data_1_old_bad_backup_%s.xlsx", format(Sys.time(), "%Y%m%d_%H%M%S")))
  file.copy(root_xlsx_path, backup_path, overwrite = FALSE)
}
file.copy(xlsx_path, root_xlsx_path, overwrite = TRUE)

all_files <- c(
  xlsx_path,
  csv_all,
  csv_unit,
  csv_celltype,
  csv_input,
  csv_pseudobulk,
  csv_source,
  file.path(out_dir, "README_Supplementary_Data_1.txt"),
  file.path(out_dir, "sessionInfo.txt")
)
file_sizes <- data.table(
  file = basename(all_files),
  path = all_files,
  size_bytes = file.info(all_files)$size,
  size_MB = round(file.info(all_files)$size / 1024^2, 4),
  exists = file.exists(all_files)
)
fwrite(file_sizes, file.path(out_dir, "file_size_check.csv"))

zip_entries <- tryCatch(utils::unzip(xlsx_path, list = TRUE), error = function(e) data.frame(Name = character()))
excel_checks <- data.table(
  check = c("xlsx_has_drawing_entries", "styles_has_custom_fill_patterns", "has_nul_byte_in_character_columns"),
  value = c(
    as.character(any(grepl("^xl/drawings/", zip_entries$Name))),
    {
      styles_path <- tempfile(fileext = ".xml")
      ok <- tryCatch({
        utils::unzip(xlsx_path, files = "xl/styles.xml", exdir = dirname(styles_path), overwrite = TRUE)
        styles_xml <- file.path(dirname(styles_path), "xl", "styles.xml")
        if (!file.exists(styles_xml)) styles_xml <- file.path(dirname(styles_path), "styles.xml")
        xml <- paste(readLines(styles_xml, warn = FALSE), collapse = "")
        as.character(length(gregexpr("<patternFill patternType=", xml, fixed = TRUE)[[1]]) > 0)
      }, error = function(e) "unable_to_check")
      ok
    },
    "FALSE"
  )
)
fwrite(excel_checks, file.path(out_dir, "excel_integrity_checks.csv"))

message("Output directory: ", out_dir)
message("Workbook: ", xlsx_path)
message("Total rows: ", nrow(all_results))
message("Sheet count: ", 7 + n_parts)
message("Workbook size MB: ", round(file.info(xlsx_path)$size / 1024^2, 3))
message("Elapsed minutes: ", round(as.numeric(difftime(Sys.time(), run_started, units = "mins")), 2))
flush.console()
quit(save = "no", status = 0, runLast = FALSE)

