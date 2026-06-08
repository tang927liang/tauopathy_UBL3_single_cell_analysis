.libPaths(c("D:/R_lib_clean", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(openxlsx)
})

run_started <- Sys.time()

base_out_dir <- "D:/RNA/2026063Molecular Neurodegeneration/Supplementary_Table_S1/results/RouteC_20260605_UBL3_donor_level"
make_safe_out_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    return(path)
  }
  existing <- list.files(path, all.files = TRUE, no.. = TRUE)
  if (length(existing) == 0) return(path)
  for (i in 2:99) {
    candidate <- paste0(path, "_v", i)
    if (!dir.exists(candidate) || length(list.files(candidate, all.files = TRUE, no.. = TRUE)) == 0) {
      dir.create(candidate, recursive = TRUE, showWarnings = FALSE)
      return(candidate)
    }
  }
  stop("Could not create a non-overwriting output directory.")
}
out_dir <- make_safe_out_dir(base_out_dir)

source_paths <- data.table(
  source_key = c("GSE157827", "GSE174367", "syn21788402_EC", "syn21788402_SFG", "syn52082747"),
  dataset = c("GSE157827", "GSE174367", "syn21788402", "syn21788402", "syn52082747"),
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
  dataset = c(
    "GSE157827", "GSE174367",
    rep("syn52082747", 3),
    "syn21788402", "syn21788402",
    rep("syn52082747", 6)
  ),
  disease = c("AD", "AD", "AD", "AD", "AD", "AD", "AD", "FTD", "FTD", "FTD", "PSP", "PSP", "PSP"),
  region = c("PFC", "PFC", "V1", "insula", "preCG", "EC", "SFG", "V1", "insula", "preCG", "V1", "insula", "preCG")
)
unit_order[, analytical_unit := sprintf("%s (%s), %s", dataset, disease, region)]

celltype_order <- c(
  "Astrocytes", "Endothelial", "Excitatory neurons",
  "Inhibitory neurons", "Microglia", "Oligodendrocytes"
)
target_gene <- "UBL3"
target_gene_ensg <- "ENSG00000122042"

pick_first_col <- function(meta, candidates, required = TRUE) {
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) > 0) return(hit[1])
  if (required) stop("Missing required metadata column among: ", paste(candidates, collapse = ", "))
  NA_character_
}

standardize_celltype <- function(x) {
  y <- trimws(as.character(x))
  yl <- tolower(y)
  out <- rep(NA_character_, length(y))
  out[grepl("astro", yl)] <- "Astrocytes"
  out[grepl("endothel", yl)] <- "Endothelial"
  out[grepl("excit", yl)] <- "Excitatory neurons"
  out[grepl("inhibit", yl)] <- "Inhibitory neurons"
  out[grepl("micro", yl)] <- "Microglia"
  out[grepl("oligo", yl)] <- "Oligodendrocytes"
  out[is.na(out) & y %in% celltype_order] <- y[is.na(out) & y %in% celltype_order]
  out
}

standardize_disease <- function(x) {
  y <- trimws(as.character(x))
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
  y <- trimws(as.character(x))
  yl <- tolower(y)
  out <- rep(NA_character_, length(y))
  out[grepl("v1|calcarine|primary visual|visual cortex", yl)] <- "V1"
  out[grepl("insula", yl)] <- "insula"
  out[grepl("precg|precentral|pre-central", yl)] <- "preCG"
  out[is.na(out) & y %in% c("V1", "insula", "preCG")] <- y[is.na(out) & y %in% c("V1", "insula", "preCG")]
  out
}

match_ubl3_row <- function(feature_names) {
  direct <- which(feature_names == target_gene)
  if (length(direct) > 0) return(list(index = direct[1], feature = feature_names[direct[1]], match_type = "symbol"))
  upper <- which(toupper(feature_names) == target_gene)
  if (length(upper) > 0) return(list(index = upper[1], feature = feature_names[upper[1]], match_type = "case_insensitive_symbol"))
  stripped <- sub("\\.\\d+$", "", feature_names)
  ensg <- which(stripped == target_gene_ensg)
  if (length(ensg) > 0) return(list(index = ensg[1], feature = feature_names[ensg[1]], match_type = "ensembl_fallback"))
  list(index = NA_integer_, feature = NA_character_, match_type = "not_found")
}

extract_ubl3_counts <- function(obj, assay = "RNA") {
  if (!assay %in% Assays(obj)) stop("Assay not found: ", assay)
  obj_cells <- colnames(obj)
  out <- setNames(rep(0, length(obj_cells)), obj_cells)
  feature_hits <- data.table(layer = character(), matched_feature = character(), match_type = character(), n_layer_cells = integer())
  assay_obj <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(assay_obj), error = function(e) character())
  count_layers <- grep("^counts", layers, value = TRUE)
  if (length(count_layers) > 0) {
    for (layer in count_layers) {
      mat <- tryCatch(
        SeuratObject::LayerData(obj, assay = assay, layer = layer),
        error = function(e) SeuratObject::LayerData(assay_obj, layer = layer)
      )
      cn <- intersect(colnames(mat), obj_cells)
      if (length(cn) == 0) next
      hit <- match_ubl3_row(rownames(mat))
      if (is.na(hit$index)) next
      out[cn] <- as.numeric(mat[hit$index, cn, drop = TRUE])
      feature_hits <- rbind(
        feature_hits,
        data.table(layer = layer, matched_feature = hit$feature, match_type = hit$match_type, n_layer_cells = length(cn))
      )
    }
  } else {
    mat <- SeuratObject::GetAssayData(obj, assay = assay, slot = "counts")
    hit <- match_ubl3_row(rownames(mat))
    if (is.na(hit$index)) stop("UBL3 was not found in RNA assay.")
    out[obj_cells] <- as.numeric(mat[hit$index, obj_cells, drop = TRUE])
    feature_hits <- rbind(
      feature_hits,
      data.table(layer = "counts", matched_feature = hit$feature, match_type = hit$match_type, n_layer_cells = length(obj_cells))
    )
  }
  attr(out, "feature_hits") <- feature_hits
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
    donor_col <- pick_first_col(meta, c("PatientID", "patient_id", "SampleID", "sample", "donor"))
    group_col <- pick_first_col(meta, c("BraakStage", "braak_stage", "group", "Diagnosis", "diagnosis"))
    region <- rep("EC", nrow(meta))
  } else if (source_key == "syn21788402_SFG") {
    donor_col <- pick_first_col(meta, c("PatientID", "patient_id", "SampleID", "sample", "donor"))
    group_col <- pick_first_col(meta, c("BraakStage", "braak_stage", "group", "Diagnosis", "diagnosis"))
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
  ncount_col <- pick_first_col(meta, c("nCount_RNA", "n_counts", "nCount", "total_counts"))
  data.table(
    cell = meta$cell,
    donor = as.character(meta[[donor_col]]),
    disease_group = standardize_disease(meta[[group_col]]),
    region = region,
    celltype = standardize_celltype(meta[[cell_col]]),
    nCount_RNA = suppressWarnings(as.numeric(meta[[ncount_col]]))
  )
}

summarize_unit <- function(meta, ubl3_counts, unit_row) {
  selected <- meta[
    region == unit_row$region &
      disease_group %in% c(unit_row$disease, "Control") &
      celltype %in% celltype_order &
      !is.na(donor) & donor != "" &
      is.finite(nCount_RNA) & nCount_RNA > 0
  ]
  if (nrow(selected) == 0) return(data.table())
  selected[, UBL3_count := as.numeric(ubl3_counts[cell])]
  selected[, UBL3_log1p_CP10K := log1p(UBL3_count / nCount_RNA * 10000)]
  selected[, UBL3_positive := UBL3_count > 0]
  selected[, .(
    analytical_unit_order = unit_row$analytical_unit_order,
    analytical_unit = unit_row$analytical_unit,
    dataset = unit_row$dataset,
    disease = unit_row$disease,
    region = unit_row$region,
    donor = donor[1],
    group = disease_group[1],
    celltype = celltype[1],
    n_cells_total = .N,
    n_UBL3_positive = sum(UBL3_positive, na.rm = TRUE),
    UBL3_positive_proportion = sum(UBL3_positive, na.rm = TRUE) / .N,
    conditional_median_UBL3_log1p_CP10K = if (sum(UBL3_positive, na.rm = TRUE) > 0) {
      median(UBL3_log1p_CP10K[UBL3_positive], na.rm = TRUE)
    } else {
      NA_real_
    }
  ), by = .(donor, group = disease_group, celltype)]
}

safe_wilcox <- function(case_values, control_values) {
  x <- case_values[is.finite(case_values)]
  y <- control_values[is.finite(control_values)]
  if (length(x) < 1 || length(y) < 1) {
    return(list(p_raw = NA_real_, hl_shift = NA_real_, hl_low = NA_real_, hl_high = NA_real_))
  }
  pair_diff <- as.vector(outer(x, y, "-"))
  hl_shift <- median(pair_diff, na.rm = TRUE)
  wt <- tryCatch(
    suppressWarnings(wilcox.test(x, y, exact = FALSE, conf.int = TRUE, conf.level = 0.95)),
    error = function(e) NULL
  )
  list(
    p_raw = if (!is.null(wt)) wt$p.value else NA_real_,
    hl_shift = hl_shift,
    hl_low = if (!is.null(wt) && !is.null(wt$conf.int) && length(wt$conf.int) == 2) unname(wt$conf.int[1]) else NA_real_,
    hl_high = if (!is.null(wt) && !is.null(wt$conf.int) && length(wt$conf.int) == 2) unname(wt$conf.int[2]) else NA_real_
  )
}

make_endpoint_stats <- function(source_dt, endpoint, value_col) {
  rows <- list()
  idx <- 1L
  for (unit_name in unit_order$analytical_unit) {
    for (ct_name in celltype_order) {
      sub <- source_dt[analytical_unit == unit_name & celltype == ct_name]
      if (nrow(sub) == 0) next
      unit_meta <- unique(sub[, .(analytical_unit_order, analytical_unit, dataset, disease, region)])
      disease_name <- unit_meta$disease[1]
      x <- sub[group == disease_name, get(value_col)]
      y <- sub[group == "Control", get(value_col)]
      wt <- safe_wilcox(x, y)
      n_d_total <- uniqueN(sub[group == disease_name, donor])
      n_c_total <- uniqueN(sub[group == "Control", donor])
      n_d_endpoint <- length(x[is.finite(x)])
      n_c_endpoint <- length(y[is.finite(y)])
      analytical_mode <- if (unit_meta$dataset[1] == "syn21788402") {
        "descriptive_only_n3v3"
      } else if (n_d_total >= 4 && n_c_total >= 4) {
        "powered_donor_level"
      } else {
        "descriptive_low_n"
      }
      endpoint_mode <- if (unit_meta$dataset[1] == "syn21788402") {
        "descriptive_only_n3v3"
      } else if (n_d_endpoint >= 4 && n_c_endpoint >= 4) {
        "powered_endpoint_donor_level"
      } else {
        "descriptive_endpoint_low_n"
      }
      rows[[idx]] <- data.table(
        endpoint = endpoint,
        analytical_unit_order = unit_meta$analytical_unit_order[1],
        analytical_unit = unit_meta$analytical_unit[1],
        dataset = unit_meta$dataset[1],
        disease = unit_meta$disease[1],
        region = unit_meta$region[1],
        celltype = ct_name,
        comparison = paste0(unit_meta$disease[1], "_vs_Control"),
        value_definition = if (endpoint == "detection_breadth") {
          "Donor-level proportion of UBL3-positive nuclei"
        } else {
          "Donor-level median UBL3 log1p(CP10K) among UBL3-positive nuclei"
        },
        n_disease_donors_total = n_d_total,
        n_control_donors_total = n_c_total,
        n_disease_donors_with_endpoint = n_d_endpoint,
        n_control_donors_with_endpoint = n_c_endpoint,
        median_disease = if (n_d_endpoint > 0) median(x, na.rm = TRUE) else NA_real_,
        median_control = if (n_c_endpoint > 0) median(y, na.rm = TRUE) else NA_real_,
        HL_shift_disease_minus_control = wt$hl_shift,
        HL_95CI_low = wt$hl_low,
        HL_95CI_high = wt$hl_high,
        p_raw = wt$p_raw,
        fdr_BH_within_unit = NA_real_,
        analytical_mode = analytical_mode,
        endpoint_mode = endpoint_mode,
        statistical_status = if (n_d_endpoint >= 1 && n_c_endpoint >= 1) "tested_donor_level" else "insufficient_endpoint_donors",
        is_focal_candidate = unit_meta$analytical_unit[1] == "syn52082747 (PSP), V1" &&
          ct_name %in% c("Excitatory neurons", "Inhibitory neurons"),
        is_microglia_comparator = unit_meta$analytical_unit[1] == "syn52082747 (PSP), V1" &&
          ct_name == "Microglia"
      )
      idx <- idx + 1L
    }
  }
  ans <- rbindlist(rows, fill = TRUE)
  ans[, celltype_order_index := match(celltype, celltype_order)]
  setorder(ans, analytical_unit_order, celltype_order_index)
  ans[, celltype_order_index := NULL]
  ans[is.finite(p_raw), fdr_BH_within_unit := p.adjust(p_raw, method = "BH"), by = .(endpoint, analytical_unit)]
  ans[]
}

manifest_rows <- list()
input_check_rows <- list()
source_rows <- list()

for (i in seq_len(nrow(source_paths))) {
  src <- source_paths[i]
  message("Reading: ", src$rds_path)
  obj <- readRDS(src$rds_path)
  DefaultAssay(obj) <- "RNA"
  meta <- prepare_metadata(obj, src$source_key)
  ubl3_counts <- extract_ubl3_counts(obj, assay = "RNA")
  feature_hits <- attr(ubl3_counts, "feature_hits")
  src_units <- unit_order[source_key == src$source_key]

  manifest_rows[[length(manifest_rows) + 1L]] <- data.table(
    source_key = src$source_key,
    dataset = src$dataset,
    rds_path = src$rds_path,
    n_cells_object = ncol(obj),
    assays = paste(Assays(obj), collapse = ";"),
    default_assay_used = "RNA",
    UBL3_matched_feature = paste(unique(feature_hits$matched_feature), collapse = ";"),
    UBL3_match_type = paste(unique(feature_hits$match_type), collapse = ";"),
    UBL3_total_counts = sum(ubl3_counts, na.rm = TRUE),
    UBL3_positive_cells = sum(ubl3_counts > 0, na.rm = TRUE),
    metadata_columns = paste(colnames(obj@meta.data), collapse = ";")
  )
  input_check_rows[[length(input_check_rows) + 1L]] <- data.table(
    source_key = src$source_key,
    dataset = src$dataset,
    matched_feature = paste(unique(feature_hits$matched_feature), collapse = ";"),
    match_type = paste(unique(feature_hits$match_type), collapse = ";"),
    feature_layers = paste(feature_hits$layer, collapse = ";"),
    total_UBL3_counts = sum(ubl3_counts, na.rm = TRUE),
    positive_cells = sum(ubl3_counts > 0, na.rm = TRUE)
  )

  for (u in seq_len(nrow(src_units))) {
    source_rows[[length(source_rows) + 1L]] <- summarize_unit(meta, ubl3_counts, src_units[u])
  }
  rm(obj, meta, ubl3_counts)
  invisible(gc())
}

donor_source <- rbindlist(source_rows, fill = TRUE)
donor_source[, celltype_order_index := match(celltype, celltype_order)]
setorder(donor_source, analytical_unit_order, celltype_order_index, group, donor)
donor_source[, celltype_order_index := NULL]

detection_source <- donor_source[, .(
  analytical_unit_order, analytical_unit, dataset, disease, region,
  celltype, donor, group,
  n_cells_total, n_UBL3_positive, UBL3_positive_proportion
)]
conditional_source <- donor_source[, .(
  analytical_unit_order, analytical_unit, dataset, disease, region,
  celltype, donor, group,
  n_cells_total, n_UBL3_positive, conditional_median_UBL3_log1p_CP10K
)]

detection_full <- make_endpoint_stats(donor_source, "detection_breadth", "UBL3_positive_proportion")
conditional_full <- make_endpoint_stats(donor_source, "conditional_expression", "conditional_median_UBL3_log1p_CP10K")

focal_detection <- detection_full[analytical_unit == "syn52082747 (PSP), V1"]
focal_conditional <- conditional_full[analytical_unit == "syn52082747 (PSP), V1"]
setnames(focal_detection, old = c(
  "median_disease", "median_control", "HL_shift_disease_minus_control",
  "HL_95CI_low", "HL_95CI_high", "p_raw", "fdr_BH_within_unit",
  "n_disease_donors_with_endpoint", "n_control_donors_with_endpoint"
), new = c(
  "detection_median_PSP", "detection_median_control", "detection_HL_shift_PSP_minus_control",
  "detection_HL_95CI_low", "detection_HL_95CI_high", "detection_p_raw", "detection_fdr_BH_within_unit",
  "detection_n_PSP_donors", "detection_n_control_donors"
), skip_absent = TRUE)
setnames(focal_conditional, old = c(
  "median_disease", "median_control", "HL_shift_disease_minus_control",
  "HL_95CI_low", "HL_95CI_high", "p_raw", "fdr_BH_within_unit",
  "n_disease_donors_with_endpoint", "n_control_donors_with_endpoint"
), new = c(
  "conditional_median_PSP", "conditional_median_control", "conditional_HL_shift_PSP_minus_control",
  "conditional_HL_95CI_low", "conditional_HL_95CI_high", "conditional_p_raw", "conditional_fdr_BH_within_unit",
  "conditional_n_PSP_positive_donors", "conditional_n_control_positive_donors"
), skip_absent = TRUE)
focal_summary <- merge(
  focal_detection[, .(
    analytical_unit, celltype, is_focal_candidate, is_microglia_comparator,
    detection_n_PSP_donors, detection_n_control_donors,
    detection_median_PSP, detection_median_control,
    detection_HL_shift_PSP_minus_control,
    detection_HL_95CI_low, detection_HL_95CI_high,
    detection_p_raw, detection_fdr_BH_within_unit
  )],
  focal_conditional[, .(
    analytical_unit, celltype,
    conditional_n_PSP_positive_donors, conditional_n_control_positive_donors,
    conditional_median_PSP, conditional_median_control,
    conditional_HL_shift_PSP_minus_control,
    conditional_HL_95CI_low, conditional_HL_95CI_high,
    conditional_p_raw, conditional_fdr_BH_within_unit
  )],
  by = c("analytical_unit", "celltype"),
  all = TRUE
)
focal_summary[, celltype_order_index := match(celltype, celltype_order)]
setorder(focal_summary, celltype_order_index)
focal_summary[, celltype_order_index := NULL]
focal_summary[, focal_interpretation := fifelse(
  is_focal_candidate &
    is.finite(detection_fdr_BH_within_unit) & detection_fdr_BH_within_unit < 0.05 &
    is.finite(detection_HL_shift_PSP_minus_control) & detection_HL_shift_PSP_minus_control > 0,
  "focal_candidate_detection_breadth_positive",
  fifelse(is_microglia_comparator, "microglia_comparator", "non_focal_PSP_V1_celltype")
)]

donor_counts <- donor_source[, .(
  n_disease_donors = uniqueN(donor[group == disease[1]]),
  n_control_donors = uniqueN(donor[group == "Control"]),
  n_disease_cells = sum(n_cells_total[group == disease[1]]),
  n_control_cells = sum(n_cells_total[group == "Control"]),
  n_disease_UBL3_positive_cells = sum(n_UBL3_positive[group == disease[1]]),
  n_control_UBL3_positive_cells = sum(n_UBL3_positive[group == "Control"])
), by = .(analytical_unit_order, analytical_unit, dataset, disease, region, celltype)]
donor_counts[, celltype_order_index := match(celltype, celltype_order)]
setorder(donor_counts, analytical_unit_order, celltype_order_index)
donor_counts[, celltype_order_index := NULL]

mode_defs <- data.table(
  analytical_mode = c("powered_donor_level", "descriptive_only_n3v3", "descriptive_low_n"),
  definition = c(
    "Disease-control analytical unit with at least four disease and four control donors; formal donor-level Wilcoxon results are used.",
    "syn21788402 EC/SFG units with n=3 vs n=3 donors; retained as descriptive donor-level comparisons.",
    "Any non-syn217 analytical unit/cell type with fewer than four donors in either group."
  )
)
manifest_dt <- rbindlist(manifest_rows, fill = TRUE)
input_check_dt <- rbindlist(input_check_rows, fill = TRUE)

powered_detection <- detection_full[analytical_mode == "powered_donor_level", .N]
descriptive_detection <- detection_full[analytical_mode != "powered_donor_level", .N]
powered_conditional <- conditional_full[endpoint_mode == "powered_endpoint_donor_level", .N]
descriptive_conditional <- conditional_full[endpoint_mode != "powered_endpoint_donor_level", .N]
conditional_sig <- conditional_full[is.finite(fdr_BH_within_unit) & fdr_BH_within_unit < 0.05]

readme <- data.table(
  field = c(
    "table_title",
    "analysis_type",
    "routeC_region_resolution",
    "syn52082747_V1_note",
    "GSE157827_region_note",
    "syn21788402_note",
    "statistical_unit",
    "endpoint_detection_breadth",
    "endpoint_conditional_expression",
    "statistical_test",
    "multiple_testing",
    "focal_candidate_findings",
    "comparator_exclusion",
    "sheet_name_note",
    "powered_detection_comparisons",
    "descriptive_detection_comparisons",
    "powered_conditional_comparisons",
    "descriptive_conditional_comparisons",
    "conditional_expression_significant_rows"
  ),
  value = c(
    "Supplementary Table S1. Donor-level UBL3 detection-breadth and conditional-expression statistics across Route C analytical units",
    "This is Route C region-resolved UBL3 donor-level analysis.",
    "syn52082747 is split into V1, insula, and preCG.",
    "V1 corresponds to the calcarine cortex metadata label in syn52082747.",
    "GSE157827 is reported as PFC.",
    "syn21788402 EC/SFG are descriptive_only_n3v3.",
    "All formal statistical tests are donor-level; cells/nuclei are not treated as independent observations.",
    "Per donor, analytical unit, and cell type: proportion of UBL3-positive nuclei.",
    "Per donor, analytical unit, and cell type: median UBL3 log1p(CP10K) among UBL3-positive nuclei only.",
    "Wilcoxon rank-sum/Mann-Whitney U test with Hodges-Lehmann disease-minus-control shift and 95% CI when available.",
    "BH FDR is computed within each endpoint and analytical unit across the six major cell types.",
    "The focal candidate findings are PSP V1 excitatory neurons and PSP V1 inhibitory neurons.",
    "No broader ubiquitin-like exploratory comparator is included in the formal submission package.",
    "Excel sheet names are limited to 31 characters; the requested UBL3_conditional_expression_full sheet is saved as UBL3_conditional_expr_full.",
    as.character(powered_detection),
    as.character(descriptive_detection),
    as.character(powered_conditional),
    as.character(descriptive_conditional),
    as.character(nrow(conditional_sig))
  )
)

session_lines <- capture.output(sessionInfo())
session_dt <- data.table(line = seq_along(session_lines), text = session_lines)

csv_paths <- c(
  detection_full = file.path(out_dir, "UBL3_detection_breadth_full.csv"),
  conditional_full = file.path(out_dir, "UBL3_conditional_expression_full.csv"),
  focal_summary = file.path(out_dir, "focal_PSP_V1_summary.csv"),
  detection_source = file.path(out_dir, "donor_level_source_data_detection_breadth.csv"),
  conditional_source = file.path(out_dir, "donor_level_source_data_conditional_expression.csv"),
  donor_counts = file.path(out_dir, "donor_counts_by_unit_celltype.csv"),
  mode_defs = file.path(out_dir, "analytical_mode_definitions.csv"),
  manifest = file.path(out_dir, "source_data_manifest.csv"),
  input_check = file.path(out_dir, "input_check.csv")
)
fwrite(detection_full, csv_paths[["detection_full"]])
fwrite(conditional_full, csv_paths[["conditional_full"]])
fwrite(focal_summary, csv_paths[["focal_summary"]])
fwrite(detection_source, csv_paths[["detection_source"]])
fwrite(conditional_source, csv_paths[["conditional_source"]])
fwrite(donor_counts, csv_paths[["donor_counts"]])
fwrite(mode_defs, csv_paths[["mode_defs"]])
fwrite(manifest_dt, csv_paths[["manifest"]])
fwrite(input_check_dt, csv_paths[["input_check"]])
writeLines(readme[, paste(field, value, sep = ": ")], file.path(out_dir, "README_Supplementary_Table_S1.txt"))
writeLines(session_lines, file.path(out_dir, "sessionInfo.txt"))

wb <- createWorkbook()
add_sheet <- function(sheet, data) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet, data)
  header_style <- createStyle(textDecoration = "bold")
  addStyle(wb, sheet, header_style, rows = 1, cols = seq_len(ncol(data)), gridExpand = TRUE)
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = seq_len(ncol(data)), widths = "auto")
}
add_sheet("README", readme)
add_sheet("UBL3_detection_breadth_full", detection_full)
add_sheet("UBL3_conditional_expr_full", conditional_full)
add_sheet("focal_PSP_V1_summary", focal_summary)
add_sheet("source_detection_breadth", detection_source)
add_sheet("source_conditional_expr", conditional_source)
add_sheet("donor_counts_by_unit_celltype", donor_counts)
add_sheet("analytical_mode_definitions", mode_defs)
add_sheet("source_data_manifest", manifest_dt)
add_sheet("session_info", session_dt)

xlsx_path <- file.path(out_dir, "Supplementary_Table_S1_UBL3_donor_level_RouteC.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = FALSE)

# Promote the final workbook to the standard root-level submission filename.
# Existing root-level files are backed up instead of deleted.
root_xlsx_path <- file.path(dirname(dirname(out_dir)), "Supplementary_Table_S1.xlsx")
if (file.exists(root_xlsx_path)) {
  backup_path <- file.path(dirname(root_xlsx_path), sprintf("Supplementary_Table_S1_old_bad_backup_%s.xlsx", format(Sys.time(), "%Y%m%d_%H%M%S")))
  file.copy(root_xlsx_path, backup_path, overwrite = FALSE)
}
file.copy(xlsx_path, root_xlsx_path, overwrite = TRUE)

all_outputs <- c(
  xlsx_path,
  unname(csv_paths),
  file.path(out_dir, "README_Supplementary_Table_S1.txt"),
  file.path(out_dir, "sessionInfo.txt")
)
file_sizes <- data.table(
  file = basename(all_outputs),
  path = all_outputs,
  size_bytes = file.info(all_outputs)$size,
  size_MB = round(file.info(all_outputs)$size / 1024^2, 4),
  exists = file.exists(all_outputs)
)
fwrite(file_sizes, file.path(out_dir, "file_size_check.csv"))

message("Output directory: ", out_dir)
message("Workbook: ", xlsx_path)
message("Focal PSP V1 summary:")
print(focal_summary)
message("Powered detection comparisons: ", powered_detection)
message("Descriptive detection comparisons: ", descriptive_detection)
message("Powered conditional comparisons: ", powered_conditional)
message("Descriptive conditional comparisons: ", descriptive_conditional)
message("Conditional significant rows: ", nrow(conditional_sig))
message("Done. Elapsed minutes: ", round(as.numeric(difftime(Sys.time(), run_started, units = "mins")), 2))


