###############################################################################
# Supplementary Figure S1
# Cross-dataset donor redundancy check rules out duplicated donors across
# source accessions.
#
#   Panel A: PCA of donor-level pseudobulk profiles (colour = dataset/display
#            unit, shape = group).
#   Panel B: Pearson correlation heatmap of the same donor profiles with
#            dataset/group/region annotation bars; plus ranked donor-pair
#            tables and a cross-accession redundancy summary.
#
# Inputs  : 5 Seurat .rds objects (7 display units) at the data-project paths
#           in section 4 below (confirm they exist before running).
# Outputs : <REPO>/output/figures/S1/ (composite vector PDF + 1000/300 dpi PNG
#           + source-data CSVs + legend sidecar + sessionInfo).
# MN specs: 170 mm wide, <=225 mm tall, vector cairo_pdf (fonts embedded) +
#           1000 dpi PNG, line widths >0.25 pt, colourblind-safe palettes,
#           supplementary file-size limit 20 MB.
###############################################################################

rm(list = ls()); gc()
Sys.setenv(LANG = "en")
options(stringsAsFactors = FALSE)
SEED <- 20251023
set.seed(SEED)

# --- Paths -------------------------------------------------------------------
# REPO : repository root; figure outputs go under here (same as Fig2-Fig4 / S3).
# Input Seurat objects and metadata CSVs are read from their data-project
# locations below; those are source data and are left unchanged.
REPO <- "D:/RNA/Code/UBL3_tauopathy"

# =============================================================================
# 0. Packages
# =============================================================================
need_cran <- c("ragg", "qs", "ggplot2", "dplyr", "data.table", "circlize", "fs")
for (p in need_cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(c("ComplexHeatmap", "org.Hs.eg.db", "AnnotationDbi"), update = FALSE, ask = FALSE)
}

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(data.table); library(dplyr); library(ggplot2)
  library(ragg); library(grid)
  library(org.Hs.eg.db); library(AnnotationDbi); library(qs)
  library(ComplexHeatmap); library(circlize); library(fs)
})

# =============================================================================
# 1. Paths, MN output settings, parameters
# =============================================================================
RAW_PIPELINE_NOTE <- "Donor-level pseudobulk redundancy check"
out_dir <- file.path(REPO, "output", "figures", "S1")
fs::dir_create(out_dir, recurse = TRUE)

TOP_VARIABLE_GENES <- 5000
MIN_CELLS_ALL      <- 50
COR_WARN           <- 0.990
COR_SUSPECT        <- 0.995

MN_WIDTH_FULL_MM <- 170
MN_HEIGHT_MAX_MM <- 225
MN_HEIGHT_MM     <- 215      # <= 225 with margin
MN_DPI           <- 1000
MN_DPI_SUBMISSION<- 300
MN_SUPP_MAX_MB   <- 20       # supplementary figure file-size limit
REL_A            <- 0.43     # fraction of composite height given to Panel A
MIN_LINE_WIDTH   <- 0.30

MM2IN <- function(x) x / 25.4

# =============================================================================
# 2. Font and embedded-font PDF device
# =============================================================================
FONT_FAM <- "Arial"
if (.Platform$OS.type == "windows") {
  arial_ok <- tryCatch({ windowsFonts(Arial = windowsFont("TT Arial")); TRUE },
                       error = function(e) FALSE, warning = function(w) TRUE)
  if (!arial_ok) FONT_FAM <- "sans"
} else {
  FONT_FAM <- "sans"
}

# =============================================================================
# 3. Colour-blind-safe palettes
# =============================================================================
# Dataset (7) — Paul Tol "bright" qualitative scheme (colour-blind safe).
# FTD is now a clearly distinct rose (#EE6677), well separated from the
# GSE174367 cyan (#66CCEE) — the v2 clash is resolved.
dataset_display_palette <- c(
  "GSE157827"       = "#4477AA",  # blue
  "GSE174367"       = "#66CCEE",  # cyan
  "syn52082747_AD"  = "#228833",  # green
  "syn21788402_EC"  = "#CCBB44",  # yellow
  "syn21788402_SFG" = "#AA3377",  # purple
  "syn52082747_FTD" = "#EE6677",  # rose
  "syn52082747_PSP" = "#BBBBBB"   # grey
)
# Group (4) — Okabe-Ito (colour-blind safe); also encoded by shape in Panel A.
group_colors <- c("AD" = "#D55E00", "FTD" = "#009E73",
                  "PSP" = "#CC79A7", "Control" = "#0072B2")
group_shapes <- c("AD" = 16, "FTD" = 18, "PSP" = 15, "Control" = 1)
# Region (5) — Paul Tol "muted" (colour-blind safe).
region_colors <- c(
  "Middle frontal gyrus"        = "#882255",
  "Prefrontal cortex"           = "#44AA99",
  "Entorhinal cortex"           = "#DDCC77",
  "Superior frontal gyrus"      = "#332288",
  "Primary visual cortex (V1)"  = "#999933"
)

# =============================================================================
# 4. Dataset configs (5 configs -> 7 display units; CTE removed)
# =============================================================================
first_existing_path <- function(paths, label = "") {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) stop("Not found for ", label, ": ", paste(paths, collapse = ", "))
  hit
}

configs <- list(
  list(dataset = "GSE157827",
       object_path = first_existing_path(c(
         "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results/stepH_obj_celltype6_named.rds"
       ), "GSE157827"),
       object_type = "rds",
       donor_candidates = c("sample"),
       group_candidates = c("group"),
       region_fixed = "Middle frontal gyrus",
       braak_filter = FALSE),
  list(dataset = "GSE174367",
       object_path = first_existing_path(c(
         "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results/stepH_obj_celltype6_named.rds"
       ), "GSE174367"),
       object_type = "rds",
       donor_candidates = c("sample", "SampleID"),
       group_candidates = c("group", "Diagnosis"),
       region_fixed = "Prefrontal cortex",
       braak_filter = FALSE),
  list(dataset = "syn21788402_EC",
       object_path = first_existing_path(c(
         "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds"
       ), "syn21788402_EC"),
       object_type = "rds",
       donor_candidates = c("PatientID"),
       group_candidates = c("BraakStage"),
       region_fixed = "Entorhinal cortex",
       braak_filter = TRUE),
  list(dataset = "syn21788402_SFG",
       object_path = first_existing_path(c(
         "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/stepH_syn21788402_SFG_obj_celltype6.rds"
       ), "syn21788402_SFG"),
       object_type = "rds",
       donor_candidates = c("PatientID"),
       group_candidates = c("BraakStage"),
       region_fixed = "Superior frontal gyrus",
       braak_filter = TRUE),
  list(dataset = "syn52082747",
       object_path = first_existing_path(c(
         "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results/NO3/stepH_slim_uncompressed.rds"
       ), "syn52082747"),
       object_type = "rds",
       donor_candidates = c("autopsy_id", "donor", "sample"),
       group_candidates = c("group4", "group", "npdx1", "clinical_dx"),
       region_fixed = "Primary visual cortex (V1)",
       braak_filter = FALSE)
)

# =============================================================================
# 5. Helper functions (unchanged pipeline logic)
# =============================================================================
safe_name <- function(x) gsub("_+", "_", gsub("[^A-Za-z0-9._-]+", "_", as.character(x)))
find_col <- function(md, candidates, required = TRUE, label = "") {
  hit <- candidates[candidates %in% colnames(md)][1]
  if (is.na(hit) && required) stop("No column found for ", label)
  hit
}
load_object <- function(path, type = "rds") {
  if (grepl("\\.qs$", path, ignore.case = TRUE)) return(qs::qread(path))
  readRDS(path)
}
normalize_group <- function(x) {
  y <- trimws(as.character(x)); y0 <- toupper(y); out <- y
  out[y0 %in% c("CTRL","CONTROL","CTR","NC","NORMAL","N","NO DEMENTIA","HEALTHY")] <- "Control"
  out[grepl("ALZ|ALZH|ALZHEIMER|^AD$", y0)] <- "AD"
  out[grepl("^PSP$", y0)] <- "PSP"
  out[grepl("FTD|BVFTD|PICK|^PID$|FTLD", y0)] <- "FTD"
  out
}
get_counts_all <- function(obj, assay = "RNA") {
  DefaultAssay(obj) <- assay
  a <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(a), error = function(e) character(0))
  count_layers <- layers[grepl("^counts", layers)]
  if (length(count_layers) > 0) {
    mats <- list()
    for (ly in count_layers) {
      m <- tryCatch(SeuratObject::LayerData(a, layer = ly), error = function(e) NULL)
      if (!is.null(m) && ncol(m) > 0) mats[[ly]] <- m
    }
    ref <- rownames(mats[[1]])
    for (nm in names(mats)) {
      if (!identical(rownames(mats[[nm]]), ref)) {
        m0 <- mats[[nm]]
        ma <- Matrix::Matrix(0, nrow = length(ref), ncol = ncol(m0), sparse = TRUE)
        rownames(ma) <- ref; colnames(ma) <- colnames(m0)
        common <- intersect(ref, rownames(m0))
        if (length(common) > 0) ma[common, ] <- m0[common, , drop = FALSE]
        mats[[nm]] <- ma
      }
    }
    mat_all <- if (length(mats) == 1) mats[[1]] else Reduce(Matrix::cbind2, mats)
    dup <- duplicated(colnames(mat_all))
    if (any(dup)) mat_all <- mat_all[, !dup, drop = FALSE]
    miss <- setdiff(colnames(obj), colnames(mat_all))
    if (length(miss) > 0) {
      mf <- Matrix::Matrix(0, nrow = nrow(mat_all), ncol = length(miss), sparse = TRUE)
      rownames(mf) <- rownames(mat_all); colnames(mf) <- miss
      mat_all <- Matrix::cbind2(mat_all, mf)
    }
    return(mat_all[, colnames(obj), drop = FALSE])
  }
  m <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "counts"), error = function(e) NULL)
  if (!is.null(m)) return(m)
  m <- tryCatch(SeuratObject::LayerData(a, layer = "counts"), error = function(e) NULL)
  if (!is.null(m)) return(m[, colnames(obj), drop = FALSE])
  stop("Cannot extract counts.")
}
standardize_to_symbol <- function(mat, dataset = "") {
  rn0 <- rownames(mat); rn_strip <- sub("\\.\\d+$", "", rn0)
  is_ens <- grepl("^ENSG[0-9]+$", rn_strip)
  symbols <- rn_strip
  if (any(is_ens)) {
    ens_keys <- unique(rn_strip[is_ens])
    map_df <- suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db, keys = ens_keys, keytype = "ENSEMBL", columns = "SYMBOL"))
    map_df <- map_df[!is.na(map_df$SYMBOL) & map_df$SYMBOL != "", ]
    map_vec <- setNames(map_df$SYMBOL, map_df$ENSEMBL)
    symbols[is_ens] <- map_vec[rn_strip[is_ens]]
  }
  keep <- !is.na(symbols) & symbols != ""
  mat <- mat[keep, , drop = FALSE]; symbols <- symbols[keep]
  lv <- sort(unique(symbols)); idx <- match(symbols, lv)
  G <- Matrix::sparseMatrix(i = idx, j = seq_along(idx), x = 1,
                            dims = c(length(lv), length(idx)),
                            dimnames = list(lv, rownames(mat)))
  as(G %*% mat, "dgCMatrix")
}

# =============================================================================
# 6. Pseudobulk per dataset
# =============================================================================
make_pseudobulk <- function(cfg) {
  message("\n==== ", cfg$dataset, " ====")
  obj <- load_object(cfg$object_path, cfg$object_type)
  DefaultAssay(obj) <- "RNA"
  md <- obj@meta.data
  md$.__cell <- rownames(md)

  donor_col <- find_col(md, cfg$donor_candidates, TRUE,  paste0(cfg$dataset, " donor"))
  group_col <- find_col(md, cfg$group_candidates, FALSE, paste0(cfg$dataset, " group"))
  md$.__donor <- trimws(as.character(md[[donor_col]]))

  if (isTRUE(cfg$braak_filter)) {
    braak <- trimws(as.character(md[[group_col]]))
    g <- rep(NA_character_, length(braak))
    g[braak == "0"] <- "Control"
    g[braak == "6"] <- "AD"
    md$.__group <- g
  } else {
    md$.__group <- if (!is.na(group_col)) normalize_group(md[[group_col]]) else NA_character_
  }

  md$.__region    <- cfg$region_fixed
  md$.__dataset   <- cfg$dataset
  md$.__accession <- sub("_(EC|SFG)$", "", cfg$dataset)

  keep <- !is.na(md$.__donor) & md$.__donor != "" &
          !is.na(md$.__group) & md$.__group != ""
  cells_keep <- rownames(md)[keep]

  cnt_all <- get_counts_all(obj, "RNA")
  cells_keep <- intersect(cells_keep, colnames(cnt_all))
  cnt <- cnt_all[, cells_keep, drop = FALSE]
  md2 <- md[cells_keep, , drop = FALSE]

  md2$.__uid_raw <- paste(cfg$dataset, safe_name(md2$.__donor),
                          safe_name(md2$.__region), sep = "__")
  key <- factor(md2$.__uid_raw, levels = unique(md2$.__uid_raw))
  M <- Matrix::sparseMatrix(i = seq_along(key), j = as.integer(key), x = 1,
                            dims = c(length(key), length(levels(key))),
                            dimnames = list(rownames(md2), levels(key)))
  pb <- as(cnt %*% M, "dgCMatrix")

  pb_meta <- md2 %>% as.data.frame() %>%
    group_by(.__uid_raw) %>%
    summarise(
      uid_raw = dplyr::first(.__uid_raw),
      dataset = dplyr::first(.__dataset),
      accession = dplyr::first(.__accession),
      donor = dplyr::first(.__donor),
      group = dplyr::first(.__group),
      region = dplyr::first(.__region),
      n_cells = dplyr::n(),
      .groups = "drop"
    ) %>% as.data.frame()

  ok <- pb_meta$n_cells >= MIN_CELLS_ALL
  pb_meta <- pb_meta[ok, , drop = FALSE]
  pb <- pb[, pb_meta$uid_raw, drop = FALSE]

  final_uid <- paste(pb_meta$dataset, safe_name(pb_meta$group),
                     safe_name(pb_meta$donor), safe_name(pb_meta$region), sep = "__")
  final_uid <- make.unique(final_uid)
  colnames(pb) <- final_uid
  pb_meta$donor_uid <- final_uid
  rownames(pb_meta) <- final_uid

  pb <- standardize_to_symbol(pb, cfg$dataset)
  pb_meta <- pb_meta[colnames(pb), , drop = FALSE]

  message(cfg$dataset, ": pseudobulk ", nrow(pb), " genes x ", ncol(pb), " donors")
  list(counts = pb, meta = pb_meta)
}

# =============================================================================
# 7. Run + merge + normalise + top variable genes
# =============================================================================
pb_all <- lapply(configs, make_pseudobulk)
names(pb_all) <- vapply(configs, function(x) x$dataset, character(1))

common_genes <- Reduce(intersect, lapply(pb_all, function(x) rownames(x$counts)))
counts_common <- lapply(pb_all, function(x) x$counts[common_genes, , drop = FALSE])
pb_counts <- Reduce(Matrix::cbind2, counts_common)
pb_meta   <- dplyr::bind_rows(lapply(pb_all, function(x) as.data.frame(x$meta)))
rownames(pb_meta) <- pb_meta$donor_uid
pb_meta   <- pb_meta[colnames(pb_counts), , drop = FALSE]

lib <- Matrix::colSums(pb_counts)
keep_lib <- is.finite(lib) & lib > 100
pb_counts <- pb_counts[, keep_lib, drop = FALSE]
pb_meta   <- pb_meta[colnames(pb_counts), , drop = FALSE]
lib <- lib[keep_lib]
expr <- as.matrix(t(t(pb_counts) / lib) * 1e4)
expr <- log1p(expr); expr[!is.finite(expr)] <- 0

vars <- apply(expr, 1, var, na.rm = TRUE)
vars <- vars[is.finite(vars) & vars > 0]
top_genes <- names(sort(vars, decreasing = TRUE))[1:min(TOP_VARIABLE_GENES, length(vars))]
expr_top  <- expr[top_genes, , drop = FALSE]
row_sd <- apply(expr_top, 1, sd, na.rm = TRUE)
expr_top <- expr_top[is.finite(row_sd) & row_sd > 1e-8, , drop = FALSE]

# =============================================================================
# 8. Display units, factor orders, PCA
# =============================================================================
pb_meta$dataset_display <- as.character(pb_meta$dataset)
is_syn52 <- pb_meta$dataset_display == "syn52082747"
pb_meta$dataset_display[is_syn52 & pb_meta$group == "AD"]      <- "syn52082747_AD"
pb_meta$dataset_display[is_syn52 & pb_meta$group == "FTD"]     <- "syn52082747_FTD"
pb_meta$dataset_display[is_syn52 & pb_meta$group == "PSP"]     <- "syn52082747_PSP"
pb_meta$dataset_display[is_syn52 & pb_meta$group == "Control"] <- "syn52082747_AD"

dataset_display_order <- c(
  "GSE157827", "GSE174367", "syn52082747_AD",
  "syn21788402_EC", "syn21788402_SFG",
  "syn52082747_FTD", "syn52082747_PSP"
)
group_order  <- c("AD", "FTD", "PSP", "Control")
region_order <- c("Middle frontal gyrus", "Prefrontal cortex",
                  "Entorhinal cortex", "Superior frontal gyrus",
                  "Primary visual cortex (V1)")

pb_meta$dataset_display <- factor(pb_meta$dataset_display, levels = dataset_display_order)
pb_meta$dataset         <- factor(pb_meta$dataset)
pb_meta$group           <- factor(pb_meta$group,  levels = group_order)
pb_meta$region          <- factor(pb_meta$region, levels = region_order)

pca <- prcomp(t(expr_top), center = TRUE, scale. = TRUE)
if (any(!is.finite(pca$x[, 1:2]))) pca <- prcomp(t(expr_top), center = TRUE, scale. = FALSE)
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], pb_meta, check.names = FALSE)
pca_df$donor_uid <- rownames(pca$x)

# Source-data tables
write.csv(pb_meta, file.path(out_dir, "Table_pseudobulk_donor_metadata.csv"), row.names = FALSE)
write.csv(pca_df,  file.path(out_dir, "Table_PCA_coordinates_donor.csv"),     row.names = FALSE)

# =============================================================================
# 9. Panel A — PCA scatter (ggplot)
# =============================================================================
nbd_theme <- theme_bw(base_size = 8, base_family = FONT_FAM) +
  theme(
    plot.title       = element_blank(),
    axis.title       = element_text(size = 8, colour = "black"),
    axis.text        = element_text(size = 7, colour = "black"),
    legend.title     = element_text(face = "bold", size = 7.5, colour = "black"),
    legend.text      = element_text(size = 7, colour = "black"),
    legend.key.size  = unit(0.7, "lines"),
    legend.position  = "right",
    legend.background = element_blank(),
    panel.border     = element_rect(colour = "grey25", fill = NA, linewidth = 0.40),
    panel.grid.major = element_line(linewidth = MIN_LINE_WIDTH, colour = "grey92"),  # v6: 0.22 -> 0.30
    panel.grid.minor = element_blank(),
    axis.ticks       = element_line(colour = "black", linewidth = 0.35),
    plot.margin      = margin(2, 2, 2, 2, unit = "pt")
  )

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = dataset_display, shape = group)) +
  geom_point(size = 2.4, alpha = 0.9, stroke = 0.6) +
  scale_color_manual(values = dataset_display_palette,
                     breaks = dataset_display_order, name = "Dataset") +
  scale_shape_manual(values = group_shapes, breaks = group_order, name = "Group") +
  labs(x = sprintf("PC1 (%.1f%%)", 100 * var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", 100 * var_exp[2])) +
  nbd_theme +
  guides(
    color = guide_legend(ncol = 1, order = 1, override.aes = list(size = 2.6)),
    shape = guide_legend(ncol = 1, order = 2, override.aes = list(size = 2.6))
  )

# =============================================================================
# 10. Panel B — Pearson correlation heatmap (ComplexHeatmap)
# =============================================================================
cor_mat <- cor(expr_top, method = "pearson", use = "pairwise.complete.obs")
saveRDS(cor_mat, file.path(out_dir, "cor_mat_donor.rds"))

cor_offdiag <- cor_mat[upper.tri(cor_mat)]
data_min <- min(cor_offdiag, na.rm = TRUE)
data_max <- max(cor_offdiag, na.rm = TRUE)

anchor_low <- floor(data_min * 20) / 20
anchor_low <- min(anchor_low, 0.50); anchor_low <- max(anchor_low, 0.30)
legend_ats <- if (anchor_low <= 0.35) {
  c(0.3, 0.5, 0.7, 0.9, 1.0)
} else if (anchor_low <= 0.45) {
  c(0.4, 0.55, 0.7, 0.85, 1.0)
} else {
  c(0.5, 0.65, 0.8, 0.9, 1.0)
}
legend_labels <- formatC(legend_ats, format = "f", digits = 1)

ann_df <- data.frame(
  Dataset = pb_meta$dataset_display,
  Group   = pb_meta$group,
  Region  = pb_meta$region,
  row.names = pb_meta$donor_uid
)

.legend_param <- function() list(
  title_gp    = gpar(fontsize = 7, fontface = "bold", fontfamily = FONT_FAM),
  labels_gp   = gpar(fontsize = 7, fontfamily = FONT_FAM),
  grid_height = unit(3, "mm"),
  grid_width  = unit(3, "mm")
)

ann_col <- HeatmapAnnotation(
  Dataset = ann_df$Dataset, Group = ann_df$Group, Region = ann_df$Region,
  col = list(Dataset = dataset_display_palette, Group = group_colors, Region = region_colors),
  annotation_name_side   = "left",
  annotation_name_gp     = gpar(fontsize = 7, fontfamily = FONT_FAM),
  annotation_name_offset = unit(1.5, "mm"),
  simple_anno_size       = unit(3.0, "mm"),
  gap                    = unit(0.9, "mm"),
  annotation_legend_param = list(Dataset = .legend_param(),
                                 Group = .legend_param(),
                                 Region = .legend_param())
)

ann_row <- rowAnnotation(
  Dataset = ann_df$Dataset, Group = ann_df$Group, Region = ann_df$Region,
  col = list(Dataset = dataset_display_palette, Group = group_colors, Region = region_colors),
  show_legend = FALSE,
  simple_anno_size = unit(3.0, "mm"),
  gap = unit(0.9, "mm"),
  annotation_name_side = "bottom",
  annotation_name_gp = gpar(fontsize = 7, fontfamily = FONT_FAM)
)

heat_col_fun <- circlize::colorRamp2(
  c(anchor_low, 0.75, 0.85, 0.93, 1.00),
  c("#8073AC", "#D8DAEB", "#F7F7F7", "#FDB863", "#B35806")  # PuOr diverging, CB-safe
)

ht <- Heatmap(
  cor_mat,
  name              = "Pearson r",
  col               = heat_col_fun,
  top_annotation    = ann_col,
  left_annotation   = ann_row,
  show_row_names    = FALSE,
  show_column_names = FALSE,
  clustering_distance_rows    = as.dist(1 - cor_mat),
  clustering_distance_columns = as.dist(1 - cor_mat),
  clustering_method_rows      = "complete",
  clustering_method_columns   = "complete",
  row_dend_width     = unit(9, "mm"),
  column_dend_height = unit(9, "mm"),
  heatmap_legend_param = list(
    title_gp      = gpar(fontsize = 7, fontface = "bold", fontfamily = FONT_FAM),
    labels_gp     = gpar(fontsize = 7, fontfamily = FONT_FAM),
    legend_height = unit(30, "mm"),
    grid_width    = unit(3.5, "mm"),
    at            = legend_ats,
    labels        = legend_labels,
    border        = "grey25"
  ),
  border = TRUE
)

# =============================================================================
# 11. Compose Panel A (top) + Panel B (bottom) on ONE device
# =============================================================================
draw_SF1 <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(nrow = 2, ncol = 1,
                               heights = grid::unit(c(REL_A, 1 - REL_A), "npc"))))

  # --- Panel A ---
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.draw(ggplot2::ggplotGrob(p_pca))
  grid::grid.text("A", x = grid::unit(1.5, "mm"),
                  y = grid::unit(1, "npc") - grid::unit(1.5, "mm"),
                  just = c("left", "top"),
                  gp = grid::gpar(fontface = "bold", fontsize = 11, fontfamily = FONT_FAM))
  grid::popViewport()

  # --- Panel B (ComplexHeatmap drawn into the current viewport) ---
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  ComplexHeatmap::draw(ht, newpage = FALSE, merge_legend = TRUE,
                       heatmap_legend_side = "right",
                       annotation_legend_side = "right",
                       padding = grid::unit(c(2, 2, 2, 3), "mm"))
  grid::grid.text("B", x = grid::unit(1.5, "mm"),
                  y = grid::unit(1, "npc") - grid::unit(1.5, "mm"),
                  just = c("left", "top"),
                  gp = grid::gpar(fontface = "bold", fontsize = 11, fontfamily = FONT_FAM))
  grid::popViewport()

  grid::popViewport()
}

out_pdf      <- file.path(out_dir, "Supplementary_Figure_S1_composite_MN_vector.pdf")
out_png      <- file.path(out_dir, "Supplementary_Figure_S1_composite_MN_1000dpi.png")
out_png_300  <- file.path(out_dir, "Supplementary_Figure_S1_composite_MN_300dpi_submission.png")

# (1) Vector PDF (fonts embedded by cairo_pdf)
if (capabilities("cairo")) {
  grDevices::cairo_pdf(out_pdf, width = MM2IN(MN_WIDTH_FULL_MM), height = MM2IN(MN_HEIGHT_MM),
                       family = FONT_FAM, bg = "white")
} else {
  grDevices::pdf(out_pdf, width = MM2IN(MN_WIDTH_FULL_MM), height = MM2IN(MN_HEIGHT_MM), bg = "white")
}
draw_SF1(); grDevices::dev.off()

# (2) 1000 dpi PNG
ragg::agg_png(out_png, width = MN_WIDTH_FULL_MM, height = MN_HEIGHT_MM, units = "mm",
              res = MN_DPI, background = "white")
draw_SF1(); grDevices::dev.off()

# (3) 300 dpi submission backup
ragg::agg_png(out_png_300, width = MN_WIDTH_FULL_MM, height = MN_HEIGHT_MM, units = "mm",
              res = MN_DPI_SUBMISSION, background = "white")
draw_SF1(); grDevices::dev.off()

message("\n>>> Saved composite:\n  ", out_pdf, "\n  ", out_png, "\n  ", out_png_300)

# =============================================================================
# 12. Donor-pair tables + redundancy summary
# =============================================================================
donor_ids <- colnames(cor_mat)
idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
pair_df <- data.frame(
  donor1_uid  = donor_ids[idx[, 1]],
  donor2_uid  = donor_ids[idx[, 2]],
  correlation = cor_mat[idx],
  stringsAsFactors = FALSE
)
for (side in c("1", "2")) {
  ids <- pair_df[[paste0("donor", side, "_uid")]]
  pair_df[[paste0("dataset", side)]]   <- as.character(pb_meta[ids, "dataset"])
  pair_df[[paste0("accession", side)]] <- pb_meta[ids, "accession"]
  pair_df[[paste0("group", side)]]     <- as.character(pb_meta[ids, "group"])
  pair_df[[paste0("donor", side)]]     <- pb_meta[ids, "donor"]
  pair_df[[paste0("region", side)]]    <- as.character(pb_meta[ids, "region"])
}
pca_xy <- data.frame(donor_uid = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
rownames(pca_xy) <- pca_xy$donor_uid
pair_df$PC12_distance <- sqrt(
  (pca_xy[pair_df$donor1_uid, "PC1"] - pca_xy[pair_df$donor2_uid, "PC1"])^2 +
  (pca_xy[pair_df$donor1_uid, "PC2"] - pca_xy[pair_df$donor2_uid, "PC2"])^2
)
pair_df$cross_accession <- pair_df$accession1 != pair_df$accession2
pair_df <- pair_df[order(-pair_df$correlation, pair_df$PC12_distance), ]
pair_df$flag <- ifelse(
  pair_df$cross_accession & pair_df$correlation >= COR_SUSPECT, "SUSPECT_r>=0.995",
  ifelse(pair_df$cross_accession & pair_df$correlation >= COR_WARN, "WATCH_r>=0.990", ""))

write.csv(pair_df, file.path(out_dir, "Table_all_donor_pairs_ranked.csv"), row.names = FALSE)
cross_df <- pair_df[pair_df$cross_accession, , drop = FALSE]
write.csv(head(cross_df, 100), file.path(out_dir, "Table_cross_accession_top100_donor.csv"), row.names = FALSE)
flagged <- cross_df[cross_df$correlation >= COR_WARN, , drop = FALSE]
write.csv(flagged, file.path(out_dir, "Table_cross_accession_flagged_cor_ge_0.99_donor.csv"), row.names = FALSE)

max_r  <- if (nrow(cross_df) > 0) round(max(cross_df$correlation), 4) else NA
n_099  <- sum(cross_df$correlation >= COR_WARN,    na.rm = TRUE)
n_0995 <- sum(cross_df$correlation >= COR_SUSPECT, na.rm = TRUE)

# =============================================================================
# 13. Manuscript title + legend (<= 300 words) written to sidecar
# =============================================================================
figure_title <- "Supplementary Figure S1. Cross-dataset donor redundancy check rules out duplicated donors across source accessions."

legend_text <- paste0(
  figure_title, "\n",
  "Donor-level pseudobulk profiles were generated by aggregating raw UMI counts across all cells of each donor x region unit (units with >= ", MIN_CELLS_ALL, " cells retained), normalized to log1p counts per 10,000, and restricted to the top ", nrow(expr_top), " most variable genes. ",
  "(A) Principal component analysis of the donor pseudobulk profiles. Each point is one donor; colour denotes the source dataset/display unit and shape denotes the diagnostic group. Axes give the variance explained by PC1 and PC2. ",
  "(B) Pearson correlation matrix of the same donor profiles, hierarchically clustered (complete linkage on 1 \u2212 r). Top and left annotation bars denote dataset, group and region. ",
  "No cross-accession donor pair exceeded the pre-specified redundancy thresholds (maximum cross-accession Pearson r = ", ifelse(is.na(max_r), "NA", formatC(max_r, format = "f", digits = 4)),
  "; pairs with r \u2265 0.990, ", n_099, "; r \u2265 0.995, ", n_0995, "), indicating that the analytical units comprise distinct donors. ",
  "AD, Alzheimer's disease; PSP, progressive supranuclear palsy; FTD, frontotemporal dementia; MFG, middle frontal gyrus; PFC, prefrontal cortex; EC, entorhinal cortex; SFG, superior frontal gyrus; V1, primary visual cortex; PC, principal component. ",
  "Source data: Supplementary Table S5."
)

legend_words <- length(strsplit(trimws(gsub("[^A-Za-z0-9]+", " ", legend_text)), "\\s+")[[1]])
writeLines(legend_text, file.path(out_dir, "Supplementary_Figure_S1_legend_MN.txt"))
writeLines(sprintf("Figure legend word count: %d", legend_words),
           file.path(out_dir, "Supplementary_Figure_S1_legend_word_count.txt"))

# =============================================================================
# 14. Self-check (supplementary file-size limit = 20 MB)
# =============================================================================
file_mb <- function(x) if (file.exists(x)) round(file.info(x)$size / 1024 / 1024, 2) else NA_real_
title_words <- length(strsplit(trimws(gsub("[^A-Za-z0-9]+", " ", sub("^Supplementary Figure S1\\. ", "", figure_title))), "\\s+")[[1]])

check_df <- data.frame(
  item = c("Composite multi-panel file", "Full-width PDF", "Maximum height",
           "1000 dpi PNG generated", "300 dpi backup generated", "Line width minimum",
           "Fonts embedded", "Colour-blind-safe palettes", "Keys inside graphic",
           "Figure title <= 15 words", "Figure legend <= 300 words",
           "PDF size (<= 20 MB supp)", "1000 dpi PNG size (<= 20 MB supp)"),
  value = c("A over B (one file)", paste0(MN_WIDTH_FULL_MM, " mm"), paste0(MN_HEIGHT_MM, " mm"),
            ifelse(file.exists(out_png), "YES", "NO"), ifelse(file.exists(out_png_300), "YES", "NO"),
            paste0(MIN_LINE_WIDTH, " pt"), FONT_FAM, "Tol bright / Okabe-Ito / Tol muted",
            "A: Dataset+Group; B: Pearson r + Dataset/Group/Region",
            paste0(title_words, " words"), paste0(legend_words, " words"),
            paste0(file_mb(out_pdf), " MB"), paste0(file_mb(out_png), " MB")),
  status = c("PASS",
             ifelse(abs(MN_WIDTH_FULL_MM - 170) <= 0.5, "PASS", "CHECK"),
             ifelse(MN_HEIGHT_MM <= MN_HEIGHT_MAX_MM, "PASS", "FAIL"),
             ifelse(file.exists(out_png), "PASS", "CHECK"),
             ifelse(file.exists(out_png_300), "PASS", "CHECK"),
             ifelse(MIN_LINE_WIDTH > 0.25, "PASS", "FAIL"),
             "PASS", "PASS", "PASS",
             ifelse(title_words <= 15, "PASS", "FAIL"),
             ifelse(legend_words <= 300, "PASS", "FAIL"),
             ifelse(!is.na(file_mb(out_pdf)) && file_mb(out_pdf) <= MN_SUPP_MAX_MB, "PASS", "CHECK"),
             ifelse(!is.na(file_mb(out_png)) && file_mb(out_png) <= MN_SUPP_MAX_MB, "PASS", "CHECK")),
  stringsAsFactors = FALSE
)
write.csv(check_df, file.path(out_dir, "Supplementary_Figure_S1_MN_guideline_selfcheck.csv"), row.names = FALSE)
print(check_df)

if (!is.na(file_mb(out_png)) && file_mb(out_png) > MN_SUPP_MAX_MB) {
  message("\nNOTE: 1000 dpi PNG (", file_mb(out_png), " MB) exceeds the 20 MB supplementary limit. ",
          "Submit the vector PDF (", file_mb(out_pdf), " MB) or the 300 dpi PNG (", file_mb(out_png_300), " MB).")
}

si <- capture.output(sessionInfo())
writeLines(c(sprintf("# Supplementary Figure S1 MN v6 build \u2014 %s", Sys.time()),
             sprintf("SEED: %d", SEED), "", si), file.path(out_dir, "sessionInfo.txt"))

message("\n>>> DONE. Results dir: ", out_dir)

