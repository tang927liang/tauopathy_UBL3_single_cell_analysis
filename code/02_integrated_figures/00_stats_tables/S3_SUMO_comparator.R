# ========================================================================
# 脚本名称: make_SuppTable_SUMO_NBDV1senseirevised.R  (FIXED)
#
# 目的:
#   汇总 7 个 disease/control 比较 (CTE 两条已删) 的 SUMO1/2/3 donor-level
#   comparator results, 输出 NBDV1senseirevised 版 Supplementary Table。
#
# 【本次唯一修改】read_donor_median_tbl() 里取"每 donor 中位数"那一列的方式：
#   原 guess_value_col(prefer=c("expr","median","value","log1p","cp10k")) 匹配不到
#   名为 `val` 的中位数列（"value"≠"val"），于是退到第一个数值列；而 syn21788402
#   的 `donor` 列是纯数字编号("1","2","8","9","10")会被 readr 当成数值列，结果
#   误取 donor 编号、算出 donor 编号的中位数（EC AD={8,9,10}→9，Control={1,2,3}→2）。
#   修复：排除 donor/id/计数类数值列后，优先取 val/median/expr/value，取不到再退第一列。
#   （powered 数据集不受影响：它们 donor 列非数字，本就不会被误选。）
# ========================================================================

suppressPackageStartupMessages({
  library(fs); library(readr); library(dplyr); library(purrr)
  library(stringr); library(stringi); library(tidyr); library(tibble)
  library(openxlsx); library(glue)
})

# -----------------------------
# 一、用户参数区
# -----------------------------
out_dir  <- "D:/RNA/supptable/Supp5SUMO/Results/NBDV1senseirevised"
out_xlsx <- path(out_dir,
  "Supplementary_Table_SUMO_donor_level_comparator_results_NBDV1senseirevised.xlsx")
dir_create(out_dir, recurse = TRUE)

## ★ 7 个比较 (CTE 的 H/I 两行已删)
comparison_cfg <- tribble(
  ~panel, ~dataset,       ~disease, ~region,                       ~comparison,       ~case_label, ~control_label, ~n_case_total, ~n_control_total, ~source_dir,
  "A",   "GSE157827",    "AD",     "Middle frontal gyrus",        "AD_vs_Control",   "AD",       "Control",      12,            9,                "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE157827/results",
  "B",   "GSE174367",    "AD",     "Prefrontal cortex",           "AD_vs_Control",   "AD",       "Control",      11,            7,                "D:/RNA/UBL3_AD_Project/data/sn_scRNA/GSE174367/results",
  "C",   "syn52082747",  "AD",     "Primary visual cortex (V1)",  "AD_vs_Control",   "AD",       "Control",      10,            10,               "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results",
  "D",   "syn21788402",  "AD",     "Entorhinal cortex",           "AD_vs_Control",   "AD",       "Control",      3,             3,                "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/NO2_overlap_hist_density_SUMO_GSE157827style_EC",
  "E",   "syn21788402",  "AD",     "Superior frontal gyrus",      "AD_vs_Control",   "AD",       "Control",      3,             3,                "D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/resultsmodify/NO2_overlap_hist_density_SUMO_GSE157827style_SFG",
  "F",   "syn52082747",  "FTD",    "Primary visual cortex (V1)",  "FTD_vs_Control",  "FTD",      "Control",      9,             10,               "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results",
  "G",   "syn52082747",  "PSP",    "Primary visual cortex (V1)",  "PSP_vs_Control",  "PSP",      "Control",      11,            10,               "D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/results"
)

genes <- c("SUMO1", "SUMO2", "SUMO3")
celltype_order <- c(
  "Astrocytes", "Endothelial", "Excitatory neurons",
  "Inhibitory neurons", "Microglia", "Oligodendrocytes"
)

# -----------------------------
# 二、基础辅助函数 (与原版一致)
# -----------------------------
to_utf8_safe <- function(x) {
  x <- as.character(x)
  x <- iconv(x, from = "", to = "UTF-8", sub = "")
  x <- stringi::stri_enc_toutf8(x, is_unknown_8bit = TRUE, validate = TRUE)
  x <- gsub("[[:cntrl:]]", " ", x)
  x <- gsub("\\s+", " ", x)
  x[is.na(x)] <- ""
  x
}

sanitize_df_for_excel <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.factor(col)) col <- as.character(col)
    if (inherits(col, c("POSIXct", "POSIXt"))) col <- as.character(col)
    if (is.character(col)) col <- to_utf8_safe(col)
    col
  })
  names(df) <- to_utf8_safe(names(df))
  df
}

check_bad_encoding_column <- function(df) {
  for (nm in names(df)) {
    col <- df[[nm]]
    if (is.factor(col)) col <- as.character(col)
    if (is.character(col)) {
      test <- try(stringi::stri_length(col), silent = TRUE)
      if (inherits(test, "try-error")) {
        cat("BAD COLUMN:", nm, "\n"); print(utils::head(col, 20))
        stop(paste("Encoding problem in column:", nm))
      }
    }
  }
  invisible(TRUE)
}

resolve_existing_dir <- function(root) {
  cand <- unique(c(
    root,
    str_replace(root, fixed("sn_scRNA"), "sn_RNA"),
    str_replace(root, fixed("sn_RNA"), "sn_scRNA")
  ))
  hit <- cand[dir_exists(cand)]
  if (length(hit) >= 1) {
    if (normalizePath(hit[1], winslash = "/", mustWork = TRUE) !=
        normalizePath(root, winslash = "/", mustWork = FALSE)) {
      message("路径自动修正: ", root, "  -->  ", hit[1])
    }
    return(hit[1])
  }
  stop("目录不存在 (已尝试 sn_scRNA / sn_RNA 替换): ", root)
}

find_one_file <- function(root, patterns, recurse = TRUE) {
  root <- resolve_existing_dir(root)
  files <- dir_ls(root, recurse = recurse, type = "file")
  for (pat in patterns) {
    hit <- files[str_detect(path_file(files), regex(pat, ignore_case = TRUE))]
    if (length(hit) == 1) return(hit)
    if (length(hit) > 1) {
      hit2 <- hit[str_detect(hit, regex("SUMO", ignore_case = TRUE))]
      if (length(hit2) == 1) return(hit2)
      hit3 <- hit[str_detect(hit, regex("donorMWU|donorMedian|cellHist", ignore_case = TRUE))]
      if (length(hit3) == 1) return(hit3)
      stop(glue("目录 {root} 中模式 {pat} 命中多个文件:\n{paste(hit, collapse = '\n')}"))
    }
  }
  stop(glue("目录 {root} 中找不到目标。尝试的模式:\n{paste(patterns, collapse = '\n')}"))
}

clean_group_label <- function(x) {
  x %>% as.character() %>% str_replace_all("[\r\n]+", " ") %>%
    str_replace_all("\\s+", " ") %>% str_trim()
}

extract_group_key <- function(x) {
  x <- clean_group_label(x)
  x <- str_replace(x, "\\s*\\(n\\s*=.*$", "")
  str_trim(x)
}

extract_total_n_from_group_label <- function(x) {
  x <- clean_group_label(x)
  out <- str_match(x, "\\(n\\s*=\\s*([0-9]+)\\)")[, 2]
  suppressWarnings(as.integer(out))
}

guess_group_col <- function(df) {
  cand <- names(df)[str_detect(names(df),
    regex("group|group_lab|group4|condition|disease", ignore_case = TRUE))]
  if (length(cand) == 0) stop("无法识别 group 列。列名: ", paste(names(df), collapse = ", "))
  cand[1]
}

guess_celltype_col <- function(df) {
  cand <- names(df)[str_detect(names(df),
    regex("celltype|cell_type|celltype6", ignore_case = TRUE))]
  if (length(cand) == 0) stop("无法识别 cell type 列。列名: ", paste(names(df), collapse = ", "))
  cand[1]
}

guess_value_col <- function(df, prefer = c("expr","median","value","cp10k","log1p","n_","count")) {
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  if (length(num_cols) == 0) stop("数据中无数值列。")
  num_cols2 <- setdiff(num_cols, c("x", "y"))
  if (length(num_cols2) == 0) num_cols2 <- num_cols
  for (p in prefer) {
    hit <- num_cols2[str_detect(num_cols2, regex(p, ignore_case = TRUE))]
    if (length(hit) >= 1) return(hit[1])
  }
  num_cols2[1]
}

# ★FIX 专用：稳健地选出"每 donor 中位数"那一列（排除会被当成数值的 donor/id/计数列）
pick_donor_median_col <- function(df) {
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  num_cols <- setdiff(num_cols, c("x", "y"))
  # 排除 donor 编号 / id / 计数 / 索引类数值列
  bad <- grepl("donor|^id$|_id$|index|^n$|count|cells|rank|^row",
               num_cols, ignore.case = TRUE)
  cand <- num_cols[!bad]
  if (length(cand) == 0) cand <- num_cols   # 兜底
  # 在剩余列里优先取像"中位数 / 表达值"的列
  hit <- cand[grepl("val|median|expr|value|cp10k|log1p", cand, ignore.case = TRUE)]
  vc <- if (length(hit) >= 1) hit[1] else cand[1]
  if (is.na(vc) || length(vc) == 0)
    stop("找不到 donor 中位数值列 (期望 val / median_expr 之类的数值列)。列名: ",
         paste(names(df), collapse = ", "))
  vc
}

read_csv_safely <- function(file) readr::read_csv(file, show_col_types = FALSE, progress = FALSE)

map_case_control <- function(group_values, case_label, control_label) {
  g <- extract_group_key(group_values)
  case_when(
    str_to_upper(g) == str_to_upper(case_label)    ~ "case",
    str_to_upper(g) == str_to_upper(control_label) ~ "control",
    TRUE ~ NA_character_
  )
}

# -----------------------------
# 三、各类输入文件的读取器
# -----------------------------
read_stats_tbl <- function(source_dir, gene, comparison) {
  f <- find_one_file(source_dir, c(
    glue("^STATS_{gene}_{comparison}_MWU_donorMedian_BH.*\\.csv$"),
    glue("^STATS_{gene}_{comparison}_donor.*\\.csv$"),
    glue("^STATS_{gene}_MWU_donorMedian_BH.*\\.csv$")
  ), recurse = TRUE)
  df <- read_csv_safely(f)
  cell_col  <- guess_celltype_col(df)
  p_col     <- names(df)[str_detect(names(df), regex("p_raw|pvalue|p_value|^p$", ignore_case = TRUE))][1]
  padj_col  <- names(df)[str_detect(names(df), regex("padj|fdr|adj", ignore_case = TRUE))][1]
  label_col <- names(df)[str_detect(names(df), regex("label", ignore_case = TRUE))][1]
  if (is.na(p_col) || is.na(padj_col)) stop("统计文件无法识别 p_raw / padj 列: ", f)
  df %>% transmute(
    cell_type = .data[[cell_col]],
    p_raw     = .data[[p_col]],
    padj      = .data[[padj_col]],
    stat_label = if (!is.na(label_col)) .data[[label_col]] else NA_character_
  )
}

read_n_donors_expr_tbl <- function(source_dir, gene, comparison, case_label, control_label) {
  f <- find_one_file(source_dir, c(
    glue("^CHECK_{gene}_n_donors_by_celltype_{comparison}.*\\.csv$"),
    glue("^CHECK_{gene}_n_by_celltype_{comparison}_donor.*\\.csv$"),
    glue("^CHECK_{gene}_n_donors_by_celltype_.*\\.csv$"),
    glue("^CHECK_{gene}_n_by_celltype_.*donor.*\\.csv$")
  ), recurse = TRUE)
  df <- read_csv_safely(f)
  cell_col  <- guess_celltype_col(df)
  group_col <- guess_group_col(df)
  value_col <- guess_value_col(df, prefer = c("n_donors", "exprGT0", "count", "n_"))
  df %>% transmute(
    cell_type = .data[[cell_col]],
    group_raw = .data[[group_col]],
    group_key = map_case_control(.data[[group_col]], case_label, control_label),
    n_donors_exprGT0 = .data[[value_col]],
    total_n_from_label = extract_total_n_from_group_label(.data[[group_col]])
  ) %>% filter(!is.na(group_key))
}

read_n_cells_expr_tbl <- function(source_dir, gene, comparison, case_label, control_label) {
  source_dir <- resolve_existing_dir(source_dir)
  files <- dir_ls(source_dir, recurse = TRUE, type = "file")
  possible_patterns <- c(
    glue("^CHECK_{gene}_n_cells_by_celltype_{comparison}.*\\.csv$"),
    glue("^CHECK_{gene}_n_by_celltype_{comparison}_cell.*\\.csv$"),
    glue("^CHECK_{gene}_n_cells_by_celltype_.*\\.csv$"),
    glue("^CHECK_{gene}_n_by_celltype_.*cell.*\\.csv$")
  )
  for (pat in possible_patterns) {
    hit <- files[str_detect(path_file(files), regex(pat, ignore_case = TRUE))]
    if (length(hit) == 1) {
      df <- read_csv_safely(hit)
      cell_col  <- guess_celltype_col(df)
      group_col <- guess_group_col(df)
      value_col <- guess_value_col(df, prefer = c("n_cells", "exprGT0", "count", "n_"))
      return(df %>% transmute(
        cell_type = .data[[cell_col]],
        group_key = map_case_control(.data[[group_col]], case_label, control_label),
        n_cells_exprGT0 = .data[[value_col]]
      ) %>% filter(!is.na(group_key)))
    }
  }
  f <- find_one_file(source_dir, c(
    glue("^INTERMEDIATE_{gene}_df2_{comparison}_exprGT0.*\\.csv$"),
    glue("^INTERMEDIATE_{gene}_df0_exprGT0.*\\.csv$"),
    glue("^INTERMEDIATE_{gene}_df2_.*exprGT0.*\\.csv$")
  ), recurse = TRUE)
  df <- read_csv_safely(f)
  cell_col  <- guess_celltype_col(df)
  group_col <- guess_group_col(df)
  df %>% transmute(
    cell_type = .data[[cell_col]],
    group_key = map_case_control(.data[[group_col]], case_label, control_label)
  ) %>% filter(!is.na(group_key)) %>%
    count(cell_type, group_key, name = "n_cells_exprGT0")
}

read_donor_median_tbl <- function(source_dir, gene, comparison, case_label, control_label) {
  f <- find_one_file(source_dir, c(
    glue("^INTERMEDIATE_{gene}_stat_input_{comparison}_donor.*\\.csv$"),
    glue("^INTERMEDIATE_{gene}_stat_input_.*donorMedian.*\\.csv$"),
    glue("^INTERMEDIATE_{gene}_stat_input_.*donor.*\\.csv$")
  ), recurse = TRUE)
  df <- read_csv_safely(f)
  cell_col  <- guess_celltype_col(df)
  group_col <- guess_group_col(df)
  # ★FIX: 用稳健的列选择，避免把数字编号的 donor 列当成中位数值列
  value_col <- pick_donor_median_col(df)
  df %>% transmute(
    cell_type = .data[[cell_col]],
    group_key = map_case_control(.data[[group_col]], case_label, control_label),
    donor_median_value = .data[[value_col]]
  ) %>% filter(!is.na(group_key)) %>%
    group_by(cell_type, group_key) %>%
    summarise(donor_median_group = median(donor_median_value, na.rm = TRUE), .groups = "drop")
}

# -----------------------------
# 四、单个 comparison × gene 汇总 (与原版一致)
# -----------------------------
summarise_one_comparison_gene <- function(cfg_row, gene) {
  src             <- cfg_row$source_dir[[1]]
  disease         <- cfg_row$disease[[1]]
  dataset         <- cfg_row$dataset[[1]]
  region          <- cfg_row$region[[1]]
  comparison      <- cfg_row$comparison[[1]]
  case_label      <- cfg_row$case_label[[1]]
  control_label   <- cfg_row$control_label[[1]]
  panel           <- cfg_row$panel[[1]]
  n_case_total    <- cfg_row$n_case_total[[1]]
  n_control_total <- cfg_row$n_control_total[[1]]

  stats_tbl  <- read_stats_tbl(src, gene, comparison)
  ndonor_tbl <- read_n_donors_expr_tbl(src, gene, comparison, case_label, control_label)
  ncell_tbl  <- read_n_cells_expr_tbl(src, gene, comparison, case_label, control_label)
  dmed_tbl   <- read_donor_median_tbl(src, gene, comparison, case_label, control_label)

  ndonor_wide <- ndonor_tbl %>%
    select(cell_type, group_key, n_donors_exprGT0) %>%
    pivot_wider(names_from = group_key, values_from = n_donors_exprGT0, names_prefix = "n_") %>%
    rename(n_case_donors_exprGT0 = n_case, n_control_donors_exprGT0 = n_control)

  ncell_wide <- ncell_tbl %>%
    select(cell_type, group_key, n_cells_exprGT0) %>%
    pivot_wider(names_from = group_key, values_from = n_cells_exprGT0, names_prefix = "n_") %>%
    rename(n_case_cells_exprGT0 = n_case, n_control_cells_exprGT0 = n_control)

  dmed_wide <- dmed_tbl %>%
    pivot_wider(names_from = group_key, values_from = donor_median_group, names_prefix = "median_") %>%
    rename(donor_median_case = median_case, donor_median_control = median_control)

  full_join(stats_tbl, ndonor_wide, by = "cell_type") %>%
    full_join(ncell_wide, by = "cell_type") %>%
    full_join(dmed_wide, by = "cell_type") %>%
    mutate(
      panel = panel, disease = disease, dataset = dataset, region = region,
      comparison = comparison, gene = gene,
      case_label = case_label, control_label = control_label,
      n_case_total_donors = n_case_total, n_control_total_donors = n_control_total,
      delta_case_minus_control = donor_median_case - donor_median_control,
      endpoint            = "conditional expression among gene-positive cells (UMI > 0)",
      donor_summary_stat  = "median expr per donor within cell type",
      statistical_test    = "Mann-Whitney U / Wilcoxon rank-sum test",
      multiple_testing    = "Benjamini-Hochberg across six cell types within each cohort"
    ) %>%
    select(
      panel, disease, dataset, region, comparison, gene, cell_type,
      case_label, control_label,
      n_case_total_donors, n_control_total_donors,
      n_case_donors_exprGT0, n_control_donors_exprGT0,
      n_case_cells_exprGT0, n_control_cells_exprGT0,
      donor_median_case, donor_median_control, delta_case_minus_control,
      p_raw, padj, stat_label,
      endpoint, donor_summary_stat, statistical_test, multiple_testing
    )
}

# -----------------------------
# 五、汇总所有 (7) comparison × 3 gene = 21 组合
# -----------------------------
message("正在处理 7 个 comparison × 3 SUMO 基因 = 21 组合...")
all_tbl <- map_dfr(seq_len(nrow(comparison_cfg)), function(i) {
  cfg_row <- comparison_cfg[i, , drop = FALSE]
  map_dfr(genes, ~ summarise_one_comparison_gene(cfg_row, .x))
}) %>%
  mutate(
    cell_type = factor(cell_type, levels = celltype_order),
    gene      = factor(gene, levels = genes),
    panel     = factor(panel, levels = comparison_cfg$panel)
  ) %>%
  arrange(panel, gene, cell_type) %>%
  mutate(
    cell_type = as.character(cell_type),
    gene      = as.character(gene),
    panel     = as.character(panel)
  )

all_tbl <- sanitize_df_for_excel(all_tbl)
if ("stat_label" %in% names(all_tbl)) {
  all_tbl$stat_label <- iconv(all_tbl$stat_label, from = "", to = "ASCII//TRANSLIT", sub = "")
  all_tbl$stat_label[is.na(all_tbl$stat_label)] <- ""
}
all_tbl <- sanitize_df_for_excel(all_tbl)
check_bad_encoding_column(all_tbl)

# ★FIX 的健全性自检：donor 中位数应在 log1p(CP10k) 量级（约 0–5），若仍出现 >5 的
#   "整数常数"，多半还是取错了列，提前报警。
.bad_med <- with(all_tbl, donor_median_case[is.finite(donor_median_case) & donor_median_case > 5])
if (length(.bad_med) > 0) {
  warning("⚠ 检测到 donor_median_case > 5 的异常值（疑似又取到 donor 编号列），请核对：\n  ",
          paste(round(unique(.bad_med), 2), collapse = ", "))
}

message(sprintf("✅ 总表生成: %d 行 × %d 列", nrow(all_tbl), ncol(all_tbl)))

# -----------------------------
# 六、写 Excel 工作簿
# -----------------------------
wb <- createWorkbook()

## README ("nine" → "seven")
addWorksheet(wb, "README")
readme_lines <- tibble(
  Item = c(
    "Workbook purpose",
    "Table title",
    "Number of comparisons",
    "Cohort exclusions",
    "Primary endpoint",
    "Unit of inference",
    "Statistical test",
    "Multiple-testing correction",
    "n_case_total_donors / n_control_total_donors",
    "n_case_donors_exprGT0 / n_control_donors_exprGT0",
    "n_case_cells_exprGT0 / n_control_cells_exprGT0",
    "donor_median_case / donor_median_control",
    "delta_case_minus_control",
    "Source of this table"
  ),
  Description = c(
    "Supplementary Table for SUMO1/2/3 donor-level comparator results (NBDV1 senseirevised, CTE excluded).",
    "SUMO1/2/3 donor-level comparator results across seven disease-control comparisons.",
    "Seven (7) disease-vs-control comparisons spanning AD, FTD, and PSP cohorts (CTE cohorts excluded in this revision).",
    "GSE155114 (CTE) and GSE261807 (CTE) cohorts were excluded from this revised analysis to focus the manuscript on AD and other tauopathy cohorts.",
    "Conditional expression among gene-positive cells (UMI > 0).",
    "Donor.",
    "Mann-Whitney U / Wilcoxon rank-sum test on donor-wise median expression within cell type.",
    "Benjamini-Hochberg correction across six major cell types within each cohort.",
    "Total donor counts used for the comparison, recorded according to the current manuscript / finalized cohort definition.",
    "Number of donors contributing expr > 0 cells for the specified gene and cell type.",
    "Number of expr > 0 cells included in the conditional-expression analysis for the specified gene and cell type.",
    "Median of donor-wise median expression values for case / control groups.",
    "donor_median_case - donor_median_control.",
    "This workbook is generated from intermediate/stat/check CSV files in the original dataset result directories."
  )
)
readme_lines <- sanitize_df_for_excel(readme_lines)
writeData(wb, "README", readme_lines)
setColWidths(wb, "README", cols = 1:2, widths = c(28, 110))
freezePane(wb, "README", firstActiveRow = 2)

## ALL_SUMO_long
addWorksheet(wb, "ALL_SUMO_long")
writeData(wb, "ALL_SUMO_long", all_tbl)
freezePane(wb, "ALL_SUMO_long", firstActiveRow = 2)
setColWidths(wb, "ALL_SUMO_long", cols = 1:ncol(all_tbl), widths = "auto")

## 每个 comparison 一个 sheet (7 sheets)
for (i in seq_len(nrow(comparison_cfg))) {
  rowi <- comparison_cfg[i, ]
  sheet_name <- glue("{rowi$panel}_{rowi$dataset}_{rowi$disease}")
  sheet_name <- to_utf8_safe(sheet_name)
  sheet_name <- str_sub(sheet_name, 1, 31)

  sub_tbl <- all_tbl %>%
    filter(panel == rowi$panel) %>%
    arrange(gene, factor(cell_type, levels = celltype_order))
  sub_tbl <- sanitize_df_for_excel(sub_tbl)
  check_bad_encoding_column(sub_tbl)

  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, sub_tbl)
  freezePane(wb, sheet_name, firstActiveRow = 2)
  setColWidths(wb, sheet_name, cols = 1:ncol(sub_tbl), widths = "auto")
}

## 表头样式
headerStyle <- createStyle(textDecoration = "bold", halign = "center", valign = "center")
sheet_names_written <- c(
  "README", "ALL_SUMO_long",
  vapply(seq_len(nrow(comparison_cfg)), function(i) {
    nm <- glue("{comparison_cfg$panel[i]}_{comparison_cfg$dataset[i]}_{comparison_cfg$disease[i]}")
    nm <- to_utf8_safe(nm); str_sub(nm, 1, 31)
  }, FUN.VALUE = character(1))
)
for (s in unique(sheet_names_written)) {
  try(addStyle(wb, s, headerStyle, rows = 1, cols = 1:50,
               gridExpand = TRUE, stack = TRUE), silent = TRUE)
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

# -----------------------------
# 七、汇总
# -----------------------------
message("\n========================================================")
message("✅ Supplementary Table (NBDV1senseirevised) 已生成")
message("========================================================")
message("文件: ", out_xlsx)
message("总比较数: ", nrow(comparison_cfg), " (CTE 已排除)")
message("总 sheet 数: ", length(unique(sheet_names_written)),
        " (README + ALL_SUMO_long + 7 个分 comparison)")
message("总数据行: ", nrow(all_tbl),
        " (= 7 comparison × 3 gene × 6 celltype)")
