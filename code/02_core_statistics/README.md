# Core statistical analysis scripts

This folder contains the stable Route C analysis code for the UBL3 tauopathy
single-nucleus RNA-seq manuscript. These scripts generate the numerical source
tables and statistical outputs that support the manuscript. They are organized
by analysis module rather than by figure number so that journal layout changes
do not require rewriting the core analysis archive.

## Scripts

- `01_ubl3_donor_level_endpoints.R`
  - Builds the donor-level UBL3 detection-breadth and conditional-expression
    result tables.
  - Supports the main donor-level claims and Supplementary Table S1.

- `02_psp_v1_robustness_sensitivity.R`
  - Runs robustness and sensitivity analyses for the PSP V1 excitatory- and
    inhibitory-neuron candidate findings.
  - Supports the robustness/sensitivity workbook now reported as
    Supplementary Table S2.

- `02b_supplementary_table_s2_xlsx_cleanup.py`
  - Performs final workbook cleanup for Supplementary Table S2.

- `03_pseudobulk_deseq2_output.R`
  - Generates the complete donor-level pseudobulk DESeq2 output deposited as
    Supplementary Data 1 / data-resource output.

## Local path configuration for external reruns

The R scripts retain the author's original Windows paths as defaults, but every
large input/output path used by `code/02_core_statistics/` can be redirected with
environment variables. Set these before running the scripts on another computer:

```powershell
$env:UBL3_GSE157827_RDS="D:/your_project/processed/GSE157827_stepH_obj_celltype6_named.rds"
$env:UBL3_GSE174367_RDS="D:/your_project/processed/GSE174367_stepH_obj_celltype6_named.rds"
$env:UBL3_SYN21788402_EC_RDS="D:/your_project/processed/stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds"
$env:UBL3_SYN21788402_SFG_RDS="D:/your_project/processed/stepH_syn21788402_SFG_obj_celltype6.rds"
$env:UBL3_SYN52082747_3REGIONS_RDS="D:/your_project/processed/syn52082747_3regions_stepH_slim_uncompressed_full_seurat.rds"
$env:UBL3_ROUTE_C_BASE_DIR="D:/your_project/RouteC_outputs"
$env:UBL3_R_LIB="D:/your_R_library"
```

Additional optional output redirection variables are available for individual
workbooks:

```powershell
$env:UBL3_SUPPTABLE_S1_OUT_DIR="D:/your_project/RouteC_outputs/Supplementary_Table_S1/results"
$env:UBL3_SUPPTABLE_S2_OUT_DIR="D:/your_project/RouteC_outputs/Supplementary_Table_S2/results"
$env:UBL3_SUPPDATA1_OUT_DIR="D:/your_project/RouteC_outputs/Supplementary_Data_1/results"
$env:UBL3_FIG4_SOURCE_DIR="D:/your_project/Figure4/source_tables"
```

Example:

```powershell
Rscript --vanilla code/02_core_statistics/01_ubl3_donor_level_endpoints.R
```

If no environment variables are set, the scripts use the archived local paths
from the submitting workstation. This is intentional for provenance, but external
users should point the variables above to their own regenerated or downloaded
processed objects.
## What is intentionally not included here

Final journal-formatted figure layout scripts are not part of this stable core
analysis archive. Figure composition, panel order, labels and page layout can
change during manuscript revision without changing the underlying donor-level
statistics. The manuscript figures are generated from the source tables and
statistical outputs produced by the scripts above.

Exploratory comparator scripts from earlier manuscript drafts are not included in the formal Route C submission archive.

## Statistical unit

All inferential UBL3 analyses use donors as the statistical units. Cell-level or
nucleus-level plots in the manuscript are descriptive visualizations only.


