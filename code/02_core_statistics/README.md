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


