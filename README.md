# UBL3 tauopathy single-nucleus analysis: core reproducibility code

This repository contains the stable core reproducibility code for a
region-resolved donor-level re-analysis of public human cortical tauopathy
single-nucleus RNA-seq datasets focused on UBL3 mRNA detectability.

The final Route C analysis uses 4 public datasets and 13 disease-region
analytical units:

- GSE157827 (AD), PFC
- GSE174367 (AD), PFC
- syn21788402 (AD), EC
- syn21788402 (AD), SFG
- syn52082747 (AD), V1
- syn52082747 (AD), insula
- syn52082747 (AD), preCG
- syn52082747 (FTD/PiD), V1
- syn52082747 (FTD/PiD), insula
- syn52082747 (FTD/PiD), preCG
- syn52082747 (PSP), V1
- syn52082747 (PSP), insula
- syn52082747 (PSP), preCG

The two syn21788402 analytical units are retained as descriptive-only n = 3 vs
3 comparisons. The focal candidate finding is limited to PSP V1 excitatory and
inhibitory neurons.

## Scope

This GitHub repository provides:

1. dataset-level upstream scripts that regenerate the harmonized six-class
   Seurat objects used by the analysis;
2. stable donor-level statistical scripts and supplementary-data generation
   scripts;
3. metadata and session information needed to trace the code-to-manuscript
   relationship.

Final journal-formatted figure layout scripts are intentionally not part of
this stable core analysis archive. Figure composition, panel order, labels,
and page layout can change during revision without changing the underlying
donor-level statistics. The manuscript figures are generated from the source
tables and statistical outputs produced by `code/02_core_statistics/`.

Large Seurat RDS objects and large generated Excel/figure outputs are not
stored in GitHub. They are reproducible from public source data and local
scripts, or supplied separately as manuscript additional files / repository data
resources where file size requires it.

## Repository structure

```text
UBL3_tauopathy/
|-- README.md
|-- LICENSE
|-- .gitignore
|-- code/
|   |-- 01_per_dataset/
|   |   |-- GSE157827/
|   |   |-- GSE174367/
|   |   |-- syn21788402/
|   |   `-- syn52082747/
|   `-- 02_core_statistics/
|       |-- 01_ubl3_donor_level_endpoints.R
|       |-- 02_psp_v1_robustness_sensitivity.R
|       |-- 02b_supplementary_table_s2_xlsx_cleanup.py
|       |-- 03_pseudobulk_deseq2_output.R
|       `-- README.md
|-- metadata/
|   |-- data_sources.tsv
|   |-- processed_objects_manifest.tsv
|   `-- code_to_manuscript_map.tsv
`-- sessionInfo/
    |-- README.md
    |-- GSE157827_sessionInfo.txt
    |-- GSE174367_sessionInfo.txt
    |-- syn21788402_sessionInfo.txt
    |-- syn52082747_sessionInfo.txt
    |-- 01_ubl3_donor_level_endpoints_sessionInfo.txt
    |-- 02_psp_v1_robustness_sensitivity_sessionInfo.txt
    `-- 03_pseudobulk_deseq2_output_sessionInfo.txt
```

## Code organization

### code/01_per_dataset

These scripts regenerate or reconstruct the harmonized `celltype6` processed
objects consumed by the integrated donor-level analyses.

- `GSE157827/`: public raw counts to harmonized AD/control PFC object.
- `GSE174367/`: public raw count matrix and metadata to harmonized AD/control
  PFC object.
- `syn21788402/`: author-provided public SCE objects to EC and SFG harmonized
  objects; EC/SFG are descriptive-only n = 3 vs 3 analytical units.
- `syn52082747/`: public author-processed object plus source metadata to the final
  region-resolved full Seurat object with V1, insula, and preCG labels.

### code/02_core_statistics

These scripts generate the stable donor-level source tables and statistical
outputs underlying the manuscript:

- UBL3 donor-level detection-breadth and conditional-expression analyses.
- PSP V1 robustness and sensitivity analyses for the focal candidate signal.
- Complete pseudobulk DESeq2 output for Supplementary Data 1.

Exploratory comparator analyses from earlier drafts are not included in the formal Route C submission archive.

## Data availability

Raw sc/snRNA-seq data are publicly available through the accessions listed in
`metadata/data_sources.tsv`. The analysis-ready Seurat RDS objects are not
committed to GitHub because they are large reproducible intermediate objects.
The final syn52082747 full region-resolved Seurat object contains 590,224 cells
and is generated directly from the public author-processed object and metadata
by the script in `code/01_per_dataset/syn52082747/`.

Supplementary Data 1 is too large for ordinary GitHub storage and should be
provided through the manuscript data resource / Zenodo record described in the
manuscript.

## Local path configuration

The scripts keep the author's Windows paths as defaults so that the archived
submission run can be reproduced on the originating workstation. On another
computer, set the environment variables below before running the scripts, or
edit the corresponding variables near the top of each script.

Required processed-object paths for `code/02_core_statistics/`:

```powershell
$env:UBL3_GSE157827_RDS="D:/your_project/processed/GSE157827_stepH_obj_celltype6_named.rds"
$env:UBL3_GSE174367_RDS="D:/your_project/processed/GSE174367_stepH_obj_celltype6_named.rds"
$env:UBL3_SYN21788402_EC_RDS="D:/your_project/processed/stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds"
$env:UBL3_SYN21788402_SFG_RDS="D:/your_project/processed/stepH_syn21788402_SFG_obj_celltype6.rds"
$env:UBL3_SYN52082747_3REGIONS_RDS="D:/your_project/processed/syn52082747_3regions_stepH_slim_uncompressed_full_seurat.rds"
```

Optional output and library paths:

```powershell
$env:UBL3_ROUTE_C_BASE_DIR="D:/your_project/RouteC_outputs"
$env:UBL3_R_LIB="D:/your_R_library"
$env:UBL3_SUPPTABLE_S1_OUT_DIR="D:/your_project/RouteC_outputs/Supplementary_Table_S1/results"
$env:UBL3_SUPPTABLE_S2_OUT_DIR="D:/your_project/RouteC_outputs/Supplementary_Table_S2/results"
$env:UBL3_SUPPDATA1_OUT_DIR="D:/your_project/RouteC_outputs/Supplementary_Data_1/results"
$env:UBL3_FIG4_SOURCE_DIR="D:/your_project/Figure4/source_tables"
```

The large processed Seurat objects are not included in GitHub. They can be
regenerated with `code/01_per_dataset/` or supplied locally by setting the RDS
path variables above. The journal-formatted figure layout scripts are separate
from this stable core archive; the core scripts generate the numerical source
tables used by the manuscript figures and supplementary files.
## Reproducibility notes

- Random seeds are declared in each script.
- R session information from successful runs is stored in `sessionInfo/`.
- Processed-object provenance is recorded in
  `metadata/processed_objects_manifest.tsv`.
- The code-to-manuscript mapping is recorded in
  `metadata/code_to_manuscript_map.tsv`.

## License

Released under the MIT License. See `LICENSE`.


