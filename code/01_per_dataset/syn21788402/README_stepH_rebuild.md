# syn21788402 StepH Rebuild

This folder contains the cleaned reproduction code for generating the two
stepH Seurat objects used by the manuscript:

- `redo/EC/stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds`
- `redo/SFG/stepH_syn21788402_SFG_obj_celltype6.rds`

The implementation is organized from the first section of
`D:/codex/21788402/redo/syn21788402 全程代码.docx`.

## Inputs

Default raw input directory:

```text
D:/RNA/UBL3_AD_Project/data/sn_scRNA/syn21788402/rawdata
```

Required files:

- `sce.EC.scAlign.assigned.rds`
- `sce.SFG.scAlign.assigned.rds`

You can override the raw-data location:

```powershell
$env:SYN21788402_RAW_DIR = "D:/path/to/rawdata"
```

## Parameters Preserved From The Original Code

- `SEED <- 42`
- `options(Seurat.object.assay.version = "v3")`
- `CreateSeuratObject(min.cells = 0, min.features = 0)`
- RNA `data` layer is set directly from SCE `logcounts`
- UMAP is rerun from `CCA.ALIGNED` using all available dimensions
- `clusterAssignment` is kept as identities and copied to `cluster_use`
- `clusterCellType` is mapped to six major cell classes:
  - `Astro -> Astrocytes`
  - `Exc -> Excitatory neurons`
  - `Micro -> Microglia`
  - `Endo -> Endothelial`
  - `Inh -> Inhibitory neurons`
  - `OPC` and `Oligo -> Oligodendrocytes`

Note: older archived objects under `D:/codex/21788402/EC` and
`D:/codex/21788402/SFG` mapped `OPC` cells to `Inhibitory neurons`. The corrected
reference objects placed in `D:/codex/21788402/redo` map `OPC` cells to
`Oligodendrocytes`, matching the manuscript methods and the first section of the
Word code. The default here follows that corrected mapping. To intentionally
reproduce the older non-method mapping, set:

```powershell
$env:SYN21788402_EC_OPC_TO = "Inhib"
$env:SYN21788402_SFG_OPC_TO = "Inhib"
```

## Run

From `D:/codex/21788402`:

```powershell
& 'C:/Program Files/R/R-4.5.3/bin/Rscript.exe' redo/scripts/01_rebuild_stepH_syn21788402.R
& 'C:/Program Files/R/R-4.5.3/bin/Rscript.exe' redo/scripts/02_verify_stepH_syn21788402.R
```

For GitHub/code-deposition use, the two target objects are also provided as
standalone one-file scripts that do not source `00_config.R`:

```powershell
& 'C:/Program Files/R/R-4.5.3/bin/Rscript.exe' redo/scripts/github_rebuild_stepH_syn21788402_EC_from_raw.R
& 'C:/Program Files/R/R-4.5.3/bin/Rscript.exe' redo/scripts/github_rebuild_stepH_syn21788402_SFG_from_raw.R
```

Each standalone script starts from the corresponding raw SCE file, copies
`counts`, `logcounts`, metadata, and `CCA.ALIGNED`, reruns UMAP with `SEED = 42`,
maps `OPC` to `Oligodendrocytes`, writes the final RDS, writes celltype QC CSVs,
and saves `sessionInfo()`.

On this workstation, the local `SeuratObject` library was built under a
different R minor version, so R can complete the work and then return a Windows
cleanup access-violation code during DLL unloading. The PowerShell wrapper
handles only that known post-completion cleanup code and still fails on normal R
errors:

```powershell
powershell -ExecutionPolicy Bypass -File redo/scripts/run_stepH_rebuild.ps1
```

The script first tries normal `SingleCellExperiment` APIs when available. On
this workstation that package is not installed, so the script uses a read-only
raw S4-slot extractor to access the same `counts`, `logcounts`, `colData`, and
`CCA.ALIGNED` matrices from the deposited SCE RDS files.

## Outputs

The rebuild script writes the target RDS files, session info, and small QC CSVs:

- `redo/sessionInfo_stepH_reconstruction.txt`
- `redo/EC/stepH_syn21788402_EC_celltype6_counts_percent.csv`
- `redo/EC/stepH_syn21788402_EC_clusterCellType_to_celltype6.csv`
- `redo/SFG/stepH_syn21788402_SFG_celltype6_counts_percent.csv`
- `redo/SFG/stepH_syn21788402_SFG_clusterCellType_to_celltype6.csv`

The verification script compares dimensions, cell order, and celltype6 counts
against the corrected reference objects in the `redo` root:

- `redo/stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds`
- `redo/stepH_syn21788402_SFG_obj_celltype6.rds`

For stricter object-level validation, run:

```powershell
& 'C:/Program Files/R/R-4.5.3/bin/Rscript.exe' redo/scripts/03_compare_rebuilt_to_correct_reference.R
```

This additionally compares key metadata columns, `clusterCellType -> celltype6`
maps, reductions, UMAP/CCA embeddings, and sparse RNA `counts`/`data` matrix
shape and sums.

The downstream Supplementary Fig. S5 SUMO-by-celltype smoke test using the
rebuilt syn21788402 objects is:

```powershell
& 'C:/Program Files/R/R-4.5.3/bin/Rscript.exe' redo/scripts/04_run_suppfigS5_SUMO_byCelltype_rebuilt_syn21788402.R
powershell -ExecutionPolicy Bypass -File redo/scripts/05_compare_suppfigS5_to_original_results2.ps1
```

Outputs are written to `redo/SuppFigS5_SUMO_byCelltype_rebuilt_test`.

The downstream Supplementary Fig. S3 reproduction using the rebuilt
syn21788402 objects is:

```powershell
& 'C:/Program Files/R/R-4.5.3/bin/Rscript.exe' redo/scripts/06_run_suppfigS3_rebuilt_syn21788402.R
powershell -ExecutionPolicy Bypass -File redo/scripts/07_compare_suppfigS3_to_original_results.ps1
& 'C:/Program Files/R/R-4.5.3/bin/Rscript.exe' redo/scripts/08_check_suppfigS3_render_equivalence.R
```

Outputs are written to `redo/SuppFigureS3_rebuilt_test`. The comparison checks
the Panel A dotplot source values, UMAP pair PNGs, legend, and final rendered
figure outputs against the manuscript `Results` folder. The original S3
1000-dpi PNG in the manuscript folder has a libpng `IDAT` checksum error on
this workstation, so strict file-hash comparison is not used as the decisive
test for that one file; the 300-dpi composite PNG, Panel A values, all UMAP
pair PNGs, and a 300-dpi render of the vector PDF are compared instead.

The Supplementary Data 1 packaging reproduction is:

```powershell
& 'C:/Program Files/R/R-4.5.3/bin/Rscript.exe' redo/scripts/09_packaging_SuppData1_rebuilt.R
```

On this workstation the R `openxlsx` package is not installed, so the same
packaging logic was also implemented as a streamed Python runner for local
verification:

```powershell
& 'C:/Users/setou/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' redo/scripts/09_packaging_suppdata1_rebuild.py --tag Submission_rebuilt_full
& 'C:/Users/setou/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' redo/scripts/10_compare_suppdata1_workbooks.py --new 'D:/codex/21788402/redo/SuppData1_rebuilt_test/Submission_rebuilt_full/Supplementary_Data_1.xlsx' --out-dir 'D:/codex/21788402/redo/SuppData1_rebuilt_test'
```

The comparison is cell-value based rather than xlsx zip-byte based, because
workbook package metadata and compression vary by writer. In the local check,
all README/cohort sheet values matched the manuscript workbook exactly, with
zero numeric difference across all sheets.
