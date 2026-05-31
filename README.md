# UBL3 expression across tauopathies — analysis code

Analysis code for the study of *UBL3* expression in single-cell / single-nucleus
RNA-seq data across tauopathies (Alzheimer's disease, progressive supranuclear
palsy, and frontotemporal lobar degeneration / Pick's disease). The pipeline
harmonizes several public datasets to six common cell types and tests two
endpoints — *UBL3* **detection breadth** (fraction of expressing nuclei) and
**conditional expression** (level among expressing nuclei) — at the donor level.

This repository contains the code to reproduce every figure and statistical
table in the manuscript. Raw data are public (see **Data availability**); large
intermediate objects are not stored here.

---

## Repository structure

```
UBL3_tauopathy/
├─ code/
│  ├─ 01_per_dataset/                # per-dataset preprocessing -> harmonized celltype6 objects
│  │  ├─ GSE157827/                  # from raw 10x; AD, middle frontal gyrus
│  │  │  ├─ GSE157827_upstream_to_stepH_obj_celltype6.named.R
│  │  │  ├─ README.md
│  │  │  └─ sessionInfo.txt
│  │  ├─ GSE174367/                  # from raw h5; AD, prefrontal cortex
│  │  │  ├─ reproduce_GSE174367_stepH_obj_celltype6.named.R
│  │  │  ├─ stepH_celltype6_done_sessionInfo.txt
│  │  │  └─ README.md
│  │  ├─ syn21788402/                # EC + SFG; author SCE -> 6 classes, OPC->Oligo
│  │  │  ├─ github_rebuild_stepH_syn21788402_part1.R
│  │  │  ├─ github_rebuild_stepH_syn21788402_part2.R
│  │  │  ├─ sessionInfo_stepH_reconstruct.txt
│  │  │  └─ README_stepH_rebuild.md
│  │  └─ syn52082747/                # AD/PSP/FTD, V1; author object -> 6 classes
│  │     ├─ make_syn52082747_stepH_slim_object.R
│  │     ├─ sessionInfo_stepH.txt
│  │     └─ README.md
│  └─ 02_integrated_figures/
│     ├─ 00_stats_tables/
│     │  Data1_pseudobulk_DEseq2.R, S1_detection_breadth.R,
│     │  S2_housekeeping_baseline.R, S3_SUMO_comparator.R, S4_robustness/
│     ├─ main/
│     │  Fig1_workflow.R, Fig2_S3_celltype_UBL3.R,
│     │  Fig3_detection_breadth.R, Fig4_robustness_sensitivity.R,
│     │  Fig5_conditional_expression.R
│     └─ supplementary/
│        S1_donor_redundancy.R, S2_processing_pipeline.R,
│        S4_conditional_all_comparisons.R, S5_SUMO_comparator.R
├─ data/
│  # empty placeholders (GSE157827/ GSE174367/ syn21788402/ syn52082747/);
│  # raw data are obtained from the public accessions below
├─ output/
│  ├─ figures/        Fig1..Fig5, S1..S5  (figure scripts write here)
│  └─ stats_tables/   Data1_pseudobulk_DEseq2, S1_detection_breadth,
│                     S2_housekeeping_baseline, S3_SUMO_comparator, S4_robustness
└─ sessionInfo/        # R session logs for reproducibility
```

Note: Supplementary Figure S3 is produced together with Figure 2 by
`main/Fig2_S3_celltype_UBL3.R`, so `supplementary/` contains S1, S2, S4 and S5
only. Figure-number references inside the scripts follow the current manuscript.

---

## How to run

The two stages are independent: stage 1 produces the harmonized Seurat objects
(written to each data project's `results/` location), stage 2 reads those
objects and produces the figures and tables.

1. **Per-dataset preprocessing — `code/01_per_dataset/`**
   - `GSE157827` and `GSE174367` start from the raw matrices and run the full
     pipeline (QC → normalization → integration → clustering → 6-cell-type
     annotation). These two also underpin Supplementary Figure S2.
   - `syn52082747` starts from the authors' processed object and harmonizes
     groups + cell types.
   - `syn21788402` (EC + SFG) starts from the authors' aligned SCE objects
     (`syn21788402_preprocess.R`), reusing the authors' annotation mapped to the
     6 classes with OPC -> Oligodendrocytes; outputs go to `resultsmodify/`. See
     `syn21788402/README_provenance.md`.
   - Output of this stage: the `stepH_*` objects listed under **Data availability**.

2. **Statistical tables — `code/02_integrated_figures/00_stats_tables/`**
   - Produce the supplementary tables and Supplementary Data 1 (full-transcriptome
     pseudobulk DESeq2). Written to `output/stats_tables/`.

3. **Figures — `code/02_integrated_figures/main/` and `supplementary/`**
   - Read the `stepH_*` objects (and the stats tables for Fig 4) and write
     composite figures to `output/figures/<Fig|S>/`.

### Paths

Scripts use absolute Windows paths from the original project
(`D:/RNA/UBL3_AD_Project/...`, `D:/RNA/UBL3_PiD_Project/...`,
`D:/RNA/Code/UBL3_tauopathy`). To run elsewhere:

- Set the `REPO` variable at the top of each integrated-figure script to your
  local clone of this repository; figure outputs then go to
  `<REPO>/output/figures/<N>/`.
- Point the per-dataset input paths at your local copies of the source objects.
- Note: `syn21788402_preprocess.R` also sets an environment-specific library
  path (`clean_lib <- "D:/Rlibs/R45_seurat_clean"` via `.libPaths()`); adjust or
  remove that line to match your R installation.

---

## Data availability

All single-cell / single-nucleus RNA-seq data are public. The harmonized,
6-cell-type Seurat objects consumed by the integrated analyses are:

| Accession (repository) | Region / groups | Harmonized object |
| --- | --- | --- |
| **GSE157827** (GEO) | middle frontal gyrus; AD, Control | `stepH_obj_celltype6_named.rds` |
| **GSE174367** (GEO) | prefrontal cortex; AD, Control | `stepH_obj_celltype6_named.rds` |
| **syn21788402** (Synapse) | entorhinal cortex (EC); AD, Control | `resultsmodify/stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds` |
| **syn21788402** (Synapse) | superior frontal gyrus (SFG); AD, Control | `resultsmodify/stepH_syn21788402_SFG_obj_celltype6.rds` |
| **syn52082747** (Synapse) | primary visual cortex (V1); AD, PSP, FTD, Control | `stepH_slim_uncompressed.rds` |

The `data/` folders in this repository are intentionally empty: raw inputs are
downloaded from the accessions above, and the derived `stepH_*` objects are
intermediate files that live with the data project, not in this repository. The
syn21788402 objects are built by `syn21788402_preprocess.R` and written to
`.../syn21788402/resultsmodify/`; the integrated scripts read them from there.
See `code/01_per_dataset/syn21788402/README_provenance.md`.

---

## Environment

- Random seed `SEED <- 20251023` is fixed throughout for reproducibility.
- `GSE157827`: R 4.5.1, Seurat v5.
- `GSE174367`: R 4.4.3, Seurat 4.3.0, SeuratObject 5.2.0, Matrix 1.7-4,
  data.table 1.17.8.
- Integrated figures/tables use, among others: `Seurat`/`SeuratObject`,
  `Matrix`, `dplyr`, `ggplot2`, `ragg`, `cowplot`, `patchwork`, `lemon`,
  `ggtext`, `ggh4x`, `scattermore`, `shadowtext`, `ComplexHeatmap`, `circlize`,
  `magick`, `pdftools`, `DESeq2`, `org.Hs.eg.db`, `AnnotationDbi`, `qs`, `fs`,
  `scales`.
- Each script writes a `sessionInfo` log next to its outputs; consult those for
  exact package versions.

### Figure specifications

Figures are built to Molecular Neurodegeneration / BMC requirements: full-width
170 mm, height ≤ 225 mm, vector PDF (`cairo_pdf`, fonts embedded) plus 1000 dpi
PNG, line widths > 0.25 pt, colour-blind-safe palettes, and figure keys inside
the graphic. Figure titles and legends are written to sidecar text files for
placement in the manuscript.

---

## Citation

If you use this code, please cite the associated paper and this repository.

- Paper: *[citation to be added on acceptance]*
- Code: *[Zenodo DOI to be added]*

---

## License

Released under the MIT License — see [`LICENSE`](LICENSE).
