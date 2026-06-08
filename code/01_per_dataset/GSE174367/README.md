# GSE174367 Reproduction to `stepH_obj_celltype6_named.rds`

This folder contains the minimal reproducible code for regenerating the upstream
GSE174367 Seurat object used in the manuscript:

```text
stepH_obj_celltype6_named.rds
```

The workflow starts from the author-deposited GSE174367 raw count matrix and
cell metadata, then performs QC, per-sample normalization, hierarchical CCA
integration, PCA, Louvain clustering, UMAP, cluster marker review, and final
six-class cell-type annotation.

## Required Input Files

Download the following files from GEO accession GSE174367 and place them here:

```text
data/GSE174367/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5
data/GSE174367/GSE174367_snRNA-seq_cell_meta.csv.gz
```

The raw data files are not included in this GitHub repository because they are
large public data files.

## Software Environment

Verified successful reproduction environment:

```text
R 4.5.3
Seurat 5.4.0
SeuratObject 5.3.0
Matrix 1.7-5
data.table 1.18.0
SEED = 20251023
```

The script checks for required R packages at runtime. Package-version differences
can change the binary hash of the `.rds` file, but the expected biological
outputs are the cell-level annotations and cell-type counts listed below.

## Run

From the repository root:

```r
source("01_GSE174367_rebuild_stepH_celltype6.R")
```

By default, outputs are written to:

```text
results/GSE174367_stepH/
```

You can override input and output paths before running:

```r
Sys.setenv(GSE174367_RAW_DIR = "D:/path/to/GSE174367")
Sys.setenv(GSE174367_OUT_DIR = "D:/path/to/output")
source("01_GSE174367_rebuild_stepH_celltype6.R")
```

The final output will be:

```text
results/GSE174367_stepH/stepH_obj_celltype6_named.rds
```

## Expected Final Cell Types

The final manuscript annotation is `celltype6`, containing six harmonized
classes:

```text
Astrocytes             4603
Endothelial             494
Excitatory neurons     2256
Inhibitory neurons     6451
Microglia              4309
Oligodendrocytes      37159
```

The script also keeps an audit column named `celltype7`. In that audit column,
`Peri` and `Endo` are shown separately. For the final `celltype6` annotation,
they are merged:

```text
Peri 342 + Endo 152 = Endothelial 494
```

## Included Provenance Files

```text
expected_tables/stepH_celltype6_counts_percent_GSE174367.csv
expected_tables/stepH_cluster_label_map_GSE174367.csv
expected_tables/verification_stepH_vs_historical.txt
sessionInfo.txt
```

`verification_stepH_vs_historical.txt` documents that the regenerated object in
the local manuscript rerun matched the original historical object in dimensions,
cell order, cluster counts, sample-by-celltype counts, and cell-level
`celltype6`/`celltype7` labels.

## Notes on Large Files

Do not commit generated `.rds` checkpoint files to ordinary Git. They are large
intermediate objects and are not planned for Zenodo deposition because they can
be regenerated from public inputs using this script. Deposit them in Zenodo,
Figshare, Synapse, or Git LFS only if reviewers specifically request direct
analysis-ready objects.

