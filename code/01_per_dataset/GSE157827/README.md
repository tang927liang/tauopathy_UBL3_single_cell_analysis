# GSE157827 upstream reproducibility pipeline

This repository contains the upstream R workflow used to regenerate the
manuscript-ready Seurat object:

`stepH_obj_celltype6_named.rds`

Analytical unit in the final Route C manuscript: `GSE157827 (AD), PFC`.

The script follows the original manuscript parameters for GSE157827 preprocessing,
CCA integration, PCA/UMAP/clustering, and final six-cell-type annotation.

## Files

- `01_GSE157827_rebuild_stepH_celltype6.R`: single-file upstream pipeline.
- `sessionInfo.txt`: R session information from the verified rerun environment.

Large generated objects (`*.rds`, `*.qs`) are intentionally excluded from GitHub.
They should be deposited in a research data repository such as Zenodo if they are
shared as analysis-ready data.

## Required input

Full rerun from the merged object:

```r
MERGED_RDS = "GSE157827_merged_with_group.rds"
```

Verified shortcut from an existing StepF checkpoint:

```r
STEPF_RDS = "stepF_cluster_umap_v2.rds"
START_FROM_STEPF = "true"
```

## Main parameters

- Seed: `20251023`
- QC: `nFeature_RNA > 200`, `nCount_RNA < 20000`, `percent.mt < 20`
- Normalization: `LogNormalize`, `scale.factor = 1e4`
- HVG: `FindVariableFeatures(selection.method = "vst", nfeatures = 1000)`
- Integration: hierarchical CCA, `anchor.features = 2000`, `dims = 1:20`
- PCA: `npcs = 50`
- JackStraw: `dims = 50`, `num.replicate = 100`
- Neighbors: `dims = 1:20`, `k.param = 20`
- Clustering: Louvain, `resolution = 1`, `algorithm = 1`, `n.start = 10`, `n.iter = 10`
- UMAP: `uwot`, cosine metric, `n.neighbors = 30`, `min.dist = 0.3`, `spread = 1`, spectral init

## Example run

From an existing StepF object:

```powershell
$env:START_FROM_STEPF='true'
$env:STEPF_RDS='D:/path/to/stepF_cluster_umap_v2.rds'
$env:OUTPUT_DIR='D:/path/to/output'
Rscript --vanilla 01_GSE157827_rebuild_stepH_celltype6.R
```

From the merged object:

```powershell
$env:MERGED_RDS='D:/path/to/GSE157827_merged_with_group.rds'
$env:OUTPUT_DIR='D:/path/to/output'
Rscript --vanilla 01_GSE157827_rebuild_stepH_celltype6.R
```

## Expected final cell counts

The script checks these counts before writing the final object:

| celltype6 | n_cells |
|---|---:|
| Astrocytes | 19157 |
| Endothelial | 2420 |
| Excitatory neurons | 60095 |
| Inhibitory neurons | 34821 |
| Microglia | 7562 |
| Oligodendrocytes | 45451 |

## Output

The default final output is:

```r
stepH_obj_celltype6_named.rds
```

Set `DRY_RUN_NO_SAVE=true` to verify the annotation and count checks without
writing the large final RDS.

