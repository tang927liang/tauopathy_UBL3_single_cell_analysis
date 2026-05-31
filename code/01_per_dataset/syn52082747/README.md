# syn52082747 StepH upstream object reproduction

This folder contains the code and runtime record needed to reproduce the
upstream Seurat object used for the syn52082747 downstream analyses.

## Contents

- `make_syn52082747_stepH_slim_uncompressed.R`  
  Single all-in-one R script that regenerates `stepH_slim_uncompressed.rds`.
- `sessionInfo_stepH.txt`  
  R session information captured from the successful reproduction run.
- `README.md`  
  This file.

## Input files expected on disk

The script expects the original raw/intermediate inputs at:

```text
D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/rawdata/NO2_step7_obj_original_backup.rds
D:/RNA/UBL3_PiD_Project/data/sn_RNA/syn52082747/rawdata/liger_subcluster_metadata_v2.csv
```

## Main output

```text
D:/codex/52082747/redo/stepH_slim_uncompressed.rds
```

The successful reproduction run generated:

```text
MD5:   3f8c89bd11b467bf3c9f3de78362817f
Size:  30071151345 bytes (28.01 GiB)
Cells: 590541
Genes: 30309
Assay: RNA
Reductions: pca, tsne, umap
```

The regenerated RDS was verified to be byte-identical to the historical
manuscript object:

```text
D:/codex/52082747/stepH_slim_uncompressed.rds
```

## Manuscript Table 2 Reporting Note

The script starts from the author-provided processed Seurat object
`NO2_step7_obj_original_backup.rds`, then adds harmonized `group4`,
`celltype6`, UBL3 log1p(CP10K), audit tables, and the slim uncompressed RDS.
It does not rerun the source publication's integration, dimensional reduction,
or clustering. Therefore Table 2 should mark those fields as source/inherited
rather than this-study clustering parameters.

Recommended Table 2 values for this reconstructed object:

```text
Cells post-QC: 590,541
Integration: source Seurat/Azimuth object, inherited
Dims/PCs: 1-100 source/inherited
Resolution: 0.1 source/inherited
Clusters: 178 inherited
R / Seurat: 4.4.3 / 5.4.0
```

## Run

```bash
Rscript make_syn52082747_stepH_slim_uncompressed.R
```

By default the script both regenerates the RDS and verifies it against the
historical object if the historical file is available.

Optional environment switches:

```text
RUN_MAKE_STEPH=true|false
RUN_VERIFY_STEPH=true|false
VERIFY_HASH_RDS=true|false
VERIFY_LOAD_HISTORICAL_OBJECT=true|false
```

## Notes for archiving

The generated `.rds` file is intentionally not included in the GitHub code
repository because it is approximately 28 GiB. It should be deposited in a data
repository such as Zenodo, with the DOI cited in the manuscript Availability of
data and materials section.
