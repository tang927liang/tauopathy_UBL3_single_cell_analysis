# syn21788402 — harmonized 6-cell-type objects (provenance)

## Files in this folder

- `syn21788402_preprocess.R` — **canonical creation script.** Builds the two
  harmonized objects end-to-end from the authors' public SingleCellExperiment
  objects (Synapse syn21788402, Leng et al. 2021), including the
  OPC -> Oligodendrocytes assignment. Run this to regenerate the objects.
- `syn21788402_celltype6_umap.R` — per-dataset UMAP / analysis step that consumes
  the harmonized objects.
- `README_provenance.md` — this file.

## Inputs (public; Synapse syn21788402)

- `sce_EC_scAlign_assigned.rds`  — entorhinal cortex (EC)
- `sce_SFG_scAlign_assigned.rds` — superior frontal gyrus (SFG)

Both are the authors' SingleCellExperiment objects with QC, scAlign integration,
clustering and **per-cell type annotation** already performed
(`colData: clusterCellType, clusterAssignment, SampleID, BraakStage`;
`reducedDims: CCA / CCA.ALIGNED`; assays `counts` + scran `logcounts`).

## What `syn21788402_preprocess.R` does

1. Builds a Seurat object from the authors' `counts`, using the authors' scran
   `logcounts` directly as the `data` layer (no re-normalization).
2. Attaches the `CCA.ALIGNED` embedding and runs UMAP (`seed.use = 42`, matching
   the authors' recorded command).
3. Maps the authors' `clusterCellType` (7 classes) to the 6 harmonized classes
   used across this study. **OPC -> Oligodendrocytes** (oligodendrocyte lineage:
   OPCs express PDGFRA/CSPG4/OLIG1/2 and are near-zero for GAD1/GAD2/SLC32A1);
   the original `clusterCellType` ("OPC") is kept on the object for traceability.
4. Saves the two objects.

## Outputs — canonical objects for the integrated analyses

- EC : `stepH_syn21788402_EC_obj_labeled_celltype7_celltype6.rds`
- SFG: `stepH_syn21788402_SFG_obj_celltype6.rds`

The script writes these to `.../syn21788402/resultsmodify/`. They are
intermediate objects derived from the public author objects and are not stored
in this repository. Disease/control grouping is defined by the
`stepP_*_matched_cells_meta.csv` files (Braak0 = Control n=3 vs Braak6 = AD n=3;
Braak2 excluded), which should accompany the objects when running downstream.

## Method note (relative to the other datasets)

GSE157827 / GSE174367 / syn52082747 assign the six cell types by **de-novo
cluster-marker scoring**. syn21788402 instead **reuses the authors' published
per-cluster annotation** (`clusterCellType`) mapped to the same six classes,
because this dataset is distributed already annotated. The OPC ->
Oligodendrocytes decision and the resulting celltype6-labelled object are the
same in spirit; the cell-type *source* differs and should be described as such
in the methods.
