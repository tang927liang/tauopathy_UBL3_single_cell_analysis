GSE157827 preprocessing

Dataset:
GSE157827

Disease:
Alzheimer's disease

Brain region:
Middle frontal gyrus

Main output:
stepH_obj_celltype6_named.rds

Pipeline:

Raw 10x matrices
→ Merge samples
→ QC filtering
→ Per-sample normalization
→ HVG selection (1000)
→ CCA integration (dims 1:20)
→ PCA (50 PCs)
→ Clustering (resolution 1.0)
→ UMAP
→ Marker-based cell-type annotation
→ Six-cell-type harmonized object

Original analysis environment:

R 4.5.1
Seurat 5.3.0
SeuratObject 5.1.99.9000
Matrix 1.7-4
future 1.67.0

Expected result:

169,506 post-QC cells
43 Louvain clusters

Important:

Re-running under newer Seurat versions may produce different cluster numbers.

The submitted manuscript used the original checkpoint objects generated in the environment listed above.
