# Pancancer_IMM
## 📁 scRNA_BulkRNA/

This directory contains scripts for integrated analysis of single-cell RNA-seq (scRNA-seq) and bulk RNA-seq (TCGA) data.  
The workflow covers scRNA-seq preprocessing and cell-type annotation, hdWGCNA co-expression network construction, and TCGA-based immune signature scoring, clustering, and visualization.

---

### 🔹 scRNA-seq analysis

- **`01_merge_scRNA.r`**  
  Merges multiple scRNA-seq cohorts and performs normalization, highly variable gene selection, and exclusion of genes not used for clustering (e.g. mitochondrial, ribosomal, stress-related genes).  
  Batch effects are corrected using Harmony, followed by UMAP dimensionality reduction, clustering, and manual annotation of major immune and stromal cell types.

- **`02_Cellular_immunity_cells.r`**  
  Performs subclustering and annotation of tumor-infiltrating lymphocytes (TILs), further resolving CD4 T, CD8 T, NK, and NKT cell subsets.  
  Processed Seurat objects for each subset are saved for downstream analyses.

- **`03_Myeloid_cells.r`**  
  Focuses on myeloid cells and performs secondary clustering and annotation, including TANs, TAMs, DCs, and MDSCs.  
  Individual Seurat objects for each myeloid subtype are generated.

---

### 🔹 hdWGCNA co-expression network analysis

- **`04_hdWGCNA.r`**  
  Constructs hdWGCNA co-expression networks for a given cell type (e.g. mast cells).  
  The analysis includes metacell construction to stabilize co-expression estimates, soft-threshold power selection, module detection, module eigengene calculation, and module connectivity analysis.  
  Module activity is quantified using UCell scores, and patient-level pseudobulk comparisons (e.g. responders vs non-responders) are performed.

- **`geom_flat_violin.R`**  
  A custom ggplot2 geometry used to generate flat violin (raincloud-style) plots for visualization of module scores.

---

### 🔹 Bulk RNA-seq (TCGA) integration

- **`05_TCGA_clustering.r`**  
  Performs pan-cancer TCGA analysis using ssGSEA to score immune-related gene signatures, including hdWGCNA-derived modules from scRNA-seq.  
  Samples are clustered into two immune states (*immune-inflamed* and *immune-depleted*) using k-means clustering, and tumor-specific heatmaps and group labels are generated.

- **`06_tSNE_3D.r`**  
  Integrates ssGSEA score matrices and immune group labels across all tumor types and performs 3D t-SNE analysis at the pan-cancer level.  
  Samples are visualized in a three-dimensional embedding space colored by immune status.

---
