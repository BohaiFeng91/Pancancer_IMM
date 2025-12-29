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
## 📁 integrated_pathology_DL/

This directory contains scripts and notebooks for **deep learning–based pathology feature extraction and integration**
The workflow focuses on whole-slide image (WSI) feature extraction using transMIL models, as well as image-text integration with clinical and pathological information.

---

### 🔹 Pathology feature extraction (WSI-level)

- **`ctran.py`**  
  Defines the **CTransPath** backbone used for pathology feature extraction.  
  A convolutional stem is implemented and integrated into a Swin Transformer architecture (via `timm`) to generate patch-level embeddings from histopathology images.
CTransPath pretrained weights (ctranspath.phth) can be downloaded here: https://drive.google.com/file/d/1dhysqcv_Ct_A96qOF8i6COTK3jLb56vx/view
- **`Ctranspath_h5.ipynb`**  
  Uses the CTransPath model to extract patch-level features from WSIs and save them in HDF5 (`.h5`) format.  
  These features serve as inputs for downstream slide-level modeling.

- **`TransMIL_feature_extraction.ipynb`**  
  Performs feature aggregation and preprocessing for **TransMIL**, a transformer-based multiple instance learning framework.  
  Extracted features are organized at the slide level for immune-related stratification.

---

### 🔹 Clinical and multimodal representation learning

- **`Bio_ClinicalBERT.ipynb`**  
  Applies **Bio/ClinicalBERT** to encode clinical text or pathology-related annotations into dense embeddings.  
  These representations can be integrated with WSI-derived features for downstream analysis.

---

### 🔹 Purpose and integration

The scripts in this folder enable:
- Transformer-based feature extraction from histopathology patches  
- TransMIL–ready representations at the slide level  
- Integration of pathology, clinical text, and molecular features for downstream predictive modeling  

This module is designed to interface with the scRNA-seq– and TCGA-based analyses, facilitating cross-modal immune and pathology-informed studies.
If you use the Ctranspath and TransMIL in `integrated_pathology_DL/`, please cite the following papers:

Wang X, Yang S, Zhang J, Wang M, Zhang J, Yang W, Huang J, Han X. Transformer-based unsupervised contrastive learning for histopathological image classification. Med Image Anal. 2022 Oct;81:102559. doi: 10.1016/j.media.2022.102559. Epub 2022 Jul 30. PMID: 35952419.

Zhuchen Shao, Hao Bian, Yang Chen, Yifeng Wang, Jian Zhang, Xiangyang Ji, and Yongbing Zhang. 2021. TransMIL: transformer based correlated multiple instance learning for whole slide image classification. In Proceedings of the 35th International Conference on Neural Information Processing Systems (NIPS '21). Curran Associates Inc., Red Hook, NY, USA, Article 164, 2136–2147.


## 📁 Model/

This directory contains trained models and prediction utilities for downstream analysis and inference.  
It includes pretrained model files, tumor-specific model checkpoints, and notebooks for generating and evaluating predictions.

---

### 🔹 Prediction and evaluation

- **`model_prediction.ipynb`**  
A notebook for loading pretrained `.pkl` models and performing inference on new samples.  
The workflow includes:
- Loading tumor-specific models  
- Reading input feature tables  
- Generating prediction scores or class labels  
- Exporting results for downstream analysis

- **`test.csv`**  
Example input feature table used for model inference.  
Rows correspond to samples and columns correspond to model input features.

- **`prediction.csv`**  
Output file containing model prediction results generated by `model_prediction.ipynb`.

---

### 🔹 Usage notes

- Ensure that the input feature format (column names and order) matches the features used during model training.
- Select the appropriate tumor-specific `.pkl` model when performing inference.
---
