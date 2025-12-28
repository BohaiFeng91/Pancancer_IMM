## Basic settings
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(666)

## Load required packages
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(scCustomize)

## Custom color palette for visualization only
my_colors <- as.character(Polychrome::palette36.colors())

## Load merged Seurat object
load("scRNA_merge.Rdata")

## Normalization & highly variable gene selection
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 5000)
scale.genes <- VariableFeatures(scRNA)

## Exclude genes not used for clustering
## (e.g. mitochondrial, ribosomal, stress, IFN-related genes)
genes_to_exclude <- read.table("genes_to_exclude.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
genes_to_exclude <- unique(as.character(genes_to_exclude$V1))
scale.genes <- setdiff(scale.genes, genes_to_exclude)

## Scaling and PCA
scRNA <- ScaleData(scRNA, features = scale.genes, verbose = FALSE)
scRNA <- RunPCA(scRNA, features = scale.genes, verbose = FALSE)
## Assess PCA dimensionality
ElbowPlot(scRNA, reduction = "pca", ndims = 30)
## Cumulative variance explained by PCs
xx <- cumsum(scRNA[["pca"]]@stdev^2)
xx <- xx / max(xx)
which(xx > 0.9)  # PCs explaining >90% variance
pc.num <- 1:30  # chosen based on ElbowPlot and cumulative variance (>90%)

## Batch correction using Harmony
## Cohort is used as the batch variable
scRNA <- RunHarmony(scRNA, reduction = "pca", group.by.vars = "Cohort", reduction.save = "harmony")

## UMAP embedding and clustering
scRNA <- RunUMAP(scRNA, dims = pc.num, reduction = "harmony", reduction.name = "umap")
scRNA <- FindNeighbors(scRNA, dims = pc.num, reduction = "harmony")
scRNA <- FindClusters(scRNA, resolution = seq(from = 0.1, 
                                              to = 1.0, 
                                              by = 0.1))

## Visualization of major lineage markers
features <- rev(c(
  "COL1A2",
  "CD3D","CD3E","CD3G","CD8A","CD4","TRAC","FOXP3",
  "NKG7","GNLY","KLRD1",
  "CD79A","IGHM","IGHG3","IGHA2",
  "LYZ","MARCO","CD68","FCGR3A",
  "KIT","MS4A2","GATA2"
))

plot1 <- FeaturePlot_scCustom(scRNA, features = features, reduction = "umap")
ggsave("plot/total_FeaturePlot.tiff", plot1, device = "tiff", width = 16, height = 20, units = "in")

Idents(scRNA)=scRNA@meta.data$RNA_snn_res.0.8
scRNA$seurat_clusters=scRNA@meta.data$RNA_snn_res.0.8
scRNA$seurat_clusters <- factor(scRNA$seurat_clusters, levels = 0:24)
scRNA@active.ident <- factor(scRNA@active.ident, levels = 0:24)
DimPlot(scRNA, reduction = "umap", cols = my_colors, label = FALSE)
DotPlot_scCustom(scRNA, features = features, scale = TRUE) + RotatedAxis()


## Manual annotation of major cell types
celltype <- c(
  "0"="Cellular_immunity_TILs","1"="Cellular_immunity_TILs","2"="Cellular_immunity_TILs","3"="Cellular_immunity_TILs",
  "4"="B_cells","5"="Cellular_immunity_TILs","6"="Cellular_immunity_TILs","7"="Myeloid_cells",
  "8"="Cellular_immunity_TILs","9"="B_cells","10"="Cellular_immunity_TILs","11"="Myeloid_cells",
  "12"="Cellular_immunity_TILs","13"="Myeloid_cells","14"="Cellular_immunity_TILs","15"="B_cells",
  "16"="Fibroblasts","17"="Myeloid_cells","18"="Mast_cells","19"="Fibroblasts",
  "20"="Myeloid_cells","21"="B_cells","22"="B_cells","23"="B_cells","24"="B_cells"
)
scRNA$celltype <- unname(celltype[as.character(Idents(scRNA))])

## Save processed object
save(scRNA, file = "scRNA_merge_nom.Rdata")

save(subset(scRNA,celltype=="Cellular_immunity_TILs"),
     file = "Cellular_immunity_TILs.Rdata")
save(subset(scRNA,celltype=="B_cells"),
     file = "B_cells.Rdata")
save(subset(scRNA,celltype=="Mast_cells"),
     file = "Mast_cells.Rdata")
save(subset(scRNA,celltype=="Myeloid_cells"),
     file = "Myeloid_cells.Rdata")