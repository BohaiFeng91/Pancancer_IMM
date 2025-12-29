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
load("Myeloid_cells.Rdata")

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
features <- c("CD14","FCGR3A","CD68","CD33",
              "CD163","CD86","MRC1","CD80",
              "XCR1","CCR7","IL3RA","ITGAX","ITGAM","HLA-DRA")
plot = FeaturePlot_scCustom(scRNA, features = features
                            ,num_columns=5,reduction = "umap")
ggsave("Cellular_immunity_TIL_FeaturePlot.tiff", plot, device = "tiff", width = 20, height = 20, units = "in")


Idents(scRNA)=scRNA@meta.data$RNA_snn_res.0.8
scRNA$seurat_clusters=scRNA@meta.data$RNA_snn_res.0.8
scRNA$seurat_clusters <- factor(scRNA$seurat_clusters, levels = 0:14)
scRNA@active.ident <- factor(scRNA@active.ident, levels = 0:14)
DimPlot(scRNA, reduction = "umap", cols = my_colors, label = FALSE)
DotPlot_scCustom(scRNA, features = features, scale = TRUE) + RotatedAxis()


## Manual annotation of major cell types
celltype <- c("0" = "TANs","1" = "TANs","2" = "TAMs",
              "3" = "TAMs","4" = "TAMs","5" = "TAMs",
              "6" = "MDSCs","7" = "DCs","8" = "DCs",
              "9" = "DCs", "10" = "MDSCs","11" = "DCs",
              "12" = "DCs","13" = "MDSCs","14" = "DCs")
scRNA$celltype <- unname(celltype[as.character(Idents(scRNA))])

## Save processed object
save(scRNA, file = "Myeloid_cells.Rdata")
save(subset(scRNA,celltype=="TANs"),
     file = "TANs.Rdata")
save(subset(scRNA,celltype=="TAMs"),
     file = "TAMs.Rdata")
save(subset(scRNA,celltype=="DCs"),
     file = "DCs.Rdata")
save(subset(scRNA,celltype=="MDSCs"),
     file = "MDSCs.Rdata")