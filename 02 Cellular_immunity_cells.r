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
load("Cellular_immunity_TILs.Rdata")

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
                                              to = 2.0, 
                                              by = 0.1))

## Visualization of major lineage markers
features <- c("CD3D","CD8A", "CD4","CCL5","PRF1",'S1PR5','FGFBP2',
              "CXCL13", "PDCD1",'LAYN','GZMB',"CCR7","XCL1", "XCL2"
              ,"FOXP3","CCR8",'RTKN2','KLRB1','CTSH','CD200',
              'BTLA',"CXCR5",'EOMES',"TCF7","NKG7","NCAM1","FCGR3A","ZNF683")

plot = FeaturePlot_scCustom(scRNA, features = features
                            ,num_columns=5,reduction = "umap")
ggsave("Cellular_immunity_TIL_FeaturePlot.tiff", plot, device = "tiff", width = 20, height = 20, units = "in")


Idents(scRNA)=scRNA@meta.data$RNA_snn_res.1.8
scRNA$seurat_clusters=scRNA@meta.data$RNA_snn_res.1.8
scRNA$seurat_clusters <- factor(scRNA$seurat_clusters, levels = 0:34)
scRNA@active.ident <- factor(scRNA@active.ident, levels = 0:34)
DimPlot(scRNA, reduction = "umap", cols = my_colors, label = FALSE)
DotPlot_scCustom(scRNA, features = features, scale = TRUE) + RotatedAxis()


## Manual annotation of major cell types
celltype <- c("0" = "CD4_T","1" = "CD4_T","2" = "CD8_T",
              "3" = "CD8_T","4" = "CD8_T","5" = "CD4_T",
              "6" = "CD8_T","7" = "CD4_T","8" = "CD4_T",
              "9" = "CD4_T","10" = "NK","11" = "CD4_T",
              "12" = "CD8_T","13" = "CD8_T","14" = "CD4_T",
              "15" = "CD8_T","16" = "CD8_T","17" = "CD8_T",
              "18" = "NKT","19" = "CD4_T","20" = "NK",
              "21" = "CD8_T","22" = "CD8_T","23" = "CD8_T",
              "24" = "CD8_T","25" = "CD8_T","26" = "CD4_T",
              "27" = "CD8_T","28" = "CD8_T","29" = "CD8_T",
              "30" = "NK","31" = "CD4_T","32" = "NK",
              "33" = "CD8_T","34" = "CD8_T")
scRNA$celltype <- unname(celltype[as.character(Idents(scRNA))])

## Save processed object
save(scRNA, file = "Cellular_immunity_TILs.Rdata")
save(subset(scRNA,celltype=="CD8_T"),
     file = "CD8_T.Rdata")
save(subset(scRNA,celltype=="TAMs"),
     file = "TAMs.Rdata")
save(subset(scRNA,celltype=="NKT"),
     file = "NKT.Rdata")
save(subset(scRNA,celltype=="NK"),
     file = "NK.Rdata")