## ------------------------------------------------------------
## hdWGCNA for celltypes (e.g. Mast cells)
## ------------------------------------------------------------
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(666)

library(Seurat)
library(harmony)
library(hdWGCNA)
library(tidyverse)
library(patchwork)
library(WGCNA)
library(UCell)
library(ggsignif)
library(scCustomize)

## Color palette (visualization only)
my_colors <- as.character(Polychrome::palette36.colors())

## If you use a custom geom, keep it in repo and source it
source("geom_flat_violin.R")

## ------------------------------------------------------------
## Input: Mast cell Seurat object
## ------------------------------------------------------------
cell <- "Mast"
load("Mast_cells.Rdata")  # should load object named scRNA (or change below accordingly)

## ------------------------------------------------------------
## Gene filtering: remove unwanted genes for network construction
## genes_to_exclude_hdWGCNA.txt should contain one gene symbol per line
## ------------------------------------------------------------
genes_to_exclude <- read.table(
  "genes_to_exclude_hdWGCNA.txt",
  header = FALSE, sep = "\t", stringsAsFactors = FALSE
)
genes_to_exclude <- unique(as.character(genes_to_exclude$V1))

## Use existing genes in object as starting set
genes_use <- setdiff(rownames(scRNA), genes_to_exclude)

## Keep only selected genes (optional but consistent with your intention)
scRNA <- subset(scRNA, features = genes_use)

## ------------------------------------------------------------
## Basic preprocessing (needed for metacells and downstream)
## ------------------------------------------------------------
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", verbose = FALSE)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", verbose = FALSE)
scRNA <- ScaleData(scRNA, verbose = FALSE)
scRNA <- RunPCA(scRNA, verbose = FALSE)

## Harmony integration for visualization/clustering (not the WGCNA network itself)
pc.num <- 1:10
scRNA <- RunHarmony(scRNA, reduction = "pca", group.by.vars = "Cohort", reduction.save = "harmony")
scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = pc.num, reduction.name = "umap")
scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = pc.num)

## Fix a single resolution for reviewer version
scRNA <- FindClusters(scRNA, resolution = 0.8)

## ------------------------------------------------------------
## hdWGCNA setup
## gene_select="fraction": keep genes expressed in >= fraction of cells
## ------------------------------------------------------------
scRNA <- SetupForWGCNA(
  scRNA,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = cell
)

## ------------------------------------------------------------
## Build metacells to stabilize co-expression estimation
## group.by: build metacells within each response x Cohort stratum
## ident.group: use response as identity grouping
## ------------------------------------------------------------
scRNA <- MetacellsByGroups(
  seurat_obj = scRNA,
  k = 30,
  max_shared = 10,
  group.by = c("response", "Cohort"),
  ident.group = "response"
)

scRNA <- NormalizeMetacells(scRNA)

## ------------------------------------------------------------
## Set expression matrix (DatExpr) for WGCNA
## Here network is constructed using RE group (as in your script)
## ------------------------------------------------------------
seurat_obj <- SetDatExpr(
  scRNA,
  group_name = "RE",
  group.by = "response",
  assay = "RNA",
  slot = "data"
)

## ------------------------------------------------------------
## Soft-threshold power selection
## ------------------------------------------------------------
seurat_obj <- TestSoftPowers(seurat_obj)
plot_list <- PlotSoftPowers(seurat_obj)

tiff("SoftPower.tiff", width = 2500, height = 1600, res = 300)
wrap_plots(plot_list, ncol = 2)
dev.off()

softpower <- 4  # chosen power in PlotSoftPowers

## ------------------------------------------------------------
## Construct network + module detection
## ------------------------------------------------------------
seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = softpower,
  setDatExpr = FALSE,
  overwrite_tom = TRUE
)

## Scale WGCNA genes + Harmony on metacells (optional; consistent with hdWGCNA workflow)
seurat_obj <- ScaleData(seurat_obj, features = GetWGCNAGenes(seurat_obj), verbose = FALSE)
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "Cohort")

seurat_obj <- ModuleEigengenes(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj)

saveRDS(seurat_obj, file = paste0("hdWGCNA_object_", cell, ".rds"))

tiff(paste0("Dendrogram_", cell, ".tiff"), width = 1500, height = 1500, res = 300)
PlotDendrogram(seurat_obj, main = paste0("hdWGCNA Dendrogram - ", cell))
dev.off()

## ------------------------------------------------------------
## Module scoring (UCell) + pseudobulk comparison (RE vs NR)
## ------------------------------------------------------------
seurat_obj <- ModuleExprScore(seurat_obj, n_genes = 100, method = "UCell")

MEs <- GetMEs(seurat_obj, harmonized = TRUE)
mods <- setdiff(colnames(MEs), "grey")

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

## Patient-level pseudobulk: average module scores per Patient_ID x response
pb_values <- seurat_obj@meta.data %>%
  group_by(Patient_ID, response) %>%
  summarise(
    across(all_of(mods), ~ mean(.x, na.rm = TRUE)),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 5)

my_comparisons <- list(c("RE", "NR"))

## Plot each module; save one file per module to avoid overwrite
for (mod in mods) {
  
  y <- pb_values[[mod]]
  
  p <- ggplot(pb_values, aes(x = response, y = y, fill = response)) +
    geom_flat_violin(position = position_nudge(x = 0.2), alpha = 0.8) +
    geom_point(aes(color = response),
               position = position_jitter(width = 0.15),
               size = 1, alpha = 0.2) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    geom_signif(comparisons = my_comparisons, step_increase = 0.1,
                map_signif_level = FALSE, test = t.test, size = 1, textsize = 6) +
    scale_fill_manual(values = my_colors[c(3, 6)]) +
    scale_colour_manual(values = my_colors[c(3, 6)]) +
    labs(
      title = paste0(cell, "_", mod),
      x = NULL,
      y = "Patient-level mean module eigengene (pseudobulk)"
    ) +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  
  ggsave(
    filename = paste0(cell, "_", mod, "_Pseudobulk.tiff"),
    plot = p,
    device = "tiff",
    width = 3.5,
    height = 4.5,
    units = "in",
    dpi = 300
  )
}