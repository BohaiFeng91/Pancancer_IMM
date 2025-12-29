# ============================================================
# TCGA pan-cancer: ssGSEA immune signatures + 2-class clustering
# Outputs per tumor: heatmap.tiff, Group.rds, Matrix.rds
# Input structure:
#   TCGA/<TumorType>/mRNA_FPKM.txt   (genes x samples, with a column "Tag")
# ============================================================

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GSVA)
  library(limma)
  library(pheatmap)
  library(Polychrome)
})

set.seed(666)
my_colors <- as.character(Polychrome::palette36.colors())

# -------- User-defined signatures (EDIT HERE) --------
# IMPORTANT: define TIL_gs / TLS_gs at the beginning (consistent naming downstream)
IMM_genelist <- list(
  CD274  = c("CD274"),
  CD8    = c("CD8A", "CD8B"),
  # TIL_gs = c("..."),
  # TLS_gs = c("..."),
  # You can append hdWGCNA modules here as additional gene sets, e.g.:
  # B_blue = c("GENE1","GENE2"), Mast_turquoise = c("GENE3","GENE4"), ...
)

# Order of signatures to show (only keep those that exist in IMM_genelist)
col_order <- c(
  "CD274","CD8","TIL_gs","TLS_gs",
  "B_blue","Mast_turquoise","TANs_turquoise",
  "B_green","B_red",
  "CD4_blue","CD8_blue","CD8_yellow",
  "NK_blue","NK_brown","NKT_blue",
  "MDSCs_yellow"
)

# Composite score features for defining immune-inflamed vs immune-depleted (EDIT if needed)
score_features <- c("CD274","CD8","TIL_gs","TLS_gs","B_blue","Mast_turquoise","TANs_turquoise")

# -------- IO settings --------
inputDir  <- "TCGA"
outputDir <- "TCGA_pancancer"
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
tumor_types <- list.dirs(inputDir, full.names = FALSE, recursive = FALSE)

# -------- Helpers --------
normalize <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (isTRUE(all.equal(rng[1], rng[2]))) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

read_tcga_fpkm <- function(path) {
  rt <- fread(path, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
  stopifnot("Tag" %in% colnames(rt))
  gene <- rt$Tag
  rt <- dplyr::select(rt, ends_with("-01"))  # primary tumor samples
  mat <- as.matrix(rt)
  rownames(mat) <- gene
  storage.mode(mat) <- "numeric"
  mat <- avereps(mat)                        # collapse duplicate genes
  mat <- mat[rowMeans(mat, na.rm = TRUE) > 0, , drop = FALSE]
  mat
}

# ============================================================
# Main loop
# ============================================================
for (tumor in tumor_types) {
  message("Processing ", tumor)
  
  inputFile <- file.path(inputDir, tumor, "mRNA_FPKM.txt")
  if (!file.exists(inputFile)) {
    message("  File not found: ", inputFile)
    next
  }
  
  # 1) Read expression (genes x samples)
  mat <- read_tcga_fpkm(inputFile)
  
  # 2) ssGSEA
  # NOTE: gsva expects gene sets as list(name -> character vector of genes)
  ssgseaScore <- gsva(
    expr = mat,
    gset.idx.list = IMM_genelist,
    method = "ssgsea",
    kcdf = "Gaussian",
    abs.ranking = TRUE
  )
  
  # 3) Normalize each signature to [0,1] across samples (row-wise in ssgseaScore)
  ssgsea <- t(apply(ssgseaScore, 1, normalize))
  data <- as.data.frame(t(ssgsea))  # samples x signatures
  
  # Keep and order columns (robust to missing gene sets)
  keep_cols <- intersect(col_order, colnames(data))
  data <- data[, keep_cols, drop = FALSE]
  
  # 4) K-means clustering (2 classes) on selected signatures
  km <- kmeans(data, centers = 2, nstart = 50, algorithm = "MacQueen")
  ann <- data.frame(cluster = km$cluster)
  rownames(ann) <- rownames(data)
  
  # 5) Decide which cluster is "immune-inflamed" using composite mean score
  feats <- intersect(score_features, colnames(data))
  if (length(feats) == 0) stop("None of score_features found in data columns.")
  
  composite <- rowMeans(data[, feats, drop = FALSE], na.rm = TRUE)  # per sample
  composite_by_cluster <- tapply(composite, ann$cluster, mean, na.rm = TRUE)
  high_cluster <- as.integer(names(which.max(composite_by_cluster)))
  
  ann$cluster <- ifelse(ann$cluster == high_cluster, "immune-inflamed", "immune-depleted")
  ann$cluster <- factor(ann$cluster, levels = c("immune-inflamed", "immune-depleted"))
  
  # Order samples by cluster for heatmap
  ann <- ann[order(ann$cluster), , drop = FALSE]
  ordered_samples <- rownames(ann)
  
  # 6) Heatmap matrix (signatures x samples)
  data_t <- t(data)[, ordered_samples, drop = FALSE]
  ann_colors <- list(cluster = c("immune-inflamed" = my_colors[3], "immune-depleted" = my_colors[6]))
  
  outTumorDir <- file.path(outputDir, tumor)
  dir.create(outTumorDir, showWarnings = FALSE, recursive = TRUE)
  
  # Gap lines after the first 4 (CD274/CD8/TIL/TLS) and after 7 (add B/Mast/TANs) if present
  gap1 <- min(4, nrow(data_t))
  gap2 <- min(7, nrow(data_t))
  gaps_row <- unique(c(gap1, gap2))
  gaps_row <- gaps_row[gaps_row > 0 & gaps_row < nrow(data_t)]
  
  tiff(file.path(outTumorDir, "heatmap.tiff"), width = 6.5, height = 5, units = "in", res = 300)
  pheatmap(
    data_t,
    annotation_col = ann,
    annotation_colors = ann_colors,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    gaps_row = gaps_row
  )
  dev.off()
  
  # 7) Save
  saveRDS(ann,    file = file.path(outTumorDir, "Group.rds"))
  saveRDS(data_t, file = file.path(outTumorDir, "Matrix.rds"))
}
