# =========================
# Pan-cancer TCGA: 3D t-SNE (final/reviewer-friendly)
# Purpose:
#   1) Collect per-tumor sample Group labels (immune-inflamed / immune-depleted)
#   2) Collect per-tumor signature score matrices (Matrix.rds)
#   3) Merge into a pan-cancer matrix (samples × signatures)
#   4) Run 3D t-SNE on signature space
#   5) Plot 3D embedding colored by immune group
# Inputs (per tumor folder under group_base_dir):
#   - Group.rds  : data.frame with rownames = sample IDs, one column "cluster" (or equivalent)
#   - Matrix.rds : signature score matrix with columns = samples (typical), rows = signatures
# Output:
#   - plot/tsne_Group.tiff
#   - tsne.rdata
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(Polychrome)
  library(Rtsne)
  library(plot3D)
})

set.seed(666)
dir.create("plot", showWarnings = FALSE)

# ---- Settings ----
group_base_dir <- "pancancer_TCGA"
tumor_types <- list.dirs(group_base_dir, full.names = FALSE, recursive = FALSE)
tumor_types <- tumor_types[!grepl("outs", tumor_types, ignore.case = TRUE)]

# Use a consistent two-color palette for groups
pal2 <- as.character(Polychrome::palette36.colors())[c(3, 6)]  # inflamed / excluded

# ---- 1) Read and merge Group labels across tumor types ----
all_group <- lapply(tumor_types, function(tt) {
  f <- file.path(group_base_dir, tt, "Group.rds")
  g <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(g)) return(NULL)
  
  # Expect rownames are sample IDs; expect one column (e.g. "cluster")
  g <- data.frame(Sample = rownames(g), Group = g[[1]], stringsAsFactors = FALSE)
  g$TumorType <- tt
  g
})
merged_group <- bind_rows(all_group) %>% distinct(Sample, .keep_all = TRUE)

# ---- 2) Read and merge signature matrices across tumor types ----
all_mat <- lapply(tumor_types, function(tt) {
  f <- file.path(group_base_dir, tt, "Matrix.rds")
  m <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  
  # Your upstream Matrix.rds appears to be signatures × samples (as used before)
  # Combine by columns (samples) across tumor types
  m
})
merged_sig <- do.call(cbind, all_mat)  # signatures × samples
merged_sig <- t(merged_sig)            # samples × signatures
merged_sig <- as.data.frame(merged_sig)
merged_sig$Sample <- rownames(merged_sig)

# ---- 3) Merge signatures + group label ----
merged_data <- merged_sig %>%
  inner_join(merged_group[, c("Sample", "Group")], by = "Sample")

# Ensure consistent group naming (match your previous logic)
merged_data$Group <- ifelse(merged_data$Group == "immune-inflamed", "inflamed", "excluded")
merged_data$Group <- factor(merged_data$Group, levels = c("inflamed", "excluded"))

# ---- 4) 3D t-SNE on signature space ----
# Remove the non-numeric columns: Sample + Group
data_matrix <- merged_data %>%
  select(-Sample, -Group) %>%
  as.matrix()

tsne_res <- Rtsne(data_matrix, dims = 3, perplexity = 50, check_duplicates = FALSE)

tsne_df <- data.frame(
  X1 = tsne_res$Y[, 1],
  X2 = tsne_res$Y[, 2],
  X3 = tsne_res$Y[, 3],
  Sample = merged_data$Sample,
  Group  = merged_data$Group
)
save(tsne_df, file = "tsne.rdata")

# ---- 5) Plot: colored by Group ----
tiff("plot/tsne_Group.tiff", width = 1700, height = 1500, res = 300)
scatter3D(
  x = tsne_df$X1, y = tsne_df$X2, z = tsne_df$X3,
  colvar = as.numeric(tsne_df$Group),
  col = pal2,
  pch = 16, cex = 0.35,
  main = "",
  xlab = "t-SNE 1", ylab = "t-SNE 2", zlab = "t-SNE 3",
  theta = 210, phi = 45,
  colkey = FALSE
)
legend(
  "right",
  title = "Group",
  legend = levels(tsne_df$Group),
  pch = 21,
  pt.bg = pal2,
  bty = "n",
  bg = "white"
)
dev.off()