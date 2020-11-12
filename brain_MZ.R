library("edgeR")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("cowplot")
library("monocle3")

rna_path <- "C:/Users/miles/Downloads/brain/"

mz.data <- Read10X(data.dir = paste(rna_path, "data/MZ/outs/filtered_feature_bc_matrix/", sep=""))

mz <- CreateSeuratObject(counts = mz.data, project = "MZ")
mz$cond <- "MZ"
mz$sample <- rep("MZ", ncol(mz))
mz <- NormalizeData(mz, normalization.method = "LogNormalize", scale.factor = 100000)

FeatureScatter(mz, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
mz <- subset(mz, subset = nFeature_RNA > 300)
mz <- subset(mz, subset = nFeature_RNA < 2000)
mz <- ScaleData(object = mz, vars.to.regress = NULL)
mz <- FindVariableFeatures(object = mz, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
mz <- RunPCA(mz, npcs = 30, verbose = FALSE, features = mz@assays$RNA@var.features)
# JackStraw Values: 1.737076e-140 1.758434e-105 1.284583e-152  5.534194e-68  2.436125e-88 1.378233e-109 9.145379e-96  2.720489e-72  1.141054e-84  2.525496e-53  3.244567e-62  1.508244e-35 1.418974e-39  3.448115e-15  1.140568e-20  4.684089e-28  5.821600e-15  1.835193e-12 4.046250e-24  5.736455e-10
# Steep dropoff after the 9th PC and ElbowPlot shows a dropoff after the 9th PC

# mz <- JackStraw(mz, num.replicate = 100)
# mz <- ScoreJackStraw(mz, dims = 1:20)
# JackStrawPlot(mz, dims = 1:20)

mz <- RunUMAP(mz, reduction = "pca", dims = 1:9) 
mz <- FindNeighbors(mz, reduction = "umap", dims = 1:2)
mz <- FindClusters(mz, resolution = 0.35)
DimPlot(mz, reduction = "umap", split.by = "cond", label = TRUE)
# saveRDS(mz, "C:/Users/miles/Downloads/brain/brain_scripts/mz_shiny/data/mz.rds")

# DEG
rna_path <- "C:/Users/miles/Downloads/brain/"
nk.markers <- FindAllMarkers(mz)
nk.markers <- nk.markers[,c(ncol(nk.markers)-1, ncol(nk.markers), 1:(ncol(nk.markers)-2)),]
nk.markers <- nk.markers[which(nk.markers$p_val_adj < 0.05),]
write.table(nk.markers, file = paste(rna_path, "/results/mz_deg.tsv", sep=""), sep="\t", quote = FALSE, row.names = FALSE)

################################
# Combined with b1, b2, and c1 #
################################
# Load old data
combined_mz <- readRDS("C:/Users/miles/Downloads/brain/brain_scripts/brain_mz_shiny/data/B1C1C2MZ_combined_031020.rds")
combined <- readRDS("C:/Users/miles/Downloads/brain/brain_scripts/brain_shiny/data/combined.rds")
mz <- readRDS("C:/Users/miles/Downloads/brain/brain_scripts/mz_shiny/data/mz.rds")

# Compare old data to new data
Idents(combined_mz) <- "sample"
combined_mz <- SetIdent(combined_mz, cells=WhichCells(combined_mz, idents = c("b1", "b2", "c1")), value="sand")
combined_mz <- SetIdent(combined_mz, cells=WhichCells(combined_mz, idents = c("c2")), value="rock")
combined_mz$rock_sand <- combined_mz@active.ident

num_combined_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
num_mz_clusters <- as.numeric(tail(levels(mz@meta.data$seurat_clusters), n=1))
for (combined_clust in 0:num_combined_clusters) {
  combined_mz <- SetIdent(combined_mz, cells=WhichCells(combined, idents = combined_clust), value=paste("sand", combined_clust, sep="_"))
}
for (mz_clust in 0:num_mz_clusters) {
  combined_mz <- SetIdent(combined_mz, cells=paste("CTRL", WhichCells(mz, idents = mz_clust), "4", sep="_"), value=paste("rock", mz_clust, sep="_"))
}

# Plot all old and new data together
combined_mz$orig.cluster <- combined_mz@active.ident
Idents(combined_mz) <- "orig.cluster"
png(filename = paste(rna_path, "results/painting/brain_MZ/all_old.png", sep=""), width = 4800, height = 2400, unit="px", res=300)
DimPlot(combined_mz, reduction = "umap", split.by = "rock_sand", pt.size = 1, label = TRUE)
dev.off()

# Plot just the rock data
png(filename = paste(rna_path, "results/painting/brain_MZ/mz_old.png", sep=""), width = 4800, height = 2400, unit="px", res=300)
rock_clusters <- unique(combined_mz$orig.cluster)[startsWith(as.character(unique(combined_mz$orig.cluster)), "rock") & unique(combined_mz$orig.cluster) != "rock"]
DimPlot(combined_mz[,WhichCells(combined_mz, idents=rock_clusters)], reduction = "umap", split.by = "rock_sand", label = TRUE)
dev.off()

# Plot the rock data with the new data
png(filename = paste(rna_path, "results/painting/brain_MZ/mz_old_comparison.png", sep=""), width = 7200, height = 2400, unit="px", res=300)
Idents(combined_mz) <- "orig.cluster"
rock_clusters <- unique(combined_mz$orig.cluster)[startsWith(as.character(unique(combined_mz$orig.cluster)), "rock") & unique(combined_mz$orig.cluster) != "rock"]
rock_cells <- WhichCells(combined_mz, idents=rock_clusters)
p1 <- DimPlot(combined_mz[,rock_cells], reduction = "umap", split.by = "rock_sand", label = TRUE)
Idents(combined_mz) <- "seurat_clusters"
p2 <- DimPlot(combined_mz, reduction = "umap", label = TRUE)
p3 <- DimPlot(combined_mz[,rock_cells], reduction = "umap", label = TRUE)
plot_grid(p1, p3, p2, ncol = 3, nrow = 1)
dev.off()

# Plot just the sand data
Idents(combined_mz) <- "orig.cluster"
png(filename = paste(rna_path, "results/painting/brain_MZ/sand_old.png", sep=""), width = 4800, height = 2400, unit="px", res=300)
sand_clusters <- unique(combined_mz$orig.cluster)[startsWith(as.character(unique(combined_mz$orig.cluster)), "sand") & unique(combined_mz$orig.cluster) != "sand"]
DimPlot(combined_mz[,WhichCells(combined_mz, idents=sand_clusters)], reduction = "umap", split.by = "rock_sand", label = TRUE)
dev.off()

# Plot the sand data next to the new data
png(filename = paste(rna_path, "results/painting/brain_MZ/sand_old_comparison.png", sep=""), width = 7200, height = 2400, unit="px", res=300)
Idents(combined_mz) <- "orig.cluster"
sand_clusters <- unique(combined_mz$orig.cluster)[startsWith(as.character(unique(combined_mz$orig.cluster)), "sand") & unique(combined_mz$orig.cluster) != "sand"]
sand_cells <- WhichCells(combined_mz, idents=sand_clusters)
p1 <- DimPlot(combined_mz[,sand_cells], reduction = "umap", split.by = "rock_sand", label = TRUE)
Idents(combined_mz) <- "seurat_clusters"
p2 <- DimPlot(combined_mz, reduction = "umap", label = TRUE)
p3 <- DimPlot(combined_mz[,sand_cells], reduction = "umap", label = TRUE)
plot_grid(p1, p3, p2, ncol = 3, nrow = 1)
dev.off()

# Count number of cells in each cluster
Idents(combined_mz) <- "sample"
b1 <- WhichCells(combined_mz, idents = "b1")
b2 <- WhichCells(combined_mz, idents = "b2")
c1 <- WhichCells(combined_mz, idents = "c1")
c2 <- WhichCells(combined_mz, idents = "c2")
Idents(combined_mz) <- "seurat_clusters"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
df <- c()
for (i in 0:num_clusters) {
  all_cells <- length(WhichCells(combined_mz, idents = i))
  b1_cells <- 0
  b2_cells <- 0
  c1_cells <- 0
  c2_cells <- 0
  
  try(b1_cells <- length(WhichCells(combined_mz[,b1], idents = i)), silent=TRUE)
  try(b2_cells <- length(WhichCells(combined_mz[,b2], idents = i)), silent=TRUE)
  try(c1_cells <- length(WhichCells(combined_mz[,c1], idents = i)), silent=TRUE)
  try(c2_cells <- length(WhichCells(combined_mz[,c2], idents = i)), silent=TRUE)
  
  df <- rbind(df, t(c(i, all_cells, b1_cells, b2_cells, c1_cells, c2_cells)))
}
colnames(df) <- c("cluster", "total_num_cells", "b1_num_cells", "b2_num_cells", "c1_num_cells", "c2_num_cells")
write.table(df, paste(rna_path, "/results/num_cells_in_brain_mz.tsv", sep=""), sep="\t", row.names = FALSE, quote=FALSE)