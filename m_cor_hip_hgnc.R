# Mouse Atlas
library("loomR")
library("scrattch.io")
library("Seurat")
library("Matrix")
library("stringr")
library("dplyr")
library("cowplot")
library("ggplot2")
library("biomaRt")
options(stringsAsFactors = FALSE)

rna_path  <- "/nv/hp10/ggruenhagen3/scratch/brain/brain_scripts/"
tome_path <- "/nv/hp10/ggruenhagen3/scratch/brain/data/m_cor_hip/"

# Load the Seurat Objects with HGNC gene names
tj <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/tj.rds")
m_cor_hip <- readRDS("/nv/hp10/ggruenhagen3/scratch/brain/data/m_cor_hip/m_cor_hip_hgnc_common.rds")
b1b2mz    <- readRDS("/nv/hp10/ggruenhagen3/scratch/brain/data/B1C1C2MZ_combined_031020_hgnc_common.rds")
m_cor_hip$org <- "Mouse"
b1b2mz$org <- "Cichlid"

# Genes that are in both datasets
all_genes <- c(rownames(m_cor_hip), rownames(b1b2mz))
common <- all_genes[duplicated(all_genes)]

combined <- merge(x=m_cor_hip, y=b1b2mz, merge.data = TRUE)
# Standard Pipeline (except excluding HVG's that aren't in common)
combined <- FindVariableFeatures(object = combined, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
combined <- ScaleData(object = combined, vars.to.regress = NULL)
var_common <- combined@assays$RNA@var.features[which(combined@assays$RNA@var.features %in% common)]
combined <- RunPCA(combined, npcs = 50, verbose = FALSE, features = var_common)
combined <- RunTSNE(object = combined)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "umap", dims = 1:2)
combined <- FindClusters(combined, resolution = 0.1)

# Add the old metadata to the new combined object
mouse_cells   <- colnames(combined)[which(colnames(combined) %in% colnames(m_cor_hip))]
cichlid_cells <- colnames(combined)[which(colnames(combined) %in% colnames(b1b2mz))]
m_metadata <- m_cor_hip@meta.data %>% slice(match(colnames(combined), mouse_cells)) # Order for the old metadata to match to order of the new cells
c_metadata <- b1b2mz@meta.data %>% slice(match(colnames(combined), cichlid_cells)) # Order for the old metadata to match to order of the new cells
for (col in colnames(m_metadata)) {
  combined@meta.data[paste("m_", col, sep="")] <- c(m_metadata[[col]], rep("cichlid", ncol(b1b2mz)))
}
for (col in colnames(c_metadata)) {
  combined@meta.data[paste("c_", col, sep="")] <- c(rep("mouse", nrow(m_metadata)), c_metadata[[col]])
}

# integrated    <- readRDS("/nv/hp10/ggruenhagen3/scratch/brain/data/m_cor_hip/integrated.rds")
filename <- paste(tome_path, "/results/w_cichlid/org_split_4.png", sep="")
png(filename, width = 4000, height = 3000, unit = "px", res = 150)
p <- DimPlot(combined, split.by="org", label=FALSE) + guides(col=guide_legend(ncol=4))
print(p)
dev.off()
system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/Brain/m_cor_hip/results/w_cichlid/", sep=""))

# Paint by mouse cell type, one old cichlid cluster at a time
Idents(combined) <- "m_cell_type_alias_label"
Idents(b1b2mz) <- "seurat_clusters"
n = length(levels(Idents(m_cor_hip)))
hues = seq(15, 375, length = n + 1)
cols = hcl(h = hues, l = 65, c = 100)[1:n]
num_old_cichlid_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
for (cluster in 0:num_old_cichlid_clusters) {
  print(paste("Old Cichlid Cluster", cluster))
  c_cluster_cells <- WhichCells(b1b2mz, idents = c(cluster))
  filename <- paste(tome_path, "/results/w_cichlid/old_cichlid_cluster/c_", cluster, ".png", sep="")
  png(filename, width = 4000, height = 3000, unit = "px", res = 150)
  p <- DimPlot(combined[, c(mouse_cells, c_cluster_cells)], reduction = "umap", label = FALSE)  + ggtitle(paste("Old cichlid cluster", cluster, "on Mouse Cortex + Hippocampus Data")) + theme_classic() + guides(col=guide_legend(ncol=4))
  print(p)
  dev.off()
  system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/Brain/m_cor_hip/results/w_cichlid/old_cichlid_cluster/", sep=""))
  # Split By Organism
  filename <- paste(tome_path, "/results/w_cichlid/old_cichlid_cluster/c_", cluster, "_split.png", sep="")
  png(filename, width = 4000, height = 3000, unit = "px", res = 150)
  p <- DimPlot(combined[, c(mouse_cells, c_cluster_cells)], reduction = "umap", split.by = "org", label = FALSE)  + ggtitle(paste("Old cichlid cluster", cluster, "on Mouse Cortex + Hippocampus Data")) + theme_classic() + guides(col=guide_legend(ncol=4))
  print(p)
  dev.off()
    system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/Brain/m_cor_hip/results/w_cichlid/old_cichlid_cluster/", sep=""))
  print("DONE")
}