library("stringr")
library("ggplot2")
library("biomaRt")
library("DropSeq.util")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("cowplot")
library("RColorBrewer")
library("dplyr")

source("/nv/hp10/ggruenhagen3/scratch/brain/brain_scripts/all_f.R")

# path <- "C:/Users/miles/Downloads/brain/data/9_brain/"
path <- "/nv/hp10/ggruenhagen3/scratch/brain/data/9_brain/region/striatum/"

# metadata <- readRDS(paste0(path, "annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS"))
dge.path <- paste0(path, "F_GRCm38.81.P60Striatum.raw.dge.txt.gz")
dge <- loadSparseDge(dge.path)
cluster_assign <- readRDS(paste0(path, "F_GRCm38.81.P60Striatum.cluster.assign.RDS"))
subcluster_assign <- readRDS(paste0(path, "F_GRCm38.81.P60Striatum.subcluster.assign.RDS"))

obj_str <- "striatum"
obj <- CreateSeuratObject(counts = dge, project = obj_str)
obj <- FindVariableFeatures(object = obj, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
obj <- ScaleData(object = obj, vars.to.regress = NULL)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE) 
# obj <- RunTSNE(object = obj)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
obj <- FindNeighbors(obj, reduction = "umap", dims = 1:2)
obj <- FindClusters(obj, resolution = 0.25)

# Add Metadata
cluster_assign <- cluster_assign[match(colnames(obj), names(cluster_assign))]
if (length(colnames(obj)) != length(cluster_assign)) {
  hi <- c(as.vector(cluster_assign), rep("NA", length(colnames(obj)) - length(cluster_assign)))
}
obj$cluster_assign <- cluster_assign

tj <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/tj.rds")
obj_hgnc <- convertToHgncObj(obj, "mouse")
b1b2mz_hgnc <- readRDS("/nv/hp10/ggruenhagen3/scratch/brain/data/B1C1C2MZ_combined_031020_hgnc.rds")
result <- keepCommonGenesObj(obj_hgnc, b1b2mz_hgnc)
obj_hgnc_common <- result[[1]]
b1b2mz_hgnc_common <- result[[2]]
obj_hgnc_common[[paste0("NormalizeData.RNA")]] <- tj[[paste0("NormalizeData.RNA")]]
b1b2mz_hgnc_common[[paste0("NormalizeData.RNA")]] <- tj[[paste0("NormalizeData.RNA")]]
b1b2mz_hgnc_common <- FindVariableFeatures(b1b2mz_hgnc_common)
obj_hgnc_common <- FindVariableFeatures(obj_hgnc_common)
test <- RunCCA(b1b2mz_hgnc_common, obj_hgnc_common, renormalize = TRUE, rescale = TRUE)
combined <- FindVariableFeatures(test)
combined <- RunPCA(combined, npcs = 20, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20) 
combined <- FindNeighbors(combined, reduction = "umap", dims = 1:2)
combined <- FindClusters(combined, resolution = 0.1)

filename <- paste0(path, "test_3.png")
png(filename, width = 1000, height = 1000)
p <- DimPlot(obj, reduction = "umap", label = TRUE)
print(p)
dev.off()
system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/Brain/9_brain/", sep=""))

filename <- paste0(path, "test_4.png")
png(filename, width = 1000, height = 1000)
p <- DimPlot(obj, reduction = "umap", split.by = "cond", label = TRUE)
print(p)
dev.off()
system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/Brain/9_brain/", sep=""))


# 
# filename <- paste0(path, "test.png")
# png(filename, width = 1000, height = 1000)
# p <- DimPlot(obj, reduction = "umap", label = TRUE)
# print(p)
# dev.off()
# system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/Brain/9_brain/", sep=""))
# 
# 
# filename <- paste0(path, "test_2.png")
# png(filename, width = 1000, height = 1000)
# Idents(obj) <- "cluster_assign"
# p <- DimPlot(obj, reduction = "umap", label = TRUE)
# print(p)
# dev.off()
# system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/Brain/9_brain/", sep=""))