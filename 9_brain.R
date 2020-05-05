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
options(future.globals.maxSize = 6000 * 1024^2)

source("/nv/hp10/ggruenhagen3/scratch/brain/brain_scripts/all_f.R")

# path <- "C:/Users/miles/Downloads/brain/data/9_brain/"
global_path <- "/nv/hp10/ggruenhagen3/scratch/brain/data/9_brain/region/"
regions <- dir(global_path, pattern =paste("*", sep=""))
regions <- c("globulus_pallidus", "hippocampus", "striatum")

# metadata <- readRDS(paste0(path, "annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS"))
# regions <- regions[2:length(regions)]
for (region in regions) {
  obj_str <- region
  print(paste0("Region:", obj_str))
  path <- paste0(global_path, region, "/")
  dge.path <- list.files(path, pattern = paste("*.gz", sep=""), full.names = TRUE)
  dge <- loadSparseDge(dge.path)
  cluster_assign_orig <- readRDS( list.files(path, pattern = paste("*.cluster.assign.RDS", sep=""), full.names = TRUE)[1] )
  subcluster_assign   <- readRDS( list.files(path, pattern = paste("*.subcluster.assign.RDS", sep=""), full.names = TRUE) )
  
  obj <- CreateSeuratObject(counts = dge, project = obj_str)
  # obj <- FindVariableFeatures(object = obj, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
  # obj <- ScaleData(object = obj, vars.to.regress = NULL)
  # obj <- RunPCA(obj, npcs = 30, verbose = FALSE) 
  # obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
  # obj <- FindNeighbors(obj, reduction = "umap", dims = 1:2)
  # obj <- FindClusters(obj, resolution = 0.25)
  # 
  # Combine the datasets
  tj <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/tj.rds")
  obj_hgnc <- convertToHgncObj(obj, "mouse")
  b1b2mz_hgnc <- readRDS("/nv/hp10/ggruenhagen3/scratch/brain/data/B1C1C2MZ_combined_031020_hgnc.rds")
  result <- keepCommonGenesObj(obj_hgnc, b1b2mz_hgnc)
  obj_hgnc_common <- result[[1]]
  b1b2mz_hgnc_common <- result[[2]]
  obj_hgnc_common[[paste0("NormalizeData.RNA")]] <- tj[[paste0("NormalizeData.RNA")]]
  b1b2mz_hgnc_common[[paste0("NormalizeData.RNA")]] <- tj[[paste0("NormalizeData.RNA")]]
  
  # Remove old objects
  rm(dge)
  rm(obj)
  rm(obj_hgnc)
  rm(b1b2mz_hgnc)
  
  # Integrate them
  obj_hgnc_common <- SCTransform(obj_hgnc_common)
  b1b2mz_hgnc_common <- SCTransform(b1b2mz_hgnc_common)
  integrate.features <- SelectIntegrationFeatures(object.list = list(obj_hgnc_common, b1b2mz_hgnc_common), nfeatures = 3000)
  integrate.list     <- PrepSCTIntegration(object.list = list(obj_hgnc_common, b1b2mz_hgnc_common), anchor.features = integrate.features, verbose = FALSE)
  integrated.anchors <- FindIntegrationAnchors(object.list = integrate.list, normalization.method = "SCT", anchor.features = integrate.features, verbose = FALSE)
  integrated <-IntegrateData(anchorset = integrated.anchors, normalization.method = "SCT", verbose = FALSE)
  integrated <- FindVariableFeatures(integrated)
  integrated <- RunPCA(integrated, npcs = 20, verbose = FALSE)
  integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:20) 
  integrated <- FindNeighbors(integrated, reduction = "umap", dims = 1:2)
  integrated <- FindClusters(integrated, resolution = 0.1)
  # # merge(obj_hgnc_common, b1b2mz_hgnc_common)
  # anchors <- FindIntegrationAnchors(object.list = list(b1b2mz_hgnc_common, obj_hgnc_common), reference = 2, dims = 1:30, k.filter = 150)
  # integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
  # DefaultAssay(integrated) <- "integrated"
  # # DefaultAssay(integrated) <- "RNA"
  # integrated <- FindVariableFeatures(integrated)
  # integrated <- ScaleData(object = integrated, vars.to.regress = NULL)
  # integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
  # integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:20)
  # integrated <- FindNeighbors(integrated, reduction = "umap", dims = 1:2)
  # integrated <- FindClusters(integrated, resolution = 0.1)
  
  integrated$cond[which(is.na(integrated$cond))] <- obj_str
  integrated$cond <- factor(integrated$cond, levels = c('BHVE', 'CTRL', obj_str))
  integrated$cond <- plyr::revalue(integrated$cond, c("BHVE" = "BHVE", "CTRL" = "CTRL", obj_str = obj_str, "NA" = obj_str))
  # Add cluster metadata
  cluster_assign <- cluster_assign_orig[match(colnames(obj_hgnc_common), names(cluster_assign_orig))]
  if (length(colnames(integrated)) != length(cluster_assign)) {
    cluster_assign <- c(as.vector(cluster_assign), rep("NA", length(colnames(obj_hgnc_common)) - length(cluster_assign)))
  }
  orig.cluster <- rep(-1, length(colnames(integrated)))
  orig.cluster[match(colnames(b1b2mz_hgnc_common), colnames(integrated))[1]:tail(match(colnames(b1b2mz_hgnc_common), colnames(integrated)), 1)] <- b1b2mz_hgnc_common$seurat_clusters
  orig.cluster[match(colnames(obj_hgnc_common),    colnames(integrated))[1]:tail(match(colnames(obj_hgnc_common),    colnames(integrated)), 1)] <- cluster_assign
  integrated$orig.cluster <- orig.cluster
  
  # Plots
  filename <- paste0(path, "combined.png")
  png(filename, width = 1000, height = 600)
  Idents(integrated) <- "seurat_clusters"
  p <- DimPlot(integrated, reduction = "umap", label = TRUE)
  print(p)
  dev.off()
  system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/Brain/9_brain/region/", region, sep=""))
  
  filename <- paste0(path, "split.png")
  png(filename, width = 1000, height = 600)
  Idents(integrated) <- "orig.cluster"
  p1 <- DimPlot(integrated[,colnames(b1b2mz_hgnc_common)], reduction = "umap", label = TRUE) + xlim(c(min(integrated@reductions$umap@cell.embeddings[,1]), max(integrated@reductions$umap@cell.embeddings[,1]))) + ylim(c(min(integrated@reductions$umap@cell.embeddings[,2]), max(integrated@reductions$umap@cell.embeddings[,1])))
  p2 <- DimPlot(integrated[,colnames(obj_hgnc_common)], reduction = "umap", label = TRUE)    + xlim(c(min(integrated@reductions$umap@cell.embeddings[,1]), max(integrated@reductions$umap@cell.embeddings[,1]))) + ylim(c(min(integrated@reductions$umap@cell.embeddings[,2]), max(integrated@reductions$umap@cell.embeddings[,1])))
  p <- plot_grid(p1, p2)
  print(p)
  dev.off()
  system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/Brain/9_brain/region/", region, sep=""))
  
  for (i in 0:40) {
    print(i)
    filename <- paste0(path, "split/cluster_", i, ".png")
    png(filename, width = 1000, height = 600)
    cells_to_plot <- c(cluster_cells, colnames(obj_hgnc_common))
    Idents(integrated) <- "cond"
    Idents(b1b2mz_hgnc_common) <- "seurat_clusters"
    cluster_cells <- WhichCells(b1b2mz_hgnc_common, idents = i)
    p <- DimPlot(integrated[,cells_to_plot], reduction = "umap", label = TRUE)
    # p1 <- DimPlot(integrated[,cluster_cells], reduction = "umap", label = TRUE) + xlim(c(min(integrated@reductions$umap@cell.embeddings[,1]), max(integrated@reductions$umap@cell.embeddings[,1]))) + ylim(c(min(integrated@reductions$umap@cell.embeddings[,2]), max(integrated@reductions$umap@cell.embeddings[,1])))
    # p2 <- DimPlot(integrated[,colnames(obj_hgnc_common)], reduction = "umap", label = TRUE) + xlim(c(min(integrated@reductions$umap@cell.embeddings[,1]), max(integrated@reductions$umap@cell.embeddings[,1]))) + ylim(c(min(integrated@reductions$umap@cell.embeddings[,2]), max(integrated@reductions$umap@cell.embeddings[,1])))
    # p <- plot_grid(p1, p2)
    print(p)
    dev.off()
    system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/Brain/9_brain/region/", region, "/split", sep=""))
  }
}

filename <- paste0(path, "vln.png")
png(filename, width = 1000, height = 1000)
Idents(integrated) <- "cond"
p <- VlnPlot(integrated, features = "PVALB", slot = "data")
print(p)
dev.off()
system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/Brain/9_brain/region/", region, sep=""))

# The Mallet
# b1b2mz_hgnc_common <- FindVariableFeatures(b1b2mz_hgnc_common)
# obj_hgnc_common <- FindVariableFeatures(obj_hgnc_common)
# test <- RunCCA(b1b2mz_hgnc_common, obj_hgnc_common, renormalize = TRUE, rescale = TRUE)
# combined <- FindVariableFeatures(test)
# combined <- RunPCA(combined, npcs = 20, verbose = FALSE)
# combined <- RunUMAP(combined, reduction = "pca", dims = 1:20) 
# combined <- FindNeighbors(combined, reduction = "umap", dims = 1:2)
# combined <- FindClusters(combined, resolution = 0.1)

# filename <- paste0(path, "test_3.png")
# png(filename, width = 1000, height = 1000)
# p <- DimPlot(combined, reduction = "umap", label = TRUE)
# print(p)
# dev.off()
# system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/Brain/9_brain/", sep=""))
# 
# filename <- paste0(path, "test_4.png")
# png(filename, width = 1000, height = 1000)
# p <- DimPlot(combined, reduction = "umap", split.by = "cond", label = TRUE)
# print(p)
# dev.off()
# system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/Brain/9_brain/", sep=""))
