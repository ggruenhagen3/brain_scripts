library(Seurat)
library(dplyr)

suppressMessages(require(Seurat))
suppressMessages(require(ggplot2))
suppressMessages(require(cowplot))

install.packages('scater')
library(scater)

#install scran
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scran")
suppressMessages(require(scran))
library(scran)

BiocManager::install("BiocParallel")
suppressMessages(require(BiocParallel))

BiocManager::install("BiocNeighbors")
suppressMessages(require(BiocNeighbors))

rna_path <-"C:/Users/miles/Downloads/brain/" 

b1.data <- Read10X(data.dir = paste(rna_path, "data/B1-from-bcl-lncRNA/outs/filtered_feature_bc_matrix/", sep=""))
b2.data <- Read10X(data.dir = paste(rna_path, "data/B2-from-bcl-lncRNA/outs/filtered_feature_bc_matrix/", sep=""))
c1.data <- Read10X(data.dir = paste(rna_path, "data/C1-from-bcl-lncRNA/outs/filtered_feature_bc_matrix/", sep=""))
c2.data <- Read10X(data.dir = paste(rna_path, "data/MZ-from-bcl-lncRNA/outs/filtered_feature_bc_matrix/", sep=""))

b1 <- CreateSeuratObject(counts = b1.data, project = "behav")
b2 <- CreateSeuratObject(counts = b2.data, project = "behav2")
c1 <- CreateSeuratObject(counts = c1.data, project = "control")
c2 <- CreateSeuratObject(counts = c2.data, project = "control2")

b1$sample <- "b1"
b2$sample <- "b2"
c1$sample <- "c1"
c2$sample <- "c2"

b1$cond <- "BHVE"
c1$cond <- "CTRL"
b2$cond <- "BHVE"
c2$cond <- "CTRL"

b1 <- NormalizeData(b1, normalization.method = "LogNormalize")
b2 <- NormalizeData(b2, normalization.method = "LogNormalize")
c1 <- NormalizeData(c1, normalization.method = "LogNormalize")
c2 <- NormalizeData(c2, normalization.method = "LogNormalize")

VlnPlot(b1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
hist(b1$nFeature_RNA, labels = TRUE, breaks = 50)
hist(b1$nCount_RNA, labels = TRUE, breaks = 50)
b1 <- subset(b1, subset = nFeature_RNA > 500)
b1 <- subset(b1, subset = nFeature_RNA < 2500)
b1 <- subset(b1, subset = nCount_RNA > 500)
b1 <- subset(b1, subset = nCount_RNA < 6000)
#b1 <- FindVariableFeatures(object = b1, selection.method = "mvp", num.bin = 50, binning.method = "equal_width", verbose = TRUE)
b1 <- FindVariableFeatures(object = b1, selection.method = "vst", verbose = TRUE)
VariableFeaturePlot(b1)

b1_var = b1@assays$RNA@var.features

# Set up behave2 object
VlnPlot(b2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
hist(b2$nFeature_RNA, labels = TRUE, breaks = 50)
hist(b2$nCount_RNA, labels = TRUE, breaks = 50)
b2 <- subset(b2, subset = nFeature_RNA > 500)
b2 <- subset(b2, subset = nFeature_RNA < 2500)
b2 <- subset(b2, subset = nCount_RNA > 500)
b2 <- subset(b2, subset = nCount_RNA < 6000)
#b2 <- FindVariableFeatures(object = b2, selection.method = "mvp", num.bin = 50, binning.method = "equal_width", verbose = TRUE)
b2 <- FindVariableFeatures(object = b2, selection.method = "vst", verbose = TRUE)
VariableFeaturePlot(b2)

b2_var = b2@assays$RNA@var.features

# Set up control object
VlnPlot(c1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
hist(c1$nFeature_RNA, labels = TRUE, breaks = 50)
hist(c1$nCount_RNA, labels = TRUE, breaks = 50)
c1 <- subset(c1, subset = nFeature_RNA > 500)
c1 <- subset(c1, subset = nFeature_RNA < 2500)
c1 <- subset(c1, subset = nCount_RNA > 500)
c1 <- subset(c1, subset = nCount_RNA < 6000)
#c1 <- FindVariableFeatures(object = c1, selection.method = "mvp", num.bin = 50, binning.method = "equal_width", verbose = TRUE)
c1 <- FindVariableFeatures(object = c1, selection.method = "vst", verbose = TRUE)
VariableFeaturePlot(c1)

c1_var = c1@assays$RNA@var.features

# Set up control object
VlnPlot(c2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
hist(c2$nFeature_RNA, labels = TRUE, breaks = 50)
hist(c2$nCount_RNA, labels = TRUE, breaks = 50)
c2 <- subset(c2, subset = nFeature_RNA > 300)
c2 <- subset(c2, subset = nFeature_RNA < 2000)
# c2 <- subset(c2, subset = nCount_RNA > 500)
# c2 <- subset(c2, subset = nCount_RNA < 6000)
#c2 <- FindVariableFeatures(object = c2, selection.method = "mvp", num.bin = 50, binning.method = "equal_width", verbose = TRUE)
c2 <- FindVariableFeatures(object = c2, selection.method = "vst", verbose = TRUE)
VariableFeaturePlot(c2)

c2_var = c2@assays$RNA@var.features

combined <- merge(x=b1, y=c(b2,c1,c2), merge.data = TRUE, add.cell.ids = c("BHVE", "BHVE", "CTRL", "CTRL"))


#combined <- readRDS("C:/Users/zjohnson37/Desktop/from-bcl/c41_scale1mil_b1b2c1_012219 (1).rds")
#combined <- readRDS("C:/Users/zjohnson37/Desktop/from-bcl/c41_scale1mil_b1b2c1_121619.rds")
DefaultAssay(combined) <- "RNA"


#looking at gene expression plots following different workflows
VlnPlot(c1, features="egr1")
VlnPlot(b1, features="egr1")
VlnPlot(b2, features="egr1")
VlnPlot(combined, features="egr1")

var_list <- unique(c(b1_var, b2_var, c1_var, c2_var))
length(unique(var_list))

combined <- FindVariableFeatures(object = combined, selection.method = "vst", nfeatures = 10000, verbose = TRUE)

DefaultAssay(combined) <- "RNA"

combined <- ScaleData(combined, verbose = TRUE)
combined <- RunPCA(combined, dim = 50, verbose = TRUE)

pca = combined@reductions$pca
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / sum(eigValues)
write.csv(varExplained,"PCA_var_explained_031020.csv", row.names = TRUE)

combined <- FindNeighbors(combined, reduction = "pca", dims = 1:50)
combined <- FindClusters(combined, resolution = 1.86)
combined <- RunUMAP(combined, dims = 1:50)

saveRDS(combined, file = "C:/Users/zjohnson37/Desktop/from-bcl/B1C1C2MZ_combined_031320.rds")
combined <- readRDS(file = "C:/Users/zjohnson37/Desktop/from-bcl/B1C1C2MZ_combined_031020.rds")
combined <- readRDS(file = "C:/Users/zjohnson37/Desktop/from-bcl/B1C1C2MZ_combined_031320.rds")

p0 <- DimPlot(combined, reduction = "umap", label = TRUE)

old_cluster <- colnames(combined)
for (i in 0:40) {
  old_cluster[which(substr(old_cluster, 6, 21) %in% substr(WhichCells(b1b2c1mz, idents = i), 6, 21))] <- i
  old_cluster[which(old_cluster %in% WhichCells(b1b2c1mz, idents = i))] <- i
}
old_cluster <- as.numeric(old_cluster)
old_cluster[is.na(old_cluster)] <- "Removed"
old_cluster <- factor(old_cluster, levels = c(0:40, "Removed"))
combined$old_cluster <- old_cluster
Idents(combined) <- combined$old_cluster
DimPlot(combined, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 1.2) + ggtitle("Zack - Master Clusters")
removed_cells <- WhichCells(combined, idents = "Removed") 
Idents(combined) <- "seurat_clusters"
DimPlot(combined[,removed_cells], reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 1.2) + ggtitle("Removed Cells")
DimPlot(combined, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 1.2)
