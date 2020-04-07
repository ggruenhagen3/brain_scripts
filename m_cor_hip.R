# Mouse Atlas
library(loomR)
library("scrattch.io")
library("Seurat")
library("Matrix")
library("stringr")
library("dplyr")
library("cowplot")
library("ggplot2")
options(stringsAsFactors = FALSE)

rna_path  <- "/nv/hp10/ggruenhagen3/scratch/brain/brain_scripts/"
tome_path <- "/nv/hp10/ggruenhagen3/scratch/brain/data/m_cor_hip/"

tome <- paste(tome_path, "transcrip.tome", sep="")
exons       <- read_tome_dgCMatrix(tome,"data/t_exon")
introns     <- read_tome_dgCMatrix(tome,"data/t_intron")
sample_name <- read_tome_sample_names(tome)
gene_name   <- read_tome_gene_names(tome)

# Plot original TSNE
tsne_coordinates <- read.csv(paste(tome_path, "2d_coordinates.csv", sep=""), header = TRUE)
metadata <- read.csv(paste(tome_path, "sample_annotations.csv", sep=""), header = TRUE)
head(tsne_coordinates)
head(metadata)
original_df <- inner_join(tsne_coordinates, metadata, by = "sample_name")
head(original_df)
# original_df <- cbind(tsne_coordinates[2:3], metadata$cluster_color)
# colnames(original_df) <- c("tsne_1", "tsne_2", "cluster_color")
png(paste(tome_path, "original_tsne.png", sep=""), width = 4800, height = 3600, unit = "px")
ggplot(original_df, aes(tsne_1, tsne_2, color=cluster_color)) + geom_point() + theme_classic()
dev.off()

# Redo Clustering in Seurat
colnames(exons) <- sample_name
colnames(introns) <- sample_name
rownames(exons) <- gene_name
rownames(introns) <- gene_name
m_cor_hip <- exons + introns
m_cor_hip <- CreateSeuratObject(counts = m_cor_hip, project = "M_COR_HIP")
m_cor_hip <- NormalizeData(m_cor_hip, normalization.method = "LogNormalize", scale.factor = 100000)
# exons <- CreateSeuratObject(counts = exons, project = "EXON")
# introns <- CreateSeuratObject(counts = introns, project = "INTRON")
# exons$cond <- "EXON"
# introns$cond <- "INTRON"
# exons <- NormalizeData(exons, normalization.method = "LogNormalize", scale.factor = 100000)
# introns <- NormalizeData(introns, normalization.method = "LogNormalize", scale.factor = 100000)
# m_cor_hip <- merge(exons, introns, merge.data = TRUE)
cat(paste("Number of rows in metadata:", nrow(metadata), "\n"))
cat(paste("Number of rows in m_cor_hip (no subset):", nrow(metadata), "\n"))
m_cor_hip$orig_clust <- metadata$cluster_order
m_cor_hip <- subset(m_cor_hip, subset = nFeature_RNA > 500)
m_cor_hip <- FindVariableFeatures(object = m_cor_hip, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
m_cor_hip <- ScaleData(object = m_cor_hip, vars.to.regress = NULL)
m_cor_hip <- RunPCA(m_cor_hip, npcs = 50, verbose = FALSE)
m_cor_hip <- RunTSNE(object = m_cor_hip)
m_cor_hip <- RunUMAP(m_cor_hip, reduction = "pca", dims = 1:20)
m_cor_hip <- FindNeighbors(m_cor_hip, reduction = "umap", dims = 1:2)
m_cor_hip <- FindClusters(m_cor_hip, resolution = 0.1)

png(paste(tome_path, "new_umap.png", sep=""), width = 4800, height = 4800, unit = "px")
Idents(m_cor_hip) <- "seurat_clusters"
DimPlot(m_cor_hip, reduction = "umap", label = TRUE)
dev.off()

png(paste(tome_path, "old_umap.png", sep=""), width = 4800, height = 4800, unit = "px")
Idents(m_cor_hip) <- "cluster_order"
DimPlot(m_cor_hip, reduction = "umap", label = TRUE)
dev.off()

png(paste(tome_path, "new_tsne.png", sep=""), width = 4800, height = 4800, unit = "px")
Idents(m_cor_hip) <- "seurat_clusters"
DimPlot(m_cor_hip, reduction = "tsne", label = TRUE)
dev.off()

png(paste(tome_path, "old_tsne.png", sep=""), width = 4800, height = 4800, unit = "px")
Idents(m_cor_hip) <- "cluster_order"
DimPlot(m_cor_hip, reduction = "tsne", label = TRUE)
dev.off()
