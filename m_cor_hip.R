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
colnames(original_df) <- c("tsne_1", "tsne_2", "cluster_color")
png(paste(tome_path, "original_tsne.png", sep=""), width = 7200, height = 7200, unit = "px")
ggplot(original_df, aes(tsne_1, tsne_2, color=cluster_color)) + geom_point()
dev.off()


# Redo Clustering in Seurat
# exons <- CreateSeuratObject(counts = exons, project = "EXON")
# introns <- CreateSeuratObject(counts = introns, project = "INTRON")
# exons$cond <- "EXON"
# introns$cond <- "INTRON"
