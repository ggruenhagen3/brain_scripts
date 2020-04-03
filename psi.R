# Load Libraries
library("pSI")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("ggplot2")
library("dplyr")

# Read in Data
combined <- readRDS("C:/Users/miles/Downloads/brain/brain_scripts/brain_mz_shiny/data/B1C1C2MZ_combined_031020.rds")
rna_path <- "C:/Users/miles/Downloads/brain/"
marker_path <- paste(rna_path, "data/markers/", sep="")
marker_files <- dir(marker_path, pattern =paste("*.txt", sep=""))
markers <- data.frame(gene <- c(), bio <- c())
for (i in 1:length(marker_files)) {
  file <- read.table(paste(marker_path, marker_files[i], sep=""), header=FALSE, sep="\t", stringsAsFactors=FALSE)
  markers <- rbind(markers, file[,1:2])
}
colnames(markers) <- c("gene", "bio")
bio <- "ROCK_UP"
markers <- markers[which(markers$bio == bio),]
valid_genes <- markers$gene

# Make a dataframe with genes as rows and clusters as columns
cluster_df <- data.frame()
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
cells_per_cluster <- c()
for (i in 0:num_clusters) {
  this_cells <- WhichCells(combined, idents = i)
  cluster_df <- rbind(cluster_df, rowSums(combined@assays$RNA@counts[,this_cells]))
  cells_per_cluster <- c(cells_per_cluster, length(this_cells))
}
cluster_df <- t(cluster_df)
colnames(cluster_df) <- 0:num_clusters
rownames(cluster_df) <- rownames(combined@assays$RNA@counts)
cluster_df_norm <- cluster_df / cells_per_cluster
write.table(cluster_df_norm, paste(rna_path, "/data/psi_input.tsv", sep=""), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Find the Specificity Index
# results <- specificity.index(cluster_df, e_min = 0.01)
results_norm <- specificity.index(cluster_df_norm, e_min = 1/3)
colnames(results) <- 0:num_clusters
colnames(results_norm) <- 0:num_clusters

# fisher_results <- fisher.iteration(results, valid_genes)
fisher_results_norm <- fisher.iteration(results_norm, valid_genes, p.adjust = FALSE)
fisher_results_norm$q <- p.adjust(fisher_results_norm$`0.05 - nominal`)
fisher_results_norm$q.sig <- fisher_results_norm$q < 0.05

# Find the Specificity per Cluster
sum_clust_norm <- colSums(results_norm[valid_genes,] > 0, na.rm = TRUE)
sum_psi_norm <- colSums(results_norm[valid_genes,], na.rm = TRUE)
sum_clust <- colSums(results[valid_genes,] > 0, na.rm = TRUE)
sum_psi <- colSums(results[valid_genes,], na.rm = TRUE)

sum_clust_norm[order(sum_clust_norm, decreasing = TRUE)]
test <- sum_clust/cells_per_cluster
test[order(test, decreasing=TRUE)]
