library("shiny")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("dplyr")
library("ggplot2")
rna_path <- "C:/Users/miles/Downloads/brain/"
combined <- readRDS("C:/Users/miles/Downloads/brain/brain_scripts/brain_shiny/data/combined.rds")

gene_names <- rownames(combined)
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))

# For HVG
Idents(combined) <- combined$sample
b1 <- subset(combined, cells = WhichCells(combined, idents = "b1"))
b2 <- subset(combined, cells = WhichCells(combined, idents = "b2"))
Idents(b1) <- b1$seurat_clusters
Idents(b2) <- b2$seurat_clusters
b1_results <- data.frame(genes <- gene_names)
b2_results <- data.frame(genes <- gene_names)
results <- data.frame(genes <- gene_names)
colnames(results)[1] <- "gene"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
for (i in 0:num_clusters) {
  b1_this_clust <- WhichCells(b1, idents = i)
  b2_this_clust <- WhichCells(b2, idents = i)
  b1_sums <- rowSums(b1@assays$RNA@counts[gene_names,b1_this_clust])
  b2_sums <- rowSums(b2@assays$RNA@counts[gene_names,b2_this_clust])
  b1_avg <- b1_sums/length(b1_this_clust)
  b2_avg <- b2_sums/length(b2_this_clust)
  sds <- sapply(1:length(b1_avg), function(x) sd(c(b1_avg[x], b2_avg[x])))
  results <- cbind(results, sds)
  colnames(results)[i+2] <- i
  
  # b1_results <- cbind(b1_results, b1_sums/length(b1_this_clust))
  # b2_results <- cbind(b2_results, b2_sums/length(b2_this_clust))
  # colnames(b1_results)[i+2] <- i
  # colnames(b2_results)[i+2] <- i
}
write.table(results, paste(rna_path, "/results/hvg_b1_b2_var_expr_lvl_per_cluster.tsv", sep=""), sep = "\t", row.names = FALSE, quote = FALSE)

# For All Genes
b1.data <- Read10X(data.dir = paste(rna_path, "data/BHVE-JTS03-B1-from-bcl/outs/filtered_feature_bc_matrix/", sep=""))
b2.data <- Read10X(data.dir = paste(rna_path, "data/BHVE-JTS02-B2-from-bcl/outs/filtered_feature_bc_matrix/", sep=""))
c1.data <- Read10X(data.dir = paste(rna_path, "data/CTRL-JTS03-C1-from-bcl/outs/filtered_feature_bc_matrix/", sep=""))

b1 <- CreateSeuratObject(counts = b1.data, project = "BHVE")
b2 <- CreateSeuratObject(counts = b2.data, project = "BHVE")
c1 <- CreateSeuratObject(counts = c1.data, project = "CTRL")

all_combined <- merge(x=c1, y=c(b1,b2), merge.data = TRUE, add.cell.ids = c("CTRL", "BHVE", "BHVE"))
gene_names <- rownames(all_combined)

# Average Expr Lvl per Cell B1 v B2
Idents(combined) <- combined$sample
b1 <- subset(combined, cells = WhichCells(combined, idents = "b1"))
b2 <- subset(combined, cells = WhichCells(combined, idents = "b2"))
Idents(b1) <- b1$seurat_clusters
Idents(b2) <- b2$seurat_clusters
b1_results <- data.frame(genes <- gene_names)
b2_results <- data.frame(genes <- gene_names)
results <- data.frame(genes <- gene_names)
colnames(results)[1] <- "gene"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
for (i in 0:num_clusters) {
  b1_this_clust <- WhichCells(b1, idents = i)
  b2_this_clust <- WhichCells(b2, idents = i)
  b1_sums <- rowSums(b1@assays$RNA@counts[gene_names,b1_this_clust])
  b2_sums <- rowSums(b2@assays$RNA@counts[gene_names,b2_this_clust])
  b1_avg <- b1_sums/length(b1_this_clust)
  b2_avg <- b2_sums/length(b2_this_clust)
  sds <- sapply(1:length(b1_avg), function(x) sd(c(b1_avg[x], b2_avg[x])))
  results <- cbind(results, sds)
  colnames(results)[i+2] <- i
}
write.table(results, paste(rna_path, "/results/all_genes_b1_b2_var_expr_lvl_per_cluster.tsv", sep=""), sep = "\t", row.names = FALSE, quote = FALSE)

# Proportion of cells displaying the transcript B1 v B2
Idents(b1) <- b1$seurat_clusters
Idents(b2) <- b2$seurat_clusters
b1_results <- data.frame(genes <- gene_names)
b2_results <- data.frame(genes <- gene_names)
results <- data.frame(genes <- gene_names)
colnames(results)[1] <- "gene"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
for (i in 0:num_clusters) {
  b1_this_clust <- WhichCells(b1, idents = i)
  b2_this_clust <- WhichCells(b2, idents = i)
  b1_sums <- rowSums(b1@assays$RNA@counts[gene_names,b1_this_clust] != 0) # num cells expr this gene
  b2_sums <- rowSums(b2@assays$RNA@counts[gene_names,b2_this_clust] != 0) # num cells expr this gene
  b1_avg <- b1_sums/length(b1_this_clust)
  b2_avg <- b2_sums/length(b2_this_clust)
  sds <- sapply(1:length(b1_avg), function(x) sd(c(b1_avg[x], b2_avg[x])))
  results <- cbind(results, sds)
  colnames(results)[i+2] <- i
}
write.table(results, paste(rna_path, "/results/all_genes_b1_b2_var_prop_expr_per_cluster.tsv", sep=""), sep = "\t", row.names = FALSE, quote = FALSE)


# Average Expr Lvl per Cell B1 v C1
Idents(combined) <- combined$sample
b1 <- subset(combined, cells = WhichCells(combined, idents = "b1"))
c1 <- subset(combined, cells = WhichCells(combined, idents = "c1"))
Idents(b1) <- b1$seurat_clusters
Idents(c1) <- c1$seurat_clusters
results <- data.frame(genes <- gene_names)
colnames(results)[1] <- "gene"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
for (i in 0:num_clusters) {
  b1_this_clust <- WhichCells(b1, idents = i)
  c1_this_clust <- WhichCells(c1, idents = i)
  b1_sums <- rowSums(b1@assays$RNA@counts[gene_names,b1_this_clust])
  c1_sums <- rowSums(as.matrix(c1@assays$RNA@counts[gene_names,c1_this_clust]))
  b1_avg <- b1_sums/length(b1_this_clust)
  c1_avg <- c1_sums/length(c1_this_clust)
  sds <- sapply(1:length(b1_avg), function(x) sd(c(b1_avg[x], c1_avg[x])))
  results <- cbind(results, sds)
  colnames(results)[i+2] <- i
}
write.table(results, paste(rna_path, "/results/all_genes_b1_c1_var_expr_lvl_per_cluster.tsv", sep=""), sep = "\t", row.names = FALSE, quote = FALSE)

# Proportion of cells displaying the transcript B1 v C1
Idents(b1) <- b1$seurat_clusters
Idents(c1) <- c1$seurat_clusters
results <- data.frame(genes <- gene_names)
colnames(results)[1] <- "gene"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
for (i in 0:num_clusters) {
  b1_this_clust <- WhichCells(b1, idents = i)
  c1_this_clust <- WhichCells(c1, idents = i)
  b1_sums <- rowSums(b1@assays$RNA@counts[gene_names,b1_this_clust] != 0) # num cells expr this gene
  c1_sums <- rowSums(as.matrix(c1@assays$RNA@counts[gene_names,c1_this_clust] != 0)) # num cells expr this gene
  b1_avg <- b1_sums/length(b1_this_clust)
  c1_avg <- c1_sums/length(c1_this_clust)
  sds <- sapply(1:length(b1_avg), function(x) sd(c(b1_avg[x], c1_avg[x])))
  results <- cbind(results, sds)
  colnames(results)[i+2] <- i
}
write.table(results, paste(rna_path, "/results/all_genes_b1_c1_var_prop_expr_per_cluster.tsv", sep=""), sep = "\t", row.names = FALSE, quote = FALSE)


# Average Expr Lvl per Cell B2 v C1
results <- data.frame(genes <- gene_names)
colnames(results)[1] <- "gene"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
for (i in 0:num_clusters) {
  b2_this_clust <- WhichCells(b2, idents = i)
  c1_this_clust <- WhichCells(c1, idents = i)
  b2_sums <- rowSums(b2@assays$RNA@counts[gene_names,b2_this_clust])
  c1_sums <- rowSums(as.matrix(c1@assays$RNA@counts[gene_names,c1_this_clust]))
  b2_avg <- b2_sums/length(b2_this_clust)
  c1_avg <- c1_sums/length(c1_this_clust)
  sds <- sapply(1:length(b2_avg), function(x) sd(c(b2_avg[x], c1_avg[x])))
  results <- cbind(results, sds)
  colnames(results)[i+2] <- i
}
write.table(results, paste(rna_path, "/results/all_genes_b2_c1_var_expr_lvl_per_cluster.tsv", sep=""), sep = "\t", row.names = FALSE, quote = FALSE)

# Proportion of cells displaying the transcript B2 v C1
Idents(b2) <- b2$seurat_clusters
Idents(c1) <- c1$seurat_clusters
results <- data.frame(genes <- gene_names)
colnames(results)[1] <- "gene"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
for (i in 0:num_clusters) {
  b2_this_clust <- WhichCells(b2, idents = i)
  c1_this_clust <- WhichCells(c1, idents = i)
  b2_sums <- rowSums(b2@assays$RNA@counts[gene_names,b2_this_clust] != 0) # num cells expr this gene
  c1_sums <- rowSums(as.matrix(c1@assays$RNA@counts[gene_names,c1_this_clust] != 0)) # num cells expr this gene
  b2_avg <- b2_sums/length(b2_this_clust)
  c1_avg <- c1_sums/length(c1_this_clust)
  sds <- sapply(1:length(b2_avg), function(x) sd(c(b2_avg[x], c1_avg[x])))
  results <- cbind(results, sds)
  colnames(results)[i+2] <- i
}
write.table(results, paste(rna_path, "/results/all_genes_b2_c1_var_prop_expr_per_cluster.tsv", sep=""), sep = "\t", row.names = FALSE, quote = FALSE)

# Unused
# results <- data.frame()
# for (gene in test_gene_names) {
#   new_row <- c(gene)
#   for (i in 0:num_clusters) {
#     this_clust <- cells_in_clust[[i+1]]
#     this_sum <- sum(combined@assays$RNA@counts[ gene, this_clust ])
#     new_row <- c(new_row, this_sum/length(this_clust))
#   }
#   results <- rbind(results, t(new_row))
# }