library("dplyr")
library("Matrix")
library("Seurat")
library("stringr")
library("ggplot2")

# rna_path <- "C:/Users/miles/Downloads/brain/"
rna_path <- "~/scratch/brain/"
source(paste0(rna_path, "/brain_scripts/all_f.R"))

# combined <- readRDS(paste(rna_path, "/brain_scripts/brain_shiny/data/combined.rds", sep = ""))
combined <- readRDS(paste(rna_path, "/data/B1C1C2MZ_combined_031020.rds", sep=""))
marker_path <- paste(rna_path, "data/markers/", sep="")
marker_files <- dir(marker_path, pattern =paste("*.txt", sep=""))

markers <- data.frame(gene <- c(), bio <- c())
for (i in 1:length(marker_files)) {
  file <- read.table(paste(marker_path, marker_files[i], sep=""), header=FALSE, sep="\t", stringsAsFactors=FALSE)
  markers <- rbind(markers, file[,1:2])
}
colnames(markers) <- c("gene", "bio")
bio <- "ROCK_SAND"
markers <- markers[which(markers$bio == bio),]
gene_names <- rownames(combined@assays$RNA)
marker_genes <- markers$gene
valid_genes <- marker_genes
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
down_avg_avg_gene <- rep(0, num_clusters+1)
run_num <- 50

# No Perm, Bootstrap
down_list <- lapply(0:num_clusters, function(x) c())
no_list <- lapply(0:num_clusters, function(x) c())
mean_down_list <- lapply(0:num_clusters, function(x) c())
mean_no_list <- lapply(0:num_clusters, function(x) c())
for (run in 1:run_num) {
  cat(paste("no_perm", run, "\n"))
  mat <- downsample(combined, marker_genes, run)
  
  cells_per_cluster <- c()
  genes_per_cluster <- c()
  for (i in 0:num_clusters) {
    this_cells <- WhichCells(combined, idents = i)
    if (run == 1) {
      mean_down_list[[i+1]] <- rep(0, length(this_cells))
      mean_no_list[[i+1]] <- rep(0, length(this_cells))
    }
    this_genes <- sapply(marker_genes, function(gene) length(rownames(mat)[which(mat[gene, this_cells] > 0)])) # number of transcripts expressed by cells in the cluster
    down_list[[i+1]] <- c(down_list[[i+1]], this_genes)
    no_list[[i+1]] <- c(no_list[[i+1]], length(this_cells) - this_genes)
    mean_down_list[[i+1]] <- mean_down_list[[i+1]] + this_genes
    mean_no_list[[i+1]] <- mean_no_list[[i+1]] + (length(this_cells) - this_genes)
  }
}
mean_down_list <- lapply(1:length(mean_down_list), function(x) mean_down_list[[x]]/run_num)
mean_no_list <- lapply(1:length(mean_no_list), function(x) mean_no_list[[x]]/run_num)

df <- data.frame()
for (i in 1:(num_clusters+1)) {
  all_other_clusters <- c(1:(num_clusters+1))[-i]
  wilcox_p <- wilcox.test(x = down_list[[i]], y  = unlist(down_list[all_other_clusters]), alternative = "greater")$p.value # wilcox to compare cluster to all other cells
  t_p      <- t.test(     x = down_list[[i]], y  = unlist(down_list[all_other_clusters]), alternative = "greater")$p.value
  
  # this_cells <- WhichCells(combined, idents = (i-1))
  this_cluster <- c(sum(mean_down_list[[i]]),                        sum(mean_no_list[[i]]))
  all_clusters <- c(sum(unlist(mean_down_list[all_other_clusters])), sum(unlist(mean_no_list[all_other_clusters])))
  contig_table <- data.frame(this_cluster <- this_cluster, all_clusters <- all_clusters)
  fisher_p <- fisher.test(contig_table, alternative = "greater")$p.value
  
  df <- rbind(df, t(c(i, wilcox_p, t_p, fisher_p)))
}
colnames(df) <- c("cluster", "wilcox_p", "t_p", "fisher_p")
df$wilcox_q <- p.adjust(df$wilcox_p, method = "bonferroni")
df$t_q <- p.adjust(df$t_p, method = "bonferroni")
df$fisher_q <- p.adjust(df$fisher_p, method = "bonferroni")
# df$cluster[which(df$wilcox_q < 0.05)]
# df$cluster[which(df$t_q < 0.05)]
df$cluster[which(df$fisher_q < 0.05)]
write.table(df, file = paste(rna_path, "/results/down_3_", bio, ".tsv", sep=""), sep = "\t", row.names = FALSE, quote=FALSE)
