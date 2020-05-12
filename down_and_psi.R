library("dplyr")
library("Matrix")
library("Seurat")
library("stringr")
library("ggplot2")
library("pSI")

# rna_path <- "C:/Users/miles/Downloads/brain/"
rna_path <- "~/scratch/brain/"
source(paste0(rna_path, "/brain_scripts/all_f.R"))

combined <- readRDS(paste(rna_path, "/brain_scripts/brain_shiny/data/combined.rds", sep = ""))
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
run_num <- 3
big_mat  <- matrix(0L, nrow=length(marker_genes), ncol = ncol(combined@assays$RNA@counts), dimnames = list(marker_genes, colnames(combined@assays$RNA@counts)))

# No Perm, Bootstrap
down_list <- lapply(0:num_clusters, function(x) c())
for (run in 1:run_num) {
  cat(paste("no_perm", run, "\n"))
  mat <- downsample(combined, marker_genes, run)
  big_mat <- big_mat + mat
}
new_obj <- CreateSeuratObject(counts = big_mat)
new_obj$seurat_clusters <- combined$seurat_clusters
Idents(new_obj) <- new_obj$seurat_clusters
cluster_df <- myAverageExpression(new_obj, slot = "counts")
results <- specificity.index(cluster_df, e_min = 0.05)
colnames(results) <- 0:num_clusters
fisher_results <- fisher.iteration(results, valid_genes, p.adjust = FALSE)
fisher_results$q <- p.adjust(fisher_results$`0.05 - nominal`)
fisher_results$q.sig <- fisher_results$q < 0.05

write.table(results, file = paste(rna_path, "/results/down_and_psi_raw_", bio, ".tsv", sep=""), sep = "\t", row.names = FALSE, quote=FALSE)
write.table(fisher_results, file = paste(rna_path, "/results/down_and_psi_", bio, ".tsv", sep=""), sep = "\t", row.names = FALSE, quote=FALSE)