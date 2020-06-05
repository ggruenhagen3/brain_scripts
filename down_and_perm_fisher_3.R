library("dplyr")
library("Matrix")
library("Seurat")
library("stringr")

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
total_genes_per_cluster <- rep(0, num_clusters+1)
run_num <- 50

# No Perm, Bootstrap
down_pos <- lapply(0:num_clusters, function(x) 0)
down_neg <- lapply(0:num_clusters, function(x) 0)
for (run in 1:run_num) {
  cat(paste("no_perm", run, "\n"))
  mat <- downsample(combined, marker_genes, run)
  
  for (i in 0:num_clusters) {
    this_cells <- WhichCells(combined, idents = i)
    this_genes <- length(which(as.vector(mat[,this_cells]) != 0))
    this_neg   <- length(which(as.vector(mat[,this_cells]) == 0))
    # genes_per_cluster <- c(genes_per_cluster, length(which(as.vector(combined@assays$RNA@counts[ran_markers,this_cells]) != 0))) # genes
    down_pos[[i+1]] <- down_pos[[i+1]] + this_genes
    down_neg[[i+1]] <- down_neg[[i+1]] + this_neg
  }
}
down_avg_pos <- lapply(1:length(down_pos), function(x) down_pos[[x]]/50)
down_avg_neg <- lapply(1:length(down_neg), function(x) down_neg[[x]]/50)
print(down_avg_pos)

# Perm, Bootstrap
backup_ids <- combined@meta.data$seurat_clusters
perm_pos <- lapply(0:num_clusters, function(x) 0)
perm_neg <- lapply(0:num_clusters, function(x) 0)
for (run in (run_num+1):(run_num+run_num)) {
  cat(paste("perm", run, "\n"))
  set.seed(run)
  mat <- downsample(combined, marker_genes, run)
  
  new_cells <- shuffleClusters(combined)
  num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
  gene_names <- rownames(combined@assays$RNA)
  for (i in 0:num_clusters) {
    this_cells <- new_cells[[i+1]]
    this_genes <- length(which(as.vector(mat[valid_genes,this_cells]) != 0))
    this_neg   <- length(which(as.vector(mat[valid_genes,this_cells]) == 0))
    perm_pos[[i+1]] <- perm_pos[[i+1]] + this_genes
    perm_neg[[i+1]] <- perm_neg[[i+1]] + this_neg
  }
}
perm_avg_pos <- lapply(1:length(perm_pos), function(x) perm_pos[[x]]/50)
perm_avg_neg <- lapply(1:length(perm_neg), function(x) perm_neg[[x]]/50)

# Compare empirical data to the permutated data on a PER CLUSTER basis
df <- data.frame()
for (i in 0:num_clusters) {
  down <- c(down_avg_pos[[i+1]], down_avg_neg[[i+1]])
  perm <- c(perm_avg_pos[[i+1]], perm_avg_neg[[i+1]])
  contig_table <- data.frame(down, perm)
  fisher_p <- fisher.test(contig_table, alternative = "greater")$p.value
  df <- rbind(df, t(c(i, down, perm, fisher_p)) )
}
df$fisher_q <- p.adjust(df[,6], method = "hochberg")
df$q_sig <- df$fisher_q < 0.05

# colnames(df) <- c("cluster", "down_avg_gene", "down_num_cells", "down")
# df <- as.data.frame(df)
# df$avg_gene_per_cell_per_cluster <- as.numeric(as.vector(df$avg_gene_per_cell_per_cluster))
# df$cluster <- factor(df$cluster, levels = 1:(num_clusters+1))
# png(paste0(rna_path, "/results/down_and_perm_2_", bio, ".png"), width = 1800, height = 1000, res = 150)
# p <- ggplot(df, aes(x = cluster, y = avg_gene_per_cell_per_cluster, fill = cond)) + geom_boxplot(alpha = 0.6) + geom_jitter(shape=16, position=position_jitterdodge(), alpha = 0.3, aes(colour = cond)) + scale_colour_manual(values=c("#999999", "#56B4E9", "#3ac9bb")) + scale_fill_manual(values=c("#999999", "#3ac9bb", "#56B4E9")) + ggtitle(paste(bio, "- Average Genes per Cell per Cluster"))
# print(p)
# dev.off()

write.table(df, file = paste(rna_path, "/results/down_and_perm_fisher_3_", bio, ".tsv", sep=""), sep = "\t", row.names = FALSE, quote=FALSE)
