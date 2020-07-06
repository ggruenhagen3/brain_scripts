library("Seurat")
library("Matrix")
library("reticulate")
library("cowplot")
library("biomaRt")
library("stringr")
library("dplyr")
library("ggplot2")

deg_file <- "/nv/hp10/ggruenhagen3/scratch/brain/data/lncRNA_lncRNA_deg.csv"
deg <- read.csv(deg_file, header = TRUE, stringsAsFactors = FALSE)
deg$cluster <- factor(deg$cluster)

lncRNA <- readRDS("/nv/hp10/ggruenhagen3/scratch/brain/data/lncRNA.RDS")
combined <- readRDS("/nv/hp10/ggruenhagen3/scratch/d_tooth/data/tj.rds")
lncRNA_genes <- rownames(lncRNA)[which(! rownames(lncRNA) %in% rownames(combined))]
gtf <- read.table("/nv/hp10/ggruenhagen3/scratch/brain/data/ens_w_ncbi_sort.gtf", sep = "\t", stringsAsFactors = FALSE)
non_zero_lncRNA = lncRNA_genes[which(rowSums(as.matrix(lncRNA@assays$RNA@counts[lncRNA_genes,])) != 0)]
num_clusters = 41

gtf_gene_name <- c()
for (i in 1:nrow(gtf)) {
  start <- gregexpr(pattern ='gene_name', gtf$V9[i])[[1]]
  stop  <- gregexpr(pattern =';', substr(gtf$V9[i], start, nchar(gtf$V9[i])))[[1]][1]
  gene_name <- substr(gtf$V9[i], start+10, start+stop-2)
  if (start == -1) {
    gene_name <- substr(gtf$V9[i], start+10, start+stop)
  }
  gtf_gene_name <- c(gtf_gene_name, gene_name)
}
gtf$gene_name <- gtf_gene_name
colnames(gtf) <- c("LG", "source", "type", "start", "stop", "idk", "idk1", "idk2", "info", "gene_name")
# gtf <- separate(gtf, col = 9, into = c("gene_id", "gene_name"), sep=";")
# gtf$gene_name[which(startsWith(gtf$gene_name, " gene_name"))] <- substr(gtf$gene_name[which(startsWith(gtf$gene_name, " gene_name"))], 12, 10000L)
deg_17 <- deg$X[which(deg$cluster == 17 & deg$avg_logFC > 0)]
neighbor_df <- data.frame()
neighbor_degree <- c(-2, -1, 1, 2)
for (deg in non_zero_lncRNA) {
  if (grepl(".", deg, fixed = TRUE)) {
    stop <- gregexpr(pattern ='\\.', deg)[[1]]
    deg <- substr(deg, 0, stop-1)
  }
  print(deg)
  deg_i <- as.numeric(rownames(gtf[which(gtf$gene_name == deg),]))
  neighbors <- list()
  for (degree in neighbor_degree) {
    if ( length(deg_i) > 0 && length(degree) > 0 && deg_i + degree >= 1 && deg_i + degree < length(gtf$gene_name)) { # weird errors on PACE, maybe it went to 0
      neighbors[as.character(degree)] = gtf$gene_name[deg_i + degree]
    }
  }
  if (all(neighbors %in% rownames(lncRNA))) {
    for (cluster in 0:num_clusters) {
      this_cells <- WhichCells(lncRNA, idents = cluster)
      # newRow <- data.frame()
      for (degree in neighbor_degree) {
        neighbor <- neighbors[[as.character(degree)]]
        this_mean <- mean(lncRNA@assays$RNA@counts[neighbor, this_cells])
        this_mean_norm <- this_mean/mean(lncRNA$nCount_RNA[this_cells])
        neighbor_df <- rbind(neighbor_df, t(c(deg, neighbor, degree, cluster, this_mean, this_mean_norm)))
      }
    } 
  } else {
    print(deg)
    for (degree in neighbor_degree) {
      print(neighbors[as.character(degree)])
    }
    print("------------------------")
  }
}
colnames(neighbor_df) <- c("deg", "neighbor", "degree", "cluster", "mean", "mean_norm")
neighbor_df$cluster <- factor(neighbor_df$cluster)
neighbor_df$degree <- factor(neighbor_df$degree)
neighbor_df$mean_norm <- as.numeric(as.vector(neighbor_df$mean_norm))

png(filename="/nv/hp10/ggruenhagen3/scratch/brain/results/lncRNA_neighbor_expr_norm_all_4.png", width = 1000, height = 600)
p = ggplot(neighbor_df, aes(cluster, mean_norm, fill = degree)) + stat_summary(fun = "mean", geom = "bar", position=position_dodge()) + ylab("Mean expression of Neighbor of All lncRNA / Mean expression of All Genes") + scale_fill_manual(values = c("#111d5e", "#c70039", "#f37121", "#ffbd69"))
print(p)
dev.off()

neighbor_df_2 <- data.frame()
for (cluster in 0:num_clusters) {
  this_rows <- neighbor_df[which(neighbor_df$cluster == cluster),]
  mean_norms <- c()
  for (degree in neighbor_degree) {
    cluster_mean_degree <- mean(this_rows$mean_norm[which(this_rows$degree == as.character(degree))])
    mean_norms <- c(mean_norms, cluster_mean_degree)
  }
  this_order <- sort(mean_norms, decreasing = TRUE, index.return = TRUE)$ix
  newRow <- data.frame( rep(cluster, length(neighbor_degree)), neighbor_degree[this_order], 1:length(neighbor_degree)  )
  # neighbor_df_2 <- rbind(neighbor_df_2, t(c(cluster, max_degree)))
  neighbor_df_2 <- rbind(neighbor_df_2, newRow)
}
colnames(neighbor_df_2) <- c("cluster", "degree", "rank")
neighbor_df_2$degree <- factor(neighbor_df_2$degree, levels=neighbor_degree)
neighbor_df_2$rank <- factor(neighbor_df_2$rank)


png(filename="/nv/hp10/ggruenhagen3/scratch/brain/results/lncRNA_all_neighbor_rank_4.png", width = 1000, height = 600)
p = ggplot(neighbor_df_2, aes(x=rank, fill = degree)) + geom_bar()  + scale_fill_manual(values = c("#111d5e", "#c70039", "#f37121", "#ffbd69")) + ggtitle("Ranking of Neighbors Across All Clusters From Highest to Lowest Expression")
print(p)
dev.off()

neighbor_df_3 <- data.frame()
for (degree in neighbor_degree) {
  this_rows <- neighbor_df[which(neighbor_df$degree == degree),]
  this_genes <- unique(this_rows$neighbor)
  this_mean <- colSums(lncRNA@assays$RNA@counts[this_genes,])
  this_mean_norm <- mean(this_mean/lncRNA$nCount_RNA)
  neighbor_df_3 <- rbind(neighbor_df_3, t(c(as.character(degree),this_mean_norm)))
}
colnames(neighbor_df_3) <- c("degree", "mean_norm")

png(filename="/nv/hp10/ggruenhagen3/scratch/brain/results/lncRNA_neighbor_expr_norm_avg_all_4.png", width = 1000, height = 600)
p = ggplot(neighbor_df_3, aes(x = degree, y=as.numeric(as.vector(mean_norm)), fill = factor(degree))) + geom_bar(stat="identity") + ylab("Mean expression of All Neighbors of lncRNA / Mean expression of All Genes") + scale_fill_manual(values = c("#111d5e", "#c70039", "#f37121", "#ffbd69")) + ggtitle("Average Expression for Each Neighbor Across All Cells")
print(p)
dev.off()

