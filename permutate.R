
findAllTrans <- function(combined, shuffled) {
  # Avg trans per cluster for ALL genes
  Idents(object = combined) <- shuffled
  num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
  results_expr <- data.frame()
  cells_per_cluster <- c()
  trans_per_cluster <- c()
  for (i in 0:num_clusters) {
    this_cells <- WhichCells(combined, idents = i)
    # genes_per_cluster <- c(genes_per_cluster, length(which(as.vector(combined@assays$RNA@counts[,this_cells]) != 0))) # genes
    trans_per_cluster <- c(trans_per_cluster, sum(colSums(as.matrix(combined@assays$RNA@counts[,this_cells])))) # trans
    cells_per_cluster <- c(cells_per_cluster, length(this_cells))
  }
  avg_trans_per_cell_per_cluster <- trans_per_cluster/cells_per_cluster
  all_trans <- as.data.frame(avg_trans_per_cell_per_cluster)
  all_trans$trans_per_cluster <- trans_per_cluster
  all_trans$cells_per_cluster <- cells_per_cluster
  rownames(all_trans) <- 0:40
  
  return(all_trans)
}


# Bootstrap
backup_ids <- combined@meta.data$seurat_clusters
shuffled_avg_gene <- c()
shuffled_avg_gene_by_all_trans <- c()
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
all_all_trans <- c()
fold_diff_all_trans <- c()
perm_down_avg_gene <- lapply(0:num_clusters, function(x) c())
for (run in 101:200) {
  if (run%%10 == 0) {
    print(run)
  }
  set.seed(run)
  shuffled <- sample(backup_ids)
  
  Idents(object = combined) <- shuffled
  num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
  gene_names <- rownames(combined@assays$RNA)
  cells_per_cluster <- c()
  genes_per_cluster <- c()
  all_trans <- findAllTrans(combined, shuffled)
  all_all_trans <- c(all_all_trans, all_trans[,1])
  for (i in 0:num_clusters) {
    this_cells <- WhichCells(combined, idents = i)
    genes_per_cluster <- c(genes_per_cluster, length(which(as.vector(combined@assays$RNA@counts[valid_genes,this_cells]) != 0))) # genes
    cells_per_cluster <- c(cells_per_cluster, length(this_cells))
    avg_gene <- length(which(as.vector(combined@assays$RNA@counts[valid_genes,this_cells]) != 0))
    perm_down_avg_gene[[i+1]] <- c(perm_down_avg_gene[[i+1]], avg_gene/all_trans[i+1,2])
  }
  avg_gene_per_cell_per_cluster <- genes_per_cluster/cells_per_cluster
  shuffled_avg_gene <- c(shuffled_avg_gene, avg_gene_per_cell_per_cluster)
  fold_diff_all_trans <- c(fold_diff_all_trans, max(all_trans[,1])/min(all_trans[,1]))
  shuffled_avg_gene_by_all_trans <- c(shuffled_avg_gene_by_all_trans, avg_gene_per_cell_per_cluster/all_trans[,1])
}
quantile(shuffled_avg_gene, c(.0275, .975))
quantile(shuffled_avg_gene_by_all_trans, c(.0275, .975))

Idents(object = combined) <- backup_ids
all_trans <- read.csv("C:/Users/miles/Downloads/brain/results/all_avg_gene_per_cluster_trans.csv", stringsAsFactors = FALSE)
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
gene_names <- rownames(combined@assays$RNA)
cells_per_cluster <- c()
genes_per_cluster <- c()
for (i in 0:num_clusters) {
  this_cells <- WhichCells(combined, idents = i)
  genes_per_cluster <- c(genes_per_cluster, length(which(as.vector(combined@assays$RNA@counts[valid_genes,this_cells]) != 0))) # genes
  cells_per_cluster <- c(cells_per_cluster, length(this_cells))
}
avg_gene_per_cell_per_cluster <- genes_per_cluster/cells_per_cluster
# avg_gene_per_cell_per_cluster <- genes_per_cluster/all_trans[,2]
df <- as.data.frame(avg_gene_per_cell_per_cluster)
df$genes_per_cluster <- genes_per_cluster
df$cells_per_cluster <- cells_per_cluster
rownames(df) <- 0:num_clusters

df2 <- data.frame()
sig_clusters <- c()
for (i in 0:num_clusters) {
  isSig <- FALSE
  sig <- quantile(perm_down_avg_gene[[i+1]], c(0.975))
  if ( df[i+1,2]/all_trans[i+1,3] > sig ) {
    # print(i)
    # print(sig)
    sig_clusters <- c(sig_clusters, i)
    isSig <- TRUE
  }
  df2 <- rbind(df2, t(c(i, df[i+1,1]/all_trans[i+1,2], sig, isSig)))
}
print(sig_clusters)