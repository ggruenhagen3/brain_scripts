
findAllTrans <- function(combined, new_cells) {
  # Avg trans per cluster for ALL genes
  num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
  results_expr <- data.frame()
  cells_per_cluster <- c()
  trans_per_cluster <- c()
  for (i in 0:num_clusters) {
    # this_cells <- WhichCells(combined, idents = i)
    this_cells <- new_cells[[i+1]]
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

pickNewCells <- function(combined, num_clusters, num_cells) {
  new_cells <- c()
  for (i in 1:num_cells) {
    ran_cluster <- sample(0:num_clusters, 1)
    this_cells <- names(combined$seurat_clusters[which(combined$seurat_clusters == ran_cluster)])
    new_cells <- c(new_cells, sample(this_cells,1))
  }
  
  return(new_cells)
}

shuffleClusters <- function(combined) {
  # The selection process for a new cluster should be as follows:
  # 1. Pick a random cluster 0-40
  # 2. Pick a random cell from that cluster to be a part of the new cluster
  # This means that the new data set would likely have duplicate cells
  new_cells <- lapply(0:num_clusters, function(x) c())
  num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
  for (i in 0:num_clusters) {
    num_cells <- length(combined$seurat_clusters[which(combined$seurat_clusters == i)])
    new_cells[[i+1]] <- pickNewCells(combined, num_clusters, num_cells)
  }
  
  return(new_cells)
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
  
  # Idents(object = combined) <- shuffled
  new_cells <- shuffleClusters(combined)
  num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
  gene_names <- rownames(combined@assays$RNA)
  cells_per_cluster <- c()
  genes_per_cluster <- c()
  all_trans <- findAllTrans(combined, new_cells)
  all_all_trans <- c(all_all_trans, all_trans[,1])
  for (i in 0:num_clusters) {
    # this_cells <- WhichCells(combined, idents = i)
    this_cells <- new_cells[[i+1]]
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