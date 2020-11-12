
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
    trans_per_cluster <- c(trans_per_cluster, sum(colSums(combined@assays$RNA@counts[,this_cells]))) # trans
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
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
valid_genes <- validGenes(markers$gene, gene_names)
num_genes_per_cluster <- c()
for (run in 101:200) {
  if (run%%10 == 0) {
    print(run)
  }
  set.seed(run)
  shuffled <- sample(backup_ids)
  
  Idents(object = combined) <- shuffled
  # gene_names_per_cluster2 <- list()
  for (i in 0:num_clusters) {
    this_cells <- WhichCells(combined, idents = i)
    threshold <- 0.25 * length(this_cells)
    this_gene_names <- rownames(combined@assays$RNA@counts[which(rowSums(as.matrix(combined@assays$RNA@counts[valid_genes, this_cells])) > threshold),])
    num_genes_per_cluster <- c(num_genes_per_cluster, length(this_gene_names))
    # gene_names_per_cluster2[[i+1]] <- this_gene_names
  }
}
quantile(num_genes_per_cluster, c(.0275, .975))