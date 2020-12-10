library("Seurat")
library("Matrix")
library("qvalue")
library("jaccard")
bb <- readRDS("~/scratch/brain/data/bb_clustered_102820.rds")
bb$seurat_clusters = bb$seuratclusters15
Idents(bb) = bb$seurat_clusters
obj <- bb
gene_names <- rownames(obj)[which(rowSums(as.matrix(obj@assays$RNA@counts)) != 0)]

# ###############
# # Transcripts #
# ###############
setToMax = function(dat, x) {
  dat[which(dat > x)] = x
  return(dat)
}

for (cluster in 0:14) {
  this_cells = WhichCells(obj, idents = cluster)
  mat_trans  <- matrix(0L, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
  mat2_trans <- matrix(0L, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
  
  for (cell in this_cells) {
    # for (col in 1:100) {
    if (col %% 100 == 0) {
      print(col)
    }
    dat = obj@assays$RNA@counts[,cell]
    non_zero_genes = names(dat[which(dat > 0)])
    tmp <- mat_trans[non_zero_genes, non_zero_genes]
    dat = obj@assays$RNA@counts[non_zero_genes, col]
    this_trans_mat = t(as.matrix(sapply(1:length(non_zero_genes), function(x) 2*setToMax(dat,dat[x]))))
    mat_trans[non_zero_genes, non_zero_genes] = tmp + this_trans_mat
  }
  
  gene_trans = rowSums(obj@assays$RNA@counts[gene_names,])
  mat2_trans = t(as.matrix(sapply(1:length(gene_names), function(x) gene_trans + gene_trans[x])))
  
  mat3_trans = mat_trans/mat2_trans
  mat3_trans_p = matrix(jaccard.rahman(as.vector(mat3_trans)), length(gene_names), length(gene_names), dimnames=list(gene_names, gene_names))
  print("Done Filling Matrices")
  
  saveRDS(mat3_trans, paste0("~/scratch/brain/data/bb_j", cluster, ".RDS"))
  saveRDS(mat3_trans_p, paste0("~/scratch/brain/data/bb_j_p", cluster, ".RDS"))
}