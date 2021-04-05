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
  
  col = 0
  for (cell in this_cells) {
    # for (col in 1:100) {
    if (col %% 100 == 0) {
      print(col)
    }
    dat = obj@assays$RNA@counts[,cell]
    non_zero_genes = names(dat[which(dat > 0)])
    tmp <- mat_trans[non_zero_genes, non_zero_genes]
    dat = obj@assays$RNA@counts[non_zero_genes, cell]
    this_trans_mat = t(as.matrix(sapply(1:length(non_zero_genes), function(x) 2*setToMax(dat,dat[x]))))
    mat_trans[non_zero_genes, non_zero_genes] = tmp + this_trans_mat
    col = col + 1
  }
  
  gene_trans = rowSums(obj@assays$RNA@counts[gene_names,])
  mat2_trans = t(as.matrix(sapply(1:length(gene_names), function(x) gene_trans + gene_trans[x])))
  
  mat3_trans = mat_trans/mat2_trans
  mat3_trans_p = matrix(jaccard.rahman(as.vector(mat3_trans)), length(gene_names), length(gene_names), dimnames=list(gene_names, gene_names))
  print("Done Filling Matrices")
  
  saveRDS(mat3_trans, paste0("~/scratch/brain/data/bb_j", cluster, ".RDS"))
  saveRDS(mat3_trans_p, paste0("~/scratch/brain/data/bb_j_p", cluster, ".RDS"))
}

############
# Analysis #
############
for (cluster in 0:14) {
  
  mat_j = readRDS(paste0("~/scratch/brain/data/bb_j", cluster, ".RDS"))
  mat_j_p = readRDS(paste0("~/scratch/brain/data/bb_j_p", cluster, ".RDS"))
  
  mat = mat_j_p
  sig_df = data.frame()
  for (row in 1:(length(gene_names)-1)) {
    if (row %% 1000 == 0) { print(row) }
    j = mat_j[row, (row+1):ncol(mat)]
    p = mat[row, (row+1):ncol(mat)]
    q = p.adjust(p, method = "bonferroni")
    
    sig_ind = which(q < 0.05)
    if (length(sig_ind) > 0) {
      newRow = data.frame(gene = gene_names[row], q = q[sig_ind], p = p[sig_ind], j = j[sig_ind], close_gene = names(q[sig_ind]))
      sig_df = rbind(sig_df, newRow)
    }
  }
  write.csv(sig_df, paste0("~/scratch/brain/results/bb_j_sig", cluster, ".csv"), row.names = F)
  
}
