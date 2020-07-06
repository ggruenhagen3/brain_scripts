library("Seurat")
library("Matrix")
library("qvalue")
library("jaccard")
lncRNA <- readRDS("/nv/hp10/ggruenhagen3/scratch/brain/data/lncRNA.RDS")
obj <- lncRNA
gene_names <- rownames(obj)[which(rowSums(as.matrix(obj@assays$RNA@counts)) != 0)]
mat_data_p   = matrix(0, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
mat_data_cor = matrix(0, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
for (col in 2:length(gene_names)) {
  if (col %% 100 == 0) {
    print(col)
  }
  gene1 = gene_names[1]
  for ( row in 1:(col-1) ) {
    gene2 = gene_names[2]
    cor_res = cor.test(obj@assays$RNA@data[gene1,], obj@assays$RNA@data[gene2,])
    mat_data_p[row,col]   = cor_res$p.value
    mat_data_cor[row,col] = cor_res$estimate
  }
}
saveRDS(mat_data_p,   "/nv/hp10/ggruenhagen3/scratch/brain/data/mat_data_p.RDS")
saveRDS(mat_data_cor, "/nv/hp10/ggruenhagen3/scratch/brain/data/mat_data_cor.RDS")