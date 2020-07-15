library("Seurat")
library("Matrix")
library("qvalue")
library("jaccard")
library("foreach")
library("doParallel")
lncRNA <- readRDS("/nv/hp10/ggruenhagen3/scratch/brain/data/lncRNA.RDS")
obj <- lncRNA
gene_names <- rownames(obj)[which(rowSums(as.matrix(obj@assays$RNA@counts)) != 0)]
mat_data_p   = matrix(0, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
mat_data_cor = matrix(0, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
# test = as.matrix(t(obj@assays$RNA@data))
test = obj@assays$RNA@data
ptm <- proc.time()
cl <- makeCluster(8)
registerDoParallel(cl)
# One for pvalues
df_p  = foreach(col = 2:500, .combine='rbind') %dopar% {
  gene1 = gene_names[col]
  thisRow = c()
  for ( row in 1:(col-1) ) {
    gene2 = gene_names[row]
    cor_res = cor.test(obj@assays$RNA@data[gene1,], obj@assays$RNA@data[gene2,])
    thisRow = c(thisRow, cor_res$p.value)
    # mat_data_p[row,col]   = cor_res$p.value
    # mat_data_cor[row,col] = cor_res$estimate
  }
  thisRow = c(thisRow, rep(0, length(gene_names) - length(thisRow))) # fill the rest with 0's
  thisRow
}
# One for correlation
df_p  = foreach(col = 2:500, .combine='rbind') %dopar% {
  gene1 = gene_names[col]
  thisRow = c()
  for ( row in 1:(col-1) ) {
    gene2 = gene_names[row]
    cor_res = cor.test(obj@assays$RNA@data[gene1,], obj@assays$RNA@data[gene2,])
    thisRow = c(thisRow, cor_res$estimate)
    # mat_data_p[row,col]   = cor_res$p.value
    # mat_data_cor[row,col] = cor_res$estimate
  }
  thisRow = c(thisRow, rep(0, length(gene_names) - length(thisRow))) # fill the rest with 0's
  thisRow
}
stopCluster(cl)
proc.time() - ptm
# saveRDS(df_p,   "/nv/hp10/ggruenhagen3/scratch/brain/data/df_p.RDS")
# saveRDS(df_cor, "/nv/hp10/ggruenhagen3/scratch/brain/data/df_cor.RDS")

ptm <- proc.time()
for (col in 2:500) {
  gene1 = gene_names[col]
  for ( row in 1:(col-1) ) {
    gene2 = gene_names[row]
    cor_res = cor.test(obj@assays$RNA@data[gene1,], obj@assays$RNA@data[gene2,])
    mat_data_p[row,col]   = cor_res$p.value
    mat_data_cor[row,col] = cor_res$estimate
  }
}
proc.time() - ptm

# for (col in 2:length(gene_names)) {
#   if (col %% 100 == 0) {
#     print(col)
#   }
#   gene1 = gene_names[col]
#   for ( row in 1:(col-1) ) {
#     gene2 = gene_names[row]
#     cor_res = cor.test(obj@assays$RNA@data[gene1,], obj@assays$RNA@data[gene2,])
#     # cor_res = cor.test(formula = ~ gene1 + gene2, data=test))
#     # mat_data_p[row,col]   = cor_res$p.value
#     # mat_data_cor[row,col] = cor_res$estimate
#   }
# }

# system.time(cor.pvalues <- future_sapply(2:100, function(x) { sapply(1:(x-1), function(y) { cor.test(test[gene_names[x],],test[gene_names[y],])$p.value })}))

saveRDS(mat_data_p,   "/nv/hp10/ggruenhagen3/scratch/brain/data/mat_data_p.RDS")
saveRDS(mat_data_cor, "/nv/hp10/ggruenhagen3/scratch/brain/data/mat_data_cor.RDS")