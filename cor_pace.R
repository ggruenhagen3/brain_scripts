library("Seurat")
library("Matrix")
library("qvalue")
library("jaccard")
library("parallel")
print(paste("Num Cors:", detectCores()))
# bb <- readRDS("~/scratch/brain/data/bb_clustered_102820.rds")
paul_epi = readRDS("~/scratch/d_tooth/data/epi_full.rds")
paul_epi$cond[is.na(paul_epi$cond)] <- "INJR"
paul_ctrl = subset(paul_epi, cells = colnames(paul_epi)[which(paul_epi$cond == "CTRL")])
obj <- paul_ctrl
gene_names <- rownames(obj)[which(rowSums(as.matrix(obj@assays$RNA@counts)) > 2)]
mat_data_p   = matrix(0, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
mat_data_cor = matrix(0, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
# test = as.matrix(t(obj@assays$RNA@data))
# ptm <- proc.time()
# cl <- makeCluster(40, outfile="/nv/hp10/ggruenhagen3/scratch/brain/brain_scripts/cor_pace.log")
# registerDoParallel(cl)
# # One for pvalues
# df_p  = foreach(col = 2:length(gene_names), .combine='rbind') %dopar% {
#   gene1 = gene_names[col]
#   thisRow = c()
#   for ( row in 1:(col-1) ) {
#     gene2 = gene_names[row]
#     cor_res = cor.test(obj@assays$RNA@data[gene1,], obj@assays$RNA@data[gene2,])
#     thisRow = c(thisRow, cor_res$p.value)
#     # mat_data_p[row,col]   = cor_res$p.value
#     # mat_data_cor[row,col] = cor_res$estimate
#   }
#   thisRow = c(thisRow, rep(NA, length(gene_names) - length(thisRow))) # fill the rest with NA's
#   print(paste("Task", col, " done"))
#   thisRow
# }
# saveRDS(df_p,   "/nv/hp10/ggruenhagen3/scratch/brain/data/df_p.RDS")
# print("Saved p-values")
# # One for correlation
# df_p  = foreach(col = 2:length(gene_names), .combine='rbind') %dopar% {
#   gene1 = gene_names[col]
#   thisRow = c()
#   for ( row in 1:(col-1) ) {
#     gene2 = gene_names[row]
#     cor_res = cor.test(obj@assays$RNA@data[gene1,], obj@assays$RNA@data[gene2,])
#     thisRow = c(thisRow, cor_res$estimate)
#     # mat_data_p[row,col]   = cor_res$p.value
#     # mat_data_cor[row,col] = cor_res$estimate
#   }
#   thisRow = c(thisRow, rep(NA, length(gene_names) - length(thisRow))) # fill the rest with NA's
#   print(paste("Task", col, " done"))
#   thisRow
# }
# stopCluster(cl)
# proc.time() - ptm
# saveRDS(df_cor, "/nv/hp10/ggruenhagen3/scratch/brain/data/df_cor.RDS")
# print("Saved correlations")

# for (col in 2:length(gene_names)) {
#   if (col %% 100 == 0) {
#     print(col)
#   }
#   gene1 = gene_names[col]
#   for ( row in 1:(col-1) ) {
#     gene2 = gene_names[row]
#     cor_res = cor.test(obj@assays$RNA@data[gene1,], obj@assays$RNA@data[gene2,])
#     # cor_res = cor.test(formula = ~ gene1 + gene2, data=test))
#     mat_data_p[row,col]   = cor_res$p.value
#     mat_data_cor[row,col] = cor_res$estimate
#   }
# }

# system.time(cor.pvalues <- future_sapply(2:100, function(x) { sapply(1:(x-1), function(y) { cor.test(test[gene_names[x],],test[gene_names[y],])$p.value })}))

# Attempt as of 02/10/2021
my_cor_t = function(r, n) (r * sqrt(n - 2))/sqrt(1 - r**2)
my_cor_p = function(t, n) 2*pt(-abs(t), df=n-2)
fx = function(x) cor(t(as.matrix(obj@assays$RNA@data[x,])), y = NULL)

r_mat = cor(t(as.matrix(obj@assays$RNA@data[,])), y = NULL)
t_mat = my_cor_t(r_mat, ncol(obj))
p_mat = my_cor_p(t_mat, ncol(obj))

saveRDS(r_mat, "~/scratch/d_tooth/data/paul_epi_ctrl_cor.RDS")
saveRDS(p_mat,   "~/scratch/d_tooth/data/paul_epi_ctrl_cor_p.RDS")

# gene_names <- rownames(bb)[which(rowSums(bb@assays$RNA@counts) > 2)]
# fx <- function(gene1) cor.test(bb@assays$RNA@data[gene1, ], bb@assays$RNA@data["egr1", ])$p.value
# print("No Parallel:")
# system.time(
#   results <- lapply(gene_names[1:200], fx)
# )
# 
# print("Yes Parallel:")
# numCores <- detectCores()
# system.time(
#   results <- mclapply(gene_names, fx, mc.cores = numCores)
# )
