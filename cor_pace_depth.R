library("Seurat")
library("Matrix")
library("qvalue")
library("jaccard")
library("parallel")
print(paste("Num Cors:", detectCores()))
bb <- readRDS("~/scratch/brain/data/bb_clustered_102820.rds")
obj <- bb
gene_names <- rownames(obj)[which(rowSums(as.matrix(obj@assays$RNA@counts)) > 2)]
# mat_data_p   = matrix(0, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
# mat_data_cor = matrix(0, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))

numCores = detectCores()
depth_cor = mclapply(1:nrow(obj), function(x) cor(obj@assays$RNA@data[x,], obj$depth), mc.cores = numCores)

numCores = detectCores()
gsi_cor = mclapply(1:nrow(obj), function(x) cor(obj@assays$RNA@data[x,], obj$gsi), mc.cores = numCores)

df = data.frame(gene = rownames(obj), depth_cor = unlist(depth_cor), gsi_cor = unlist(gsi_cor))

# Attempt as of 02/10/2021
my_cor_t = function(r, n) (r * sqrt(n - 2))/sqrt(1 - r**2)
my_cor_p = function(t, n) 2*pt(-abs(t), df=n-2)
fx = function(x) cor(t(as.matrix(obj@assays$RNA@data[x,])), y = NULL)

r_mat = cor(t(as.matrix(obj@assays$RNA@data[,])), y = NULL)
t_mat = my_cor_t(r_mat, ncol(obj))
p_mat = my_cor_p(t_mat, ncol(obj))

saveRDS(r_mat, "~/scratch/d_tooth/data/hm_adult_cor.RDS")
saveRDS(p_mat,   "~/scratch/d_tooth/data/hm_adult_cor_p.RDS")
