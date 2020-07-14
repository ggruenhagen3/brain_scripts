library("future.apply")
library("Seurat")
library("Matrix")
library("qvalue")
library("jaccard")
lncRNA <- readRDS("/nv/hp10/ggruenhagen3/scratch/brain/data/lncRNA.RDS")
obj <- lncRNA
gene_names <- rownames(obj)[which(rowSums(as.matrix(obj@assays$RNA@counts)) != 0)]
test = obj@assays$RNA@data
print("Parallelized:")
system.time(cor.pvalues <- future_sapply(2:100, function(x) { sapply(1:(x-1), function(y) { cor.test(test[gene_names[x],],test[gene_names[y],])$p.value })}))
print("Normal:")
system.time(cor.pvalues <- sapply(2:100, function(x) { sapply(1:(x-1), function(y) { cor.test(test[gene_names[x],],test[gene_names[y],])$p.value })}))
