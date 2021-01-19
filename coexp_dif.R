# Read Packages
library("Seurat")
library("Matrix")

# Read in Seurat Object
bb <- readRDS("~/scratch/brain/data/bb_clustered_102820.rds")
obj <- bb
gene_names <- rownames(obj)[which(rowSums(as.matrix(obj@assays$RNA@counts)) > 1)]
exp = GetAssayData(obj, assay = "RNA", slot="data")

# Find the positive cells for each gene
gene_cells = lapply(gene_names, function(x) c())
names(gene_cells) = gene_names
all_cells = colnames(obj)

for (gene in gene_names) {
  gene_cells[[gene]] = colnames(obj)[which(obj@assays$RNA@counts[gene,] != 0)]
}

#
mat_pos = matrix(0L, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
mat_neg = matrix(0L, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
mat_dif = matrix(0L, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
for (i in c(2)) {
  gene1 = gene_names[i]
  gene1_cells = gene_cells[[gene1]]
  pos = rowMeans(exp[,which(  colnames(obj) %in% gene1_cells )])
  neg = rowMeans(exp[,which(! colnames(obj) %in% gene1_cells )])
  mat_pos[i,] = pos
  mat_neg[i,] = neg
  mat_dif[i,] = pos - neg
}

saveRDS(mat_pos, "~/scratch/brain/data/bb_co_pos.RDS")
saveRDS(mat_neg, "~/scratch/brain/data/bb_co_neg.RDS")
saveRDS(mat_dif, "~/scratch/brain/data/bb_co_dif.RDS")