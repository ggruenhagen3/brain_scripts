# Read Packages
library("Seurat")
library("Matrix")
library("matrixTests")

# Read in Seurat Object
bb <- readRDS("~/scratch/brain/data/bb_clustered_102820.rds")
obj <- bb
gene_names <- rownames(obj)[which(rowSums(as.matrix(obj@assays$RNA@counts)) > 1)]
exp = GetAssayData(obj, assay = "RNA", slot="data")
exp = as.matrix(exp)

# Find the positive cells for each gene
gene_cells = lapply(gene_names, function(x) c())
names(gene_cells) = gene_names
all_cells = colnames(obj)

for (gene in gene_names) {
  gene_cells[[gene]] = colnames(obj)[which(obj@assays$RNA@counts[gene,] != 0)]
}

# Find the Expression of Genes in Gene1+ cells
# Find the Expression of Genes in Gene1- cells
# Find the Difference in expression between those two
mat_pos = matrix(0L, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
mat_neg = matrix(0L, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
mat_dif = matrix(0L, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
mat_p   = matrix(0L, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))

for (i in 1:length(gene_names)) {
  print(i)
  gene1 = gene_names[i]
  gene1_cells = gene_cells[[gene1]]
  pos = rowMeans(exp[gene_names, colnames(obj)[which(  colnames(obj) %in% gene1_cells )]]) # expression of all genes in gene1 positive cells
  neg = rowMeans(exp[gene_names, colnames(obj)[which(! colnames(obj) %in% gene1_cells )]]) # expression of all genes in gene1 negative cells
  mat_pos[i,] = pos
  mat_neg[i,] = neg
  mat_dif[i,] = pos - neg # difference in expression of all genes in gene1 positive vs negative cells
  # mat_p[i,] = row_kruskalwallis(exp[gene_names,], g = colnames(obj) %in% gene1_cells)$pvalue
}

saveRDS(mat_pos, "~/scratch/brain/data/bb_co_pos.RDS")
saveRDS(mat_neg, "~/scratch/brain/data/bb_co_neg.RDS")
saveRDS(mat_dif, "~/scratch/brain/data/bb_co_dif.RDS")
saveRDS(mat_p, "~/scratch/brain/data/bb_co_dif_p.RDS")