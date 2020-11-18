# BiocManager::install("qvalue")
library("Seurat")
library("Matrix")
library("qvalue")
library("jaccard")
# iegs <- read.csv("C:/Users/miles/Downloads/zack_IEG_list_061720.csv", header = FALSE, stringsAsFactors = F)
# iegs <- iegs$V1
bb <- readRDS("~/scratch/brain/data/bb_clustered_102820.rds")
obj <- bb
gene_names <- rownames(obj)[which(rowSums(as.matrix(obj@assays$RNA@counts)) != 0)]


sig_df = read.csv("~/scratch/brain/results/bb_j_sig.csv", stringsAsFactors = F)

print("Finding Cells per Gene")
for (gene in sig_genes) {
  gene_cells[[gene]] = colnames(obj)[which(obj@assays$RNA@counts[gene,] != 0)]
}

print("Finding Correlations")
int = c()
my_cor = c()
my_cor_p = c()
for (i in 1:nrow(sig_df)) {
  if (row %% 10000 == 0) { print(row) }
  gene1 = sig_df$gene[i]
  gene2 = sig_df$close_gene[i]
  
  this_cor = cor.test(obj@assays$RNA@data[gene1,], obj@assays$RNA@data[gene2,])
  my_cor = c(my_cor, this_cor$estimate)
  my_cor_p = c(my_cor_p, this_cor$p.value)
  
  int = c(int, length(which(gene_cells[[gene1]] %in% gene_cells[[gene2]])))
}
sig_df$int = int
sig_df$cor = my_cor
sig_df$cor_p = my_cor_p

write.csv(sig_df, "~/scratch/brain/results/bb_j_sig_cor.csv", row.names = F)
print("Done")