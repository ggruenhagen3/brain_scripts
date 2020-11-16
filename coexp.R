# BiocManager::install("qvalue")
library("Seurat")
library("Matrix")
library("qvalue")
library("jaccard")
# iegs <- read.csv("C:/Users/miles/Downloads/zack_IEG_list_061720.csv", header = FALSE, stringsAsFactors = F)
# iegs <- iegs$V1
bb <- readRDS("~/scratch/brain/data/bb.RDS")
obj <- bb
# combined <- readRDS("C:/Users/miles/Downloads/brain/brain_scripts/brain_mz_shiny/data/B1C1C2MZ_combined_031020.rds")
# obj <- combined
gene_names <- rownames(obj)[which(rowSums(as.matrix(obj@assays$RNA@counts)) != 0)]

# First try
# mat  <- matrix(0, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
# 
# ptm <- proc.time()
# for (col in 1:ncol(obj)) {
# # for (col in 1:100) {
#   if (col %% 100 == 0) {
#     print(col)
#   }
#   dat = obj@assays$RNA@counts[,col]
#   non_zero_genes = names(dat[which(dat > 0)])
#   tmp <- mat[non_zero_genes, non_zero_genes]
#   mat[non_zero_genes, non_zero_genes] = tmp + 1
# }
# proc.time() - ptm
# saveRDS(mat, "C:/Users/miles/Downloads/brain/data/coexp_mat.RDS")
# write.csv(mat, "C:/Users/miles/Downloads/brain/data/coexp_mat.csv", quote = FALSE)

# mat2  <- matrix(0, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
# ptm <- proc.time()
# obj_mat = obj@assays$RNA@counts
# for (col in 2:length(gene_names)) {
#   # for (col in 1:100) {
#   if (col %% 100 == 0) {
#     print(col)
#   }
#   gene1 <- gene_names[col]
#   for ( row in 1:(col-1) ) {
#     gene2 <- gene_names[row]
#     tmp = colSums(obj_mat[c(gene1, gene2),])
#     mat2[row, col] = length(tmp[which(tmp > 0)])
#   }
# }
# gene_cell = list()
# for (gene in gene_names) {
#   gene_cell[[gene]] = colnames(obj)[which(obj@assays$RNA@counts[gene,] != 0)]
# }
# for (col in 2:length(gene_names)) {
#   # for (col in 1:100) {
#   if (col %% 100 == 0) {
#     print(col)
#   }
#   gene1 <- gene_names[col]
#   gene1_cells = gene_cell[[gene1]]
#   for ( row in 1:(col-1) ) {
#     gene2 <- gene_names[row]
#     gene2_cells = gene_cell[[gene2]]
#     mat2[row, col] = length(unique(c(gene1_cells, gene2_cells)))
#   }
# }
# proc.time() - ptm
# mat3 = mat/mat2

# Normalization bc some genes are more common than others
# for (col in 1:length(gene_names)) {
#   if (col %% 1000 == 0) {
#     print(col)
#   }
#   col_vector <- mat[,col]
#   mat[col,] = col_vector/col_vector[col]
# }

# Keep only the upper triangle
# mat[lower.tri(mat,diag=TRUE)] <- 0

# maxs <- data.frame()
# for (ieg in iegs) {
#   ieg_row <- mat3[ieg,]
#   my_max <- names(which.max(ieg_row[which(names(ieg_row) != ieg)]))
#   newRow <- data.frame(ieg, my_max, ieg_row[my_max], my_max %in% iegs)
#   maxs <- rbind(maxs, newRow)
# }
# colnames(maxs) <- c("ieg", "ieg_max_gene", "ieg_max_num", "isIEG")
# ggplot(maxs, aes())

######
# Bi #
######
# mat_bi  = matrix(0, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
# mat_j   = matrix(0, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
# gene_bi = lapply(gene_names, function(x) rep(0, ncol(obj)))
# names(gene_bi) = gene_names
# 
# for (col in 1:ncol(obj)) {
#   if (col %% 100 == 0) {
#     print(col)
#   }
#   dat = obj@assays$RNA@counts[,col]
#   non_zero_genes = names(dat[which(dat > 0)])
#   for (gene in non_zero_genes) {
#     gene_bi[[gene]][col] = 1
#   }
#   # gene_bi[[non_zero_genes]][col] = 1
# }
# for (col in 2:length(gene_names)) {
#   if (col %% 100 == 0) {
#     print(col)
#   }
#   gene1 <- gene_names[col]
#   for ( row in 1:(col-1) ) {
#     gene2 <- gene_names[row]
#     mat_bi[row, col] = jaccard.test(gene_bi[[gene1]], gene_bi[[gene2]], method = "mca", accuracy=1e-5)$pvalue
#     mat_j[row, col]  = jaccard(gene_bi[[gene1]], gene_bi[[gene2]])
#   }
# }
# saveRDS(mat_bi,   "/nv/hp10/ggruenhagen3/scratch/brain/data/mat_bi.RDS")
# saveRDS(mat_j,   "/nv/hp10/ggruenhagen3/scratch/brain/data/mat_j.RDS")

# ###############
# # Transcripts #
# ###############
setToMax = function(dat, x) {
  dat[which(dat > x)] = x
  return(dat)
}
mat_trans  <- matrix(0L, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
mat2_trans <- matrix(0L, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))

for (col in 1:ncol(obj)) {
  # for (col in 1:100) {
  if (col %% 100 == 0) {
    print(col)
  }
  dat = obj@assays$RNA@counts[,col]
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

saveRDS(mat3_trans, "~/scratch/brain/data/bb_j.RDS")
saveRDS(mat3_trans_p, "~/scratch/brain/data/bb_j_p.RDS")
# egr1_top = sort(mat3_trans_p["egr1",])
# 
# non_zero_lncRNA = lncRNA_genes[which(lncRNA_genes %in% gene_names)]
# maxs <- data.frame()
# gene_list = iegs
# for (gene in gene_list) {
#   row <- mat3_trans[gene,]
#   my_max <- names(which.max(row[which(names(row) != gene)]))
#   newRow <- data.frame(gene, my_max, row[my_max], my_max %in% gene_list)
#   maxs <- rbind(maxs, newRow)
# }
# colnames(maxs) <- c("gene", "coexp_max_gene", "coexp_max_num", "inList")
# 
# all_maxs = data.frame()
# for (gene in gene_names) {
#   row <- mat3_trans[gene,]
#   my_max <- names(which.max(row[which(names(row) != gene)]))
#   newRow <- data.frame(gene, my_max, row[my_max])
#   all_maxs <- rbind(all_maxs, newRow)
# }
# colnames(all_maxs) <- c("gene", "coexp_max_gene", "coexp_max_num")
# isReciprocal = c()
# for (gene in gene_names) {
#   df_row = all_maxs[which(all_maxs$gene == gene),]
#   max_gene = as.vector(df_row$coexp_max_gene)
#   isReciprocal = c(isReciprocal, all_maxs$coexp_max_gene[which(all_maxs$gene == max_gene)] == gene)
# }
# all_maxs$isReciprocal = isReciprocal
# colnames(all_maxs)[ncol(all_maxs)] = "isReci"


# ptm <- proc.time()
# for (col in 2:length(gene_names)) {
#   # for (col in 1:100) {
#   if (col %% 100 == 0) {
#     print(col)
#   }
#   gene1 <- gene_names[col]
#   gene1_cells = gene_cell[[gene1]]
#   for ( row in 1:(col-1) ) {
#     gene2 <- gene_names[row]
#     gene2_cells = gene_cell[[gene2]]
#     int_cells = gene1_cells[which(gene1_cells %in% gene2_cells)]
#     mat[row,col]         = sum(rowSums(obj@assays$RNA@counts[c(gene1, gene2), int_cells]))
#     mat2_trans[row, col] = sum(rowSums(obj@assays$RNA@counts[c(gene1, gene2), ]))
#   }
# }
# proc.time() - ptm
# mat3_trans = mat_trans/mat2_trans

#
# mat3_trans = readRDS("C:/Users/miles/Downloads/brain/data/int_uni_lncRNA_mat.RDS")
# mat3_trans_p = matrix(jaccard.rahman(as.vector(mat3_trans)), nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
# egr1_val = matr3_trans

#########################
# Spearman Correlations #
#########################
# mat_data_p   = matrix(0, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
# mat_data_cor = matrix(0, nrow=length(gene_names), ncol = length(gene_names), dimnames = list(gene_names, gene_names))
# for (col in 2:length(gene_names)) {
#     # if (col %% 100 == 0) {
#       print(col)
#     # }
#   gene1 = gene_names[1]
#   for ( row in 1:(col-1) ) {
#     gene2 = gene_names[2]
#     cor_res = cor.test(obj@assays$RNA@data[gene1,], obj@assays$RNA@data[gene2,])
#     mat_data_p[row,col]   = cor_res$p.value
#     mat_data_cor[row,col] = cor_res$estimate
#   }
# }

#################
# Visualization #
#################
# Local
# test <- mat[iegs$V1[1:200], iegs$V1[1:200]]
# test[lower.tri(test,diag=TRUE)] <- 0
# heatmap(test, Rowv=NA, Colv=NA, revC=TRUE, scale="none")
# 
# library("igraph")
# ig_obj = graph_from_incidence_matrix(test[1:3,1:3], weighted = TRUE)
# plot(ig_obj, vertex.label="")
# deg <- degree(ig_obj, mode="all")
# V(ig_obj)$size <- deg*3
# plot(ig_obj, vertex.label="")

############
# Analysis #
############
# mat_fish = readRDS("C:/Users/miles/Downloads/brain/data/mat_fish.RDS")
# 
# mat = mat_fish
# sig_df = data.frame()
# for (row in 1:length(gene_names)) {
#   if (row %% 1000 == 0) { print(row) }
#   q = p.adjust(mat[row, (row+1):ncol(mat)], method = "bonferroni")
#   q_sig = q[which(q < 0.05)]
#   if (length(q_sig) > 0) {
#     newRow = data.frame(gene = rep(gene_names[row], length(q_sig)), q = q_sig, close_gene = names(q_sig))
#     sig_df = rbind(sig_df, newRow)
#   }
# }