library("dplyr")
library("Matrix")
library("Seurat")
library("stringr")

# rna_path <- "C:/Users/miles/Downloads/brain/"
rna_path <- "~/scratch/brain/"
source(paste0(rna_path, "/brain_scripts/all_f.R"))

# combined <- readRDS(paste(rna_path, "/brain_scripts/brain_shiny/data/combined.rds", sep = ""))
combined <- readRDS(paste(rna_path, "/data/B1C1C2MZ_combined_031020.rds", sep=""))
marker_path <- paste(rna_path, "data/markers/", sep="")
marker_files <- dir(marker_path, pattern =paste("*.txt", sep=""))

markers <- data.frame(gene <- c(), bio <- c())
for (i in 1:length(marker_files)) {
  file <- read.table(paste(marker_path, marker_files[i], sep=""), header=FALSE, sep="\t", stringsAsFactors=FALSE)
  markers <- rbind(markers, file[,1:2])
}
colnames(markers) <- c("gene", "bio")
bio <- "ROCK_SAND"
markers <- markers[which(markers$bio == bio),]
gene_names <- rownames(combined@assays$RNA)
marker_genes <- markers$gene
valid_genes <- marker_genes
other_genes = rownames(combined)[which(! rownames(combined) %in% valid_genes)]
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
total_genes_per_cluster <- rep(0, num_clusters+1)
run_num <- 1000
test_clusters = c(1,5,12)

mat = combined@assays$RNA@counts
perm_df = data.frame()

# Real Data
run = 1
Idents(combined) = combined$seurat_clusters
for (i in 0:num_clusters) {
  this_cells <- WhichCells(combined, idents = i)
  this_pos <- sum(rowSums(mat[valid_genes,this_cells]))
  this_neg <- sum(rowSums(mat[other_genes,this_cells]))
  
  test_ratio = sum(colSums(mat[valid_genes,this_cells]))/sum(colSums(mat[,this_cells]))
  
  other_cells = colnames(combined)[which(! colnames(combined) %in% this_cells)]
  other_pos = sum(rowSums(mat[valid_genes, other_cells]))
  other_neg = sum(rowSums(mat[other_genes, other_cells]))
  
  this_pos_pct = (this_pos)/(this_pos + this_neg)
  this_neg_pct = (this_neg)/(this_pos + this_neg)
  
  other_pos_pct = (other_pos)/(other_pos + other_neg)
  other_neg_pct = (other_neg)/(other_pos + other_neg)
  
  contig_table = matrix(c(this_pos, other_pos, this_neg, other_neg), nrow=2)
  # contig_table = matrix(c(log2(this_pos), log2(other_pos), log2(this_neg), log2(other_neg)), nrow=2)
  p = fisher.test(contig_table, alternative = "greater")[[1]]
  
  newRow = data.frame("Real", run, i, this_pos/this_neg, other_pos/other_neg, (this_pos/this_neg)/(other_pos/other_neg), test_ratio)
  colnames(newRow) = c("isReal", "Run", "Cluster", "this_ratio", "other_ratio", "this_to_other", "test_ratio")
  perm_df = rbind(perm_df, newRow)
}
# perm_df$q = p.adjust(perm_df$p, method="bonferroni")

# Permutation
backup_ids <- combined@meta.data$seurat_clusters
for (run in 1:run_num) {
  cat(paste0(run, "."))
  # new_cells <- shuffleClusters(combined)
  set.seed(run)
  combined$shuffled = sample(as.numeric(as.vector(combined$seurat_clusters)))
  Idents(combined) = combined$shuffled
  for (i in 0:num_clusters) {
    this_cells = WhichCells(combined, idents=i)
    # this_cells <- new_cells[[i+1]]
    this_pos <- sum(rowSums(mat[valid_genes,this_cells]))
    this_neg <- sum(rowSums(mat[other_genes,this_cells]))
    
    test_ratio = sum(colSums(mat[valid_genes,this_cells]))/sum(colSums(mat[,this_cells]))
    
    other_cells = colnames(combined)[which(! colnames(combined) %in% this_cells)]
    other_pos = sum(rowSums(mat[valid_genes, other_cells]))
    other_neg = sum(rowSums(mat[other_genes, other_cells]))
    
    # this_pos_pct = (this_pos)/(this_pos + this_neg)
    # this_neg_pct = (this_neg)/(this_pos + this_neg)
    # 
    # other_pos_pct = (other_pos)/(other_pos + other_neg)
    # other_neg_pct = (other_neg)/(other_pos + other_neg)
    # 
    # contig_table = matrix(c(this_pos, other_pos, this_neg, other_neg), nrow=2)
    # p = fisher.test(contig_table, alternative = "greater")[[1]]
    
    newRow = data.frame("Boot", run, i, this_pos/this_neg, other_pos/other_neg, (this_pos/this_neg)/(other_pos/other_neg), test_ratio)
    colnames(newRow) = c("isReal", "Run", "Cluster", "this_ratio", "other_ratio", "this_to_other", "test_ratio")
    perm_df = rbind(perm_df, newRow)
  }
}
cat("\n")
write.table(perm_df, "~/scratch/brain/results/perm_raw.tsv", sep="\t", quote = F)

# res_df = data.frame()
# for (i in test_clusters) {
#   this_df = perm_df[which(perm_df$isReal == "Boot" & perm_df$Cluster == i),]
#   real_value = perm_df$this_ratio[which(perm_df$isReal == "Real" & perm_df$Cluster == i)]
# 
#   this_df$above = this_df$this_ratio >= real_value
#   print(ggplot(this_df, aes(this_ratio, alpha=.7, fill=above)) + geom_histogram(alpha=0.5, color = "purple") + geom_vline(aes(xintercept = real_value)) + geom_text(aes(x=real_value, label="Real Value"), y = Inf, hjust=0, vjust=1, color = "black") + guides(color=F, alpha=F, fill=F) + ggtitle(paste("Cluster", i)))
#   res_df = rbind(res_df, t(c(i, length(which(this_df$above)), length(which(this_df$above))/run_num )))
# }

# Idents(combined) = combined$seurat_clusters
# markerExpPerCellPerCluster(combined, valid_genes, n_markers = F, correct = T)
# colnames()