# Load Libraries
library("pSI")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("ggplot2")
library("dplyr")
library("tidyverse")
# Read in Data
combined <- readRDS("C:/Users/miles/Downloads/brain/brain_scripts/brain_mz_shiny/data/B1C1C2MZ_combined_031020.rds")
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
rna_path <- "C:/Users/miles/Downloads/brain/"
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
valid_genes <- markers$gene

# Plot number of DEGs per cluster and DEGs in our gene list per cluster
# is pSI necessary if this does the same?
degs = read.csv("C:/Users/miles/Downloads/brain/results/all_markers_B1C1C2MZ_042820.csv")
degs_our = degs[which(degs$gene %in% valid_genes & degs$avg_logFC > 0.5),]
ggplot(degs_our, aes(cluster)) + geom_bar() + labs(title="DEGs in our gene list per cluster")
ggplot(degs,     aes(cluster)) + geom_bar() + labs(title="DEGs per cluster")

# Make a dataframe with genes as rows and clusters as columns
# cluster_df <- data.frame()
# num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
# cells_per_cluster <- c()
# for (i in 0:num_clusters) {
#   this_cells <- WhichCells(combined, idents = i)
#   cluster_df <- rbind(cluster_df, rowSums(combined@assays$RNA@counts[,this_cells]))
#   cells_per_cluster <- c(cells_per_cluster, length(this_cells))
# }
# cluster_df <- t(cluster_df)
# colnames(cluster_df) <- 0:num_clusters
# rownames(cluster_df) <- rownames(combined@assays$RNA@counts)
# cluster_df_norm <- cluster_df / cells_per_cluster
# write.table(cluster_df_norm, paste(rna_path, "/data/psi_input.tsv", sep=""), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
cluster_df_norm <- AverageExpression(combined, slot = "counts", assay = "RNA")[[1]]
cluster_df_norm = myAverageExpression(combined)

# Find the Specificity Index
# results <- specificity.index(cluster_df, e_min = 0.01)
# colnames(results) <- 0:num_clusters
results_norm <- specificity.index(cluster_df_norm, e_min = 1/38)
colnames(results_norm) <- 0:num_clusters

# pSI actually uses BH correction
# fisher_results <- fisher.iteration(results, valid_genes)
fisher_results_norm <- fisher.iteration(results_norm, valid_genes, p.adjust = FALSE)
fisher_results_norm$q <- p.adjust(fisher_results_norm$`0.05 - nominal`)
fisher_results_norm$q.sig <- fisher_results_norm$q < 0.05

# Find the genes in each cluster that are below a p-value threshold
pSI_list <- data.frame()
for (i in 0:num_clusters) {
  this_genes <- rownames(results_norm)[which(results_norm[,i+1] < 0.05)]
  pSI_list <- rbind(pSI_list, data.frame(this_genes, cluster <- rep(i, length(this_genes))))
}
colnames(pSI_list) = c("gene", "cluster")
pSI_list_our = pSI_list[which(pSI_list$gene %in% valid_genes),]
ggplot(pSI_list_our, aes(cluster)) + geom_bar()
write.table(pSI_list, paste0(rna_path, "/results/pSI_genes_ROCK_SAND.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Plot the number of pSI specific genes per cluster
my_list = my.pSI.list(results_norm)
test = my_list$pSi_0.05
df = sapply(1:ncol(test), function(x) length(test[which(! is.na(test[,x])),x]))
df = data.frame(cluster = 0:40, count = df)
ggplot(df, aes(cluster, y=count )) + geom_bar(stat="identity")

# Plot the number of pSI specific genes in our list per cluster
df = sapply(1:ncol(test), function(x) length(test[which(test[,x] %in% valid_genes),x]))
df = data.frame(cluster = 0:40, count = df)
ggplot(df, aes(cluster, y=count )) + geom_bar(stat="identity")

# Find the Specificity per Cluster
sum_clust_norm <- colSums(results_norm[valid_genes,] > 0, na.rm = TRUE)
sum_psi_norm <- colSums(results_norm[valid_genes,], na.rm = TRUE)
sum_clust <- colSums(results[valid_genes,] > 0, na.rm = TRUE)
sum_psi <- colSums(results[valid_genes,], na.rm = TRUE)

sum_clust_norm[order(sum_clust_norm, decreasing = TRUE)]
test <- sum_clust/cells_per_cluster
test[order(test, decreasing=TRUE)]


## Testing Bias ##
Idents(combined) <- "seurat_clusters"
mat <- downsample(combined, rownames(combined@assays$RNA@counts), run = 1)
cluster_0_cells <- WhichCells(combined, idents = 0)
df <- as.data.frame(rowSums(mat[,cluster_0_cells]))
df$gene <- rownames(df)
colnames(df)[1] <- "all_trans"
df$all_trans_per_cell <- df$all_trans/length(cluster_0_cells)

mat2 <- downsample(combined, rownames(combined@assays$RNA@counts), run = 2)
ran_12 <- sample(cluster_0_cells, 12)
df <- cbind(df, rowSums(mat2[,ran_12]))
colnames(df)[4] <- "ran_12_trans"
df$ran_12_trans_per_cell <- df$ran_12_trans/12
df2 <- df %>% pivot_longer(c("all_trans_per_cell", "ran_12_trans_per_cell"), names_to = "variable", values_to = "value")
ggplot(df2, aes(value, fill = variable, color = variable)) + geom_histogram(position = "identity", alpha = 0.5) + ylab("Number of Genes") + xlab("Transcripts_per_cell")
ggplot(df2, aes(value, fill = variable, color = variable)) + geom_histogram(position = "identity", alpha = 0.5) + ylab("Number of Genes") + xlab("Transcripts_per_cell") + xlim(-0.05,0.25) + ggtitle("Zoom")
ggplot(df2[which(df2$value > 0),], aes(value, fill = variable, color = variable)) + geom_histogram(position = "identity", alpha = 0.5) + ylab("Number of Genes") + xlab("Transcripts_per_cell") + xlim(0,0.25) + ggtitle("Non-zero Zoom")

# Test how pSI works
test <- cluster_df_norm[,1:3]
test2 <- specificity.index(cluster_df_norm, e_min = 1/12)
fc12 <- test[,1]/test[,2]
fc23 <- test[,2]/test[,3]
fc13 <- test[,1]/test[,3]
gene_rank_12 <- rownames(test)[order(fc12, decreasing = TRUE)]
gene_rank_23 <- rownames(test)[order(fc23, decreasing = TRUE)]
gene_rank_13 <- rownames(test)[order(fc13, decreasing = TRUE)]
(which(gene_rank_13 == "nup93") + which(gene_rank_23 == "nup93"))/(3-1)

overall <- c()
overall_q <- c()
isSmall <- c()
for (i in 1:1000) {
  ran_list <- sample(rownames(combined), clust_17_deg_size)
  names(data) <- ran_list
  res <- fgsea(pathways, data, nperm = 1000)
  overall <- c(overall, res$pval)
  other_p <- fisher_results_norm_1_12$`0.05 - nominal`
  other_p[18] <- res$pval
  other_q <- p.adjust(other_p)
  overall_q <- c(overall_q, other_q[18])
  isSmall <- c(isSmall, res$pval < 0.004601227)
}
barplot(isSmall)
hist(overall, breaks = seq(0,1, by = 0.05), col = "lightgrey", main = "pSI pvalue for 1000 random lists of same size as cluster 17 DEGs")
hist(overall_q, breaks = seq(0,1, by = 0.05), col = "cadetblue", main = "pSI qvalue for 1000 random lists of same size as cluster 17 DEGs")

all_res <- list()
for (bio in unique(backup_markers$bio)) {
  valid_genes <- backup_markers$gene[which(backup_markers$bio == bio)]
  this_res <- fisher.iteration(results_norm, valid_genes, p.adjust = FALSE)
  this_res$q <- p.adjust(this_res$`0.05 - nominal`)
  this_res$q.sig <- this_res$q < 0.05
  all_res[[bio]] <- rownames(this_res)[which(this_res$q.sig == TRUE)]
}

my.pSI.list <- function(pSIs,write.csv=TRUE){
  
  #Reassign gene names as uppercase for continuity in downstream analysis
  
  #Assigning variables used in loop, ie number of genes, number of samples etc.
  gene_l <- length(pSIs[,2])
  n_samples <- length(pSIs[1,])
  
  
  for(b in c(0.05, 0.01, 0.001, 0.0001)){
    
    for(i in 1:n_samples){
      #Remove genes with pSI values of NA, ie ignore genes whose SI values did not fall within the top 10% of all genes
      sig_tmp <- pSIs[!is.na(pSIs[,i]),]
      #Keep genes with pSI values below threshold, ie ignore genes whose pSI values were above set pSI threshold value
      sig_tmp <- sig_tmp[sig_tmp[,i] < b,]
      sig_tmp <- data.frame(rownames(sig_tmp), stringsAsFactors=FALSE)
      
      #Binds together the lists of genes specific to each sample type for a given threshold
      if(!i==1){
        gene_bind <- cbindX(gene_bind, sig_tmp)
        
      }     
      else{
        gene_bind <- sig_tmp
      }  
    }
    #Assigns columns names to temporary loop data frame as well assign temporary loop data frame a descriptive name
    #colnames(gene_bind) <- paste(colnames(pSIs), b, sep="_")
    if(b<0.001){
      colnames(gene_bind) <- paste(colnames(pSIs), "0.0001", sep="_")
    }else{
      colnames(gene_bind) <- paste(colnames(pSIs), b, sep="_")
    }
    if(b<0.001){
      assign(paste("pSi", "0.0001", sep="_"),gene_bind)
    }else{
      assign(paste("pSi", b, sep="_"),gene_bind)
    }
    #assign(paste("pSi", b, sep="_"),gene_bind)
    
    #Boolean used to indicate if CSV files of the output should be writted to the current working directory
    if(write.csv==TRUE){
      write.table(t(gene_bind), paste("pSi", paste(b,".csv",sep=""), sep="_"), col.names=FALSE, sep=",", na="")
    }
    
    if(b<0.05){
      out <- c(out,gene_bind) 
    }else{
      out <-gene_bind
    }  
    
  }
  
  #create list of 6 dataframes, one for each pSI threshold
  pSi_list <- list(get("pSi_0.0001"), get("pSi_0.001"), get("pSi_0.01"), get("pSi_0.05"))
  pSI_rang <- c("pSi_0.0001", "pSi_0.001", "pSi_0.01", "pSi_0.05")
  names(pSi_list)<- pSI_rang
  
  return(pSi_list)
  
}