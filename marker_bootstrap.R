library("stringr")
library("dplyr")
library("Seurat")

rna_path <- "C:/Users/miles/Downloads/brain/"
# rna_path <- "~/scratch/brain/"
combined <- readRDS(paste(rna_path, "/brain_scripts/brain_shiny/data/combined.rds", sep = ""))
marker_path <- paste(rna_path, "data/markers/", sep="")
marker_files <- dir(marker_path, pattern =paste("*.txt", sep=""))


geneCap <- function(gene, gene_names) {
  # Gene the gene name in the right format
  gene_lower <- tolower(gene)
  gene_upper <- toupper(gene)
  gene_title <- str_to_title(gene)
  error <- FALSE
  if (gene_lower %in% gene_names) {
    gene <- gene_lower
  } else if (gene_upper %in% gene_names) {
    gene <- gene_upper
  } else if (gene_title %in% gene_names) {
    gene <- gene_title
  } else {
    error <- TRUE
  }
  
  return(c(gene, error))
}

validGenes <- function(genes, gene_names) {
  valid_genes <- c()
  for (gene in genes) {
    result <- geneCap(gene, gene_names)
    gene <- result[1]
    error <- as.logical(result[2])
    if (! error) {
      valid_genes <- c(valid_genes, gene)
    }
  } # end gene for
  unique(valid_genes)
  return(valid_genes)
} # end validGenes function

markers <- data.frame(gene <- c(), bio <- c())
for (i in 1:length(marker_files)) {
  file <- read.table(paste(marker_path, marker_files[i], sep=""), header=FALSE, sep="\t", stringsAsFactors=FALSE)
  file[,1] <- toupper(file[,1])
  markers <- rbind(markers, file[,1:2])
}
colnames(markers) <- c("gene", "bio")
markers <- markers[which(markers$bio == "DISC_ASE"),]
markers <- markers[which(! startsWith(markers$gene, "LOC")),]
gene_names <- rownames(combined@assays$RNA)
valid_genes <- validGenes(markers$gene, gene_names)

# Avg genes in bootstrap
big_df <- data.frame()
for (j in 1:50) {
  print(j)
  set.seed(j)
  ran_markers <- sample(valid_genes, 50)
  Idents(object = combined) <- "seurat_clusters"
  num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
  gene_names <- rownames(combined@assays$RNA)
  cells_per_cluster <- c()
  genes_per_cluster <- c()
  for (i in 0:num_clusters) {
    this_cells <- WhichCells(combined, idents = i)
    genes_per_cluster <- c(genes_per_cluster, length(which(as.vector(combined@assays$RNA@counts[valid_genes,this_cells]) != 0))) # genes
    cells_per_cluster <- c(cells_per_cluster, length(this_cells))
  }
  avg_gene_per_cell_per_cluster <- genes_per_cluster/cells_per_cluster
  df <- as.data.frame(avg_gene_per_cell_per_cluster)
  df$genes_per_cluster <- genes_per_cluster
  df$cells_per_cluster <- cells_per_cluster
  rownames(df) <- 0:40
  # df <- df[order(df$avg_gene_per_cell_per_cluster, decreasing = TRUE),]
  big_df <- rbind(big_df, t(order(df$avg_gene_per_cell_per_cluster, decreasing = TRUE)))
}
rownames(big_df) <- rep(paste("run_", 1:50, sep =""))
colnames(big_df) <- rep(paste("rank", 1:41, sep =""))

rankSums <- c()
for (cluster in 1:41) {
  rankSums <- c(rankSums, sum(which(big_df[,] == cluster, arr.ind = TRUE)[,2]))
}
order(rankSums)
