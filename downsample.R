library("dplyr")
library("Matrix")
library("Seurat")
library("stringr")

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
  
  return(valid_genes)
} # end validGenes function

downsample <- function(combined, marker_genes, run) {
  set.seed(run)
  min_trans <- min(combined$nCount_RNA)
  gene_names <- rownames(combined@assays$RNA@counts)
  # new_matrix <- matrix(, nrow = nrow(combined@assays$RNA@counts), ncol = ncol(combined@assays$RNA@counts), dimnames = list(gene_names, colnames(combined@assays$RNA@counts)))
  # new_new_matrix <- matrix(, nrow=nrow(combined@assays$RNA@counts))
  marker_matrix  <- matrix(, nrow=length(marker_genes), ncol = ncol(combined@assays$RNA@counts), dimnames = list(marker_genes, colnames(combined@assays$RNA@counts)))
  i <- 0
  for (cell in colnames(combined@assays$RNA@counts)) {
    # i <- i + 1
    # if (i%%500 == 1) {
    #   print(cell)
    # }
    # start.time <- Sys.time()
    
    trans_names <- rep(gene_names, combined@assays$RNA@counts[,cell])
    ran_trans_names <- sample(trans_names, min_trans)
    ran_trans_names <- ran_trans_names[which(ran_trans_names %in% marker_genes)]
    ran_df <- as.data.frame(table(ran_trans_names))
    zero_gene_names <- marker_genes[which(! marker_genes %in% ran_trans_names)]
    zero_df <- setNames(data.frame(zero_gene_names <- zero_gene_names, Freq <- rep(0, length(zero_gene_names))), c("ran_trans_names", "Freq"))
    ran_df <- rbind(ran_df, zero_df)
    rownames(ran_df) <- ran_df$ran_trans_names
    ran_df <- ran_df[marker_genes,2]
    # new_matrix[,cell] <-as.matrix(ran_df)
 
    
    marker_matrix[,cell] <- as.matrix(ran_df)
    # new_new_matrix <- cbind(new_new_matrix, as.matrix(ran_df))
    
    # end.time <- Sys.time()
    # time.taken <- end.time - start.time
    # print(time.taken)
  }
  
  return(marker_matrix)
}
## END FUNCTIONS ##
rna_path <- "~/scratch/brain/"
combined <- readRDS(paste(rna_path, "/brain_scripts/brain_shiny/data/combined.rds", sep = ""))
marker_path <- paste(rna_path, "data/markers/", sep="")
marker_files <- dir(marker_path, pattern =paste("*.txt", sep=""))

markers <- data.frame(gene <- c(), bio <- c())
for (i in 1:length(marker_files)) {
  file <- read.table(paste(marker_path, marker_files[i], sep=""), header=FALSE, sep="\t", stringsAsFactors=FALSE)
  file[,1] <- toupper(file[,1])
  markers <- rbind(markers, file[,1:2])
}
colnames(markers) <- c("gene", "bio")
markers <- markers[which(markers$bio == "DISC_ASE"),]
gene_names <- rownames(combined@assays$RNA)
marker_genes <- validGenes(markers$gene, gene_names)
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
big_df <- data.frame()

# Avg trans per cluster for ALL genes
Idents(object = combined) <- "seurat_clusters"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
results_expr <- data.frame()
cells_per_cluster <- c()
trans_per_cluster <- c()
for (i in 0:num_clusters) {
  print(i)
  this_cells <- WhichCells(combined, idents = i)
  # genes_per_cluster <- c(genes_per_cluster, length(which(as.vector(combined@assays$RNA@counts[,this_cells]) != 0))) # genes
  trans_per_cluster <- c(trans_per_cluster, sum(colSums(combined@assays$RNA@counts[,this_cells]))) # trans
  cells_per_cluster <- c(cells_per_cluster, length(this_cells))
}
avg_trans_per_cell_per_cluster <- trans_per_cluster/cells_per_cluster
all_trans <- as.data.frame(avg_trans_per_cell_per_cluster)
all_trans$trans_per_cluster <- trans_per_cluster
all_trans$cells_per_cluster <- cells_per_cluster
rownames(all_trans) <- 0:40


# Bootstrap
for (run in 1:50) {
  print(run)
  mat <- downsample(combined, marker_genes, run)
  
  cells_per_cluster <- c()
  genes_per_cluster <- c()
  for (i in 0:num_clusters) {
    this_cells <- WhichCells(combined, idents = i)
    # genes_per_cluster <- c(genes_per_cluster, length(which(as.vector(combined@assays$RNA@counts[ran_markers,this_cells]) != 0))) # genes
    genes_per_cluster <- c(genes_per_cluster, length(which(as.vector(mat[,this_cells]) != 0))) # genes
    cells_per_cluster <- c(cells_per_cluster, length(this_cells))
  }
  avg_gene_per_cell_per_cluster <- genes_per_cluster/cells_per_cluster
  df <- as.data.frame(avg_gene_per_cell_per_cluster / all_trans[,1]) # addition
  rownames(df) <- 0:40
  big_df <- rbind(big_df, t(order(df[,1], decreasing = TRUE)-1))
}

rownames(big_df) <- rep(paste("run_", 1:50, sep =""))
colnames(big_df) <- rep(paste("rank", 1:41, sep =""))

rankSums <- c()
for (cluster in 0:40) {
  rankSums <- c(rankSums, sum(which(big_df[,] == cluster, arr.ind = TRUE)[,2]))
}
order(rankSums)-1
write.csv(order(rankSums)-1, file = paste(rna_path, "/results/down_boot_rank_all_trans.csv", sep=""), row.names = FALSE)
          