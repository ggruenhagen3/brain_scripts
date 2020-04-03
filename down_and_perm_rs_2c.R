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

pickNewCells <- function(combined, num_clusters, num_cells) {
  new_cells <- c()
  for (i in 1:num_cells) {
    ran_cluster <- sample(0:num_clusters, 1)
    this_cells <- names(combined$seurat_clusters[which(combined$seurat_clusters == ran_cluster)])
    new_cells <- c(new_cells, sample(this_cells,1))
  }
  
  return(new_cells)
}

shuffleClusters <- function(combined) {
  # The selection process for a new cluster should be as follows:
  # 1. Pick a random cluster 0-40
  # 2. Pick a random cell from that cluster to be a part of the new cluster
  # This means that the new data set would likely have duplicate cells
  new_cells <- lapply(0:num_clusters, function(x) c())
  num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
  for (i in 0:num_clusters) {
    num_cells <- length(combined$seurat_clusters[which(combined$seurat_clusters == i)])
    new_cells[[i+1]] <- pickNewCells(combined, num_clusters, num_cells)
  }
  
  return(new_cells)
}
## END FUNCTIONS ##
# rna_path <- "C:/Users/miles/Downloads/brain/"
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
bio <- "RAN_2"
markers <- markers[which(markers$bio == bio),]
print("Before gene_names")
gene_names <- rownames(combined@assays$RNA)
print("After gene_names")
marker_genes <- unique(validGenes(markers$gene, gene_names))
valid_genes <- marker_genes
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
down_avg_avg_gene <- rep(0, num_clusters+1)
total_genes_per_cluster <- rep(0, num_clusters+1)
run_num <- 50

# No Perm, Bootstrap
down_avg_gene <- lapply(0:num_clusters, function(x) c())
for (run in 1:run_num) {
  cat(paste("no_perm", run, "\n"))
  mat <- downsample(combined, marker_genes, run)
  
  cells_per_cluster <- c()
  genes_per_cluster <- c()
  for (i in 0:num_clusters) {
    this_cells <- WhichCells(combined, idents = i)
    # genes_per_cluster <- c(genes_per_cluster, length(which(as.vector(combined@assays$RNA@counts[ran_markers,this_cells]) != 0))) # genes
    genes_per_cluster <- c(genes_per_cluster, length(which(as.vector(mat[,this_cells]) != 0))) # genes
    cells_per_cluster <- c(cells_per_cluster, length(this_cells))
    down_avg_gene[[i+1]] <- c(down_avg_gene[[i+1]], length(which(as.vector(mat[,this_cells]) != 0)))
  }
  # avg_gene_per_cell_per_cluster <- genes_per_cluster/cells_per_cluster
  # down_avg_avg_gene <- down_avg_avg_gene + avg_gene_per_cell_per_cluster
  total_genes_per_cluster <- total_genes_per_cluster + genes_per_cluster
}
down_avg_avg_gene <- total_genes_per_cluster / run_num
print(down_avg_avg_gene)

# Perm, Bootstrap
backup_ids <- combined@meta.data$seurat_clusters
perm_down_avg_gene <- lapply(0:num_clusters, function(x) c())
for (run in (run_num+1):(run_num+run_num)) {
  cat(paste("perm", run, "\n"))
  set.seed(run)
  mat <- downsample(combined, marker_genes, run)
  
  new_cells <- shuffleClusters(combined)
  num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
  gene_names <- rownames(combined@assays$RNA)
  cells_per_cluster <- c()
  genes_per_cluster <- c()
  for (i in 0:num_clusters) {
    this_cells <- new_cells[[i+1]]
    this_genes <- length(which(as.vector(mat[valid_genes,this_cells]) != 0))
    genes_per_cluster <- c(genes_per_cluster, this_genes) # genes
    cells_per_cluster <- c(cells_per_cluster, length(this_cells))
    perm_down_avg_gene[[i+1]] <- c(perm_down_avg_gene[[i+1]], this_genes)
  }
  # avg_gene_per_cell_per_cluster <- genes_per_cluster/cells_per_cluster
  # perm_down_avg_gene <- c(perm_down_avg_gene, avg_gene_per_cell_per_cluster)
}

# Compare empirical data to the permutated data on a PER CLUSTER basis
df <- data.frame()
for (i in 0:num_clusters) {
  p <- wilcox.test(down_avg_gene[[i+1]], perm_down_avg_gene[[i+1]])$p.value
  
  df <- rbind(df, t(c(i, p, mean(down_avg_gene[[i+1]]) > mean(perm_down_avg_gene[[i+1]]))) )
}
df$q <- p.adjust(df[,2], method = "hochberg")
df$q_sig_enrich <- df$q < 0.05 & df[,3] == TRUE
df$q_sig_unenrich <- df$q < 0.05 & df[,3] == FALSE
write.table(df, file = paste(rna_path, "/results/down_and_perm_rs_2c_", bio, ".tsv", sep=""), sep = "\t", row.names = FALSE, quote=FALSE)
