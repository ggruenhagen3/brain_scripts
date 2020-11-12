library("dplyr")
library("Matrix")
library("Seurat")
library("stringr")
library("ggplot2")

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
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
down_avg_avg_gene <- rep(0, num_clusters+1)
run_num <- 50

# No Perm, Bootstrap
down_list <- lapply(0:num_clusters, function(x) c())
for (run in 1:run_num) {
  cat(paste("no_perm", run, "\n"))
  mat <- downsample(combined, marker_genes, run)
  
  cells_per_cluster <- c()
  genes_per_cluster <- c()
  for (i in 0:num_clusters) {
    this_cells <- WhichCells(combined, idents = i)
    this_genes <- length(which(as.vector(mat[,this_cells]) != 0))
    # genes_per_cluster <- c(genes_per_cluster, length(which(as.vector(combined@assays$RNA@counts[ran_markers,this_cells]) != 0))) # genes
    genes_per_cluster <- c(genes_per_cluster, this_genes) # genes
    cells_per_cluster <- c(cells_per_cluster, length(this_cells))
    down_list[[i+1]] <- c(down_list[[i+1]], this_genes/ length(this_cells))
  }
  avg_gene_per_cell_per_cluster <- genes_per_cluster/cells_per_cluster
  down_avg_avg_gene <- down_avg_avg_gene + avg_gene_per_cell_per_cluster
}
down_avg_avg_gene <- down_avg_avg_gene / run_num
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
    perm_down_avg_gene[[i+1]] <- c(perm_down_avg_gene[[i+1]], this_genes/length(this_cells))
  }
  avg_gene_per_cell_per_cluster <- genes_per_cluster/cells_per_cluster
  # perm_down_avg_gene <- c(perm_down_avg_gene, avg_gene_per_cell_per_cluster)
}

# Compare empirical data to 97.5th percentile of the permutated data on a PER CLUSTER basis
sig_clusters  <- c()
alpha <- 0.01
bonferroni_alpha <- (alpha/num_clusters)/length(gene_names)
# upper_tail <- bonferroni_alpha/2
upper_tail <- alpha/2
for (i in 0:num_clusters) {
  sig <- quantile(perm_down_avg_gene[[i+1]], c(1-upper_tail))
  print(paste("Max value for perm cluster ", i, ":", max(perm_down_avg_gene[[i+1]])))
  print(paste("Mean Real value for cluster", i, ":", down_avg_avg_gene[i+1]))
  print(paste("Min Real Value for cluster", i, ":", min(down_list[[i+1]])))
  bot_25 <- quantile(down_list[[i+1]], c(0.25))
  print(paste("Bot 25 Value for cluster", i, ":", bot_25))
  if ( min(down_list[[i+1]]) > max(perm_down_avg_gene[[i+1]]) ) {
    sig_clusters <- c(sig_clusters, i)
  }
}
# sig <- quantile(perm_down_avg_gene, c(.975))
# print(sig)
# sig_clusters <- which(down_avg_avg_gene > sig)-1

df <- t(as.data.frame(lapply(1:length(down_list), function(x) t(cbind(down_list[[x]], rep(x, length(num_clusters)+1), rep("Down", length(num_clusters)+1))))))
df <- rbind(df, t(as.data.frame(lapply(1:length(perm_down_avg_gene), function(x) t(cbind(perm_down_avg_gene[[x]], rep(x, length(num_clusters)+1), rep("Down+Perm", length(num_clusters)+1)))))))
colnames(df) <- c("avg_gene_per_cell_per_cluster", "cluster", "cond")
df <- as.data.frame(df)
df$avg_gene_per_cell_per_cluster <- as.numeric(as.vector(df$avg_gene_per_cell_per_cluster))
df$cluster <- factor(df$cluster, levels = 1:(num_clusters+1))
png(paste0(rna_path, "/results/down_and_perm_2_", bio, ".png"), width = 1800, height = 1000, res = 150)
p <- ggplot(df, aes(x = cluster, y = avg_gene_per_cell_per_cluster, fill = cond)) + geom_boxplot(alpha = 0.6) + geom_jitter(shape=16, position=position_jitterdodge(), alpha = 0.3, aes(colour = cond)) + scale_colour_manual(values=c("#999999", "#56B4E9", "#3ac9bb")) + scale_fill_manual(values=c("#999999", "#3ac9bb", "#56B4E9")) + ggtitle(paste(bio, "- Average Genes per Cell per Cluster"))
print(p)
dev.off()

print(sig_clusters)
write.csv(sig_clusters, file = paste(rna_path, "/results/down_and_perm_2_", bio, ".csv", sep=""), row.names = FALSE)
