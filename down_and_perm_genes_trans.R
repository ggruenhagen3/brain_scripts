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
# Find what genes from a gene list are enriched in each cluster
# Downsample and downsample+permutate. Compare the pos cells for each gene for each cluster
# to the downsampled+permutated data.
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
marker_genes <- unique(validGenes(markers$gene, gene_names))
valid_genes <- marker_genes
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
num_runs <- 50

sort(valid_genes)
gene_df <- data.frame(gene <- rep(valid_genes, num_clusters+1), cluster <- rep(0:num_clusters,each=length(valid_genes)), sum <- rep(0, length(valid_genes)*(num_clusters+1)), p <-  rep(0, length(valid_genes)*(num_clusters+1)), q <-  rep(0, length(valid_genes)*(num_clusters+1)), up <-  rep(FALSE, length(valid_genes)*(num_clusters+1)))
perm_gene_df <- gene_df
colnames(gene_df) <- c("gene","cluster", "sum", "p", "q", "up")
colnames(perm_gene_df) <- c("gene","cluster", "sum", "p", "q", "up")

# No Perm, Bootstrap
for (run in 1:num_runs) {
  cat(paste("no_perm", run, "\n"))
  mat <- downsample(combined, marker_genes, run)
  # mat[mat != 0] <- 1 # allows for rowSums to give the number of positive cells
  
  cells_per_cluster <- c()
  for (i in 0:num_clusters) {
    start_i <- 1 + i*length(valid_genes)
    end_i <- start_i + length(valid_genes)
    this_cells <- WhichCells(combined, idents = i)
    gene_df[start_i:end_i,3] <- gene_df[start_i:end_i,3] + rowSums(mat[,this_cells])
    cells_per_cluster <- c(cells_per_cluster, length(this_cells))
  }
}
# Take the average of the 50 runs
gene_df[,3] <- gene_df[,3] / num_runs

# Perm, Bootstrap
backup_ids <- combined@meta.data$seurat_clusters
for (run in (num_runs+1):(num_runs+num_runs)) {
  cat(paste("perm", run, "\n"))
  set.seed(run)
  shuffled <- sample(backup_ids)
  mat <- downsample(combined, marker_genes, run)
  # mat[mat != 0] <- 1 # allows for rowSums to give the number of positive cells
  
  Idents(object = combined) <- shuffled
  num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
  gene_names <- rownames(combined@assays$RNA)
  for (i in 0:num_clusters) {
    start_i <- 1 + i*length(valid_genes)
    end_i <- start_i + length(valid_genes)
    this_cells <- WhichCells(combined, idents = i)
    perm_gene_df[start_i:end_i,3] <- perm_gene_df[start_i:end_i,3] + rowSums(mat[,this_cells])
  }
}
# Take the average of the 50 runs
perm_gene_df[,3] <- perm_gene_df[,3] / num_runs

# Do a Fisher's Exact Test on the number of positive cells in Down data vs Down+perm data
# on a per gene per cluster basis
gene_df <- gene_df[complete.cases(gene_df), ]
perm_gene_df <- perm_gene_df[complete.cases(perm_gene_df), ]
for (i in 1:nrow(gene_df)) {
  cluster <- gene_df[i,2] + 1
  no_perm <- c(round(gene_df[i,3]), cells_per_cluster[cluster])
  perm <- c(round(perm_gene_df[i,3]), cells_per_cluster[cluster])
  contig_table <- data.frame(no_perm <- no_perm, perm <- perm)
  fisher_p <- fisher.test(contig_table)$p.value
  gene_df[i,4] <- fisher_p
  if ( (no_perm[1]/no_perm[2]) > (perm[1]/perm[2]) ) {
    gene_df$up <- TRUE
  }
}
gene_df$q <- p.adjust(gene_df$p, method = "hochberg")
gene_df$p_sig <- gene_df$p < 0.05
gene_df$q_sig <- gene_df$q < 0.05
write.csv(gene_df[which(gene_df$p_sig),], file = paste(rna_path, "/results/down_perm_sig_gene_trans_clusters.csv", sep=""), row.names = FALSE)