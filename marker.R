library("stringr")
library("dplyr")
library("Seurat")
library("MASS")
options(warn=-1)

rna_path <- "C:/Users/miles/Downloads/brain/"
# rna_path <- "C:/Users/miles/Downloads/d_tooth/"
# rna_path <- "~/scratch/brain/"
combined <- readRDS(paste(rna_path, "/brain_scripts/brain_shiny/data/combined.rds", sep = ""))
marker_path <- paste(rna_path, "data/markers/", sep="")
marker_files <- dir(marker_path, pattern =paste("*.txt", sep=""))
# marker_files <- dir(marker_path, pattern =paste("*.tsv", sep=""))

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

markers <- data.frame(gene <- c(), bio <- c())
for (i in 1:length(marker_files)) {
  file <- read.table(paste(marker_path, marker_files[i], sep=""), header=FALSE, sep="\t", stringsAsFactors=FALSE)
  # file[,1] <- toupper(file[,1])
  markers <- rbind(markers, file[,1:2])
}
colnames(markers) <- c("gene", "bio")
backup_markers <- markers

# Paint the markers
combined@active.assay <- "RNA"
gene_names <- rownames(combined@assays$RNA)
markers <- markers[which(markers$bio == "DISC_ASE"),]
markers <- unique(markers)
for (i in 1:nrow(markers)) {
  gene <- markers[i, 1]
  bio <- markers[i, 2]
  
  result <- geneCap(gene, gene_names)
  gene <- result[1]
  error <- as.logical(result[2])
  
  if (i %% 50 == 0) {
        print(i)
        print(gene)
  }
    
  if (! error) {
    png(filename = paste(rna_path, "results/painting/", bio, "/", gene, "_by_sample_umap.png", sep=""), width = 1200, height = 500, unit="px")
    p <- FeaturePlot(combined, features = c(gene), split.by = "sample", reduction = "umap", pt.size = 2, label=TRUE, order = TRUE)
    print(p)
    dev.off()
  }
}

# DEG
big_results <- data.frame()
markers_bio <- unique(backup_markers$bio)
for (bio in markers_bio) {
  markers <- backup_markers[which(backup_markers$bio == bio),]
  # markers <- markers[which(markers$bio == "DISC_ASE"),]
  results <- data.frame(gene <- c(), cluster <- c())
  # deg <- read.csv(paste(rna_path, "results/c41_conserved_and_all_markers_by_cluster_121319.csv", sep=""))
  num_clusters <- max(deg$cluster)
  for (i in 1:nrow(markers)) {
    
    gene <- markers[i,1]
    result <- geneCap(gene, gene_names)
    gene <- result[1]
    error <- as.logical(result[2])
    
    # if (! error) {
    this_rows <- deg[which(deg$gene == gene),]
    # this_rows <- deg[which(deg[,2] == gene),]
    if (nrow(this_rows) > 0) {
      results <- rbind(results, data.frame(rep(gene, nrow(this_rows)), this_rows$cluster))
    }
    # }
  }
  if (nrow(results) > 0) {
    colnames(results) <- c("gene", "cluster")
    results$cluster <- factor(results$cluster, levels = 0:num_clusters)
    # results_table <- table(factor(results$cluster, levels = 0:num_clusters))
    # big_results <- rbind(big_results, results_table)
    # rownames(big_results)[nrow(big_results)] <- bio
  }
}
colnames(big_results) <- 0:num_clusters
heatmap(as.matrix(big_results), Rowv=NA, Colv=NA, revC = TRUE, xlab = "Jaw Cluster", main = "Markers in Jaw DEG", scale = "none")
ggplot(results, aes(cluster)) + geom_bar() + ggtitle(paste(bio, "in Cichlid Jaw")) + scale_x_discrete(drop=FALSE)
# 
## Expression per cluster
# markers <- markers[which(markers$bio == "MC_UP"),]
# gene_names <- rownames(combined@assays$RNA)
# num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
# results_expr <- data.frame()
# colnames(results_expr) <- c("gene", "bio", rep(0:40))
# for (i in 1:nrow(markers)) {
# 
#   gene <- markers[i, 1]
#   bio <- markers[i, 2]
#   if (i %% 50 == 0) {
#     print(i)
#     print(gene)
#   }
# 
#   # Gene the gene name in the right format
#   result <- geneCap(gene, gene_names)
#   gene <- result[1]
#   error <- as.logical(result[2])
# 
#   if (! error) {
#     Idents(combined) <- combined$seurat_clusters
#     expr <- FetchData(object = combined, vars = gene)
#     if (nrow(unique(expr)) > 1) {
#       pos <- combined[, which(x = expr > 1)]
#       num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
#       cluster_expr <- c()
#       for (j in 0:num_clusters) {
#         num_cells_in_cluster <- length(WhichCells(combined, idents = j))
#         num_pos_cells_in_cluster <- 0
#         try(num_pos_cells_in_cluster <- length(WhichCells(pos, idents = j)), silent = TRUE)
#         cluster_expr <- c(cluster_expr, num_pos_cells_in_cluster/num_cells_in_cluster)
#       }
#       new_row <- c(gene, bio, as.character(cluster_expr))
#     } else {
#       # Gene has 0 expression
#       print(paste(gene, "has 0 expression"))
#       new_row <- c(gene, bio, rep(0, num_clusters+1))
#     }
#   } else {
#     print(paste(gene, "not found"))
#     new_row <- c(gene, bio, rep("NOT_FOUND", num_clusters+1))
#   }
#   new_row <- t(data.frame(new_row))
#   colnames(new_row) <- c("gene", "bio", rep(0:40))
#   results_expr <- rbind(results_expr, new_row)
#   colnames(results_expr) <- c("gene", "bio", rep(0:40))
# 
# } # end marker gene for
# colnames(results_expr) <- c("gene", "bio", rep(0:40))
# write.csv(results_expr, file = "C:/Users/miles/Downloads/brain/results/mc_up_per_cluster.csv")

## Report the proportion of cells that express the gene and the mean transcript count per cell
# markers <- markers[which(markers$bio == "MC_UP"),]
# gene_names <- rownames(combined@assays$RNA)
# uniq_markers <- unique(markers$gene)
# results_expr <- data.frame()
# genes_found <- c()
# for (i in 1:length(uniq_markers)) {
# 
#   gene <- uniq_markers[i]
#   # gene <- markers[i, 1]
#   # bio <- markers[i, 2]
#   if (i %% 50 == 0) {
#     print(i)
#     print(gene)
#   }
# 
#   result <- geneCap(gene, gene_names)
#   gene <- result[1]
#   error <- as.logical(result[2])
# 
#   if (! error) {
#     genes_found <- c(genes_found, gene)
#     Idents(combined) <- combined$seurat_clusters
#     expr <- FetchData(object = combined, vars = gene)
#     if (nrow(unique(expr)) > 1) {
#       pos <- combined[, which(x = expr > 1)]
#       mean_expr <- mean(pos@assays$RNA@counts[gene,])
#       pct_expr  <- (ncol(pos) / ncol(combined)) * 100
#       new_row <- c(gene, pct_expr, mean_expr)
#       results_expr <- rbind(results_expr, t(as.data.frame(new_row)))
#     } else {
#       new_row <- c(gene, 0, 0)
#       results_expr <- rbind(results_expr, t(as.data.frame(new_row)))
#     }
#   } else {
#     new_row <- c(gene, "NOT_FOUND", "NOT_FOUND")
#     results_expr <- rbind(results_expr, t(as.data.frame(new_row)))
#   }
# }
# colnames(results_expr) <- c("gene", "pct_cells", "mean_trans")
# genes_found[which(!genes_found %in% results_expr$gene)]
# write.csv(results_expr, file = "C:/Users/miles/Downloads/brain/results/mc_up_per_cluster.csv", row.names = FALSE)
# results_expr <- results_expr[which(results_expr$pct_cells != "NOT_FOUND"),]
# write.csv(results_expr, file = "C:/Users/miles/Downloads/brain/results/mc_up_per_cluster_found.csv", row.names = FALSE)

# Find the Average Number of Genes in the List Expressed per cell in that Cluster
gene_names <- rownames(combined@assays$RNA)
uniq_markers <- unique(markers$gene)
Idents(object = combined) <- "seurat_clusters"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
results_expr <- data.frame()
combined@active.assay <- "RNA"
for (i in 1:nrow(markers)) {

  gene <- uniq_markers[i]
  result <- geneCap(gene, gene_names)
  gene <- result[1]
  error <- as.logical(result[2])

  if (! error) {
    expr <- FetchData(object = combined, vars = gene)
    if (nrow(unique(expr)) > 1) {
      pos <- combined[, which(x = expr > 0)]
      cells_per_cluster <- c(gene)
      for (i in 0:num_clusters) {
        if (i %in% pos@active.ident) { # There might not be pos cells in some clusters
          cells_per_cluster <- c(cells_per_cluster, length(WhichCells(pos, idents = i)))
        } else {
          cells_per_cluster <- c(cells_per_cluster, 0)
        }
      }
    }
    results_expr <- rbind(results_expr, t(cells_per_cluster))
  } # end error if
} # end gene for
colnames(results_expr) <- c("gene", 0:40)
cells_per_cluster <- c()
for (i in 0:num_clusters) {
  cells_per_cluster <- c(cells_per_cluster, length(WhichCells(combined, idents = i)))
}
tmp <- results_expr[2:ncol(results_expr)] %>% mutate_all(as.character)
tmp2 <- tmp %>% mutate_all(as.numeric)
genes_per_cluster <- colSums(tmp2)
avg_genes_per_cell_per_cluster <- genes_per_cluster/cells_per_cluster
df <- as.data.frame(avg_genes_per_cell_per_cluster)
df$genes_per_cluster <- genes_per_cluster
df$cells_per_cluster <- cells_per_cluster
write.csv(df, file = "C:/Users/miles/Downloads/brain/results/disc_ase.csv", row.names = TRUE)
background_ratio <- c(sum(df$genes_per_cluster), sum(df$cells_per_cluster))
back_rat_chi <- c(sum(df$genes_per_cluster), sum(nrow(tmp2) * cells_per_cluster) - sum(df$genes_per_cluster))
results <- data.frame()
results_disc_to_local <- data.frame()
results_chi <- data.frame()
num_non_zero_genes <- length(which(rowSums(as.matrix(combined@assays$RNA[rownames(combined),])) > 1))
for (i in 1:(num_clusters+1)) {
  num_genes <- nrow(tmp2)
  total_num_cells <- num_genes * cells_per_cluster[i]
  pos_cells <- df[i,2]
  local_back_rat_chi <- c(df2$genes_per_cluster[i], 22739*cells_per_cluster[i])
  # cont_table <- data.frame(cluster <- c(pos_cells, total_num_cells-pos_cells), back_rat_chi <- back_rat_chi) # global background
  cont_table <- data.frame(cluster <- c(pos_cells, total_num_cells-pos_cells), back_rat_chi <- local_back_rat_chi) # local background
  # chi_p <- chisq.test(cont_table)$p.value # global background
  chi_p <- chisq.test(cont_table)$p.value # local background
  
  # pairwise <- c()
  # for (j in 1:num_clusters+1) {
  #   total_num_cells_j <- num_genes * cells_per_cluster[j]
  #   pos_cells_j <- df[j,2]
  #   local_back_rat_chi_j <- c(df2$genes_per_cluster[j], (num_non_zero_genes * cells_per_cluster[j]) - df2$genes_per_cluster[j])
  #   cont_table <- data.frame(cluster_1 <- c(total_num_cells-2*pos_cells, local_back_rat_chi[2] - local_back_rat_chi[1]), cluster_2 <- c(total_num_cells_j-2*pos_cells_j, local_back_rat_chi_j[2] - local_back_rat_chi_j[1]))
  #   chi_p_j <- chisq.test(cont_table)$p.value # local background
  #   pairwise <- c(pairwise, chi_p_j)
  # }
  disc <- cont_table[1,1]/cont_table[2,1]
  local <- cont_table[1,2]/cont_table[2,2]
  disc_to_local <- disc/local
  results_disc_to_local <- rbind(results_disc_to_local, t(c(i-1, disc, local, disc_to_local)))
  fisher_p <- fisher.test(data.frame(cluster <- t(df[i,2:3]), background_ratio <- background_ratio))$p.value
  results <- rbind(results, t(c(i-1, df[i,2]/df[i,3], fisher_p)))
  results_chi <- rbind(results_chi, t(c(i-1, pos_cells/(total_num_cells-pos_cells), local_back_rat_chi[1]/local_back_rat_chi[2], chi_p)))
}
colnames(results) <- c("cluster", "avg_gene_per_cell_per_cluster", "fisher_p")
colnames(results_chi) <- c("cluster", "disc_ase_pos/neg", "all_pos/neg", "chi_p")
colnames(results_disc_to_local) <- c("cluster", "disc_ase_pos_neg_ratio", "local_background_pos_neg_ratio", "disc_local_ratio")
results$q <- p.adjust(results$fisher_p, method = "bonferroni")
results_chi$chi_q <- p.adjust(results_chi$chi_p, method = "bonferroni")
results_chi$sig_high <- results_chi$chi_q < 0.05 & results_chi$avg_gene_per_cell_per_cluster > back_rat_chi[1] / back_rat_chi[2]
write.csv(results_disc_to_local, file = "C:/Users/miles/Downloads/brain/results/disc_local_ratio.csv", row.names = FALSE)

# Find the Average Number of Genes in the List Expressed per cell in that Cluster - Efficient
markers <- backup_markers[which(backup_markers$bio == "DISC_ASE"),]
valid_genes <- markers$gene
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
rownames(df) <- 0:num_clusters
write.csv(df, file = "C:/Users/miles/Downloads/brain/results/enrichment/b1b2c1mz/disc_ase.csv", row.names = TRUE)

# Find the Average Number of TRANSCRIPTS in the List Expressed per cell in that Cluster
gene_names <- rownames(combined@assays$RNA)
uniq_markers <- unique(markers$gene)
Idents(object = combined) <- "seurat_clusters"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
results_expr <- data.frame()
valid_genes <- c()
combined@active.assay <- "RNA"
for (i in 1:nrow(markers)) {
  
  gene <- uniq_markers[i]
  result <- geneCap(gene, gene_names)
  gene <- result[1]
  error <- as.logical(result[2])
  
  if (! error) {
    valid_genes <- c(valid_genes, gene)
    expr <- FetchData(object = combined, vars = gene)
    if (nrow(unique(expr)) > 1) {
      pos <- combined[, which(x = expr > 0)]
      cells_per_cluster <- c(gene)
      for (i in 0:num_clusters) {
        if (i %in% pos@active.ident) { # There might not be pos cells in some clusters
          trans <- sum(as.matrix(combined@assays$RNA@counts[gene,]))
          cells_per_cluster <- c(cells_per_cluster, trans)
        } else {
          cells_per_cluster <- c(cells_per_cluster, 0)
        }
      }
    }
    results_expr <- rbind(results_expr, t(cells_per_cluster))
  } # end error if
  # else {
  #   cat(paste("Can't find gene", gene, "\n"))
  #   sim_genes <- gene_names[which(startsWith(tolower(gene_names), tolower(gene)))]
  #   cat(paste("Found these instead", sim_genes, "\n"))
  # 
  # 
  #   for (sim_gene in sim_genes) {
  #     expr <- FetchData(object = combined, vars = sim_gene)
  #     if (nrow(unique(expr)) > 1) {
  #       pos <- combined[, which(x = expr > 1)]
  #       cells_per_cluster <- c(sim_gene)
  #       for (i in 0:num_clusters) {
  #         if (i %in% pos@active.ident) { # There might not be pos cells in some clusters
  #           cells_per_cluster <- c(cells_per_cluster, length(WhichCells(pos, idents = i)))
  #         } else {
  #           cells_per_cluster <- c(cells_per_cluster, 0)
  #         }
  #       }
  #     }
  #     results_expr <- rbind(results_expr, t(cells_per_cluster))
  #   }# end sim_gene for
  # }
} # end gene for
colnames(results_expr) <- c("gene", 0:40)
cells_per_cluster <- c()
for (i in 0:num_clusters) {
  cells_per_cluster <- c(cells_per_cluster, length(WhichCells(combined, idents = i)))
}
tmp <- results_expr[2:ncol(results_expr)] %>% mutate_all(as.character)
tmp2 <- tmp %>% mutate_all(as.numeric)
genes_per_cluster <- colSums(tmp2)
avg_trans_per_cell_per_cluster <- genes_per_cluster/cells_per_cluster
df <- as.data.frame(avg_trans_per_cell_per_cluster)
df$trans_per_cluster <- genes_per_cluster
df$cells_per_cluster <- cells_per_cluster
write.csv(df, file = "C:/Users/miles/Downloads/brain/results/disc_ase_trans.csv", row.names = TRUE)

# Avg gene per cluster for ALL genes
Idents(object = combined) <- "seurat_clusters"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
results_expr <- data.frame()
cells_per_cluster <- c()
genes_per_cluster <- c()
box_results <- list()
for (i in 0:num_clusters) {
  print(i)
  this_cells <- WhichCells(combined, idents = i)
  genes_per_cluster <- c(genes_per_cluster, length(which(as.vector(combined@assays$RNA@counts[,this_cells]) != 0)))
  cells_per_cluster <- c(cells_per_cluster, length(this_cells))
  box_results[[i+1]] <- colSums(as.matrix(combined@assays$RNA@counts[,this_cells]) != 0)
}
avg_gene_per_cell_per_cluster <- genes_per_cluster/cells_per_cluster
df <- as.data.frame(avg_gene_per_cell_per_cluster)
df$genes_per_cluster <- genes_per_cluster
df$cells_per_cluster <- cells_per_cluster
rownames(df) <- 0:num_clusters
boxplot(box_results, ylab="Genes per Cell", names = 0:num_clusters)
box_p <- data.frame()
for (i in 1:(length(box_results)-1)) {
  for (j in (i+1):length(box_results)) {
    p <- wilcox.test(box_results[[i]], box_results[[j]])$p.value
    box_p <- rbind(box_p, t(c(i, j, p)))
  }
}
box_p$q <- p.adjust(box_p$V3, method="bonferroni")
box_p$sig <- box_p$q < 0.05
colnames(box_p) <- c("cluster_1", "cluster_2", "p", "q", "sig")
write.csv(box_p, file = "C:/Users/miles/Downloads/brain/results/genes_per_cell_sig.csv", row.names = FALSE)
write.csv(df, file = "C:/Users/miles/Downloads/brain/results/all_avg_gene_per_cluster.csv", row.names = TRUE)

# Avg trans per cluster for ALL genes
Idents(object = combined) <- "seurat_clusters"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
results_expr <- data.frame()
cells_per_cluster <- c()
genes_per_cluster <- c()
sd_per_cluster <- c()
box_results <- list()
for (i in 0:num_clusters) {
  print(i)
  this_cells <- WhichCells(combined, idents = i)
  # genes_per_cluster <- c(genes_per_cluster, length(which(as.vector(combined@assays$RNA@counts[,this_cells]) != 0))) # genes
  genes_per_cluster <- c(genes_per_cluster, sum(colSums(combined@assays$RNA@counts[,this_cells]))) # trans
  cells_per_cluster <- c(cells_per_cluster, length(this_cells))
  sd_per_cluster <- c(sd_per_cluster, sd(colSums(combined@assays$RNA@counts[,this_cells])))
  box_results[[i+1]] <- colSums(combined@assays$RNA@counts[,this_cells])
}
avg_gene_per_cell_per_cluster <- genes_per_cluster/cells_per_cluster
df <- as.data.frame(avg_gene_per_cell_per_cluster)
df$genes_per_cluster <- genes_per_cluster
df$cells_per_cluster <- cells_per_cluster
rownames(df) <- 0:num_clusters
boxplot(box_results, ylab="Transcripts per Cell", names = 0:num_clusters)
box_p <- data.frame()
for (i in 1:(length(box_results)-1)) {
  for (j in (i+1):length(box_results)) {
    p <- wilcox.test(box_results[[i]], box_results[[j]])$p.value
    box_p <- rbind(box_p, t(c(i, j, p)))
  }
}
box_p$q <- p.adjust(box_p$V3, method="bonferroni")
box_p$sig <- box_p$q < 0.05
colnames(box_p) <- c("cluster_1", "cluster_2", "p", "q", "sig")
write.csv(box_p, file = "C:/Users/miles/Downloads/brain/results/trans_per_cell_sig.csv", row.names = FALSE)
write.csv(df, file = "C:/Users/miles/Downloads/brain/results/all_avg_trans_per_cluster_w_sd.csv", row.names = TRUE)

# No thresholding
valid_genes <- validGenes(markers$gene, gene_names)
gene_names_per_cluster <- data.frame()
gene_names_per_cluster2 <- list()
for (i in 0:num_clusters) {
  print(i)
  this_cells <- WhichCells(combined, idents = i)
  this_gene_names <- rownames(combined@assays$RNA@counts[which(rowSums(as.matrix(combined@assays$RNA@counts[valid_genes, this_cells])) != 0),])
  gene_names_per_cluster <- rbind(gene_names_per_cluster, t(c(i, length(this_gene_names), paste(this_gene_names, collapse=","))))
  gene_names_per_cluster2[[i+1]] <- this_gene_names
}
common_genes <- Reduce(intersect, gene_names_per_cluster2)
colnames(gene_names_per_cluster) <- c("cluster", "number_of_disc_ase_genes_expressed", "disc_ase_genes_expressed")
write.table(gene_names_per_cluster, file = "C:/Users/miles/Downloads/brain/results/disc_ase_genes_exp_per_cluster.tsv", sep="\t", row.names = FALSE, quote=FALSE)
write.table(common_genes, file = "C:/Users/miles/Downloads/brain/results/common_disc_ase_genes.tsv", sep="\t", row.names = FALSE, quote=FALSE, col.names = FALSE)

# Thresholding
valid_genes <- validGenes(markers$gene, gene_names)
gene_names_per_cluster <- data.frame()
gene_names_per_cluster2 <- list()
for (i in 0:num_clusters) {
  print(i)
  this_cells <- WhichCells(combined, idents = i)
  threshold <- 0.25 * length(this_cells)
  this_gene_names <- rownames(combined@assays$RNA@counts[which(rowSums(as.matrix(combined@assays$RNA@counts[valid_genes, this_cells])) > threshold),])
  gene_names_per_cluster <- rbind(gene_names_per_cluster, t(c(i, length(this_gene_names), paste(this_gene_names, collapse=","))))
  gene_names_per_cluster2[[i+1]] <- this_gene_names
}
common_genes <- Reduce(intersect, gene_names_per_cluster2)
colnames(gene_names_per_cluster) <- c("cluster", "number_of_disc_ase_genes_expressed", "disc_ase_genes_expressed")
write.table(gene_names_per_cluster, file = "C:/Users/miles/Downloads/brain/results/disc_ase_genes_exp_per_cluster_25_thresh.tsv", sep="\t", row.names = FALSE, quote=FALSE)
write.table(common_genes, file = "C:/Users/miles/Downloads/brain/results/common_disc_ase_genes.tsv", sep="\t", row.names = FALSE, quote=FALSE, col.names = FALSE)

# Number of unique genes expressed per cluster
genes_per_cluster <- c()
for (i in 0:num_clusters) {
  this_cells <- WhichCells(combined, idents = i)
  genes_per_cluster <- c(genes_per_cluster, length(which(rowSums(as.matrix(combined@assays$RNA@counts[,this_cells])) != 0))) # genes
}

# Histograms - num cells
library("ggplot2")
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
ran <- backup_markers[which(backup_markers$bio == "RAN"),1]
rs <- backup_markers[which(backup_markers$bio == "ROCK_SAND"),1]
num_cells <- data.frame()
for (gene in ran) {
  num_cells <- rbind(num_cells, t(c("RAN", length(which(as.vector(combined@assays$RNA@counts[gene,]) != 0)))))
}
for (gene in rs) {
  num_cells <- rbind(num_cells, t(c("ROCK_SAND", length(which(as.vector(combined@assays$RNA@counts[gene,]) != 0)))))
}
colnames(num_cells) <- c("type", "num_cells")
num_cells$num_cells <- as.numeric(as.vector(num_cells$num_cells)) + 0.1
ggplot(num_cells, aes(x=num_cells, fill = type, color=type)) + geom_histogram(stat = "bin", alpha=0.5, position="identity") + ggtitle("Number of Cells Expressing a Random List of Genes vs Rock-Sand Brain DEGs") + xlab("Number of Cells") + ylab("Number of Genes")
# p+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

# Histograms - num cluster
ran <- backup_markers[which(backup_markers$bio == "RAN"),1]
rs <- backup_markers[which(backup_markers$bio == "ROCK_SAND"),1]
num_clust <- data.frame()
threshold <- 0.05
for (gene in ran) {
  num_clust_gene <- 0
  for (i in 0:num_clusters) {
    this_cells <- WhichCells(combined, idents = i)
    num_cells <- length(which(as.vector(combined@assays$RNA@counts[gene,this_cells]) != 0))
    if( num_cells/length(this_cells) > threshold ) {
      num_clust_gene <- num_clust_gene + 1
    }
  }
  num_clust <- rbind(num_clust, t(c("RAN", num_clust_gene)))
  # unique(combined$seurat_clusters[which(as.vector(combined@assays$RNA@counts[gene,]) != 0)])
  # num_cells <- rbind(num_cells, t(c("RAN", length(which(as.vector(combined@assays$RNA@counts[gene,]) != 0)))))
}
for (gene in rs) {
  num_clust_gene <- 0
  for (i in 0:num_clusters) {
    this_cells <- WhichCells(combined, idents = i)
    num_cells <- length(which(as.vector(combined@assays$RNA@counts[gene,this_cells]) != 0))
    if( num_cells/length(this_cells) > threshold ) {
      num_clust_gene <- num_clust_gene + 1
    }
  }
  num_clust <- rbind(num_clust, t(c("ROCK_SAND", num_clust_gene)))
  # num_cells <- rbind(num_cells, t(c("ROCK_SAND", length(which(as.vector(combined@assays$RNA@counts[gene,]) != 0)))))
}
colnames(num_clust) <- c("type", "num_clusters")
# num_clust$num_clusters <- as.numeric(as.vector(num_clust$num_clusters)) + 0.1
# num_clust$num_clusters <- as.factor(num_clust$num_clusters)
num_clust$num_clusters <- as.numeric(as.vector(num_clust$num_clusters))
# ggplot(num_clust, aes(x=num_clusters, fill = type, color=type)) + geom_histogram(stat = "bin", alpha=0.5, position="identity") + ggtitle(paste("Number of Clusters Expressing RAN vs Rock-Sand at ", threshold*100, "%", sep="")) + xlab("Number of Cells") + ylab("Number of Genes")
ggplot(num_clust, aes(x=num_clusters, fill = type, color=type)) + geom_bar(alpha=0.5, position="identity") + ggtitle(paste("Number of Clusters Expressing RAN vs Rock-Sand at ", threshold*100, "%", sep="")) + xlab("Number of Clusters") + ylab("Number of Genes")
# p+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

# p1 <- hist(as.numeric(num_cells[which(num_cells$type == "RAN"), 2]), col=rgb(1,0,0,0.5), xlim = c(0,800), ylim = c(0,300))
# p2 <- hist(as.numeric(num_cells[which(num_cells$type == "ROCK_SAND"), 2]), col=rgb(0,0,1,0.5), xlim = c(0,800), ylim = c(0,300))
# plot(p1, col=rgb(0,0,1,1/4), xlim=c(0,800), ylim = c(0,350))
# plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,800), ylim = c(0,350), add=T)
# legend("topright", legend=c("RAN", "ROCK_SAND"), fill=c("blue", "red"))

# Paint the average number of genes from the list expressed by each cell
bio <- "DISC_ASE"
markers <- backup_markers[which(backup_markers$bio == bio),]
valid_genes <- markers$gene
Idents(object = combined) <- "seurat_clusters"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
gene_names <- rownames(combined@assays$RNA)
genes_per_cell <- c()
trans_per_cell <- c()
mat <- as.matrix(combined@assays$RNA@counts[valid_genes,])
trans_per_cell <- colSums(mat)
mat[mat != 0] <- 1
genes_per_cell <- colSums(mat)
combined$genes_per_cell <- genes_per_cell
combined$trans_per_cell <- trans_per_cell
combined$genes_per_cell_norm_all_trans <- genes_per_cell/combined$nCount_RNA
combined$trans_per_cell_norm_all_trans <- trans_per_cell/combined$nCount_RNA
# matrix2 <- matrix(genes_per_cell, nrow = 1, ncol = ncol(combined), dimnames = list("geneList", colnames(combined)))
# combined[["geneList"]] <- CreateAssayObject(counts = matrix2)
# matrix3 <- matrix(trans_per_cell, nrow = 1, ncol = ncol(combined), dimnames = list("transList", colnames(combined)))
# combined[["transList"]] <- CreateAssayObject(counts = matrix3)

png(paste0(rna_path, "/results/enrichment/b1b2c1mz/background_trans_all_genes.png"), width = 2000, height = 1200, res = 200)
p <- FeaturePlot(combined, features = "nCount_RNA", reduction = "umap", label = TRUE, pt.size = 2, order = TRUE) + ggtitle(paste0("Transcripts per Cell from All Genes"))
print(p)
dev.off()
png(paste0(rna_path, "/results/enrichment/b1b2c1mz/genes_per_cell_", bio, ".png"), width = 2000, height = 1200, res = 200)
p <- FeaturePlot(combined, features = "genes_per_cell", reduction = "umap", label = TRUE, pt.size = 2, order = TRUE) + ggtitle(paste0("Genes per Cell from Gene List: ", bio))
print(p)
dev.off()
png(paste0(rna_path, "/results/enrichment/b1b2c1mz/trans_per_cell_", bio, ".png"), width = 2000, height = 1200, res = 200)
p <- FeaturePlot(combined, features = "trans_per_cell", reduction = "umap", label = TRUE, pt.size = 2, order = TRUE) + ggtitle(paste0("Transcripts per Cell from Gene List: ", bio))
print(p)
dev.off()
png(paste0(rna_path, "/results/enrichment/b1b2c1mz/genes_per_cell_norm_", bio, ".png"), width = 2000, height = 1200, res = 200)
p <- FeaturePlot(combined, features = "genes_per_cell_norm_all_trans", reduction = "umap", label = TRUE, pt.size = 2, order = TRUE) + ggtitle(paste0("Genes per Cell from Gene List: ", bio, ", Normalized by All Transcripts"))
print(p)
dev.off()
png(paste0(rna_path, "/results/enrichment/b1b2c1mz/trans_per_cell_norm_", bio, ".png"), width = 2000, height = 1200, res = 200)
p <- FeaturePlot(combined, features = "trans_per_cell_norm_all_trans", reduction = "umap", label = TRUE, pt.size = 2, order = TRUE) + ggtitle(paste0("Transcripts per Cell from Gene List: ", bio, ", Normalized by All Transcripts"))
print(p)
dev.off()

# Histogram of total #transcripts for all genes in set (x axis) and number of nuclei on y axis (number of nuclei that fall into that expression bin)
iegs <- read.csv("C:/Users/miles/Downloads/zack_IEG_list_061720.csv", header = FALSE, stringsAsFactors = F)
valid_genes <- iegs$V1 # TODO change this to whatever list Zack wants
tops <- c(10, 20, 50, 100)
for (top in tops) {
  combined$thisSum <- colSums(combined@assays$RNA@counts[valid_genes[1:top],])
  df <- data.frame(thisSum<-combined$thisSum)
  # combined$this
  p <- ggplot(df, aes(thisSum)) + geom_bar() + xlab("# transcripts for all genes in set") + ylab("Number of Nuclei With that # Transcripts") + ggtitle(paste0("Top ", top, " IEGs"))
  print(p)
}

# Do cells with higher expression also have a greater transcriptional diversity (more genes with >1 transcripts)?
mat = combined@assays$RNA@counts
mat[which(mat > 0)]=1
df = data.frame(counts = combined$nCount_RNA, data = colSums(combined@assays$RNA@data), trans_div = colSums(mat))
ggplot(df, aes(x=trans_div, y=counts)) + geom_point() + labs(x="Transcriptional Diversity", y="Total Transcripts",title="Transcriptional Diversity vs Total Transcripts in a Cell")
ggplot(df, aes(x=trans_div, y=data)) + geom_point() + labs(x="Transcriptional Diversity", y="Total Normalized Transcripts", title="Transcriptional Diversity vs Total Normalized Transcripts in a Cell")
