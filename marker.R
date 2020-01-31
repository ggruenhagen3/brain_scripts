library("stringr")
library("dplyr")
library("Seurat")
options(warn=-1)

rna_path <- "C:/Users/miles/Downloads/brain/"
rna_path <- "~/scratch/brain/"
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

markers <- data.frame(gene <- c(), bio <- c())
for (i in 1:length(marker_files)) {
  file <- read.table(paste(marker_path, marker_files[i], sep=""), header=FALSE, sep="\t", stringsAsFactors=FALSE)
  file[,1] <- toupper(file[,1])
  markers <- rbind(markers, file[,1:2])
}
colnames(markers) <- c("gene", "bio")
# markers[564,] <- c("SST", "test")

# Paint the markers
gene_names <- rownames(combined@assays$RNA)
markers <- markers[which(markers$bio == "ASE"),]
markers <- unique(markers)
for (i in 1:nrow(markers)) {
  gene <- markers[i, 1]
  bio <- markers[i, 2]
  
  result <- geneCap(gene, gene_names)
  gene <- result[1]
  error <- as.logical(result[2])
  print(gene)
  print(error)
  if (! error) {
    png(filename = paste(rna_path, "results/painting/", bio, "/", gene, "_umap.png", sep=""), width = 900, height = 500, unit="px")
    p <- FeaturePlot(combined, features = c(gene), split.by = "cond", reduction = "umap", pt.size = 2, label=TRUE, order = TRUE)
    print(p)
    dev.off()
  }
}
# 
# # DEG
# results <- data.frame(cluster <- c(), gene <- c(), bio <- c())
# for (i in 1:length(cluster_files)) {
#   tryCatch({
#     cluster_genes <- read.table(paste(cluster_path, cluster_files[i], sep=""), sep ="\n")
#     cluster_genes <- cluster_genes[,1]
#     cluster_genes <- toupper(cluster_genes)
#   },  error=function(e) NULL)
# 
#   new_results <- markers[which(markers$gene %in% cluster_genes),]
#   if (nrow(new_results) > 0) {
#     new_results <- cbind(new_results, rep(i, nrow(new_results)))
#     results <- rbind(results, new_results)
#   }
# }
# colnames(results) <- c("gene", "bio", "cluster")
# 
# ## Expression per cluster
# # markers <- markers[which(markers$bio != "LG11_HIGH_FST"),]
# gene_names <- rownames(combined@assays$RNA)
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
#   error <- result[2]
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
# write.csv(results_expr, file = "C:/Users/miles/Downloads/brain/results/results_expr_z_n_b.csv")
# 
# 
# ## Report the proportion of cells that express the gene and the mean transcript count per cell
# # markers <- markers[which(markers$bio == "LG11_HIGH_FST"),]
# markers <- markers[which(markers$bio == "ASE"),]
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
# write.csv(results_expr, file = "C:/Users/miles/Downloads/brain/results/results_expr_ASE.csv", row.names = FALSE)
# results_expr <- results_expr[which(results_expr$pct_cells != "NOT_FOUND"),]
# write.csv(results_expr, file = "C:/Users/miles/Downloads/brain/results/results_expr_ASE_found.csv", row.names = FALSE)
# 
# 
# # Find the Average Number of Genes in the List Expressed per cell in that Cluster
# gene_names <- rownames(combined@assays$RNA)
# uniq_markers <- unique(markers$gene)
# Idents(object = combined) <- "seurat_clusters"
# num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
# results_expr <- data.frame()
# for (i in 1:nrow(markers)) {
#   
#   gene <- uniq_markers[i]
#   result <- geneCap(gene, gene_names)
#   gene <- result[1]
#   error <- as.logical(result[2])
#   
#   if (! error) {
#     expr <- FetchData(object = combined, vars = gene)
#     if (nrow(unique(expr)) > 1) {
#       pos <- combined[, which(x = expr > 1)]
#       cells_per_cluster <- c(gene)
#       for (i in 0:num_clusters) {
#         if (i %in% pos@active.ident) { # There might not be pos cells in some clusters
#           cells_per_cluster <- c(cells_per_cluster, length(WhichCells(pos, idents = i)))
#         } else {
#           cells_per_cluster <- c(cells_per_cluster, 0)
#         }
#       }
#     }
#     results_expr <- rbind(results_expr, t(cells_per_cluster))
#   } # end error if
#   else {
#     cat(paste("Can't find gene", gene, "\n"))
#     sim_genes <- gene_names[which(startsWith(tolower(gene_names), tolower(gene)))]
#     cat(paste("Found these instead", sim_genes, "\n"))
#     
#     
#     for (sim_gene in sim_genes) {
#       expr <- FetchData(object = combined, vars = sim_gene)
#       if (nrow(unique(expr)) > 1) {
#         pos <- combined[, which(x = expr > 1)]
#         cells_per_cluster <- c(sim_gene)
#         for (i in 0:num_clusters) {
#           if (i %in% pos@active.ident) { # There might not be pos cells in some clusters
#             cells_per_cluster <- c(cells_per_cluster, length(WhichCells(pos, idents = i)))
#           } else {
#             cells_per_cluster <- c(cells_per_cluster, 0)
#           }
#         }
#       }
#       results_expr <- rbind(results_expr, t(cells_per_cluster))
#     }# end sim_gene
#   }
# } # end gene for
# colnames(results_expr) <- c("gene", 0:40)
# cells_per_cluster <- c()
# for (i in 0:num_clusters) {
#   cells_per_cluster <- c(cells_per_cluster, length(WhichCells(combined, idents = i)))
# }
# tmp <- results_expr[2:ncol(results_expr)] %>% mutate_all(as.character)
# tmp2 <- tmp %>% mutate_all(as.numeric)
# genes_per_cluster <- colSums(tmp2)
# avg_genes_per_cell_per_cluster <- genes_per_cluster/cells_per_cluster
# df <- as.data.frame(avg_genes_per_cell_per_cluster)
# df$genes_per_cluster <- genes_per_cluster
# df$cells_per_cluster <- cells_per_cluster
# write.csv(df, file = "C:/Users/miles/Downloads/brain/results/avg_genes_per_cell_per_cluster_converted_names.csv", row.names = TRUE)
