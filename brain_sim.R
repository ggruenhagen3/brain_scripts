library("edgeR")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("dplyr")
library("cowplot")
library("ggplot2")
# set.seed(1)
# rm(.Random.seed)

rna_path <- "C:/Users/miles/Downloads/brain/"
combined <- readRDS("C:/Users/miles/Downloads/brain/brain_scripts/brain_mz_shiny/data/B1C1C2MZ_combined_031020.rds")

Idents(combined) <- "sample"
b1 <- combined[, WhichCells(combined, idents = "b1")]
b2 <- combined[, WhichCells(combined, idents = "b2")]
c1 <- combined[, WhichCells(combined, idents = "c1")]
mz <- combined[, WhichCells(combined, idents = "c2")]
Idents(b1) <- "seurat_clusters"
Idents(b2) <- "seurat_clusters"
Idents(c1) <- "seurat_clusters"
Idents(mz) <- "seurat_clusters"
Idents(combined) <- "seurat_clusters"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
rock_unique <- c(13, 29)
sand_unique <- c(1, 2)
unique_clusters <- c(rock_unique, sand_unique)
num_common_clusters <- (num_clusters+1) - length(unique_clusters)

num_pools <- 4
num_cells <- 3000

# Sample BHVE
bhve_pools <- list()
clusters_to_use <- c(40:30, 28:14, 12:0)
for (p in 1:num_pools) {
  # Sample Cells Proportionally from each cluster
  new_cells_b1 <- c()
  new_cells_b2 <- c()
  for (i in clusters_to_use) {
    this_cluster <- WhichCells(combined, idents = i)
    pro_real <- length(this_cluster)/ncol(combined)
    num_cluster_sim_cells <- pro_real * (num_cells/2)
    if (i == 0) { # fill the rest of the cells with all cluster 0 cells
      num_cluster_sim_cells <- (num_cells/2) - length(new_cells_b1)
    }
    new_cells_b1 <- c(new_cells_b1, WhichCells(b1, idents = i, downsample = num_cluster_sim_cells, seed = p))
    new_cells_b2 <- c(new_cells_b2, WhichCells(b2, idents = i, downsample = num_cluster_sim_cells, seed = p))
  }
  
  # new_cells_b1 <- sample(colnames(b1@assays$RNA@counts), num_cells/2, replace = TRUE)
  # new_cells_b2 <- sample(colnames(b2@assays$RNA@counts), num_cells/2, replace = TRUE)
  
  new_matrix <- cbind(b1@assays$RNA@counts[,new_cells_b1], b2@assays$RNA@counts[,new_cells_b2])
  # new_data_matrix <- cbind(b1@assays$RNA@data[,new_cells_b1], b2@assays$RNA@data[,new_cells_b2])
  
  bhve_pools[[p]] <- CreateSeuratObject(new_matrix, project = paste0("pool", p))
  bhve_pools[[p]] <- NormalizeData(bhve_pools[[p]], normalization.method = "LogNormalize", scale.factor = 100000)
  bhve_pools[[p]]$cond <- paste0("bhve_", p)
  bhve_pools[[p]]$orig.cluster <- c(as.vector(b1$seurat_clusters[new_cells_b1]), as.vector(b2$seurat_clusters[new_cells_b2]))
  # saveRDS(bhve_pools[[p]], paste0(rna_path, "/data/brain_sim/bhve_", p, ".rds"))
} # end pools for

# Sample CTRL - Method 1
ctrl_pools_1 <- list()
for (p in 1:num_pools) {
  new_cells_c1 <- c()
  for (i in clusters_to_use) {
    this_cluster <- WhichCells(combined, idents = i)
    pro_real <- length(this_cluster)/ncol(combined)
    num_cluster_sim_cells <- pro_real * (num_cells)
    if (i == 0) { # fill the rest of the cells with all cluster 0 cells
      num_cluster_sim_cells <- (num_cells) - length(new_cells_c1)
    }
    set.seed(i)
    # new_cells_c1 <- c(new_cells_c1, sample(WhichCells(c1, idents = i), num_cluster_sim_cells, replace=TRUE))
    new_cells_c1 <- c(new_cells_c1, WhichCells(c1, idents = i, downsample = num_cluster_sim_cells, seed = p))
  }
  
  # Randomly sample 3,000 cells
  ctrl_pools_1[[p]] <- CreateSeuratObject(c1@assays$RNA@counts[,c(new_cells_c1)], project = paste0("pool", i))
  ctrl_pools_1[[p]] <- NormalizeData(ctrl_pools_1[[p]], normalization.method = "LogNormalize", scale.factor = 100000)
  ctrl_pools_1[[p]]$cond <- paste0("ctrl_1_", p)
  ctrl_pools_1[[p]]$orig.cluster <- as.vector(c1$seurat_clusters[new_cells_c1])
  # saveRDS(ctrl_pools_1[[p]], paste0(rna_path, "/data/brain_sim/ctrl_1_", p, ".rds"))
} # end pools for

# Sample CTRL - Method 2
# ctrl_pools_2 <- list()
# for (i in 1:num_pools) {
#   obj <- c1
#   
#   # Randomly pick c1/mz 
#   if ( sample(c(1,2), 1) == 1 ) {
#     obj <- c1
#   } else {
#     obj <- mz
#   } # end if
#   
#   # Randomly sample 3,000 cells
#   new_cells <- sample(colnames(obj@assays$RNA@counts), num_cells, replace = TRUE)
#   
#   ctrl_pools_2[[i]] <- CreateSeuratObject(obj@assays$RNA@counts[,c(new_cells)], project = paste0("pool", i))
#   ctrl_pools_2[[i]]$cond <- paste0("ctrl_2_", i)
#   saveRDS(ctrl_pools_2[[i]], paste0(rna_path, "/data/brain_sim/ctrl_2_", i, ".rds"))
# } # end pools for

# Cluster
bhve <- merge(x = bhve_pools[[1]], y = unlist(bhve_pools[2:num_pools]))
ctrl <- merge(x = ctrl_pools_1[[1]], y = unlist(ctrl_pools_1[2:num_pools]))
sim <- merge(x = bhve, y = ctrl)
sim <- NormalizeData(sim, normalization.method = "LogNormalize", scale.factor = 100000)
sim <- FindVariableFeatures(object = sim, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
sim <- ScaleData(object = sim, vars.to.regress = NULL)
sim <- RunPCA(sim, npcs = 50, verbose = FALSE)
sim <- RunUMAP(sim, reduction = "pca", dims = 1:30)
sim <- FindNeighbors(sim, reduction = "umap", dims = 1:2)
sim <- FindClusters(sim, resolution = 0.1)
DimPlot(sim)
DimPlot(sim, split.by = "orig.ident") + NoLegend()

num_clusters <- as.numeric(tail(levels(sim@meta.data$seurat_clusters), n=1))
cells_cluster <- c()
for (i in 0:num_clusters ) {
  cells_cluster <- c(cells_cluster, length(WhichCells(sim, idents = i)))
}

sim$cond_cluster <- paste0(sim$orig.ident, sim$seurat_clusters) # it's actually orig.ident, but who cares
alt_markers <- data.frame()
Idents(sim) <- sim$cond_cluster
for (i in 0:num_clusters) {
  cat(paste0(i, ": "))
  new_row <- FindMarkers(sim, ident.1 = paste0("BHVE", i), ident.2 = paste0("CTRL", i))
  new_row$cluster <- i
  new_row$gene <- rownames(new_row)
  rownames(new_row) <- NULL
  alt_markers <- rbind(alt_markers, new_row)
}
alt_markers$q <- p.adjust(alt_markers$p_val, method = "bonferroni")
nrow(alt_markers[which(alt_markers$q < 0.05),])
length(unique(alt_markers$gene[which(alt_markers$q < 0.05)]))

sim$cond_cluster <- paste0(sim$orig.ident, sim$orig.cluster) # it's actually orig.ident, but who cares
alt_markers_orig <- data.frame()
Idents(sim) <- sim$cond_cluster
for (i in clusters_to_use) {
  cat(paste0(i, ": "))
  new_row <- FindMarkers(sim, ident.1 = paste0("BHVE", i), ident.2 = paste0("CTRL", i))
  new_row$cluster <- i
  new_row$gene <- rownames(new_row)
  rownames(new_row) <- NULL
  alt_markers_orig <- rbind(alt_markers_orig, new_row)
}
alt_markers_orig$q <- p.adjust(alt_markers_orig$p_val, method = "bonferroni")
nrow(alt_markers_orig[which(alt_markers_orig$q < 0.05),])
length(unique(alt_markers_orig$gene[which(alt_markers_orig$q < 0.05)]))

# alt_markers$gene <- rownames(alt_markers)
# alt_markers_orig$gene <- rownames(alt_markers_orig)
write.table(alt_markers, "C:/Users/miles/Downloads/brain/results/sim_new_cluster_sig_deg.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(alt_markers_orig, "C:/Users/miles/Downloads/brain/results/sim_old_cluster_sig_deg.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
###################
# Null Hypothesis #
###################
null_bhve_pools <- list()
clusters_to_use <- c(39:30, 28:14, 12:0)
for (p in 1:num_pools) {
  
  new_cells <- list(c(), c(), c(), c())
  samples <- list(b1, b2, c1, mz)
  num_samples <- length(samples)
  for (j in 1:num_samples) {
    sample <- samples[[j]]
    for (i in clusters_to_use) {
      this_cluster <- WhichCells(combined, idents = i)
      pro_real <- length(this_cluster)/ncol(combined) + 0.0005
      num_cluster_sim_cells <- ceiling(pro_real * (num_cells/num_samples))
      if (i == 0) { # fill the rest of the cells with all cluster 0 cells
        num_cluster_sim_cells <- (num_cells/num_samples) - length(new_cells[[j]])
      }
      # print(pro_real)
      # print(num_cluster_sim_cells)
      new_cells[[j]] <- c(new_cells[[j]], WhichCells(sample, idents = i, downsample = num_cluster_sim_cells, seed = p))
    }
  }

  
  # Randomly sample 750 cells from each individual
  # new_cells_b1 <- sample(colnames(b1@assays$RNA@counts), num_cells/num_pools, replace = TRUE)
  # new_cells_b2 <- sample(colnames(b2@assays$RNA@counts), num_cells/num_pools, replace = TRUE)
  # new_cells_c1 <- sample(colnames(c1@assays$RNA@counts), num_cells/num_pools, replace = TRUE)
  # new_cells_mz <- sample(colnames(mz@assays$RNA@counts), num_cells/num_pools, replace = TRUE)
  
  new_matrix <- cbind(b1@assays$RNA@counts[,new_cells[[1]]], b2@assays$RNA@counts[,new_cells[[2]]], c1@assays$RNA@counts[,new_cells[[3]]], mz@assays$RNA@counts[,new_cells[[4]]])
  # new_data_matrix <- cbind(b1@assays$RNA@data[,new_cells[[1]]], b2@assays$RNA@data[,new_cells[[2]]], c1@assays$RNA@data[,new_cells[[3]]], mz@assays$RNA@data[,new_cells[[4]]])
  
  null_bhve_pools[[p]] <- CreateSeuratObject(new_matrix, project = paste0("pool", p))
  # null_bhve_pools[[i]] <- SetAssayData(object = obj_b_2, slot = 'data', new.data = new_data_matrix)
  null_bhve_pools[[p]] <- NormalizeData(null_bhve_pools[[p]], normalization.method = "LogNormalize", scale.factor = 100000)
  null_bhve_pools[[p]]$cond <- paste0("bhve_", p)
  null_bhve_pools[[p]]$sim.orig.ident <- "BHVE"
  null_bhve_pools[[p]]$orig.cluster <- c(b1$seurat_clusters[new_cells[[1]]], b2$seurat_clusters[new_cells[[2]]], c1$seurat_clusters[new_cells[[3]]], mz$seurat_clusters[new_cells[[4]]])
  # saveRDS(null_bhve_pools[[p]], paste0(rna_path, "/data/brain_sim/null_bhve_", p, ".rds"))
} # end pools for

null_ctrl_pools <- list()
for (p in 1:num_pools) {
  
  new_cells <- list(c(), c(), c(), c())
  samples <- c(b1, b2, c1, mz)
  num_samples <- length(samples)
  for (j in 1:num_samples) {
    sample <- samples[[j]]
    for (i in clusters_to_use) {
      this_cluster <- WhichCells(combined, idents = i)
      pro_real <- length(this_cluster)/ncol(combined) + 0.0005
      num_cluster_sim_cells <- ceiling(pro_real * (num_cells/num_samples))
      if (i == 0) { # fill the rest of the cells with all cluster 0 cells
        num_cluster_sim_cells <- (num_cells/num_samples) - length(new_cells[[j]])
      }
      new_cells[[j]] <- c(new_cells[[j]], WhichCells(sample, idents = i, downsample = num_cluster_sim_cells, seed = p))
    }
  }
  
  # Randomly sample 750 cells from each individual
  # new_cells_b1 <- sample(colnames(b1@assays$RNA@counts), num_cells/num_pools, replace = TRUE)
  # new_cells_b2 <- sample(colnames(b2@assays$RNA@counts), num_cells/num_pools, replace = TRUE)
  # new_cells_c1 <- sample(colnames(c1@assays$RNA@counts), num_cells/num_pools, replace = TRUE)
  # new_cells_mz <- sample(colnames(mz@assays$RNA@counts), num_cells/num_pools, replace = TRUE)
  
  new_matrix <- cbind(b1@assays$RNA@counts[,new_cells[[1]]], b2@assays$RNA@counts[,new_cells[[2]]], c1@assays$RNA@counts[,new_cells[[3]]], mz@assays$RNA@counts[,new_cells[[4]]])
  
  null_ctrl_pools[[p]] <- CreateSeuratObject(new_matrix, project = paste0("pool", p))
  null_ctrl_pools[[p]] <- NormalizeData(null_ctrl_pools[[p]], normalization.method = "LogNormalize", scale.factor = 100000)
  null_ctrl_pools[[p]]$cond <- paste0("ctrl_", p)
  null_ctrl_pools[[p]]$sim.orig.ident <- "CTRL"
  null_ctrl_pools[[p]]$orig.cluster <- c(b1$seurat_clusters[new_cells[[1]]], b2$seurat_clusters[new_cells[[2]]], c1$seurat_clusters[new_cells[[3]]], mz$seurat_clusters[new_cells[[4]]])
  # saveRDS(null_ctrl_pools[[p]], paste0(rna_path, "/data/brain_sim/null_ctrl_", p, ".rds"))
} # end pools for

null_bhve <- merge(x = null_bhve_pools[[1]], y = unlist(null_bhve_pools[2:num_pools]))
null_ctrl <- merge(x = null_ctrl_pools[[1]], y = unlist(null_ctrl_pools[2:num_pools]))
null_sim <- merge(x = null_bhve, y = null_ctrl, merge.data = TRUE)
# null_sim <- NormalizeData(null_sim, normalization.method = "LogNormalize", scale.factor = 100000)
null_sim <- FindVariableFeatures(object = null_sim, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
null_sim <- ScaleData(object = null_sim, vars.to.regress = NULL)
null_sim <- RunPCA(null_sim, npcs = 50, verbose = FALSE)
null_sim <- RunUMAP(null_sim, reduction = "pca", dims = 1:30)
null_sim <- FindNeighbors(null_sim, reduction = "umap", dims = 1:2)
null_sim <- FindClusters(null_sim, resolution = 0.1)
DimPlot(null_sim)
DimPlot(null_sim) + theme(legend.position = "none")

Idents(null_sim) <- "sim.orig.ident"
null_sim$cond_cluster <- paste0(null_sim$sim.orig.ident, null_sim$seurat_clusters) # it's actually orig.ident, but who cares
null_markers <- data.frame()
Idents(null_sim) <- null_sim$cond_cluster
for (i in 0:num_clusters) {
  cat(paste0(i, ": "))
  new_row <- NULL
  try(new_row <- FindMarkers(null_sim, ident.1 = paste0("BHVE", i), ident.2 = paste0("CTRL", i)), silent = FALSE)
  new_row$cluster <- i
  null_markers <- rbind(null_markers, new_row)
}
null_markers$q <- p.adjust(null_markers$p_val, method = "bonferroni")
nrow(null_markers[which(null_markers$q < 0.05),])

Idents(null_sim) <- "sim.orig.ident"
null_sim$cond_cluster <- paste0(null_sim$sim.orig.ident, null_sim$orig.cluster) # it's actually orig.ident, but who cares
null_markers_orig <- data.frame()
Idents(null_sim) <- null_sim$cond_cluster
for (i in clusters_to_use) {
  cat(paste0(i, ": "))
  new_row <- NULL
  try(new_row <- FindMarkers(null_sim, ident.1 = paste0("BHVE", i), ident.2 = paste0("CTRL", i)), silent = FALSE)
  new_row$cluster <- i
  null_markers_orig <- rbind(null_markers_orig, new_row)
}
null_markers_orig$q <- p.adjust(null_markers_orig$p_val, method = "bonferroni")
nrow(null_markers_orig[which(null_markers_orig$q < 0.05),])

# Plots
identical(colnames(ctrl_pools_1[[1]]), colnames(ctrl_pools_1[[2]]))
sim$orig.cluster <- as.numeric(sim$orig.cluster)
ggplot(as.data.frame(sim$orig.cluster), aes(sim$orig.cluster)) + geom_bar() + ggtitle("Proportion of Old Clusters in the Alt Simulated Data")
null_sim$orig.cluster <- as.numeric(null_sim$orig.cluster)
ggplot(as.data.frame(null_sim$orig.cluster), aes(null_sim$orig.cluster)) + geom_bar() + ggtitle("Proportion of Old Clusters in the Null Simulated Data")
ggplot(as.data.frame(combined$seurat_clusters), aes(combined$seurat_clusters)) + geom_bar() + ggtitle("Proportion of Old Clusters in the Original Data")
