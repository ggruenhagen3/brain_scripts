library("stringr")
library("ggplot2")
library("biomaRt")
library("DropSeq.util")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("cowplot")
library("RColorBrewer")
library("dplyr")
library("pSI")

source("/nv/hp10/ggruenhagen3/scratch/brain/brain_scripts/all_f.R")

# path <- "C:/Users/miles/Downloads/brain/data/9_brain/"
top_global_path <- "/nv/hp10/ggruenhagen3/scratch/brain/data/9_brain/"
global_path <- paste0(top_global_path, "region/")
regions <- dir(global_path, pattern =paste("*", sep=""))
regions <- c("globulus_pallidus", "hippocampus", "striatum")
regions_abr <- c("GP", "HC", "STR")

cell_types <- dir(paste0(top_global_path, "/cell_type/"), pattern =paste("*", sep=""))

# print("Finding pSI for our object")
combined <- readRDS("/nv/hp10/ggruenhagen3/scratch/brain/data/B1C1C2MZ_combined_031020.rds")
# our_psi <- AverageExpression(combined, slot = counts)
all_deg <- read.csv("/nv/hp10/ggruenhagen3/scratch/brain/data/all_markers_B1C1C2MZ_042820.csv", stringsAsFactors = FALSE)
all_deg <- convertMzebraDFToMouse(all_deg, 8)
ann <- readRDS(paste0(top_global_path, "annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS"))

# ## By CLUSTER ##
# all_fisher <- data.frame()
# for (i in 1:length(regions)) {
#   region <- regions[i]
#   obj_str <- region
#   obj_abr <- regions_abr[i]
#   region_ann <- ann[which(ann$tissue == obj_abr),]
#   region_ann$cluster <- as.numeric(lapply(region_ann$subcluster, function(x) str_split(x, "-")[[1]][1]))
#   print(paste("Loading data for", obj_str))
#   path <- paste0(global_path, region, "/")
#   dge.path <- list.files(path, pattern = paste("*.gz", sep=""), full.names = TRUE)
#   dge <- loadSparseDge(dge.path)
#   obj <- CreateSeuratObject(counts = dge, project = obj_str)
#   cluster_assign      <- readRDS( list.files(path, pattern = paste("*.cluster.assign.RDS", sep=""), full.names = TRUE)[1] )
#   subcluster_assign   <- readRDS( list.files(path, pattern = paste("*.subcluster.assign.RDS", sep=""), full.names = TRUE) )
#   cluster_assign     <- cluster_assign[match(colnames(obj), names(cluster_assign))]
#   subcluster_assign  <- subcluster_assign[match(colnames(obj), names(subcluster_assign))]
#   obj$cluster_assign    <- cluster_assign
#   obj$subcluster_assign <- subcluster_assign
#   obj$cluster_assign    <- as.numeric(lapply(subcluster_assign, function(x) str_split(x, "-")[[1]][1]))
#   Idents(obj) <- obj$cluster_assign
#   cluster_avg <- AverageExpression(obj)[[1]]
# 
#   print(paste0("Finding pSI for ", obj_str))
#   results <- specificity.index(cluster_avg)
#   num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
#   for ( i in 0:num_clusters) {
#     cluster_deg <- all_deg$gene[which(all_deg$cluster == i)]
#     # cluster_deg <- convertMzebraGeneListToMouse(cluster_deg)# gene_names are converted to mouse
#     fisher_results <- fisher.iteration(results, cluster_deg, p.adjust = TRUE)
#     fisher_results$cluster <- i
#     fisher_results$region  <- region
#     fisher_results$region_cluster        <- 1:nrow(fisher_results)
#     fisher_results$region_cluster_label  <- sapply(levels(Idents(obj)), function(x) unique(region_ann$class[which(region_ann$cluster == x)])[1])
#     fisher_results$region_cluster_marker <- sapply(levels(Idents(obj)), function(x) unique(region_ann$class_marker[which(region_ann$cluster == x)])[1])
#     all_fisher <- rbind(all_fisher, fisher_results)
#   }
# 
# }
# all_fisher$q <- p.adjust(all_fisher[,2], method = "bonferroni")
# # all_fisher$brain_cluster <- rownames(all_fisher)
# write.table(all_fisher, "/nv/hp10/ggruenhagen3/scratch/brain/results/3_brain_region_all.tsv", sep="\t", quote = FALSE, row.names = FALSE)
# write.table(all_fisher[which(all_fisher$q < 0.05),], "/nv/hp10/ggruenhagen3/scratch/brain/results/3_brain_region_sig.tsv", sep="\t", quote = FALSE, row.names = FALSE)

## BY SUBCLUSTER ##
all_fisher <- data.frame()
for (i in 1:length(regions)) {
  region <- regions[i]
  obj_str <- region
  obj_abr <- regions_abr[i]
  region_ann <- ann[which(ann$tissue == obj_abr),]
  region_ann$cluster <- as.numeric(lapply(region_ann$subcluster, function(x) str_split(x, "-")[[1]][1]))
  print(paste("Loading data for", obj_str))
  path <- paste0(global_path, region, "/")
  dge.path <- list.files(path, pattern = paste("*.gz", sep=""), full.names = TRUE)
  dge <- loadSparseDge(dge.path)
  obj <- CreateSeuratObject(counts = dge, project = obj_str)
  cluster_assign      <- readRDS( list.files(path, pattern = paste("*.cluster.assign.RDS", sep=""), full.names = TRUE)[1] )
  subcluster_assign   <- readRDS( list.files(path, pattern = paste("*.subcluster.assign.RDS", sep=""), full.names = TRUE) )
  cluster_assign     <- cluster_assign[match(colnames(obj), names(cluster_assign))]
  subcluster_assign  <- subcluster_assign[match(colnames(obj), names(subcluster_assign))]
  obj$cluster_assign    <- cluster_assign
  obj$subcluster_assign <- subcluster_assign
  Idents(obj) <- subcluster_assign
  cluster_avg <- AverageExpression(obj)[[1]]

  print(paste0("Finding pSI for ", obj_str))
  results <- specificity.index(cluster_avg)
  num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
  for ( i in 0:num_clusters) {
    cluster_deg <- all_deg$gene[which(all_deg$cluster == i)]
    # cluster_deg <- convertMzebraGeneListToMouse(cluster_deg)# gene_names are converted to mouse
    fisher_results <- fisher.iteration(results, cluster_deg, p.adjust = TRUE)
    fisher_results$cluster <- i
    fisher_results$region  <- region
    fisher_results$region_subcluster  <- levels(Idents(obj))
    fisher_results$region_full_name   <- sapply(levels(Idents(obj)), function(x) region_ann$full_name[which(region_ann$subcluster == x)])
    fisher_results$region_common_name <- sapply(levels(Idents(obj)), function(x) region_ann$common_name[which(region_ann$subcluster == x)])
    all_fisher <- rbind(all_fisher, fisher_results)
  }

}
all_fisher$q <- p.adjust(all_fisher[,2], method = "bonferroni")
write.table(all_fisher, "/nv/hp10/ggruenhagen3/scratch/brain/results/3_brain_region_sub_all.tsv", sep="\t", quote = FALSE, row.names = FALSE)
write.table(all_fisher[which(all_fisher$q < 0.05),], "/nv/hp10/ggruenhagen3/scratch/brain/results/3_brain_region_sub_sig.tsv", sep="\t", quote = FALSE, row.names = FALSE)

## BY CELL_TYPE ##
obj_list <- list()
for (i in 1:length(cell_types)) {
  cell_type <- cell_types[i]
  print(paste("Loading data for", cell_type))
  path <- paste0(top_global_path, "cell_type/", cell_type, "/")
  dge.path <- list.files(path, pattern = paste("*.gz", sep=""), full.names = TRUE)
  dge <- loadSparseDge(dge.path)
  obj <- CreateSeuratObject(counts = dge, project = cell_type)
  obj$cell_type <- cell_type
  obj_list[[cell_type]] <- obj
  # subcluster_assign   <- readRDS( list.files(path, pattern = paste("*.subcluster.assign.RDS", sep=""), full.names = TRUE) )
}
big <- merge(x = obj_list[["astrocytes"]], y = c(obj_list[["endothelial"]], obj_list[["fibroblast"]], obj_list[["micro_macro"]], obj_list[["mural"]], obj_list[["oligo"]], obj_list[["poly"]]))
Idents(big) <- big$cell_type
cluster_avg <- AverageExpression(big)[[1]]
print(paste0("Finding pSI for all the cell types"))

all_fisher_cell <- data.frame()
results <- specificity.index(cluster_avg)
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
for ( i in 0:num_clusters) {
  cluster_deg <- all_deg$gene[which(all_deg$cluster == i)]
  # cluster_deg <- convertMzebraGeneListToMouse(cluster_deg)# gene_names are converted to mouse
  fisher_results <- fisher.iteration(results, cluster_deg, p.adjust = TRUE)
  fisher_results$cluster <- i
  fisher_results$cell_type <- levels(Idents(big))
  all_fisher_cell <- rbind(all_fisher_cell, fisher_results)
}
all_fisher_cell$q <- p.adjust(all_fisher_cell[,1], method = "bonferroni")
write.table(all_fisher, "/nv/hp10/ggruenhagen3/scratch/brain/results/9_brain_region_cell_all.tsv", sep="\t", quote = FALSE, row.names = FALSE)
write.table(all_fisher[which(all_fisher$q < 0.05),], "/nv/hp10/ggruenhagen3/scratch/brain/results/9_brain_region_cell_sig.tsv", sep="\t", quote = FALSE, row.names = FALSE)

## By MetaCells (Local) ##
regions_abr <- c("GP", "HC", "STR")
meta_data <- readRDS("C:/Users/miles/Downloads/brain/data/9_brain/metacells.BrainCellAtlas_Saunders_version_2018.04.01.rds")
meta <- CreateSeuratObject(meta_data)
meta <- meta[,which(startsWith(colnames(meta), regions_abr))]
ann <- readRDS("C:/Users/miles/Downloads/brain/data/9_brain/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS")
combined <- readRDS("C:/Users/miles/Downloads/brain/brain_scripts/brain_mz_shiny/data/B1C1C2MZ_combined_031020.rds")
combined_mouse <- convertToMouseObj(combined)
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
Idents(combined_mouse) <- "seurat_clusters"
cluster_avg <- AverageExpression(combined_mouse)[[1]]
results <- specificity.index(cluster_avg)

meta_deg <- FindAllMarkers(meta)

all_fisher <- data.frame()
for (i in 1:length(regions_abr)) {
  region <- regions_abr[i]
  meta_deg_region <- meta_deg$gene[which(meta_deg$cluster == region)]
  fisher_results <- fisher.iteration(results, meta_deg_region, p.adjust = TRUE)
  fisher_results$cluster <- 0:num_clusters
  fisher_results$region <- region
  all_fisher <- rbind(all_fisher, fisher_results)
}
sig <- all_fisher[which(all_fisher[,2] < 0.05),]

## By MetaCells (PACE) ##
meta_data <- readRDS(paste0(top_global_path, "metacells.BrainCellAtlas_Saunders_version_2018.04.01.RDS"))
meta <- CreateSeuratObject(meta_data)
meta <- meta[,which(startsWith(colnames(meta), regions_abr))]
meta$subregion <- colnames(meta)
Idents(meta) <- "subregion"
combined_mouse <- convertToMouseObj(combined)
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
Idents(combined_mouse) <- "seurat_clusters"
cluster_avg <- AverageExpression(combined_mouse)[[1]]
results <- specificity.index(cluster_avg)

meta_deg <- FindAllMarkers(meta)

all_fisher <- data.frame()
for (i in 1:ncol(meta)) {
  subregion <- colnames(meta)[i]
  meta_deg_subregion <- meta_deg$gene[which(meta_deg$cluster == subregion)]
  fisher_results <- fisher.iteration(results, meta_deg_subregion, p.adjust = TRUE)
  fisher_results$cluster <- 0:num_clusters
  fisher_results$subregion <- subregion
  all_fisher <- rbind(all_fisher, fisher_results)
}
nrow(all_fisher[which(all_fisher[,2] < 0.05),])
write.table(all_fisher, "/nv/hp10/ggruenhagen3/scratch/brain/results/3_brain_reciprocal.tsv", sep = "\t", quote = FALSE, row.names = FALSE)