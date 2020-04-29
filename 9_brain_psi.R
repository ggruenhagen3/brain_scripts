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
global_path <- "/nv/hp10/ggruenhagen3/scratch/brain/data/9_brain/region/"
regions <- dir(global_path, pattern =paste("*", sep=""))
regions <- c("cerebellum", "ento_peduncular", "globulus_pallidus", "hippocampus", "nigra", "striatum", "thalamus")

# print("Finding pSI for our object")
combined <- readRDS("/nv/hp10/ggruenhagen3/scratch/brain/data/B1C1C2MZ_combined_031020.rds")
# our_psi <- AverageExpression(combined, slot = counts)
all_deg <- read.csv("/nv/hp10/ggruenhagen3/scratch/brain/data/all_markers_B1C1C2MZ_042820.csv", stringsAsFactors = FALSE)
all_deg <- convertMzebraDFToMouse(all_deg, 8)

all_fisher <- data.frame()
for (region in regions) {
  obj_str <- region
  print(paste("Loading data for", obj_str))
  path <- paste0(global_path, region, "/")
  dge.path <- list.files(path, pattern = paste("*.gz", sep=""), full.names = TRUE)
  dge <- loadSparseDge(dge.path)
  cluster_assign_orig <- readRDS( list.files(path, pattern = paste("*.cluster.assign.RDS", sep=""), full.names = TRUE)[1] )
  subcluster_assign   <- readRDS( list.files(path, pattern = paste("*.subcluster.assign.RDS", sep=""), full.names = TRUE) )
  obj <- CreateSeuratObject(counts = dge, project = obj_str)
  cluster_assign <- cluster_assign_orig[match(colnames(obj), names(cluster_assign_orig))]
  obj$cluster_assign <- cluster_assign
  Idents(obj) <- cluster_assign
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
    all_fisher <- rbind(all_fisher, fisher_results)
  }
  
}
all_fisher$q <- p.adjust(all_fisher[,1], method = "bonferroni")
all_fisher$brain_cluster <- rownames(all_fisher)
write.table(all_fisher, "/nv/hp10/ggruenhagen3/scratch/brain/results/9_brain_region_all.tsv", sep="\t", quote = FALSE, row.names = FALSE)
write.table(all_fisher[which(all_fisher$q < 0.05),], "/nv/hp10/ggruenhagen3/scratch/brain/results/9_brain_region_sig.tsv", sep="\t", quote = FALSE, row.names = FALSE)