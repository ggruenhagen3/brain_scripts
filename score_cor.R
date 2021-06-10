# Load Genes
score_genes = read.csv("~/scratch/brain/data/ieg_like_fos_egr1_npas4_detected_011521.csv")[,1]

# Load BB
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("SeuratObject")
bb = readRDS(paste0(rna_path, "data/bb_cc_04072021.RDS"))
Idents(bb) = bb$seurat_clusters

# Deal with package issues
library(pacman)
p_unload(SeuratDisk)
p_unload(Seurat)
p_load(Seurat)

# Calculate Score
mat = bb@assays$RNA@counts
mat[which(mat > 1)] = 1
bb$score <- colSums(mat[score_genes, ]) / bb$nFeature_RNA

library(parallel)
bulk_res = unlist(mclapply(rownames(bb), function(x) score_cor(bb$score, x) , mc.cores = detectCores()))

score_cor = function(score, gene, cells.use = NULL) {
  if (length(cells.use) > 0) {
    return(cor(score, bb@assays$RNA@data[gene, cells.use]))
  } else {
    return(cor(score, bb@assays$RNA@data[gene, ]))
  }
}
