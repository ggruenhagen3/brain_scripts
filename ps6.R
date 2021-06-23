# Load BB and Functions
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("SeuratObject")
bb = readRDS(paste0(rna_path, "data/bb_cc_04072021.RDS"))
Idents(bb) = bb$seurat_clusters
library(pacman)
p_unload(SeuratDisk)
p_unload(Seurat)
p_load(Seurat)
library("plyr")
library("parallel")

# 15 Cluster conversion
library("scales")
cols15 = gc.ramp <- hue_pal()(15)
convert15 = data.frame(old = 0:14, new = c("8_Ex", "9_Ex", "4_In", "15_In/Ex", "1_Astro/MG", "10_Ex", "5_In", "11_Ex", "6_In", "2_OPC/Oligo", "12_Ex", "13_Ex", "14_Ex", "3_Peri", "7_In"))
convert15$new = gsub("In", "GABA", convert15$new)
convert15$new = gsub("Ex", "Glut", convert15$new)
convert15$new.full = convert15$new
convert15 = separate(data = convert15, col = new, into = c("new.num", "new.junk"), sep = "_")
convert15 = convert15[order(as.numeric(convert15$new.num), decreasing = F),]
convert15$col = cols15

# 53 Cluster Conversion
convert53 = data.frame(old = 0:52, new = c("4.1_In", "10.1_Ex", "15.1_In/Ex", "9.1_Ex", "8.1_Ex", "1.1_Astro", "6_In", "5.1_In", "9.2_Ex", "8.2_Ex", "15.2_In", "11.1_Ex", "8.3_Ex", "8.4_Ex", "9.3_Ex", "4.2_In", "8.5_Ex", "5.2_In", "8.6_Ex", "8.7_Ex", "1.2_Astro", "4.3_In", "4.4_In", "9.4_Ex", "9.5_Ex", "8.8_Ex", "9.6_Ex", "4.5_In", "12_Ex", "8.9_Ex", "10.2_Ex", "2.1_OPC", "15.3_In", "11.2_Ex", "15.4_In", "4.6_In", "9.7_Ex", "13_Ex", "14_Ex", "4.7_In", "11.3_Ex", "9.8_Ex", "8-9_Ex", "15.5_In/Ex", "4.8_In", "1.3_MG", "2.2_Oligo", "15.6_Ex", "8.10_Ex", "8.11_Ex", "3_Peri", "15.7_Ex", "7_In"))
convert53$new = gsub("In", "GABA", convert53$new)
convert53$new = gsub("Ex", "Glut", convert53$new)
convert53$new2 = convert53$new
convert53 = separate(data = convert53, col = new2, into = c("new.id", "new.gaba"), sep = "_")
convert53 = separate(data = convert53, col = new.id, into = c("new.parent", "new.sub"), sep = "\\.")
convert53$new.parent[which(convert53$old == 42)] = 8
convert53$new.sub[which(convert53$old == 42)] = 12

# Load Gene Lists
ieg_genes = read.csv("~/scratch/brain/data/ieg_like_fos_egr1_npas4_detected_011521.csv", stringsAsFactors = F)[,1]
prog_genes = read.csv("~/scratch/brain/data/progenitor_sox2_nes_coexpress_051721.csv", stringsAsFactors = F)[,1]
neurogen_genes = read.csv("~/scratch/brain/data/neurogen_genes_final_050621.csv", stringsAsFactors = F)[,1]

# Calculate Score
mat = bb@assays$RNA@counts
mat[which(mat > 1)] = 1

# IEG
bb$ieg_score = colSums(mat[ieg_genes, ])
bb$ieg_score_norm = bb$ieg_score / bb$nFeature_RNA
# ieg_score_95 = colnames(bb)[which(bb$ieg_score >= quantile(bb$ieg_score, 0.95))]
ieg_score_95 = colnames(bb)[which(bb$ieg_score >= 2)]

# Function to find BHVE vs CTRL DEGs within a cluster, within high scoring cells
FindMarkersPair = function(pair, cluster_cells) {
  bb$pos_cluster = bb$pos & colnames(bb) %in% cluster_cells
  bb$pos_cluster_bhve = bb$pos_cluster
  bb$pos_cluster_bhve[which(bb$pos_cluster & bb$sample == paste0("b", pair))] = "BHVE"
  bb$pos_cluster_bhve[which(bb$pos_cluster & bb$sample == paste0("c", pair))] = "CTRL"
  
  Idents(bb) = bb$pos_cluster_bhve
  ieg_ps6 = FindMarkers(bb, ident.1 = "BHVE", ident.2 = "CTRL", min.pct = 0.01, logfc.threshold = 0)
  ieg_ps6$gene = rownames(ieg_ps6)
  ieg_ps6$n_bhve = length(which(bb$pos_cluster_bhve == "BHVE"))
  ieg_ps6$n_ctrl = length(which(bb$pos_cluster_bhve == "CTRL"))
  return(ieg_ps6)
}

bb$pos = F
bb$pos[ieg_score_95] = T

# IEG 15 cluster level
numCores = detectCores()
ieg_ps6_clust15_full = data.frame()
for (clust15 in c(1, 2, 3)) {
  print(paste0("IEG, Cluster15 = ", clust15))
  this_cluster_cells = colnames(bb)[which(bb$seuratclusters15 == clust15)]
  clust_all_pairs_list = mclapply(1:5, function(x) FindMarkersPair(x, this_cluster_cells), mc.cores = numCores)
  
  # Concatenate dataframes and format the data
  multi_inner = join_all(clust_all_pairs_list, by='gene', type='inner')
  rownames(multi_inner) = multi_inner$gene
  multi_inner$gene = NULL
  colnames(multi_inner) = c( paste0(colnames(multi_inner)[1:7], ".1"), paste0(colnames(multi_inner)[8:14], ".2"), paste0(colnames(multi_inner)[15:21], ".3"), paste0(colnames(multi_inner)[22:28], ".4"), paste0(colnames(multi_inner)[29:35], ".5"))
  multi_inner$gene = rownames(multi_inner)
  multi_inner$cluster = clust15
  multi_inner$cluster_new = convert15$new.full[match(multi_inner$cluster, convert15$old)]
  multi_inner$cat = "ieg"
  ieg_ps6_clust15_full = rbind(ieg_ps6_clust15_full, multi_inner)
}
write.csv(ieg_ps6_clust15_full, "~/scratch/brain/results/ieg_ps6_clust15_full_061521.csv")

# IEG 53 cluster level
numCores = detectCores()
ieg_ps6_clust53_full = data.frame()
for (clust53 in c(3, 14, 27, 2)) {
  print(paste0("IEG, Cluster53 = ", clust53))
  this_cluster_cells = colnames(bb)[which(bb$seuratclusters53 == clust53)]
  clust_all_pairs_list = mclapply(1:5, function(x) FindMarkersPair(x, this_cluster_cells), mc.cores = numCores)
  
  # Concatenate dataframes and format the data
  multi_inner = join_all(clust_all_pairs_list, by='gene', type='inner')
  rownames(multi_inner) = multi_inner$gene
  multi_inner$gene = NULL
  colnames(multi_inner) = c( paste0(colnames(multi_inner)[1:7], ".1"), paste0(colnames(multi_inner)[8:14], ".2"), paste0(colnames(multi_inner)[15:21], ".3"), paste0(colnames(multi_inner)[22:28], ".4"), paste0(colnames(multi_inner)[29:35], ".5"))
  multi_inner$gene = rownames(multi_inner)
  multi_inner$cluster = clust53
  multi_inner$cluster_new = convert53$new[match(multi_inner$cluster, convert53$old)]
  multi_inner$cat = "ieg"
  ieg_ps6_clust53_full = rbind(ieg_ps6_clust53_full, multi_inner)
}
write.csv(ieg_ps6_clust53_full, "~/scratch/brain/results/ieg_ps6_clust53_full_061521.csv")

# Progenitor Score
bb$prog_score = colSums(mat[prog_genes, ])
bb$prog_score_norm = bb$prog_score / bb$nFeature_RNA
prog_score_95 = colnames(bb)[which(bb$prog_score >= 10)]

bb$pos = F
bb$pos[prog_score_95] = T

# Prog 15 cluster level
numCores = detectCores()
prog_ps6_clust15_full = data.frame()
for (clust15 in c(0)) {
  print(paste0("Prog, Cluster15 = ", clust15))
  this_cluster_cells = colnames(bb)[which(bb$seuratclusters15 == clust15)]
  clust_all_pairs_list = mclapply(1:5, function(x) FindMarkersPair(x, this_cluster_cells), mc.cores = numCores)
  
  # Concatenate dataframes and format the data
  multi_inner = join_all(clust_all_pairs_list, by='gene', type='inner')
  rownames(multi_inner) = multi_inner$gene
  multi_inner$gene = NULL
  colnames(multi_inner) = c( paste0(colnames(multi_inner)[1:7], ".1"), paste0(colnames(multi_inner)[8:14], ".2"), paste0(colnames(multi_inner)[15:21], ".3"), paste0(colnames(multi_inner)[22:28], ".4"), paste0(colnames(multi_inner)[29:35], ".5"))
  multi_inner$gene = rownames(multi_inner)
  multi_inner$cluster = clust15
  multi_inner$cluster_new = convert15$new.full[match(multi_inner$cluster, convert15$old)]
  multi_inner$cat = "prog"
  prog_ps6_clust15_full = rbind(prog_ps6_clust15_full, multi_inner)
}
write.csv(prog_ps6_clust15_full, "~/scratch/brain/results/prog_ps6_clust15_full_061521.csv")

# Prog 53 cluster level
numCores = detectCores()
prog_ps6_clust53_full = data.frame()
for (clust53 in c(1, 27)) {
  print(paste0("Prog, Cluster53 = ", clust53))
  this_cluster_cells = colnames(bb)[which(bb$seuratclusters53 == clust53)]
  clust_all_pairs_list = mclapply(1:5, function(x) FindMarkersPair(x, this_cluster_cells), mc.cores = numCores)
  
  # Concatenate dataframes and format the data
  multi_inner = join_all(clust_all_pairs_list, by='gene', type='inner')
  rownames(multi_inner) = multi_inner$gene
  multi_inner$gene = NULL
  colnames(multi_inner) = c( paste0(colnames(multi_inner)[1:7], ".1"), paste0(colnames(multi_inner)[8:14], ".2"), paste0(colnames(multi_inner)[15:21], ".3"), paste0(colnames(multi_inner)[22:28], ".4"), paste0(colnames(multi_inner)[29:35], ".5"))
  multi_inner$gene = rownames(multi_inner)
  multi_inner$cluster = clust53
  multi_inner$cluster_new = convert53$new[match(multi_inner$cluster, convert53$old)]
  multi_inner$cat = "prog"
  prog_ps6_clust53_full = rbind(prog_ps6_clust53_full, multi_inner)
}
write.csv(prog_ps6_clust53_full, "~/scratch/brain/results/prog_ps6_clust53_full_065321.csv")

# Neurogenesis Score
bb$neuro_score = colSums(mat[neurogen_genes, ])
bb$neuro_score_norm = bb$neuro_score / bb$nFeature_RNA
neuro_score_95 = colnames(bb)[which(bb$neuro_score >= 1)]
bb$pos = F
bb$pos[neuro_score_95] = T


# Neuro 53 cluster level
numCores = detectCores()
neuro_ps6_clust53_full = data.frame()
for (clust53 in c(12)) {
  print(paste0("neuro, Cluster53 = ", clust53))
  this_cluster_cells = colnames(bb)[which(bb$seuratclusters53 == clust53)]
  clust_all_pairs_list = mclapply(1:5, function(x) FindMarkersPair(x, this_cluster_cells), mc.cores = numCores)
  
  # Concatenate dataframes and format the data
  multi_inner = join_all(clust_all_pairs_list, by='gene', type='inner')
  rownames(multi_inner) = multi_inner$gene
  multi_inner$gene = NULL
  colnames(multi_inner) = c( paste0(colnames(multi_inner)[1:7], ".1"), paste0(colnames(multi_inner)[8:14], ".2"), paste0(colnames(multi_inner)[15:21], ".3"), paste0(colnames(multi_inner)[22:28], ".4"), paste0(colnames(multi_inner)[29:35], ".5"))
  multi_inner$gene = rownames(multi_inner)
  multi_inner$cluster = clust53
  multi_inner$cluster_new = convert53$new[match(multi_inner$cluster, convert53$old)]
  multi_inner$cat = "neuro"
  neuro_ps6_clust53_full = rbind(neuro_ps6_clust53_full, multi_inner)
}
write.csv(neuro_ps6_clust53_full, "~/scratch/brain/results/neurogen_ps6_clust53_full_065321.csv")
