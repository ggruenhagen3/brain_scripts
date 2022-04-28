#*******************************************************************************
# Helper Functinos =============================================================
#*******************************************************************************
randomCombo = function(x) {
  set.seed(x)
  random_pcrc     = sample(not_pcrc_clust,     n_pcrc_clust)
  random_neurogen = sample(not_neurogen_clust, n_neurogen_clust)
  return(mean(cor_mat[random_pcrc, random_neurogen]))
}
createPermDataMat = function(x) {
  set.seed(x)
  this_data_mat = do.call('rbind', lapply(1:nrow(data_mat), function(y) sample(data_mat[y,])))
  rownames(this_data_mat) = rownames(data_mat)
  return(this_data_mat)
}
permModuleWoPermMat = function(x) {
  this_data_mat = perm_data_mats[[x]]
  this_cor_mat = cor(x = t(as.matrix(this_data_mat)))
  this_cor_mat = this_cor_mat[pcrc_clust, neurogen_clust]
  this_cor_mat[is.na(this_cor_mat)] = 0
  return(mean(this_cor_mat))
}
permModule = function(x) {
  set.seed(x)
  this_data_mat = do.call('rbind', lapply(1:nrow(data_mat), function(y) sample(data_mat[y,])))
  rownames(this_data_mat) = rownames(data_mat)
  this_cor_mat = cor(x = t(as.matrix(this_data_mat)))
  this_cor_mat = this_cor_mat[pcrc_clust, neurogen_clust]
  this_cor_mat[is.na(this_cor_mat)] = 0
  return(mean(this_cor_mat))
}
permModuleRGCCor = function(x) {
  set.seed(x)
  this_counts_mat = do.call('rbind', lapply(1:nrow(full_counts_mat), function(y) sample(full_counts_mat[y,])))
  rownames(this_counts_mat) = rownames(full_counts_mat)
  rgc_sub$perm_mod_score_pcrc = colSums(this_counts_mat[c(pcrc_clust), ] > 0)
  rgc_sub$perm_mod_score_neurogen = colSums(this_counts_mat[c(neurogen_clust), ] > 0)
  rgc_sub$perm_mod_score = colSums(this_counts_mat > 0)
  
  perm_mod_qui = cor(rgc_sub$perm_mod_score, rgc_sub$quiescent_score)
  perm_mod_cyc = cor(rgc_sub$perm_mod_score, rgc_sub$cycling_score)
  
  perm_mod_pcrc_qui = cor(rgc_sub$perm_mod_score_pcrc, rgc_sub$quiescent_score)
  perm_mod_pcrc_cyc = cor(rgc_sub$perm_mod_score_pcrc, rgc_sub$cycling_score)
  
  perm_mod_neurogen_qui = cor(rgc_sub$perm_mod_score_neurogen, rgc_sub$quiescent_score)
  perm_mod_neurogen_cyc = cor(rgc_sub$perm_mod_score_neurogen, rgc_sub$cycling_score)
  
  return(c(x, perm_mod_qui, perm_mod_cyc, perm_mod_pcrc_qui, perm_mod_pcrc_cyc, perm_mod_neurogen_qui, perm_mod_neurogen_cyc))
}
#*******************************************************************************
# Main =========================================================================
#*******************************************************************************
# Load Data
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("SeuratObject")
bb = readRDS(paste0(rna_path, "data/bb_demux_102021.rds"))
Idents(bb) = bb$seurat_clusters

print(paste0("Loading Data: ", format(Sys.time(), "%X ")))
neurogen = read.csv("~/scratch/brain/data/conserved_neurogenesis_positive88_zfish_mouse_cichlid.csv")[,3]
pcrc = read.csv("~/scratch/brain/fst/pc_20_rc_20_10kb_bins_25kb_genes_on_lg_11_peak_by_bin.csv")[,2]
neurogen = neurogen[which(!duplicated(neurogen))]
data_mat = bb@assays$RNA@data[c(pcrc, neurogen),]
full_counts_mat = bb@assays$RNA@counts[c(pcrc, neurogen),]
cor_mat = cor(x = t(as.matrix(data_mat)))
cor_mat = cor_mat[pcrc, neurogen]
cor_mat[is.na(cor_mat)] = 0

pcrc_clust     = c("cobl", "LOC101479283", "wdr73", "plekhg4b", "grik5", "LOC101476487", "LOC101476914", "ddr1", "LOC101477204", "plekhf2")
neurogen_clust = c("csf1r", "LOC101480727", "vegfa", "LOC101484715", "arhgef10", "stat3", "erbb2", "smo", "epha3", "LOC101469419", "LOC101487687", "boc", "pax6", "metrn", "LOC101469466")
n_pcrc_clust     = length(pcrc_clust)
n_neurogen_clust = length(neurogen_clust)
not_pcrc_clust     = rownames(cor_mat)[which(!rownames(cor_mat) %in% pcrc_clust)]
not_neurogen_clust = colnames(cor_mat)[which(!colnames(cor_mat) %in% neurogen_clust)]
real_mean = mean(cor_mat[pcrc_clust, neurogen_clust])

# 1
# Do other combos of PCRC and Neurogen of the same size as the real module (10x15),
# have as strong of an average correlation coefficient?
print(paste0("Starting Random Combos: ", format(Sys.time(), "%X ")))
library(parallel)
nperm = 10000
random_means = unlist(mclapply(1:nperm, function(x) randomCombo(x), mc.cores = detectCores()))
print(paste0("Number of Random Combo Average Correlation Greater Than Real: ", length(which(random_means >= real_mean))))
print(paste0("Random Combos Done: ", format(Sys.time(), "%X ")))

# 2
# Is the average correlation coefficient of the module greater than permutations of the module?
# Where a permutation is when the expression values for each are shuffled separately.
## perm_data_mats = mclapply(1:nperm, function(x) createPermDataMat(x), mc.cores = detectCores())
# print("")
# print(paste0("Starting Permutations of the Module: ", format(Sys.time(), "%X ")))
# perm_mod_means = unlist(mclapply(1:nperm, function(x) permModule(x), mc.cores = detectCores()))
# print(paste0("Number of Permutations of the Module Greater Than Real: ", length(which(perm_mod_means >= real_mean))))
# print(paste0("Permutations of the Module Done: ", format(Sys.time(), "%X ")))
# 
# print("")
# print(paste0("Writing Output to csv: ", format(Sys.time(), "%X ")))
# write.csv(data.frame(random_combo_mean = random_means, perm_mod_mean = perm_mod_means, real_mean = real_mean), "~/scratch/brain/results/pcrc_neurogen_perm_042822.csv")
# print(paste0("All Done: ", format(Sys.time(), "%X ")))

print("")
print(paste0("Loading RGC Object: ", format(Sys.time(), "%X ")))
rgc_sub = readRDS("~/scratch/brain/data//rgc_subclusters_reclustered_q_c_nb_scores.rds")

# Module Relationship to Quiescent and Cycling RGCs
# 3a. Is the correlation between module score and quiescent score in RGCs greater than permutations?
# 3b. Is the correlation between module score for CDGs and quiescent score in RGCs greater than permutations?
# 3c. Is the correlation between module score for pNGs and quiescent score in RGCs greater than permutations?
# 3d. Is the correlation between module score and cycling score in RGCs more negative than permutations?
# 3e. Is the correlation between module score for CDGs and cycling score in RGCs more negative than permutations?
# 3f. Is the correlation between module score for pNGs and cycling score in RGCs more negative than permutations?
rgc_sub$mod_score_pcrc = colSums(rgc_sub@assays$RNA@counts[c(pcrc_clust), ] > 0)
rgc_sub$mod_score_neurogen = colSums(rgc_sub@assays$RNA@counts[c(neurogen_clust), ] > 0)
rgc_sub$mod_score = colSums(rgc_sub@assays$RNA@counts[c(pcrc_clust, neurogen_clust), ] > 0)

real_mod_qui = cor(rgc_sub$mod_score, rgc_sub$quiescent_score)
real_mod_cyc = cor(rgc_sub$mod_score, rgc_sub$cycling_score)
real_mod_pcrc_qui = cor(rgc_sub$mod_score_pcrc, rgc_sub$quiescent_score)
real_mod_pcrc_cyc = cor(rgc_sub$mod_score_pcrc, rgc_sub$cycling_score)
real_mod_neurogen_qui = cor(rgc_sub$mod_score_neurogen, rgc_sub$quiescent_score)
real_mod_neurogen_cyc = cor(rgc_sub$mod_score_neurogen, rgc_sub$cycling_score)

print("")
print(paste0("Starting Permutations Examing the Relationship b/w the Module and Quiescent/Cycling RGCs: ", format(Sys.time(), "%X ")))
perm_mod_rgc = mclapply(1:nperm, function(x) permModuleRGCCor(x), mc.cores = detectCores())
print(paste0("Permutations of the Module Done: ", format(Sys.time(), "%X ")))
perm_mod_rgc_df = as.data.frame(t(as.data.frame(perm_mod_rgc)))
colnames(perm_mod_rgc_df) = c("perm_num", "perm_mod_qui", "perm_mod_cyc", "perm_mod_pcrc_qui", "perm_mod_pcrc_cyc", "perm_mod_neurogen_qui", "perm_mod_neurogen_cyc")
rownames(perm_mod_rgc_df) = perm_mod_rgc_df$perm_num
print(paste0("Number of Permutations of the Module Stronger Than Real - Mod Quiescent: ", length(which(perm_mod_rgc_df$perm_mod_qui          >= real_mod_qui))))
print(paste0("Number of Permutations of the Module Stronger Than Real - Mod Cycling: ",   length(which(perm_mod_rgc_df$perm_mod_cyc          <= real_mod_cyc))))
print(paste0("Number of Permutations of the Module Stronger Than Real - Mod Quiescent: ", length(which(perm_mod_rgc_df$perm_mod_pcrc_qui     >= real_mod_pcrc_qui))))
print(paste0("Number of Permutations of the Module Stronger Than Real - Mod Quiescent: ", length(which(perm_mod_rgc_df$perm_mod_pcrc_cyc     <= real_mod_pcrc_cyc))))
print(paste0("Number of Permutations of the Module Stronger Than Real - Mod Quiescent: ", length(which(perm_mod_rgc_df$perm_mod_neurogen_qui >= real_mod_neurogen_qui))))
print(paste0("Number of Permutations of the Module Stronger Than Real - Mod Quiescent: ", length(which(perm_mod_rgc_df$perm_mod_neurogen_cyc <= real_mod_neurogen_cyc))))

write.csv(perm_mod_rgc_df, "~/scratch/brain/results/pcrc_neurogen_perm_mod_rgc_042822.csv")

