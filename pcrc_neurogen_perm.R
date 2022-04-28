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
#*******************************************************************************
# Main =========================================================================
#*******************************************************************************
# Load Data
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("SeuratObject")
bb = readRDS(paste0(rna_path, "data/bb_demux_102021.rds"))
Idents(bb) = bb$seurat_clusters

neurogen = read.csv("~/scratch/brain/data/conserved_neurogenesis_positive88_zfish_mouse_cichlid.csv")[,3]
pcrc = read.csv("~/scratch/brain/fst/pc_20_rc_20_10kb_bins_25kb_genes_on_lg_11_peak_by_bin.csv")[,2]
neurogen = neurogen[which(!duplicated(neurogen))]
data_mat = bb@assays$RNA@data[c(pcrc, neurogen),]
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
cat(paste0("Starting Random Combos: ", format(Sys.time(), "%X ")))
library(parallel)
nperm = 10000
random_means = unlist(mclapply(1:nperm, function(x) randomCombo(x), mc.cores = detectCores()))
cat(paste0("Number of Random Combo Average Correlation Greater Than Real: ", length(which(random_means >= real_mean))))
cat(paste0("Random Combos Done: ", format(Sys.time(), "%X ")))

# 2
# Is the average correlation coefficient of the module greater than permutations of the module?
# Where a permutation is when the expression values for each are shuffled separately.
# perm_data_mats = mclapply(1:nperm, function(x) createPermDataMat(x), mc.cores = detectCores())
cat("")
cat(paste0("Starting Permutations of the Module: ", format(Sys.time(), "%X ")))
perm_mod_means = unlist(mclapply(1:nperm, function(x) permModule(x), mc.cores = detectCores()))
cat(paste0("Number of Permutations of the Module Greater Than Real: ", length(which(perm_mod_means >= real_mean))))
cat(paste0("Permutations of the Module Done: ", format(Sys.time(), "%X ")))

cat("")
cat(paste0("Writing Output to csv: ", format(Sys.time(), "%X ")))
write.csv(data.frame(random_combo_mean = random_means, perm_mod_mean = perm_mod_means, real_mean = real_mean), "~/scratch/brain/results/pcrc_neurogen_perm_042822.csv")
cat(paste0("All Done: ", format(Sys.time(), "%X ")))
