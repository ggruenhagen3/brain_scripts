#==================================================================================================
# Helper Functions ================================================================================
#==================================================================================================
# Single Run Function
combosRes = function(perm, cluster_level) {
  # Set the random samples
  bb$subsample = perm_labels[[perm]]
  
  # Find the mean ieg_like_score in each cluster
  # ptm <- proc.time()
  df = bb@meta.data[, c("ieg_like_score", "seuratclusters15", "subsample")]
  ieg_long = aggregate(ieg_like_score ~ seuratclusters15 + subsample, df, mean)  # finds mean per groups
  ieg_mat = acast(ieg_long, subsample ~ seuratclusters15, value.var = "ieg_like_score")
  # mean_time = proc.time() - ptm
  
  # Find the Correlation between Cluster Combos in Behave Samples
  # ptm2 <- proc.time()
  b_r_mat = cor(ieg_mat[which(startsWith(rownames(ieg_mat), "b")),], y = NULL, use = "pairwise.complete.obs") # use = "pairwise.complete.obs" removes the NA subsamples
  b_r_mat[upper.tri(b_r_mat, diag = T)] = NA
  # b_r_long = na.omit(melt(b_r_mat))
  # colnames(b_r_long) = c("cluster1", "cluster2", "r")
  
  # Find the Correlation between Cluster Combos in Control Samples
  c_r_mat = cor(ieg_mat[which(startsWith(rownames(ieg_mat), "c")),], y = NULL, use = "pairwise.complete.obs")
  c_r_mat[upper.tri(c_r_mat, diag = T)] = NA
  # c_r_long = na.omit(melt(c_r_mat))
  # colnames(c_r_long) = c("cluster1", "cluster2", "r")
  
  bvc_mat = b_r_mat - c_r_mat
  bvc_long = na.omit(melt(bvc_mat))
  bvc = bvc_long[,3]
  names(bvc) = paste0(bvc_long[,1], "_", bvc_long[,2])
  
  # end_time = proc.time()
  # print(paste0("Mean time: ", mean_time["elapsed"]))
  # print(paste0("R^2 time: ", end_time["elapsed"] - ptm2["elapsed"]))
  # print(paste0("Total time: ", end_time["elapsed"] - ptm["elapsed"]))
  
  return(bvc)
}

#==================================================================================================
# Main Body =======================================================================================
#==================================================================================================
# Load BB Dataset
# rna_path = "C:/Users/miles/Downloads/brain/"
rna_path = "~/scratch/brain/"
bb = readRDS(paste0(rna_path, "data/bb_subsample_02222021.RDS"))
source(paste0(rna_path, "brain_scripts/all_f.R"))
set.seed(156)

# Load IEG and IEG Like Genes
ieg_cons = c("LOC101487312", "egr1", "npas4", "jun", "homer1")
# ieg_like = read.csv("C:/Users/miles/Downloads/ieg_like_fos_egr1_npas4_detected_011521.csv", stringsAsFactors = F)[,2]
ieg_like = read.csv(paste0(rna_path, "/results/ieg_like_fos_egr1_npas4_detected_011521.csv"), stringsAsFactors = F)[,"ieg_like"]
ieg_like = c(ieg_like, "jun")
bb$ieg_like_score <- colSums(bb@assays$RNA@data[ieg_like,])

# Setup Permutations
n_perm = 100001
perm_labels = lapply(1:n_perm, function(x) sample(unname(as.vector(bb$subsample))))
bb$backup_subsample = bb$subsample

# Parallelize Finding Difference in Behave and Control for Cluster Combos
library("parallel")
numCores = detectCores()
cluster_level = "seuratclusters15"
# perm_labels = mclapply(1:n_perm, function(x) sample(unname(as.vector(bb$subsample))), mc.cores = numCores)
# length(which(sapply(1:n_perm, function(x) is.null(perm_labels[[x]]))))
all_combos = mclapply(1:n_perm, function(perm) combosRes(perm), mc.cores = numCores)

# For some reason, one of the 100k returned a vector of length 1377 instead of 1378.
# Double check that all returned vectors are the right size
bad_idx = mclapply(1:n_perm, function(perm) length(all_combos[[perm]]) != 105, mc.cores = numCores)
bad_idx = which(unlist(bad_idx))
all_combos[bad_idx] = NULL

# Merge List of Combos into One Dataframe
perm_bvc_df = data.frame(t(plyr::ldply(all_combos, rbind))) # merge by name
perm_bvc_df = perm_bvc_df[,1:n_perm]
colnames(perm_bvc_df) = 1:n_perm
combos = colsplit(rownames(perm_bvc_df), pattern = "\\_", names = c('cluster1', 'cluster2')) # split combos vector into two columns
perm_bvc_df = cbind(combos, perm_bvc_df)

# Save the Results
write.csv(perm_bvc_df, "~/scratch/brain/results/ieg_covar_c15_p100k_bvc.csv")

# # After permutations are done:
# # Load perm results
# perm_bvc_df = read.csv("~/scratch/brain/results/ieg_covar_c53_p100k_bvc.csv")
# rownames(perm_bvc_df) = perm_bvc_df$X
# perm_bvc_df$X = NULL
# perm_bvc_df$X100001 = NULL
# colnames(perm_bvc_df) = c("cluster1", "cluster2", 1:n_perm)
# 
# # Real Results
# perm_labels = list()
# perm_labels[[n_perm+1]] = bb$backup_subsample
# real_combos = combosRes(n_perm+1)

# Graph the results
# df_bvc_plot = data.frame(cluster1 = rep(perm_df_r_b$Cluster1, n_perm), cluster2 = rep(perm_df_r_b$Cluster2, n_perm),
#                          bvc = unlist( perm_df_r_b[,3:(n_perm+2)] - perm_df_r_c[,3:(n_perm+2)] ), isReal = rep(FALSE, n_perm))
# real_combos$bvc = as.numeric(as.vector(real_combos$r_behave)) - as.numeric(as.vector(real_combos$r_control))
# real_combos$abs_bvc = abs(real_combos$bvc)
# real_combos$isReal = TRUE
# df_bvc_plot = rbind(df_bvc_plot, real_combos[,c("cluster1", "cluster2", "bvc", "isReal")])
# df_bvc_plot$combo = paste0(df_bvc_plot$cluster1, "_", df_bvc_plot$cluster2)
# df_bvc_plot$abs_bvc = abs(df_bvc_plot$bvc)
# 
# png("~/scratch/brain/results/ieg_covar_c53_p100k_r_bvc.png", width = 2500, height = 400, res = 90)
# ggplot(df_bvc_plot, aes(combo, bvc, fill = isReal, color = isReal)) + geom_boxplot(alpha = 0.7)
# dev.off()
# png("~/scratch/brain/results/ieg_covar_c53_p100k_r_abs_bvc.png", width = 2500, height = 400, res = 90)
# ggplot(df_bvc_plot, aes(combo, abs_bvc, fill = isReal, color = isReal)) + geom_boxplot(alpha = 0.7)
# dev.off()

# library("parallel")
# numCores = detectCores()
# df_bvc_plot3 = perm_bvc_df
# df_bvc_plot3$bvc = 0
# df_bvc_plot3[names(real_combos),"bvc"] = real_combos
# perm_greater_boolean = df_bvc_plot3[,as.character(c(1:n_perm))] > df_bvc_plot3$bvc
# df_bvc_plot3$n_perm_greater = rowSums(perm_greater_boolean)
# 
# df_bvc_plot3$cluster1 = factor(df_bvc_plot3$cluster1, levels = 0:14)
# df_bvc_plot3$cluster2 = factor(df_bvc_plot3$cluster2, levels = 0:14)
# 
# perm_greater_boolean_abs = abs(df_bvc_plot3[,as.character(c(1:n_perm))]) > abs(df_bvc_plot3$bvc)
# df_bvc_plot3$abs_n_perm_greater = rowSums(perm_greater_boolean_abs)
# 
# png("~/scratch/brain/results/ieg_covar_c53_p100k_r_bvc_perm_greater_raw.png", width = 850, height = 800, res = 90)
# ggplot(df_bvc_plot3, aes(cluster1, cluster2, fill = n_perm_greater)) + geom_tile() + scale_fill_viridis(discrete=F, limits=c(0, n_perm), begin = 1, end = 0) + ggtitle("Behave - Control Correation")
# dev.off()
# 
# png("~/scratch/brain/results/ieg_covar_c53_p100k_r_abs_bvc_perm_greater.png", width = 850, height = 800, res = 90)
# ggplot(df_bvc_plot3, aes(cluster1, cluster2, fill = abs_n_perm_greater)) + geom_tile() + scale_fill_viridis(discrete=F, limits=c(0, n_perm), begin = 1, end = 0) + labs(fill ="n_perm_greater") + ggtitle("Absolute Value Behave - Control Correation")
# dev.off()
# 
# df_bvc_plot3$p = df_bvc_plot3$n_perm_greater / n_perm
# df_bvc_plot3$q = p.adjust(df_bvc_plot3$p, method = "BH")
#   
# df_bvc_plot3$abs_p = df_bvc_plot3$abs_n_perm_greater / n_perm
# df_bvc_plot3$abs_q = p.adjust(df_bvc_plot3$abs_p, method = "BH")
# length(which(df_bvc_plot3$abs_p < 0.05))
# length(which(df_bvc_plot3$abs_q < 0.05))
# # df_bvc_plot3[which(df_bvc_plot3$abs_q < 0.05),c(1,2,"abs_bvc", "abs_n_perm_greater", "abs_p", "abs_q")]
# 
# # P value per combo
# z_scores = t(scale(t(df_bvc_plot3[x, c("bvc", as.character(1:n_perm))]))) # scale combos per row 
# p_from_z = lapply(1:nrow(z_scores), function(x) 2*pnorm(-abs(z_scores[x, "bvc"])) )
# # p_from_z = lapply(1:ncol(z_scores_mat), function(x) 2*pnorm(-abs(z_scores_mat["real",x]), mean = mean(z_scores_mat[,x]), sd = sd(z_scores_mat[,x])) )
# p = unlist(p_from_z)
# names(p) = rownames(z_scores)
# q = p.adjust(p, method = "BH")
# names(q) = rownames(z_scores)
# length(which(p < 0.05))
# length(which(q < 0.05))
# q[which.min(q)]
