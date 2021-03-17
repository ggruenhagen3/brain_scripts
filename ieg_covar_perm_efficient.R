#==================================================================================================
# Helper Functions ================================================================================
#==================================================================================================
# Single Run Function
combosRes = function(perm, cluster_level) {
  # Set the random samples
  bb$subsample = perm_labels[[perm]]
  
  # Find the mean ieg_like_score in each cluster
  # ptm <- proc.time()
  df = bb@meta.data[, c("ieg_like_score", "seuratclusters53", "subsample")]
  ieg_long = aggregate(ieg_like_score ~ seuratclusters53 + subsample, df, mean)  # finds mean per groups
  ieg_mat = acast(ieg_long, subsample ~ seuratclusters53, value.var = "ieg_like_score")
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
set.seed(155)

# Load IEG and IEG Like Genes
ieg_cons = c("LOC101487312", "egr1", "npas4", "jun", "homer1")
# ieg_like = read.csv("C:/Users/miles/Downloads/ieg_like_fos_egr1_npas4_detected_011521.csv", stringsAsFactors = F)[,2]
ieg_like = read.csv(paste0(rna_path, "/results/ieg_like_fos_egr1_npas4_detected_011521.csv"), stringsAsFactors = F)[,"ieg_like"]
ieg_like = c(ieg_like, "jun")
bb$ieg_like_score <- colSums(bb@assays$RNA@data[ieg_like,])

# Setup Permutations
n_perm = 100000
perm_labels = lapply(1:n_perm, function(x) sample(unname(as.vector(bb$subsample))))
bb$backup_subsample = bb$subsample

# Parallelize Finding Difference in Behave and Control for Cluster Combos
library("parallel")
numCores = detectCores()
cluster_level = "seuratclusters53"
# perm_labels = mclapply(1:n_perm, function(x) sample(unname(as.vector(bb$subsample))), mc.cores = numCores)
# length(which(sapply(1:n_perm, function(x) is.null(perm_labels[[x]]))))
all_combos = mclapply(1:n_perm, function(perm) combosRes(perm), mc.cores = numCores)

# Generate the Cluster Combos Once Outside of the For Loop
combos = data.frame()
num_clusters = max(as.numeric(as.vector(levels(bb@meta.data[,c(cluster_level)]))))
for (cluster2 in 0:(num_clusters-1)) {
  for (cluster1 in (cluster2+1):num_clusters) {
    combos = rbind(combos, t(c(cluster1, cluster2)))
  }
}

# Reformat Data
perm_bvc_mat = matrix(unlist(all_combos), nrow = nrow(combos), dimnames = list(1:nrow(combos), 1:n_perm))
perm_bvc_df = data.frame(perm_bvc_mat)
perm_bvc_df = cbind(combos, perm_bvc_df)

# Save the Results
write.csv(perm_bvc_df, "~/scratch/brain/results/ieg_covar_c53_p100k_bvc.csv")

# # After permutations are done:
# # Load perm results
# perm_df = read.csv("~/scratch/brain/results/ieg_covar_c53_p1000_summary.csv")
# perm_df_r_b = read.csv("~/scratch/brain/results/ieg_covar_c53_p1000_r_b.csv")
# perm_df_r_c = read.csv("~/scratch/brain/results/ieg_covar_c53_p1000_r_c.csv")
# perm_df_p_b = read.csv("~/scratch/brain/results/ieg_covar_c53_p1000_p_b.csv")
# perm_df_p_c = read.csv("~/scratch/brain/results/ieg_covar_c53_p1000_p_c.csv")
# perm_df_p_bvc = read.csv("~/scratch/brain/results/ieg_covar_c53_p1000_p_bvc.csv")
# perm_df$X = perm_df_r_b$X = perm_df_r_c$X = perm_df_p_b$X = perm_df_p_c$X = perm_df_p_bvc$X = NULL
# 
# # Real Results
# perm_labels[[n_perm+1]] = bb$backup_subsample
# real_combos = combosRes(n_perm+1)
# 
# # Graph the results
# df_bvc_plot = data.frame(cluster1 = rep(perm_df_r_b$Cluster1, n_perm), cluster2 = rep(perm_df_r_b$Cluster2, n_perm),
#                          bvc = unlist( perm_df_r_b[,3:(n_perm+2)] - perm_df_r_c[,3:(n_perm+2)] ), isReal = rep(FALSE, n_perm))
# real_combos$bvc = as.numeric(as.vector(real_combos$r_behave)) - as.numeric(as.vector(real_combos$r_control))
# real_combos$abs_bvc = abs(real_combos$bvc)
# real_combos$isReal = TRUE
# df_bvc_plot = rbind(df_bvc_plot, real_combos[,c("cluster1", "cluster2", "bvc", "isReal")])
# df_bvc_plot$combo = paste0(df_bvc_plot$cluster1, "_", df_bvc_plot$cluster2)
# df_bvc_plot$abs_bvc = abs(df_bvc_plot$bvc)
# 
# png("~/scratch/brain/results/ieg_covar_c53_p1000_r_bvc.png", width = 2500, height = 400, res = 90)
# ggplot(df_bvc_plot, aes(combo, bvc, fill = isReal, color = isReal)) + geom_boxplot(alpha = 0.7)
# dev.off()
# png("~/scratch/brain/results/ieg_covar_c53_p1000_r_abs_bvc.png", width = 2500, height = 400, res = 90)
# ggplot(df_bvc_plot, aes(combo, abs_bvc, fill = isReal, color = isReal)) + geom_boxplot(alpha = 0.7)
# dev.off()
# 
# df_bvc_plot3 = real_combos[,c("cluster1", "cluster2", "bvc")]
# df_bvc_plot3 = cbind(df_bvc_plot3, perm_df_r_b[,3:(n_perm+2)] - perm_df_r_c[,3:(n_perm+2)])
# colnames(df_bvc_plot3) = c("cluster1", "cluster2", "bvc", 1:n_perm)
# df_bvc_plot3$n_perm_greater = sapply(1:nrow(df_bvc_plot3), function(x) length(which(df_bvc_plot3[x,as.character(c(1:n_perm))] > df_bvc_plot3[x,c("bvc")])) )
# # df_bvc_plot3 = df_bvc_plot3[order(as.numeric(as.vector(df_bvc_plot3$cluster1)), as.numeric(as.vector(df_bvc_plot3$cluster2))), ]
# df_bvc_plot3 = df_bvc_plot3[order(as.numeric(as.vector(df_bvc_plot3$cluster2))), ]
# 
# png("~/scratch/brain/results/ieg_covar_c53_p1000_r_bvc_perm_greater_raw.png", width = 850, height = 800, res = 90)
# ggplot(df_bvc_plot3, aes(cluster1, cluster2, fill = n_perm_greater)) + geom_tile() + scale_fill_viridis(discrete=F, limits=c(0, max(df_bvc_plot2$n_perm_greater)), begin = 1, end = 0)
# dev.off()
# 
# df_bvc_plot2 = real_combos[,c("cluster1", "cluster2", "abs_bvc")]
# df_bvc_plot2 = cbind(df_bvc_plot2, abs(perm_df_r_b[,3:(n_perm+2)] - perm_df_r_c[,3:(n_perm+2)]))
# colnames(df_bvc_plot2) = c("cluster1", "cluster2", "abs_bvc", 1:n_perm)
# df_bvc_plot2$n_perm_greater = sapply(1:nrow(df_bvc_plot2), function(x) length(which(df_bvc_plot2[x,as.character(c(1:n_perm))] > df_bvc_plot2[x,c("abs_bvc")])) )
# # df_bvc_plot2 = df_bvc_plot2[order(as.numeric(as.vector(df_bvc_plot2$cluster1)), as.numeric(as.vector(df_bvc_plot2$cluster2))), ]
# df_bvc_plot2 = df_bvc_plot2[order(as.numeric(as.vector(df_bvc_plot2$cluster2))), ]
# 
# png("~/scratch/brain/results/ieg_covar_c53_p1000_r_bvc_perm_greater.png", width = 850, height = 800, res = 90)
# ggplot(df_bvc_plot2, aes(cluster1, cluster2, fill = n_perm_greater)) + geom_tile() + scale_fill_viridis(discrete=F, limits=c(0, max(df_bvc_plot2$n_perm_greater)), begin = 1, end = 0)
# dev.off()
# 
# # P value per combo
# z_scores = lapply(1:nrow(df_bvc_plot3), function(x) scale(as.numeric(as.vector(df_bvc_plot3[x, c("bvc", as.character(1:n_perm))]))) )
# # z_scores = lapply(1:nrow(df_bvc_plot3), function(x) scale(as.numeric(as.vector(df_bvc_plot3[x, c("bvc", paste0("X",1:n_perm))]))) )
# z_scores_mat = as.matrix(data.frame(z_scores))
# rownames(z_scores_mat) = c("real", as.character(1:n_perm))
# colnames(z_scores_mat) = paste0(df_bvc_plot3$cluster1, "_", df_bvc_plot3$cluster2)
# 
# p_from_z = lapply(1:ncol(z_scores_mat), function(x) 2*pnorm(-abs(z_scores_mat["real",x])) )
# # p_from_z = lapply(1:ncol(z_scores_mat), function(x) 2*pnorm(-abs(z_scores_mat["real",x]), mean = mean(z_scores_mat[,x]), sd = sd(z_scores_mat[,x])) )
# p = unlist(p_from_z)
# names(p) = colnames(z_scores_mat)
# q = p.adjust(p, method = "BH")
# names(q) = colnames(z_scores_mat)
# length(which(p < 0.05))
# length(which(q < 0.05))
# q[which.min(q)]
