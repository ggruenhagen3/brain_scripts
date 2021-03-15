#==================================================================================================
# Helper Functions ================================================================================
#==================================================================================================
# Single Run Function
combosRes = function(perm) {
  # Set the random samples
  bb$subsample = perm_labels[[perm]]
  
  # Find the mean ieg_like_score in each cluster
  cluster15_ieg = data.frame()
  cluster_level = "seuratclusters53"
  for (cluster in levels(bb@meta.data[,c(cluster_level)])) {
    cluster_cells = colnames(bb)[which(bb@meta.data[,c(cluster_level)] == cluster)]
    for (subsample in sort(unique(bb$subsample))) {
      subsample_cells = colnames(bb)[which(bb$subsample == subsample)]
      this_cells = cluster_cells[which(cluster_cells %in% subsample_cells)]
      num_cells = length(this_cells)
      ieg_like_mean = mean(bb$ieg_like_score[this_cells])
      cluster15_ieg = rbind(cluster15_ieg, t(c(cluster, subsample, num_cells, ieg_like_mean)))
    }
  }
  colnames(cluster15_ieg) = c("cluster", "sample", "num_cells", "mean_ieg_like_score")
  cluster15_ieg$mean_ieg_like_score = as.numeric(as.vector(cluster15_ieg$mean_ieg_like_score))
  cluster15_ieg$num_cells = as.numeric(as.vector(cluster15_ieg$num_cells))
  
  # Find the Correlation and R^2 of the cluster combos
  cluster15_ieg_combos = data.frame()
  for (cluster1 in levels(bb@meta.data[,c(cluster_level)])) {
    for (cluster2 in levels(bb@meta.data[,c(cluster_level)])) {
      cluster1 = as.numeric(cluster1)
      cluster2 = as.numeric(cluster2)
      if (cluster1 > cluster2) {
        this_df = cbind( cluster15_ieg[which(cluster15_ieg$cluster == cluster1),], cluster15_ieg[which(cluster15_ieg$cluster == cluster2),] )
        colnames(this_df) = c("cluster1", "sample", "num_cells_1", "mean_ieg_like_score_1", "cluster2", "sample2", "num_cells_2", "mean_ieg_like_score_2")
        this_df = this_df[which( this_df$num_cells_1 > 0 & this_df$num_cells_2 > 0 ),]
        this_df$sample = factor(this_df$sample, levels = unique(sort(bb$subsample)))
        this_r2 = cor(this_df$mean_ieg_like_score_1, this_df$mean_ieg_like_score_2) ^ 2
        this_df$isBehave = startsWith(as.character(as.vector(this_df$sample)), "b")
        
        # Behave
        this_df_b = this_df[which( this_df$isBehave ),]
        this_r_b = cor(this_df_b$mean_ieg_like_score_1, this_df_b$mean_ieg_like_score_2)
        this_r2_b = this_r_b ^ 2
        num_b = nrow(this_df_b)
        
        # Control
        this_df_c = this_df[which( ! this_df$isBehave ),]
        this_r_c = cor(this_df_c$mean_ieg_like_score_1, this_df_c$mean_ieg_like_score_2)
        this_r2_c = this_r_c ^ 2
        num_c = nrow(this_df_c)
        
        cluster15_ieg_combos = rbind(cluster15_ieg_combos, t(c(cluster1, cluster2, this_r2, this_r2_b, this_r2_c, this_r_b, this_r_c, num_b, num_c)))
      }
    }
  }
  colnames(cluster15_ieg_combos) = c("cluster1", "cluster2", "this_r2", "r2_behave", "r2_control", "r_behave", "r_control", "num_b", "num_c")
  cluster15_ieg_combos[c("this_r2", "r2_behave", "r2_control", "r_behave", "r_control", "num_b", "num_c")] <- sapply(cluster15_ieg_combos[c("this_r2", "r2_behave", "r2_control", "r_behave", "r_control", "num_b", "num_c")],as.vector)
  cluster15_ieg_combos[c("this_r2", "r2_behave", "r2_control", "r_behave", "r_control", "num_b", "num_c")] <- sapply(cluster15_ieg_combos[c("this_r2", "r2_behave", "r2_control", "r_behave", "r_control", "num_b", "num_c")],as.numeric)
  cluster15_ieg_combos$cluster1 = factor(cluster15_ieg_combos$cluster1, levels = levels(bb@meta.data[,c(cluster_level)]))
  cluster15_ieg_combos$cluster2 = factor(cluster15_ieg_combos$cluster2, levels = levels(bb@meta.data[,c(cluster_level)]))
  
  return(cluster15_ieg_combos)
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
bb$ieg_like_score = colSums(bb@assays$RNA@data[ieg_like,])

# Setup Permutations
n_perm = 1000
perm_labels = lapply(1:n_perm, function(x) sample(unname(as.vector(bb$subsample))))
bb$backup_subsample = bb$subsample
perm_df_r_b = data.frame()
perm_df_r_c = data.frame()
perm_df_p_b = data.frame()
perm_df_p_c = data.frame()
perm_df = data.frame()

# Parallelize Finding Cluster Combos
library("parallel")
numCores = detectCores()
all_combos = mclapply(1:n_perm, function(perm) combosRes(perm), mc.cores = numCores)

for (perm in 1:n_perm) {
  cat(paste0(perm, "."))
  
  # Load Parallelized Results
  cluster15_ieg_combos = all_combos[[perm]]
  
  # Find Significant Correlations
  my_cor_t = function(r, n) (r * sqrt(n - 2))/sqrt(1 - r**2)
  my_cor_p = function(t, n) 2*pt(-abs(t), df=n-2)
  
  # Convert Behave correlations to p values
  r_behave_mat = acast(cluster15_ieg_combos, cluster1 ~ cluster2, value.var="r_behave")
  n_behave_mat = acast(cluster15_ieg_combos, cluster1 ~ cluster2, value.var="num_b")
  t_behave_mat = my_cor_t(r_behave_mat, n_behave_mat)
  p_behave_mat = my_cor_p(t_behave_mat, n_behave_mat)
  
  # Convert Control correlations to p values
  r_control_mat = acast(cluster15_ieg_combos, cluster1 ~ cluster2, value.var="r_control")
  n_control_mat = acast(cluster15_ieg_combos, cluster1 ~ cluster2, value.var="num_c")
  t_control_mat = my_cor_t(r_control_mat, n_control_mat)
  p_control_mat = my_cor_p(t_control_mat, n_control_mat)
  
  # Find Number of Significant Behave Correlations
  # p_behave_mat[upper.tri(p_behave_mat)] = NA
  p_behave_df = na.omit(melt(p_behave_mat))
  colnames(p_behave_df) = c("Cluster1", "Cluster2", "p")
  p_behave_df$bh = p.adjust(p_behave_df$p, method = "BH")
  num_p_sig_b = length(which(p_behave_df$p < 0.05))
  num_q_sig_b = length(which(p_behave_df$bh < 0.05))
  
  # Find Number of Significant Control Correlations
  # p_control_mat[upper.tri(p_control_mat)] = NA
  p_control_df = na.omit(melt(p_control_mat))
  colnames(p_control_df) = c("Cluster1", "Cluster2", "p")
  p_control_df$bh = p.adjust(p_control_df$p, method = "BH")
  num_p_sig_c = length(which(p_control_df$p < 0.05))
  num_q_sig_c = length(which(p_control_df$bh < 0.05))
  
  # P-value of Behave vs Control Correlation Differences - Two-tailed
  p_bvc_mat = r_to_p(r_behave_mat, r_control_mat, n_behave_mat, n_control_mat)
  p_bvc_mat[upper.tri(p_bvc_mat)] = NA
  p_bvc_df = na.omit(melt(p_bvc_mat))
  colnames(p_bvc_df) = c("Cluster1", "Cluster2", "p")
  p_bvc_df$bh = p.adjust(p_bvc_df$p, method = "BH")
  num_p_sig_bvc = length(which(p_bvc_df$p < 0.05))
  num_q_sig_bvc = length(which(p_bvc_df$bh < 0.05))
  
  # Store permutation data
  if (perm == 1) {
    perm_df_r_b = cluster15_ieg_combos[,c("cluster1", "cluster2")]
    perm_df_r_c = cluster15_ieg_combos[,c("cluster1", "cluster2")]
    perm_df_p_b = p_behave_df[,c("Cluster1", "Cluster2")]
    perm_df_p_c = p_control_df[,c("Cluster1", "Cluster2")]
    perm_df_p_bvc = p_bvc_df[,c("Cluster1", "Cluster2")]
  }
  perm_df_r_b = cbind(perm_df_r_b, cluster15_ieg_combos[,c("r_behave")])
  perm_df_r_c = cbind(perm_df_r_c, cluster15_ieg_combos[,c("r_control")])
  perm_df_p_b = cbind(perm_df_p_b, p_behave_df[,c("bh")])
  perm_df_p_c = cbind(perm_df_p_c, p_control_df[,c("bh")])
  perm_df_p_bvc = cbind(perm_df_p_bvc, p_bvc_df[,c("bh")])
  perm_df = rbind(perm_df, t(c(perm, num_p_sig_bvc, num_q_sig_bvc, num_p_sig_b, num_q_sig_b, num_p_sig_c, num_q_sig_c)))
}
colnames(perm_df) = c("perm", "num_p_sig_bvc", "num_q_sig_bvc", "num_p_sig_b", "num_q_sig_b", "num_p_sig_c", "num_q_sig_c")
colnames(perm_df_r_b) = colnames(perm_df_r_c) = colnames(perm_df_p_b) = colnames(perm_df_p_c) = c("Cluster1", "Cluster2", 1:n_perm)

# Save the Results
write.csv(perm_df, "~/scratch/brain/results/ieg_covar_c53_p1000_summary.csv")
write.csv(perm_df_r_b, "~/scratch/brain/results/ieg_covar_c53_p1000_r_b.csv")
write.csv(perm_df_r_c, "~/scratch/brain/results/ieg_covar_c53_p1000_r_c.csv")
write.csv(perm_df_p_b, "~/scratch/brain/results/ieg_covar_c53_p1000_p_b.csv")
write.csv(perm_df_p_c, "~/scratch/brain/results/ieg_covar_c53_p1000_p_c.csv")
write.csv(perm_df_p_bvc, "~/scratch/brain/results/ieg_covar_c53_p1000_p_bvc.csv")

# perm_labels[[n_perm+1]] = bb$backup_subsample
# real_combos = combosRes(perm_labels[[n_perm+1]])

# Graph the results
# df_bvc_plot = data.frame(cluster1 = rep(perm_df_r_b$Cluster1, n_perm), cluster2 = rep(perm_df_r_b$Cluster2, n_perm),
#                          bvc = unlist( perm_df_r_b[,3:102] - perm_df_r_c[,3:102] ), isReal = rep(FALSE, n_perm))
# real_combos$bvc = as.numeric(as.vector(real_combos$r_behave)) - as.numeric(as.vector(real_combos$r_control))
# real_combos$abs_bvc = abs(real_combos$bvc)
# real_combos$isReal = TRUE
# df_bvc_plot = rbind(df_bvc_plot, real_combos[,c("cluster1", "cluster2", "bvc", "isReal")])
# df_bvc_plot$combo = paste0(df_bvc_plot$cluster1, "_", df_bvc_plot$cluster2)
# df_bvc_plot$abs_bvc = abs(df_bvc_plot$bvc)
# 
# png("~/scratch/brain/results/ieg_covar_c15_p100_r_bvc.png", width = 1500, height = 400, res = 90)
# ggplot(df_bvc_plot, aes(combo, bvc, fill = isReal, color = isReal)) + geom_boxplot(alpha = 0.7)
# dev.off()
# png("~/scratch/brain/results/ieg_covar_c15_p100_r_abs_bvc.png", width = 1500, height = 400, res = 90)
# ggplot(df_bvc_plot, aes(combo, abs_bvc, fill = isReal, color = isReal)) + geom_boxplot(alpha = 0.7)
# dev.off()
# 
# df_bvc_plot3 = real_combos[,c("cluster1", "cluster2", "bvc")]
# df_bvc_plot3 = cbind(df_bvc_plot3, perm_df_r_b[,3:102] - perm_df_r_c[,3:102])
# df_bvc_plot3$n_perm_greater = sapply(1:nrow(df_bvc_plot3), function(x) length(which(df_bvc_plot3[x,as.character(c(1:100))] > df_bvc_plot3[x,c("bvc")])) )
# # df_bvc_plot3 = df_bvc_plot3[order(as.numeric(as.vector(df_bvc_plot3$cluster1)), as.numeric(as.vector(df_bvc_plot3$cluster2))), ]
# df_bvc_plot3 = df_bvc_plot3[order(as.numeric(as.vector(df_bvc_plot3$cluster2))), ]
# 
# png("~/scratch/brain/results/ieg_covar_c15_p100_r_bvc_perm_greater_raw.png", width = 650, height = 600, res = 90)
# ggplot(df_bvc_plot3, aes(cluster1, cluster2, fill = n_perm_greater)) + geom_tile() + scale_fill_viridis(discrete=F, limits=c(0, max(df_bvc_plot2$n_perm_greater)))
# dev.off()
# 
# df_bvc_plot2 = real_combos[,c("cluster1", "cluster2", "abs_bvc")]
# df_bvc_plot2 = cbind(df_bvc_plot2, abs(perm_df_r_b[,3:102] - perm_df_r_c[,3:102]))
# df_bvc_plot2$n_perm_greater = sapply(1:nrow(df_bvc_plot2), function(x) length(which(df_bvc_plot2[x,as.character(c(1:100))] > df_bvc_plot2[x,c("abs_bvc")])) )
# # df_bvc_plot2 = df_bvc_plot2[order(as.numeric(as.vector(df_bvc_plot2$cluster1)), as.numeric(as.vector(df_bvc_plot2$cluster2))), ]
# df_bvc_plot2 = df_bvc_plot2[order(as.numeric(as.vector(df_bvc_plot2$cluster2))), ]
# 
# png("~/scratch/brain/results/ieg_covar_c15_p100_r_bvc_perm_greater.png", width = 650, height = 600, res = 90)
# ggplot(df_bvc_plot2, aes(cluster1, cluster2, fill = n_perm_greater)) + geom_tile() + scale_fill_viridis(discrete=F, limits=c(0, max(df_bvc_plot2$n_perm_greater)))
# dev.off()

