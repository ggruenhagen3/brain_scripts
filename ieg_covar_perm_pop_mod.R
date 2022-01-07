#==================================================================================================
# Helper Functions ================================================================================
#==================================================================================================
# 10 and 9
permSubsamples = function(x) {
  isEven = x %% 2 == 0
  if (isEven) { num1 = 10; } else { num1 = 9; }
  num2 = 19 - num1
  real_subs = unique(bb$subsample)
  b_subs = real_subs[1:19]
  c_subs = real_subs[20:38]
  new_b = c(sample(b_subs, num1), sample(c_subs, num2))
  new_c = c(b_subs[which(! b_subs %in% new_b)], c_subs[which(! c_subs %in% new_b)])
  my_replace = c(new_b, new_c)
  names(my_replace) = c(b_subs, c_subs)
  return(plyr::revalue(bb$subsample, my_replace))
}
# Single Run Function
combosRes = function(perm, cluster_level) {
  # Set the random samples
  bb$subsample = factor(perm_labels[[perm]])
  
  # Find the mean ieg_like_score in each cluster
  ptm <- proc.time()
  ieg_long_list = list()
  for (i in 1:length(pop_list)) {
    df = bb@meta.data[pop_list[[i]], c("ieg_like_score", "subsample")]
    ieg_long_list[[i]] = aggregate(ieg_like_score ~ subsample, df, mean, drop = F)[,2]
  }
  ieg_mat = as.matrix(do.call(cbind, ieg_long_list))
  rownames(ieg_mat) = levels(bb$subsample)
  
  # Find the Correlation between Cluster Combos in Behave Samples
  b_r_mat = cor(ieg_mat[which(startsWith(rownames(ieg_mat), "b")),], y = NULL, use = "pairwise.complete.obs") # use = "pairwise.complete.obs" removes the NA subsamples
  b_r_mat[upper.tri(b_r_mat, diag = T)] = NA
  
  # Find the Correlation between Cluster Combos in Control Samples
  c_r_mat = cor(ieg_mat[which(startsWith(rownames(ieg_mat), "c")),], y = NULL, use = "pairwise.complete.obs")
  c_r_mat[upper.tri(c_r_mat, diag = T)] = NA
  
  bvc_mat = b_r_mat - c_r_mat
  bvc_long = na.omit(melt(bvc_mat))
  bvc = bvc_long[,3]
  names(bvc) = paste0(bvc_long[,1], "_", bvc_long[,2])
  
  # end_time = proc.time()
  # print(paste0("Total time: ", end_time["elapsed"] - ptm["elapsed"]))
  cat(paste0("- Single Perm End Time: ", format(Sys.time(), "%X ")))
  
  return(bvc)
}

#==================================================================================================
# Main Body =======================================================================================
#==================================================================================================
# Load BB Dataset
# rna_path = "C:/Users/miles/Downloads/brain/"
rna_path = "~/scratch/brain/"
bb = readRDS(paste0(rna_path, "data/bb_demux_102021.rds"))
source(paste0(rna_path, "brain_scripts/all_f.R"))
set.seed(156) # seed for the script

# Load IEG and IEG Like Genes
# ieg_like = read.csv("C:/Users/miles/Downloads/ieg_like_fos_egr1_npas4_detected_011521.csv", stringsAsFactors = F)[,1]
# ieg_like = read.csv(paste0(rna_path, "/results/ieg_like_011521.txt"), stringsAsFactors = F)[,1]
bb$ieg_like_score <- bb$ieg_score

# Load Populations
# pop_df = read.csv("C:/Users/miles/Downloads/ieg_pops.csv")
pop_df = read.csv(paste0(rna_path, "/data/ieg_pops.csv"))

# Create Populations
pop_list = list()
for (i in 1:nrow(pop_df)) {
  print(i)
  my.level   = pop_df[i, "level"]
  my.cluster = pop_df[i, "cluster"]
  my.gp      = pop_df[i, "mzebra"]
  if (my.level == "primary")    { bb$cluster = bb$seuratclusters15 }
  if (my.level == "secondary")  { bb$cluster = bb$seuratclusters53 }
  if (my.level == "all")        { bb$cluster = "All"               }
  if (my.gp == FALSE) { bb$gp = TRUE } else { bb$gp = bb@assays$RNA@counts[my.gp,] > 0 }
  bb@meta.data[, paste0("pop", i)] = bb$cluster == my.cluster & bb$gp
  pop_list[[i]] = colnames(bb)[which( bb@meta.data[, paste0("pop", i)] )]
}

# Setup Permutations
n_perm = 100001
# perm_labels = lapply(1:n_perm, function(x) sample(unname(as.vector(bb$subsample))))
perm_labels = lapply(1:n_perm, function(x) permSubsamples(x))
bb$backup_subsample = bb$subsample

# Parallelize Finding Difference in Behave and Control for Cluster Combos
library("parallel")
numCores = detectCores()
all_combos = mclapply(1:n_perm, function(perm) combosRes(perm), mc.cores = numCores)

# For some reason, one of the 100k returned a vector of length 1377 instead of 1378.
# Double check that all returned vectors are the right size
bad_idx = mclapply(1:n_perm, function(perm) length(all_combos[[perm]]) != 231, mc.cores = numCores)
bad_idx = which(unlist(bad_idx))
print(paste0("Found ", length(bad_idx), " permutations that didn't have the required number of permutations."))
all_combos[bad_idx] = NULL

# Merge List of Combos into One Dataframe
perm_bvc_df = data.frame(t(plyr::ldply(all_combos, rbind))) # merge by name
perm_bvc_df = perm_bvc_df[,1:ncol(perm_bvc_df)]
colnames(perm_bvc_df) = 1:ncol(perm_bvc_df)
combos = colsplit(rownames(perm_bvc_df), pattern = "\\_", names = c('cluster1', 'cluster2')) # split combos vector into two columns
perm_bvc_df = cbind(combos, perm_bvc_df)

# Save the Results
write.csv(perm_bvc_df, "~/scratch/brain/results/ieg_covar_pop_mod_p100k_bvc.csv")

# # After permutations are done:
# # Load perm results
# perm_bvc_df = as.data.frame(data.table::fread("~/scratch/brain/results/ieg_covar_pop_mod_p100k_bvc.csv"))
# rownames(perm_bvc_df) = perm_bvc_df[,1]
# perm_bvc_df[,1] = NULL
# perm_bvc_df[,"100001"] = NULL
# 
# # Real Results
# perm_labels = list()
# perm_labels[[1]] = bb$subsample
# real_combos = combosRes(1)
# 
# df_bvc_plot3 = perm_bvc_df
# df_bvc_plot3$bvc = 0
# df_bvc_plot3[names(real_combos),"bvc"] = real_combos
# perm_greater_boolean = df_bvc_plot3[,as.character(c(1:(n_perm-1)))] > df_bvc_plot3$bvc
# df_bvc_plot3$n_perm_greater = rowSums(perm_greater_boolean)
# 
# df_bvc_plot3$cluster1 = factor(df_bvc_plot3$cluster1, levels = 0:14)
# df_bvc_plot3$cluster2 = factor(df_bvc_plot3$cluster2, levels = 0:14)
# 
# perm_greater_boolean_abs = abs(df_bvc_plot3[,as.character(c(1:(n_perm-1)))]) > abs(df_bvc_plot3$bvc)
# df_bvc_plot3$abs_n_perm_greater = rowSums(perm_greater_boolean_abs)
# 
# png("~/scratch/brain/results/ieg_covar_c15_p100k_r_bvc_perm_greater_raw_pop.png", width = 850, height = 800, res = 90)
# ggplot(df_bvc_plot3, aes(cluster1, cluster2, fill = n_perm_greater)) + geom_tile() + scale_fill_viridis(discrete=F, limits=c(0, n_perm), begin = 1, end = 0) + ggtitle("Behave - Control Correation")
# dev.off()
# 
# png("~/scratch/brain/results/ieg_covar_c15_p100k_r_abs_bvc_perm_greater_pop.png", width = 850, height = 800, res = 90)
# ggplot(df_bvc_plot3, aes(cluster1, cluster2, fill = abs_n_perm_greater)) + geom_tile() + scale_fill_viridis(discrete=F, limits=c(0, n_perm), begin = 1, end = 0) + labs(fill ="n_perm_greater") + ggtitle("Absolute Value Behave - Control Correation")
# dev.off()
# 
# df_bvc_plot3$p = df_bvc_plot3$n_perm_greater / n_perm
# df_bvc_plot3$q = p.adjust(df_bvc_plot3$p, method = "BH")
# length(which(df_bvc_plot3$p < 0.05))
# length(which(df_bvc_plot3$q < 0.05))
# 
# df_bvc_plot3$abs_p = df_bvc_plot3$abs_n_perm_greater / n_perm
# df_bvc_plot3$abs_q = p.adjust(df_bvc_plot3$abs_p, method = "BH")
# length(which(df_bvc_plot3$abs_p < 0.05))
# length(which(df_bvc_plot3$abs_q < 0.05))
# df_bvc_plot3$dup = sapply(1:nrow(df_bvc_plot3), function(x) length(which(duplicated(df_bvc_plot3[x,as.character(1:100000)]))))
# # df_bvc_plot3[which(df_bvc_plot3$abs_q < 0.05),c(1,2,"abs_bvc", "abs_n_perm_greater", "abs_p", "abs_q")]
# 
# # P value per combo
# z_scores = t(scale(t(df_bvc_plot3[, c("bvc", as.character(1:n_perm))]))) # scale combos per row
# p_from_z = lapply(1:nrow(z_scores), function(x) 2*pnorm(-abs(z_scores[x, "bvc"])) )
# # p_from_z = lapply(1:ncol(z_scores_mat), function(x) 2*pnorm(-abs(z_scores_mat["real",x]), mean = mean(z_scores_mat[,x]), sd = sd(z_scores_mat[,x])) )
# p = unlist(p_from_z)
# names(p) = rownames(z_scores)
# q = p.adjust(p, method = "BH")
# names(q) = rownames(z_scores)
# length(which(p < 0.05))
# length(which(q < 0.05))
# q[which.min(q)]
