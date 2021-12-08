#**********************************************************************
# Helper Functions ====================================================
#**********************************************************************
myBBmm = function(x) {
  #' Expects that data has already been subset by cluster.
  #' Finds p value of bower_activity_index by score.
  #' 
  #' @param x score column
  test_var = "bower_activity_index"
  ff = as.formula(paste0(x, " ~ bower_activity_index + gsi"))
  rf = as.formula(" ~ (subject %in% sample %in% run) + (subject %in% pair)")
  bbmm <- BBmm(fixed.formula = ff, random.formula = rf, m=88, data = df, show = F)
  this_res = data.frame(summary(bbmm)$fixed.coefficients)
  
  print(paste0("Done ", x, ". ", proc.time()[3] - clust_start_time, "s elapsed."))
  return(this_res$p.value[which( rownames(this_res) == test_var)])
}

myBBmmLocal = function(x, clust) {
  #' Expects that data has already been subset by cluster.
  #' Finds p value of bower_activity_index by score.
  #' 
  #' @param x score column
  local_df = fread(paste0("~/scratch/brain/results/real_and_rand_meta/meta15/meta_", i, ".csv"))
  ff = as.formula(paste0(x, " ~ ", ff_str))
  rf = as.formula(" ~ (subject %in% sample %in% run) + (subject %in% pair)")
  bbmm <- BBmm(fixed.formula = ff, random.formula = rf, m=88, data = local_df, show = F)
  this_res = data.frame(summary(bbmm)$fixed.coefficients)
  this_res = this_res$p.value[2:nrow(this_res)]
  print(paste0("Done ", x))
  return(this_res)
}

myBBmmVector = function(x) {
  #' Expects that data has already been subset by cluster.
  #' Finds p value of bower_activity_index by score.
  #' 
  #' @param x score column
  test_var = "as.numeric(bower_activity_index)"
  ff = as.formula(paste0(x, " ~ as.numeric(bower_activity_index) + as.numeric(gsi)"))
  rf = as.formula(" ~ (subject %in% sample %in% run) + (subject %in% pair)")
  # x2 = df[,x]
  # bai = df[, "bower_activity_index"]
  bbmm <- BBmm(fixed.formula = df[, x] ~ df[, "bower_activity_index"] + df[, "gsi"], random.formula = ~ (df[,"subject"] %in% df[,"sample"] %in% df[, "run"]) + (df[,"subject"] %in% df[,"pair"]), m=88, show = TRUE)
  this_res = data.frame(summary(bbmm)$fixed.coefficients)
  print(paste0("Done ", x))
  return(this_res$p.value[which( rownames(this_res) == test_var)])
}

#**********************************************************************
# Body ================================================================
#**********************************************************************
# Read Input
args = commandArgs(trailingOnly=TRUE)
gene_list = args[1]
this_model = as.numeric(args[2])
print(paste0("Performing BootStraps on ", gene_list))
print(paste0("Using Model ", this_model))

# Load Libraries
library(parallel)
library(qvalue)
library(tidyverse)
library(lme4)
library(PROreg)
library(data.table)

# # Load Data
# bbmm_start_time <- proc.time()[[3]]
# rna_path = "~/scratch/brain/"
# # source(paste0(rna_path, "brain_scripts/all_f.R"))
# # library("SeuratObject")
# bb = readRDS(paste0(rna_path, "data/bb_demux_102021.rds"))
# # bb_metadata = read.csv(paste0(rna_path, "data/bb_meta_data.csv"))
# bb_metadata = fread(paste0(rna_path, "data/bb_meta_data.csv"))
# bbmm_stop_time = proc.time()[[3]]
# print(paste0("Loading took: ", bbmm_stop_time-bbmm_start_time))


# # Clusters to Test
# skip_clusters=c(12,13,14)
# k = 0

# # Random Scores
# neurogen = read.csv("~/scratch/brain/data/conserved_neurogenesis_positive_zfish_mouse_cichlid.csv")[,4]
# neurogen = read.csv("~/scratch/brain/data/ieg_like_fos_egr1_npas4_detected_011521.csv")[,1]
# gene_counts = data.frame(rowSums(bb@assays$RNA@counts))
# gene_counts$gene = rownames(gene_counts)
# gene_counts = gene_counts[order(gene_counts[,1]),]
# neurogen_idx = which(gene_counts[,2] %in% neurogen)
# ran_pools = list()
# search_space = seq(-200, 200)
# search_space = search_space[order(abs(search_space))][2:length(search_space)]
# for (gene in neurogen) {
#   print(gene)
#   gene_neurogen_idx = which(gene_counts[,2] == gene)
#   ran_pools[[gene]] = c()
#   search_space_i = 1
#   while(length(ran_pools[[gene]]) < 100) {
#     idx_to_try = gene_neurogen_idx + search_space[search_space_i]
#     if (idx_to_try > 0 & idx_to_try <= nrow(bb) & (! idx_to_try %in% neurogen_idx) ) {
#       ran_pools[[gene]] = c(ran_pools[[gene]], gene_counts[idx_to_try, 2])
#     }
#     search_space_i = search_space_i + 1
#   }
# }
# # Create Random Lists of Equal Size to the real
# ran_lists = lapply(1:nperm, function(x) {
#   this_ran_list = c()
#   for (gene in neurogen) { this_ran_list = c(this_ran_list, sample(ran_pools[[gene]], 1)) }
#   return(this_ran_list)
# })
# # Create Scores for random lists
# mat = bb@assays$RNA@counts
# mat[which(mat > 1)] = 1
# for (x in 1:100) { bb@meta.data[,paste0("neurogen", x)] = colSums(mat[ran_lists[[x]],]) }
# bb$subject = bb$trial_id
# df = bb@meta.data[,c("subject", "sample", "run", "pair", "neurogen_score", "ieg_score", "bower_activity_index", "gsi", paste0("neurogen", 1:100), paste0("ieg", 1:100))]
# write.csv(df, "~/scratch/brain/results/bb_real_and_rand_meta.csv")
# for (i in 0:11) {
#   clust_df = df[colnames(bb)[which(bb$seuratclusters15 == i)],]
#   write.csv(clust_df, paste0("~/scratch/brain/results/real_and_rand_meta/meta15/meta_", i, ".csv"))
# }

# # Subset Data by Cluster
# df = bb_metadata[which(bb_metadata$seuratclusters15 == k),]
# df$log_spawn_events = as.numeric(df$log_spawn_events)
# df$bower_activity_index = as.numeric(df$bower_activity_index)
# df$gsi = as.numeric(df$gsi)
# df$neurogen_score = as.numeric(df$neurogen_score)
# df$pair = as.factor(df$pair)
# df$subject = as.factor(df$trial_id)
# df$cond = as.factor(df$cond)

ff_str = switch(this_model, 
                "cond + log_spawn_events + gsi",                 # Model 1
                "cond + log_spawn_events",                       # Model 2
                "cond + gsi",                                    # Model 3
                "bower_activity_index + log_spawn_events + gsi", # Model 4
                "bower_activity_index + log_spawn_events",       # Model 5
                "bower_activity_index + gsi",                    # Model 6
                "log_spawn_events + gsi")                        # Model 7
print(paste0("Model Fixed Effects formula: score ~ ", ff_str))

# Parameters and Print Messages
nperm = 2
run_vars = paste0(gene_list, 1:nperm)
num.cores = detectCores()
print(paste0("Performing ", nperm, " Bootstraps"))
print(paste0("Number of Cores: ", num.cores))
print(paste0("Start Time: ", format(Sys.time(), "%X")))
very_start_time <- proc.time()[[3]]

# Calculate p's by cluster
all_p = data.frame()
for (i in 0:11) {
  print(paste0("--- Cluster ", i, " --- "))
  print(paste0("Cluster Start Time: ", format(Sys.time(), "%X")))
  clust_start_time <- proc.time()[[3]]
  clust_p = mclapply(run_vars, function(x) myBBmmLocal(x, i), mc.cores = num.cores)
  all_p = rbind(all_p, do.call(rbind, clust_p))
  print(paste0("Cluster End Time: ", format(Sys.time(), "%X")))
  clust_end_time <- proc.time()[[3]]
  print(paste0("Number of Seconds Elapsed for Cluster ", i, ": ", clust_end_time-clust_start_time))
}

# Write Data and Finish Up
print("All Done")
write.csv(all_p, paste0("~/scratch/brain/results/bbmm15_", gene_list, "_model", this_model, ".csv"))
very_stop_time = proc.time()[[3]]
print(paste0("Total End Time: ", format(Sys.time(), "%X")))
print(paste0("Number of Seconds Elapsed: ", very_stop_time-very_start_time))
