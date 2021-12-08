# *** MODEL 4 ***
# 1. bbmm <- BBmm(fixed.formula = neurogen_score ~ as.factor(cond) + as.numeric(log_spawn_events) + as.numeric(gsi), random.formula = ~ (subject %in% sample %in% run) + (subject %in% pair) , m=88, data = df, show = TRUE)
# 2. bbmm <- BBmm(fixed.formula = neurogen_score ~ as.factor(cond) + as.numeric(gsi), random.formula = ~ (subject %in% sample %in% run) + (subject %in% pair) , m=88, data = df, show = TRUE)
# 3. bbmm <- BBmm(fixed.formula = neurogen_score ~ as.factor(cond) + as.numeric(gsi), random.formula = ~ (subject %in% sample %in% run) + (subject %in% pair) , m=88, data = df, show = TRUE)
# 4. bbmm <- BBmm(fixed.formula = neurogen_score ~ as.numeric(bower_activity_index) + as.numeric(log_spawn_events) + as.numeric(gsi), random.formula = ~ (subject %in% sample %in% run) + (subject %in% pair) , m=88, data = df, show = TRUE) (edited) 
# 5. bbmm <- BBmm(fixed.formula = neurogen_score ~ as.numeric(bower_activity_index) + as.numeric(log_spawn_events), random.formula = ~ (subject %in% sample %in% run) + (subject %in% pair) , m=88, data = df, show = TRUE) (edited) 
# 6. bbmm <- BBmm(fixed.formula = neurogen_score ~ as.numeric(bower_activity_index) + as.numeric(gsi), random.formula = ~ (subject %in% sample %in% run) + (subject %in% pair) , m=88, data = df, show = TRUE)
# 7. bbmm <- BBmm(fixed.formula = neurogen_score ~ as.numeric(log_spawn_events) + as.numeric(gsi), random.formula = ~ (subject %in% sample %in% run) + (subject %in% pair) , m=88, data = df, show = TRUE) (edited) 

#**********************************************************************
# Helper Functions ====================================================
#**********************************************************************
myBBmm = function(x) {
  #' Expects that data has already been subset by cluster.
  #' Finds p value of bower_activity_index by score.
  #' 
  #' @param x score column
  test_var = "as.numeric(bower_activity_index)"
  ff = as.formula(paste0(x, " ~ as.numeric(bower_activity_index) + as.numeric(gsi)"))
  rf = as.formula(" ~ (subject %in% sample %in% run) + (subject %in% pair)")
  bbmm <- BBmm(fixed.formula = ff, random.formula = rf, m=88, data = df2, show = TRUE)
  this_res = data.frame(summary(bbmm)$fixed.coefficients)
  print(paste0("Done ", x))
  return(this_res$p.value[which( rownames(this_res) == test_var)])
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
# Load Libraries
library(parallel)
library(qvalue)
library(tidyverse)
library(lme4)
library(PROreg)
library(data.table)

# Load Data
bbmm_start_time <- proc.time()[[3]]
rna_path = "~/scratch/brain/"
# source(paste0(rna_path, "brain_scripts/all_f.R"))
# library("SeuratObject")
# bb = readRDS(paste0(rna_path, "data/bb_demux_102021.rds"))
# bb_metadata = read.csv(paste0(rna_path, "data/bb_meta_data.csv"))
bb_metadata = fread(paste0(rna_path, "data/bb_meta_data.csv"))
bbmm_stop_time = proc.time()[[3]]
print(paste0("Loading took: ", bbmm_stop_time-bbmm_start_time))


# Clusters to Test
skip_clusters=c(12,13,14)
k = 0

# Temporary Random Scores (Same Mean as Neurogen Score -> 10)
run_vars = c("neurogen_score")
for (i in 1:10) {
  this_var = paste0("ran", i)
  run_vars = c(run_vars, this_var)
  bb_metadata[, this_var] = abs(bb_metadata$neurogen_score + round(rnorm(n = nrow(bb_metadata))))
}

# Subset Data by Cluster
df = bb_metadata[which(bb_metadata$seuratclusters15 == k),]
df$log_spawn_events = as.numeric(df$log_spawn_events)
df$bower_activity_index = as.numeric(df$bower_activity_index)
df$gsi = as.numeric(df$gsi)
df$neurogen_score = as.numeric(df$neurogen_score)
df$pair = as.factor(df$pair)
df$subject = as.factor(df$trial_id)
df$cond = as.factor(df$cond)

# Do smaller dataframes run faster?
# df = df[,c("subject", "sample", "run", "pair", "neurogen_score", "bower_activity_index", "gsi")]

num.cores = detectCores()
print(paste0("Number of Cores: ", num.cores))
print(paste0("BBmm Start Time: ", format(Sys.time(), "%X")))
bbmm_start_time <- proc.time()[[3]]
# res = myBBmm("neurogen_score")
# res2 = myBBmmVector("neurogen_score")
# res = unlist(mclapply(run_vars, function(x) myBBmm(x), mc.cores = num.cores))
# names(res) = run_vars
# print(res)
bbmm <- BBmm(fixed.formula = neurogen_score ~ bower_activity_index + gsi, random.formula = ~ (subject %in% sample %in% run) + (subject %in% pair) , m=88, data = df, show = TRUE)
bbmm_stop_time = proc.time()[[3]]
print(paste0("BBmm on Cluster 0 (15 level) on Real took: ", bbmm_stop_time-bbmm_start_time))
print(paste0("BBmm End Time: ", format(Sys.time(), "%X")))
