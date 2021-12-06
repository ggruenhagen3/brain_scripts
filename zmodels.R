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
  bbmm <- BBmm(fixed.formula = ff, random.formula = rf, m=88, data = df, show = TRUE)
  this_res = data.frame(summary(bbmm)$fixed.coefficients)
  return(this_res$p.value[which( rownames(this_res) == test_var)])
}


#**********************************************************************
# Body ================================================================
#**********************************************************************
# Load Data + Libraries
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("SeuratObject")
bb = readRDS(paste0(rna_path, "data/bb_demux_102021.rds"))

# Load Other Libraries
library(parallel)
library(qvalue)
library(tidyverse)
library(lme4)
library(PROreg)

# Clusters to Test
skip_clusters=c(12,13,14)
k = 0

# Temporary Random Scores (Same Mean as Neurogen Score -> 10)
run_vars = c("neurogen_score")
for (i in 1:100) {
  this_var = paste0("ran", i)
  run_vars = c(run_vars, this_var)
  bb@meta.data[, this_var] = abs(bb$neurogen_score + round(rnorm(n = ncol(bb))))
}

# Subset Data by Cluster
df = bb@meta.data[which(bb$seuratclusters15 == k),]
df$pair = as.factor(df$pair)
df$subject = as.factor(df$trial_id)

bbmm_start_time <- proc.time()[[3]]
res = unlist(mclapply(run_vars, function(x) myBBmm(x), mc.cores = detectCores()))
names(res) = run_vars
print(res)
# bbmm <- BBmm(fixed.formula = neurogen_score ~ as.numeric(bower_activity_index) + as.numeric(gsi), random.formula = ~ (subject %in% sample %in% run) + (subject %in% pair) , m=88, data = df, show = TRUE)
bbmm_stop_time = proc.time()[[3]]
print(paste0("BBmm on Cluster 0 (15 level) w/ 100 Randoms took: ", bbmm_stop_time-bbmm_start_time))

# saveRDS(bbmm, "~/scratch/brain/results/bbmm_test.rds")