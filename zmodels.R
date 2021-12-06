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
bb$ran1 = round(rnorm(n = ncol(bb), mean = 10))
bb$ran2 = round(rnorm(n = ncol(bb), mean = 10))
bb$ran3 = round(rnorm(n = ncol(bb), mean = 10))
bb$ran4 = round(rnorm(n = ncol(bb), mean = 10))
bb$ran5 = round(rnorm(n = ncol(bb), mean = 10))
bb$ran6 = round(rnorm(n = ncol(bb), mean = 10))
bb$ran7 = round(rnorm(n = ncol(bb), mean = 10))
bb$ran8 = round(rnorm(n = ncol(bb), mean = 10))
bb$ran9 = round(rnorm(n = ncol(bb), mean = 10))
bb$ran0 = round(rnorm(n = ncol(bb), mean = 10))
run_vars = c("neurogen_score", paste0("ran", 0:9))

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
print(paste0("BBmm on Cluster 0 (15 level) w/ 10 Randoms took: ", bbmm_stop_time-bbmm_start_time))

# saveRDS(bbmm, "~/scratch/brain/results/bbmm_test.rds")