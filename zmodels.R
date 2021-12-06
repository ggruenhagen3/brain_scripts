# Load Data + Libraries
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("SeuratObject")
bb = readRDS(paste0(rna_path, "data/bb_demux_102021.rds"))

# Load Other Libraries
library(Seurat)
library(qvalue)
library(tidyverse)
library(lme4)
library(PROreg)

skip_clusters=c(12,13,14)

k = 0

df = bb@meta.data[which(bb$seuratclusters15 == k),]
df$pair = as.factor(df$pair)

bbmm_start_time <- proc.time()
bbmm <- BBmm(fixed.formula = neurogen_score ~ as.numeric(bower_activity_index) + as.numeric(gsi), random.formula = ~ (subject %in% sample %in% run) + (subject %in% pair) , m=88, data = df, show = TRUE)
bbmm_stop_time = proc.time()
print(paste0("BBmm on Cluster 0 (15 level) took: ", bbmm_stop_time-bbmm_start_time))

saveRDS(bbmm, "~/scratch/brain/results/bbmm_test.rds")