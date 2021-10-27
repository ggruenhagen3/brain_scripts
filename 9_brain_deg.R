# Load Libraries
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("parallel")

# Set Paths
global_path <- "/storage/coda1/p-js585/0/ggruenhagen3/George/rich_project_pb1/data/9_brain/all/"
setwd(global_path)
regions = c("Pallidus", "Hippocampus", "Striatum")
tissue = c("GP", "HC", "STR")

# Read Objects
obj_list = lapply(regions, function(x) readRDS(paste0(x, "_clean.RDS")))
names(obj_list) = regions

# Find Cluster DEGs
for (region in regions) {
  Idents(obj_list[[region]]) = obj_list[[region]]$cluster
}

myFindAllMarkersParallel = function(obj) {
  deg = FindAllMarkers(obj)
  deg_sig = deg[which(deg$p_val_adj < 0.05),]
  return(deg_sig)
}

print("Finding DEGs")
all_deg_list = mclapply(obj_list, function(x) myFindAllMarkersParallel(obj), mc.cores = detectCores())
print("Done")

print("Writing DEGs")
for (i in 1:3) {
  this_deg = all_deg_list[i]
  region = regions[i]
  write.csv(this_deg, paste0(region, "_cluster_deg.csv"))
}
print("Done")