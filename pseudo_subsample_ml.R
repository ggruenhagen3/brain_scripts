rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("SeuratObject")
bb = readRDS(paste0(rna_path, "data/bb_cc_04072021.RDS"))
Idents(bb) = bb$seurat_clusters


library(pacman)
p_unload(SeuratDisk)
p_unload(Seurat)
p_load(Seurat)


# 100 Permutations
for (i in 1:100) {
  print(i)
  subsample_stats_df = data.frame()
  bb$pseudo_subsample = "NA"
  set.seed(i)
  for (sub_sample in levels(bb$subsample)) {
    sub_sample_cells = colnames(bb)[which(bb$subsample == sub_sample)]
    train_num_cells = floor(length(sub_sample_cells) * 0.8)
    train_cells = sample(sub_sample_cells, train_num_cells)
    test_cells = sub_sample_cells[which(!sub_sample_cells %in% train_cells)]
    bb$pseudo_subsample[train_cells] = paste0(sub_sample, "_train")
    bb$pseudo_subsample[test_cells] = paste0(sub_sample, "_test")
    # subsample_stats_df = rbind(subsample_stats_df, t(c(length(sub_sample_cells), train_num_cells, length(test_cells))))
  }
  Idents(bb) = bb$pseudo_subsample
  pseudo_avg = myAverageExpression(bb, slot = "counts")
  pseudo_avg = t(pseudo_avg)
  write.csv(pseudo_avg, paste0("~/scratch/brain/data/pseudo_subsample_ml/pseudo_", i, ".csv"))
}
