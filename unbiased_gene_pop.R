#**********************************************************************
# Helper Functions ====================================================
#**********************************************************************
pct_FC_in_GP = function(gp) {
  gp_cells = colnames(bmat)[which(bmat[gp,] > 0)]
  gp_cells_1 = gp_cells[which(gp_cells %in% group1_cells)]
  gp_cells_2 = gp_cells[which(gp_cells %in% group2_cells)]
  pct_FC = pct_dif_avg_logFC(obj, cells.1 = gp_cells_1, cells.2 = gp_cells_2)
  pct_FC$mym = (pct_FC$pct_dif / 100) * pct_FC$avg_logFC * log(bmat_row_sums[gp])
  pct_FC = pct_FC[order(pct_FC$mym, decreasing = T),]
  top_50  = sum(pct_FC$mym[1:50])
  top_100 = sum(pct_FC$mym[1:100])
  top_500 = sum(pct_FC$mym[1:500])
  num_above_2  = length(which(pct_FC$mym >= 2))
  num_above_25 = length(which(pct_FC$mym >= 25))
  num_above_3  = length(which(pct_FC$mym >= 3))
  return(c(top_50, top_100, top_500, num_above_2, num_above_25, num_above_3))
}

#**********************************************************************
# Body ================================================================
#**********************************************************************
# Read Input
args = commandArgs(trailingOnly=TRUE)
obj_str = args[1]

# Load Libraries
library("parallel")
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("SeuratObject")

# Load Data
if (obj_str == 'bb') {
  obj = readRDS(paste0(rna_path, "data/bb_demux_102021.rds"))
  obj$group = obj$cond
  obj$group = plyr::revalue(obj$group, replace = c("BHVE" = 'group1', "CTRL" = 'group2'))
} else if (obj_str == 'clown') {
  obj = readRDS(paste0(rna_path, "data/anenomefish_clustered_061821.rds"))
  obj$subsample = paste0(obj$sample, "_", obj$subsample)
  obj$group = obj$sex
  obj$group = plyr::revalue(obj$group, replace = c("f" = 'group1', "m" = 'group2'))
}
# obj = switch(obj_str, 'bb' = readRDS(paste0(rna_path, "data/bb_demux_102021.rds")))
obj_subsamples = sort(unique(obj$subsample))
group1_cells = colnames(obj)[which(obj$group == "group1")]
group2_cells = colnames(obj)[which(obj$group == "group2")]

# Create a Binary Matrix
bmat = obj@assays$RNA@counts
bmat[which(bmat > 1)] = 1
bmat_row_sums = rowSums(bmat)

# Find Genes w/ At Least X Cells in Every Fish
sub_thresh = 5
sub_gpops = rownames(obj)
for (sub in obj_subsamples) {
  sub_sums = rowSums(bmat[,which(obj$subsample == sub)])
  sub_sums = sub_sums[which(sub_sums >= sub_thresh)]
  sub_gpops = sub_gpops[which(sub_gpops %in% names(sub_sums))]
}

# Find Gene Pops
gpop_thresh = 1000
gpops = names(bmat_row_sums)[which(bmat_row_sums >= gpop_thresh)]
gpops = gpops[which(gpops %in% sub_gpops)]

# Calculate My Metric
gpop_res = mclapply(gpops, function(gp) pct_FC_in_GP(gp), mc.cores = detectCores())
gp_df = setNames(as.data.frame(t(as.data.frame(gpop_res)), row.names = 1:length(gpops)), c("top_50", "top_100", "top_500", "num_above_2", "num_above_25", "num_above_3"))
gp_df$gene = gpops
gp_df$num_cells = bmat_row_sums[gpops]
gp_df = gp_df %>% relocate(gene:num_cells, .before = 1)

# Write data
write.csv(gp_df, paste0("~/scratch/brain/results/unbiased_gene_pop_", obj_str, ".csv"))