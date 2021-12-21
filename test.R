args = commandArgs(trailingOnly=TRUE)
cur_level = as.character(args[1])

rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("SeuratObject")
bb = readRDS(paste0(rna_path, "data/bb_demux_102021.rds"))
Idents(bb) = bb$seurat_clusters

z15 = read.csv("~/scratch/brain/results/out_bb15_bbmm_demux_deg_all_tests_for_volcano_plotting_121321.csv")
z53 = read.csv("~/scratch/brain/results/out_bb53_glmmseq_demux_deg_all_tests_for_volcano_plotting.csv")
gcm15 = readRDS("~/scratch/brain/results/adjusted_glmmseq_ffm_15.rds")
gcm53 = readRDS("~/scratch/brain/results/adjusted_glmmseq_ffm_53.rds")

if  (cur_level == "53") { cz = z53 }                   else { cz = z15 }
if  (cur_level == "53") { cmeta = "seuratclusters53" } else { cmeta = "seuratclusters15" }
if  (cur_level == "53") { clusters = 0:52 }            else { clusters = 0:14 }
if  (cur_level == "53") { gcm = gcm53 }                else { gcm = gcm15 }

cz$zgenes = str_replace(cz$mzebra, pattern = "\\.", "-")
cz = cz[which(cz$zgenes %in% rownames(bb)),]
adj_sub_mean    = setNames(data.frame(matrix(0, nrow = nrow(cz), ncol = 38)), sort(unique(bb$subsample)))
data_sub_mean   = setNames(data.frame(matrix(0, nrow = nrow(cz), ncol = 38)), sort(unique(bb$subsample)))
counts_sub_mean = setNames(data.frame(matrix(0, nrow = nrow(cz), ncol = 38)), sort(unique(bb$subsample)))
cz$adj_mean_of_mean_b = cz$adj_mean_of_mean_c = cz$adj_mean_b = cz$adj_sign_pair = cz$adj_mean_c = 0
cz$data_mean_of_mean_b = cz$data_mean_of_mean_c = cz$data_mean_b = cz$data_sign_pair = cz$data_mean_c = 0
cz$counts_mean_of_mean_b = cz$counts_mean_of_mean_c = cz$counts_mean_b = cz$counts_sign_pair = cz$counts_mean_c = 0

for (i in 1:nrow(cz)) {
  # if (i %% 1000 == 0) { print(i) }
  print(i)
  czgene = cz$zgenes[i]
  cluster = cz$cluster[i]
  clust_idx = which(bb@meta.data[,cmeta] == i)
  gcm_df = data.frame(adj = gcm[czgene, clust_idx], data = bb@assays$RNA@data[czgene, clust_idx], counts = bb@assays$RNA@counts[czgene, clust_idx], sample = bb$sample[clust_idx], cond = bb$cond[clust_idx], subsample = bb$subsample[clust_idx], pair = bb$pair[clust_idx], bai = bb$bower_activity_index[clust_idx], gsi = bb$gsi[clust_idx], spawn = bb$log_spawn_events[clust_idx])
  if (nrow(gcm_df) > 0) {
    adj_agr    = aggregate(adj    ~ subsample + cond + pair + bai, gcm_df, mean)
    data_agr   = aggregate(data   ~ subsample + cond + pair + bai, gcm_df, mean)
    counts_agr = aggregate(counts ~ subsample + cond + pair + bai, gcm_df, mean)
    
    adj_pair_sign_vector = adj_agr$adj[order(adj_agr$pair)[c(FALSE, TRUE)]] - adj_agr$adj[order(adj_agr$pair)[c(TRUE, FALSE)]]
    data_pair_sign_vector = data_agr$data[order(data_agr$pair)[c(FALSE, TRUE)]] - data_agr$data[order(data_agr$pair)[c(TRUE, FALSE)]]
    counts_pair_sign_vector = counts_agr$counts[order(counts_agr$pair)[c(FALSE, TRUE)]] - counts_agr$counts[order(counts_agr$pair)[c(TRUE, FALSE)]]
    
    adj_sub_mean[i, ] = adj_agr$adj[match(colnames(adj_sub_mean), adj_agr$subsample)]
    data_sub_mean[i, ] = data_agr$data[match(colnames(data_sub_mean), data_agr$subsample)]
    counts_sub_mean[i, ] = counts_agr$counts[match(colnames(counts_sub_mean), counts_agr$subsample)]
    
    cz$adj_mean_of_mean_b = mean(adj_agr$adj[which(adj_agr$cond == "BHVE")])
    cz$adj_mean_of_mean_c = mean(adj_agr$adj[which(adj_agr$cond == "CTRL")])
    cz$adj_mean_b = mean(gcm_df$adj[which(gcm_df$cond == "BHVE")])
    cz$adj_mean_c = mean(gcm_df$adj[which(gcm_df$cond == "CTRL")])
    cz$adj_sign_pair = length(which(adj_pair_sign_vector > 0))
    cz$data_sign_pair = length(which(data_pair_sign_vector > 0))
    cz$counts_sign_pair = length(which(counts_pair_sign_vector > 0))
  }
}

write.csv(adj_sub_mean,    paste0("~/scratch/brain/results/zdeg", cur_level, "_adj_subsample_means.csv"))
write.csv(data_sub_mean,   paste0("~/scratch/brain/results/zdeg", cur_level, "_data_subsample_means.csv"))
write.csv(counts_sub_mean, paste0("~/scratch/brain/results/zdeg", cur_level, "_counts_subsample_means.csv"))
write.csv(cz,              paste0("~/scratch/brain/results/zdeg", cur_level, "_summary.csv"))
