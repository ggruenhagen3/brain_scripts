args = commandArgs(trailingOnly=TRUE)
cur_level = as.character(args[1])

rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("SeuratObject")
bb = readRDS(paste0(rna_path, "data/bb_demux_102021.rds"))
Idents(bb) = bb$seurat_clusters
sub_meta = aggregate(nCount_RNA ~ subsample + cond + pair, bb@meta.data, mean)
sub_meta = sub_meta[order(sub_meta$subsample),]

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
all_pairs = sort(unique(bb$pair))
adj_sub_mean    = setNames(data.frame(matrix(0, nrow = nrow(cz), ncol = 38)), sort(unique(bb$subsample)))
data_sub_mean   = setNames(data.frame(matrix(0, nrow = nrow(cz), ncol = 38)), sort(unique(bb$subsample)))
counts_sub_mean = setNames(data.frame(matrix(0, nrow = nrow(cz), ncol = 38)), sort(unique(bb$subsample)))
cz$adj_mean_of_mean_b = cz$adj_mean_of_mean_c = cz$adj_mean_b = cz$adj_sign_pair = cz$adj_mean_c = 0
cz$data_mean_of_mean_b = cz$data_mean_of_mean_c = cz$data_mean_b = cz$data_sign_pair = cz$data_mean_c = 0
cz$counts_mean_of_mean_b = cz$counts_mean_of_mean_c = cz$counts_mean_b = cz$counts_sign_pair = cz$counts_mean_c = 0
cz$num_non_zero_pair = 0

top_hits = order(cz$P_bower_activity_index, decreasing = F)[1:100]
# top_hits = c("30_LOC106674892", "11_LOC106674892", "29_LOC101483255", "2_arhgef18", "3_hs3st5", "36_LOC101470924", "20_LOC101470091")
# top_hits = which(cz$cluster_genes %in% top_hits)
for (i in top_hits) {
# for (i in 1:nrow(cz)) {
  if (i %% 1000 == 0) { print(i) }
  # print(i)
  czgene = cz$zgenes[i]
  cluster = cz$cluster[i]
  clust_idx = which(bb@meta.data[,cmeta] == cluster)
  gcm_df = data.frame(adj = gcm[czgene, clust_idx], data = bb@assays$RNA@data[czgene, clust_idx], counts = bb@assays$RNA@counts[czgene, clust_idx], subsample = bb$subsample[clust_idx], cond = bb$cond[clust_idx])
  gcm_df$subsample = factor(gcm_df$subsample, levels = sub_meta$subsample)
  if (nrow(gcm_df) > 0) {
    adj_agr    = aggregate(adj    ~ subsample, gcm_df, mean, drop = F)
    data_agr   = aggregate(data   ~ subsample, gcm_df, mean, drop = F)
    counts_agr = aggregate(counts ~ subsample, gcm_df, mean, drop = F)
    
    adj_agr[, c("pair", "cond")] = sub_meta[match(adj_agr$subsample, sub_meta$subsample), c("pair", "cond")]
    data_agr[, c("pair", "cond")] = sub_meta[match(data_agr$subsample, sub_meta$subsample), c("pair", "cond")]
    counts_agr[, c("pair", "cond")] = sub_meta[match(counts_agr$subsample, sub_meta$subsample), c("pair", "cond")]
    
    # Aggregate removes 0 subsamples and messes up BHVE vs CTRL
    # clust_idx = which(bb@meta.data[,cmeta] == cluster & bb$subsample == "b1.1")
    # gcm_df %>% mutate(b = factor(b, letters[1:5])) %>%
    #   group_by(b) %>%
    #   summarise(out = n()) %>%
    #   complete(b, fill = list(out = 0))
    # for (j in all_pairs) {
    #   
    # }
    
    adj_pair_sign_vector = adj_agr$adj[order(adj_agr$pair)[c(FALSE, TRUE)]] - adj_agr$adj[order(adj_agr$pair)[c(TRUE, FALSE)]]
    data_pair_sign_vector = data_agr$data[order(data_agr$pair)[c(FALSE, TRUE)]] - data_agr$data[order(data_agr$pair)[c(TRUE, FALSE)]]
    counts_pair_sign_vector = counts_agr$counts[order(counts_agr$pair)[c(FALSE, TRUE)]] - counts_agr$counts[order(counts_agr$pair)[c(TRUE, FALSE)]]

    
    adj_sub_mean[i, ] = adj_agr$adj[match(colnames(adj_sub_mean), adj_agr$subsample)]
    data_sub_mean[i, ] = data_agr$data[match(colnames(data_sub_mean), data_agr$subsample)]
    counts_sub_mean[i, ] = counts_agr$counts[match(colnames(counts_sub_mean), counts_agr$subsample)]
    
    cz$adj_mean_of_mean_b[i] = mean(adj_agr$adj[which(adj_agr$cond == "BHVE")], na.rm = T)
    cz$adj_mean_of_mean_c[i] = mean(adj_agr$adj[which(adj_agr$cond == "CTRL")], na.rm = T)
    cz$adj_mean_b[i] = mean(gcm_df$adj[which(gcm_df$cond == "BHVE")])
    cz$adj_mean_c [i]= mean(gcm_df$adj[which(gcm_df$cond == "CTRL")])
    cz$adj_sign_pair[i] = length(which(adj_pair_sign_vector > 0))
    
    cz$data_mean_of_mean_b[i] = mean(data_agr$data[which(data_agr$cond == "BHVE")], na.rm = T)
    cz$data_mean_of_mean_c[i] = mean(data_agr$data[which(data_agr$cond == "CTRL")], na.rm = T)
    cz$data_mean_b[i] = mean(gcm_df$data[which(gcm_df$cond == "BHVE")])
    cz$data_mean_c [i]= mean(gcm_df$data[which(gcm_df$cond == "CTRL")])
    cz$data_sign_pair[i] = length(which(data_pair_sign_vector > 0))
    
    cz$counts_mean_of_mean_b[i] = mean(counts_agr$counts[which(counts_agr$cond == "BHVE")], na.rm = T)
    cz$counts_mean_of_mean_c[i] = mean(counts_agr$counts[which(counts_agr$cond == "CTRL")], na.rm = T)
    cz$counts_mean_b[i] = mean(gcm_df$counts[which(gcm_df$cond == "BHVE")])
    cz$counts_mean_c [i]= mean(gcm_df$counts[which(gcm_df$cond == "CTRL")])
    cz$counts_sign_pair[i] = length(which(counts_pair_sign_vector > 0))
    cz$num_non_zero_pair[i] = length(which(! is.na(adj_pair_sign_vector) ))
  }
}
cz2 = cz[top_hits,]

write.csv(adj_sub_mean,    paste0("~/scratch/brain/results/zdeg", cur_level, "_adj_subsample_means.csv"))
write.csv(data_sub_mean,   paste0("~/scratch/brain/results/zdeg", cur_level, "_data_subsample_means.csv"))
write.csv(counts_sub_mean, paste0("~/scratch/brain/results/zdeg", cur_level, "_counts_subsample_means.csv"))
write.csv(cz,              paste0("~/scratch/brain/results/zdeg", cur_level, "_summary.csv"))

# Goes to bcs
library(qs)
qs::qsave(bb, "~/Downloads/bb_012122.qs")
qs::qsave()