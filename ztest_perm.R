# Helper Functions ***********************************************************************
singleRun = function(markers, returnP = T) {
  avg_exp = colSums(exp[markers,])
  avg_exp = avg_exp/bb$nFeature_RNA
  cluster_p = c()
  cluster_d = c()
  cluster_up = c()
  cluster_down = c()
  for (cluster in clusters) {
    cluster_cells <- WhichCells(bb, idents = cluster)
    clust_exp = avg_exp[cluster_cells]
    other_exp = avg_exp[which(! colnames(bb) %in% cluster_cells)]
    
    p = z.test(clust_exp, other_exp, sigma.x = sd(clust_exp), sigma.y = sd(other_exp), alternative = "greater")$p.value
    
    all_exp = c(clust_exp, other_exp)
    test= effsize::cohen.d(all_exp, c(rep("cluster", length(clust_exp)), rep("other",   length(other_exp))))
    d=test$estimate
    up=test$conf.int[2]
    down = test$conf.int[1]
    
    cluster_p = c(cluster_p, p)
    cluster_d = c(cluster_d, d)
    cluster_up = c(cluster_up, up)
    cluster_down = c(cluster_down, down)
  }
  if (returnP) 
    return(cluster_p)
  else
    return(cluster_d)
}

singleRunGeneDefined = function(markers, genePops = zGenePops, returnP = T) {
  avg_exp = colSums(exp[markers,])
  avg_exp = avg_exp/bb$nFeature_RNA
  cluster_p = c()
  cluster_d = c()
  cluster_up = c()
  cluster_down = c()
  for (gene in genePops) {
    gene_pop_cells <- colnames(bb)[which(exp[gene,] > 0)]
    gene_pop_exp = avg_exp[gene_pop_cells]
    other_exp = avg_exp[which(! colnames(bb) %in% gene_pop_cells)]
    
    p = z.test(gene_pop_exp, other_exp, sigma.x = sd(gene_pop_exp), sigma.y = sd(other_exp), alternative = "greater")$p.value
    
    all_exp = c(gene_pop_exp, other_exp)
    test= effsize::cohen.d(all_exp, c(rep("cluster", length(gene_pop_exp)), rep("other",   length(other_exp))))
    d=test$estimate
    up=test$conf.int[2]
    down = test$conf.int[1]
    
    cluster_p = c(cluster_p, p)
    cluster_d = c(cluster_d, d)    
    cluster_up = c(cluster_up, up)
    cluster_down = c(cluster_down, down)
  }
  
  cat(paste0("- Single Perm End Time: ", format(Sys.time(), "%X ")))
  
  if (returnP) 
    return(cluster_p)
  else
    return(cluster_down)
}

# Body *************************************************************************************
# Load BB
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("SeuratObject")
bb = readRDS(paste0(rna_path, "data/bb_demux_102021.rds"))
Idents(bb) = bb$seuratclusters53

# Set Number of Permutations
nperm = 10000

# Load in Real PCRC List
# pcrc = read.csv("~/scratch/brain/fst/pc_20_rc_20_10kb_bins_25kb_genes_on_lg_11_peak_by_bin.csv")[,2]
# zGenePops = read.csv("~/scratch/brain/data/goi_1plus_by_trial_id_122121.csv")[,1]
# pcrc = c("cobl", "ddr1", "fhod3", "LOC101476487", "LOC101476914", "LOC101477204", "LOC101479283", "plekhf2", "plekhg4b", "wdr73", "LOC101476922", "boc", "LOC101487687", "epha3", "erbb2", "met", "metrn", "pax6", "smo", "LOC101480727", "LOC101469466", "vegfa", "LOC101469419")
# pcrc = c("cobl", "ddr1", "fhod3", "grik5", "LOC101476914", "LOC101477204", "LOC101479283", "LOC105941351", "nbeal2", "plekhf2", "plekhg4b", "wdr73")
pcrc = c("LOC101463785", "hes1", "LOC101470264", "s100b", "LOC101484507", "LOC101464395", "id1", "fezf2", "jag1", "fabp7", "fgfbp3", "slc1a2", "LOC101463816", "selenop")

# Sort genes by their # of UMIs
gene_counts = data.frame(rowSums(bb@assays$RNA@counts))
gene_counts$gene = rownames(gene_counts)
gene_counts = gene_counts[order(gene_counts[,1]),]
pcrc_idx = which(gene_counts[,2] %in% pcrc)

# Find pools of genes with comparable expression levels as the real list
print(paste0("Finding Gene Pools Start Time: ", format(Sys.time(), "%X")))
ran_pools = list()
search_space = seq(-200, 200)
search_space = search_space[order(abs(search_space))][2:length(search_space)]
for (gene in pcrc) {
  gene_pcrc_idx = which(gene_counts[,2] == gene)
  ran_pools[[gene]] = c()
  search_space_i = 1
  while(length(ran_pools[[gene]]) < 100) {
    idx_to_try = gene_pcrc_idx + search_space[search_space_i]
    if (idx_to_try > 0 & idx_to_try <= nrow(bb) & (! idx_to_try %in% pcrc_idx) ) {
      ran_pools[[gene]] = c(ran_pools[[gene]], gene_counts[idx_to_try, 2])
    }
    search_space_i = search_space_i + 1
  }
}

# Create Random Lists of Equal Size to the real
ran_lists = lapply(1:nperm, function(x) {
  this_ran_list = c()
  for (gene in pcrc) { this_ran_list = c(this_ran_list, sample(ran_pools[[gene]], 1)) }
  return(this_ran_list)
})

# Prepare the Data
exp = GetAssayData(bb, assay = "RNA", slot='counts')
exp[which(exp > 0)] = 1
clusters = sort(unique(as.numeric(as.vector(Idents(bb)))))

# Find Results for the Random Gene Lists
# perm_res = lapply(1:nperm, function(x) singleRun(ran_lists[[x]], returnP = F))
# real_res = singleRun(pcrc, returnP = F)
# real_res_log = -log10(real_res)

library("parallel")
print(paste0("Doing Perms Start Time: ", format(Sys.time(), "%X")))
print("")
# perm_res = mclapply(1:nperm, function(x) singleRunGeneDefined(ran_lists[[x]], genePops = zGenePops, returnP = F), mc.cores = detectCores())
perm_res = mclapply(1:nperm, function(x) singleRun(ran_lists[[x]], returnP = F), mc.cores = detectCores())
perm_df = as.data.frame(t(as.data.frame(perm_res)))
rownames(perm_df) = 1:nperm
colnames(perm_df) = clusters

# Visualize Results
# perm_df = as.data.frame(t(as.data.frame(perm_res)))
# rownames(perm_df) = 1:nperm
# colnames(perm_df) = clusters
# perm_df_melt = melt(perm_df)
# perm_df_melt$neg_log_p = -log10(perm_df_melt$value)
# perm_df_melt$above = perm_df_melt$neg_log_p > real_res_log[as.numeric(as.vector(perm_df_melt$variable)) + 1]
  
# ggplot(perm_df_melt, aes(x = value, fill = above, color = above)) + geom_histogram() + facet_wrap(~ variable)
write.csv(perm_df, "~/scratch/brain/results/ztest_perm_q_10k_53_060722.csv")
# write.csv(perm_df, "~/scratch/brain/results/ztest_perm_10k_all_dgene_120321.csv")

# p_df = data.frame()
# perm_df_log = -log10(perm_df)
# for (cluster in clusters) {
#   # neg = length(which(perm_df_log[,as.character(cluster)] <= real_res_log))
#   neg = length(which(perm_df[,as.character(cluster)] <= real_res[cluster + 1]))
#   p_df = rbind(p_df, data.frame(cluster, neg))
# }
# p_df$cluster = factor(p_df$cluster, levels = clusters)
# p_df$p = ((nperm - p_df$neg) / nperm) * 100
# p_df$bon = p.adjust(p_df$p/100, method = "BH")

# ggplot(p_df, aes(x = cluster, y = neg)) + geom_bar(stat = 'identity') + geom_text(aes(label=neg),hjust=0.5, vjust=1, color = 'white') + ggtitle("Number of Perms Less Than Or Equal to Real")
# ggplot(p_df, aes(x = cluster, y = p))   + geom_bar(stat = 'identity', fill = 'gray60') + geom_text(aes(label=p),hjust=0.5, vjust=1, color = 'black')   + ggtitle("p per cluster") + theme_bw()

# # p_df = data.frame()
# # # perm_df_log = -log10(perm_df)
# # # for (gene in zGenePops) {
# # for (cluster in clusters) {
# #   # print(gene)
# #   # pgene = str_replace(gene, "-", ".")
# #   neg = length(which(perm_df[,as.character(cluster)] <= real_res[cluster+1]))
# #   # neg = length(which( perm_df[,pgene] <= real_res[which(zGenePops == gene)] ))
# #   p_df = rbind(p_df, data.frame(cluster, neg))
# # }
# # p_df$cluster = factor(p_df$cluster, levels = clusters)
# # # p_df$gene = factor(p_df$gene, levels = unique(zGenePops))
# # p_df$p = ((nperm - p_df$neg) / nperm)
# # p_df$bh = p.adjust(p_df$p, method = "BH")
# # p_df$bon = p.adjust(p_df$p, method = "bonferroni")
# # ggplot(p_df, aes(x = cluster, y = neg)) + geom_bar(stat = 'identity') + geom_text(aes(label=neg),hjust=0.5, vjust=1, color = 'white') + ggtitle("Number of Perms Less Than Or Equal to Real")
# # ggplot(p_df, aes(x = cluster, y = p))   + geom_bar(stat = 'identity', fill = 'gray60') + geom_text(aes(label=p),hjust=0.5, vjust=1, color = 'black')   + ggtitle("p per cluster") + theme_bw()
