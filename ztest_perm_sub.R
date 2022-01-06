# Helper Functions ***********************************************************************
singleRunGeneDefined = function(markers, returnP = T) {
  avg_exp = colSums(exp[markers,])
  avg_exp = avg_exp/bb$nFeature_RNA
  cluster_p = c()
  cluster_d = c()
  for (i in 1:nrow(zdf)) {
    gene = zdf$gene[i]
    cluster_cells = colnames(bb)[which(bb$cluster == zdf$cluster[i])]
    gene_pop_cells <- colnames(bb)[which(exp[gene, cluster_cells] > 0)]
    gene_pop_exp = avg_exp[gene_pop_cells]
    other_exp = avg_exp[which(! colnames(bb) %in% gene_pop_cells)]
    
    p = z.test(gene_pop_exp, other_exp, sigma.x = sd(gene_pop_exp), sigma.y = sd(other_exp), alternative = "greater")$p.value
    
    all_exp = c(gene_pop_exp, other_exp)
    test= effsize::cohen.d(all_exp, c(rep("cluster", length(gene_pop_exp)), rep("other",   length(other_exp))))
    d=test$estimate
    
    cluster_p = c(cluster_p, p)
    cluster_d = c(cluster_d, d)
  }
  
  cat(paste0("- Single Perm End Time: ", format(Sys.time(), "%X ")))
  
  if (returnP) 
    return(cluster_p)
  else
    return(cluster_d)
}

# Body *************************************************************************************
# Load BB
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("SeuratObject")
bb = readRDS(paste0(rna_path, "data/bb_demux_102021.rds"))

# Set Number of Permutations
nperm = 10000

# Set Cluster Level
cluster_level = "53"

# Load in Real PCRC List
pcrc = read.csv("~/scratch/brain/fst/pc_20_rc_20_10kb_bins_25kb_genes_on_lg_11_peak_by_bin.csv")[,2]

# Load Gene Cluster Combos
zdf15 = read.csv("~/scratch/brain/data/goi_by_15clusters_1plus_by_trial_id_122021.csv")
zdf53 = read.csv("~/scratch/brain/data/goi_by_53clusters_1plus_by_trial_id_122221.csv")
if (cluster_level == "15") { zdf = zdf15; bb$cluster = bb$seuratclusters15; }
if (cluster_level == "53") { zdf = zdf53; bb$cluster = bb$seuratclusters53; }
zdf$gene = str_replace(zdf$gene, "\\.1", "")
Idents(bb) = bb$cluster

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

# Find Results for the Real PCRC
# real_res = singleRunGeneDefined(pcrc, returnP = F)

library("parallel")
print(paste0("Doing Perms Start Time: ", format(Sys.time(), "%X")))
print("")
perm_res = mclapply(1:nperm, function(x) singleRunGeneDefined(ran_lists[[x]], returnP = F), mc.cores = detectCores())
perm_df = as.data.frame(t(as.data.frame(perm_res)))
rownames(perm_df) = 1:nperm
colnames(perm_df) = clusters

write.csv(perm_df, paste0("~/scratch/brain/results/ztest_perm_10k_", cluster_level, "_by_goi_010522.csv"))

# p_df = data.frame()
# # perm_df_log = -log10(perm_df)
# for (gene in zGenePops) {
#   print(gene)
#   pgene = str_replace(gene, "-", ".")
#   neg = length(which( perm_df[,pgene] <= real_res[which(zGenePops == gene)] ))
#   p_df = rbind(p_df, data.frame(cluster, neg))
# }
# p_df$gene = factor(p_df$gene, levels = unique(zGenePops))
# p_df$p = ((nperm - p_df$neg) / nperm)
# p_df$bh = p.adjust(p_df$p, method = "BH")
# p_df$bon = p.adjust(p_df$p, method = "bonferroni")
# # ggplot(p_df, aes(x = cluster, y = neg)) + geom_bar(stat = 'identity') + geom_text(aes(label=neg),hjust=0.5, vjust=1, color = 'white') + ggtitle("Number of Perms Less Than Or Equal to Real")
# # ggplot(p_df, aes(x = cluster, y = p))   + geom_bar(stat = 'identity', fill = 'gray60') + geom_text(aes(label=p),hjust=0.5, vjust=1, color = 'black')   + ggtitle("p per cluster") + theme_bw()
