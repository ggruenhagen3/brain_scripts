rna_path = "~/scratch/brain/"
bb = readRDS(paste0(rna_path, "data/bb_subsample_02222021.RDS"))
source(paste0(rna_path, "brain_scripts/all_f.R"))

bsub_ctrl = as.vector(bb$subsample)
bsub_ctrl[which(startsWith(bsub_ctrl, "c"))] = "CTRL"

Idents(bb) = bsub_ctrl
bsubs = sort(unique(bsub_ctrl))[1:(length(sort(unique(bsub_ctrl))))-1]
all_res = data.frame()
p_df = data.frame()
for (bsub in bsubs) {
  print(bsub)
  this_res = FindMarkers(bb, ident.1 = bsub, ident.2 = "CTRL", min.pct = 0.01, logfc.threshold = 0.10)
  this_res$gene = rownames(this_res)
  this_res$subsample = bsub
  this_res = this_res[which(this_res$p_val_adj < 0.05),]
  all_res = rbind(all_res, this_res)
  
  newRow = data.frame(subsample = bsub, n_deg = nrow(this_res), avg_avg_logFC = mean(abs(this_res$avg_logFC)), bsub_size = length(which(Idents(bb) == bsub)))
  p_df = rbind(p_df, newRow)
}

p_df = data.frame()
for (bsub in bsubs) {
  this_res = all_res[which(all_res$subsample == bsub),]
  this_res$pct_dif = abs(this_res$pct.1 - this_res$pct.2)
  newRow = data.frame(subsample = bsub, n_deg = nrow(this_res), avg_avg_logFC = mean(abs(this_res$avg_logFC)), bsub_size = length(which(Idents(bb) == bsub)), avg_pct_dif = mean(this_res$pct_dif))
  p_df = rbind(p_df, newRow)
}

p_df$sample = substr(p_df$subsample, 1, 2)
ggplot(p_df, aes(subsample, n_deg)) + geom_bar(stat="identity") + ylab("Number of DEGs")
ggplot(p_df, aes(subsample, avg_avg_logFC)) + geom_bar(stat="identity") + ylab("Average avg_logFC")
ggplot(p_df, aes(subsample, avg_pct_dif)) + geom_bar(stat="identity") + ylab("Average pct_dif")

hist(table(all_res$gene), breaks = 50, xlab = "Number of Subsamples Found As DEG", main="Overlap of Genes in Subsample DEGs")