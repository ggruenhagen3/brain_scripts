#rna_path = "C:/Users/miles/Downloads/brain/"
rna_path = "~/research/brain/"
setwd(rna_path)
# rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("SeuratObject")
# bb = readRDS(paste0(rna_path, "data/bb_subsample_02222021.RDS"))
bb = readRDS(paste0(rna_path, "data/bb_demux_102021.rds"))
Idents(bb) = bb$seurat_clusters


library(pacman)
p_unload(SeuratDisk)
p_unload(Seurat)
p_load(Seurat)

#==========================================================================================
# Final Object ============================================================================
#==========================================================================================
rna_path = "C:/Users/miles/Downloads/brain/"
bb = readRDS(paste0(rna_path, "data/bb_clustered_102820.rds"))

loop_samples =  c("b1", "b2", "b3", "b4", "b5", "c1", "c2", "c3", "c4", "c5")
loop_samples = c("c5")
for (s in loop_samples) {
  print(s)
  soup = read.table(paste0("~/scratch/brain/ffm/JTS07-", toupper(s), "/outs/soup_end.tsv"), header = F, sep = "\t")
  colnames(soup)[1:2] = c("cell", "best")
  soup$demux = bb$demux[which(bb$sample == s)]
  soup$best_demux = paste0(soup$best, "_", soup$demux)
  sd_df = as.data.frame(table(soup$best_demux))
  sd_df[,c("soup", "demux")] = reshape2::colsplit(sd_df[,1], "_", c("1", "2"))
  print(sd_df)
  
  all_soup_max_num = 0
  all_soup_other_num = 0
  for (soup_s in unique(sd_df$soup)) {
    print(soup_s)
    max_num = max(sd_df$Freq[which(sd_df$soup == soup_s)])
    other_num = sum(sd_df$Freq[which(sd_df$soup == soup_s & sd_df$Freq != max_num)])
    print(max_num/(max_num+other_num))
    all_soup_max_num = all_soup_max_num + max_num
    all_soup_other_num = all_soup_other_num + other_num
  }
  print(paste("All Soup #:", all_soup_max_num))
  print(paste("All Soup %:", all_soup_max_num/(all_soup_max_num+all_soup_other_num)))
  
  # png(paste0("~/scratch/brain/results/soup_to_demux", s, ".png"), width = 500, height = 400)
  # print(ggplot(sd_df, aes(x=soup, y = Freq, color = demux, fill = demux)) + geom_bar() + ggtitle(s))
  # dev.off()
}

# Add Depth Change
bb$depth = 0
bb$depth[which(bb$sample == "b1")] = 70.7531187
bb$depth[which(bb$sample == "b2")] = 78.17879535
bb$depth[which(bb$sample == "b3")] = 36.58969697
bb$depth[which(bb$sample == "b4")] = 170.4355415
bb$depth[which(bb$sample == "b5")] = 118.5065215

# Add GSI
bb$gsi = 0
bb$gsi[which(bb$sample == "b1")] = 0.4225
bb$gsi[which(bb$sample == "b2")] = 0.625
bb$gsi[which(bb$sample == "b3")] = 0.545
bb$gsi[which(bb$sample == "b4")] = 0.5466667
bb$gsi[which(bb$sample == "b5")] = 0.715
bb$gsi[which(bb$sample == "c1")] = 0.4175
bb$gsi[which(bb$sample == "c2")] = 0.5325
bb$gsi[which(bb$sample == "c3")] = 0.37825
bb$gsi[which(bb$sample == "c4")] = 0.4133333
bb$gsi[which(bb$sample == "c5")] = 0.365

# Add RNA velocity stuff to bb
all_rna = read.csv("~/research/brain/data/rna_velo_all_data.csv", stringsAsFactors = F)
all_rna$good_names = paste0(all_rna$cond, "_", all_rna$X, "_1-1")
df = data.frame()
function()
for (barcode in all_rna$X) {
  new_cell = colnames(bb)[which(grepl(barcode, colnames(bb)))]
  if (length(new_cell) > 0) {
    newRow = data.frame(barcode = barcode, bb_cell = new_cell)
    df = rbind(df, newRow)
  }
}

# Find DEGs between each sample
all_deg = list()
sig_deg = list()
Idents(bb) = bb$sample
bb_samples = unique(bb$sample)
sample_df = data.frame()
for (i in 2:length(bb_samples)) {
  print(paste0("i = ", i))
  for (j in 1:(i-1)) {
    this_deg = FindMarkers(bb, ident.1 = bb_samples[i], ident.2 = bb_samples[j], logfc.threshold = 0.1, min.pct = 0.01)
    this_sig_deg = this_deg[which(this_deg$p_val_adj < 0.05),]
    all_deg[[paste0(bb_samples[i], "_", bb_samples[j])]] = this_deg
    sig_deg[[paste0(bb_samples[i], "_", bb_samples[j])]] = this_sig_deg
    sample_df = rbind(sample_df, t(c(i, j, bb_samples[i], bb_samples[j], nrow(this_deg), nrow(this_sig_deg))))
  }
}
sample_df[c(1,2,5,6)] <- lapply(sample_df[c(1,2,5,6)], as.numeric)
colnames(sample_df) = c("i", "j", "sample1", "sample2", "all_deg", "sig_deg")

png("~/scratch/brain/results/bb_sample_sig_deg.png", width = 600, height = 600)
print(ggplot(sample_df, aes(x=sample1, y = sample2, fill = sig_deg)) + geom_tile() + scale_fill_viridis() + ggtitle("Significant DEGs b/w Pairwise Samples"))
dev.off()

png("~/scratch/brain/results/bb_sample_all_deg.png", width = 600, height = 600)
print(ggplot(sample_df, aes(x=sample1, y = sample2, fill = all_deg)) + geom_tile() + scale_fill_viridis() + ggtitle("All DEGs b/w Pairwise Samples"))
dev.off()

png("~/scratch/brain/results/bb_sample_sig_deg_red_blue.png", width = 600, height = 600)
print(ggplot(sample_df, aes(x=sample1, y = sample2, fill = sig_deg)) + geom_tile() + scale_fill_gradientn(colors = pal(50)) + ggtitle("Significant DEGs b/w Pairwise Samples"))
dev.off()

png("~/scratch/brain/results/bb_sample_all_deg_red_blue.png", width = 600, height = 600)
print(ggplot(sample_df, aes(x=sample1, y = sample2, fill = all_deg)) + geom_tile() + scale_fill_gradientn(colors = pal(50)) + ggtitle("All DEGs b/w Pairwise Samples"))
dev.off()

# Find DEGs between each subsample
all_deg = list()
sig_deg = list()
Idents(bb) = bb$subsample
bb_susamples = unique(bb$subsample)
subsample_df = data.frame()
for (i in 2:length(bb_susamples)) {
  print(paste0("i = ", i))
  for (j in 1:(i-1)) {
    this_deg = FindMarkers(bb, ident.1 = bb_susamples[i], ident.2 = bb_susamples[j], logfc.threshold = 0.1, min.pct = 0.01)
    this_sig_deg = this_deg[which(this_deg$p_val_adj < 0.05),]
    all_deg[[paste0(bb_susamples[i], "_", bb_susamples[j])]] = this_deg
    sig_deg[[paste0(bb_susamples[i], "_", bb_susamples[j])]] = this_sig_deg
    subsample_df = rbind(subsample_df, t(c(i, j, bb_susamples[i], bb_susamples[j], nrow(this_deg), nrow(this_sig_deg))))
  }
}
subsample_df[c(1,2,5,6)] <- lapply(subsample_df[c(1,2,5,6)], as.numeric)
colnames(subsample_df) = c("i", "j", "subsample1", "subsample2", "all_deg", "sig_deg")
subsample_df$subsample1 = as.vector(bb_susamples[subsample_df$i])
subsample_df$subsample2 = as.vector(bb_susamples[subsample_df$j])

png("~/scratch/brain/results/bb_subsample_sig_deg.png", width = 900, height = 900)
print(ggplot(subsample_df, aes(x=subsample1, y = subsample2, fill = sig_deg)) + geom_tile() + scale_fill_viridis() + ggtitle("Significant DEGs b/w Pairwise Subsamples"))
dev.off()

png("~/scratch/brain/results/bb_subsample_all_deg.png", width = 900, height = 900)
print(ggplot(subsample_df, aes(x=subsample1, y = subsample2, fill = all_deg)) + geom_tile() + scale_fill_viridis() + ggtitle("All DEGs b/w Pairwise Subsamples"))
dev.off()

png("~/scratch/brain/results/bb_subsample_sig_deg_red_blue.png", width = 900, height = 900)
print(ggplot(subsample_df, aes(x=subsample1, y = subsample2, fill = sig_deg)) + geom_tile() + scale_fill_gradientn(colors = pal(50)) + ggtitle("Significant DEGs b/w Pairwise Subsamples"))
dev.off()

png("~/scratch/brain/results/bb_subsample_all_deg_red_blue.png", width = 900, height = 900)
print(ggplot(subsample_df, aes(x=subsample1, y = subsample2, fill = all_deg)) + geom_tile() + scale_fill_gradientn(colors = pal(50)) + ggtitle("All DEGs b/w Pairwise Subsamples"))
dev.off()

# PCA of samples
library("factoextra")
sample_avg = read.table("~/research/brain/data/t_sample_avg.txt")
test = sample_avg[,which(colSums(sample_avg) != 0)]
res.pca <- prcomp(test, scale = T)
groups = factor(c(rep("Behave", 5), rep("Control", 5)))
fviz_pca_ind(res.pca, col.ind = groups, palette = c("#00AFBB",  "#FC4E07"), addEllipses = T, ellipse.type = "confidence", legend.title = "Groups", repel = T)

# PCA of subsamples
library("factoextra")
subsample_avg = read.table("~/research/brain/data/bb_subsample_data.txt")
test = subsample_avg[,which(colSums(subsample_avg) != 0)]
res.pca <- prcomp(test, scale = T)
groups = factor(c(rep("Behave", 19), rep("Control", 19)))
groups2 = factor(colsplit(as.vector(rownames(res.pca$x)), pattern = "\\.", names = c('1', "2"))[,1])
groups3 = factor(substr(groups2, 2, 2))
fviz_pca_ind(res.pca, col.ind = groups, palette = c("#00AFBB",  "#FC4E07"), addEllipses = T, ellipse.type = "confidence", legend.title = "Groups", repel = T)
fviz_pca_ind(res.pca, col.ind = groups2, addEllipses = T, ellipse.type = "confidence", legend.title = "Groups", repel = T)
fviz_pca_ind(res.pca, col.ind = groups3, addEllipses = T, ellipse.type = "confidence", legend.title = "Groups", repel = T)

lg11 = read.table("~/research/brain/data/markers/LG11_highFST_genes_umd2a_120920.txt", stringsAsFactors = F)[,1]
mat = as.matrix(bb@assays$RNA@counts[lg11,])
mat[which(mat > 0)] = 1
mat_scale = scale(colSums(mat)/bb$nFeature_RNA)
bb$enrich_z = mat_scale
abs_max = max(abs(bb$enrich_z))
floor(seq(-abs_max, abs_max))
myFeaturePlot(bb, "enrich_z") + scale_color_gradientn(colors = pal2(100), limits = c(-abs_max, abs_max))
myFeaturePlot(bb, "enrich_z") + scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(50), limits = c(-abs_max, abs_max))
myFeaturePlot(bb, "enrich_z") + scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11, "RdBu")))(50), limits = c(-abs_max, abs_max))
myFeaturePlot(bb, "enrich_z") + scale_color_gradientn(colors = colorRampPalette(brewer.pal(11, "BrBG"))(50), limits = c(-abs_max, abs_max))
myFeaturePlot(bb, "enrich_z", my.pt.size = 0.01) + scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11, "PuOr")))(50), limits = c(-abs_max, abs_max))
myFeaturePlot(bb, "enrich_z") + scale_colour_continuous_diverging( limits = c(-abs_max, abs_max))
myFeaturePlot(bb, "enrich_z") + scale_colour_continuous_divergingx( limits = c(-abs_max, abs_max))
myFeaturePlot(bb, "enrich_z", my.pt.size = 0.20) + scale_colour_viridis_c( limits = c(-abs_max, abs_max))
bb$enrich_z2 = bb$enrich_z ** 2
myFeaturePlot(bb, "enrich_z2") + scale_colour_viridis_c() + ggtitle("Z Score Squared")
myFeaturePlot(bb, "enrich_z2") + scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11, "RdBu")))(50)) + ggtitle("Z Score Squared")
myFeaturePlot(bb, "enrich_z") + scale_color_gradientn(colors = c(rep("lightgray", 750), rep("red", 250)), limits = c(-abs_max, abs_max), breaks = c(-abs_max, -abs_max/2, 0, abs_max/2, abs_max), labels = c(0, 25, 50, 75, 100))

# Sample/Subsample Stats
library("DropletUtils")
all_res = lapply(unique(bb$sample), findMeanReadsPerNuc)
all_res_df = bind_rows(all_res)
bb$reads = all_res_df$reads[match(colnames(bb), all_res_df$real)]
ggplot(bb@meta.data, aes(x=nCount_RNA, y = reads, color = sample)) + geom_point()

# Stats by Sample
sample_df = setNames(aggregate(nFeature_RNA ~ sample, bb@meta.data, mean), list("sample", "mean_nFeature_RNA"))
sample_df$median_nFeature_RNA = aggregate(nFeature_RNA ~ sample, bb@meta.data, median)[,2]
sample_df$mean_nCount_RNA = aggregate(nCount_RNA ~ sample, bb@meta.data, mean)[,2]
sample_df$mean_reads = aggregate(reads ~ sample, bb@meta.data, mean)[,2]*2
sample_df$num_reads = aggregate(reads ~ sample, bb@meta.data, sum)[,2]*2
sample_df$num_nuc = table(bb@meta.data$sample)
sample_df$genes_detected = unlist(lapply(unique(bb$sample), function(x) length(which(rowSums(bb@assays$RNA@counts[,which(bb$sample == x)]) > 0))))
write.csv(sample_df, "~/research/brain/results/sample_stat.csv")

# Stats by Subsample
subsample_df = setNames(aggregate(nFeature_RNA ~ subsample, bb@meta.data, mean), list("subsample", "mean_nFeature_RNA"))
subsample_df$median_nFeature_RNA = aggregate(nFeature_RNA ~ subsample, bb@meta.data, median)[,2]
subsample_df$mean_nCount_RNA = aggregate(nCount_RNA ~ subsample, bb@meta.data, mean)[,2]
subsample_df$mean_reads = aggregate(reads ~ subsample, bb@meta.data, mean)[,2]*2
subsample_df$num_reads = aggregate(reads ~ subsample, bb@meta.data, sum)[,2]*2
subsample_df$num_nuc = table(bb@meta.data$subsample)
subsample_df$genes_detected = unlist(lapply(unique(bb$subsample), function(x) length(which(rowSums(bb@assays$RNA@counts[,which(bb$subsample == x)]) > 0))))
write.csv(subsample_df, "~/research/brain/results/subsample_stat.csv")

ggplot(sample_df, aes(x=median_nFeature_RNA, y = mean_nCount_RNA, color = sample)) + geom_point()
ggplot(sample_df, aes(x=median_nFeature_RNA, y = genes_detected, color = sample)) + geom_point()
ggplot(subsample_df, aes(x=median_nFeature_RNA, y = mean_nCount_RNA, color = subsample)) + geom_point()
ggplot(subsample_df, aes(x=median_nFeature_RNA, y = genes_detected, color = subsample)) + geom_point()
ggplot(subsample_df, aes(x=mean_reads, y = genes_detected, color = subsample)) + geom_point()
ggplot(subsample_df, aes(x=num_nuc, y = genes_detected, color = subsample)) + geom_point()

mat = bb@assays$RNA@counts
mat[which(mat > 0)] = 1
for (sample in bb$sample) {
  num_sample_cells = length(which(bb$sample == sample))
  this_sample_min_cells = min_pct * num_sample_cells
  idx = which(rowSums(mat) >= this_sample_min_cells)
}

findMeanReadsPerNuc = function(sample) {
  print(sample)
  this_reads = read10xMolInfo(paste0("~/research/brain/data/", sample, "_molecule_info.h5"))
  this_convert = data.frame(real = colnames(bb)[which(bb$sample == sample)])
  this_convert$sub = colsplit(this_convert$real, pattern = "-", names = c("1", "2"))[,1]
  this_convert$sub_sub = colsplit(this_convert$sub, pattern = "_", names = c("1", "2"))[,2]
  this_reads_good = this_reads$data[which(this_reads$data$cell %in% this_convert$sub_sub),]
  reads_per_cell = aggregate(reads ~ cell, this_reads_good, sum)
  this_convert$reads = reads_per_cell$reads[match(this_convert$sub_sub, reads_per_cell$cell)]
  
  # num_reads = sum(this_reads_good$reads)*2
  # mean_reads = (sum(this_reads_good$reads)/nrow(this_convert))*2
  
  return(this_convert)
}

findCell = function(cells_strip, obj, sample) {
  df = data.frame()
  cells_to_search = colnames(obj)[which(obj$sample == sample)]
  
  findCellMini = function(old_cell) {
    new_cell = cells_to_search[which(grepl(old_cell, cells_to_search))]
    return(new_cell)
  }
  
  numCores = detectCores()
  new_cells = unlist(mclapply(cells_strip, findCellMini, mc.cores = numCores))
  df = data.frame(old_cell = cells_strip, new_cell = new_cells)
  
  return(df)
}

# DEG Enrichment
int_genes = read.csv("C:/Users/miles/Downloads/brain/data/markers/bb_all_interesting_genes_final_3count_051021.csv", stringsAsFactors = F)
int_genes$gene[which(int_genes$gene == "slac17a6")] = "slc17a6"
deg15 = read.csv("C:/Users/miles/Downloads/brain/results/bb/bb_all_cluster_15_degs.csv", stringsAsFactors = F)
deg53 = read.csv("C:/Users/miles/Downloads/brain/results/bb/bb_all_cluster_53_degs.csv", stringsAsFactors = F)
deg15 = deg15[which(deg15$p_val_adj < 0.05),]
deg53 = deg53[which(deg53$p_val_adj < 0.05),]

num_deg = findMarkersInDEG(deg = deg15, markers = int_genes$gene)[[2]]
pct_deg = findMarkersInDEG(deg = deg15, markers = int_genes$gene, pct = T)[[2]]
test15 = markersInDEGUpDown(deg15, int_genes$gene)
test53 = markersInDEGUpDown(deg53, int_genes$gene)

non_zero_genes = rownames(bb)[which(rowSums(bb@assays$RNA@counts > 0))]
all_clust15_deg = unique(deg15$gene[which(deg15$avg_logFC > 0)])
all_class_genes = unique(int_genes$gene)

int_clust15 = data.frame()
for (clust in 0:52) {
  clust_pos_deg = deg53[which(deg53$avg_logFC > 0 & deg53$cluster == clust),]
  # clust_pos_deg = clust_pos_deg[which(clust_pos_deg$gene %in% all_class_genes),]
  num_clust_deg = nrow(clust_pos_deg)
  for (class in unique(int_genes$class)) {
    class_genes = int_genes$gene[which(int_genes$class == class | int_genes$class2 == class)]
    class_genes = class_genes[which(class_genes %in% non_zero_genes)]
    # class_genes = class_genes[which(class_genes %in% all_clust15_deg)]
    num_class_genes = length(class_genes)
    
    # Overlap Values
    raw_ovlp = length(which(class_genes %in% clust_pos_deg$gene))
    ovlp_norm_clust = raw_ovlp / num_clust_deg
    ovlp_norm_list = raw_ovlp / num_class_genes
    ovlp_norm_clust_list = raw_ovlp / (num_clust_deg + num_class_genes)
    
    # Statistical Test
    # contig_table = data.frame(isClass = c(raw_ovlp, num_class_genes - length(raw_ovlp)), notClass = c(num_clust_deg - length(raw_ovlp), length(non_zero_genes) - length(unique(c(class_genes, clust_pos_deg$gene)))))
    # contig_table = data.frame(isClass = c(raw_ovlp, num_class_genes - length(raw_ovlp)), notClass = c(num_clust_deg - length(raw_ovlp), length(all_clust15_deg) - length(unique(c(class_genes, clust_pos_deg$gene))) ))
    contig_table = data.frame(isClass = c(raw_ovlp, num_class_genes - length(raw_ovlp)), 
                              notClass = c(num_clust_deg - length(raw_ovlp), length(all_clust15_deg) - length(unique(c(class_genes, clust_pos_deg$gene))) ))
    rownames(contig_table) = c("isDEG", "notDEG")
    p = fisher.test(contig_table, alternative = "greater")$p.value
    
    int_clust15 = rbind(int_clust15, data.frame(cluster = clust, class = class, raw_ovlp = raw_ovlp, ovlp_norm_clust = ovlp_norm_clust, ovlp_norm_list = ovlp_norm_list, ovlp_norm_clust_list = ovlp_norm_clust_list, p = p, this_clust = contig_table[1,1]/contig_table[1,2], other = contig_table[2,1]/contig_table[2,2]))
  }
}
int_clust15$cluster = factor(int_clust15$cluster, levels = 0:52)
int_clust15$bh = p.adjust(int_clust15$p, method = "BH")
int_clust15$bon = p.adjust(int_clust15$p, method = "bonferroni")

ggplot(int_clust15, aes(x = cluster, y = ovlp_norm_clust_list, fill = class)) + geom_bar(stat="identity") + ylab("Normalized Overlap") + xlab("Cluster") + ggtitle("Overlap with Interesting Genes and Cluster DEGs Normalized for Length of Both")
ggplot(int_clust15, aes(x = cluster, y = raw_ovlp, fill = class)) + geom_bar(stat="identity") + ylab("Raw Overlap") + xlab("Cluster") + ggtitle("Raw Overlap with Interesting Genes and Cluster DEGs")

big_res = data.frame()
for (i in 0:14) {
  clust_pos = deg15[which(deg15$cluster == i & deg15$avg_logFC > 0),]
  clust_pos_ranks = 1:nrow(clust_pos)
  names(clust_pos_ranks) = clust_pos$gene
  fgseaRes_pos = fgsea(setNames(list(lg11), c("lg11")), clust_pos_ranks, nperm = 1000)
  # fgseaRes_pos = fgsea(setNames(list(int_genes$gene), c("int")), clust_pos_ranks, nperm = 1000)
  fgseaRes_pos$cluster = i
  fgseaRes_pos$isPos = T
  # plotEnrichment(int_genes$gene, clust0_pos_ranks)
  
  clust_neg = deg15[which(deg15$cluster == i & deg15$avg_logFC < 0),]
  clust_neg_ranks = 1:nrow(clust_neg)
  names(clust_neg_ranks) = clust_neg$gene
  fgseaRes_neg = fgsea(setNames(list(lg11), c("lg11")), clust_neg_ranks, nperm = 1000)
  # fgseaRes_neg = fgsea(setNames(list(int_genes$gene), c("int")), clust_neg_ranks, nperm = 1000)
  fgseaRes_neg$cluster = i
  fgseaRes_neg$isPos = F
  
  big_res = rbind(big_res, fgseaRes_pos)
  big_res = rbind(big_res, fgseaRes_neg)
}


# Enrichment Testing
res = separateEnrichTest(bb, this_list)
test = res[[2]]
bb_all_enrich = xlsx::read.xlsx("C:/Users/miles/Downloads/brain/results/bb/bb_all_enrich_results.xlsx", sheetIndex = 1)
test$new = test$p_adj
test$old = bb_all_enrich$ieg
test = test[,c("cluster", "new", "old")]
test2 = melt(test, id.vars = c("cluster"))
ggplot(test2, aes(cluster, value, color = variable)) + geom_line(size = 1) + ylab("Bonferroni Adjusted p-value") + ggtitle("Enrichment Test on IEGs")

test = res[[1]]
test$p = as.numeric(as.vector(test$p))
ggplot(test, aes(cluster, p, color = gene, group = gene)) + geom_line()

# Subset Data by neuronal clusters only
bb_neuron_15 = subset(bb, cells = WhichCells(bb, idents = c(4, 9, 13), invert = T))
bb_neuron_53 = subset(bb, cells = WhichCells(bb, idents = c(5, 20, 31, 45, 46, 50), invert = T))

# Enrichment Lists
disc_ase_dig_v_build_pit_up_in_dig_111820 = read.table(paste0(rna_path, "data/markers/disc_ase_dig_v_build_pit_up_in_dig_111820.txt"), stringsAsFactors = F)[,1]
ase_down_all_111820 = read.table(paste0(rna_path, "data/markers/ase_down_all_111820.txt"), stringsAsFactors = F)[,1]
ase_up_all_111820 = read.table(paste0(rna_path, "data/markers/ase_up_all_111820.txt"), stringsAsFactors = F)[,1]
disc_ase_bhve_v_ctrl_111820 = read.table(paste0(rna_path, "data/markers/disc_ase_bhve_v_ctrl_111820.txt"), stringsAsFactors = F)[,1]
ase_same_dir_all_111820 = read.table(paste0(rna_path, "data/markers/ase_same_dir_all_111820.txt"), stringsAsFactors = F)[,1]
disc_ase_dig_v_build_111820 = read.table(paste0(rna_path, "data/markers/disc_ase_dig_v_build_111820.txt"), stringsAsFactors = F)[,1]
pc_20_rc_30_25kb_genes_binned = read.table(paste0(rna_path, "data/markers/pc_20_rc_30_25kb_genes_binned.txt"), stringsAsFactors = F)[,1]
pc_20_rc_40_25kb_genes_binned = read.table(paste0(rna_path, "data/markers/pc_20_rc_40_25kb_genes_binned.txt"), stringsAsFactors = F)[,1]
LG11_highFST_genes_umd2a_120920 =  read.table(paste0(rna_path, "data/markers/LG11_highFST_genes_umd2a_120920.txt"), stringsAsFactors = F)[,1]

my_list = setNames(list(disc_ase_dig_v_build_pit_up_in_dig_111820, ase_down_all_111820, ase_up_all_111820, disc_ase_bhve_v_ctrl_111820, ase_same_dir_all_111820, disc_ase_dig_v_build_111820, pc_20_rc_30_25kb_genes_binned, pc_20_rc_40_25kb_genes_binned, LG11_highFST_genes_umd2a_120920), 
                   c("disc_ase_dig_v_build_pit_up_in_dig_111820", "ase_down_all_111820", "ase_up_all_111820", "disc_ase_bhve_v_ctrl_111820", "ase_same_dir_all_111820", "disc_ase_dig_v_build_111820", "pc_20_rc_30_25kb_genes", "pc_20_rc_40_25kb_genes", "LG11_highFST_genes_umd2a_120920"))

# Plot UMAPs of lists colored by positive effect size
for (i in 1:length(my_list)) {
  base_name = names(my_list[i])
  Idents(bb) = bb$seuratclusters15
  res15 = markerExpPerCellPerClusterQuick(bb, my_list[[i]])
  
  res15[[3]]$color = res15[[3]]$mag_pos
  v_cols = viridis(4)
  res15[[3]]$color = plyr::revalue(as.character(res15[[3]]$color), replace = c("negligible" = v_cols[1], "small" = v_cols[2], "medium" = v_cols[3], large = v_cols[4]))
  p15 = DimPlot(bb, label = T, pt.size = 0.7) + scale_color_manual(values = res15[[3]]$color) + NoLegend()
  
  Idents(bb) = bb$seuratclusters53
  res53 = markerExpPerCellPerClusterQuick(bb, my_list[[i]])
  
  res53[[3]]$color = res53[[3]]$mag_pos
  v_cols = viridis(4)
  res53[[3]]$color = plyr::revalue(as.character(res53[[3]]$color), replace = c("negligible" = v_cols[1], "small" = v_cols[2], "medium" = v_cols[3], large = v_cols[4]))
  p53 = DimPlot(bb, label = T, pt.size = 0.7) + scale_color_manual(values = res53[[3]]$color) + NoLegend()
  
  png_name = paste0("C:/Users/miles/Downloads/brain/results/bb/enrichment/", base_name, "_15_umap.png")
  png(file = png_name, width = 1000, height = 1000, res = 200)
  print(p15)
  dev.off()
  
  png_name = paste0("C:/Users/miles/Downloads/brain/results/bb/enrichment/", base_name, "_53_umap.png")
  png(file = png_name, width = 1000, height = 1000, res = 200)
  print(p53)
  dev.off()
}

# Create BarPlots of # and % of Lists in DEGs per cluster
for (i in 1:length(my_list)) {
  base_name = names(my_list[i])
  # res = markerExpPerCellPerCluster(bb_neuron_15, my_list[[i]], correct = T)
  res15 = markersInDEGUpDown(deg15, my_list[[i]])
  res15_pct = markersInDEGUpDown(deg15, my_list[[i]], pct=T)
  res53 = markersInDEGUpDown(deg53, my_list[[i]])
  res53_pct = markersInDEGUpDown(deg53, my_list[[i]], pct=T)
  
  png_name = paste0("C:/Users/miles/Downloads/brain/results/bb/enrichment/", base_name, "_15_deg.png")
  png(file = png_name, width = 1500, height = 1000, res = 200)
  print(res15[[2]])
  dev.off()
  
  png_name = paste0("C:/Users/miles/Downloads/brain/results/bb/enrichment/", base_name, "_15_deg_logFC.png")
  png(file = png_name, width = 1500, height = 1000, res = 200)
  print(res15[[3]])
  dev.off()
  
  png_name = paste0("C:/Users/miles/Downloads/brain/results/bb/enrichment/", base_name, "_15_deg_pct.png")
  png(file = png_name, width = 1500, height = 1000, res = 200)
  print(res15_pct[[2]])
  dev.off()
  
  # res = markerExpPerCellPerCluster(bb_neuron_53, my_list[[i]], correct = T)
  png_name = paste0("C:/Users/miles/Downloads/brain/results/bb/enrichment/", base_name, "_53_deg.png")
  png(file = png_name, width = 2500, height = 1000, res = 200)
  print(res53[[2]])
  dev.off()
  
  png_name = paste0("C:/Users/miles/Downloads/brain/results/bb/enrichment/", base_name, "_53_deg_logFC.png")
  png(file = png_name, width = 2500, height = 1000, res = 200)
  print(res53[[3]])
  dev.off()
  
  png_name = paste0("C:/Users/miles/Downloads/brain/results/bb/enrichment/", base_name, "_53_deg_pct.png")
  png(file = png_name, width = 2500, height = 1000, res = 200)
  print(res53_pct[[2]])
  dev.off()
}

# For the Master Excel Sheet
levels15 = apply(rev(expand.grid(unique(bb$cond), levels(bb$seuratclusters15))), 1, paste, collapse=" ")
levels53 = apply(rev(expand.grid(unique(bb$cond), levels(bb$seuratclusters53))), 1, paste, collapse=" ")
df15_p = data.frame(cluster_cond = levels15)
df15_e = data.frame(cluster_cond = levels15)
df53_p = data.frame(cluster_cond = levels53)
df53_e = data.frame(cluster_cond = levels53)
df15_bon = data.frame(cluster_cond = levels15)
df53_bon = data.frame(cluster_cond = levels53)
my_list = c("ase_down_all_111820", "ase_up_all_111820", "ase_same_dir_all_111820", "disc_ase_dig_v_build_111820", "disc_ase_bhve_v_ctrl_111820", "disc_ase_dig_v_build_pit_up_in_dig_111820", "LG11_highFST_genes_umd2a_120920", "pc_20_rc_20_25kb_genes", "pc_20_rc_30_25kb_genes", "pc_20_rc_40_25kb_genes", "pc_30_rc_30_25kb_genes", "pc_70_rc_70", "pc_80_rc_80", "ieg_like_011521")
for (i in 1:length(my_list)) {
  base_name = my_list[[i]]
  this_list = read.table(paste0("C:/Users/miles/Downloads/brain/data/markers/", base_name, ".txt"), header = F, stringsAsFactors = F)[,1]
  
  bb$seurat_clusters = bb$seuratclusters15
  Idents(bb) = bb$seurat_clusters
  res15 = markerExpPerCellPerClusterBVC(bb, this_list, correct =T)[[3]]
  df15_p[,c(paste(i))] = res15$p
  df15_e[,c(paste(i))] = res15$d
  df15_bon[,c(paste(i))] = p.adjust(res15$p, method = "bonferroni")
  
  bb$seurat_clusters = bb$seuratclusters53
  Idents(bb) = bb$seurat_clusters
  res53 = markerExpPerCellPerClusterBVC(bb, this_list, correct =T)[[3]]
  df53_p[,c(paste(i))] = res53$p
  df53_e[,c(paste(i))] = res53$d
  df53_bon[,c(paste(i))] = p.adjust(res53$p, method = "bonferroni")
}
colnames(df15_p)[2:ncol(df15_p)] = my_list
colnames(df15_e)[2:ncol(df15_e)] = my_list
colnames(df53_p)[2:ncol(df53_p)] = my_list
colnames(df53_e)[2:ncol(df53_e)] = my_list
colnames(df15_bon)[2:ncol(df15_bon)] = my_list
colnames(df53_bon)[2:ncol(df53_bon)] = my_list
df15_e_dif = df15_e[c(TRUE, FALSE), 2:ncol(df53_e)] - df15_e[c(FALSE, TRUE), 2:ncol(df53_e)]
df53_e_dif = df53_e[c(TRUE, FALSE), 2:ncol(df53_e)] - df53_e[c(FALSE, TRUE), 2:ncol(df53_e)]
write.csv(df15_p, "C:/Users/miles/Downloads/df15_p.csv")
write.csv(df15_e, "C:/Users/miles/Downloads/df15_e.csv")
write.csv(df53_p, "C:/Users/miles/Downloads/df53_p.csv")
write.csv(df53_e, "C:/Users/miles/Downloads/df53_e.csv")
write.csv(df53_e_dif, "C:/Users/miles/Downloads/df53_e_dif.csv")
write.csv(df15_e_dif, "C:/Users/miles/Downloads/df15_e_dif.csv")
write.csv(df15_bon, "C:/Users/miles/Downloads/df15_bon.csv")
write.csv(df53_bon, "C:/Users/miles/Downloads/df53_bon.csv")

# Add the bonferronni correction to the excel sheets
test = readxl::read_excel("C:/Users/miles/Downloads/brain/results/bb/bb_all_enrich_results.xlsx", "pvalue 15")
test2 = data.frame(lapply(2:ncol(test), function(x) p.adjust(test[[x]], method = "bonferroni")))
write.csv(test2, "C:/Users/miles/Downloads/df15_bon.csv")

# Are Samples Disproportionately Represented by cluster?
df = data.frame(Cluster = bb$seuratclusters53, Replicate=bb$sample)
col_pal = c("#9d0208", "#d00000", "#dc2f02", "#e85d04", "#f48c06", "#03045e", "#023e8a", "#0077b6", "#0096c7", "#00b4d8")
ggplot(df, aes(x = Cluster, fill = Replicate, cond = Replicate)) + geom_bar(alpha = 0.8) + scale_color_manual(values = col_pal) + scale_fill_manual(values = col_pal) + ggtitle("Number of Cells per Cluster per Replicate")

# Are Samples Disproportionately Represented by cluster?
df2 = data.frame()
df3 = data.frame()
df4 = data.frame(Replicate = unique(bb$sample), row.names = unique(bb$sample))
for (cluster in levels(bb$seuratclusters53)) {
  cluster_cells = colnames(bb)[which(bb$seuratclusters53 == cluster)]
  for (sample in unique(bb$sample)) {
    sample_cells = colnames(bb)[which(bb$sample == sample)]
    cluster_sample = cluster_cells[which( cluster_cells %in% sample_cells )]
    pct = (length(cluster_sample)/length(cluster_cells)) * 100
    df2 = rbind(df2, t(c(cluster, sample, length(cluster_sample), pct)))
    df3 = rbind(df3, t(c(sample, cluster, length(cluster_sample), length(cluster_sample)/length(sample_cells) * 100)))
  }
  df4 = cbind(df4, c( df3$V4[which(df3$V2 == cluster)] ))
}
colnames(df2) = c("Cluster", "Replicate", "Num_cells", "Percent")
colnames(df3) = c("Replicate", "Cluster", "Num_cells", "Percent")
df4$Replicate = NULL
colnames(df4) = levels(bb$seuratclusters53)
rownames(df4) = unique(bb$sample)
df2$Percent = as.numeric(as.vector(df2$Percent))
ggplot(df2, aes(x = Cluster, y = Percent, fill = Replicate, color = Replicate)) + geom_bar(stat = "identity") + scale_color_manual(values = col_pal) + scale_fill_manual(values = col_pal) + ggtitle("Composition of Clusters by Replicate")

# Are Samples Disproportionately Represented by cluster Across pairs?
df2$pair = str_sub(as.character(df2$Replicate),-1,-1)
pdf("C:/Users/miles/Downloads/bvc_cluster_num_cells.pdf", width = 20, height = 16)
print(ggpaired(df2, x = "cond", y = "Percent", id = "pair", line.color = "gray", palette = "jco", color = "cond", facet.by = "Cluster", line.size = 0.4) + stat_compare_means(label = "p.format", paired = TRUE, method = "t.test"))
dev.off()

test = rowSums(mat)/ncol(bb)
hist(test, breaks  = 100, xlab = "Percent of Cells a Gene is Expressed In", main = "Percent of Cells a Gene is Expressed In", xaxt="n")
axis(side=1, at=seq(0, 1, .02), labels=seq(0,1,.02))

Idents(bb) = bb$seuratclusters53
clust13_cells = WhichCells(bb, idents = 13)
# zack_0 = zack15$mzebra[which(zack15$cluster == 0)]
gene_df2 = gene_df[order(gene_df$avg_logFC, decreasing = T),]
zack_0 = as.vector(gene_df2$clust_non_zero[which(gene_df2$pct_dif > 5)][1:20])
for (gene in zack_0) {
  pdf(paste0("C:/Users/miles/Downloads/brain/results/bb/paintings/53_13/", gene, "_cluster13.pdf"), width = 8, height = 3.5)
  print(bvcVis(bb, gene, cells.use = clust13_cells, only.pos = T, mode = "violin_split", cell.alpha = 0.3))
  dev.off()
}
test = zack15[which(zack15$cluster == 0),]
test = test[order(test$abs_avg_logFC, decreasing = T),]

res = data.frame()
for (cluster in 0:14) {
  df4 = df3[which(df3$Cluster == cluster),]
  clust13_counts = as.numeric(as.vector(df4$Num_cells))
  df4$pair = str_sub(as.character(df4$Replicate),-1,-1)
  df4$cond = substr(as.character(df4$Replicate),1, 1)
  contingency_table <- data.frame(count = c(clust13_counts, as.vector(table(bb$sample)) - clust13_counts ),
                                  target = c(rep(1, 10), rep(0, 10)),
                                  cond = c(rep(c(rep(1, 5), rep(0, 5)), 2)),
                                  pair = c(rep(1:5, 4)))
  t2 = xtabs(count ~ target + cond + pair, data=contingency_table)
  test_res = mantelhaen.test(t2)
  res = rbind(res, t(c(cluster, df4$Num_cells, df4$Percent, test_res$statistic, test_res$p.value)))
}
colnames(res) = c("cluster", paste0(unique(bb$sample), "_num"), paste0(unique(bb$sample), "_pct"), "test_stat", "p")
res$bh = p.adjust(as.numeric(as.vector(res$p)), method = "BH")
res$bon = p.adjust(as.numeric(as.vector(res$p)), method = "bonferroni")
write.csv(res, "C:/Users/miles/Downloads/nuclei_count_by_cluster15_and_condition_12821.csv", row.names = F)

df3$pair = str_sub(as.character(df3$Replicate),-1,-1)
df3$cond = substr(as.character(df3$Replicate),1, 1)
t3 = xtabs(Num_cells ~ cond + pair + Cluster, data = df3)
mantelhaen.test(t3)

t.test(df2_13$Percent[1:5], df2_13$Percent[6:10], paired = T)

# Create a Dataframe for ML that has samples as rows and metadata as columns
df_meta = data.frame()
for ( sample in unique(bb$sample) ) {
  sample_cells = colnames(bb)[which( bb$sample == sample )]
  metasums = colSums(bb@meta.data[sample_cells,c("nCount_RNA", "nFeature_RNA", "pct_mt")])
  df_meta = rbind(df_meta, metasums)
}
rownames(df_meta) = unique(bb$sample)
colnames(df_meta) = c("nCount_RNA", "nFeature_RNA", "pct_mt")

pdf("C:/Users/miles/Downloads/zack_num.pdf", width = 20, height = 16)
print(ggpaired(df3, x = "cond", y = "Num_cells", id = "pair", line.color = "gray", palette = "jco", color = "cond", facet.by = "Cluster", line.size = 0.4) + stat_compare_means(label = "p.format", paired = TRUE, method = "t.test"))
dev.off()

# Enrichment Test on List that Zack made that's a combination of high FST LG 11, pc_rc, and ASE
Idents(bb) = bb$seuratclusters15
res = markerExpPerCellPerCluster(bb, z_list2)
png_name = paste0("C:/Users/miles/Downloads/brain/results/bb/enrichment/z2_pcrc_lg11_15.png")
png(file = png_name, width = 1500, height = 1000, res = 200)
print(res[[1]])
dev.off()
write.csv(res[[3]], "C:/Users/miles/Downloads/brain/results/bb/enrichment/z2_pcrc_lg11_15.csv")

png_name = paste0("C:/Users/miles/Downloads/brain/results/bb/enrichment/z2_pcrc_lg11_53.png")
png(file = png_name, width = 2500, height = 1000, res = 200)
print(res[[1]])
dev.off()
write.csv(res[[3]], "C:/Users/miles/Downloads/brain/results/bb/enrichment/z2_pcrc_lg11_53.csv")

# BHVE vs CTRL plots
Idents(bb) = bb$seuratclusters15
pdf(paste0("C:/Users/miles/Downloads/brain/results/bb/paintings/LOC101468629_box.pdf"), width = 8, height = 3.5)
print(bvcVis(bb, "LOC101468629", only.pos = T, mode = "box", cell.alpha = 0.3, add_title = "- Only Positive Cells"))
dev.off()

df = data.frame(cell = colnames(bb), score = bb$ieg_like_score, cond = bb$cond, cluster15 = bb$seuratclusters15, cluster53 = bb$seuratclusters53)
mean15 = data.frame(cluster15 = rep(levels(bb$seuratclusters15),2),
                    means15 = c(sapply(levels(bb$seuratclusters15), function(x) mean(df$score[which(df$cluster15 == x & df$cond == "BHVE")])), sapply(levels(bb$seuratclusters15), function(x) mean(df$score[which(df$cluster15 == x & df$cond == "CTRL")]))),
                    cond = c(rep("BHVE", length(levels(bb$seuratclusters15))), rep("CTRL", length(levels(bb$seuratclusters15)))))
mean53 = data.frame(cluster53 = levels(bb$seuratclusters53),
                    means53 = sapply(levels(bb$seuratclusters53), function(x) mean(df$score[which(df$cluster53 == x)])))
ggplot(df, aes(x = cluster15, y = score, fill = cond, color = cond)) + geom_violin(alpha = 0.6) + geom_jitter(position=position_dodge2(width = 0.6), alpha = 0.01)
ggplot(df, aes(x = cluster15, y = score, fill = cond, color = cond)) + geom_point(position = position_dodge2(width = 0.8), alpha = 0.05) + geom_bar(data = mean15, aes(x = cluster15, y = means15, fill = cond, color = cond), stat = "identity", alpha = 0.3, position = position_dodge2(width = 0.5))

png_name = paste0("C:/Users/miles/Downloads/brain/results/bb/enrichment/neurogenesis_53_bvc.png")
png(file = png_name, width = 3000, height = 1000, res = 200)
print(res[[1]])
dev.off()

Idents(bb) = bb$seuratclusters15
png_name = paste0("C:/Users/miles/Downloads/brain/results/bb/paintings/rsrp1_15_bvc.png")
png(file = png_name, width = 3000, height = 1400, res = 200)
myFeaturePlot(bb, feature = "rsrp1", my.split.by = "cond", my.col.pal = pal)
dev.off()

Idents(bb) = bb$seuratclusters15
png_name = paste0("C:/Users/miles/Downloads/brain/results/bb/paintings/rsrp1_15_0_bvc.png")
png(file = png_name, width = 2000, height = 700, res = 200)
myFeaturePlot(bb, feature = "rsrp1", my.split.by = "cond", my.col.pal = pal, cells.use = colnames(bb)[which(bb$seuratclusters15 == 0)])
dev.off()

zack15 = read.csv("C:/Users/miles/Downloads/bb15_combined_all_pairs_strict_hits_120920_pct_logFC.csv", stringsAsFactors = F)
for (i in 1:nrow(zack15)) {
  gene = zack15$mzebra[i]
  human = zack15$human[i]
  cluster = zack15$cluster[i]
  upBhve = zack15$z[i] < 0
  subtitle_str = paste0(gene, " (", human, ") in Cluster ", cluster, "\n",
                        "Expression Significantly Greater in ")
  subtitle_str = ifelse(upBhve, paste0(subtitle_str, "BHVE"), paste0(subtitle_str, "CTRL"))
  
  Idents(bb) = bb$seuratclusters15
  png_name = paste0("C:/Users/miles/Downloads/brain/results/bb/paintings/bvc15/", gene, "_", cluster, ".png")
  png(file = png_name, width = 1800, height = 1500, res = 200)
  p_list = list()
  p_list[[1]] = myFeaturePlot(bb, feature = gene, my.pt.size = 1.75, my.split.by = "cond", my.col.pal = pal, cells.use = colnames(bb)[which(bb$seuratclusters15 == cluster)], na.blank = T)
  # p_list[[2]] = bvcVis(bb, gene, cells.use = WhichCells(bb, idents = cluster), meta = "cond", mode = "violin") + labs(title = NULL)
  p_list[[2]] = bvcVis(bb, gene, cells.use = WhichCells(bb, idents = cluster), mode = "violin_split") + labs(title = NULL) + xlab("")
  p = plot_grid(plotlist=p_list, ncol = 1)
  my.title.obj = ggdraw() + draw_label(subtitle_str)
  p1 = plot_grid(my.title.obj, p, ncol = 1, rel_heights=c(0.1, 1))
  print(p1)
  dev.off()
}

# Old Attempt at Specific Markers
deg15 = read.csv("C:/Users/miles/Downloads/brain/results/bb/bb_all_sig_cluster_15_degs.csv")
deg15$pct_dif = deg15$pct.1 - deg15$pct.2
rank = c()
for (i in 0:14) {
  this_clust = deg15[which(deg15$cluster == i),]
  rank = c(rank, 1:nrow(this_clust))
}
deg15$rank = rank
deg15_order = deg15[order(deg15$pct.2),]
deg15_thresh = deg15[which(deg15$pct.1 > 0.1 & deg15$pct.2 < 0.02),]
deg15_thresh$hgnc = gene_info$human[match(as.character(as.vector(deg15_thresh$gene)), gene_info$mzebra)]
# deg15_order = deg15[which(deg15$avg_logFC > 0.5),]
# top10 = t(as.data.frame(lapply(0:14, function(x) t(deg15_order[which(deg15_order$cluster == x),][1:10,]))))
top10 = t(as.data.frame(lapply(0:14, function(x) t(deg15_thresh[which(deg15_thresh$cluster == x),][1:10,]))))
top10 = as.data.frame(top10)
top10$gene = as.character(as.vector(top10$gene))
top10 = top10[which( ! is.na(top10$gene) ),]
for (i in 1:nrow(top10)) {
  clust = as.numeric(as.vector(top10$cluster[i]))
  gene = top10$gene[i]
  this_rank = as.numeric(as.vector(top10$rank[i]))
  png(paste0(rna_path, "/results/bb/paintings/15_markers/", clust, "_", gene, ".png"), width = 750, height = 750, res = 100)
  print(FeaturePlot(bb, gene, order = T, pt.size = 1, label = T) + ggtitle(paste(gene, "Hit For", clust, "- Rank", this_rank)))
  dev.off()
}
write.csv(deg15_thresh, "C:/Users/miles/Downloads/specific_markers_15_10_02_12821.csv")

bb$my_clust = bb$seuratclusters15
bb$my_clust[which(! bb$my_clust %in% c(0, 1, 2, 5, 7, 10))] = NA
# bb$my_clust = bb$seuratclusters53
# bb$my_clust[which(! bb$my_clust %in% c(0, 1, 2, 7, 8, 10, 12, 13, 14, 16, 17, 18, 22, 27, 28, 29, 30, 31, 33, 36, 39, 40, 44, 45, 46, 48))] = NA
Idents(bb) = bb$my_clust
DimPlot(bb, label = T) + ggtitle("Clusters that Predicted All Samples")

replace_vect = c(63, 320, 67, NA, NA, NA, NA, 51, 109, NA, 3, NA, 582, 68, 174, NA, 820, 3, 867, NA, NA, NA, 7, NA, NA, NA, NA, 479, 67, 3, 854, 75, NA, 16, NA, NA, 364, NA, NA, 717, 48, NA, NA, NA, 296, 4, 10, NA, 9, NA, NA, NA, NA)
names(replace_vect) = levels(bb$seuratclusters53)
bb$my_clust = plyr::revalue(as.character(bb$seuratclusters53), replace = replace_vect)
# bb$my_clust = plyr::revalue(as.character(bb$seuratclusters53), replace = c("0" = 63, "1" = 320, "2" = 67, "3" = NA, "4" = NA, "5" = NA, "6" = NA, "7" = 51, "8" = 109, "9" = NA, "10" = 3, "11" = NA, "12" = 582, "13" = 68, "14" = 174, "15" = NA, "16" = 820, "17" = 3, "18" = 867, "19" = NA, "20" = NA, "21" = NA, "22" = 7, "23" = NA, "24" = NA, "25" = NA, "26" = NA, "27" = 479, "28" = 67, "29" = 3, "30" = 854, "31" = 75, "32" = NA, "33" = 16, "34" = NA, "35" = NA, "36" = 364, "37" = NA, "38" = NA, "39" = 717, "40" = 48, "41" = NA, "42" = NA, "43" = NA, 44 = 296, "45" = 4, "46" = 10, "47" = NA, "48" = 9, "49" = NA, "50" = NA, "51" = NA, "52" = NA))
bb$my_clust = as.numeric(as.vector(bb$my_clust))
FeaturePlot(bb, "my_clust", label = T) + scale_colour_viridis(direction = -1, name = "Min. # of Genes") + ggtitle("Clusters that Predicted All Samples")

# 
Idents(bb) = paste0(bb$subsample, ".", bb$seuratclusters15)
res = myAverageExpression(bb, slot="counts")
res2 = data.frame()
for (subsample in levels(bb$subsample)) {
  subsample_df = res[,which( startsWith(colnames(res), subsample) )]
  subsample_df = unlist(subsample_df[,paste0(subsample, ".", 0:12)])
  res2 = rbind(res2, subsample_df)
}
rownames(res2) = levels(bb$subsample)
colnames(res2) = as.vector(outer(rownames(res), 0:12, paste, sep="."))
write.table(res2, "~/scratch/brain/data/bb_subsample_counts_15.txt")

# Boxplots of IEG Module Expression by Condition per Cluster
ieg_module = as.vector(read.csv("C:/Users/miles/Downloads/wgcna_module_15_ieg_03012021.csv", header = T)[,3])
bb$ieg_module = colSums(bb@assays$RNA@data[ieg_module,])
ieg_df = data.frame()
for (cluster in 0:52) {
  for (cond in c("BHVE", "CTRL")) {
    score = bb$ieg_module[which(bb$seuratclusters53 == cluster & bb$cond == cond)]
    newRow = data.frame(cluster=cluster, cond=cond, score = score)
    ieg_df = rbind(ieg_df, newRow)
  }
}

pal = colorRampPalette(rev(brewer.pal(11,"Spectral")))
ieg_df$cluster = factor(ieg_df$cluster, levels=0:52)
ggplot(ieg_df, aes(cluster, score, fill = cond, color = score, group = cluster)) + geom_boxplot(alpha = 0.5) + geom_jitter(position=position_dodge2(width = 0.6), alpha = 0.01) + scale_fill_manual(values = c("white", "gray47")) + scale_color_gradientn(colors = pal(50)) + ylab("IEG Module Score per Cell") + ggtitle("IEG Module Score in Cluster 21 BHVE vs CTRL")
# ggplot(ieg_df, aes(x=Cluster, y=Score, color = Score, fill = Condition)) + geom_boxplot(alpha=0.6) + geom_jitter(position=position_dodge2(width = 0.6), alpha = 0.01) + scale_fill_manual(values = c("white", "gray47")) + scale_color_gradientn(colors = pal(50))

png("C:/Users/miles/Downloads/brain/results/bb/ieg_module_score_by_cluster15.png", width = 1000, height = 600)
print(ggplot(ieg_df, aes(cluster, score, fill = cond, color = score)) + geom_boxplot(alpha = 0.5) + geom_jitter(position=position_dodge2(width = 0.6), alpha = 0.05) + scale_fill_manual(values = c("white", "gray47")) + scale_color_gradientn(colors = pal(50)) + ylab("IEG Module Score per Cell") + ggtitle("IEG Module Score in All Clusters BHVE vs CTRL"))
dev.off()

png("C:/Users/miles/Downloads/brain/results/bb/ieg_like_score_by_cluster15_2.png", width = 400, height = 600)
print(ggplot(test2, aes(x = Condition, y=Score, color = Score, fill = Condition)) + geom_boxplot(alpha=0.6) + geom_jitter(position=position_dodge2(width = 0.6), alpha = 0.05) + scale_fill_manual(values = c("white", "gray47")) + scale_color_gradientn(colors = pal(50)))
dev.off()

test = ieg_df[which(ieg_df$cluster == 21),]
ggplot(test, aes(cond, score, fill = cond, color = cond)) + geom_boxplot(alpha = 0.5) + geom_jitter(position=position_dodge2(width = 0.6), alpha = 0.15) + ylab("IEG Module Score per Cell") + ggtitle("IEG Module Score in Cluster 21 (53 Level) BHVE vs CTRL")

# Testing Multiple Feature Plots
specific_2 = read.csv("C:/Users/miles/Downloads/brain/results/bb/specific_markers_unique_in_subclusters_from_2_53_top5_03192021.csv")
png("C:/Users/miles/Downloads/specific_2_double.png", width = 1000, height = 1000, res = 150)
multiFeaturePlot(bb, c("rgs11", "itgb8", "gli3", "LOC101481488", "LOC101472921", "LOC101474403", "LOC101484334", "LOC101471543"))
dev.off()
png("C:/Users/miles/Downloads/specific_2_mix.png", width = 1000, height = 1000, res = 150)
multiFeaturePlot2(bb, c("rgs11", "itgb8", "gli3", "LOC101481488", "LOC101472921", "LOC101474403", "LOC101484334", "LOC101471543"))
dev.off()
png("C:/Users/miles/Downloads/test_2_double.png", width = 1000, height = 1000, res = 150)
multiFeaturePlot(bb,  c("egr1", "bdnf"))
dev.off()
png("C:/Users/miles/Downloads/test_2_mix.png", width = 1000, height = 1000, res = 150)
multiFeaturePlot2(bb,  c("egr1", "bdnf"))
dev.off()
png("C:/Users/miles/Downloads/ex_in_2_double.png", width = 1000, height = 1000, res = 150)
multiFeaturePlot(bb,  c("slc17a6", "gad2"))
dev.off()
png("C:/Users/miles/Downloads/ex_in_2_mix.png", width = 1000, height = 1000, res = 150)
multiFeaturePlot2(bb,  c("slc17a6", "gad2"))
dev.off()

# Make a Matrix of Sample by Cluster Proportion for ML
prop_df = data.frame(sample = as.vector(unique(sort(bb$sample))))
for (cluster in levels(bb$seuratclusters53)) {
  cluster_sizes = c()
  for (sample in prop_df$sample) {
    cluster_sizes = c(cluster_sizes, length(which(bb$seuratclusters53 == cluster & bb$sample == sample)))
  }
  prop_df[,cluster] = cluster_sizes
}
prop_df_pct = lapply(1:nrow(prop_df), function(x) prop_df[x,2:ncol(prop_df)]/sum(prop_df[x,2:ncol(prop_df)]))
prop_df_pct = matrix(unlist(prop_df_pct), nrow = nrow(prop_df), byrow = T, dimnames = list(prop_df$sample, colnames(prop_df)[2:ncol(prop_df)]))
write.table(prop_df_pct, "~/scratch/brain/data/bb_sample_cluster_pct.txt", row.names = T, quote = F)

# Make a Matrix of Subsample by Cluster Proportion for ML
prop_df2 = data.frame(subsample = as.vector(unique(sort(bb$subsample))))
for (cluster in levels(bb$seuratclusters53)) {
  cluster_sizes = c()
  for (subsample in prop_df2$subsample) {
    cluster_sizes = c(cluster_sizes, length(which(bb$seuratclusters53 == cluster & bb$subsample == subsample)))
  }
  prop_df2[,cluster] = cluster_sizes
}
prop_df_pct2 = lapply(1:nrow(prop_df2), function(x) prop_df2[x,2:ncol(prop_df2)]/sum(prop_df2[x,2:ncol(prop_df2)]))
prop_df_pct2 = matrix(unlist(prop_df_pct2), nrow = nrow(prop_df2), byrow = T, dimnames = list(prop_df2$subsample, colnames(prop_df2)[2:ncol(prop_df2)]))
write.table(prop_df_pct2, "~/scratch/brain/data/bb_subsample_cluster_pct.txt", row.names = T, quote = F)

# Heatmap of Cor in Interesting genes.
all_cor = readRDS("~/scratch/brain/data/mat_data_cor.RDS")
b_cor = readRDS("~/scratch/brain/data/bb_b_cor.RDS")
c_cor = readRDS("~/scratch/brain/data/bb_c_cor.RDS")
genes_df = read.csv("~/scratch/brain/results/bb_all_interesting_genes_final_3count_051021.csv", stringsAsFactors = F)
genes_df$gene[which(genes_df$gene == "slac17a6")] = "slc17a6"

p_all_cor = all_cor[genes_df$gene, genes_df$gene]
diag(p_all_cor) = 0
p_b_cor = b_cor[genes_df$gene, genes_df$gene]
p_c_cor = c_cor[genes_df$gene, genes_df$gene]
p_dif_cor = p_b_cor - p_c_cor
# pheatmap::pheatmap(p_all_cor, cellwidth = 15, cellheight = 15, cluster_rows = T, cluster_cols = F, border_color = NA, angle_col = "315", filename = "~/scratch/brain/results/bb_interesting_genes_cor_all.pdf")
png("~/scratch/brain/results/bb_interesting_genes_cor_all.png", width = 20000, height = 20000)
pheatmap::pheatmap(p_all_cor, cellwidth = 15, cellheight = 15, cluster_rows = T, cluster_cols = F, border_color = NA, angle_col = "315")
dev.off()
system(paste0("rclone copy ~/scratch/brain/results/bb_interesting_genes_cor_all.png dropbox:BioSci-Streelman/George/Brain/bb/results/coexp/"))

df = aggregate(seuratclusters15 ~ subsample, bb@meta.data, length)
df$sample = substr(df$subsample, 1, 2)
df$num = substr(df$subsample, 4, 5)
pal = c("#9d0208", "#dc2f02", "#e85d04", "#f48c06")
ggplot(df, aes(x = sample, y = seuratclusters15, fill = num, color = num)) + geom_bar(stat = "identity") + xlab("") + ylab("") + scale_fill_manual(values = pal) + scale_color_manual(values = pal) + theme_classic() + scale_y_continuous(expand = c(0,0))

clust_nums = aggregate(nCount_RNA ~ seuratclusters15, bb@meta.data, length)
sample_nums = aggregate(nCount_RNA ~ sample, bb@meta.data, length)
df = aggregate(nCount_RNA ~ subsample + seuratclusters15, bb@meta.data, length)
df$sample = substr(df$subsample, 1, 2)
df$num = substr(df$subsample, 4, 5)
df$clust_nums = clust_nums$nCount_RNA[match(df$seuratclusters15, clust_nums$seuratclusters15)]
df$sample_nums = sample_nums$nCount_RNA[match(df$sample, sample_nums$sample)]
df$prop = df$nCount_RNA / df$clust_nums * 100
df$prop2 = df$nCount_RNA / df$sample_nums * 100
df$good_names = factor(convert15$new.full[match(df$seuratclusters15, convert15$old)], levels = convert15$new.full)
df$col = factor(convert15$col[match(df$seuratclusters15, convert15$old)], levels = convert15$col)
pal = c("#9d0208", "#d00000", "#dc2f02", "#e85d04", "#f48c06", "#03045e", "#023e8a", "#0077b6", "#0096c7", "#00b4d8")
ggplot(df, aes(x = seuratclusters15, y = nCount_RNA, fill = sample, color = sample)) + geom_bar(stat = "identity") + xlab("") + ylab("") + scale_fill_manual(values = pal) + scale_color_manual(values = pal) + theme_classic() + scale_y_continuous(expand = c(0,0))
ggplot(df, aes(x = good_names, y = nCount_RNA, fill = sample, color = sample)) + geom_bar(stat = "identity") + xlab("") + ylab("") + scale_fill_manual(values = pal) + scale_color_manual(values = pal) + theme_classic() + scale_y_continuous(expand = c(0,0))

ggplot(df, aes(x = seuratclusters15, y = prop, fill = sample, color = sample)) + geom_bar(stat = "identity") + xlab("") + ylab("Cluster Proportion") + scale_fill_manual(values = pal) + scale_color_manual(values = pal) + theme_classic() + scale_y_continuous(expand = c(0,0))
ggplot(df, aes(x = good_names, y = prop, fill = sample, color = sample)) + geom_bar(stat = "identity") + xlab("") + ylab("Cluster Proportion") + scale_fill_manual(values = pal) + scale_color_manual(values = pal) + theme_classic() + scale_y_continuous(expand = c(0,0))
ggplot(df, aes(x = sample, y = prop2, fill = col, color = col)) + geom_bar(stat = "identity") + xlab("") + ylab("Cluster Proportion") + theme_classic() + scale_y_continuous(expand = c(0,0)) + scale_color_identity() + scale_fill_identity()


ggplot(df, aes(x = good_names, y = prop, fill = subsample, color = subsample)) + geom_bar(stat = "identity") + xlab("") + ylab("Cluster Proportion") + theme_classic() + scale_y_continuous(expand = c(0,0))
ggplot(df, aes(x = good_names, y = nCount_RNA, fill = subsample, color = subsample)) + geom_bar(stat = "identity") + xlab("") + ylab("") + theme_classic() + scale_y_continuous(expand = c(0,0))



#==========================================================================================
# BHVE v CTRL
#==========================================================================================
gene_names <- rownames(bb)[which(rowSums(as.matrix(bb@assays$RNA@counts)) > 0)]
df = data.frame()
bhve_cells = colnames(bb)[which(bb$cond == "BHVE")]
ctrl_cells = colnames(bb)[which(bb$cond == "CTRL")]
for (gene in gene_names) {
  pos_cells = colnames(bb)[which(bb@assays$RNA@counts[gene,] != 0)]
  this_b_cells = length(which(pos_cells %in% bhve_cells))
  this_c_cells = length(which(pos_cells %in% ctrl_cells))
  df = rbind(df, t(c(gene, this_b_cells, this_c_cells)))
}
colnames(df) = c("gene", "n_bhve_cells", "n_ctrl_cells")
df$pct_bhve_cells = df$n_bhve_cells / length(bhve_cells)
df$pct_ctrl_cells = df$n_ctrl_cells / length(ctrl_cells)

bb@meta.data$my_velo = bb$initial_size_spliced/bb$initial_size_unspliced
bb@meta.data$cyto = bb@meta.data$my_velo

png("~/scratch/brain/results/rna_velocity_color_my_box.png", width = 800, height = 1000, res = 120)
# print(FeaturePlot(bb, "velocity_self_transition", order = T, pt.size = 0.8) + scale_color_viridis())
# print(FeaturePlot(bb, "my_velo", order = T, pt.size = 0.8) + scale_color_viridis())
print(cytoScoreByIdent(bb, my_idents = c(4, 9, 12, 13, 16, 18, 19, 25, 29, 42, 48, 49)) + ylab("Spliced/Unspliced"))
# print(myFeaturePlot(bb, "cyto",my.col.pal = pal,  cells.use = WhichCells(bb, idents = c(4, 9, 12, 13, 16, 18, 19, 25, 29, 42, 48, 49))))
dev.off()

#==========================================================================================
# Enrichment Testing on GeneX+ vs GeneX- for Many Lists ===================================
#==========================================================================================
geneX_list = read.csv("C:/Users/miles/Downloads/brain/data/markers/bb_all_interesting_genes_final_3count_051021.csv")
geneX_list$gene[which(geneX_list$gene == "slac17a6")] = "slc17a6"
# enrich_list = read.csv("C:/Users/miles/Downloads/brain/data/markers/LG11_highFST_genes_umd2a_120920.txt", he)
# enrich_list = read.table("C:/Users/miles/Downloads/brain/data/markers/LG11_highFST_genes_umd2a_120920.txt", header = F)[,1]
enrich_list = read.table("C:/Users/miles/Downloads/brain/data/markers/pcrc_FST20_30_LG11_evolution_genes_031821.csv", header = T)[,1]
geneX_list_unique = unique(geneX_list$gene)

full_res = data.frame()
for (i in 1:length(geneX_list_unique)) {
# for (i in 1:10) {
  gene = geneX_list_unique[i]
  if (i %% 10 == 0) { print(i) }
  this_res = markerExpPerCellPerClusterQuickGeneX(bb, enrich_list, gene)
  if (! is.na(this_res) ) {
    this_res$gene = gene
    this_res = this_res[,c(ncol(this_res), 1:(ncol(this_res)-1))]
    full_res = rbind(full_res, this_res)
  }
}
rownames(full_res) = paste0(full_res$gene, "_", full_res$cluster)
full_res$bh = p.adjust(full_res$p, method = "BH")
full_res$bon = p.adjust(full_res$p, method = "bonferroni")
full_res$hgnc = gene_info$human[match(full_res$gene, gene_info$mzebra)]
write.csv(full_res, "~/research/brain/results/bb_all_interesting_pcrc_lg11_040520.csv")

# full_res2 = full_res[which(full_res$bon < 0.05,),]
full_res2 = full_res
full_res2$neg_log_bon = -log10(full_res2$bon)
ggplot(full_res2, aes(x = num_cells, y = neg_log_bon, color = mag_pos, size = mag_pos)) + geom_point() + scale_color_manual(values = viridis(4)) + scale_size_manual(values = c(0.5, 1, 2, 3)) + ylab("Negative log10 Bonferroni Adjusted p-value") + xlab("Number of Cells") + geom_text(data=subset(full_res2, mag_pos == "large"), aes(num_cells,neg_log_bon,label=hgnc), color = "black")

#==========================================================================================
# Homolog =================================================================================
#==========================================================================================
gene_info$mzebra_description_clean = as.vector(data.frame(do.call('rbind', strsplit(as.character(gene_info$mzebra_description),'[Source',fixed=TRUE)))[,1])
all_desc = gene_info$mzebra_description_clean[which( ! is.na(gene_info$mzebra_description_clean) & gene_info$biotype == "protein_coding")]
homo_df = data.frame()
homo_df_small = data.frame()
for (desc in unique(all_desc)) {
  this_df = gene_info[which(gene_info$mzebra_description_clean == desc),]
  homo_df = rbind(homo_df, data.frame( rep(desc, nrow(this_df)), this_df$mzebra ))
  homo_df_small = rbind(homo_df_small, data.frame(desc, paste0(this_df$mzebra, collapse = ", ")))
}
colnames(homo_df) = c("description", "gene")
colnames(homo_df_small) = c("description", "genes")
homo_df[which(homo_df$description == "5S ribosomal RNA")]
hist(table(homo_df$description), breaks = 35)

#==========================================================================================
# IEG =====================================================================================
#==========================================================================================
ieg_cons = c("LOC101487312", "egr1", "npas4", "jun", "homer1")
ieg_like = read.csv("C:/Users/miles/Downloads/ieg_like_fos_egr1_npas4_detected_011521.csv", stringsAsFactors = F)[,2]
ieg_like = c(ieg_like, "jun")
# mat = bb@assays$RNA@counts
# mat[which(mat > 1)] = 1
bb$ieg_score = colSums(mat[ieg_cons, ])
# bb$ieg_like_score = colSums(mat[ieg_like, ])
bb$ieg_like_score = colSums(bb@assays$RNA@data[ieg_like,])
temp = rev(brewer.pal(11,"Spectral"))
temp[6] = "gold" # this is what plotCytoTRACE uses
pal = colorRampPalette(temp)
p = FeaturePlot(bb, "ieg_like_score", label = T, order = T, pt.size = 1, split.by = "cond")
p1 = p[[1]] + scale_color_gradientn(limits = c(min(bb$ieg_like_score), max(bb$ieg_like_score)), colors = rev(brewer.pal(11, "Spectral")))
p2 = p[[2]] + scale_color_gradientn(limits = c(min(bb$ieg_like_score), max(bb$ieg_like_score)), colors = rev(brewer.pal(11, "Spectral")))
plot_grid(plotlist=list(p1, p2), ncol = 2)

z_ieg15_res = read.csv("C:/Users/miles/Downloads/bb15_bbmm_demux_combined_diff_ieg_score_cond_hmp_q_by_cluster_102221.csv")
z_ieg15_res$neg_log_q = -log10(z_ieg15_res$q_hmp)
Idents(bb) = bb$seuratclusters15
bb$tmp = z_ieg15_res$neg_log_q[match(bb$seuratclusters15, z_ieg15_res$cluster)]
# Idents(bb) = plyr::revalue(Idents(bb), replace = c())
FeaturePlot(bb, feature = "tmp", label = T) + scale_color_viridis_c()

ggplot(bb@meta.data, aes(x = seuratclusters15, y = ieg_like_score, color = cond)) + geom_boxplot(alpha = 0.6) + geom_point(position = position_jitterdodge(), alpha = 0.01)
ggplot(bb@meta.data, aes(x = good15_names, y = ieg_like_score, color = cond, fill = cond)) + geom_split_violin(alpha = 0.3, scale = "width") + stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.25, position = position_dodge(width = .25))

ieg_like = read.csv(paste0(rna_path, "/data/markers/ieg_like_011521.txt"), stringsAsFactors = F)[,1]
p_df = aggregate(ieg_like_score ~ subsample + seuratclusters15, data = bb@meta.data, mean)
p_df$cond = startsWith(p_df$subsample, "b")
ggplot(p_df, aes(x = seuratclusters15, y = ieg_like_score, color = cond, fill = cond)) + geom_boxplot(alpha = 0.6) + geom_point(position=position_jitterdodge(jitter.width = 0.2), alpha = 0.4)

bb$cyto = read.table("C:/Users/miles/Downloads/brain/results/bb/bb_cyto.txt")[,1]
Idents(bb) = bb$seuratclusters15
bb_0 = subset(bb, idents = 0)
Idents(bb_0) = bb$seuratclusters53
cytoScoreByIdent(bb_0)
cytoScoreByIdent(bb, my_idents = c(4, 9, 12, 13, 16, 18, 19, 25, 29, 48, 49))

pdf("C:/Users/miles/Downloads/brain/results/bb/ieg_like_score_data15.pdf", width = 30, height = 5)
p = FeaturePlot(bb, "ieg_like_score", label = T, order = T, pt.size = 1, split.by = "sample", cols = pal(50))
print(p)
dev.off()

ieg_cells = lapply(ieg_cons, function(x) colnames(bb)[which( bb@assays$RNA@counts[x,] > 0 )])
all_ieg_cells = unique(unlist(ieg_cells))
length(all_ieg_cells[which(all_ieg_cells %in% colnames(bb)[which(bb$cond == "BHVE")])])/length(colnames(bb)[which(bb$cond == "BHVE")])
length(all_ieg_cells[which(all_ieg_cells %in% colnames(bb)[which(bb$cond == "CTRL")])])/length(colnames(bb)[which(bb$cond == "CTRL")])

# Convert to z-score and find behave vs control differences
# bb$ieg_z = scale(as.numeric(as.vector(bb$ieg_like_score)))
bb$ieg_z = as.numeric(as.vector(bb$ieg_like_score))
bhve_samples = c("b1", "b2", "b3", "b4", "b5")
ctrl_samples = c("c1", "c2", "c3", "c4", "c5")
bhve_cells = colnames(bb)[which(bb$sample %in% bhve_samples)]
ctrl_cells = colnames(bb)[which(bb$sample %in% ctrl_samples)]
Idents(bb) = bb$seuratclusters15
ieg_df = data.frame()
ieg_df_stats = data.frame()
for (clust in levels(Idents(bb))) {
  clust_cells = WhichCells(bb, idents = clust)
  clust_b_cells = clust_cells[which(clust_cells %in% bhve_cells)]
  clust_c_cells = clust_cells[which(clust_cells %in% ctrl_cells)]
  clust_b = bb$ieg_like_score[clust_b_cells]
  clust_c = bb$ieg_like_score[clust_c_cells]
  scores = c(clust_b, clust_c)
  ieg_df = rbind(ieg_df, data.frame(scores, clust, c(rep("BHVE", length(clust_b_cells)), rep("CTRL", length(clust_c_cells))) ))
  ieg_df_stats = rbind(ieg_df_stats, data.frame(rep(clust, 2), c(mean(clust_b), mean(clust_c)), c( median(clust_b), median(clust_c)), c("BHVE", "CTRL"), c(length(clust_b_cells), length(clust_c_cells)) ))
}
colnames(ieg_df) = c("Score", "Cluster", "Condition")
colnames(ieg_df_stats) = c("Cluster", "Mean", "Median", "Condition", "Num_Cells")
ieg_df$Score = as.numeric(as.vector(ieg_df$Score))
pal = colorRampPalette(rev(brewer.pal(11,"Spectral")))
ggplot(ieg_df, aes(x=Cluster, y=Score, color = Score, fill = Condition)) + geom_boxplot(alpha=0.6) + geom_jitter(position=position_dodge2(width = 0.6), alpha = 0.01) + scale_fill_manual(values = c("white", "gray47")) + scale_color_gradientn(colors = pal(50))
ggplot(ieg_df_stats, aes(x=Cluster, y=Mean, color = Condition, fill = Condition)) + geom_bar(alpha = 0.8, stat = "identity", position = "dodge")

# Check to see how zack's cluster behave vs control compares to all gene-cluster combos
# gene_df = data.frame()
Idents(bb) = bb$seuratclusters53
mode = "sample" # mode can equal sample or pair
for (i in 1:1) {
  # cluster = 13
  # cluster_cells <- WhichCells(bb, idents = cluster)
  cluster_cells = colnames(bb)
  clust_non_zero = rownames(bb)[which(rowSums(bb@assays$RNA@counts[,cluster_cells]) > 0)]
  clust_b_cells = cluster_cells[which(cluster_cells %in% bhve_cells)]
  cluster_b_mean_exp = rowMeans(expm1(bb@assays$RNA@data[clust_non_zero, clust_b_cells]))
  
  clust_c_cells = cluster_cells[which(cluster_cells %in% ctrl_cells)]
  cluster_c_mean_exp = rowMeans(expm1(bb@assays$RNA@data[clust_non_zero, clust_c_cells]))
  
  avg_logFC = log(cluster_b_mean_exp + 1) - log(cluster_c_mean_exp + 1)
  avg_logFC = unname(avg_logFC)
  
  # Num Cells
  # in_zack = clust_non_zero %in% zack$mzebra[which(zack$cluster == 0)]
  # in_zack = clust_non_zero %in% zack_0
  mat_b = bb@assays$RNA@counts[clust_non_zero, clust_b_cells]
  mat_b[which(mat_b> 1)] = 1
  num_b_cells = rowSums(mat_b)
  mat_c = bb@assays$RNA@counts[clust_non_zero, clust_c_cells]
  mat_c[which(mat_c> 1)] = 1
  num_c_cells = rowSums(mat_c)
  pct_b_cells = (num_b_cells/length(clust_b_cells)) * 100
  pct_c_cells = (num_c_cells/length(clust_c_cells)) * 100
  pct_dif = pct_b_cells - pct_c_cells
  
  if (mode == "pair") {
    for (pair in 1:5) {
      print(samples)
      sample_cells = colnames(bb)[which(bb$sample == paste0("b", pair) | bb$sample == paste0("c", pair))]
      sample_cells_b = clust_b_cells[which(clust_b_cells %in% sample_cells)]
      sample_cells_c = clust_c_cells[which(clust_c_cells %in% sample_cells)]
      sample_b_mean = rowMeans(expm1(bb@assays$RNA@data[clust_non_zero, sample_cells_b]))
      sample_c_mean = rowMeans(expm1(bb@assays$RNA@data[clust_non_zero, sample_cells_c]))
      sample_avg_logFC = log(sample_b_mean + 1) - log(sample_c_mean + 1)
      assign(paste0("avg_logFC_pair", pair), sample_avg_logFC)
    }
    same_dir = sign(avg_logFC_pair1) == sign(avg_logFC_pair2) & 
      sign(avg_logFC_pair1) == sign(avg_logFC_pair3) & 
      sign(avg_logFC_pair1) == sign(avg_logFC_pair4) & 
      sign(avg_logFC_pair1) == sign(avg_logFC_pair5)
    gene_df = data.frame(clust_non_zero, num_b_cells, num_c_cells, pct_b_cells, pct_c_cells, pct_dif, abs(pct_dif), avg_logFC, abs(avg_logFC), avg_logFC_pair1, avg_logFC_pair2, avg_logFC_pair3, avg_logFC_pair4, avg_logFC_pair5, same_dir)
  } else if (mode == "sample") {
    sample_means = list()
    for (sample in unique(bb$sample)) {
      sample_cells = colnames(bb)[which(bb$sample == sample)]
      sample_mean = rowMeans(bb@assays$RNA@data[clust_non_zero, sample_cells])
      sample_means[[sample]] = sample_mean
    }
    gene_df = data.frame(clust_non_zero, num_b_cells, num_c_cells, pct_b_cells, pct_c_cells, pct_dif, abs(pct_dif), avg_logFC, abs(avg_logFC), as.data.frame(sample_means))
  }

  
  # gene_df = data.frame(clust_non_zero, in_zack, num_b_cells, num_c_cells, pct_b_cells, pct_c_cells, pct_dif, abs(pct_dif), avg_logFC, abs(avg_logFC), avg_logFC_pair1, avg_logFC_pair2, avg_logFC_pair3, avg_logFC_pair4, avg_logFC_pair5, same_dir)
  
  
}

# Sample Expression for IEGs
ieg_plot_df = data.frame(sample = unique(bb$sample), ieg_like_mean = colMeans(ieg_gene_df_sample[,c(10:19)]), ieg_mean = colMeans(ieg_gene_df_sample[which(ieg_gene_df_sample$clust_non_zero %in% ieg_cons),c(10:19)]))
ieg_plot_df_2 = melt(ieg_plot_df)
ggplot(ieg_plot_df_2, aes(sample, value, fill = variable, color = variable)) + geom_bar(stat = "identity", position = position_dodge2())


# Investigate whether Pair 4 is different than the others
same_df = data.frame()
for (pair1 in c(1, 2, 3, 4, 5)) {
  for (pair2 in c(1, 2, 3, 4, 5)) {
    num_same = length(which( sign(gene_df[,9+pair1]) == sign(gene_df[,9+pair2]) ))
    same_df = rbind(same_df, t(c( pair1, pair2, num_same )))
  }
}
colnames(same_df) = c("pair1", "pair2", "num_same")
same_df$num_same = as.numeric(as.vector(same_df$num_same))
same_df$pair2 = factor(same_df$pair2, levels = 1:5)
ggplot(same_df, aes(pair1, num_same, color = pair2, fill = pair2)) + geom_bar(stat="identity", position = position_dodge2()) + ylab("Number of Genes With Consistent Direction in Pair 1 and 2")

zack = read.csv("C:/Users/miles/Downloads/bb15_combined_all_pairs_strict_hits_120920_pct_logFC.csv")
hist(abs(avg_logFC), breaks = 50)
zack_0 = zack$abs_avg_logFC[which(zack$cluster == 0)]
hist(zack_0, breaks = 50)
length(abs(avg_logFC)[which(abs(avg_logFC) > mean(zack_0))])

comp = data.frame(avg_logFC=c(abs(avg_logFC), zack_0), 
                  variable=c(rep("All", length(avg_logFC)), rep("Zack", length(zack_0))))
# ggplot(comp, aes(avg_logFC, color = variable, fill = variable)) + geom_histogram(alpha = 0.6) + scale_y_continuous(name = "Freq (All)", sec.axis = sec_axis(~./1000, name="Freq (Zack)"))
ggplot(comp, aes(avg_logFC, color = variable, fill = variable)) + geom_histogram(alpha = 0.6)
gene_df = gene_df[order(in_zack),]
ggplot(gene_df, aes(x = avg_logFC, y = pct_dif, color = in_zack, fill = in_zack, size = in_zack, alpha = in_zack)) + geom_point() + scale_size_manual(values = c(1,2), name = "Genes Selected") + scale_alpha_manual(values=c(0.4, 0.6), name = "Genes Selected") + ggtitle("BHVE vs CTRL for All Genes in Cluster 13") + labs(fill='Genes Selected', color = "Genes Selected") 

gene_df$potential = gene_df$abs.avg_logFC. > 0.02 & gene_df$abs.pct_dif. > 1 & gene_df$same_dir & !gene_df$in_zack
gene_df$potential[which(gene_df$in_zack)] = "true_hit"
gene_df = gene_df[order(gene_df$potential),]
ggplot(gene_df, aes(x = avg_logFC, y = pct_dif, color = potential, fill = potential, size = potential)) + geom_point(alpha = 0.5) + scale_size_manual(values = c(1,2,2)) + ggtitle("Potential Hits in Cluster 0")
ggplot(gene_df, aes(x = avg_logFC, y = pct_dif, color = potential, fill = potential, size = potential)) + geom_point(alpha = 0.5) + scale_size_manual(values = c(1,2,2)) + ggtitle("Potential Hits")

a = rnorm(1000)
hist(a, breaks = 50)
b = scale(a)
hist(b, breaks = 50)

a = colSums(bb@assays$RNA@scale.data[ieg_like,])
b = colSums(scale(bb@assays$RNA@data[ieg_like,]))
hist(a, breaks = 50)
hist(b, breaks = 50)

#==========================================================================================
# IEG Covariance ==========================================================================
#==========================================================================================
library("ggpmisc")
library("hrbrthemes")
ieg_cons = c("LOC101487312", "egr1", "npas4", "jun", "homer1")
ieg_like = read.csv("C:/Users/miles/Downloads/brain/data/markers/ieg_like_011521.txt", stringsAsFactors = F)[,1]
# ieg_like = c(ieg_like, "jun")
bb$ieg_like_score = colSums(bb@assays$RNA@data[ieg_like,])

# Find the mean ieg_like_score in each cluster
cluster15_ieg = data.frame()
cluster_level = "seuratclusters15"
for (cluster in levels(bb@meta.data[,c(cluster_level)])) {
  cluster_cells = colnames(bb)[which(bb@meta.data[,c(cluster_level)] == cluster)]
  for (subsample in sort(unique(bb$subsample))) {
    subsample_cells = colnames(bb)[which(bb$subsample == subsample)]
    this_cells = cluster_cells[which(cluster_cells %in% subsample_cells)]
    num_cells = length(this_cells)
    ieg_like_mean = mean(bb$ieg_like_score[this_cells])
    cluster15_ieg = rbind(cluster15_ieg, t(c(cluster, subsample, num_cells, ieg_like_mean)))
  }
}
colnames(cluster15_ieg) = c("cluster", "sample", "num_cells", "mean_ieg_like_score")
cluster15_ieg$mean_ieg_like_score = as.numeric(as.vector(cluster15_ieg$mean_ieg_like_score))
cluster15_ieg$num_cells = as.numeric(as.vector(cluster15_ieg$num_cells))

# Find the Correlation and R^2 of the cluster combos
cluster15_ieg_combos = data.frame()
for (cluster1 in levels(bb@meta.data[,c(cluster_level)])) {
  for (cluster2 in levels(bb@meta.data[,c(cluster_level)])) {
    if (cluster1 != cluster2) {
      print(paste("Cluster 1,2: ", cluster1, cluster2))
      # this_df = cluster15_ieg[which(cluster15_ieg$cluster %in% c(cluster1, cluster2)),]
      this_df = cbind( cluster15_ieg[which(cluster15_ieg$cluster == cluster1),], cluster15_ieg[which(cluster15_ieg$cluster == cluster2),] )
      colnames(this_df) = c("cluster1", "sample", "num_cells_1", "mean_ieg_like_score_1", "cluster2", "sample2", "num_cells_2", "mean_ieg_like_score_2")
      this_df = this_df[which( this_df$num_cells_1 > 0 & this_df$num_cells_2 > 0 ),]
      this_df$sample = factor(this_df$sample, levels = unique(sort(bb$subsample)))
      this_r2 = cor(this_df$mean_ieg_like_score_1, this_df$mean_ieg_like_score_2) ^ 2
      this_r2_str = format(round(this_r2, 4), nsmall = 4)
      this_df$isBehave = startsWith(as.character(as.vector(this_df$sample)), "b")
      
      # All Plot (ie r^2 for behave and control together)
      # png(paste0(rna_path, "/results/bb/ieg_covariance_15/", cluster1, "_", cluster2, ".png"), width = 600, height = 500, res = 90)
      # print(ggplot(this_df, aes(mean_ieg_like_score_1, mean_ieg_like_score_2, color = sample)) + geom_point() + geom_smooth(method=lm , color="#569688", fill="#69b3a2", se=TRUE) + scale_colour_discrete(drop=TRUE, limits = levels(this_df$sample)) + labs(title= paste("Cluster", cluster1, "vs Cluster", cluster2, "- R^2 =", this_r2_str)))
      # dev.off()
      
      # Behave
      this_df_b = this_df[which( this_df$isBehave ),]
      this_r_b = cor(this_df_b$mean_ieg_like_score_1, this_df_b$mean_ieg_like_score_2)
      this_r2_b = this_r_b ^ 2
      this_r2_b_str = format(round(this_r2_b, 4), nsmall = 4)
      num_b = nrow(this_df_b)
      
      # Control
      this_df_c = this_df[which( ! this_df$isBehave ),]
      this_r_c = cor(this_df_c$mean_ieg_like_score_1, this_df_c$mean_ieg_like_score_2)
      this_r2_c = this_r_c ^ 2
      this_r2_c_str = format(round(this_r2_c, 4), nsmall = 4)
      num_c = nrow(this_df_c)
      
      # Behave and Control as Separate lines
      # png(paste0(rna_path, "/results/bb/ieg_covariance_15_b_c/", cluster1, "_", cluster2, ".png"), width = 600, height = 500, res = 90)
      # print(ggplot(this_df, aes(mean_ieg_like_score_1, mean_ieg_like_score_2, color = isBehave)) + geom_point() + geom_smooth(method=lm, se=FALSE) + scale_colour_discrete(drop=TRUE, limits = levels(this_df$sample)) + labs(title= paste("Cluster", cluster1, "vs Cluster", cluster2), subtitle = paste("Behave R^2 =", this_r2_b_str, "and Control R^2 =", this_r2_c_str)))
      # dev.off()
      
      cluster15_ieg_combos = rbind(cluster15_ieg_combos, t(c(cluster1, cluster2, this_r2, this_r2_b, this_r2_c, this_r_b, this_r_c, num_b, num_c)))
    }
  }
}
colnames(cluster15_ieg_combos) = c("cluster1", "cluster2", "this_r2", "r2_behave", "r2_control", "r_behave", "r_control", "num_b", "num_c")
cluster15_ieg_combos$cluster1 = factor(cluster15_ieg_combos$cluster1, levels = levels(bb@meta.data[,c(cluster_level)]))
cluster15_ieg_combos$cluster2 = factor(cluster15_ieg_combos$cluster2, levels = levels(bb@meta.data[,c(cluster_level)]))
cluster15_ieg_combos$this_r2 = as.numeric(as.vector(cluster15_ieg_combos$this_r2))
cluster15_ieg_combos$r2_behave = as.numeric(as.vector(cluster15_ieg_combos$r2_behave))
cluster15_ieg_combos$r2_control = as.numeric(as.vector(cluster15_ieg_combos$r2_control))
cluster15_ieg_combos$r_behave = as.numeric(as.vector(cluster15_ieg_combos$r_behave))
cluster15_ieg_combos$r_control = as.numeric(as.vector(cluster15_ieg_combos$r_control))
cluster15_ieg_combos$num_b = as.numeric(as.vector(cluster15_ieg_combos$num_b))
cluster15_ieg_combos$num_c = as.numeric(as.vector(cluster15_ieg_combos$num_c))

# R^2 Plots
p_r2_a = ggplot(cluster15_ieg_combos, aes(cluster1, cluster2, fill = this_r2)) + geom_tile() + scale_fill_viridis(discrete=F, limits=c(0, max(cluster15_ieg_combos$this_r2))) + theme_classic() + ggtitle("All Samples IEG Covariance") + labs(fill = "R^2")
p_r2_b = ggplot(cluster15_ieg_combos, aes(cluster1, cluster2, fill = r2_behave)) + geom_tile() + scale_fill_viridis(discrete=F, limits=c(0, max(cluster15_ieg_combos$r2_behave))) + theme_classic() + ggtitle("Behave IEG Covariance") + labs(fill = "R^2")
p_r2_c = ggplot(cluster15_ieg_combos, aes(cluster1, cluster2, fill = r2_control)) + geom_tile() + scale_fill_viridis(discrete=F, limits=c(0, max(cluster15_ieg_combos$r2_control))) + theme_classic() + ggtitle("Control IEG Covariance") + labs(fill = "R^2")

png("C:/Users/miles/Downloads/brain/results/bb/ieg_covariance_all.png", width = 800, height = 750, res = 90)
print(p_r2_a)
dev.off()
png("C:/Users/miles/Downloads/brain/results/bb/ieg_covariance_behave.png", width = 800, height = 750, res = 90)
print(p_r2_b)
dev.off()
png("C:/Users/miles/Downloads/brain/results/bb/ieg_covariance_control.png", width = 800, height = 750, res = 90)
print(p_r2_c)
dev.off()

# R Plots
cluster15_ieg_combos_ind = sapply(1:nrow(cluster15_ieg_combos), function(x) as.numeric(as.vector(cluster15_ieg_combos$cluster1[x])) > as.numeric(as.vector(cluster15_ieg_combos$cluster2[x])))
cluster15_ieg_combos2 = cluster15_ieg_combos[which(cluster15_ieg_combos_ind),]
p_r_b = ggplot(cluster15_ieg_combos2, aes(cluster1, cluster2, fill = r_behave)) + geom_tile() + scale_fill_distiller(palette = "Spectral") + theme_classic() + ggtitle("Behave IEG Covariance") + labs(fill = "R")
p_r_c = ggplot(cluster15_ieg_combos2, aes(cluster1, cluster2, fill = r_control)) + geom_tile() + scale_fill_distiller(palette = "Spectral") + theme_classic() + ggtitle("Control IEG Covariance") + labs(fill = "R")

png("C:/Users/miles/Downloads/brain/results/bb/ieg_covariance_behave_r.png", width = 800, height = 750, res = 90)
print(p_r_b)
dev.off()
png("C:/Users/miles/Downloads/brain/results/bb/ieg_covariance_control_r.png", width = 800, height = 750, res = 90)
print(p_r_c)
dev.off()

# Find Significant Correlations
my_cor_t = function(r, n) (r * sqrt(n - 2))/sqrt(1 - r**2)
my_cor_p = function(t, n) 2*pt(-abs(t), df=n-2)

# Convert Behave correlations to p values
r_behave_mat = acast(cluster15_ieg_combos, cluster1 ~ cluster2, value.var="r_behave")
n_behave_mat = acast(cluster15_ieg_combos, cluster1 ~ cluster2, value.var="num_b")
t_behave_mat = my_cor_t(r_behave_mat, n_behave_mat)
p_behave_mat = my_cor_p(t_behave_mat, n_behave_mat)

# Convert Control correlations to p values
r_control_mat = acast(cluster15_ieg_combos, cluster1 ~ cluster2, value.var="r_control")
n_control_mat = acast(cluster15_ieg_combos, cluster1 ~ cluster2, value.var="num_c")
t_control_mat = my_cor_t(r_control_mat, n_control_mat)
p_control_mat = my_cor_p(t_control_mat, n_control_mat)

# Find Number of Significant Behave Correlations
p_behave_mat[upper.tri(p_behave_mat)] = NA
p_behave_df = na.omit(melt(p_behave_mat))
colnames(p_behave_df) = c("Cluster1", "Cluster2", "p")
p_behave_df$bh = p.adjust(p_behave_df$p, method = "BH")
length(which(p_behave_df$p < 0.05))
length(which(p_behave_df$bh < 0.05))

# Find Number of Significant Control Correlations
p_control_mat[upper.tri(p_control_mat)] = NA
p_control_df = na.omit(melt(p_control_mat))
colnames(p_control_df) = c("Cluster1", "Cluster2", "p")
p_control_df$bh = p.adjust(p_control_df$p, method = "BH")
length(which(p_control_df$p < 0.05))
length(which(p_control_df$bh < 0.05))

# Plot a Single: Behave and Control ScatterPlot with Separate Trendlines
singleBVCScatter = function(cluster1, cluster2) {
  this_df = cbind( cluster15_ieg[which(cluster15_ieg$cluster == cluster1),], cluster15_ieg[which(cluster15_ieg$cluster == cluster2),] )
  colnames(this_df) = c("cluster1", "sample", "num_cells_1", "mean_ieg_like_score_1", "cluster2", "sample2", "num_cells_2", "mean_ieg_like_score_2")
  this_df = this_df[which( this_df$num_cells_1 > 0 & this_df$num_cells_2 > 0 ),]
  this_df$sample = factor(this_df$sample, levels = unique(sort(bb$subsample)))
  this_r2 = cor(this_df$mean_ieg_like_score_1, this_df$mean_ieg_like_score_2) ^ 2
  this_r2_str = format(round(this_r2, 4), nsmall = 4)
  this_df$isBehave = startsWith(as.character(as.vector(this_df$sample)), "b")
  # Behave
  this_df_b = this_df[which( this_df$isBehave ),]
  this_r2_b_str = format(round(cor(this_df_b$mean_ieg_like_score_1, this_df_b$mean_ieg_like_score_2) ^ 2, 4), nsmall = 4)
  # Control
  this_df_c = this_df[which( ! this_df$isBehave ),]
  this_r2_c_str = format(round(cor(this_df_c$mean_ieg_like_score_1, this_df_c$mean_ieg_like_score_2)^2, 4), nsmall = 4)
  # Plot
  png(paste0(rna_path, "/results/bb/ieg_covariance_b_c/", cluster1, "_", cluster2, ".png"), width = 600, height = 500, res = 90)
  print(ggplot(this_df, aes(mean_ieg_like_score_1, mean_ieg_like_score_2, color = isBehave)) + geom_point() + geom_smooth(method=lm, se=FALSE) + labs(title= paste("Cluster", cluster1, "vs Cluster", cluster2), subtitle = paste("Behave R^2 =", this_r2_b_str, "and Control R^2 =", this_r2_c_str)))
  dev.off()
}

# singleBVCScatter(15, 2)
# singleBVCScatter(23, 0)

# Behave vs Control
p_bvc_mat = r_to_p(r_behave_mat, r_control_mat, n_behave_mat, n_control_mat)
p_bvc_mat[upper.tri(p_bvc_mat)] = NA
p_bvc_df = na.omit(melt(p_bvc_mat))
colnames(p_bvc_df) = c("Cluster1", "Cluster2", "p")
p_bvc_df$bh = p.adjust(p_bvc_df$p, method = "BH")
length(which(p_bvc_df$p < 0.05))
length(which(p_bvc_df$bh < 0.05))
# singleBVCScatter(18, 13)
singleBVCScatter(1, 2)
singleBVCScatter(4, 9)

colnames(p_bvc_df)[c(3,4)] = c("p_bvc", "bh_bvc")
p_bvc_df[, c('p_b', 'bh_b')] = p_behave_df[,c('p', 'bh')]
p_bvc_df[, c('p_c', 'bh_c')] = p_control_df[,c('p', 'bh')]
p_bvc_df[, c('r_b', 'r_c', 'num_b', 'num_c')] = cluster15_ieg_combos[match(paste0(p_bvc_df$Cluster1, p_bvc_df$Cluster2), paste0(cluster15_ieg_combos$cluster1, cluster15_ieg_combos$cluster2)), c('r_behave', 'r_control', 'num_b', 'num_c')]

p_bvc_df$Cluster1 = factor(p_bvc_df$Cluster1, levels = levels(bb@meta.data[,c(cluster_level)]))
p_bvc_df$Cluster2 = factor(p_bvc_df$Cluster2, levels = levels(bb@meta.data[,c(cluster_level)]))
png("C:/Users/miles/Downloads/brain/results/bb/ieg_covariance_bvc_bh_15.png", width = 800, height = 750, res = 90)
print(ggplot(p_bvc_df, aes(Cluster1, Cluster2, fill = bh)) + geom_tile() + scale_fill_viridis(begin = 1, end = 0) + theme_classic() + ggtitle("Behave vs Control Fisher's Z Transformation") + labs(fill = "BH"))
dev.off()

r_dif_mat = r_behave_mat - r_control_mat
r_abs_dif_mat = abs(r_behave_mat - r_control_mat)
r_dif_mat[upper.tri(r_dif_mat)] = NA
r_abs_dif_mat[upper.tri(r_abs_dif_mat)] = NA
r_dif_df = na.omit( melt(r_dif_mat) )
r_abs_dif_df = na.omit( melt(r_abs_dif_mat) )
r_dif_df$Var1 = factor(r_dif_df$Var1, levels = 0:14)
r_dif_df$Var2 = factor(r_dif_df$Var2, levels = 0:14)
r_abs_dif_df$Var1 = factor(r_abs_dif_df$Var1, levels = 0:14)
r_abs_dif_df$Var2 = factor(r_abs_dif_df$Var2, levels = 0:14)
p_r_bvc = ggplot(r_dif_df, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(11,"Spectral")))(50), limits=c(-max(abs(r_dif_df$value), na.rm = T), max(abs(r_dif_df$value), na.rm = T))) + theme_classic() + ggtitle("Behave vs Control Correlation Difference") + labs(fill = "R Dif") + xlab("Cluster 1") + ylab("Cluster 2")
p_r_abs_bvc = ggplot(r_abs_dif_df, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_viridis(discrete = F) + theme_classic() + ggtitle("Absolute Behave vs Control Correlation Difference") + labs(fill = "Abs R Dif") + xlab("Cluster 1") + ylab("Cluster 2")
png("C:/Users/miles/Downloads/brain/results/bb/ieg_covariance_bvc_r.png", width = 800, height = 750, res = 90)
print(p_r_bvc)
dev.off()
png("C:/Users/miles/Downloads/brain/results/bb/ieg_covariance_abs_bvc_r.png", width = 800, height = 750, res = 90)
print(p_r_abs_bvc)
dev.off()

#******************************************************************************************
# Demux Matching===========================================================================
#******************************************************************************************
demux = read.table("~/scratch/brain/ffm/JTS07-C1/outs/demux/c1_demux.best", header = T)
df_b1 = df[which(df$sample == "c1"),]
df_b1$demux = demux$SNG.1ST[match(df_b1$sub_sub, demux$BARCODE)]
all_combos = expand.grid(unique(df_b1$demux), unique(df_b1$sub))
all_combos = paste0(all_combos[,1], "-", all_combos[,2])
all_combos = sort(all_combos)
real_paste = paste0(df_b1$demux, "-", df_b1$sub)
real_paste = factor(real_paste, levels = all_combos)
df_b1_p = data.frame(table(real_paste))
df_b1_p[,c("demux", "scsplit")] = colsplit(df_b1_p$real_paste, pattern = "-", names = c("demux", "scsplit"))
png("~/scratch/brain/results/c1_scsplit_to_demux.png", width = 700, height = 500, res = 90)
ggplot(df_b1_p, aes(x = scsplit, y = Freq, fill = demux, color = demux)) + geom_bar(stat = "identity", position = "dodge")
dev.off()

df = bb@meta.data
df$sub = colsplit(rownames(df), pattern = "_", names = c("1", "2"))[,2]
df$sub_sub = colsplit(df$sub, pattern = "_", names = c("1", "2"))[,1]
df$demux = "idk"
for (sample in c("b1", "b2", "b3", "b4", "b5", "c1", "c2", "c3", "c4", "c5")) {
  sample_caps = toupper(sample)
  demux = read.table(paste0("~/scratch/brain/ffm/JTS07-", sample_caps, "/outs/demux/", sample, "_demux.best"), header = T)
  df$demux[which(df$sample == sample)] = demux$SNG.1ST[match(df$sub_sub[which(df$sample == sample)], demux$BARCODE)]
}

subsample_demux = c( "1B18" = "b1.1", "1B4"  = "b1.2", "1B5"  = "b1.3", "1B11" = "b1.4",
                     "2B10" = "b2.1", "2B13" = "b2.2", "2B17" = "b2.3", "2B19" = "b2.4",
                     "3B7"  = "b3.1", "3B6"  = "b3.2", "3B2"  = "b3.3", "3B9"  = "b3.4",
                     "4B12" = "b4.1", "4B14" = "b4.2", "4B25" = "b4.3",
                     "5B26" = "b5.1", "5B24" = "b5.2", "5B22" = "b5.3", "5B23" = "b5.4",
                     "1C5"  = "c1.1", "1C11" = "c1.2", "1C4"  = "c1.3", "1C17" = "c1.4",
                     "2C19" = "c2.1", "2C10" = "c2.2", "2C14" = "c2.3", "2C18" = "c2.4",
                     "3C7"  = "c3.1", "3C9"  = "c3.2", "3C6"  = "c3.3", "3C2"  = "c3.4",
                     "4C25" = "c4.1", "4C12" = "c4.2", "4C13" = "c4.3",
                     "5C26" = "c5.1", "5C22" = "c5.2", "5C23" = "c5.3", "5C24" = "c5.4" )
df$subsample = plyr::revalue(df$demux, replace = subsample_demux)
bb$demux = df$demux
bb$subsample = df$subsample

demux_res
demux_sing2 = read.delim("~/Downloads/b2_demux_merge.sing2")
demux_res2$S1 = demux_res2$S2 = demux_res2$S3 = demux_res2$S4 = ""
for (barcode in demux_res2$BARCODE) {
  this_rows = demux_sing2[which(demux_sing2$BARCODE == barcode),]
  demux_res2[which(demux_res2$BARCODE == barcode), c("S1", "S2", "S3", "S4")] = this_rows$LLK1
}
demux_res2$BS1 = demux_res2$BS2 = demux_res2$BS3 = demux_res2$BS4 = ""
for (barcode in demux_res2$BARCODE) {
  demux_res2[which(demux_res2$BARCODE == barcode), c("BS1", "BS2", "BS3", "BS4")] = names(sort(demux_res2[which(demux_res2$BARCODE == barcode), c("S1", "S2", "S3", "S4")]))
}
demux_res2$BS1_BS2 = sapply(1:nrow(demux_res2), function(x) as.numeric(demux_res2[x,demux_res2$BS1[x]]) / as.numeric(demux_res2[x,demux_res2$BS2[x]]))
demux_res2$BS1_BS3 = sapply(1:nrow(demux_res2), function(x) as.numeric(demux_res2[x,demux_res2$BS1[x]]) / as.numeric(demux_res2[x,demux_res2$BS3[x]]))
demux_res2$BS1_BS4 = sapply(1:nrow(demux_res2), function(x) as.numeric(demux_res2[x,demux_res2$BS1[x]]) / as.numeric(demux_res2[x,demux_res2$BS4[x]]))

ggplot(demux_res2, aes(x = is1B4, y = BS1_BS4, color = is1B4, fill = is1B4)) + geom_boxplot(alpha = 0.6)

demux_res3 = melt(demux_res2, id.vars = colnames(demux_res2)[which(! colnames(demux_res2) %in% c("BS1_BS2", "BS1_BS3", "BS1_BS4"))])
demux_res3$variable = paste(demux_res3$variable, demux_res3$is1B4)
demux_res3$variable = factor(demux_res3$variable, levels = c("BS1_BS2 FALSE", "BS1_BS3 FALSE", "BS1_BS4 FALSE", "BS1_BS2 TRUE", "BS1_BS3 TRUE", "BS1_BS4 TRUE"))
demux_res3$variable = plyr::revalue(demux_res3$variable, replace = c("BS1_BS2 FALSE" = "BS1/BS2 Real", "BS1_BS3 FALSE" = "BS1/BS3 Real", "BS1_BS4 FALSE" = "BS1/BS4 Real", "BS1_BS2 TRUE" = "BS1/BS2 B1", "BS1_BS3 TRUE" = "BS1/BS3 B1", "BS1_BS4 TRUE" = "BS1/BS4 B1"))
my_pal = c("#57CC99", "#38A3A5", "#22577A", "#FCD2D1", "#FE8F8F", "#FF5C58")
ggplot(demux_res3, aes(x = value, color = variable, fill = variable)) + geom_density(alpha = 0.6) + scale_color_manual(values = my_pal) + scale_fill_manual(values = my_pal)
library(ggridges)
ggplot(demux_res3, aes(x = value, y = variable, color = variable, fill = variable)) + geom_density_ridges(alpha = 0.9) + scale_color_manual(values = my_pal) + scale_fill_manual(values = my_pal) + theme_ridges() + theme(legend.position = "none")

#******************************************************************************************
# Receptor Diff ===========================================================================
#******************************************************************************************
int_df = read.csv("~/research/brain/data/markers/bb_all_interesting_genes_final_3count_051021.csv")
int_df$gene[which(int_df$gene == "slac17a6")] = "slc17a6"
intBHVECTRL = function(gene) {
  gene_pos_cells = colnames(bb)[which(bb@assays$RNA@counts[gene,] > 0)]
  bhve_cells = colnames(bb)[which(bb$cond == "BHVE")]
  ctrl_cells = colnames(bb)[which(bb$cond == "CTRL")]
  gene_pos_bhve_cells = gene_pos_cells[which(gene_pos_cells %in% bhve_cells)]
  gene_pos_ctrl_cells = gene_pos_cells[which(gene_pos_cells %in% ctrl_cells)]
  res = pct_dif_avg_logFC(bb, cells.1 = gene_pos_bhve_cells, cells.2 = gene_pos_ctrl_cells)
  res$pct_dif = res$pct.1 - res$pct.2
  res$abs_pct_dif = abs(res$pct_dif)
  res$abs_avg_logFC = abs(res$avg_logFC)
  res = res[order(res$abs_avg_logFC, decreasing = T),]
  logFC_sum = sum(res$abs_avg_logFC[1:100])
  pct_sum = sum(res$abs_pct_dif[1:100])
  return(list(logFC_sum, pct_sum))
}



library(parallel)
test = mclapply(int_df$gene, intBHVECTRL, mc.cores = 24)
test = c()
for (i in 1:nrow(int_df)) {
  print(i)
  this_sum = intBHVECTRL(int_df$gene[i])
  test = c(test, this_sum)
}

zfp36_deg = read.csv("~/research/brain/results/zfp36_bvc.csv")
FeaturePlot(bb, "LOC101486937", order = T)

prokr2_deg = read.csv("~/research/brain/results/prokr2_bvc.csv")
FeaturePlot(bb, "LOC101465959", order = T)

#******************************************************************************************
# BHVE DEG Module==========================================================================
#******************************************************************************************
library(WGCNA)
bb15_deg = read.csv("~/research/brain/results/bb15_glmmseq_cond_gsi_control_and_cond_spawn_control_sig_genes_q_by_cluster_100521.csv")
deg_clusters = sort(unique(bb15_deg$cluster))

mat = bb@assays$RNA@counts
mat[which(mat > 1)] = 1

for (this_clust in deg_clusters) {
  pos_clust_deg = bb15_deg$X[which(bb15_deg$cluster == this_clust & bb15_deg$condCTRL < 0 & bb15_deg$condCTRL.1 < 0)]
  neg_clust_deg = bb15_deg$X[which(bb15_deg$cluster == this_clust & bb15_deg$condCTRL > 0 & bb15_deg$condCTRL.1 > 0)]
  
  if (length(pos_clust_deg) > 0) {
    pos_score = colSums(mat[pos_clust_deg, which(bb$seuratclusters15 == this_clust & bb$cond == "BHVE")])
    
    exp = bb@assays$RNA@data[pos_clust_deg, which(bb$seuratclusters15 == this_clust & bb$cond == "BHVE")]
    pheatmap::pheatmap(exp, show_rownames = F, scale = "row", show_colnames = F, filename = "test.pdf", width = 10, height = 10)
    
    r_mat = cor(x = as.matrix(t(exp)))
    diag(r_mat) = 0
    # r_mat[which(upper.tri(r_mat, diag = T))] = NA
    pheatmap::pheatmap(r_mat, show_rownames = F, show_colnames = F, filename = "test2.pdf", width = 10, height = 10)
    
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    sft = pickSoftThreshold(r_mat, powerVector = powers, verbose = 5)
    sft = min(sft$fitIndices$Power[which(sft$fitIndices$SFT.R.sq >= 0.9)])
    # sft = 6
    # r_mat = abs(r_mat) ** sft
    adj_mat = adjacency.fromSimilarity(r_mat)
    TOM = TOMsimilarity(adj_mat)
    dissTOM = 1-TOM
    geneTree = hclust(as.dist(dissTOM), method = "average")
    plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
    minModuleSize = 10;
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize);
    table(dynamicMods)
    df = data.frame(module = dynamicMods, row.names = rownames(r_mat))
    my_cols = brewer.pal(max(df$module)+1, "Set2")
    names(my_cols) = sort(unique(df$module))
    my_cols = list(module = my_cols)
    # hcl = "none"
    # my_callback = function(hcl, mat) { print(hcl); hcl <<- hcl; return(hcl) }
    test =pheatmap::pheatmap(r_mat, show_rownames = F, show_colnames = F, cutree_rows = 4, annotation_row = df, filename = "test4.pdf", width = 10, height = 10, annotation_colors = my_cols)
    hcl = test$tree_col
    dend = as.dendrogram(hcl)
    clust1 = dend[[1]] %>% labels
    df = data.frame(module = rep("High", length(clust1)), row.names = clust1)
    df = rbind(df, data.frame(module = rep("none", length(which( !rownames(r_mat) %in% clust1 ))), row.names = rownames(r_mat)[which( !rownames(r_mat) %in% clust1 )]))
    test = pheatmap::pheatmap(r_mat, show_rownames = F, show_colnames = F, cutree_rows = 2, cutree_cols = 2, annotation_row = df, annotation_col = df, filename = "test5.pdf", width = 10, height = 10)
    
    pos_score2 = colSums(mat[clust1, which(bb$seuratclusters15 == this_clust & bb$cond == "CTRL")])
    hist(pos_score2, breaks = 50, main = "CTRL Cells, High Module Score, Cluster 0")
    pos_score2 = colSums(mat[clust1, which(bb$seuratclusters15 == this_clust & bb$cond == "BHVE")])
    hist(pos_score2, breaks = 50, main = "BHVE Cells, High Module Score, Cluster 0")
    exp = bb@assays$RNA@data[clust1, which(bb$seuratclusters15 == this_clust & bb$cond == "BHVE")]
    pheatmap::pheatmap(exp, show_rownames = F, show_colnames = F, filename = "test6.pdf", width = 10, height = 10)
    myFeaturePlot(bb, "tmp", my.col.pal = pal, cells.use = colnames(bb)[which(bb$seuratclusters15 == 0 & bb$cond == "BHVE")])
    df = bb@meta.data[which(bb$seuratclusters15 == 0 & bb$cond == "BHVE"),]
    ggplot(df, aes(x = seuratclusters53, y = tmp, color = tmp)) + geom_violin() + geom_boxplot(width = 0.2) + geom_point(position = position_jitter(), alpha = 0.1) + scale_color_gradientn(colors = pal(50)) + ylab("High Module Score") + ggtitle("High Module in BHVE Cells in Cluster 0")
  }
}

# baDEG
badeg = read.csv("~/research/brain/results/deg_depth_build_badeg_glmmseq_demux_all_clusters_all_tests_pair_subjectinpair_pool_subjectinpool_sig_all_genes_100821_q_hgnc.csv")
badeg_tg = read.csv("~/Downloads/badeg_most_significantly_enriched_pathways_results_100821.csv")
badeg_gene_converter = read.csv("~/Downloads/all_clusters_cond_deg_list_alt_names.csv")

badeg$cluster_hgnc = paste0(badeg$cluster, "_", badeg$hgnc)
badeg_gene_converter$cluster[which(badeg_gene_converter$cluster == "8_Glut")] = 0
# badeg_gene_converter$cluster_num = convert15$old[match(badeg_gene_converter$cluster, convert15$new.full)]
badeg_gene_converter$cluster_hgnc = paste0(badeg_gene_converter$cluster, "_", badeg_gene_converter$hgnc)
badeg$hgnc2 = badeg_gene_converter$final[match(badeg$cluster_hgnc, badeg_gene_converter$cluster_hgnc)]

tg1_hgnc = unlist(strsplit(badeg_tg$Hit.in.Query.List[1], ","))
badeg_tg1 = badeg$mzebra[which( badeg$cluster == badeg_tg$cluster[1] & badeg$hgnc2 %in% tg1_hgnc )]
exp = bb@assays$RNA@data[badeg_tg1, which(bb$seuratclusters15 == badeg_tg$cluster[1] & bb$cond == "BHVE")]
pheatmap::pheatmap(exp, show_rownames = F, show_colnames = F, filename = "test.pdf", width = 10, height = 10)

bb$tmp = colSums(bb@assays$RNA@data[badeg_tg1,])
df = bb@meta.data[which(bb$seuratclusters15 == 0 & bb$cond == "BHVE"),]
ggplot(df, aes(x = seuratclusters53, y = tmp, color = tmp)) + geom_violin() + geom_boxplot(width = 0.2) + geom_point(position = position_jitter(), alpha = 0.1) + scale_color_gradientn(colors = pal(50)) + ylab("High Module Score") + ggtitle("Neurexin Module in BHVE Cells in Cluster 0")
df = bb@meta.data[which(bb$seuratclusters15 == 0),]
ggplot(df, aes(x = seuratclusters53, y = tmp, color = cond)) + geom_violin() + geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), alpha = 0.6, color = "black", outlier.size = 0, aes(fill = cond)) + geom_point(position = position_jitterdodge(), alpha = 0.1) + ylab("Neurexin Module Score") + ggtitle("Neurexin Module in Cluster 0")

#******************************************************************************************
# Neuron Time==============================================================================
#******************************************************************************************
library("SeuratWrappers")
library("monocle3")
ieg_like25 = read.csv("~/research/brain//data/ieg_like_fos_egr1_npas4_detected_011521.csv")[,1]
ieg_like61 = read.csv("~/research/brain/results/IEG61_list_for_GG_081021.csv")[,1]
ieg_like450 = read.csv("~/research/brain/results/IEG450_list_for_GG_081221.csv")[,1]
# ieg_like = read.csv("~/research/brain//data/ieg_like_fos_egr1_npas4_detected_011521.csv")[,1]
# ieg_like = read.csv("~/research/brain/results/IEG61_list_for_GG_081021.csv")[,1]
# ieg_like = read.csv("~/research/brain/results/IEG450_list_for_GG_081221.csv")[,1]
bb_neuron = subset(bb, cells = colnames(bb)[which(! bb$seuratclusters15 %in% c(4, 9,13) )])
ieg_like_mat = bb_neuron@assays$RNA@counts[ieg_like25,]
ieg_like_mat[which( ieg_like_mat > 1 )] = 1
# bb.cds = as.cell_data_set(bb_neuron)
Seurat_Object_Diet <- DietSeurat(bb_neuron, graphs = "umap")
bb.cds = as.cell_data_set(Seurat_Object_Diet)
bb.cds = preprocess_cds(bb.cds, num_dim = 100, use_genes = ieg_like61) # "You're computing too large a percentage of total singular values, use a standard svd instead."
# bb.cds = align_cds(bb.cds, alignment_group = "batch")
bb.cds <- reduce_dimension(bb.cds)
bb.cds <- cluster_cells(bb.cds)
bb.cds <- learn_graph(bb.cds, use_partition = F)
bb.cds <- order_cells(bb.cds, root_cells=colnames(bb_neuron)[which(colSums(ieg_like_mat) == 0)])
# # Programtically pick root
# get_earliest_principal_node <- function(cds, time_bin="130-170"){
#   cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
# 
#   closest_vertex <-
#     cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
#   closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
#   root_pr_nodes <-
#     igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
#                                                               (which.max(table(closest_vertex[cell_ids,]))))]
# 
#   root_pr_nodes
# }
# bb.cds <- order_cells(bb.cds, root_pr_nodes=get_earliest_principal_node(bb.cds))

png("~/scratch/brain/results/bb_neuron_monocle450.png", width = 1000, height = 900, res = 100)
plot_cells(bb.cds, color_cells_by = "pseudotime", label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=0)
dev.off()
write.csv(bb.cds@principal_graph_aux@listData$UMAP$pseudotime, "~/scratch/brain/results/bb_neuron_monocle450.csv")

# bb <- Seurat::AddMetaData(object = bb, metadata = bb.cds@principal_graph_aux@listData$UMAP$pseudotime, col.name = "n_time")
# bb$n_time[which(! is.finite(bb$n_time) )] = NA
# myFeaturePlot(bb, feature = "n_time", na.blank = T, my.col.pal = colorRampPalette(inferno(100)))
n_time = read.csv("~/research/brain/results/bb_neuron_monocle450.csv")[,2]
# bb_neuron = Seurat::AddMetaData(object = bb_neuron, metadata = bb.cds@principal_graph_aux@listData$UMAP$pseudotime, col.name = "n_time")
bb_neuron = Seurat::AddMetaData(object = bb_neuron, metadata = n_time, col.name = "n_time")
bb_neuron$n_time[which(! is.finite(bb_neuron$n_time) )] = NA
myFeaturePlot(bb_neuron, feature = "n_time", na.blank = T, my.col.pal = colorRampPalette(inferno(100)))

# write.csv(bb_neuron$n_time, "~/research/brain/results/bb_neuron_ntime_ieg61_root025.csv")

# Create a dataframe with a column for cell, the pseudotime, and the data for each IEG
time_df = data.frame(cell = colnames(bb_neuron), time = bb_neuron$n_time)
for (ieg in ieg_like25) {
  time_df[,ieg] = bb_neuron@assays$RNA@data[ieg,]
}
write.csv(time_df, "~/research/brain/results/bb_neuron_ntime_df_ieg61_root025.csv")
write.csv(time_df, "~/research/brain/results/bb_neuron_ntime_df_ieg450_root025.csv")

# time_df = read.csv("~/research/brain/results/bb_neuron_ntime_61_1.csv")
time_df$sum = rowSums(time_df[,3:ncol(time_df)])

ggplot(time_df, aes(time, npas4, color = time)) + geom_point() + geom_smooth(formula = y ~ x, method='loess') + theme_bw()

# Do Loess Smoothing
time_df_p = data.frame()
for (ieg in ieg_like25) {
  print(which(ieg_like25 == ieg))
  loessMod75 <- loess(time_df[,ieg] ~ time_df[,"time"], span=0.75)
  smoothed75 <- predict(loessMod75)
  this_df = time_df[,c("cell", "time")]
  this_df$ieg = ieg
  this_df$smooth = smoothed75
  time_df_p = rbind(time_df_p, this_df)
}
time_df_p$label = gene_info$human[match(time_df_p$ieg, gene_info$mzebra)]
time_df_p$label[which(time_df_p$ieg == "egr1")] = "EGR1A"
time_df_p$label[which(time_df_p$ieg == "LOC101475150")] = "EGR1B"
time_df_p$label[which(time_df_p$ieg == "LOC101487312")] = "FOSA"
time_df_p$label[which(time_df_p$ieg == "LOC101487601")] = "FOSB"
time_df_p$label[which(time_df_p$ieg == "LOC101486758")] = "RTN4RL2A"
time_df_p$label[which(time_df_p$ieg == "LOC101465305")] = "RTN4RL2B"
time_df_p = time_df_p[order(time_df_p$time),]
ggplot(time_df_p, aes(time, smooth, color = label)) + geom_line() + geom_label_repel(data=time_df_p[tail(which(time_df_p$label == "JUND"), 1),], aes(label = label), nudge_x = 1, na.rm = TRUE)

# Find time of peak
ieg_peaks = unlist(lapply(ieg_like25, function(x) { 
  this_time_df_p = time_df_p[which(time_df_p$ieg == x),]
  this_time_df_p$time[which.max(this_time_df_p$smooth)] 
}))
ieg_peaks = data.frame(ieg = ieg_like25, hgnc = gene_info$human[match(ieg_like25, gene_info$mzebra)], time_of_peak = ieg_peaks)
write.csv(ieg_peaks, "~/research/brain/results/bb_neuron_time_peaks_ieg450_root025")

# mono_umap = as.data.frame(reducedDim(bb.cds, "UMAP"))
# mono_umap$cluster = bb_neuron$seuratclusters15
# ggplot(mono_umap, aes(x = V1, y = V2, color = cluster)) + geom_point(size = 0.5) + theme_classic()
# mono_umap$ieg_sum = colSums(ieg_like_mat)
# mono_umap$ieg_sum_nonzero = mono_umap$ieg_sum > 0
# ggplot(mono_umap, aes(x = V1, y = V2, color = ieg_sum)) + geom_point() + theme_classic()
# ggplot(mono_umap, aes(x = V1, y = V2, color = ieg_sum_nonzero)) + geom_point() + theme_classic()

#******************************************************************************************
# BHVE Time================================================================================
#******************************************************************************************
subsample_time = c("b1.1" = 1,  "b1.2" = 2,  "b1.3" = 3,  "b1.4" = 4,
                   "b2.1" = 5,  "b2.2" = 6,  "b2.3" = 7,  "b2.4" = 8,
                   "b3.1" = 9,  "b3.2" = 10, "b3.3" = 11, "b3.4" = 12,
                   "b4.1" = 13, "b4.2" = 14, "b4.3" = 15,
                   "b5.1" = 16, "b5.2" = 17, "b5.3" = 18, "b5.4" = 19,
                   "c1.1" = 20, "c1.2" = 21, "c1.3" = 22, "c1.4" = 23,
                   "c2.1" = 24, "c2.2" = 25, "c2.3" = 26, "c2.4" = 27,
                   "c3.1" = 28, "c3.2" = 29, "c3.3" = 30, "c3.4" = 31,
                   "c4.1" = 32, "c4.2" = 33, "c4.3" = 34,
                   "c5.1" = 35, "c5.2" = 36, "c5.3" = 37, "c5.4" = 38)
subsample_time = c("b1.1" = 1,  "b1.2" = 2,  "b1.3" = 3,  "b1.4" = 4,
                   "b2.1" = 5,  "b2.2" = 6,  "b2.3" = 7,  "b2.4" = 8,
                   "b3.1" = 9,  "b3.2" = 10, "b3.3" = 11, "b3.4" = 12,
                   "b4.1" = 13, "b4.2" = 14, "b4.3" = 15,
                   "b5.1" = 16, "b5.2" = 17, "b5.3" = 18, "b5.4" = 00,
                   "c1.1" = 00, "c1.2" = 00, "c1.3" = 00, "c1.4" = 00,
                   "c2.1" = 00, "c2.2" = 00, "c2.3" = 00, "c2.4" = 00,
                   "c3.1" = 00, "c3.2" = 00, "c3.3" = 00, "c3.4" = 00,
                   "c4.1" = 00, "c4.2" = 00, "c4.3" = 00,
                   "c5.1" = 00, "c5.2" = 00, "c5.3" = 00, "c5.4" = 00)
bb$time = as.numeric(as.vector(plyr::revalue(bb$subsample, subsample_time)))

red_pal = c("#431313", "#f75656")
yellow_pal = c("#434313","#f7f456")
green_pal = c("#13431b", "#56f76e")
gene_time_df = data.frame(data = bb@assays$RNA@data["fosb",], time = bb$time)
# loess_mod <- loess(data ~ time, gene_time_df)
# pred <- predict(loess_mod, gene_time_df, se=TRUE)
# gene_time_df$lwl <- pred$fit-1.96*pred$se.fit
# gene_time_df$upl <- pred$fit+1.96*pred$se.fit
ggplot(gene_time_df, aes(time, data, color = time)) + geom_point() + geom_smooth(formula = y ~ x, method='loess', color = green_pal[1]) + theme_bw() + scale_color_gradientn(colors = green_pal) + ggtitle("fosb")
ggplot(gene_time_df, aes(time, data, color = time)) + geom_point() + geom_smooth(formula = y ~ x, method='loess') + theme_bw() + ggtitle("EGR1")
ggplot(gene_time_df, aes(time, data)) + geom_boxplot() + geom_jitter(position=position_dodge2(width = 0.6), alpha = 0.01)

#
#
#
all_r = read.table("~/Downloads/all_no_cr_relate.relatedness2", header = T)
all_r$isFirst = all_r$RELATEDNESS_PHI > 0.225 & all_r$RELATEDNESS_PHI < 0.48
ggplot(all_r, aes(INDV1, INDV2, fill = RELATEDNESS_PHI)) + geom_raster() + scale_fill_viridis_c() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
ggplot(all_r, aes(INDV1, INDV2, fill = isFirst)) + geom_raster() + scale_fill_viridis_d() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

#*******************************************************************************
# Add Hulls Around Cluster =====================================================
#*******************************************************************************
library("ggforce")
library("concaveman")
df = data.frame(cluster15 = bb$seuratclusters15, cluster53 = bb$seuratclusters53, col15 = convert15$col[match(bb$seuratclusters15, convert15$old)], col53 = convert53$col[match(bb$seuratclusters53, convert53$old)])
df[,c("UMAP_1", "UMAP_2")] = as.data.frame(bb@reductions$umap@cell.embeddings)
ggplot(df, aes(UMAP_1, UMAP_2, color = col53)) + geom_point(size = 0.6) + theme_void() + geom_mark_hull(aes(group = col15, fill = col15, color = col15), concavity = 20, expand = unit(2.5, "mm")) + scale_color_identity() + scale_fill_identity()

#******************************************************************************************
# Sierra ==================================================================================
#******************************************************************************************
library("Sierra")
FindPeaks("~/scratch/brain/bs/JTS07/pbs/b4_peaks.txt", gtf.file = "~/scratch/m_zebra_ref/GCF_000238955.4_M_zebra_UMD2a_genomic.gff", bamfile = "~/scratch/brain/bs/JTS07/pbs/B1_intron/outs/possorted_genome_bam.bam", junctions.file = "~/scratch/brain/bs/JTS07/pbs/b4_intron_junctions.bed", ncores = 24)

my.peak.data.set.table = data.frame(Peak_file = c(sapply(1:5, function(x) paste0("~/scratch/brain/bs/JTS07/pbs/b", x, "_peaks.txt")), sapply(1:5, function(x) paste0("~/scratch/brain/bs/JTS07/pbs/c", x, "_peaks.txt"))), 
                                    Identifier = c(sapply(1:5, function(x) paste0("b", x)), sapply(1:5, function(x) paste0("c", x))))
MergePeakCoordinates(peak.dataset.table = my.peak.data.set.table, output.file = "~/scratch/brain/bs/JTS07/pbs/peak_coord_tmp.txt", ncores = 24)
CountPeaks(peak.sites.file = "~/scratch/brain/bs/JTS07/pbs/peak_coord_tmp.txt", gtf.file = "~/scratch/m_zebra_ref/GCF_000238955.4_M_zebra_UMD2a_genomic.gff", whitelist.file = "~/scratch/brain/bs/JTS07/pbs/B1_intron/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", bamfile = "~/scratch/brain/bs/JTS07/pbs/B1_intron/outs/possorted_genome_bam.bam", output.dir = "~/scratch/brain/bs/JTS07/pbs/peak_counts/")
AnnotatePeaksFromGTF(peak.sites.file = "~/scratch/brain/bs/JTS07/pbs/peak_coord.txt", 
                     gtf.file = "~/scratch/m_zebra_ref/GCF_000238955.4_M_zebra_UMD2a_genomic.gff",
                     output.file = "~/scratch/brain/bs/JTS07/pbs/peak_annot.txt", 
                     genome = NULL)

my.count.dirs = list("b1" = "~/scratch/brain/bs/JTS07/pbs/b1_peak_counts/", "b2" = "~/scratch/brain/bs/JTS07/pbs/b2_peak_counts/", 
                     "b3" = "~/scratch/brain/bs/JTS07/pbs/b3_peak_counts/", "b4" = "~/scratch/brain/bs/JTS07/pbs/b4_peak_counts/", "b5" = "~/scratch/brain/bs/JTS07/pbs/b5_peak_counts/", 
                     "c1" = "~/scratch/brain/bs/JTS07/pbs/c1_peak_counts/", "c2" = "~/scratch/brain/bs/JTS07/pbs/c2_peak_counts/", 
                     "c3" = "~/scratch/brain/bs/JTS07/pbs/c3_peak_counts/", "c4" = "~/scratch/brain/bs/JTS07/pbs/c4_peak_counts/", "c5" = "~/scratch/brain/bs/JTS07/pbs/c5_peak_counts/")
for (pair in 4:5) {
  print(paste("Pair:", pair))
  print("Aggregating Counts")
  AggregatePeakCounts(peak.sites.file = "~/scratch/brain/bs/JTS07/pbs/peak_coord.txt",
                      count.dirs = unlist(my.count.dirs[c(paste0("b", pair), paste0("c", pair))]),
                      exp.labels = c(paste0("b", pair), paste0("c", pair)),
                      output.dir = paste0("~/scratch/brain/bs/JTS07/pbs/counts_pair", pair,"/"))
  print("Creating Seurat Object")
  peak.counts = ReadPeakCounts(data.dir = paste0("~/scratch/brain/bs/JTS07/pbs/counts_pair", pair,"/"))
  peak.counts.meta.data = data.frame(barcode = reshape2::colsplit(colnames(peak.counts), "-", c(1,2))[,1], sample = reshape2::colsplit(colnames(peak.counts), "-", c(1,2))[,2])
  bb_barcodes = data.frame( sample = bb$sample, real = colnames(bb), barcode = reshape2::colsplit(reshape2::colsplit(reshape2::colsplit(colnames(bb), "_", c(1,2))[,2], "_", c(1,2))[,1], "-", c(1,2))[,1] )
  for (sample in unique(bb$sample)) {
    if (sample %in% unique(peak.counts.meta.data$sample)) {
      sample_bb_barcodes = bb_barcodes[which(bb_barcodes$sample == sample),]
      colnames(peak.counts)[which(peak.counts.meta.data$sample == sample)] = sample_bb_barcodes$real[match(peak.counts.meta.data$barcode[which(peak.counts.meta.data$sample == sample)], sample_bb_barcodes$barcode)] 
    }
  }
  peaks.seurat2 = NewPeakSeurat(peak.data = peak.counts, annot.info = peak.annotations, cell.idents = bb$seuratclusters15[which(bb$sample %in% c(paste0("b", pair), paste0("c", pair)))])
  saveRDS(peaks.seurat2, paste0("~/scratch/brain/data/pair", pair, "_sierra.rds"))
}

AggregatePeakCounts(peak.sites.file = "~/scratch/brain/bs/JTS07/pbs/peak_coord.txt",
                    count.dirs = unlist(my.count.dirs),
                    exp.labels = names(my.count.dirs),
                    output.dir = paste0("~/scratch/brain/bs/JTS07/pbs/all_counts/"))
peak.counts = ReadPeakCounts(data.dir = paste0("~/scratch/brain/bs/JTS07/pbs/all_counts/"))
peak.counts.meta.data = data.frame(barcode = reshape2::colsplit(colnames(peak.counts), "-", c(1,2))[,1], sample = reshape2::colsplit(colnames(peak.counts), "-", c(1,2))[,2])
for (sample in unique(bb$sample)) {
  if (sample %in% unique(peak.counts.meta.data$sample)) {
    sample_bb_barcodes = bb_barcodes[which(bb_barcodes$sample == sample),]
    colnames(peak.counts)[which(peak.counts.meta.data$sample == sample)] = sample_bb_barcodes$real[match(peak.counts.meta.data$barcode[which(peak.counts.meta.data$sample == sample)], sample_bb_barcodes$barcode)] 
  }
}
peaks.seurat2 = NewPeakSeurat(peak.data = peak.counts, annot.info = peak.annotations, cell.idents = bb$seuratclusters15)

# Load BB
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("SeuratObject")
bb = readRDS(paste0(rna_path, "data/bb_cc_04072021.RDS"))
Idents(bb) = bb$seurat_clusters
library(pacman)
p_unload(SeuratDisk)
p_unload(Seurat)
p_load(Seurat)

# peak.counts.b1 <- ReadPeakCounts(data.dir = "~/scratch/brain/bs/JTS07/pbs/b1_peak_counts/")
# peak.counts.c1 <- ReadPeakCounts(data.dir = "~/scratch/brain/bs/JTS07/pbs/c1_peak_counts/")
peak.counts = ReadPeakCounts(data.dir = "~/scratch/brain/bs/JTS07/pbs/all_counts_tmp/")
peak.annotations <- read.table("~/scratch/brain/bs/JTS07/pbs/peak_annot.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
peak.annotations$gene_id = sub("(.*).*:.*:.*-.*:.*", "\\1", rownames(peak.annotations))
# rownames(peak.annotations) = str_replace_all(rownames(peak.annotations), "_", "-")

peak.counts.meta.data = data.frame(barcode = reshape2::colsplit(colnames(peak.counts), "-", c(1,2))[,1], sample = reshape2::colsplit(colnames(peak.counts), "-", c(1,2))[,2])
bb_barcodes = data.frame( sample = bb$sample, real = colnames(bb), barcode = reshape2::colsplit(reshape2::colsplit(reshape2::colsplit(colnames(bb), "_", c(1,2))[,2], "_", c(1,2))[,1], "-", c(1,2))[,1] )
for (sample in unique(bb$sample)) {
  if (sample %in% unique(peak.counts.meta.data$sample)) {
    sample_bb_barcodes = bb_barcodes[which(bb_barcodes$sample == sample),]
    colnames(peak.counts)[which(peak.counts.meta.data$sample == sample)] = sample_bb_barcodes$real[match(peak.counts.meta.data$barcode[which(peak.counts.meta.data$sample == sample)], sample_bb_barcodes$barcode)] 
  }
}
# colnames(peak.counts) = bb_barcodes$real[match(colnames(peak.counts), bb_barcodes$barcode)]

peaks.seurat2 = NewPeakSeurat(peak.data = peak.counts, 
                              annot.info = peak.annotations, 
                              cell.idents = bb$seuratclusters15[which(bb$sample %in% c("b1", "c1"))])
saveRDS(peaks.seurat2, "~/scratch/brain/data/pair1_sierra.rds")

# peaks.seurat <- PeakSeuratFromTransfer(peak.data = peak.counts, genes.seurat = bb, annot.info = peak.annotations)
# saveRDS(peaks.seurat, "~/scratch/brain/data/pair1_sierra.rds")

# Local
# bbs = readRDS("~/scratch/brain/data/pair4_sierra.rds")
bbs = readRDS("~/research/brain/data/pair1_sierra.rds")
colnames_to_transfer = colnames(bb@meta.data)[-which(colnames(bb@meta.data) %in% c("orig.ident", "nCount_RNA", "nFeature_RNA", "geneLvlID"))]
bbs@meta.data[, colnames_to_transfer] = bb@meta.data[match(colnames(bbs), colnames(bb)), colnames_to_transfer]
bbs@reductions$umap = Seurat::CreateDimReducObject(embeddings = bb@reductions$umap@cell.embeddings[match(colnames(bbs), colnames(bb)),], key = "UMAP_", assay = "RNA")
rownames(bbs@reductions$umap@cell.embeddings) = str_replace_all(rownames(bbs@reductions$umap@cell.embeddings), "\\.", "-")
bbs$cond.clust15 = paste0(bbs$cond, ".", bbs$seuratclusters15)
bbs$cond.clust53 = paste0(bbs$cond, ".", bbs$seuratclusters53)

Idents(bbs) = bbs$cond
bulk.res = DUTest(bbs, population.1 = "BHVE", population.2 = "CTRL", exp.thresh = 0.01, fc.thresh = 0)
bulk.res$peak = rownames(bulk.res)

Idents(bbs) = bbs$cond.clust15
res15 = data.frame()
for (clust15 in 11:14) {
  print(clust15)
  this.res15 = DUTest(bbs, population.1 = paste0("BHVE.", clust15), population.2 = paste0("CTRL.", clust15), exp.thresh = 0.01, fc.thresh = 0)
  if (nrow(this.res15) > 0) {
    this.res15$peak = rownames(this.res15)
    this.res15$cluster = clust15
    res15 = rbind(res15, this.res15)
  }
}

bulk.res = list()
bulk.res[[1]] = read.csv("~/research/brain/results/pair1_dtu_bulk.csv")
bulk.res[[2]] = read.csv("~/research/brain/results/pair2_dtu_bulk.csv")
bulk.res[[3]] = read.csv("~/research/brain/results/pair3_dtu_bulk.csv")
bulk.res[[4]] = read.csv("~/research/brain/results/pair4_dtu_bulk.csv")
bulk.res[[5]] = read.csv("~/research/brain/results/pair5_dtu_bulk.csv")
common_res = bulk.res[[1]]$peak[which(bulk.res[[1]]$peak %in% bulk.res[[2]]$peak
                                      & bulk.res[[1]]$peak %in% bulk.res[[3]]$peak
                                      & bulk.res[[1]]$peak %in% bulk.res[[4]]$peak
                                      & bulk.res[[1]]$peak %in% bulk.res[[5]]$peak)]
table(table(c( bulk.res[[1]]$peak, bulk.res[[2]]$peak, bulk.res[[3]]$peak, bulk.res[[4]]$peak, bulk.res[[5]]$peak)))

bulk.res.full = data.frame(i = 1:length(common_res))
for (i in c(1, 5)) {
  print(i)
  colnames(bulk.res[[i]]) = paste0(colnames(bulk.res[[i]]), ".", i)
  bulk.res[[i]] = bulk.res[[i]][which(bulk.res[[i]]$peak %in% common_res),]
  bulk.res.full = cbind(bulk.res.full, bulk.res[[i]])
}
bulk.res.full = bulk.res.full
bulk.res.full$consistent = unlist(sapply(1:nrow(bulk.res.full), function(x) sign(bulk.res.full$Log2_fold_change.1[x]) == sign(bulk.res.full$Log2_fold_change.5[x]) ))
# bulk.res.full$consistent = unlist(sapply(1:nrow(bulk.res.full), function(x) sign(bulk.res.full$Log2_fold_change.1[x]) == sign(bulk.res.full$Log2_fold_change.2[x]) & sign(bulk.res.full$Log2_fold_change.1[x]) == sign(bulk.res.full$Log2_fold_change.2[x]) & sign(bulk.res.full$Log2_fold_change.3[x]) == sign(bulk.res.full$Log2_fold_change.4[x]) & sign(bulk.res.full$Log2_fold_change.1[x]) == sign(bulk.res.full$Log2_fold_change.5[x])))

clust15.res = list()
clust15.res[[1]] = read.csv("~/research/brain/results/pair1_dtu_15.csv")
clust15.res[[2]] = read.csv("~/research/brain/results/pair2_dtu_15.csv")
clust15.res[[3]] = read.csv("~/research/brain/results/pair3_dtu_15.csv")
clust15.res[[4]] = read.csv("~/research/brain/results/pair4_dtu_15.csv")
clust15.res[[5]] = read.csv("~/research/brain/results/pair5_dtu_15.csv")
for (i in 0:14) {
  print(i)
  common_res = clust15.res[[1]]$peak[which(clust15.res[[1]]$peak[which(clust15.res[[1]]$cluster == i)] %in% clust15.res[[2]]$peak[which(clust15.res[[2]]$cluster == i)]
                                        & clust15.res[[1]]$peak[which(clust15.res[[1]]$cluster == i)] %in% clust15.res[[3]]$peak[which(clust15.res[[3]]$cluster == i)]
                                        & clust15.res[[1]]$peak[which(clust15.res[[1]]$cluster == i)] %in% clust15.res[[4]]$peak[which(clust15.res[[4]]$cluster == i)]
                                        & clust15.res[[1]]$peak[which(clust15.res[[1]]$cluster == i)] %in% clust15.res[[5]]$peak[which(clust15.res[[5]]$cluster == i)])]
  # print(length(common_res))
  print(table(table(c( clust15.res[[1]]$peak[which(clust15.res[[1]]$cluster == i)], clust15.res[[2]]$peak[which(clust15.res[[2]]$cluster == i)], clust15.res[[3]]$peak[which(clust15.res[[3]]$cluster == i)], clust15.res[[4]]$peak[which(clust15.res[[4]]$cluster == i)], clust15.res[[5]]$peak[which(clust15.res[[5]]$cluster == i)]))))
}

clust53.res = list()
clust53.res[[1]] = read.csv("~/research/brain/results/pair1_dtu_53.csv")
clust53.res[[2]] = read.csv("~/research/brain/results/pair2_dtu_53.csv")
clust53.res[[3]] = read.csv("~/research/brain/results/pair3_dtu_53.csv")
clust53.res[[4]] = read.csv("~/research/brain/results/pair4_dtu_53.csv")
clust53.res[[5]] = read.csv("~/research/brain/results/pair5_dtu_53.csv")
for (i in 0:52) {
  print(i)
  common_res = clust53.res[[1]]$peak[which(clust53.res[[1]]$peak[which(clust53.res[[1]]$cluster == i)] %in% clust53.res[[2]]$peak[which(clust53.res[[2]]$cluster == i)]
                                           & clust53.res[[1]]$peak[which(clust53.res[[1]]$cluster == i)] %in% clust53.res[[3]]$peak[which(clust53.res[[3]]$cluster == i)]
                                           & clust53.res[[1]]$peak[which(clust53.res[[1]]$cluster == i)] %in% clust53.res[[4]]$peak[which(clust53.res[[4]]$cluster == i)]
                                           & clust53.res[[1]]$peak[which(clust53.res[[1]]$cluster == i)] %in% clust53.res[[5]]$peak[which(clust53.res[[5]]$cluster == i)])]
  # print(length(common_res))
  print(table(table(c( clust53.res[[1]]$peak[which(clust53.res[[1]]$cluster == i)], clust53.res[[2]]$peak[which(clust53.res[[2]]$cluster == i)], clust53.res[[3]]$peak[which(clust53.res[[3]]$cluster == i)], clust53.res[[4]]$peak[which(clust53.res[[4]]$cluster == i)], clust53.res[[5]]$peak[which(clust53.res[[5]]$cluster == i)]))))
}

bbs = readRDS("~/research/brain/data/all_sierra.rds")
peak_info = Tool(bbs, "Sierra")
peak_info$Peak_number_num = reshape2::colsplit(peak_info$Peak_number, "Peak", c(1,2))[,2]
Peak_number_num_max = aggregate(Peak_number_num ~ Gene_name, peak_info, max)
peak_info$Peak_number_num_max = Peak_number_num_max[match(peak_info$Gene_name, Peak_number_num_max[,1]),2]
all15 = read.csv("~/research/brain/results/all_dtu_15.csv")
all53 = read.csv("~/research/brain/results/all_dtu_53.csv")
all15 = all15[order(all15$padj, decreasing = F),]
all53 = all53[order(all53$padj, decreasing = F),]
all15$pct_dif = all15$population1_pct - all15$population2_pct
all15 = all15[which(all15$Log2_fold_change > 0 & all15$pct_dif > 0),]
all53$pct_dif = all53$population1_pct - all53$population2_pct
all53 = all53[which(all53$Log2_fold_change > 0 & all53$pct_dif > 0),]
FeaturePlot(bbs, all15$peak[2], label = T, order = T)
FeaturePlot(bb, all15$gene_name[2], label = T, order = T)

# Find DTU Hits where Peaks have different expression profiles
deg15 = read.csv("~/research/brain/results/bb_all_markers_15clusters_102820_more_info.csv")
deg53 = read.csv("~/research/brain/results/bb_all_markers_53clusters_102720_more_info.csv")
int_dtu15 = data.frame()
int_dtu53 = data.frame()
for (clust15 in 0:14) {
  this_dtu = all15[which(all15$cluster == clust15),]
  this_deg = deg15[which(deg15$cluster == clust15),]
  int_dtu15 = rbind(int_dtu15, this_dtu[which(! this_dtu$gene_name %in% this_deg$gene & this_dtu$Log2_fold_change > 0 & this_dtu$population1_pct > this_dtu$population2_pct ),])
}
for (clust53 in 0:52) {
  this_dtu = all53[which(all53$cluster == clust53),]
  this_deg = deg53[which(deg53$cluster == clust53),]
  int_dtu53 = rbind(int_dtu53, this_dtu[which(! this_dtu$gene_name %in% this_deg$gene & this_dtu$Log2_fold_change > 0 & this_dtu$population1_pct > this_dtu$population2_pct ),])
}
int_dtu15$hgnc = gene_info$human[match(int_dtu15$gene_name, gene_info$mzebra)]
int_dtu53$hgnc = gene_info$human[match(int_dtu53$gene_name, gene_info$mzebra)]

int_dtu15_hgnc = unique(int_dtu15$hgnc[which(! is.na(int_dtu15$hgnc) )])
clipboard(int_dtu15_hgnc)
clipboard(unique(int_dtu53$hgnc[which(! is.na(int_dtu53$hgnc) )]))
FeaturePlot(bbs, int_dtu15$peak[1], label = T, order = T)
FeaturePlot(bb, int_dtu15$gene_name[1], label = T, order = T)

total_gene_counts_test = data.frame(peak = rownames(peak_info), gene = peak_info$Gene_name, count = unname(bbs@assays$RNA@counts[,1]))
total_gene_counts_test = aggregate(count ~ gene, total_gene_counts_test, sum)
bb@assays$RNA@counts[total_gene_counts_test$gene[1:5],1]

# The DTUTest from Sierra doesn't do exactly what I'm looking for.
# A peak could be a DTU even if it matches closesly with the gene expression profile.
# For example, a peak could be a DTU for cluster 3, but if the parent gene is
# a DEG for cluster 3, then that peak isn't so interesting anymore.
# I'm interested in cases where the expression profiles of the peaks are distinct
# from one another. Like where one peak is expressed only in cluster 3 and another
# peak from the same parent gene is expressed only in cluster 7. That would be
# an extreme example, but that's what I'm looking for.
# Helper Functions
peaksSimilarRange = function(gene) {
  #' For a gene, find peaks that are less than standard deviations away from the
  #' mean number of counts per peak.
  idx_remove = 1 # required to enter while loop
  idx = which(min_peak_info$Gene_name == gene)
  while(length(idx_remove) > 0) {
    idx_remove = c()
    mean_value = mean(min_peak_info$count[idx])
    values = abs(min_peak_info$count[idx] - mean_value) / mean_value
    idx_remove = which(values >= 0.75 | values <= -0.75)
    # values = scale(min_peak_info$count[idx])
    # idx_remove = which(values >= 1 | values <= -1)
    if (length(idx_remove) > 0) {
      idx = idx[-idx_remove]
    }
  }
  rownames(min_peak_info)[idx]
}
myDTU = function(gene, cluster) {
  # names = clust15_peak_gene_ratio$peak[which(clust15_peak_gene_ratio$gene == gene)]
  values = scale(clust15_peak_gene_ratio[which(clust15_peak_gene_ratio$gene == gene),cluster])[,1]
  values
}
myDTU2 = function(peak) {
  # names = clust15_peak_gene_ratio$peak[which(clust15_peak_gene_ratio$gene == gene)]
  values = scale(t(clust15_peak_gene_ratio[which(clust15_peak_gene_ratio$peak == peak), as.character(0:14)]))[,1]
  values
}

# Ensure that peaks express a minimum number of counts
min_count = 50
min_count_sample = 5
min_peaks = rownames(bbs)[which(rowSums(bbs@assays$RNA@counts) >= min_count)]

sample_above_min = list()
for (this_sample in unique(bbs$sample)) {
  sample_above_min[[this_sample]] = rownames(bbs)[which(rowSums(bbs@assays$RNA@counts[,which(bb$sample == this_sample)]) >= min_count_sample)]
}
sample_above_min_df = as.data.frame(table(unlist(sample_above_min)))
min_peaks = min_peaks[which(min_peaks %in% sample_above_min_df$Var1[which(sample_above_min_df$Freq == 10)])]
min_peak_info = peak_info[which(rownames(peak_info) %in% min_peaks),]
min_peak_info$count = rowSums(bbs@assays$RNA@counts[min_peaks,])

# Ensure that peaks have similar levels of expression
library(parallel)
num.cores = detectCores()
similarPeaks = unlist(mclapply(unique(min_peak_info$Gene_name), peaksSimilarRange, mc.cores = num.cores/3))
min_peak_info = min_peak_info[similarPeaks,]
min_Peak_number_num_max = aggregate(Peak_number ~ Gene_name, min_peak_info, length)
min_peak_info$Num_parts = min_Peak_number_num_max[match(min_peak_info$Gene_name, min_Peak_number_num_max[,1]),2]

# Calculate the percent of cells expressing each peak per cluster
clust15_cells = list()
clust15_peak_pct = data.frame(peak = rownames(min_peak_info), gene = min_peak_info$Gene_name)
for (clust15 in 0:14) {
  print(clust15)
  clust_cells = colnames(bbs)[which(bbs$seuratclusters15 == clust15)]
  
  mat1 = bbs@assays$RNA@counts[rownames(min_peak_info), clust_cells]
  mat1[which(mat1 > 1)] = 1
  num_pos_cells1 = rowSums(mat1)
  pct_pos_cells1 = (num_pos_cells1/length(clust_cells)) * 100
  
  clust15_peak_pct = cbind(clust15_peak_pct, pct_pos_cells1)
}
colnames(clust15_peak_pct) = c("peak", "gene", 0:14)

# Divide that number by the total percent for all peaks belonging to the same parent gene
clust15_gene_pct = as.data.frame(lapply(0:14, function(x) { 
  this_df = clust15_peak_pct[, c("gene", as.character(x))]
  colnames(this_df) = c("gene", "cluster")
  this_res = aggregate(cluster ~ gene, this_df, sum)[,2]
}))
colnames(clust15_gene_pct) = 0:14
rownames(clust15_gene_pct) = sort(unique(clust15_peak_pct$gene))
clust15_peak_gene_ratio = as.data.frame(lapply(as.character(0:14), function(x) clust15_peak_pct[,x] / clust15_gene_pct[match(clust15_peak_pct$gene, rownames(clust15_gene_pct)),x]))
colnames(clust15_peak_gene_ratio) = 0:14
rownames(clust15_peak_gene_ratio) = clust15_peak_pct$peak
clust15_peak_gene_ratio$peak = clust15_peak_pct$peak
clust15_peak_gene_ratio$gene = clust15_peak_pct$gene

# Scale values by cluster to find how many standard deviations a peak is from the other peaks
parent_genes = unique(clust15_peak_gene_ratio$gene)
clust15_scale = data.frame(peak = unlist(lapply(parent_genes, function(gene) clust15_peak_gene_ratio$peak[which(clust15_peak_gene_ratio$gene == gene)])))
clust15_scale$gene = min_peak_info$Gene_name[match(clust15_scale$peak, rownames(min_peak_info))]
for (clust15 in as.character(0:14)) {
  print(clust15)
  scaled_values = unlist(mcmapply(myDTU, parent_genes, rep(clust15, length(parent_genes)), mc.cores = num.cores))
  clust15_scale = cbind(clust15_scale, scaled_values)
}
colnames(clust15_scale) = c("peak", "gene", 0:14)
clust15_scale$num_parts = min_peak_info$Num_parts[match(clust15_scale$peak, rownames(min_peak_info))]
clust15_scale$gene_part = min_peak_info$Gene_part[match(clust15_scale$peak, rownames(min_peak_info))]

# Experimenting
test = unlist(mclapply(clust15_peak_pct$peak, myDTU2, mc.cores = num.cores/3))
test2 = as.data.frame(matrix(test, ncol = 15))
colnames(test2) = 0:14
test2$peak = clust15_peak_pct$peak
test2$gene = clust15_peak_pct$gene
test2$sum = unlist(mclapply(test2$peak, scaleDif, mc.cores = num.cores/3))
test2$num_parts = min_peak_info$Num_parts[match(test2$peak, rownames(min_peak_info))]
test2$gene_part = min_peak_info$Gene_part[match(test2$peak, rownames(min_peak_info))]
scaleDif = function(peak) {
  gene = test2$gene[which(test2$peak == peak)]
  dif = sum(abs(test2[which(test2$peak == peak), as.character(0:14)] - colMeans(test2[which(test2$gene == gene), as.character(0:14)], na.rm = T)), na.rm = T)
  dif
}

neuronal_clusters = as.character(0:14)[which(! as.character(0:14) %in% c(4, 9,13) )]
clust15_scale2 = clust15_scale[,as.character(0:14)]
clust15_scale2[clust15_scale2 < 0] = 0
clust15_scale_abs = abs(clust15_scale2[,as.character(0:14) ])
clust15_scale_abs$sum = rowSums(clust15_scale_abs, na.rm = T)
clust15_scale_abs$peak = clust15_scale$peak
clust15_scale_abs$gene = clust15_scale$gene
clust15_scale_abs$num_parts = clust15_scale$num_parts
clust15_scale_abs$gene_part = clust15_scale$gene_part
clust15_scale_abs$count = min_peak_info$count[match(clust15_scale_abs$peak, rownames(min_peak_info))]
clust15_scale_abs = clust15_scale_abs[order(clust15_scale_abs$sum, decreasing = T)[1:10000],]
ggplot(clust15_scale_abs, aes(x = count, y = sum, color = sum)) + geom_point() + scale_color_viridis_c()

# plot
clust15_scale[which(clust15_scale[,"0"] > 3 & clust15_peak_pct[match(clust15_scale$peak, clust15_peak_pct$peak),"0"] > 3)[1:5], c("peak", "gene", "0", "num_parts", "gene_part")]
for (i in as.character(0:14)) {
  df0 = clust15_scale[which(clust15_scale[,i] > 0.5), c("peak", "gene", i, "num_parts", "gene_part")]
  colnames(df0)[which(colnames(df0) == i)] = "cluster_z_score"
  df0$pct = clust15_peak_pct[match(df0$peak, clust15_peak_pct$peak),i]
  pdf(paste0("~/research/brain/results/sierra_cluster/", i, ".pdf"), width = 5, height = 5)
  print(ggplot(df0, aes(x = cluster_z_score, y = pct)) + geom_point(alpha = 0.4) + geom_text_repel(data=df0[which(df0$cluster_z_score >= 2.5 & df0$pct >= 5),], aes(label = gene_part), color = "red") + ggtitle(paste0("Cluster ", i, ": Scaled(Peak/Total) vs Percent of cells [Only Keep Peaks in Similar Range Per Gene]")))
  dev.off()
}

clust15_scale[which(clust15_scale[,"14"] > 8 & clust15_peak_pct[match(clust15_scale$peak, clust15_peak_pct$peak),"14"] > 3)[1:5], c("peak", "gene", "14", "num_parts", "gene_part")]

# DotPlot
gene = "LOC101477058"
DotPlot(bbs, features = rownames(min_peak_info)[which(min_peak_info$Gene_name == gene)])  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + scale_x_discrete(labels = min_peak_info$Gene_part[which(min_peak_info$Gene_name == gene)])
FeaturePlot(bbs, rownames(min_peak_info)[which(min_peak_info$Gene_name == gene)], order = T, label = T, pt.size = 1)

temp = rev(brewer.pal(11,"Spectral"))
temp[6] = "gold" # this is what plotCytoTRACE uses
pal = colorRampPalette(temp)
# DimPlot(bbs)
# PlotRelativeExpressionUMAP(bbs, row.names(this.res15)[c(1,3)])
FeaturePlot(bbs, row.names(this.res15)[c(1,3)])
myFeaturePlot(bbs, row.names(this.res15)[1], my.split.by = "cond", my.col.pal = pal, na.blank = T)

#==========================================================================================
# PS6 =====================================================================================
#==========================================================================================
ieg_genes = read.csv("C:/Users/miles/Downloads/ieg_like_fos_egr1_npas4_detected_011521.csv", stringsAsFactors = F)[,1]
prog_genes = read.csv("C:/Users/miles/Downloads/progenitor_sox2_nes_coexpress_051721.csv", stringsAsFactors = F)[,1]
neurogen_genes = read.csv("C:/Users/miles/Downloads/neurogen_genes_final_050621.csv", stringsAsFactors = F)[,1]

# Calculate Score
mat = bb@assays$RNA@counts
mat[which(mat > 1)] = 1

bb$ieg_score = colSums(mat[ieg_genes, ])
bb$ieg_score_norm = bb$ieg_score / bb$nFeature_RNA

# ieg_score_norm_95 = colnames(bb)[which(bb$ieg_score_norm >= quantile(bb$ieg_score_norm, 0.95))]
# DimPlot(bb, cells.highlight = ieg_score_norm_95, order = T)
ieg_score_95 = colnames(bb)[which(bb$ieg_score >= quantile(bb$ieg_score, 0.95))]
DimPlot(bb, cells.highlight = ieg_score_9, order = T)

ieg_ps6_full = data.frame()
for (clust15 in c(1, 2, 3)) {
  print(paste0("IEG, Cluster15 = ", clust15))
  bb$pos = F
  bb$pos[ieg_score_95] = T
  bb$pos[which(bb$seuratclusters15 != clust15)] = F
  Idents(bb) = bb$pos
  ieg_ps6 = FindMarkers(bb, ident.1 = T, ident.2 = F, min.pct = 0.01, logfc.threshold = 0)
  ieg_ps6$gene = rownames(ieg_ps6)
  ieg_ps6$cluster = clust15
  ieg_ps6$cat = "ieg"
  ieg_ps6$n_1 = length(which(bb$pos == T))
  ieg_ps6_full = rbind(ieg_ps6_full, ieg_ps6)
}
write.csv(ieg_ps6_full, "C:/Users/miles/Downloads/brain/results/bb/ps6_ieg_15.csv")

ieg_ps6_full_clust53 = data.frame()
for (clust53 in c()) {
  print(paste0("IEG, Cluster53 = ", clust53))
  bb$pos = F
  bb$pos[ieg_score_95] = T
  bb$pos[which(bb$seuratclusters53 != clust53)] = F
  Idents(bb) = bb$pos
  ieg_ps6 = FindMarkers(bb, ident.1 = T, ident.2 = F, min.pct = 0.01, logfc.threshold = 0)
  ieg_ps6$gene = rownames(ieg_ps6)
  ieg_ps6$cluster = clust53
  ieg_ps6$cat = "ieg"
  ieg_ps6$n_1 = length(which(bb$pos == T))
  ieg_ps6_full_clust53 = rbind(ieg_ps6_full_clust53, ieg_ps6)
}



bb$prog_score = colSums(mat[prog_genes, ])
bb$prog_score_norm = bb$prog_score / bb$nFeature_RNA

# prog_score_norm_95 = colnames(bb)[which(bb$prog_score_norm >= quantile(bb$prog_score_norm, 0.95))]
# DimPlot(bb, cells.highlight = prog_score_norm_95, order = T)
prog_score_95 = colnames(bb)[which(bb$prog_score >= quantile(bb$prog_score, 0.95))]
DimPlot(bb, cells.highlight = prog_score_95, order = T)

#==========================================================================================
# Depth and GSI Cor =======================================================================
#==========================================================================================
# 04/02/2021
df = read.csv("C:/Users/miles/Downloads/bb_depth_gsi_cor.csv")
df = df[order(abs(df$depth_cor), decreasing = T),]
df$depth_cor_rank = 1:nrow(df)
df = df[order(abs(df$gsi_cor), decreasing = T),]
df$gsi_cor_rank = 1:nrow(df)
df = df[order(abs(df$depth_cor), decreasing = T),]
df$human = gene_info$human[match(df$gene, gene_info$mzebra)]

my_cor_t = function(r, n) (r * sqrt(n - 2))/sqrt(1 - r**2)
my_cor_p = function(t, n) 2*pt(-abs(t), df=n-2)

df$depth_cor_t = my_cor_t(df$depth_cor, ncol(bb))
df$depth_cor_p = my_cor_p(df$depth_cor_t, ncol(bb))
df$depth_cor_bh = p.adjust(df$depth_cor_p, method = "BH")
df$depth_cor_bon = p.adjust(df$depth_cor_p, method = "bonferroni")

df$gsi_cor_t = my_cor_t(df$gsi_cor, ncol(bb))
df$gsi_cor_p = my_cor_p(df$gsi_cor_t, ncol(bb))
df$gsi_cor_bh = p.adjust(df$gsi_cor_p, method = "BH")
df$gsi_cor_bon = p.adjust(df$gsi_cor_p, method = "bonferroni")
write.csv(df, "C:/Users/miles/Downloads/brain/results/bb/bb_depth_gsi_cor_more_info.csv")

bvcVis(bb, "LOC101469196", mode = "violin_split")
png("C:/Users/miles/Downloads/brain/results/bb/strongest_cor_painting.png", width = 2000, height = 5000, res = 100)
myFeaturePlot(bb, "LOC101469196", my.split.by = "sample", my.pt.size = 1.5, my.col.pal = pal, na.blank = T)
dev.off()

df_pair_dif = data.frame()
df_pair_clust_dif = data.frame()
for (pair in 1:5) {
  print("***")
  print(pair)
  print("***")
  sample_b = colnames(bb)[which(bb$sample == paste0("b", pair))]
  sample_c = colnames(bb)[which(bb$sample == paste0("c", pair))]
  this_df = pct_dif_avg_logFC(bb, cells.1 = sample_b, cells.2 = sample_c)
  this_df$pair = pair
  df_pair_dif = rbind(df_pair_dif, this_df)
  
  for (cluster in levels(bb$seuratclusters53)) {
    print(cluster)
    cluster_cells = colnames(bb)[which(bb$seuratclusters53 == cluster)]
    sample_b_cluster = sample_b[which(sample_b %in% cluster_cells)]
    sample_c_cluster = sample_c[which(sample_c %in% cluster_cells)]
    if (length(sample_b_cluster) > 1 & length(sample_c_cluster) > 1) {
      this_df = pct_dif_avg_logFC(bb, cells.1 = sample_b_cluster, cells.2 = sample_c_cluster)
      this_df$pair = pair
      this_df$cluster = cluster
      df_pair_clust_dif = rbind(df_pair_clust_dif, this_df)
    }
  }
}

#==========================================================================================
# B vs C Network ==========================================================================
#==========================================================================================
# real = read.csv("~/scratch/brain/results/cor_pr_real_strength.csv")
# colnames(real)[1] = "gene"
# real = left_join(data.frame(gene = rownames(bb)), real, by = "gene")
# rownames(real) = real$gene
# real = read.csv("~/scratch/brain/results/py_ns_real.csv")
real = read.csv("~/scratch/brain/results/py_ns_real_abs.csv")
rownames(real) = real$X
real$X = NULL
real = real[which(real$Dif != 0),]

perm = data.frame(gene = rownames(bb))
for (i in 1:20) {
  print(i)
  # perm_small = read.csv(paste0("~/scratch/brain/results/cor_pr_perm/perm_", i, ".csv"))
  # colnames(perm_small)[1] = "gene"
  # perm = left_join(perm, perm_small, by = "gene")
  # perm_small = read.csv(paste0("~/scratch/brain/results/py_ns/perm_", i, ".csv"))
  # perm_small = read.csv(paste0("~/scratch/brain/results/py_ns/perm_abs___", i, ".csv"))
  perm_small = read.csv(paste0("~/scratch/brain/results/py_ns_ieg5/perm_abs_ieg5_", i, ".csv"))
  perm = cbind(perm, perm_small[,-c(1)])
}
rownames(perm) = perm$gene
perm$gene = NULL
perm = data.matrix(perm)
colnames(perm) = c(1:1000)

gene = "LOC101470450"
png(paste0("~/scratch/brain/results/py_ns_bulk_", gene, ".png"), width = 600, height = 400)
# hist_df = data.frame(perm_dif = perm[gene,])
# print(ggplot(hist_df, aes(x=perm_dif)) + geom_histogram(color="black", fill="lightgray") + geom_vline(xintercept = real[gene,"Dif"]) + annotate(x=real[gene,"Dif"],y=+Inf,label="Real",vjust=2,geom="label") + ggtitle(paste0("NS Difference in Permutations Compared to Real for ", gene)) + xlab("NS Difference"))
hist_df = data.frame(idx = 1:1000, perm_dif = unname(t(gene_res)))
print(ggplot(hist_df, aes(x=perm_dif)) + geom_histogram(color="black", fill="lightgray") + geom_vline(xintercept = 273894.59160636) + annotate(x=273894.59160636,y=+Inf,label="Real",vjust=2,geom="label") + ggtitle(paste0("Sum of Absolute NS Difference in Permutations Compared to Real for ", gene)) + xlab("Sum of Abosulte NS Difference"))
dev.off()

test_stat = data.frame(gene = rownames(real), num_greater = 0, pct_greater = 0)
rownames(test_stat) = rownames(real)
for (gene in rownames(real)) {
  this_num_greater = length(which( perm[gene, ] > real[gene, "Dif"] ))
  test_stat[gene,] = c(gene, this_num_greater, this_num_greater / ncol(perm))
}
test_stat[which( is.na(real$Dif) ),c("num_greater", "pct_greater")] = NA

num_greater_table = table(factor(test_stat$num_greater, levels = 0:1000))
test_stat_95 = test_stat[which(test_stat$num_greater > 0.95*1000),]
test_stat_100 = test_stat[which(test_stat$num_greater > 1*1000),]
ieg_cons = c("LOC101487312", "egr1", "npas4", "jun", "homer1")
ieg_like = read.csv(paste0(rna_path, "/results/ieg_like_fos_egr1_npas4_detected_011521.csv"))[,"ieg_like"]
test_stat[ieg_cons,]
test_stat[ieg_like,]
test_stat[which( test_stat$gene %in% ieg_like & (test_stat$pct_greater < 0.05 | test_stat$pct_greater > 0.95) ),]

df = data.frame()
for (gene in rownames(bb)) {
  df = rbind(df, t(c(gene, length(which(bb@assays$RNA@counts[gene,] > 0)))))
}
colnames(df) = c("gene", "num_cells")
df$num_cells = as.numeric(as.vector(df$num_cells))
write.csv(df, "~/scratch/brain/data/gene_num_cells.csv")

png(paste0(rna_path, "/results/py_ns_10k.png"))
hist(as.numeric(test_stat$num_greater[which(! is.na(test_stat$num_greater) )]), xlab = "Number of Permutations that are Greater than the Real B vs C", breaks = 50, col = "lightgray")
dev.off()

png(paste0(rna_path, "/results/py_ns_real_10k.png"))
hist(as.numeric(real$Dif), xlab = "BHVE vs CTRL Node Strength Difference", breaks = 50, col = "lightgray")
dev.off()

png(paste0(rna_path, "/results/py_ns_perm_10k.png"))
hist(as.numeric(rowMeans(perm)), xlab = "Average Permuted BHVE vs CTRL Node Strength Difference", breaks = 50, col = "lightgray")
dev.off()

non_zero = rownames(bb)[which( rowSums(bb@assays$RNA@counts) > 0 )]
mean_df = data.frame(gene = non_zero, real_dif = real[non_zero, "Dif"], perm_mean = rowMeans(perm[non_zero,]))
mean_df2 = melt(mean_df)

png(paste0(rna_path, "/results/cor_res_means_10k.png"))
ggplot(mean_df2, aes(variable, value, fill = variable, color = variable)) + geom_boxplot(alpha = 0.6) + geom_jitter(position=position_dodge2(width = 0.6), alpha = 0.01)
dev.off()


test = rowMeans(perm)
png(paste0(rna_path, "/results/cor_res_means_perm.png"))
hist(test, breaks = 50)
dev.off()

real_means = real$Dif
png(paste0(rna_path, "/results/cor_res_means_real.png"))
hist(real_means, breaks = 50)
dev.off()

#====================#
# 15 Cluster Results #
#====================#
for (i in 0:14) {
  print(i)
  # real = read.csv(paste0("~/scratch/brain/results/py_ns_cluster15_real/cluster15_", i, ".csv"))
  real = read.csv(paste0("~/scratch/brain/results/py_ns_cluster15_real_abs/real_abs_cluster15_", i, "_1.csv"))
  rownames(real) = real$X
  real$X = NULL
  real = real[which(real$Dif != 0),]
  
  print("Reading in Permutations")
  perm = data.frame(gene = rownames(bb))
  for (j in 1:2) {
    # perm_small = read.csv(paste0("~/scratch/brain/results/py_ns_cluster15/cluster15_", i, "_", j, ".csv"))
    perm_small = read.csv(paste0("~/scratch/brain/results/py_ns_cluster15_abs/perm_abs_cluster15_", i, "_", j, ".csv"))
    perm = cbind(perm, perm_small[,-c(1)])
  }
  rownames(perm) = perm$gene
  perm$gene = NULL
  perm = data.matrix(perm)
  colnames(perm) = c(1:1000)
  
  print("Finding the test statistic")
  test_stat = data.frame(gene = rownames(real), num_greater = 0, pct_greater = 0)
  rownames(test_stat) = rownames(real)
  for (gene in rownames(real)) {
    this_num_greater = length(which( perm[gene, ] > real[gene, "Dif"] ))
    test_stat[gene,] = c(gene, this_num_greater, this_num_greater / ncol(perm))
  }
  test_stat[which( is.na(real$Dif) ),c("num_greater", "pct_greater")] = NA
  real = cbind(real, test_stat[,c("num_greater", "pct_greater")])
  
  # Calculating the number of cells in the cluster
  print("Calculating the number of cells in the cluster")
  b_mat = bb@assays$RNA@counts[,which(bb$seuratclusters15 == i & bb$cond == "BHVE")]
  b_mat[which(b_mat > 0)] = 1
  c_mat = bb@assays$RNA@counts[,which(bb$seuratclusters15 == i & bb$cond == "CTRL")]
  c_mat[which(c_mat > 0)] = 1
  num_cells = data.frame(B_n = rowSums(b_mat), C_n = rowSums(c_mat))
  num_cells$num_cells = num_cells$B_n + num_cells$C_n
  real[,c("B_n", "C_n", "n")] = num_cells[match(rownames(real), rownames(num_cells)),]
  real$num_greater = as.numeric(as.vector(real$num_greater))
  p_df = data.frame(low = 2*(abs(1000 - real$num_greater)/1000), high = 2* (1 - (abs(1000-real$num_greater)/1000)))
  real$two.tail.p = apply(p_df, 1, FUN=min)
  real$bh = p.adjust(real$two.tail.p, method = "BH")
  # write.csv(real, paste0("~/scratch/brain/results/py_ns_cluster15_real/cluster15_", i, "_res.csv"))
  write.csv(real, paste0("~/scratch/brain/results/py_ns_cluster15_real_abs/cluster15_", i, "_res.csv"))
  
  print("Subsetting by significant hits")
  num_cells_in_cluster = length(which(bb$seuratclusters15 == i))
  sig = real[which(real$bh <= 0.05 & real$B_n >= 5 & real$C_n >= 5),]
  # write.csv(sig, paste0("~/scratch/brain/results/py_ns_cluster15_real/cluster15_", i, "_res_sig.csv"))
  write.csv(sig, paste0("~/scratch/brain/results/py_ns_cluster15_real_abs/cluster15_", i, "_res_sig.csv"))
  system(paste0("rclone copy ~/scratch/brain/results/py_ns_cluster15_real_abs/cluster15_", i, "_res_sig.csv dropbox:BioSci-Streelman/George/Brain/bb/results/py_ns/abs/cluster15/"))
  system(paste0("rclone copy ~/scratch/brain/results/py_ns_cluster15_real_abs/cluster15_", i, "_res.csv dropbox:BioSci-Streelman/George/Brain/bb/results/py_ns/abs/cluster15/"))
}

big_df = data.frame()
clust_stats = data.frame()
for (i in 0:14) {
  this_df = read.csv(paste0("~/scratch/brain/results/py_ns_cluster15_real_abs/cluster15_", i, "_res_sig.csv"))
  colnames(this_df)[1] = "gene"
  this_df$cluster = i
  print(paste0("Cluster ", i, ": ", dim(this_df)))
  big_df = rbind(big_df, this_df)
}
big_df$hgnc = gene_info$human[match(big_df$gene, gene_info$mzebra)]
write.csv(big_df, "~/scratch/brain/results/py_ns_cluster15_all.csv")
system(paste0("rclone copy ~/scratch/brain/results/py_ns_cluster15_all.csv dropbox:BioSci-Streelman/George/Brain/bb/results/py_ns/abs/cluster15"))

big_df$cat = factor(big_df$cat)
big_df_b_up = big_df[which(big_df$Dif > 0),]
big_df_c_up = big_df[which(big_df$Dif < 0),]
png("~/scratch/brain/results/py_ns_cluster15_stats.png", width = 1500, height = 500, res = 90)
print(ggplot(big_df, aes(x=cluster, fill = Dif > 0, color = Dif > 0)) + geom_bar(alpha = 0.6) + guides(fill=guide_legend(title="BHVE Up"), color=guide_legend(title="BHVE Up")))
dev.off()

png("~/scratch/brain/results/py_ns_cluster15_b_up_stats.png", width = 1500, height = 500, res = 90)
print(ggplot(big_df_b_up, aes(x=cluster, fill = cat, color = cat)) + geom_bar(alpha = 0.6) + ggtitle("BHVE Up") + scale_color_discrete(drop=FALSE, limits=levels(big_df_b_up$cat), name = "Category") + scale_fill_discrete(drop=FALSE, limits=levels(big_df_b_up$cat), name = "Category"))
dev.off()

png("~/scratch/brain/results/py_ns_cluster15_c_up_stats.png", width = 1500, height = 500, res = 90)
print(ggplot(big_df_c_up, aes(x=cluster, fill = cat, color = cat)) + geom_bar(alpha = 0.6) + ggtitle("CTRL Up") + scale_color_discrete(drop=FALSE, limits=levels(big_df_c_up$cat), name = "Category") + scale_fill_discrete(drop=FALSE, limits=levels(big_df_c_up$cat), name = "Category"))
dev.off()


# Wrap a for loop around this
library("rhdf5")
# h5f = H5Fopen("~/scratch/brain/results/py_cor15/cluster15_1_cor_bhve.h5")
h5f = H5Fopen("~/scratch/brain/results/test/perm_abs.h5")
b_cor = h5f$name
h5closeAll()
h5f = H5Fopen("~/scratch/brain/results/py_cor15/cluslster15_1_cor_ctrl.h5")
c_cor = h5f$name
h5closeAll()
tf = read.csv("~/scratch/brain/data/tf_mouse_human_mzebra.csv")
num_cells = read.csv("~/scratch/brain/data/gene_num_cells.csv")
clust1_res = read.csv("~/scratch/brain/results/py_ns_cluster15_real/cluster15_1_res.csv")
tf$num_cells = clust1_res$num_cells[match(tf$mzebra, clust1_res$X)]
tf_mz = tf$mzebra[which(tf$mzebra != "NA" & tf$num_cells >= 10)]

colnames(b_cor) = num_cells$gene
rownames(b_cor) = num_cells$gene
colnames(c_cor) = num_cells$gene
rownames(c_cor) = num_cells$gene

clust1_sig = read.csv("~/scratch/brain/results/py_ns_cluster15_real/cluster15_1_res_sig.csv")

tf_mz_idx = which(num_cells$gene %in% tf_mz)
clust_idx = which(num_cells$gene %in% clust1_sig$X)
tf_clust_idx = unique(c(tf_mz_idx, clust_idx))
b_tf = b_cor[clust_idx, clust_idx]
c_tf = c_cor[clust_idx, clust_idx]
b_tf[which(upper.tri(b_tf, diag = T))] = NA
c_tf[which(upper.tri(c_tf, diag = T))] = NA
b_tf_melt = reshape2::melt(b_tf)
c_tf_melt = reshape2::melt(c_tf)
b_tf_melt = b_tf_melt[which(! is.na(b_tf_melt$value) ),]
c_tf_melt = c_tf_melt[which(! is.na(c_tf_melt$value) ),]

# Change the Gene Names to HGNC
gene_info = read.table("~/scratch/m_zebra_ref/gene_info.txt", sep="\t", header = T, stringsAsFactors = F)
b_tf_melt_hgnc = b_tf_melt
b_tf_melt_hgnc$Var1 = gene_info$human[match(b_tf_melt_hgnc$Var1, gene_info$mzebra)]
b_tf_melt_hgnc$Var2 = gene_info$human[match(b_tf_melt_hgnc$Var2, gene_info$mzebra)]

# If there's no HGNC name, use the LOC's
b_tf_melt_hgnc$Var1[which(is.na(b_tf_melt_hgnc$Var1))] = as.vector(b_tf_melt$Var1[which(is.na(b_tf_melt_hgnc$Var1))])
b_tf_melt_hgnc$Var2[which(is.na(b_tf_melt_hgnc$Var2))] = as.vector(b_tf_melt$Var2[which(is.na(b_tf_melt_hgnc$Var2))])

# Do the same for control
c_tf_melt_hgnc = c_tf_melt
c_tf_melt_hgnc$Var1 = gene_info$human[match(c_tf_melt_hgnc$Var1, gene_info$mzebra)]
c_tf_melt_hgnc$Var2 = gene_info$human[match(c_tf_melt_hgnc$Var2, gene_info$mzebra)]
c_tf_melt_hgnc$Var1[which(is.na(c_tf_melt_hgnc$Var1))] = as.vector(c_tf_melt$Var1[which(is.na(c_tf_melt_hgnc$Var1))])
c_tf_melt_hgnc$Var2[which(is.na(c_tf_melt_hgnc$Var2))] = as.vector(c_tf_melt$Var2[which(is.na(c_tf_melt_hgnc$Var2))])

# Create a dataframe for node information
node_df = data.frame(idx = clust_idx, gene = num_cells$gene[clust_idx], isTF = T, col = "#d3d3d3")
node_df$isTF[which(! node_df$gene %in% tf_mz )] = F
node_df$col[which(node_df$gene %in% clust1_sig$X & clust1_sig$Dif > 0)] = "#e8eb34"
node_df$col[which(node_df$gene %in% clust1_sig$X & clust1_sig$Dif < 0)] = "#3295a8" 
node_df$b_size = clust1_sig$B[match(node_df$gene, clust1_sig$X)]
node_df$c_size = clust1_sig$C[match(node_df$gene, clust1_sig$X)]
# node_df$b_size = clust1_res$B[match(node_df$gene, clust1_res$X)]
# node_df$c_size = clust1_res$C[match(node_df$gene, clust1_res$X)]
node_df$Id = gene_info$human[match(node_df$gene, gene_info$mzebra)]
node_df$Id[which(is.na(node_df$Id))] = node_df$gene[which(is.na(node_df$Id))]
node_df$Label = node_df$Id
node_df = node_df[,c("Id", "Label", colnames(node_df)[1:(ncol(node_df)-2)])]

# Write Nodes and Edges to a file
colnames(b_tf_melt_hgnc) = c("Source", "Target", "Value")
colnames(c_tf_melt_hgnc) = c("Source", "Target", "Value")
write.csv(b_tf_melt_hgnc, "~/scratch/brain/grn/cluster15_1_network/b_no_tf_edges.csv", row.names = F)
write.csv(c_tf_melt_hgnc, "~/scratch/brain/grn/cluster15_1_network/c_no_tf_edges.csv", row.names = F)
write.csv(node_df, "~/scratch/brain/grn/cluster15_1_network/nodes_no_tf.csv", row.names = F)
system(paste0("rclone copy ~/scratch/brain/grn/cluster15_1_network/b_no_tf_edges.csv dropbox:BioSci-Streelman/George/Brain/bb/results/py_ns/cluster15/cluster1_network"))
system(paste0("rclone copy ~/scratch/brain/grn/cluster15_1_network/c_no_tf_edges.csv dropbox:BioSci-Streelman/George/Brain/bb/results/py_ns/cluster15/cluster1_network"))
system(paste0("rclone copy ~/scratch/brain/grn/cluster15_1_network/nodes_no_tf.csv dropbox:BioSci-Streelman/George/Brain/bb/results/py_ns/cluster15/cluster1_network"))

# Write Nodes and Edges to a file - Bulk
colnames(b_tf_melt_hgnc) = c("Source", "Target", "Value")
colnames(c_tf_melt_hgnc) = c("Source", "Target", "Value")
write.csv(b_tf_melt_hgnc, "~/scratch/brain/grn/bulk_network/b_no_tf_edges.csv", row.names = F)
write.csv(c_tf_melt_hgnc, "~/scratch/brain/grn/bulk_network/c_no_tf_edges.csv", row.names = F)
write.csv(node_df, "~/scratch/brain/grn/bulk_network/nodes_no_tf.csv", row.names = F)
system(paste0("rclone copy ~/scratch/brain/grn/bulk_network/b_no_tf_edges.csv dropbox:BioSci-Streelman/George/Brain/bb/results/py_ns/abs/bulk/bulk_network"))
system(paste0("rclone copy ~/scratch/brain/grn/bulk_network/c_no_tf_edges.csv dropbox:BioSci-Streelman/George/Brain/bb/results/py_ns/abs/bulk/bulk_network"))
system(paste0("rclone copy ~/scratch/brain/grn/bulk_network/nodes_no_tf.csv dropbox:BioSci-Streelman/George/Brain/bb/results/py_ns/abs/bulk/bulk_network"))

bvc_tf = bvc[clust_idx, clust_idx]
bvc_tf[which(upper.tri(bvc_tf, diag = T))] = NA
bvc_tf_melt = reshape2::melt(bvc_tf)
bvc_tf_melt = bvc_tf_melt[which(! is.na(bvc_tf_melt$value) ),]

bvc_p_tf = bvc_p[clust_idx, clust_idx]
bvc_p_tf[which(upper.tri(bvc_p_tf, diag = T))] = NA
bvc_p_tf_melt = reshape2::melt(bvc_p_tf)
bvc_p_tf_melt = bvc_p_tf_melt[which(! is.na(bvc_p_tf_melt$value) ),]

# Make a big heatmap
png("~/scratch/brain/grn/cluster15_1_network/b_tf_heatmap.png", width = nrow(b_tf)*33, height = ncol(b_tf)*33)
pheatmap::pheatmap(b_tf, cellwidth = 30, cellheight = 30, cluster_rows = F, cluster_cols = F)
dev.off()
system(paste0("rclone copy ~/scratch/brain/grn/cluster15_1_network/b_tf_heatmap.png dropbox:BioSci-Streelman/George/Brain/bb/results/py_ns/cluster15/cluster1_network"))

sig_tf_df =  read.csv("C:/Users/miles/Downloads/brain/results/bb/bb_all_interesting_genes_in_ns15.csv", stringsAsFactors = F)
sig_tf = sig_tf_df$gene[which(sig_tf_df$class == "tf" | sig_tf_df$class2 == "tf")]

in_path = "~/scratch/brain/results/py_ns_cluster15_real_abs/"
# in_path = "C:/Users/miles/Downloads/brain/results/bb/cluster15/"
sig_tf_heat = data.frame()
for (i in convert15$new.full) {
  j = convert15$old[which(convert15$new.full == i)]
  ns_df = read.csv(paste0(in_path, "cluster15_", j, "_res.csv"), stringsAsFactors = F)
  ns_df = ns_df[which(ns_df$X %in% sig_tf),]
  # ns_df = ns_df[which(ns_df$num_cells >= 5),]
  ns_df = ns_df[which(ns_df$n >= 5),]
  ns_df$cluster = i
  ns_df$nMoreExtreme = 500-abs(ns_df$num_greater-500)
  sig_tf_heat = rbind(sig_tf_heat, ns_df)
}
sig_tf_mat = acast(sig_tf_heat, cluster ~ X, value.var = "num_greater")
sig_tf_mat = sig_tf_mat[order(as.numeric(convert15$new.num), decreasing = T), ]
colnames(sig_tf_mat) = sig_tf_df$label[match(colnames(sig_tf_mat), sig_tf_df$gene)]
sig_tf_mat[which(is.na(sig_tf_mat))] = 500
isSig = matrix("", nrow = nrow(sig_tf_mat), ncol = ncol(sig_tf_mat))
isSig[which(sig_tf_mat <= 2 | sig_tf_mat > 998)] = "*"
# pheatmap::pheatmap(sig_tf_mat, cellwidth = 20, cellheight = 20, cluster_rows = F, cluster_cols = T, angle_col = "315", display_numbers = isSig, number_color = "black", fontsize_number = 14, filename = "C:/Users/miles/Downloads/brain/results/bb/py_ns_cluster15_tf_sig_heat.pdf")
pheatmap::pheatmap(sig_tf_mat, show_colnames = T, cellwidth = 20, cellheight = 20, cluster_rows = F, cluster_cols = T, angle_col = "315", display_numbers = isSig, number_color = "black", fontsize_number = 14, color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100), filename = "~/scratch/brain/results/py_ns_abs_cluster15_tf_sig_heat.pdf")

# Plot num_greater vs Dif
p_list = list()
for (i in 0:14) {
  ns_df = read.csv(paste0(in_path, "cluster15_", i, "_res.csv"), stringsAsFactors = F)
  ns_df = ns_df[which(ns_df$num_cells >= 5),]
  ns_df$cluster = i
  ns_df$nMoreExtreme = 500-abs(ns_df$num_greater-500)
  ns_df$sig = ns_df$bh < 0.05
  # print(ggplot(ns_df, aes(x = Dif, y = num_greater, color = sig)) + geom_point() + scale_colour_viridis_d(begin = 1, end = 0, name = "Sig", limits = c(T, F)) + geom_text_repel(data=ns_df[which(ns_df$sig),], aes(label=X),hjust=0, vjust=0, color  = "black"))
  p_list[[i+1]] = ggplot(ns_df, aes(x = Dif, y = num_greater, color = sig)) + geom_point() + scale_colour_viridis_d(begin = 1, end = 0, name = "Significant", limits = c(T, F)) + ylab("Number of Permutations Greater Than Real") + xlab("BHVE NS - CTRL NS") + ggtitle(paste0("Node Strength (NS) Difference in Cluster ", i))
}
p = plot_grid(plotlist=p_list, ncol = 5)
pdf("C:/Users/miles/Downloads/brain/results/bb/py_ns_cluster15_dif.pdf", width = 50, height = 20)
print(plot_grid(plotlist=p_list, ncol = 5))
dev.off()

# Make DotPlots
big_df2_top = data.frame()
for ( i in 0:10 ) {
  this_row = big_df2[which(big_df2$cluster == i),]
  this_row = this_row[which(this_row$Dif < 0),]
  if (nrow(this_row) >= 5) {
    this_row = this_row[order(this_row$Dif, decreasing = F),]
    newRow = data.frame(cluster = rep(this_row$cluster[1],10), cond = rep(c("Behave", "Control"),5), 
                        node_strength = as.vector(sapply(1:5, function(x) c(this_row$B[x], this_row$C[x]) )),
                        dif = as.vector(sapply(1:5, function(x) rep(this_row$Dif[x], 2))),
                        gene = as.vector(sapply(1:5, function(x) rep(this_row$label[x], 2))),
                        num = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5))
    big_df2_top = rbind(big_df2_top, newRow)
  } else {
    print(paste("Cluster", i, "didn't have 5 genes where BHVE is UP.")) 
  }
  
}
big_df2_top$cluster = factor(big_df2_top$cluster, levels = 0:14)
ggplot(big_df2_top, aes(x = num, y=cluster, size = node_strength, color = dif, label = gene)) + geom_point() + geom_text(aes(label=gene),hjust=0, vjust=0) + facet_wrap( ~ cond)
ggplot(big_df2_top, aes(x = cond, y=cluster, size = node_strength, color = dif)) + geom_point() + geom_text(aes(label=gene),hjust=1.2, vjust=0.5, size = 3, color = "black") + facet_wrap( ~ num) + ggtitle("Top 5 Genes with Significantly Greater Control NS") + guides(color=guide_legend(title="NS Dif"), size=guide_legend(title="NS")) + scale_color_continuous(high = "#132B43", low = "#56B1F7")

# Score Cor Results - Bulk
# num_cells = write.csv("~/scratch/brain/data/gene_num_cells.csv", stringsAsFactors = F)
print("Calculating number of cells")
mat = bb@assays$RNA@counts[,which(bb$cond == "BHVE")]
mat[which(mat > 0)] = 1
bhve_num_cells = rowSums(mat)
mat = bb@assays$RNA@counts[,which(bb$cond == "CTRL")]
mat[which(mat > 0)] = 1
ctrl_num_cells = rowSums(mat)
# gene_lists = c("ieg", "prog", "neurogen")
gene_lists = c("pcrclg11")
for (gene_list in gene_lists) {
  real = read.csv(paste0("~/scratch/brain/results/bulk_real_", gene_list, "_score_cor_bvc.csv"), stringsAsFactors = F)
  rownames(real) = real$X
  real$X = NULL
  real = real[which(real$Dif != 0),]
  perm = read.csv(paste0("~/scratch/brain/results/bulk_perm_1000_", gene_list, "_score_cor_bvc.csv"), stringsAsFactors = F)
  rownames(perm) = perm$X
  perm$X = NULL
  
  real$nMoreExtreme = perm$nMoreExtreme[match(rownames(real), rownames(perm))]
  # real$num_cells = num_cells$num_cells[match(rownames(real), num_cells$gene)]
  real$bhve_num_cells = bhve_num_cells[match(rownames(real), rownames(bb))]
  real$ctrl_num_cells = ctrl_num_cells[match(rownames(real), rownames(bb))]
  
  real$perm.two.tail.p = 2 * (real$nMoreExtreme/1000)
  real$perm.bh = p.adjust(real$perm.two.tail.p, method = "BH")
  real$fisher.p = r_to_p(real$BHVE, real$CTRL, ncol(bb), ncol(bb))
  real$fisher.bh = p.adjust(real$fisher.p, method = "BH") 
  write.csv(real, paste0("~/scratch/brain/results/bulk_real_", gene_list, "_score_cor_bvc_w_stats.csv"))
  system(paste0("rclone copy ~/scratch/brain/results/bulk_real_", gene_list, "_score_cor_bvc_w_stats.csv dropbox:BioSci-Streelman/George/Brain/bb/results/score_cor/bulk/"))
  
  # real_sig = real[which(real$bh < 0.05 & real$num_cells >= 10),]
  real_sig = real[which(real$perm.bh < 0.05 & real$fisher.bh < 0.05 & real$bhve_num_cells >= 5 & real$ctrl_num_cells >= 5),]
  write.csv(real_sig, paste0("~/scratch/brain/results/bulk_real_", gene_list, "_score_cor_bvc_w_stats_sig.csv"))
  system(paste0("rclone copy ~/scratch/brain/results/bulk_real_", gene_list, "_score_cor_bvc_w_stats_sig.csv dropbox:BioSci-Streelman/George/Brain/bb/results/score_cor/bulk/"))
  
  png(paste0("~/scratch/brain/results/bulk_score_cor_", gene_list, "_extreme.png"))
  hist(real$nMoreExtreme, breaks = 50, xlab = "Number of Correlations More Extreme than Permutations", main = "")
  dev.off()
  system(paste0("rclone copy ~/scratch/brain/results/bulk_score_cor_", gene_list, "_extreme.png dropbox:BioSci-Streelman/George/Brain/bb/results/score_cor/bulk/"))
  
  print(gene_list)
  print(paste0("Number of significant bulk results: ", nrow(real_sig)))
  print(paste0("Max Significant Difference: ", rownames(real_sig)[which.max(abs(real_sig$Dif))], ", ", max(abs(real_sig$Dif))))
}


# Score Cor Results - 15 Cluster Level
# gene_lists = c("ieg", "prog", "neurogen")
gene_lists = c("pcrclg11")
for (gene_list in gene_lists) {
  print(gene_list)
  real = read.csv(paste0("~/scratch/brain/results/clust15_real_", gene_list, "_score_cor_bvc.csv"), stringsAsFactors = F)
  rownames(real) = real$X
  real$X = NULL
  colnames(real) = 0:14
  perm = read.csv(paste0("~/scratch/brain/results/clust15_perm_1000_", gene_list, "_score_cor_bvc.csv"), stringsAsFactors = F)
  rownames(perm) = perm$X
  perm$X = NULL
  colnames(perm) = 0:14
  
  real_full = read.csv(paste0("~/scratch/brain/results/clust15_real_full_", gene_list, "_score_cor_bvc.csv"), stringsAsFactors = F)
  rownames(real_full) = real_full$X
  real_full$X = NULL
  
  clust15_sig = data.frame()
  clust15_full = data.frame()
  for (i in 0:14) {
    # print(i)
    # print("Calculating the number of cells in the cluster")
    num_clust_cells = length(which(bb$seuratclusters15 == i))
    mat = bb@assays$RNA@counts[,which(bb$seuratclusters15 == i & bb$cond == "BHVE")]
    mat[which(mat > 0)] = 1
    bhve_num_cells = rowSums(mat)
    mat = bb@assays$RNA@counts[,which(bb$seuratclusters15 == i & bb$cond == "CTRL")]
    mat[which(mat > 0)] = 1
    ctrl_num_cells = rowSums(mat)
    
    # this_real = data.frame(Dif = real[,as.character(i)], num_cells = num_cells)
    this_real = data.frame(BHVE = real_full[,paste0("BHVE_",i)], CTRL = real_full[,paste0("CTRL_",i)], Dif = real[,as.character(i)], bhve_num_cells = bhve_num_cells, ctrl_num_cells = ctrl_num_cells)
    rownames(this_real) = rownames(real)
    this_real = this_real[which(this_real$Dif != 0),]
    
    this_real$nMoreExtreme = perm[match(rownames(this_real), rownames(perm)), as.character(i)]
    this_real$two.tail.p = 2 * (this_real$nMoreExtreme/1000)
    this_real$perm.bh = p.adjust(this_real$two.tail.p, method = "BH")
    this_real$fisher.p = r_to_p(this_real$BHVE, this_real$CTRL, num_clust_cells, num_clust_cells)
    this_real$fisher.bh = p.adjust(this_real$fisher.p, method = "BH")
    
    this_real$gene = rownames(this_real)
    this_real$cluster = i
    rownames(this_real) = paste0(rownames(this_real), "_", i)
    clust15_full = rbind(clust15_full, this_real)
    
    # this_real_sig = this_real[which(this_real$bh < 0.05 & this_real$num_cells >= 10),]
    this_real_sig = this_real[which(this_real$perm.bh < 0.05 & this_real$fisher.bh < 0.05 & this_real$bhve_num_cells >= 5 & this_real$ctrl_num_cells),]
    if ( nrow(this_real_sig) > 0 ) {
      clust15_sig = rbind(clust15_sig, this_real_sig)
    }
    
    print(paste0("Number of significant results for cluster ", i, ": ", nrow(this_real_sig)))
    # if ( nrow(this_real_sig) > 0 )
    #   print(paste0("Max Significant Difference: ", rownames(this_real_sig)[which.max(abs(this_real_sig$Dif))], ", ", max(abs(this_real_sig$Dif))))
  } # end cluster for
  print(paste0("Total Number of 15 Cluster Significant results for ", gene_list, ": ", nrow(clust15_sig)))
  write.csv(clust15_sig, paste0("~/scratch/brain/results/clust15_real_", gene_list, "_score_cor_bvc_w_stats_sig.csv"))
  system(paste0("rclone copy ~/scratch/brain/results/clust15_real_", gene_list, "_score_cor_bvc_w_stats_sig.csv dropbox:BioSci-Streelman/George/Brain/bb/results/score_cor/clust15/"))
  
  # Across Cluster BH Correction
  clust15_full$full.fisher.bh = p.adjust(clust15_full$fisher.p, method = "BH")
  clust15_full$full.perm.bh = p.adjust(clust15_full$two.tail.p, method = "BH")
  clust15_full_sig = clust15_full[which(clust15_full$full.fisher.bh < 0.05 & clust15_full$full.perm.bh < 0.05 & clust15_full$bhve_num_cells >= 5 & clust15_full$ctrl_num_cells >= 5),]
  print(paste0("Total Number of FULL 15 Cluster Significant results for ", gene_list, ": ", nrow(clust15_full_sig)))
  write.csv(clust15_sig, paste0("~/scratch/brain/results/clust15_real_full_", gene_list, "_score_cor_bvc_w_stats_sig.csv"))
  system(paste0("rclone copy ~/scratch/brain/results/clust15_real_full_", gene_list, "_score_cor_bvc_w_stats_sig.csv dropbox:BioSci-Streelman/George/Brain/bb/results/score_cor/clust15/"))
} # end gene_list for

# Score Cor Results - 53 Cluster Level
# gene_lists = c("ieg", "prog", "neurogen")
gene_lists = c("pcrclg11")
for (gene_list in gene_lists) {
  print(gene_list)
  real = read.csv(paste0("~/scratch/brain/results/clust53_real_", gene_list, "_score_cor_bvc.csv"), stringsAsFactors = F)
  rownames(real) = real$X
  real$X = NULL
  colnames(real) = 0:52
  perm = read.csv(paste0("~/scratch/brain/results/clust53_perm_1000_", gene_list, "_score_cor_bvc.csv"), stringsAsFactors = F)
  rownames(perm) = perm$X
  perm$X = NULL
  colnames(perm) = 0:52
  
  real_full = read.csv(paste0("~/scratch/brain/results/clust53_real_full_", gene_list, "_score_cor_bvc.csv"), stringsAsFactors = F)
  rownames(real_full) = real_full$X
  real_full$X = NULL
  
  clust53_sig = data.frame()
  clust53_full = data.frame()
  for (i in 0:52) {
    # print(i)
    # print("Calculating the number of cells in the cluster")
    num_clust_cells = length(which(bb$seuratclusters53 == i))
    mat = bb@assays$RNA@counts[,which(bb$seuratclusters53 == i & bb$cond == "BHVE")]
    mat[which(mat > 0)] = 1
    bhve_num_cells = rowSums(mat)
    mat = bb@assays$RNA@counts[,which(bb$seuratclusters53 == i & bb$cond == "CTRL")]
    mat[which(mat > 0)] = 1
    ctrl_num_cells = rowSums(mat)
    
    # this_real = data.frame(Dif = real[,as.character(i)], num_cells = num_cells)
    this_real = data.frame(BHVE = real_full[,paste0("BHVE_",i)], CTRL = real_full[,paste0("CTRL_",i)], Dif = real[,as.character(i)], bhve_num_cells = bhve_num_cells, ctrl_num_cells = ctrl_num_cells)
    rownames(this_real) = rownames(real)
    this_real = this_real[which(this_real$Dif != 0),]
    
    this_real$nMoreExtreme = perm[match(rownames(this_real), rownames(perm)), as.character(i)]
    this_real$two.tail.p = 2 * (this_real$nMoreExtreme/1000)
    this_real$perm.bh = p.adjust(this_real$two.tail.p, method = "BH")
    this_real$fisher.p = r_to_p(this_real$BHVE, this_real$CTRL, num_clust_cells, num_clust_cells)
    this_real$fisher.bh = p.adjust(this_real$fisher.p, method = "BH")
    
    this_real$gene = rownames(this_real)
    this_real$cluster = i
    rownames(this_real) = paste0(rownames(this_real), "_", i)
    clust53_full = rbind(clust53_full, this_real)
    
    # this_real_sig = this_real[which(this_real$bh < 0.05 & this_real$num_cells >= 10),]
    this_real_sig = this_real[which(this_real$perm.bh < 0.05 & this_real$fisher.bh < 0.05 & this_real$bhve_num_cells >= 5 & this_real$ctrl_num_cells),]
    if ( nrow(this_real_sig) > 0 ) {
      clust53_sig = rbind(clust53_sig, this_real_sig)
    }
    
    # print(paste0("Number of significant results for cluster ", i, ": ", nrow(this_real_sig)))
    # if ( nrow(this_real_sig) > 0 )
    #   print(paste0("Max Significant Difference: ", rownames(this_real_sig)[which.max(abs(this_real_sig$Dif))], ", ", max(abs(this_real_sig$Dif))))
  } # end cluster for
  print(paste0("Total Number of 53 Cluster Significant results for ", gene_list, ": ", nrow(clust53_sig)))
  write.csv(clust53_sig, paste0("~/scratch/brain/results/clust53_real_", gene_list, "_score_cor_bvc_w_stats_sig.csv"))
  system(paste0("rclone copy ~/scratch/brain/results/clust53_real_", gene_list, "_score_cor_bvc_w_stats_sig.csv dropbox:BioSci-Streelman/George/Brain/bb/results/score_cor/clust53/"))
  
  # Across Cluster BH Correction
  clust53_full$full.fisher.bh = p.adjust(clust53_full$fisher.p, method = "BH")
  clust53_full$full.perm.bh = p.adjust(clust53_full$two.tail.p, method = "BH")
  clust53_full_sig = clust53_full[which(clust53_full$full.fisher.bh < 0.05 & clust53_full$full.perm.bh < 0.05 & clust53_full$bhve_num_cells >= 5 & clust53_full$ctrl_num_cells >= 5),]
  print(paste0("Total Number of FULL 53 Cluster Significant results for ", gene_list, ": ", nrow(clust53_full_sig)))
  write.csv(clust53_sig, paste0("~/scratch/brain/results/clust53_real_full_", gene_list, "_score_cor_bvc_w_stats_sig.csv"))
  system(paste0("rclone copy ~/scratch/brain/results/clust53_real_full_", gene_list, "_score_cor_bvc_w_stats_sig.csv dropbox:BioSci-Streelman/George/Brain/bb/results/score_cor/clust53/"))
} # end gene_list for


# png("~/scratch/brain/results/score_cor_ieg_sig.png")
# ggplot(real_sig, aes(Dif, nMoreExtreme, color = num_cells)) + geom_point()
# dev.off()

#==========================================================================================
# Split Into Individuals ==================================================================
#==========================================================================================
all_p = read.csv("C:/Users/miles/Downloads/all_p.csv", header = F, stringsAsFactors = F)
colnames(all_p) = c("barcode", "1", "2", "3", "4", "sample")
# substr_cells = substr(colnames(bb),6,23)
all_p$cell = sapply(1:nrow(all_p), function(x) colnames(bb)[which(bb$sample == all_p$sample[x] & substr_cells == all_p$barcode[x])])
all_p$best = sapply(1:nrow(all_p), function(x) colnames(all_p)[which.max(as.numeric(as.vector(all_p[x,2:5])))+1] )
all_p$best_p = sapply(1:nrow(all_p), function(x) max(as.numeric(as.vector(all_p[x,2:5]))) )

bb$sub = all_p$best[match(colnames(bb), all_p$cell)]
bb$p = all_p$best_p[match(colnames(bb), all_p$cell)]
bb$subsample = paste0(bb$sample, ".", bb$sub)
ind_list =  setNames(lapply(sort(unique(bb$subsample)), function(x) colnames(bb)[which(bb$subsample == x)]), sort(unique(bb$subsample)))
bb$subsample = factor(bb$subsample, levels = sort(unique(bb$subsample)))
Idents(bb) = bb$subsample
png("C:/Users/miles/Downloads/bb_subsample.png", width = 6000, height = 3000, res = 100)
DimPlot(bb, pt.size = 1, split.by = "subsample", ncol=10)
dev.off()

Idents(bb) = bb$seuratclusters15
FeaturePlot(bb, "p", order = T, label = T, pt.size = 1)

#==========================================================================================
# Sample ID from Unique SNPs ==============================================================
#==========================================================================================
sample_score = data.frame(sample = factor(bb$sample, levels = unique(bb$sample)))
sample_score$b1 = read.table(paste0(rna_path, "/data/cell_b1.txt"), skip=1, header = F, stringsAsFactors = F)[,2]
sample_score$b2 = read.table(paste0(rna_path, "/data/cell_b2.txt"), skip=1, header = F, stringsAsFactors = F)[,2]
sample_score$b3 = read.table(paste0(rna_path, "/data/cell_b3.txt"), skip=1, header = F, stringsAsFactors = F)[,2]
sample_score$b4 = read.table(paste0(rna_path, "/data/cell_b4.txt"), skip=1, header = F, stringsAsFactors = F)[,2]
sample_score$b5 = read.table(paste0(rna_path, "/data/cell_b5.txt"), skip=1, header = F, stringsAsFactors = F)[,2]
sample_score$c1 = read.table(paste0(rna_path, "/data/cell_c1.txt"), skip=1, header = F, stringsAsFactors = F)[,2]
sample_score$c2 = read.table(paste0(rna_path, "/data/cell_c2.txt"), skip=1, header = F, stringsAsFactors = F)[,2]
sample_score$c3 = read.table(paste0(rna_path, "/data/cell_c3.txt"), skip=1, header = F, stringsAsFactors = F)[,2]
sample_score$c4 = read.table(paste0(rna_path, "/data/cell_c4.txt"), skip=1, header = F, stringsAsFactors = F)[,2]
sample_score$c5 = read.table(paste0(rna_path, "/data/cell_c5.txt"), skip=1, header = F, stringsAsFactors = F)[,2]

# UMAP Plot
bb$b1_score = sample_score$b1
png("C:/Users/miles/Downloads/test.png", width = 4000, height = 750, res = 100)
print(FeaturePlot(bb, "b1_score", order = T, pt.size = 1, label = T, split.by = "sample"))
dev.off()

for (i in 2:ncol(sample_score)) {
  this_sample = colnames(sample_score)[i]
  # temp = c(brewer.pal(10, "Paired")[10], rep(brewer.pal(10, "Paired")[9], 9))
  print(ggplot(sample_score, aes(x=sample_score[,1], y=sample_score[,i], fill = sample_score[,1], color = sample_score[,1])) + geom_boxplot(alpha = 0.6) + NoLegend() + ggtitle(paste("Number of", this_sample, "SNPs")) + ylab(paste("Number of", this_sample, "SNPs")) + xlab("Sample"))
}
# temp = c(brewer.pal(10, "Paired")[c(TRUE, FALSE)], brewer.pal(10, "Paired")[c(FALSE, TRUE)])

sample_score$best = NULL
sample_score$best = sapply(1:nrow(sample_score), function(x) colnames(sample_score)[which.max(as.numeric(as.vector(sample_score[x,2:ncol(sample_score)])))+1] )
sample_score = sample_score[which(sample_score$sample %in% colnames(sample_score)[2:ncol(sample_score)]),]
ggplot(sample_score, aes(x=sample_score[,1], fill = best)) + geom_bar(stat = "count", position = "stack") + xlab("Real Sample") + guides(fill=guide_legend(title="Pred. Sample"))

sample_score$correct = TRUE
sample_score$correct[which(sample_score$sample != sample_score$best)] = FALSE
ggplot(sample_score, aes(x = sample, fill = correct)) + geom_bar(stat = "count", position = "dodge") + xlab("Sample") + guides(fill=guide_legend(title="Pred. Correct"))

b1 = subset(bb, idents = "b1")
b1$pred = "b1"
pred_cells = list()
for (pred_sample in sort(unique(sample_score$best))) {
  pred_cells[[pred_sample]] = rownames(sample_score)[which(sample_score$best == pred_sample)]
  b1$pred[which(colnames(b1) %in% pred_cells[[pred_sample]])] = pred_sample
}

Idents(bb) = bb$sample
DimPlot(subset(bb, idents = "b1"), cells.highlight = pred_cells, cols.highlight = rainbow(length(unique(sample_score$best)) + 1)) + ggtitle("B1") + guides(color=guide_legend(title="Pred. Sample"))

# temp = c(brewer.pal(10, "Paired")[c(TRUE, FALSE)], brewer.pal(10, "Paired")[c(FALSE, TRUE)])
DimPlot(subset(bb, idents = "b1"), cells.highlight = pred_cells, cols.highlight = brewer.pal(10, "Spectral"))
DimPlot(subset(bb, idents = "b1"), cells.highlight = rev(pred_cells), cols.highlight = rev(brewer.pal(11, "Spectral")), order = T, sizes.highlight = 1.2)
DimPlot(subset(bb, idents = "b1"), cells.highlight = pred_cells, cols.highlight = c("gray", rep("blue", 9)), order = T, sizes.highlight = 1.2)

Idents(b1) = b1$pred
DimPlot(b1, order = T, pt.size = 1.2) + ggtitle("Predicted Labels for Real B1 Cells")

svg(paste0(rna_path, "/results/bb/cluster_15_small.svg"), width = 6, height = 5)
new15 = changeClusterID(bb$seuratclusters15, clust_level = 15)
new15 = factor(new15, levels = c("1_Astro/MG", "2_OPC/Oligo", "3_Peri", "4_In", "5_In", "6_In", "7_In", "8_Ex", "9_Ex", "10_Ex", "11_Ex", "12_Ex", "13_Ex", "14_Ex", "15_InEx"))
Idents(bb) = new15
DimPlot(bb, label = T)
dev.off()

#==========================================================================================
# MS Figures ==============================================================================
#==========================================================================================
# DEGs DotPlots
new_clust_15_order = c("1_Astro/MG", "2_OPC/Oligo", "3_Peri", "4_In", "5_In", "6_In", "7_In", "8_Ex", "9_Ex", "10_Ex", "11_Ex", "12_Ex", "13_Ex", "14_Ex", "15_InEx")
deg15 = read.csv("C:/Users/miles/Downloads/bb_all_markers_15clusters_102820_more_info.csv")
deg15$cluster_new = changeClusterID(deg15$cluster, clust_level = 15)
num_x = 5
topX = as.vector(sapply(new_clust_15_order, function(x) deg15$gene[which(deg15$cluster_new == x)[1:num_x]]))
topX = unique(topX)

bb = ScaleData(bb, features = topX)
new15 = changeClusterID(bb$seuratclusters15, clust_level = 15)
new15 = factor(new15, levels = c("1_Astro/MG", "2_OPC/Oligo", "3_Peri", "4_In", "5_In", "6_In", "7_In", "8_Ex", "9_Ex", "10_Ex", "11_Ex", "12_Ex", "13_Ex", "14_Ex", "15_InEx"))
Idents(bb) = new15

library("scales")
xtext_unique_cols = gc.ramp <- hue_pal()(15)
xtext_cols = sapply(xtext_unique_cols, function(x) rep(x, num_x))
xtext_cols = as.vector(xtext_cols)

jungle_cols = rev(c("#f94144", "#f3722c", "#f8961e", "#f9c74f", "#79bf43", "#43aa8b", "#577590"))
cyto_cols = rev(brewer.pal(11,"Spectral"))
cyto_cols[6] = "gold" # this is what plotCytoTRACE uses

pdf(paste0(rna_path, "/results/bb/top5_default.pdf"), width = 18, height = 6)
print(DotPlot(bb, features = rev(topX)) + ylab("Cluster") + xlab("Gene") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = xtext_cols), axis.text.y = element_text(colour = xtext_unique_cols)))
dev.off()

pdf(paste0(rna_path, "/results/bb/top5_my_jungle.pdf"), width = 18, height = 6)
pal = colorRampPalette(jungle_cols)
print(DotPlot(bb, features = rev(topX)) + scale_color_gradientn(colors = pal(50)) + ylab("Cluster") + xlab("Gene") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = xtext_cols), axis.text.y = element_text(colour = xtext_unique_cols)))
dev.off()

pdf(paste0(rna_path, "/results/bb/top5_my_cyto.pdf"), width = 18, height = 6)
pal = colorRampPalette(cyto_cols)
print(DotPlot(bb, features = rev(topX)) + scale_color_gradientn(colors = pal(50)) + ylab("Cluster") + xlab("Gene") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = xtext_cols), axis.text.y = element_text(colour = xtext_unique_cols)))
dev.off()

# DEG Heatmaps
pdf(paste0(rna_path, "/results/bb/heatmap_top5_default.pdf"), width = 10, height = 10)
print(DoHeatmap(bb, features = rev(topX)))
dev.off()

pdf(paste0(rna_path, "/results/bb/heatmap_top5_my_v.pdf"), width = 18, height = 6)
pal = colorRampPalette(jungle_cols)
print(markerHeatmap(bb, markers = topX) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = xtext_cols), axis.text.y = element_text(colour = xtext_unique_cols)))
dev.off()

# Violin Plots
library(Seurat)
library(patchwork)
library(ggplot2)
modify_vlnplot<- function(obj, 
                          feature, 
                          xtext_col,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0, colour = xtext_col), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          xtext_unique_cols = NULL,
                          xtext_cols = NULL,
                          ...) {
  
  plot_list<- purrr::map2(features, xtext_cols, function(x, y) modify_vlnplot(obj = obj,feature = x, xtext_col = y, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1, colour = xtext_unique_cols), 
          axis.ticks.x = element_line(),
          axis.text.y = element_text(colour = xtext_cols)
          )
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  
  # # change y-axis title text colors
  # plot_list<- purrr::map2(plot_list, xtext_cols, function(x,y) {
  #   x + theme()
  # })
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

pdf(paste0(rna_path, "/results/bb/vlnplot_top5.pdf"), width = 10, height = 18)
print(StackedVlnPlot(bb, rev(topX), xtext_cols = rev(xtext_cols), xtext_unique_cols = xtext_unique_cols))
dev.off()

# Specific Markers
deg15$pct_dif = as.numeric(as.vector(deg15$pct.1)) - as.numeric(as.vector(deg15$pct.2))
rank = c()
for (i in 0:14) {
  this_clust = deg15[which(deg15$cluster == i),]
  rank = c(rank, 1:nrow(this_clust))
}
deg15$rank = rank
deg15_order = deg15[order(deg15$pct.2),]
num_x = 5
topX_s = as.vector(sapply(new_clust_15_order, function(x) deg15_order$gene[which(deg15_order$cluster_new == x & deg15_order$n_gene_appears == 1)[1:num_x]]))
topX_s2 = as.vector(sapply(new_clust_15_order, function(x) deg15$gene[which(deg15$cluster_new == x & deg15$n_gene_appears == 1)[1:num_x]]))


pdf(paste0(rna_path, "/results/bb/specific2_top5_my_jungle.pdf"), width = 18, height = 6)
pal = colorRampPalette(jungle_cols)
print(DotPlot(bb, features = rev(topX_s2)) + scale_color_gradientn(colors = pal(50)) + ylab("Cluster") + xlab("Gene") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = xtext_cols), axis.text.y = element_text(colour = xtext_unique_cols)))
dev.off()

pdf(paste0(rna_path, "/results/bb/specific2_vlnplot_top5.pdf"), width = 10, height = 18)
print(StackedVlnPlot(bb, rev(topX_s2), xtext_cols = rev(xtext_cols), xtext_unique_cols = xtext_unique_cols))
dev.off()

deg53 = read.csv("C:/Users/miles/Downloads/bb_all_markers_53clusters_102720_more_info.csv")
deg53$pct_dif = as.numeric(as.vector(deg53$pct.1)) - as.numeric(as.vector(deg53$pct.2))
rank = c()
for (i in 0:52) {
  this_clust = deg53[which(deg53$cluster == i),]
  rank = c(rank, 1:nrow(this_clust))
}
deg53$rank = rank
deg53$cluster_new = changeClusterID(deg53$cluster)
deg53_order = deg53[order(deg53$pct.2),]
num_x = 5
# topX_s = as.vector(sapply(new_clust_53_order, function(x) deg53_order$gene[which(deg53_order$cluster_new == x & deg53_order$n_gene_appears == 1)[1:num_x]]))
# topX_s2 = as.vector(sapply(new_clust_53_order, function(x) deg53$gene[which(deg53$cluster_new == x & deg53$n_gene_appears == 1)[1:num_x]]))
topX_s = as.vector(sapply(0:52, function(x) deg53_order$gene[which(deg53_order$cluster == x & deg53_order$n_gene_appears == 1)[1:num_x]]))
topX_s2 = as.vector(sapply(0:52, function(x) deg53$gene[which(deg53$cluster == x & deg53$n_gene_appears == 1)[1:num_x]]))
selectTop5 = function(x) {
  df = deg53[which(deg53$cluster == x & deg53$n_gene_appears == 1)[1:num_x],]
  df = df[complete.cases(df),]
  return(df)
}
topX_s3 = lapply(0:52, function(x) selectTop5(x) )
deg53_topX = data.frame(plyr::ldply(topX_s3, rbind))
topX_s4 = lapply(0:52, function(x) deg53[which(deg53$cluster == x & deg53$n_gene_appears == 1),])
deg53_unique = data.frame(plyr::ldply(topX_s4, rbind))
write.csv(deg53_unique, "C:/Users/miles/Downloads/brain/results/bb/specific_markers_53_03192021.csv")

cluster2_subs = c(0, 15, 21, 22, 27, 35, 39, 44)
deg53_2 = deg53[which(deg53$cluster %in% cluster2_subs),]
deg_table = table(deg53_2$gene)
deg53_2$n_gene_appears = deg_table[match(deg53_2$gene, names(deg_table))]
topX_s4 = lapply(cluster2_subs, function(x) deg53_2[which(deg53_2$cluster == x & deg53_2$n_gene_appears == 1),])
deg53_2_unique = data.frame(plyr::ldply(topX_s4, rbind))

clust2.pct.2 = c()
for (i in 1:nrow(deg53_2_unique)) {
  row = deg53_2_unique[i,]
  gene = as.character(row$gene)
  cluster = row$cluster
  other_cluster_cells = colnames(bb)[which(bb$seuratclusters53 %in% cluster2_subs & bb$seuratclusters53 != cluster)]
  gene_cells = colnames(bb)[which(bb@assays$RNA@counts[gene,] != 0)]
  other_cluster_gene_cells = other_cluster_cells[which(other_cluster_cells %in% gene_cells)]
  clust2.pct.2 = c(clust2.pct.2, length(other_cluster_gene_cells)/length(other_cluster_cells))
}
deg53_2_unique$clust2.pct.2 = clust2.pct.2
write.csv(deg53_2_unique, "C:/Users/miles/Downloads/brain/results/bb/specific_markers_unique_in_subclusters_from_2_53_03192021.csv")

selectTop5 = function(x) {
  deg = deg53_2_unique[order(deg53_2_unique$clust2.pct.2),]
  df = deg[which(deg$cluster == x)[1:num_x],]
  df = df[complete.cases(df$gene),]
  return(df)
}
topX_s6 = lapply(cluster2_subs, function(x) selectTop5(x) )
deg53_2_topX = data.frame(plyr::ldply(topX_s6, rbind))
write.csv(deg53_2_topX, "C:/Users/miles/Downloads/brain/results/bb/specific_markers_unique_in_subclusters_from_2_53_top5_03192021.csv")

#==========================================================================================
# Markers =================================================================================
#==========================================================================================
svg(paste0(rna_path, "/results/bb/paintings/dotplot_bhve_v_ctrl.svg"), width = 7, height = 4)
myDotPlot(bb, c("rsrp1", "smarce1", "LOC101475628"))[[1]] + coord_flip()
dev.off()

svg(paste0(rna_path, "/results/bb/paintings/rsrp1_bvc.svg"), width = 3, height = 6)
bvcVis(bb, "rsrp1", only.pos = T, meta = "cond", mode = 'box', cell.alpha = 0.02)
dev.off()
svg(paste0(rna_path, "/results/bb/paintings/smarce1_bvc.svg"), width = 3, height = 6)
bvcVis(bb, "smarce1", only.pos = T, meta = "cond", mode = 'box')
dev.off()
svg(paste0(rna_path, "/results/bb/paintings/lrch3_bvc.svg"), width = 3, height = 6)
bvcVis(bb, "lrch3", only.pos = T, meta = "cond", mode = 'box')
dev.off()

pdf(paste0(rna_path, "/results/bb/paintings/rsrp1_sample.pdf"), width = 7, height = 4)
bvcVis(bb, "rsrp1", only.pos = T, meta = "sample", mode = 'violin_split')
dev.off()
pdf(paste0(rna_path, "/results/bb/paintings/rsrp1_bvc_split.pdf"), width = 4, height = 6)
bvcVis(bb, "rsrp1", only.pos = T, meta = "cond", mode = 'violin_split')
dev.off()

bri = list()
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "Glia", skip=1, col_types = "text"))
bri[["glia"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp[,2])),"\\s")))
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "Vc, Vd", skip=6, col_types = "text"))
bri[["vc_vd"]] = gsub("[\r\n;]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "Vs-Vp, Vv-Vl", skip=5, col_types = "text"))
bri[["vs_vp_vv_vl"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "OB (gcl)", skip=5, col_types = "text"))
bri[["ob_gcl"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "OB (mcl, gl)", skip=9, col_types = "text"))
bri[["ob_mcl_gl"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "BNSM", skip=0, col_types = "text"))
bri[["bnsm"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "Dp", skip=11, col_types = "text"))
bri[["dp"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "Dm", skip=9, col_types = "text"))
bri[["dm"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "Dl", skip=9, col_types = "text"))
bri[["dl"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))

bri_clean = sapply(1:length(bri), function(x) unique(bri[[x]][which(bri[[x]] %in% rownames(bb))]))
names(bri_clean) = names(bri)
for (i in 1:length(bri_clean)) {
  markers = bri_clean[[i]]
  name = names(bri_clean)[i]
  
  # UMAP Expression
  png(file=paste0(rna_path, "/results/bb/paintings/", name, ".png"), width = 750, height = 500)
  print(markerExpPerCell(bb, markers) + ggtitle(paste0("Expression of ", str_to_title(name))))
  dev.off()
  
  # DotPlot
  dot_res = myDotPlot(bb, markers)
  png(file=paste0(rna_path, "/results/bb/paintings/", name, "_dot.png"), width = 750, height = 750)
  print(dot_res[[1]] + ggtitle(paste0("Expression of ", str_to_title(name))))
  dev.off()
}

#==========================================================================================
# BHVE vs CTRL Prediction =================================================================
#==========================================================================================
library(glmnet)
mat = bb@assays$RNA@counts
mat[which(mat > 1)] = 1
mat_rowsums = rowSums(mat)
# high_count_genes = names(mat_rowsums)[which(mat_rowsums >= 100)]
# bhve_ctrl_PA_diff = abs(rowSums(mat[,which(bb$cond == "BHVE")]) - rowSums(mat[,which(bb$cond == "CTRL")]))
# high_count_genes = names(bhve_ctrl_PA_diff)[which(bhve_ctrl_PA_diff >= 500)]
# bb = ScaleData(bb, features = high_count_genes)
Idents(bb) = bb$subsample
# subsample_exp_avgs = myAverageExpression(bb, features = high_count_genes, slot = "data")
# colnames(subsample_exp_avgs) = levels(Idents(bb))
# all_exp_avgs = myAverageExpression(bb, cells = colnames(bb)[which(bb$seuratclusters15 == 0)], slot = "data")
all_exp_avgs = as.data.frame(lapply(levels(Idents(bb)), function(x) rowSums(mat[,which(bb$subsample == x & bb$seuratclusters15 == 0)])/length(which(bb$subsample == x & bb$seuratclusters15 == 0))))
colnames(all_exp_avgs) = levels(Idents(bb))
# subsample_meta = unique(bb@meta.data[,which(!colnames(bb@meta.data) %in% c("nCount_RNA", "nFeature_RNA", "pct_mt", "RNA_snn_res.1.3", "seurat_clusters", "seuratclusters53", "RNA_snn_res.0.02", "seuratclusters15", "RNA_snn_res.0.018", "RNA_snn_res.0.015", "seuratclusters14"))])
subsample_meta = unique(bb@meta.data[,c("subsample", "pair", "pool", "gsi", "standard_length", "cond", colnames(bb@meta.data)[which(startsWith(colnames(bb@meta.data), "build"))], colnames(bb@meta.data)[which(startsWith(colnames(bb@meta.data), "depth"))], colnames(bb@meta.data)[which(startsWith(colnames(bb@meta.data), "spawn"))] )])
rownames(subsample_meta) = subsample_meta$subsample
subsample_meta$depth_adj_gsi = subsample_meta$depth/(subsample_meta$gsi)

nums_fine = c(15, 20, 25, 30, 35, 40, 50)
nums_coarse = c(10, 50, 100, 200, 500, 1000)
for (j in nums_coarse) {
  clust0Pred(j)
}

clust0Pred = function(num_genes, doScale = T, doPCA = F, this_model = "rf", bhve_metric = "depth") {
  if (! this_model %in% c("rf", "lasso", "cv")) { print("Not a valid model. Exiting."); return(NULL) }
  big_df = data.frame()
  this_cluster = 0
  # for (i in 1:5) {
  for (i in subsample_meta$subsample) {
    print(i)
    # test_samples = unique(bb$subsample)[which( substr(unique(bb$subsample), 2, 2) == i )]
    test_samples = i
    train = all_exp_avgs[,which( ! colnames(all_exp_avgs) %in% test_samples )]
    test = all_exp_avgs[,test_samples]
    
    train = as.data.frame(t(train))
    test = as.data.frame(t(test))
    colnames(test) = colnames(train)
    rownames(test) = test_samples
    
    train_bvc_mean_dif = colMeans(train[rownames(train)[which( startsWith(rownames(train), "b") )],]) - colMeans(train[rownames(train)[which( startsWith(rownames(train), "c") )],])
    # bhve_up = unlist(lapply(colnames(train), function(x) max(train[rownames(train)[which( startsWith(rownames(train), "c") )],x]) - min(train[rownames(train)[which( startsWith(rownames(train), "b") )],x]) ))
    # non_zero_i = which(sapply(colnames(train), function(x) min(train[,x]) != 0))
    # non_zero_values = bhve_up[non_zero_i]
    # non_zero_values = names(non_zero_i)[order(non_zero_values)][1:25]
    # train = train[,non_zero_values]
    # test = test[,non_zero_values]
    # bhve_up = unlist(lapply(colnames(train), function(x) min(train[rownames(train)[which( startsWith(rownames(train), "b") )],x]) > max(train[rownames(train)[which( startsWith(rownames(train), "c") )],x])))
    # ctrl_up = unlist(lapply(colnames(train), function(x) min(train[rownames(train)[which( startsWith(rownames(train), "c") )],x]) > max(train[rownames(train)[which( startsWith(rownames(train), "b") )],x])))
    high_count_genes = names(train_bvc_mean_dif[order(abs(train_bvc_mean_dif), decreasing = T)])[1:num_genes]
    train = train[,high_count_genes]
    test = test[,high_count_genes]
    # non_zero_i = which(colSums(train) != 0 & colSums(test) != 0)
    # train = train[,non_zero_i]
    # test = test[,non_zero_i]
    colnames(train) = str_replace(colnames(train), "-", "_")
    colnames(test) = str_replace(colnames(test), "-", "_")
    
    if (doPCA) {
      train.pca <- prcomp(train, center = F, scale = F)
      test.pca = predict(train.pca, test)
      train = as.data.frame(train.pca$x)
      test = as.data.frame(test.pca)
    }
    
    if(doScale) {
      print("Not Scaling bc Predicting 1 Ind.")
      # test = as.data.frame(scale(test))
      # train = as.data.frame(scale(train))
    }
  
    lambdas <- 10^seq(4, -3, by = -.1)
    if (this_model == "cv") {
      ridge_reg = glmnet(train, subsample_meta[match(rownames(train), subsample_meta$subsample), bhve_metric], nlambda = 25, alpha = 0, family = 'gaussian', lambda = lambdas)
      cv_ridge <- cv.glmnet(as.matrix(train), subsample_meta[match(rownames(train), subsample_meta$subsample), bhve_metric], alpha = 0, lambda = lambdas)
      optimal_lambda <- cv_ridge$lambda.min
      predictions_test <- predict(cv_ridge, s = optimal_lambda, newx = as.matrix(test))
    } 
    if (this_model == "lasso") {
      lasso_reg <- cv.glmnet(as.matrix(train), subsample_meta[match(rownames(train), subsample_meta$subsample), bhve_metric], alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
      lambda_best <- lasso_reg$lambda.min
      lasso_model <- glmnet(as.matrix(train), subsample_meta[match(rownames(train), subsample_meta$subsample), bhve_metric], alpha = 1, lambda = lambda_best, standardize = TRUE)
      predictions_test <- predict(lasso_model, s = lambda_best, newx = as.matrix(test))  
    } 
    if (this_model == "rf") {
      train[, bhve_metric] = subsample_meta[match(rownames(train), subsample_meta$subsample), bhve_metric]
      test[, bhve_metric] = subsample_meta[match(rownames(test), subsample_meta$subsample), bhve_metric]
      myFormula = as.formula(paste(bhve_metric, "~ ."))
      m <- caret::train(myFormula, data = train, method = "ranger")
      predictions_test <- predict(m, test)
    }
    big_df = rbind(big_df, data.frame(subsample = rownames(test), pred = unname(predictions_test), real = subsample_meta[match(rownames(test), subsample_meta$subsample), bhve_metric]))
  } # end pair for
  big_df$cond = startsWith(big_df$subsample, "b")
  print(ggplot(big_df, aes(pred, real, color = abs(cond))) + geom_point() + ggtitle("Test") + geom_text_repel(aes(label = subsample)) + ggtitle(paste("Number of Genes:", num_genes)))
  print(ggplot(big_df, aes(x = cond, y = pred, color = abs(cond))) + geom_boxplot() + geom_point(stat = "identity", position = position_jitterdodge()) + ggtitle(paste("Number of Genes:", num_genes)))
}


input_thresh = c(400, 450, 500, 550, 600, 650, 700, 750, 800)
input_thresh = seq(700, 800, by = 5)
all_runs_df = as.data.frame(mclapply(input_thresh, singleRun, mc.cores = 24))
colnames(all_runs_df) = input_thresh
rownames(all_runs_df) = c("mean_b", "mean_c", "sd_b", "sd_c", "min_b", "min_c", "max_b", "max_c")
all_runs_df
test = mclapply(input_thresh, singleRun, mc.cores = 24)
all_runs_df = ldply(test, data.frame)
png(paste0("~/scratch/brain/results/bvc_pred_demux_scale_std_depth_small2.png"), width = 2000, height = 1600, res = 120)
print(ggplot(all_runs_df, aes(x = cond, y = pred, color = abs(cond))) + geom_boxplot() + geom_point(stat = "identity", position = position_jitterdodge()) + ggtitle("Test") + facet_wrap(unique(all_runs_df$thresh)))
dev.off()

input_thresh = seq(5, 100, by = 5)
test = mclapply(input_thresh, singleRun, mc.cores = 24)
all_runs_df = as.data.frame(mclapply(input_thresh, singleRun, mc.cores = 24))

other_method_res = singleRun(760, bhve_metric = "depth")
png(paste0("~/scratch/brain/results/bvc_pred_demux_test.png"), width = 800, height = 600, res = 120)
print(ggplot(other_method_res, aes(x = cond, y = pred, color = abs(cond))) + geom_boxplot() + geom_point(stat = "identity", position = position_jitterdodge()) + ggtitle("Test"))
dev.off()
png(paste0("~/scratch/brain/results/bvc_pred_demux_test1.png"), width = 1000, height = 800, res = 120)
print(ggplot(other_method_res, aes(x = real, y = pred, color = abs(cond))) + geom_point() + ggtitle("Test")) + geom_text_repel(aes(label = subsample))
dev.off()

test = optimize(singleRun, interval = c(5, 100))

singleRun = function(thresh, doScale = T, doPlot = F, bhve_metric = "std_depth") {
  # high_count_genes = names(bhve_ctrl_PA_diff)[which(bhve_ctrl_PA_diff >= thresh)]
  Idents(bb) = bb$subsample
  # subsample_exp_avgs = as.data.frame(lapply(levels(Idents(bb)), function(x) rowSums(mat[high_count_genes,which(bb$subsample == x)])/length(which(bb$subsample == x))))
  # colnames(subsample_exp_avgs) = levels(Idents(bb))
  # all_exp_avgs = as.data.frame(lapply(levels(Idents(bb)), function(x) rowSums(mat[,which(bb$subsample == x)])/length(which(bb$subsample == x))))
  all_exp_avgs = myAverageExpression(bb, slot = "data")
  colnames(all_exp_avgs) = levels(Idents(bb))
  
  all_pair_df = data.frame()
  for (i in 1:5) {
    # print("Finding PA diff")
    test_samples = unique(bb$subsample)[which( substr(unique(bb$subsample), 2, 2) == i )]
    this_bhve_ctrl_PA_diff = abs(rowSums(mat[,which(bb$cond == "BHVE" & ! bb$subsample %in% test_samples)]) - rowSums(mat[,which(bb$cond == "CTRL" & ! bb$subsample %in% test_samples)]))
    # high_count_genes = names(this_bhve_ctrl_PA_diff)[which(this_bhve_ctrl_PA_diff >= thresh)]
    this_bhve_ctrl_PA_diff = this_bhve_ctrl_PA_diff[order(this_bhve_ctrl_PA_diff, decreasing = T)[1:thresh]]
    high_count_genes = names(this_bhve_ctrl_PA_diff)
    subsample_exp_avgs = all_exp_avgs[high_count_genes,]
    # print(length(high_count_genes))
    # print("Done")
    train = subsample_exp_avgs[,which( ! colnames(subsample_exp_avgs) %in% test_samples )]
    test = subsample_exp_avgs[,test_samples]
    train = as.data.frame(t(train))
    test = as.data.frame(t(test))
    colnames(test) = colnames(train)
    rownames(test) = test_samples
    
    if(doScale) {
      test = as.data.frame(scale(test))
      train = as.data.frame(scale(train))
    }
    train[, bhve_metric] = subsample_meta[match(rownames(train), subsample_meta$subsample), bhve_metric]
    test[, bhve_metric] = subsample_meta[match(rownames(test), subsample_meta$subsample), bhve_metric]
    
    # lambdas <- 10^seq(4, -3, by = -.1)
    # ridge_reg = glmnet(train, subsample_meta[match(rownames(train), subsample_meta$subsample), bhve_metric], nlambda = 25, alpha = 0, family = 'gaussian', lambda = lambdas)
    # cv_ridge <- cv.glmnet(as.matrix(train), subsample_meta[match(rownames(train), subsample_meta$subsample), bhve_metric], alpha = 0, lambda = lambdas)
    # optimal_lambda <- cv_ridge$lambda.min
    # predictions_test <- predict(ridge_reg, s = optimal_lambda, newx = as.matrix(test))
    # library("rpart")
    myFormula = as.formula(paste(bhve_metric, "~ ."))
    m <- caret::train(myFormula, data = train, method = "ranger")
    predictions_test <- predict(m, test)
    
    # print(eval_results(subsample_meta[match(rownames(test), subsample_meta$subsample), "depth"], predictions_test, test))
    all_pair_df = rbind(all_pair_df, data.frame(subsample = rownames(test), pred = unname(predictions_test), real = subsample_meta[match(rownames(test), subsample_meta$subsample), bhve_metric]))
  }
  all_pair_df$cond = startsWith(all_pair_df$subsample, "b")
  all_pair_df$thresh = thresh
  
  if (doPlot) {
    print("plotting")
    png(paste0("~/scratch/brain/results/bvc_pred_demux_scale", thresh, ".png"), width = 1000, height = 800, res = 120)
    print(ggplot(all_pair_df, aes(pred, real, color = abs(pred - real))) + geom_point() + ggtitle("Test") + geom_text_repel(aes(label = subsample)))
    dev.off()
  }
  
  mean_b = mean(all_pair_df$pred[which(all_pair_df$cond)])
  mean_c = mean(all_pair_df$pred[which(all_pair_df$cond == F)])
  sd_b = sd(all_pair_df$pred[which(all_pair_df$cond)])
  sd_c = sd(all_pair_df$pred[which(all_pair_df$cond == F)])
  min_b = min(all_pair_df$pred[which(all_pair_df$cond)])
  min_c = min(all_pair_df$pred[which(all_pair_df$cond == F)])
  max_b = max(all_pair_df$pred[which(all_pair_df$cond)])
  max_c = max(all_pair_df$pred[which(all_pair_df$cond == F)])
  # rss = sum((all_pair_df$pred - all_pair_df$real)^2)
  # return(rss)
  return(all_pair_df)
  # return(c(mean_b, mean_c, sd_b, sd_c, min_b, min_c, max_b, max_c))
}

rPartMod <- caret::train(depth ~ ., data=train, method="RRF")
rpartImp <- varImp(rPartMod)

test_samples = c("b1.1", "b1.2", "b1.3", "b1.4", "c1.1", "c1.2", "c1.3", "c1.4")
Idents(bb) = bb$subsample
subsample_exp_avgs = myAverageExpression(bb, features = high_count_genes, slot = "data")
train = subsample_exp_avgs[,which( ! colnames(subsample_exp_avgs) %in% test_samples )]
test = subsample_exp_avgs[,test_samples]

train = as.data.frame(t(train))
test = as.data.frame(t(test))
colnames(test) = colnames(train)
rownames(test) = test_samples
train[, "build_events"] = subsample_meta[match(rownames(train), subsample_meta$subsample), "build_events"]
test[, "build_events"] = subsample_meta[match(rownames(test), subsample_meta$subsample), "build_events"]
# train[, colnames(subsample_meta)[2:ncol(subsample_meta)]] = subsample_meta[match(rownames(train), subsample_meta$subsample), colnames(subsample_meta)[2:ncol(subsample_meta)]]
# test[, colnames(subsample_meta)[2:ncol(subsample_meta)]] = subsample_meta[match(rownames(test), subsample_meta$subsample), colnames(subsample_meta)[2:ncol(subsample_meta)]]
my_model = glm(build_events ~ ., data = train)
predict(my_model, newdata=test, type="response")

lambdas <- 10^seq(4, -3, by = -.1)
ridge_reg = glmnet(train, subsample_meta[match(rownames(train), subsample_meta$subsample), "build_events"], nlambda = 25, alpha = 0, family = 'gaussian', lambda = lambdas)
cv_ridge <- cv.glmnet(as.matrix(train), subsample_meta[match(rownames(train), subsample_meta$subsample), "build_events"], alpha = 0, lambda = lambdas)
optimal_lambda <- cv_ridge$lambda.min
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  
  
  # Model performance metrics
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
  
}
# Prediction and evaluation on train data
predictions_train <- predict(ridge_reg, s = optimal_lambda, newx = as.matrix(train))
eval_results(subsample_meta[match(rownames(train), subsample_meta$subsample), "build_events"], predictions_train, train)
# Prediction and evaluation on test data
predictions_test <- predict(ridge_reg, s = optimal_lambda, newx = as.matrix(test))
eval_results(subsample_meta[match(rownames(test), subsample_meta$subsample), "build_events"], predictions_test, test)

p_df = data.frame(pred = unname(predictions_train), real = subsample_meta[match(rownames(train), subsample_meta$subsample), "build_events"], subsample = rownames(train))
ggplot(p_df, aes(pred, real, color = abs(pred - real))) + geom_point() + ggtitle("Train") + geom_text_repel(aes(label = subsample))
p_df = data.frame(pred = unname(predictions_test), real = subsample_meta[match(rownames(test), subsample_meta$subsample), "build_events"], subsample = rownames(test))
ggplot(p_df, aes(pred, real, color = abs(pred - real))) + geom_point() + ggtitle("Test") + geom_text_repel(aes(label = subsample))

# cv.lasso <- cv.glmnet(as.matrix(train), subsample_meta[match(rownames(train), subsample_meta$subsample), "build_events"], alpha=1, parallel=TRUE, standardize=TRUE, type.measure='auc')
# plot(cv.lasso)
lasso_reg <- cv.glmnet(as.matrix(train), subsample_meta[match(rownames(train), subsample_meta$subsample), "build_events"], alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
lambda_best <- lasso_reg$lambda.min 
lasso_model <- glmnet(as.matrix(train), subsample_meta[match(rownames(train), subsample_meta$subsample), "build_events"], alpha = 1, lambda = lambda_best, standardize = TRUE)
predictions_train <- predict(lasso_model, s = lambda_best, newx = as.matrix(train))
eval_results(subsample_meta[match(rownames(train), subsample_meta$subsample), "build_events"], predictions_train, train)
predictions_test <- predict(lasso_model, s = lambda_best, newx = as.matrix(test))
eval_results(subsample_meta[match(rownames(test), subsample_meta$subsample), "build_events"], predictions_test, test)



train = t(bb@assays$RNA@data[high_count_genes, which(bb$subsample != "b1.1")])
test = t(bb@assays$RNA@data[high_count_genes, which(bb$subsample == "b1.1")])
# train_sparse = sparse.model.matrix(~.,train)
fit = glmnet(train, bb$build_events[which(bb$subsample != "b1.1")])
cv <- cv.glmnet(train,bb$build_events[which(bb$subsample != "b1.1")],nfolds=3)
pred_train <- predict(fit, train, type="response",s=cv$lambda.min)
pred_test <- predict(fit, test, type="response",s=cv$lambda.min)
mean(abs(pred_train - bb$build_events[which(bb$subsample != "b1.1")]))
mean(abs(pred_test - bb$build_events[which(bb$subsample == "b1.1")]))

p_df = data.frame(pred = pred_train, real = bb$build_events[which(bb$subsample != "b1.1")])
ggplot(p_df, aes(pred, real, color = abs(pred - real))) + geom_point() + ggtitle("Train")

bhve_samples = c("b1", "b2", "b3", "b4", "b5")
ctrl_samples = c("c1", "c2", "c3", "c4", "c5")

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

# Scrap Code for Python/PACE
bb$sample.clust = paste0(bb$sample, ".", bb$seuratclusters53)
Idents(bb) = bb$sample.clust
avgs = myAverageExpression(bb)

small_mat = as.matrix(avgs)
# small_mat = avgs[which(rownames(avgs) %in% zack2$mzebra),]
# small_mat = small_mat[match(zack2$mzebra[which(! duplicated(zack2$mzebra) )], rownames(small_mat)),]
new_avg = list()
samples = unique(as.vector(bb$sample))
clusters = sort(unique(as.numeric(as.vector(bb$seuratclusters53))))
for (sample in samples) {
  for (cluster in clusters) {
    n_cells = length(which(bb$sample.clust == paste0(sample, ".", cluster)))
    if (n_cells == 0) {
      new_avg[[sample]] = c(new_avg[[sample]], rep(0, nrow(small_mat)))
    } else {
      new_avg[[sample]] = c(new_avg[[sample]], small_mat[,which(colnames(small_mat) == paste0(sample, ".", cluster))])
    }
  }
}
test = t(as.data.frame(new_avg))
colnames(test) = paste0( rep(clusters, each=nrow(small_mat)) , ".", rownames(small_mat))

# PCA
library("factoextra")
test = read.table("C:/Users/miles/Downloads/c_hit_cluster_order_data_small.txt")
test2 = test[,1:8]
colnames(test2) = unlist(strsplit(colnames(test2), ".", fixed = T))[c(FALSE, TRUE)]
colnames(test2) = sapply(colnames(test2), function(x) paste0(x, " (", gene_info$human[which(gene_info$mzebra == x)], ")") )
res.pca <- prcomp(test2, scale = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, col.ind = "cos2", repel = TRUE)
fviz_pca_ind(res.pca, col.ind = "x", gradient.cols = brewer.pal(n = 11, name = "RdYlBu"), repel = TRUE)
fviz_pca_var(res.pca, col.var = "contrib", gradient.cols = rev(grey.colors(6))[3:6], repel = TRUE)
fviz_pca_var(res.pca, col.var = "x", gradient.cols = brewer.pal(n = 11, name = "RdYlBu"), repel = TRUE)
fviz_pca_biplot(res.pca, repel = TRUE, col.var = "#2E9FDF", col.ind = "#696969")
fviz_pca_biplot(res.pca, repel = TRUE, col.var = "x", col.ind = "x", gradient.cols = brewer.pal(n = 11, name = "RdYlBu"))

groups = factor(c(rep("Behave", 5), rep("Control", 5)))
fviz_pca_ind(res.pca, col.ind = groups, palette = c("#00AFBB",  "#FC4E07"), addEllipses = T, ellipse.type = "confidence", legend.title = "Groups", repel = T)
fviz_pca_var(res.pca, col.var = "contrib", gradient.cols = rev(grey.colors(6))[3:6], repel = TRUE)

# Plot value of coefficients
coef = read.csv("C:/Users/miles/Downloads/coef_cluster.txt")
coef_df = data.frame()
for (i in 1:nrow(coef)) {
  coef_df = rbind(coef_df, unname(cbind(rownames(coef)[i], t(coef[i,]))))
}
coef_df$V2 = as.numeric(as.vector(coef_df$V2))
coef_df$Predicts = coef_df$V2 > 0
coef_df$Predicts = plyr::revalue(as.character(coef_df$Predicts), replace = c("TRUE" = "Behave", "FALSE" = "Control"))
ggplot(coef_df, aes(x = V1, y = V2, fill = Predicts, color = Predicts)) + geom_bar(stat = "identity", alpha = 0.6) + xlab("Gene") + ylab("Predictive Power")

# PCA and subset by a set of genes
test = scan("C:/Users/miles/Downloads/sample_clust53_avg.txt", what = "character")
if (bulk)
  test[14074] = "my_var%"
test = test[which( ! test %in% c(bhve_samples, ctrl_samples) )]
test[which(test == "my_var%")] = "c5"
test = t(matrix(test, ncol = 11, byrow = F))
colnames(test) = test[1,]
test = data.matrix(test[-1,])
rownames(test) = c(bhve_samples, ctrl_samples)
class(test) = "numeric"
coef = read.csv("C:/Users/miles/Downloads/coef_53_4_172021.txt", stringsAsFactors = F)
all_combos = c(coef$X0, coef$X1, coef$X2, coef$X3, coef$X4)
test2 = test[, unique(all_combos)]
res.pca <- prcomp(test2)
groups = factor(c(rep("Behave", 5), rep("Control", 5)))
fviz_pca_ind(res.pca, col.ind = groups, palette = c("#00AFBB",  "#FC4E07"), addEllipses = T, ellipse.type = "confidence", legend.title = "Groups", repel = T)
genes = unlist(strsplit(all_combos, ".", fixed = T))[c(FALSE, TRUE)]
genes_df = as.data.frame(table(genes))
clust = unlist(strsplit(all_combos, ".", fixed = T))[c(TRUE, FALSE)]
clust_df = as.data.frame(table(clust))
clust_df$clust = factor(clust_df$clust, levels = 0:14)
ggplot(clust_df, aes(x=clust, y=Freq)) + geom_bar(stat="identity") + ggtitle("Overlap of Clusters")
ggplot(genes_df, aes(x=Freq)) + geom_bar() + ggtitle("Overlap of Genes")
ggplot(combo_df, aes(x=Freq)) + geom_bar() + ggtitle("Overlap of Combos")

# SubSample 
coef_df = read.csv("C:/Users/miles/Downloads/ml_subsample_counts_bulk.txt")
coef_df$X = NULL
Idents(bb) = bb$subsample
avgs = myAverageExpression(bb, slot = "counts")
small_mat = as.matrix(avgs)
test2 = t(small_mat[which(rownames(small_mat) %in% as.vector(t(coef_df)))])
groups = factor(c(rep("Behave", 19), rep("Control", 19)))
pdf("C:/Users/miles/Downloads/ml_subsample_counts_bulk.pdf", width = 8, height = 8)
print(fviz_pca_ind(res.pca, col.ind = groups, palette = c("#00AFBB",  "#FC4E07"), addEllipses = T, ellipse.type = "confidence", legend.title = "Groups", repel = T))
dev.off()
genes_df = as.data.frame(table( as.vector(t(coef_df)) ))
ggplot(genes_df, aes(x=Freq)) + geom_bar() + ggtitle("Overlap of Genes in Models") + xlab("Number of Pairs a Gene is Used In") + ylab("Number of Genes")

# Biased Approach
Idents(bb) = bb$seuratclusters53
t_test = t(test)
all_same = data.frame()
for (i in 1:nrow(t_test)) {
  if ( ! any(t_test[i,c(1,2,3,4,5)] == 0) & ! any(t_test[i,c(6,7,8,9,10)] == 0) ) {
    b_max = max(t_test[i,c(1,2,3,4,5)])
    b_min = min(t_test[i,c(1,2,3,4,5)])
    c_max = max(t_test[i,c(6,7,8,9,10)])
    c_min = min(t_test[i,c(6,7,8,9,10)])
    
    this_gene = rownames(t_test)[i]
    if (b_min > c_max) {
      # num_cells = length(which(bb@assays$RNA@counts[this_gene,] > 0))
      this_cells = WhichCells(bb, idents = str_split(this_gene, "\\.")[[1]][1])
      num_cells = length(which(bb@assays$RNA@counts[str_split(this_gene, "\\.")[[1]][2], this_cells] > 0))
      tot_dif = sum(t_test[i,c(1,2,3,4,5)]) - sum(t_test[i,c(6,7,8,9,10)])
      max_dif = b_min - c_max
      all_same = rbind(all_same, t(c( this_gene, "BHVE", num_cells, tot_dif, max_dif )))
    }
    else if (c_min > b_max) {
      # num_cells = length(which(bb@assays$RNA@counts[this_gene,] > 0))
      this_cells = WhichCells(bb, idents = str_split(this_gene, "\\.")[[1]][1])
      num_cells = length(which(bb@assays$RNA@counts[str_split(this_gene, "\\.")[[1]][2], this_cells] > 0))
      tot_dif = sum(t_test[i,c(6,7,8,9,10)]) - sum(t_test[i,c(1,2,3,4,5)])
      max_dif = c_min - b_max
      all_same = rbind(all_same, t(c( this_gene, "CTRL", num_cells, tot_dif, max_dif )))
    }
  } # non-zero if
}
colnames(all_same) = c("feature", "up_cond", "num_cells", "tot_dif", "max_dif")
all_same$cluster = unlist(strsplit(as.vector(all_same$feature), ".", fixed = T))[c(TRUE, FALSE)]
all_same$gene = unlist(strsplit(as.vector(all_same$feature), ".", fixed = T))[c(FALSE, TRUE)]
all_same$num_cells = as.numeric(as.vector(all_same$num_cells))
all_same$tot_dif = as.numeric(as.vector(all_same$tot_dif))
all_same$max_dif = as.numeric(as.vector(all_same$max_dif))
all_same$hgnc = gene_info$human[match(all_same$gene, gene_info$mzebra)]
all_same = cbind(all_same, t_test[as.vector(all_same$feature),])
write.table(all_same$feature[order(all_same$max_dif, decreasing = T)], "C:/Users/miles/Downloads/53_biased_genes.txt", quote = F, row.names = F, col.names = F)
write.csv(all_same, "C:/Users/miles/Downloads/53_biased_data_table.csv")

convert15$new.full = as.vector(convert15$new.full)
convert15$new.full = factor(convert15$new.full, levels = convert15$new.full[order(as.numeric(convert15$new.num))])
ggplot(convert15, aes(x=new.full, y = res_pct, color = col, fill = col)) + geom_bar(stat = "identity") + scale_color_identity() + scale_fill_identity()

convert15_2 = data.frame()
for (i in 0:14) {
  n = convert15$res[which(convert15$old == i)]
  # if (is.na(n)) { newRows = data.frame(old = i, new = convert15$new.full[which(convert15$old == i)], col = convert15$col[which(convert15$old == i)], num = 1:n) }
  for ( j in 1:38 ) {
    this_col = convert15$col[which(convert15$old == i)]
    if ( is.na(n) ) { this_col = "lightgray"} else { 
      if ( j > n ) { this_col = "gray" }
    }
    convert15_2 = rbind(convert15_2, data.frame(old = i, new = convert15$new.full[which(convert15$old == i)], col = this_col, num = j))
  }
  # newRows = data.frame(old = i, new = convert15$new.full[which(convert15$old == i)], col = convert15$col[which(convert15$old == i)], num = 1:n)
  # convert15_2 = rbind(convert15_2, newRows)
}
pdf("C:/Users/miles/Downloads/clust15_ml_061021.pdf", width = 12, height = 5)
ggplot(convert15_2, aes(x = new, y = num, color = col, fill = col)) + geom_tile(width = 0.6, linetype = "solid") + scale_color_identity() + scale_fill_identity() + theme_classic() + scale_y_continuous(expand = c(0, 0)) + xlab("") + ylab("Number of Correctly Predicted Individuals")
dev.off()

ggplot(convert15_2, aes(x = new, y = num, color = col, fill = col)) + geom_point(size = 3) + scale_color_identity() + scale_fill_identity() + theme_classic()
ggplot(convert15_2, aes(xmin = old, xmax = old+1, ymin = num, ymax = num+0.5, color = col, fill = col)) + geom_rect(size = 2) + scale_color_identity() + scale_fill_identity() + theme_classic()

#=======================================================================================
# BHVE vs CTRL Correlation =============================================================
#=======================================================================================
library("parallel")
numCores <- detectCores()
bhve_binary = bb$cond
bhve_binary = as.numeric(as.vector( plyr::revalue(as.character(bhve_binary), replace = c("BHVE" = 1, "CTRL" = 0)) ))
mat2 = rbind(bhve_binary, bhve_binary)
bhve_mat = cor(t(as.matrix(bb@assays$RNA@data[,])), y = t(mat2)) 
gene_cells = read.csv("~/scratch/brain/results/cells_b_v_c.csv")
gene_cells$X = NULL
bhve_mat = data.frame(bhve_mat)
bhve_mat$bhve_binary.1 = NULL
test = merge(bhve_mat, gene_cells)
test$above = abs(test$bhve_binary) > quantile(abs(bhve_mat$bhve_binary), 0.99, na.rm = T)
test$sig = test$q < 0.05

png(file = "~/scratch/brain/results/bhve_cor.png", width = 1000, height = 1000, res = 100)
ggplot(test, aes(pct_dif, bhve_binary, color = sig)) + geom_point() + ylab("Correlation With Behve") + xlab("% Difference in Number of Cells")
dev.off()
system(paste0("rclone copy ~/scratch/brain/results/bhve_cor.png dropbox:BioSci-Streelman/George/Brain/bb/results/coexp/"))


#==========================================================================================
# BHVE vs CTRL Bootstrap ==================================================================
#==========================================================================================
# Purpose: Determine if there's a real BHVE vs CTRL signal by
#          comparing the real BHVE vs CTRL to simulated ones.

# Load object called obj
obj = bb

# Constants
n_boot = 50 # number of bootstraps
do_replace = F # sample with replacement?
bhve_samples = c("b1", "b2", "b3", "b4", "b5")
ctrl_samples = c("c1", "c2", "c3", "c4", "c5")
clusters = sort(unique(obj$seurat_clusters))
all_samples = c(bhve_samples, ctrl_samples)
bhve_cells = colnames(obj)[which(obj$sample %in% bhve_samples)]
ctrl_cells = colnames(obj)[which(obj$sample %in% ctrl_samples)]
obj$seurat_clusters_vect = as.vector(obj$seurat_clusters)
names(obj$seurat_clusters_vect) = colnames(obj)

# Generate random samples
print("Generating random samples and Finding DEGs")
sim_degs = data.frame()
for (i in 1:n_boot) {
  set.seed(i+100)
  # Some Useful Print Statements
  if(i == n_boot) {
    cat(paste(i, "\n"))
  } else if (i %% (n_boot/10) == 0 || i == 1) {
    cat(i)
  } else {
    cat(".")
  }

  # Shuffle Clusters while keeping the # cells per sample per cluster the same
  obj$sim = obj$seurat_clusters_vect
  Idents(obj) = obj$seurat_clusters
  for (cluster in clusters) {
    this_cells = WhichCells(obj, idents = cluster, seed = NULL)
    sim = sample(obj$cond[this_cells], length(this_cells))
    obj$sim[this_cells] = unname(sim)
  }
  
  # Find DEGs
  obj$sim_cond = paste0(obj$seurat_clusters, obj$sim)
  Idents(obj) = obj$sim_cond
  for (cluster in clusters) {
    degs = FindMarkers(obj, ident.1 = paste0(cluster, "BHVE"), ident.2 = paste0(cluster, "CTRL"), only.pos=F, verbose = F, min.pct = 0.01, logfc.threshold = 0.1)
    if ( ! is.null(degs) && nrow(degs) > 0 ) {
      degs$boot = i
      degs$cluster = cluster
      sim_degs = rbind(sim_degs, degs)
    }
  }
} # end boot for
# sim_degs = sim_degs[which(sim_degs$boot <= 50),]
write.table(sim_degs, paste0(rna_path, "results/bb/sim_bhve_v_ctrl.tsv"), sep = "\t", quote = F, row.names = F)

perm_df2 = data.frame()
for (i in 1:n_boot) {
  degs = sim_degs[which(sim_degs$boot == i),]
  qobj = qvalue(degs$p_val)
  degs$q = qobj$qvalues
  degs$bh = p.adjust(degs$p_val, method = "BH")
  
  n_p = length(which(degs$p_val < 0.05))
  n_q = length(which(degs$q < 0.05))
  n_bh = length(which(degs$bh < 0.05))
  n_bon = length(which(degs$p_val_adj < 0.05))
  
  perm_df2 = rbind(perm_df2, t(c(i, n_p, n_q, n_bh, n_bon)))
}

real_value_p = length(which(zack$p_val < 0.05))
real_value_q = length(which(zack$q < 0.05))
real_value_bh = length(which(zack$bh < 0.05))
real_value_bon = length(which(zack$p_val_adj < 0.05))
colnames(perm_df2) = c("Boot", "n_p", "n_q", "n_bh", "n_bon")

perm_df2$above = perm_df2$n_bon > real_value_bon
ggplot(perm_df2, aes(n_bon, alpha=.7, fill=above)) + geom_histogram(alpha=0.5) + geom_vline(aes(xintercept = real_value_bon)) + geom_text(aes(x=real_value_bon, label="Real Value"), y = Inf, hjust=0, vjust=1, color = "black") + xlab("# of Gene w/ adjusted p-value < 0.05 (Bonferroni)") + ggtitle("Comparison Between Bootstrap Values and Real Value") + guides(color=F, alpha=F, fill=F)

perm_df2$above = perm_df2$n_bh > real_value_bh
ggplot(perm_df2, aes(n_bh, alpha=.7, fill=above)) + geom_histogram(alpha=0.5) + geom_vline(aes(xintercept = real_value_bh)) + geom_text(aes(x=real_value_bh, label="Real Value"), y = Inf, hjust=0, vjust=1, color = "black") + xlab("# of Gene w/ adjusted FDR adjusted pvalue < 0.05 (BH)") + ggtitle("Comparison Between Bootstrap Values and Real Value") + guides(color=F, alpha=F, fill=F)

perm_df2$above = perm_df2$n_q > real_value_q
ggplot(perm_df2, aes(n_q, alpha=.7, fill=above)) + geom_histogram(alpha=0.5) + geom_vline(aes(xintercept = real_value_q)) + geom_text(aes(x=real_value_q, label="Real Value"), y = Inf, hjust=0, vjust=1, color = "black") + xlab("# of Gene w/ qvalue < 0.05") + ggtitle("Comparison Between Bootstrap Values and Real Value") + guides(color=F, alpha=F, fill=F)


perm_df2$above = perm_df2$n_p > real_value_p
ggplot(perm_df2, aes(n_p, alpha=.7, fill=above)) + geom_histogram(alpha=0.5) + geom_vline(aes(xintercept = real_value_p)) + geom_text(aes(x=real_value_p, label="Real Value"), y = Inf, hjust=0, vjust=1, color = "black") + xlab("# of Gene w/ pvalue < 0.05") + ggtitle("Comparison Between Bootstrap Values and Real Value") + guides(color=F, alpha=F, fill=F)


# Test Cluster DEG vs Bulk DEG for bhve v ctrl
obj$seurat_clusters_cond = paste0(obj$seurat_clusters, "_", obj$cond)
Idents(obj) = obj$seurat_clusters_cond
# bhve_cells = names(obj$cond[which(obj$cond == "BHVE")])
cluster_deg = data.frame()
for (cluster in sort(unique(obj$seurat_clusters))) {
  print(cluster)
  tryCatch ({
    this_deg = FindMarkers(obj, paste0(cluster, "_BHVE"), paste0(cluster, "_CTRL"), min.pct = 0.01, logfc.threshold = 0.1)
    this_deg$cluster = cluster
    this_deg$gene = rownames(this_deg)
    # this_deg = this_deg[which(this_deg$p_val_adj < 0.05),]
    cluster_deg = rbind(cluster_deg, this_deg)
  }, error = function(err) {
    print("Failed on this cluster")
  })


}
nrow(cluster_deg)
length(unique(rownames(cluster_deg)))

Idents(obj) = obj$cond
bulk_deg = FindAllMarkers(obj, only.pos = F)
# bulk_deg = bulk_deg[which(bulk_deg$p_val_adj < 0.05),]
nrow(bulk_deg)/2
length(unique(bulk_deg$gene))

cluster_deg$bonferroni = p.adjust(cluster_deg$p_val, method="bonferroni")
cluster_deg$bh = p.adjust(cluster_deg$p_val, method="BH")
qobj = qvalue(cluster_deg$p_val)
cluster_deg$qvalue = qobj$qvalues
nrow(cluster_deg[which(cluster_deg$q < 0.05),])

all_coef = read.csv("C:/Users/miles/Downloads/scgnn_bin_pred_coef.csv")
all_coef_melt = melt(all_coef)
all_coef$mean = rowMeans(all_coef[2:ncol(all_coef)])
all_coef$stderr = sapply(1:nrow(all_coef), function(x) sd(all_coef[x, 2:ncol(all_coef)]/sqrt(length(ncol(all_coef)-1))) )
# all_coef_p = sapply(all_coef$X, function(x) t.test(all_coef_melt$value[which(all_coef_melt$X == x)], all_coef_melt$value[which(all_coef_melt$X != x)])$p.value )
all_coef_p = 2*pnorm(-abs(scale(all_coef$mean)))
all_coef_bh = p.adjust(all_coef_p, method = "BH")
names(all_coef_p) = names(all_coef_bh) = all_coef$X
sort(all_coef_p[which(all_coef_p < 0.05)])
sort(all_coef_bh[which(all_coef_bh < 0.05)])
head(all_coef[order(all_coef$mean, decreasing = T), c("X", "mean")])
all_coef[match(names(sort(all_coef_bh[which(all_coef_bh < 0.05)])), all_coef$X), c("X", "mean", "stderr")]
hist(all_coef$mean, breaks = 50)

all_coef_melt = all_coef_melt[which(all_coef_melt$X %in% names(all_coef_bh[which(all_coef_bh < 0.05)])),]
all_coef_melt$hgnc = gene_info$human[match(all_coef_melt$X, gene_info$mzebra)]
all_coef_melt$hgnc[which(is.na(all_coef_melt$hgnc))] = as.vector(all_coef_melt$X[which(is.na(all_coef_melt$hgnc))])
all_coef_melt$hgnc[which(all_coef_melt$hgnc != "LOC101476487")] = tolower(all_coef_melt$hgnc[which(all_coef_melt$hgnc != "LOC101476487")])
all_coef_melt$X = factor(all_coef_melt$X, levels = all_coef$X[order(abs(all_coef$mean), decreasing = T)])
all_coef_melt$hgnc = factor(all_coef_melt$hgnc, levels = unique(all_coef_melt$hgnc[order(all_coef_melt$X)]))
all_coef_melt$abs_value = abs(all_coef_melt$value)
# ggplot(all_coef_melt[which(all_coef_melt$X %in% all_coef$X[order(abs(all_coef$mean), decreasing = T)[1:50]]),], aes(x = X, y = abs_value, color = X, fill = X)) + geom_boxplot(alpha = 0.7) + geom_point(alpha = 0.2, position = position_jitter()) + theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("Importance to Model") + scale_y_continuous(expand = c(0,0))
pdf("~/research/brain/results/ml_genes.pdf", width = 6, height = 3)
ggplot(all_coef_melt[which(all_coef_melt$X %in% all_coef$X[order(abs(all_coef$mean), decreasing = T)[1:50]]),], aes(x = hgnc, y = abs_value, color = hgnc, fill = hgnc)) + geom_boxplot(alpha = 0.7, outlier.shape = NA) + geom_point(alpha = 0.5, size = 0.5, position = position_jitter()) + theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic")) + xlab("") + ylab("Importance to Model") + scale_y_continuous(expand = c(0,0))
dev.off()

#*****************************************************************************************
# Neurogenesis ===========================================================================
#*****************************************************************************************
ieg = read.csv("~/research/brain/data/ieg_like_fos_egr1_npas4_detected_011521.csv")[,1]
neurogen = read.csv("~/research/brain/data/conserved_neurogenesis_positive_zfish_mouse_cichlid.csv")[,4]
fst = read.csv("~/research/brain/data/pcrc_FST20_30_LG11_evolution_genes_031821.csv")[,1]
non_zero_genes = rownames(bb)[which(rowSums(bb@assays$RNA@counts) > 0)]

mat = bb@assays$RNA@counts
mat[which(mat > 1)] = 1

bb$ieg = colSums(mat[ieg,])
bb$neurogen = colSums(mat[neurogen,])
bb$fst = colSums(mat[fst,])
real_df = data.frame(ieg = bb$ieg, neurogen = bb$neurogen, fst = bb$fst, subsample = bb$subsample)
real_df_agr = aggregate(. ~ subsample, data = real_df, mean)
cor(real_df_agr$ieg, real_df_agr$neurogen)
cor(real_df_agr$fst, real_df_agr$neurogen)
cor(bb$ieg, bb$neurogen)
cor(bb$fst, bb$neurogen)

real_df$cluster15 = bb$seuratclusters15
real_df$cluster53 = bb$seuratclusters53
real_df_agr = aggregate(. ~ subsample + cluster15, data = real_df, mean)
real_clust15_res = data.frame()
for (i in 0:14) {
  real_clust15_res = rbind(real_clust15_res, t(c(i, cor(real_df_agr$fst[which(real_df_agr$cluster15 == i)], real_df_agr$neurogen[which(real_df_agr$cluster15 == i)]))))
}
real_clust15_res$V1 = factor(real_clust15_res$V1, levels = real_clust15_res$V1)
ggplot(real_clust15_res, aes(x = V1, y = V2, color = V1, fill = V1)) + geom_bar(stat = 'identity', alpha = 0.8) + xlab("Cluster 15 Level") + ylab("Real Correlation") + scale_y_continuous(limits = c(0,1), expand = c(0, 0)) + theme_bw()
real_df_agr = aggregate(. ~ subsample + cluster53, data = real_df, mean)
real_clust53_res = data.frame()
for (i in 0:52) {
  real_clust53_res = rbind(real_clust53_res, t(c(i, cor(real_df_agr$fst[which(real_df_agr$cluster53 == i)], real_df_agr$neurogen[which(real_df_agr$cluster53 == i)]))))
}
real_clust53_res$V1 = factor(real_clust53_res$V1, levels = real_clust53_res$V1)
ggplot(real_clust53_res, aes(x = V1, y = V2, color = V1, fill = V1)) + geom_bar(stat = 'identity', alpha = 0.8) + xlab("Cluster 53 Level") + ylab("Real Correlation") + scale_y_continuous(limits = c(0,1), expand = c(0, 0)) + theme_bw()

permNeurogenCor = function(x, clust_level = "bulk", useNuc = T, num_genes = 51) {
  # Set Seed
  set.seed(x)
  
  # Check Input
  if (useNuc != T & useNuc != F) { print("useNuc must be T or F."); return(NULL); }
  if (clust_level != "bulk" & clust_level != "15" & clust_level != "53") { print("clust_level must be 'bulk', '15', or '53'."); return(NULL); }
  if (! is.numeric(num_genes) )          { print("num_genes must be numeric");  return(NULL); }

  # Create Permutation Gene Set
  perm_genes = ran_lists[[x]]
  
  # Calculate Score
  perm_score = colSums(mat[perm_genes,])

  # Find Correlation
  if (! useNuc ) {
    # Finding Correlations b/w Subsamples
    perm_meta = data.frame(neurogen = bb$neurogen, perm = perm_score, subsample = bb$subsample)
    num_clust = 0
    if (clust_level == "15") {
      perm_meta$cluster = bb$seuratclusters15
      num_clust = 14
    } else if (clust_level == "53") {
      perm_meta$cluster = bb$seuratclusters53
      num_clust = 52
    }

    if (num_clust != 0) {
      perm_meta_agr = aggregate(. ~ subsample + cluster, data = perm_meta, mean)
      perm_clust_res = data.frame()
      for ( i in 0:num_clust ) {
        perm_clust_res = rbind(perm_clust_res, t(c(i, cor(perm_meta_agr$perm[which(perm_meta_agr$cluster == i)], perm_meta_agr$neurogen[which(perm_meta_agr$cluster == i)]) )))
      } 
      perm_cor = perm_clust_res[,2]
    } else {
      perm_meta_agr = aggregate(. ~ subsample, data = perm_meta, mean)
      perm_cor = cor(perm_meta_agr$perm, perm_meta_agr$neurogen)
    }
  } else {
    # Finding Correlations in Nuclei
    perm_cor = cor(bb$neurogen, perm_score)
  }
  
  return(perm_cor)
}

library("parallel")
num_perms = 10000

# 10k Perm FST w/ Neurogen: Bulk, Nuc, NonZero
perm_cors = unlist(mclapply( 1:num_perms, function(x) permNeurogenCor(x, clust_level = 'bulk', useNuc = T, num_genes = 51), mc.cores = detectCores() ))
perm_df = data.frame(perm_num = 1:num_perms, perm_cor = perm_cors)
perm_df$isAbove = cor(bb$fst, bb$neurogen) > perm_df$perm_cor
ggplot(perm_df, aes(perm_cor, color = isAbove, fill = isAbove)) + geom_histogram(alpha = 0.5) + scale_fill_manual(values = c("gray40", "goldenrod1"), name = "Real Greater than Perm") + scale_color_manual(values = c("gray40", "goldenrod1"), name = "Real Greater than Perm") + theme_bw() + scale_y_continuous(expand = c(0,0), name = "Number of Perms") + xlab("Correlation w/ Neurogenesis") + ggtitle("10k Perms of 51 Random Nonzero Genes R w/ Neurogenesis Compared to FST R", subtitle = paste0("p-value = ", (num_perms - length(which(perm_df$isAbove)))/ num_perms))
perm_df_bulk_nuc_nonzero = perm_df

pdf("~/scratch/brain/results/pcrc_neurogen_10k.pdf", width = 8, height = 5)
ggplot(perm_df, aes(perm_cor, color = isAbove, fill = isAbove)) + geom_histogram(alpha = 0.5) + scale_fill_manual(values = c("gray40", "goldenrod1"), name = "Real Greater than Perm") + scale_color_manual(values = c("gray40", "goldenrod1"), name = "Real Greater than Perm") + theme_bw() + scale_y_continuous(expand = c(0,0), name = "Number of Perms") + xlab("Correlation w/ Neurogenesis") + ggtitle("10k Perms of 51 Random Nonzero Genes R w/ Neurogenesis Compared to FST R", subtitle = paste0("p-value = ", (num_perms - length(which(perm_df$isAbove)))/ num_perms))
dev.off()

# 10k Perm FST w/ Neurogen: Bulk, Subsample, NonZero
perm_cors = unlist(mclapply( 1:num_perms, function(x) permNeurogenCor(x, clust_level = 'bulk', useNuc = F, num_genes = 51), mc.cores = detectCores() ))
perm_df = data.frame(perm_num = 1:num_perms, perm_cor = perm_cors)
perm_df$isAbove = cor(real_df_agr$fst, real_df_agr$neurogen) > perm_df$perm_cor
perm_df_bulk_sub_nonzero = perm_df
pdf("~/scratch/brain/results/pcrc_neurogen_sub_10k.pdf", width = 8, height = 5)
ggplot(perm_df_bulk_sub_nonzero, aes(perm_cor, color = isAbove, fill = isAbove)) + geom_histogram(alpha = 0.5) + scale_fill_manual(values = c("gray40", "goldenrod1"), name = "Real Greater than Perm") + scale_color_manual(values = c("gray40", "goldenrod1"), name = "Real Greater than Perm") + theme_bw() + scale_y_continuous(expand = c(0,0), name = "Number of Perms") + xlab("Correlation w/ Neurogenesis") + ggtitle("10k Perms of 51 Random Genes R w/ Neurogenesis Compared to FST R", subtitle = paste0("p-value = ", (num_perms - length(which(perm_df$isAbove)))/ num_perms))
dev.off()

# 10k Perm FST w/ Neurogen: 15, Subsample, NonZero
perm_cors = mclapply( 1:num_perms, function(x) permNeurogenCor(x, clust_level = '15', useNuc = F, num_genes = 51), mc.cores = detectCores() )
perm_df = as.data.frame(t(as.data.frame(perm_cors)))
real_clust15_res$num_above = sapply(1:ncol(perm_df), function(x) num_perms - length(which(real_clust15_res[x,2] > perm_df[,x])) )
real_clust15_res$isSig = real_clust15_res$num_above < 50
perm_df_15_sub_nonzero = perm_df
pdf("~/scratch/brain/results/pcrc_neurogen_sub_15_10k.pdf", width = 8, height = 5)
ggplot(real_clust15_res, aes(x=V1, y = num_above, fill = isSig, color = isSig)) + geom_bar(alpha = 0.5, stat = 'identity') + theme_bw() + scale_y_continuous(expand = c(0,0), name = "Number of Perms Greater than Real") + xlab("Cluster") + ggtitle("10k Perms of 51 Random Genes R w/ Neurogenesis Compared to FST R by 15 Cluster") + geom_hline(yintercept = 50, lty = 2) + scale_color_manual(values = c('gray40', "goldenrod1"), name = 'Significant') + scale_fill_manual(values = c('gray40', "goldenrod1"), name = 'Significant') + geom_text(aes(label = num_above), vjust = -0.2)
dev.off()

# 10k Perm FST w/ Neurogen: 53, Subsample, NonZero
perm_cors = mclapply( 1:num_perms, function(x) permNeurogenCor(x, clust_level = '53', useNuc = F, num_genes = 51), mc.cores = detectCores() )
perm_df = as.data.frame(t(as.data.frame(perm_cors)))
real_clust53_res$num_above = sapply(1:ncol(perm_df), function(x) num_perms - length(which(real_clust53_res[x,2] > perm_df[,x])) )
real_clust53_res$isSig = real_clust53_res$num_above < 50
perm_df_53_sub_nonzero = perm_df
pdf("~/scratch/brain/results/pcrc_neurogen_sub_53_10k.pdf", width = 8, height = 5)
ggplot(real_clust53_res, aes(x=V1, y = num_above, fill = isSig, color = isSig)) + geom_bar(alpha = 0.5, stat = 'identity') + theme_bw() + scale_y_continuous(expand = c(0,0), name = "Number of Perms Greater than Real") + xlab("Cluster") + ggtitle("10k Perms of 51 Random Genes R w/ Neurogenesis Compared to FST R by 53 Cluster") + geom_hline(yintercept = 50, lty = 2) + scale_color_manual(values = c('gray40', "goldenrod1"), name = 'Significant') + scale_fill_manual(values = c('gray40', "goldenrod1"), name = 'Significant') + geom_text(aes(label = num_above), vjust = -0.2)
dev.off()

# Couple of sig Effects
perm_df = data.frame(perm_num = 1:num_perms, perm_cor = perm_df_53_sub_nonzero[,34])
perm_df$isAbove = real_clust53_res[34,2] > perm_df$perm_cor
ggplot(perm_df, aes(perm_cor, color = isAbove, fill = isAbove)) + geom_histogram(alpha = 0.5) + scale_fill_manual(values = c("gray40", "goldenrod1"), name = "Real Greater than Perm") + scale_color_manual(values = c("gray40", "goldenrod1"), name = "Real Greater than Perm") + theme_bw() + scale_y_continuous(expand = c(0,0), name = "Number of Perms") + xlab("Correlation w/ Neurogenesis") + ggtitle("10k Perms of 51 Random Nonzero Genes R w/ Neurogenesis Compared to FST R - Cluster 33 (53)", subtitle = paste0("p-value = ", (num_perms - length(which(perm_df$isAbove)))/ num_perms))

#==========================================================================================
# Initial Clustering ======================================================================
#==========================================================================================
# Load Filtered Feature Matrix
b1.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-B1/"))
b2.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-B2/"))
b3.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-B3/"))
b4.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-B4/"))
b5.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-B5/"))
c1.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-C1/"))
c2.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-C2/"))
c3.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-C3/"))
c4.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-C4/"))
c5.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-C5/"))

b1 = CreateSeuratObject(counts = b1.data)
b2 = CreateSeuratObject(counts = b2.data)
b3 = CreateSeuratObject(counts = b3.data)
b4 = CreateSeuratObject(counts = b4.data)
b5 = CreateSeuratObject(counts = b5.data)
c1 = CreateSeuratObject(counts = c1.data)
c2 = CreateSeuratObject(counts = c2.data)
c3 = CreateSeuratObject(counts = c3.data)
c4 = CreateSeuratObject(counts = c4.data)
c5 = CreateSeuratObject(counts = c5.data)

# Add Metdata
b1$sample = "b1"
b2$sample = "b2"
b3$sample = "b3"
b4$sample = "b4"
b5$sample = "b5"
c1$sample = "c1"
c2$sample = "c2"
c3$sample = "c3"
c4$sample = "c4"
c5$sample = "c5"

b1$cond = "bhve"
b2$cond = "bhve"
b3$cond = "bhve"
b4$cond = "bhve"
b5$cond = "bhve"
c1$cond = "ctrl"
c2$cond = "ctrl"
c3$cond = "ctrl"
c4$cond = "ctrl"
c5$cond = "ctrl"

bb = merge(b1, list(b2,b3,b4,b5,c1,c2,c3,c4,c5), add.cell.ids = c("b1", "b2", "b3", "b4", "b5", "c1", "c2", "c3", "c4", "c5"))

rm(b1.data); rm(b1)
rm(b2.data); rm(b2)
rm(b3.data); rm(b3)
rm(b4.data); rm(b4)
rm(b5.data); rm(b5)
rm(c1.data); rm(c1)
rm(c2.data); rm(c2)
rm(c3.data); rm(c3)
rm(c4.data); rm(c4)
rm(c5.data); rm(c5)

# Quality Control: Mitochondrial DNA, Dead Cells, and Doublets
# mt_genes = scan("~/scratch/m_zebra_ref/mt.txt", what = "character")
mt_genes = scan("C:/Users/miles/Downloads/all_research/mt.txt", what = "character")
mt_genes = str_replace(mt_genes, "_", "-")
mt_genes = mt_genes[which(mt_genes %in% rownames(bb))]
bb$pct_mt = PercentageFeatureSet(bb, features = mt_genes)
FeatureScatter(bb, feature1 = "nCount_RNA", feature2 = "pct_mt") + NoLegend()
FeatureScatter(bb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
print(paste("Number of Cells in BB Before Filterning:", ncol(bb)))

bb = subset(bb, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & pct_mt < 5)
print(paste("Number of Cells in BB After Filterning:", ncol(bb)))
bb = NormalizeData(bb, normalization.method = "LogNormalize", scale.factor = 10000)
bb = FindVariableFeatures(bb, selection.method = "vst", nfeatures = 2000)
bb = ScaleData(bb, features = rownames(bb))
bb = RunPCA(bb, features = VariableFeatures(bb))
bb = RunUMAP(bb, reduction = "pca", dims = 1:50)
bb = FindNeighbors(bb, reduction="umap", dims = 1:2)
bb = FindClusters(bb, resolution = .30)
# DimPlot(bb, reduction = "umap", split.by = "cond", label = TRUE)

png_name = paste0(rna_path, "/results/bb/umap_10k_50_30.png")
png(file = png_name, width = 2000, height = 1500, res = 200)
print(DimPlot(bb, reduction = "umap", label = TRUE, pt.size = 0.7))
dev.off()
system(paste0("rclone copy ", png_name, " dropbox:BioSci-Streelman/George/Brain/bb/tmp/"))

#=======================================================================================
# GO Enrichment ========================================================================
#=======================================================================================
library("parallel")
numCores = detectCores()
all_go = read.csv("~/research/brain/data/cluster_compare_all_levels_MF_042021.csv", stringsAsFactors = F)
# all_go = read.csv("~/research/brain/data/53cluster_compare_all_levels_MF_042021.csv", stringsAsFactors = F)
all_go = read.csv("~/Downloads/cluster15_compare_all_levels_CC_042121.csv", stringsAsFactors = F)
all_go$p = -1
all_go$bh = -1
all_go$bon = -1
levels_categories = data.frame()
for (level_cur in 2:13) {
  print(level_cur)
  go_cur_idx = which(all_go$level == level_cur)
  go_cur = all_go[go_cur_idx,]
  levels_categories = rbind(levels_categories, data.frame(Description = unique(go_cur$Description), level = level_cur))
  go_cur = cbind(go_cur, colsplit(go_cur$GeneRatio, pattern = "/", names = c("Hit", "Total")))
  go_cur_tots = aggregate(Hit ~ Description, go_cur, sum)
  go_cur_tots$Total = aggregate(Total ~ Description, go_cur, sum)[,2]
  go_cur$Hit_Total = go_cur_tots$Hit[match(go_cur$Description, go_cur_tots$Description)]
  go_cur$Total_Total = go_cur_tots$Total[match(go_cur$Description, go_cur_tots$Description)]
  
  # Set up contingency table values
  # a b
  # c d
  go_cur$a = go_cur$Hit
  go_cur$b = go_cur$Hit_Total - go_cur$Hit
  go_cur$c = go_cur$Total - go_cur$Hit
  go_cur$d = go_cur$Total_Total - go_cur$a - go_cur$b - go_cur$c
  
  # Calculate p value
  findP = function(x) {
    contig_table = matrix(c(go_cur[x,"a"], go_cur[x,"b"], go_cur[x,"c"], go_cur[x,"d"]), nrow = 2, byrow = T)
    fisher.test(contig_table, alternative = "greater")$p.value
  }
  go_cur$p = unlist(mclapply(1:nrow(go_cur), findP, mc.cores = numCores))
  go_cur_adj = mclapply(unique(go_cur$ID), function(x) { data.frame(X = go_cur$X[which(go_cur$ID == x)], bh = p.adjust(go_cur$p[which(go_cur$ID == x)], method = "BH"), bon = p.adjust(go_cur$p[which(go_cur$ID == x)], method = "bonferroni") ) }, mc.cores = numCores)
  go_cur_adj = do.call(rbind, go_cur_adj)
  all_go$p[go_cur_idx] = go_cur$p[match(all_go$X[go_cur_idx], go_cur$X)]
  all_go[go_cur_idx, c("bh", "bon")] = go_cur_adj[match(all_go$X[go_cur_idx], go_cur_adj$X), c("bh", "bon")]
}
# all_go$bh = p.adjust(all_go$p, method = "BH")
# all_go$bon = p.adjust(all_go$p, method = "bonferroni")
all_go$new_cluster = changeClusterID(all_go$Cluster, returnFactor = T)
all_go$Cluster_ID = paste0(all_go$Cluster, "_", all_go$ID)
all_go$level_ID = paste0(all_go$level, "_", all_go$ID)
all_go$neg_log_p = -log10(all_go$p)
all_go_sig = all_go[which(all_go$bh < 0.05),]
all_go_sig = all_go_sig[which(! duplicated(all_go_sig$Cluster_ID) ),]

top.x = 5
# all_go_sig_top = mclapply(rev(levels(all_go$new_cluster)), function(x) { tmp = all_go_sig[which(all_go_sig$new_cluster == x),]; tmp[order(tmp$level, rev(tmp$p), decreasing = T)[1:top.x],] }, mc.cores = numCores)
all_go_sig_top = mclapply(rev(levels(all_go$new_cluster)), function(x) { tmp = all_go_sig[which(all_go_sig$new_cluster == x),]; tmp[order(tmp$p, decreasing = F)[1:top.x],] }, mc.cores = numCores)
all_go_sig_top = do.call(rbind, all_go_sig_top)
all_go_sig_top = all_go_sig_top[which(! is.na(all_go_sig_top$ID) ),]
all_go_sig_top$Description_level = paste0(all_go_sig_top$Description, " (", all_go_sig_top$level, ")")
all_go_sig_top$col = convert15$col[match(all_go_sig_top$Cluster, convert15$old)]

all_go_sig_top_p = all_go[which(all_go$level_ID %in% all_go_sig_top$level_ID),]
all_go_sig_top_p$Description_level = paste0(all_go_sig_top_p$Description, " (", all_go_sig_top_p$level, ")")
all_go_sig_top_p$Description_level = factor(all_go_sig_top_p$Description_level, levels = unique(all_go_sig_top$Description_level))
all_go_sig_top_p$col = all_go_sig_top$col[match(all_go_sig_top_p$X, all_go_sig_top$X)]
scaled_nlp = mclapply(levels(all_go_sig_top_p$Description_level), function(x) data.frame(X = all_go_sig_top_p$X[which(all_go_sig_top_p$Description_level == x)], scaled_nlp = scale(all_go_sig_top_p$neg_log_p[which(all_go_sig_top_p$Description_level == x)])), mc.cores = numCores)
scaled_nlp = do.call(rbind, scaled_nlp)
all_go_sig_top_p$scaled_nlp = scaled_nlp$scaled_nlp[match(all_go_sig_top_p$X, scaled_nlp$X)]

# Matrix GO (rows) by cluster (columns). Colored by -log10(Diff Enrich P)
pdf("~/research/brain/results/supplement/cluster_deg_diff_enrich_top_sig.pdf", width = 10, height = 9)
ggplot(all_go_sig_top_p, aes(x = new_cluster, y = Description_level, fill = scaled_nlp)) + geom_tile() + scale_fill_viridis_c(name = expression("Scaled -"*Log["10"]*"(p)")) + coord_fixed() + xlab("") + ylab("")  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = convert15$col, size = 10), axis.text.y = element_text(color = all_go_sig_top$col[which(! duplicated(all_go_sig_top$Description_level))]))
dev.off()

# Matrix GO (rows) by cluster (columns). Colored by -log10(Diff Enrich P) and anything that is sig is automatically yellow
pdf("~/research/brain/results/supplement/cluster_deg_diff_enrich_cap_top_sig.pdf", width = 10, height = 9)
all_go_sig_top_p2 = all_go_sig_top_p
all_go_sig_top_p2$neg_log_p[which(all_go_sig_top_p2$neg_log_p >= min(all_go_sig_top$neg_log_p))] = NA
print(ggplot(all_go_sig_top_p2, aes(x = new_cluster, y = Description_level, fill = neg_log_p)) + geom_tile() + scale_fill_viridis_c(name = expression("-"*Log["10"]*"(p)"), na.value = "#FDE725") + coord_fixed() + xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = convert15$col, size = 10), axis.text.y = element_text(color = all_go_sig_top$col[which(! duplicated(all_go_sig_top$Description_level))])))
dev.off()

# Boxplot of -log10(Diff Enrich P)
all_go_sig_top_p$isSig = all_go_sig_top_p$Cluster_ID %in% all_go_sig_top$Cluster_ID
ggplot(all_go_sig_top_p, aes(x = neg_log_p, y = Description_level)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(color = new_cluster, size = isSig)) + theme_bw() + scale_size_manual(values = c(1, 3)) + scale_color_manual(values = convert15$col) + scale_x_continuous(expand = c(0,0.1)) + xlab(expression("-"*Log["10"]*"(p)")) + ylab("") + theme(axis.text.y = element_text(color = all_go_sig_top$col[which(! duplicated(all_go_sig_top$Description_level))]))

# all_go$bh_sig = all_go$bh < 0.05
# for (level_cur in 2:13) {
#   go_cur = all_go[which(all_go$level == level_cur),]
#   
#   # Loose
#   go_cur$bh = p.adjust(go_cur$p, method = "BH")
#   go_cur$bh_sig = go_cur$bh < 0.05
#   
#   go_cur_good_categories = aggregate(bh_sig ~ Description, go_cur, sum)
#   go_cur_good_categories = go_cur_good_categories$Description[which(go_cur_good_categories$bh_sig > 0)]
#   go_cur = go_cur[which(go_cur$Description %in% go_cur_good_categories),]
#   png(paste0("~/research/brain/results/go/dif_enrich_loose_53/", level_cur, ".png"), width = 1000, height = 1000, res = 90)
#   print(ggplot(go_cur, aes(x = Description, y = new_cluster, fill = bh_sig)) + geom_tile() + scale_fill_viridis_d() + ggtitle(paste0("Significantly Differentially Enriched MF Categories by BH at Level ", level_cur)) + xlab("Category") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + coord_fixed())
#   dev.off()
# }

# Top Hits
# top_go_15 = read.csv("~/research/brain/data/15cluster_enrichGO_compareCluster_MF_042021.csv", stringsAsFactors = F)
# top_go_15 = read.csv("~/research/brain/data/53cluster_enrichGO_compareCluster_MF_042021.csv", stringsAsFactors = F)
# top_go_15 = read.csv("~/research/brain/data/15cluster_enrichGO_compareCluster_MF_all_pvalues_042121.csv", stringsAsFactors = F)
# top_go_15 = read.csv("~/research/brain/data/53cluster_enrichGO_compareCluster_MF_all_pvalues_all_set_sizes_test_042121.csv", stringsAsFactors = F)
# top_go_15 = read.csv("~/research/brain/data/15cluster_enrichGO_compareCluster_CC_all_pvalues_all_set_sizes_test_042121.csv", stringsAsFactors = F)
# top_go_15 = read.csv("~/research/brain/data/53cluster_enrichGO_compareCluster_CC_all_pvalues_all_set_sizes_test_042121.csv", stringsAsFactors = F)
# top_go_15 = read.csv("~/research/brain/data/53cluster_enrichGO_compareCluster_BP_all_pvalues_all_set_sizes_test_042121.csv", stringsAsFactors = F)
top_go_15$new_cluster = changeClusterID(top_go_15$Cluster, returnFactor = T)
top_go_15$neg_log_p_adj = -log10(top_go_15$p.adjust)
# top_go_15$level = levels_categories$level[match(top_go_15$Description, levels_categories$Description)]
for (level_cur in 1) {
# for (level_cur in 2:13) {
  
  go_cur = top_go_15
  # go_cur = top_go_15[which(top_go_15$level == level_cur),]
  go_cur2 = data.frame()
  sig_categories = unique(go_cur$Description[which(go_cur$p.adjust < 0.05)])
  if ( length(sig_categories) > 1 ) {
    for (cat_cur in sig_categories) {
      # For sig/non-sig
      # newRows = as.data.frame(table(go_cur$new_cluster[which(go_cur$Description == cat_cur)]))
      # newRows$Cat = cat_cur
      # For continuous
      newRows = go_cur[which(go_cur$Description == cat_cur), c("new_cluster", "neg_log_p_adj", "Description")]
      go_cur2 = rbind(go_cur2, newRows)
    }
    colnames(go_cur2) = c("new_cluster", "sig", "Description")
    
    fullRows = expand.grid(unique(top_go_15$new_cluster),sig_categories)
    fullRows$pvalue = 1
    colnames(fullRows) = c("new_cluster", "Description", "sig")
    
    go_cur2 = rbind(go_cur2, fullRows[which(! paste(fullRows$new_cluster, fullRows$Description) %in% paste(go_cur2$new_cluster, go_cur2$Description) ), c("new_cluster", "sig", "Description")])
    
    this_height = 100+33*length(sig_categories)
    if (this_height < 400)
      this_height = 400
    
    go_cur2_mat = acast(go_cur2, Description ~ new_cluster, value.var = "sig")
    old_rownames = rownames(go_cur2_mat)
    rownames(go_cur2_mat) = substr(rownames(go_cur2_mat), 1, 50L)
    rownames(go_cur2_mat)[which(old_rownames != rownames(go_cur2_mat))] = paste0(str_trim(rownames(go_cur2_mat)[which(old_rownames != rownames(go_cur2_mat))]), "...")
    
    # All
    png(paste0("~/research/brain/results/go/enrich_53_BP_scale.png"), width = 2500,  height = this_height, res = 90)
    pheatmap::pheatmap(go_cur2_mat, color = colorRampPalette(viridis(100))(100), na_col = viridis(1), cluster_cols = F, cluster_rows = T, angle_col = 45, scale = "row", cellwidth = 25, cellheight = 25, border_color = "black")
    # print(ggplot(go_cur2, aes(x = Description, y = new_cluster, fill = sig)) + geom_tile() + scale_fill_viridis_c() + ggtitle(paste0("Significantly Enriched MF Categories by BH at Level ", level_cur)) + xlab("Category") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + coord_fixed() + theme(panel.background = element_rect(fill = viridis(1)), panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
    dev.off()
    
    # By level
    # png(paste0("~/research/brain/results/go/enrich_15/", level_cur, "_enrich.png"), width = 2500,  height = this_height, res = 90)
    # pheatmap::pheatmap(go_cur2_mat, color = colorRampPalette(viridis(100))(100), na_col = viridis(1), cluster_cols = F, cluster_rows = T, angle_col = 45, scale = "none", cellwidth = 25, cellheight = 25, border_color = "black")
    # # print(ggplot(go_cur2, aes(x = Description, y = new_cluster, fill = sig)) + geom_tile() + scale_fill_viridis_c() + ggtitle(paste0("Significantly Enriched MF Categories by BH at Level ", level_cur)) + xlab("Category") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + coord_fixed() + theme(panel.background = element_rect(fill = viridis(1)), panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
    # dev.off()
  }
}

#=======================================================================================
# Cluster Size on DEG ==================================================================
#=======================================================================================
my_data = matrix(0L, nrow = 6, ncol = 1000, dimnames = list(paste0("gene", 1:6), paste0("cell", 1:1000)))

# Pseudo BHVE Cluster 1 is 1:200
my_data[1, sample(1:200, 40)] = 2 
my_data[6, sample(1:200, 40)] = rnorm(40, mean = 2, sd = 1)
my_data[2, 1:200] = 1
my_data[3, 1:200] = 0

# Pseudo CTRL Cluster 1 is 201:400
my_data[1, sample(201:400, 40)] = 1
my_data[6, sample(201:400, 40)] = rnorm(40, mean = 1, sd = 1)
my_data[2, 201:400] = 1
my_data[3, 201:400] = 0

# Pseduo BHVE Cluster 2 is 401:700
my_data[1, sample(401:700, 60)] = 2
my_data[6, sample(401:700, 60)] = rnorm(60, mean = 2, sd = 1)
my_data[2, 401:700] = 0
my_data[3, 401:700] = 1

# Pseudo CTRL Cluster 2 is 701: 1000
my_data[1, sample(701:1000, 60)] = 1
my_data[6, sample(701:1000, 60)] = rnorm(60, mean = 1, sd = 1)
my_data[2, 701:1000] = 0
my_data[3, 701:1000] = 1

my_data[4,] = my_data[3,]*2
my_data[5,] = my_data[2,]*2
pseduo_obj = CreateSeuratObject(my_data)
pseduo_obj$seurat_clusters = c(rep(1, 400), rep(2, 600))
pseduo_obj$cond = c(rep("BHVE", 200), rep("CTRL", 200), rep("BHVE", 300), rep("CTRL", 300))
pseduo_obj$clust_cond = paste0(pseduo_obj$seurat_clusters, ".", pseduo_obj$cond)
Idents(pseduo_obj) = pseduo_obj$clust_cond
pseduo_obj = NormalizeData(pseduo_obj)
pseduo_deg = FindMarkers(pseduo_obj, ident.1 = "1.BHVE", ident.2 = "1.CTRL")
pseduo_deg = rbind(pseduo_deg, FindMarkers(pseduo_obj, ident.1 = "2.BHVE", ident.2 = "2.CTRL"))

pseduo_obj = ScaleData(pseduo_obj, features = rownames(pseduo_obj))
pseduo_obj = RunPCA(pseduo_obj, features = c("gene2", "gene3", "gene4", "gene5"))
pseduo_obj = RunUMAP(pseduo_obj, reduction = "pca", dims = 1:2)
# pseduo_obj = FindNeighbors(pseduo_obj, reduction="umap", dims = 1:2)
# pseduo_obj = FindClusters(pseduo_obj, resolution = .30)
Idents(pseduo_obj) = pseduo_obj$seurat_clusters
DimPlot(pseduo_obj, reduction = "umap", split.by = "cond", label = TRUE)
FeaturePlot(pseduo_obj, "gene1", order = T, label = T, split.by = "cond")
FeaturePlot(pseduo_obj, "gene6", order = T, label = T, split.by = "cond")
myFeaturePlot(pseduo_obj, "gene1", my.split.by = "cond")
myFeaturePlot(pseduo_obj, "gene6", my.split.by = "cond")


df = data.frame()
for (old_cluster in unique(old$seurat_clusters)) {
  
}

#======================================================================================
# Gene Set ============================================================================
#======================================================================================
neuro_df = as.data.frame(readxl::read_excel("~/research/brain/data/aau5324_Moffitt_Table-S3.xlsx"))
colnames(neuro_df) = c("Neuropeptides", "TF", "Neuromodulators")
neuro_df = neuro_df[-1,]
neuro_df[1:87,1] = neuro_df[2:88,1]
neuro_df[77:85,1] = neuro_df[78:86,1]
neuro_df[86,1] = NA
cluster15_deg = read.csv("~/research/brain/results/bb_all_markers_15clusters_102820_more_info.csv")

neuro_deg = data.frame()
for (cluster15 in levels(bb$seuratclusters15)) {
  deg_rows = cluster15_deg[which(cluster15_deg$cluster == cluster15 & cluster15_deg$p_val_adj < 0.05 & cluster15_deg$avg_logFC > 0),]
  deg_genes = deg_rows$human[which(! is.na(deg_rows$human) )]
  ranks_logFC = deg_rows$avg_logFC[which(! is.na(deg_rows$human) )]
  ranks_p = deg_rows$p_val[which(! is.na(deg_rows$human) )]
  names(ranks_logFC) = deg_genes
  names(ranks_p) = deg_genes
  for (gene_list in colnames(neuro_df)) {
    this_list = neuro_df[,c(gene_list)]
    this_list = this_list[which(!is.na(this_list))]
    this_list = toupper(this_list)
    num_ovlp = length(which(deg_genes %in% this_list))
    
    this_list_list = list(this_list)
    names(this_list_list) = c(gene_list)
    
    if (num_ovlp > 0) {
      p_logFC = fgsea(this_list_list, ranks_logFC)$pval
      p_p = fgsea(this_list_list, ranks_p)$pval
      png(paste0("~/research/brain/results/fgsea/fgsea_logFC_", cluster15, "_", gene_list, ".png"), width = 400, height = 400)
      print(plotEnrichment(this_list_list[[gene_list]], ranks_logFC))
      dev.off()
      png(paste0("~/research/brain/results/fgsea/fgsea_p_", cluster15, "_", gene_list, ".png"), width = 400, height = 400)
      print(plotEnrichment(this_list_list[[gene_list]], ranks_p))
      dev.off()
    } else { p_logFC = 1; p_p = 1 }
    
    neuro_deg = rbind(neuro_deg, t(c(cluster15, gene_list, p_logFC, p_p, num_ovlp, length(deg_genes), length(this_list))))
  }
}
colnames(neuro_deg) = c("cluster", "gene_list", "p_logFC", "p_p", "num_ovlp", "num_degs", "num_list")
neuro_deg[,c("cluster", "p_logFC", "p_p", "num_ovlp", "num_degs", "num_list")] <- lapply(neuro_deg[,c("cluster", "p_logFC", "p_p", "num_ovlp", "num_degs", "num_list")], as.numeric)
neuro_deg$neg_logp_logFC = -log10(neuro_deg$p_logFC)
neuro_deg$neg_logp_p = -log10(neuro_deg$p_p)
neuro_deg$new = convert15$new.full[match(neuro_deg$cluster, convert15$old)]
neuro_deg$new = factor(neuro_deg$new, levels = convert15$new.full)
neuro_deg$corrected_ovlp = neuro_deg$num_ovlp / neuro_deg$num_degs / neuro_deg$num_list

ggplot(neuro_deg, aes(x=new,y=gene_list, fill = neg_logp_logFC)) + geom_tile() + scale_fill_viridis_c(name = "-log(p)") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("") + ggtitle("GSEA on Avg LogFC Ranks")
ggplot(neuro_deg, aes(x=new,y=gene_list, fill = neg_logp_p)) + geom_tile() + scale_fill_viridis_c(name = "-log(p)") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("") + ggtitle("GSEA on p-value Ranks")
ggplot(neuro_deg, aes(x=new,y=gene_list, fill = num_ovlp)) + geom_tile() + scale_fill_viridis_c(name = "Ovlp") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("") + ggtitle("Raw Overlap")
ggplot(neuro_deg, aes(x=new,y=gene_list, fill = corrected_ovlp)) + geom_tile() + scale_fill_viridis_c(name = "Corr. Ovlp") + coord_fixed() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("") + ggtitle("Corrected Overlap")
write.csv(neuro_deg, "~/research/brain/results/fgsea_results_04132021.csv")

#======================================================================================
# For Katie ===========================================================================
#======================================================================================
# old = readRDS("C:/Users/miles/Downloads/brain/brain_scripts/brain_mz_shiny/data/B1C1C2MZ_combined_031020.rds")
# old = readRDS("~/Downloads/B1C1C2MZ_combined_katie.rds")
old = readRDS("C:/Users/miles/Downloads/combined.sct_ensemble.rds")
# old = readRDS("C:/Users/miles/Downloads/brain/data/bb_subsample_02222021.RDS")
new = readRDS("C:/Users/miles/Downloads/combined.sct_ncbi.rds")
# new = readRDS("C:/Users/miles/Downloads/B1C1C2MZ_combinedv3.rds")
# new = readRDS("~/Downloads/mcmzbb_combined_vst_1.8res.rds")
# new = readRDS("C:/Users/miles/Downloads/mcmz_bb_combinedv2.rds")
# new = readRDS("~/Downloads/combined_sctransform1.8res.rds")

new_cells = sapply(colnames(new), function(x) substr(x, 0, 21))
old_cells = sapply(colnames(old), function(x) substr(x, 0, 21))

new_df = data.frame(new_cell = colnames(new), new_sample = substr(new$sample, 0, 2), new_cluster = new$seurat_clusters)
old_df = data.frame(old_cell = colnames(old), old_sample = old$sample, old_cluster = old$seurat_clusters)

new_df$cell_strip = str_remove(new_df$new_cell, "MCMZ_")
new_df$cell_strip = colsplit(new_df$cell_strip, "_", names = c("begin", "barcode"))[,2]
new_df$cell_strip = colsplit(new_df$cell_strip, "-", names = c("barcode", "end"))[,1]

old_df$old_sample = str_replace(old_df$old_sample, "MC_BHVE", "MC")
old_df$old_sample = str_replace(old_df$old_sample, "MC_CTRL", "MC")
old_df$old_sample = str_replace(old_df$old_sample, "MZ_CTRL", "MZ")
old_df$cell_strip = str_remove(old_df$old_cell, "MCMZ_")
old_df$cell_strip = colsplit(old_df$cell_strip, "_", names = c("begin", "barcode"))[,2]
old_df$cell_strip = colsplit(old_df$cell_strip, "-", names = c("barcode", "end"))[,1]

new_df[,c("old_cell", "old_cluster")] = old_df[match( paste0(new_df$cell_strip, new_df$new_sample), paste0(old_df$cell_strip, old_df$old_sample) ), c("old_cell", "old_cluster")]
df = new_df
df = df[which(df$new_sample == "MZ"),]

df2_mat = matrix(0L, nrow = length(unique(old$seurat_clusters)), ncol = length(unique(new$seurat_clusters)), dimnames = list(levels(old$seurat_clusters), levels(new$seurat_clusters)))
df2 = data.frame()
for (old_cluster in unique(old$seurat_clusters)) {
  for (new_cluster in unique(new$seurat_clusters)) {
    combo = length(which(df$old_cluster == old_cluster & df$new_cluster == new_cluster))
    df2 = rbind(df2, t(c(old_cluster, new_cluster, combo)))
    df2_mat[old_cluster, new_cluster] = combo
  }
}
colnames(df2) = c("old_cluster", "new_cluster", "value")

df2$value = as.numeric(as.vector(df2$value))
df2$old_cluster = factor(df2$old_cluster, levels = levels(old$seurat_clusters))
df2$new_cluster = factor(df2$new_cluster, levels = levels(new$seurat_clusters))

png("~/Downloads/katie_to_sc_big_bar_rock.png", width = 1000, height = 600)
ggplot(df2, aes(old_cluster, value, fill = new_cluster, color = new_cluster)) + geom_bar(stat = "identity") + xlab("Katie Clusters") + ylab("Number of Cells from New Cluster in Katie Cluster") + ggtitle("Map of New Clusters to Katie's Clusters")
dev.off()
png("~/Downloads/katie_ens_v_ncbi.png", width = 750, height = 750)
heatmap.2(df2_mat, scale = "row", Rowv=FALSE, Colv=FALSE, cexRow = 1, cexCol = 1, dendrogram = "none", trace = "none", xlab = "NCBI Clusters", ylab = "ENSEMBL Clusters")
dev.off()
df2_mat_melt = melt(df2_mat)
png("~/Downloads/katie_to_sc_big_heat_3_rock.png", width = 750, height = 750)
ggplot(df2_mat_melt, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + scale_fill_viridis() + xlab("Kaite's Clusters") + ylab("New Clusters")
dev.off()
write.csv(df2_mat, "~/Downloads/katie_to_sc_rock.csv")

png("C:/Users/miles/Downloads/zackbb_to_katie_big_bar.png", width = 1000, height = 600)
ggplot(df2, aes(old_cluster, value, fill = new_cluster, color = new_cluster)) + geom_bar(stat = "identity") + xlab("Zack BB Clusters") + ylab("Number of Cells from Katie Cluster in BB Cluster") + ggtitle("Map of Katie's Clusters to Zack's BB Clusters")
dev.off()

png("C:/Users/miles/Downloads/zack_bb_to_katie_big_heat_2.png", width = 750, height = 750)
heatmap.2(df2_mat, scale = "row", Rowv=FALSE, Colv=FALSE, cexRow = 1, cexCol = 1, dendrogram = "none", trace = "none", xlab = "Katie's Clusters", ylab = "Zack's BB Clusters")
dev.off()

df2_mat_melt = melt(df2_mat)
colnames(df2_mat_melt) = c("Zack", "Katie", "Value")
df2_mat_melt$Zack = factor(df2_mat_melt$Zack, levels = 0:max(df2_mat_melt$Zack))
df2_mat_melt$Katie = factor(df2_mat_melt$Katie, levels = 0:max(df2_mat_melt$Katie))
png("C:/Users/miles/Downloads/zack_bb_to_katie_big_heat_3.png", width = 750, height = 750)
ggplot(df2_mat_melt, aes(x = Zack, y = Katie, fill = Value)) + geom_tile() + scale_fill_viridis() + xlab("Zack's BB Clusters") + ylab("Katie's Clusters")
dev.off()
write.csv(df2_mat, "C:/Users/miles/Downloads/zack_bb_to_katie.csv")


# Look for Batch Effects
pdf(paste0(rna_path, "/results/b1b2c1mz_ncount_RNA.pdf"), width = 7, height = 5)
col_pal = col_pal = c("#9d0208", "#dc2f02", "#023e8a", "#0077b6")
bvcVis(old, "nCount_RNA", mode = "violin") + xlab("Number of Transcripts per Cell") + scale_fill_manual(values=col_pal) + scale_color_manual(values=col_pal)
dev.off()

pdf(paste0(rna_path, "/results/b1b2c1mz_nfeature_RNA.pdf"), width = 7, height = 5)
col_pal = col_pal = c("#9d0208", "#dc2f02", "#023e8a", "#0077b6")
bvcVis(old, "nFeature_RNA", mode = "violin") + xlab("Number of Features per Cell") + scale_fill_manual(values=col_pal) + scale_color_manual(values=col_pal)
dev.off()

df_count_stats = data.frame(sample = c("b1", "b2", "c1", "c2"), mean = rep(0, 4), median = rep(0,4))
df_feature_stats = data.frame(sample = c("b1", "b2", "c1", "c2"), mean = rep(0, 4), median = rep(0,4))
rownames(df_count_stats) = df_count_stats$sample
rownames(df_feature_stats) = df_feature_stats$sample
for (sample_i in df_count_stats$sample) {
  print(sample_i)
  df_count_stats[sample_i,c(2,3)] = c( mean(old$nCount_RNA[WhichCells(old, idents = sample_i)]), median(old$nCount_RNA[WhichCells(old, idents = sample_i)]) )
  df_feature_stats[sample_i,c(2,3)] = c( mean(old$nFeature_RNA[WhichCells(old, idents = sample_i)]), median(old$nFeature_RNA[WhichCells(old, idents = sample_i)]) )
}
df_count_stats_2 = melt(df_count_stats)
df_feature_stats_2 = melt(df_feature_stats)

pdf(paste0(rna_path, "/results/b1b2c1mz_ncount_RNA_bar.pdf"), width = 7, height = 5)
col_pal = col_pal = c("#9d0208", "#dc2f02", "#023e8a", "#0077b6")
ggplot(df_count_stats_2, aes(sample, value, fill = variable, color = variable)) + geom_bar(stat = "identity", position = position_dodge2()) + xlab("Sample")
dev.off()

pdf(paste0(rna_path, "/results/b1b2c1mz_nfeature_RNA_bar.pdf"), width = 7, height = 5)
col_pal = col_pal = c("#9d0208", "#dc2f02", "#023e8a", "#0077b6")
ggplot(df_feature_stats_2, aes(sample, value, fill = variable, color = variable)) + geom_bar(stat = "identity", position = position_dodge2()) + xlab("Sample")
dev.off()

Idents(old) = old$seurat_clusters
for (cluster in levels(old$seurat_clusters)) {
  this_cells = WhichCells(old, idents = cluster)
  # Look for Batch Effects
  pdf(paste0(rna_path, "/results/b1b2c1mz_batch/", cluster, "_ncount_vln.pdf"), width = 7, height = 5)
  col_pal = col_pal = c("#9d0208", "#dc2f02", "#023e8a", "#0077b6")
  print(bvcVis(old, "nCount_RNA", mode = "violin", cells.use = this_cells) + xlab("Number of Transcripts per Cell") + scale_fill_manual(values=col_pal) + scale_color_manual(values=col_pal) + labs(subtitle = paste0("Cluster ", cluster)))
  dev.off()
  
  pdf(paste0(rna_path, "/results/b1b2c1mz_batch/", cluster, "_nfeature_vln.pdf"), width = 7, height = 5)
  col_pal = col_pal = c("#9d0208", "#dc2f02", "#023e8a", "#0077b6")
  print(bvcVis(old, "nFeature_RNA", mode = "violin", cells.use = this_cells) + xlab("Number of Features per Cell") + scale_fill_manual(values=col_pal) + scale_color_manual(values=col_pal) + labs(subtitle = paste0("Cluster ", cluster)))
  dev.off()
}

df_batch = data.frame()
df_batch2 = data.frame()
Idents(old) = old$seurat_clusters
for (cluster in levels(old$seurat_clusters)) {
  cluster_cells = WhichCells(old, idents = cluster)
  mc_cells = colnames(old)[which(old$sample %in% c("b1", "b2", "c1"))]
  this_cells = cluster_cells[which( cluster_cells %in% mc_cells )]
  mc_pct = length(this_cells)/length(mc_cells) * 100
  mc_mean_nCount = mean(old$nCount_RNA[this_cells])

  sample = "c2"
  sample_cells = colnames(old)[which(old$sample == sample)]
  this_cells = cluster_cells[which( cluster_cells %in% sample_cells )]
  mz_pct = length(this_cells)/length(sample_cells) * 100
  mz_mean_nCount = mean(old$nCount_RNA[this_cells])
  df_batch = rbind(df_batch, t(c(cluster, mz_pct - mc_pct, mz_mean_nCount - mc_mean_nCount)))
  for (sample in c("b1", "b2", "c1", "c2")) {
    sample_cells = colnames(old)[which(old$sample == sample)]
    this_cells = cluster_cells[which( cluster_cells %in% sample_cells )]
    df_batch2 = rbind(df_batch2, t(c(cluster, sample, length(this_cells)/length(sample_cells), mean(old$nCount_RNA[this_cells]), mean(old$nFeature_RNA[this_cells]))))
  }
}
colnames(df_batch) = c("cluster", "mz_mc_pct_dif", "mz_mc_mean_nCount_dif")
colnames(df_batch2) = c("cluster", "sample", "pct_cluster_of_sample", "mean_nCount", "mean_nFeature")
df_batch$unique = "non-unique"
df_batch$unique[which( df_batch$cluster %in% c(1,2,3,6,17,24) )] = "sand"
df_batch$unique[which( df_batch$cluster %in% c(13, 29) )] = "rock"
df_batch2$unique = "non-unique"
df_batch2$unique[which( df_batch2$cluster %in% c(1,2,3,6,17,24) )] = "sand"
df_batch2$unique[which( df_batch2$cluster %in% c(13, 29) )] = "rock"

df_batch$mz_mc_pct_dif = as.numeric(as.vector(df_batch$mz_mc_pct_dif))
df_batch$mz_mc_mean_nCount_dif = as.numeric(as.vector(df_batch$mz_mc_mean_nCount_dif))
df_batch$abs_nCount_dif = abs(df_batch$mz_mc_mean_nCount_dif)
df_batch2$mean_nCount = as.numeric(as.vector(df_batch2$mean_nCount))

png("C:/Users/miles/Downloads/reads_per_nuc_rock_sand.png", width = 750, height = 450)
ggplot(df_batch, aes(mz_mc_mean_nCount_dif, mz_mc_pct_dif, fill = unique, color = unique, size = 0.6)) + geom_point() + ylab("% Difference in MZ vs MC Nuclei") + xlab("Difference in # Reads per Nuclei") + ggtitle("Reads per Nuclei Differences in Rock and Sand Unique Clusters")
dev.off()

ggplot(df_batch2, aes(cluster, mean_nCount, fill = sample, color = sample)) + geom_bar(stat="identity", position = position_dodge2())


old_down = downsampleObj(old)
saveRDS(old_down, "C:/Users/miles/Downloads/b1b2c1mz_down_03092021.rds")

# Z Score Learning Exercise
df_bvc_plot3 = read.csv("C:/Users/miles/Downloads/ieg_covar_c53_p1000_bvc_summary.csv")
raw = as.numeric(as.vector(df_bvc_plot3[which(df_bvc_plot3$cluster1 == 18 & df_bvc_plot3$cluster2 == 13), c(4:1004)]))
scaled = scale(raw)
hist(scaled, breaks = 50)
mean(scaled)
z_scores = lapply(1:nrow(df_bvc_plot3), function(x) scale(as.numeric(as.vector(df_bvc_plot3[x, c("bvc", paste0("X",1:1000))]))) )

# PCA
library("factoextra")
Idents(bb) = bb$subsample
subsample_avg = t(myAverageExpression(bb, slot = "data"))
test = subsample_avg[,which(colSums(subsample_avg) != 0)]
test = test[, names(sum_dif[order(sum_dif, decreasing = T)[1:200]])]
res.pca <- prcomp(test, scale = T)
groups = factor(c(rep("Behave", 19), rep("Control", 19)))
groups2 = factor(colsplit(as.vector(unique(bb$subsample)), pattern = "\\.", names = c('1', "2"))[,1])
fviz_pca_ind(res.pca, col.ind = groups, palette = c("#00AFBB",  "#FC4E07"), addEllipses = T, ellipse.type = "confidence", legend.title = "Groups", repel = T)
# fviz_pca_ind(res.pca, col.ind = groups2, addEllipses = T, ellipse.type = "confidence", legend.title = "Groups", repel = T)

pdf("C:/Users/miles/Downloads/test3.pdf", width = 8, height = 8)
print(test)
dev.off()

library(sp)
# test = mySingleGeneHull(bb, "dlx5")
ggplot(test$data, aes(UMAP_1, UMAP_2)) + geom_point(size = 0.9, colour = test$data$col, alpha = 0.9) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + stat_density2d(data = test$data[which(test$data$gene == 'nkx2-1'),], alpha = 0.5, geom="polygon", contour_var = "count", position="identity", bins = 2, size = 0.9, contour = T, show.legend = T)
ggbld <- ggplot_build(test)
gdata <- ggbld$data[[3]]
gdata %>% ggplot() + geom_point(aes(x, y, color = level))
test1 = gdata[which(gdata$level == levels(gdata$level)[2]),]
SpatialPolygons(lapply(levels(gdata$level)[2], function(x) {
  pts <- test1[test1$level == x,]
  Polygons(list(Polygon(as.matrix(data.frame(x=pts$x, y=pts$y)))), as.character(x))
})) -> polys
polys_dat <- SpatialPolygonsDataFrame(polys, 
                                      data.frame(id=sapply(slot(polys, "polygons"), slot, "ID"),
                                                 row.names=sapply(slot(polys, "polygons"), slot, "ID"),
                                                 stringsAsFactors=FALSE))
plot(polys)
# test1 = test1[order( test1$y ),]
# ggplot(test1, aes(x, y, color = level)) + geom_point()

#***********************************************************************************
# scGNN ============================================================================
#***********************************************************************************
emb = read.csv("C:/Users/miles/Downloads/CichlidDataFolder_embedding.csv")
res = read.csv("C:/Users/miles/Downloads/CichlidDataFolder_results.txt")
ggplot(emb, aes(x = embedding0, y = embedding2)) + geom_point() + theme_classic()

bb$scgnn = factor(res$Celltype, levels = sort(unique(res$Celltype)))
Idents(bb) = bb$scgnn
DimPlot(bb, label = T)

#*******************************************************************************
# Vole and Clown ===============================================================
#*******************************************************************************
vole = readRDS("~/Downloads/vole_101421.rds")
clown = readRDS("~/Downloads/anenomefish_clustered_061821.rds")

# Vole Bulk OXTR DEGs
Idents(vole) = vole$oxtr
oxtr.deg = FindMarkers(vole, ident.1 = "high", ident.2 = "low")
oxtr.deg$gene = rownames(oxtr.deg)
oxtr.deg$isSig = oxtr.deg$p_val_adj < 0.05
oxtr.deg$neg_log_p = -log10(oxtr.deg$p_val)
this_thresh = min(oxtr.deg$neg_log_p[which(oxtr.deg$isSig)])
ggplot(oxtr.deg, aes(x = avg_log2FC, y = neg_log_p, color = isSig)) + geom_point(alpha = 0.5) + geom_hline(yintercept = this_thresh, lty = 2) + geom_text_repel(data = oxtr.deg[which(oxtr.deg$isSig & (oxtr.deg$neg_log_p > 90 | abs(oxtr.deg$avg_log2FC) > 200)),], aes(label = gene)) + scale_color_manual(values = c("gray40", "#6a040f")) + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression("-"*Log["10"]*" P")) + ggtitle("Vole Oxytocin Receptor DEG Volcano Plot") + theme_bw()

# Vole Bulk Sex DEGs
Idents(vole) = vole$sex
vole.sex.deg = FindMarkers(vole, ident.1 = "f", ident.2 = "m")
vole.sex.deg$gene = rownames(vole.sex.deg)
vole.sex.deg$isSig = vole.sex.deg$p_val_adj < 0.05
vole.sex.deg$neg_log_p = -log10(vole.sex.deg$p_val)
this_thresh = min(vole.sex.deg$neg_log_p[which(vole.sex.deg$isSig)])
ggplot(vole.sex.deg, aes(x = avg_log2FC, y = neg_log_p, color = isSig)) + geom_point(alpha = 0.5) + geom_hline(yintercept = this_thresh, lty = 2) + geom_text_repel(data = vole.sex.deg[which(vole.sex.deg$isSig & (vole.sex.deg$neg_log_p > 75 | abs(vole.sex.deg$avg_log2FC) > 200)),], aes(label = gene)) + scale_color_manual(values = c("gray40", "#1a759f")) + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression("-"*Log["10"]*" P")) + ggtitle("Vole Sex DEG Volcano Plot") + theme_bw()

# Comparison of Vole Bulk DEGs
vole.deg.combined = oxtr.deg[which(oxtr.deg$gene %in% vole.sex.deg$gene),]
colnames(vole.deg.combined) = paste0(colnames(vole.deg.combined), "_oxtr")
vole.deg.combined[,paste0(colnames(vole.sex.deg), "_sex")] = vole.sex.deg[match(vole.deg.combined$gene, vole.sex.deg$gene),]
vole.deg.combined$Significant = "None"
vole.deg.combined$Significant[which(vole.deg.combined$p_val_adj_oxtr < 0.05)] = "Oxtr"
vole.deg.combined$Significant[which(vole.deg.combined$p_val_adj_sex < 0.05)] = "Sex"
vole.deg.combined$Significant[which(vole.deg.combined$p_val_adj_oxtr < 0.05 & vole.deg.combined$p_val_adj_sex < 0.05)] = "Both"
vole.deg.combined$Significant = factor(vole.deg.combined$Significant, levels = c("Both", "Oxtr", "Sex", "None"))
ggplot(vole.deg.combined, aes(x = neg_log_p_oxtr, y = neg_log_p_sex, color = Significant)) + geom_point(alpha = 0.5) + scale_color_manual(values = c("goldenrod3", "#6a040f", "#1a759f", "gray40")) + xlab(expression("Oxytocin Receptor -"*Log["10"]*" P")) + ylab(expression("Sex -"*Log["10"]*" P")) + ggtitle("Correlation of Oxtr DEGs and Sex DEGs", subtitle = paste0("R = ", cor(vole.deg.combined$neg_log_p_oxtr, vole.deg.combined$neg_log_p_sex))) + geom_text_repel(data = vole.deg.combined[which(vole.deg.combined$neg_log_p_oxtr > 100 | vole.deg.combined$neg_log_p_sex > 100),], aes(label = gene_oxtr)) + theme_bw()
ggplot(vole.deg.combined, aes(x = avg_log2FC_oxtr, y = avg_log2FC_sex, color = Significant)) + geom_point(alpha = 0.5) + scale_color_manual(values = c("goldenrod3", "#6a040f", "#1a759f", "gray40")) + xlab(expression("Oxytocin Receptor "*Log["2"]*" Fold Change")) + ylab(expression("Sex "*Log["2"]*" Fold Change")) + ggtitle("Correlation of Oxtr DEGs and Sex DEGs", subtitle = paste0("R = ", cor(vole.deg.combined$avg_log2FC_oxtr[which(is.finite(vole.deg.combined$avg_log2FC_oxtr) & is.finite(vole.deg.combined$avg_log2FC_sex))], vole.deg.combined$avg_log2FC_sex[which(is.finite(vole.deg.combined$avg_log2FC_oxtr) & is.finite(vole.deg.combined$avg_log2FC_sex))]))) + geom_text_repel(data = vole.deg.combined[which(abs(vole.deg.combined$avg_log2FC_oxtr) > 225 | abs(vole.deg.combined$avg_log2FC_sex) > 225),], aes(label = gene_oxtr)) + theme_bw()

# Violin / BoxPlot for Voles by Subsample
vgene = "Lsamp"
cond = "oxtr"
vole$vgene = vole@assays$RNA@data[vgene,]
vdf = aggregate(vgene ~ subsample + sample + sex + oxtr, data = vole@meta.data, mean)
valid_combos = expand.grid(unique(vole$sample), c(0,1,2,3))
valid_combos = paste0(valid_combos[,1], "_", valid_combos[,2])
vdf = vdf[which(vdf$subsample %in% valid_combos),]
vdf$sex = plyr::revalue(vdf$sex, replace = c("f" = "Female", "m" = "Male"))
vdf$oxtr = plyr::revalue(vdf$oxtr, replace = c("high" = "High", "low" = "Low"))
vdf$subsample = factor(vdf$subsample, levels = vdf$subsample[order(vdf[, cond])])
vdf$cond = vdf[, cond]
if (cond == "oxtr") { my_pal = c("#6a040f", "gray40") }
if (cond == "sex")  { my_pal = c("#1a759f", "gray40") }
ggplot(vdf, aes(x = cond, y = vgene, color = cond, fill = cond)) + geom_violin(alpha = 0.3) + geom_boxplot(width = 0.3, alpha =0.75) + geom_point(position = position_jitter(width = 0.1)) + ylab("Normalized Expression") + xlab("") + ggtitle(vgene) + scale_color_manual(values = my_pal, guide = 'none') + scale_fill_manual(values = my_pal, guide = 'none') + theme_bw()

# Clown Bulk Sex DEGs
Idents(clown) = clown$sex
clown.sex.deg = FindMarkers(clown, ident.1 = "f", ident.2 = "m")
clown.sex.deg$gene = rownames(clown.sex.deg)
clown.sex.deg$isSig = clown.sex.deg$p_val_adj < 0.05
clown.sex.deg$neg_log_p = -log10(clown.sex.deg$p_val)
this_thresh = min(clown.sex.deg$neg_log_p[which(clown.sex.deg$isSig)])
ggplot(clown.sex.deg, aes(x = avg_log2FC, y = neg_log_p, color = isSig)) + geom_point(alpha = 0.5) + geom_hline(yintercept = this_thresh, lty = 2) + geom_text_repel(data = clown.sex.deg[which(clown.sex.deg$isSig & (clown.sex.deg$neg_log_p > 75 | abs(clown.sex.deg$avg_log2FC) > 200)),], aes(label = gene)) + scale_color_manual(values = c("gray40", "goldenrod3")) + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression("-"*Log["10"]*" P")) + ggtitle("Clown Sex DEG Volcano Plot")

# UMAP Plot of Clown Gene
p = FeaturePlot(clown, "clu", split.by = 'subsample2', order = T, ncol = 3, by.col = T)
p_list = list(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], p[[7]], p[[8]], p[[9]], p[[10]], p[[11]])
cowplot::plot_grid(plotlist=p_list, nrow = 3)

# Violin / BoxPlot for Clown by Subsample
clown$subsample2 = paste0(clown$sample, "_", clown$subsample)
cgene = "adcy5"
clown$cgene = clown@assays$RNA@data[cgene,]
cdf = aggregate(cgene ~ subsample2 + sample + sex, data = clown@meta.data, mean)
cdf$sex = plyr::revalue(cdf$sex, replace = c("f" = "Female", "m" = "Male"))
cdf$subsample = factor(cdf$subsample, levels = cdf$subsample[order(cdf$sex)])
my_pal = c("#6a040f", "#1a759f")
ggplot(cdf, aes(x = sex, y = cgene, color = sex, fill = sex)) + geom_violin(alpha = 0.3) + geom_boxplot(width = 0.3, alpha =0.75) + geom_point(position = position_jitter(width = 0.1)) + ylab("Normalized Expression") + xlab("") + ggtitle(cgene) + scale_color_manual(values = my_pal, guide = 'none') + scale_fill_manual(values = my_pal, guide = 'none') + theme_bw()

# Clown Cluster Sex DEGs
library("scales")
clust.cols = gc.ramp <- hue_pal()(max(as.numeric(clown$seurat_clusters))+1)
names(clust.cols) = 0:max(as.numeric(clown$seurat_clusters))
clown$clust_sex = paste0(clown$seurat_clusters, "_", clown$sex)
Idents(clown) = clown$clust_sex
clown.sex.clust.deg = data.frame()
for (this_clust in 0:(max(as.numeric(clown$seurat_clusters))-1)) {
  print(this_clust)
  this_deg = FindMarkers(clown, ident.1 = paste0(this_clust, "_f"), ident.2 = paste0(this_clust, "_m"))
  this_deg$gene = rownames(this_deg)
  this_deg$cluster = this_clust
  clown.sex.clust.deg = rbind(clown.sex.clust.deg, this_deg)
}
clown.sex.clust.deg$isSig = clown.sex.clust.deg$p_val_adj < 0.05
clown.sex.clust.deg$neg_log_p = -log10(clown.sex.clust.deg$p_val)
clown.sex.clust.deg$clust_sig = as.numeric(clown.sex.clust.deg$cluster)
clown.sex.clust.deg$clust_sig[which(clown.sex.clust.deg$p_val_adj > 0.05)] = 'none'
clown.sex.clust.deg$clust_sig = factor(clown.sex.clust.deg$clust_sig, levels = c(0:max(as.numeric(clown$seurat_clusters)), 'none'))
# clown.sex.clust.deg$col = clust.cols[as.character(clown.sex.clust.deg$cluster)]
# clown.sex.clust.deg$col[which(clown.sex.clust.deg$p_val_adj > 0.05)] = 'gray60'
this_thresh = min(clown.sex.clust.deg$neg_log_p[which(clown.sex.clust.deg$isSig)])
ggplot(clown.sex.clust.deg, aes(x = avg_log2FC, y = neg_log_p, color = isSig)) + geom_point(alpha = 0.5) + geom_hline(yintercept = this_thresh, lty = 2) + geom_text_repel(data = clown.sex.clust.deg[which(clown.sex.clust.deg$isSig & (clown.sex.clust.deg$neg_log_p > 15 | abs(clown.sex.clust.deg$avg_log2FC) > 2)),], aes(label = gene)) + scale_color_manual(values = c("gray40", "goldenrod3")) + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression("-"*Log["10"]*" P")) + ggtitle("Clown Sex Cluster DEG Volcano Plot")
# ggplot(clown.sex.clust.deg, aes(x = avg_log2FC, y = neg_log_p, color = col)) + geom_point(alpha = 0.5) + geom_hline(yintercept = this_thresh, lty = 2) + geom_text_repel(data = clown.sex.clust.deg[which(clown.sex.clust.deg$isSig & (clown.sex.clust.deg$neg_log_p > 15 | abs(clown.sex.clust.deg$avg_log2FC) > 2)),], aes(label = gene)) + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression("-"*Log["10"]*" P")) + ggtitle("Clown Sex Cluster DEG Volcano Plot")
ggplot(clown.sex.clust.deg, aes(x = avg_log2FC, y = neg_log_p, color = clust_sig)) + geom_point(alpha = 0.5) + geom_hline(yintercept = this_thresh, lty = 2) + geom_text_repel(data = clown.sex.clust.deg[which(clown.sex.clust.deg$isSig & (clown.sex.clust.deg$neg_log_p > 15 | abs(clown.sex.clust.deg$avg_log2FC) > 2)),], aes(label = gene)) + scale_color_manual(values = clust.cols) + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression("-"*Log["10"]*" P")) + ggtitle("Clown Sex Cluster DEG Volcano Plot") 

ss_means = list()
for (ss in sort(unique(clown$subsample2))) {
  ss_means[[ss]] = rowMeans(clown@assays$RNA@data[,which(clown$subsample2 == ss)])
}
ss_means = as.data.frame(ss_means)
sep_cgenes = sapply(1:nrow(ss_means), function(x) all(outer(X = as.numeric(ss_means[x,1:6]), Y = as.numeric(ss_means[x,7:12]), FUN = '>')) | all(outer(X = as.numeric(ss_means[x,1:6]), Y = as.numeric(ss_means[x,7:12]), FUN = '<')) )
up_cgenes = sapply(rownames(ss_means)[which(sep_cgenes)], function(x) ss_means[x, 1] > ss_means[x, 7] )
big_cgenes = sapply(rownames(ss_means)[which(sep_cgenes)], function(x) abs(sum(ss_means[x, 1:6]) - sum(ss_means[x, 7:12])) )
names(big_cgenes) = rownames(ss_means)[which(sep_cgenes)]
# big_cgenes = sort(big_cgenes, decreasing = T)
sd_cgenes = sapply(rownames(ss_means)[which(sep_cgenes)], function(x) sum(sd(ss_means[x, 1:6]), sd(ss_means[x, 7:12])) )
names(sd_cgenes) = rownames(ss_means)[which(sep_cgenes)]
cgenes_df = data.frame(gene = names(big_cgenes), mean_dif = big_cgenes, sd_sum = sd_cgenes, up_in_f = up_cgenes)
cgenes_df$up_in_f = plyr::revalue(as.character(cgenes_df$up_in_f), replace = c("TRUE" = "Female", "FALSE" = "Male"))
ggplot(cgenes_df, aes(x = mean_dif, y = sd_sum, color = mean_dif)) + geom_point() + geom_text_repel(data = cgenes_df[which(cgenes_df$mean_dif > 1),], aes(label = gene), color = "black") + scale_color_viridis_c() + xlab("Difference in Female Mean of Avg Expression - Male Mean of Avg Expression") + ylab("Sum of Female and Male Standard Deviation") + theme_bw()

clust_sub_num = aggregate(nCount_RNA ~ subsample2 + seurat_clusters, data = clown@meta.data, length)
clust_big = sapply(0:max(as.numeric(clown$seurat_clusters)), function(x) all(clust_sub_num$nCount_RNA[which(clust_sub_num$seurat_clusters == x)] > 10) )
clust_big = seq(0, max(as.numeric(clown$seurat_clusters)))[which(clust_big)]

clust_sep_cgenes = c()
clust_cgenes_df = data.frame()
for (this_clust in clust_big) {
  print(this_clust)
  ss_means = list()
  for (ss in sort(unique(clown$subsample2))) {
    ss_means[[ss]] = rowMeans(clown@assays$RNA@data[,which(clown$subsample2 == ss & clown$seurat_clusters == this_clust)])
  }
  ss_means = as.data.frame(ss_means)
  this_sep_cgenes = sapply(1:nrow(ss_means), function(x) all(outer(X = as.numeric(ss_means[x,1:6]), Y = as.numeric(ss_means[x,7:12]), FUN = '>')) | all(outer(X = as.numeric(ss_means[x,1:6]), Y = as.numeric(ss_means[x,7:12]), FUN = '<')) )
  this_up_cgenes = sapply(rownames(ss_means)[which(this_sep_cgenes)], function(x) ss_means[x, 1] > ss_means[x, 7] )
  this_big_cgenes = sapply(rownames(ss_means)[which(this_sep_cgenes)], function(x) abs(sum(ss_means[x, 1:6]) - sum(ss_means[x, 7:12])) )
  names(this_big_cgenes) = rownames(ss_means)[which(this_sep_cgenes)]
  sd_cgenes = sapply(rownames(ss_means)[which(this_sep_cgenes)], function(x) sum(sd(ss_means[x, 1:6]), sd(ss_means[x, 7:12])) )
  names(sd_cgenes) = rownames(ss_means)[which(this_sep_cgenes)]
  this_cgenes_df = data.frame(gene = names(this_big_cgenes), mean_dif = this_big_cgenes, sd_sum = sd_cgenes, up_in_f = this_up_cgenes, cluster = rep(this_clust, length(this_big_cgenes)))
  clust_cgenes_df = rbind(clust_cgenes_df, this_cgenes_df)
}
clust_cgenes_df$human = clown_coverter$human[match(clust_cgenes_df$gene, clown_coverter$cgene)]
clust_cgenes_df$cluster = factor(clust_cgenes_df$cluster, levels = 0:max(as.numeric(clown$seurat_clusters)))
ggplot(clust_cgenes_df, aes(x = mean_dif, y = sd_sum, color = up_in_f)) + geom_point(alpha = 0.75) + geom_text_repel(data = clust_cgenes_df[which(clust_cgenes_df$mean_dif > 3.4),], aes(label = gene), color = "black") + scale_color_manual(values = c("#6a040f", "#1a759f")) + xlab("Difference in Female Mean of Avg Expression - Male Mean of Avg Expression") + ylab("Sum of Female and Male Standard Deviation") + theme_bw() + theme(legend.title = element_blank())
ggplot(clust_cgenes_df, aes(x = mean_dif, y = sd_sum, color = cluster)) + geom_point(alpha = 0.75) + geom_text_repel(data = clust_cgenes_df[which(clust_cgenes_df$mean_dif > 3.4),], aes(label = gene), color = "black") + scale_color_manual(values = clust.cols)              + xlab("Difference in Female Mean of Avg Expression - Male Mean of Avg Expression") + ylab("Sum of Female and Male Standard Deviation") + theme_bw() + theme(legend.title = element_blank())
write.csv(clust_cgenes_df, "~/research/brain/results/clown/clust_sex_segregating.csv")

# Number of Sex Segregating Genes per Cluster
ggplot(as.data.frame(table(clust_cgenes_df$cluster)), aes(Var1, Freq, color = Var1, fill = Var1)) + geom_bar(stat = "identity") + scale_color_manual(values = clust.cols) + scale_fill_manual(values = clust.cols) + theme_bw() + xlab("Cluster") + ylab("Number of Sex DEGs") + scale_y_continuous(expand = c(0,0)) + theme(legend.title=element_blank())


# Clown Converter
censembl = useEnsembl("ensembl", mirror = "uswest", dataset = "apercula_gene_ensembl")
human =  useEnsembl("ensembl", mirror = "uswest", dataset = "hsapiens_gene_ensembl")
ex_gene_converter = getLDS(attributes = c("external_gene_name", 'description'), filters = "external_gene_name", values = rownames(clown) , mart = censembl, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
ens_converter = getLDS(attributes = c("ensembl_gene_id", 'description'), filters = "ensembl_gene_id", values = rownames(clown) , mart = censembl, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
colnames(ens_converter)[1] = "Gene.name"
clown_coverter = rbind(ex_gene_converter, ens_converter)
clown_coverter = data.frame(cgene = rownames(clown), human = clown_coverter$HGNC.symbol[match(rownames(clown), clown_coverter$Gene.name)], description = clown_coverter$Gene.description[match(rownames(clown), clown_coverter$Gene.name)])
write.csv(clown_coverter, "~/research/brain/results/clown/clown_gene_converter.csv")

# Number of Sex DEGs per Cluster
ggplot(as.data.frame(table(clown.sex.clust.deg$cluster[which(clown.sex.clust.deg$isSig)])), aes(Var1, Freq, color = Var1, fill = Var1)) + geom_bar(stat = "identity") + scale_color_manual(values = clust.cols) + scale_fill_manual(values = clust.cols) + theme_bw() + xlab("Cluster") + ylab("Number of Sex DEGs") + scale_y_continuous(expand = c(0,0)) + theme(legend.title=element_blank())

# Violin / BoxPlot for Clown by Cluster Subsample
clown$subsample2 = paste0(clown$sample, "_", clown$subsample)
cgene = "whrna"
this_clust = 2
clown$cgene = clown@assays$RNA@data[cgene,which(clown$seurat_clusters == this_clust)]
cdf = aggregate(cgene ~ subsample2 + sample + sex, data = clown@meta.data, mean)
cdf$sex = plyr::revalue(cdf$sex, replace = c("f" = "Female", "m" = "Male"))
cdf$subsample = factor(cdf$subsample, levels = cdf$subsample[order(cdf$sex)])
my_pal = c("#6a040f", "#1a759f")
ggplot(cdf, aes(x = sex, y = cgene, color = sex, fill = sex)) + geom_violin(alpha = 0.3) + geom_boxplot(width = 0.3, alpha =0.75) + geom_point(position = position_jitter(width = 0.1)) + ylab("Normalized Expression") + xlab("") + ggtitle(paste(cgene, "in Cluster", this_clust)) + scale_color_manual(values = my_pal, guide = 'none') + scale_fill_manual(values = my_pal, guide = 'none') + theme_bw()

clown_meta = clown@meta.data
ggplot(clown_meta, aes(x = subsample2, y = nCount_RNA, color = subsample2, fill = subsample2)) + geom_violin(alpha = 0.4) + geom_boxplot(width = 0.2, alpha = 0.2) + xlab("") + ylab("Number of UMIs")
ggplot(clown_meta, aes(x = subsample2, y = nFeature_RNA, color = subsample2, fill = subsample2)) + geom_violin(alpha = 0.4) + geom_boxplot(width = 0.2, alpha = 0.2) + xlab("") + ylab("Number of Features")
ggplot(clown_meta, aes(x = nCount_RNA, y = nFeature_RNA, color = subsample2)) + geom_point()
FeaturePlot(clown, 'nCount_RNA', order = T, label = T)
clown_meta_agr = aggregate(nCount_RNA ~ seurat_clusters + subsample2, data = clown_meta, length)
clown_meta_agr$sex = startsWith(clown_meta_agr$subsample2, "f")
clown_meta_agr$sex = plyr::revalue(as.character(clown_meta_agr$sex), replace = c("TRUE" = "Female", "FALSE" = "Male"))
ggplot(clown_meta_agr, aes(x = seurat_clusters, y = nCount_RNA, fill = sex, color = sex)) + geom_boxplot(alpha = 0.3) + geom_point(position = position_jitterdodge(), width = 0.05, alpha = 0.6) + scale_color_manual(values = my_pal) + scale_fill_manual(values = my_pal) + ylab("Number of Cells")

#*******************************************************************************
# Brianna Markers ==============================================================
#*******************************************************************************
# *** BRIANNA 15 ***
brianna15 = xlsx::read.xlsx("~/Downloads/heatmapmarkerlist_george1.xlsx", sheetIndex = 1, startRow = 1)
colnames(brianna15)[which(colnames(brianna15) == "Gene.")] = "Gene"
brianna15$Category = str_replace_all(trimws(brianna15$Category, which = "both"), "[^[:alnum:]\\s]", "")
brianna15$Category[which( startsWith(brianna15$Category, "Neuroanat") )]  = "NeuroanatNeurodev TF"
brianna15$LOCID = str_replace_all(trimws(brianna15$LOCID, which = "both"), "[^[:alnum:]\\s]", "")
brianna15$gene_name = gtf$gene_name[match(brianna15$LOCID, gtf$loc)]
brianna15$col = plyr::revalue(brianna15$Category, replace = c("Neuromodulator" = "#00E7EC", "Neuromodulatory Receptor" = "#FDD615", "NeuroanatNeurodev TF" = "#FE04FF"))

# Create a dataframe that has all the combos of clusters and genes
all_combos = expand.grid(unique(brianna15$gene_name), convert15$new.full)
colnames(all_combos) = c("gene_name", "cluster")
all_combos[, colnames(brianna15)[which(! colnames(brianna15) %in% colnames(all_combos))]] = brianna15[match(all_combos$gene_name, brianna15$gene_name), colnames(brianna15)[which(! colnames(brianna15) %in% colnames(all_combos))]] 
all_combos[, c("avg_logFC", "pct.1", "pct.2", "pct_dif", "num.1", "num.2")] = 0

# Find the expression levels of Brianna's genes
for (cluster in convert15$new.full) {
  print(cluster)
  this_idx = which(all_combos$cluster == cluster)
  this_cluster_cells = colnames(bb)[which(bb$seuratclusters15 == convert15$old[which(convert15$new.full == cluster)])]
  clust_pct_fc = pct_dif_avg_logFC(bb, this_cluster_cells, colnames(bb)[which(! colnames(bb) %in% this_cluster_cells)], features = brianna15$gene_name)
  all_combos[this_idx, c("avg_logFC", "pct.1", "pct.2", "pct_dif", "num.1", "num.2")] = clust_pct_fc[match(all_combos$gene_name[this_idx], clust_pct_fc$genes), c("avg_logFC", "pct.1", "pct.2", "pct_dif", "num.1", "num.2")]
}

# Transparancy theme
all_combos$col4 = "gray98"
all_combos$col4[which(all_combos$pct.1 >= 5 )]  = paste0(all_combos$col[which(all_combos$pct.1 >= 5 )],  "70")
all_combos$col4[which(all_combos$pct.1 >= 10 )] = paste0(all_combos$col[which(all_combos$pct.1 >= 10 )], "80")
all_combos$col4[which(all_combos$pct.1 >= 20 )] = paste0(all_combos$col[which(all_combos$pct.1 >= 20 )], "90")
all_combos$col4[which(all_combos$pct.1 >= 30 )] = paste0(all_combos$col[which(all_combos$pct.1 >= 30 )], "ff")

# Transparancy theme continuous
all_combos$col5 = "gray98"
OldRange = (25 - 5)  
NewRange = (100 - 70)  
some.alpha.values = floor( (((all_combos$pct.1[which(all_combos$pct.1 >= 5)] - 5) * NewRange) / OldRange) + 70 )
some.alpha.values[which(some.alpha.values >= 100)] = 100
some.alpha.values = as.character(some.alpha.values)
some.alpha.values[which(some.alpha.values == "100")] = "ff"
all_combos$col5[which(all_combos$pct.1 >= 5 )]  = paste0(all_combos$col[which(all_combos$pct.1 >= 5 )], some.alpha.values)

# Plot the matrix
all_combos15 = all_combos
brianna_order = xlsx::read.xlsx("~/Downloads/heatmapmarkerlist_george2.xlsx", sheetIndex = 1, startRow = 1)
all_combos$Gene = factor(all_combos$Gene, levels = brianna_order$Gene)
all_combos$cluster = factor(all_combos$cluster, levels = convert15$new.full[order(as.numeric(convert15$new.num), decreasing = T)])
pdf("~/research/brain/results/bri15_markers_heatmap_8.pdf", height = 3.5, width = 12)
ggplot(all_combos[which(! is.na(all_combos$Gene)),], aes(x = Gene, y = cluster, fill = col4)) + geom_tile(color = "gray40") + scale_fill_identity() + coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic")) + xlab("") + ylab("") + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0, 0))
dev.off()

# *** BRIANNA 53 ***
brianna53 = xlsx::read.xlsx("~/Downloads/53heatmapmarkerlist_george6.xlsx", sheetIndex = 1, startRow = 1)
colnames(brianna53)[1] = "Category"
colnames(brianna53)[4] = "Rank.Within"
colnames(brianna53)[5] = "Rank.Overall"
brianna53[, 6:ncol(brianna53)] = NULL
colnames(brianna53)[which(colnames(brianna53) == "Gene.")] = "Gene"
brianna53$Category = str_replace_all(trimws(brianna53$Category, which = "both"), "[^[:alnum:]\\s]", "")
brianna53$Category[which( startsWith(brianna53$Category, "Neuroanat") )]  = "NeuroanatNeurodev TF"
brianna53$LOCID = str_replace_all(trimws(brianna53$LOCID, which = "both"), "[^[:alnum:]\\s]", "")
brianna53$gene_name = gtf$gene_name[match(brianna53$LOCID, gtf$loc)]
brianna53$col = plyr::revalue(brianna53$Category, replace = c("Neuromodulator" = "#00E7EC", "Neuromodulatory Receptor" = "#FDD615", "NeuroanatNeurodev TF" = "#FE04FF"))

# Last minute changes the gene order of just a few genes
drd1_idx = which(brianna53$Gene == "drd1")
nr4_idx = which(brianna53$Gene == "nr4a2b (nurr1)")
neurod1_idx = which(brianna53$Gene == "neurod1")
neurod6b_idx = which(brianna53$Gene == "neurod6b")
cck_idx = which(brianna53$Gene == "cck")
nos1_idx = which(brianna53$Gene == "nos1")
brianna53_idx = c(1:(neurod6b_idx-1), neurod1_idx, neurod6b_idx, (neurod1_idx+1):(cck_idx-1), nos1_idx, cck_idx, (cck_idx+1):(drd1_idx-1), (drd1_idx+1):(nr4_idx-1), drd1_idx, nr4_idx:nrow(brianna53))

brianna53[which( grepl("LOC101487266", brianna53$LOCID) ), c("LOCID", "Gene", "gene_name")] = c("LOC101487266", "nr4a2b (nurr1)", "LOC101487266")

# Create a dataframe that has all the combos of clusters and genes
all_combos = expand.grid(unique(brianna53$gene_name), convert53$new)
colnames(all_combos) = c("gene_name", "cluster")
all_combos[, colnames(brianna53)[which(! colnames(brianna53) %in% colnames(all_combos))]] = brianna53[match(all_combos$gene_name, brianna53$gene_name), colnames(brianna53)[which(! colnames(brianna53) %in% colnames(all_combos))]] 
all_combos[, c("avg_logFC", "pct.1", "pct.2", "pct_dif", "num.1", "num.2")] = 0

# Find the expression levels of Brianna's genes
for (cluster in convert53$new) {
  print(cluster)
  this_idx = which(all_combos$cluster == cluster)
  # this_cluster_cells = colnames(bb)[which(bb$seuratclusters15 == convert15$old[which(convert15$new.full == cluster)])]
  this_cluster_cells = colnames(bb)[which(bb$seuratclusters53 == convert53$old[which(convert53$new == cluster)])]
  clust_pct_fc = pct_dif_avg_logFC(bb, this_cluster_cells, colnames(bb)[which(! colnames(bb) %in% this_cluster_cells)], features = brianna53$gene_name)
  all_combos[this_idx, c("avg_logFC", "pct.1", "pct.2", "pct_dif", "num.1", "num.2")] = clust_pct_fc[match(all_combos$gene_name[this_idx], clust_pct_fc$genes), c("avg_logFC", "pct.1", "pct.2", "pct_dif", "num.1", "num.2")]
}

# Transparancy theme
all_combos$col4 = "gray98"
all_combos$col4[which(all_combos$pct.1 >= 5 )]  = paste0(all_combos$col[which(all_combos$pct.1 >= 5 )],  "70")
all_combos$col4[which(all_combos$pct.1 >= 10 )] = paste0(all_combos$col[which(all_combos$pct.1 >= 10 )], "80")
all_combos$col4[which(all_combos$pct.1 >= 20 )] = paste0(all_combos$col[which(all_combos$pct.1 >= 20 )], "90")
all_combos$col4[which(all_combos$pct.1 >= 30 )] = paste0(all_combos$col[which(all_combos$pct.1 >= 30 )], "ff")

# Transparancy theme continuous
all_combos$col5 = "gray98"
OldRange = (25 - 5)  
NewRange = (100 - 70)  
some.alpha.values = floor( (((all_combos$pct.1[which(all_combos$pct.1 >= 5)] - 5) * NewRange) / OldRange) + 70 )
some.alpha.values[which(some.alpha.values >= 100)] = 100
some.alpha.values = as.character(some.alpha.values)
some.alpha.values[which(some.alpha.values == "100")] = "ff"
all_combos$col5[which(all_combos$pct.1 >= 5 )]  = paste0(all_combos$col[which(all_combos$pct.1 >= 5 )], some.alpha.values)

# Plot the matrix
all_combos$cluster = factor(all_combos$cluster, levels = rev(unique(convert53$new)))
all_combos$Gene = factor(all_combos$Gene, levels = unique(brianna53$Gene[brianna53_idx]))
test = acast(Gene ~ cluster, data = all_combos, value.var = 'pct.1')
pheatmap::pheatmap(test)
pdf("~/research/brain/results/bri53_markers_heatmap14.pdf", height = 10, width = 12)
ggplot(all_combos[which(! is.na(all_combos$Gene)),], aes(x = Gene, y = cluster, fill = col4)) + geom_tile(color = "gray40") + scale_fill_identity() + coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic")) + xlab("") + ylab("") + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0, 0))
dev.off()
all_combos53 = all_combos

#*** Brianna All ***
all_combos15$Rank.Within = 1:nrow(all_combos15)
all_combos15$level = "15"
all_combos53$level = "53"
ball_common_genes = brianna15$LOCID[which(brianna15$LOCID %in% brianna53$LOCID)]
ball = rbind(all_combos15[which( all_combos15$LOCID %in% ball_common_genes ),], all_combos53[which( all_combos53$LOCID %in% ball_common_genes ),])

# Zack last minute changes
other = read.csv("~/Downloads/goi_less_10k_122821_by_cat_hgnc.csv")
other = other[which(other$cat == "other"),]
ball$col6 = ball$col5
ball$col6[which( ball$gene_name %in% other$mzebra & startsWith(ball$col5, "#")  )] = paste0("#002DD1", substr(ball$col6[which( ball$gene_name %in% other$mzebra & startsWith(ball$col5, "#")  )], 8, 9))
new.order = unique(brianna53$Gene[brianna53_idx])
new.order = c(new.order[which(!new.order %in% ball$Gene[which(ball$gene_name %in% other$mzebra)])], new.order[which(new.order %in% ball$Gene[which(ball$gene_name %in% other$mzebra)])])

ball$Gene[which(ball$LOCID == "LOC101464700")] = "emx3 (emx1)"
ball$Gene = factor(ball$Gene, levels = new.order)
ball$cluster = factor(ball$cluster, levels = rev(convert_all$cluster))
pdf("~/research/brain/results/bri_all_markers_heatmap4.pdf", height = 12, width = 12)
ggplot(ball, aes(x = Gene, y = cluster, fill = col6)) + geom_tile(color = "gray40") + scale_fill_identity() + coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic", size = 10), axis.text.y = element_text(size = ifelse(rev(convert_all$level) == "primary", 10, 8), face = ifelse(rev(convert_all$level) == "primary", "bold", "plain"), color = rev(convert_all$color))) + xlab("") + ylab("") + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0, 0))
dev.off()

ball.vln = unique(brianna15[,c("LOCID", "Gene", "gene_name")])
# ball.vln = unique(all_combos15$Gene[which(! is.na(all_combos15$Gene) ),c("LOCID", "Gene")])
ball.vln$Gene = factor(ball.vln$Gene, levels = brianna_order$Gene)
pdf("~/research/brain/results/supplement/vlnplot_bri.pdf", height = 17, width = 8)
StackedVlnPlot(obj = bb, features = ball.vln$gene_name, xcols = convert15$col, sec.axis.names = ball.vln$Gene )
dev.off()

# Read in the GTF (necessary for converting LOC to NCBI gene name)
gtf = read.delim("~/research/all_research/GCF_000238955.4_M_zebra_UMD2a_genomic.gtf", header = F)
gtf = gtf[which(gtf$V3 == "gene"),]
# gtf = gtf[5:nrow(gtf),]
# gtf = gtf[which(gtf$V1 != "###"),]
gtf$gene_name = unlist(str_split(as.character(gtf$V9),';', n = 2))[c(TRUE, FALSE)]
gtf$gene_name = unlist(str_split(as.character(gtf$gene_name),' ', n = 2))[c(FALSE, TRUE)]
gtf$loc = unlist(str_split(as.character(gtf$V9),'db_xref'))[c(FALSE, TRUE)]
gtf$loc = unlist(str_split(as.character(gtf$loc),';', n = 2))[c(TRUE, FALSE)]
gtf$loc = paste0("LOC", unlist(str_split(as.character(gtf$loc),'GeneID:', n = 2))[c(FALSE, TRUE)])
gtf$loc = trimws(gtf$loc, which = "both")
gene_info = read.table(paste0("~/research/all_research/gene_info.txt"), sep="\t", header = T, stringsAsFactors = F)
brainna_gene_key = xlsx::read.xlsx("~/Downloads/brain_gene_marker_key_draft2.xlsx", sheetIndex = 3, endRow = 89)
brianna15 = xlsx::read.xlsx("~/Downloads/brain_gene_marker_key_draft2.xlsx", sheetIndex = 1, startRow = 4, endRow = 19)
bb$good15_names = factor(convert15$new.full[match(bb$seuratclusters15, convert15$old)], levels = convert15$new.full)
Idents(bb) = bb$good15_names
gene15_df = data.frame(unlist(brianna15[,2:5]), c(rep('id', 15), rep('anat', 15), rep('mod', 15), rep('modr', 15)), rep(convert15$new.full, 4))
gene15_df_list = lapply(1:nrow(gene15_df), function(x) { tmp = unlist(strsplit(gene15_df[x,1], ",")); tmp_df = data.frame(bname = tmp, cat = rep(gene15_df[x,2], length(tmp)), cluster = rep(gene15_df[x,3], length(tmp))); return(tmp_df)} )
gene15_df <- do.call("rbind", gene15_df_list)
gene15_df$bname = str_replace(gene15_df$bname, "\\*", "")
gene15_df$bname = str_replace_all(trimws(gene15_df$bname, which = "both"), "[^[:alnum:]]", "")
# gene15_df = gene15_df[which( !is.na(gene15_df$bname) & ! duplicated(gene15_df$bname) ),]
gene15_df = gene15_df[which( !is.na(gene15_df$bname) ),]

# Special Cases
gene15_df$bname[which(gene15_df$bname == "lhx2b")] = "lhx2"
gene15_df$bname[which(gene15_df$bname == "dlx5")] = "dlx5a"
gene15_df$bname[which(gene15_df$bname == "pax6b")] = "pax6"
gene15_df$bname[which(gene15_df$bname == "zic2")] = "zic2a"
gene15_df$bname[which(gene15_df$bname == "her42")] = "her4.2"
gene15_df$bname[which(gene15_df$bname == "nkx21")] = "nkx2-1"
gene15_df$bname[which(gene15_df$bname == "gpr54")] = "GPR54"  # same LOC as slc18a2
gene15_df$bname[which(gene15_df$bname == "crhbp")] = "crh-bp"
gene15_df = gene15_df[which(gene15_df$bname != "lhx1"),] # not in brianna's key!

# Get LOC
gene15_df$loc1 = brainna_gene_key$M..zebra.ID[match(gene15_df$bname, brainna_gene_key$Gene)]
gene15_df$loc1 = str_replace_all(trimws(gene15_df$loc1, which = "both"), "[^[:alnum:]]", "")
gene15_df$loc2 = brainna_gene_key$M..zebra.ID[match(gene15_df$bname, brainna_gene_key$Alias.es.)]
gene15_df$loc2 = str_replace_all(trimws(gene15_df$loc2, which = "both"), "[^[:alnum:]]", "")
gene15_df$loc_final = gene15_df$loc1
gene15_df$loc_final[which(is.na(gene15_df$loc1))] = gene15_df$loc2[which(is.na(gene15_df$loc1))]
gene15_df$clean = gtf$gene_name[match(gene15_df$loc_final, gtf$loc)]
gene15_df$use = ! is.na(gene15_df$clean) & ! duplicated(gene15_df$clean)
gene15_df$col = plyr::revalue(gene15_df$cat, replace = c("id" = "#00E7EC", "anat" = "#FDD615", "mod" = "#FE04FF", "modr" = "#002DD1"))
# DotPlot(bb, features = gene15_df$clean[which(gene15_df$use)]) + scale_x_discrete(labels = gene15_df$bname[which(gene15_df$use)]) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = gene15_df$col[which(gene15_df$use)], face = "italic")) + ylab("")

all_combos = expand.grid(unique(gene15_df$clean), unique(gene15_df$cluster))
colnames(all_combos) = c("clean", "cluster")
all_combos[, colnames(gene15_df)[which(! colnames(gene15_df) %in% colnames(all_combos))]] = "na"
all_combos = all_combos[,colnames(gene15_df)]
all_combos$bname = gene15_df$bname[match(all_combos$clean, gene15_df$clean)]
gene15_df$hit = T
all_combos$hit = F
gene15_df = rbind(gene15_df, all_combos[which(! paste0(all_combos$clean, all_combos$cluster) %in%  paste0(gene15_df$clean, gene15_df$cluster) ),])
# gene15_df$hgnc = gene_info$human[match(gene15_df$clean, gene_info$mzebra)]
# gene15_df$hgnc[which(gene15_df$clean == "LOC101483038")] = "SOX10"
gene15_df$col[which(gene15_df$col == "na")] = "gray95"
gene15_df$col2 = "gray95"
gene15_df$col2[which(gene15_df$hit)] = convert15$col[match(gene15_df$cluster[which(gene15_df$hit)], convert15$new.full)]
gene15_df$cluster = factor(gene15_df$cluster, levels = convert15$new.full)
gene15_df = gene15_df[which(gene15_df$bname != "GPR54"),]
gene15_df$label = factor(tolower(gene15_df$bname), levels = unique(tolower(gene15_df$bname[order(gene15_df$cat, gene15_df$bname)])))
gene15_df[, c("avg_logFC", "pct.1", "pct.2", "pct_dif", "num.1", "num.2")] = 0

for (cluster in convert15$new.full) {
  print(cluster)
  this_idx = which(gene15_df$cluster == cluster)
  this_cluster_cells = colnames(bb)[which(bb$seuratclusters15 == convert15$old[which(convert15$new.full == cluster)])]
  clust_pct_fc = pct_dif_avg_logFC(bb, this_cluster_cells, colnames(bb)[which(! colnames(bb) %in% this_cluster_cells)])
  gene15_df[this_idx, c("avg_logFC", "pct.1", "pct.2", "pct_dif", "num.1", "num.2")] = clust_pct_fc[match(gene15_df$clean[this_idx], clust_pct_fc$genes), c("avg_logFC", "pct.1", "pct.2", "pct_dif", "num.1", "num.2")]
}

gene15_df$col3 = gene15_df$col
gene15_df$col3[which(gene15_df$use != "na" & gene15_df$pct.1 >= 5 )]  = paste0(gene15_df$col[which(gene15_df$use != "na" & gene15_df$pct.1 >= 5 )],  "40")
gene15_df$col3[which(gene15_df$use != "na" & gene15_df$pct.1 >= 25 )] = paste0(gene15_df$col[which(gene15_df$use != "na" & gene15_df$pct.1 >= 25 )], "60")
gene15_df$col3[which(gene15_df$use != "na" & gene15_df$pct.1 >= 50 )] = paste0(gene15_df$col[which(gene15_df$use != "na" & gene15_df$pct.1 >= 50 )], "80")
gene15_df$col3[which(gene15_df$use != "na" & gene15_df$pct.1 >= 75 )] = paste0(gene15_df$col[which(gene15_df$use != "na" & gene15_df$pct.1 >= 75 )], "ff")

gene15_df$col4 = gene15_df$col
gene15_df$col4[which(gene15_df$use != "na" & gene15_df$pct.1 >= 5 )]  = paste0(gene15_df$col[which(gene15_df$use != "na" & gene15_df$pct.1 >= 5 )],  "40")
gene15_df$col4[which(gene15_df$use != "na" & gene15_df$pct.1 >= 20 )] = paste0(gene15_df$col[which(gene15_df$use != "na" & gene15_df$pct.1 >= 20 )], "60")
gene15_df$col4[which(gene15_df$use != "na" & gene15_df$pct.1 >= 40 )] = paste0(gene15_df$col[which(gene15_df$use != "na" & gene15_df$pct.1 >= 40 )], "80")
gene15_df$col4[which(gene15_df$use != "na" & gene15_df$pct.1 >= 60 )] = paste0(gene15_df$col[which(gene15_df$use != "na" & gene15_df$pct.1 >= 60 )], "ff")


pdf("~/research/brain/results/bri15_markers_heatmap_6.pdf", height = 3.1, width = 12)
ggplot(gene15_df[which(! is.na(gene15_df$clean) ),], aes(x = label, y = cluster, fill = col4)) + geom_tile(color = "gray40") + scale_fill_identity() + coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic")) + xlab("") + ylab("") + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0, 0))
dev.off()

pdf("~/research/brain/results/bri15_markers_heatmap_5.pdf", height = 3.1, width = 12)
ggplot(gene15_df[which(! is.na(gene15_df$clean) ),], aes(y = label, x = cluster, fill = col3)) + geom_tile(color = "gray40") + scale_fill_identity() + coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic")) + xlab("") + ylab("") + coord_flip() + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0, 0))
dev.off()

pdf("~/research/brain/results/bri15_markers_heatmap_4.pdf", height = 3.1, width = 12)
ggplot(gene15_df[which(! is.na(gene15_df$clean) ),], aes(y = label, x = cluster, fill = col)) + geom_tile(color = "gray40") + scale_fill_identity() + coord_fixed() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("") + coord_flip() + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0, 0))
dev.off()

pdf("C:/Users/miles/Downloads/bri15_markers_heatmap_3.pdf", height = 13, width = 6)
ggplot(gene15_df[which(! is.na(gene15_df$clean) ),], aes(y = label, x = cluster, fill = col)) + geom_tile(color = "gray40") + scale_fill_identity() + coord_fixed() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(face = "italic")) + xlab("") + ylab("")
dev.off()

pdf("C:/Users/miles/Downloads/bri15_markers_heatmap_2.pdf", height = 13, width = 6)
ggplot(gene15_df[which(! is.na(gene15_df$clean) ),], aes(y = label, x = cluster, fill = col2)) + geom_tile(color = "gray40") + scale_fill_identity() + coord_fixed() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(face = "italic")) + xlab("") + ylab("")
dev.off()

pdf("C:/Users/miles/Downloads/bri15_markers_heatmap_1.pdf", height = 13, width = 6)
ggplot(gene15_df[which(! is.na(gene15_df$clean) ),], aes(y = label, x = cluster, fill = hit)) + geom_tile(color = "black") + coord_fixed() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(face = "italic")) + xlab("") + ylab("") + scale_fill_manual(values = c(viridis(3)[2], "gold"), guide = 'none')
dev.off()

BriDotPlot(gene15_df$clean[which(gene15_df$use)])

BriDotPlot = function(features) {
  exp_mat = myAverageExpression(bb, features = features)
  scale_exp_mat = scale(t(exp_mat))
  dp_df = melt(scale_exp_mat)
  # dp_df = aggregate(exp ~ feature, exp_df, mean)
  colnames(dp_df) = c("Cluster", "Feature", "Expression")
  ggplot(dp_df, aes(x = Feature, y = Cluster, fill = Expression)) + geom_point() + scale_fill_gradientn(colors = c("lightgrey", "blue")) + theme_classic()
}

#****************************************************************************
# scGNN =====================================================================
#****************************************************************************
library("data.table")
imp_mat = read.csv("C:/Users/miles/Downloads/CichlidDataFolder_recon.csv")
t_imp_mat <- transpose(imp_mat)
t_imp_mat[,c("seuratclusters15", "seuratclusters53", "subsample", "sample")] = bb@meta.data[,c("seuratclusters15", "seuratclusters53", "subsample", "sample")]
colnames(t_imp_mat) = t_imp_mat[1,]
t_imp_mat = t_imp_mat[-c(1),]
t_imp_mat2 <- mutate_all(t_imp_mat, function(x) as.numeric(as.character(x)))
t_imp_mat2[,c("seuratclusters15", "seuratclusters53", "subsample", "sample")] = bb@meta.data[,c("seuratclusters15", "seuratclusters53", "subsample", "sample")]
ggplot(t_imp_mat2, aes(satb1, mettl15, color = subsample)) + geom_point() + NoLegend()
ggplot(t_imp_mat2, aes(x = subsample, y = ccdc71, color = subsample)) + geom_boxplot() + geom_point(alpha = 0.01, position = position_jitterdodge()) + NoLegend()
# bb <- SetAssayData(object = bb, slot = "imp", new.data = imp_mat2)

saveRDS(t_imp_mat2, "C:/Users/miles/Downloads/brain/data/scgnn_imputed.rds")

all_gene_df = aggregate(. ~ subsample, t_imp_mat2[,c(colnames(t_imp_mat2)[1:2000], "subsample")], mean)
sep_cgenes = sapply(2:ncol(all_gene_df), function(x) all(outer(X = as.numeric(all_gene_df[1:19,x]), Y = as.numeric(all_gene_df[20:38,x]), FUN = '>')) | all(outer(X = as.numeric(all_gene_df[1:19,x]), Y = as.numeric(all_gene_df[20:38,x]), FUN = '<')) )
up_cgenes = sapply(colnames(all_gene_df)[which(sep_cgenes) + 1], function(x) all_gene_df[1,x] > all_gene_df[20,x] )
big_cgenes = sapply(colnames(all_gene_df)[which(sep_cgenes) + 1], function(x) abs(sum(all_gene_df[1:19,x]) - sum(all_gene_df[20:38,x])) )
names(big_cgenes) = colnames(all_gene_df)[which(sep_cgenes) + 1]
sd_cgenes = sapply(colnames(all_gene_df)[which(sep_cgenes) + 1], function(x) sum(sd(all_gene_df[1:19,x]), sd(all_gene_df[20:38,x])) )
names(sd_cgenes) = colnames(all_gene_df)[which(sep_cgenes) + 1]
cgenes_df = data.frame(gene = names(big_cgenes), mean_dif = big_cgenes, sd_sum = sd_cgenes, up_in_b = up_cgenes)
cgenes_df$up_in_b = plyr::revalue(as.character(cgenes_df$up_in_b), replace = c("TRUE" = "BHVE", "FALSE" = "CTRL"))
ggplot(cgenes_df, aes(x = mean_dif, y = sd_sum, color = mean_dif)) + geom_point() + geom_text_repel(data = cgenes_df[which(cgenes_df$mean_dif > 1),], aes(label = gene), color = "black") + scale_color_viridis_c() + xlab("Difference in BHVE Mean of Avg Expression - CTRL Mean of Avg Expression") + ylab("Sum of BHVE and CTRL Standard Deviation") + theme_bw()

clust15_sep_df = data.frame()
for (this_clust in 0:52) {
  print(this_clust)
  all_gene_df = aggregate(. ~ subsample, t_imp_mat2[which(bb$seuratclusters53 == this_clust),c(colnames(t_imp_mat2)[1:2000], "subsample")], mean)
  sep_cgenes = sapply(2:ncol(all_gene_df), function(x) all(outer(X = as.numeric(all_gene_df[1:19,x]), Y = as.numeric(all_gene_df[20:38,x]), FUN = '>')) | all(outer(X = as.numeric(all_gene_df[1:19,x]), Y = as.numeric(all_gene_df[20:38,x]), FUN = '<')) )
  up_cgenes = sapply(colnames(all_gene_df)[which(sep_cgenes) + 1], function(x) all_gene_df[1,x] > all_gene_df[20,x] )
  big_cgenes = sapply(colnames(all_gene_df)[which(sep_cgenes) + 1], function(x) abs(sum(all_gene_df[1:19,x]) - sum(all_gene_df[20:38,x])) )
  names(big_cgenes) = colnames(all_gene_df)[which(sep_cgenes) + 1]
  sd_cgenes = sapply(colnames(all_gene_df)[which(sep_cgenes) + 1], function(x) sum(sd(all_gene_df[1:19,x]), sd(all_gene_df[20:38,x])) )
  names(sd_cgenes) = colnames(all_gene_df)[which(sep_cgenes) + 1]
  if (length(which(sep_cgenes)) > 0) {
    this_df = data.frame(gene = names(big_cgenes), mean_dif = big_cgenes, sd_sum = sd_cgenes, up_in_b = up_cgenes, cluster = this_clust)
    this_df$up_in_b = plyr::revalue(as.character(this_df$up_in_b), replace = c("TRUE" = "BHVE", "FALSE" = "CTRL"))
    clust15_sep_df = rbind(clust15_sep_df, this_df)
  }
}

all_gene_df$cond = c(rep("BHVE", 19), rep("CTRL", 19))
ggplot(all_gene_df, aes(x = cond, y = npas3, fill = cond, color = cond)) + geom_boxplot(alpha = 0.6) + geom_point(position = position_jitter()) + scale_fill_manual(values = rev(viridis(2))) + scale_color_manual(values = rev(viridis(2)))

test_df = data.frame(satb1 = as.vector(bb@assays$RNA@data['satb1',]))
test_df[,c("seuratclusters15", "seuratclusters53", "subsample", "sample")] = bb@meta.data[,c("seuratclusters15", "seuratclusters53", "subsample", "sample")]
ggplot(test_df, aes(x = subsample, y = asic4, color = subsample)) + geom_boxplot() + geom_point(alpha = 0.01, position = position_jitterdodge()) + NoLegend() + ggtitle("Non-Imputed Data")

cor_df = data.frame(gene = colnames(imp_mat2), cor_counts = 0, cor_data = 0, row.names = colnames(imp_mat2))
for (gene in colnames(imp_mat2)) {
  print(gene)
  cor_df[gene, "cor_counts"] = cor(bb@assays$RNA@counts[gene,], imp_mat2[,gene])
  cor_df[gene, "cor_data"] = cor(bb@assays$RNA@data[gene,], imp_mat2[,gene])
}
cor_df[which(is.na(cor_df$cor_counts) | is.na(cor_df$cor_data)), c("cor_counts", "cor_data")] = 0
ggplot(cor_df, aes(cor_counts, cor_data, color = cor_counts)) + geom_point() + scale_color_gradientn(colors = plasma(100))

sum_df = data.frame(gene = colnames(imp_mat2), count_sum = 0, data_sum = 0, imp_sum = 0, row.names = colnames(imp_mat2))
sum_df$count_sum = rowSums(bb@assays$RNA@counts[colnames(imp_mat2),])
sum_df$data_sum = rowSums(bb@assays$RNA@data[colnames(imp_mat2),])
sum_df$imp_sum = colSums(imp_mat2[,colnames(imp_mat2)])
ggplot(sum_df, aes(count_sum, imp_sum, color = count_sum)) + geom_point() + scale_color_gradientn(colors = plasma(100)) + xlab("Sum Counts") + ylab("Sum Imputed")
ggplot(sum_df, aes(data_sum, imp_sum, color = data_sum)) + geom_point() + scale_color_gradientn(colors = plasma(100)) + xlab("Sum Data") + ylab("Sum Imputed")

# Plot the gene with the highest correlation with bower_activity_index
p_df = data.frame(value = imp_mat2["unc13b",], bai = bb$bower_activity_index)
ggplot(p_df, aes(x = bai, y = value)) + geom_point() + xlab("Bower Avitvity Index")

# Plot the gene with the highest correlation with bower_activity_index
p_df = data.frame(value1 = imp_mat2["unc13b",], value2 = imp_mat2["LOC101473063",], value3 = -imp_mat2["LOC101484787",], bai = bb$bower_activity_index)
p_df = p_df[order(p_df$bai, decreasing = F),]
# ggplot(p_df, aes(x = value1, y = value2, color = bai)) + geom_point() + xlab("Bower Avitvity Index") + scale_color_gradientn(colors = viridis(100))
plot_ly(p_df, x = ~value1, y = ~value2, z = ~value3, color = ~bai, size = I(10), colors = viridis(38))

imp_cor = sapply(rownames(imp_mat2), function(x) cor(imp_mat2[x,], bb$bower_activity_index))
imp_cor_df = data.frame(gene = rownames(imp_mat2), imp_cor = imp_cor, imp_cor_abs = abs(imp_cor))
imp_cor_df$sub_cor = 0
sub_meta = aggregate(bower_activity_index ~ subsample, bb@meta.data, mean)
for (i in 1:nrow(imp_mat2)) {
  gene = rownames(imp_mat2)[i]
  gdf = data.frame(value = imp_mat2[gene,], subsample = bb$subsample)
  gdf$subsample = factor(gdf$subsample, levels = sub_meta$subsample)
  gdf_agr = aggregate(value ~ subsample, gdf, mean, drop = F)
  gdf_agr$value[which(is.na(gdf_agr$value))] = 0
  gdf_agr$bai = sub_meta$bower_activity_index[match(gdf_agr$subsample, sub_meta$subsample)]
  imp_cor_df$sub_cor[i] = cor(gdf_agr$value, gdf_agr$bai)
}
imp_cor_df$sub_cor_abs = abs(imp_cor_df$sub_cor)


gdf = data.frame(value1 = imp_mat2["LOC101484787",], value2 = imp_mat2["LOC101468154",], subsample = bb$subsample)
gdf$subsample = factor(gdf$subsample, levels = sub_meta$subsample)
gdf_agr = aggregate(. ~ subsample, gdf, mean, drop = F)
gdf_agr$value1[which(is.na(gdf_agr$value1))] = 0
gdf_agr$value2[which(is.na(gdf_agr$value2))] = 0
gdf_agr$bai = sub_meta$bower_activity_index[match(gdf_agr$subsample, sub_meta$subsample)]
# ggplot(gdf_agr, aes(x = value, y = bai)) + geom_point()
ggplot(gdf_agr, aes(x = value1, y = value2, color = bai)) + geom_point() + scale_color_gradientn(colors = viridis(100))

#****************************************************************************
# Z-Test Random Genes =======================================================
#****************************************************************************
pcrc = read.csv("C:/Users/miles/Downloads/brain/data/markers/pcrc_FST20_30_LG11_evolution_genes_031821.csv")[,1]
gene_counts = data.frame(rowSums(bb@assays$RNA@counts))
gene_counts$gene = rownames(gene_counts)
gene_counts = gene_counts[order(gene_counts[,1]),]
pcrc_idx = which(gene_counts[,2] %in% pcrc)
ran_genes_1 = gene_counts[pcrc_idx + 1, 2]
ran_genes_2 = gene_counts[pcrc_idx + 2, 2]
ran_genes_3 = gene_counts[pcrc_idx + 3, 2]
ran_genes_4 = gene_counts[pcrc_idx + 4, 2]
ran_genes_5 = gene_counts[pcrc_idx + 5, 2]

real_res = markerExpPerCellPerClusterQuick(bb, pcrc)[[4]]
ran_1_df = markerExpPerCellPerClusterQuick(bb, ran_genes_1)[[4]]
ran_2_df = markerExpPerCellPerClusterQuick(bb, ran_genes_2)[[4]]
ran_3_df = markerExpPerCellPerClusterQuick(bb, ran_genes_3)[[4]]
ran_4_df = markerExpPerCellPerClusterQuick(bb, ran_genes_4)[[4]]
ran_5_df = markerExpPerCellPerClusterQuick(bb, ran_genes_5)[[4]]

big_df = rbind(real_res, ran_1_df, ran_2_df, ran_3_df, ran_4_df, ran_5_df)
big_df$real = c(rep("real", nrow(real_res)), rep(1, nrow(ran_1_df)), rep(2, nrow(ran_2_df)), rep(3, nrow(ran_3_df)), rep(4, nrow(ran_4_df)), rep(5, nrow(ran_5_df)))
ggplot(big_df, aes(x = real, y = avg_cluster_exp)) + geom_boxplot() + facet_wrap(~ cluster)

big_df2 = real_res[, c(1, 2)]
colnames(big_df2) = c("cluster", "real")
big_df2$ran1 = ran_1_df$avg_cluster_exp
big_df2$ran2 = ran_2_df$avg_cluster_exp
big_df2$ran3 = ran_3_df$avg_cluster_exp
big_df2$ran4 = ran_4_df$avg_cluster_exp
big_df2$ran5 = ran_5_df$avg_cluster_exp
write.csv(big_df2, "C:/Users/miles/Downloads/bb_pcrc_real_and_ran_53.csv")

cor(colSums(bb@assays$RNA@data[pcrc,]), colSums(bb@assays$RNA@data[ran_genes_1,]))
cor(colSums(bb@assays$RNA@data[pcrc,]), colSums(bb@assays$RNA@data[ran_genes_2,]))
cor(colSums(bb@assays$RNA@data[pcrc,]), colSums(bb@assays$RNA@data[ran_genes_3,]))
cor(colSums(bb@assays$RNA@data[pcrc,]), colSums(bb@assays$RNA@data[ran_genes_4,]))
cor(colSums(bb@assays$RNA@data[pcrc,]), colSums(bb@assays$RNA@data[ran_genes_5,]))

big_exp = GetAssayData(bb, assay = "RNA", slot="counts")
big_exp[which(big_exp > 1)] = 1
clusters = sort(unique(as.numeric(as.vector(Idents(bb)))))

ran_pools = list()
search_space = seq(-200, 200)
search_space = search_space[order(abs(search_space))][2:length(search_space)]
for (gene in pcrc) {
  gene_pcrc_idx = pcrc_idx[which(pcrc == gene)]
  ran_pools[[gene]] = c()
  search_space_i = 1
  while(length(ran_pools[[gene]]) < 100) {
    idx_to_try = gene_pcrc_idx + search_space[search_space_i]
    if (idx_to_try > 0 & idx_to_try <= nrow(bb)) {
      ran_pools[[gene]] = c(ran_pools[[gene]], gene_counts[idx_to_try, 2])
    }
    search_space_i = search_space_i + 1
  }
}

mySingleRun = function(markers) {
  exp = big_exp[markers,]
  
  per_cluster_df = data.frame()
  for (cluster in clusters) {
    cluster_cells <- WhichCells(bb, idents = cluster)
    avg_cluster_exp = colSums(exp[markers, cluster_cells])
    avg_cluster_exp = avg_cluster_exp/bb$nFeature_RNA[cluster_cells]
    mean_avg_cluster_exp = mean(avg_cluster_exp)
    per_cluster_df = rbind(per_cluster_df, t(c(cluster, mean_avg_cluster_exp)))
  }
  return(per_cluster_df[,2])
}

library("rhdf5")
# h5f = H5Fopen("~/scratch/brain/results/py_cor15/cluster15_1_cor_bhve.h5")
h5f = H5Fopen("~/scratch/brain/results/py_ns/real_no_bvc__cor_all.h5")
r_mat = h5f$name
h5closeAll()
rownames(r_mat) = colnames(r_mat) = rownames(bb)
gene_info = read.table("~/scratch/m_zebra_ref/gene_info.txt", sep="\t", header = T, stringsAsFactors = F)
cons_gene_info = gene_info[which(gene_info$human_mart == gene_info$human_pat),]
dup_gene_df = data.frame(table(cons_gene_info$human))
cons_gene_info_dup = cons_gene_info[which( cons_gene_info$human %in% dup_gene_df[which(dup_gene_df[,2] == 2),1] ),]
cons_gene_info_dup = cons_gene_info_dup[order(cons_gene_info_dup$human),]
dup_gene_pair = data.frame(gene1 = cons_gene_info_dup$mzebra[which(duplicated(cons_gene_info_dup$human))], gene2 = cons_gene_info_dup$mzebra[which(! duplicated(cons_gene_info_dup$human))], ens1 = cons_gene_info_dup$ens[which(duplicated(cons_gene_info_dup$human))], ens2 = cons_gene_info_dup$ens[which(! duplicated(cons_gene_info_dup$human))], human = cons_gene_info_dup$human[which(! duplicated(cons_gene_info_dup$human))])
dup_gene_pair$cor = diag(r_mat[dup_gene_pair$gene1, dup_gene_pair$gene2])
dup_gene_pair = dup_gene_pair[order(abs(dup_gene_pair$cor), decreasing = T),]
dup_gene_pair$n1 = gene_counts[match(dup_gene_pair$gene1, gene_counts[,2]),1]
dup_gene_pair$n2 = gene_counts[match(dup_gene_pair$gene2, gene_counts[,2]),1]

cons_gene_info = gene_info
cons_gene_info = cons_gene_info[which(!duplicated(cons_gene_info$mzebra)),]
dup_gene_df = data.frame(table(cons_gene_info$human_pat))
cons_gene_info_dup = cons_gene_info[which( cons_gene_info$human_pat %in% dup_gene_df[which(dup_gene_df[,2] == 2),1] ),]
cons_gene_info_dup = cons_gene_info_dup[order(cons_gene_info_dup$human_pat),]
dup_gene_pair2 = data.frame(gene1 = cons_gene_info_dup$mzebra[which(duplicated(cons_gene_info_dup$human_pat))], gene2 = cons_gene_info_dup$mzebra[which(! duplicated(cons_gene_info_dup$human_pat))], ens1 = cons_gene_info_dup$ens[which(duplicated(cons_gene_info_dup$human_pat))], ens2 = cons_gene_info_dup$ens[which(! duplicated(cons_gene_info_dup$human_pat))], human = cons_gene_info_dup$human_pat[which(! duplicated(cons_gene_info_dup$human_pat))])
dup_gene_pair2$cor = diag(r_mat[dup_gene_pair2$gene1, dup_gene_pair2$gene2])
dup_gene_pair2 = dup_gene_pair2[order(abs(dup_gene_pair2$cor), decreasing = T),]
dup_gene_pair2$n1 = gene_counts[match(dup_gene_pair2$gene1, gene_counts[,2]),1]
dup_gene_pair2$n2 = gene_counts[match(dup_gene_pair2$gene2, gene_counts[,2]),1]

goi = read.csv("~/scratch/brain/data/markers/goi_1plus_by_trial_id_and_cat_120121_hgnc.csv")
goi$X.1 = goi$X = NULL
goi_mat = r_mat[unique(goi$mzebra), unique(goi$mzebra)]
diag(goi_mat) = 0
global_hc = 0
callback = function(hc, mat){
  global_hc <<- hc
}
# pheatmap::pheatmap(goi_mat, clustering_callback = callback, cellwidth = 10, cellheight = 10, file = "~/scratch/brain/results/goi_heatmap_big.pdf")
pheatmap::pheatmap(goi_mat, clustering_callback = callback, file = "~/scratch/brain/results/goi_heatmap_small.pdf")

pop17 = c("th", "tac1", "LOC101464862", "LOC101480131", "npy", "igf2", "LOC101477361", "gad1", "adrb1", "gad2", "LOC101484392", "LOC101467991", "LOC101470177", "LOC101466282", "LOC101469465", "LOC101463985", "calcr")
pop17_mat = r_mat[unique(pop17), unique(pop17)]
diag(pop17_mat) = 0
pheatmap::pheatmap(pop17_mat, clustering_callback = callback, file = "~/scratch/brain/results/pop17_heatmap_small.pdf")

ran_df = read.csv("~/Downloads/pcrc_matched_exp_level_ran_lists_10k.csv")
ran_df$X = NULL
mat = as.matrix(bb@assays$RNA@counts)
mat[which(mat > 1)] = 1
pcrc = read.csv("~/research/brain/data/pcrc_FST20_30_LG11_evolution_genes_031821.csv")[,1]

library("parallel")
ran_score_sum = sapply(1:10000, function(x) sum(colSums(mat[as.character(ran_df[x,]),])) )
# ran_score_sum = unlist(mclapply(1:10000, function(x) sum(colSums(mat[as.character(ran_df[x,]),])), mc.cores = detectCores() ))
ran_counts_sum = sapply(1:10000, function(x) sum(colSums(bb@assays$RNA@counts[as.character(ran_df[x,]),])) )
ran_sum_df = data.frame(score = ran_score_sum, counts = ran_counts_sum, isReal = F)
real_sum_df = data.frame(score = sum(colSums(mat[pcrc,])), counts = sum(colSums(bb@assays$RNA@counts[as.character(pcrc),])), isReal = T)
all_sum_df = rbind(real_sum_df, ran_sum_df)
ggplot(all_sum_df, aes(x = isReal, y = score, fill = isReal, color = isReal)) + geom_boxplot(alpha = 0.2, outlier.shape = NA) + geom_point(position = position_jitter(), alpha = 0.1) + ggtitle("Matched Random Genes Have Significantly Higher Scores (p = 0.03)") + NoLegend()
ggplot(all_sum_df, aes(x = isReal, y = counts, fill = isReal, color = isReal)) + geom_boxplot(alpha = 0.2, outlier.shape = NA) + geom_point(position = position_jitter(), alpha = 0.1) + ggtitle("Matched Random Genes Have Similar Counts (By Creation)") + NoLegend()

full_ran_counts = as.vector(sapply(1:10000, function(x) colSums(bb@assays$RNA@counts[as.character(ran_df[x,]),]) ))
full_ran_score = as.vector(sapply(1:10000, function(x) colSums(mat[as.character(ran_df[x,]),]) ))
ran_score_mean = sapply(1:10000, function(x) mean(colSums(mat[as.character(ran_df[x,]),])) )
ran_counts_mean = sapply(1:10000, function(x) mean(colSums(bb@assays$RNA@counts[as.character(ran_df[x,]),])) )

full_ran_counts = readRDS("~/research/brain/data/full_ran_counts.RDS")
real_mean_df = data.frame(score = mean(colSums(mat[as.character(pcrc),])), counts = mean(colSums(bb@assays$RNA@counts[as.character(pcrc),])), isReal = T)
real_counts = colSums(bb@assays$RNA@counts[as.character(pcrc),])
all_counts_df = rbind(data.frame(counts = real_counts, isReal = T), data.frame(counts = full_ran_counts, isReal = F))
pdf("~/scratch/brain/results/real_ran_counts_density.pdf", width = 7, height = 6)
ggplot(all_counts_df, aes(counts, color = isReal, fill = isReal)) + geom_density()
dev.off()

rc = read.delim("rc_fst_closest_bedtools.bed", header = F)
rc = rc[which( abs(rc$V16) < 25000 ),]
rc$gene_name = colsplit(rc$V15, "; ", c('1', '2'))[, 1]
rc$gene_name = colsplit(rc$gene_name, "gene_id ", c('1', '2'))[, 2]
colnames(rc) = c("CHROM", "BIN_START", "BIN_END", "N_VAR", "WEIGHTED_FST", "MEAN_FST", "CHROM_1", "REF", "BIOTYPE", "GENE_START", "GENE_END", "GTF_1", "DIR", "GTF_2", "GENE_INFO", "DIST_TO_GENE", "GENE_NAME" )
rc$BIN_ID = paste0(rc$CHROM, "_", rc$BIN_START)

pc = read.delim("pc_fst_closest_bedtools.bed", header = F)
pc = pc[which( abs(pc$V16) < 25000 & pc$V1 == "NC_036790.1" & pc$V2 > 5950001 & pc$V3 < 25280000),]
pc$gene_name = colsplit(pc$V15, "; ", c('1', '2'))[, 1]
pc$gene_name = colsplit(pc$gene_name, "gene_id ", c('1', '2'))[, 2]
colnames(pc) = c("CHROM", "BIN_START", "BIN_END", "N_VAR", "WEIGHTED_FST", "MEAN_FST", "CHROM_1", "REF", "BIOTYPE", "GENE_START", "GENE_END", "GTF_1", "DIR", "GTF_2", "GENE_INFO", "DIST_TO_GENE", "GENE_NAME" )
pc$BIN_ID = paste0(pc$CHROM, "_", pc$BIN_START)

pcrc_merged = merge(pc, rc, by = "BIN_ID", suffixes = c("_PC", "_RC"))
pcrc_thresh = pcrc_merged[which(pcrc_merged$WEIGHTED_FST_PC >= 0.2 & pcrc_merged$WEIGHTED_FST_RC >= 0.3),]
pcrc_bin = sort(unique(pcrc_thresh$GENE_NAME))
write.csv(pcrc_bin, "pc_20_rc_30_10kb_bins_25kb_genes_on_lg_11_peak_by_bin.csv")
write.csv(pcrc_thresh, "pc_20_rc_30_10kb_bins_25kb_genes_on_lg_11_peak_by_bin_df.csv")

pcrc = sort(unique(rc$GENE_NAME[which(rc$GENE_NAME %in% pc$GENE_NAME)]))
write.csv(pcrc, "pc_20_rc_30_10kb_bins_25kb_genes_on_lg_11_peak.csv")

#*******************************************************************************
# Zack Models for IEG, Neurogen, and PCRC ======================================
#*******************************************************************************
ieg_dir = "C:/Users/miles/Downloads/ieg_real_vs_random_final/ieg_real_vs_random_final/"
ieg_dir15  = paste0(ieg_dir, "bb15/")
ieg_dir53  = paste0(ieg_dir, "bb53/")
ieg_dirgoi = paste0(ieg_dir, "by_goi/")

neurogen_dir = "C:/Users/miles/Downloads/neurogen_real_vs_random_final/"
neurogen_dir15  = paste0(neurogen_dir, "bb15/")
neurogen_dir53  = paste0(neurogen_dir, "bb53/")
neurogen_dirgoi = paste0(neurogen_dir, "by_goi/")

dirgoi = paste0(zdir, "by_goi/")
pcrc_dir = "C:/Users/miles/Downloads/pcrc_real_vs_random_final/"
pcrc_dir15  = paste0(dirgoi, "bb15/")
pcrc_dir53  = paste0(dirgoi, "bb53/")
pcrc_dirgoi = paste0(dirgoi, "by_goi/")

ieg_zdf15      = analyze15Dir(ieg_dir15)
neurogen_zdf15 = analyze15Dir(neurogen_dir15)
pcrc_zdf15     = analyze15Dir(pcrc_dir15)

ieg_zdf53      = analyze53Dir(ieg_dir53)
neurogen_zdf53 = analyze53Dir(neurogen_dir53)
pcrc_zdf53     = analyze53Dir(pcrc_dir53)

ieg_zdfgoi      = analyzegoiDir(ieg_dirgoi)
neurogen_zdfgoi = analyzegoiDir(neurogen_dirgoi)
pcrc_zdfgoi     = analyzegoiDir(pcrc_dirgoi)

ieg_zdf15$list      = ieg_zdf53$list      = ieg_zdfgoi$list      = "ieg"
neurogen_zdf15$list = neurogen_zdf53$list = neurogen_zdfgoi$list = "neurogen"
pcrc_zdf15$list     = pcrc_zdf53$list     = pcrc_zdfgoi$list     = "pcrc"

zdf15 = rbind(ieg_zdf15, neurogen_zdf15, pcrc_zdf15)
zdf15 = reshape2::melt(zdf15)
zdf15 = zdf15[which(zdf15$isReal == "real" & zdf15$variable %in% c("pct_sig", "neg_log_p_max", "neg_log_hmp")),]
zdf15$thresh = plyr::revalue(zdf15$variable, replace = c("pct_sig" = "100", "neg_log_p_max" = as.character(-log10(0.05)), "neg_log_hmp" = as.character(-log10(0.05))))
zdf15$thresh = as.numeric(as.vector(zdf15$thresh))
zdf15$list_cluster_exp_cond = paste(zdf15$list, zdf15$cluster, zdf15$exp_cond)
list_cluster_exp_cond_sig = sapply(unique(zdf15$list_cluster_exp_cond), function(x) { 
  this_df = zdf15[which(zdf15$list_cluster_exp_cond == x),]
  return(all(this_df$value >= this_df$thresh))
})
names(list_cluster_exp_cond_sig) = unique(zdf15$list_cluster_exp_cond)
zdf15$isSig = list_cluster_exp_cond_sig[match(zdf15$list_cluster_exp_cond, names(list_cluster_exp_cond_sig))]
my_pal = viridis(n = 3)
pdf("C:/Users/miles/Downloads/bb15_all_lists_15.pdf", width = 15, height = 7)
ggplot(zdf15, aes(x = cluster, y = value, shape = isSig, color = list, group = list)) + geom_point(aes(alpha = isSig, size = isSig)) + geom_line(alpha = 0.8) + scale_alpha_manual(values = c(0.6, 0.9)) + scale_size_manual(values = c(2, 2.3)) + scale_shape_manual(values = c(1, 19)) + scale_color_manual(values = my_pal) + xlab("") + facet_grid(variable ~ exp_cond, scales = "free") + geom_hline(aes(yintercept = thresh), linetype = "dashed", color = "gray40") + theme_bw()
dev.off()

zdf53 = rbind(ieg_zdf53, neurogen_zdf53, pcrc_zdf53)
zdf53 = reshape2::melt(zdf53)
zdf53 = zdf53[which(zdf53$isReal == "real" & zdf53$variable %in% c("pct_sig", "neg_log_p_max", "neg_log_hmp")),]
zdf53$thresh = plyr::revalue(zdf53$variable, replace = c("pct_sig" = "100", "neg_log_p_max" = as.character(-log10(0.05)), "neg_log_hmp" = as.character(-log10(0.05))))
zdf53$thresh = as.numeric(as.vector(zdf53$thresh))
zdf53$list_cluster_exp_cond = paste(zdf53$list, zdf53$cluster, zdf53$exp_cond)
list_cluster_exp_cond_sig = sapply(unique(zdf53$list_cluster_exp_cond), function(x) { 
  this_df = zdf53[which(zdf53$list_cluster_exp_cond == x),]
  return(all(this_df$value >= this_df$thresh))
})
names(list_cluster_exp_cond_sig) = unique(zdf53$list_cluster_exp_cond)
zdf53$isSig = list_cluster_exp_cond_sig[match(zdf53$list_cluster_exp_cond, names(list_cluster_exp_cond_sig))]
zdf53$group = paste(zdf53$list, zdf53$exp_cond)
my_pal = viridis(n = 3)
pdf("C:/Users/miles/Downloads/bb53_all_lists_bower_53.pdf", width = 15, height = 7)
ggplot(zdf53[which(zdf53$exp_cond == "bower"),], aes(x = cluster, y = value, shape = isSig, color = list, group = group)) + geom_point(aes(alpha = isSig, size = isSig)) + geom_line(alpha = 0.8) + scale_alpha_manual(values = c(0.6, 0.9)) + scale_size_manual(values = c(2, 2.3)) + scale_shape_manual(values = c(1, 19)) + scale_color_manual(values = my_pal) + xlab("") + facet_wrap(~ variable, scales = "free", ncol = 1) + geom_hline(aes(yintercept = thresh), linetype = "dashed", color = "gray40") + theme_bw()
dev.off()
pdf("C:/Users/miles/Downloads/bb53_all_lists_cond_53.pdf", width = 15, height = 7)
ggplot(zdf53[which(zdf53$exp_cond == "condition"),], aes(x = cluster, y = value, shape = isSig, color = list, group = group)) + geom_point(aes(alpha = isSig, size = isSig)) + geom_line(alpha = 0.8) + scale_alpha_manual(values = c(0.6, 0.9)) + scale_size_manual(values = c(2, 2.3)) + scale_shape_manual(values = c(1, 19)) + scale_color_manual(values = my_pal) + xlab("") + facet_wrap(~ variable, scales = "free", ncol = 1) + geom_hline(aes(yintercept = thresh), linetype = "dashed", color = "gray40") + theme_bw()
dev.off()
pdf("C:/Users/miles/Downloads/bb53_all_lists_bai_53.pdf", width = 15, height = 7)
ggplot(zdf53[which(zdf53$exp_cond == "bower activity index"),], aes(x = cluster, y = value, shape = isSig, color = list, group = group)) + geom_point(aes(alpha = isSig, size = isSig)) + geom_line(alpha = 0.8) + scale_alpha_manual(values = c(0.6, 0.9)) + scale_size_manual(values = c(2, 2.3)) + scale_shape_manual(values = c(1, 19)) + scale_color_manual(values = my_pal) + xlab("") + facet_wrap(~ variable, scales = "free", ncol = 1) + geom_hline(aes(yintercept = thresh), linetype = "dashed", color = "gray40") + theme_bw()
dev.off()
pdf("C:/Users/miles/Downloads/bb53_all_lists_gsi_53.pdf", width = 15, height = 7)
ggplot(zdf53[which(zdf53$exp_cond == "gsi"),], aes(x = cluster, y = value, shape = isSig, color = list, group = group)) + geom_point(aes(alpha = isSig, size = isSig)) + geom_line(alpha = 0.8) + scale_alpha_manual(values = c(0.6, 0.9)) + scale_size_manual(values = c(2, 2.3)) + scale_shape_manual(values = c(1, 19)) + scale_color_manual(values = my_pal) + xlab("") + facet_wrap(~ variable, scales = "free", ncol = 1) + geom_hline(aes(yintercept = thresh), linetype = "dashed", color = "gray40") + theme_bw()
dev.off()
pdf("C:/Users/miles/Downloads/bb53_all_lists_spawn_53.pdf", width = 15, height = 7)
ggplot(zdf53[which(zdf53$exp_cond == "spawn"),], aes(x = cluster, y = value, shape = isSig, color = list, group = group)) + geom_point(aes(alpha = isSig, size = isSig)) + geom_line(alpha = 0.8) + scale_alpha_manual(values = c(0.6, 0.9)) + scale_size_manual(values = c(2, 2.3)) + scale_shape_manual(values = c(1, 19)) + scale_color_manual(values = my_pal) + xlab("") + facet_wrap(~ variable, scales = "free", ncol = 1) + geom_hline(aes(yintercept = thresh), linetype = "dashed", color = "gray40") + theme_bw()
dev.off()

zdfgoi = rbind(ieg_zdfgoi, neurogen_zdfgoi, pcrc_zdfgoi)
zdfgoi = reshape2::melt(zdfgoi)
zdfgoi = zdfgoi[which(zdfgoi$isReal == "real" & zdfgoi$variable %in% c("pct_sig", "neg_log_p_max", "neg_log_hmp")),]
zdfgoi$thresh = plyr::revalue(zdfgoi$variable, replace = c("pct_sig" = "100", "neg_log_p_max" = as.character(-log10(0.05)), "neg_log_hmp" = as.character(-log10(0.05))))
zdfgoi$thresh = as.numeric(as.vector(zdfgoi$thresh))
zdfgoi$list_mzebra_exp_cond = paste(zdfgoi$list, zdfgoi$mzebra, zdfgoi$exp_cond)
list_mzebra_exp_cond_sig = sapply(unique(zdfgoi$list_mzebra_exp_cond), function(x) { 
  this_df = zdfgoi[which(zdfgoi$list_mzebra_exp_cond == x),]
  return(all(this_df$value >= this_df$thresh))
})
names(list_mzebra_exp_cond_sig) = unique(zdfgoi$list_mzebra_exp_cond)
zdfgoi$isSig = list_mzebra_exp_cond_sig[match(zdfgoi$list_mzebra_exp_cond, names(list_mzebra_exp_cond_sig))]
zdfgoi$group = paste(zdfgoi$list, zdfgoi$exp_cond)
my_pal = viridis(n = 3)
pdf("C:/Users/miles/Downloads/bbgoi_all_lists_bower_goi.pdf", width = 25, height = 7)
ggplot(zdfgoi[which(zdfgoi$exp_cond == "bower"),], aes(x = mzebra, y = value, shape = isSig, color = list, group = group)) + geom_point(aes(alpha = isSig, size = isSig)) + geom_line(alpha = 0.8) + scale_alpha_manual(values = c(0.6, 0.9)) + scale_size_manual(values = c(2, 2.3)) + scale_shape_manual(values = c(1, 19)) + scale_color_manual(values = my_pal) + xlab("") + facet_wrap(~ variable, scales = "free", ncol = 1) + geom_hline(aes(yintercept = thresh), linetype = "dashed", color = "gray40") + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()
pdf("C:/Users/miles/Downloads/bbgoi_all_lists_cond_goi.pdf", width = 25, height = 7)
ggplot(zdfgoi[which(zdfgoi$exp_cond == "condition"),], aes(x = mzebra, y = value, shape = isSig, color = list, group = group)) + geom_point(aes(alpha = isSig, size = isSig)) + geom_line(alpha = 0.8) + scale_alpha_manual(values = c(0.6, 0.9)) + scale_size_manual(values = c(2, 2.3)) + scale_shape_manual(values = c(1, 19)) + scale_color_manual(values = my_pal) + xlab("") + facet_wrap(~ variable, scales = "free", ncol = 1) + geom_hline(aes(yintercept = thresh), linetype = "dashed", color = "gray40") + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()
pdf("C:/Users/miles/Downloads/bbgoi_all_lists_bai_goi.pdf", width = 25, height = 7)
ggplot(zdfgoi[which(zdfgoi$exp_cond == "bower activity index"),], aes(x = mzebra, y = value, shape = isSig, color = list, group = group)) + geom_point(aes(alpha = isSig, size = isSig)) + geom_line(alpha = 0.8) + scale_alpha_manual(values = c(0.6, 0.9)) + scale_size_manual(values = c(2, 2.3)) + scale_shape_manual(values = c(1, 19)) + scale_color_manual(values = my_pal) + xlab("") + facet_wrap(~ variable, scales = "free", ncol = 1) + geom_hline(aes(yintercept = thresh), linetype = "dashed", color = "gray40") + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()
pdf("C:/Users/miles/Downloads/bbgoi_all_lists_gsi_goi.pdf", width = 25, height = 7)
ggplot(zdfgoi[which(zdfgoi$exp_cond == "gsi"),], aes(x = mzebra, y = value, shape = isSig, color = list, group = group)) + geom_point(aes(alpha = isSig, size = isSig)) + geom_line(alpha = 0.8) + scale_alpha_manual(values = c(0.6, 0.9)) + scale_size_manual(values = c(2, 2.3)) + scale_shape_manual(values = c(1, 19)) + scale_color_manual(values = my_pal) + xlab("") + facet_wrap(~ variable, scales = "free", ncol = 1) + geom_hline(aes(yintercept = thresh), linetype = "dashed", color = "gray40") + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()
pdf("C:/Users/miles/Downloads/bbgoi_all_lists_spawn_goi.pdf", width = 25, height = 7)
ggplot(zdfgoi[which(zdfgoi$exp_cond == "spawn"),], aes(x = mzebra, y = value, shape = isSig, color = list, group = group)) + geom_point(aes(alpha = isSig, size = isSig)) + geom_line(alpha = 0.8) + scale_alpha_manual(values = c(0.6, 0.9)) + scale_size_manual(values = c(2, 2.3)) + scale_shape_manual(values = c(1, 19)) + scale_color_manual(values = my_pal) + xlab("") + facet_wrap(~ variable, scales = "free", ncol = 1) + geom_hline(aes(yintercept = thresh), linetype = "dashed", color = "gray40") + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

analyze15Dir = function(dir15) {
  zdf15 = data.frame()
  for (f in list.files(dir15)) {
    zdf = read.csv(paste0(dir15, f))
    zdf_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_bbmm"))]
    zdf$num_sig = rowSums( zdf[, zdf_p_names] < 0.05 )
    zdf$pct_sig = (zdf$num_sig / length(zdf_p_names)) * 100
    zdf$p_max = sapply(1:nrow(zdf), function(x) max(zdf[x, zdf_p_names]) )
    zdf$isReal = "none"
    zdf$exp_cond = "none"
    
    # Real v Ran 1 v Ran 2
    if      (grepl("rand1", f)) { zdf$isReal = "rand1" } 
    else if (grepl("rand2", f)) { zdf$isReal = "rand2" } 
    else                        { zdf$isReal = "real"  }
    
    # Experimental Condition
    if      (grepl("bower", f)) { zdf$exp_cond = "bower" }
    else if (grepl("gsi", f))   { zdf$exp_cond = "gsi"   }
    else if (grepl("spawn", f)) { zdf$exp_cond = "spawn" }
    else                        { print("ERROR")         }
    
    if(zdf$exp_cond[1] == "bower") { zdf$hmp = zdf$hmp_bower_all }
    
    zdf15 = rbind(zdf15, zdf[,c('exp_cond', 'isReal', 'cluster', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
    
    # Split Bower into Cond and BAI
    if (zdf$exp_cond[1] == "bower") {
      zdf_cond_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_bbmm_cond"))]
      zdf_bai_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_bbmm_bower_activity"))]
      
      zdf_cond = zdf
      zdf_cond$exp_cond = "condition"
      zdf_bai = zdf
      zdf_bai$exp_cond = "bower activity index"
      
      zdf_cond$num_sig = rowSums( zdf[, zdf_cond_p_names] < 0.05 )
      zdf_cond$pct_sig = (zdf_cond$num_sig / length(zdf_cond_p_names)) * 100
      zdf_cond$p_max = sapply(1:nrow(zdf_cond), function(x) max(zdf_cond[x, zdf_cond_p_names]) )
      zdf_cond$hmp = zdf_cond$hmp_cond
      
      zdf_bai$num_sig = rowSums( zdf[, zdf_bai_p_names] < 0.05 )
      zdf_bai$pct_sig = (zdf_bai$num_sig / length(zdf_bai_p_names)) * 100
      zdf_bai$p_max = sapply(1:nrow(zdf_bai), function(x) max(zdf_bai[x, zdf_bai_p_names]) )
      zdf_bai$hmp = zdf_bai$hmp_cond
      
      zdf15 = rbind(zdf15, zdf_cond[,c('exp_cond', 'isReal', 'cluster', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
      zdf15 = rbind(zdf15, zdf_bai[ ,c('exp_cond', 'isReal', 'cluster', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
    }
  }
  
  zdf15$neg_log_p_max = -log10(zdf15$p_max)
  zdf15$neg_log_hmp = -log10(zdf15$hmp)
  zdf15$cluster = factor(zdf15$cluster, levels = unique(zdf15$cluster))
  zdf15$isReal = factor(zdf15$isReal, levels = c("real", "rand1", "rand2"))
  zdf15$group = paste(zdf15$isReal, zdf15$exp_cond)
  return(zdf15)
}

analyze53Dir = function(dir53) {
  
  zdf53 = data.frame()
  for (f in list.files(dir53)) {
    zdf = read.csv(paste0(dir53, f))
    zdf_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_bbmm"))]
    zdf$num_sig = rowSums( zdf[, zdf_p_names] < 0.05 )
    zdf$pct_sig = (zdf$num_sig / length(zdf_p_names)) * 100
    zdf$p_max = sapply(1:nrow(zdf), function(x) max(zdf[x, zdf_p_names]) )
    zdf$isReal = "none"
    zdf$exp_cond = "none"
    
    # Real v Ran 1 v Ran 2
    if      (grepl("rand1", f)) { zdf$isReal = "rand1" } 
    else if (grepl("rand2", f)) { zdf$isReal = "rand2" } 
    else                        { zdf$isReal = "real"  }
    
    # Experimental Condition
    if      (grepl("bower", f)) { zdf$exp_cond = "bower" }
    else if (grepl("gsi", f))   { zdf$exp_cond = "gsi"   }
    else if (grepl("spawn", f)) { zdf$exp_cond = "spawn" }
    else                        { print("ERROR")         }
    
    # zdf$hmp = zdf[, which(startsWith(colnames(zdf), "hmp")[1])]
    if(zdf$exp_cond[1] == "bower") { zdf$hmp = zdf$hmp_bower_all }
    
    zdf53 = rbind(zdf53, zdf[,c('exp_cond', 'isReal', 'cluster', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
    
    # Split Bower into Cond and BAI
    if (zdf$exp_cond[1] == "bower") {
      zdf_cond_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_bbmm_cond"))]
      zdf_bai_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_bbmm_bower_activity"))]
      
      zdf_cond = zdf
      zdf_cond$exp_cond = "condition"
      zdf_bai = zdf
      zdf_bai$exp_cond = "bower activity index"
      
      zdf_cond$num_sig = rowSums( zdf[, zdf_cond_p_names] < 0.05 )
      zdf_cond$pct_sig = (zdf_cond$num_sig / length(zdf_cond_p_names)) * 100
      zdf_cond$p_max = sapply(1:nrow(zdf_cond), function(x) max(zdf_cond[x, zdf_cond_p_names]) )
      zdf_cond$hmp = zdf_cond$hmp_cond
      
      zdf_bai$num_sig = rowSums( zdf[, zdf_bai_p_names] < 0.05 )
      zdf_bai$pct_sig = (zdf_bai$num_sig / length(zdf_bai_p_names)) * 100
      zdf_bai$p_max = sapply(1:nrow(zdf_bai), function(x) max(zdf_bai[x, zdf_bai_p_names]) )
      zdf_bai$hmp = zdf_bai$hmp_cond
      
      zdf53 = rbind(zdf53, zdf_cond[,c('exp_cond', 'isReal', 'cluster', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
      zdf53 = rbind(zdf53, zdf_bai[ ,c('exp_cond', 'isReal', 'cluster', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
    }
  }
  
  zdf53$neg_log_p_max = -log10(zdf53$p_max)
  zdf53$neg_log_hmp = -log10(zdf53$hmp)
  zdf53$cluster = factor(zdf53$cluster, levels = unique(zdf53$cluster))
  zdf53$isReal = factor(zdf53$isReal, levels = c("real", "rand1", "rand2"))
  zdf53$group = paste(zdf53$isReal, zdf53$exp_cond)
  return(zdf53)
}

analyzegoiDir = function(dirgoi) {
  zdfgoi = data.frame()
  for (f in list.files(dirgoi)) {
    zdf = read.csv(paste0(dirgoi, f))
    zdf_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_bbmm"))]
    zdf$num_sig = rowSums( zdf[, zdf_p_names] < 0.05 )
    zdf$pct_sig = (zdf$num_sig / length(zdf_p_names)) * 100
    zdf$p_max = sapply(1:nrow(zdf), function(x) max(zdf[x, zdf_p_names]) )
    zdf$isReal = "none"
    zdf$exp_cond = "none"
    zdf$human = zdf$goi
    
    # Real v Ran 1 v Ran 2
    if      (grepl("rand1", f)) { zdf$isReal = "rand1" }
    else if (grepl("rand2", f)) { zdf$isReal = "rand2" }
    else                        { zdf$isReal = "real"  }
    
    # Experimental Condition
    if      (grepl("bower", f)) { zdf$exp_cond = "bower" }
    else if (grepl("gsi", f))   { zdf$exp_cond = "gsi"   }
    else if (grepl("spawn", f)) { zdf$exp_cond = "spawn" }
    else                        { print("ERROR")         }
    
    # zdf$hmp = zdf[, which(startsWith(colnames(zdf), "hmp")[1])]
    if(zdf$exp_cond[1] == "bower") { zdf$hmp = zdf$hmp_bower_all }
    
    zdfgoi = rbind(zdfgoi, zdf[,c('exp_cond', 'isReal', 'mzebra', 'human', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
    
    # Split Bower into Cond and BAI
    if (zdf$exp_cond[1] == "bower") {
      zdf_cond_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_bbmm_cond"))]
      zdf_bai_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_bbmm_bower_activity"))]
      
      zdf_cond = zdf
      zdf_cond$exp_cond = "condition"
      zdf_bai = zdf
      zdf_bai$exp_cond = "bower activity index"
      
      zdf_cond$num_sig = rowSums( zdf[, zdf_cond_p_names] < 0.05 )
      zdf_cond$pct_sig = (zdf_cond$num_sig / length(zdf_cond_p_names)) * 100
      zdf_cond$p_max = sapply(1:nrow(zdf_cond), function(x) max(zdf_cond[x, zdf_cond_p_names]) )
      zdf_cond$hmp = zdf_cond$hmp_cond
      
      zdf_bai$num_sig = rowSums( zdf[, zdf_bai_p_names] < 0.05 )
      zdf_bai$pct_sig = (zdf_bai$num_sig / length(zdf_bai_p_names)) * 100
      zdf_bai$p_max = sapply(1:nrow(zdf_bai), function(x) max(zdf_bai[x, zdf_bai_p_names]) )
      zdf_bai$hmp = zdf_bai$hmp_cond
      
      zdfgoi = rbind(zdfgoi, zdf_cond[,c('exp_cond', 'isReal', 'mzebra', 'human', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
      zdfgoi = rbind(zdfgoi, zdf_bai[ ,c('exp_cond', 'isReal', 'mzebra', 'human', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
    }
  }
  
  zdfgoi$goi = zdfgoi$mzebra
  zdfgoi$neg_log_p_max = -log10(zdfgoi$p_max)
  zdfgoi$neg_log_hmp = -log10(zdfgoi$hmp)
  zdfgoi$goi = factor(zdfgoi$goi, levels = sort(unique(zdfgoi$goi)))
  zdfgoi$isReal = factor(zdfgoi$isReal, levels = c("real", "rand1", "rand2"))
  zdfgoi$group = paste(zdfgoi$isReal, zdfgoi$exp_cond)
  zdfgoi = zdfgoi[order(zdfgoi$isReal, decreasing = T),]
  return(zdfgoi)
}

# Help Zack
# zdir = "C:/Users/miles/Downloads/ieg_real_vs_random_final/ieg_real_vs_random_final/"
# zdir = "C:/Users/miles/Downloads/neurogen_real_vs_random_final/"
zdir = "C:/Users/miles/Downloads/pcrc_real_vs_random_final/"
dir15  = paste0(zdir, "bb15/")
dir53  = paste0(zdir, "bb53/")
dirgoi = paste0(zdir, "by_goi/")

zdf15 = data.frame()
for (f in list.files(dir15)) {
  zdf = read.csv(paste0(dir15, f))
  zdf_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_bbmm"))]
  zdf$num_sig = rowSums( zdf[, zdf_p_names] < 0.05 )
  zdf$pct_sig = (zdf$num_sig / length(zdf_p_names)) * 100
  zdf$p_max = sapply(1:nrow(zdf), function(x) max(zdf[x, zdf_p_names]) )
  zdf$isReal = "none"
  zdf$exp_cond = "none"
  
  # Real v Ran 1 v Ran 2
  if      (grepl("rand1", f)) { zdf$isReal = "rand1" } 
  else if (grepl("rand2", f)) { zdf$isReal = "rand2" } 
  else                        { zdf$isReal = "real"  }
  
  # Experimental Condition
  if      (grepl("bower", f)) { zdf$exp_cond = "bower" }
  else if (grepl("gsi", f))   { zdf$exp_cond = "gsi"   }
  else if (grepl("spawn", f)) { zdf$exp_cond = "spawn" }
  else                        { print("ERROR")         }
  
  # zdf$hmp = zdf[, which(startsWith(colnames(zdf), "hmp")[1])]
  if(zdf$exp_cond[1] == "bower") { zdf$hmp = zdf$hmp_bower_all }
  
  zdf15 = rbind(zdf15, zdf[,c('exp_cond', 'isReal', 'cluster', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
  
  # Split Bower into Cond and BAI
  if (zdf$exp_cond[1] == "bower") {
    zdf_cond_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_bbmm_cond"))]
    zdf_bai_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_bbmm_bower_activity"))]
    
    zdf_cond = zdf
    zdf_cond$exp_cond = "condition"
    zdf_bai = zdf
    zdf_bai$exp_cond = "bower activity index"
    
    zdf_cond$num_sig = rowSums( zdf[, zdf_cond_p_names] < 0.05 )
    zdf_cond$pct_sig = (zdf_cond$num_sig / length(zdf_cond_p_names)) * 100
    zdf_cond$p_max = sapply(1:nrow(zdf_cond), function(x) max(zdf_cond[x, zdf_cond_p_names]) )
    zdf_cond$hmp = zdf_cond$hmp_cond
    
    zdf_bai$num_sig = rowSums( zdf[, zdf_bai_p_names] < 0.05 )
    zdf_bai$pct_sig = (zdf_bai$num_sig / length(zdf_bai_p_names)) * 100
    zdf_bai$p_max = sapply(1:nrow(zdf_bai), function(x) max(zdf_bai[x, zdf_bai_p_names]) )
    zdf_bai$hmp = zdf_bai$hmp_cond
    
    zdf15 = rbind(zdf15, zdf_cond[,c('exp_cond', 'isReal', 'cluster', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
    zdf15 = rbind(zdf15, zdf_bai[ ,c('exp_cond', 'isReal', 'cluster', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
  }
}

zdf15$neg_log_p_max = -log10(zdf15$p_max)
zdf15$neg_log_hmp = -log10(zdf15$hmp)
# zdf15$realSig = zdf15$isReal == "real" & 
zdf15$cluster = factor(zdf15$cluster, levels = unique(zdf15$cluster))
zdf15$isReal = factor(zdf15$isReal, levels = c("real", "rand1", "rand2"))
zdf15$group = paste(zdf15$isReal, zdf15$exp_cond)
pdf("C:/Users/miles/Downloads/pcrc_bb15_p_max.pdf", width = 12, height = 10)
ggplot(zdf15, aes(x = cluster, y = neg_log_p_max, shape = exp_cond, color = isReal, group = group)) + geom_point(size = 2) + geom_line() + scale_color_manual(values = c("blue", "gray50", "gray70")) + xlab("") + ylab(expression(-Log["10"]*" Max P")) + facet_wrap(~ exp_cond, ncol = 2) + theme_bw()
dev.off()
pdf("C:/Users/miles/Downloads/pcrc_bb15_hmp.pdf", width = 12, height = 10)
ggplot(zdf15, aes(x = cluster, y = neg_log_hmp,   shape = exp_cond, color = isReal, group = group)) + geom_point(size = 2) + geom_line() + scale_color_manual(values = c("blue", "gray50", "gray70")) + xlab("") + ylab(expression(-Log["10"]*" HMP")) + facet_wrap(~ exp_cond, ncol = 2) + theme_bw()
dev.off()
pdf("C:/Users/miles/Downloads/pcrc_bb15_pct_sig.pdf", width = 12, height = 10)
ggplot(zdf15, aes(x = cluster, y = pct_sig,       shape = exp_cond, color = isReal, group = group)) + geom_point(size = 2) + geom_line() + scale_color_manual(values = c("blue", "gray50", "gray70")) + xlab("") + ylab("Percent of Significant Models") + facet_wrap(~ exp_cond, ncol = 2) + theme_bw() + scale_y_continuous(limits = c(0, 100))
dev.off()

zdf15_heat = data.frame()
for (ec in c("bower", "condition", "bower activity index", "gsi", "spawn")) {
  for (cluster in levels(zdf15$cluster)) {
    this_df = zdf15[which(zdf15$exp_cond == ec & zdf15$cluster == cluster),]
    
    zdf15_heat = rbind(zdf15_heat, data.frame(ec = ec, cluster = cluster, sig = this_df$hmp[which(this_df$isReal == "real")] < 0.05 & this_df$pct_sig[which(this_df$isReal == "real")] == 100 & this_df$pct_sig[which(this_df$isReal == "rand1")] < 100 & this_df$pct_sig[which(this_df$isReal == "rand2")] < 100 ))
  }
}
zdf15_heat$cluster = as.numeric(zdf15_heat$cluster)
zdf15_heat = acast(data = zdf15_heat, formula = ec ~ cluster) * 1
pheatmap::pheatmap(zdf15_heat, cluster_rows = F, cluster_cols = F, color = viridis(100), angle_col = 0, filename = "C:/Users/miles/Downloads/pcrc_bb15_binary.pdf", cellwidth = 25, cellheight = 25, legend = F)


zdf53 = data.frame()
for (f in list.files(dir53)) {
  zdf = read.csv(paste0(dir53, f))
  zdf_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_"))]
  zdf$num_sig = rowSums( zdf[, zdf_p_names] < 0.05 )
  zdf$pct_sig = (zdf$num_sig / length(zdf_p_names)) * 100
  zdf$p_max = sapply(1:nrow(zdf), function(x) max(zdf[x, zdf_p_names]) )
  zdf$isReal = "none"
  zdf$exp_cond = "none"
  
  # Real v Ran 1 v Ran 2
  if      (grepl("rand1", f)) { zdf$isReal = "rand1" } 
  else if (grepl("rand2", f)) { zdf$isReal = "rand2" } 
  else                        { zdf$isReal = "real"  }
  
  # Experimental Condition
  if      (grepl("bower", f)) { zdf$exp_cond = "bower" }
  else if (grepl("gsi", f))   { zdf$exp_cond = "gsi"   }
  else if (grepl("spawn", f)) { zdf$exp_cond = "spawn" }
  else                        { print("ERROR")         }
  
  # zdf$hmp = zdf[, which(startsWith(colnames(zdf), "hmp")[1])]
  if(zdf$exp_cond[1] == "bower") { zdf$hmp = zdf$hmp_bower_all }
  
  zdf53 = rbind(zdf53, zdf[,c('exp_cond', 'isReal', 'cluster', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
  
  # Split Bower into Cond and BAI
  if (zdf$exp_cond[1] == "bower") {
    zdf_cond_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_bbmm_cond"))]
    zdf_bai_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_bbmm_bower_activity"))]
    
    zdf_cond = zdf
    zdf_cond$exp_cond = "condition"
    zdf_bai = zdf
    zdf_bai$exp_cond = "bower activity index"
    
    zdf_cond$num_sig = rowSums( zdf[, zdf_cond_p_names] < 0.05 )
    zdf_cond$pct_sig = (zdf_cond$num_sig / length(zdf_cond_p_names)) * 100
    zdf_cond$p_max = sapply(1:nrow(zdf_cond), function(x) max(zdf_cond[x, zdf_cond_p_names]) )
    zdf_cond$hmp = zdf_cond$hmp_cond
    
    zdf_bai$num_sig = rowSums( zdf[, zdf_bai_p_names] < 0.05 )
    zdf_bai$pct_sig = (zdf_bai$num_sig / length(zdf_bai_p_names)) * 100
    zdf_bai$p_max = sapply(1:nrow(zdf_bai), function(x) max(zdf_bai[x, zdf_bai_p_names]) )
    zdf_bai$hmp = zdf_bai$hmp_cond
    
    zdf53 = rbind(zdf53, zdf_cond[,c('exp_cond', 'isReal', 'cluster', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
    zdf53 = rbind(zdf53, zdf_bai[ ,c('exp_cond', 'isReal', 'cluster', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
  }
}

zdf53$neg_log_p_max = -log10(zdf53$p_max)
zdf53$neg_log_hmp = -log10(zdf53$hmp)
zdf53$cluster = factor(zdf53$cluster, levels = unique(zdf53$cluster))
zdf53$isReal = factor(zdf53$isReal, levels = c("real", "rand1", "rand2"))
zdf53$group = paste(zdf53$isReal, zdf53$exp_cond)
pdf("C:/Users/miles/Downloads/pcrc_bb53_p_max.pdf", width = 20, height = 10)
ggplot(zdf53, aes(x = cluster, y = neg_log_p_max, shape = exp_cond, color = isReal, group = group)) + geom_point(size = 2) + geom_line() + scale_color_manual(values = c("blue", "gray50", "gray70")) + xlab("") + ylab(expression(-Log["10"]*" Max P")) + facet_wrap(~ exp_cond, ncol = 2) + theme_bw()
dev.off()
pdf("C:/Users/miles/Downloads/pcrc_bb53_hmp.pdf", width = 20, height = 10)
ggplot(zdf53, aes(x = cluster, y = neg_log_hmp,   shape = exp_cond, color = isReal, group = group)) + geom_point(size = 2) + geom_line() + scale_color_manual(values = c("blue", "gray50", "gray70")) + xlab("") + ylab(expression(-Log["10"]*" HMP")) + facet_wrap(~ exp_cond, ncol = 2) + theme_bw()
dev.off()
pdf("C:/Users/miles/Downloads/pcrc_bb53_pct_sig.pdf", width = 20, height = 10)
ggplot(zdf53, aes(x = cluster, y = pct_sig,       shape = exp_cond, color = isReal, group = group)) + geom_point(size = 2) + geom_line() + scale_color_manual(values = c("blue", "gray50", "gray70")) + xlab("") + ylab("Percent of Significant Models") + facet_wrap(~ exp_cond, ncol = 2) + theme_bw() + scale_y_continuous(limits = c(0, 100))
dev.off()

zdf53_heat = data.frame()
for (ec in c("bower", "condition", "bower activity index", "gsi", "spawn")) {
  for (cluster in levels(zdf53$cluster)) {
    this_df = zdf53[which(zdf53$exp_cond == ec & zdf53$cluster == cluster),]
    
    zdf53_heat = rbind(zdf53_heat, data.frame(ec = ec, cluster = cluster, sig = this_df$hmp[which(this_df$isReal == "real")] < 0.05 & this_df$pct_sig[which(this_df$isReal == "real")] == 100 & this_df$pct_sig[which(this_df$isReal == "rand1")] < 100 & this_df$pct_sig[which(this_df$isReal == "rand2")] < 100 ))
  }
}
zdf53_heat$cluster = as.numeric(zdf53_heat$cluster)
zdf53_heat = acast(data = zdf53_heat, formula = ec ~ cluster) * 1
pheatmap::pheatmap(zdf53_heat, cluster_rows = F, cluster_cols = F, color = viridis(100), angle_col = 0, filename = "C:/Users/miles/Downloads/pcrc_bb53_binary.pdf", cellwidth = 25, cellheight = 25, legend = F)


zdfgoi = data.frame()
for (f in list.files(dirgoi)) {
  zdf = read.csv(paste0(dirgoi, f))
  zdf_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_"))]
  zdf$num_sig = rowSums( zdf[, zdf_p_names] < 0.05 )
  zdf$pct_sig = (zdf$num_sig / length(zdf_p_names)) * 100
  zdf$p_max = sapply(1:nrow(zdf), function(x) max(zdf[x, zdf_p_names]) )
  zdf$isReal = "none"
  zdf$exp_cond = "none"
  zdf$human = zdf$goi

  # Real v Ran 1 v Ran 2
  if      (grepl("rand1", f)) { zdf$isReal = "rand1" }
  else if (grepl("rand2", f)) { zdf$isReal = "rand2" }
  else                        { zdf$isReal = "real"  }

  # Experimental Condition
  if      (grepl("bower", f)) { zdf$exp_cond = "bower" }
  else if (grepl("gsi", f))   { zdf$exp_cond = "gsi"   }
  else if (grepl("spawn", f)) { zdf$exp_cond = "spawn" }
  else                        { print("ERROR")         }

  # zdf$hmp = zdf[, which(startsWith(colnames(zdf), "hmp")[1])]
  if(zdf$exp_cond[1] == "bower") { zdf$hmp = zdf$hmp_bower_all }

  zdfgoi = rbind(zdfgoi, zdf[,c('exp_cond', 'isReal', 'mzebra', 'human', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
  
  # Split Bower into Cond and BAI
  if (zdf$exp_cond[1] == "bower") {
    zdf_cond_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_bbmm_cond"))]
    zdf_bai_p_names = colnames(zdf)[which(startsWith(colnames(zdf), "p_bbmm_bower_activity"))]
    
    zdf_cond = zdf
    zdf_cond$exp_cond = "condition"
    zdf_bai = zdf
    zdf_bai$exp_cond = "bower activity index"
    
    zdf_cond$num_sig = rowSums( zdf[, zdf_cond_p_names] < 0.05 )
    zdf_cond$pct_sig = (zdf_cond$num_sig / length(zdf_cond_p_names)) * 100
    zdf_cond$p_max = sapply(1:nrow(zdf_cond), function(x) max(zdf_cond[x, zdf_cond_p_names]) )
    zdf_cond$hmp = zdf_cond$hmp_cond
    
    zdf_bai$num_sig = rowSums( zdf[, zdf_bai_p_names] < 0.05 )
    zdf_bai$pct_sig = (zdf_bai$num_sig / length(zdf_bai_p_names)) * 100
    zdf_bai$p_max = sapply(1:nrow(zdf_bai), function(x) max(zdf_bai[x, zdf_bai_p_names]) )
    zdf_bai$hmp = zdf_bai$hmp_cond
    
    zdfgoi = rbind(zdfgoi, zdf_cond[,c('exp_cond', 'isReal', 'mzebra', 'human', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
    zdfgoi = rbind(zdfgoi, zdf_bai[ ,c('exp_cond', 'isReal', 'mzebra', 'human', 'num_sig', 'pct_sig', 'p_max', 'hmp')])
  }
}

zdfgoi$goi = zdfgoi$mzebra
zdfgoi$neg_log_p_max = -log10(zdfgoi$p_max)
zdfgoi$neg_log_hmp = -log10(zdfgoi$hmp)
zdfgoi$goi = factor(zdfgoi$goi, levels = sort(unique(zdfgoi$goi)))
zdfgoi$isReal = factor(zdfgoi$isReal, levels = c("real", "rand1", "rand2"))
zdfgoi$group = paste(zdfgoi$isReal, zdfgoi$exp_cond)
zdfgoi = zdfgoi[order(zdfgoi$isReal, decreasing = T),]
pdf("C:/Users/miles/Downloads/pcrc_bbgoi_p_max.pdf", width = 25, height = 10)
ggplot(zdfgoi, aes(x = goi, y = neg_log_p_max, shape = exp_cond, color = isReal, group = group)) + geom_point(size = 2) + geom_line() + scale_color_manual(values = c("blue", "gray50", "gray70")) + xlab("") + ylab(expression(-Log["10"]*" Max P")) + facet_wrap(~ exp_cond, ncol = 1) + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()
pdf("C:/Users/miles/Downloads/pcrc_bbgoi_hmp.pdf", width = 25, height = 10)
ggplot(zdfgoi, aes(x = goi, y = neg_log_hmp,   shape = exp_cond, color = isReal, group = group)) + geom_point(size = 2) + geom_line() + scale_color_manual(values = c("blue", "gray50", "gray70")) + xlab("") + ylab(expression(-Log["10"]*" HMP")) + facet_wrap(~ exp_cond, ncol = 1) + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()
pdf("C:/Users/miles/Downloads/pcrc_bbgoi_pct_sig.pdf", width = 25, height = 10)
ggplot(zdfgoi, aes(x = goi, y = pct_sig,       shape = exp_cond, color = isReal, group = group)) + geom_point(size = 2) + geom_line() + scale_color_manual(values = c("blue", "gray50", "gray70")) + xlab("") + ylab("Percent of Significant Models") + facet_wrap(~ exp_cond, ncol = 1) + theme_bw() + scale_y_continuous(limits = c(0, 100)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

zdfgoi_heat = data.frame()
incomplete_goi = c()
for (ec in c("bower", "condition", "bower activity index", "gsi", "spawn")) {
  for (goi in levels(zdfgoi$goi)) {
    this_df = zdfgoi[which(zdfgoi$exp_cond == ec & zdfgoi$goi == goi),]
    if (nrow(this_df) == 3)
      zdfgoi_heat = rbind(zdfgoi_heat, data.frame(ec = ec, goi = goi, sig = this_df$hmp[which(this_df$isReal == "real")] < 0.05 & this_df$pct_sig[which(this_df$isReal == "real")] == 100 & this_df$pct_sig[which(this_df$isReal == "rand1")] < 100 & this_df$pct_sig[which(this_df$isReal == "rand2")] < 100 ))
    else
      incomplete_goi = c(incomplete_goi, goi)
  }
}
zdfgoi_heat = acast(data = zdfgoi_heat, formula = ec ~ goi) * 1
zdfgoi_heat = zdfgoi_heat[, which(colSums(zdfgoi_heat) > 0)]
colnames(zdfgoi_heat) = paste0(colnames(zdfgoi_heat), " (", gene_info$human[match(colnames(zdfgoi_heat), gene_info$mzebra)], ")")
pheatmap::pheatmap(zdfgoi_heat, cluster_rows = F, cluster_cols = F, color = viridis(100), angle_col = 45, filename = "C:/Users/miles/Downloads/pcrc_bbgoi_binary.pdf", cellwidth = 25, cellheight = 25, legend = F)


zdf15_heat_ieg = read.csv("C:/Users/miles/Downloads/ieg_sum15.csv")
zdf53_heat_ieg = read.csv("C:/Users/miles/Downloads/ieg_sum53.csv")
zdfgoi_heat_ieg = read.csv("C:/Users/miles/Downloads/ieg_sumgoi.csv")
zdf15_heat_ieg$X = paste("IEG: ", zdf15_heat_ieg$X)
zdf53_heat_ieg$X = paste("IEG: ", zdf53_heat_ieg$X)
zdfgoi_heat_ieg$X = paste("IEG: ", zdfgoi_heat_ieg$X)

zdf15_heat_neurogen = read.csv("C:/Users/miles/Downloads/neurogen_sum15.csv")
zdf53_heat_neurogen = read.csv("C:/Users/miles/Downloads/neurogen_sum53.csv")
zdfgoi_heat_neurogen = read.csv("C:/Users/miles/Downloads/neurogen_sumgoi.csv")
zdf15_heat_neurogen$X = paste("Neurogenesis: ", zdf15_heat_neurogen$X)
zdf53_heat_neurogen$X = paste("Neurogenesis: ", zdf53_heat_neurogen$X)
zdfgoi_heat_neurogen$X = paste("Neurogenesis: ", zdfgoi_heat_neurogen$X)

zdf15_heat_pcrc = read.csv("C:/Users/miles/Downloads/pcrc_sum15.csv")
zdf53_heat_pcrc = read.csv("C:/Users/miles/Downloads/pcrc_sum53.csv")
zdfgoi_heat_pcrc = read.csv("C:/Users/miles/Downloads/pcrc_sumgoi.csv")
zdf15_heat_pcrc$X = paste("FST: ", zdf15_heat_pcrc$X)
zdf53_heat_pcrc$X = paste("FST: ", zdf53_heat_pcrc$X)
zdfgoi_heat_pcrc$X = paste("FST: ", zdfgoi_heat_pcrc$X)

sum15 = rbind(zdf15_heat_ieg, zdf15_heat_neurogen, zdf15_heat_pcrc)
rownames(sum15) = sum15$X
sum15$X = NULL
colnames(sum15) = substr(colnames(sum15), 2, 100)
pheatmap::pheatmap(sum15, cluster_rows = F, cluster_cols = F, color = viridis(100), angle_col = 0, filename = "C:/Users/miles/Downloads/sum_bb15_binary.pdf", cellwidth = 25, cellheight = 25, legend = F)

sum53 = rbind(zdf53_heat_ieg, zdf53_heat_neurogen, zdf53_heat_pcrc)
rownames(sum53) = sum53$X
sum53$X = NULL
colnames(sum53) = substr(colnames(sum53), 2, 100)
pheatmap::pheatmap(sum53, cluster_rows = F, cluster_cols = F, color = viridis(100), angle_col = 0, filename = "C:/Users/miles/Downloads/sum_bb53_binary.pdf", cellwidth = 25, cellheight = 25, legend = F)

sumgoi = rbind(zdfgoi_heat_ieg, zdfgoi_heat_neurogen, zdfgoi_heat_pcrc)
rownames(sumgoi) = sumgoi$X
sumgoi$X = NULL
colnames(sumgoi) = paste0(colnames(sumgoi), " (", gene_info$human[match(colnames(sumgoi), gene_info$mzebra)], ")")
# colnames(sumgoi) = str_replace(colnames(sumgoi), "\\..", " (")
# colnames(sumgoi) = str_replace(colnames(sumgoi), "\\.", ")")
sumgoi = sumgoi[, which(colSums(sumgoi) > 0)]
# colnames(sumgoi)[ncol(sumgoi)] = "nkx2.1 (NKX2.1)"
pheatmap::pheatmap(sumgoi, cluster_rows = F, cluster_cols = T, color = viridis(100), angle_col = 45, filename = "C:/Users/miles/Downloads/sum_bbgoi_binary.pdf", cellwidth = 25, cellheight = 25, legend = F)


test = read.csv("C:/Users/miles/Downloads/doi_10.5061_dryad.8t8s248__v1/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv")
ran_lists = read.csv("C:/Users/miles/Downloads/pcrc_matched_exp_level_ran_lists_10k.csv")

bb$rand1 = colSums(mat[as.character(ran_lists[2,]),])
df = data.frame(value = c(bb$pcrc, bb$rand1), real = c(rep(T, ncol(bb)), rep(F, ncol(bb))))
# ggplot(df, aes(x = value, color = real, fill = real)) + geom_density(alpha = 0.2)
ggplot(df, aes(x = value, color = real, fill = real)) + geom_histogram(alpha = 0.2, position = 'identity', binwidth = 1) + ggtitle("Score")

df = data.frame(value = c(colSums(bb@assays$RNA@counts[pcrc,]), colSums(bb@assays$RNA@counts[as.character(ran_lists[2,]),])), real = c(rep(T, ncol(bb)), rep(F, ncol(bb))))
ggplot(df, aes(x = value, color = real, fill = real)) + geom_histogram(alpha = 0.2, position = 'identity', binwidth = 1) + ggtitle("Counts")

library(ggExtra)
i = 2
df = data.frame(Counts = c(colSums(bb@assays$RNA@counts[pcrc,]), colSums(bb@assays$RNA@counts[as.character(ran_lists[i,]),])),
                Score  = c(bb$pcrc, colSums(mat[as.character(ran_lists[i,]),])),
                real = c(rep(T, ncol(bb)), rep(F, ncol(bb))))
ggplot(df, aes(x = Counts, y = Score, color = real, fill = real)) + geom_jitter(alpha = 0.5, width = 0.2, height = 0.2) + ggtitle(paste0("Score vs Counts (i = ", i, ")"))
# p = ggplot(df, aes(x = Counts, y = Score, color = real, fill = real)) + geom_point() + ggtitle(paste0("Score vs Counts (i = ", i, ")"))
# ggExtra::ggMarginal(p, type = "histogram")

fst = read.csv("C:/Users/miles/Downloads/pcrc_fst_file_for_george.csv")
pcrc = read.csv("C:/Users/miles/Downloads/pc_20_rc_30_10kb_bins_25kb_genes_on_lg_11_peak.csv")[,2]
pcrc_bin = read.csv("C:/Users/miles/Downloads/pc_20_rc_30_10kb_bins_25kb_genes_on_lg_11_peak_by_bin.csv")[,2]
pcrc_bin_bin = read.csv("C:/Users/miles/Downloads/pc_20_rc_30_10kb_bins_25kb_genes_on_lg_11_peak_by_bin_df.csv")
pcrc_bin_bin$lg = lgConverter(pcrc_bin_bin$CHROM_PC)
pcrc_bin_bin$bin = paste0(pcrc_bin_bin$lg, "_", pcrc_bin_bin$BIN_START_PC)
fst$hasGene = F
fst$hasPCRC = F
fst$hasPCRC_BinGene = F
fst$hasPCRC_BinOf_BinGene = fst$bin %in% pcrc_bin_bin$bin
fst$pcrcGene = ""
fst$pcrcBinGene = ""
bin_size = 10000
gtf$lg = lgConverter(gtf$V1)
gtf$gene_start_round = (floor(gtf$V4  / bin_size) * bin_size) + 1
gtf$gene_end_round   = ceiling(gtf$V5 / bin_size) * bin_size
for (i in 1:nrow(gtf)) {
  if (i %% 1000 == 0)
    cat(paste0(i, "."))
  gene_start_round = gtf$gene_start_round[i]
  gene_end_round = gtf$gene_end_round[i]
  gene_bins = seq(gene_start_round, gene_end_round, by = bin_size)
  gene_bins = paste0(gtf$lg[i], "_", gene_bins)
  fst$hasGene[match(gene_bins, fst$bin)] = T
  if(gtf$gene_name[i] %in% pcrc) {
    fst$hasPCRC[match(gene_bins, fst$bin)] = T
    fst$pcrcGene[match(gene_bins, fst$bin)] = gtf$gene_name[i]
  }
  if(gtf$gene_name[i] %in% pcrc_bin) {
    fst$hasPCRC_BinGene[match(gene_bins, fst$bin)] = T
    fst$pcrcBinGene[match(gene_bins, fst$bin)] = gtf$gene_name[i]
  }
}
write.csv(fst, "C:/Users/miles/Downloads/pcrc_fst_file_for_zack_120221.csv")

# org_list = human_ids$NCBI.gene..formerly.Entrezgene..ID[match(org_list, human_ids$Gene.name)]
# Urn = All Genes In Category. White Ball = Term Genes. Black Ball = Non Term Genes
# k = unname(unlist(df_merge[1, "GenesInQuery1"])) # Number of balls drawn from urn
# m = unname(unlist(df_merge[1, "GenesInTerm1"])) # Number of White Balls
# n = unname(unlist(df_merge[1, "TotalGenes1"])) - unname(unlist(df_merge[1, "GenesInTerm1"])) # Number of black balls
# x = unname(unlist(df_merge[1, "GenesInTermInQuery1"])) # Number of white balls drawn from urn
# phyper(q = x, m = m, n = n, k = k)
# signif(phyper(x, m, n, k) - cumsum(dhyper(x, m, n, k)), digits = 8)

#************************************************************************
# Zack Models for Neurogen - 100 Perms ==================================
#************************************************************************
library(lazyeval)
gene_list_name = "neurogen"
cluster_level = 15
zmodel = lapply(1:7, function(x) { 
  this_zm = read.csv(paste0(rna_path, "bbmm", cluster_level, "_", gene_list_name, "_model", x, ".csv"))
  colnames(this_zm) = paste0(colnames(this_zm), "_", x)
  return(this_zm)})

big_zm = do.call(cbind, zmodel)
big_zm[, which( startsWith(colnames(big_zm), "X") | startsWith(colnames(big_zm), "perm") | startsWith(colnames(big_zm), "cluster") )] = NULL
big_zm$perm = zmodel[[1]]$perm
big_zm$cluster = zmodel[[1]]$cluster
big_zm = big_zm %>% relocate(perm:cluster, .before = 1)

# Find Maximum P Value per Variable
big_zm = big_zm %>% mutate(cond_max_p  = pmax(!!!rlang::syms(colnames(big_zm)[which(startsWith(colnames(big_zm), "cond_p"))])))
big_zm = big_zm %>% mutate(bai_max_p   = pmax(!!!rlang::syms(colnames(big_zm)[which(startsWith(colnames(big_zm), "bower_activity_index_p"))])))
big_zm = big_zm %>% mutate(gsi_max_p   = pmax(!!!rlang::syms(colnames(big_zm)[which(startsWith(colnames(big_zm), "gsi_p"))])))
big_zm = big_zm %>% mutate(spawn_max_p = pmax(!!!rlang::syms(colnames(big_zm)[which(startsWith(colnames(big_zm), "log_spawn_events_p"))])))
big_zm = big_zm %>% mutate(bower_max_p = pmax(!!!rlang::syms(colnames(big_zm)[which(startsWith(colnames(big_zm), "cond_p") | startsWith(colnames(big_zm), "bower_activity_index_p"))])))

# Find Number of Significant Models per Variable
big_zm$cond_num_sig  = big_zm %>% select(!!!rlang::syms(colnames(big_zm)[which(startsWith(colnames(big_zm), "cond_p"))])) %>% mutate(cond_num_sig  = rowSums(. < 0.05)) %>% pull('cond_num_sig')
big_zm$bai_num_sig   = big_zm %>% select(!!!rlang::syms(colnames(big_zm)[which(startsWith(colnames(big_zm), "bower_activity_index_p"))])) %>% mutate(bai_num_sig  = rowSums(. < 0.05)) %>% pull('bai_num_sig')
big_zm$gsi_num_sig   = big_zm %>% select(!!!rlang::syms(colnames(big_zm)[which(startsWith(colnames(big_zm), "gsi_p"))])) %>% mutate(gsi_num_sig  = rowSums(. < 0.05)) %>% pull('gsi_num_sig')
big_zm$spawn_num_sig = big_zm %>% select(!!!rlang::syms(colnames(big_zm)[which(startsWith(colnames(big_zm), "log_spawn_events_p"))])) %>% mutate(spawn_num_sig  = rowSums(. < 0.05)) %>% pull('spawn_num_sig')
big_zm$bower_num_sig = big_zm %>% select(!!!rlang::syms(colnames(big_zm)[which(startsWith(colnames(big_zm), "cond_p") | startsWith(colnames(big_zm), "bower_activity_index_p"))])) %>% mutate(bower_num_sig  = rowSums(. < 0.05)) %>% pull('bower_num_sig')

# Plot Num Sig
zm_p = rbindlist(lapply(0:11, function(x)  {this_zm = big_zm[which(big_zm$cluster == x),]; data.frame(cond = length(which(this_zm$cond_num_sig == 3)), bai = length(which(this_zm$bai_num_sig == 3)), bower = length(which(this_zm$bower_num_sig == 6)), gsi = length(which(this_zm$gsi_num_sig == 5)), spawn = length(which(this_zm$spawn_num_sig == 5)) ) }))
zm_p$cluster = 0:11
zm_p = melt(zm_p, id.var = c('cluster'))
zm_p$cluster = factor(zm_p$cluster, levels = 0:11)
ggplot(zm_p, aes(x = cluster, y = value, fill = variable)) + geom_bar(stat = 'identity', position = position_dodge()) + scale_y_continuous(expand = c(0, 0)) + theme_classic() + ylab("# of Perms Meeting # of Sig Model Threshold") + xlab("Cluster") + ggtitle("Neurogenesis 15 Cluster Level")

# Plot P Max
zm_p = melt(big_zm[, c('cluster', 'cond_max_p', 'bai_max_p', 'bower_max_p', 'gsi_max_p', 'spawn_max_p')], id.var = c('cluster'))
zm_p$cluster = factor(zm_p$cluster, levels = 0:11)
zm_p$value = -log10(zm_p$value)
ggplot(zm_p, aes(x = cluster, y = value, fill = variable)) + geom_boxplot(position = position_dodge()) + scale_y_continuous(expand = c(0, 0)) + theme_classic() + ylab("-Log10(Max P)") + xlab("Cluster") + ggtitle("Neurogenesis 15 Cluster Level") + geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'gray40')

tj.deg.sig.pos$neg_log_bon = -log10(tj.deg.sig.pos$p_val_adj)
paste0("BON v 9: ", cor(tj.deg.sig.pos$p_val_adj[which(is.finite(tj.deg.sig.pos$all_metric9))], tj.deg.sig.pos$all_metric9[which(is.finite(tj.deg.sig.pos$all_metric9))]))
paste0("BON v 8: ", cor(tj.deg.sig.pos$p_val_adj[which(is.finite(tj.deg.sig.pos$all_metric8))], tj.deg.sig.pos$all_metric8[which(is.finite(tj.deg.sig.pos$all_metric8))]))
paste0("BON v 7: ", cor(tj.deg.sig.pos$p_val_adj[which(is.finite(tj.deg.sig.pos$all_metric7))], tj.deg.sig.pos$all_metric7[which(is.finite(tj.deg.sig.pos$all_metric7))]))
paste0("BON v 6: ", cor(tj.deg.sig.pos$p_val_adj[which(is.finite(tj.deg.sig.pos$all_metric6))], tj.deg.sig.pos$all_metric6[which(is.finite(tj.deg.sig.pos$all_metric6))]))
paste0("BON v 5: ", cor(tj.deg.sig.pos$p_val_adj[which(is.finite(tj.deg.sig.pos$all_metric5))], tj.deg.sig.pos$all_metric5[which(is.finite(tj.deg.sig.pos$all_metric5))]))
paste0("BON v 4: ", cor(tj.deg.sig.pos$p_val_adj[which(is.finite(tj.deg.sig.pos$all_metric4))], tj.deg.sig.pos$all_metric4[which(is.finite(tj.deg.sig.pos$all_metric4))]))
paste0("BON v 3: ", cor(tj.deg.sig.pos$p_val_adj[which(is.finite(tj.deg.sig.pos$all_metric3))], tj.deg.sig.pos$all_metric3[which(is.finite(tj.deg.sig.pos$all_metric3))]))
paste0("BON v 2: ", cor(tj.deg.sig.pos$p_val_adj[which(is.finite(tj.deg.sig.pos$all_metric2))], tj.deg.sig.pos$all_metric2[which(is.finite(tj.deg.sig.pos$all_metric2))]))
paste0("BON v 1: ", cor(tj.deg.sig.pos$p_val_adj[which(is.finite(tj.deg.sig.pos$all_metric))], tj.deg.sig.pos$all_metric[which(is.finite(tj.deg.sig.pos$all_metric))]))

tj.deg.sig.pos$num_1 = tj.deg.sig.pos$num.cells * tj.deg.sig.pos$pct.1
tj.deg.sig.pos$num_2 = (ncol(tj) - tj.deg.sig.pos$num.cells) * tj.deg.sig.pos$pct.2

tj.deg.sig.pos$all_metric_tmp = tj.deg.sig.pos$pct_dif * tj.deg.sig.pos$avg_log2FC * log(tj.deg.sig.pos$num.cells)
cor(tj.deg.sig.pos$neg_log_bon[which(is.finite(tj.deg.sig.pos$all_metric_tmp))], tj.deg.sig.pos$all_metric_tmp[which(is.finite(tj.deg.sig.pos$all_metric_tmp))])

tj.deg$num.cells = tj.num.cells$Freq[match(tj.deg$cluster, tj.num.cells$cluster)]
tj.deg$pct_dif = tj.deg$pct.1 - tj.deg$pct.2
tj.deg$abs_pct_dif = abs(tj.deg$pct_dif)
tj.deg$mym = tj.deg$pct_dif * tj.deg$avg_log2FC * log(tj.deg$num.cells)
tj.deg$neg_log_bon = -log10(tj.deg$p_val_adj)
tj.deg$isSig = tj.deg$p_val_adj < 0.05
cor(tj.deg$neg_log_bon, tj.deg$mym)
ggplot(tj.deg, aes(x = neg_log_bon, y = mym, color = as.numeric(isSig))) + geom_point() + geom_smooth(method='lm', formula= y~x) 
ggplot(tj.deg, aes(x = isSig, y = mym, color = as.numeric(isSig))) + geom_boxplot() + geom_point(alpha = 0.1, position = position_jitterdodge())
ggplot(tj.deg, aes(x = isSig, y = neg_log_bon, color = as.numeric(isSig))) + geom_boxplot() + geom_point(alpha = 0.1, position = position_jitterdodge())

bb15$all_metric_tmp = bb15$pct_dif * bb15$avg_logFC * log(bb15$num.cells)
cor(bb15$neg_log_bon[which(is.finite(bb15$all_metric_tmp))], bb15$all_metric_tmp[which(is.finite(bb15$all_metric_tmp))])

testdf = rbind(data.frame(value = p_df15[,172+1], level = "15"), data.frame(value = p_df53[,213+1], level = "53"))
ggplot(testdf, aes(x = value, color = level, fill = level)) + geom_histogram(alpha = 0.25) + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) + theme_bw()
10000 - length(which(testdf$value[which(testdf$level == "15")] > 0.73955740316691))
10000 - length(which(testdf$value[which(testdf$level == "53")] > 0.754682436447482))

# bvc = read.delim("C:/Users/miles/Downloads/brain/results/bb/sim_bhve_v_ctrl.tsv")
bb15 = read.csv("C:/Users/miles/Downloads/brain/results/bb/bb_all_cluster_15_degs.csv")
# bb15 = read.csv("C:/Users/miles/Downloads/brain/results/bb/bb_all_cluster_53_degs.csv")
bb.num.cells = as.data.frame(table(bb$seuratclusters15))
bb15$num.cells = bb.num.cells$Freq[match(bb15$cluster, bb.num.cells$Var1)]
bb15$pct_dif = bb15$pct.1 - bb15$pct.2
bb15$abs_pct_dif = abs(bb15$pct_dif)
bb15$mym = bb15$pct_dif * bb15$avg_logFC * log(bb15$num.cells)
bb15$neg_log_bon = -log10(bb15$p_val_adj)
bb15$isSig = bb15$p_val_adj < 0.05
bb15$neg_log_bon[which(! is.finite(bb15$neg_log_bon) )] = 500
cor(bb15$neg_log_bon[which(bb15$neg_log_bon != 500)], bb15$mym[which(bb15$neg_log_bon != 500)])
ggplot(bb15, aes(x = neg_log_bon, y = mym, color = as.numeric(isSig))) + geom_point() + geom_smooth(method='lm', formula= y~x) 
ggplot(bb15, aes(x = isSig, y = mym, color = as.numeric(isSig))) + geom_boxplot() + geom_point(alpha = 0.1, position = position_jitterdodge())
ggplot(bb15, aes(x = isSig, y = neg_log_bon, color = as.numeric(isSig))) + geom_boxplot() + geom_point(alpha = 0.1, position = position_jitterdodge())

# Best so far: pct_dif * avg_log2FC * log(num.cells)

tj.deg$all_metric9 = tj.deg$pct_dif * tj.deg$avg_log2FC * log(tj.deg$num.cells)
tj.deg$all_metric8 = tj.deg$pct_dif * tj.deg$avg_log2FC
paste0("LOG BON v 9: ", cor(tj.deg$neg_log_bon[which(is.finite(tj.deg$all_metric9))], tj.deg$all_metric9[which(is.finite(tj.deg$all_metric9))]))
paste0("LOG BON v 8: ", cor(tj.deg$neg_log_bon[which(is.finite(tj.deg$all_metric8))], tj.deg$all_metric8[which(is.finite(tj.deg$all_metric8))]))

paste0("LOG BON v 9: ", cor(tj.deg.sig.pos$neg_log_bon[which(is.finite(tj.deg.sig.pos$all_metric9))], tj.deg.sig.pos$all_metric9[which(is.finite(tj.deg.sig.pos$all_metric9))]))
paste0("LOG BON v 8: ", cor(tj.deg.sig.pos$neg_log_bon[which(is.finite(tj.deg.sig.pos$all_metric8))], tj.deg.sig.pos$all_metric8[which(is.finite(tj.deg.sig.pos$all_metric8))]))
paste0("LOG BON v 7: ", cor(tj.deg.sig.pos$neg_log_bon[which(is.finite(tj.deg.sig.pos$all_metric7))], tj.deg.sig.pos$all_metric7[which(is.finite(tj.deg.sig.pos$all_metric7))]))
paste0("LOG BON v 6: ", cor(tj.deg.sig.pos$neg_log_bon[which(is.finite(tj.deg.sig.pos$all_metric6))], tj.deg.sig.pos$all_metric6[which(is.finite(tj.deg.sig.pos$all_metric6))]))
paste0("LOG BON v 5: ", cor(tj.deg.sig.pos$neg_log_bon[which(is.finite(tj.deg.sig.pos$all_metric5))], tj.deg.sig.pos$all_metric5[which(is.finite(tj.deg.sig.pos$all_metric5))]))
paste0("LOG BON v 4: ", cor(tj.deg.sig.pos$neg_log_bon[which(is.finite(tj.deg.sig.pos$all_metric4))], tj.deg.sig.pos$all_metric4[which(is.finite(tj.deg.sig.pos$all_metric4))]))
paste0("LOG BON v 3: ", cor(tj.deg.sig.pos$neg_log_bon[which(is.finite(tj.deg.sig.pos$all_metric3))], tj.deg.sig.pos$all_metric3[which(is.finite(tj.deg.sig.pos$all_metric3))]))
paste0("LOG BON v 2: ", cor(tj.deg.sig.pos$neg_log_bon[which(is.finite(tj.deg.sig.pos$all_metric2))], tj.deg.sig.pos$all_metric2[which(is.finite(tj.deg.sig.pos$all_metric2))]))
# paste0("LOG BON v 1: ", cor(tj.deg.sig.pos$neg_log_bon[which(is.finite(tj.deg.sig.pos$all_metric))], tj.deg.sig.pos$all_metric[which(is.finite(tj.deg.sig.pos$all_metric))]))

tj.deg.sig.pos$cluster_gene = paste(tj.deg.sig.pos$cluster, tj.deg.sig.pos$gene)
my_f = mean(tj.deg.sig.pos$neg_log_bon) / mean(tj.deg.sig.pos$all_metric8)
tj.deg.sig.pos$tmp = abs(tj.deg.sig.pos$all_metric8 * my_f - tj.deg.sig.pos$neg_log_bon)
ggplot(tj.deg.sig.pos, aes(x = p_val_adj, y = all_metric3)) + geom_point()
ggplot(tj.deg.sig.pos, aes(x = neg_log_bon, y = all_metric3, color = avg_log2FC)) + geom_point() + geom_smooth(method='lm', formula= y~x) + geom_text_repel(data = tj.deg.sig.pos[which(tj.deg.sig.pos$tmp > 25),], aes(label = cluster_gene))

glmmseq_counts = readRDS("~/scratch/brain/results/deg_lmer_demux/53_clusters/glmmseq_counts.rds")
test_names = str_replace(rownames(glmmseq_counts), pattern = "\\.", "-")
rownames(glmmseq_counts) = test_names
colnames(glmmseq_counts) = str_replace(colnames(glmmseq_counts), pattern = "\\.", "-")
a = as.matrix(glmmseq_counts)
b = as.matrix(bb@assays$RNA@counts[test_names, which(bb$seuratclusters53 == 0)])
identical(a, b)

convert53 = read.csv("/path/to/convert53.csv")
clust53_new_col_list2 = convert53$col
names(clust53_new_col_list2) = colsplit(convert53$new, "_", c("num", "ex"))[,1]

all_pairs = 1:length(unique(bb$pair))
names(all_pairs) = unique(bb$pair)
bb$pair2 = plyr::revalue(bb$pair, replace = all_pairs)
p_list = list()
p_list[["empty"]] = ggplot() + geom_point() + theme_void()
for (pair in unique(bb$pair2)) {
  print(pair)
  pair_cells = colnames(bb)[which(bb$pair2 == pair)]
  b_cells = colnames(bb)[which(bb$pair2 == pair & bb$cond == "BHVE")]
  c_cells = colnames(bb)[which(bb$pair2 == pair & bb$cond == "CTRL")]
  # p_list[[pair]] = DimPlot(bb, cells = pair_cells)
  p1 = DimPlot(bb, cells = b_cells) + theme_void() + scale_color_manual(values = clust53_new_col_list2) + NoLegend()
  p2 = DimPlot(bb, cells = c_cells) + theme_void() + scale_color_manual(values = clust53_new_col_list2) + NoLegend()
  p3 = plot_grid(plotlist=list(p1, p2), ncol = 2)
  title <- ggdraw() + draw_label(pair)
  p4 = plot_grid(title, p3, ncol=1, rel_heights=c(0.1, 1))
  p_list[[pair]] = p4
}
pdf("C:/Users/miles/Downloads/test.pdf", width = 6, height = 15)
p = plot_grid(plotlist=p_list, ncol = 2)
print(p)
dev.off()

gcm = glmmseq_counts
res = results@stats
i = which(rownames(res) == "wdr17")

my_int = res[i, "(Intercept)"]
my_bower = res[i, "bower_activity_index"]
my_gsi   = res[i, "gsi"]
my_spawn = res[i, "log_spawn_events"]

bai_cell = bb$bower_activity_index[which(bb$seuratclusters53 == 0)]
gsi_cell = bb$gsi[which(bb$seuratclusters53 == 0)]
spawn_cell = bb$log_spawn_events[which(bb$seuratclusters53 == 0)]

gcm[i, ] = glmmseq_counts[i, ] - my_gsi * gsi_cell - my_spawn * spawn_cell
gcm_df = data.frame(adj_counts = t(gcm[i,]), sample = bb$sample[which(bb$seuratclusters53 == 0)], cond = bb$cond[which(bb$seuratclusters53 == 0)], subsample = bb$subsample[which(bb$seuratclusters53 == 0)], pair = bb$pair[which(bb$seuratclusters53 == 0)], bai = bai_cell, gsi = gsi_cell, spawn = spawn_cell)
colnames(gcm_df)[1] = "adj_counts"
gcm_df_agr = aggregate(adj_counts ~ subsample + cond + pair + bai, gcm_df, mean)

pdf("~/scratch/brain/results/bb53_0_wdr17_glmmseq.pdf", width = 4, height = 4)
print(ggplot(gcm_df_agr, aes(x = cond, y = adj_counts, color = cond)) + geom_boxplot(outlier.shape = NA) + geom_line(linetype = "dashed", color = "gray40", alpha = 0.4, aes(group=pair)) + geom_point(alpha = 0.75, position = position_jitterdodge()))
dev.off()

pdf("~/scratch/brain/results/bb53_0_wdr17_glmmseq_w_bai.pdf", width = 5, height = 5)
print(ggplot(gcm_df_agr, aes(x = adj_counts, y = bai, color = cond, label = subsample)) + geom_point(alpha = 0.75))
dev.off()

cz = read.csv("C:/Users/miles/Downloads/zdeg15_summary.csv")
adj_sub_mean = read.csv("C:/Users/miles/Downloads/zdeg15_adj_subsample_means.csv")
cz$cluster_gene = paste0(cz$cluster, "_", cz$mzebra)
rownames(cz) = cz$cluster_gene
rownames(adj_sub_mean) = cz$cluster_gene
cz$adj_sign_pair_abs = 9.5 + abs(9.5 - cz$adj_sign_pair)
cz$data_sign_pair_abs = 9.5 + abs(9.5 - cz$data_sign_pair)
cz$counts_sign_pair_abs = 9.5 + abs(9.5 - cz$counts_sign_pair)
cz$adj_mean_of_mean_dif = cz$adj_mean_of_mean_b - cz$adj_mean_of_mean_c
cz$data_mean_of_mean_dif = cz$data_mean_of_mean_b - cz$data_mean_of_mean_c
cz$counts_mean_of_mean_dif = cz$counts_mean_of_mean_b - cz$counts_mean_of_mean_c
cz$adj_mean_of_mean_dif_abs = abs(cz$adj_mean_of_mean_dif)
cz$data_mean_of_mean_dif_abs = abs(cz$data_mean_of_mean_dif)
cz$counts_mean_of_mean_dif_abs = abs(cz$counts_mean_of_mean_dif)
cz$adj_mean_dif = cz$adj_mean_b - cz$adj_mean_c
cz$data_mean_dif = cz$data_mean_b - cz$data_mean_c
cz$counts_mean_dif = cz$counts_mean_b - cz$counts_mean_c
cz$adj_mean_dif_abs = abs(cz$adj_mean_dif)
cz$data_mean_dif_abs = abs(cz$data_mean_dif)
cz$counts_mean_dif_abs = abs(cz$counts_mean_dif)
cz$neg_log_p_bower = -log10(cz$P_bower_activity_index)
cz$neg_log_p_gsi = -log10(cz$P_gsi)
cz$neg_log_p_spawn = -log10(cz$P_log_spawn_events)
cz$bower_p_sig = cz$P_bower_activity_index < 0.05

cz_sig = cz[which(cz$sig_bower_behavior == 1),]
adj_sub_mean_sig = adj_sub_mean[which(cz$sig_bower_behavior == 1),]
hist(cz_sig$adj_sign_pair_abs, breaks = 9, main = "BVC DEGs on 15 level - Adjusted", xlab = "Number of Pairs w/ Consistent Direction")
hist(cz_sig$data_sign_pair_abs, breaks = 9, main = "BVC DEGs on 15 level - Data", xlab = "Number of Pairs w/ Consistent Direction")
hist(cz_sig$counts_sign_pair_abs, breaks = 9, main = "BVC DEGs on 15 level - Counts", xlab = "Number of Pairs w/ Consistent Direction")

p_df = as.data.frame(table(paste0(cz$bower_p_sig, "_", cz$adj_sign_pair_abs)))
p_df[, c("bower_p_sig", "adj_sign_pair_abs")] = reshape2::colsplit(p_df$Var1, "_", c("1", "2"))
p_df$pct = ( p_df$Freq / length(which(cz$bower_p_sig)) ) * 100
p_df$pct[which(!p_df$bower_p_sig)] = ( p_df$Freq[which(!p_df$bower_p_sig)] / length(which(! cz$bower_p_sig)) ) * 100
# ggplot(p_df, aes(x = adj_sign_pair_abs, y = Freq, fill = bower_p_sig)) + geom_bar(stat = 'identity', position = position_dodge2()) + ggtitle("Adjusted")
ggplot(p_df, aes(x = adj_sign_pair_abs, y = pct, fill = bower_p_sig)) + geom_bar(stat = 'identity', position = position_dodge2()) + xlab("Number of Pairs w/ Consistent Direction") + ylab("Percent of Rows") + ggtitle("adj")

ggplot(cz_sig, aes(x = adj_mean_of_mean_b, y = adj_mean_of_mean_c)) + geom_point()
ggplot(cz_sig, aes(x = adj_mean_of_mean_dif, y = data_mean_of_mean_dif)) + geom_point()
ggplot(cz_sig, aes(x = adj_mean_of_mean_dif, y = counts_mean_of_mean_dif)) + geom_point()
ggplot(cz_sig, aes(x = adj_mean_dif, y = counts_mean_dif)) + geom_point()

ggplot(cz_sig, aes(x = adj_mean_of_mean_dif, y = bower_activity_index)) + geom_point()
ggplot(cz_sig, aes(x = adj_mean_dif, y = bower_activity_index)) + geom_point()
ggplot(cz_sig, aes(x = adj_mean_of_mean_dif_abs, y = neg_log_p_bower)) + geom_point()
ggplot(cz_sig, aes(x = adj_mean_dif_abs, y = neg_log_p_bower)) + geom_point()

ggplot(cz_sig, aes(x = adj_mean_dif_abs, y = neg_log_p_bower, color = adj_sign_pair_abs)) + geom_point() + scale_color_gradientn(colors = viridis(100))
ggplot(cz_sig, aes(x = neg_log_p_gsi, y = neg_log_p_bower)) + geom_point()
ggplot(cz_sig, aes(x = neg_log_p_spawn, y = neg_log_p_bower)) + geom_point()
ggplot(cz_sig, aes(x = adj_mean_of_mean_dif, y = bower_activity_index, color = neg_log_p_bower)) + geom_point() + scale_color_gradientn(colors = viridis(100))

# ClownFish ERE
cdf = read.delim("~/scratch/brain/results/clown_ere_closest.bed", header = F)
cdf$dist = cdf[,ncol(cdf)]
cdf$loc = reshape2::colsplit(cdf[,12], ";", c("1", "2"))[,1]
cdf$loc = reshape2::colsplit(cdf$loc, " ", c("1", "2"))[,2]
cdf$gene_name = reshape2::colsplit(cdf[,12], "; gene_source", c("1", "2"))[,1]
cdf$gene_name = reshape2::colsplit(cdf$gene_name, "gene_name ", c("1", "2"))[,2]
cdf$gene = cdf$gene_name
cdf$gene[which(cdf$gene == "")] = cdf$loc[which(cdf$gene == "")]
# cdf = cdf[which(cdf[,12] != "."),]
# cdf$gene = cdf$gene_name = cdf$loc
cdf2 = cdf[which( abs(cdf$dist) < 25000 ),]
cdf2$class = "distal"
cdf2$class[which(cdf2$dist <= 5000 & cdf2$dist > 0)] = "promoter"
cdf2$class[which(cdf2$dist == 0)] = "intragenic"
cdf3 = cdf2[, c("gene", "class")]
write.csv(cdf3, "~/scratch/brain/results/clown_ere_class.csv")

cz_sig = z53[which(z53$sig_bower_behavior == 1),]
cz_sig = cz_sig[order(cz_sig$P_bower_activity_index, decreasing = F),]
cz_sig$zgenes = str_replace(cz_sig$mzebra, "\\.", "-")
for (i in 1:5) {
  print(pairedBoxViolin(bb, gcm53, cz_sig$zgenes[i], 53, cz_sig$cluster[i]) + ggtitle(i))
}

#************************************************************************
# IEG Summary Figure ====================================================
#************************************************************************
# ieg_sum = read.csv("C:/Users/miles/Downloads/summary_sig_and_trend_IEG_hits_all_analyses_010321.csv")
# ieg_sum = read.csv("C:/Users/miles/Downloads/ieg_summary_data_by_cat_cluster_goi_010321.csv")
ieg_sum = read.csv("C:/Users/miles/Downloads/ieg_bower_quiver_gsi_cluster_goi_beta_effects_for_dotplot.csv")
ieg_sum = ieg_sum[which(ieg_sum$plot),]
ieg_sum$hgnc = ieg_sum$gene
ieg_sum$gene_pop = ieg_sum$goi
ieg_sum$cluster[which(ieg_sum$cluster == FALSE)] = "All"
ieg_sum$gene_pop[which(ieg_sum$gene_pop == FALSE)] = "All"
ieg_sum$level_old = paste0(ieg_sum$level, "_", ieg_sum$cluster)
ieg_sum$level_old_gp = paste0(ieg_sum$level_old, "_", ieg_sum$gene_pop)
ieg_sum$cat_level_old_gp = paste0(ieg_sum$cat, "_", ieg_sum$level_old_gp)

convert_all = data.frame(cluster = c(convert15$new.full, convert53$new), color = c(convert15$col, convert53$col), old = c(convert15$old, convert53$old))
convert_all = convert_all[which(! duplicated(convert_all$cluster) ),]
convert_all[, c("new.id", "new.gaba")]    = colsplit(convert_all$cluster, "_", names = c("new.id", "new.gaba"))
convert_all[, c("new.parent", "new.sub")] = colsplit(convert_all$new.id, "\\.", names = c("new.id", "new.gaba"))
convert_all$new.sub[which(is.na(convert_all$new.sub))] = 0
convert_all[which(convert_all$cluster == "8-9_Glut"), c("new.parent", "new.sub")] = c(8, 12) # Assign a parent and subcluster to the 8-9_Glut special case
convert_all = convert_all[order(as.numeric(convert_all$new.parent), as.numeric(convert_all$new.sub), decreasing = F),]
convert_all$level = plyr::revalue(as.character(convert_all$new.sub == 0), replace = c("TRUE" = "primary", "FALSE" = "secondary"))
convert_all = rbind(data.frame(cluster = "All", color = viridis(1), new.id = 0, new.gaba = "all", new.parent = 0, new.sub = 0, level = "all", old = "All"), convert_all)
convert_all$level_old = paste0(convert_all$level, "_", convert_all$old)

all_combos = expand.grid(convert_all$level_old, unique(ieg_sum$gene_pop[which(ieg_sum$gene_pop != FALSE)]) )
colnames(all_combos) = c("level_old", "gene_pop")
all_combos$level_old_gp = paste0(all_combos$level_old, "_", all_combos$gene_pop)
all_combos[,colnames(convert_all)] = convert_all[match(all_combos$level_old, convert_all$level_old),]
all_combos$hgnc = gene_info$human[match(all_combos$gene_pop, gene_info$mzebra)]
all_combos$hgnc[which(is.na(all_combos$hgnc))] = "All"
all_combos$bsig = NA
all_combos$gsig = NA
all_combos$qsig = NA

for (i in 1:nrow(ieg_sum)) {
  my.sig = TRUE
  my.cat = ieg_sum$cat[i]
  my.x = ieg_sum$level_old_gp[i]
  my.idx = which(all_combos$level_old_gp == my.x)
  my.cat.col = plyr::revalue(my.cat, replace = c("bower" = "bsig", "gsi" = "gsig", "quiver" = "qsig"))
  all_combos[my.idx, my.cat.col] = my.sig
}

all_combos$pcol = "white"
non_na_rows = which( ! is.na(all_combos$bsig) | ! is.na(all_combos$gsig) | ! is.na(all_combos$qsig) )
all_combos$pcol[non_na_rows] = unlist(lapply(non_na_rows, function(x) {
  this.col = "error"
  if ( ! is.na(all_combos[x,"bsig"]) ) { this.col = "#FDE725" }
  if ( ! is.na(all_combos[x,"gsig"]) ) { this.col = "#2AB07F" }
  if ( ! is.na(all_combos[x,"qsig"]) ) { this.col = "#433E85" }
  if ( ! is.na(all_combos[x,"bsig"]) & ! is.na(all_combos[x,"gsig"]) ) { this.col = "#94CC52" }
  if ( ! is.na(all_combos[x,"bsig"]) & ! is.na(all_combos[x,"qsig"]) ) { this.col = "#A09355" }
  if ( ! is.na(all_combos[x,"gsig"]) & ! is.na(all_combos[x,"qsig"]) ) { this.col = "#377782" }
  if ( ! is.na(all_combos[x,"bsig"]) & ! is.na(all_combos[x,"gsig"]) & ! is.na(all_combos[x,"qsig"]) ) { this.col = "#799c63" }
  return(this.col)
}))

all_combos$tran = "FALSE"
non_na_rows = which( ! is.na(all_combos$bsig) | ! is.na(all_combos$gsig) | ! is.na(all_combos$qsig) )
all_combos$tran[non_na_rows] = unlist(lapply(non_na_rows, function(x) {
  this.tran = "FF"
  if ( ! is.na(all_combos[x,"bsig"]) &  ! all_combos[x,"bsig"] ) { this.tran = "85" }
  if ( ! is.na(all_combos[x,"gsig"]) &  ! all_combos[x,"gsig"] ) { this.tran = "85" }
  if ( ! is.na(all_combos[x,"qsig"]) &  ! all_combos[x,"qsig"] ) { this.tran = "85" }
  return(this.tran)
}))
all_combos$pcol_tran = "white"
all_combos$pcol_tran[non_na_rows] = paste0(all_combos$pcol[non_na_rows], all_combos$tran[non_na_rows])
all_combos$sig_level = plyr::revalue(as.character(all_combos$tran), replace = c("FALSE" = "FALSE", "FF" = "TRUE", "85" = "trending"))
all_combos$sig = F
all_combos$sig[which(all_combos$sig_level == "TRUE")] = T
all_combos$trending = F
all_combos$trending[which(all_combos$sig_level == "trending")] = T
all_combos$pcol_tran2 = "white"
all_combos$pcol_tran2[non_na_rows] = paste0(all_combos$pcol[non_na_rows], "85")

all_combos$cluster = factor(all_combos$cluster, levels = rev(convert_all$cluster))
# all_combos$hgnc = factor(all_combos$hgnc, levels = c("All", sort(unique(all_combos$hgnc[which(all_combos$hgnc != "All")]))))

ggplot(all_combos, aes(x = hgnc, y = cluster, fill = pcol)) + geom_tile(color = "gray60") + scale_fill_identity() + coord_fixed() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic"), axis.text.y = element_text(colour = rev(convert_all$color), face=ifelse(rev(convert_all$level) =="secondary","plain","bold"), size=ifelse(rev(convert_all$level) =="secondary", 8, 10))) + xlab("") + ylab("") + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0, 0))

num_hits_by_clust = data.frame(table(all_combos$cluster[which(all_combos$pcol != "white" & all_combos$cluster != "All")]))
num_hits_by_gp    = data.frame(table(all_combos$gene_pop[which(all_combos$pcol != "white" & all_combos$cluster != "All")]))
convert_all$num_hits = num_hits_by_clust$Freq[match(convert_all$cluster, num_hits_by_clust$Var1)]
convert_all_small = convert_all[which(convert_all$num_hits != 0 & convert_all$cluster != "All"),]
all_combos_small = all_combos[which(all_combos$cluster %in% convert_all_small$cluster),]
all_combos_small$cluster = factor(all_combos_small$cluster, levels = rev(convert_all_small$cluster))
all_combos_small$gp_hits = num_hits_by_gp$Freq[match(all_combos_small$gene_pop, num_hits_by_gp$Var1)]
all_combos_small = all_combos_small[which(all_combos_small$gp_hits > 0),]
length(which(! all_combos_small$pcol %in% c("white", "#FDE725", "#2AB07F", "#433E85")))
all_combos_small[which(! all_combos_small$pcol %in% c("white", "#FDE725", "#2AB07F", "#433E85")),]

# # all_combos_small$num = as.numeric(as.vector(plyr::revalue(all_combos_small$pcol, replace = c("white" = 0, "#FDE725" = 1, "#2AB07F" = 2, "#433E85" = 3, "#94CC52" = 1.5, "#A09355" = 1, "#377782" = 2.5, "#799c63" = 4))))
# all_combos_small$num = as.numeric(as.vector(plyr::revalue(all_combos_small$pcol, replace = c("white" = 0, "#FDE725" = 1.1, "#2AB07F" = 1.2, "#433E85" = 1.3, "#94CC52" = 1.15, "#A09355" = 1.2, "#377782" = 1.25, "#799c63" = 1.4))))
# combos_mat = reshape(all_combos_small[,c("cluster", "gene_pop", "num")], idvar = "gene_pop", timevar = "cluster", direction = "wide")
# rownames(combos_mat) = combos_mat[,1]
# combos_mat[,1] = NULL
# colnames(combos_mat) = str_sub(colnames(combos_mat), 5, 10000L)
# combos_mat = scale(combos_mat)
# dist_mat <- dist(combos_mat, method = 'euclidean')
# hclust_avg <- hclust(dist_mat, method = 'average')
# gene_order = rownames(combos_mat)[hclust_avg$order]
# all_combos_small$gene_pop = factor(all_combos_small$gene_pop, levels = gene_order)

order_combos = data.frame(table(all_combos_small$gene_pop[which(all_combos_small$pcol != "white")]))
order_combos$border = data.frame(table(all_combos_small$gene_pop[which( ! is.na(all_combos_small$bsig) )]))[,2]
order_combos$gorder = data.frame(table(all_combos_small$gene_pop[which( ! is.na(all_combos_small$gsig) )]))[,2]
order_combos$qorder = data.frame(table(all_combos_small$gene_pop[which( ! is.na(all_combos_small$qsig) )]))[,2]
order_combos = order_combos[order(order_combos$qorder, order_combos$gorder, order_combos$border, decreasing = T),]
all_combos_small$gene_pop = factor(all_combos_small$gene_pop, levels = order_combos$Var1)

all_combos_small$hgnc = tolower(gene_info$human[match(all_combos_small$gene_pop, gene_info$mzebra)])
all_combos_small$hgnc[which(all_combos_small$gene_pop == "LOC101474236")] = "LOC101474236"
all_combos_small$hgnc[which(all_combos_small$gene_pop == "All")] = "All"
all_combos_small$hgnc = factor(all_combos_small$hgnc, levels = all_combos_small$hgnc[match(order_combos$Var1, all_combos_small$gene_pop)])

pdf("C:/Users/miles/Downloads/ieg_summary_3_sig.pdf", width = 7, height = 5)
ggplot(all_combos_small, aes(x = hgnc, y = cluster, fill = pcol_tran)) + geom_tile(color = "gray60") + scale_fill_identity() + coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic"), axis.text.y = element_text(colour = rev(convert_all_small$color), face=ifelse(rev(convert_all_small$level) =="secondary","plain","bold"), size=ifelse(rev(convert_all_small$level) =="secondary", 8, 10))) + xlab("") + ylab("") + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0, 0))
dev.off()
# all_combos_small$pcol[which(all_combos_small$pcol == "white")] = "gray60"
# ggplot(all_combos_small, aes(x = gene_pop, y = cluster, fill = pcol_tran, color = pcol)) + geom_tile(aes(size = sig), width = 0.8, height = 0.8) + scale_color_identity() + scale_fill_identity() + scale_size_manual(values = c(0.8, 1.4)) + coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic"), axis.text.y = element_text(colour = rev(convert_all_small$color), face=ifelse(rev(convert_all_small$level) =="secondary","plain","bold"), size=ifelse(rev(convert_all_small$level) =="secondary", 8, 10))) + xlab("") + ylab("") + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0, 0))
# all_combos_small$pcol[which(all_combos_small$pcol == "white")] = "gray60"
# all_combos_small$pcol_tran[which(all_combos_small$pcol_tran == "white")] = "gray90"
# all_combos_small = all_combos_small[order(all_combos_small$trending, decreasing = F),]
pdf("C:/Users/miles/Downloads/ieg_summary_1.pdf", width = 7, height = 5)
ggplot(all_combos_small, aes(x = hgnc, y = cluster, fill = pcol_tran, color = pcol)) + geom_tile(width = 0.8, height = 0.8, size = 0.725) + scale_color_identity() + scale_fill_identity() + coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic"), axis.text.y = element_text(colour = rev(convert_all_small$color), face=ifelse(rev(convert_all_small$level) =="secondary","plain","bold"), size=ifelse(rev(convert_all_small$level) =="secondary", 8, 10))) + xlab("") + ylab("")
dev.off()
# all_combos_small$pcol[which(all_combos_small$pcol == "white")] = "gray60"
# all_combos_small = all_combos_small[order(all_combos_small$sig, decreasing = F),]
# ggplot(all_combos_small, aes(x = gene_pop, y = cluster, fill = pcol_tran2)) + geom_tile(color = "gray60") + geom_tile(data = all_combos_small[which(all_combos_small$sig),], aes(color = pcol), size = 1.5) + scale_color_identity() + scale_fill_identity() + coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic"), axis.text.y = element_text(colour = rev(convert_all_small$color), face=ifelse(rev(convert_all_small$level) =="secondary","plain","bold"), size=ifelse(rev(convert_all_small$level) =="secondary", 8, 10))) + xlab("") + ylab("")
all_combos_small = all_combos_small[order(all_combos_small$trending, decreasing = F),]
pdf("C:/Users/miles/Downloads/ieg_summary_3.pdf", width = 7, height = 5)
ggplot(all_combos_small, aes(x = hgnc, y = cluster, fill = pcol_tran, color = pcol)) + geom_tile(aes(size = trending)) + scale_color_identity() + scale_fill_identity() + scale_size_manual(values = c(0.6, 1), guide = F) + coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic"), axis.text.y = element_text(colour = rev(convert_all_small$color), face=ifelse(rev(convert_all_small$level) =="secondary","plain","bold"), size=ifelse(rev(convert_all_small$level) =="secondary", 8, 10))) + xlab("") + ylab("")
dev.off()
# ggplot(all_combos_small, aes(x = gene_pop, y = cluster, fill = pcol_tran)) + geom_tile() + geom_tile(data = all_combos_small[which(all_combos_small$trending),], aes(color = pcol), size = 0.8) + scale_color_identity() + scale_fill_identity() + coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic"), axis.text.y = element_text(colour = rev(convert_all_small$color), face=ifelse(rev(convert_all_small$level) =="secondary","plain","bold"), size=ifelse(rev(convert_all_small$level) =="secondary", 8, 10))) + xlab("") + ylab("")

# Clownfish Count - Fisher's Exact Test
Idents(clown) = clown$clusters49
fpdf = data.frame()
for (clust in sort(unique(Idents(clown)))) {
  f_cells = length(which(clown$clusters49 == clust & clown$sex == "f"))
  m_cells = length(which(clown$clusters49 == clust & clown$sex == "m"))
  other_f = length(which(clown$clusters49 != clust & clown$sex == "f"))
  other_m = length(which(clown$clusters49 != clust & clown$sex == "m"))
  contig_table = data.frame(m = c(m_cells, other_m), f = c(f_cells, other_f))
  fpdf = rbind(fpdf, data.frame(cluster = clust, fp = fisher.test(contig_table)$p.value))
}


z15 = read.csv("~/scratch/brain/results/out_bb15_bbmm_demux_deg_all_tests_for_volcano_plotting_121321.csv")
z53 = read.csv("~/scratch/brain/results/out_bb53_glmmseq_demux_deg_all_tests_for_volcano_plotting.csv")
cur_level = "53"
if  (cur_level == "53") { cz = z53 }                   else { cz = z15 }
if  (cur_level == "53") { cmeta = "seuratclusters53" } else { cmeta = "seuratclusters15" }
if  (cur_level == "53") { clusters = 0:52 }            else { clusters = 0:14 }
all_zgenes = unique(cz$mzebra)
all_zgenes = str_replace(all_zgenes, pattern = "\\.", "-")
all_zgenes = all_zgenes[which(all_zgenes %in% rownames(bb))]
# gcm = sparseMatrix(i = integer(0), j = integer(0), dims = c(length(all_zgenes), ncol(bb)), dimnames = list(all_zgenes, colnames(bb)))
# gcm = as(gcm, "dgCMatrix")
gcm_list = list()
for (i in clusters) {
  print(i)
  # zgenes = cz$mzebra[which(cz$cluster == i)]
  clust_idx = which(bb@meta.data[,cmeta] == i)
  this_cz = cz[which(cz$cluster == i),]
  this_cz$zgenes = this_cz$mzebra
  print(paste0(length(which(! this_cz$zgenes %in% rownames(bb))), " genes not in bb without correction."))
  this_cz$zgenes = str_replace(this_cz$zgenes, pattern = "\\.", "-")
  print(paste0(length(which(! this_cz$zgenes %in% rownames(bb))), " genes not in bb at first pass."))
  print(head(this_cz$zgenes)[which(! this_cz$zgenes %in% rownames(bb))])
  this_cz = this_cz[which(this_cz$zgenes %in% rownames(bb)),]
  glmmseq_counts = bb@assays$RNA@counts[this_cz$zgenes, clust_idx]
  
  if (nrow(this_cz) != 0) {
    gsi_cell = bb$gsi[clust_idx]
    spawn_cell = bb$log_spawn_events[clust_idx]
    new_mods = mclapply(1:nrow(this_cz), function(x) adjGlmmseqCounts(x), mc.cores = 24)
    new_mods_mtx = do.call(rbind, new_mods)
    rownames(new_mods_mtx) = this_cz$zgenes
    zgenes_not_in_cluster = all_zgenes[which(! all_zgenes %in% this_cz$zgenes )]
    extra_rows = matrix(data = 0L, nrow = length(zgenes_not_in_cluster), ncol = length(clust_idx))
    rownames(extra_rows) = zgenes_not_in_cluster
    colnames(extra_rows) = colnames(bb)[clust_idx]
    new_mods_mtx = rbind(new_mods_mtx, extra_rows)
    new_mods_mtx = new_mods_mtx[ order(row.names(new_mods_mtx)), ]
    new_mods_mtx = as(new_mods_mtx, "dgCMatrix")
  }
  else {
    zgenes_not_in_cluster = all_zgenes[which(! all_zgenes %in% this_cz$zgenes )]
    extra_rows = matrix(data = 0L, nrow = length(zgenes_not_in_cluster), ncol = length(clust_idx))
    rownames(extra_rows) = zgenes_not_in_cluster
    colnames(extra_rows) = colnames(bb)[clust_idx]
    new_mods_mtx = extra_rows
    new_mods_mtx = new_mods_mtx[ order(row.names(new_mods_mtx)), ]
  }
  gcm_list[[paste0("cluster_", i)]] = new_mods_mtx
}

gcm = do.call(cbind, gcm_list)
gcm = gcm[,colnames(bb)]
saveRDS(gcm, "~/scratch/brain/results/adjusted_glmmseq_ffm_53.rds")
# saveRDS(gcm, "~/scratch/brain/results/adjusted_glmmseq_ffm_15.rds")

adjGlmmseqCounts = function(x) {
  zg = this_cz$zgenes[x]
  bg = this_cz$mzebra[x]
  my_gsi   = this_cz$gsi[x]
  my_spawn = this_cz$log_spawn_events[x]
  
  return(glmmseq_counts[zg, ] - my_gsi * gsi_cell - my_spawn * spawn_cell)
}

gene = "wdr17"
gcm_df = data.frame(adj_counts = gcm53[gene,which(bb$seuratclusters53 == 0)], sample = bb$sample[which(bb$seuratclusters53 == 0)], cond = bb$cond[which(bb$seuratclusters53 == 0)], subsample = bb$subsample[which(bb$seuratclusters53 == 0)], pair = bb$pair[which(bb$seuratclusters53 == 0)], bai = bb$bower_activity_index[which(bb$seuratclusters53 == 0)], gsi = bb$gsi[which(bb$seuratclusters53 == 0)], spawn = bb$log_spawn_events[which(bb$seuratclusters53 == 0)])
colnames(gcm_df)[1] = "adj_counts"
gcm_df_agr = aggregate(adj_counts ~ subsample + cond + pair + bai, gcm_df, mean)

pdf("~/scratch/brain/results/bb53_0_wdr17_glmmseq_test.pdf", width = 4, height = 4)
print(ggplot(gcm_df_agr, aes(x = cond, y = adj_counts, color = cond)) + geom_boxplot(outlier.shape = NA) + geom_line(linetype = "dashed", color = "gray40", alpha = 0.4, aes(group=pair)) + geom_point(alpha = 0.75, position = position_jitterdodge()))
dev.off()

if  (cur_level == "53") { gcm = gcm53 } else { gcm = gcm15 }
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


# 9 Plot DEGs =====================================================
bb15 = read.csv("C:/Users/miles/Downloads/bb15_deg_all_split_by_up_or_down_121621.csv")
bb53 = read.csv("C:/Users/miles/Downloads/bb53_deg_all_split_by_up_or_down_121621.csv")
bb15$sig_any = bb15$sig_bower_behavior == 1 | bb15$sig_gsi == 1 | bb15$sig_log_spawn_events == 1
bb53$sig_any = bb53$sig_bower_behavior == 1 | bb53$sig_gsi == 1 | bb53$sig_log_spawn_events == 1
bb15 = bb15[which(bb15$sig_any & bb15$mzebra %in% rownames(bb)),]
bb53 = bb53[which(bb53$sig_any & bb53$mzebra %in% rownames(bb)),]
# all_bdeg = unique(c(bb15$mzebra[which(bb15$sig_bower_behavior == 1)], bb53$mzebra[which(bb53$sig_bower_behavior == 1)]))
# all_gdeg = unique(c(bb15$mzebra[which(bb15$sig_gsi == 1)], bb53$mzebra[which(bb53$sig_gsi == 1)]))
# all_qdeg = unique(c(bb15$mzebra[which(bb15$sig_log_spawn_events == 1)], bb53$mzebra[which(bb53$sig_log_spawn_events == 1)]))
# all_deg = data.frame(table(c(all_bdeg, all_gdeg, all_qdeg)))
# all_deg = as.vector(all_deg$Var1[which(all_deg$Freq == 3)])
all_deg = read.csv("C:/Users/miles/Downloads/overlap_of_bDEGs_qDEGs_gDEGs_012021.csv")[,5]
ere = read.csv("C:/Users/miles/Downloads/ere_class.csv")[,2]
my.pt.size = 0.6
p_list = list()
for (cat in c("bower_behavior", "gsi", "log_spawn_events")) {
  print(cat)
  if (cat == "bower_behavior") { cat_mod = "bower_activity_index" } else { cat_mod = cat }
  this_bb15 = bb15[which(bb15[, paste0("sig_", cat)] == 1),]
  this_bb53 = bb53[which(bb53[, paste0("sig_", cat)] == 1),]
  # this_bb15_up = bb15[which(bb15[, cat_mod] > 0)]
  # this_bb15_down = bb15[which(bb15[, cat_mod] < 0)]
  score_df = data.frame(cell = colnames(bb), UMAP_1 = bb@reductions$umap@cell.embeddings[,"UMAP_1"], UMAP_2 = bb@reductions$umap@cell.embeddings[,"UMAP_2"], up_score = 0, down_score = 0, all_score = 0, ovlp_score = 0, ere_score = 0)
  for (i in 1:nrow(this_bb15)) {
    this_gene = this_bb15$mzebra[i]
    this_cluster = this_bb15$cluster[i]
    isUp = this_bb15[i, cat_mod] > 0
    this_score = as.numeric(bb@assays$RNA@counts[this_gene,] > 0)
    score_df$all_score[which(bb$seuratclusters15 == this_cluster)] = score_df$all_score[which(bb$seuratclusters15 == this_cluster)] + this_score[which(bb$seuratclusters15 == this_cluster)]
    if (isUp) { score_df$up_score[which(bb$seuratclusters15 == this_cluster)] = score_df$up_score[which(bb$seuratclusters15 == this_cluster)] + this_score[which(bb$seuratclusters15 == this_cluster)] } else { score_df$down_score[which(bb$seuratclusters15 == this_cluster)] = score_df$down_score[which(bb$seuratclusters15 == this_cluster)] + this_score[which(bb$seuratclusters15 == this_cluster)] }
    if (this_gene %in% all_deg) { score_df$ovlp_score[which(bb$seuratclusters15 == this_cluster)] = score_df$ovlp_score[which(bb$seuratclusters15 == this_cluster)] + this_score[which(bb$seuratclusters15 == this_cluster)] }
    if (this_gene %in% ere) { score_df$ere_score[which(bb$seuratclusters15 == this_cluster)] = score_df$ere_score[which(bb$seuratclusters15 == this_cluster)] + this_score[which(bb$seuratclusters15 == this_cluster)] }
  }
  for (i in 1:nrow(this_bb53)) {
    this_gene = this_bb53$mzebra[i]
    this_cluster = this_bb53$cluster[i]
    isUp = this_bb53[i, cat_mod] > 0
    this_score = as.numeric(bb@assays$RNA@counts[this_gene,] > 0)
    score_df$all_score[which(bb$seuratclusters53 == this_cluster)] = score_df$all_score[which(bb$seuratclusters53 == this_cluster)] + this_score[which(bb$seuratclusters53 == this_cluster)]
    if (isUp) { score_df$up_score[which(bb$seuratclusters53 == this_cluster)] = score_df$up_score[which(bb$seuratclusters53 == this_cluster)] + this_score[which(bb$seuratclusters53 == this_cluster)] } else { score_df$down_score[which(bb$seuratclusters53 == this_cluster)] = score_df$down_score[which(bb$seuratclusters53 == this_cluster)] + this_score[which(bb$seuratclusters53 == this_cluster)] }
    if (this_gene %in% all_deg) { score_df$ovlp_score[which(bb$seuratclusters53 == this_cluster)] = score_df$ovlp_score[which(bb$seuratclusters53 == this_cluster)] + this_score[which(bb$seuratclusters53 == this_cluster)] }
    if (this_gene %in% ere) { score_df$ere_score[which(bb$seuratclusters15 == this_cluster)] = score_df$ere_score[which(bb$seuratclusters15 == this_cluster)] + this_score[which(bb$seuratclusters15 == this_cluster)] }
  }
  score_df = score_df[order(score_df$up_score, decreasing = F),]
  # score_df$up_score[which(score_df$up_score == 0)] = NA
  p_list[[paste0(cat, "_", "up")]] = ggplot(score_df, aes(x = UMAP_1, y = UMAP_2, color = up_score)) + geom_point(size = my.pt.size) + scale_color_gradientn(colors = viridis(100), na.value = "grey80") + theme_void() + NoLegend()
  score_df = score_df[order(score_df$down_score, decreasing = F),]
  # score_df$down_score[which(score_df$down_score == 0)] = NA
  p_list[[paste0(cat, "_", "down")]] = ggplot(score_df, aes(x = UMAP_1, y = UMAP_2, color = all_score)) + geom_point(size = my.pt.size) + scale_color_gradientn(colors = viridis(100), na.value = "grey80") + theme_void() + NoLegend()
  score_df = score_df[order(score_df$all_score, decreasing = F),]
  # score_df$all_score[which(score_df$all_score == 0)] = NA
  p_list[[paste0(cat, "_", "all")]] = ggplot(score_df, aes(x = UMAP_1, y = UMAP_2, color = all_score)) + geom_point(size = my.pt.size) + scale_color_gradientn(colors = viridis(100), na.value = "grey80") + theme_void() + NoLegend()
  score_df = score_df[order(score_df$ovlp_score, decreasing = F),]
  # score_df$ovlp_score[which(score_df$ovlp_score == 0)] = NA
  p_list[[paste0(cat, "_", "ovlp")]] = ggplot(score_df, aes(x = UMAP_1, y = UMAP_2, color = ovlp_score)) + geom_point(size = my.pt.size) + scale_color_gradientn(colors = viridis(100), na.value = "grey80") + theme_void() + NoLegend()
  p_list[[paste0(cat, "_", "ere")]] = ggplot(score_df, aes(x = UMAP_1, y = UMAP_2, color = ere_score)) + geom_point(size = my.pt.size) + scale_color_gradientn(colors = viridis(100), na.value = "grey80") + theme_void() + NoLegend()
}

pdf("C:/Users/miles/Downloads/zack_15_deg_plot.pdf", width = 6*5, height = 6*3)
p = plot_grid(plotlist=p_list, ncol = 5)
print(p)
dev.off()

png("C:/Users/miles/Downloads/zack_15_deg_plot.png", width = 350*6, height = 350*3)
p = plot_grid(plotlist=p_list, ncol = 5)
print(p)
dev.off()

# print(ggplot(pdf, aes(x = time, y = value)) + geom_point(data = pdf[which(!pdf$isMean),], size = 2.5, alpha = 0.2, aes(color = variable)) + geom_point(data = pdf[which(pdf$isMean),], size = 2.5, color = "black") + geom_smooth(data = pdf[which(pdf$isMean),], method = "loess", se = F, color = "gray40") + theme_classic() + ylab("R2") + xlab("Time (min to flash freeze)") + ggtitle(paste0("bDEG Hits Up at ", i_clean, ". Depth_adj R2 w/ Adjusted")) + scale_x_continuous(breaks = rev(unique(pdf$time)), labels = rev(unique(pdf$time))) + NoLegend())

sub_meta = aggregate(depth_5_35 + depth_15_45 + depth_25_55 + depth_35_65 + depth_45_75 + depth_55_85 + depth_65_95 + build_5_35 + build_15_45 + build_25_55 + build_35_65 + build_45_75 + build_55_85 + build_65_95 ~ subsample, bb@meta.data, mean)
mean_df = read.csv("~/Downloads/ieg_summary_subsample_means_bower.csv")

#***************************************************
# Neurogen Summary Figure ==========================
#***************************************************
neurogen15_s = read.csv("C:/Users/miles/Downloads/neurogen_all_analyses_hmp_output_for_george_022522/neurogen_all_analyses_hmp_output_for_george_022522/out_neurogen_by_goi_by_cluster_bbmm_demux_log_spawn_events_hmp_calculated_across_all_goi_020122_hgnc.csv")
neurogen15_g = read.csv("C:/Users/miles/Downloads/neurogen_all_analyses_hmp_output_for_george_022522/neurogen_all_analyses_hmp_output_for_george_022522/out_neurogen_by_goi_by_cluster_bbmm_demux_gsi_hmp_calculated_across_all_goi_020122_hgnc.csv")
neurogen15_b = read.csv("C:/Users/miles/Downloads/neurogen_all_analyses_hmp_output_for_george_022522/neurogen_all_analyses_hmp_output_for_george_022522/out_neurogen_by_goi_by_cluster_bbmm_demux_bower_behavior_hmp_calculated_across_all_goi_020122_hgnc.csv")
neurogen53_s = read.csv("C:/Users/miles/Downloads/neurogen_all_analyses_hmp_output_for_george_022522/neurogen_all_analyses_hmp_output_for_george_022522/out_neurogen_by_goi_by_53cluster_bbmm_demux_log_spawn_events_hmp_calculated_across_all_goi_020122_hgnc.csv")
neurogen53_g = read.csv("C:/Users/miles/Downloads/neurogen_all_analyses_hmp_output_for_george_022522/neurogen_all_analyses_hmp_output_for_george_022522/out_neurogen_by_goi_by_53cluster_bbmm_demux_gsi_hmp_calculated_across_all_goi_020122_hgnc.csv")
neurogen53_b = read.csv("C:/Users/miles/Downloads/neurogen_all_analyses_hmp_output_for_george_022522/neurogen_all_analyses_hmp_output_for_george_022522/out_neurogen_by_goi_by_53cluster_bbmm_demux_bower_behavior_hmp_calculated_across_all_goi_020122_hgnc.csv")

neurogen15_s$cluster_gene = paste0(neurogen15_s$cluster, "_", neurogen15_s$goi)
neurogen15_g$cluster_gene = paste0(neurogen15_g$cluster, "_", neurogen15_g$goi)
neurogen15_b$cluster_gene = paste0(neurogen15_b$cluster, "_", neurogen15_b$goi)
neurogen53_s$cluster_gene = paste0(neurogen53_s$cluster, "_", neurogen53_s$goi)
neurogen53_g$cluster_gene = paste0(neurogen53_g$cluster, "_", neurogen53_g$goi)
neurogen53_b$cluster_gene = paste0(neurogen53_b$cluster, "_", neurogen53_b$goi)

neurogen15_s$sig_s = neurogen15_s$sig_log_spawn_events == 5 & neurogen15_s$hmp < 0.05
neurogen15_g$sig_g = neurogen15_g$sig_gsi == 5 & neurogen15_g$hmp < 0.05
neurogen15_b$sig_b = neurogen15_b$sig_cond == 3 & neurogen15_b$hmp_cond < 0.05
neurogen53_s$sig_s = neurogen53_s$sig_log_spawn_events == 5 & neurogen53_s$hmp < 0.05
neurogen53_g$sig_g = neurogen53_g$sig_gsi == 5 & neurogen53_g$hmp < 0.05
neurogen53_b$sig_b = neurogen53_b$sig_cond == 3 & neurogen53_b$hmp_cond < 0.05

neurogen53_b$sig_g = neurogen53_g$sig_g[match(neurogen53_b$cluster_gene, neurogen53_g$cluster_gene)]
neurogen53_b$sig_s = neurogen53_s$sig_s[match(neurogen53_b$cluster_gene, neurogen53_s$cluster_gene)]
neurogen15_b$sig_g = neurogen15_g$sig_g[match(neurogen15_b$cluster_gene, neurogen15_g$cluster_gene)]
neurogen15_b$sig_s = neurogen15_s$sig_s[match(neurogen15_b$cluster_gene, neurogen15_s$cluster_gene)]

neurogen15_b$level_old = paste0("primary_", neurogen15_b$cluster)
neurogen15_b$cluster_all = convert_all$cluster[match(neurogen15_b$level_old, convert_all$level_old)]
neurogen53_b$level_old = paste0("secondary_", neurogen53_b$cluster)
neurogen53_b$cluster_all = convert_all$cluster[match(neurogen53_b$level_old, convert_all$level_old)]

neurogen53_b$pcol = unlist(lapply(1:nrow(neurogen53_b), function(x) {
  this.col = "white"
  if ( neurogen53_b[x,"sig_b"] ) { this.col = "#FDE725" }
  if ( neurogen53_b[x,"sig_g"] ) { this.col = "#2AB07F" }
  if ( neurogen53_b[x,"sig_s"] ) { this.col = "#433E85" }
  if ( neurogen53_b[x,"sig_b"] & neurogen53_b[x,"sig_g"] ) { this.col = "#94CC52" }
  if ( neurogen53_b[x,"sig_b"] & neurogen53_b[x,"sig_s"] ) { this.col = "#A09355" }
  if ( neurogen53_b[x,"sig_g"] & neurogen53_b[x,"sig_s"] ) { this.col = "#377782" }
  if ( neurogen53_b[x,"sig_b"] & neurogen53_b[x,"sig_g"] & neurogen53_b[x,"sig_s"] ) { this.col = "#799c63" }
  return(this.col)
}))
neurogen15_b$pcol = unlist(lapply(1:nrow(neurogen15_b), function(x) {
  this.col = "white"
  if ( neurogen15_b[x,"sig_b"] ) { this.col = "#FDE725" }
  if ( neurogen15_b[x,"sig_g"] ) { this.col = "#2AB07F" }
  if ( neurogen15_b[x,"sig_s"] ) { this.col = "#433E85" }
  if ( neurogen15_b[x,"sig_b"] & neurogen15_b[x,"sig_g"] ) { this.col = "#94CC52" }
  if ( neurogen15_b[x,"sig_b"] & neurogen15_b[x,"sig_s"] ) { this.col = "#A09355" }
  if ( neurogen15_b[x,"sig_g"] & neurogen15_b[x,"sig_s"] ) { this.col = "#377782" }
  if ( neurogen15_b[x,"sig_b"] & neurogen15_b[x,"sig_g"] & neurogen15_b[x,"sig_s"] ) { this.col = "#799c63" }
  return(this.col)
}))

cols_to_keep = c("pcol", "level_old", "cluster_all", "mzebra", "sig_b", "sig_g", "sig_s")
neurogen_all = rbind(neurogen15_b[,cols_to_keep], neurogen53_b[,cols_to_keep])
neurogen_all = neurogen_all[which(neurogen_all$pcol != "white"),]

all <- neurogen_all %>% expand(cluster_all, mzebra)
all$cluster_all_mzebra = paste0(all$cluster_all, "_", all$mzebra)
all$pcol = "white"
match_idx = which(all$cluster_all_mzebra %in% paste0(neurogen_all$cluster_all, "_", neurogen_all$mzebra))
all[match_idx,cols_to_keep] = neurogen_all[match(all$cluster_all_mzebra[match_idx], paste0(neurogen_all$cluster_all, "_", neurogen_all$mzebra)), cols_to_keep]

order_combos = data.frame(table(all$mzebra[which(all$pcol != "white")]))
order_combos$border = data.frame(table(all$mzebra[which( all$sig_b )]))[,2]
order_combos$gorder = data.frame(table(all$mzebra[which( all$sig_g )]))[,2]
order_combos$qorder = data.frame(table(all$mzebra[which( all$sig_s )]))[,2]
order_combos = order_combos[order(order_combos$qorder, order_combos$gorder, order_combos$border, decreasing = T),]

all$mzebra = factor(all$mzebra, levels = order_combos$Var1)
all$hgnc = gene_info$human[match(all$mzebra, gene_info$mzebra)]
all$label = as.vector(all$mzebra)
all$label[which(startsWith(all$label, "LOC"))] = paste0(all$mzebra[which(startsWith(all$label, "LOC"))], " (", tolower(all$hgnc[which(startsWith(all$label, "LOC"))]), ")")
all$label = factor(all$label, levels = unique(all$label[order(all$mzebra)]))
all = all[which(! is.na(all$cluster_all) ),]
all$cluster_all = factor(all$cluster_all, levels = convert_all$cluster)

pdf("C:/Users/miles/Downloads/neurogen_sum.pdf", width = 14, height  = 6)
# ggplot(all, aes(x = label, y = cluster_all, fill = pcol)) + geom_tile(color = "gray60") + scale_fill_identity() + coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic"), axis.text.y = element_text(colour = rev(convert_all$color[which(convert_all$cluster %in% all$cluster_all)]), face=ifelse(rev(convert_all$level[which(convert_all$cluster %in% all$cluster_all)]) =="secondary","plain","bold"), size=ifelse(rev(convert_all$level[which(convert_all$cluster %in% all$cluster_all)]) =="secondary", 8, 10))) + xlab("") + ylab("") + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0, 0)) + coord_fixed()
ggplot(all, aes(x = label, y = cluster_all, fill = pcol)) + geom_tile(color = "gray60") + scale_fill_identity() + coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic"), axis.text.y = element_text(colour = convert_all$color[which(convert_all$cluster %in% all$cluster_all)], face=ifelse(convert_all$level[which(convert_all$cluster %in% all$cluster_all)] =="secondary","plain","bold"), size=ifelse(convert_all$level[which(convert_all$cluster %in% all$cluster_all)] =="secondary", 8, 10))) + xlab("") + ylab("") + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0, 0)) + coord_fixed()
dev.off()

# 2
neurogen_sum = read.csv("C:/Users/miles/Downloads/bower_quiver_gsi_cluster_goi_beta_effects_for_dotplot_030122.csv")
neurogen_sum$is_sig = neurogen_sum$plot
neurogen_sum = neurogen_sum[which(neurogen_sum$is_sig),]
neurogen_sum$hgnc = neurogen_sum$gene
neurogen_sum$gene_pop = neurogen_sum$goi
neurogen_sum$cluster[which(neurogen_sum$cluster == FALSE)] = "All"
neurogen_sum$gene_pop[which(neurogen_sum$gene_pop == FALSE)] = "All"
neurogen_sum$level_old = paste0(neurogen_sum$level, "_", neurogen_sum$cluster)
neurogen_sum$level_old_gp = paste0(neurogen_sum$level_old, "_", neurogen_sum$gene_pop)
neurogen_sum$cat_level_old_gp = paste0(neurogen_sum$cat, "_", neurogen_sum$level_old_gp)

convert_all = data.frame(cluster = c(convert15$new.full, convert53$new), color = c(convert15$col, convert53$col), old = c(convert15$old, convert53$old))
convert_all = convert_all[which(! duplicated(convert_all$cluster) ),]
convert_all[, c("new.id", "new.gaba")]    = reshape2::colsplit(convert_all$cluster, "_", names = c("new.id", "new.gaba"))
convert_all[, c("new.parent", "new.sub")] = reshape2::colsplit(convert_all$new.id, "\\.", names = c("new.id", "new.gaba"))
convert_all$new.sub[which(is.na(convert_all$new.sub))] = 0
convert_all[which(convert_all$cluster == "8-9_Glut"), c("new.parent", "new.sub")] = c(8, 12) # Assign a parent and subcluster to the 8-9_Glut special case
convert_all = convert_all[order(as.numeric(convert_all$new.parent), as.numeric(convert_all$new.sub), decreasing = F),]
convert_all$level = plyr::revalue(as.character(convert_all$new.sub == 0), replace = c("TRUE" = "primary", "FALSE" = "secondary"))
convert_all = rbind(data.frame(cluster = "All", color = viridis::viridis(1), new.id = 0, new.gaba = "all", new.parent = 0, new.sub = 0, level = "all", old = "All"), convert_all)
convert_all$level_old = paste0(convert_all$level, "_", convert_all$old)
convert_all$cluster = stringr::str_replace(convert_all$cluster, "Astro", "RGC")

all_combos = expand.grid(convert_all$level_old, unique(neurogen_sum$gene_pop[which(neurogen_sum$gene_pop != FALSE)]) )
colnames(all_combos) = c("level_old", "gene_pop")
all_combos$level_old_gp = paste0(all_combos$level_old, "_", all_combos$gene_pop)
all_combos[,colnames(convert_all)] = convert_all[match(all_combos$level_old, convert_all$level_old),]
all_combos$hgnc = gene_info$human[match(all_combos$gene_pop, gene_info$mzebra)]
all_combos$hgnc[which(is.na(all_combos$hgnc))] = "All"
all_combos$bsig = NA
all_combos$gsig = NA
all_combos$qsig = NA

for (i in 1:nrow(neurogen_sum)) {
  my.sig = neurogen_sum$is_sig[i]
  my.cat = neurogen_sum$cat[i]
  my.x = neurogen_sum$level_old_gp[i]
  my.idx = which(all_combos$level_old_gp == my.x)
  my.cat.col = plyr::revalue(my.cat, replace = c("bower" = "bsig", "gsi" = "gsig", "quiver" = "qsig"))
  if (! is.na(my.cat.col) )
    all_combos[my.idx, my.cat.col] = my.sig
  if (my.x == "secondary_36_bhlhe22")
    all_combos[my.idx, "bsig"] = T
}

all_combos$pcol = "white"
non_na_rows = which( ! is.na(all_combos$bsig) | ! is.na(all_combos$gsig) | ! is.na(all_combos$qsig) )
all_combos$pcol[non_na_rows] = unlist(lapply(non_na_rows, function(x) {
  this.col = "error"
  if ( ! is.na(all_combos[x,"bsig"]) ) { this.col = "#FDE725" }
  if ( ! is.na(all_combos[x,"gsig"]) ) { this.col = "#2AB07F" }
  if ( ! is.na(all_combos[x,"qsig"]) ) { this.col = "#433E85" }
  if ( ! is.na(all_combos[x,"bsig"]) & ! is.na(all_combos[x,"gsig"]) ) { this.col = "#94CC52" }
  if ( ! is.na(all_combos[x,"bsig"]) & ! is.na(all_combos[x,"qsig"]) ) { this.col = "#A09355" }
  if ( ! is.na(all_combos[x,"gsig"]) & ! is.na(all_combos[x,"qsig"]) ) { this.col = "#377782" }
  if ( ! is.na(all_combos[x,"bsig"]) & ! is.na(all_combos[x,"gsig"]) & ! is.na(all_combos[x,"qsig"]) ) { this.col = "#799c63" }
  return(this.col)
}))

all_combos$tran = "FALSE"
non_na_rows = which( ! is.na(all_combos$bsig) | ! is.na(all_combos$gsig) | ! is.na(all_combos$qsig) )
all_combos$tran[non_na_rows] = unlist(lapply(non_na_rows, function(x) {
  this.tran = "FF"
  if ( ! is.na(all_combos[x,"bsig"]) &  ! all_combos[x,"bsig"] ) { this.tran = "85" }
  if ( ! is.na(all_combos[x,"gsig"]) &  ! all_combos[x,"gsig"] ) { this.tran = "85" }
  if ( ! is.na(all_combos[x,"qsig"]) &  ! all_combos[x,"qsig"] ) { this.tran = "85" }
  return(this.tran)
}))
all_combos$pcol_tran = "white"
all_combos$pcol_tran[non_na_rows] = paste0(all_combos$pcol[non_na_rows], all_combos$tran[non_na_rows])
all_combos$sig_level = plyr::revalue(as.character(all_combos$tran), replace = c("FALSE" = "FALSE", "FF" = "TRUE", "85" = "trending"))
all_combos$sig = F
all_combos$sig[which(all_combos$sig_level == "TRUE")] = T
all_combos$trending = F
all_combos$trending[which(all_combos$sig_level == "trending")] = T
all_combos$pcol_tran2 = "white"
all_combos$pcol_tran2[non_na_rows] = paste0(all_combos$pcol[non_na_rows], "85")

all_combos$cluster = factor(all_combos$cluster, levels = rev(convert_all$cluster))
all_combos$hgnc = factor(all_combos$hgnc, levels = c("All", sort(unique(all_combos$hgnc[which(all_combos$hgnc != "All")]))))

ggplot(all_combos, aes(x = hgnc, y = cluster, fill = pcol)) + geom_tile(color = "gray60") + scale_fill_identity() + coord_fixed() + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic"), axis.text.y = element_text(colour = rev(convert_all$color), face=ifelse(rev(convert_all$level) =="secondary","plain","bold"), size=ifelse(rev(convert_all$level) =="secondary", 8, 10))) + xlab("") + ylab("") + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0, 0))

num_hits_by_clust = data.frame(table(all_combos$cluster[which(all_combos$pcol != "white" & all_combos$cluster != "All")]))
num_hits_by_gp    = data.frame(table(all_combos$gene_pop[which(all_combos$pcol != "white" & all_combos$cluster != "All")]))
convert_all$num_hits = num_hits_by_clust$Freq[match(convert_all$cluster, num_hits_by_clust$Var1)]
convert_all_small = convert_all[which(convert_all$num_hits != 0 & convert_all$cluster != "All"),]
all_combos_small = all_combos[which(all_combos$cluster %in% convert_all_small$cluster),]
all_combos_small$cluster = factor(all_combos_small$cluster, levels = rev(convert_all_small$cluster))
all_combos_small$gp_hits = num_hits_by_gp$Freq[match(all_combos_small$gene_pop, num_hits_by_gp$Var1)]
all_combos_small = all_combos_small[which(all_combos_small$gp_hits > 0),]
length(which(! all_combos_small$pcol %in% c("white", "#FDE725", "#2AB07F", "#433E85")))
all_combos_small[which(! all_combos_small$pcol %in% c("white", "#FDE725", "#2AB07F", "#433E85")),]

order_combos = data.frame(table(all_combos_small$gene_pop[which(all_combos_small$pcol != "white")]))
order_combos$border = data.frame(table(all_combos_small$gene_pop[which( ! is.na(all_combos_small$bsig) )]))[,2]
order_combos$gorder = data.frame(table(all_combos_small$gene_pop[which( ! is.na(all_combos_small$gsig) )]))[,2]
order_combos$qorder = data.frame(table(all_combos_small$gene_pop[which( ! is.na(all_combos_small$qsig) )]))[,2]
order_combos = order_combos[order(order_combos$qorder, order_combos$gorder, order_combos$border, decreasing = T),]
all_combos_small$gene_pop = factor(all_combos_small$gene_pop, levels = order_combos$Var1)

all_combos_small$hgnc = tolower(gene_info$human[match(all_combos_small$gene_pop, gene_info$mzebra)])
all_combos_small$hgnc[which(all_combos_small$gene_pop == "LOC101474236")] = "LOC101474236"
all_combos_small$hgnc[which(all_combos_small$gene_pop == "All")] = "All"
all_combos_small$hgnc = factor(all_combos_small$hgnc, levels = unique(all_combos_small$hgnc[match(order_combos$Var1, all_combos_small$gene_pop)]))

pdf("C:/Users/miles/Downloads/neurogen_summary_3_sig.pdf", width = 10, height = 5)
ggplot(all_combos_small, aes(x = hgnc, y = cluster, fill = pcol_tran)) + geom_tile(color = "gray60") + scale_fill_identity() + coord_fixed() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic"), axis.text.y = element_text(colour = rev(convert_all_small$color), face=ifelse(rev(convert_all_small$level) =="secondary","plain","bold"), size=ifelse(rev(convert_all_small$level) =="secondary", 8, 10))) + xlab("") + ylab("") + scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0, 0))
dev.off()

#**********************************************************************
# Temporal Data =======================================================
#**********************************************************************
adj = readRDS("~/scratch/brain/results/adjusted_glmmseq_ffm_15.rds")
# tdf_long = read.csv("C:/Users/miles/Downloads/build_spawn_depth_rolling30_010622.csv")
# tdf_long = read.csv("~/scratch/brain/results/build_spawn_depth_rolling30_010622.csv")
# tdf_wide = reshape(tdf_long, idvar = 'pid', timevar = 'colname', direction = 'wide')
# tdf_wide = read.csv("~/scratch/brain/results/ieg_behavior_time_series_data_011722.csv")
tdf_wide = read.csv("~/scratch/brain/results/ieg_behavior_time_series_data_011822.csv")
tdf_wide = tdf_wide[,which( startsWith(colnames(tdf_wide), "trial_id") | startsWith(colnames(tdf_wide), "build") | startsWith(colnames(tdf_wide), "depth") )]
zdf = read.csv("~/scratch/brain/results/bb15_deg_all_split_by_up_or_down_121621.csv")
# zdf = read.csv("C:/Users/miles/Downloads/bb15_deg_all_split_by_up_or_down_121621.csv")
zdf = zdf[which(zdf$mzebra %in% rownames(bb)),]
zdf$sig_any = zdf$sig_bower_behavior == 1 | zdf$sig_gsi == 1 | zdf$sig_log_spawn_events == 1
zdf = zdf[which(zdf$sig_any),]
bb$cluster = bb$seuratclusters15
bb$trial_id = factor(bb$trial_id)

exp_mean = matrix(0L, nrow = 38, ncol = nrow(zdf))
for (i in 1:nrow(zdf)) {
  if (i %% 100 == 0) { cat(paste0(i, "."))}
  zgene = zdf$mzebra[i]
  zcluster = zdf$cluster[i]
  # exp_df = data.frame(value = bb@assays$RNA@data[zgene, which(bb$cluster == zcluster)], trial_id = bb$trial_id[which(bb$cluster == zcluster)])
  exp_df = data.frame(value = adj[zgene, which(bb$cluster == zcluster)], trial_id = bb$trial_id[which(bb$cluster == zcluster)])
  exp_agr = aggregate(value ~ trial_id, exp_df, mean, drop = F)
  exp_mean[,i] = exp_agr[, 2]
}
write.csv(exp_mean, "C:/Users/miles/Downloads/bb15_deg_data_means_sig.csv")

# big_cor_df = cor_df = data.frame(cat = c(rep("build", 7), rep("depth", 7), rep("spawn", 7)), time = rep(seq(5, 65, by = 10), 3) )
big_cor_df = cor_df = data.frame(cat = c(rep("build", 7), rep("depth", 7), rep("depth_adj", 7)), time = rep(seq(5, 65, by = 10), 3) )
big_cor_df$time = factor(big_cor_df$time)
big_cor_df$value = 0
for( i in 1:nrow(zdf) ) {
  if (i %% 10 == 0) { cat(paste0(i, "."))}
  cat_list = c()
  if (zdf[i, "sig_bower_behavior"] == 1)   { cat_list = c(cat_list, "bdeg15_cors") }
  if (zdf[i, "sig_gsi"] == 1)              { cat_list = c(cat_list, "gdeg15_cors") }
  if (zdf[i, "sig_log_spawn_events"] == 1) { cat_list = c(cat_list, "qdeg15_cors") }
  # tdf_wide$value = exp_mean[match(levels(bb$trial_id), tdf_wide$pid), i]
  tdf_wide$value = exp_mean[match(levels(bb$trial_id), tdf_wide$trial_id), i]
  my_cors = sapply(2:(ncol(tdf_wide)-1), function(x) cor(tdf_wide[,x], tdf_wide[, ncol(tdf_wide)]) )
  big_cor_df[, paste0("res", i)] = my_cors
  big_cor_df$value = big_cor_df[, paste0("res", i)]
  
  # for (this_cat in cat_list) {
  #   png(paste0("~/scratch/brain/results/", this_cat, "/res", i, ".png"), width = 400, height = 300)
  #   print(ggplot(big_cor_df, aes(x = time, y = value, color = cat)) + geom_point(size = 2.5) + theme_classic() + ylab("Pearson r") + xlab("Time") + ggtitle("bDEGs w/ BHVE Up at the 15 Cluster Level - Adjusted"))
  #   dev.off()
  # }
}
# big_cor_df$value = NULL

countChangeDir = function(x) {
  value_dif = unname(cor_mat[,x])[1:6] - unname(cor_mat[,x])[2:7]
  value_dif_sign = sign(value_dif)
  return(sum(value_dif_sign[1:5] != value_dif_sign[2:6]))
}

# big_cor_df = read.csv("~/scratch/brain/results/bb15_deg_cors.csv")
big_cor_df = read.csv("~/scratch/brain/results/bb15_deg_cors_011922.csv")
big_cor_df$X = NULL
big_0 = sapply(colnames(big_cor_df), function(x) length(which(is.na(big_cor_df[,x])))) 
big_cor_df = big_cor_df[,which(big_0 < nrow(big_cor_df))]
big_cor_df[,which(startsWith(colnames(big_cor_df), "res"))] = big_cor_df[,which(startsWith(colnames(big_cor_df), "res"))] ^ 2
# cor_mat = as.matrix(big_cor_df[which(big_cor_df$cat == "build"), which(startsWith(colnames(big_cor_df), "res") & colnames(big_cor_df) %in% paste0("res", which(zdf$sig_bower_behavior == 1)) )])
# rownames(cor_mat) = big_cor_df$time[which(big_cor_df$cat == "build")]
cor_mat = as.matrix(big_cor_df[which(big_cor_df$cat == "depth_adj"), which(startsWith(colnames(big_cor_df), "res") & colnames(big_cor_df) %in% paste0("res", which(zdf$sig_bower_behavior == 1)) )])
rownames(cor_mat) = big_cor_df$time[which(big_cor_df$cat == "depth_adj")]
cor_mat = scale(cor_mat)
hcl = "hi"
my_callback = function(hcl, mat) { print(hcl); hcl <<- hcl; return(hcl) }
pheatmap::pheatmap(cor_mat, cluster_rows = F, clustering_callback = my_callback, filename = "~/scratch/brain/results/depth_adj_bower_deg53_scaled_r2.pdf")
dend = as.dendrogram(hcl)
late_clust = c(dend[[1]][[1]][[1]] %>% labels, dend[[1]][[1]][[2]][[1]] %>% labels)
cor_0    = sapply(1:ncol(cor_mat), function(x) countChangeDir(x))
cor_mat = cor_mat[, which(cor_0 < 2)]
# cor_mat1 = cor_mat[, which(cor_0 < 1)]
cor_maxs = sapply(1:ncol(cor_mat), function(x) which.max(cor_mat[,x]))
# cor_sort = sapply(1:ncol(cor_mat), function(x) which.max(cor_mat[,x]))
clust5  = which(cor_maxs == 1)
clust15 = which(cor_maxs == 2)
clust25 = which(cor_maxs == 3)
clust35 = which(cor_maxs == 4)
clust45 = which(cor_maxs == 5)
clust55 = which(cor_maxs == 6)
clust65 = which(cor_maxs == 7)
late = c(clust5, clust15, clust25)
middle = c(clust25, clust35, clust45)
early = c(clust45, clust55, clust65)
pheatmap::pheatmap(cor_mat, cluster_rows = F, clustering_callback = my_callback, filename = "~/scratch/brain/results/depth_adj_bower_deg53_scaled_r2.pdf")
pheatmap::pheatmap(cor_mat[,clust5], cluster_rows = F, clustering_callback = my_callback, filename = "~/scratch/brain/results/depth_adj_bower_deg53_scaled_r2_5.pdf")
pheatmap::pheatmap(cor_mat[,clust15], cluster_rows = F, clustering_callback = my_callback, filename = "~/scratch/brain/results/depth_adj_bower_deg53_scaled_r2_15.pdf")
pheatmap::pheatmap(cor_mat[,clust25], cluster_rows = F, clustering_callback = my_callback, filename = "~/scratch/brain/results/depth_adj_bower_deg53_scaled_r2_25.pdf")
pheatmap::pheatmap(cor_mat[,clust35], cluster_rows = F, clustering_callback = my_callback, filename = "~/scratch/brain/results/depth_adj_bower_deg53_scaled_r2_35.pdf")
pheatmap::pheatmap(cor_mat[,clust45], cluster_rows = F, clustering_callback = my_callback, filename = "~/scratch/brain/results/depth_adj_bower_deg53_scaled_r2_45.pdf")
pheatmap::pheatmap(cor_mat[,clust55], cluster_rows = F, clustering_callback = my_callback, filename = "~/scratch/brain/results/depth_adj_bower_deg53_scaled_r2_55.pdf")
pheatmap::pheatmap(cor_mat[,clust65], cluster_rows = F, clustering_callback = my_callback, filename = "~/scratch/brain/results/depth_adj_bower_deg53_scaled_r2_65.pdf")
pheatmap::pheatmap(cor_mat[,late]   , cluster_rows = F, clustering_callback = my_callback, filename = "~/scratch/brain/results/depth_adj_bower_deg53_scaled_r2_late.pdf")
pheatmap::pheatmap(cor_mat[,middle] , cluster_rows = F, clustering_callback = my_callback, filename = "~/scratch/brain/results/depth_adj_bower_deg53_scaled_r2_middle.pdf")
pheatmap::pheatmap(cor_mat[,early]  , cluster_rows = F, clustering_callback = my_callback, filename = "~/scratch/brain/results/depth_adj_bower_deg53_scaled_r2_early.pdf")

for (i in list(5, 15, 25, 35, 45, 55, 65, "early", "middle", "late")) {
  print(i)
  if (class(i) == "numeric") {
    i_clean = i + 25
    this_clust = get(paste0("clust", i))
  }
  if (class(i) == "character") {
    i_clean = str_to_title(i)
    this_clust = get(i)
  }
  if (length(this_clust) > 0) {
    png(paste0(paste0("~/scratch/brain/results/depth_adj_bower_deg53_", i, ".png")), width = 600, height = 500)
    this_df = big_cor_df[which(big_cor_df$cat == "depth_adj"),c("time", colnames(cor_mat)[this_clust])]
    if (length(this_clust) > 1) { this_df$mean = rowMeans(this_df[,colnames(cor_mat)[this_clust]]) }
    else                        { this_df$mean = this_df[,colnames(cor_mat)[this_clust]] }
    pdf = melt(this_df, id.var = "time")
    pdf$time = as.numeric(as.vector(pdf$time)) + 25
    pdf$isMean = pdf$variable == "mean"
    # print(ggplot(pdf, aes(x = time, y = value, color = variable)) + geom_point(size = 2.5, alpha = 0.4) + theme_classic() + ylab("R2") + xlab("Time (min to flash freeze)") + ggtitle("bDEGs w/ BHVE Up at the 15 Cluster Level - Adjusted") + geom_smooth(method = "loess", se = F, alpha = 0.4) + scale_x_reverse() + NoLegend())
    print(ggplot(pdf, aes(x = time, y = value)) + geom_point(data = pdf[which(!pdf$isMean),], size = 2.5, alpha = 0.2, aes(color = variable)) + geom_point(data = pdf[which(pdf$isMean),], size = 2.5, color = "black") + geom_smooth(data = pdf[which(pdf$isMean),], method = "loess", se = F, color = "gray40") + theme_classic() + ylab("R2") + xlab("Time (min to flash freeze)") + ggtitle(paste0("bDEG Hits Up at ", i_clean, ". Depth_adj R2 w/ Adjusted")) + scale_x_continuous(breaks = rev(unique(pdf$time)), labels = rev(unique(pdf$time))) + NoLegend())
    dev.off()
  }
}

zdf$peak = NA
zdf$peak[as.numeric(substr(colnames(cor_mat), 4, 10))] = plyr::revalue(as.character(cor_maxs), replace = c("1" = "30", "2" = "40", "3" = '50', '4' = '60', '5' = '70', '6' = '80', '7' = '90'))
zdf$n_change_dir = NA
zdf$n_change_dir[as.numeric(substr(colnames(cor_mat), 4, 10))] = cor_0

early_late_de = diffEnrichTG(zdf$hgnc[which(zdf$peak %in% c("30", "40", "50") & zdf$n_change_dir < 2)], zdf$hgnc[which(zdf$peak %in% c("70", "80", "90") & zdf$n_change_dir < 2)], path_to_gene_info = "~/scratch/m_zebra_ref/gene_info.txt")

library(parallel)
shuffled_res = unlist(mclapply(1:10000, function(x) shuffleCors(x), mc.cores=detectCores()))

shuffleCors = function(x) {
  set.seed(x)
  this_cor_mat = do.call('cbind', lapply(1:ncol(cor_mat), function(col) unname(sample(cor_mat[,col]))))
  this_cor_0    = c()
  for(col in 1:ncol(this_cor_mat)) {
    value_dif = unname(this_cor_mat[,col])[1:6] - unname(this_cor_mat[,col])[2:7]
    value_dif_sign = sign(value_dif)
    this_cor_0 = c(this_cor_0, sum(value_dif_sign[1:5] != value_dif_sign[2:6]))
  }
  return(length(which(this_cor_0 < 2)))
}

#******************************************************************
# IEG FeaturePlot =================================================
#******************************************************************
ieg_cons = c("LOC101487312", "egr1", "npas4")
ieg_like = read.csv("C:/Users/miles/Downloads/ieg_like_fos_egr1_npas4_detected_011521.csv", stringsAsFactors = F)[,1]
ieg_like = c(ieg_cons, ieg_like[which(! ieg_like %in% ieg_cons)])
mat = bb@assays$RNA@counts
mat[which(mat > 1)] = 1
ieg_query_col = "#21918c"
ieg_target_col = "purple"
ieg_both_col = "yellow"
none_col = "gray90"
my.pt.size = 0.4
p_list = list()
for (ieg_target in ieg_like) {
  for (ieg_query in ieg_cons) {
    codf = data.frame(ieg_query = mat[ieg_query,], ieg_target = mat[ieg_target,], UMAP_1 = bb@reductions$umap@cell.embeddings[,"UMAP_1"], UMAP_2 = bb@reductions$umap@cell.embeddings[,"UMAP_2"], col = none_col)
    codf$col[which(codf$ieg_query == 1 & codf$ieg_target == 1)] = ieg_both_col
    codf$col[which(codf$ieg_query == 1 & codf$ieg_target != 1)] = ieg_query_col
    codf$col[which(codf$ieg_query != 1 & codf$ieg_target == 1)] = ieg_target_col
    codf$sum = rowSums(codf[, c("ieg_query", "ieg_target")])
    codf = codf[order(codf$sum, decreasing = F),]
    p = ggplot(codf, aes(UMAP_1, UMAP_2, color = col)) + geom_point(size = my.pt.size) + theme_void() + scale_color_identity()
    p_list[[length(p_list)+1]] = p
  }
}
ppp = 3
pdf("C:/Users/miles/Downloads/ieg_25_by_3_scatter_plot.pdf", width = 3*ppp, height = 25*ppp)
# ppp = 100
# png("C:/Users/miles/Downloads/ieg_25_by_3_scatter_plot.png", width = 3*ppp, height = 25*ppp)
print(plot_grid(plotlist=p_list, ncol = 3))
dev.off()

# value_dif = unname(cor_mat[,3])[1:6] - unname(cor_mat[,3])[2:7]
# value_dif_sign = sign(value_dif)
# sum(value_dif_sign[1:5] == value_dif_sign[2:6])
# testdf = data.frame(time = rep(unique(big_cor_df$time),2), value = c(unname(cor_mat[,3]), predict(sm.spline(unique(big_cor_df$time), unname(cor_mat[,3])), unique(big_cor_df$time), 1)), isSlope = c(rep(F, 7), rep(T, 7)))
# png("~/scratch/brain/results/test3.png", width = 400, height = 300)
# ggplot(testdf, aes(x = time, y = value, color = isSlope)) + geom_point() + geom_smooth(data = testdf[which(!testdf$isSlope),], method = "loess")
# dev.off()

# exp_mean_bower = exp_mean[,which(zdf$sig_bower_behavior == 1)]
# mean_exp_mean_bower = rowMeans(exp_mean_bower)
# exp_mean_bower_pos = exp_mean[,which(zdf$sig_bower_behavior == 1 & zdf$bower_activity_index > 0)]
# mean_exp_mean_bower_pos = rowMeans(exp_mean_bower_pos)
# tdf_wide$value = mean_exp_mean_bower_pos[match(levels(bb$trial_id), tdf_wide$pid)]
# 
# my_cors = sapply(2:(ncol(tdf_wide)-1), function(x) cor(tdf_wide[,x], tdf_wide[, ncol(tdf_wide)]) )
# cor_df = data.frame(value = my_cors, cat = c(rep("build", 7), rep("depth", 7), rep("spawn", 7)), time = rep(seq(5, 65, by = 10), 3) )
# # cor_df$time = factor(cor_df$time)
# ggplot(cor_df, aes(x = time, y = value, color = cat)) + geom_point(size = 2.5) + theme_classic() + ylab("Pearson r") + xlab("Time") + ggtitle("bDEGs w/ BHVE Up at the 15 Cluster Level - Data Slot") + geom_smooth(method = "loess")
# 5-35 > 30
# 90 to 30

bb15_b = read.csv("C:/Users/miles/Downloads/bb15_b_extra_columns.csv")
bb15_d = read.csv("C:/Users/miles/Downloads/bb15_d_extra_columns.csv")
bb53_b = read.csv("C:/Users/miles/Downloads/bb53_b_extra_columns.csv")
bb53_d = read.csv("C:/Users/miles/Downloads/bb53_d_extra_columns.csv")

bb15_b$peak = factor(bb15_b$peak, unique(bb15_b$peak[which(!is.na(bb15_b$peak))]))

peak_df = data.frame(time = rep(unique(bb15_b$peak[which(!is.na(bb15_b$peak))]), 4), cat = c(rep("build_15", 7), rep("depth_adj_15", 7), rep("build_53",7), rep("depth_adj_53",7)), num_peak = c(as.numeric(table(bb15_b$peak[which(bb15_b$n_change_dir < 2)])), as.numeric(table(bb15_d$peak[which(bb15_d$n_change_dir < 2)])), as.numeric(table(bb53_b$peak[which(bb53_b$n_change_dir < 2)])), as.numeric(table(bb53_d$peak[which(bb53_d$n_change_dir < 2)]))))
peak_df$time = as.numeric(as.vector(peak_df$time))
peak_df$time = factor(peak_df$time, levels = sort(unique(peak_df$time)))
ggplot(peak_df, aes(x = time, fill = cat, y = num_peak)) + geom_bar(position = position_dodge2(), stat = 'identity') + ylab("Number of DEGs that Peak at the Time Bin") + xlab("Time Bin")

bb15_b_genes = bb15_b$mzebra[which(bb15_b$n_change_dir < 2)]
bb15_b_genes = bb15_b_genes[which(duplicated(bb15_b_genes))]
bb15_b_dup = bb15_b[which(bb15_b$mzebra %in% bb15_b_genes & bb15_b$n_change_dir < 2), c("mzebra", "peak")]
bb15_b_dup[order(bb15_b_dup$mzebra),]

bb15_d_genes = bb15_d$mzebra[which(bb15_d$n_change_dir < 2)]
bb15_d_genes = bb15_d_genes[which(duplicated(bb15_d_genes))]
bb15_d_dup = bb15_d[which(bb15_d$mzebra %in% bb15_d_genes & bb15_d$n_change_dir < 2), c("mzebra", "peak")]
bb15_d_dup[order(bb15_d_dup$mzebra),]

bb53_b_genes = bb53_b$mzebra[which(bb53_b$n_change_dir < 2)]
bb53_b_genes = bb53_b_genes[which(duplicated(bb53_b_genes))]
bb53_b_dup = bb53_b[which(bb53_b$mzebra %in% bb53_b_genes & bb53_b$n_change_dir < 2), c("mzebra", "peak", "X", "n_change_dir")]
bb53_b_dup[order(bb53_b_dup$mzebra),]

# IEG Temporal =====================================================
ieg_like = read.csv("C:/Users/miles/Downloads/ieg_like_fos_egr1_npas4_detected_011521.csv", stringsAsFactors = F)[,1]
ieg_sum = read.csv("C:/Users/miles/Downloads/ieg_summary_data_by_cat_cluster_goi_010321.csv")
ieg_sum = ieg_sum[which(ieg_sum$is_sig),]
ieg_sum$hgnc = ieg_sum$gene
ieg_sum$gene_pop = ieg_sum$mzebra
ieg_sum$cluster[which(ieg_sum$cluster == FALSE)] = "All"
ieg_sum$gene_pop[which(ieg_sum$gene_pop == FALSE)] = "All"
ieg_sum$level_old = paste0(ieg_sum$level, "_", ieg_sum$cluster)
ieg_sum$level_old_gp = paste0(ieg_sum$level_old, "_", ieg_sum$gene_pop)
ieg_sum$cat_level_old_gp = paste0(ieg_sum$cat, "_", ieg_sum$level_old_gp)
ieg_sum = ieg_sum[which(ieg_sum$cat != "gsi"),]
bb$trial_id = factor(bb$trial_id)

tdf = read.csv("C:/Users/miles/Downloads/ieg_behavior_time_series_data_011822.csv")
tdf = tdf[,which( startsWith(colnames(tdf), "trial_id") | startsWith(colnames(tdf), "build") | startsWith(colnames(tdf), "depth") | startsWith(colnames(tdf), "spawn") )]
tdf = tdf[,which(! endsWith(colnames(tdf), ".1") ),]
tdf = tdf[,which(! (startsWith(colnames(tdf), "depth") & ! endsWith(colnames(tdf), "adj")) )]
tdf$trial_id = factor(tdf$trial_id, levels = levels(bb$trial_id))

num_cell_df = data.frame()
big_mean_df = data.frame(matrix(0L, nrow = 38, ncol = nrow(ieg_sum)*length(ieg_like)), row.names = levels(bb$trial_id))
big_cor_df = data.frame(matrix(0L, nrow = 7, ncol = nrow(ieg_sum)*length(ieg_like)), row.names = c(30, 40, 50, 60, 70, 80, 90))
colnames(big_mean_df) = colnames(big_cor_df) = sapply(ieg_sum$cat_level_old_gp, function(x) paste0(ieg_like, "_", x))
for (i in 1:nrow(ieg_sum)) {
  cat(paste0(i, "."))
  for (ieg_like_gene in ieg_like) {
    my.cat = ieg_sum$cat[i]
    my.gene = ieg_sum$mzebra[i]
    my.level = ieg_sum$level[i]
    my.cluster = ieg_sum$cluster[i]
    if (my.level == "primary")   { bb$cluster = bb$seuratclusters15 }
    if (my.level == "secondary") { bb$cluster = bb$seuratclusters53 }
    if (my.level == "all")       { bb$cluster = "All"               }
    if (my.cat == "bower") { this_tdf = tdf[,which( startsWith(colnames(tdf), "depth") )] }
    if (my.cat == "quiver") { this_tdf = tdf[,which( startsWith(colnames(tdf), "spawn") )] }
    if (my.gene == F) { my.gene.cells = rep(T, ncol(bb)) } else { my.gene.cells = bb@assays$RNA@counts[my.gene,] > 0 }
    df = data.frame(value = bb@assays$RNA@data[ieg_like_gene, which(my.gene.cells & bb$cluster == my.cluster)], trial_id = bb$trial_id[which(my.gene.cells & bb$cluster == my.cluster)])
    
    my_x = bb@assays$RNA@counts[ieg_like_gene, ] > 0 & my.gene.cells & bb$cluster == my.cluster
    num_cell_df = rbind( num_cell_df, data.frame(ieg_pop = ieg_sum$cat_level_old_gp[i], ieg_like = ieg_like_gene, num_cell = length(which(my_x)), num_pair = length(unique(bb$pair[which(my_x)])), num_ind = length(unique(bb$subsample[which(my_x)])), num_sample = length(unique(bb$sample[which(my_x)])), num_sample_paired = length(unique(substr(bb$sample[which(my_x)], 2, 2))), inds = paste0(unique(bb$subsample[which(my_x)]), collapse = ", ") ))
    if (num_cell_df$num_pair[nrow(num_cell_df)] == 19) {
      mean_df = aggregate(value ~ trial_id, df, mean, drop = F)
      big_mean_df[,paste0(ieg_like_gene, "_", ieg_sum$cat_level_old_gp[i])] = mean_df[,2]
      big_cor_df[,paste0(ieg_like_gene, "_", ieg_sum$cat_level_old_gp[i])] = sapply(1:ncol(this_tdf), function(x) cor(mean_df[,2], this_tdf[,x]) ^ 2 ) 
    }
    
  }
}
print("Done")

cor_wide = as.data.frame(pivot_longer(big_cor_df, cols = colnames(big_cor_df)))
cor_wide$time = unlist(lapply(rownames(big_cor_df), function(x) rep(x, ncol(big_cor_df))))
cor_wide$time = as.numeric(cor_wide$time)
cor_wide[, c("ieg_like", "cat_level_old_gp")] = reshape2::colsplit(cor_wide$name, "_", c("1", "2"))
for (ieg_like_gene in ieg_like) {
  pdf(paste0("C:/Users/miles/Downloads/brain/results/bb/ieg_time/", ieg_like_gene, ".pdf") , width = 5, height = 5)
  print(ggplot(cor_wide[which(cor_wide$ieg_like == ieg_like_gene),], aes(x = time, y = value)) + geom_point(alpha = 0.4) + geom_smooth(method = "loess", se = T) + NoLegend() + ylab("R2") + xlab("Time (min to flash freeze)") + ggtitle(paste0("All IEG Results w/ IEG: ", ieg_like_gene)))
  dev.off()
}

# cor_mean = aggregate(value ~ time + ieg_like, cor_wide, mean)
# pdf(paste0("C:/Users/miles/Downloads/ieg_mean_cor.pdf") , width = 5, height = 5)
# print(ggplot(cor_mean, aes(x = time, y = value, color = ieg_like)) + geom_point() + geom_smooth(method = "loess", se = F) + NoLegend() + ylab("R2") + xlab("Time (min to flash freeze)") + scale_x_continuous(breaks = rev(unique(cor_wide$time)), labels = rev(unique(cor_wide$time))) + ggtitle("All IEG Means"))
# dev.off()
pdf(paste0("C:/Users/miles/Downloads/ieg_mean_cor_conf_int.pdf") , width = 5, height = 5)
print(ggplot(cor_wide, aes(x = time, y = value, color = ieg_like)) + geom_smooth(method = "loess", se = T) + NoLegend() + ylab("R2") + xlab("Time (min to flash freeze)") + scale_x_continuous(breaks = rev(unique(cor_wide$time)), labels = rev(unique(cor_wide$time))) + ggtitle("All IEG Means"))
dev.off()

ieg_mean_bower = read.csv("C:/Users/miles/Downloads/ieg_summary_subsample_means_bower.csv")
ieg_mean = ieg_mean_bower
ieg_mean[,1] = NULL
tdf$subsample = bb$subsample[match(tdf$trial_id, bb$trial_id)]
tdf = tdf[order(tdf$subsample),]
my_r2 = data.frame()
for (i in 1:ncol(ieg_mean)) {
  this_col = colnames(tdf)[which(! colnames(tdf) %in% c("trial_id", "subsample"))]
  this_r2 = sapply(this_col, function(x) cor(ieg_mean[,i], tdf[,x]) ^ 2 ) 
  my_r2 = rbind(my_r2, data.frame(r2 = this_r2, cat = reshape2::colsplit(this_col, "_", c("1", "2"))[,1], time = rep(c(30, 40, 50, 60, 70, 80, 90), 3) ))
}

my_cor_mean = aggregate(r2 ~ cat + time, my_r2, mean)
ggplot(my_cor_mean, aes(x = time, y = r2, color = cat)) + geom_point() + geom_smooth(method = "loess", se = F)

# 9 Plot DEGs =====================================================
bb15 = read.csv("C:/Users/miles/Downloads/bb15_deg_all_split_by_up_or_down_121621.csv")
bb53 = read.csv("C:/Users/miles/Downloads/bb53_deg_all_split_by_up_or_down_121621.csv")
bb15$sig_any = bb15$sig_bower_behavior == 1 | bb15$sig_gsi == 1 | bb15$sig_log_spawn_events == 1
bb53$sig_any = bb53$sig_bower_behavior == 1 | bb53$sig_gsi == 1 | bb53$sig_log_spawn_events == 1
bb15 = bb15[which(bb15$sig_any & bb15$mzebra %in% rownames(bb)),]
bb53 = bb53[which(bb53$sig_any & bb53$mzebra %in% rownames(bb)),]
# all_bdeg = unique(c(bb15$mzebra[which(bb15$sig_bower_behavior == 1)], bb53$mzebra[which(bb53$sig_bower_behavior == 1)]))
# all_gdeg = unique(c(bb15$mzebra[which(bb15$sig_gsi == 1)], bb53$mzebra[which(bb53$sig_gsi == 1)]))
# all_qdeg = unique(c(bb15$mzebra[which(bb15$sig_log_spawn_events == 1)], bb53$mzebra[which(bb53$sig_log_spawn_events == 1)]))
# all_deg = data.frame(table(c(all_bdeg, all_gdeg, all_qdeg)))
# all_deg = as.vector(all_deg$Var1[which(all_deg$Freq == 3)])
all_deg = read.csv("C:/Users/miles/Downloads/overlap_of_bDEGs_qDEGs_gDEGs_012021.csv")[,5]
my.pt.size = 0.6
p_list = list()
for (cat in c("bower_behavior", "gsi", "log_spawn_events")) {
  print(cat)
  if (cat == "bower_behavior") { cat_mod = "bower_activity_index" } else { cat_mod = cat }
  this_bb15 = bb15[which(bb15[, paste0("sig_", cat)] == 1),]
  this_bb53 = bb53[which(bb53[, paste0("sig_", cat)] == 1),]
  # this_bb15_up = bb15[which(bb15[, cat_mod] > 0)]
  # this_bb15_down = bb15[which(bb15[, cat_mod] < 0)]
  score_df = data.frame(cell = colnames(bb), UMAP_1 = bb@reductions$umap@cell.embeddings[,"UMAP_1"], UMAP_2 = bb@reductions$umap@cell.embeddings[,"UMAP_2"], up_score = 0, down_score = 0, all_score = 0, ovlp_score = 0)
  for (i in 1:nrow(this_bb15)) {
    this_gene = this_bb15$mzebra[i]
    this_cluster = this_bb15$cluster[i]
    isUp = this_bb15[i, cat_mod] > 0
    this_score = as.numeric(bb@assays$RNA@counts[this_gene,] > 0)
    score_df$all_score[which(bb$seuratclusters15 == this_cluster)] = score_df$all_score[which(bb$seuratclusters15 == this_cluster)] + this_score[which(bb$seuratclusters15 == this_cluster)]
    if (isUp) { score_df$up_score[which(bb$seuratclusters15 == this_cluster)] = score_df$up_score[which(bb$seuratclusters15 == this_cluster)] + this_score[which(bb$seuratclusters15 == this_cluster)] } else { score_df$down_score[which(bb$seuratclusters15 == this_cluster)] = score_df$down_score[which(bb$seuratclusters15 == this_cluster)] + this_score[which(bb$seuratclusters15 == this_cluster)] }
    if (this_gene %in% all_deg) { score_df$ovlp_score[which(bb$seuratclusters15 == this_cluster)] = score_df$ovlp_score[which(bb$seuratclusters15 == this_cluster)] + this_score[which(bb$seuratclusters15 == this_cluster)] }
  }
  for (i in 1:nrow(this_bb53)) {
    this_gene = this_bb53$mzebra[i]
    this_cluster = this_bb53$cluster[i]
    isUp = this_bb53[i, cat_mod] > 0
    this_score = as.numeric(bb@assays$RNA@counts[this_gene,] > 0)
    score_df$all_score[which(bb$seuratclusters53 == this_cluster)] = score_df$all_score[which(bb$seuratclusters53 == this_cluster)] + this_score[which(bb$seuratclusters53 == this_cluster)]
    if (isUp) { score_df$up_score[which(bb$seuratclusters53 == this_cluster)] = score_df$up_score[which(bb$seuratclusters53 == this_cluster)] + this_score[which(bb$seuratclusters53 == this_cluster)] } else { score_df$down_score[which(bb$seuratclusters53 == this_cluster)] = score_df$down_score[which(bb$seuratclusters53 == this_cluster)] + this_score[which(bb$seuratclusters53 == this_cluster)] }
    if (this_gene %in% all_deg) { score_df$ovlp_score[which(bb$seuratclusters53 == this_cluster)] = score_df$ovlp_score[which(bb$seuratclusters53 == this_cluster)] + this_score[which(bb$seuratclusters53 == this_cluster)] }
  }
  score_df = score_df[order(score_df$up_score, decreasing = F),]
  # score_df$up_score[which(score_df$up_score == 0)] = NA
  p_list[[paste0(cat, "_", "up")]] = ggplot(score_df, aes(x = UMAP_1, y = UMAP_2, color = up_score)) + geom_point(size = my.pt.size) + scale_color_gradientn(colors = viridis(100), na.value = "grey80") + theme_void() + NoLegend()
  score_df = score_df[order(score_df$down_score, decreasing = F),]
  # score_df$down_score[which(score_df$down_score == 0)] = NA
  p_list[[paste0(cat, "_", "down")]] = ggplot(score_df, aes(x = UMAP_1, y = UMAP_2, color = all_score)) + geom_point(size = my.pt.size) + scale_color_gradientn(colors = viridis(100), na.value = "grey80") + theme_void() + NoLegend()
  score_df = score_df[order(score_df$all_score, decreasing = F),]
  # score_df$all_score[which(score_df$all_score == 0)] = NA
  p_list[[paste0(cat, "_", "all")]] = ggplot(score_df, aes(x = UMAP_1, y = UMAP_2, color = all_score)) + geom_point(size = my.pt.size) + scale_color_gradientn(colors = viridis(100), na.value = "grey80") + theme_void() + NoLegend()
  score_df = score_df[order(score_df$ovlp_score, decreasing = F),]
  # score_df$ovlp_score[which(score_df$ovlp_score == 0)] = NA
  p_list[[paste0(cat, "_", "ovlp")]] = ggplot(score_df, aes(x = UMAP_1, y = UMAP_2, color = ovlp_score)) + geom_point(size = my.pt.size) + scale_color_gradientn(colors = viridis(100), na.value = "grey80") + theme_void() + NoLegend()
}

pdf("C:/Users/miles/Downloads/zack_9_deg_plot.pdf", width = 6*4, height = 6*3)
p = plot_grid(plotlist=p_list, ncol = 4)
print(p)
dev.off()

png("C:/Users/miles/Downloads/zack_9_deg_plot.png", width = 1300, height = 1000)
p = plot_grid(plotlist=p_list, ncol = 4)
print(p)
dev.off()

# print(ggplot(pdf, aes(x = time, y = value)) + geom_point(data = pdf[which(!pdf$isMean),], size = 2.5, alpha = 0.2, aes(color = variable)) + geom_point(data = pdf[which(pdf$isMean),], size = 2.5, color = "black") + geom_smooth(data = pdf[which(pdf$isMean),], method = "loess", se = F, color = "gray40") + theme_classic() + ylab("R2") + xlab("Time (min to flash freeze)") + ggtitle(paste0("bDEG Hits Up at ", i_clean, ". Depth_adj R2 w/ Adjusted")) + scale_x_continuous(breaks = rev(unique(pdf$time)), labels = rev(unique(pdf$time))) + NoLegend())

sub_meta = aggregate(depth_5_35 + depth_15_45 + depth_25_55 + depth_35_65 + depth_45_75 + depth_55_85 + depth_65_95 + build_5_35 + build_15_45 + build_25_55 + build_35_65 + build_45_75 + build_55_85 + build_65_95 ~ subsample, bb@meta.data, mean)
mean_df = read.csv("~/Downloads/ieg_summary_subsample_means_bower.csv")

perm_sub = lapply(unique(bb$subsample), function(x) sample(1:38, length(which(bb$subsample == x)), replace = T))
names(perm_sub) = unique(bb$subsample)
bb$perm = sapply(bb$subsample, function(x) { new = perm_sub[[x]][length(perm_sub[[x]])]; perm_sub[[x]] <<- perm_sub[[x]][-length(perm_sub[[x]])]; return(new); })

perm_sub = lapply(unique(bb$subsample), function(x) { num_b = sample(9:10, 1); sample(1:38, length(which(bb$subsample == x)), replace = T); } )

bb$perm1 = as.numeric(plyr::revalue(bb$subsample, replace = setNames(sample(1:38, 38), unique(bb$subsample))))
bb$perm2 = permSubsamples(1) # the number you enter in the function is the seed for randomness

df = unique(bb@meta.data[, c("subsample", "perm1", "perm2")])
rownames(df) = 1:38
length(which( startsWith(df$subsample, "b") & df$perm1_b ))
length(which( startsWith(df$subsample, "b") & df$perm2_b ))

permSubsamples = function(x) {
  set.seed(x)
  isEven = x %% 2 == 0
  if (isEven) { num1 = 10; } else { num1 = 9; }
  num2 = 19 - num1
  b_subs = 1:19
  c_subs = 20:38
  real_subs = as.vector(unique(bb$subsample))
  b_subs_names = real_subs[1:19]
  c_subs_names = real_subs[20:38]
  new_b = c(sample(b_subs, num1), sample(c_subs, num2))
  new_c = c(b_subs[which(! b_subs %in% new_b)], c_subs[which(! c_subs %in% new_b)])
  my_replace = c(new_b, new_c)
  names(my_replace) = c(sample(b_subs_names), sample(c_subs_names))
  return(as.numeric(plyr::revalue(bb$subsample, my_replace)))
}

#*******************************************************************************
# Supplement ===================================================================
#*******************************************************************************
subsample_df$pair.new = subsample_df$pair
subsample_df$pair.new[1:19] = 1:19
bsubdf = subsample_df[which( startsWith(subsample_df$subsample, "b") ),]
csubdf = subsample_df[which( startsWith(subsample_df$subsample, "c") ),]
subsample_df$pair.new[20:38] = bsubdf$pair.new[match(csubdf$pair, bsubdf$pair)]
bb$pair.new = subsample_df$pair.new[match(bb$subsample, subsample_df$subsample)]
bb$subsample.new = paste0( tolower(substr(bb$cond, 1, 1)), bb$pair.new  )
bb$subsample.new = factor(bb$subsample.new, levels = unique(bb$subsample.new[order(bb$cond, rev(bb$sample), decreasing = T)]))
my_pal = rainbow(19)
my_pal = darken(my_pal, 0.15)
show_col(my_pal)

pdf("~/research/brain/results/supplement/ncountrna.pdf", width = 10, height = 4)
# ggplot(bb@meta.data, aes(x = subsample, y = nCount_RNA, color = subsample)) + geom_boxplot(outlier.shape = NA) + geom_point(alpha = 0.05, position = position_jitterdodge()) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + scale_y_continuous(expand = c(0,0)) + NoLegend() + xlab("")
ggplot(bb@meta.data, aes(x = subsample.new, y = nCount_RNA, color = pair.new)) + geom_violin() + geom_boxplot(outlier.shape = NA, width=0.25) + geom_point(alpha = 0.025, position = position_jitterdodge()) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + scale_y_continuous(expand = c(0,0)) + NoLegend() + xlab("")
dev.off()

pdf("~/research/brain/results/supplement/nfeaturerna.pdf", width = 10, height = 4)
# ggplot(bb@meta.data, aes(x = subsample, y = nCount_RNA, color = subsample)) + geom_boxplot(outlier.shape = NA) + geom_point(alpha = 0.05, position = position_jitterdodge()) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + scale_y_continuous(expand = c(0,0)) + NoLegend() + xlab("")
ggplot(bb@meta.data, aes(x = subsample.new, y = nFeature_RNA, color = pair.new)) + geom_violin() + geom_boxplot(outlier.shape = NA, width=0.25) + geom_point(alpha = 0.025, position = position_jitterdodge()) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + scale_y_continuous(expand = c(0,0)) + NoLegend() + xlab("")
dev.off()

pdf("~/research/brain/results/supplement/pctmt.pdf", width = 10, height = 4)
# ggplot(bb@meta.data, aes(x = subsample, y = nCount_RNA, color = subsample)) + geom_boxplot(outlier.shape = NA) + geom_point(alpha = 0.05, position = position_jitterdodge()) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + scale_y_continuous(expand = c(0,0)) + NoLegend() + xlab("")
ggplot(bb@meta.data, aes(x = subsample.new, y = pct_mt, color = pair.new)) + geom_violin() + geom_boxplot(outlier.shape = NA, width=0.25) + geom_point(alpha = 0.025, position = position_jitterdodge()) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + scale_y_continuous(expand = c(0,0)) + NoLegend() + xlab("")
dev.off()

# Reads (see line number 200)

pdf("~/research/brain/results/supplement/reads_subsample.pdf", width = 10, height = 4)
# ggplot(bb@meta.data, aes(x = subsample, y = nCount_RNA, color = subsample)) + geom_boxplot(outlier.shape = NA) + geom_point(alpha = 0.05, position = position_jitterdodge()) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + scale_y_continuous(expand = c(0,0)) + NoLegend() + xlab("")
ggplot(bb@meta.data, aes(x = subsample.new, y = reads, color = pair.new)) + geom_violin() + geom_boxplot(outlier.shape = NA, width=0.25) + geom_point(alpha = 0.01, position = position_jitterdodge()) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + scale_y_continuous(expand = c(0,0)) + NoLegend() + xlab("") + ylab("Number of Reads")
dev.off()

# PCRC
pcrc = read.csv("~/research/brain/data/pc_20_rc_20_10kb_bins_25kb_genes_on_lg_11_peak_by_bin.csv")[,2]
convert15$new.full = str_replace(convert15$new.full, "Astro", "RGC")
convert53$new = str_replace(convert53$new, "Astro", "RGC")
Idents(bb) = factor(convert53$new[match(bb$seuratclusters53, convert53$old)], levels = convert53$new)
real_res = markerExpPerCellPerClusterQuick(bb, pcrc)
pdf("~/research/brain/results/pcrc2020_53_020922.pdf", width = 12, height = 5)
real_res[[1]] + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + ggtitle("") + ylab("Mean Normalized Expression of Castle-Divergent Genes") + xlab("") + ylab("")
dev.off()

Idents(bb) = factor(convert15$new.full[match(bb$seuratclusters15, convert15$old)], levels = convert15$new.full)
real_res = markerExpPerCellPerClusterQuick(bb, pcrc, pt.alpha = 0.02)
pdf("~/research/brain/results/pcrc2020_15_020922.pdf", width = 6, height = 5)
real_res[[1]] + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) + ggtitle("") + ylab("Mean Normalized Expression of Castle-Divergent Genes") + xlab("") + ylab("")
dev.off()

mat = bb@assays$RNA@counts[pcrc,]
mat[which(mat>1)]=1
bb$pcrc = colSums(mat)
bb$pcrc.norm = bb$pcrc / bb$nFeature_RNA
pdf("~/research/brain/results/pcrc2020_umap.pdf", width = 6, height = 6)
print(FeaturePlot(bb, "pcrc", order = T, pt.size = 0.9)+ NoLegend() + ggtitle("") + scale_color_viridis_c()) 
dev.off()

pdf("~/research/brain/results/pcrc2020_norm_umap.pdf", width = 6, height = 6)
print(FeaturePlot(bb, "pcrc.norm", order = T, pt.size = 0.9)+ NoLegend() + ggtitle("") + scale_color_viridis_c()) 
dev.off()

sum15 = data.frame(readxl::read_excel("~/Downloads/pcrc2020_enrichment.xlsx", sheet = "15 Summary"))
sum53 = data.frame(readxl::read_excel("~/Downloads/pcrc2020_enrichment.xlsx", sheet = "53 Summary"))
gp = data.frame(readxl::read_excel("~/Downloads/pcrc2020_enrichment.xlsx", sheet = "Gene Pop Summary"))
gp15 = data.frame(readxl::read_excel("~/Downloads/pcrc2020_enrichment.xlsx", sheet = "Gene Pop By 15 Cluster Summary"))
gp53 = data.frame(readxl::read_excel("~/Downloads/pcrc2020_enrichment.xlsx", sheet = "Gene Pop By 53 Cluster Summary"))

sum15$level = "15"
sum53$level = "53"
gp15$level = "gp15"
gp53$level = "gp53"

sum15$new = convert15$new.full[match(sum15$cluster, convert15$old)]
sum53$new = convert53$new[match(sum53$cluster, convert53$old)]

sum.cluster = rbind(sum15[, c("cluster", "new", "d", "mag_pos", "p_perm")], sum53[, c("cluster", "new", "d", "mag_pos", "p_perm")])
sum.cluster = sum.cluster[which(!duplicated(sum.cluster$new)),]
sum.cluster$new = factor(sum.cluster$new, levels = convert_all$cluster)
sum.cluster$neg_log_perm = -log10(sum.cluster$p_perm)
sum.cluster$sig = sum.cluster$p_perm < 0.05
sum.cluster$mag_pos = factor(sum.cluster$mag_pos, levels = c("negligible", "small", "medium", "large"))
pdf("~/research/brain/results/pcrc_perm_and_d.pdf", width = 12, height = 5)
ggplot(sum.cluster, aes(x = new, y = neg_log_perm, fill = mag_pos, color = sig)) + geom_bar(stat = "identity") + geom_hline(yintercept = -log10(0.05), color = "gray60", linetype = "dashed") + xlab("") + ylab(expression("-"*Log["10"]*" P")) + scale_y_continuous(expand = c(0,0)) + scale_fill_viridis(discrete = T,drop=TRUE, limits = levels(sum.cluster$mag_pos), name = "Effect Size") + scale_color_manual(values = c("white", "black"), guide="none") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + labs(fill="Cohen's d") + ylab("")
dev.off()

zgp = read.csv("~/Downloads/goi_less_10k_122821_by_cat_hgnc.csv")
# zgp$gp.col = plyr::revalue(zgp$cat, replace = c("receptor" = "#0CECDD", "tf" = "goldenrod1", "ligand" = "#FF67E7", "other" = "#C400FF"))
# zgp$gp.col = plyr::revalue(zgp$cat, replace = c("receptor" = "#542E71", "tf" = "#FB3640", "ligand" = "#FDCA40", "other" = "#A799B7"))
zgp$gp.col = plyr::revalue(zgp$cat, replace = c("receptor" = "goldenrod2", "tf" = "#FE04FF", "ligand" = "#002DD1", "other" = "#00E7EC"))

convert_all2 = rbind(convert_all[2:nrow(convert_all),], convert_all[1,])

gp$perm_p = gp$p
gp$cluster = "All"
gp$cluster_name = "All"
gp$level = "all"
gp$hgnc = gp$human
# gp15$perm_p = gp15$p
gp15$level = "primary"
gp53$level = "secondary"
gp15$gene = gp15$mzebra
gp53$gene = gp53$mzebra
gp15$cluster_name = gp15$cell_type
gp53$cluster_name = gp53$cell_type
cols.to.keep = c("level", "cluster", "cluster_name", "gene", "hgnc", "perm_p", "neg")
gpsig = rbind(gp[which(gp$perm_p < 0.05), cols.to.keep], gp15[which(gp15$perm_p < 0.05), cols.to.keep], gp53[which(gp53$perm_p < 0.05), cols.to.keep])
gpsig$level_old = paste0(gpsig$level, "_", gpsig$cluster)
gpsig$level_old_gp = paste0(gpsig$level_old, "_", gpsig$gene)
convert_all2$id = 1:nrow(convert_all2)
gpsig$id = convert_all2$id[match(gpsig$level_old, convert_all2$level_old)]
gpsig = gpsig[which(gpsig$gene %in% zgp$mzebra),]

all_combos = expand.grid(convert_all2$level_old, unique(gpsig$gene))
colnames(all_combos) = c("level_old", "gene")
all_combos$level_old_gp = paste0(all_combos$level_old, "_", all_combos$gene)
all_combos$isSig = all_combos$level_old_gp %in% gpsig$level_old_gp
all_combos$neg = gpsig$neg[match(all_combos$level_old_gp, gpsig$level_old_gp)]
all_combos$pos = 10000 - all_combos$neg
all_combos$pct.pos = (all_combos$pos / 10000) * 100
all_combos$pct.neg = (all_combos$neg / 10000) * 100
all_combos$perm_p = gpsig$perm_p[match(all_combos$level_old_gp, gpsig$level_old_gp)]
all_combos$neg_log_p = -log10(all_combos$perm_p)
all_combos$col = "white"
all_combos$col[which(all_combos$isSig)] = convert_all2$color[match(all_combos$level_old[which(all_combos$isSig)], convert_all2$level_old)]
all_combos$cluster = convert_all2$cluster[match(all_combos$level_old, convert_all2$level_old)]
all_combos$cluster = factor(all_combos$cluster, levels = convert_all2$cluster[which(convert_all2$cluster %in% all_combos$cluster[which(all_combos$isSig)])])
all_combos = all_combos[which( !is.na(all_combos$cluster) ),]
all_combos$id = convert_all2$id[match(all_combos$level_old, convert_all2$level_old)]
all_combos$gene = factor(all_combos$gene, levels = unique(gpsig$gene[order(gpsig$id, decreasing = F)]))
all_combos$hgnc = gene_info$human[match(all_combos$gene, gene_info$mzebra)]
all_combos$gene_hgnc = paste0(all_combos$gene, " (", tolower(all_combos$hgnc), ")")
all_combos$label = as.vector(all_combos$gene)
all_combos$label[which(startsWith(all_combos$label, "LOC"))] = all_combos$gene_hgnc[which(startsWith(all_combos$label, "LOC"))]
all_combos$label = factor(all_combos$label, levels = unique(all_combos$label[order(all_combos$gene)]))

# ggplot(all_combos, aes(x = label, y = cluster, color = col, size = neg)) + geom_point() + theme_classic() + scale_color_identity() + scale_size_continuous(limits = c(9500, 10000), name = "# of Permutations Less Than Real") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = zgp$gp.col[match(levels(all_combos$gene), zgp$mzebra)]), axis.text.y = element_text(color = convert_all2$color[match(levels(all_combos$cluster), convert_all2$cluster)], face = ifelse(convert_all2$level[match(levels(all_combos$cluster), convert_all2$cluster)] == "secondary", "plain", "bold"), size = ifelse(convert_all2$level[match(levels(all_combos$cluster), convert_all2$cluster)] == "secondary", 8, 10))) + xlab("") + ylab("")

# ggplot(all_combos, aes(x = label, y = cluster, color = col, size = neg_log_p)) + geom_point() + theme_classic() + scale_color_identity() + scale_size_continuous(name = expression("-"*Log["10"]*"(p)")) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = zgp$gp.col[match(levels(all_combos$gene), zgp$mzebra)]), axis.text.y = element_text(color = convert_all2$color[match(levels(all_combos$cluster), convert_all2$cluster)], face = ifelse(convert_all2$level[match(levels(all_combos$cluster), convert_all2$cluster)] == "secondary", "plain", "bold"), size = ifelse(convert_all2$level[match(levels(all_combos$cluster), convert_all2$cluster)] == "secondary", 8, 10))) + xlab("") + ylab("")

pdf("~/research/brain/results/supplement/pcrc_gp.pdf", width = 20, height = 6)
# Circle
ggplot(all_combos, aes(x = label, y = cluster, color = col, size = pct.neg)) + geom_point() + theme_classic() + scale_color_identity() + scale_size_continuous(name = "% of Permutations Less Than Real") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = zgp$gp.col[match(levels(all_combos$gene), zgp$mzebra)]), axis.text.y = element_text(color = convert_all2$color[match(levels(all_combos$cluster), convert_all2$cluster)], face = ifelse(convert_all2$level[match(levels(all_combos$cluster), convert_all2$cluster)] == "secondary", "plain", "bold"), size = ifelse(convert_all2$level[match(levels(all_combos$cluster), convert_all2$cluster)] == "secondary", 8, 10))) + xlab("") + ylab("")
# Square
# ggplot(all_combos, aes(x = label, y = cluster, fill = col)) + geom_tile() + scale_fill_identity() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(color = convert_all$color[match(levels(all_combos$cluster), convert_all$cluster)], face = ifelse(convert_all$level[match(levels(all_combos$cluster), convert_all$cluster)] == "secondary", "plain", "bold"), size = ifelse(convert_all$level[match(levels(all_combos$cluster), convert_all$cluster)] == "secondary", 8, 10))) + xlab("") + ylab("")
dev.off()

# IEG
pdf("~/research/brain/results/ieg_score_umap.pdf", width = 6, height = 6)
print(FeaturePlot(bb, "ieg_score", order = T)+ NoLegend() + ggtitle("") + scale_color_viridis_c()) 
dev.off()

xmin = 0.02
xmax = 0.04
ieg_df15 = data.frame(score = bb$ieg_score, cond = bb$cond, old = bb$seuratclusters15, new = convert15$new.full[match(bb$seuratclusters15, convert15$old)], level = "15")
ieg_df53 = data.frame(score = bb$ieg_score, cond = bb$cond, old = bb$seuratclusters53, new = convert53$new[match(bb$seuratclusters53, convert53$old)], level = "53")

ieg_df15$new = factor(ieg_df15$new, levels = convert15$new.full)
ieg_df15$cond = factor(ieg_df15$cond, levels = c("CTRL", "BHVE"))
pdf("~/research/brain/results/ieg_by_cluster15_cond.pdf", width = 8, height = 5)
ggplot(ieg_df15, aes(x = new, y = score, color = cond, fill = cond)) + geom_split_violin(alpha = 0.3, adjust = 8) + theme_classic() + stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.75, position = position_dodge(width = .75)) + scale_fill_manual(values = c("#21918C", "goldenrod1")) + scale_color_manual(values = c("#21918C", "goldenrod1")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = convert15$col, size = 15), axis.text.y = element_text(size = 15)) + xlab("") + ylab("IEG Score")
dev.off()

ieg_df53$new = factor(ieg_df53$new, levels = convert53$new)
ieg_df53$cond = factor(ieg_df53$cond, levels = c("CTRL", "BHVE"))
pdf("~/research/brain/results/ieg_by_cluster53_cond.pdf", width = 14, height = 5)
ggplot(ieg_df53, aes(x = new, y = score, color = cond, fill = cond)) + geom_split_violin(alpha = 0.3, adjust = 8) + theme_classic() + stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.75, position = position_dodge(width = .75)) + scale_fill_manual(values = c("#21918C", "goldenrod1")) + scale_color_manual(values = c("#21918C", "goldenrod1")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = convert53$col, size = 12), axis.text.y = element_text(size = 15)) + xlab("") + ylab("IEG Score")
dev.off()

# Neurogen
pdf("~/research/brain/results/neurogen_score_umap.pdf", width = 6, height = 6)
print(FeaturePlot(bb, "neurogen_score", order = T)+ NoLegend() + ggtitle("") + scale_color_viridis_c()) 
dev.off()

neurogen_df15 = data.frame(score = bb$neurogen_score, cond = bb$cond, old = bb$seuratclusters15, new = convert15$new.full[match(bb$seuratclusters15, convert15$old)], level = "15")
neurogen_df53 = data.frame(score = bb$neurogen_score, cond = bb$cond, old = bb$seuratclusters53, new = convert53$new[match(bb$seuratclusters53, convert53$old)], level = "53")

neurogen_df = rbind(neurogen_df15, neurogen_df53)
neurogen_df$new = factor(neurogen_df$new, levels = rev(convert_all$cluster))
neurogen_df$cond = factor(neurogen_df$cond, levels = c("CTRL", "BHVE"))
xmin = 0.2
xmax = 0.5
pdf("~/research/brain/results/neurogen_by_cluster_cond.pdf", width = 5, height = 15)
ggplot(neurogen_df, aes(x = new, y = score, color = cond, fill = cond)) + geom_split_violin(alpha = 0.3) + theme_classic() + stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.75, position = position_dodge(width = .75)) + coord_flip() + scale_fill_manual(values = c("#21918C", "goldenrod1")) + scale_color_manual(values = c("#21918C", "goldenrod1")) + theme(axis.text.y = element_text(color = rev(convert_all$color)))
dev.off()

neurogen_df15$new = factor(neurogen_df15$new, levels = convert15$new.full)
neurogen_df15$cond = factor(neurogen_df15$cond, levels = c("CTRL", "BHVE"))
pdf("~/research/brain/results/neurogen_by_cluster15_cond.pdf", width = 8, height = 5)
ggplot(neurogen_df15, aes(x = new, y = score, color = cond, fill = cond)) + geom_split_violin(alpha = 0.3) + theme_classic() + stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.75, position = position_dodge(width = .75)) + scale_fill_manual(values = c("#21918C", "goldenrod1")) + scale_color_manual(values = c("#21918C", "goldenrod1")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = convert15$col, size = 15), axis.text.y = element_text(size = 15)) + xlab("") + ylab("Neurogenesis Score")
dev.off()

neurogen_df53$new = factor(neurogen_df53$new, levels = convert53$new)
neurogen_df53$cond = factor(neurogen_df53$cond, levels = c("CTRL", "BHVE"))
pdf("~/research/brain/results/neurogen_by_cluster53_cond.pdf", width = 14, height = 5)
ggplot(neurogen_df53, aes(x = new, y = score, color = cond, fill = cond)) + geom_split_violin(alpha = 0.3) + theme_classic() + stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.75, position = position_dodge(width = .75)) + scale_fill_manual(values = c("#21918C", "goldenrod1")) + scale_color_manual(values = c("#21918C", "goldenrod1")) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = convert53$col, size = 12), axis.text.y = element_text(size = 15)) + xlab("") + ylab("Neurogenesis Score")
dev.off()


# PCA
pca.df = data.frame(x = bb@reductions$pca@cell.embeddings[,1], y = bb@reductions$pca@cell.embeddings[,2], reads = bb$reads, ncount = bb$nCount_RNA, nfeature = bb$nFeature_RNA, subject = bb$subsample.new, sample = bb$sample, run = bb$run, cond = bb$cond)
pca.sum.df = data.frame(subject = unique(bb$subsample.new), x = 0, y = 0)
pca.sum.df$sample = pca.df$sample[match(pca.sum.df$subject, pca.df$subject)]
for (i in 1:nrow(pca.sum.df)) {
  pca.sum.df$x[i] = median(pca.df$x[which(pca.df$subject == pca.sum.df$subject[i])])
  pca.sum.df$y[i] = median(pca.df$y[which(pca.df$subject == pca.sum.df$subject[i])])
}

pdf("~/research/brain/results/pca_subject.pdf", width = 6, height = 6)
print(ggplot(pca.df, aes(x = x, y = y, color = subject)) + geom_point() + xlab("PCA 1") + ylab("PCA 2") + theme_classic() + NoLegend()) 
dev.off()

pdf("~/research/brain/results/pca_median.pdf", width = 6, height = 6)
# + stat_ellipse() 
print(ggplot(pca.sum.df, aes(x = x, y = y, color = sample, group = sample)) + geom_point() + xlab("PCA 1") + ylab("PCA 2") + theme_classic() + NoLegend() + geom_text_repel(aes(label = subject), size = 6)) 
dev.off()

pdf("~/research/brain/results/ncount_nfeature_subject.pdf", width = 6, height = 6)
print(ggplot(pca.df, aes(x = ncount, y = nfeature, color = subject)) + geom_point() + xlab("Number of UMIs") + ylab("Number of Genes") + theme_classic() + NoLegend()) 
dev.off()

pdf("~/research/brain/results/ncount_read_sample.pdf", width = 6, height = 6)
print(ggplot(pca.df, aes(x = ncount, y = reads, color = sample)) + geom_point() + xlab("Number of UMIs") + ylab("Number of Reads") + theme_classic() + NoLegend()) 
dev.off()

# ViolinPlot
convert15 = convert15[order(as.numeric(convert15$new.num), decreasing = F),]
Idents(bb) = factor(convert15$new.full[match(bb$seuratclusters15, convert15$old)], levels = convert15$new.full)
bb15 = read.csv("~/research/brain/data/bb_all_sig_cluster_15_degs.csv")
bb15$pct.dif = bb15$pct.1 - bb15$pct.2
bb15 = bb15[which(bb15$avg_logFC > 0 & bb15$pct.dif > 0),]
bb15$fc.pct = bb15$avg_logFC * bb15$pct.dif
bb15.m4 = data.frame(cluster = unlist(lapply(0:14, function(x) rep(x, 4))), new = "", col = "", gene = "", reason = rep(c("Top DEG", "Top logFC", "Top % Dif", "Top FC * % Dif")), num.j = 0, logFC = 0, pct.dif = 0, rank = 0)
bb15.m4[, c("new", "col")] = convert15[match(bb15.m4$cluster, convert15$old), c("new.full", "col")]
for (i in 0:14) {
  bb15.i = bb15[which(bb15$cluster == i),]
  # bb15.i = bb15.i[order(bb15.i$p_val_adj, decreasing = F),]
  bb15.i$rank = 1:nrow(bb15.i)
  genes.to.add = c("", "", "", "")
  num.j = c(1, 1, 1, 1)
  
  # 1. Top Seurat DEG
  genes.to.add[1] = bb15.i$gene[1]
  
  # 2. Top logFC
  j = 1
  bb15.i = bb15.i[order(bb15.i$avg_logFC, decreasing = T),]
  while (j == 1 || genes.to.add[2] %in% genes.to.add[1]) {
    genes.to.add[2] = bb15.i$gene[j]
    j = j + 1
  }
  num.j[2] = j -1
  
  # 3. Highest % Diff
  j = 1
  bb15.i = bb15.i[order(bb15.i$pct.dif, decreasing = T),]
  while (j == 1 || genes.to.add[3] %in% genes.to.add[c(1,2)]) {
    genes.to.add[3] = bb15.i$gene[j]
    j = j + 1
  }
  num.j[3] = j -1

  # 4. LogFC * % Diff
  j = 1
  bb15.i = bb15.i[order(bb15.i$fc.pct, decreasing = T),]
  while (j == 1 || genes.to.add[4] %in% genes.to.add[c(1,2,3)]) {
    genes.to.add[4] = bb15.i$gene[j]
    j = j + 1
  }
  num.j[4] = j -1
  
  bb15.m4$gene[which(bb15.m4$cluster == i)] = genes.to.add
  bb15.m4$num.j[which(bb15.m4$cluster == i)] = num.j
  bb15.m4[which(bb15.m4$cluster == i), c("rank", "logFC", "pct.dif")] = bb15.i[match(genes.to.add, bb15.i$gene), c("rank", "avg_logFC", "pct.dif")]
}

bb15.m4$new = factor(bb15.m4$new, levels = convert15$new.full)
bb15.m4 = bb15.m4[order(bb15.m4$new),]
my_cols = unique(bb15.m4$col)
names(my_cols) = convert15$new.full
bb15.m4$hgnc = tolower(gene_info$human[match(bb15.m4$gene, gene_info$mzebra)])
bb15.m4$hgnc[which( is.na(bb15.m4$hgnc) )] = bb15.m4$gene[which( is.na(bb15.m4$hgnc) )]
pdf("~/research/brain/results/supplement/vlnplot_deg.pdf", height = 16, width = 8)
StackedVlnPlot(obj = bb, features = bb15.m4$gene, xcols = unique(bb15.m4$col), ycols = bb15.m4$col, cols = my_cols, sec.axis.names = bb15.m4$hgnc)
dev.off()

# DEGs by LG
deg.list = list()
deg.list[["15"]] = read.csv("~/research/brain/data/bb_all_cluster_15_degs.csv")
deg.list[["53"]] = read.csv("~/research/brain/data/bb_all_cluster_53_degs.csv")
gtf$lg = lgConverter(gtf$V1)
gtf$lg = factor(gtf$lg, levels = unique(gtf$lg))

gtf.meta = aggregate(gene_name ~ lg, gtf, length)
gtf$lg.num.genes = gtf.meta$gene_name[match(gtf$lg, gtf.meta$lg)]

for (clevel in c("15", "53")) {
  deg = deg.list[[clevel]] 
  gtf.cols = c("lg", "lg.num.genes", "V4", "V5")
  deg[, gtf.cols] = gtf[match(deg$X, gtf$gene_name), gtf.cols]
  deg = deg[which(deg$p_val_adj < 0.05),]
  
  cluster.df = data.frame()
  for (tc in unique(deg$cluster)) {
    cdeg = deg[which(deg$cluster == tc),]
    num.deg = nrow(cdeg)
    num.deg.lg = aggregate(X ~ lg, cdeg, length, drop = F)[,2]
    cluster.df = rbind(cluster.df, data.frame(clevel = clevel, cluster = tc, lg = gtf.meta$lg, num.deg = num.deg, num.deg.lg = num.deg.lg, num.lg.genes = gtf.meta$gene_name))
  }
  
  cluster.df$num.deg.lg[which(is.na(cluster.df$num.deg.lg))] = 0
  cluster.df$cluster = factor(cluster.df$cluster, levels = unique(deg$cluster))
  cluster.df$num.deg.lg.norm.c = cluster.df$num.deg.lg / cluster.df$num.deg # Normalize by number of cluster DEGs
  cluster.df$num.deg.lg.norm.lg = cluster.df$num.deg.lg / cluster.df$num.lg.genes # Normalize by number of LG genes
  cluster.df$num.deg.lg.norm.c.lg = cluster.df$num.deg.lg.norm.c / cluster.df$num.lg.genes # Normalize by both
  
  ggplot(cluster.df, aes(x = lg, y = cluster, fill = num.deg.lg))           + geom_tile() + scale_fill_viridis_c() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + coord_fixed() + ggtitle("Raw")
  ggplot(cluster.df, aes(x = lg, y = cluster, fill = num.deg.lg.norm.c))    + geom_tile() + scale_fill_viridis_c() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + coord_fixed() + ggtitle("Normalize by # of Cluster DEGs")
  ggplot(cluster.df, aes(x = lg, y = cluster, fill = num.deg.lg.norm.lg))   + geom_tile() + scale_fill_viridis_c() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + coord_fixed() + ggtitle("Normalize by # of LG Genes")
  ggplot(cluster.df, aes(x = lg, y = cluster, fill = num.deg.lg.norm.c.lg)) + geom_tile() + scale_fill_viridis_c() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + coord_fixed() + ggtitle("Normalize by Both")
  
}

# GOI in primary and secondary clusters
goi = read.csv("~/Downloads/goi_less_10k_122821_by_cat_hgnc.csv")
goi$cat.col = plyr::revalue(goi$cat, replace = c("receptor" = "goldenrod2", "tf" = "#FE04FF", "ligand" = "#002DD1", "other" = "#00E7EC"))
goi.cluster.df = data.frame()
for (i in 2:nrow(convert_all)) {
  if (convert_all$level[i] == "primary")   { bb$cluster = bb$seuratclusters15 }
  if (convert_all$level[i] == "secondary") { bb$cluster = bb$seuratclusters53 }
  this.cells = colnames(bb)[which(bb$cluster == convert_all$old[i])]
  other.cells = colnames(bb)[which(bb$cluster != convert_all$old[i])]
  this.df = pct_dif_avg_logFC(bb, cells.1 = this.cells, cells.2 = other.cells, features = goi$mzebra)
  this.df$new = convert_all$cluster[i]
  this.df$level = convert_all$level[i]
  this.df$old = convert_all$old[i]
  goi.cluster.df = rbind(goi.cluster.df, this.df)
}

goi.cluster.df$new = factor(goi.cluster.df$new, levels = convert_all$cluster[2:nrow(convert_all)])
goi.cluster.df$genes = factor(goi.cluster.df$genes, levels = goi$mzebra[order(goi$cat)])
goi.cluster.df$cat = goi$cat[match(goi.cluster.df$genes, goi$mzebra)]

goi.cluster.df$avg_logFC[which(goi.cluster.df$avg_logFC < 0)] = 0
goi.cluster.df$pct_dif[which(goi.cluster.df$pct_dif < 0)] = 0

pdf("~/research/brain/results/goi_genes_in_clusters.pdf", width = 15, height = 20)
ggplot(goi.cluster.df, aes(x = new, y = genes, size = pct_dif, color = avg_logFC)) + geom_point() + scale_color_viridis_c() + theme(axis.text.y = element_text(color = goi$cat.col[order(goi$cat)]), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = convert_all$color[2:nrow(convert_all)]))
dev.off()

# Local Ztest investigation
zpops = data.frame(level = c("secondary", "secondary", "secondary", "secondary", "primary", "secondary", "secondary", "all", "all", "all", "all", "all", "all"), 
                   cluster = c(31, 5, 20, 17, 1, 21, 22, "all", "all", "all", "all", "all", "all"), 
                   label = c("2.1_OPC", "1.1_RGC", "1.2_RGC", "5.2_GABA th+", "9_Glut hrh3+ (LOC101478323)", "4.3_GABA", "4.4_GABA", "adrb1+", "ghr+", "nr3c2+ (LOC101466282)", "ghrhr+", "esr2+", "gad2+"), 
                   gene = c("all", "all", "all", "all", "all", "all", "all", "adrb1", "LOC101464862", "LOC101466282", "ghrhr", "esr2", "gad2"))
df = data.frame()
for (i in 1:nrow(zpops)) {
  # for (i in 1) {
  print(i)
  level.i = zpops$level[i]
  cluster.i = zpops$cluster[i]
  label.i = zpops$label[i]
  gene.i = zpops$gene[i]
  for (pcrc.gene in pcrc) {
    bb$cluster = bb$seuratclusters15
    if (level.i == "secondary") { bb$cluster = bb$seuratclusters53 }
    if (level.i == "all")       { bb$cluster = "all"               }
    bb$gene = T
    if (gene.i != "all") { bb$gene = bb@assays$RNA@counts[gene.i,] > 0 }
    bb$pcrc.gene = bb@assays$RNA@counts[pcrc.gene,] > 0
    zpop = colnames(bb)[which(bb$cluster == cluster.i & bb$gene & bb$pcrc.gene)]
    df = rbind(df, data.frame( i = i, level = level.i, cluster = cluster.i, gene = gene.i, pcrc.gene = pcrc.gene, num.cells = length(zpop), this.z = sum(1/bb$nFeature_RNA[zpop]), other.z = sum(1/bb$nFeature_RNA[! colnames(bb) %in% zpop]) ))
  }
}

df$this.v.other = df$this.z / df$other.z
# ggplot(df, aes(x = label.i, y = this.v.other, fill = pcrc.gene)) + geom_bar(stat = 'identity') + ylab("[PCRC Gene Present / # of Genes] in Pop vs [PCRC Gene Present / # of Genes] NOT in Pop") + NoLegend()
write.csv(df, "C:/Users/miles/Downloads/pcrc_genes_driving_effects.csv")

mat = readRDS("~/scratch/brain/data/scgnn_imputed.rds")
test = as.data.frame(mat)
big_mat = data.frame()
for (i in 0:14) {
  print(i)
  test[, c("subsample", "inc")] = c(bb$subsample, bb$seuratclusters15 == i)
  test2 = aggregate(. ~ subsample + inc, test, mean)
  big_mat = rbind(big_mat, test2[which(test2$inc == F),])
}
