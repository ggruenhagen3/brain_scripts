rna_path = "C:/Users/miles/Downloads/brain/"
# rna_path = "~/research/brain/"
bb = readRDS(paste0(rna_path, "data/bb_subsample_02222021.RDS"))
source(paste0(rna_path, "brain_scripts/all_f.R"))

#=====================================================================================
# Figure 1 ===========================================================================
#=====================================================================================

# *** 1B. 15 cluster and 53 cluster level UMAP *** #
# Define 15 level colors
for (k in 1:15) {
print(k)
library("scales")
cols15 = gc.ramp <- hue_pal()(15)
cols15 = cols15[c(10:15, 1:9)]

# names(cols15) = 0:14

# Change Cluster IDs on 53 Level
convert53 = data.frame(old = 0:52, new = c("4.1_In", "10.1_Ex", "15.1_In/Ex", "9.1_Ex", "8.1_Ex", "1.1_Astro", "6_In", "5.1_In", "9.2_Ex", "8.2_Ex", "15.2_In", "11.1_Ex", "8.3_Ex", "8.4_Ex", "9.3_Ex", "4.2_In", "8.5_Ex", "5.2_In", "8.6_Ex", "8.7_Ex", "1.2_Astro", "4.3_In", "4.4_In", "9.4_Ex", "9.5_Ex", "8.8_Ex", "9.6_Ex", "4.5_In", "12_Ex", "8.9_Ex", "10.2_Ex", "2.1_OPC", "15.3_In", "11.2_Ex", "15.4_In", "4.6_In", "9.7_Ex", "13_Ex", "14_Ex", "4.7_In", "11.3_Ex", "9.8_Ex", "8-9_Ex", "15.5_In/Ex", "4.8_In", "1.3_MG", "2.2_Oligo", "15.6_Ex", "8.10_Ex", "8.11_Ex", "3_Peri", "15.7_Ex", "7_In"))
convert53$new = gsub("In", "GABA", convert53$new)
convert53$new = gsub("Ex", "Glut", convert53$new)
convert53$new2 = convert53$new

# Split the parent and subcluster apart
convert53 = separate(data = convert53, col = new2, into = c("new.id", "new.gaba"), sep = "_")
convert53 = separate(data = convert53, col = new.id, into = c("new.parent", "new.sub"), sep = "\\.")
# Assign a parent and subcluster to the 8-9_Glut special case
convert53$new.parent[which(convert53$old == 42)] = 8
convert53$new.sub[which(convert53$old == 42)] = 12
# Change Cluster IDs on 15 Level
convert15 = data.frame(old = 0:14, new = c("8_Ex", "9_Ex", "4_In", "15_In/Ex", "1_Astro/MG", "10_Ex", "5_In", "11_Ex", "6_In", "2_OPC/Oligo", "12_Ex", "13_Ex", "14_Ex", "3_Peri", "7_In"))
convert15$new = gsub("In", "GABA", convert15$new)
convert15$new = gsub("Ex", "Glut", convert15$new)
convert15$new.full = convert15$new
convert15 = separate(data = convert15, col = new, into = c("new.num", "new.junk"), sep = "_")
convert15 = convert15[order(as.numeric(convert15$new.num), decreasing = T),]
convert15$col = cols15
# Get the parent and color of the 53 cluster 
convert53$new.parent.old = convert15$old[match(convert53$new.parent, convert15$new.num)]
convert53$new.parent.old.col = convert15$col[match(convert53$new.parent, convert15$new.num)]
convert53$row = rownames(convert53)
convert53$col = convert53$new.parent.old.col

# Lighten and Darken Parent Color
for (clust15_new_num in 1:15) {
  this_rows = convert53[which(convert53$new.parent == clust15_new_num),]
  max_sub = max(as.numeric(this_rows$new.sub))
  if (! is.na(max_sub)) {
    this_col = this_rows$new.parent.old.col[1]
    mod_values = seq(0.1, 1, by = 0.1)
    mod_values = sort(rep(mod_values, 2))
    mod_values = mod_values[1:max_sub]
    new_cols = c()
    for (i in 1:max_sub) {
      mod_col = lighten(this_col, mod_values[i])
      if (i %% 2 == 0) {
        #darken half the colors
        mod_col = darken(this_col, mod_values[i])
      } 
      new_cols = c(new_cols, mod_col)
    } # end max_sub for
    convert53[as.numeric(this_rows$row), "col"] = new_cols
  }
} # end outer for
convert53 = convert53[order(as.numeric(convert53$new.parent), decreasing = F),]
tmp = convert53[which(convert53$new == "8-9_Glut"),]
convert53 = convert53[-which(convert53$new == "8-9_Glut"),]
convert53 = rbind(convert53[1:29,], tmp, convert53[30:nrow(convert53),])

bb$good_names = convert53$new[match(bb$seuratclusters53, convert53$old)]
bb$good_names_num = colsplit(bb$good_names, "_", c("num", "ex"))[,1]
Idents(bb) = bb$good_names_num
clust53_new_col_list2 = convert53$col
names(clust53_new_col_list2) = colsplit(convert53$new, "_", c("num", "ex"))[,1]

pdf (paste0("~/research/brain/results/cols_", k, "_rev.pdf"), width = 5, height = 5)
print(DimPlot(bb, label = T, pt.size = 1) + scale_color_manual(values = clust53_new_col_list2) + NoLegend())
dev.off()
}

pdf("C:/Users/miles/Downloads/brain/results/bb/53_clusters_on_15_light_dark_label_04012021.pdf", width = 12, height = 12)
print(DimPlot(bb, label = T, pt.size = 1) + scale_color_manual(values = clust53_new_col_list2) + NoLegend())
dev.off()

pdf("C:/Users/miles/Downloads/brain/results/bb/53_clusters_on_15_light_dark_no_label_04012021.pdf", width = 12, height = 12)
print(DimPlot(bb, label = F, pt.size = 1) + scale_color_manual(values = clust53_new_col_list2) + NoLegend())
dev.off()

# *** 1E. DotPlot *** #

# Read in Cluster Specific Markers
specific = read.csv(paste0(rna_path, "results/specific_markers_15_10_02_12821.csv"))
specific = specific[-which(specific$gene == "LOC101464395"),]
specific = specific[-which(specific$gene == "LOC101478779"),]
specific = specific[-which(specific$gene == "zic1"),]
specific = specific[-which(specific$gene == "LOC101464181"),] # remove duplicated gene
specific = specific[-which(specific$gene == "LOC101475710"),] # remove duplicated gene
specific = specific[-which(specific$gene == "LOC112435918"),] # trying to pick genes without LOC from cluster 5
specific = specific[-which(specific$gene == "LOC112435917"),] # trying to pick genes without LOC from cluster 5
deg_table = table(specific$gene)
specific$n_gene_appears = deg_table[match(specific$gene, names(deg_table))]

# Select Top 2 Markers
specific = specific[which(specific$n_gene_appears == 1),]
specific2 = as.data.frame(specific %>% group_by(cluster) %>% slice(head(row_number(), 2)))
specific2$col = convert15$col[match(specific2$cluster, convert15$old)]
specific2$new.full = convert15$new.full[match(specific2$cluster, convert15$old)]
specific2$new.num = as.numeric(convert15$new.num[match(specific2$cluster, convert15$old)])
specific2 = specific2[order(specific2$new.num, decreasing = F),]
specific2$labels = specific2$gene
specific2$labels[which(startsWith(specific2$gene, "LOC") & !is.na(specific2$hgnc))] = specific2$hgnc[which(startsWith(specific2$gene, "LOC") & !is.na(specific2$hgnc))]
specific2$labels = tolower(specific2$labels)

# Get the List of Genes to use in the plot
top2 = specific2$gene
top2_labels = specific2$labels
known_markers = rev(c("fabp7", "LOC101468392", "LOC101483038", "olig2", "lama2", "LOC101478779", "LOC101469680", "syp", "gad2", "gad1", "slc17a6", "LOC101484681"))
known_labels = rev(c("fabp7", "gfap", "sox10", "olig2", "lama2", "slc6a13", "map2", "syp", "gad2", "gad1", "slc17a6", "slc17a7"))
# known_markers = rev(c("fabp7", "LOC101468392", "LOC101483038", "olig2", "lama2", "LOC101476633", "syp", "gad2", "gad1", "slc17a6", "LOC101484681"))
# known_labels = rev(c("fabp7", "GFAP", "sox10", "olig2", "lama2", "map2", "syp", "gad2", "gad1", "slc17a6", "slc17a7"))
all_genes = c(top2, known_markers)
all_labels = c(top2_labels, known_labels)
bb = ScaleData(bb, features = all_genes)

# Set colors for plot
xtext_unique_cols = unique(specific2$col)
xtext_cols = c( as.vector(specific2$col), rep("black", length(known_markers)) )
pal = rev(brewer.pal(7, "RdYlBu"))
pal = colorRampPalette(pal)

# Convert cluster IDs
bb$good15_names = convert15$new.full[match(bb$seuratclusters15, convert15$old)]
bb$good15_names = factor(bb$good15_names, levels = convert15$new.full)
Idents(bb) = bb$good15_names
dotp = DotPlot(bb, features = all_genes) + scale_color_gradientn(colors = pal(50)) + ylab("Cluster") + xlab("Gene") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = xtext_cols, face = "italic"), axis.text.y = element_text(colour = xtext_unique_cols)) + scale_x_discrete(labels = all_labels)
print(dotp)

# DotPlot
pdf(paste0(rna_path, "/results/known_markers_w_specific_dotplot.pdf"), width = 12, height = 5)
print(dotp)
dev.off()

pal = rev(brewer.pal(7, "RdYlBu"))
pal = colorRampPalette(pal)
# names(zpcrc) = c("4_GABA", "6_GABA", "14_Glut", "13_Glut", "4_GABA", "12_Glut", "3_Peri", "6_GABA", "9_Glut", "1_Astro/MG")
zpcrc_deg = deg15[which(deg15$gene %in% zpcrc),]
names(zpcrc) = sapply(zpcrc, function(x) { 
  this_df = zpcrc_deg[which(zpcrc_deg$gene == x),]
  return(this_df$cluster[which.min(this_df$p_val_adj)])
})
names(zpcrc) = convert15$new.full[match(names(zpcrc), convert15$old)]
names(zpcrc)[which(is.na(names(zpcrc)))] = "2_OPC/Oligo"
zpcrc = zpcrc[order(match(names(zpcrc), convert15$new.full))]
xtext_cols = convert15$col[match(names(zpcrc), convert15$new.full)]
names(zpcrc) = NULL
pdf("C:/Users/miles/Downloads/brain/results/pcrc_cluster_markers.pdf", width = 12, height = 5)
print(DotPlot(bb, features = zpcrc) + scale_color_gradientn(colors = pal(50)) + ylab("Cluster") + xlab("Gene") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic", colour = xtext_cols), axis.text.y = element_text(colour = convert15$col)) + scale_x_discrete(labels = zpcrc))
dev.off()

# Volcano Plots
library("EnhancedVolcano")
bvc15 = read.csv("C:/Users/miles/Downloads/bb15_combined_diff_markers_by_sample_111720_pqz.csv", stringsAsFactors = F)
bvc53 = read.csv("C:/Users/miles/Downloads/bb53_combined_diff_markers_by_sample_012620_pqz.csv", stringsAsFactors = F)
bvcbulk = read.csv("C:/Users/miles/Downloads/real_combined_diff_markers_by_condition_by_pair1235_121820_pqz.csv", stringsAsFactors = F)

# bvc53_sig = read.csv("C:/Users/miles/Downloads/bb53_combined_diff_markers_by_sample_012620_pqz.csv")
bvc15_sig = read.csv("C:/Users/miles/Downloads/brain/results/bb/bb15_combined_all_pairs_strict_hits_120920_pct_logFC.csv", stringsAsFactors = F)
bulk_sig = read.csv("C:/Users/miles/Downloads/combined_all_all_replicates_bulk_hits_121620.csv")

# Volcano plot for all 15 level results
bvc15_vp = data.frame()
for (i in unique(bvc15_sig$cluster)) {
  cluster_i = i
  bvc15_0 = bvc15[which(bvc15$cluster == cluster_i),]
  bvc15_0_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$seuratclusters15 == cluster_i & bb$cond == "BHVE")], cells.2 = colnames(bb)[which(bb$seuratclusters15 == cluster_i & bb$cond == "CTRL")])
  bvc15_0$avg_logFC = bvc15_0_stats$avg_logFC[match(bvc15_0$gene_id, bvc15_0_stats$genes)]
  bvc15_0$pct_dif = bvc15_0_stats$pct_dif[match(bvc15_0$gene_id, bvc15_0_stats$genes)]
  bvc15_0$same_dir = sign(bvc15_0$p1) == sign(bvc15_0$p2) & sign(bvc15_0$p1) == sign(bvc15_0$p3) & sign(bvc15_0$p1) == sign(bvc15_0$p4) & sign(bvc15_0$p1) == sign(bvc15_0$p5)
  bvc15_0$Significant = "NA"
  bvc15_0$Significant[which(bvc15_0$same_dir)] = "False"
  bvc15_0$Significant[which(bvc15_0$same_dir & bvc15_0$p1 < 0 & bvc15_0$gene_id %in% bvc15_sig$mzebra[which(bvc15_sig$cluster == cluster_i)])] = "BHVE"
  bvc15_0$Significant[which(bvc15_0$same_dir & bvc15_0$p1 > 0 & bvc15_0$gene_id %in% bvc15_sig$mzebra[which(bvc15_sig$cluster == cluster_i)])] = "CTRL"
  bvc15_0$Significant = factor(bvc15_0$Significant, levels=c("BHVE", "CTRL", "False", "NA"))
  bvc15_0$hgnc = gene_info$human[match(bvc15_0$gene_id, gene_info$mzebra)]
  bvc15_0$label = bvc15_0$hgnc
  bvc15_0$label[which(is.na(bvc15_0$label))] = bvc15_0$gene_id[which(is.na(bvc15_0$label))]
  bvc15_0$label[which(bvc15_0$gene_id == "rsrp1")] = "rsrp1"
  
  bvc15_0 = bvc15_0[order(bvc15_0$Significant, decreasing = T),]
  bvc15_0$top_quantile = bvc15_0$Significant %in% c("BHVE", "CTRL") & bvc15_0$q %in% head(sort(bvc15_0$q[which(bvc15_0$Significant %in% c("BHVE", "CTRL"))]),5)
  bvc15_0$cluster = cluster_i
  bvc15_0$new = convert15$new.full[match(bvc15_0$cluster, convert15$old)]
  bvc15_0$col = convert15$col[match(bvc15_0$cluster, convert15$old)]
  bvc15_vp = rbind(bvc15_vp, bvc15_0)
}

pdf(paste0("C:/Users/miles/Downloads/bvc15.pdf"), width = 7, height = 6)
print(ggplot(bvc15_vp, aes(x = avg_logFC, y = -log10(p), color = Significant, size = abs(pct_dif), alpha = Significant)) + geom_point() + scale_color_manual(values = c(viridis(3)[c(3, 1, 2)], "gray90"), drop = F) + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" P")) + scale_size(range = c(0.05, 5), name = "% Difference") + scale_alpha_manual(values = c(0.9, 0.9, 0.5, 0.5), drop = F) + geom_text_repel(data=bvc15_vp[which(bvc15_vp$top_quantile),], aes(label = label), color = "black", size = 4, show.legend = FALSE) + scale_y_sqrt() + theme_light())
dev.off()

# bvc15_vp$Significant = as.vector(bvc15_vp$Significant)
# bvc15_vp$Significant[which(bvc15_vp$Significant %in% c("BHVE", "CTRL"))] = bvc15_vp$new[which(bvc15_vp$Significant %in% c("BHVE", "CTRL"))]
# bvc15_vp$Significant[which( ! bvc15_vp$Significant %in% convert15$new.full)] = NA
# bvc15_vp$col[which(is.na(bvc15_vp$Significant))] = "gray90"
bvc15_vp$col_sig = "gray90"
bvc15_vp$col_sig[which( bvc15_vp$Significant != "NA" )] = bvc15_vp$col[which( bvc15_vp$Significant != "NA" )]
bvc15_vp$alpha = 0.58
bvc15_vp$alpha[which(bvc15_vp$Significant == "False")] = 0.6
bvc15_vp$alpha[which(bvc15_vp$Significant %in% c("BHVE", "CTRL"))] = 0.9
bvc15_vp = bvc15_vp[order(bvc15_vp$Significant, decreasing = T),]
pdf(paste0("C:/Users/miles/Downloads/bvc15_color.pdf"), width = 7, height = 6)
# print(ggplot(bvc15_vp, aes(x = avg_logFC, y = -log10(p), color = col, size = abs(pct_dif), alpha = Significant)) + geom_point() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" P")) + scale_size(range = c(0.05, 5), name = "% Difference") + scale_alpha_manual(values = c(0.9, 0.9, 0.5, 0.5), drop = F) + geom_text_repel(data=bvc15_vp[which(bvc15_vp$top_quantile),], aes(label = label), color = "black", size = 4, show.legend = FALSE) + scale_y_sqrt() + theme_light() + scale_color_identity())
print(ggplot(bvc15_vp, aes(x = avg_logFC, y = -log10(p), color = col_sig, size = abs(pct_dif), alpha = alpha)) + geom_point() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" P")) + scale_size(range = c(0.05, 5), name = "% Difference") + geom_text_repel(data=bvc15_vp[which(bvc15_vp$top_quantile),], aes(label = label), color = "black", size = 4, show.legend = FALSE) + scale_y_sqrt() + theme_light() + scale_color_identity())
dev.off()

pdf(paste0("C:/Users/miles/Downloads/bvc15_", cluster_i, "_s.pdf"), width = 7, height = 6)
# ggplot(bvc15_0, aes(x = avg_logFC, y = -log10(p), color = Significant, size = abs(pct_dif), alpha = Significant)) + geom_point(alpha = 0.5) + scale_color_manual(values = c("gold", "#006769")) + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(Log["10"]*" P")) + scale_size(range = c(0.05, 5), name = "% Difference") + scale_alpha_manual(values = c(0.9, 0.5), name="") + geom_text_repel(data=top_quantile, aes(label = label), color = "black", size = 4, show.legend = FALSE) + scale_y_sqrt() + theme_light()
print(ggplot(bvc15_0, aes(x = avg_logFC, y = -log10(p), color = Significant, size = abs(pct_dif), alpha = Significant)) + geom_point() + scale_color_manual(values = c(viridis(3)[c(3, 1, 2)], "gray90"), drop = F) + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" P")) + scale_size(range = c(0.05, 5), name = "% Difference") + scale_alpha_manual(values = c(0.9, 0.9, 0.5, 0.5), drop = F) + geom_text_repel(data=top_quantile, aes(label = label), color = "black", size = 4, show.legend = FALSE) + scale_y_sqrt() + theme_light())
dev.off()

bulk_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$cond == "BHVE")], cells.2 = colnames(bb)[which(bb$cond == "CTRL")])
bvcbulk$same_dir = sign(bvcbulk$p1) == sign(bvcbulk$p2) & sign(bvcbulk$p1) == sign(bvcbulk$p3) & sign(bvcbulk$p1) == sign(bvcbulk$p4) & sign(bvcbulk$p1) == sign(bvcbulk$p5)
# bvcbulk$avg_logFC = bulk_stats_s$avg_logFC[match(bvcbulk$gene_id, bulk_stats_s$genes)]
# bvcbulk$pct_dif = bulk_stats_s$pct_dif[match(bvcbulk$gene_id, bulk_stats_s$genes)]
bvcbulk$avg_logFC = bulk_stats$avg_logFC[match(bvcbulk$gene_id, bulk_stats$genes)]
bvcbulk$pct_dif = bulk_stats$pct_dif[match(bvcbulk$gene_id, bulk_stats$genes)]
bvcbulk$Significant = "NA"
bvcbulk$Significant[which(bvcbulk$same_dir)] = "False"
bvcbulk$Significant[which(bvcbulk$same_dir & bvcbulk$p1 < 0 & bvcbulk$gene_id %in% bulk_sig$gene_id)] = "BHVE"
bvcbulk$Significant[which(bvcbulk$same_dir & bvcbulk$p1 > 0 & bvcbulk$gene_id %in% bulk_sig$gene_id)] = "CTRL"
bvcbulk$Significant = factor(bvcbulk$Significant, levels=c("BHVE", "CTRL", "False", "NA"))
bvcbulk$hgnc = gene_info$human[match(bvcbulk$gene_id, gene_info$mzebra)]
bvcbulk$label = bvcbulk$hgnc
bvcbulk$label[which(is.na(bvcbulk$label))] = bvcbulk$gene_id[which(is.na(bvcbulk$label))]
bvcbulk$label[which(bvcbulk$gene_id == "rsrp1")] = "rsrp1"

bvcbulk = bvcbulk[order(bvcbulk$Significant, decreasing = T),]
top_quantile = bvcbulk[which(bvcbulk$Significant %in% c("BHVE", "CTRL") & bvcbulk$bh %in% head(sort(bvcbulk$bh[which(bvcbulk$Significant %in% c("BHVE", "CTRL"))]),5)),]

pdf(paste0("C:/Users/miles/Downloads/bvcbulk.pdf"), width = 7, height = 6)
# ggplot(bvcbulk, aes(x = avg_logFC, y = -log10(p), color = Significant, size = abs(pct_dif), alpha = Significant)) + geom_point(alpha = 0.5) + scale_color_manual(values = c("gold", "#006769")) + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(Log["10"]*" P")) + scale_size(range = c(0.05, 5), name = "% Difference") + scale_alpha_manual(values = c(0.9, 0.5), name="") + geom_text_repel(data=top_quantile, aes(label = label), color = "black", size = 4, show.legend = FALSE) + scale_y_sqrt() + theme_light()
ggplot(bvcbulk, aes(x = avg_logFC, y = -log10(p), color = Significant, size = abs(pct_dif), alpha = Significant)) + geom_point() + scale_color_manual(values = c(viridis(3)[c(3, 1, 2)], "gray90")) + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" P")) + scale_size(range = c(0.05, 5), name = "% Difference") + scale_alpha_manual(values = c(0.9, 0.9, 0.5, 0.5)) + geom_text_repel(data=top_quantile, aes(label = label), color = "black", size = 4, show.legend = FALSE) + scale_y_sqrt() + theme_light()
dev.off()


# Heatmap with Hierarchical Clustering for Bulk BVC DEGs
pair1_bulk_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$sample == "b1")], cells.2 = colnames(bb)[which(bb$sample == "c1")])
pair2_bulk_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$sample == "b2")], cells.2 = colnames(bb)[which(bb$sample == "c2")])
pair3_bulk_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$sample == "b3")], cells.2 = colnames(bb)[which(bb$sample == "c3")])
pair4_bulk_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$sample == "b4")], cells.2 = colnames(bb)[which(bb$sample == "c4")])
pair5_bulk_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$sample == "b5")], cells.2 = colnames(bb)[which(bb$sample == "c5")])

bulk_common_genes = pair1_bulk_stats$genes[which( pair1_bulk_stats$genes %in% pair2_bulk_stats$genes & pair1_bulk_stats$genes %in% pair3_bulk_stats$genes & pair1_bulk_stats$genes %in% pair4_bulk_stats$genes & pair1_bulk_stats$genes %in% pair5_bulk_stats$genes )]
all_pairs_bulk_stats = data.frame(genes = bulk_common_genes,
                                  pair1 = pair1_bulk_stats$avg_logFC[match(bulk_common_genes, pair1_bulk_stats$genes)],
                                  pair2 = pair2_bulk_stats$avg_logFC[match(bulk_common_genes, pair2_bulk_stats$genes)],
                                  pair3 = pair3_bulk_stats$avg_logFC[match(bulk_common_genes, pair3_bulk_stats$genes)],
                                  pair4 = pair4_bulk_stats$avg_logFC[match(bulk_common_genes, pair4_bulk_stats$genes)],
                                  pair5 = pair5_bulk_stats$avg_logFC[match(bulk_common_genes, pair5_bulk_stats$genes)])
all_sample_logFC = data.frame(genes = all_pairs_bulk_stats$genes)
all_sample_logFC[,c("b1", "c1")] = c(all_pairs_bulk_stats$pair1, -all_pairs_bulk_stats$pair1)
all_sample_logFC[,c("b2", "c2")] = c(all_pairs_bulk_stats$pair2, -all_pairs_bulk_stats$pair2)
all_sample_logFC[,c("b3", "c3")] = c(all_pairs_bulk_stats$pair3, -all_pairs_bulk_stats$pair3)
all_sample_logFC[,c("b4", "c4")] = c(all_pairs_bulk_stats$pair4, -all_pairs_bulk_stats$pair4)
all_sample_logFC[,c("b5", "c5")] = c(all_pairs_bulk_stats$pair5, -all_pairs_bulk_stats$pair5)

rownames(all_sample_logFC) = all_sample_logFC$genes
all_sample_logFC$genes = NULL
all_sample_logFC_mat = as.matrix(all_sample_logFC)
all_sample_logFC_mat_test = all_sample_logFC_mat[bulk_sig$gene_id,]

Idents(bb) = bb$sample
bulk_avg_exp_sample = myAverageExpression(bb)
bulk_avg_exp_sample_mat = as.matrix(bulk_avg_exp_sample)
bulk_avg_exp_sample_mat_test = bulk_avg_exp_sample_mat[bulk_sig$gene_id,]

col_annot = data.frame(cond = rep(c("BHVE", "CTRL"), 5), row.names = colnames(all_sample_logFC_mat_test))

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

pdf("C:/Users/miles/Downloads/bvcbulk_sample_heatmap_logFC_viridis_test1.pdf")
# pheatmap::pheatmap(bulk_avg_exp_sample_mat_test, show_rownames = F)
mat_breaks <- quantile_breaks(all_sample_logFC_mat_test, n = 11)
# pheatmap::pheatmap(all_sample_logFC_mat_test, border_color = NA, angle_col = 0, show_rownames = F, color = viridis(length(mat_breaks) - 1), breaks = mat_breaks, legend = F, annotation_col = col_annot, annotation_colors = list(cond = c(CTRL = "#007C7D", BHVE = "gold")))
pheatmap::pheatmap(all_sample_logFC_mat_test, border_color = NA, angle_col = 0, show_rownames = F, color = "white", legend = F, annotation_col = col_annot, annotation_colors = list(cond = c(CTRL = "#007C7D", BHVE = "gold")))
dev.off()

# Modify ordering of the clusters using clustering callback option
callback = function(hc, mat){
  if (length(hc$order) > 10) {
    sv = unlist(lapply(1:nrow(mat), function(x) max(mat[x,c("b1", "b2", "b3", "b4", "b5")])))
    dend = reorder(as.dendrogram(hc), wts = sv)
    tmp = dend[[2]][[2]][[2]][[1]]
    dend[[2]][[2]][[2]][[1]] = dend[[2]][[2]][[2]][[2]]
    dend[[2]][[2]][[2]][[2]] = tmp
    tmp = dend[[1]]
    dend[[1]] = dend[[2]]
    dend[[2]] = tmp
    as.hclust(dend)
  } else {
    dendsort(hc, isReverse = T)
  }
}

# pdf("C:/Users/miles/Downloads/bvcbulk_sample_heatmap_logFC_viridis2.pdf")
png("C:/Users/miles/Downloads/bvcbulk_sample_heatmap_logFC_viridis2.png", width = 1500, height = 2000, res = 300)
# pheatmap::pheatmap(bulk_avg_exp_sample_mat_test, show_rownames = F)
mat_breaks <- quantile_breaks(all_sample_logFC_mat_test, n = 11)
# pheatmap::pheatmap(all_sample_logFC_mat_test, border_color = NA, angle_col = 0, show_rownames = F, color = viridis(length(mat_breaks) - 1), breaks = mat_breaks, legend = F, annotation_col = col_annot, annotation_colors = list(cond = c(CTRL = "#007C7D", BHVE = "gold")))
pheatmap::pheatmap(all_sample_logFC_mat_test, clustering_callback = callback, border_color = "black", angle_col = 0, show_rownames = F, color = viridis(length(mat_breaks) - 1), breaks = mat_breaks, legend = F, annotation_col = col_annot, annotation_colors = list(cond = c(CTRL = "#007C7D", BHVE = "gold")))
dev.off()

# Heatmap with Hierarchical Clustering for 15 Cluster BVC DEGs
heat15 = data.frame()
for (i in 0:14) {
  if ( length(which(bvc15_sig$cluster == i)) > 0 ) {
    print(i)
    pair1_bulk_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$sample == "b1" & bb$seuratclusters15 == i)], cells.2 = colnames(bb)[which(bb$sample == "c1" & bb$seuratclusters15 == i)])
    pair2_bulk_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$sample == "b2" & bb$seuratclusters15 == i)], cells.2 = colnames(bb)[which(bb$sample == "c2" & bb$seuratclusters15 == i)])
    pair3_bulk_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$sample == "b3" & bb$seuratclusters15 == i)], cells.2 = colnames(bb)[which(bb$sample == "c3" & bb$seuratclusters15 == i)])
    pair4_bulk_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$sample == "b4" & bb$seuratclusters15 == i)], cells.2 = colnames(bb)[which(bb$sample == "c4" & bb$seuratclusters15 == i)])
    pair5_bulk_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$sample == "b5" & bb$seuratclusters15 == i)], cells.2 = colnames(bb)[which(bb$sample == "c5" & bb$seuratclusters15 == i)])
    
    bulk_common_genes = pair1_bulk_stats$genes[which( pair1_bulk_stats$genes %in% pair2_bulk_stats$genes & pair1_bulk_stats$genes %in% pair3_bulk_stats$genes & pair1_bulk_stats$genes %in% pair4_bulk_stats$genes & pair1_bulk_stats$genes %in% pair5_bulk_stats$genes )]
    all_pairs_bulk_stats = data.frame(genes = bulk_common_genes,
                                      pair1 = pair1_bulk_stats$avg_logFC[match(bulk_common_genes, pair1_bulk_stats$genes)],
                                      pair2 = pair2_bulk_stats$avg_logFC[match(bulk_common_genes, pair2_bulk_stats$genes)],
                                      pair3 = pair3_bulk_stats$avg_logFC[match(bulk_common_genes, pair3_bulk_stats$genes)],
                                      pair4 = pair4_bulk_stats$avg_logFC[match(bulk_common_genes, pair4_bulk_stats$genes)],
                                      pair5 = pair5_bulk_stats$avg_logFC[match(bulk_common_genes, pair5_bulk_stats$genes)])
    all_sample_logFC = data.frame(genes = all_pairs_bulk_stats$genes)
    all_sample_logFC[,c("b1", "c1")] = c(all_pairs_bulk_stats$pair1, -all_pairs_bulk_stats$pair1)
    all_sample_logFC[,c("b2", "c2")] = c(all_pairs_bulk_stats$pair2, -all_pairs_bulk_stats$pair2)
    all_sample_logFC[,c("b3", "c3")] = c(all_pairs_bulk_stats$pair3, -all_pairs_bulk_stats$pair3)
    all_sample_logFC[,c("b4", "c4")] = c(all_pairs_bulk_stats$pair4, -all_pairs_bulk_stats$pair4)
    all_sample_logFC[,c("b5", "c5")] = c(all_pairs_bulk_stats$pair5, -all_pairs_bulk_stats$pair5)
    
    rownames(all_sample_logFC) = all_sample_logFC$genes
    all_sample_logFC$genes = NULL
    all_sample_logFC_mat = as.matrix(all_sample_logFC)
    all_sample_logFC_mat_test = all_sample_logFC_mat[unique(bvc15_sig$mzebra[which(bvc15_sig$cluster == i)]),] 
    
    if ( length(which(bvc15_sig$cluster == i)) == 1 ) {
      heat15 = rbind(heat15, all_sample_logFC_mat_test)
      rownames(heat15)[nrow(heat15)] = paste0(bvc15_sig$mzebra[which(bvc15_sig$cluster == i)], "_", i)
    } else {
      rownames(all_sample_logFC_mat_test) = paste0(rownames(all_sample_logFC_mat_test), "_", i)
      heat15 = rbind(heat15, all_sample_logFC_mat_test)
    }
  
  } # end if cluster has a sig hit
} # end cluster for
heat15 = as.matrix(heat15)
# png("C:/Users/miles/Downloads/bvc15_sample_heatmap_logFC_viridis2.png", width = 1500, height = 2000, res = 300)
heat15_myc = heat15
heat15_myc = cbind(heat15_myc, rowSums(heat15_myc[,c(T,F)]) > rowSums(heat15_myc[, c(F,T)]))
heat15_myc = cbind(heat15_myc, reshape2::colsplit(rownames(heat15), "_", c('a', 'b'))[,2])
heat15_myc2 = data.frame()
heat15_myc2_row_names = c()
for (new.num in 1:15) {
  i = convert15$old[which(convert15$new.num == new.num)]
  this_rows = heat15_myc[which(heat15_myc[,11] == 1 & heat15_myc[,12] == i), 1:10]
  heat15_myc2 = rbind(heat15_myc2, this_rows)
  colnames(heat15_myc2) = c("b1", "c1", "b2", "c2", "b3", "c3", "b4", "c4", "b5", "c5")
  heat15_myc2_row_names = c(heat15_myc2_row_names, rownames(heat15_myc)[which(heat15_myc[,11] == 1 & heat15_myc[,12] == i)])
}
for (new.num in 1:15) {
  i = convert15$old[which(convert15$new.num == new.num)]
  this_rows = heat15_myc[which(heat15_myc[,11] == 0 & heat15_myc[,12] == i), 1:10]
  heat15_myc2 = rbind(heat15_myc2, this_rows)
  heat15_myc2_row_names = c(heat15_myc2_row_names, rownames(heat15_myc)[which(heat15_myc[,11] == 0 & heat15_myc[,12] == i)])
}
rownames(heat15_myc2) = heat15_myc2_row_names

callback = function(hc, mat){
  if (length(hc$order) == 10) {
    # dend = as.dendrogram(hc)
    # tmp = dend[[2]][[2]][[2]][[1]]
    # dend[[2]][[2]][[2]][[1]] = dend[[2]][[2]][[2]][[2]]
    # dend[[2]][[2]][[2]][[2]] = tmp
    # tmp = dend[[1]]
    # dend[[1]] = dend[[2]]
    # dend[[2]] = tmp
    # as.hclust(dend)
    dendsort(hc, isReverse = T)
  }
}

png("C:/Users/miles/Downloads/bvc15_sample_heatmap_logFC_viridis2.png", width = 1000, height = 1500, res = 200)
mat_breaks <- quantile_breaks(heat15, n = 11)
col_annot = data.frame(cond = rep(c("BHVE", "CTRL"), 5), row.names = colnames(heat15))
row_annot = data.frame(cluster = factor(reshape2::colsplit(rownames(heat15_myc2), "_", c('a', 'b'))[,2]), row.names = rownames(heat15_myc2))
row_annot$cluster = convert15$new.full[match(row_annot$cluster, convert15$old)]
row_annot$cluster = factor(row_annot$cluster, levels = c("1_Astro/MG", "4_GABA", "5_GABA", "6_GABA", "8_Glut", "9_Glut", "10_Glut", "11_Glut", "12_Glut", "15_GABA/Glut"))
row_annot_color = setNames(convert15$col, convert15$new.full)
pheatmap::pheatmap(heat15_myc2, clustering_callback = callback, border_color = NA, angle_col = 0, show_rownames = F, color = viridis(length(mat_breaks) - 1), breaks = mat_breaks, legend = F, annotation_col = col_annot, annotation_colors = list(cond = c(CTRL = "#007C7D", BHVE = "gold"), cluster = row_annot_color), annotation_row = row_annot, cluster_rows = F, annotation_names_row = F, annotation_names_col = T)
dev.off()

pdf("C:/Users/miles/Downloads/bvc15_sample_heatmap_logFC_viridis2.pdf")
mat_breaks <- quantile_breaks(heat15, n = 11)
col_annot = data.frame(cond = rep(c("BHVE", "CTRL"), 5), row.names = colnames(heat15))
row_annot = data.frame(cluster = factor(reshape2::colsplit(rownames(heat15_myc2), "_", c('a', 'b'))[,2]), row.names = rownames(heat15_myc2))
row_annot$cluster = convert15$new.full[match(row_annot$cluster, convert15$old)]
row_annot$cluster = factor(row_annot$cluster, levels = c("1_Astro/MG", "4_GABA", "5_GABA", "6_GABA", "8_Glut", "9_Glut", "10_Glut", "11_Glut", "12_Glut", "15_GABA/Glut"))
row_annot_color = setNames(convert15$col, convert15$new.full)
pheatmap::pheatmap(heat15_myc2, clustering_callback = callback, border_color = NA, angle_col = 0, show_rownames = F, color = viridis(length(mat_breaks) - 1), breaks = mat_breaks, legend = F, annotation_col = col_annot, annotation_colors = list(cond = c(CTRL = "#007C7D", BHVE = "gold"), cluster = row_annot_color), annotation_row = row_annot, cluster_rows = F, annotation_names_row = F, annotation_names_col = T)
dev.off()

# 15 Cluster BVC DEGs UMAP
df = as.data.frame(bb@reductions$umap@cell.embeddings)
df$cond = bb$cond
df$seuratclusters15 = bb$seuratclusters15
df$seuratclusters53 = bb$seuratclusters53
df$new15 = convert15$new.full[match(df$seuratclusters15, convert15$old)]
df$col15 = convert15$col[match(df$seuratclusters15, convert15$old)]

topDEG = bvc15_sig %>% group_by(cluster) %>% slice_min(order_by = q, n = 1)

df$gene = topDEG$mzebra[match(df$seuratclusters15, topDEG$cluster)]
df$hgnc = gene_info$human[match(df$gene, gene_info$mzebra)]
df$label = df$hgnc
df$label[which(is.na(df$label))] = df$gene[which(is.na(df$label))]

gene_coord = match(df$gene[which(!is.na(df$gene))], rownames(bb))
cell_coord = match(rownames(df)[which(!is.na(df$gene))], colnames(bb))
df$value = NA
df$value[which( !is.na(df$gene) )] = bb@assays$RNA@data[cbind(gene_coord, cell_coord)]
df$value_count = 0
df$value_count[which( !is.na(df$gene) )] = bb@assays$RNA@counts[cbind(gene_coord, cell_coord)]

# lightgray_hex = "#D3D3D3FF"
# lightgray_hex = "#000000"
lightgray_hex = "#1a1a1a"
df$col15_blend = lightgray_hex
for (i in 0:14) {
  this_df = df[which(df$seuratclusters15 == i),]
  if (! is.na(this_df$gene[1])) {
    n = nrow(this_df)
    this_values = this_df$values
    # m <- grDevices::colorRamp(c(lightgray_hex, this_df$col15[1]))( (1:n)/n )
    m <- grDevices::colorRamp(c("blue", "yellow", "red"))( (1:n)/n )
    this_df$col15_blend = lightgray_hex
    this_df$col15_blend[which(this_df$value_count != 0)] = colour_values(this_df$value[which(this_df$value_count != 0)], palette = m)
    df[rownames(this_df), "col15_blend"] = this_df$col15_blend
  }
}

my.pt.size = 0.2
p_list = list()
df = df[order(df$value, decreasing = F),]
p_list[["BHVE"]] = ggplot(df[which(df$cond == "BHVE"),], aes(UMAP_1, UMAP_2, col = col15_blend)) + geom_point(size = my.pt.size) + theme(plot.title = element_text(hjust = 0.5)) + scale_color_identity() + ggtitle("BHVE") + theme_dark() + theme(panel.background = element_rect(fill = "black", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p_list[["CTRL"]] = ggplot(df[which(df$cond == "CTRL"),], aes(UMAP_1, UMAP_2, col = col15_blend)) + geom_point(size = my.pt.size) + theme(plot.title = element_text(hjust = 0.5)) + scale_color_identity() + ggtitle("CTRL") + theme_dark() + theme(panel.background = element_rect(fill = "black", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("C:/Users/miles/Downloads/bvc15_umap_rdylbu.png", width = 2000, height = 1000, res = 120)
plot_grid(plotlist=p_list, ncol = 2)
dev.off()

#************************************************************************
# Prediction ============================================================
#************************************************************************
all_gene_df = read.csv("C:/Users/miles/Downloads/best_pred_clust0_genes.csv")
all_genes = unname(unlist(all_gene_df[,3:ncol(all_gene_df)]))
all_genes_usage = as.data.frame(table(all_genes))
all_genes_usage$all_genes = as.vector(all_genes_usage$all_genes)
all_genes_usage$pct = all_genes_usage$Freq / 19 * 100
all_genes_usage$label = gene_info$human[match(all_genes_usage$all_genes, gene_info$mzebra)]
all_genes_usage$label[which(is.na(all_genes_usage$label))] = all_genes_usage$all_genes[which(is.na(all_genes_usage$label))]
all_genes_usage = all_genes_usage[order(all_genes_usage$label),]
ggplot(all_genes_usage, aes(x = label, y = pct, color = pct, fill = pct)) + geom_bar(stat = 'identity', alpha = 0.6, position = position_dodge2()) + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + ylab("Percent of Models") + scale_color_continuous(guide = 'none') + scale_fill_continuous(guide = 'none')

# PCA
library("factoextra")
Idents(bb) = factor(bb$subsample, levels = sort(unique(bb$subsample)))
# avg_exp = myAverageExpression(bb, features = unique(all_genes))
avg_exp = myAverageExpression(bb, features = unique(all_genes_usage$all_genes[which(all_genes_usage$Freq == 19)]))
res.pca <- prcomp(t(as.matrix(avg_exp)), scale = T, center = T)
groups = factor(c(rep("Behave", 19), rep("Control", 19)))
# pdf("C:/Users/miles/Downloads/ml_subsample_counts_bulk.pdf", width = 8, height = 8)
print(fviz_pca_ind(res.pca, col.ind = groups, palette = c("#00AFBB",  "#FC4E07"), addEllipses = T, ellipse.type = "confidence", legend.title = "Groups", repel = T))
# dev.off()

# Conf
conf_df = read.csv("C:/Users/miles/Downloads/best_pred_clust0_conf.csv")
ggplot(conf_df, aes(x = cond, y = conf, color = cond, fill = cond)) + geom_boxplot(alpha = 0.5) + geom_jitter(position = position_jitter()) + scale_color_manual(values = c("#00AFBB",  "#FC4E07")) + scale_fill_manual(values = c("#00AFBB",  "#FC4E07")) + theme_bw() + ylab("Likelihood of Behavior") + xlab("")

# Clusters
clust_acc_empty = data.frame(X = 0:52, cluster = 0, mean_b = 0, mean_c = 0, mean_dif = 0, num_b_in_c = 0, num_c_in_b = 0, all_cor = 0, cluster15 = 0:52)
clust_acc = read.csv("C:/Users/miles/Downloads/scgnn_pred_cluster53.csv")
clust_acc$cluster15 = clust_acc$X
clust_acc = as.data.frame(rbindlist(list(clust_acc, data.frame("all", 0, 0, 0 ,0, 0, 0, 1, "all"))))
clust_acc = rbind(clust_acc, clust_acc_empty[which(! clust_acc_empty$X %in% clust_acc$X ),])
clust_acc$accuracy = clust_acc$all_cor
# clust_acc = xlsx::read.xlsx("C:/Users/miles/Downloads/params.xlsx", sheetName = "Final")
# clust_acc$new = convert15$new.full[match(clust_acc$cluster15, convert15$old)]
# clust_acc$col = convert15$col[match(clust_acc$cluster15, convert15$old)]
clust_acc$new = convert53$new[match(clust_acc$cluster15, convert53$old)]
clust_acc$col = convert53$col[match(clust_acc$cluster15, convert53$old)]
clust_acc$accuracy = as.numeric(as.vector(clust_acc$accuracy))
clust_acc$accuracy[which( clust_acc$accuracy == "NA" )] = 0
clust_acc[which(clust_acc$cluster15 == "all"), c("new", "col")] = c("All", "#377863")
# clust_acc$new = factor(clust_acc$new, levels = c(convert15$new.full, "All"))
clust_acc$new = factor(clust_acc$new, levels = c(convert53$new, "All"))
ggplot(clust_acc, aes(x = new, y = accuracy * 100, color = col, fill = col)) + geom_bar(stat = 'identity', alpha = 0.75) + scale_color_identity() + scale_fill_identity() + theme_bw() + theme(panel.grid.major.x = element_blank()) + scale_y_continuous(limits = c(0,100), expand = c(0, 0)) + ylab("Accuracy") + xlab("") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# bb$accuracy = clust_acc$accuracy[match(bb$seuratclusters15, clust_acc$cluster15)]
# Idents(bb) = convert15$new.full[match(bb$seuratclusters15, convert15$old)]
bb$accuracy = clust_acc$accuracy[match(bb$seuratclusters53, clust_acc$cluster15)]
Idents(bb) = convert53$new[match(bb$seuratclusters53, convert53$old)]
FeaturePlot(bb, "accuracy", order = T, label = T) + scale_color_viridis_c()

#************************************************************************
# PCRC Plot =============================================================
#************************************************************************
pcrc15  = as.data.frame(readxl::read_xlsx("C:/Users/miles/Downloads/pcrc_enrichment.xlsx", sheet = "15 Summary"))
pcrc53  = as.data.frame(readxl::read_xlsx("C:/Users/miles/Downloads/pcrc_enrichment.xlsx", sheet = "53 Summary"))
pcrcgoi = as.data.frame(readxl::read_xlsx("C:/Users/miles/Downloads/pcrc_enrichment.xlsx", sheet = "Gene Pops Summary"))
colnames(pcrcgoi)[which(colnames(pcrcgoi) == "Cohen's d")] = 'd'
pcrcgoi$cluster = tolower( pcrcgoi$human_gene )

pcrc15$num_cells = as.vector(table(bb$seuratclusters15))
pcrc53$num_cells = as.vector(table(bb$seuratclusters53))
pcrcgoi$num_cells = sapply(pcrcgoi$gene, function(x) length(which(bb@assays$RNA@counts[x,] > 0)))

pcrc15$level  = "15"
pcrc53$level  = "53"
pcrcgoi$level = "goi"

pcrc15$neg_log_p1k  = -log10(pcrc15$p_1k)
pcrc53$neg_log_p1k  = -log10(pcrc53$p_1k)
pcrcgoi$neg_log_p1k = -log10(pcrcgoi$p_1k)

pcrc15$neg_log_p10k  = -log10(pcrc15$p_10k)
pcrc53$neg_log_p10k  = -log10(pcrc53$p_10k)
pcrcgoi$neg_log_p10k = -log10(pcrcgoi$p_10k)

pcrc15$neg_log_pz  = -log10(pcrc15$p_val_adj)
pcrc53$neg_log_pz  = -log10(pcrc53$p_val_adj)
pcrcgoi$neg_log_pz = -log10(pcrcgoi$p_val_adj)

pcrc15$pz_sig  = as.numeric(pcrc15$p_val_adj) < 0.05
pcrc53$pz_sig  = as.numeric(pcrc53$p_val_adj) < 0.05
pcrcgoi$pz_sig = as.numeric(pcrcgoi$p_val_adj) < 0.05

pcrc15$p10k_sig  = as.numeric(pcrc15$p_10k) < 0.05
pcrc53$p10k_sig  = as.numeric(pcrc53$p_10k) < 0.05
pcrcgoi$p10k_sig = as.numeric(pcrcgoi$p_10k) < 0.05
pcrcgoi$p10k_sig[which(is.na(pcrcgoi$p10k_sig))] = F

pcrc15$both_sig  = pcrc15$pz_sig & pcrc15$p10k_sig
pcrc53$both_sig  = pcrc53$pz_sig & pcrc53$p10k_sig
pcrcgoi$both_sig = pcrcgoi$pz_sig & pcrcgoi$p10k_sig

pcrc15$neg_log_cells  = -log10(pcrc15$num_cells)
pcrc53$neg_log_cells  = -log10(pcrc53$num_cells)
pcrcgoi$neg_log_cells = -log10(pcrcgoi$num_cells)

pcrc15$neg_log_plot  = pcrc15$neg_log_p10k
pcrc53$neg_log_plot  = pcrc53$neg_log_p10k
pcrcgoi$neg_log_plot = pcrcgoi$neg_log_p1k

pcrcall = rbind(pcrc15[, c("d", "neg_log_plot", "level", "both_sig", "num_cells", "cluster")], pcrc53[, c("d", "neg_log_plot", "level", "both_sig", "num_cells", "cluster")], pcrcgoi[, c("d", "neg_log_plot", "level", "both_sig", "num_cells", "cluster")])
my_pal = c("#0a9396", "#ee9b00", "#ae2012")
pdf("C:/Users/miles/Downloads/pcrcall_w_label_tiny.pdf", width = 6, height = 4)
ggplot(pcrcall, aes(x = d, y = neg_log_plot, shape = both_sig, alpha = both_sig, color = level, size = num_cells, group = 1)) + geom_point(stroke = 1) + xlab("Cohen's d") + ylab(expression(-Log["10"]*" P 10k")) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60") + scale_color_manual(values = my_pal) + scale_shape_manual(values = c(1, 19)) + scale_alpha_manual(values = c(0.4, 0.8)) + theme_classic() + geom_text_repel(data = pcrcall[which(pcrcall$both_sig & pcrcall$neg_log_plot > -log10(0.05)),], aes(label = cluster))
dev.off()
# pdf("C:/Users/miles/Downloads/pcrcall_w_curve.pdf", width = 6, height = 4)
# ggplot(pcrcall, aes(x = d, y = neg_log_plot, shape = both_sig, color = level, group = 1)) + geom_point(size = 2) + xlab("Cohen's d") + ylab(expression(-Log["10"]*" P 10k")) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60") + scale_color_manual(values = my_pal) + scale_shape_manual(values = c(3, 19)) + theme_classic() + stat_smooth(method = "loess", size = 1, color = "purple") 
# dev.off()

pdf("C:/Users/miles/Downloads/pcrc15.pdf", width = 6, height = 4)
ggplot(pcrc15, aes(x = d, y = neg_log_p10k, shape = pz_sig, color = both_sig, group = 1)) + geom_point(size = 2) + ggtitle("15 Level") + xlab("Cohen's d") + ylab(expression(-Log["10"]*" P 10k")) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60") + scale_color_manual(values = c(viridis(3)[2], "gold"), guide = 'none') + scale_shape_manual(values = c(3, 19)) + theme_classic()
dev.off()
pdf("C:/Users/miles/Downloads/pcrc53.pdf", width = 6, height = 4)
ggplot(pcrc53, aes(x = d, y = neg_log_p10k, shape = pz_sig, color = both_sig, group = 1)) + geom_point(size = 1.75) + ggtitle("53 Level") + xlab("Cohen's d") + ylab(expression(-Log["10"]*" P 10k")) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60") + scale_color_manual(values = c(viridis(3)[2], "gold"), guide = 'none') + scale_shape_manual(values = c(3, 19)) + theme_classic()
dev.off()
pdf("C:/Users/miles/Downloads/pcrcgoi.pdf", width = 6, height = 4)
ggplot(pcrcgoi, aes(x = d, y = neg_log_p1k, shape = pz_sig, color = both_sig, group = 1)) + geom_point(size = 1.75) + ggtitle("GOI") + xlab("Cohen's d") + ylab(expression(-Log["10"]*" P 1k")) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60") + scale_color_manual(values = c(viridis(3)[2], "gold"), guide = 'none') + scale_shape_manual(values = c(3, 19)) + theme_classic()
dev.off()

# pdf("C:/Users/miles/Downloads/pcrc15_w_curve.pdf", width = 6, height = 4)
# ggplot(pcrc15, aes(x = d, y = neg_log_p10k, shape = pz_sig, color = both_sig, group = 1)) + geom_point(size = 2) + ggtitle("15 Level") + xlab("Cohen's d") + ylab(expression(-Log["10"]*" P 10k")) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60") + scale_color_manual(values = c(viridis(3)[2], "gold"), guide = 'none') + scale_shape_manual(values = c(3, 19)) + stat_smooth(method = "loess", size = 1, color = "purple") + theme_classic()
# dev.off()
# pdf("C:/Users/miles/Downloads/pcrc53_w_curve.pdf", width = 6, height = 4)
# ggplot(pcrc53, aes(x = d, y = neg_log_p10k, shape = pz_sig, color = both_sig, group = 1)) + geom_point(size = 1.75) + ggtitle("53 Level") + xlab("Cohen's d") + ylab(expression(-Log["10"]*" P 10k")) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60") + scale_color_manual(values = c(viridis(3)[2], "gold"), guide = 'none') + scale_shape_manual(values = c(3, 19)) + stat_smooth(method = "loess", size = 1, color = "purple") + theme_classic()
# dev.off()
# pdf("C:/Users/miles/Downloads/pcrcgoi_w_curve.pdf", width = 6, height = 4)
# ggplot(pcrcgoi, aes(x = d, y = neg_log_p1k, shape = pz_sig, color = both_sig, group = 1)) + geom_point(size = 1.75) + ggtitle("GOI") + xlab("Cohen's d") + ylab(expression(-Log["10"]*" P 1k")) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60") + scale_color_manual(values = c(viridis(3)[2], "gold"), guide = 'none') + scale_shape_manual(values = c(3, 19)) + stat_smooth(method = "loess", size = 1, color = "purple") + theme_classic()
# dev.off()

#************************************************************************
# Demux DEG =============================================================
#************************************************************************
ddeg = read.csv("C:/Users/miles/Downloads/deg_glmmseq_demux_all_clusters_cond_gsi_control_pair_subject_pool_subject_random_091821_q.csv")
ddeg[, c("new.full", "col")] = convert15[match(ddeg$cluster, convert15$old), c("new.full", "col")]
ddeg$neg_log_cond_p = -log10(ddeg$P_cond)
ddeg$neg_log_gsi_p = -log10(ddeg$P_gsi)
ddeg$neg_log_cond_q = -log10(ddeg$q_cond)
ddeg$neg_log_gsi_q = -log10(ddeg$q_gsi)
ddeg$col2 = ddeg$col
ddeg$col2[which(ddeg$q_cond >= 0.05)] = "darkgray"
ddeg$cond_is_sig = ddeg$q_cond < 0.05
ddeg$hgnc = gene_info$human[match(ddeg$X, gene_info$mzebra)]
ddeg$label = ddeg$hgnc
ddeg$label[which( is.na(ddeg$label) )] = ddeg$X[which( is.na(ddeg$label) )]

# Cond -Log10P vs GSI -Log10P
ggplot(ddeg, aes(x = neg_log_gsi_p, y = neg_log_cond_p, color = col)) + geom_point() + scale_color_identity()
ggplot(ddeg, aes(x = neg_log_gsi_p, y = neg_log_cond_p, color = col)) + geom_point() + scale_color_identity() + ylab(expression(-Log["10"]*" P of Condition")) + xlab(expression(-Log["10"]*" P of GSI")) + scale_y_sqrt() + scale_x_sqrt() + theme_light() 
pdf("C:/Users/miles/Downloads/bb_demux_clust15_deg.pdf", width = 8, height = 6)
ggplot(ddeg[which(ddeg$q_gsi > 0.05),], aes(x = neg_log_gsi_p, y = neg_log_cond_p, color = col2, alpha = cond_is_sig)) + geom_point(size = 3) + scale_color_identity() + ylab(expression(-Log["10"]*" P of Condition")) + xlab(expression(-Log["10"]*" P of GSI")) + scale_y_sqrt() + scale_x_sqrt() + scale_alpha_manual(values = c(0.3, 0.8), guide = F) + theme_light()
dev.off()

# Calculate pct and avg_logFC dif for each cluster
ddeg[,c("avg_logFC", "pct.1", "pct.2", "pct_dif", "num.1", "num.2")] = 0
pct_fc_df = data.frame()
for (this_clust in 0:14) {
  print(this_clust)
  this_res = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$seuratclusters15 == this_clust & bb$cond == "BHVE")], cells.2 = colnames(bb)[which(bb$seuratclusters15 == this_clust & bb$cond == "CTRL")])
  this_res$old = this_clust
  ddeg[which(ddeg$cluster == this_clust),c("avg_logFC", "pct.1", "pct.2", "pct_dif", "num.1", "num.2")] = this_res[match(ddeg$X[which(ddeg$cluster == this_clust)], this_res$genes), c("avg_logFC", "pct.1", "pct.2", "pct_dif", "num.1", "num.2")]
  pct_fc_df = rbind(pct_fc_df, this_res)
}

# Volcano Plot
pdf("C:/Users/miles/Downloads/bb_demux_clust15_deg_volcano.pdf", height = 6, width = 5)
ggplot(ddeg, aes(x = avg_logFC, y = neg_log_cond_p, color = col2)) + geom_point() + scale_y_sqrt() + scale_color_identity() + geom_text_repel(data=ddeg[which(ddeg$cond_is_sig & (ddeg$neg_log_cond_q > 25 | ddeg$avg_logFC > 0.26)),], aes(label =  label)) + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" P")) + theme_light()
dev.off()
ggplot(ddeg, aes(x = avg_logFC, y = neg_log_cond_p, color = col2)) + geom_point() + scale_y_log10() + scale_color_identity() + geom_text_repel(data=ddeg[which(ddeg$cond_is_sig & (ddeg$neg_log_cond_q > 25 | ddeg$avg_logFC > 0.26)),], aes(label =  label)) + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" P")) + theme_light()

# *************************************************************************************************************
# Circular Figure =============================================================================================
# *************************************************************************************************************
library(circlize)
mat1 = as.matrix(bb$nCount_RNA)
split = bb$seurat_clusters

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
scaleValues = function(values, set.min = NULL, set.max = NULL) {
  if (is.null(set.min)) { set.min = min(values) }
  if (is.null(set.max)) { set.max = max(values) }
  values_norm = (values - set.min) / (set.max - set.min)
  col_fun <- colorRamp(viridis(100))
  # col_fun = colorRamp(rev(brewer.pal(11, "RdYlBu")))
  cols <- col_fun(values_norm)
  rownames(cols) = 1:length(values)
  all.cols = data.frame(cols)
  cols = cols[which(! is.na(cols[,1])),]
  hex.cols = rgb(cols[,1], cols[,2], cols[,3], maxColorValue = 255)
  all.cols[, 4] = 0
  all.cols[rownames(cols), 4] = hex.cols
  all.cols[which(is.na(all.cols[,1])), 4] = "#FDE725FF"
  return(all.cols[,4])
}

mat = bb@assays$RNA@counts
mat[which(mat > 1)] = 1

# Calculate Differences in IEG Score, Neurgoen Score, and Cluster Proportion between BHVE and Control
bb$annot15 = factor(convert15$new.full[match(bb$seuratclusters15, convert15$old)])
bb$annot53 = factor(convert53$new[match(bb$seuratclusters53, convert53$old)])
bb$subsample = factor(bb$subsample, levels = unique(sort(bb$subsample)))
neuro_gen = read.csv("~/research/brain/data/conserved_neurogenesis_positive_zfish_mouse_cichlid.csv")[,3]
ieg_like = read.csv("~/research/brain/data/ieg_like_fos_egr1_npas4_detected_011521.csv")[,1]
bhve_evo= read.csv("~/research/brain/data/pcrc_FST20_30_LG11_evolution_genes_031821.csv")[,1]
neuro_gen = neuro_gen[which(neuro_gen %in% rownames(bb))]
neuro_gen_scores = list()
ieg_like_scores = list()
neuro_gen_scores53 = list()
neuro_gen_scores53_b = list()
neuro_gen_scores53_c = list()
ieg_like_scores53 = list()
num_cells53 = list()
print("Finding 15 Level Differences")
for (cluster in convert15$new.full) {
  # for (sample in unique(bb$sample)) {
  #   neuro_gen_scores[paste0(cluster, "_", sample)] = mean(mat[neuro_gen, which(bb$annot15 == cluster & bb$sample == sample)])
  # }
  neuro_gen_scores[[cluster]] = mean(mat[neuro_gen, which(bb$annot15 == cluster & bb$cond == "BHVE")]) - mean(mat[neuro_gen, which(bb$annot15 == cluster & bb$cond == "CTRL")])
  ieg_like_scores[[cluster]] = mean(mat[ieg_like, which(bb$annot15 == cluster & bb$cond == "BHVE")]) - mean(mat[ieg_like, which(bb$annot15 == cluster & bb$cond == "CTRL")])
}
print("Finding 53 Level Differences")
for (cluster in convert53$new) {
  neuro_gen_scores53[[cluster]] = mean(mat[neuro_gen, which(bb$annot53 == cluster & bb$cond == "BHVE")]) - mean(mat[neuro_gen, which(bb$annot53 == cluster & bb$cond == "CTRL")])
  ieg_like_scores53[[cluster]] = mean(mat[ieg_like, which(bb$annot53 == cluster & bb$cond == "BHVE")]) - mean(mat[ieg_like, which(bb$annot53 == cluster & bb$cond == "CTRL")])
  
  neuro_gen_scores53_b[[cluster]] = mean(mat[neuro_gen, which(bb$annot53 == cluster & bb$cond == "BHVE")])
  neuro_gen_scores53_c[[cluster]] = mean(mat[neuro_gen, which(bb$annot53 == cluster & bb$cond == "CTRL")])
  
  # num_cells53 = aggregate(nCount_RNA ~ sample, bb@meta.data[which(bb$annot == cluster),], length)[,2] / aggregate(nCount_RNA ~ sample, bb@meta.data, length)[,2]
  this_num_cells53 = aggregate(nCount_RNA ~ subsample, bb@meta.data[which(bb$annot53 == cluster),], length, drop = F)
  this_num_cells53[,2] = this_num_cells53[,2] / aggregate(nCount_RNA ~ subsample, bb@meta.data[which(bb$annot15 == "8_Glut"),], length)[,2]
  this_num_cells53[which( is.na(this_num_cells53[,2]) ),2] = 0
  this_num_cells53$sample = substr(this_num_cells53$subsample, 1, 2)
  this_pool = aggregate(nCount_RNA ~ sample, this_num_cells53, mean)
  this_r = cor(this_pool[,2], aggregate(depth ~ sample, bb@meta.data, median)[match(this_pool$sample, unique(bb$sample)),2])
  num_cells53[[cluster]] = this_r
}
Idents(bb) = bb$seuratclusters15
bhve_evo15_res = markerExpPerCellPerCluster(bb, bhve_evo)[[3]]
bhve_evo15_res$new_clust = convert15$new.full[match(bhve_evo15_res$cluster, convert15$old)]
Idents(bb) = bb$seuratclusters53
bhve_evo53_res = markerExpPerCellPerCluster(bb, bhve_evo)[[3]]
bhve_evo53_res$new_clust = convert53$new[match(bhve_evo53_res$cluster, convert53$old)]

# Find Overlap of Cluster DEGs for Chord Diagram
deg15 = read.csv("~/research/brain/results/bb_all_markers_15clusters_102820_more_info.csv")
deg_num = data.frame(old = 0:14, new = convert15$new.full[match(0:14, convert15$old)], num_deg = aggregate(gene ~ cluster, deg15, length)[,2])
deg_ovlp = data.frame()
dev_ovlp_big = data.frame()
for (i in 1:15) {
  i_new = convert15$new.full[i]
  i_old = convert15$old[i]
  for (j in (i+1):15) {
    j_new = convert15$new.full[j]
    j_old = convert15$old[j]
    deg_ovlp = rbind(deg_ovlp, data.frame( i_new, j_new, length(which(deg15$gene[which(deg15$cluster == i_old)] %in% deg15$gene[which(deg15$cluster == j_old)])) ))
  }
  for (j in 1:15) {
    j_new = convert15$new.full[j]
    j_old = convert15$old[j]
    dev_ovlp_big = rbind(dev_ovlp_big, data.frame( i_new, j_new, length(which(deg15$gene[which(deg15$cluster == i_old)] %in% deg15$gene[which(deg15$cluster == j_old)])) ))
  }
}
colnames(deg_ovlp) = c("cluster1", "cluster2", "ovlp")
deg_ovlp$pct = deg_ovlp$ovlp / deg_num$num_deg[match(deg_ovlp$cluster1, deg_num$new)] / deg_num$num_deg[match(deg_ovlp$cluster2, deg_num$new)]
deg_mat = acast(deg_ovlp, cluster2 ~ cluster1, value.var = "pct")
deg_mat = deg_mat[which( rownames(deg_mat) != "NA" ),]
deg_mat = deg_mat[match(convert15$new.full, rownames(deg_mat))[2:15], colnames(deg_mat)[match(convert15$new.full, colnames(deg_mat))]]

colnames(dev_ovlp_big) = c("cluster1", "cluster2", "ovlp")
dev_ovlp_big$pct = dev_ovlp_big$ovlp / deg_num$num_deg[match(dev_ovlp_big$cluster1, deg_num$new)] / deg_num$num_deg[match(dev_ovlp_big$cluster2, deg_num$new)]
deg_mat_big = acast(dev_ovlp_big, cluster2 ~ cluster1, value.var = "pct")
deg_mat_big = deg_mat_big[which( rownames(deg_mat_big) != "NA" ),]
deg_mat_big = deg_mat_big[match(convert15$new.full, rownames(deg_mat_big)), colnames(deg_mat_big)[match(convert15$new.full, colnames(deg_mat_big))]]
test = hclust(dist(deg_mat_big), method = "average")
dend = as.dendrogram(test)
circos.dendrogram(dend, facing = "inside") # Can't get this dendrogram to work, I think it needs sectors initiliazed.
# col_mat = t(as.matrix(sapply(1:nrow(deg_mat), function(x) rep(convert15$col[which( convert15$new.full == rownames(deg_mat)[x] )], ncol(deg_mat)) )))
# rownames(col_mat) = rownames(deg_mat)
# colnames(col_mat) = colnames(deg_mat)
my_col = convert15$col
names(my_col) = convert15$new.full
# chordDiagram(deg_mat, order = convert15$new.full, row.col = convert15$col[2:15], column.col = convert15$col, directional = 0)
chordDiagram(deg_mat, order = convert15$new.full, grid.col = my_col, directional = 0)

# Read in Zack's Results
neurogen15 = read.csv("~/research/brain/results/bb15_combined_new_pos88_neurogen_module_score_by_sample_070721_pqz.csv")
neurogen53 = read.csv("~/research/brain/results/bb53_combined_new_neurogen_pos88_module_score_by_sample_070121_pqz.csv")
ieg15 = read.csv("~/research/brain/results/bb15_combined_ieg_genes_half_detected_score_gg_method_051321_pqz.csv")
ieg53 = read.csv("~/research/brain/results/bb53_combined_ieg_genes_half_detected_score_gg_method_051321_pqz.csv")
ieg15$new_clust = neurogen15$new_clust = convert15$new.full[match(ieg15$cluster, convert15$old)]
ieg53$new_clust = convert53$new[match(ieg53$cluster, convert53$old)]
neurogen53$new_clust = convert53$new[match(neurogen53$cluster, convert53$old)]
ieg15$neg_log_p = -log10(ieg15$p.1)
ieg53$neg_log_p = -log10(ieg53$p.1)
neurogen15$neg_log_p = -log10(neurogen15$p.1)
neurogen53$neg_log_p = -log10(neurogen53$p.1)

# Prepare dataframes for Circos - 15 Level
bb$annot = factor(convert53$new[match(bb$seurat_clusters, convert53$old)], levels = convert53$new)
# bb_df15 = data.frame(cluster = factor(convert15$new.full), col = convert15$col, 
#                      new_sub = as.vector(table(factor(convert53$new.parent, levels = 1:15))), 
#                      neuro_gen = unlist(neuro_gen_scores), 
#                      neuro_gen_col = scaleValues(as.numeric(unlist(neuro_gen_scores))), 
#                      ieg_like = unlist(ieg_like_scores), 
#                      ieg_like_col = scaleValues(as.numeric(unlist(ieg_like_scores))) )
bb_df15 = data.frame(cluster = factor(convert15$new.full), col = convert15$col, 
                     num = reshape2::colsplit(factor(convert15$new.full), "_", c(1,2))[,1],
                     new_sub = as.vector(table(factor(convert53$new.parent, levels = 1:15))), 
                     neuro_gen = neurogen15$neg_log_p[match(convert15$new.full, neurogen15$new_clust)],
                     neuro_gen_sig = neurogen15$bh[match(convert15$new.full, neurogen15$new_clust)] < 0.05,
                     ieg_like = ieg15$neg_log_p[match(convert15$new.full, ieg15$new_clust)],
                     ieg_like_sig = ieg15$bh[match(convert15$new.full, ieg15$new_clust)] < 0.05,
                     bhve_evo = bhve_evo15_res$d[match(convert15$new.full, bhve_evo15_res$new_clust)],
                     bhve_evo_sig = bhve_evo15_res$p_val_adj[match(convert15$new.full, bhve_evo15_res$new_clust)] < 0.05,
                     bhve_evo_col = scaleValues(bhve_evo15_res$d[match(convert15$new.full, bhve_evo15_res$new_clust)]))
bb_df15$cluster = factor(bb_df15$cluster, levels = bb_df15$cluster)
bb_df15$neuro_gen[which( is.na(bb_df15$neuro_gen) )] = 0
bb_df15$neuro_gen_sig[which( is.na(bb_df15$neuro_gen_sig) )] = F
bb_df15$ieg_like[which( is.na(bb_df15$ieg_like) )] = 0
bb_df15$ieg_like_sig[which( is.na(bb_df15$ieg_like_sig) )] = F
bb_df15$neuro_gen_col = scaleValues(bb_df15$neuro_gen)
bb_df15$ieg_like_col = scaleValues(bb_df15$ieg_like)

# Prepare dataframes for Circos - 53 Level
bb_df53 = data.frame(cluster = levels(bb$annot), col = convert53$col, num_cells = aggregate(nCount_RNA ~ annot, bb@meta.data, length)[,2], 
                     neuro_gen = neurogen53$neg_log_p[match(levels(bb$annot), neurogen53$new_clust)],
                     neuro_gen_sig = neurogen53$bh[match(levels(bb$annot), neurogen53$new_clust)] < 0.05, 
                     ieg_like = ieg53$neg_log_p[match(levels(bb$annot), ieg53$new_clust)], 
                     ieg_like_sig = ieg53$bh[match(levels(bb$annot), ieg53$new_clust)] < 0.05, 
                     prop = unlist(num_cells53), neuro_gen_b = unlist(neuro_gen_scores53_b),
                     bhve_evo = bhve_evo53_res$d[match(convert53$new, bhve_evo53_res$new_clust)],
                     bhve_evo_sig = bhve_evo53_res$p_val_adj[match(convert53$new, bhve_evo53_res$new_clust)] < 0.05,
                     bhve_evo_col = scaleValues(bhve_evo53_res$d[match(convert53$new, bhve_evo53_res$new_clust)]))
bb_df53$cluster = factor(bb_df53$cluster, levels = bb_df53$cluster)
bb_df53$cluster_label = reshape2::colsplit(bb_df53$cluster, "_", c(1,2))[,1]
bb_df53$parent = reshape2::colsplit(bb_df53$cluster_label, "\\.", c(1, 2))[,1]
bb_df53$parent[which(bb_df53$cluster == "8-9_Glut")] = "8"
bb_df53$sub = reshape2::colsplit( reshape2::colsplit(bb_df53$cluster, "\\.", c(1, 2))[,2], "_", c(1,2) )[,1]
bb_df53$sub[which( is.na(bb_df53$sub) )] = "1"
bb_df53$sub[which( bb_df53$cluster == "8-9_Glut" )] = "12"
bb_df53$sub = as.numeric(bb_df53$sub)
bb_df53$neuro_gen[which( is.na(bb_df53$neuro_gen) )] = 0
bb_df53$neuro_gen_sig[which( is.na(bb_df53$neuro_gen_sig) )] = F
bb_df53$ieg_like[which( is.na(bb_df53$ieg_like) )] = 0
bb_df53$ieg_like_sig[which( is.na(bb_df53$ieg_like_sig) )] = F
bb_df53$neuro_gen_col = scaleValues(bb_df53$neuro_gen)
bb_df53$ieg_like_col = scaleValues(bb_df53$ieg_like)

# 15 Level
pdf("~/research/brain/results/sum_fig.pdf", width = 11, height = 11)
my_track_height = 0.1
# circos.par("canvas.xlim" = c(-2.4, 2.4), "canvas.ylim" = c(-2.4, 2.4))
# circos.par("canvas.xlim" = c(-3, 3), "canvas.ylim" = c(-3, 3))
circos.par("gap.after" = c(rep(0, 14), 10), cell.padding = c(0, 0, 0, 0), "start.degree" = 90, track.margin = c(0.01, 0.01), track.height = my_track_height)
circos.initialize(bb_df15$cluster, xlim = cbind(rep(0, nrow(bb_df15)), bb_df15$new_sub))
border_buffer = 0.01

# BHVE EVO
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  # 15 Level
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df15$bhve_evo_col[CELL_META$sector.numeric.index], border = NA)

  # 53 Level
  num_new_sub = bb_df15$new_sub[CELL_META$sector.numeric.index]
  for (i in 1:num_new_sub) {
    idx53 = which(bb_df53$parent == CELL_META$sector.numeric.index & bb_df53$sub == i)
    xstart = 0 + (i-1)*1
    xend = xstart + 1
    circos.rect(xstart, 0.5, xend, 1, col = bb_df53$bhve_evo_col[idx53], border = NA)
    circos.text((xstart + xend)/2, CELL_META$ycenter + 0.25, bb_df53$cluster_label[idx53], facing = "downward", niceFacing = TRUE, adj = c(0.5, 0.5), cex = 0.6)
  }
})
circos.text(-0.5, 0.5, sector.index = "1_Astro/MG", "B_EVO", track.index = 1)
for (clust15 in bb_df15$cluster[which(bb_df15$bhve_evo_sig)]) {
  xstart = 0
  xend = bb_df15$new_sub[which(bb_df15$cluster == clust15)]
  circos.rect(xstart, 0, xend, 0.5, col = NA, border = "#FDE725FF", lwd = 2, sector.index = clust15, track.index = 1)
} # Denote significance with yellow border
for (clust53 in bb_df53$cluster[which(bb_df53$bhve_evo_sig)]) {
  clust15 = bb_df15$cluster[which(bb_df15$num == bb_df53$parent[which(bb_df53$cluster == clust53)])]
  xstart = bb_df53$sub[which(bb_df53$cluster == clust53)]-1
  xend = bb_df53$sub[which(bb_df53$cluster == clust53)]
  circos.rect(xstart, 0.5, xend, 1, col = NA, border = "#FDE725FF", lwd = 2, sector.index = clust15, track.index = 1)
} # Denote significance with yellow border

# IEG
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  # 15 Level
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df15$ieg_like_col[CELL_META$sector.numeric.index], border = NA)
  xstart = CELL_META$xlim[1]
  xend = CELL_META$xlim[2]
  
  # 53 Level
  num_new_sub = bb_df15$new_sub[CELL_META$sector.numeric.index]
  for (i in 1:num_new_sub) {
    idx53 = which(bb_df53$parent == CELL_META$sector.numeric.index & bb_df53$sub == i)
    xstart = 0 + (i-1)*1
    xend = xstart + 1
    circos.rect(xstart, 0.5, xend, 1, col = bb_df53$ieg_like_col[idx53], border = NA)
    circos.text((xstart + xend)/2, CELL_META$ycenter + 0.25, bb_df53$cluster_label[idx53], facing = "downward", niceFacing = TRUE, adj = c(0.5, 0.5), cex = 0.6)
  }
})
for (clust15 in bb_df15$cluster[which(bb_df15$ieg_like_sig)]) {
  xstart = 0
  xend = bb_df15$new_sub[which(bb_df15$cluster == clust15)]
  circos.rect(xstart, 0, xend, 0.5, col = NA, border = "#FDE725FF", lwd = 2, sector.index = clust15, track.index = 2)
} # Denote significance with yellow border
for (clust53 in bb_df53$cluster[which(bb_df53$ieg_like_sig)]) {
  clust15 = bb_df15$cluster[which(bb_df15$num == bb_df53$parent[which(bb_df53$cluster == clust53)])]
  xstart = bb_df53$sub[which(bb_df53$cluster == clust53)]-1
  xend = bb_df53$sub[which(bb_df53$cluster == clust53)]
  circos.rect(xstart, 0.5, xend, 1, col = NA, border = "#FDE725FF", lwd = 2, sector.index = clust15, track.index = 2)
} # Denote significance with yellow border
circos.text(-0.5, 0.5, sector.index = "1_Astro/MG", "IEG", track.index = 2)

# Neurogenesis
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  # 15 Level
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df15$neuro_gen_col[CELL_META$sector.numeric.index], border = NA)
  xstart = CELL_META$xlim[1]
  xend = CELL_META$xlim[2]
  
  # 53 Level
  num_new_sub = bb_df15$new_sub[CELL_META$sector.numeric.index]
  for (i in 1:num_new_sub) {
    idx53 = which(bb_df53$parent == CELL_META$sector.numeric.index & bb_df53$sub == i)
    xstart = 0 + (i-1)*1
    xend = xstart + 1
    circos.rect(xstart, 0.5, xend, 1, col = bb_df53$neuro_gen_col[idx53], border = NA)
    circos.text((xstart + xend)/2, CELL_META$ycenter + 0.25,
                bb_df53$cluster_label[idx53], facing = "downward", niceFacing = TRUE,
                adj = c(0.5, 0.5), cex = 0.6)
  }
})
for (clust15 in bb_df15$cluster[which(bb_df15$neuro_gen_sig)]) {
  xstart = 0
  xend = bb_df15$new_sub[which(bb_df15$cluster == clust15)]
  circos.rect(xstart, 0, xend, 0.5, col = NA, border = "#FDE725FF", lwd = 2, sector.index = clust15, track.index = 3)
} # Denote significance with yellow border
for (clust53 in bb_df53$cluster[which(bb_df53$neuro_gen_sig)]) {
  clust15 = bb_df15$cluster[which(bb_df15$num == bb_df53$parent[which(bb_df53$cluster == clust53)])]
  xstart = bb_df53$sub[which(bb_df53$cluster == clust53)]-1
  xend = bb_df53$sub[which(bb_df53$cluster == clust53)]
  circos.rect(xstart, 0.5, xend, 1, col = NA, border = "#FDE725FF", lwd = 2, sector.index = clust15, track.index = 3)
} # Denote significance with yellow border
circos.text(-.5, 0.5, sector.index = "1_Astro/MG", "Neurogen", track.index = 3)

# Cluster Colors
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  bg.col = bb_df15$col[CELL_META$sector.numeric.index]
  
  # 15 Level
  circos.rect(CELL_META$cell.xlim[1], -0.5, CELL_META$cell.xlim[2], 0.5, col = bb_df15$col[CELL_META$sector.numeric.index], border = NA)
  circos.text(CELL_META$xcenter, 0.5 - mm_y(2),
              CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE,
              adj = c(1, 0.5), cex = 0.6, col = darken(bb_df15$col[CELL_META$sector.numeric.index], amount = 0.5))
  
  # 53 Level
  num_new_sub = bb_df15$new_sub[CELL_META$sector.numeric.index]
  for (i in 1:num_new_sub) {
    idx53 = which(bb_df53$parent == CELL_META$sector.numeric.index & bb_df53$sub == i)
    xstart = 0 + (i-1)*1
    xend = xstart + 1
    circos.rect(xstart, 0.5, xend, 1, col = bb_df53$col[idx53], border = NA)
    circos.text((xstart + xend)/2, CELL_META$ycenter + 0.25,
                bb_df53$cluster_label[idx53], facing = "downward", niceFacing = TRUE,
                adj = c(0.5, 0.5), cex = 0.6)
  }
})
dev.off()


#*******************************************************************************
# Circle Figure 2 ==============================================================
#*******************************************************************************
bb$annot15 = factor(convert15$new.full[match(bb$seuratclusters15, convert15$old)])
bb$annot53 = factor(convert53$new[match(bb$seuratclusters53, convert53$old)])
bb$subsample = factor(bb$subsample, levels = unique(sort(bb$subsample)))

zdf15 = read.csv("~/research/brain/results/bb15_all_lists_bower_all.csv")
zdf53 = read.csv("~/research/brain/results/bb53_all_lists_bower_all.csv")
bdeg15 = read.csv("~/research/brain/results/bb15_all_data_for_all_sig_cluster_genes.csv")
bdeg53 = read.csv("~/research/brain/results/bb53_all_data_for_all_sig_cluster_genes.csv")
pcrc15 = as.data.frame(readxl::read_xlsx("~/research/brain/results/pcrc_enrichment.xlsx", sheet = "15 Summary"))
pcrc53 = as.data.frame(readxl::read_xlsx("~/research/brain/results/pcrc_enrichment.xlsx", sheet = "53 Summary"))

zdf15$ieg_isSig   = zdf15$pct_sig_ieg   == 100 & zdf15$hmp_ieg   < 0.05
zdf15$neuro_isSig = zdf15$pct_sig_neuro == 100 & zdf15$hmp_neuro < 0.05
zdf15$pcrc_isSig  = zdf15$pct_sig_pcrc  == 100 & zdf15$hmp_pcrc  < 0.05

zdf53$ieg_isSig   = zdf53$pct_sig_ieg   == 100 & zdf53$hmp_ieg   < 0.05
zdf53$neuro_isSig = zdf53$pct_sig_neuro == 100 & zdf53$hmp_neuro < 0.05
zdf53$pcrc_isSig  = zdf53$pct_sig_pcrc  == 100 & zdf53$hmp_pcrc  < 0.05

pcrc15$pcrc_enrich_isSig = pcrc15$p_val_adj < 0.05 & pcrc15$p_10k < 0.05
pcrc53$pcrc_enrich_isSig = pcrc53$p_val_adj < 0.05 & pcrc53$p_10k < 0.05

pcrc15$neg_log_p1k = -log10(pcrc15$p_1k)
pcrc15$neg_log_p10k = -log10(pcrc15$p_10k)
pcrc53$neg_log_p10k = -log10(pcrc53$p_10k)

bdeg15_num = aggregate(formula = X ~ cluster.1, data = bdeg15, FUN = length)
bdeg53_num = aggregate(formula = X ~ cluster.1, data = bdeg53, FUN = length)

num_cells_15 = aggregate(nCount_RNA ~ seuratclusters15, bb@meta.data, length)

bb_df15 = convert15
bb_df15$num = as.character(bb_df15$new.num)
bb_df15$cluster = bb_df15$new.full
bb_df15$num_cells = num_cells_15$nCount_RNA[match(bb_df15$old, num_cells_15$seuratclusters15)]
bb_df15$new_sub = as.vector(table(factor(convert53$new.parent, levels = 1:15)))
bb_df15$bdeg_num = bdeg15_num$X[match(bb_df15$old, bdeg15_num$cluster.1)]
bb_df15$zieg   = zdf15$neg_log_p_max_ieg[match(bb_df15$old,      zdf15$cluster)]
bb_df15$zneuro = zdf15$neg_log_p_max_neurogen[match(bb_df15$old, zdf15$cluster)]
bb_df15$zpcrc  = zdf15$neg_log_p_max_pcrc[match(bb_df15$old,   zdf15$cluster)]
bb_df15[is.na(bb_df15)] = 0
bb_df15$zieg_isSig   = zdf15$ieg_isSig[match(bb_df15$old,   zdf15$cluster)]
bb_df15$zneuro_isSig = zdf15$neuro_isSig[match(bb_df15$old, zdf15$cluster)]
bb_df15$zpcrc_isSig  = zdf15$pcrc_isSig[match(bb_df15$old,  zdf15$cluster)]
bb_df15$pcrc_enrich = pcrc15$neg_log_p10k[match(bb_df15$old, pcrc15$cluster)]
bb_df15$pcrc_enrich_isSig = pcrc15$pcrc_enrich_isSig[match(bb_df15$old, pcrc15$cluster)]
bb_df15[is.na(bb_df15)] = FALSE
bb_df15$bdeg_num_col = scaleValues(bb_df15$bdeg_num)
bb_df15$zieg_col     = scaleValues(bb_df15$zieg)
bb_df15$zneuro_col   = scaleValues(bb_df15$zneuro)
bb_df15$zpcrc_col    = scaleValues(bb_df15$zpcrc)
bb_df15$pcrc_enrich_col = scaleValues(bb_df15$pcrc_enrich)

bb_df53 = convert53
bb_df53$cluster = bb_df53$new
bb_df53$cluster_label = reshape2::colsplit(bb_df53$cluster, "_", c(1,2))[,1]
bb_df53$parent = reshape2::colsplit(bb_df53$cluster_label, "\\.", c(1, 2))[,1]
bb_df53$parent[which(bb_df53$cluster == "8-9_Glut")] = "8"
bb_df53$sub = reshape2::colsplit( reshape2::colsplit(bb_df53$cluster, "\\.", c(1, 2))[,2], "_", c(1,2) )[,1]
bb_df53$sub[which( is.na(bb_df53$sub) )] = "1"
bb_df53$sub[which( bb_df53$cluster == "8-9_Glut" )] = "12"
bb_df53$sub = as.numeric(bb_df53$sub)
bb_df53$bdeg_num = bdeg53_num$X[match(bb_df53$old, bdeg53_num$cluster.1)]
bb_df53$zieg   = zdf53$neg_log_p_max_ieg[match(bb_df53$old,      zdf53$cluster)]
bb_df53$zneuro = zdf53$neg_log_p_max_neurogen[match(bb_df53$old, zdf53$cluster)]
bb_df53$zpcrc  = zdf53$neg_log_p_max_pcrc[match(bb_df53$old,   zdf53$cluster)]
bb_df53[is.na(bb_df53)] = 0
bb_df53$pcrc_enrich = pcrc53$neg_log_p10k[match(bb_df53$old, pcrc53$cluster)]
bb_df53$pcrc_enrich_isSig = pcrc53$pcrc_enrich_isSig[match(bb_df53$old, pcrc53$cluster)]
bb_df53$zieg_isSig   = zdf53$ieg_isSig[match(bb_df53$old,   zdf53$cluster)]
bb_df53$zneuro_isSig = zdf53$neuro_isSig[match(bb_df53$old, zdf53$cluster)]
bb_df53$zpcrc_isSig  = zdf53$pcrc_isSig[match(bb_df53$old,  zdf53$cluster)]
bb_df53[is.na(bb_df53)] = FALSE
bb_df53$bdeg_num_col    = scaleValues(bb_df53$bdeg_num)
bb_df53$zieg_col        = scaleValues(bb_df53$zieg)
bb_df53$zneuro_col      = scaleValues(bb_df53$zneuro)
bb_df53$zpcrc_col       = scaleValues(bb_df53$zpcrc)
bb_df53$pcrc_enrich_col = scaleValues(bb_df53$pcrc_enrich)


pdf("~/research/brain/results/sum_fig2.pdf", width = 11, height = 11)
my_track_height = 0.1
circos.par("gap.after" = c(rep(0, 14), 10), cell.padding = c(0, 0, 0, 0), "start.degree" = 90, track.margin = c(0.01, 0.01), track.height = my_track_height)
circos.initialize(bb_df15$cluster, xlim = cbind(rep(0, nrow(bb_df15)), bb_df15$new_sub))
border_buffer = 0.01
sigBorderThickness = 3
sigBorderCol = "goldenrod"
this_track_index = 0

# BDEG Num
this_track_index = this_track_index + 1
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  # 15 Level
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df15$bdeg_num_col[CELL_META$sector.numeric.index], border = NA)
  
  # 53 Level
  num_new_sub = bb_df15$new_sub[CELL_META$sector.numeric.index]
  for (i in 1:num_new_sub) {
    idx53 = which(bb_df53$parent == CELL_META$sector.numeric.index & bb_df53$sub == i)
    xstart = 0 + (i-1)*1
    xend = xstart + 1
    circos.rect(xstart, 0.5, xend, 1, col = bb_df53$bdeg_num_col[idx53], border = NA)
    circos.text((xstart + xend)/2, CELL_META$ycenter + 0.25, bb_df53$cluster_label[idx53], facing = "downward", niceFacing = TRUE, adj = c(0.5, 0.5), cex = 0.6)
  }
})
circos.text(-0.5, 0.5, sector.index = "1_Astro/MG", "# DEG", track.index = this_track_index)

# BHVE EVO Enrichment
this_track_index = this_track_index + 1
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  # 15 Level
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df15$pcrc_enrich_col[CELL_META$sector.numeric.index], border = NA)
  
  # 53 Level
  num_new_sub = bb_df15$new_sub[CELL_META$sector.numeric.index]
  for (i in 1:num_new_sub) {
    idx53 = which(bb_df53$parent == CELL_META$sector.numeric.index & bb_df53$sub == i)
    xstart = 0 + (i-1)*1
    xend = xstart + 1
    circos.rect(xstart, 0.5, xend, 1, col = bb_df53$pcrc_enrich_col[idx53], border = NA)
    circos.text((xstart + xend)/2, CELL_META$ycenter + 0.25, bb_df53$cluster_label[idx53], facing = "downward", niceFacing = TRUE, adj = c(0.5, 0.5), cex = 0.6)
  }
})
circos.text(-0.5, 0.5, sector.index = "1_Astro/MG", "BEVO 10k", track.index = this_track_index)
for (clust15 in bb_df15$cluster[which(bb_df15$pcrc_enrich_isSig)]) {
  xstart = 0
  xend = bb_df15$new_sub[which(bb_df15$cluster == clust15)]
  circos.rect(xstart, 0, xend, 0.5, col = NA, border = sigBorderCol, lwd = sigBorderThickness, sector.index = clust15, track.index = this_track_index)
} # Denote significance with yellow border
for (clust53 in bb_df53$cluster[which(bb_df53$pcrc_enrich_isSig)]) {
  clust15 = bb_df15$cluster[which(bb_df15$num == bb_df53$parent[which(bb_df53$cluster == clust53)])]
  xstart = bb_df53$sub[which(bb_df53$cluster == clust53)]-1
  xend = bb_df53$sub[which(bb_df53$cluster == clust53)]
  circos.rect(xstart, 0.5, xend, 1, col = NA, border = sigBorderCol, lwd = sigBorderThickness, sector.index = clust15, track.index = this_track_index)
} # Denote significance with yellow border


# Z BHVE EVO 
this_track_index = this_track_index + 1
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  # 15 Level
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df15$zpcrc_col[CELL_META$sector.numeric.index], border = NA)
  
  # 53 Level
  num_new_sub = bb_df15$new_sub[CELL_META$sector.numeric.index]
  for (i in 1:num_new_sub) {
    idx53 = which(bb_df53$parent == CELL_META$sector.numeric.index & bb_df53$sub == i)
    xstart = 0 + (i-1)*1
    xend = xstart + 1
    circos.rect(xstart, 0.5, xend, 1, col = bb_df53$zpcrc_col[idx53], border = NA)
    circos.text((xstart + xend)/2, CELL_META$ycenter + 0.25, bb_df53$cluster_label[idx53], facing = "downward", niceFacing = TRUE, adj = c(0.5, 0.5), cex = 0.6)
  }
})
circos.text(-0.5, 0.5, sector.index = "1_Astro/MG", "BEVO", track.index = this_track_index)
for (clust15 in bb_df15$cluster[which(bb_df15$zpcrc_isSig)]) {
  xstart = 0
  xend = bb_df15$new_sub[which(bb_df15$cluster == clust15)]
  circos.rect(xstart, 0, xend, 0.5, col = NA, border = sigBorderCol, lwd = sigBorderThickness, sector.index = clust15, track.index = this_track_index)
} # Denote significance with yellow border
for (clust53 in bb_df53$cluster[which(bb_df53$zpcrc_isSig)]) {
  clust15 = bb_df15$cluster[which(bb_df15$num == bb_df53$parent[which(bb_df53$cluster == clust53)])]
  xstart = bb_df53$sub[which(bb_df53$cluster == clust53)]-1
  xend = bb_df53$sub[which(bb_df53$cluster == clust53)]
  circos.rect(xstart, 0.5, xend, 1, col = NA, border = sigBorderCol, lwd = sigBorderThickness, sector.index = clust15, track.index = this_track_index)
} # Denote significance with yellow border

# IEG
this_track_index = this_track_index + 1
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  # 15 Level
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df15$zieg_col[CELL_META$sector.numeric.index], border = NA)
  xstart = CELL_META$xlim[1]
  xend = CELL_META$xlim[2]
  
  # 53 Level
  num_new_sub = bb_df15$new_sub[CELL_META$sector.numeric.index]
  for (i in 1:num_new_sub) {
    idx53 = which(bb_df53$parent == CELL_META$sector.numeric.index & bb_df53$sub == i)
    xstart = 0 + (i-1)*1
    xend = xstart + 1
    circos.rect(xstart, 0.5, xend, 1, col = bb_df53$zieg_col[idx53], border = NA)
    circos.text((xstart + xend)/2, CELL_META$ycenter + 0.25, bb_df53$cluster_label[idx53], facing = "downward", niceFacing = TRUE, adj = c(0.5, 0.5), cex = 0.6)
  }
})
for (clust15 in bb_df15$cluster[which(bb_df15$zieg_isSig)]) {
  xstart = 0
  xend = bb_df15$new_sub[which(bb_df15$cluster == clust15)]
  circos.rect(xstart, 0, xend, 0.5, col = NA, border = sigBorderCol, lwd = sigBorderThickness, sector.index = clust15, track.index = this_track_index)
} # Denote significance with yellow border
for (clust53 in bb_df53$cluster[which(bb_df53$zieg_isSig)]) {
  clust15 = bb_df15$cluster[which(bb_df15$num == bb_df53$parent[which(bb_df53$cluster == clust53)])]
  xstart = bb_df53$sub[which(bb_df53$cluster == clust53)]-1
  xend = bb_df53$sub[which(bb_df53$cluster == clust53)]
  circos.rect(xstart, 0.5, xend, 1, col = NA, border = sigBorderCol, lwd = sigBorderThickness, sector.index = clust15, track.index = this_track_index)
} # Denote significance with yellow border
circos.text(-0.5, 0.5, sector.index = "1_Astro/MG", "IEG", track.index = this_track_index)

# Neurogenesis
this_track_index = this_track_index + 1
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  # 15 Level
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df15$zneuro_col[CELL_META$sector.numeric.index], border = NA)
  xstart = CELL_META$xlim[1]
  xend = CELL_META$xlim[2]
  
  # 53 Level
  num_new_sub = bb_df15$new_sub[CELL_META$sector.numeric.index]
  for (i in 1:num_new_sub) {
    idx53 = which(bb_df53$parent == CELL_META$sector.numeric.index & bb_df53$sub == i)
    xstart = 0 + (i-1)*1
    xend = xstart + 1
    circos.rect(xstart, 0.5, xend, 1, col = bb_df53$zneuro_col[idx53], border = NA)
    circos.text((xstart + xend)/2, CELL_META$ycenter + 0.25,
                bb_df53$cluster_label[idx53], facing = "downward", niceFacing = TRUE,
                adj = c(0.5, 0.5), cex = 0.6)
  }
})
for (clust15 in bb_df15$cluster[which(bb_df15$zneuro_isSig)]) {
  xstart = 0
  xend = bb_df15$new_sub[which(bb_df15$cluster == clust15)]
  circos.rect(xstart, 0, xend, 0.5, col = NA, border = sigBorderCol, lwd = sigBorderThickness, sector.index = clust15, track.index = this_track_index)
} # Denote significance with yellow border
for (clust53 in bb_df53$cluster[which(bb_df53$zneuro_isSig)]) {
  clust15 = bb_df15$cluster[which(bb_df15$num == bb_df53$parent[which(bb_df53$cluster == clust53)])]
  xstart = bb_df53$sub[which(bb_df53$cluster == clust53)]-1
  xend = bb_df53$sub[which(bb_df53$cluster == clust53)]
  circos.rect(xstart, 0.5, xend, 1, col = NA, border = sigBorderCol, lwd = sigBorderThickness, sector.index = clust15, track.index = this_track_index)
} # Denote significance with yellow border
circos.text(-.5, 0.5, sector.index = "1_Astro/MG", "NG", track.index = this_track_index)

# Cluster Colors
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  bg.col = bb_df15$col[CELL_META$sector.numeric.index]
  
  # 15 Level
  circos.rect(CELL_META$cell.xlim[1], -0.5, CELL_META$cell.xlim[2], 0.5, col = bb_df15$col[CELL_META$sector.numeric.index], border = NA)
  this_facing = "clockwise"
  this_adj = c(1, 0.5)
  this_y = 0.5 - mm_y(2)
  if (CELL_META$sector.index %in% c("15_GABA/Glut", "1_Astro/MG", "8_Glut", "9_Glut")) { this_facing = "downward"; this_adj = c(0.5, 0.5); this_y = CELL_META$ycenter - 0.5; }
  circos.text(CELL_META$xcenter, this_y,
              CELL_META$sector.index, facing = this_facing, niceFacing = TRUE,
              adj = this_adj, cex = 0.6, col = darken(bb_df15$col[CELL_META$sector.numeric.index], amount = 0.5))
  
  # 53 Level
  num_new_sub = bb_df15$new_sub[CELL_META$sector.numeric.index]
  for (i in 1:num_new_sub) {
    idx53 = which(bb_df53$parent == CELL_META$sector.numeric.index & bb_df53$sub == i)
    xstart = 0 + (i-1)*1
    xend = xstart + 1
    circos.rect(xstart, 0.5, xend, 1, col = bb_df53$col[idx53], border = NA)
    circos.text((xstart + xend)/2, CELL_META$ycenter + 0.25,
                bb_df53$cluster_label[idx53], facing = "downward", niceFacing = TRUE,
                adj = c(0.5, 0.5), cex = 0.6)
  }
})
highlight.sector("2_OPC/Oligo", col = "#FF000020", border = "red", lwd = 1, padding = c(-0.073, 0, 0.15, -0.5))
highlight.sector("4_GABA", col = "#FF000020", border = "red", lwd = 1, padding = c(-0.073, 0.086, 0.15, 0.065))
highlight.sector("9_Glut", col = "#FF000020", border = "red", lwd = 1, padding = c(-0.073, -0.625, 0.15, -0.25))
highlight.sector("10_Glut", col = "#FF000020", border = "red", lwd = 1, padding = c(-0.073, 0, 0.15, -0.5))
dev.off()

#*******************************************************************************
# Circle Figure 3 ==============================================================
#*******************************************************************************
bb$annot15 = factor(convert15$new.full[match(bb$seuratclusters15, convert15$old)])
bb$annot53 = factor(convert53$new[match(bb$seuratclusters53, convert53$old)])
bb$subsample = factor(bb$subsample, levels = unique(sort(bb$subsample)))

# Generate random bDEG results
bdeg.res = read.csv("~/Downloads/degs_by_primary_seecondary_cluster_for_dot_heatmap_122921.csv")
bdeg.res$cluster_new = str_replace(bdeg.res$cluster_new, "Astro", "RGC")
bdeg.res$old_new = paste0(bdeg.res$cluster, " - ", bdeg.res$cluster_new)
bdeg15 = bdeg.res[which(bdeg.res$cat == "bower" & bdeg.res$old_new %in% paste0(convert15$old, " - ", convert15$new.full)),]
bdeg53 = bdeg.res[which(bdeg.res$cat == "bower" & bdeg.res$old_new %in% paste0(convert53$old, " - ", convert53$new)),]

# Generate random Neurogen results
ng = xlsx::read.xlsx("~/Downloads/out_neurogen_summary_for_dot_plotting_012722_cond_not_BAI.xlsx", sheetIndex = 1)
ng = ng[which(ng$cat == "bower"),]
ng$cluster_new = str_replace(ng$cluster_new, "Astro", "RGC")
ng$ng = ng$prop_sig
ng15 = ng[which(ng$level == "primary"),]
ng53 = ng[which(ng$level == "secondary"),]

# Load PCRC Results
pcrc15 = xlsx::read.xlsx("~/Downloads/pcrc2020_enrichment.xlsx", sheetName = "15 Summary")
pcrc53 = xlsx::read.xlsx("~/Downloads/pcrc2020_enrichment.xlsx", sheetName = "53 Summary")
pcrc15$neg.prop = pcrc15$neg_perm / 10000
pcrc53$neg.prop = pcrc53$neg_10k / 10000
pcrc15$isSig = pcrc15$p_perm < 0.05
pcrc53$isSig = pcrc53$p_10k < 0.05

# Load IEG Results
ieg.res = read.csv("~/Downloads/ieg_summary_for_dotplot_heatmap.csv")
ieg.res$num.model = as.numeric(plyr::revalue(ieg.res$cat, replace = c("bower" = 6, "gsi" = 5, "quiver" = 5)))
ieg.res$prop.sig = ieg.res$n_sig / ieg.res$num.model
ieg.res$cluster_new = str_replace(ieg.res$cluster_new, "Astro", "RGC")
ieg15 = ieg.res[which(ieg.res$level == "primary"),]
ieg53 = ieg.res[which(ieg.res$level == "secondary"),]

# Find # of cells per cluster
num_cells_15 = aggregate(nCount_RNA ~ seuratclusters15, bb@meta.data, length)
num_cells_53 = aggregate(nCount_RNA ~ seuratclusters53, bb@meta.data, length)

bb_df15 = convert15
bb_df15$num = as.character(bb_df15$new.num)
bb_df15$cluster = bb_df15$new.full
bb_df15$num_cells = num_cells_15$nCount_RNA[match(bb_df15$old, num_cells_15$seuratclusters15)]
bb_df15$new_sub = as.vector(table(factor(convert53$new.parent, levels = 1:15)))
bb_df15$new.junk = str_replace(bb_df15$new.junk, "Astro", "RGC")
bb_df15$ieg = ieg15$prop.sig[match(bb_df15$old, ieg15$cluster_old)]
bb_df15$num.deg = bdeg15$freq[match(bb_df15$old, bdeg15$cluster)]
bb_df15$ng = ng15$ng[match(bb_df15$old, ng15$cluster)]
bb_df15$pcrc = pcrc15$neg.prop[match(bb_df15$old, pcrc15$cluster)]
bb_df15[is.na(bb_df15)] = 0
bb_df15$pcrc_isSig  = pcrc15$isSig[match(bb_df15$old, pcrc15$cluster)]
bb_df15$num.deg.col = scaleValues(bb_df15$num.deg, set.min = 0, set.max = 60)
bb_df15$ieg.col     = scaleValues(bb_df15$ieg)
bb_df15$ng.col   = scaleValues(bb_df15$ng)
bb_df15$pcrc.col    = scaleValues(bb_df15$pcrc)

bb_df53 = convert53
bb_df53$cluster = bb_df53$new
bb_df53$cluster_label = reshape2::colsplit(bb_df53$cluster, "_", c(1,2))[,1]
bb_df53$parent = reshape2::colsplit(bb_df53$cluster_label, "\\.", c(1, 2))[,1]
bb_df53$parent[which(bb_df53$cluster == "8-9_Glut")] = "8"
bb_df53$sub = reshape2::colsplit( reshape2::colsplit(bb_df53$cluster, "\\.", c(1, 2))[,2], "_", c(1,2) )[,1]
bb_df53$sub[which( is.na(bb_df53$sub) )] = "1"
bb_df53$sub[which( bb_df53$cluster == "8-9_Glut" )] = "12"
bb_df53$sub = as.numeric(bb_df53$sub)
bb_df53$ieg = ieg53$prop.sig[match(bb_df53$old, ieg53$cluster_old)]
bb_df53$num.deg = bdeg53$freq[match(bb_df53$old, bdeg53$cluster)]
bb_df53$ng = ng53$ng[match(bb_df53$old, ng53$cluster)]
bb_df53$pcrc = pcrc53$neg.prop[match(bb_df53$old, pcrc53$cluster)]
bb_df53[is.na(bb_df53)] = 0
bb_df53$pcrc_isSig  = pcrc53$isSig[match(bb_df53$old, pcrc53$cluster)]
bb_df53$num.deg.col = scaleValues(bb_df53$num.deg, set.min = 0, set.max = 60)
bb_df53$ieg.col     = scaleValues(bb_df53$ieg)
bb_df53$ng.col   = scaleValues(bb_df53$ng)
bb_df53$pcrc.col    = scaleValues(bb_df53$pcrc, set.min = 0, set.max = 9950)

library(circlize)
pdf("~/research/brain/results/sum_fig3.pdf", width = 11, height = 11)
my_track_height = 0.1
circos.par("gap.after" = c(rep(0, 14), 10), cell.padding = c(0, 0, 0, 0), "start.degree" = 90, track.margin = c(0.01, 0.01), track.height = my_track_height)
circos.initialize(bb_df15$cluster, xlim = cbind(rep(0, nrow(bb_df15)), bb_df15$new_sub))
border_buffer = 0.01
sigBorderThickness = 3
sigBorderCol = "goldenrod"
this_track_index = 0

# BDEG Num
this_track_index = this_track_index + 1
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  # 15 Level
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df15$num.deg.col[CELL_META$sector.numeric.index], border = NA)
  
  # 53 Level
  num_new_sub = bb_df15$new_sub[CELL_META$sector.numeric.index]
  for (i in 1:num_new_sub) {
    idx53 = which(bb_df53$parent == CELL_META$sector.numeric.index & bb_df53$sub == i)
    xstart = 0 + (i-1)*1
    xend = xstart + 1
    circos.rect(xstart, 0.5, xend, 1, col = bb_df53$num.deg.col[idx53], border = NA)
    circos.text((xstart + xend)/2, CELL_META$ycenter + 0.25, bb_df53$cluster_label[idx53], facing = "downward", niceFacing = TRUE, adj = c(0.5, 0.5), cex = 0.6)
  }
})
circos.text(-0.5, 0.5, sector.index = "1_RGC/MG", "# DEG", track.index = this_track_index)

# IEG
this_track_index = this_track_index + 1
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  # 15 Level
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df15$ieg.col[CELL_META$sector.numeric.index], border = NA)
  xstart = CELL_META$xlim[1]
  xend = CELL_META$xlim[2]
  
  # 53 Level
  num_new_sub = bb_df15$new_sub[CELL_META$sector.numeric.index]
  for (i in 1:num_new_sub) {
    idx53 = which(bb_df53$parent == CELL_META$sector.numeric.index & bb_df53$sub == i)
    xstart = 0 + (i-1)*1
    xend = xstart + 1
    circos.rect(xstart, 0.5, xend, 1, col = bb_df53$ieg.col[idx53], border = NA)
    circos.text((xstart + xend)/2, CELL_META$ycenter + 0.25, bb_df53$cluster_label[idx53], facing = "downward", niceFacing = TRUE, adj = c(0.5, 0.5), cex = 0.6)
  }
})
circos.text(-0.5, 0.5, sector.index = "1_RGC/MG", "IEG", track.index = this_track_index)

# Neurogenesis
this_track_index = this_track_index + 1
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  # 15 Level
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df15$ng.col[CELL_META$sector.numeric.index], border = NA)
  xstart = CELL_META$xlim[1]
  xend = CELL_META$xlim[2]
  
  # 53 Level
  num_new_sub = bb_df15$new_sub[CELL_META$sector.numeric.index]
  for (i in 1:num_new_sub) {
    idx53 = which(bb_df53$parent == CELL_META$sector.numeric.index & bb_df53$sub == i)
    xstart = 0 + (i-1)*1
    xend = xstart + 1
    circos.rect(xstart, 0.5, xend, 1, col = bb_df53$ng.col[idx53], border = NA)
    circos.text((xstart + xend)/2, CELL_META$ycenter + 0.25,
                bb_df53$cluster_label[idx53], facing = "downward", niceFacing = TRUE,
                adj = c(0.5, 0.5), cex = 0.6)
  }
})
circos.text(-.5, 0.5, sector.index = "1_RGC/MG", "NG", track.index = this_track_index)

# BHVE EVO Enrichment
this_track_index = this_track_index + 1
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  # 15 Level
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df15$pcrc.col[CELL_META$sector.numeric.index], border = NA)
  
  # 53 Level
  num_new_sub = bb_df15$new_sub[CELL_META$sector.numeric.index]
  for (i in 1:num_new_sub) {
    idx53 = which(bb_df53$parent == CELL_META$sector.numeric.index & bb_df53$sub == i)
    xstart = 0 + (i-1)*1
    xend = xstart + 1
    circos.rect(xstart, 0.5, xend, 1, col = bb_df53$pcrc.col[idx53], border = NA)
    circos.text((xstart + xend)/2, CELL_META$ycenter + 0.25, bb_df53$cluster_label[idx53], facing = "downward", niceFacing = TRUE, adj = c(0.5, 0.5), cex = 0.6)
  }
})
circos.text(-0.5, 0.5, sector.index = "1_RGC/MG", "PCRC", track.index = this_track_index)
# for (clust15 in bb_df15$cluster[which(bb_df15$pcrc_isSig)]) {
#   xstart = 0
#   xend = bb_df15$new_sub[which(bb_df15$cluster == clust15)]
#   circos.rect(xstart, 0, xend, 0.5, col = NA, border = sigBorderCol, lwd = sigBorderThickness, sector.index = clust15, track.index = this_track_index)
# } # Denote significance with yellow border
# for (clust53 in bb_df53$cluster[which(bb_df53$pcrc_isSig)]) {
#   clust15 = bb_df15$cluster[which(bb_df15$num == bb_df53$parent[which(bb_df53$cluster == clust53)])]
#   xstart = bb_df53$sub[which(bb_df53$cluster == clust53)]-1
#   xend = bb_df53$sub[which(bb_df53$cluster == clust53)]
#   circos.rect(xstart, 0.5, xend, 1, col = NA, border = sigBorderCol, lwd = sigBorderThickness, sector.index = clust15, track.index = this_track_index)
# } # Denote significance with yellow border

# Cluster Colors
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  bg.col = bb_df15$col[CELL_META$sector.numeric.index]
  
  # 15 Level
  circos.rect(CELL_META$cell.xlim[1], -0.5, CELL_META$cell.xlim[2], 0.5, col = bb_df15$col[CELL_META$sector.numeric.index], border = NA)
  this_facing = "clockwise"
  this_adj = c(1, 0.5)
  this_y = 0.5 - mm_y(2)
  this.cex = 0.8
  if (CELL_META$sector.index %in% c("2_OPC/Oligo")) { this.cex = 0.6 }
  if (CELL_META$sector.index %in% c("15_GABA/Glut", "1_RGC/MG", "8_Glut", "9_Glut")) { this_facing = "downward"; this_adj = c(0.5, 0.5); this_y = CELL_META$ycenter - 0.5; }
  circos.text(CELL_META$xcenter, 0,
              bb_df15$new.junk[CELL_META$sector.numeric.index], facing = this_facing, niceFacing = TRUE,
              adj = c(0.5, 0.5), cex = this.cex, col = darken(bb_df15$col[CELL_META$sector.numeric.index], amount = 0.5))
  circos.text(CELL_META$xcenter, -0.8,
              bb_df15$new.num[CELL_META$sector.numeric.index], facing = "downward", niceFacing = TRUE,
              adj = c(0.5, 0.5), cex = 1, col = darken(bb_df15$col[CELL_META$sector.numeric.index], amount = 0.3))
  
  # 53 Level
  num_new_sub = bb_df15$new_sub[CELL_META$sector.numeric.index]
  for (i in 1:num_new_sub) {
    idx53 = which(bb_df53$parent == CELL_META$sector.numeric.index & bb_df53$sub == i)
    xstart = 0 + (i-1)*1
    xend = xstart + 1
    circos.rect(xstart, 0.5, xend, 1, col = bb_df53$col[idx53], border = NA)
    circos.text((xstart + xend)/2, CELL_META$ycenter + 0.25,
                bb_df53$cluster_label[idx53], facing = "downward", niceFacing = TRUE,
                adj = c(0.5, 0.5), cex = 0.6)
  }
})
dev.off()

# Old code:
# track1_breaks = pretty(unlist(neuro_gen_scores), n = 2)
# circos.track(ylim = c(min(track1_breaks), max(track1_breaks)), track.height = 0.30, bg.border = "black", panel.fun = function(x, y) {
#   for (tb in track1_breaks) {
#     if (tb != track1_breaks[1] & tb != track1_breaks[length(track1_breaks)])
#       circos.segments(CELL_META$cell.xlim[1], tb, CELL_META$cell.xlim[2], tb, col = "gray60", lty = 2)
#   }
#   
#   xrange = CELL_META$cell.xlim[2] - CELL_META$cell.xlim[1]
#   this_gap = 0.15*xrange
#   this_gap = ifelse(this_gap > 1.5, 1.5, this_gap)
#   for (sample in unique(bb$sample)) {
#     if ( startsWith(sample, "b") ) {
#       this_col = "#d0000090"
#       this_x = myJitter(CELL_META$xcenter-this_gap, 0.1*xrange)
#     } else {
#       this_col = "#023e8a90"
#       this_x = myJitter(CELL_META$xcenter+this_gap, 0.1*xrange)
#     }
#     print(sample)
#     this_score = as.numeric(neuro_gen_scores[paste0(CELL_META$sector.index, "_", sample)])
#     print(this_score)
#     circos.points(this_x, this_score, col = this_col, pch = 20, cex = 1.5)
#   }
# })
# circos.yaxis(at = track1_breaks, sector.index = "1_Astro/MG", track.index = 1, side = "left")
# # 53 Level
# par(new = TRUE)
# circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1, 1))
# # circos.par("canvas.xlim" = c(-0.75, 0.75), "canvas.ylim" = c(-0.75, 0.75))
# circos.par("gap.after" = c(rep(0, 52), 10), cell.padding = c(0, 0, 0, 0), "start.degree" = 90, track.margin = c(0.02, 0.02))
# circos.initialize(bb_df53$cluster, x = bb_df53$cluster, matrix(rep(c(0,1), 53), nrow = 53, byrow = T))
# # Cluster Proportion
# track1_breaks = pretty(bb_df53$prop, n = 2)
# # track1_breaks[1] = -0.25
# circos.track(ylim = c(min(track1_breaks), max(track1_breaks)), track.height = 0.15, bg.border = NA, panel.fun = function(x, y) {
#   for (tb in track1_breaks) {
#     # if (tb != track1_breaks[1] & tb != track1_breaks[length(track1_breaks)])
#     if (tb != track1_breaks[1] & tb != track1_breaks[length(track1_breaks)])
#       circos.segments(CELL_META$cell.xlim[1], tb, CELL_META$cell.xlim[2], tb, col = "gray60", lty = 2)
#     else
#       circos.segments(CELL_META$cell.xlim[1], tb, CELL_META$cell.xlim[2], tb, col = "black", lty = 1)
#   }
#   if (CELL_META$sector.numeric.index == nrow(bb_df53))
#     circos.segments(CELL_META$cell.xlim[2], CELL_META$cell.ylim[1], CELL_META$cell.xlim[2], CELL_META$cell.ylim[2], col = "black", lty = 1)
#   
#   pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
#   # circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df53$ieg_like_col[CELL_META$sector.numeric.index], border = NA)
#   circos.barplot(bb_df53$prop[CELL_META$sector.numeric.index], CELL_META$xcenter, col = bb_df53$col[CELL_META$sector.numeric.index])
# })
# circos.yaxis(at = track1_breaks, sector.index = "1.1_Astro", track.index = 1, side = "left")
# circos.text(-1, 0.5, sector.index = "1.1_Astro", "R", track.index = 1)
# # IEG
# circos.track(ylim = c(0, 1), track.height = 0.10, bg.border = NA, panel.fun = function(x, y) {
#   pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
#   circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df53$ieg_like_col[CELL_META$sector.numeric.index], border = NA)
# })
# circos.text(-0.5, 0.5, sector.index = "1.1_Astro", "IEG", track.index = 2)
# # Neurogenesis
# circos.track(ylim = c(0, 1), track.height = 0.10, bg.border = NA, panel.fun = function(x, y) {
#   pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
#   # circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df53$neuro_gen_col[CELL_META$sector.numeric.index], border = NA)
#   
#   neuro_gen_sum = bb_df53$neuro_gen_b[CELL_META$sector.numeric.index] + bb_df53$neuro_gen_c[CELL_META$sector.numeric.index]
#   circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], bb_df53$neuro_gen_b[CELL_META$sector.numeric.index]/neuro_gen_sum, col = bb_df53$neuro_gen_col[CELL_META$sector.numeric.index], border = NA)
#   circos.rect(CELL_META$cell.xlim[1], bb_df53$neuro_gen_b[CELL_META$sector.numeric.index]/neuro_gen_sum, CELL_META$cell.xlim[2], 1, col = "gray80", border = NA)
# })
# circos.text(-.5, 0.5, sector.index = "1.1_Astro", "Neurogen", track.index = 3)
# # Cluster Colors
# circos.track(ylim = c(0, 1), track.height = 0.15, bg.border = NA, panel.fun = function(x, y) {
#   pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
#   circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df53$col[CELL_META$sector.numeric.index], border = NA)
#   circos.text(CELL_META$xcenter, CELL_META$ycenter,
#               bb_df53$cluster_label[CELL_META$sector.numeric.index], facing = "downward", niceFacing = TRUE,
#               adj = c(0.5, 0.5), cex = 0.6)
# })
# dev.off()

#*******************************************************************************
# bDEGS ========================================================================
#*******************************************************************************
# badeg = read.csv("~/research/brain/results/deg_depth_build_badeg_glmmseq_demux_all_clusters_all_tests_pair_subjectinpair_pool_subjectinpool_sig_all_genes_100821_q_hgnc.csv")
# all_deg_gsi = read.csv("C:/Users/miles/Downloads/bb15_deg_glmmseq_demux_all_clusters_cond_gsi_control_pair_subjectinpair_pool_subjectinpool_good_genes_by_pair_101321_q_by_cluster.csv")
# all_deg_spawn = read.csv("C:/Users/miles/Downloads/bb15_deg_glmmseq_demux_all_clusters_cond_spawn_control_pair_subjectinpair_pool_subjectinpool_good_genes_by_pair_101321_q_by_cluster.csv")
# bdeg = read.csv("C:/Users/miles/Downloads/bb15_glmmseq_cond_gsi_control_and_cond_spawn_control_sig_genes_q_by_cluster_100521.csv")
bdeg = read.csv("C:/Users/miles/Downloads/out_glmmseq_bb15_demux_deg_bower_behavior_hmp_calculated_across_clusters_111821_hgnc.csv")
all_deg_gsi = read.csv("C:/Users/miles/Downloads/deg_glmmseq_demux_all_clusters_cond_gsi_control_pair_subjectinpair_pool_subjectinpool_good_genes_by_pair_101121_q_by_cluster.csv")
all_deg_spawn = read.csv("C:/Users/miles/Downloads/deg_glmmseq_demux_all_clusters_cond_spawn_control_pair_subjectinpair_pool_subjectinpool_good_genes_by_pair_101121_q_by_cluster.csv")
all_deg_spawn = read.csv("C:/Users/miles/Downloads/deg_glmmseq_demux_all_clusters_cond_spawn_control_pair_subjectinpair_pool_subjectinpool_good_genes_by_pair_101121_q_by_cluster (1).csv")
bdeg = read.csv("C:/Users/miles/Downloads/bb15_glmmseq_cond_gsi_control_and_cond_spawn_control_sig_genes_q_by_cluster_100521.csv")
bdeg = read.csv("~/research/brain/results/bb15_glmmseq_cond_gsi_control_and_cond_spawn_control_sig_genes_q_by_cluster_100521.csv")
all_deg_gsi$cluster_gene = paste0(all_deg_gsi$cluster, "_", all_deg_gsi$mzebra)
all_deg_spawn$cluster_gene = paste0(all_deg_spawn$cluster, "_", all_deg_spawn$mzebra)
bdeg$cluster_gene = paste0(bdeg$cluster, "_", bdeg$mzebra)
bdeg53 = read.csv("C:/Users/miles/Downloads/out_glmmseq_bb53_demux_deg_bower_behavior_hmp_calculated_across_clusters_111821_hgnc.csv")
bdeg53$cluster_gene = paste0(bdeg53$cluster, "_", bdeg53$mzebra)
sdeg = read.csv("C:/Users/miles/Downloads/out_glmmseq_bb15_demux_deg_log_spawn_events_hmp_calculated_across_clusters_111821_hgnc.csv")
sdeg$cluster_gene = paste0(sdeg$cluster, "_", sdeg$mzebra)
gdeg = read.csv("C:/Users/miles/Downloads/out_glmmseq_bb15_demux_deg_gsi_hmp_calculated_across_clusters_11182_hgnc.csv")
gdeg$cluster_gene = paste0(gdeg$cluster, "_", gdeg$mzebra)

bdeg_post = read.csv("C:/Users/miles/Downloads/bdeg_for_volcano_11_30_21.csv")
bdeg$avg_logFC_pool = bdeg_post$avg_logFC

all_pairs = sort(unique(bb$pair))
bdeg[,c(paste0("b", all_pairs), paste0("c", all_pairs))] = 0
for (this_clust in sort(unique(bdeg$cluster))) {
  this_clust_degs = bdeg$mzebra[which(bdeg$cluster == this_clust)]
  print(paste0("Cluster: ", this_clust, ". Num DEGs BAD = ", length(which(! this_clust_degs %in% rownames(bb))), ". Pct DEGs BAD = ", length(which(! this_clust_degs %in% rownames(bb))) / length(this_clust_degs) * 100))
  this_clust_degs = this_clust_degs[which(this_clust_degs %in% rownames(bb))]
  # print("")
  for (pair in all_pairs) {
    # cat(paste0(pair, "."))
    clust_pair_b_cells = colnames(bb)[which(bb$seuratclusters15 == this_clust & bb$pair == pair & bb$cond == "BHVE")]
    clust_pair_c_cells = colnames(bb)[which(bb$seuratclusters15 == this_clust & bb$pair == pair & bb$cond == "CTRL")]
    this_res = pct_dif_avg_logFC(bb, cells.1 = clust_pair_b_cells, cells.2 = clust_pair_c_cells, features = this_clust_degs)
    this_res$cluster = this_clust
    this_res$cluster_gene = paste0(this_res$cluster, "_", this_res$genes)
    bdeg[match(this_res$cluster_gene, bdeg$cluster_gene),paste0("b", pair)] = this_res$avg_logFC
    bdeg[match(this_res$cluster_gene, bdeg$cluster_gene),paste0("c", pair)] = -this_res$avg_logFC
  }
  # print("")
}
bdeg$bvc_fc_dif = sapply(1:nrow(bdeg), function(x) sum(bdeg[x, paste0("b", all_pairs)]) - sum(bdeg[x, paste0("c", all_pairs)]) )
bdeg$abs_bvc_fc_dif = abs(bdeg$bvc_fc_dif)
bdeg$avg_b_fc = rowMeans(bdeg[, paste0("b", all_pairs)])
bdeg$num_b_up = rowSums(sign(bdeg[, paste0("b", all_pairs)]) == 1)
bdeg$neg_log_hmp_cond = -log10(bdeg$hmp_cond)
bdeg_cond = bdeg[which(bdeg$sig_cond == 3),]
bdeg_melt = melt(bdeg_cond[, c("cluster_gene", "cluster", "mzebra", "hmp_cond", "neg_log_hmp_cond", "avg_b_fc", "avg_logFC_pool", paste0("b", all_pairs), paste0("c", all_pairs))], id.var = c('cluster_gene', 'mzebra', 'cluster', 'hmp_cond', 'neg_log_hmp_cond', 'avg_b_fc', 'avg_logFC_pool'))
bdeg_melt$cluster_gene = factor(bdeg_melt$cluster_gene, levels = bdeg_cond$cluster_gene[order(bdeg_cond$bvc_fc_dif, decreasing = T)])
bdeg_melt$cond = startsWith(as.vector(bdeg_melt$variable), "b")
pdf("C:/Users/miles/Downloads/volcano_pool.pdf", width = 10, height = 8)
ggplot(bdeg_melt, aes(x = avg_logFC_pool, y = neg_log_hmp_cond)) + geom_point() + xlab("Gene Cluster Combo") + theme_light() + scale_y_sqrt() + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab("-Log10(HMP COND)") 
dev.off()

ggplot(bdeg_melt, aes(x = avg_b_fc, avg_logFC_pool)) + geom_point()
cor(bdeg_melt$avg_b_fc[which(! is.na(bdeg_melt$avg_b_fc) & ! is.na(bdeg_melt$avg_logFC_pool))], bdeg_melt$avg_logFC_pool[which(! is.na(bdeg_melt$avg_b_fc) & ! is.na(bdeg_melt$avg_logFC_pool))])
# bdeg_melt$cluster = factor(bdeg_melt$cluster)
# pdf("C:/Users/miles/Downloads/volcano_test_cluster.pdf", width = 10, height = 8)
# ggplot(bdeg_melt, aes(x = cluster_gene, y = value, color = cluster)) + geom_point() + xlab("Gene Cluster Combo") + ylab("LOG FC") + theme_light() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# dev.off()

bhve_cells = colnames(bb)[which(bb$cond == "BHVE")]
ctrl_cells = colnames(bb)[which(bb$cond == "CTRL")]
all_pct_fc = data.frame()
for (this_clust in sort(unique(bdeg$cluster))) {
  this_clust_cells = colnames(bb)[which(bb$seuratclusters15 == this_clust)]
  this_res = pct_dif_avg_logFC(bb, cells.1 = this_clust_cells[which(this_clust_cells %in% bhve_cells)], cells.2 = this_clust_cells[which(this_clust_cells %in% ctrl_cells)])
  this_res$cluster = this_clust
  all_pct_fc = rbind(all_pct_fc, this_res)
}
all_pct_fc$cluster_gene = paste0(all_pct_fc$cluster, "_", all_pct_fc$genes)
bdeg[,colnames(all_pct_fc)] = all_pct_fc[match(bdeg$cluster_gene, all_pct_fc$cluster_gene),]
bdeg$cond_sig = bdeg$hmp_cond < 0.05 & bdeg$sig_cond == 3
bdeg$bai_sig =  bdeg$hmp_bower_activity_index < 0.05 & bdeg$sig_bower_activity_index == 3
bdeg$bower_sig =  bdeg$hmp_bower_all < 0.05 & bdeg$sig_bower_all == 6
bdeg$col = convert15$col[match(bdeg$cluster, convert15$old)]
bdeg$col[which(! bdeg$cond_sig )] = "gray80"
bdeg$col_bai = convert15$col[match(bdeg$cluster, convert15$old)]
bdeg$col_bai[which(! bdeg$bai_sig )] = "gray80"
bdeg$col_bower = convert15$col[match(bdeg$cluster, convert15$old)]
bdeg$col_bower[which(! bdeg$bower_sig )] = "gray80"
bdeg$neg_log_hmp_cond = -log10(bdeg$hmp_cond)
bdeg$neg_log_hmp_bai = -log10(bdeg$hmp_bower_activity_index)
bdeg$neg_log_hmp_bower = -log10(bdeg$hmp_bower_all)

bhve_cells = colnames(bb)[which(bb$cond == "BHVE")]
ctrl_cells = colnames(bb)[which(bb$cond == "CTRL")]
all_pct_fc = data.frame()
for (this_clust in 0:52) {
  print(this_clust)
  this_clust_cells = colnames(bb)[which(bb$seuratclusters53 == this_clust)]
  this_res = pct_dif_avg_logFC(bb, cells.1 = this_clust_cells[which(this_clust_cells %in% bhve_cells)], cells.2 = this_clust_cells[which(this_clust_cells %in% ctrl_cells)])
  this_res$cluster = this_clust
  all_pct_fc = rbind(all_pct_fc, this_res)
}
all_pct_fc$cluster_gene = paste0(all_pct_fc$cluster, "_", all_pct_fc$genes)
bdeg53[,colnames(all_pct_fc)] = all_pct_fc[match(bdeg53$cluster_gene, all_pct_fc$cluster_gene),]
bdeg53$cond_sig = bdeg53$hmp_cond < 0.05 & bdeg53$sig_cond == 3
bdeg53$bai_sig =  bdeg53$hmp_bower_activity_index < 0.05 & bdeg53$sig_bower_activity_index == 3
bdeg53$bower_sig =  bdeg53$hmp_bower_all < 0.05 & bdeg53$sig_bower_all == 6
bdeg53$col = convert15$col[match(bdeg53$cluster, convert15$old)]
bdeg53$col[which(! bdeg53$cond_sig )] = "gray80"
bdeg53$col_bai = convert15$col[match(bdeg53$cluster, convert15$old)]
bdeg53$col_bai[which(! bdeg53$bai_sig )] = "gray80"
bdeg53$col_bower = convert15$col[match(bdeg53$cluster, convert15$old)]
bdeg53$col_bower[which(! bdeg53$bower_sig )] = "gray80"
bdeg53$neg_log_hmp_cond = -log10(bdeg53$hmp_cond)
bdeg53$neg_log_hmp_bai = -log10(bdeg53$hmp_bower_activity_index)
bdeg53$neg_log_hmp_bower = -log10(bdeg53$hmp_bower_all)

# Spawn DEGs
sub_spawn = aggregate(log_spawn_events ~ subsample, bb@meta.data, mean)
sdeg$spawn_cor = 0
for (this_clust in 0:14) {
  print(this_clust)
  this_clust_colnames = colnames(all_avg)[which( startsWith(colnames(all_avg), as.character(this_clust)) )]
  this_spawn = sub_spawn$log_spawn_events[match(reshape2::colsplit(this_clust_colnames, "_", c('1', '2'))[,2], sub_spawn$subsample)]
  clust_avg = all_avg[, this_clust_colnames]
  sdeg_clust = sdeg[which(sdeg$cluster == this_clust),]
  if (nrow(sdeg_clust) > 0) {
    spawn_cor = sapply(sdeg_clust$mzebra, function(x) cor(as.numeric(clust_avg[x,]), this_spawn) )
    sdeg$spawn_cor[match(sdeg_clust$cluster_gene, sdeg$cluster_gene)] = spawn_cor
  }
}
sdeg$isSig = sdeg$hmp < 0.05 & sdeg$sig_log_spawn_events == 5
sdeg$col = convert15$col[match(sdeg$cluster, convert15$old)]
sdeg$col[which(! sdeg$isSig )] = "gray80"
sdeg$neg_log_hmp = -log10(sdeg$hmp)
sdeg = sdeg[order(sdeg$isSig, decreasing = F),]

# GSI DEGs
Idents(bb) = paste0(bb$seuratclusters15, "_", bb$subsample)
sub_gsi = aggregate(gsi ~ subsample, bb@meta.data, mean)
all_avg = myAverageExpression(bb)
# all_avg[,c("cluster", "subsample")] = reshape2::colsplit(all_avg, "_", c('1', '2'))
gdeg$gsi_cor = 0
for (this_clust in 0:14) {
  print(this_clust)
  this_clust_colnames = colnames(all_avg)[which( startsWith(colnames(all_avg), as.character(this_clust)) )]
  this_gsi = sub_gsi$gsi[match(reshape2::colsplit(this_clust_colnames, "_", c('1', '2'))[,2], sub_gsi$subsample)]
  clust_avg = all_avg[, this_clust_colnames]
  gdeg_clust = gdeg[which(gdeg$cluster == this_clust),]
  if (nrow(gdeg_clust) > 0) {
    gsi_cor = sapply(gdeg_clust$mzebra, function(x) cor(as.numeric(clust_avg[x,]), this_gsi) )
    gdeg$gsi_cor[match(gdeg_clust$cluster_gene, gdeg$cluster_gene)] = gsi_cor
  }
}
gdeg$isSig = gdeg$hmp < 0.05 & gdeg$sig_gsi == 5
gdeg$col = convert15$col[match(gdeg$cluster, convert15$old)]
gdeg$col[which(! gdeg$isSig )] = "gray80"
gdeg$neg_log_hmp = -log10(gdeg$hmp)
gdeg = gdeg[order(gdeg$isSig, decreasing = F),]

clust_nums = data.frame(cluster = 0:14, clust_name = convert15$new.full[match(0:14, convert15$old)], col = convert15$col[match(0:14, convert15$old)],
                        COND = as.numeric(table(factor(bdeg$cluster[which(bdeg$cond_sig)], levels = 0:14))),
                        BAI = as.numeric(table(factor(bdeg$cluster[which(bdeg$bai_sig)], levels = 0:14))),
                        GSI = as.numeric(table(factor(gdeg$cluster[which(gdeg$isSig)], levels = 0:14))),
                        SPAWN = as.numeric(table(factor(sdeg$cluster[which(sdeg$isSig)], levels = 0:14))))
clust_nums_melt = melt(clust_nums)
clust_nums_melt = clust_nums_melt[which(clust_nums_melt$variable != "cluster"),]
clust_nums_melt$clust_name = factor(clust_nums_melt$clust_name, levels = convert15$new.full)
# clust_nums_melt$variable = plyr::revalue(clust_nums_melt$variable, replace = c("bdeg"))
pdf("C:/Users/miles/Downloads/deg_nums.pdf", width = 13, height = 6)
ggplot(clust_nums_melt, aes(x = clust_name, y = value, fill = variable, color = variable, group = variable)) + geom_point(size = 2) + geom_line() + theme_classic() + xlab("") + ylab("Number of DEGs")
dev.off()

# Volcano Plot
bdeg = bdeg[order(bdeg$cond_sig, decreasing = F),]
pdf("C:/Users/miles/Downloads/bb_cond_bdeg_volcano_15.pdf", height = 6, width = 5)
ggplot(bdeg, aes(x = avg_logFC, y = neg_log_hmp_cond, color = col, alpha = cond_sig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none')
dev.off()
pdf("C:/Users/miles/Downloads/bb_cond_bdeg_volcano_15_label.pdf", height = 6, width = 5)
ggplot(bdeg, aes(x = avg_logFC, y = neg_log_hmp_cond, color = col, alpha = cond_sig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none') + geom_text_repel(data = bdeg[which(bdeg$cond_sig & bdeg$neg_log_hmp_cond > 100),] ,aes(label = cluster_gene))
dev.off()
bdeg = bdeg[order(bdeg$bai_sig, decreasing = F),]
pdf("C:/Users/miles/Downloads/bb_bai_deg_volcano_15.pdf", height = 6, width = 5)
ggplot(bdeg, aes(x = avg_logFC, y = neg_log_hmp_bai, color = col_bai, alpha = bai_sig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none')
dev.off()
pdf("C:/Users/miles/Downloads/bb_bai_deg_volcano_15_label.pdf", height = 6, width = 5)
ggplot(bdeg, aes(x = avg_logFC, y = neg_log_hmp_bai, color = col_bai, alpha = bai_sig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none') + geom_text_repel(data = bdeg[which(bdeg$bai_sig & bdeg$neg_log_hmp_bai > 10),] ,aes(label = cluster_gene))
dev.off()
bdeg = bdeg[order(bdeg$bower_sig, decreasing = F),]
pdf("C:/Users/miles/Downloads/bb_bower_deg_volcano_15.pdf", height = 6, width = 5)
ggplot(bdeg, aes(x = avg_logFC, y = neg_log_hmp_bower, color = col_bower, alpha = bower_sig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none')
dev.off()
pdf("C:/Users/miles/Downloads/bb_bower_deg_volcano_15_label.pdf", height = 6, width = 5)
ggplot(bdeg, aes(x = avg_logFC, y = neg_log_hmp_bower, color = col_bower, alpha = bower_sig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none') + geom_text_repel(data = bdeg[which(bdeg$bower_sig & bdeg$neg_log_hmp_bower > 20),] ,aes(label = cluster_gene))
dev.off()

pdf("C:/Users/miles/Downloads/bb_gsi_deg_volcano_15.pdf", height = 6, width = 5)
ggplot(gdeg, aes(x = gsi_cor, y = neg_log_hmp, color = col, alpha = isSig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab("Correlation with GSI") + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none')
dev.off()
pdf("C:/Users/miles/Downloads/bb_gsi_deg_volcano_15_label.pdf", height = 6, width = 5)
ggplot(gdeg, aes(x = gsi_cor, y = neg_log_hmp, color = col, alpha = isSig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab("Correlation with GSI") + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none') + geom_text_repel(data = gdeg[which(gdeg$isSig & gdeg$neg_log_hmp > 200),] ,aes(label = cluster_gene))
dev.off()
pdf("C:/Users/miles/Downloads/bb_spawn_deg_volcano_15.pdf", height = 6, width = 5)
ggplot(sdeg, aes(x = spawn_cor, y = neg_log_hmp, color = col, alpha = isSig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab("Correlation with Log Spawn Events") + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none')
dev.off()
pdf("C:/Users/miles/Downloads/bb_spawn_deg_volcano_15_label.pdf", height = 6, width = 5)
ggplot(sdeg, aes(x = spawn_cor, y = neg_log_hmp, color = col, alpha = isSig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab("Correlation with Log Spawn Events") + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none') + geom_text_repel(data = sdeg[which(sdeg$isSig & sdeg$neg_log_hmp > 75),] ,aes(label = cluster_gene))
dev.off()
# pdf("C:/Users/miles/Downloads/bb_bdeg_volcano_gsi_spawn_max_53.pdf", height = 6, width = 5)
# ggplot(all_deg_join, aes(x = avg_logFC.gsi, y = neg_log_q_max, color = col2, alpha = isSig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none')
# dev.off()
# pdf("C:/Users/miles/Downloads/bb_bdeg_volcano_gsi_53.pdf", height = 6, width = 5)
# ggplot(all_deg_join, aes(x = avg_logFC.gsi, y = neg_log_q_gsi, color = col2, alpha = isSig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none')
# dev.off()
# pdf("C:/Users/miles/Downloads/bb_bdeg_volcano_spawn_53.pdf", height = 6, width = 5)
# ggplot(all_deg_join, aes(x = avg_logFC.gsi, y = neg_log_q_spawn, color = col2, alpha = isSig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none')
# dev.off()
# pdf("C:/Users/miles/Downloads/bb_bdeg_gsi_spawn_53.pdf", height = 6, width = 5)
# ggplot(all_deg_join, aes(x = neg_log_q_spawn, y = neg_log_q_gsi, color = col2, alpha = isSig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab(expression(-Log["10"]*" Adjusted P (Spawn)")) + ylab(expression(-Log["10"]*" Adjusted P (GSI)")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none')
# dev.off()

#53
bdeg5353 = bdeg53[order(bdeg53$cond_sig, decreasing = F),]
pdf("C:/Users/miles/Downloads/bb_cond_bdeg53_volcano_53.pdf", height = 6, width = 5)
ggplot(bdeg53, aes(x = avg_logFC, y = neg_log_hmp_cond, color = col, alpha = cond_sig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none')
dev.off()
pdf("C:/Users/miles/Downloads/bb_cond_bdeg53_volcano_53_label.pdf", height = 6, width = 5)
ggplot(bdeg53, aes(x = avg_logFC, y = neg_log_hmp_cond, color = col, alpha = cond_sig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none') + geom_text_repel(data = bdeg53[which(bdeg53$cond_sig & bdeg53$neg_log_hmp_cond > 100),] ,aes(label = cluster_gene))
dev.off()
bdeg53 = bdeg53[order(bdeg53$bai_sig, decreasing = F),]
pdf("C:/Users/miles/Downloads/bb_bai_deg_volcano_53.pdf", height = 6, width = 5)
ggplot(bdeg53, aes(x = avg_logFC, y = neg_log_hmp_bai, color = col_bai, alpha = bai_sig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none')
dev.off()
pdf("C:/Users/miles/Downloads/bb_bai_deg_volcano_53_label.pdf", height = 6, width = 5)
ggplot(bdeg53, aes(x = avg_logFC, y = neg_log_hmp_bai, color = col_bai, alpha = bai_sig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none') + geom_text_repel(data = bdeg53[which(bdeg53$bai_sig & bdeg53$neg_log_hmp_bai > 10),] ,aes(label = cluster_gene))
dev.off()
bdeg53 = bdeg53[order(bdeg53$bower_sig, decreasing = F),]
pdf("C:/Users/miles/Downloads/bb_bower_deg_volcano_53.pdf", height = 6, width = 5)
ggplot(bdeg53, aes(x = avg_logFC, y = neg_log_hmp_bower, color = col_bower, alpha = bower_sig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none')
dev.off()
pdf("C:/Users/miles/Downloads/bb_bower_deg_volcano_53_label.pdf", height = 6, width = 5)
ggplot(bdeg53, aes(x = avg_logFC, y = neg_log_hmp_bower, color = col_bower, alpha = bower_sig)) + geom_point() + scale_y_sqrt() + scale_color_identity() + xlab(expression(Log["2"]*" Fold Change")) + ylab(expression(-Log["10"]*" Adjusted P")) + theme_light() + scale_alpha_manual(values = c(0.5, 0.9), drop = F, guide = 'none') + geom_text_repel(data = bdeg53[which(bdeg53$bower_sig & bdeg53$neg_log_hmp_bower > 20),] ,aes(label = cluster_gene))
dev.off()
# all_deg_gsi = read.csv("C:/Users/miles/Downloads/deg_glmmseq_demux_all_clusters_cond_gsi_control_pair_subjectinpair_pool_subjectinpool_good_genes_by_pair_101121_q_by_cluster (1).csv")
# all_deg_spawn = read.csv("C:/Users/miles/Downloads/deg_glmmseq_demux_all_clusters_cond_spawn_control_pair_subjectinpair_pool_subjectinpool_good_genes_by_pair_101121_q_by_cluster (2).csv")
# bdeg = read.csv("C:/Users/miles/Downloads/deg_glmmseq_demux_all_cond_sig_53clusters_pair_subjectinpair_pool_subjectinpool_good_genes_by_pair_101121_q_by_cluster_hgnc.csv")
# 
# all_deg_gsi$cluster_gene = paste0(all_deg_gsi$cluster, "_", all_deg_gsi$mzebra)
# all_deg_spawn$cluster_gene = paste0(all_deg_spawn$cluster, "_", all_deg_spawn$mzebra)
# all_deg_gsi = all_deg_gsi[which(! is.na(all_deg_gsi$cluster_gene) ),]
# all_deg_spawn = all_deg_spawn[which(! is.na(all_deg_spawn$cluster_gene) ),]
# bdeg$cluster_gene = paste0(bdeg$cluster, "_", bdeg$mzebra)
# all_deg_gsi[,colnames(all_pct_fc)] = all_pct_fc[match(all_deg_gsi$cluster_gene, all_pct_fc$cluster_gene),]
# all_deg_spawn[,colnames(all_pct_fc)] = all_pct_fc[match(all_deg_spawn$cluster_gene, all_pct_fc$cluster_gene),]
#
# all_deg_gsi = all_deg_gsi[which(! is.na(all_deg_gsi$cluster_gene) ),]
# all_deg_spawn = all_deg_spawn[which(! is.na(all_deg_spawn$cluster_gene) ),]
# all_deg_join = inner_join(all_deg_gsi, all_deg_spawn, by = "cluster_gene", suffix = c(".gsi", ".spawn"))
# all_deg_join$neg_log_q_gsi = -log10(all_deg_join$q_cond.gsi)
# all_deg_join$neg_log_q_spawn = -log10(all_deg_join$q_cond.spawn)
# # all_deg_join$col = convert15$col[match(all_deg_join$cluster.gsi, convert15$old)]
# all_deg_join$col = convert53$col[match(all_deg_join$cluster.gsi, convert53$old)]
# all_deg_join$isSig = all_deg_join$cluster_gene %in% bdeg$cluster_gene
# all_deg_join$col2 = all_deg_join$col
# all_deg_join$col2[which(!all_deg_join$isSig)] = "gray"
# all_deg_join = all_deg_join[order(all_deg_join$isSig, decreasing = F),]
# all_deg_join$q_combined = sapply(1:nrow(all_deg_join), function(x) mean(c(all_deg_join$q_cond.gsi[x], all_deg_join$q_cond.spawn[x])))
# all_deg_join$q_max = sapply(1:nrow(all_deg_join), function(x) max(c(all_deg_join$q_cond.gsi[x], all_deg_join$q_cond.spawn[x])))
# all_deg_join$neg_log_q_max = -log10(all_deg_join$q_max)
# all_deg_join$neg_log_q_combined2 = -log10(all_deg_join$q_combined)
# all_deg_join$neg_log_q_combined = sapply(1:nrow(all_deg_join), function(x) mean(c(all_deg_join$neg_log_q_gsi[x], all_deg_join$neg_log_q_spawn[x])))
# head(all_deg_join[,c("neg_log_q_gsi", "neg_log_q_spawn", "neg_log_q_combined", "neg_log_q_combined2", "neg_log_q_max")])


# subsample_meta = unique(bb@meta.data[,c("subsample", "pair", "pool", "gsi", "standard_length", "cond", colnames(bb@meta.data)[which(startsWith(colnames(bb@meta.data), "build"))], colnames(bb@meta.data)[which(startsWith(colnames(bb@meta.data), "depth"))], colnames(bb@meta.data)[which(startsWith(colnames(bb@meta.data), "spawn"))] )])
subsample_meta = unique(bb@meta.data[,c("subsample", "pair", "run", "gsi", "standard_length", "cond", "log_spawn_events", "bower_activity_index", "build_events", "depth", "depth_adj", "spawn_events" )])
rownames(subsample_meta) = subsample_meta$subsample
subsample_meta = subsample_meta[order(subsample_meta$subsample),]
subsample_pairs = list()
for (i in 1:18) {
  subsample_pairs[[i]] = subsample_meta$subsample[which(subsample_meta$pair == i)]
}

# Heatmap with Hierarchical Clustering for 15 Cluster BVC DEGs
heat15 = data.frame()
# for (i in 0:14) {
for (i in 0:52) {  
  this_clust_bdeg = bdeg[which(bdeg$cluster == i),]
  if ( nrow(this_clust_bdeg) > 0 ) {
    print(i)
    pairs_res = data.frame()
    for (p in unique(subsample_meta$pair)) {
      # this_pair_res = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$seuratclusters15 == i & bb$pair == p & bb$cond == "BHVE")], cells.2 = colnames(bb)[which(bb$seuratclusters15 == i & bb$pair == p & bb$cond == "CTRL")])
      this_pair_res = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$seuratclusters53 == i & bb$pair == p & bb$cond == "BHVE")], cells.2 = colnames(bb)[which(bb$seuratclusters53 == i & bb$pair == p & bb$cond == "CTRL")])
      this_pair_res = this_pair_res[this_clust_bdeg$genes,]
      this_pair_res$cluster = i
      this_pair_res$pair = p
      this_pair_res = rbind(this_pair_res, this_pair_res)
      this_pair_res$cond_pair = c(rep(paste0("b", p), nrow(this_pair_res)/2), rep(paste0("c", p), nrow(this_pair_res)/2))
      this_pair_res$avg_logFC[(nrow(this_pair_res)/2):nrow(this_pair_res)] = -this_pair_res$avg_logFC[(nrow(this_pair_res)/2):nrow(this_pair_res)]
      pairs_res = rbind(pairs_res, this_pair_res)
    }
    # for (p in 1:5) {
    #   # this_pair_res = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$seuratclusters15 == i & bb$pool == paste0("pair", p) & bb$cond == "BHVE")], cells.2 = colnames(bb)[which(bb$seuratclusters15 == i & bb$pool == paste0("pair", p) & bb$cond == "CTRL")])
    #   this_pair_res = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$seuratclusters53 == i & bb$pool == paste0("pair", p) & bb$cond == "BHVE")], cells.2 = colnames(bb)[which(bb$seuratclusters53 == i & bb$pool == paste0("pair", p) & bb$cond == "CTRL")])
    #   this_pair_res = this_pair_res[this_clust_bdeg$genes,]
    #   this_pair_res$cluster = i
    #   this_pair_res$pool = p
    #   this_pair_res = rbind(this_pair_res, this_pair_res)
    #   this_pair_res$cond_pool = c(rep(paste0("b_", p), nrow(this_pair_res)/2), rep(paste0("c_", p), nrow(this_pair_res)/2))
    #   this_pair_res$avg_logFC[(nrow(this_pair_res)/2):nrow(this_pair_res)] = -this_pair_res$avg_logFC[(nrow(this_pair_res)/2):nrow(this_pair_res)]
    #   pairs_res = rbind(pairs_res, this_pair_res)
    # }
    heat15 = rbind(heat15, pairs_res)
  }
}
heat15$cluster_gene = paste0(heat15$cluster, "_", heat15$genes)
heat15_wide = reshape2::acast(cluster_gene ~ cond_pair, data = heat15, value.var = "avg_logFC")
# heat15_wide[which(is.na(heat15_wide))] = 0
# heat15_wide = heat15_wide[, c(1, 9:12, 14:18)]
heat15_wide = abs(heat15_wide > 0)

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}


mat_breaks <- quantile_breaks(heat15_wide, n = 11)
pheatmap::pheatmap(heat15_wide, border_color = NA, angle_col = 0, show_rownames = F, color = viridis(length(mat_breaks) - 1), breaks = mat_breaks)


#******************************************************************
# Volcano 2.0 =====================================================
#******************************************************************
bb15 = read.csv("~/Downloads/out_bb15_bbmm_demux_deg_all_tests_for_volcano_plotting_121321.csv")
bb53 = read.csv("~/Downloads/out_bb53_glmmseq_demux_deg_all_tests_for_volcano_plotting.csv")

bb15$neg_log_bower = -log10(bb15$P_bower_activity_index)
bb15$neg_log_gsi   = -log10(bb15$P_gsi)
bb15$neg_log_spawn = -log10(bb15$P_log_spawn_events)

bb53$neg_log_bower = -log10(bb53$P_bower_activity_index)
bb53$neg_log_gsi   = -log10(bb53$P_gsi)
bb53$neg_log_spawn = -log10(bb53$P_log_spawn_events)

non_sig_col = "gray80"
bb15$col = convert15$col[match(bb15$cluster, convert15$old)]
bb15$col_bower = bb15$col
bb15$col_gsi   = bb15$col
bb15$col_spawn = bb15$col
bb15$col_bower[which(bb15$sig_bower_behavior == 0)]   = non_sig_col
bb15$col_gsi[which(bb15$sig_gsi == 0)]                = non_sig_col
bb15$col_spawn[which(bb15$sig_log_spawn_events == 0)] = non_sig_col

bb53$col = convert53$col[match(bb53$cluster, convert53$old)]
bb53$col_bower = bb53$col
bb53$col_gsi   = bb53$col
bb53$col_spawn = bb53$col
bb53$col_bower[which(bb53$sig_bower_behavior == 0)]   = 
bb53$col_gsi[which(bb53$sig_gsi == 0)]                = non_sig_col
bb53$col_spawn[which(bb53$sig_log_spawn_events == 0)] = non_sig_col

# 15 Cluster Level Plots
bb15 = bb15[order(bb15$sig_bower_behavior),]
bb15$neg_log_bower_10 = bb15$neg_log_bower ^ (1/10)
brk_f = pretty_breaks(n = 3)
bower_labels = brk_f(bb15$neg_log_bower_10, n=3)
pdf("~/research/brain/results/bb15_bower_volcano_like_10th.pdf", width = 6, height = 5)
ggplot(bb15, aes(x = bower_activity_index, y = neg_log_bower_10, color = col_bower, alpha = as.logical(sig_bower_behavior))) + geom_point() + scale_color_identity() + xlab("Beta Estimate") + ylab(expression(-Log["10"]*" P")) + scale_alpha_manual(values = c(0.2, 0.75), guide = 'none') + scale_y_continuous(labels = paste0(as.character(bower_labels), "^10"), breaks = bower_labels) + theme_bw()
dev.off()

bb15 = bb15[order(bb15$sig_gsi),]
bb15$neg_log_gsi_10 = bb15$neg_log_gsi ^ (1/10)
brk_f = pretty_breaks(n = 3)
gsi_labels = brk_f(bb15$neg_log_gsi_10, n=3)
pdf("~/research/brain/results/bb15_gsi_volcano_like_10th.pdf", width = 6, height = 5)
ggplot(bb15, aes(x = gsi, y = neg_log_gsi_10, color = col_gsi, alpha = as.logical(sig_gsi))) + geom_point() + scale_color_identity() + xlab("Beta Estimate") + ylab(expression(-Log["10"]*" P")) + scale_alpha_manual(values = c(0.2, 0.75), guide = 'none') + scale_y_continuous(labels = paste0(as.character(gsi_labels), "^10"), breaks = gsi_labels) + theme_bw()
dev.off()

bb15 = bb15[order(bb15$sig_log_spawn_events),]
bb15$neg_log_spawn_10 = bb15$neg_log_spawn ^ (1/10)
brk_f = pretty_breaks(n = 3)
spawn_labels = brk_f(bb15$neg_log_spawn_10, n=3)
pdf("~/research/brain/results/bb15_spawn_volcano_like_10th.pdf", width = 6, height = 5)
ggplot(bb15, aes(x = log_spawn_events, y = neg_log_spawn_10, color = col_spawn, alpha = as.logical(sig_log_spawn_events))) + geom_point() + scale_color_identity() + xlab("Beta Estimate") + ylab(expression(-Log["10"]*" P")) + scale_alpha_manual(values = c(0.2, 0.75), guide = 'none') + scale_y_continuous(labels = paste0(as.character(spawn_labels), "^10"), breaks = spawn_labels) + theme_bw()
dev.off()

# 53 Cluster Level Plots
bb53 = bb53[order(bb53$sig_bower_behavior),]
bb53$neg_log_bower_10 = bb53$neg_log_bower ^ (1/10)
brk_f = pretty_breaks(n = 3)
bower_labels = brk_f(bb53$neg_log_bower_10, n=3)
pdf("~/research/brain/results/bb53_bower_volcano_like_10th.pdf", width = 6, height = 5)
ggplot(bb53, aes(x = bower_activity_index, y = neg_log_bower_10, color = col_bower, alpha = as.logical(sig_bower_behavior))) + geom_point() + scale_color_identity() + xlab("Beta Estimate") + ylab(expression(-Log["10"]*" P")) + scale_alpha_manual(values = c(0.2, 0.75), guide = 'none') + scale_y_continuous(labels = paste0(as.character(bower_labels), "^10"), breaks = bower_labels) + theme_bw()
dev.off()

bb53 = bb53[order(bb53$sig_gsi),]
bb53$neg_log_gsi_10 = bb53$neg_log_gsi ^ (1/10)
brk_f = pretty_breaks(n = 3)
gsi_labels = brk_f(bb53$neg_log_gsi_10, n=3)
pdf("~/research/brain/results/bb53_gsi_volcano_like_10th.pdf", width = 6, height = 5)
ggplot(bb53, aes(x = gsi, y = neg_log_gsi_10, color = col_gsi, alpha = as.logical(sig_gsi))) + geom_point() + scale_color_identity() + xlab("Beta Estimate") + ylab(expression(-Log["10"]*" P")) + scale_alpha_manual(values = c(0.2, 0.75), guide = 'none') + scale_y_continuous(labels = paste0(as.character(gsi_labels), "^10"), breaks = gsi_labels) + theme_bw()
dev.off()

bb53 = bb53[order(bb53$sig_log_spawn_events),]
bb53$neg_log_spawn_10 = bb53$neg_log_spawn ^ (1/10)
brk_f = pretty_breaks(n = 3)
spawn_labels = brk_f(bb53$neg_log_spawn_10, n=3)
pdf("~/research/brain/results/bb53_spawn_volcano_like_10th.pdf", width = 6, height = 5)
ggplot(bb53, aes(x = log_spawn_events, y = neg_log_spawn_10, color = col_spawn, alpha = as.logical(sig_log_spawn_events))) + geom_point() + scale_color_identity() + xlab("Beta Estimate") + ylab(expression(-Log["10"]*" P")) + scale_alpha_manual(values = c(0.2, 0.75), guide = 'none') + scale_y_continuous(labels = paste0(as.character(spawn_labels), "^10"), breaks = spawn_labels) + theme_bw()
dev.off()

###################################################################
# Replicate #######################################################
###################################################################
bb = ScaleData(bb, features = c("LOC106675461", "rsrp1"))
bb$aroma = bb@assays$RNA@scale.data["LOC106675461",]
bb_df = data.frame(sample = unique(bb$sample), pair = rep(1:5, 2), col = c("#9d020890", "#d0000090", "#dc2f0290", "#e85d0490", "#f4b90690", "#03045e90", "#023e8a90", "#0077b690", "#0096c790", "#00b4d890"), cond = c(rep("BHVE", 5), rep("CTRL", 5)), num_cells = aggregate(nCount_RNA ~ sample, bb@meta.data, length)[,2], depth = aggregate(depth ~ sample, bb@meta.data, mean)[,2], aroma = aggregate(aroma ~ sample, bb@meta.data, mean)[,2] )
bb_df$prop13 = unlist(sapply(bb_df$sample, function(x) length(which(bb$seuratclusters53 == 13 & bb$sample == x)) ))
bb_df$prop13 = unlist(sapply(bb_df$sample, function(x) bb_df$prop13[which(bb_df$sample == x)] /length(which(bb$seuratclusters15 == 0 & bb$sample == x)) ))
# bb_df$prop13 = bb_df$prop13 / bb_df$num_cells
bb_df$sample = factor(bb_df$sample, levels = c("c4", "c5", "b5", "b4", "b3", "b2", "b1", "c1", "c2", "c3"))
bb_df$aroma_col = scaleValues(bb_df$aroma)
# bb_df$depth[which(bb_df$cond == "CTRL")] = F

myJitter = function(x, jitter_width = 0.15) {
  return(x + runif(1, min = -jitter_width, max = jitter_width))
}

pdf("~/research/brain/results/test2.pdf", width =  10, height = 10)
circos.par("gap.degree" = 10, cell.padding = c(0, 0, 0, 0), "start.degree" = -28, track.margin = c(0.02, 0.02))
circos.initialize(bb_df$pair, xlim = c(0, 1))

# Depth
track1_breaks = pretty(bb_df$depth, n = 3)
circos.track(ylim = c(0, max(track1_breaks)), bg.col = NA, bg.border = "black", track.height = 0.15, panel.fun = function(x, y) {
  for (tb in track1_breaks) {
    if (tb != track1_breaks[1] & tb != track1_breaks[length(track1_breaks)])
      circos.segments(0, tb, 1, tb, col = "gray90")
  }
  
  b_sample = paste0("b", CELL_META$sector.index)
  value = bb_df$depth[which(bb_df$sample == b_sample)]
  circos.barplot(value, CELL_META$xcenter, col = bb_df$col[which(bb_df$sample == b_sample)], bar_width = 0.2, border = NA)
  # c_sample = paste0("c", CELL_META$sector.index)
  # value = bb_df$depth[which(bb_df$sample == c_sample)]
  # circos.barplot(value, CELL_META$xcenter, col = bb_df$col[which(bb_df$sample == c_sample)], bar_width = 0.2, border = NA)
})
circos.yaxis(at = track1_breaks, sector.index = "1", track.index = 1, side = "right")

# Proportion
track2_breaks = c(0, pretty(bb_df$prop13, n = 1))
circos.track(ylim = c(0, max(track2_breaks)), bg.col = NA, bg.border = "black", track.height = 0.15, panel.fun = function(x, y) {
  for (tb in track2_breaks) {
    if (tb != track2_breaks[1] & tb != track2_breaks[length(track2_breaks)])
      circos.segments(0, tb, 1, tb, col = "gray90")
  }
  
  b_sample = paste0("b", CELL_META$sector.index)
  c_sample = paste0("c", CELL_META$sector.index)
  value = bb_df$prop13[which(bb_df$sample == b_sample)]
  circos.barplot(value, CELL_META$xcenter-0.05, col = bb_df$col[which(bb_df$sample == b_sample)], bar_width = 0.05, border = "black")
  value = bb_df$prop13[which(bb_df$sample == c_sample)]
  circos.barplot(value, CELL_META$xcenter+0.05, col = bb_df$col[which(bb_df$sample == c_sample)], bar_width = 0.05, border = "black")
  
})
circos.yaxis(at=track2_breaks, sector.index = "1", track.index = 2, side = "right")

# Aromatase
# hist(bb@assays$RNA@data["rsrp1",])
aroma_mid_bool = bb@assays$RNA@data["rsrp1",] < quantile(bb@assays$RNA@data["rsrp1",], 0.95)
track3_breaks = pretty(bb@assays$RNA@data["rsrp1", which(aroma_mid_bool)], n = 2)
circos.track(ylim = c(min(track3_breaks), max(track3_breaks)), bg.col = NA, bg.border = "black", track.height = 0.15, panel.fun = function(x, y) {
  # Negative zone
  # circos.rect(0, min(track3_breaks), 1, 0, col = "gray70", border = NULL)
  
  for (tb in track3_breaks) {
    if (tb != track3_breaks[1] & tb != track3_breaks[length(track3_breaks)])
      circos.segments(0, tb, 1, tb, col = "gray90")
  }
  
  # value = bb_df$aroma[which(bb_df$sample == CELL_META$sector.index)]
  # circos.barplot(value, CELL_META$xcenter, col = bb_df$aroma_col[which(bb_df$sample == CELL_META$sector.index)], bar_width = 0.2, border = NA)
  b_sample = paste0("b", CELL_META$sector.index)
  c_sample = paste0("c", CELL_META$sector.index)
  b_value = bb@assays$RNA@data["rsrp1", which(bb$sample == b_sample & aroma_mid_bool)]
  c_value = bb@assays$RNA@data["rsrp1", which(bb$sample == c_sample & aroma_mid_bool)]
  circos.violin(b_value, CELL_META$xcenter-0.15, col = bb_df$col[which(bb_df$sample == b_sample)], violin_width = 0.2)
  circos.violin(c_value, CELL_META$xcenter+0.15, col = bb_df$col[which(bb_df$sample == c_sample)], violin_width = 0.2)
  
  for (i in 1:4) {
    subsample = paste0(b_sample, ".", i)
    if (subsample %in% unique(bb$subsample))
      circos.points(myJitter(CELL_META$xcenter-0.15, 0.1), mean(bb@assays$RNA@data["rsrp1", which(bb$subsample == subsample & aroma_mid_bool)]), col = bb_df$col[which(bb_df$sample == b_sample)], pch = 20)
  }
  for (i in 1:4) {
    subsample = paste0(c_sample, ".", i)
    if (subsample %in% unique(bb$subsample))
      circos.points(myJitter(CELL_META$xcenter+0.15, 0.1), mean(bb@assays$RNA@data["rsrp1", which(bb$subsample == subsample & aroma_mid_bool)]), col = bb_df$col[which(bb_df$sample == c_sample)], pch = 20)
  }
})
circos.yaxis(at=track3_breaks, sector.index = "1", track.index = 3, side = "right")

# Sample
# circos.track(ylim = c(0, 1), bg.col = NA, bg.border = NA, track.height = 0.15, panel.fun = function(x, y) {
#   pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
#   circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] - mm_y(2),
#               CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE,
#               adj = c(1, 0.5), cex = 0.6)
#   circos.rect(0, 0, 1, 1, col = bb_df$col[which(bb_df$sample == CELL_META$sector.index)], border = NA)
# })
dev.off()
circos.clear()

#*******************************************************************************
# Supplemental =================================================================
#*******************************************************************************
bb$good15 = convert15$new.full[match(bb$seuratclusters15, convert15$old)]
bb$col15 = convert15$col[match(bb$seuratclusters15, convert15$old)]
df = bb@meta.data
pdf("~/research/brain/results/bb_clust15_nCount.pdf", width = 12, height = 5)
ggplot(df, aes(x = good15_names, y = nCount_RNA, color = col15)) + geom_violin() + geom_boxplot(width = 0.1) + scale_color_identity() + xlab("")
dev.off()
pdf("~/research/brain/results/bb_clust15_nFeature.pdf", width = 12, height = 5)
ggplot(df, aes(x = good15_names, y = nFeature_RNA, color = col15)) + geom_violin() + geom_boxplot(width = 0.1) + scale_color_identity() + xlab("")
dev.off()
pdf("~/research/brain/results/bb_pool_nCount.pdf", width = 10, height = 5)
ggplot(df, aes(x = sample, y = nCount_RNA, color = sample)) + geom_violin() + geom_boxplot(width = 0.1) + xlab("") + scale_color_manual(values = hue_pal()(10), guide = "none")
dev.off()
pdf("~/research/brain/results/bb_pool_nFeature.pdf", width = 10, height = 5)
ggplot(df, aes(x = sample, y = nFeature_RNA, color = sample)) + geom_violin() + geom_boxplot(width = 0.1) + xlab("") + scale_color_manual(values = hue_pal()(10), guide = "none")
dev.off()

subsample_df[,c("pool", "num")] = colsplit(subsample_df$subsample, "\\.", c("1", "2"))
subsample_df$pair = bb$pair[match(subsample_df$subsample, bb$subsample)]
subsample_df$num_nuc = as.numeric(subsample_df$num_nuc)
subsample_df$num = factor(subsample_df$num, levels = unique(subsample_df$num))
pdf("~/research/brain/results/bb_subsample_num_nuc.pdf", width = 3, height = 6)
ggplot(subsample_df, aes(x = num, y = pool, size = num_nuc, color = num)) + geom_point() + theme_classic() + xlab("") + ylab("") + scale_color_manual(values = c(pal(4)), guide = "none") + theme(axis.ticks.x = element_blank(), axis.line.x.bottom = element_blank(), axis.text.x = element_blank())
dev.off()

test = colorRampPalette(brewer.pal(8, "Set3"))(19)
pdf("~/research/brain/results/bb_subsample_num_nuc_2.pdf", width = 3, height = 6)
ggplot(subsample_df, aes(x = num, y = pool, size = num_nuc, color = pair)) + geom_point() + theme_classic() + xlab("") + ylab("") + scale_color_manual(values = test, guide = "none") + theme(axis.ticks.x = element_blank(), axis.line.x.bottom = element_blank(), axis.text.x = element_blank())
dev.off()

# local comment
