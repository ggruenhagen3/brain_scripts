rna_path = "C:/Users/miles/Downloads/brain/"
# rna_path = "~/research/brain/"
bb = readRDS(paste0(rna_path, "data/bb_subsample_02222021.RDS"))
source(paste0(rna_path, "brain_scripts/all_f.R"))

#=====================================================================================
# Figure 1 ===========================================================================
#=====================================================================================

# *** 1B. 15 cluster and 53 cluster level UMAP *** #
# Define 15 level colors
library("scales")
cols15 = gc.ramp <- hue_pal()(15)
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
convert15 = convert15[order(as.numeric(convert15$new.num), decreasing = F),]
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

bb$good_names = convert53$new[match(bb$seuratclusters53, convert53$old)]
bb$good_names_num = colsplit(bb$good_names, "_", c("num", "ex"))[,1]
Idents(bb) = bb$good_names_num
clust53_new_col_list2 = convert53$col
names(clust53_new_col_list2) = colsplit(convert53$new, "_", c("num", "ex"))[,1]
DimPlot(bb, label = T, pt.size = 1) + scale_color_manual(values = clust53_new_col_list2) + NoLegend()

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
