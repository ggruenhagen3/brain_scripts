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
convert53 = convert53[order(as.numeric(convert53$new.parent), decreasing = F),]
tmp = convert53[which(convert53$new == "8-9_Glut"),]
convert53 = convert53[-which(convert53$new == "8-9_Glut"),]
convert53 = rbind(convert53[1:29,], tmp, convert53[30:nrow(convert53),])

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



# pair1_bulk_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$sample == "b1")], cells.2 = colnames(bb)[which(bb$sample == "c1")])
# pair2_bulk_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$sample == "b2")], cells.2 = colnames(bb)[which(bb$sample == "c2")])
# pair3_bulk_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$sample == "b3")], cells.2 = colnames(bb)[which(bb$sample == "c3")])
# pair4_bulk_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$sample == "b4")], cells.2 = colnames(bb)[which(bb$sample == "c4")])
# pair5_bulk_stats = pct_dif_avg_logFC(bb, cells.1 = colnames(bb)[which(bb$sample == "b5")], cells.2 = colnames(bb)[which(bb$sample == "c5")])
# 
# bulk_common_genes = pair1_bulk_stats$genes[which( pair1_bulk_stats$genes %in% pair2_bulk_stats$genes & pair1_bulk_stats$genes %in% pair3_bulk_stats$genes & pair1_bulk_stats$genes %in% pair4_bulk_stats$genes & pair1_bulk_stats$genes %in% pair5_bulk_stats$genes )]
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

# *************************************************************************************************************
# Circular Figure =============================================================================================
# *************************************************************************************************************
set.seed(999)
n = 1000
df = data.frame(sectors = bb$seurat_clusters,
                x = bb$nFeature_RNA, y = bb$nCount_RNA)

library(circlize)
library(ComplexHeatmap)
circos.par("track.height" = 0.1)
circos.initialize(df$sectors, x = df$x)
circos.track(df$sectors, y = df$y,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(5), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })
# col = rep(c("#FF0000", "#00FF00"), 4)
# circos.trackPoints(df$sectors, df$x, df$y, col = col, pch = 16, cex = 0.5)
circos.trackPoints(df$sectors, df$x, df$y)
circos.text(-1, 0.5, "text", sector.index = "a", track.index = 1)
circos.clear()

mat1 = as.matrix(bb$nCount_RNA)
split = bb$seurat_clusters

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
# col_fun1 = colorRampPalette(viridis(100))
# col_fun1 = colorRamp2(c(min(bb$nCount_RNA), median(bb$nCount_RNA), max(bb$nCount_RNA)), c("blue", "white", "red"))
ii <- cut(bb$nCount_RNA, breaks = seq(min(bb$nCount_RNA), max(bb$nCount_RNA), len = 100), include.lowest = TRUE)
col_fun1 = colorRampPalette(viridis(100))(99)[ii]
col_fun1 = colorRamp2(quantile_breaks(bb$nCount_RNA, 100), colors = viridis(100))
circos.heatmap(mat1, split = split, col = col_fun1, show.sector.labels = TRUE)
circos.clear()

mat = bb@assays$RNA@counts
mat[which(mat > 1)] = 1

bb$annot15 = factor(convert15$new.full[match(bb$seurat_clusters, convert15$old)])
neuro_gen = read.csv("~/research/brain/data/neurogen_genes_final_050621.csv")[,1]
neuro_gen_scores = list()
for (cluster in convert15$new.full) {
  for (sample in unique(bb$sample)) {
    neuro_gen_scores[paste0(cluster, "_", sample)] = mean(mat[neuro_gen, which(bb$annot15 == cluster & bb$sample == sample)])
  }
}

bb$annot = factor(convert53$new[match(bb$seurat_clusters, convert53$old)], levels = convert53$new)
bb_df15 = data.frame(cluster = factor(convert15$new.full), col = convert15$col, new_sub = as.vector(table(factor(convert53$new.parent, levels = 1:15))))
# bb_df15 = rbind(bb_df15, data.frame(cluster = "16_NA", col = NA, new_sub = 3))
bb_df15$cluster = factor(bb_df15$cluster, levels = bb_df15$cluster)
  
bb_df53 = data.frame(cluster = levels(bb$annot), col = convert53$col, num_cells = aggregate(nCount_RNA ~ annot, bb@meta.data, length)[,2])

circos.par("gap.after" = c(rep(0, 14), 10), cell.padding = c(0, 0, 0, 0), "start.degree" = 90, track.margin = c(0.02, 0.02))
circos.initialize(bb_df15$cluster, xlim = cbind(rep(0, nrow(bb_df15)), bb_df15$new_sub))
track1_breaks = pretty(unlist(neuro_gen_scores), n = 2)
circos.track(ylim = c(min(track1_breaks), max(track1_breaks)), track.height = 0.30, bg.border = "black", panel.fun = function(x, y) {
  for (tb in track1_breaks) {
    if (tb != track1_breaks[1] & tb != track1_breaks[length(track1_breaks)])
      circos.segments(CELL_META$cell.xlim[1], tb, CELL_META$cell.xlim[2], tb, col = "gray60", lty = 2)
  }
  
  xrange = CELL_META$cell.xlim[2] - CELL_META$cell.xlim[1]
  this_gap = 0.15*xrange
  this_gap = ifelse(this_gap > 1.5, 1.5, this_gap)
  for (sample in unique(bb$sample)) {
    if ( startsWith(sample, "b") ) {
      this_col = "#d0000090"
      this_x = myJitter(CELL_META$xcenter-this_gap, 0.1*xrange)
    } else {
      this_col = "#023e8a90"
      this_x = myJitter(CELL_META$xcenter+this_gap, 0.1*xrange)
    }
    print(sample)
    this_score = as.numeric(neuro_gen_scores[paste0(CELL_META$sector.index, "_", sample)])
    print(this_score)
    circos.points(this_x, this_score, col = this_col, pch = 20, cex = 1.5)
  }
})
circos.yaxis(at = track1_breaks, sector.index = "1_Astro/MG", track.index = 1, side = "left")

circos.track(ylim = c(0, 1), track.height = 0.15, bg.border = NA, panel.fun = function(x, y) {
  pos = circlize:::polar2Cartesian(circlize(CELL_META$xcenter, CELL_META$ycenter))
  bg.col = bb_df15$col[CELL_META$sector.numeric.index]
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[1] - mm_y(2),
              CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE,
              adj = c(1, 0.5), cex = 0.6)
  circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1, col = bb_df15$col[CELL_META$sector.numeric.index], border = NA)
})


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


