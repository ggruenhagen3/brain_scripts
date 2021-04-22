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
