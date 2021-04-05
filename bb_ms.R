rna_path = "C:/Users/miles/Downloads/brain/"
# rna_path = "~/scratch/brain/"
bb = readRDS(paste0(rna_path, "data/bb_subsample_02222021.RDS"))
source(paste0(rna_path, "brain_scripts/all_f.R"))

#=====================================================================================
# Figure 1 ===========================================================================
#=====================================================================================

# *** 1B. 15 cluster and 53 cluster level UMAP *** #
# Define 15 level colors
library("scales")
cols15 = gc.ramp <- hue_pal()(15)
names(cols15) = 0:14

# Define Colors for Mixing
library("gplots")
mix_cols = c("gray85", "gray80", "gray75", "gray70", "gray65", "gray60", "gray55", "gray50", "gray45", "gray40", "gray35", "gray30", "gray25", "gray20")
mix_cols = col2hex(mix_cols)
mix_pal = colorRampPalette(mix_cols)
  
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
convert15 = separate(data = convert15, col = new, into = c("new.num", "new.junk"), sep = "_")
# Get the parent and color of the 53 cluster 
convert53$new.parent.old = convert15$old[match(convert53$new.parent, convert15$new.num)]
convert53$new.parent.old.col = cols15[as.character(convert53$new.parent.old)]
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

clust53_new_col_list2 = convert53$col
names(clust53_new_col_list2) = convert53$new
DimPlot(bb, label = T, pt.size = 1) + scale_color_manual(values = clust53_new_col_list2) + NoLegend()

pdf("C:/Users/miles/Downloads/brain/results/bb/53_clusters_on_15_light_dark_label_04012021.pdf", width = 12, height = 12)
print(DimPlot(bb, label = T, pt.size = 1) + scale_color_manual(values = clust53_new_col_list2) + NoLegend())
dev.off()

pdf("C:/Users/miles/Downloads/brain/results/bb/53_clusters_on_15_light_dark_no_label_04012021.pdf", width = 12, height = 12)
print(DimPlot(bb, label = F, pt.size = 1) + scale_color_manual(values = clust53_new_col_list2) + NoLegend())
dev.off()

