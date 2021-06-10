library("scales")
cols15 = gc.ramp <- hue_pal()(15)

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

lightgray_hex = "#D3D3D3FF"
df$col15_blend = lightgray_hex
for (i in 0:14) {
  this_df = df[which(df$seuratclusters15 == i),]
  if (! is.na(this_df$gene[1])) {
    n = nrow(this_df)
    this_values = this_df$values
    m <- grDevices::colorRamp(c("lightgray", this_df$col15[1]))( (1:n)/n )
    this_cols = colour_values(this_df$value, palette = m)
    df[rownames(this_df), "col15_blend"] = this_cols
  }
}

my.pt.size = 0.8
p_list = list()
df = df[order(df$value, decreasing = F),]
p_list[["BHVE"]] = ggplot(df[which(df$cond == "BHVE"),], aes(UMAP_1, UMAP_2, col = col15_blend)) + geom_point(size = my.pt.size) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + scale_color_identity() + ggtitle("BHVE")
p_list[["CTRL"]] = ggplot(df[which(df$cond == "CTRL"),], aes(UMAP_1, UMAP_2, col = col15_blend)) + geom_point(size = my.pt.size) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + scale_color_identity() + ggtitle("CTRL")
plot_grid(plotlist=p_list, ncol = 2)