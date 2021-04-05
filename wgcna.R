
r_mat = readRDS("~/scratch/brain/data/mat_data_cor.RDS")
r_mat =  r_mat[which( ! is.na(r_mat[2,]) ), which( ! is.na(r_mat[2,]) )]  # removes NAs

# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# sft = pickSoftThreshold(r_mat, powerVector = powers, verbose = 5)
# softPower = 6 # because it's the lowest at which it hits 90%

# From Simulated-00-Background.pdf
# "Optionally, the user can also specify a signed co-expression network where the
# adjacency is defined as..."

library("WGCNA")
# r_mat = abs(r_mat) ** softPower # this converts the correlation matrix to an adjacency matrix
# TOM = TOMsimilarity(r_mat, TOMType = "signed")
adj_mat = adjacency.fromSimilarity(r_mat)
TOM = TOMsimilarity(adj_mat)
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
png_name = "~/scratch/brain/results/wgcna.png"
png(png_name, width = 750, height = 750, res = 90)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04);
dev.off()

# Find Modules
minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

df = data.frame(gene = rownames(r_mat), module = dynamicMods)
df = df[which(df$module > 0),]
df[which(df$gene == "egr1"),]
write.csv(df, "~/scratch/brain/results/bb_wgcna_modules.csv")
write.csv(df[which(df$module == 25),], "~/scratch/brain/results/bb_wgcna_modules_25_ieg.csv")

# Split IEG Module into those that are positively correlated and those that are negative
df_ieg = read.csv("~/scratch/brain/results/bb_wgcna_modules_25_ieg.csv")

df_ieg$egr1_cor = r_mat["egr1", df_ieg$gene]
df_ieg$npas4_cor = r_mat["npas4", df_ieg$gene]
df_ieg$homer1_cor = r_mat["homer1", df_ieg$gene]
df_ieg$LOC101487312_cor = r_mat["LOC101487312", df_ieg$gene]

df_ieg_pos = sapply(1:nrow(df_ieg), function(x) all(df_ieg[x,4:ncol(df_ieg)] > 0) )
df_ieg_pos = df_ieg[which(df_ieg_pos),]
write.csv(df_ieg_pos, "~/scratch/brain/results/bb_wgcna_modules_25_ieg_pos.csv")

rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"), width = 800, height = 800)
png("~/scratch/brain/results/bb_wgcna_modules_ieg_03232021.png")
heatmap.2(r_mat[df_ieg_pos$gene, df_ieg_pos$gene], scale = "none", trace = "none")
dev.off()

# This is meant to merge similar modules. This changed nothing for me.
# dynamicColors = labels2colors(dynamicMods[which(dynamicMods > 0)])
# my_data = t(as.matrix(bb@assays$RNA@data[which(dynamicMods > 0),]))
# MEList = moduleEigengenes(my_data, colors = dynamicColors)
# MEs = MEList$eigengenes
# MEDiss = 1-cor(MEs)
# METree = hclust(as.dist(MEDiss), method = "average")
# MEDissThres = 0.25
# merge = mergeCloseModules(my_data, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# mergedColors = merge$colors
# colorOrder = c("grey", standardColors(50))
# moduleLabels = match(mergedColors, colorOrder)-1

df6 = read.csv("C:/Users/miles/Downloads/bb_wgcna_modules.csv", stringsAsFactors = F)
for (module in sort(unique(df6$module))) {
  print(module)
  this_genes = df6$gene[which(df6$module == module)]
  bb$this_module = colSums(bb@assays$RNA@data[this_genes,])
  png(paste0("C:/Users/miles/Downloads/brain/results/bb/modules6/module_", module, ".png"), height = 550, width = 600)
  print(FeaturePlot(bb, "this_module", order = T, label = T, pt.size = 1.5) + ggtitle(paste("Module", module)))
  dev.off()
}
