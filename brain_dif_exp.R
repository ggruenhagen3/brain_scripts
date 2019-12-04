# Install Packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("edgeR")
# BiocManager::install("Seurat")
# devtools::install_github("dynverse/dyno", force = TRUE)

# Load Packages
library("edgeR")
library("Seurat")
library("Matrix")
library("reticulate")
library("cowplot")
library("monocle")
library("monocle3")
# library(dyno)
# library(tidyverse)

rna_path <- "C:/Users/miles/Downloads/brain/"
rna_path <- "/nv/hp10/ggruenhagen3/scratch/brain/"

b1.data <- Read10X(data.dir = paste(rna_path, "data/BHVE-JTS03-B1/outs/filtered_feature_bc_matrix/", sep=""))
b2.data <- Read10X(data.dir = paste(rna_path, "data/BHVE-JTS02-B2/outs/filtered_feature_bc_matrix/", sep=""))
c1.data <- Read10X(data.dir = paste(rna_path, "data/CTRL-JTS03-C1/outs/filtered_feature_bc_matrix/", sep=""))

b1 <- CreateSeuratObject(counts = b1.data, project = "BHVE")
b2 <- CreateSeuratObject(counts = b2.data, project = "BHVE")
c1 <- CreateSeuratObject(counts = c1.data, project = "CTRL")

b1$cond <- "BHVE"
b2$cond <- "BHVE"
c1$cond <- "CTRL"

# ctrl <- RenameCells(ctrl, new.names = paste("CTRL_", colnames(ctrl@assays$RNA@counts), sep=""))
# clpp <- RenameCells(clpp, new.names = paste("CLIPP_", colnames(clpp@assays$RNA@counts), sep=""))
# 
# # Add metadata for cells expressing a gene
# gene_name <- "Celsr1"
# clpp.celsr1 <- clpp.data[which(row.names(clpp.data) == gene_name),]
# clpp.celsr1_not_0 <- clpp.celsr1[which(clpp.celsr1 != 0)]
# clpp$celsr1_label <- clpp.celsr1 != 0
# clpp$celsr1_group_label <- paste("clpp_", as.character(clpp$celsr1_label), sep="")
# 
# ctrl.celsr1 <- ctrl.data[which(row.names(ctrl.data) == gene_name),]
# ctrl.celsr1_not_0 <- ctrl.celsr1[which(ctrl.celsr1 != 0)]
# ctrl$celsr1_label <- ctrl.celsr1 != 0
# ctrl$celsr1_group_label <- paste("ctrl_", as.character(ctrl$celsr1_label), sep="")

# Normalize Data and then merge
b1 <- NormalizeData(b1, normalization.method = "LogNormalize", scale.factor = 100000)
b2 <- NormalizeData(b2, normalization.method = "LogNormalize", scale.factor = 100000)
c1 <- NormalizeData(c1, normalization.method = "LogNormalize", scale.factor = 100000)

combined <- merge(x=c1, y=c(b1,b2), merge.data = TRUE, add.cell.ids = c("CTRL", "BHVE", "BHVE"))
combined <- subset(combined, subset = nFeature_RNA > 500)
combined <- FindVariableFeatures(object = combined, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)

# Run the standard workflow for visualization and clustering
combined <- ScaleData(object = combined, vars.to.regress = NULL)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "umap", dims = 1:2)
combined <- FindClusters(combined, resolution = 0.28)
# Visualization
DimPlot(combined, reduction = "umap", split.by = "cond", label = TRUE)
png(filename = paste(rna_path, "results/umap.png", sep=""), width = 900, height = 500, unit="px")
print(p2)
dev.off()

## Monocle
# data <- as(as.matrix(combined@assays$RNA@data), 'sparseMatrix')
# pd <- new('AnnotatedDataFrame', data = combined@meta.data)
# fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
# fd <- new('AnnotatedDataFrame', data = fData)
# cds <- newCellDataSet(data, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
# cth <- newCellTypeHierarchy()
# cds <- classifyCells(cds, cth, 0.1)
# 
# cds <- estimateSizeFactors(cds)
# cds <- estimateDispersions(cds)
# disp_table <- dispersionTable(cds)
# ordering_genes <- subset(disp_table, mean_expression >= 0.1)
# cds <- setOrderingFilter(cds, ordering_genes)
# cds <- reduceDimension(cds)
# cds <- orderCells(cds)
# diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~cond")

## Monocle 3
clpp <- load_cellranger_data(paste(rna_path, "/CLIPP/", sep=""))
ctrl <- load_cellranger_data(paste(rna_path, "/CTRL/", sep=""))
cds <- combine_cds(list(clpp, ctrl))
cds <- preprocess_cds(cds, num_dim = 100)
# cds <- align_cds(cds, alignment_group = "batch")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds)

# Cluster analysis
# dist.matrix <- dist(x = Embeddings(object = combined[["pca"]])[, 0:8])
# s <- silhouette(as.vector(as.numeric(combined@meta.data$seurat_clusters)), dist = dist.matrix)
# plot(s)

# Autmatic Cell Type Annotation?
  # SCMAP was a bust because you need a reference dataset to project to.
  # I don't have a reference dataset.
  # combined.sce <- as.SingleCellExperiment(combined)
  # sce <- SingleCellExperiment(assays = list(counts = as.matrix(combined@assays$RNA@counts)), colData = combined@meta.data$seurat_clusters)
  # rowData(sce)$feature_symbol <- rownames(sce)
  # logcounts(sce) <- log2(counts(sce) + 1)
  # sce <- sce[!duplicated(rownames(sce)), ]
  # sce <- selectFeatures(sce, suppress_plot = FALSE)
  # sce@colData$cell_type1 <- as.numeric(as.vector(sce@colData$X))
  # sce <- indexCluster(sce)
  # results <- scmapCluster(sce, index_list = list(yan = metadata(sce)$scmap_cluster_index))
  ## Garnett was a bust because its pre-trained classifiers were only for lung or brain.

# Find clusters
## mCLIPP_cluster <- list()
## mCTRL_cluster <- list()
## Idents(object=combined) <- "seurat_clusters"
## colnames(mCLIPP) <- paste("CLIPP", colnames(mCLIPP), sep="_")
## colnames(mCTRL)  <- paste("CTRL",  colnames(mCTRL),  sep="_")
Idents(object = combined) <- "seurat_clusters"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
for (i in 0:num_clusters) {
  print(i)
  nk.markers <- FindMarkers(combined, ident.1 = i, verbose = FALSE)
  nk.markers$gene_name <- row.names(nk.markers)
  sig_nk.markers <- nk.markers[which(nk.markers$p_val_adj < 0.05 & abs(nk.markers$avg_logFC) > 2),]
  write.table(nk.markers, file = paste(rna_path, "/results/clusters/all_cluster_", i, ".tsv", sep=""), quote = FALSE, row.names = FALSE)
  write.table(sig_nk.markers, file = paste(rna_path, "/results/clusters/sig_cluster_", i, ".tsv", sep=""), quote = FALSE, row.names = FALSE)
  write.table(sig_nk.markers$gene_name, file = paste(rna_path, "/results/clusters/genes_sig_cluster_", i, ".tsv", sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE)
  sig_nk_pos.markers <- nk.markers[which(nk.markers$p_val_adj < 0.05 & nk.markers$avg_logFC > 2),]
  write.table(sig_nk_pos.markers, file = paste(rna_path, "/results/clusters/sig_pos_cluster_", i, ".tsv", sep=""), quote = FALSE, row.names = FALSE)
  write.table(sig_nk_pos.markers$gene_name, file = paste(rna_path, "/results/clusters/genes_sig_pos_cluster_", i, ".tsv", sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE)
}
print("Done finding clusters")

## EdgeR DEG
# names <- colnames(combined@assays$RNA@counts)
# group <- c(rep("CTRL", length(grep("CTRL*", names))), rep("CLIPP", length(grep("CLIPP*", names))))
# y <- DGEList(counts = combined@assays$RNA@counts, group = group)
# y <- calcNormFactors(y)
# # design <- model.matrix(~group) 
# design <- model.matrix(~0+group, data=y$samples)
# colnames(design) <- levels(y$samples$group)
# y <- estimateDisp(y,design)
# fit <- glmQLFit(y,design) 
# qlf <- glmQLFTest(fit,coef=2) 
# topTags(qlf)
# 
# my.contrasts <- makeContrasts(CLIPPvsCTRL=CLIPP-CTRL, levels=design)
# lrt.AvsB <- glmQLFTest(fit, contrast=my.contrasts[,"CLIPPvsCTRL"])
# topTags(lrt.AvsB)
# results_edgeR_AvsB <- topTags(lrt.AvsB, n = nrow(combined@assays$RNA@counts), sort.by = "none")
# DEsum = sum(results_edgeR_AvsB$table$FDR < .05)
# R = results_edgeR_AvsB$table[results_edgeR_AvsB$table$FDR < .05,]
# R <- R[order(R$FDR),]
# write.csv(R, file = "Differential_expression_CLIPP_CTRL.csv")

# Paint specific genes
genes <- c("Nes", "Emcn", "Mki67", "Aurkb", "Cspg4", "Sox10", "Ptprc", "Thy1", "Celsr1", "Mcam", "Rspo4", "Pax9", "Krt14")
for (gene in genes) {
  png(filename = paste(rna_path, "results/", gene, "_umap.png", sep=""), width = 900, height = 500, unit="px")
  p <- FeaturePlot(combined, features = c(gene), split.by = "cond", reduction = "umap", pt.size = 2, label=TRUE, order = TRUE)
  val <- FeaturePlot(combined, features = c("Cspg4"), reduction = "umap", label=TRUE, order = TRUE)
  print(p)
  dev.off()
}

# Plot Val's combo of genes
png(filename = paste(rna_path, "results/", "Mki67_Aurkb_umap.png", sep=""), width = 1200, height = 600, unit="px")
p3<-FeaturePlot(combined, features = c("Mki67", "Aurkb"), split.by = "cond", reduction = "umap", pt.size = 1, label=TRUE, order = TRUE)
print(p3)
dev.off()
png(filename = paste(rna_path, "results/", "Ng2_CD146_umap.png", sep=""), width = 1200, height = 600, unit="px")
FeaturePlot(combined, features = c("Cspg4", "Mcam"), split.by = "cond", reduction = "umap", pt.size = 1,label=TRUE, order = TRUE)
dev.off()
png(filename = paste(rna_path, "results/", "Rspo4_Pax9_umap.png", sep=""), width = 1200, height = 600, unit="px")
FeaturePlot(combined, features = c("Rspo4", "Pax9"), split.by = "cond", reduction = "umap", pt.size = 1,label=TRUE, order = TRUE)
dev.off()
png(filename = paste(rna_path, "results/", "Sox10_CD45_umap.png", sep=""), width = 1200, height = 600, unit="px")
FeaturePlot(combined, features = c("Sox10", "Ptprc"), split.by = "cond", reduction = "umap", pt.size = 1,label=TRUE, order = TRUE)
dev.off()

# Find Dif Exp Genes Across Conditions
# combined$celltype.cond <- paste(Idents(combined), combined$cond, sep = "_")
# combined$celltype <- Idents(combined)
# Idents(combined) <- "celltype.cond"
# dif_exp <- data.frame()
# for (i in 0:20) {
#   print(i)
#   dif_exp_i <- FindMarkers(combined, ident.1 = paste(i, "CTRL", sep="_"), ident.2 = paste(i, "CLIPP", sep="_"), verbose = FALSE)
#   dif_exp_i <- dif_exp_i[which(dif_exp_i$p_val_adj < 0.05 | abs(dif_exp_i$avg_logFC) > 2),]
#   dif_exp_i$cluster <- rep(i, nrow(dif_exp_i))
#   dif_exp <- rbind(dif_exp, dif_exp_i)
# }
# Idents(object = combined) <- "orig.ident"
# dif_exp <- FindMarkers(object = combined, ident.1 = "CTRL", ident.2 = "CLIPP", group.by = "orig.ident")
# dif_exp$gene_name <- row.names(dif_exp)
# dif_exp <- dif_exp[which(dif_exp$p_val_adj < 0.05),]
# write.table(dif_exp, file = paste(rna_path, "results/dif_exp_cond.tsv", sep=""), quote = FALSE, row.names = FALSE)
# 
# # Plot Top Features
# for (i in 1:20) {
#   feature <- top_features[i]
#   png(filename = paste(rna_path, "plots/", feature, "_tsne_norm.png", sep=""), width = 700, height = 500, unit="px")
#   p <- FeaturePlot(combined, features = feature, min.cutoff = "q9", label=TRUE)
#   print(p)
#   dev.off()
# }
# 
# Celsr1 Analysis
# gene_name <- "Celsr1"
# clpp.celsr1 <- clpp.data[which(row.names(clpp.data) == gene_name),]
# clpp.celsr1_not_0 <- clpp.celsr1[which(clpp.celsr1 != 0)]
# clpp.celsr1_label <- clpp.celsr1 != 0
# 
# ctrl.celsr1 <- ctrl.data[which(row.names(ctrl.data) == gene_name),]
# ctrl.celsr1_not_0 <- ctrl.celsr1[which(ctrl.celsr1 != 0)]
# clpp.celsr1_label <- clpp.celsr1 != 0
# 
# rna <- combined[['RNA']]
# rna <- rna@counts
# celsr1 <- rna[which(row.names(rna) == gene_name),]
# celsr1_not_0 <- celsr1[which(celsr1 != 0)]
# 
# cat(paste("CLPP", gene_name, ":"), length(clpp.celsr1_not_0), "\n")
# cat(paste("CLPP", gene_name, ">1:"), length(clpp.celsr1_not_0[which(clpp.celsr1_not_0 != 1)]), "\n")
# cat(paste("CTRL", gene_name, ":"), length(ctrl.celsr1_not_0), "\n")
# cat(paste("CTRL", gene_name, ">1:"), length(ctrl.celsr1_not_0[which(ctrl.celsr1_not_0 != 1)]), "\n")
# cat(paste("Combined", gene_name, ":"), length(celsr1_not_0), "\n")
# cat(paste("Combined", gene_name, ">1:"), length(celsr1_not_0[which(celsr1_not_0 != 1)]), "\n")

# # Differential Expression in Cells expressing a gene
# gene_name <- "Celsr1"
# Idents(object = combined) <- "celsr1_label"
# dif_exp_gene <- FindMarkers(object = combined, ident.1 = "TRUE", ident.2 = "FALSE")
# dif_exp_gene$gene_name <- row.names(dif_exp_gene)
# dif_exp_gene <- dif_exp_gene[which(dif_exp_gene$p_val_adj < 0.05 & abs(dif_exp_gene$avg_logFC) > 2),]
# write.table(dif_exp_gene, file = paste(rna_path, "/results/dif_exp_celsr1_1.tsv", sep=""), quote = FALSE, row.names = FALSE)
# 
# Idents(object = combined) <- "celsr1_group_label"
# dif_gene <- FindMarkers(object = combined, ident.1 = "clpp_TRUE", ident.2 = "ctrl_TRUE")
# sig_gene <- dif_gene[which(dif_gene$p_val_adj < 0.05),]
# 
# write.table(sig_clpp, file = paste(rna_path, "results/dif_exp_celsr_clpp.tsv", sep=""), quote = FALSE, row.names = FALSE)
# write.table(sig_ctrl, file = paste(rna_path, "results/dif_exp_celsr_ctrl.tsv", sep=""), quote = FALSE, row.names = FALSE)

# Exploration of PCA
# pca <- combined@reductions$pca
# eigValues = (pca@stdev)^2  ## EigenValues
# varExplained = eigValues / sum(eigValues)
# library(reshape2)
# varExplained$row <- seq_len(1)
# dat2 <- melt(varExplained, id.vars = "row")
# library(ggplot2)
# numbers <- 1:20
# ggplot(dat2, aes(x = 1, y = varExplained, fill = numbers)) + 
#   geom_bar(stat = "identity") +
#   xlab("\nType") +
#   ylab("Time\n") +
#   guides(fill = FALSE) +
#   theme_bw()
# 
# charts.data <- data.frame(percentage <- varExplained[1:5], year <- rep(1, 5), dim <- as.character(1:5))
# p <- ggplot() + geom_bar(aes(y = percentage, x = year, fill = dim), data = charts.data,
#                           stat="identity") + ggtitle("Variance explained by first 5 dimenstions") + ylim(0,1) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
# p + scale_fill_manual(values=wes_palette(n=5, name="Rushmore"))
# 
# cells <- pca@cell.embeddings
# weird <- cells[which(cells[,1] > 3 & cells[,2] < -3),]
# main <- cells[which(!cells[,1] > 3 & !cells[,2] < -3),]
# weird_names <- rownames(weird)
# main_names <- rownames(main)
# weird_counts <- counts[,weird_names]
# rna <- FetchData(combined, vars="nCount_RNA")
# feature <- FetchData(combined, vars="nFeature_RNA")
# avg_rna <- mean(rna$nCount_RNA)
# avg_feature <- mean(feature$nFeature_RNA)
# avg_weird_rna <- mean(rna[weird_names,])
# avg_weird_feature <- mean(feature[weird_names,])
# 
# weird_obj <- subset(combined, cells = weird_names)
# main_obj <- subset(combined, cells = main_names)
# main_obj <- RunTSNE(main_obj, reduction = "pca", dims = 1:15)
# main_obj <- FindNeighbors(main_obj, reduction = "pca", dims = 1:15)
# main_obj <- FindClusters(main_obj, resolution = 0.7)
# DimPlot(main_obj, reduction = "tsne", split.by = "cond", label = TRUE)
