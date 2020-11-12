# Install Packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("edgeR")
# BiocManager::install("Seurat")
# devtools::install_github("dynverse/dyno", force = TRUE)
args = commandArgs(trailingOnly=TRUE)

# Load Packages
library("edgeR")
library("Seurat")
library("Matrix")
# if(.Platform$OS.type == "windows") Sys.setenv(PATH= paste("C:/Anaconda3/Library/bin",Sys.getenv()["PATH"],sep=";"))
library("reticulate")
library("stringr")
library("dplyr")
# use_python("C:/Users/miles/AppData/Local/Programs/Python/Python37/", required = TRUE)
library("cowplot")
# library("monocle")
# library("monocle3")
# library(dyno)
# library(tidyverse)

rna_path <- "C:/Users/miles/Downloads/brain/"
# rna_path <- "/nv/hp10/ggruenhagen3/scratch/brain/"

b1.data <- Read10X(data.dir = paste(rna_path, "data/BHVE-JTS03-B1-from-bcl/outs/filtered_feature_bc_matrix/", sep=""))
b2.data <- Read10X(data.dir = paste(rna_path, "data/BHVE-JTS02-B2-from-bcl/outs/filtered_feature_bc_matrix/", sep=""))
c1.data <- Read10X(data.dir = paste(rna_path, "data/CTRL-JTS03-C1-from-bcl/outs/filtered_feature_bc_matrix/", sep=""))

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
combined <- RunPCA(combined, npcs = 30, verbose = FALSE, features = combined@assays$RNA@var.features)
# t-SNE and Clustering
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "umap", dims = 1:2)
if(length(args) > 0) {
  resolution <- args[1]
} else {
  resolution <- 0.28
}
combined <- FindClusters(combined, resolution = 0.28)
# Visualization
p2 <- DimPlot(combined, reduction = "umap", split.by = "cond", label = TRUE)
png(filename = paste(rna_path, "results/umap.png", sep=""), width = 900, height = 500, unit="px")
print(p2)
dev.off()

## Monocle 3
library(org.Dr.eg.db)
library(org.Mm.eg.db)
library(garnett)
classifier <- readRDS(paste(rna_path, "/data/mmBrain_20191017.RDS", sep=""))
b1 <- load_cellranger_data(paste(rna_path, "/data/BHVE-JTS03-B1-from-bcl/", sep=""))
b2 <- load_cellranger_data(paste(rna_path, "/data/BHVE-JTS02-B2-from-bcl/", sep=""))
c1 <- load_cellranger_data(paste(rna_path, "/data/CTRL-JTS03-C1-from-bcl/", sep=""))
cds <- combine_cds(list(b1, b2, c1))
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)
plot_cells(cds)
classify_cells(cds, classifier, cluster_extend = TRUE, db = org.Dr.eg.db, cds_gene_id_type = "org.Dr.egZFIN")
classify_cells(cds, classifier, cluster_extend = TRUE, db = org.Mm.eg.db, cds_gene_id_type = "SYMBOL")

output<- evaluate("test(1)")

backup_options <- options()
# options(show.error.messages= FALSE)
# options(error=traceback)
# options(error = expression(NULL))
# options()

test <- mzebraToMouseInPlace2(combined)
seurat_object <- test
data <- as(as.matrix(seurat_object@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat_object@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
cds@phenoData@data$Size_Factor[1:length(cds@phenoData@data$Size_Factor)] <- 0

# cds <- combine_cds(list(clpp, ctrl))
# cds <- preprocess_cds(cds, num_dim = 100)
# # cds <- align_cds(cds, alignment_group = "batch")
# cds <- reduce_dimension(cds)
# cds <- cluster_cells(cds)
# cds <- learn_graph(cds)
# cds <- order_cells(cds)
# plot_cells(cds)

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
print("Finding DEG between clusters")
Idents(object = combined) <- "seurat_clusters"
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
for (i in 4:num_clusters) {
  print(i)
  nk.markers <- FindMarkers(combined, ident.1 = i, verbose = FALSE)
  nk.markers$gene_name <- row.names(nk.markers)
  sig_nk.markers <- nk.markers[which(nk.markers$p_val_adj < 0.05 & abs(nk.markers$avg_logFC) > 2),]
  write.table(nk.markers, file = paste(rna_path, "/results/clusters/", num_clusters+1, "/all_cluster_", i, ".tsv", sep=""), quote = FALSE, row.names = FALSE)
  write.table(sig_nk.markers, file = paste(rna_path, "/results/clusters/", num_clusters+1, "/sig_cluster_", i, ".tsv", sep=""), quote = FALSE, row.names = FALSE)
  write.table(sig_nk.markers$gene_name, file = paste(rna_path, "/results/clusters/", num_clusters+1, "/genes_sig_cluster_", i, ".tsv", sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE)
  sig_nk_pos.markers <- nk.markers[which(nk.markers$p_val_adj < 0.05 & nk.markers$avg_logFC > 2),]
  write.table(sig_nk_pos.markers, file = paste(rna_path, "/results/clusters/", num_clusters+1, "/sig_pos_cluster_", i, ".tsv", sep=""), quote = FALSE, row.names = FALSE)
  write.table(sig_nk_pos.markers$gene_name, file = paste(rna_path, "/results/clusters/", num_clusters+1, "/genes_sig_pos_cluster_", i, ".tsv", sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE)
}
print("Done finding clusters")

# DEG between conditions/samples
Idents(combined) <- "cond"
combined@active.assay <- "RNA"
b_v_sc_b <- FindMarkers(combined, ident.1 = "BHVE", ident.2 = "CTRL", only.pos = TRUE)
b_v_sc_b <- b_v_sc_b[which(b_v_sc_b$p_val_adj < 0.05),]
b_v_sc_sc <- FindMarkers(combined, ident.1 = "CTRL", ident.2 = "BHVE", only.pos = TRUE)
b_v_sc_sc <- b_v_sc_sc[which(b_v_sc_sc$p_val_adj < 0.05),]

write.table(b_v_sc_b, paste(rna_path, "/results/bhve_v_sand_ctrl_deg_bhve.tsv", sep=""), sep="\t", quote=FALSE)
write.table(b_v_sc_sc, paste(rna_path, "/results/bhve_v_sand_ctrl_deg_sand_ctrl.tsv", sep=""), sep="\t", quote=FALSE)

# Add MZ Data to the Object
mz <- NormalizeData(mz, normalization.method = "LogNormalize", scale.factor = 1000000)
combined_mz <- merge(combined, mz, merge.data = TRUE)
# anchors <- FindIntegrationAnchors(object.list = c(combined, mz), dims = 1:50)
# combined_mz <- IntegrateData(anchorset = anchors, dims = 1:50)
# combined_mz <- FindVariableFeatures(object = combined_mz, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
combined_mz <- ScaleData(object = combined_mz, features = combined_mz@assays$RNA@var.features)
# combined_mz <- RunPCA(combined_mz, npcs = 30, verbose = FALSE, features = combined_mz@assays$RNA@var.features)
# combined_mz <- RunUMAP(combined_mz, reduction = "pca", dims = 1:50)
# combined_mz <- FindNeighbors(combined_mz, reduction = "umap", dims = 1:2)
# combined_mz <- FindClusters(combined_mz, resolution = 0.28)
# DimPlot(combined_mz, reduction = "umap", split.by = "cond", label = TRUE)
# p2 <- DimPlot(combined, reduction = "umap", split.by = "cond", label = TRUE)

# All Behave vs Rock Control
Idents(combined_mz) <- "cond"
combined_mz@active.assay <- "RNA"
b_v_rc_b <- FindMarkers(combined_mz, ident.1 = "BHVE", ident.2 = "MZ", only.pos = TRUE)
b_v_rc_b <- b_v_rc_b[which(b_v_rc_b$p_val_adj < 0.05),]
b_v_rc_rc <- FindMarkers(combined_mz, ident.1 = "MZ", ident.2 = "BHVE", logfc.threshold=0.01, only.pos = TRUE)
b_v_rc_rc <- b_v_rc_rc[which(b_v_rc_rc$p_val_adj < 0.05),]

write.table(b_v_rc_b, paste(rna_path, "/results/bhve_v_rock_ctrl_deg_bhve.tsv", sep=""), sep="\t", quote=FALSE)
write.table(b_v_rc_rc, paste(rna_path, "/results/bhve_v_rock_ctrl_deg_rock_ctrl.tsv", sep=""), sep="\t", quote=FALSE)

# All Behave vs All Control
combined_mz$isBHVE <- combined_mz$cond
combined_mz$isBHVE[which(combined_mz$isBHVE == "MZ")] <- "CTRL"
Idents(combined_mz) <- "isBHVE"
b_v_ac_b <- FindMarkers(combined_mz, ident.1 = "BHVE", ident.2 = "CTRL", only.pos = TRUE)
b_v_ac_b <- b_v_ac_b[which(b_v_ac_b$p_val_adj < 0.05),]
b_v_ac_ac <- FindMarkers(combined_mz, ident.1 = "CTRL", ident.2 = "BHVE", only.pos = TRUE)
b_v_ac_ac <- b_v_ac_ac[which(b_v_ac_ac$p_val_adj < 0.05),]

write.table(b_v_ac_b, paste(rna_path, "/results/bhve_v_all_ctrl_deg_bhve.tsv", sep=""), sep="\t", quote=FALSE)
write.table(b_v_ac_ac, paste(rna_path, "/results/bhve_v_all_ctrl_deg_all_ctrl.tsv", sep=""), sep="\t", quote=FALSE)

# Behave 1 vs All Control
combined_mz$isBHVE1 <- combined_mz$sample
combined_mz$isBHVE1[which(combined_mz$isBHVE1 == "MZ" | combined_mz$isBHVE1 == "c1")] <- "CTRL"
Idents(combined_mz) <- "isBHVE1"
b1_v_ac_b1 <- FindMarkers(combined_mz, ident.1 = "b1", ident.2 = "CTRL", only.pos = TRUE)
b1_v_ac_b1 <- b1_v_ac_b1[which(b1_v_ac_b1$p_val_adj < 0.05),]
b1_v_ac_ac <- FindMarkers(combined_mz, ident.1 = "CTRL", ident.2 = "b1", only.pos = TRUE)
b1_v_ac_ac <- b1_v_ac_ac[which(b1_v_ac_ac$p_val_adj < 0.05),]

write.table(b1_v_ac_b1, paste(rna_path, "/results/bhve1_v_all_ctrl_deg_bhve1.tsv", sep=""), sep="\t", quote=FALSE)
write.table(b1_v_ac_ac, paste(rna_path, "/results/bhve1_v_all_ctrl_deg_all_ctrl.tsv", sep=""), sep="\t", quote=FALSE)

# Behave 2 vs All Control
b2_v_ac_b2 <- FindMarkers(combined_mz, ident.1 = "b2", ident.2 = "CTRL", only.pos = TRUE)
b2_v_ac_b2 <- b2_v_ac_b2[which(b2_v_ac_b2$p_val_adj < 0.05),]
b2_v_ac_ac <- FindMarkers(combined_mz, ident.1 = "CTRL", ident.2 = "b2", only.pos = TRUE)
b2_v_ac_ac <- b2_v_ac_ac[which(b2_v_ac_ac$p_val_adj < 0.05),]

write.table(b2_v_ac_b2, paste(rna_path, "/results/bhve2_v_all_ctrl_deg_bhve2.tsv", sep=""), sep="\t", quote=FALSE)
write.table(b2_v_ac_ac, paste(rna_path, "/results/bhve2_v_all_ctrl_deg_all_ctrl.tsv", sep=""), sep="\t", quote=FALSE)

# Behave 1 vs Sand Control
Idents(combined_mz) <- "sample"
b1_v_sc_b1 <- FindMarkers(combined_mz, ident.1 = "b1", ident.2 = "c1", only.pos = TRUE)
b1_v_sc_b1 <- b1_v_sc_b1[which(b1_v_sc_b1$p_val_adj < 0.05),]
b1_v_sc_sc <- FindMarkers(combined_mz, ident.1 = "c1", ident.2 = "b1", only.pos = TRUE)
b1_v_sc_sc <- b1_v_sc_sc[which(b1_v_sc_sc$p_val_adj < 0.05),]

write.table(b1_v_sc_b1, paste(rna_path, "/results/bhve1_v_sand_ctrl_deg_bhve1.tsv", sep=""), sep="\t", quote=FALSE)
write.table(b1_v_sc_sc, paste(rna_path, "/results/bhve1_v_sand_ctrl_deg_sand_ctrl.tsv", sep=""), sep="\t", quote=FALSE)

# Behave 2 vs Sand Control
Idents(combined_mz) <- "sample"
b2_v_sc_b2 <- FindMarkers(combined_mz, ident.1 = "b2", ident.2 = "c1", only.pos = TRUE)
b2_v_sc_b2 <- b2_v_sc_b2[which(b2_v_sc_b2$p_val_adj < 0.05),]
b2_v_sc_sc <- FindMarkers(combined_mz, ident.1 = "c1", ident.2 = "b2", only.pos = TRUE)
b2_v_sc_sc <- b2_v_sc_sc[which(b2_v_sc_sc$p_val_adj < 0.05),]

write.table(b2_v_sc_b2, paste(rna_path, "/results/bhve2_v_sand_ctrl_deg_bhve2.tsv", sep=""), sep="\t", quote=FALSE)
write.table(b2_v_sc_sc, paste(rna_path, "/results/bhve2_v_sand_ctrl_deg_sand_ctrl.tsv", sep=""), sep="\t", quote=FALSE)

# Behave 1 vs Sand Control
Idents(combined_mz) <- "sample"
b1_v_rc_b1 <- FindMarkers(combined_mz, ident.1 = "b1", ident.2 = "MZ", only.pos = TRUE)
b1_v_rc_b1 <- b1_v_rc_b1[which(b1_v_rc_b1$p_val_adj < 0.05),]
b1_v_rc_rc <- FindMarkers(combined_mz, ident.1 = "MZ", ident.2 = "b1", only.pos = TRUE)
b1_v_rc_rc <- b1_v_rc_rc[which(b1_v_rc_rc$p_val_adj < 0.05),]

write.table(b1_v_rc_b1, paste(rna_path, "/results/bhve1_v_rock_ctrl_deg_bhve1.tsv", sep=""), sep="\t", quote=FALSE)
write.table(b1_v_rc_rc, paste(rna_path, "/results/bhve1_v_rock_ctrl_deg_rock_ctrl.tsv", sep=""), sep="\t", quote=FALSE)

# Behave 2 vs Sand Control
Idents(combined_mz) <- "sample"
b2_v_rc_b2 <- FindMarkers(combined_mz, ident.1 = "b2", ident.2 = "MZ", only.pos = TRUE)
b2_v_rc_b2 <- b2_v_rc_b2[which(b2_v_rc_b2$p_val_adj < 0.05),]
b2_v_rc_rc <- FindMarkers(combined_mz, ident.1 = "MZ", ident.2 = "b2", only.pos = TRUE)
b2_v_rc_rc <- b2_v_rc_rc[which(b2_v_rc_rc$p_val_adj < 0.05),]

write.table(b2_v_rc_b2, paste(rna_path, "/results/bhve2_v_rock_ctrl_deg_bhve2.tsv", sep=""), sep="\t", quote=FALSE)
write.table(b2_v_rc_rc, paste(rna_path, "/results/bhve2_v_rock_ctrl_deg_rock_ctrl.tsv", sep=""), sep="\t", quote=FALSE)

# Behave 1 vs Behave 2
Idents(combined_mz) <- "sample"
b1_v_b2_b1 <- FindMarkers(combined_mz, ident.1 = "b1", ident.2 = "b2", only.pos = TRUE)
b1_v_b2_b1 <- b1_v_b2_b1[which(b1_v_b2_b1$p_val_adj < 0.05),]
b1_v_b2_b2 <- FindMarkers(combined_mz, ident.1 = "b2", ident.2 = "b1", only.pos = TRUE)
b1_v_b2_b2 <- b1_v_b2_b2[which(b1_v_b2_b2$p_val_adj < 0.05),]

write.table(b1_v_b2_b1, paste(rna_path, "/results/bhve1_v_bhve2_deg_bhve1.tsv", sep=""), sep="\t", quote=FALSE)
write.table(b1_v_b2_b2, paste(rna_path, "/results/bhve1_v_bhve2_deg_bhve2.tsv", sep=""), sep="\t", quote=FALSE)

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
marker_path <- paste(rna_path, "data/markers/", sep="")
marker_files <- dir(marker_path, pattern =paste("*.txt", sep=""))
markers <- data.frame(gene <- c(), bio <- c())
for (i in 1:length(marker_files)) {
  file <- read.table(paste(marker_path, marker_files[i], sep=""), header=FALSE, sep="\t", stringsAsFactors=FALSE)
  # file[,1] <- str_to_title(file[,1])
  file[,1] <- tolower(file[,1])
  markers <- rbind(markers, file[,1:2])
}
colnames(markers) <- c("gene", "bio")
markers <- markers[which(markers$bio != "LG11_HIGH_FST"),]

for (i in 1:nrow(markers)) {
  gene <- markers[i,1]
  gene_lower <- tolower(gene)
  gene_upper <- toupper(gene)
  gene_title <- str_to_title(gene)
  bio <- markers[i,2]
  dir.create(paste(rna_path, "results/painting/", sep=""), showWarnings = FALSE)
  dir.create(paste(rna_path, "results/painting/", bio, sep=""), showWarnings = FALSE)
  
  if (gene_lower %in% rownames(combined@assays$RNA)) {
    gene <- gene_lower
  } else if (gene_upper %in% rownames(combined@assays$RNA)) {
    gene <- gene_upper
  } else if (gene_title %in% rownames(combined@assays$RNA)) {
    gene <- gene_title
  } else {
    gene <- paste("Did not find gene:", gene)
  }
  
  if (startsWith(gene, "Did not find gene:")) {
    print(gene)
  } else {
    print(gene)
    png(filename = paste(rna_path, "results/painting/", bio, "/", gene, "_umap.png", sep=""), width = 900, height = 500, unit="px")
    p <- FeaturePlot(combined, features = c(gene), split.by = "cond", reduction = "umap", pt.size = 2, label=TRUE, order = TRUE)
    print(p)
    dev.off()
  }
}

gene1 <- "sox9a"
gene2 <- "slc1a3a"
# expr1 <- FetchData(object = combined, vars = gene1)
# expr2 <- FetchData(object = combined, vars = gene2)
# duo <- colnames(combined[, which(x = expr1 > 1 & expr2 > 1)])
# non_duo <- colnames(combined[, which(x = expr1 < 1 | expr2 < 1)])
# duo <- WhichCells(object = combined, expression = gene1 > 1 & gene2 > 1)
# non_duo <- WhichCells(object = combined, expression = gene1 < 1 | gene2 < 1)
combined <- SetIdent(combined, cells=duo, value=paste(gene1, gene2, sep="_"))
combined <- SetIdent(combined, cells=non_duo, value="non_duo")
DimPlot(combined, reduction="umap", group.by = "ident", split.by="cond", pt.size=1.5, order=TRUE)


gene3 <- "fosb"
expr1 <- FetchData(object = combined, vars = gene3)
# cells_in_cluster <- WhichCells(combined, idents = "5")
# ctrl_cells_in_cluster <- cells_in_cluster[startsWith(WhichCells(combined, idents = "5"), "CTRL")]
# bhve_cells_in_cluster <- cells_in_cluster[startsWith(WhichCells(combined, idents = "5"), "BHVE")]
test <- combined[, which(x = expr1 > 1)]
cells_in_cluster <- WhichCells(test, idents = "5")
ctrl_cells_in_cluster <- cells_in_cluster[startsWith(WhichCells(test, idents = "5"), "CTRL")]
bhve_cells_in_cluster <- cells_in_cluster[startsWith(WhichCells(test, idents = "5"), "BHVE")]
combined@assays$RNA@counts[gene3,cells_in_cluster]
boxplot(combined@assays$RNA@counts[gene3,cells_in_cluster])
boxplot(combined@assays$RNA@counts[gene3,ctrl_cells_in_cluster])
boxplot(combined@assays$RNA@counts[gene3,bhve_cells_in_cluster])



gene3 <- "fosb"
expr1 <- FetchData(object = combined, vars = gene3)
test <- combined[, which(x = expr1 > 1)]
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
total_b1 <- c()
total_b2 <- c()
total_c1 <- c()
Idents(test) <- test$sample
test$cond.cluster <- paste(Idents(test), test$seurat_clusters, sep = "_")
Idents(test) <- test$cond.cluster
for (i in 0:num_clusters) {
  b1_cells_in_cluster <- 0
  b2_cells_in_cluster <- 0
  c1_cells_in_cluster <- 0
  
  try(b1_cells_in_cluster <- length(WhichCells(test, idents = paste("b1", i, sep="_"))), silent=TRUE)
  try(b2_cells_in_cluster <- length(WhichCells(test, idents = paste("b2", i, sep="_"))), silent=TRUE)
  try(c1_cells_in_cluster <- length(WhichCells(test, idents = paste("c1", i, sep="_"))), silent=TRUE)
  
  total_b1 <- c(total_b1, b1_cells_in_cluster)
  total_b2 <- c(total_b2, b2_cells_in_cluster)
  total_c1 <- c(total_c1, c1_cells_in_cluster)
}
df <- data.frame(condition <- c(rep("b1", length(total_b1)), rep("b2", length(total_b2)), rep("c1", length(total_c1))),
                 cluster_num <- c(0:40, 0:40, 0:40),
                 value <- c(total_b1, total_b2, total_c1))
df2 <- data.frame(condition <- c(rep("b1", length(total_b1)), rep("c1", length(total_c1))),
                 cluster_num <- c(0:40, 0:40),
                 value <- c(total_b1, total_c1))
p <- ggplot(df2, aes(fill=condition, x=cluster_num, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  ggtitle(paste("Number of Cells Expressing", gene3, "per Cluster")) +
  xlab("Cluster") +
  ylab("Number of Cells") +
  scale_x_continuous(breaks = 0:40)
p

# bp <- text(bp, x, labels = symnum(x, cutpoints = c(1, 5, 10), symbols = c("*", "**")), pos = 3)
#   1 2 3 4
# C
# B

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

# Mouse Atlas
library(loomR)
library("scrattch.io")
options(stringsAsFactors = FALSE)
combined <- readRDS("C:/Users/miles/Downloads/brain/brain_scripts/brain_shiny/data/combined.rds")
# lfile <- connect(filename = "C:/Users/miles/Downloads/transcriptome/transcrip.tome", mode = "r+")
# tome_mat <- read_tome_dgCMatrix("C:/Users/miles/Downloads/transcriptome/transcrip.tome", )
exons       <- read_tome_dgCMatrix("C:/Users/miles/Downloads/brain/data/transcrip.tome","data/t_exon")    # (or data/exon)
introns     <- read_tome_dgCMatrix("C:/Users/miles/Downloads/transcriptome/transcrip.tome","data/t_intron")  # (or data/intron)
sample_name <- read_tome_sample_names("C:/Users/miles/Downloads/transcriptome/transcrip.tome")  
gene_name   <- read_tome_gene_names("C:/Users/miles/Downloads/transcriptome/transcrip.tome")

# read_tome_dgCMatrix("/nv/hp10/ggruenhagen3/scratch/brain/data/transcriptome.zip", "/data/t_exon")
#

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

