library(Seurat)
library(dplyr)

suppressMessages(require(Seurat))
suppressMessages(require(ggplot2))
suppressMessages(require(cowplot))

install.packages('scater')
library(scater)

#install scran
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scran")
suppressMessages(require(scran))
library(scran)

BiocManager::install("BiocParallel")
suppressMessages(require(BiocParallel))

BiocManager::install("BiocNeighbors")
suppressMessages(require(BiocNeighbors))

rna_path <-"C:/Users/miles/Downloads/brain/" 

b1.data <- Read10X(data.dir = paste(rna_path, "data/B1-from-bcl-lncRNA-fake/outs/filtered_feature_bc_matrix/", sep=""))
b2.data <- Read10X(data.dir = paste(rna_path, "data/B2-from-bcl-lncRNA-fake/outs/filtered_feature_bc_matrix/", sep=""))
c1.data <- Read10X(data.dir = paste(rna_path, "data/C1-from-bcl-lncRNA-fake/outs/filtered_feature_bc_matrix/", sep=""))
c2.data <- Read10X(data.dir = paste(rna_path, "data/MZ-from-bcl-lncRNA-fake/outs/filtered_feature_bc_matrix/", sep=""))

b1 <- CreateSeuratObject(counts = b1.data, project = "behav")
b2 <- CreateSeuratObject(counts = b2.data, project = "behav2")
c1 <- CreateSeuratObject(counts = c1.data, project = "control")
c2 <- CreateSeuratObject(counts = c2.data, project = "control2")

b1$sample <- "b1"
b2$sample <- "b2"
c1$sample <- "c1"
c2$sample <- "c2"

b1$cond <- "BHVE"
c1$cond <- "CTRL"
b2$cond <- "BHVE"
c2$cond <- "CTRL"

b1 <- NormalizeData(b1, normalization.method = "LogNormalize")
b2 <- NormalizeData(b2, normalization.method = "LogNormalize")
c1 <- NormalizeData(c1, normalization.method = "LogNormalize")
c2 <- NormalizeData(c2, normalization.method = "LogNormalize")

b1 <- subset(b1, subset = nFeature_RNA > 500)
b1 <- subset(b1, subset = nFeature_RNA < 2500)
b1 <- subset(b1, subset = nCount_RNA > 500)
b1 <- subset(b1, subset = nCount_RNA < 6000)
b1 <- FindVariableFeatures(object = b1, selection.method = "vst", verbose = TRUE)

b1_var = b1@assays$RNA@var.features

# Set up behave2 object
b2 <- subset(b2, subset = nFeature_RNA > 500)
b2 <- subset(b2, subset = nFeature_RNA < 2500)
b2 <- subset(b2, subset = nCount_RNA > 500)
b2 <- subset(b2, subset = nCount_RNA < 6000)
b2 <- FindVariableFeatures(object = b2, selection.method = "vst", verbose = TRUE)

b2_var = b2@assays$RNA@var.features

# Set up control object
c1 <- subset(c1, subset = nFeature_RNA > 500)
c1 <- subset(c1, subset = nFeature_RNA < 2500)
c1 <- subset(c1, subset = nCount_RNA > 500)
c1 <- subset(c1, subset = nCount_RNA < 6000)
c1 <- FindVariableFeatures(object = c1, selection.method = "vst", verbose = TRUE)

c1_var = c1@assays$RNA@var.features

# Set up control object
c2 <- subset(c2, subset = nFeature_RNA > 500)
c2 <- subset(c2, subset = nFeature_RNA < 2500)
c2 <- subset(c2, subset = nCount_RNA > 500)
c2 <- subset(c2, subset = nCount_RNA < 6000)
c2 <- FindVariableFeatures(object = c2, selection.method = "vst", verbose = TRUE)

c2_var = c2@assays$RNA@var.features

combined <- merge(x=b1, y=c(b2,c1,c2), merge.data = TRUE, add.cell.ids = c("BHVE", "BHVE", "CTRL", "CTRL"))

#combined <- readRDS("C:/Users/zjohnson37/Desktop/from-bcl/c41_scale1mil_b1b2c1_012219 (1).rds")
#combined <- readRDS("C:/Users/zjohnson37/Desktop/from-bcl/c41_scale1mil_b1b2c1_121619.rds")
DefaultAssay(combined) <- "RNA"

brain = readRDS("C:/Users/miles/Downloads/brain/brain_scripts/brain_mz_shiny/data/B1C1C2MZ_combined_031020.rds")
combined@assays$RNA@var.features = brain@assays$RNA@var.features
# combined <- FindVariableFeatures(object = combined, selection.method = "vst", nfeatures = 10000, verbose = TRUE)

DefaultAssay(combined) <- "RNA"

combined <- ScaleData(combined, verbose = TRUE)
combined <- RunPCA(combined, dim = 50, verbose = TRUE)

combined <- FindNeighbors(combined, reduction = "pca", dims = 1:50)
combined <- FindClusters(combined, resolution = 1.86)
combined <- RunUMAP(combined, dims = 1:50)

saveRDS(combined, file = "C:/Users/zjohnson37/Desktop/from-bcl/B1C1C2MZ_combined_031320.rds")
combined <- readRDS(file = "C:/Users/zjohnson37/Desktop/from-bcl/B1C1C2MZ_combined_031020.rds")
combined <- readRDS(file = "C:/Users/zjohnson37/Desktop/from-bcl/B1C1C2MZ_combined_031320.rds")

p0 <- DimPlot(combined, reduction = "umap", label = TRUE)
saveRDS(combined, file = "C:/Users/miles/Downloads/brain/data/lncRNA_2.RDS")
deg = FindAllMarkers(lncRNA, only.pos = T)
deg = deg[which(deg$p_val_adj < 0.05),] # adjust p val is  less tha 0.05
ggplot(deg[which(deg$gene %in% lncRNA_genes),], aes(cluster)) + geom_bar() + labs(title="DEG per Cluster")
norm_df = data.frame()
for ( cluster in unique(deg$cluster) ) {
  this_rows = deg[which(deg$cluster == cluster),]
  lncRNA_rows = this_rows[which(this_rows$gene %in% lncRNA_genes),]
  norm_df = rbind(norm_df, t(c(cluster, (nrow(lncRNA_rows)/nrow(this_rows))*100 )))
}
colnames(norm_df) = c("cluster", "lncRNA_deg_pct")
norm_df$cluster = factor(norm_df$cluster)
norm_df$lncRNA_deg_pct = as.numeric(as.vector(norm_df$lncRNA_deg_pct))
ggplot(norm_df, aes(cluster, lncRNA_deg_pct)) + geom_bar(stat ="identity") 

###############
# ncRNA Plots #
###############
lncRNA <- readRDS(file = "C:/Users/miles/Downloads/brain/data/lncRNA_2.RDS")

type = "lncRNA"
mat <- as.matrix(lncRNA@assays$RNA@counts[lncRNA_genes,])
mat[mat > 0] <- 1
lncRNA$thisSum <- colSums(mat)
mat2 <- as.matrix(lncRNA@assays$RNA@counts)
mat2[mat2 > 0] <- 1
lncRNA$thisPct <- (lncRNA$thisSum/colSums(mat2)) * 100
png(paste0(rna_path, "/results/ncRNA/", type, "_unique_per_cell.png"), width = 1000, height = 600)
p <- FeaturePlot(lncRNA, features = "thisPct", order = TRUE, label = TRUE, pt.size = 1.5) + ggtitle(paste(type, "Sum"))
print(p)
dev.off()
df <- data.frame()
num_clusters <- as.numeric(tail(levels(lncRNA@meta.data$seurat_clusters), n=1))
for (i in 0:num_clusters) {
  this_cells <- WhichCells(lncRNA, idents = i)
  # thisSums <- rowSums(lncRNA@assays$RNA@counts[lncRNA_genes,this_cells])
  # thisSums <- length(thisSums[which(thisSums > 0)])
  # df <- rbind(df, t(c(i, thisSums)))
  df <- rbind(df, t(c(i, mean(lncRNA$thisSum[this_cells]))))
}
df$V1 <- factor(df$V1)
png(paste0(rna_path, "/results/ncRNA/", type, "_lncRNA_present_pct_cluster.png"), width = 900, height = 450)
p <- ggplot(df, aes(x = V1, y = V2)) + geom_bar(stat = "identity") + xlab("Cluster") + ylab(paste0("Avg Pct lncRNA Genes of Genes Present")) + ggtitle(paste(type))
print(p)
dev.off()

type = "miRNA"
lncRNA$thisSum <- colSums(as.matrix(lncRNA@assays$RNA@counts[miRNA,]))
lncRNA$thisPct <- (lncRNA$thisSum/lncRNA$nCount_RNA) * 100
png(paste0(rna_path, "/results/ncRNA/", type, ".png"), width = 1000, height = 600)
p <- FeaturePlot(lncRNA, features = "thisPct", order = TRUE, label = TRUE, pt.size = 1.5) + ggtitle(type)
print(p)
dev.off()
df <- data.frame()
num_clusters <- as.numeric(tail(levels(lncRNA@meta.data$seurat_clusters), n=1))
for (i in 0:num_clusters) {
  this_cells <- WhichCells(lncRNA, idents = i)
  df <- rbind(df, t(c(i, mean(lncRNA$thisPct[this_cells]) * 100)))
}
df$V1 <- factor(df$V1)
png(paste0(rna_path, "/results/ncRNA/", type, "_cluster.png"), width = 900, height = 450)
p <- ggplot(df, aes(x = V1, y = V2)) + geom_bar(stat = "identity") + xlab("Cluster") + ylab(paste0("Avg Pct ", type)) + ggtitle(type)
print(p)
dev.off()

type = "all_ncRNA"
lncRNA$thisSum <- colSums(as.matrix(lncRNA@assays$RNA@counts[c(lncRNA_genes, miRNA, Mt_rRNA, Mt_tRNA, scaRNA, snoRNA),]))
lncRNA$thisPct <- (lncRNA$thisSum/lncRNA$nCount_RNA) * 100
png(paste0(rna_path, "/results/ncRNA/", type, ".png"), width = 1000, height = 600)
p <- FeaturePlot(lncRNA, features = "thisPct", order = TRUE, label = TRUE, pt.size = 1.5) + ggtitle(type)
print(p)
dev.off()
df <- data.frame()
num_clusters <- as.numeric(tail(levels(lncRNA@meta.data$seurat_clusters), n=1))
for (i in 0:num_clusters) {
  this_cells <- WhichCells(lncRNA, idents = i)
  df <- rbind(df, t(c(i, mean(lncRNA$thisPct[this_cells]))))
}
df$V1 <- factor(df$V1)
png(paste0(rna_path, "/results/ncRNA/", type, "_cluster.png"), width = 900, height = 450)
p <- ggplot(df, aes(x = V1, y = V2)) + geom_bar(stat = "identity") + xlab("Cluster") + ylab(paste0("Avg Pct ", type)) + ggtitle(type)
print(p)
dev.off()

###################################
# Which lncRNA are driving 17&39? #
###################################
# Method 1: Fisher's Test
lncRNA <- readRDS("C:/Users/miles/Downloads/brain/data/lncRNA_2.RDS")
tj <- readRDS("C:/Users/miles/Downloads/d_tooth/tooth_scripts/tj_shiny/data/tj.rds")
lncRNA_genes <- rownames(lncRNA)[which(! rownames(lncRNA) %in% rownames(tj))]
non_zero_lncRNA = lncRNA_genes[which(rowSums(as.matrix(lncRNA@assays$RNA@counts[lncRNA_genes,])) != 0)]
num_clusters <- as.numeric(tail(levels(lncRNA@meta.data$seurat_clusters), n=1))

fisher_df <- data.frame()
ptm <- proc.time()
for (i in 1:length(lncRNA_genes)) {
  if (i %% 100 == 0) {
    print(i)
    proc.time() - ptm
  }
# for (gene in test) {
  # ptm <- proc.time()
  gene <- lncRNA_genes[i]
  pos_cells <- colnames(lncRNA)[which(lncRNA@assays$RNA@counts[gene,] != 0)]
  for (cluster in 0:num_clusters) {
    cluster_cells <- WhichCells(lncRNA, idents = cluster)
    not_cluster_cells <- colnames(lncRNA)[which(! colnames(lncRNA) %in% cluster_cells)]
    cluster_pos <- length(pos_cells[which(pos_cells %in% cluster_cells)])
    cluster_neg <- length(cluster_cells[which(!cluster_cells %in% pos_cells)])
    not_cluster_pos <- length(pos_cells[which(pos_cells %in% not_cluster_cells)])
    not_cluster_neg <- length(not_cluster_cells[which(!not_cluster_cells %in% pos_cells)])
    fisher_p <- fisher.test(data.frame(c(cluster_pos, cluster_neg), c(not_cluster_pos, not_cluster_neg)))$p.value
    fisher_df <- rbind(fisher_df, t(c(gene, cluster, fisher_p, cluster_pos, cluster_neg, not_cluster_pos, not_cluster_neg)))
  }
  # proc.time() - ptm
}
proc.time() - ptm
colnames(fisher_df) <- c("gene", "cluster", "fisher_p", "cluster_pos", "cluster_neg", "other_pos", "other_neg")
fisher_df$fisher_q <- p.adjust(as.numeric(as.vector(fisher_df$fisher_p)), method = "bonferroni")
nrow(fisher_df[which(fisher_df$fisher_q < 0.05),])
write.csv(fisher_df[which(fisher_df$fisher_q < 0.05),], paste0("C:/Users/miles/Downloads/brain/results/ncRNA/fisher_lncRNA.csv"), quote = FALSE, row.names = FALSE)
ggplot(fisher_df_sig, aes(factor(cluster, levels=0:num_clusters))) + geom_bar() + xlab("Cluster") + ylab("Number of Sig lncRNA genes") + ggtitle("Number of Significant lncRNA genes (from Fisher's Test) per Cluster")

# Method 2: Above/Below Avg Exp
gene_df <- data.frame()
ptm <- proc.time()
tot_cells <- ncol(lncRNA)
for (i in 1:length(lncRNA_genes)) {
  if (i %% 100 == 0) {
    print(i)
    print(proc.time() - ptm)
  }
  # for (gene in test) {
  # ptm <- proc.time()
  gene <- lncRNA_genes[i]
  pos_cells <- colnames(lncRNA)[which(lncRNA@assays$RNA@counts[gene,] != 0)]
  gene_sum <- sum(lncRNA@assays$RNA@counts[gene,])
  gene_avg_exp <- gene_sum/tot_cells
  for (cluster in 0:num_clusters) {
    cluster_cells <- WhichCells(lncRNA, idents = cluster)
    cluster_pos <- sum(lncRNA@assays$RNA@counts[gene,cluster_cells])
    cluster_gene_avg_exp <- cluster_pos/length(cluster_cells)
    gene_df <- rbind(gene_df, t(c(gene, cluster, cluster_gene_avg_exp, gene_avg_exp)))
  }
}
colnames(gene_df) <- c("gene", "cluster", "cluster_exp", "all_exp")
gene_df$isHigher <- as.numeric(as.vector(gene_df$cluster_exp)) > as.numeric(as.vector(gene_df$all_exp))
ggplot(gene_df, aes(cluster)) + geom_bar(aes(fill = isHigher)) + ylab("Number of lncRNA Genes") + ggtitle("Number of lncRNA Genes With Expression Greater Than Average")

# Find background for each cluster
df <- data.frame()
for (cluster in 0:num_clusters) {
  cluster_cells <- WhichCells(lncRNA, idents = cluster)
  df <- rbind(df, t(c(cluster, mean(lncRNA$nCount_RNA[cluster_cells]))))
}
ggplot(df, aes(factor(V1), V2)) + geom_bar(stat = "identity") + xlab("Cluster") + ylab("# Transcripts per Cell") + ggtitle("# Transcripts per Cell per Cluster")

per_cluster_df <- data.frame()
for (i in 1:length(lncRNA_genes)) {
  if (i %% 100 == 0) {
    print(i)
    print(proc.time() - ptm)
  }
  gene <- lncRNA_genes[i]
  for (cluster in 0:num_clusters) {
    cluster_cells <- WhichCells(lncRNA, idents = cluster)
    cluster_sum <- sum(lncRNA@assays$RNA@counts[gene, cluster_cells])
    per_cluster_df <- rbind(per_cluster_df, t(c(gene, cluster, cluster_sum)))
  }
}
colnames(per_cluster_df) <- c("gene", "cluster", "cluster_sum")
ggplot(per_cluster_df, aes(cluster, as.numeric(as.vector(cluster_sum)), fill=gene)) + geom_bar(stat = "identity") + NoLegend() + ylab("Number of Transcripts") + ggtitle("Number of Transcripts for each lncRNA gene")

# Sum per cluster
lncRNA$thisSum <- colSums(as.matrix(lncRNA@assays$RNA@counts[lncRNA_genes,]))
# lncRNA$thisPct <- (lncRNA$thisSum/lncRNA$nCount_RNA) * 100
df <- data.frame()
for (cluster in 0:num_clusters) {
  cluster_cells <- WhichCells(lncRNA, idents = cluster)
  df <- rbind(df, t(c(cluster, sum(lncRNA$thisSum[cluster_cells]))))
}
ggplot(df, aes(factor(V1), V2)) + geom_bar(stat = "identity") + xlab("Cluster") + ylab("# lncRNA Transcripts") + ggtitle("# lncRNA Transcripts per Cluster")

############
# DEG Plot #
############
# deg_file <- "/nv/hp10/ggruenhagen3/scratch/brain/data/lncRNA_lncRNA_deg.csv"
deg_file <- "C:/Users/miles/Downloads/brain/results/lncRNA_lncRNA_deg.csv"
deg <- read.csv(deg_file, header = TRUE, stringsAsFactors = FALSE)
deg$cluster <- factor(deg$cluster)
ggplot(deg, aes(cluster)) + geom_bar() + xlab("Cluster") + ylab("Number of DEGs") + ggtitle("lncRNA genes in DEGs")

##################################
# lncRNA DEG Neighbor Expression #
##################################
# lncRNA <- readRDS("/nv/hp10/ggruenhagen3/scratch/brain/data/lncRNA.RDS")
# combined <- readRDS("/nv/hp10/ggruenhagen3/scratch/brain/brain_scripts/brain_shiny/data/combined.rds")
# lncRNA_genes <- rownames(lncRNA)[which(! rownames(lncRNA) %in% rownames(combined))]
# gtf <- read.table("/nv/hp10/ggruenhagen3/scratch/brain/data/ens_w_ncbi_sort.gtf", sep = "\t", stringsAsFactors = FALSE)
library(tidyverse)
gtf <- read.table(file = "C:/Users/miles/Downloads/brain/data/ens_w_ncbi_2_sort.gtf", sep = "\t", stringsAsFactors = FALSE)
gtf_gene_name <- c()
for (i in 1:nrow(gtf)) {
  start <- gregexpr(pattern ='gene_name', gtf$V9[i])[[1]]
  stop  <- gregexpr(pattern =';', substr(gtf$V9[i], start, nchar(gtf$V9[i])))[[1]][1]
  gene_name <- substr(gtf$V9[i], start+10, start+stop-2)
  if (start == -1) {
    gene_name <- substr(gtf$V9[i], start+10, start+stop)
  }
  gtf_gene_name <- c(gtf_gene_name, gene_name)
}
gtf$gene_name <- gtf_gene_name
colnames(gtf) <- c("LG", "source", "type", "start", "stop", "idk", "idk1", "idk2", "info", "gene_name")
previous_gene_name = ""
gtf_mod = data.frame()
for ( i in 1:nrow(gtf)) {
  cur_gene_name = gtf$gene_name[i]
  if (cur_gene_name != previous_gene_name) {
    gtf_mod = rbind(gtf_mod, gtf[i,])
  }
  previous_gene_name = cur_gene_name
}
gtf = gtf_mod
# gtf <- separate(gtf, col = 9, into = c("gene_id", "gene_name"), sep=";")
# gtf$gene_name[which(startsWith(gtf$gene_name, " gene_name"))] <- substr(gtf$gene_name[which(startsWith(gtf$gene_name, " gene_name"))], 12, 10000L)
deg_17 <- deg$X[which(deg$cluster == 17 & deg$avg_logFC > 0)]
neighbor_df <- data.frame()
neighbor_degree <- c(-2, -1, 1, 2)
for (deg in deg_17) {
  if (grepl(".", deg, fixed = TRUE)) {
    stop <- gregexpr(pattern ='\\.', deg)[[1]]
    deg <- substr(deg, 0, stop-1)
  }  
  deg_i <- as.numeric(rownames(gtf[which(gtf$gene_name == deg),]))
  neighbors <- list()
  for (degree in neighbor_degree) {
    if (deg_i + degree >= 1 && deg_i + degree < length(gtf$gene_name)) { # weird errors on PACE, maybe it went to 0
      neighbors[as.character(degree)] = gtf$gene_name[deg_i + degree]
    }
  }
  if (all(neighbors %in% rownames(lncRNA))) {
    for (cluster in 0:num_clusters) {
      this_cells <- WhichCells(lncRNA, idents = cluster)
      # newRow <- data.frame()
      for (degree in neighbor_degree) {
        neighbor <- neighbors[[as.character(degree)]]
        this_mean <- mean(lncRNA@assays$RNA@counts[neighbor, this_cells])
        this_mean_norm <- this_mean/mean(lncRNA$nCount_RNA[this_cells])
        neighbor_df <- rbind(neighbor_df, t(c(deg, neighbor, degree, cluster, this_mean, this_mean_norm)))
      }
    } 
  } else {
    print(deg)
    for (degree in neighbor_degree) {
      print(neighbors[as.character(degree)])
    }
    print("------------------------")
  }
}
colnames(neighbor_df) <- c("deg", "neighbor", "degree", "cluster", "mean", "mean_norm")
neighbor_df$cluster <- factor(neighbor_df$cluster)
neighbor_df$degree <- factor(neighbor_df$degree)
neighbor_df$mean_norm <- as.numeric(as.vector(neighbor_df$mean_norm))
# ggplot(neighbor_df, aes(cluster, mean_norm, fill = degree)) + geom_bar(stat = "identity", position = position_dodge()) + ylab("Mean expression of Neighbor of all lncRNA / Mean expression of All Genes") + scale_fill_manual(values = c("gray", "gray47"))
# ggplot(neighbor_df, aes(cluster, mean, fill = degree)) + geom_bar(stat = "identity", position = position_dodge()) + ylab("Mean expression of Neighbor of all lncRNA") + scale_fill_manual(values = c("gray", "gray47"))
ggplot(neighbor_df, aes(cluster, mean_norm, fill = degree)) + stat_summary(fun = "mean", geom = "bar", position=position_dodge()) + ylab("Mean expression of Neighbor of lncRNA DEG in Cluster 17 / Mean expression of All Genes") + scale_fill_manual(values = c("#111d5e", "#c70039", "#f37121", "#ffbd69"))
# ggplot(neighbor_df, aes(cluster, mean, fill = degree)) + geom_bar(stat = "identity") + ylab("Mean expression of Neighbor of lncRNA DEG in Cluster 17") + scale_fill_manual(values = c("#ff9c71", "#654062"))
# saveRDS(neighbor_df, file = "C:/Users/miles/Downloads/brain/data/all_lncRNA_genes.RDS")
# test_df <- data.frame( c("a", "a", "a", "a", "b", "b", "b", "b"), c(4, 5, 5, 6, 8, 9, 9, 10), c("1", "1", "2", "2", "1", "1", "2", "2") )
# colnames(test_df) <- c("letters", "numbers", "replicate")
# ggplot(test_df, aes(letters, numbers, fill = replicate)) + geom_bar(stat = "identity", position=position_dodge()) + scale_y_continuous(breaks = c(0:20))
neighbor_df_2 <- data.frame()
for (cluster in 0:num_clusters) {
  this_rows <- neighbor_df[which(neighbor_df$cluster == cluster),]
  mean_norms <- c()
  for (degree in neighbor_degree) {
    cluster_mean_degree <- mean(this_rows$mean_norm[which(this_rows$degree == as.character(degree))])
    mean_norms <- c(mean_norms, cluster_mean_degree)
  }
  this_order <- sort(mean_norms, decreasing = TRUE, index.return = TRUE)$ix
  newRow <- data.frame( rep(cluster, length(neighbor_degree)), neighbor_degree[this_order], 1:length(neighbor_degree)  )
  # neighbor_df_2 <- rbind(neighbor_df_2, t(c(cluster, max_degree)))
  neighbor_df_2 <- rbind(neighbor_df_2, newRow)
}
colnames(neighbor_df_2) <- c("cluster", "degree", "rank")
neighbor_df_2$degree <- factor(neighbor_df_2$degree, levels=neighbor_degree)
neighbor_df_2$rank <- factor(neighbor_df_2$rank)
ggplot(neighbor_df_2, aes(x=rank, fill = degree)) + geom_bar()  + scale_fill_manual(values = c("#111d5e", "#c70039", "#f37121", "#ffbd69")) + ggtitle("Ranking of Neighbors Across All Clusters From Highest to Lowest Expression")

neighbor_df_3 <- data.frame()
for (degree in neighbor_degree) {
  this_rows <- neighbor_df[which(neighbor_df$degree == degree),]
  this_genes <- unique(this_rows$neighbor)
  this_mean <- colSums(lncRNA@assays$RNA@counts[this_genes,])
  this_mean_norm <- mean(this_mean/lncRNA$nCount_RNA)
  neighbor_df_3 <- rbind(neighbor_df_3, t(c(as.character(degree),this_mean_norm)))
}
colnames(neighbor_df_3) <- c("degree", "mean_norm")
ggplot(neighbor_df_3, aes(x = degree, y=as.numeric(as.vector(mean_norm)), fill = factor(degree))) + geom_bar(stat="identity") + ylab("Mean expression of Neighbor of lncRNA DEG in Cluster 17 / Mean expression of All Genes") + scale_fill_manual(values = c("#111d5e", "#c70039", "#f37121", "#ffbd69")) + ggtitle("Average Expression for Each Neighbor Across All Cells")


#########################################################################
# Separate lncRNA with 0 transcripts from lncRNA with non 0 transcripts #
#########################################################################
# non_zero lncRNAs
neighbor_df_3_1 = data.frame()
neighbor_degree = -20:20
neighbor_degree = neighbor_degree[which(neighbor_degree != 0)]
n_thrown = 0
non_zero_lncRNA = lncRNA_genes[which(rowSums(as.matrix(lncRNA@assays$RNA@counts[lncRNA_genes,])) != 0)]
for (deg in non_zero_lncRNA) {
  if (grepl(".", deg, fixed = TRUE)) {
    stop <- gregexpr(pattern ='\\.', deg)[[1]]
    deg <- substr(deg, 0, stop-1)
  }  
  deg_i <- as.numeric(rownames(gtf[which(gtf$gene_name == deg),]))
  neighbors <- list()
  for (degree in neighbor_degree) {
    if (deg_i + degree >= 1 && deg_i + degree < length(gtf$gene_name)) { # weird errors on PACE, maybe it went to 0
      neighbors[as.character(degree)] = gtf$gene_name[deg_i + degree]
    }
  }
  if (all(neighbors %in% rownames(lncRNA))) {
    for (degree in neighbor_degree) {
      neighbor_df_3_1 = rbind(neighbor_df_3_1, t(c(deg, degree, neighbors[[as.character(degree)]])))
    }
  } else {
        n_thrown = n_thrown+1
        print(deg)
        for (degree in neighbor_degree) {
          print(neighbors[[as.character(degree)]])
        }
        print("------------------------")
  }
}
print(n_thrown)
colnames(neighbor_df_3_1) = c("deg", "degree", "neighbor")
neighbor_df_3_2 = data.frame()
for (degree in neighbor_degree) {
  all_degree_neighbors = unique(neighbor_df_3_1$neighbor[which(neighbor_df_3_1$degree == degree)])
  newRow = data.frame(rep(degree, ncol(lncRNA)), colSums(lncRNA@assays$RNA@counts[all_degree_neighbors,])/lncRNA$nCount_RNA)
  neighbor_df_3_2 = rbind(neighbor_df_3_2, newRow)
  # neighbor_df_3_2 = rbind(neighbor_df_3_2, t(c(degree, mean(colSums(lncRNA@assays$RNA@counts[all_degree_neighbors,])/lncRNA$nCount_RNA))))
  # neighbor_df_3_2 = rbind(neighbor_df_3_2, t(c(degree, mean(colSums(lncRNA@assays$RNA@data[all_degree_neighbors,])))))
}
colnames(neighbor_df_3_2) <- c("degree", "mean")
neighbor_df_3_2$degree = factor(neighbor_df_3_2$degree, levels = neighbor_degree)
# ggplot(neighbor_df_3_2, aes(x=degree, y=mean, group=degree, fill=degree, color=degree)) + geom_boxplot(alpha=0.9) + xlab("Neighbor to lncRNA Gene") + geom_jitter(position=position_jitter(), alpha=0.01) + NoLegend() + ylab("# Neighbor Transcripts/# All Transcripts")

neighbor_df_3_2$cell = rownames(neighbor_df_3_2)
neighbor_df_3_2_non_0 = neighbor_df_3_2

# zero lncRNAs
neighbor_df_3_1 = data.frame()
neighbor_degree = -20:20
neighbor_degree = neighbor_degree[which(neighbor_degree != 0)]
n_thrown = 0
non_zero_lncRNA = lncRNA_genes[which(rowSums(as.matrix(lncRNA@assays$RNA@counts[lncRNA_genes,])) == 0)]
for (deg in non_zero_lncRNA) {
  if (grepl(".", deg, fixed = TRUE)) {
    stop <- gregexpr(pattern ='\\.', deg)[[1]]
    deg <- substr(deg, 0, stop-1)
  }  
  deg_i <- as.numeric(rownames(gtf[which(gtf$gene_name == deg),]))
  neighbors <- list()
  for (degree in neighbor_degree) {
    if (deg_i + degree >= 1 && deg_i + degree < length(gtf$gene_name)) { # weird errors on PACE, maybe it went to 0
      neighbors[as.character(degree)] = gtf$gene_name[deg_i + degree]
    }
  }
  if (all(neighbors %in% rownames(lncRNA))) {
    for (degree in neighbor_degree) {
      neighbor_df_3_1 = rbind(neighbor_df_3_1, t(c(deg, degree, neighbors[[as.character(degree)]])))
    }
  } else {
    n_thrown = n_thrown+1
    print(deg)
    for (degree in neighbor_degree) {
      print(neighbors[[as.character(degree)]])
    }
    print("------------------------")
  }
}
print(n_thrown)
colnames(neighbor_df_3_1) = c("deg", "degree", "neighbor")
neighbor_df_3_2 = data.frame()
for (degree in neighbor_degree) {
  all_degree_neighbors = unique(neighbor_df_3_1$neighbor[which(neighbor_df_3_1$degree == degree)])
  newRow = data.frame(rep(degree, ncol(lncRNA)), colSums(lncRNA@assays$RNA@counts[all_degree_neighbors,])/lncRNA$nCount_RNA)
  neighbor_df_3_2 = rbind(neighbor_df_3_2, newRow)
  # neighbor_df_3_2 = rbind(neighbor_df_3_2, t(c(degree, mean(colSums(lncRNA@assays$RNA@counts[all_degree_neighbors,])/lncRNA$nCount_RNA))))
  # neighbor_df_3_2 = rbind(neighbor_df_3_2, t(c(degree, mean(colSums(lncRNA@assays$RNA@data[all_degree_neighbors,])))))
}
colnames(neighbor_df_3_2) <- c("degree", "mean")
neighbor_df_3_2$degree = factor(neighbor_df_3_2$degree, levels = neighbor_degree)
# ggplot(neighbor_df_3_2, aes(x=degree, y=mean, group=degree, fill=degree, color=degree)) + geom_boxplot(alpha=0.9) + xlab("Neighbor to lncRNA Gene") + geom_jitter(position=position_jitter(), alpha=0.01) + NoLegend() + ylab("# Neighbor Transcripts/# All Transcripts")

neighbor_df_3_2$cell = rownames(neighbor_df_3_2)
neighbor_df_3_2_0 = neighbor_df_3_2

neighbor_df_3_2_diff = neighbor_df_3_2_non_0
neighbor_df_3_2_diff$mean = neighbor_df_3_2_non_0$mean - neighbor_df_3_2_0$mean
ggplot(neighbor_df_3_2_diff, aes(x=degree, y=mean, group=degree, fill=degree, color=degree)) + geom_boxplot(alpha=0.9) + xlab("Neighbor to lncRNA Gene") + geom_jitter(position=position_jitter(), alpha=0.01) + NoLegend() + ylab("Difference in # Neighbor Transcripts/# All Transcripts")

###################################################
# Only sum the neighbors if the lncRNA is present #
###################################################
neighbor_df_3_1 = data.frame()
neighbor_degree = -20:20
neighbor_degree = neighbor_degree[which(neighbor_degree != 0)]
n_thrown = 0
n_success = 0
i = 0
lncRNA_cells = list()
non_zero_lncRNA = lncRNA_genes[which(rowSums(as.matrix(lncRNA@assays$RNA@counts[lncRNA_genes,])) != 0)]
non_zero_norm_genes = rownames(lncRNA)[which(rowSums(lncRNA@assays$RNA@counts) > 0 & ! rownames(lncRNA) %in% lncRNA_genes)]
non_zero_norm_genes = non_zero_norm_genes[which(! grepl(".1", non_zero_norm_genes))]
# for (deg in non_zero_lncRNA) {
while (n_success < 697) {
  # Use this block for normal genes
  i = i + 1
  deg = non_zero_norm_genes[i]
  
  # Use this block for lncRNA
  # if (grepl(".", deg, fixed = TRUE)) {
  #   stop <- gregexpr(pattern ='\\.', deg)[[1]]
  #   deg <- substr(deg, 0, stop-1)
  # }
  
  deg_i <- as.numeric(which(gtf$gene_name == deg))
  lncRNA_cells[[deg]] = colnames(lncRNA)[which(lncRNA@assays$RNA@counts[deg,] != 0)]
  deg_lg = gtf$LG[deg_i]
  deg_start = gtf$start[deg_i]
  neighbors <- list()
  neighbors_dist = list()

  for (degree in neighbor_degree) {
    new_i = deg_i + degree
    if (new_i >= 1 && new_i < length(gtf$gene_name) && gtf$LG[new_i] == deg_lg ) { # weird errors on PACE, maybe it went to 0
      neighbors[[as.character(degree)]] = gtf$gene_name[new_i]
      neighbors_dist[[as.character(degree)]] = gtf$start[new_i]-deg_start
    }
  }
  
  if (length(neighbors) == length(neighbor_degree) && all(neighbors %in% rownames(lncRNA))) {
    n_success = n_success + 1
    for (degree in neighbor_degree) {
      neighbor_df_3_1 = rbind(neighbor_df_3_1, t(c(deg, degree, neighbors[[as.character(degree)]], neighbors_dist[[as.character(degree)]])))
    }
  } else {
    n_thrown = n_thrown+1
    print(deg)
    for (degree in neighbor_degree) {
      print(neighbors[[as.character(degree)]])
    }
    print("------------------------")
  }
}
colnames(neighbor_df_3_1) = c("deg", "degree", "neighbor", "distance")
used_lncRNA = unique(neighbor_df_3_1$deg)
print(length(as.vector(used_lncRNA)))

# Only sum neighbor genes in pos_cells for that lncRNA
# Not all genes will be summed in the same cell and not all cells will be summed for the same gene
tot_sum_tot  = lapply(1:length(neighbor_degree), function(x) setNames(rep(0, ncol(lncRNA)), colnames(lncRNA)))
tot_sum_list = lapply(1:length(neighbor_degree), function(x) setNames(rep(0, ncol(lncRNA)), colnames(lncRNA)))
tot_sum_neg  = lapply(1:length(neighbor_degree), function(x) setNames(rep(0, ncol(lncRNA)), colnames(lncRNA)))
w_wo_diff = data.frame()
i=0
for (gene in used_lncRNA) {
  i = i+1
  if (i %% 100 == 0) print(i)
  # get all the neighbors
  this_neighbors = neighbor_df_3_1$neighbor[which(neighbor_df_3_1$deg == gene)]
  
  pos_cells = lncRNA_cells[[gene]]
  not_cells = colnames(lncRNA)[which(! colnames(lncRNA) %in% pos_cells)]
  
  for (x in 1:length(neighbor_degree)) {
    tot_sum_tot[[x]]  = tot_sum_neg[[x]] + lncRNA@assays$RNA@counts[this_neighbors[x], ]
    tot_sum_list[[x]][pos_cells] = tot_sum_list[[x]][pos_cells] + lncRNA@assays$RNA@counts[this_neighbors[x], pos_cells]
    tot_sum_neg[[x]][not_cells]  = tot_sum_neg[[x]][not_cells] + lncRNA@assays$RNA@counts[this_neighbors[x], not_cells]
    
    w =  mean(lncRNA@assays$RNA@counts[this_neighbors[x], pos_cells]/lncRNA$nCount_RNA[pos_cells])
    wo = mean(lncRNA@assays$RNA@counts[this_neighbors[x], not_cells]/lncRNA$nCount_RNA[not_cells])
    w_wo_diff = rbind(w_wo_diff, t(c(neighbor_degree[x], w, wo, w-wo)))
  }
}
colnames(w_wo_diff) = c("degree", "w", "wo", "w-wo")
neighbor_df_3_2 = as.data.frame(lapply(1:length(neighbor_degree), function(x) tot_sum_list[[x]]/lncRNA$nCount_RNA))
colnames(neighbor_df_3_2) = neighbor_degree
neighbor_df_w = neighbor_df_3_2 %>% pivot_longer(colnames(neighbor_df_3_2), names_to = "degree", values_to = "mean")
neighbor_df_w$degree = factor(neighbor_df_w$degree, levels = neighbor_degree)
ggplot(neighbor_df_w, aes(x=degree, y=mean, group=degree, fill=degree, color=degree)) + geom_boxplot(alpha=0.9) + xlab("Neighbor to lncRNA Gene") + geom_jitter(position=position_jitter(), alpha=0.01) + NoLegend() + ylab("# Neighbor Transcripts/# All Transcripts") + labs(title="Non-zero Transcripts - Neighbor Summed if lncRNA Present")

neighbor_df_3_2 = as.data.frame(lapply(1:length(neighbor_degree), function(x) tot_sum_neg[[x]]/lncRNA$nCount_RNA))
colnames(neighbor_df_3_2) = neighbor_degree
neighbor_df_wo = neighbor_df_3_2 %>% pivot_longer(colnames(neighbor_df_3_2), names_to = "degree", values_to = "mean")
neighbor_df_wo$degree = factor(neighbor_df_wo$degree, levels = neighbor_degree)
ggplot(neighbor_df_wo, aes(x=degree, y=mean, group=degree, fill=degree, color=degree)) + geom_boxplot(alpha=0.9) + xlab("Neighbor to lncRNA Gene") + geom_jitter(position=position_jitter(), alpha=0.01) + NoLegend() + ylab("# Neighbor Transcripts/# All Transcripts") + labs(title="Non-zero Transcripts - Neighbor Summed if lncRNA NOT Present")

# for (degree in neighbor_degree) {
#   print(degree)
#   all_degree_neighbors = unique(neighbor_df_3_1$neighbor[which(neighbor_df_3_1$degree == degree)])
#   tot_sum = rep(0, ncol(lncRNA))
#   names(tot_sum) = colnames(lncRNA)
#   for (gene in used_lncRNA) {
#     not_cells = colnames(lncRNA)[which(! colnames(lncRNA) %in% lncRNA_cells[[gene]])]
#     pos_cells = lncRNA_cells[[gene]]
#     if (length(pos_cells) == 1) {
#       tot_sum[pos_cells] = sum(lncRNA@assays$RNA@counts[all_degree_neighbors,pos_cells])
#     } else {
#       cur_sum = colSums(lncRNA@assays$RNA@counts[all_degree_neighbors,])
#       cur_sum[which( !names(cur_sum) %in% pos_cells )] = 0
#       tot_sum = tot_sum + cur_sum
#     }
# 
#   }
#   newRow = data.frame(rep(degree, ncol(lncRNA)), tot_sum/lncRNA$nCount_RNA)
#   neighbor_df_3_2 = rbind(neighbor_df_3_2, newRow)
#   # neighbor_df_3_2 = rbind(neighbor_df_3_2, t(c(degree, mean(colSums(lncRNA@assays$RNA@counts[all_degree_neighbors,])/lncRNA$nCount_RNA))))
#   # neighbor_df_3_2 = rbind(neighbor_df_3_2, t(c(degree, mean(colSums(lncRNA@assays$RNA@data[all_degree_neighbors,])))))
# }
colnames(neighbor_df_3_2) <- c("degree", "mean")
neighbor_df_3_2$degree = factor(neighbor_df_3_2$degree, levels = neighbor_degree)
ggplot(neighbor_df_3_2, aes(x=degree, y=mean, group=degree, fill=degree, color=degree)) + geom_boxplot(alpha=0.9) + xlab("Neighbor to lncRNA Gene") + geom_jitter(position=position_jitter(), alpha=0.01) + NoLegend() + ylab("# Neighbor Transcripts/# All Transcripts") + labs(title="Non-zero Transcripts - Neighbor Summed if Normal Gene Present")

##############################################################
# See Corrlation between lncRNA expression and its Neighbors #
##############################################################
library("ggpmisc")
cor_p = c()
cor_cor = c()
for (n in 1:length(neighbor_degree)) {
# for (degree in neighbor_degree) {
  degree = neighbor_degree[n]
  print(degree)
  small_neighbor_df_3_2 = data.frame()
  
  # lncRNA Expression
  x = colSums(lncRNA@assays$RNA@counts[used_lncRNA,])/lncRNA$nCount_RNA
  
  # tot_sum = rep(0, ncol(lncRNA))
  # names(tot_sum) = colnames(lncRNA)
  # for (gene in used_lncRNA) {
  #   if (length(lncRNA_cells[[gene]]) == 1) {
  #     tot_sum[lncRNA_cells[[gene]]] = sum(lncRNA@assays$RNA@counts[all_degree_neighbors,lncRNA_cells[[gene]]])
  #   } else {
  #     cur_sum = colSums(lncRNA@assays$RNA@counts[all_degree_neighbors,])
  #     cur_sum[which( !names(cur_sum) %in% lncRNA_cells[[gene]] )] = 0
  #     tot_sum = tot_sum + cur_sum
  #   }
  # }
  
  # Neighbor Expression
  tot_sum = tot_sum_list[[n]]/lncRNA$nCount_RNA
  # all_degree_neighbors = unique(neighbor_df_3_1$neighbor[which(neighbor_df_3_1$degree == degree)])
  # y = colSums(lncRNA@assays$RNA@counts[all_degree_neighbors,])/lncRNA$nCount_RNA
  small_neighbor_df_3_2 = data.frame(x=x,y=tot_sum)
  
  cor_res = cor.test(x,tot_sum)
  cor_p   = c(cor_p,   cor_res$p.value)
  cor_cor = c(cor_cor, cor_res$estimate)
  
  # Plot
  my.formula <- y ~ x
  png(paste0("C:/Users/miles/Downloads/brain/results/ncRNA/cor_norm/", degree, ".png"), width = 800, height = 583)
  p = ggplot(small_neighbor_df_3_2, aes(x,y)) + geom_point(alpha=0.3) + xlab("lncRNA Expression") + ylab(paste(degree, "Neighbor Expression")) + labs(title=paste0("Corrlation: ", cor_res$estimate, ", P-value: ", cor_res$p.value)) + geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) + stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) 
  print(p)
  dev.off()
}
names(cor_p) = neighbor_degree
names(cor_cor) = neighbor_degree
sort(cor_p)
sort(cor_cor)

############################
# Observe Distance Effects #
############################
# My way
mean_exp = c()
mean_exp_lncRNA = c()
for (i in 1:nrow(neighbor_df_3_1)) {
  if (i %% 1000 == 0) { print(i) }
  this_lncRNA = neighbor_df_3_1$deg[i]
  this_neighbor = neighbor_df_3_1$neighbor[i]
  mean_exp = c(mean_exp, mean(lncRNA@assays$RNA@counts[this_neighbor,]/lncRNA$nCount_RNA))
  mean_exp_lncRNA = c(mean_exp_lncRNA, mean(lncRNA@assays$RNA@counts[this_lncRNA,]/lncRNA$nCount_RNA))
}

# Zack's Way
z_mean_exp = c()
z_mean_exp_lncRNA = c()
for (i in 1:nrow(neighbor_df_3_1)) {
  if (i %% 1000 == 0) { print(i) }
  this_lncRNA = neighbor_df_3_1$deg[i]
  this_neighbor = neighbor_df_3_1$neighbor[i]
  this_cells = lncRNA_cells[[this_lncRNA]]
  not_cells = colnames(lncRNA)[which(! colnames(lncRNA) %in% this_cells)]
  z_with           = mean(lncRNA@assays$RNA@counts[this_neighbor,this_cells]/lncRNA$nCount_RNA[this_cells])
  z_without        = mean(lncRNA@assays$RNA@counts[this_neighbor, not_cells]/lncRNA$nCount_RNA[not_cells ])
  z_with_lncRNA    = mean(lncRNA@assays$RNA@counts[this_lncRNA,  this_cells]/lncRNA$nCount_RNA[this_cells])
  z_without_lncRNA = mean(lncRNA@assays$RNA@counts[this_lncRNA,   not_cells]/lncRNA$nCount_RNA[not_cells ])
  z_mean_exp = c(z_mean_exp, z_with-z_without)
  z_mean_exp_lncRNA = c(z_mean_exp_lncRNA, z_with_lncRNA - z_without_lncRNA)
}

neighbor_df_3_1$mean_exp = mean_exp
neighbor_df_3_1$mean_exp_lncRNA = mean_exp_lncRNA
neighbor_df_3_1$distance = as.numeric(as.vector(neighbor_df_3_1$distance))
neighbor_df_3_1$abs_dist = abs(neighbor_df_3_1$distance)
neighbor_df_3_1$log_dist = log(neighbor_df_3_1$abs_dist, base = 2)
neighbor_df_3_1$my_log_dist = log(abs(neighbor_df_3_1$distance), base = 2) * sign(neighbor_df_3_1$distance)
neighbor_df_3_1$z_mean_exp = z_mean_exp
neighbor_df_3_1$z_mean_exp_lncRNA = z_mean_exp_lncRNA

ggplot(neighbor_df_3_1, aes(mean_exp_lncRNA, mean_exp, color = log_dist)) + geom_point(alpha=0.3) + ylab("Neighbor Mean Expression Across All Cells") + xlab("lncRNA Mean Expression Across All Cells") + labs(title="Distance Effects") +scale_color_gradientn(colours = rainbow(5))
ggplot(neighbor_df_3_1, aes(abs_dist, mean_exp)) + geom_point(alpha=0.3) + ylab("Neighbor Mean Expression Across All Cells") + xlab("Absolute Distance From lncRNA") + labs(title="Distance Effects")+ geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) + stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE)
ggplot(neighbor_df_3_1, aes(abs_dist, z_mean_exp, color = degree)) + geom_point(alpha=0.3) + ylab("Neighbor Mean Expression Across All Cells") + xlab("Absolute Distance From lncRNA") + labs(title="Distance Effects")

for (abs_degree in unique(abs(neighbor_degree))) {
  new_range = -abs_degree:abs_degree
  new_range = new_range[which(new_range != 0)]
  png(paste0("C:/Users/miles/Downloads/brain/results/ncRNA/dist/z_", abs_degree, ".png"))
  print(ggplot(neighbor_df_3_1[which(neighbor_df_3_1$degree %in% new_range),], aes(abs_dist, z_mean_exp)) + geom_point(alpha=0.3) + ylab("Neighbor Mean Expression Across All Cells") + xlab("Absolute Distance From lncRNA") + labs(title=paste("Distance Effects From", -abs_degree, "to", abs_degree)) + geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) + stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE))
  dev.off()
}

for (degree in neighbor_degree) {
  png(paste0("C:/Users/miles/Downloads/brain/results/ncRNA/dist/", degree, ".png"))
  print(ggplot(neighbor_df_3_1[which(neighbor_df_3_1$degree == degree),], aes(mean_exp_lncRNA, mean_exp, color = abs_dist)) + geom_point(alpha=0.3) + ylab("Neighbor Mean Expression Across All Cells") + xlab("lncRNA Mean Expression Across All Cells") + labs(title="Distance Effects") +scale_color_gradientn(colours = rainbow(5)), limits = c(0,max(neighbor_df_3_1$abs_dist)))
  dev.off()
}

#################
# Left vs Right #
#################
neighbor_df_3_3 = w_wo_diff
# for (degree in neighbor_degree) {
#   print(degree)
#   all_degree_neighbors = unique(neighbor_df_3_1$neighbor[which(neighbor_df_3_1$degree == degree)])
#   for (gene in used_lncRNA) {
#     not_cells = colnames(lncRNA)[which(! colnames(lncRNA) %in% lncRNA_cells[[gene]])]
#     if (length(lncRNA_cells[[gene]]) == 1) {
#       mean_with    = mean(sum(lncRNA@assays$RNA@counts[all_degree_neighbors, lncRNA_cells[[gene]]])/lncRNA$nCount_RNA[lncRNA_cells[[gene]]])
#       mean_without = mean(colSums(lncRNA@assays$RNA@counts[all_degree_neighbors, not_cells])/lncRNA$nCount_RNA[not_cells])
#       neighbor_df_3_3 = rbind(neighbor_df_3_3, t(c(gene, degree, mean_with, mean_without, mean_with - mean_without)))
#     } else {
#       mean_with    = mean(colSums(lncRNA@assays$RNA@counts[all_degree_neighbors, lncRNA_cells[[gene]]])/lncRNA$nCount_RNA[lncRNA_cells[[gene]]])
#       mean_without = mean(colSums(lncRNA@assays$RNA@counts[all_degree_neighbors, not_cells])/lncRNA$nCount_RNA[not_cells])
#       neighbor_df_3_3 = rbind(neighbor_df_3_3, t(c(gene, degree, mean_with, mean_without, mean_with - mean_without)))
#     }
#   }
# }
colnames(neighbor_df_3_3) <- c("degree", "with", "without",  "diff_mean")
neighbor_df_3_3$degree = factor(neighbor_df_3_3$degree, levels = neighbor_degree)
neighbor_df_3_3$diff_mean = as.numeric(as.vector(neighbor_df_3_3$diff_mean))
ggplot(neighbor_df_3_3, aes(x=degree, y=diff_mean, group=degree, fill=degree, color=degree)) + geom_boxplot(alpha=0.9) + xlab("Neighbor to lncRNA Gene") + geom_jitter(position=position_jitter(), alpha=0.03) + NoLegend() + ylab("Difference in Means b/w Cells With and Without the lncRNA") + labs(title="Non-zero Transcripts - Difference in Means")

# Plot w/wo seperately
neighbor_df_3_3$degree = as.numeric(as.vector(neighbor_df_3_3$degree))
df2 = neighbor_df_3_3 %>% pivot_longer(c(with, without), names_to = "has_lncRNA", values_to = "mean_expression")
df2$mean_expression = as.numeric(as.vector(df2$mean_expression))
df2$degree = factor(df2$degree, levels=neighbor_degree)
ggplot(df2, aes(x=degree, y=mean_expression, fill=has_lncRNA, color=has_lncRNA)) + geom_jitter(position=position_jitterdodge(), alpha=0.03) + geom_boxplot(alpha=0.9) + xlab("Neighbor to lncRNA Gene") + ylab("Difference in Means b/w Cells With and Without the lncRNA") + labs(title="Non-zero Transcripts - Difference in Means - With and Without Seperated")
