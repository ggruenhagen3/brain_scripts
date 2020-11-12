library("dplyr")
library("Matrix")
library("Seurat")
library("stringr")
library("ggplot2")
library("fgsea")

degs <- read.table(paste0(rna_path, "/results/b1b2c1mz_degs_full.tsv"), sep="\t", header = TRUE, stringsAsFactors = FALSE)

marker_path <- paste(rna_path, "data/markers/", sep="")
marker_files <- dir(marker_path, pattern =paste("*.txt", sep=""))
markers <- data.frame(gene <- c(), bio <- c())
for (i in 1:length(marker_files)) {
  file <- read.table(paste(marker_path, marker_files[i], sep=""), header=FALSE, sep="\t", stringsAsFactors=FALSE)
  markers <- rbind(markers, file[,1:2])
}
colnames(markers) <- c("gene", "bio")
bio <- "ROCK_SAND"
markers <- markers[which(markers$bio == bio),]

fgseaRes <- data.frame()
num_clusters <- max(as.numeric(as.vector(degs$cluster)))
for (i in 0:num_clusters) {
  this_degs <- degs[which(degs$cluster == i),c("gene", "avg_logFC")]
  rank <- setNames(this_degs$avg_logFC, this_degs$gene)
  pathway <- list()
  pathway[["ROCK_SAND"]] <- markers$gene
  this_res <- fgsea(pathway, rank, nperm = 1000)
  this_res$cluster <- i
  fgseaRes <- rbind(fgseaRes, this_res)
}
fgseaRes$q <- p.adjust(fgseaRes$padj)

# ClustifyR
remotes::install_github("rnabioco/clustifyr")
library("clustifyr")
res <- run_gsea(combined@assays$RNA@data, markers$gene, cluster_ids = combined$seurat_clusters)