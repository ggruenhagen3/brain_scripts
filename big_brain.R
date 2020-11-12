# Load libraries
library("Seurat")
library("Matrix")
library("reticulate")
library("cowplot")
library("biomaRt")
library("stringr")
library("dplyr")
library("ggplot2")
rna_path = "C:/Users/miles/Downloads/brain/"
# rna_path = "~/scratch/brain/"

#==========================================================================================
# Initial Clustering ======================================================================
#==========================================================================================
# Load Filtered Feature Matrix
b1.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-B1/"))
b2.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-B2/"))
b3.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-B3/"))
b4.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-B4/"))
b5.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-B5/"))
c1.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-C1/"))
c2.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-C2/"))
c3.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-C3/"))
c4.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-C4/"))
c5.data = Read10X(data.dir = paste0(rna_path, "data/bb_ffm/JTS07-C5/"))

b1 = CreateSeuratObject(counts = b1.data)
b2 = CreateSeuratObject(counts = b2.data)
b3 = CreateSeuratObject(counts = b3.data)
b4 = CreateSeuratObject(counts = b4.data)
b5 = CreateSeuratObject(counts = b5.data)
c1 = CreateSeuratObject(counts = c1.data)
c2 = CreateSeuratObject(counts = c2.data)
c3 = CreateSeuratObject(counts = c3.data)
c4 = CreateSeuratObject(counts = c4.data)
c5 = CreateSeuratObject(counts = c5.data)

# Add Metdata
b1$sample = "b1"
b2$sample = "b2"
b3$sample = "b3"
b4$sample = "b4"
b5$sample = "b5"
c1$sample = "c1"
c2$sample = "c2"
c3$sample = "c3"
c4$sample = "c4"
c5$sample = "c5"

b1$cond = "bhve"
b2$cond = "bhve"
b3$cond = "bhve"
b4$cond = "bhve"
b5$cond = "bhve"
c1$cond = "ctrl"
c2$cond = "ctrl"
c3$cond = "ctrl"
c4$cond = "ctrl"
c5$cond = "ctrl"

bb = merge(b1, list(b2,b3,b4,b5,c1,c2,c3,c4,c5), add.cell.ids = c("b1", "b2", "b3", "b4", "b5", "c1", "c2", "c3", "c4", "c5"))

rm(b1.data); rm(b1)
rm(b2.data); rm(b2)
rm(b3.data); rm(b3)
rm(b4.data); rm(b4)
rm(b5.data); rm(b5)
rm(c1.data); rm(c1)
rm(c2.data); rm(c2)
rm(c3.data); rm(c3)
rm(c4.data); rm(c4)
rm(c5.data); rm(c5)

# Quality Control: Mitochondrial DNA, Dead Cells, and Doublets
# mt_genes = scan("~/scratch/m_zebra_ref/mt.txt", what = "character")
mt_genes = scan("C:/Users/miles/Downloads/all_research/mt.txt", what = "character")
mt_genes = str_replace(mt_genes, "_", "-")
mt_genes = mt_genes[which(mt_genes %in% rownames(bb))]
bb$pct_mt = PercentageFeatureSet(bb, features = mt_genes)
FeatureScatter(bb, feature1 = "nCount_RNA", feature2 = "pct_mt") + NoLegend()
FeatureScatter(bb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
print(paste("Number of Cells in BB Before Filterning:", ncol(bb)))

bb = subset(bb, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & pct_mt < 5)
print(paste("Number of Cells in BB After Filterning:", ncol(bb)))
bb = NormalizeData(bb, normalization.method = "LogNormalize", scale.factor = 10000)
bb = FindVariableFeatures(bb, selection.method = "vst", nfeatures = 2000)
bb = ScaleData(bb, features = rownames(bb))
bb = RunPCA(bb, features = VariableFeatures(bb))
bb = RunUMAP(bb, reduction = "pca", dims = 1:50)
bb = FindNeighbors(bb, reduction="umap", dims = 1:2)
bb = FindClusters(bb, resolution = .30)
# DimPlot(bb, reduction = "umap", split.by = "cond", label = TRUE)

png_name = paste0(rna_path, "/results/bb/umap_10k_50_30.png")
png(file = png_name, width = 2000, height = 1500, res = 200)
print(DimPlot(bb, reduction = "umap", label = TRUE, pt.size = 0.7))
dev.off()
system(paste0("rclone copy ", png_name, " dropbox:BioSci-Streelman/George/Brain/bb/tmp/"))

#==========================================================================================
# Final Object ============================================================================
#==========================================================================================
rna_path = "C:/Users/miles/Downloads/brain/"
bb = readRDS(paste0(rna_path, "data/bb_clustered_102820.rds"))

png_name = "C:/Users/miles/Downloads/bb_pc_fst_90_umap.png"
png(file = png_name, width = 1000, height = 1000, res = 120)
print(p)
dev.off()

png_name = "C:/Users/miles/Downloads/90_15.png"
png(file = png_name, width = 1500, height = 1000, res = 200)
print(res[[1]])
dev.off()

#==========================================================================================
# Markers =================================================================================
#==========================================================================================
bri = list()
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "Glia", skip=1, col_types = "text"))
bri[["glia"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp[,2])),"\\s")))
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "Vc, Vd", skip=6, col_types = "text"))
bri[["vc_vd"]] = gsub("[\r\n;]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "Vs-Vp, Vv-Vl", skip=5, col_types = "text"))
bri[["vs_vp_vv_vl"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "OB (gcl)", skip=5, col_types = "text"))
bri[["ob_gcl"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "OB (mcl, gl)", skip=9, col_types = "text"))
bri[["ob_mcl_gl"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "BNSM", skip=0, col_types = "text"))
bri[["bnsm"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "Dp", skip=11, col_types = "text"))
bri[["dp"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "Dm", skip=9, col_types = "text"))
bri[["dm"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx(paste0(rna_path, "data/MC,MZ Merged Gene, Cluster Info.xlsx"), sheet = "Dl", skip=9, col_types = "text"))
bri[["dl"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))

bri_clean = sapply(1:length(bri), function(x) unique(bri[[x]][which(bri[[x]] %in% rownames(bb))]))
names(bri_clean) = names(bri)
for (i in 1:length(bri_clean)) {
  markers = bri_clean[[i]]
  name = names(bri_clean)[i]
  
  # UMAP Expression
  png(file=paste0(rna_path, "/results/bb/paintings/", name, ".png"), width = 750, height = 500)
  print(markerExpPerCell(bb, markers) + ggtitle(paste0("Expression of ", str_to_title(name))))
  dev.off()
  
  # DotPlot
  dot_res = myDotPlot(bb, markers)
  png(file=paste0(rna_path, "/results/bb/paintings/", name, "_dot.png"), width = 750, height = 750)
  print(dot_res[[1]] + ggtitle(paste0("Expression of ", str_to_title(name))))
  dev.off()
}

#==========================================================================================
# BHVE vs CTRL Bootstrap ==================================================================
#==========================================================================================
# Purpose: Determine if there's a real BHVE vs CTRL signal by
#          comparing the real BHVE vs CTRL to simulated ones.

# Load object called obj
obj = bb

# Constants
n_boot = 50 # number of bootstraps
do_replace = F # sample with replacement?
bhve_samples = c("b1", "b2", "b3", "b4", "b5")
ctrl_samples = c("c1", "c2", "c3", "c4", "c5")
clusters = sort(unique(obj$seurat_clusters))
all_samples = c(bhve_samples, ctrl_samples)
bhve_cells = colnames(obj)[which(obj$sample %in% bhve_samples)]
ctrl_cells = colnames(obj)[which(obj$sample %in% ctrl_samples)]
obj$seurat_clusters_vect = as.vector(obj$seurat_clusters)
names(obj$seurat_clusters_vect) = colnames(obj)

# Generate random samples
print("Generating random samples and Finding DEGs")
sim_degs = data.frame()
for (i in 1:n_boot) {
  set.seed(i+100)
  # Some Useful Print Statements
  if(i == n_boot) {
    cat(paste(i, "\n"))
  } else if (i %% (n_boot/10) == 0 || i == 1) {
    cat(i)
  } else {
    cat(".")
  }

  # Shuffle Clusters while keeping the # cells per sample per cluster the same
  obj$sim = obj$seurat_clusters_vect
  Idents(obj) = obj$seurat_clusters
  for (cluster in clusters) {
    this_cells = WhichCells(obj, idents = cluster, seed = NULL)
    sim = sample(obj$cond[this_cells], length(this_cells))
    obj$sim[this_cells] = unname(sim)
  }
  
  # Find DEGs
  obj$sim_cond = paste0(obj$seurat_clusters, obj$sim)
  Idents(obj) = obj$sim_cond
  for (cluster in clusters) {
    degs = FindMarkers(obj, ident.1 = paste0(cluster, "BHVE"), ident.2 = paste0(cluster, "CTRL"), only.pos=F, verbose = F, min.pct = 0.01, logfc.threshold = 0.1)
    if ( ! is.null(degs) && nrow(degs) > 0 ) {
      degs$boot = i
      degs$cluster = cluster
      sim_degs = rbind(sim_degs, degs)
    }
  }
} # end boot for
# sim_degs = sim_degs[which(sim_degs$boot <= 50),]
write.table(sim_degs, paste0(rna_path, "results/bb/sim_bhve_v_ctrl.tsv"), sep = "\t", quote = F, row.names = F)

perm_df2 = data.frame()
for (i in 1:n_boot) {
  degs = sim_degs[which(sim_degs$boot == i),]
  qobj = qvalue(degs$p_val)
  degs$q = qobj$qvalues
  degs$bh = p.adjust(degs$p_val, method = "BH")
  
  n_p = length(which(degs$p_val < 0.05))
  n_q = length(which(degs$q < 0.05))
  n_bh = length(which(degs$bh < 0.05))
  n_bon = length(which(degs$p_val_adj < 0.05))
  
  perm_df2 = rbind(perm_df2, t(c(i, n_p, n_q, n_bh, n_bon)))
}

real_value_p = length(which(zack$p_val < 0.05))
real_value_q = length(which(zack$q < 0.05))
real_value_bh = length(which(zack$bh < 0.05))
real_value_bon = length(which(zack$p_val_adj < 0.05))
colnames(perm_df2) = c("Boot", "n_p", "n_q", "n_bh", "n_bon")

perm_df2$above = perm_df2$n_bon > real_value_bon
ggplot(perm_df2, aes(n_bon, alpha=.7, fill=above)) + geom_histogram(alpha=0.5) + geom_vline(aes(xintercept = real_value_bon)) + geom_text(aes(x=real_value_bon, label="Real Value"), y = Inf, hjust=0, vjust=1, color = "black") + xlab("# of Gene w/ adjusted p-value < 0.05 (Bonferroni)") + ggtitle("Comparison Between Bootstrap Values and Real Value") + guides(color=F, alpha=F, fill=F)

perm_df2$above = perm_df2$n_bh > real_value_bh
ggplot(perm_df2, aes(n_bh, alpha=.7, fill=above)) + geom_histogram(alpha=0.5) + geom_vline(aes(xintercept = real_value_bh)) + geom_text(aes(x=real_value_bh, label="Real Value"), y = Inf, hjust=0, vjust=1, color = "black") + xlab("# of Gene w/ adjusted FDR adjusted pvalue < 0.05 (BH)") + ggtitle("Comparison Between Bootstrap Values and Real Value") + guides(color=F, alpha=F, fill=F)

perm_df2$above = perm_df2$n_q > real_value_q
ggplot(perm_df2, aes(n_q, alpha=.7, fill=above)) + geom_histogram(alpha=0.5) + geom_vline(aes(xintercept = real_value_q)) + geom_text(aes(x=real_value_q, label="Real Value"), y = Inf, hjust=0, vjust=1, color = "black") + xlab("# of Gene w/ qvalue < 0.05") + ggtitle("Comparison Between Bootstrap Values and Real Value") + guides(color=F, alpha=F, fill=F)


perm_df2$above = perm_df2$n_p > real_value_p
ggplot(perm_df2, aes(n_p, alpha=.7, fill=above)) + geom_histogram(alpha=0.5) + geom_vline(aes(xintercept = real_value_p)) + geom_text(aes(x=real_value_p, label="Real Value"), y = Inf, hjust=0, vjust=1, color = "black") + xlab("# of Gene w/ pvalue < 0.05") + ggtitle("Comparison Between Bootstrap Values and Real Value") + guides(color=F, alpha=F, fill=F)


# Test Cluster DEG vs Bulk DEG for bhve v ctrl
obj$seurat_clusters_cond = paste0(obj$seurat_clusters, "_", obj$cond)
Idents(obj) = obj$seurat_clusters_cond
# bhve_cells = names(obj$cond[which(obj$cond == "BHVE")])
cluster_deg = data.frame()
for (cluster in sort(unique(obj$seurat_clusters))) {
  print(cluster)
  tryCatch ({
    this_deg = FindMarkers(obj, paste0(cluster, "_BHVE"), paste0(cluster, "_CTRL"), min.pct = 0.01, logfc.threshold = 0.1)
    this_deg$cluster = cluster
    this_deg$gene = rownames(this_deg)
    # this_deg = this_deg[which(this_deg$p_val_adj < 0.05),]
    cluster_deg = rbind(cluster_deg, this_deg)
  }, error = function(err) {
    print("Failed on this cluster")
  })


}
nrow(cluster_deg)
length(unique(rownames(cluster_deg)))

Idents(obj) = obj$cond
bulk_deg = FindAllMarkers(obj, only.pos = F)
# bulk_deg = bulk_deg[which(bulk_deg$p_val_adj < 0.05),]
nrow(bulk_deg)/2
length(unique(bulk_deg$gene))

cluster_deg$bonferroni = p.adjust(cluster_deg$p_val, method="bonferroni")
cluster_deg$bh = p.adjust(cluster_deg$p_val, method="BH")
qobj = qvalue(cluster_deg$p_val)
cluster_deg$qvalue = qobj$qvalues
nrow(cluster_deg[which(cluster_deg$q < 0.05),])
