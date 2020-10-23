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
mt_genes = scan("C:/Users/miles/Downloads/all_research/mt.txt", what = "character")
mt_genes = str_replace(mt_genes, "_", "-")
mt_genes = mt_genes[which(mt_genes %in% rownames(bb))]
bb$pct_mt = PercentageFeatureSet(bb, features = mt_genes)
FeatureScatter(bb, feature1 = "nCount_RNA", feature2 = "pct_mt") + NoLegend()
FeatureScatter(bb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
print(paste("Number of Cells in BB Before Filterning:", ncol(bb)))

bb = subset(bb, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & pct_mt < 5)
print(paste("Number of Cells in BB After Filterning:", ncol(bb)))
bb = NormalizeData(bb, normalization.method = "LogNormalize", scale.factor = 100000)
bb = FindVariableFeatures(bb, selection.method = "vst", nfeatures = 2000)
bb = ScaleData(bb, features = rownames(bb))
bb = RunPCA(bb, features = VariableFeatures(bb))
bb = RunUMAP(bb, reduction = "pca", dims = 1:20)
bb = FindNeighbors(bb, reduction="umap", dims = 1:2)
bb = FindClusters(bb, resolution = .15)
DimPlot(bb, reduction = "umap", split.by = "cond", label = TRUE)

Cairo(file = paste0(rna_path, "/results/bb/umap.png"), width = 2000, height = 1500, res = 200)
print(DimPlot(bb, reduction = "umap", split.by = "cond", label = TRUE, pt.size = 0.7))
dev.off()

# bb = JackStraw(bb, num.replicate = 100)
# bb = ScoreJackStraw(bb, dims = 1:20)

#==========================================================================================
# Markers =================================================================================
#==========================================================================================
bri = list()
tmp = data.frame(readxl::read_xlsx("C:/Users/miles/Downloads/brain/data/MC,MZ Merged Gene, Cluster Info.xlsx", sheet = "Glia", skip=1, col_types = "text"))
bri[["glia"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp[,2])),"\\s")))
tmp = data.frame(readxl::read_xlsx("C:/Users/miles/Downloads/brain/data/MC,MZ Merged Gene, Cluster Info.xlsx", sheet = "Vc, Vd", skip=6, col_types = "text"))
bri[["vc_vd"]] = gsub("[\r\n;]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx("C:/Users/miles/Downloads/brain/data/MC,MZ Merged Gene, Cluster Info.xlsx", sheet = "Vs-Vp, Vv-Vl", skip=5, col_types = "text"))
bri[["vs_vp_vv_vl"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx("C:/Users/miles/Downloads/brain/data/MC,MZ Merged Gene, Cluster Info.xlsx", sheet = "OB (gcl)", skip=5, col_types = "text"))
bri[["ob_gcl"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx("C:/Users/miles/Downloads/brain/data/MC,MZ Merged Gene, Cluster Info.xlsx", sheet = "OB (mcl, gl)", skip=9, col_types = "text"))
bri[["ob_mcl_gl"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx("C:/Users/miles/Downloads/brain/data/MC,MZ Merged Gene, Cluster Info.xlsx", sheet = "BNSM", skip=0, col_types = "text"))
bri[["bnsm"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx("C:/Users/miles/Downloads/brain/data/MC,MZ Merged Gene, Cluster Info.xlsx", sheet = "Dp", skip=11, col_types = "text"))
bri[["dp"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx("C:/Users/miles/Downloads/brain/data/MC,MZ Merged Gene, Cluster Info.xlsx", sheet = "Dm", skip=9, col_types = "text"))
bri[["dm"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))
tmp = data.frame(readxl::read_xlsx("C:/Users/miles/Downloads/brain/data/MC,MZ Merged Gene, Cluster Info.xlsx", sheet = "Dl", skip=9, col_types = "text"))
bri[["dl"]] = gsub("[\r\n]", "", unlist(strsplit(as.character(as.vector(tmp$Gene)),"\\s")))

bri_clean = sapply(1:length(bri), function(x) unique(bri[[x]][which(bri[[x]] %in% rownames(bb))]))
names(bri_clean) = names(bri)
for (i in 1:length(bri_clean)) {
  markers = bri_clean[[i]]
  name = names(bri_clean)[i]
  
  # UMAP Expression
  Cairo(file=paste0(rna_path, "/results/bb/paintings/", name, ".png"), width = 750, height = 500)
  print(markerExpPerCell(bb, markers) + ggtitle(paste0("Expression of ", str_to_title(name))))
  dev.off()
  
  # DotPlot
  dot_res = myDotPlot(bb, markers)
  Cairo(file=paste0(rna_path, "/results/bb/paintings/", name, "_dot.png"), width = 750, height = 500)
  print(dot_res[[1]] + ggtitle(paste0("Expression of ", str_to_title(name))))
  dev.off()
}

#==========================================================================================
# BHVE vs CTRL Bootstrap ==================================================================
#==========================================================================================
# Purpose: Determine if there's a real BHVE vs CTRL signal by
#          comparing the real BHVE vs CTRL to simulated ones.

# Load object called obj
tj = readRDS("~/scratch/d_tooth/data/tj.rds") # for now load tj
tj$sample = rep(1:10, ncol(tj))[1:ncol(tj)]
obj = tj

# Constants
n_boot = 100 # number of bootstraps
do_replace = F # sample with replacement?
bhve_samples = c(1, 2, 3, 4, 5)
ctrl_sampels = c(6, 7, 8, 9, 10)

# Find the sizes of the samples
print("Finding the sizes of the samples")
samples = unique(as.vector(obj$sample))
sample_cells = list()
sizes = list()
for (sample in samples) {
  sample_cells[[sample]] = obj$sample[which(obj$sample == sample)]
  sizes[[sample]] = length(sample_cells[[sample]])
}

# # Testing a non-complicated sampling without replacement
# ind = match(boot_samples$cell[which(boot_samples$boot == 1)], colnames(obj))
# obj$test_sample = boot_samples$sample[ind]
# obj$cond = obj$test_sample %in% bhve_samples
# obj$cond = factor(obj$cond, levels = c(TRUE, FALSE))
# obj$cond = plyr::revalue(obj$cond, replace = c("TRUE" = "bhve", "FALSE" = "ctrl"))
# Idents(obj) = obj$cond
# test = FindAllMarkers(obj_2, assay = "RNA", only.pos = F, verbose = F)

# Generate random samples
print("Generating random samples")
boot_samples = data.frame()
for (i in 1:n_boot) {
  # Some Useful Print Statements
  if(i == n_boot) {
    cat(paste(i, "\n"))
  } else if (i %% (n_boot/10) == 0 || i == 1) {
    cat(i)
  } else {
    cat(".")
  }
  
  # ran_samples_cells = list() # keys: sim samples 1-10 , values: cells in that sim sample
  cells_left = colnames(obj) # cells left to sample
  
  for (sample in samples) {
    ran_samples_cells = sample(cells_left, sizes[[sample]], replace = do_replace)
    
    # if sample w/o replacement, take the used cells out
    if (! do_replace) { 
      cells_left = cells_left[which(! cells_left %in% ran_samples_cells)]
    } # end replacement if
    
    new_rows = data.frame(rep(i, sizes[[sample]]), rep(sample, sizes[[sample]]), ran_samples_cells)
    boot_samples = rbind(boot_samples, new_rows)
  } # end sample for
  
  # Error check for sampling without replacement
  if (length(cells_left > 0) && do_replace == F) {
    print(paste("Error: cells left after sampling without replacement.",
                "There should be 0 cells left, but there are", length(cells_left), "left."))
  } # end error check for
} # end boot for
cat("\n")
colnames(boot_samples) = c("boot", "sample", "cell")

# Find DEGs in between sim BHVE and CTRL
sim_degs = data.frame()
for (i in 1:n_boot) {
  # Some Useful Print Statements
  if(i == n_boot) {
    cat(paste(i, "\n"))
  } else if (i %% (n_boot/10) == 0 || i == 1) {
    cat(i)
  } else {
    cat(".")
  }
  
  all_cells = boot_samples$cell[which(boot_samples$boot == i)]
  
  # Create a new Seurat object 
  obj_2 = CreateSeuratObject(counts = obj@assays$RNA@counts[,all_cells], project = obj@project.name)
  data_matrix = obj@assays$RNA@data[,all_cells]
  data_matrix@Dimnames[[2]] = colnames(obj_2)
  obj_2 = SetAssayData(object = obj_2, slot = 'data', assay = "RNA", new.data = data_matrix)
  
  # Add metadata
  obj_2$sim_sample = boot_samples$sample[which(boot_samples$boot ==i)]
  obj_2$cond = obj_2$sim_sample %in% bhve_samples
  obj_2$cond = factor(obj_2$cond, levels = c(TRUE, FALSE))
  obj_2$cond = plyr::revalue(obj_2$cond, replace = c("TRUE" = "bhve", "FALSE" = "ctrl"))
  
  # Find DEGs
  Idents(obj_2) = obj_2$cond
  obj_2_deg = FindAllMarkers(obj_2, assay = "RNA", only.pos = F, verbose = F)
  obj_2_deg = obj_2_deg[which(obj_2_deg$p_val_adj < 0.05),]
  if (nrow(obj_2_deg) > 0) obj_2_deg$boot = i
  
  rbind(sim_degs, obj_2_deg)
} # end boot for


# Test Cluster DEG vs Bulk DEG for bhve v ctrl
obj$seurat_clusters_cond = paste0(obj$seurat_clusters, "_", obj$cond)
Idents(obj) = obj$seurat_clusters_cond
# bhve_cells = names(obj$cond[which(obj$cond == "BHVE")])
cluster_deg = data.frame()
for (cluster in sort(unique(obj$seurat_clusters))) {
  print(cluster)
  # this_cluster = names(obj$seurat_clusters[which(obj$seurat_clusters == cluster)])
  # this_cluster_bhve = this_cluster[which(this_cluster %in% bhve_cells)]
  # this_cluster_ctrl = this_cluster[which(! this_cluster %in% bhve_cells)]
  
  # this_deg = FindMarkers(obj, cells.1 = this_cluster_bhve, cells.2 = this_cluster_ctrl, ident.1 = NULL, ident.2 = NULL)
  # cells_1 = WhichCells(obj, idents = paste0(cluster, "_BHVE"))
  # cells_2 = WhichCells(obj, idents = paste0(cluster, "_CTRL"))
  tryCatch ({
    this_deg = FindMarkers(obj, paste0(cluster, "_BHVE"), paste0(cluster, "_CTRL"), only.pos = F)
    this_deg = this_deg[which(this_deg$p_val_adj < 0.05),]
    cluster_deg = rbind(cluster_deg, this_deg)
  }, error = function(err) {
    print("Failed on this cluster")
  })


}
nrow(cluster_deg)
length(unique(rownames(cluster_deg)))

Idents(obj) = obj$cond
bulk_deg = FindAllMarkers(obj, only.pos = F)
bulk_deg = bulk_deg[which(bulk_deg$p_val_adj < 0.05),]
nrow(bulk_deg)/2
length(unique(bulk_deg$gene))

cluster_deg$q = p.adjust(cluster_deg$p_val, method="bonferroni")
nrow(cluster_deg[which(cluster_deg$q < 0.05),])
