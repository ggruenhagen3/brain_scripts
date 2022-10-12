# Read Input ===================================================================
# boot.num = 1; min_dist = 0.1; n_neighbors = 10; res = 0.2; nboot = 1;
args = commandArgs(trailingOnly=TRUE)
boot.num = as.numeric(args[1])
min_dist = as.numeric(args[2])
n_neighbors = as.numeric(args[3])
res   = as.numeric(args[4])
nboot = as.numeric(args[5])
boot.pct = 0.8
num.cores = 10
message(paste0("Running w/ Parameters: boot.num=", boot.num, ", min_dist=", min_dist, ", n_neighbors=", n_neighbors, ", res=", res, ", nboot=", nboot))

# Load Libraries ===============================================================
suppressMessages(library('parallel',  quietly = T, warn.conflicts = F, verbose = F))
suppressMessages(library('Seurat',  quietly = T, warn.conflicts = F, verbose = F))
suppressMessages(library('ggplot2',  quietly = T, warn.conflicts = F, verbose = F))

# Load Object ==================================================================
obj = readRDS("~/scratch/brain/data/hb_100322_diet.rds")
num.boot.cells = round(ncol(obj@assays$RNA@counts) * boot.pct)
all.orig.labels = read.csv("~/scratch/brain/results/hb/hb_clustering_options101222.csv")
orig.labels = all.orig.labels[,paste0("mindist", min_dist, "_nneighbors", n_neighbors, "_res", res)]

# Helper Functions =============================================================
sc.cluster.boot = function(x) {
  message("Subsampling")
  this.cells = sample(colnames(obj@assays$RNA@counts), num.boot.cells, replace = F)
  this.obj = CreateSeuratObject(counts = obj@assays$RNA@counts[,this.cells], meta.data = obj@meta.data[this.cells,])
  
  # Clustering Parameters
  message("SCTransform")
  this.obj = SCTransform(this.obj, vars.to.regress = "sample", verbose = F)
  message("RunPCA")
  this.obj = RunPCA(this.obj, dim = 50, verbose = F)
  message("RunUMAP")
  this.obj = RunUMAP(this.obj, dims=1:50, min.dist=min_dist, spread=1, n.neighbors=n_neighbors, n.epochs=1000, metric="euclidean", verbose = F)
  message("FindNeighbors")
  this.obj = FindNeighbors(this.obj, reduction="umap", k.param=n_neighbors, dims=1:2, n.trees=500, prune.SNN=0, verbose = FALSE)
  message("FindClusters")
  this.obj = FindClusters(this.obj, resolution = res, algorithm = 2, verbose = F)
  this.labels = this.obj$seurat_clusters
  
  return(as.numeric(as.vector(this.labels)))
}

evaluate.boot = function(x) {
  orig.and.boot = paste0(orig.labels, "_", boot.labels[,x])
  orig.and.boot.count = data.frame(table(orig.and.boot))
  orig.and.boot.count[,c("orig", "boot")] = reshape2::colsplit(orig.and.boot.count[,1], "_", c('1', '2'))
  jcrd.mat = reshape2::acast(orig.and.boot.count, orig ~ boot, value.var = "Freq")
  jcrd.mat[is.na(jcrd.mat)] = 0
  cluster.max.jcrd = unlist(lapply(1:nrow(jcrd.mat), function(x) max(jcrd.mat[x,])))
  cluster.max.jcrd = cluster.max.jcrd / as.numeric(as.vector(table(orig.labels)))
  return(cluster.max.jcrd)
}

# Main Body ====================================================================
# boot.labels is a data.frame with rows corresponding to cells and columns 
# corresponding to bootstraps. The value of each element is the cluster the 
# cell was assigned to in that bootstap.
boot.labels = as.data.frame(mclapply(1:nboot, function(x) sc.cluster.boot(x), mc.cores = num.cores))
colnames(boot.labels) = 1:nboot
write.csv(boot.labels, paste0("~/scratch/brain/results/hb/cluster_boot/boots/mindist", min_dist, "_nneighbors", n_neighbors, "_res", res, "_labels.csv"))

# cluster.jaccard is a data.frame with rows corresponding to clusters in the original
# clustering. While the columns are the bootstraps. Each element is the highest jaccard index
# of the original cluster in any of the bootstrap clusters from that bootstrap.
cluster.jaccard.list = mclapply(1:nboot, function(x) evaluate.boot(x), mc.cores = num.cores)
cluster.jaccard = do.call('cbind', cluster.jaccard.list)
write.csv(cluster.jaccard, paste0("~/scratch/brain/results/hb/cluster_boot/boots/mindist", min_dist, "_nneighbors", n_neighbors, "_res", res, "_jcrd.csv"))

# Visualize the Results
cluster.jaccard.melt = reshape2::melt(cluster.jaccard)
cluster.jaccard.melt[,1] = cluster.jaccard.melt[,1] - 1
colnames(cluster.jaccard.melt) = c("cluster", "boot", "value") # TODO is this correct name order
pdf(paste0("~/scratch/brain/results/hb/cluster_boot/boots/mindist", min_dist, "_nneighbors", n_neighbors, "_res", res, "_jcrd.pdf"), width = 12, height = 6)
print(ggplot(cluster.jaccard.melt, aes(x = cluster, y = value, color = cluster, fill = cluster, group = cluster)) + geom_point() + NoLegend())
dev.off()
