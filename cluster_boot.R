# Read Input ===================================================================
args = commandArgs(trailingOnly=TRUE)
boot.num = as.numeric(args[1])
min_dist = as.numeric(args[2])
n_neighbors = as.numeric(args[3])
res   = as.numeric(args[4])
nboot = as.numeric(args[5])
boot.pct = 0.8
num.cores = 10

# Load Libraries ===============================================================
suppressMessages(library('parallel',  quietly = T, warn.conflicts = F, verbose = F))
suppressMessages(library('Seurat',  quietly = T, warn.conflicts = F, verbose = F))

# Load Object ==================================================================
obj = readRDS("~/scratch/brain/data/hb_100322_diet.rds")
num.boot.cells = ncol(obj) * boot.pct
orig.labels = read.csv(paste0("~/scratch/brain/results/hb/cluster_boot/hb_mindist", min.dist, "_nneighbors", n_neighbors, "_res", res, ".csv"))

# Helper Functions =============================================================
sc.cluster.boot = function(x) {
  this.obj = obj
  this.obj = subset(this.obj, cells = sample(colnames(this.obj), num.boot.cells, replace = F))
  
  # Clustering Parameters
  this.obj = SCTransform(this.obj, vars.to.regress = "sample", verbose = F)
  this.obj = RunPCA(this.obj, dim = 50, verbose = F)
  this.obj = RunUMAP(this.obj, dims=1:50, min.dist=min_dist, spread=1, n.neighbors=n_neighbors, n.epochs=1000, metric="euclidean")
  this.obj = FindNeighbors(this.obj, reduction="umap", k.param=n_neighbors, dims=1:2, n.trees=500, prune.SNN=0)
  this.obj = FindClusters(this.obj, resolution = res, algorithm = 2)
  this.labels = this.obj$seurat_clusters
  
  return(as.numeric(as.vector(this.labels)))
}

evaluate.boot = function(x) {
  orig.and.boot = paste0(orig.labels, "_", boot.labels[,x])
  orig.and.boot.count = data.frame(table(orig.and.boot))
  orig.and.boot.count[,c("orig", "boot")] = reshape2::colsplit(orig.and.boot.count[,1], "_", c('1', '2'))
  jcrd.mat = reshape2::acast(orig.and.boot.count, orig ~ boot, value.var = "Freq")
  cluster.max.jcrd = unlist(lapply(1:nrow(jcrd.mat), function(x) max(jcrd.mat)))
  return(cluster.max.jcrd)
}

# Main Body ====================================================================
# boot.labels is a data.frame with rows corresponding to cells and columns 
# corresponding to bootstraps. The value of each element is the cluster the 
# cell was assigned to in that bootstap.
boot.labels = as.data.frame(mclapply(1:nboot, function(x) sc.cluster.boot(x), mc.cores = num.cores))
write.csv(boot.labels, "~/scratch/brain/results/my_cluster_stability/test_boot_labels.csv") # TODO

# cluster.jaccard is a data.frame with rows corresponding to clusters in the original
# clustering. While the columns are the 
cluster.jaccard.list = mclapply(1:nboot, function(x) evaluate.boot(x), mc.cores = num.cores)
cluster.jaccard = do.call('cbind', cluster.jaccard.list)
write.csv(cluster.jaccard, "~/scratch/brain/results/my_cluster_stability/test_cluster_jaccard.csv") # TODO
cluster.jaccard.melt = reshape2::melt(cluster.jaccard) # TODO make sure melt works
colnames(cluster.jaccard.melt) = c("cluster", "boot", "value") # TODO is this correct name order
pdf(paste0("~/scratch/brain/results/my_cluster_stability/test_cluster_jaccard.pdf"), width = 12, height = 6)
print(ggplot(cluster.jaccard.melt, aes(x = cluster, y = value, color = cluster, fill = cluster)) + geom_boxplot(alpha = 0.5) + geom_point() + NoLegend())
dev.off()
