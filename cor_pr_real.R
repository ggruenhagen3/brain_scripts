# Load Packages
library("Seurat")
library("Matrix")
library("reshape2")
library("igraph")
library("parallel")

# Load Data
bb = readRDS("~/scratch/brain/data/bb_subsample_02222021.RDS")

# Split the Data
b_mat = t(as.matrix(bb@assays$RNA@data[,which(bb$cond == "BHVE")]))
c_mat = t(as.matrix(bb@assays$RNA@data[,which(bb$cond == "CTRL")]))
all_mats = list(b_mat, c_mat)

# Find the correlations
print("Finding Correlations.")
numCores = detectCores()
cor_mats <- mclapply(all_mats, function(mat) cor(mat, y = NULL), mc.cores = numCores)

# Prepare the Data for Graph Creation and Pagerank
cor_mats <- mclapply(cor_mats, function(mat) mat[which( ! is.na(mat[2,]) ), which( ! is.na(mat[2,]) )], mc.cores = numCores)  # removes NAs
print("Dimensions of Clean Correlation Matrices:")
junk = sapply(cor_mats, function(mat) print(dim(mat)))
melt_mats <-  mclapply(cor_mats, function(mat) setNames(melt(mat), c("Node1", "Node2", "weight")), mc.cores = numCores)

# Create the graphs
print("Creating the graphs.")
graph_objs <- mclapply(melt_mats, function(x) graph_from_data_frame(x), mc.cores = numCores)

# Find the PageRank
print("Finding pagerank of each node in the graph.")
pageranks = mclapply(graph_objs, function(g) page.rank(g)$vector, mc.cores = numCores)

# Write the PageRank BHVE and CTRL matrices to file
print("Writing to File.")
pr_df = t(plyr::ldply(pageranks, rbind))
colnames(pr_df) = c(paste0("b_real"), paste0("c_real"))
write.csv(pr_df, paste0("~/scratch/brain/results/cor_pr_real.csv"))
print("Done.")