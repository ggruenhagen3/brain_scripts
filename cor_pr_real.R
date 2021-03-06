# Load Packages
library("Seurat")
library("Matrix")
library("reshape2")
library("igraph")
library("parallel")
numCores = detectCores()

# Load Data
bb = readRDS("~/scratch/brain/data/bb_subsample_02222021.RDS")

# Split the Data
b_mat = t(as.matrix(bb@assays$RNA@data[,which(bb$cond == "BHVE")]))
c_mat = t(as.matrix(bb@assays$RNA@data[,which(bb$cond == "CTRL")]))
all_mats = list(b_mat, c_mat)
b_mat = NULL # clear memory
c_mat = NULL # clear memory

# Slightly Parallelized BHVE vs CTRL Permutation Pageranks
print(paste0("Finding Correlations in Pair."))
cor_mats <- mclapply(all_mats[c(1,2)], function(mat) cor(mat, y = NULL), mc.cores = numCores, mc.preschedule = TRUE)
all_mats = NULL # clear memory

# Save Data
saveRDS(cor_mats[[1]], "~/scratch/brain/data/bb_b_cor.RDS")
saveRDS(cor_mats[[2]], "~/scratch/brain/data/bb_c_cor.RDS")
# Find p values
my_cor_t = function(r, n) (r * sqrt(n - 2))/sqrt(1 - r**2)
my_cor_p = function(t, n) 2*pt(-abs(t), df=n-2)
b_t_mat = my_cor_t(cor_mats[[1]], ncol(cor_mats[[1]]))
b_p_mat = my_cor_p(b_t_mat, ncol(cor_mats[[1]]))
c_t_mat = my_cor_t(cor_mats[[2]], ncol(cor_mats[[2]]))
c_p_mat = my_cor_p(c_t_mat, ncol(cor_mats[[2]]))
saveRDS(b_p_mat,   "~/scratch/brain/data/bb_b_cor_p.RDS")
saveRDS(c_p_mat,   "~/scratch/brain/data/bb_c_cor_p.RDS")
b_t_mat = b_p_mat = c_t_mat = c_p_mat = NULL

# B Prep
print(paste("Cleaning B in Pair"))
cor_mats[[1]] = cor_mats[[1]][which( ! is.na(cor_mats[[1]][2,]) ), which( ! is.na(cor_mats[[1]][2,]) )]
print(paste0("Dimensions of Clean Matrix B in Pair: ", dim(cor_mats[[1]])))
# b_melt = setNames(melt(cor_mats[[1]]), c("Node1", "Node2", "weight"))

# C Prep
print(paste("Cleaning C in Pair"))
cor_mats[[2]] = cor_mats[[2]][which( ! is.na(cor_mats[[2]][2,]) ), which( ! is.na(cor_mats[[2]][2,]) )]
print(paste0("Dimensions of Clean Matrix C in Pair: ", dim(cor_mats[[2]])))
# c_melt = setNames(melt(cor_mats[[2]]), c("Node1", "Node2", "weight"))

# cor_mats = NULL # clear Memory

# # B Graph + Pagerank
# print(paste("Creating graph B in Pair"))
# graph_obj = graph_from_data_frame(b_melt)
# print("Finding pagerank of each node in graph B.")
# pr_b = page.rank(graph_obj)$vector
# b_melt = NULL # clear memory
# graph_obj = NULL # clear memory
# 
# # C Graph + Pagerank
# print(paste("Creating graph C in Pair"))
# graph_obj = graph_from_data_frame(c_melt)
# print("Finding pagerank of each node in graph C.")
# pr_c = page.rank(graph_obj)$vector
# c_melt = NULL # clear memory
# graph_obj = NULL # clear memory

# Node Strength
print(paste("Finding Node Strength of B"))
pr_b = colSums(abs(cor_mats[[1]]))
print(paste("Finding Node Strength of C"))
pr_c = colSums(abs(cor_mats[[2]]))

pr_df = t(plyr::ldply(list(pr_b, pr_c), rbind))
pr_df = as.data.frame(pr_df)
pr_df[is.na(pr_df)] = 0
pr_df$Dif = pr_df$V1 - pr_df$V2
b_c_dif = pr_df$Dif
names(b_c_dif) = rownames(pr_df)

print(paste0("Finished Hard Part."))

# Write the PageRank for all the permutated matrices to file
print("Writing to File.")
colnames(pr_df) = c("B", "C", "Dif")
write.csv(pr_df, paste0("~/scratch/brain/results/cor_pr_real_strength.csv"))
print("Done.")


# # Find the correlations
# print("Finding Correlations.")
# numCores = detectCores()
# cor_mats <- mclapply(all_mats, function(mat) cor(mat, y = NULL), mc.cores = numCores)
# 
# # Prepare the Data for Graph Creation and Pagerank
# cor_mats <- mclapply(cor_mats, function(mat) mat[which( ! is.na(mat[2,]) ), which( ! is.na(mat[2,]) )], mc.cores = numCores)  # removes NAs
# print("Dimensions of Clean Correlation Matrices:")
# junk = sapply(cor_mats, function(mat) print(dim(mat)))
# melt_mats <-  mclapply(cor_mats, function(mat) setNames(melt(mat), c("Node1", "Node2", "weight")), mc.cores = numCores)
# 
# # Create the graphs
# print("Creating the graphs.")
# graph_objs <- mclapply(melt_mats, function(x) graph_from_data_frame(x), mc.cores = numCores)
# 
# # Find the PageRank
# print("Finding pagerank of each node in the graph.")
# pageranks = mclapply(graph_objs, function(g) page.rank(g)$vector, mc.cores = numCores)
# 
# # Write the PageRank BHVE and CTRL matrices to file
# print("Writing to File.")
# pr_df = t(plyr::ldply(pageranks, rbind))
# colnames(pr_df) = c(paste0("b_real"), paste0("c_real"))
# write.csv(pr_df, paste0("~/scratch/brain/results/cor_pr_real.csv"))
# print("Done.")