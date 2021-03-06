#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Test for correct input arguments: if not, return an error
if (length(args)==0) {
  stop(paste("No arguments supplied.", 
             "Please supply the permutation number and number of permutation to complete",
             "Or at least just the permutation number."), call.=FALSE)
} else if (length(args)==1) {
  # default number of permutations
  args[2] = 10
}

# Read Inputs
perm_num = as.numeric(args[1]) # The current permutation number. This is used for the seed.
num_perm = as.numeric(args[2]) # The number of permutations to complete.
print(paste("Permutation Number:", perm_num))
output_folder = "~/scratch/brain/results/cor_pr_perm/"

# Load Packages
library("Seurat")
library("Matrix")
library("reshape2")
library("igraph")
library("parallel")
numCores = detectCores()

# Load Data
bb = readRDS("~/scratch/brain/data/bb_subsample_02222021.RDS")
mat = bb@assays$RNA@data

# Set random seed so all the permutations are different
set.seed(perm_num)

# Permute BHVE and CTRL. Split into 2 matrices.
print(paste("Permuting Data", num_perm, "times."))
perm_labels = lapply(1:num_perm, function(x) sample(bb$cond))
# b_mats = lapply(1:num_perm, function(x) t(as.matrix(bb@assays$RNA@data[,which(perm_labels[[x]] == "BHVE")])))
# c_mats = lapply(1:num_perm, function(x) t(as.matrix(bb@assays$RNA@data[,which(perm_labels[[x]] == "CTRL")])))
# all_mats = append(b_mats, c_mats)
# b_mats = NULL # clear memory
# c_mats = NULL # clear memory

bb = NULL #clear memory

# Slightly Parallelized BHVE vs CTRL Permutation Pageranks
pagerank_dif = list() # differences between behave and control in all pairs
pageranks = list() 
for (i in 1:num_perm) {
  print(paste0("Finding Correlations in Pair ", i, "."))
  all_mats = list(t(as.matrix(mat[,which(perm_labels[[i]] == "BHVE")])),
                  t(as.matrix(mat[,which(perm_labels[[i]] == "CTRL")])))
  cor_mats <- mclapply(all_mats[c(1,2)], function(mat) cor(mat, y = NULL), mc.cores = numCores, mc.preschedule = TRUE)
  all_mats = NULL # clear memory
  
  # B Prep
  print(paste("Cleaning B in Pair", i))
  cor_mats[[1]] = cor_mats[[1]][which( ! is.na(cor_mats[[1]][2,]) ), which( ! is.na(cor_mats[[1]][2,]) )]
  print(paste0("Dimensions of Clean Matrix B in Pair ", i, ": ", dim(cor_mats[[1]])))
  # b_melt = setNames(melt(cor_mats[[1]]), c("Node1", "Node2", "weight"))
  
  # C Prep
  print(paste("Cleaning C in Pair", i))
  cor_mats[[2]] = cor_mats[[2]][which( ! is.na(cor_mats[[2]][2,]) ), which( ! is.na(cor_mats[[2]][2,]) )]
  print(paste0("Dimensions of Clean Matrix C in Pair ", i, ": ", dim(cor_mats[[2]])))
  # c_melt = setNames(melt(cor_mats[[2]]), c("Node1", "Node2", "weight"))
  
  # cor_mats = NULL # clear Memory
  
  # # B Graph + Pagerank
  # print(paste("Creating graph B in Pair", i))
  # graph_obj = graph_from_data_frame(b_melt)
  # print("Finding pagerank of each node in graph B.")
  # pr_b = page.rank(graph_obj)$vector
  # b_melt = NULL # clear memory
  # graph_obj = NULL # clear memory
  # 
  # # C Graph + Pagerank
  # print(paste("Creating graph C in Pair", i))
  # graph_obj = graph_from_data_frame(c_melt)
  # print("Finding pagerank of each node in graph C.")
  # pr_c = page.rank(graph_obj)$vector
  # c_melt = NULL # clear memory
  # graph_obj = NULL #clear memory
  
  # Node Strength
  print(paste("Finding Node Strength of B in Pair", i))
  pr_b = colSums(abs(cor_mats[[1]]))
  print(paste("Finding Node Strength of C in Pair", i))
  pr_c = colSums(abs(cor_mats[[2]]))
  
  pr_df = t(plyr::ldply(list(pr_b, pr_c), rbind))
  pr_df = as.data.frame(pr_df)
  pr_df[is.na(pr_df)] = 0
  b_c_dif = pr_df$V1 - pr_df$V2
  names(b_c_dif) = rownames(pr_df)
  pageranks[[i]] = b_c_dif
  
  print(paste0("Finished ", i))
}

# Write the PageRank for all the permutated matrices to file
print("Writing to File.")
pr_df = t(plyr::ldply(pageranks, rbind))
colnames(pr_df) = c(paste0("dif_", perm_num, "_", 1:num_perm))
write.csv(pr_df, paste0("~/scratch/brain/results/cor_pr_perm/perm_", perm_num, ".csv"))
print("Done.")

# Parallelized Version - Too much memory.
# Permute BHVE and CTRL. Split into 2 matrices.
# print(paste("Permuting Data", num_perm, "times."))
# perm_labels = lapply(1:num_perm, function(x) sample(bb$cond))
# b_mats = lapply(1:num_perm, function(x) t(as.matrix(bb@assays$RNA@data[,which(perm_labels[[x]] == "BHVE")])))
# c_mats = lapply(1:num_perm, function(x) t(as.matrix(bb@assays$RNA@data[,which(perm_labels[[x]] == "CTRL")])))
# all_mats = append(b_mats, c_mats)
# b_mats = NULL # clear memory
# c_mats = NULL # clear memory
# Find the correlations for all permutations
# print("Finding Correlations.")
# numCores = detectCores()
# cor_mats <- mclapply(all_mats, function(mat) cor(mat, y = NULL), mc.cores = numCores)
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
# # Write the PageRank for all the permutated matrices to file
# print("Writing to File.")
# pr_df = t(plyr::ldply(pageranks, rbind))
# colnames(pr_df) = c(paste0("b", perm_num, ".", 1:10), paste0("c", perm_num, ".", 1:10))
# write.csv(pr_df, paste0("~/scratch/brain/results/cor_pr_perm/perm_", perm_num, ".csv"))
# print("Done.")