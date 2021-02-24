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
perm_num = args[1] # The current permutation number. This is used for the seed.
num_perm = args[2] # The number of permutations to complete.
print(paste("Permutation Number:", perm_num))
output_folder = "~/scratch/brain/results/cor_pr_perm/"

# Load Packages
library("Seurat")
library("Matrix")
library("reshape2")
library("igraph")
library("parallel")

# Load Data
bb = readRDS("~/scratch/brain/data/bb_subsample_02222021.RDS")

# Set random seed so all the permutations are different
set.seed(perm_num)

# Permute BHVE and CTRL. Split into 2 matrices.
print(paste("Permuting Data", num_perm, "times."))
perm_labels = lapply(1:num_perm, function(x) sample(c("BHVE", "CTRL"), nrow(bb), replace = T))
b_mats = lapply(1:num_perm, function(x) as.matrix(bb@assays$RNA@data[,which(perm_labels[[x]] == "BHVE")]))
c_mats = lapply(1:num_perm, function(x) as.matrix(bb@assays$RNA@data[,which(perm_labels[[x]] == "CTRL")]))
all_mats = append(b_mats, c_mats)
b_mats = NULL # clear memory
c_mats = NULL # clear memory

# Find the correlations for all permutations
print("Finding Correlations.")
numCores = detectCores()
cor_mats <- mclapply(all_mats, function(mat) cor(mat, y = NULL), mc.cores = numCores)

# Prepare the Data for Graph Creation and Pagerank
cor_mats <- mclapply(cor_mats, function(mat) mat[which( ! is.na(mat[1,]) ), which( ! is.na(mat[1,]) )], mc.cores = numCores)  # removes NAs
junk = sapply(cor_mats, function(mat) print(dim(mat)))
melt_mats <-  mclapply(cor_mats, function(mat) setNames(melt(mat), c("Node1", "Node2", "weight")), mc.cores = numCores)

# Create the graphs
print("Creating the graphs.")
graph_objs <- mclapply(melt_mats, function(x) graph_from_data_frame(x), mc.cores = numCores)

# Find the PageRank
print("Finding pagerank of each node in the graph.")
pageranks = mclapply(graph_objs, function(g) page.rank(g)$vector, mc.cores = numCores)

# Write the PageRank for all the permutated matrices to file
print("Writing to File.")
pr_df = setNames(as.data.frame(pageranks), c(paste0("b", perm_num, ".", 1:10), paste0("c", perm_num, ".", 1:10)))
write.csv(pr_df, paste0("~/scratch/brain/results/cor_pr_perm/perm_", perm_num, ".csv"))
print("Done.")