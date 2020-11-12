######################
# Enrichment Testing #
######################
# This is the final script for enrichment testing (testing for significant 
# enrichment of a list of genes of interest in single cell clusters).
# This script uses a consensus of two appproaches: down_and_perm_2 and psi. 
#
# 1. down_and_perm_2
# First, the data is downsampled to the minimum number of transcripts per cell.
# This means that transcripts for every cell is randomly discarded until it reaches
# the same number of transcripts in the cell that has the minimum number of transcripts.
# This ensures that every cell has the same expression level, thereby bypassing the problem
# of some cells/clusters have more transcripts on average than another cell/cluster.
# After downsampling, store the average number of genes in the list of interest expressed 
# per cell in each cluster (this is a list of a size equal to the number of clusters).
# Then take the average of this number for each cluster after 50 downsamples.
#
# Repeat this process for downsampled and permutated data. To permutate the data,
# new clusters of equal size as the original clusters are formed by randomly choosing a
# cluster and then randomly choosing a cell within that cluster to be in the new cluster.
#
# Finally, the average of the average number of genes in the list of interest expressed 
# per cell in each cluster in the downsampled emperical data is compared to 99.5th 
# percentile of the average number of genes in the list of interest expressed per 
# cell in each cluster in the downsampled permuted data.
#
# 2. psi
# https://cran.r-project.org/web/packages/pSI/pSI.pdf
# "This package contains functions to calculate the Specificity Index statistic, 
# which can be used for comparative quantitative analysis to identify genes enriched 
# in specific cell populations across a large number of profiles." 
#
# 3. Finding consensus of the methods
# Finally the clusters that are found significant by both methods are returned and
# are considered significantly enriched for the list of genes of interest.