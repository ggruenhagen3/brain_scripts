# Count num cells in clusters
cells_in_cluster <- c()
for (i in 0:num_clusters) {
  this_cells <- WhichCells(combined, idents = i)
  cells_in_cluster <- c(cells_in_cluster, length(this_cells))
}

# Find the number of cells expressing that gene in each cluster
combined@active.assay <- "RNA"
pct_df <- data.frame()
for (gene in ran) {
  expr1 <- FetchData(object = combined, vars = gene, slot = "counts")
  gene_cells <- colnames(combined@assays$RNA@counts[, which(x = expr1 > 0)])
  gene_table  <- table(combined$seurat_clusters[gene_cells])
  new_row <- c(gene, "RAN", as.numeric(as.vector(gene_table)))
  names(new_row) <- c("gene", "type", 0:40)
  pct_df <- rbind(pct_df, t(new_row))
  colnames(pct_df) <- c("gene", "type", 0:40)
}
for (gene in rs) {
  expr1 <- FetchData(object = combined, vars = gene, slot = "counts")
  gene_cells <- colnames(combined@assays$RNA@counts[, which(x = expr1 > 0)])
  gene_table  <- table(combined$seurat_clusters[gene_cells])
  new_row <- c(gene, "ROCK_SAND", as.numeric(as.vector(gene_table)))
  names(new_row) <- c("gene", "type", 0:40)
  pct_df <- rbind(pct_df, t(new_row))
  colnames(pct_df) <- c("gene", "type", 0:40)
}

# Find the percentage of cells expressing that gene in each cluster
pct_df[,3:43] <- sapply(3:43, function(x) as.numeric(pct_df[,x])/cells_in_cluster)

# Find the number of genes from ran vs rock-sand expressed in each cluster
pct_df_ran <- pct_df[which(pct_df$type == "RAN"),]
pct_df_ran$num_clusters <- sapply(1:nrow(pct_df_ran), function(x) length(which(pct_df_ran[x,3:43] > 0) ))
gene_names_per_cluster3 <- sapply(1:nrow(pct_df_ran), function(x) which(pct_df_ran[x,3:43] > 0 ))
cluster_table_ran <- table(unlist(gene_names_per_cluster3, recursive = FALSE))

pct_df_rs <- pct_df[which(pct_df$type == "ROCK_SAND"),]
pct_df_rs$num_clusters <- sapply(1:nrow(pct_df_rs), function(x) length(which(pct_df_rs[x,3:43] > 0) ))
gene_names_per_cluster3 <- sapply(1:nrow(pct_df_rs), function(x) which(pct_df_rs[x,3:43] > 0 ))
cluster_table_rs <- table(unlist(gene_names_per_cluster3, recursive = FALSE))

cluster_df <- rbind(as.data.frame(cluster_table_ran), as.data.frame(cluster_table_rs))
cluster_df$type <- cbind(c(rep("RAN", 41), rep("ROCK_SAND", 41)))
colnames(cluster_df) <- c("Cluster", "Number_of_Genes", "Gene_List")
ggplot(cluster_df, aes(x=Cluster, y=Number_of_Genes, color=Gene_List, fill=Gene_List)) + geom_bar(stat = "identity", position=position_dodge(), alpha=0.5) + ggtitle("Number of Genes from RAN vs Rock-Sand Expressed in Each Cluster")

# Find the number of genes from ran vs rock-sand expressed above a threshold in each cluster
pct_df_ran <- pct_df[which(pct_df$type == "RAN"),]
pct_df_ran$num_clusters <- sapply(1:nrow(pct_df_ran), function(x) length(which(pct_df_ran[x,3:43] > 0.05) ))
gene_names_per_cluster3 <- sapply(1:nrow(pct_df_ran), function(x) which(pct_df_ran[x,3:43] > 0.05 ))
cluster_table_ran <- table(unlist(gene_names_per_cluster3, recursive = FALSE))

pct_df_rs <- pct_df[which(pct_df$type == "ROCK_SAND"),]
pct_df_rs$num_clusters <- sapply(1:nrow(pct_df_rs), function(x) length(which(pct_df_rs[x,3:43] > 0.05) ))
gene_names_per_cluster3 <- sapply(1:nrow(pct_df_rs), function(x) which(pct_df_rs[x,3:43] > 0.05 ))
cluster_table_rs <- table(factor(unlist(gene_names_per_cluster3, recursive = FALSE), levels=1:41))

cluster_df <- rbind(as.data.frame(cluster_table_ran), as.data.frame(cluster_table_rs))
cluster_df$type <- cbind(c(rep("RAN", 41), rep("ROCK_SAND", 41)))
colnames(cluster_df) <- c("Cluster", "Number_of_Genes", "Gene_List")
ggplot(cluster_df, aes(x=Cluster, y=Number_of_Genes, color=Gene_List, fill=Gene_List)) + geom_bar(stat = "identity", position=position_dodge(), alpha=0.5) + ggtitle("Number of Genes from RAN vs Rock-Sand Expressed Above 5% in Each Cluster")

# Keep genes that are at thresh in 1-10 clusters
threshold <- 0.05
pct_df_ran <- pct_df[which(pct_df$type == "RAN"),]
pct_df_ran$num_clusters <- sapply(1:nrow(pct_df_ran), function(x) length(which(pct_df_ran[x,3:43] > threshold) ))
pct_df_ran <- pct_df_ran[which(pct_df_ran$num_clusters > 0 & pct_df_ran$num_clusters < 11),]
gene_names_per_cluster3 <- sapply(1:nrow(pct_df_ran), function(x) which(pct_df_ran[x,3:43] > threshold ))
cluster_table_ran <- table(factor(unlist(gene_names_per_cluster3, recursive = FALSE), levels = 1:41))

pct_df_rs <- pct_df[which(pct_df$type == "ROCK_SAND"),]
pct_df_rs$num_clusters <- sapply(1:nrow(pct_df_rs), function(x) length(which(pct_df_rs[x,3:43] > threshold) ))
pct_df_rs <- pct_df_rs[which(pct_df_rs$num_clusters > 0 & pct_df_rs$num_clusters < 11),]
gene_names_per_cluster3 <- sapply(1:nrow(pct_df_rs), function(x) which(pct_df_rs[x,3:43] > threshold ))
cluster_table_rs <- table(factor(unlist(gene_names_per_cluster3, recursive = FALSE), levels=1:41))

cluster_df <- rbind(as.data.frame(cluster_table_ran), as.data.frame(cluster_table_rs))
cluster_df$type <- cbind(c(rep("RAN", 41), rep("ROCK_SAND", 41)))
colnames(cluster_df) <- c("Cluster", "Number_of_Genes", "Gene_List")
ggplot(cluster_df, aes(x=Cluster, y=Number_of_Genes, color=Gene_List, fill=Gene_List)) + geom_bar(stat = "identity", position=position_dodge(), alpha=0.5) + ggtitle("Number of Genes from RAN vs Rock-Sand Expressed Above 5% in Each Cluster and is in 1-10 Clusters")

###############
# Transcripts #
###############
ran_sums <- colSums(as.matrix(combined@assays$RNA@counts[ran,]))
ran_sums <- as.data.frame(ran_sums)
ran_sums$cluster <- combined$seurat_clusters
cluster_trans_ran <- as.vector(sapply(0:40, function(x) sum(ran_sums$ran_sums[which(ran_sums$cluster == x)])))
cluster_trans_ran <- data.frame(cluster_trans_ran, rep("RAN", 41))
cluster_trans_ran$cluster <- 0:40
colnames(cluster_trans_ran) <- c("trans", "type", "cluster")

ran_sums <- colSums(as.matrix(combined@assays$RNA@counts[rs,]))
ran_sums <- as.data.frame(ran_sums)
ran_sums$cluster <- combined$seurat_clusters
cluster_trans_rs <- as.vector(sapply(0:40, function(x) sum(ran_sums$ran_sums[which(ran_sums$cluster == x)])))
cluster_trans_rs <- data.frame(cluster_trans_rs, rep("ROCK_SAND", 41))
cluster_trans_rs$cluster <- 0:40
colnames(cluster_trans_rs) <- c("trans", "type", "cluster")

cluster_trans_df <- rbind(cluster_trans_ran, cluster_trans_rs)
ggplot(cluster_trans_df, aes(x=cluster, y=trans, color=type, fill=type)) + geom_bar(stat = "identity", position=position_dodge(), alpha=0.5) + ggtitle("Number of Transcripts for RAN vs Rock-Sand Expressed in Each Cluster")

# Trans above a threshold
trans_df <- data.frame()
for (gene in ran) {
  this_trans <- setNames(data.frame(combined@assays$RNA@counts[gene,]), "trans")
  this_trans$cluster <- combined$seurat_clusters
  this_trans_cluster <- aggregate(this_trans$trans, by=list(Category=this_trans$cluster), FUN=sum)
  trans_df <- rbind(trans_df, t(c(gene, "RAN", as.numeric(as.vector(this_trans_cluster$x)))))
}
for (gene in rs) {
  this_trans <- setNames(data.frame(combined@assays$RNA@counts[gene,]), "trans")
  this_trans$cluster <- combined$seurat_clusters
  this_trans_cluster <- aggregate(this_trans$trans, by=list(Category=this_trans$cluster), FUN=sum)
  trans_df <- rbind(trans_df, t(c(gene, "ROCK_SAND", as.numeric(as.vector(this_trans_cluster$x)))))
}
colnames(trans_df) <- c("gene", "type", 0:40)
trans_df[,3:43] <- t(sapply(1:nrow(trans_df), function(x) as.numeric(as.vector(trans_df[x,3:43]))))

pct_df_ran <- pct_df[which(pct_df$type == "RAN"),]
pct_df_ran$num_clusters <- sapply(1:nrow(pct_df_ran), function(x) length(which(pct_df_ran[x,3:43] > 0.05) ))
pct_df_ran_genes <- pct_df_ran$gene[which(pct_df_ran$num_clusters > 0)]
trans_df_ran <- trans_df[which(trans_df$type == "RAN"),]
# trans_df_ran <- trans_df[which(trans_df$type == "RAN" & trans_df$gene %in% pct_df_ran_genes),]
sum_trans_df_ran <- data.frame(0:40, colSums(trans_df_ran[,3:43]), rep("RAN", 41))
colnames(sum_trans_df_ran) <- c("cluster", "trans", "type")

pct_df_rs <- pct_df[which(pct_df$type == "ROCK_SAND"),]
pct_df_rs$num_clusters <- sapply(1:nrow(pct_df_rs), function(x) length(which(pct_df_rs[x,3:43] > 0.05) ))
pct_df_rs_genes <- pct_df_rs$gene[which(pct_df_rs$num_clusters > 0)]
trans_df_rs <- trans_df[which(trans_df$type == "ROCK_SAND"),]
# trans_df_rs <- trans_df[which(trans_df$type == "ROCK_SAND" & trans_df$gene %in% pct_df_rs_genes),]
sum_trans_df_rs <- data.frame(0:40, colSums(trans_df_rs[,3:43]), rep("ROCK_SAND", 41))
colnames(sum_trans_df_rs) <- c("cluster", "trans", "type")

cluster_trans_df <- rbind(sum_trans_df_ran, sum_trans_df_rs)
ggplot(cluster_trans_df, aes(x=cluster, y=trans, color=type, fill=type)) + geom_bar(stat = "identity", position=position_dodge(), alpha=0.5) + ggtitle("Number of Transcripts for RAN vs Rock-Sand Expressed Above 5% in Each Cluster")

# Trans above threshold plus 1-10
pct_df_ran <- pct_df[which(pct_df$type == "RAN"),]
pct_df_ran$num_clusters <- sapply(1:nrow(pct_df_ran), function(x) length(which(pct_df_ran[x,3:43] > 0.05) ))
pct_df_ran_genes <- pct_df_ran$gene[which(pct_df_ran$num_clusters > 0 & pct_df_ran$num_clusters < 11)]
trans_df_ran <- trans_df[which(trans_df$type == "RAN" & trans_df$gene %in% pct_df_ran_genes),]
sum_trans_df_ran <- data.frame(0:40, colSums(trans_df_ran[,3:43]), rep("RAN", 41))
colnames(sum_trans_df_ran) <- c("cluster", "trans", "type")

pct_df_rs <- pct_df[which(pct_df$type == "ROCK_SAND"),]
pct_df_rs$num_clusters <- sapply(1:nrow(pct_df_rs), function(x) length(which(pct_df_rs[x,3:43] > 0.05) ))
pct_df_rs_genes <- pct_df_rs$gene[which(pct_df_rs$num_clusters > 0 & pct_df_rs$num_clusters < 11)]
trans_df_rs <- trans_df[which(trans_df$type == "ROCK_SAND" & trans_df$gene %in% pct_df_rs_genes),]
sum_trans_df_rs <- data.frame(0:40, colSums(trans_df_rs[,3:43]), rep("ROCK_SAND", 41))
colnames(sum_trans_df_rs) <- c("cluster", "trans", "type")

cluster_trans_df <- rbind(sum_trans_df_ran, sum_trans_df_rs)
ggplot(cluster_trans_df, aes(x=cluster, y=trans, color=type, fill=type)) + geom_bar(stat = "identity", position=position_dodge(), alpha=0.5) + ggtitle("Number of Transcripts for RAN vs Rock-Sand Expressed Above 5% in Each Cluster and in 1-10 Clusters")


ran <- backup_markers[which(backup_markers$bio == "RAN"),1]
rs <- backup_markers[which(backup_markers$bio == "ROCK_SAND"),1]
num_clust <- data.frame()
for (gene in ran) {
  num_clust_gene <- 0
  for (i in 0:num_clusters) {
    this_cells <- WhichCells(combined, idents = i)
    num_cells <- length(which(as.vector(combined@assays$RNA@counts[gene,this_cells]) != 0))
    if( num_cells/length(this_cells) > threshold ) {
      num_clust_gene <- num_clust_gene + 1
    }
  }
  num_clust <- rbind(num_clust, t(c(gene, "RAN", num_clust_gene)))
  # unique(combined$seurat_clusters[which(as.vector(combined@assays$RNA@counts[gene,]) != 0)])
  # num_cells <- rbind(num_cells, t(c("RAN", length(which(as.vector(combined@assays$RNA@counts[gene,]) != 0)))))
}
for (gene in rs) {
  num_clust_gene <- 0
  for (i in 0:num_clusters) {
    this_cells <- WhichCells(combined, idents = i)
    num_cells <- length(which(as.vector(combined@assays$RNA@counts[gene,this_cells]) != 0))
    if( num_cells/length(this_cells) > threshold ) {
      num_clust_gene <- num_clust_gene + 1
    }
  }
  num_clust <- rbind(num_clust, t(c(gene, "ROCK_SAND", num_clust_gene)))
  # num_cells <- rbind(num_cells, t(c("ROCK_SAND", length(which(as.vector(combined@assays$RNA@counts[gene,]) != 0)))))
}
colnames(num_clust) <- c("gene", "type", "num_clusters")
num_clust$num_clusters <- as.numeric(as.vector(num_clust$num_clusters))
ggplot(num_clust, aes(x=num_clusters, fill = type, color=type)) + geom_bar(alpha=0.5, position="identity") + ggtitle(paste("Number of Clusters Expressing RAN vs Rock-Sand at ", threshold*100, "%", sep="")) + xlab("Number of Clusters") + ylab("Number of Genes")

