
rna_path = "C:/Users/miles/Downloads/brain/"
marker_path <- paste(rna_path, "data/markers/", sep="")
marker_files <- dir(marker_path, pattern =paste("*.txt", sep=""))
markers <- data.frame(gene <- c(), bio <- c())
for (i in 1:length(marker_files)) {
  file <- read.table(paste(marker_path, marker_files[i], sep=""), header=FALSE, sep="\t", stringsAsFactors=FALSE)
  markers <- rbind(markers, file[,1:2])
}
colnames(markers) <- c("gene", "bio")
bio <- "ROCK_SAND"
markers <- markers[which(markers$bio == bio),]
valid_genes = markers$gene

# Overall
df = data.frame()
num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
for (cluster in 0:num_clusters) {
  this_cells = WhichCells(combined, idents = cluster)
  other_cells = colnames(combined)[which(! colnames(combined) %in% this_cells)]
  cluster_trans = colSums(combined@assays$RNA@data[valid_genes, this_cells ])
  other_trans   = colSums(combined@assays$RNA@data[valid_genes, other_cells])
  # this_p = wilcox.test(cluster_trans, other_trans, alternative = "greater")$p.value
  combined$this = colnames(combined) %in% this_cells
  # this_p = wilcox.test(colSums(combined@assays$RNA@data[valid_genes,]) ~ combined$this, alternative = "less")$p.value
  this_p = ks.test(cluster_trans, other_trans, alternative="less")$p.value
  # tmp_df = data.frame(exp = c(cluster_trans,other_trans), group=c(rep("cluster", length(cluster_trans)), rep("other", length(other_trans))))
  # res.aov = aov(exp ~ group , data=tmp_df)
  # sum.aov = summary(res.aov)
  # this_p = sum.aov[[1]]$`Pr(>F)`[1]
  newRow = t(c(cluster, mean(cluster_trans), mean(other_trans), this_p))
  df = rbind(df, newRow)
}
colnames(df) = c("cluster", "cluster_mean", "other_mean", "p")
df$q = p.adjust(df$p, method="bonferroni")
length(which(df$q < 0.05 & df$cluster_mean > df$other_mean))
ggplot(df, aes(cluster))

# See what gene drives the effects
multi_df = data.frame()
for (gene in valid_genes) {
  df = data.frame()
  for (cluster in 0:num_clusters) {
    this_cells = WhichCells(combined, idents = cluster)
    other_cells = colnames(combined)[which(! colnames(combined) %in% this_cells)]
    cluster_trans = combined@assays$RNA@data[gene, this_cells ]
    other_trans   = combined@assays$RNA@data[gene, other_cells]
    this_p = wilcox.test(cluster_trans, other_trans, alternative = "greater")$p.value
    newRow = t(c(cluster, mean(cluster_trans), mean(other_trans), this_p))
    df = rbind(df, newRow)
  }
  colnames(df) = c("cluster", "cluster_mean", "other_mean", "p")
  df$q = p.adjust(df$p, method="bonferroni")
  multi_df = rbind(multi_df, df[which(df$q < 0.05),])
}
multi_df$cluster = factor(multi_df$cluster, levels=0:num_clusters)
ggplot(multi_df, aes(cluster)) + geom_bar()

# Helper Functions
# findTransMean = function(obj) {
#   
# }
# Permutation Testing - Using difference of means
# For each cluster, do a permutation test. Shuffle all the clusters and store
# the difference of means between each cluster and the rest of the data.
# Then find 40 p-values.
obj = lncRNA
num_clusters <- as.numeric(tail(levels(obj@meta.data$seurat_clusters), n=1))
backup_ids <- obj@meta.data$seurat_clusters
num_perm = 200
perm_df = data.frame()
orig_test_stat = c()
for (i in 0:num_perm) {
  if (i==0) {
    for (cluster in 0:num_clusters) {
      # if(cluster == 0) { print(head(obj$seurat_clusters)) }
      this_cells  = WhichCells(obj, idents = cluster)
      other_cells = colnames(obj)[which(! colnames(obj) %in% this_cells)]
      this_trans_mean  = mean(colSums(obj@assays$RNA@data[valid_genes, this_cells ]))
      other_trans_mean = mean(colSums(obj@assays$RNA@data[valid_genes, other_cells]))
      orig_test_stat = c(orig_test_stat, this_trans_mean-other_trans_mean)
    }
  } else {
    set.seed(i)
    shuffled <- sample(backup_ids)
    obj$seurat_clusters = shuffled
    Idents(obj) = obj$seurat_clusters
    for (cluster in 0:num_clusters) {
      # if(cluster == 0) { print(head(obj$seurat_clusters))}
      this_cells  = WhichCells(obj, idents = cluster)
      other_cells = colnames(obj)[which(! colnames(obj) %in% this_cells)]
      this_trans_mean  = mean(colSums(obj@assays$RNA@data[valid_genes, this_cells ]))
      other_trans_mean = mean(colSums(obj@assays$RNA@data[valid_genes, other_cells]))
      perm_df = rbind(perm_df, t(c(i, cluster, this_trans_mean, other_trans_mean, this_trans_mean-other_trans_mean)))
    }
  }
}
colnames(perm_df) = c("perm", "cluster", "this_mean", "other_mean", "test_stat")
tmp_df = perm_df[which(perm_df$cluster == 0),]
tmp_df$above = tmp_df$test_stat > orig_test_stat[0+1]
ggplot(tmp_df, aes(test_stat, alpha=.7, fill=above, color = "maroon")) + geom_histogram() + geom_density() + geom_vline(aes(xintercept = quantile(test_stat, 0.95))) + geom_text(aes(x=quantile(test_stat, 0.95), label="95th Percentile"), y = Inf, vjust=1, color = "black")
# ggplot(tmp_df, aes(test_stat, alpha=.7, fill=above, color = "maroon")) + geom_histogram() + geom_density() + geom_vline(aes(xintercept = quantile(test_stat, 0.95))) + annotate("text", x=5, y = Inf, vjust=1, label ="text")

fig_path = "C:/Users/miles/Downloads/brain/results/ncRNA/perm/"
perm_df_p = data.frame()
for (cluster in 0:num_clusters) {
  # dir.create(file.path(fig_path, paste0("perm_", cluster)))
  tmp_df = perm_df[which(perm_df$cluster == cluster),]
  tmp_df$above = tmp_df$test_stat > orig_test_stat[cluster+1]
  # tmp_df$above[which(tmp_df$above)] = "maroon"
  # tmp_df$above[which(tmp_df$above == F)] = "blue"
  tmp_df$above = factor(tmp_df$above, levels = c(T,F))
  perm_df_p = rbind(perm_df_p, t(c(cluster, length(tmp_df$above[which(tmp_df$above == F)]) / num_perm)))
  png(paste0(fig_path, cluster, ".png"), width = 800, height=400)
  # p = ggplot(tmp_df, aes(test_stat, fill = above, color = above, alpha=0.7)) + geom_histogram() + geom_vline(aes(xintercept = quantile(test_stat, 0.95))) + scale_alpha(guide = 'none') + geom_text(aes(x=quantile(test_stat, 0.95), label="95th Percentile"), y = Inf, vjust=1, color = "black") + labs(fill = "Is Above Original Cluster Value?")
  p = ggplot(tmp_df, aes(test_stat, fill = above, color = above, alpha=0.7)) + geom_histogram() + 
    geom_vline(aes(xintercept = quantile(test_stat, 0.95))) + scale_alpha(guide = 'none') + geom_text(aes(x=quantile(test_stat, 0.95), label="95th Percentile"), y = Inf, vjust=1, color = "black") + 
    geom_vline(aes(xintercept = orig_test_stat[cluster+1])) + scale_alpha(guide = 'none') + geom_text(aes(x=orig_test_stat[cluster+1], label="True Mean"), y = Inf, vjust=1, color = "black") +
    scale_fill_discrete(limits = c('TRUE', 'FALSE')) + scale_color_discrete(limits = c('TRUE', 'FALSE'))
  print(p)
  dev.off()
}
colnames(perm_df_p) = c("cluster", "p")
perm_df_p$q = p.adjust(perm_df_p$p)
nrow(perm_df_p[which(perm_df_p$q < 0.05),])

# 
even_df = data.frame()
for (cluster in 0:num_clusters) {
  this_cells  = WhichCells(lncRNA, idents = cluster)
  this_trans_mean  = colSums(lncRNA@assays$RNA@data[, this_cells ])
  newRow = data.frame(rep(cluster, length(this_cells)), this_trans_mean)
  even_df = rbind(even_df, newRow)
}
colnames(even_df) = c("cluster", "mean_trans")
# ggplot(even_df, aes(cluster, mean_trans)) + stat_summary(fun = "mean", geom = "bar", position=position_dodge())
even_df$cluster = factor(even_df$cluster, levels=0:num_clusters)
ggplot(even_df, aes(x=cluster, y=mean_trans, group=cluster, fill=cluster, color=cluster)) + geom_boxplot(alpha=0.7) + geom_jitter(position=position_jitter(), alpha=0.1) + NoLegend() + labs(title="Data")

# Seurat Functions
WilcoxDETest <- function(
  data.use,
  cells.1,
  cells.2,
  verbose = TRUE,
  ...
) {
  data.use <- data.use[, c(cells.1, cells.2), drop = FALSE]
  j <- seq_len(length.out = length(x = cells.1))
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  limma.check <- PackageCheck("limma", error = FALSE)
  if (limma.check[1]) {
    p_val <- my.sapply(
      X = 1:nrow(x = data.use),
      FUN = function(x) {
        return(min(2 * min(limma::rankSumTestWithCorrelation(index = j, statistics = data.use[x, ])), 1))
      }
    )
  } else {
    if (getOption('Seurat.limma.wilcox.msg', TRUE)) {
      message(
        "For a more efficient implementation of the Wilcoxon Rank Sum Test,", 
        "\n(default method for FindMarkers) please install the limma package",
        "\n--------------------------------------------",
        "\ninstall.packages('BiocManager')",
        "\nBiocManager::install('limma')",
        "\n--------------------------------------------",
        "\nAfter installation of limma, Seurat will automatically use the more ",
        "\nefficient implementation (no further action necessary).",
        "\nThis message will be shown once per session"
      )
      options(Seurat.limma.wilcox.msg = FALSE)
    }
    group.info <- data.frame(row.names = c(cells.1, cells.2))
    group.info[cells.1, "group"] <- "Group1"
    group.info[cells.2, "group"] <- "Group2"
    group.info[, "group"] <- factor(x = group.info[, "group"])
    data.use <- data.use[, rownames(x = group.info), drop = FALSE]
    p_val <- sapply( 1:nrow(x = data.use), function(x) {
        return(wilcox.test(data.use[x, ] ~ group.info[, "group"], ...)$p.value)
      }
    )
  }
  return(data.frame(p_val, row.names = rownames(x = data.use)))
}