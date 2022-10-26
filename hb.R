#*******************************************************************************
# Load Libraries ===============================================================
#*******************************************************************************
wdstr = substr(getwd(), 1, 12)
switch(wdstr,
       "C:/Users/mil" = { main_path = "C:/Users/miles/Downloads/";        },
       "/home/george" = { main_path = "~/research/"                       },
       "/storage/scr" = { main_path = "/storage/scratch1/6/ggruenhagen3/" },
       "/storage/hom" = { main_path = "/storage/scratch1/6/ggruenhagen3/" },
       "/storage/cod" = { main_path = "/storage/scratch1/6/ggruenhagen3/" })
brain_dir = paste0(main_path, "brain/")
data_dir  = paste0(main_path, "st/data/")
out_dir   = paste0(main_path, "st/results/")
if (main_path == "/storage/scratch1/6/ggruenhagen3/") { data_dir = "/storage/coda1/p-js585/0/ggruenhagen3/George/rich_project_pb1/data/st/data/" }
source(paste0(brain_dir, "/brain_scripts/all_f.R"))
setwd(out_dir)

# hb = readRDS("~/research/brain/data/hb_mi")

#*******************************************************************************
# Clustering Robustness ========================================================
#*******************************************************************************
source("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/chooser/my/cluster.stability/R/pipeline.R")
batch_num = 1
results_path = paste0("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/chooser/my/cluster.stability/examples/batch_data_", batch_num, "/")
min_dists = c(0.001, 0.1, 0.2, 0.3, 0.4, 0.5)
# min_dists = c(0.001, 0.1)
n_neighbors = c(5, 10, 20, 30, 40, 50)
resolutions = seq(0.2, 2.0, by = 0.2)
scores = list()
for (min_dist in min_dists) {
  for (n_neighbor in n_neighbors) {
    for (res in resolutions) {
      this.str = paste0("X", batch_num, "_", min_dist, "_", n_neighbor, "_res.", res)
      this.sil =  readRDS(paste0(results_path, "silhouette_grouped_", this.str, ".rds"))
      this.sil$min_dist = min_dist
      this.sil$n_neighbor = n_neighbor
      scores[[length(scores)+1]] = this.sil
    }
  }
}

# scores <- purrr::map(
#   paste0(results_path, "silhouette_grouped_", paste0("X", batch_num, "_", min_dists, "_", n_neighbors, "_res.", resolutions), ".rds"),
#   readRDS
# )

library("dplyr")
scores <- dplyr::bind_rows(scores) %>%
  dplyr::group_by(min_dist, n_neighbor, res) %>%
  dplyr::mutate("n_clusters" = dplyr::n()) %>%
  dplyr::ungroup()

unique(data.frame(scores)[,c("min_dist", "n_neighbor", "res", "n_clusters")])
range(unique(data.frame(scores)[,c("min_dist", "n_neighbor", "res", "n_clusters")])$n_clusters)

meds <- scores %>%
  dplyr::group_by(min_dist, n_neighbor, res) %>%
  dplyr::summarise(
    "boot" = list(boot_median(avg_sil)),
    "n_clusters" = mean(n_clusters)
  ) %>%
  tidyr::unnest_wider(boot)

threshold <- max(meds$low_med)
meds %>% dplyr::filter(med >= threshold) %>% dplyr::arrange(n_clusters) %>% tail(n = 1)
meds$id = paste(meds$min_dist, meds$n_neighbor, meds$res, sep = "_")
scores$id = paste(scores$min_dist, scores$n_neighbor, scores$res, sep = "_")
meds = meds %>% dplyr::arrange(n_clusters)
choice <- as.character(meds %>% dplyr::filter(med >= threshold) %>% dplyr::arrange(n_clusters) %>% tail(n = 1) %>% dplyr::pull(id))

p.meds = data.frame(meds)
p.meds$id = factor(p.meds$id, levels = unique(p.meds$id))
ggplot(p.meds, aes(id, med)) + geom_crossbar(aes(ymin = low_med, ymax = high_med), fill = "grey", size = 0.25) +
                               geom_hline(aes(yintercept = threshold), colour = "blue") +
                               geom_vline(aes(xintercept = choice), colour = "red") +
                               scale_x_discrete("Min Distance _ N Neighbors _ Resolution") +
                               scale_y_continuous("Silhouette Score", expand = c(0, 0), limits = c(-1, 1), breaks = seq(-1, 1, 0.25), oob = scales::squish) +
                               cowplot::theme_minimal_hgrid() +
                               theme(axis.title = element_text(size = 8), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7, angle = 45, vjust = 1, hjust = 1), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.ticks = element_line(colour = "black"),)
ggsave(filename = paste0("~/scratch/brain/results/hb/test_sdp.png"), dpi = 300, height = 3.5, width = 7, units = "in")
# ggplot(data.frame(meds), aes(id, med)) + geom_crossbar(aes(ymin = low_med, ymax = high_med), fill = "grey", size = 0.25) +
#   geom_hline(aes(yintercept = threshold), colour = "blue") +
#   geom_vline(aes(xintercept = choice), colour = "red") +
#   geom_jitter(data = scores, aes(id, avg_sil), size = 0.35, width = 0.15) +
#   scale_x_discrete("Min Distance _ N Neighbors _ Resolution") +
#   scale_y_continuous("Silhouette Score", expand = c(0, 0), limits = c(-1, 1), breaks = seq(-1, 1, 0.25), oob = scales::squish) +
#   cowplot::theme_minimal_hgrid() +
#   theme(axis.title = element_text(size = 8), axis.text = element_text(size = 7, angle = 45, vjust = 1, hjust = 1), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.ticks = element_line(colour = "black"),)

# Real Run
my.cluster = function(res) {
        this.obj = RunUMAP(this.obj, dims=1:50, min.dist=min_dist, spread=1, n.neighbors=n_neighbors, n.epochs=1000, metric="euclidean")
        this.obj = FindNeighbors(this.obj, reduction="umap", k.param=n_neighbors, dims=1:2, n.trees=500, prune.SNN=0)
        this.obj = FindClusters(this.obj, resolution = res, algorithm = 2)
        this.labels = this.obj$seurat_clusters
        return(as.numeric(as.vector(this.labels)))
}

this.obj = readRDS("~/scratch/brain/data/hb_minimized_101122.rds")
this.obj = SCTransform(this.obj, vars.to.regress = "sample", verbose = F)
this.obj = RunPCA(this.obj, dim = 50, verbose = F)
all.clustering.df = data.frame(cell = colnames(this.obj), row.names = colnames(this.obj))
for (min_dist in c(0.001, seq(0.1, 0.5, by = 0.1))) {
        for (n_neighbors in c(5, seq(10, 50, by = 5))) {
                print(paste0("Params: mindist", min_dist, "_nneighbors", n_neighbors))
                my.cluster.out = parallel::mclapply(seq(0.2, 2.0, by = 0.2), function(x) my.cluster(x), mc.cores = 9)
                this.clustering.df = do.call(cbind, my.cluster.out)
                colnames(this.clustering.df) = paste0("mindist", min_dist, "_nneighbors", n_neighbors, "_res", seq(0.2, 2.0, by = 0.2))
                all.clustering.df = cbind(all.clustering.df, this.clustering.df)
        }
}
all.clustering.df[,1] = NULL
write.csv(all.clustering.df, "~/scratch/brain/results/hb/hb_clustering_options101222.csv")

all.clustering.df = read.csv("~/scratch/brain/results/hb/hb_clustering_options101222.csv")
all.clustering.df$X = NULL
# all.rand = data.frame()
# for (i in 1:(ncol(all.clustering.df)-1)) {
#   print(paste0("i = ", i))
#   for (j in (i+1):ncol(all.clustering.df)) {
#     all.rand = rbind(all.rand, data.frame(i = i, j = j, name.i = colnames(all.clustering.df)[i], name.j = colnames(all.clustering.df)[j], rand = fossil::rand.index(all.clustering.df[,i], all.clustering.df[,j])))
#   }
# }
rand_pair = function(this.str) {
  i = as.numeric(strsplit(this.str, "_")[[1]][1])
  j = as.numeric(strsplit(this.str, "_")[[1]][2])
  this.rand = fossil::rand.index(all.clustering.df[,i], all.clustering.df[,j])
  return(this.rand)
}

# strs.i.j = c()
# for (i in 1:(ncol(all.clustering.df)-1)) { for (j in (i+1):ncol(all.clustering.df)) { strs.i.j = c(strs.i.j, paste0(i, "_", j)) } }
strs.i.j = expand.grid(1:600, 1:600)
strs.i.j <- data.table::setDT(data.table::CJ(strs.i.j[,1], strs.i.j[,2], unique = TRUE))
strs.i.j = as.data.frame(strs.i.j[!duplicated(t(apply(strs.i.j, 1, sort))),])
strs.i.j = strs.i.j[which(strs.i.j[,1] != strs.i.j[,2]),]
strs.i.j = paste0(as.vector(strs.i.j[,1]), "_", as.vector(strs.i.j[,2]))
# rand.list = unlist(parallel::mclapply(strs.i.j, function(x) rand_pair(x), mc.cores = 24))

all.rand = data.frame()
for (i in 1:10) {
  this.rand = read.csv(paste0("~/scratch/brain/results/hb/cluster_boot/rand/rand_", i, ".csv"))
  this.rand$X = NULL
  this.rand$str = strs.i.j[this.rand$combo.idx]
  all.rand = rbind(all.rand, this.rand)
}
all.rand[,c('i','j')] = reshape2::colsplit(all.rand$str, "_", c('1', '2'))
rownames(all.rand) = all.rand$str

all.n.clust = unlist(lapply(1:ncol(all.clustering.df), function(x) max(all.clustering.df[,x])+1 ))
all.rand$i.n.clust = all.n.clust[all.rand$i]
all.rand$j.n.clust = all.n.clust[all.rand$j]
# pdf("~/scratch/brain/results/hb/test.pdf", width = 7, height = 7)
all.rand = all.rand[order(all.rand$rand, decreasing = F),]
Cairo::Cairo(file = "~/scratch/brain/results/hb/test.png", width = 2400, height = 2400, dpi = 200)
print(ggplot(all.rand, aes(x = i.n.clust, y = j.n.clust, alpha = rand, color = rand)) + geom_point() + viridis::scale_color_viridis())
dev.off()

clust.max = unlist(lapply(1:ncol(all.clustering.df), function(x) max(all.clustering.df[,x])))
small.bin.start = seq(10, 240, by = 5)
small.bin.stop  = seq(15, 245, by = 5)-1
big.bin.start = seq(0,  225, by = 25)
big.bin.stop  = seq(25, 250, by = 25)-1
small.bin.df = data.frame(start = small.bin.start, stop = small.bin.stop)
small.bin.df$Freq = unlist(lapply(1:nrow(small.bin.df), function(x) length(which(clust.max >= small.bin.df$start[x] & clust.max <= small.bin.df$stop[x])) ))

small.bin.rand = data.frame()
for (i in 17:nrow(small.bin.df)) {
  print(i)
  bin.idx = which(clust.max >= small.bin.df$start[i] & clust.max <= small.bin.df$stop[i])
  if (length(bin.idx) > 1) {
    this.strs.i.j = expand.grid(bin.idx, bin.idx)
    this.strs.i.j = data.table::setDT(data.table::CJ(this.strs.i.j[,1], this.strs.i.j[,2], unique = TRUE))
    this.strs.i.j = as.data.frame(this.strs.i.j[!duplicated(t(apply(this.strs.i.j, 1, sort))),])
    this.strs.i.j = this.strs.i.j[which(this.strs.i.j[,1] != this.strs.i.j[,2]),]
    this.strs.i.j = paste0(as.vector(this.strs.i.j[,1]), "_", as.vector(this.strs.i.j[,2]))
    
    already.calc.strs = this.strs.i.j[which(this.strs.i.j %in% all.rand$str)]
    if (i < 17) {
      already.calc.df = data.frame(bin = i, str = already.calc.strs, rand = all.rand[already.calc.strs, "rand"])
      not.calc.strs = this.strs.i.j[which(!this.strs.i.j %in% all.rand$str)]
      not.calc.df = data.frame(bin = i, str = not.calc.strs, rand = NA)
      this.rand = unlist(parallel::mclapply(not.calc.strs, function(x) rand_pair(x), mc.cores = 3))
      not.calc.df$rand = this.rand
      this.df = rbind(already.calc.df, not.calc.df)
    } else {
      not.calc.strs = this.strs.i.j
      not.calc.df = data.frame(bin = i, str = not.calc.strs, rand = NA)
      this.rand = unlist(parallel::mclapply(not.calc.strs, function(x) rand_pair(x), mc.cores = 3))
      not.calc.df$rand = this.rand
      this.df = not.calc.df
    }
    
    small.bin.rand = rbind(small.bin.rand, this.df)
  }
}

small.bin.rand$start = small.bin.df$start[small.bin.rand$bin]
pdf("~/scratch/brain/results/hb/cluster_boot/rand/small_bin_rand.pdf", width = 10, height = 5)
ggplot(small.bin.rand, aes(x = start, y = rand, color = start, fill = start, group = start)) + geom_boxplot(alpha = 0.25) + geom_point(alpha = 0.25, position = position_jitterdodge())
dev.off()

#*******************************************************************************
# Demultiplexing ===============================================================
#*******************************************************************************
demux_summary = data.frame()
for (this.str in c("b1", "b2", "b3", "c1", "c2", "c3")) {
        this.demux = data.table::fread(paste0("~/research/brain/results/demux/qcc/", this.str, "_qc_cell_demux.best"), data.table = F)
        this.demux$sample = this.str
        this.demux$barcode = this.demux$BARCODE
        this.demux$cell = paste0(this.demux$sample, "_", this.demux$barcode)
        this.demux$demux = this.demux$SNG.1ST
        this.demux$demux_n_snp = this.demux$N.SNP
        this.demux$demux_llik = this.demux$SNG.LLK1
        demux_summary = rbind(demux_summary, this.demux[,c("sample", "barcode", "cell", "demux", "demux_n_snp", "demux_llik")])
}
write.csv(demux_summary, "~/research/brain/results/hb_demux_qcc_101222.csv")
hb$demux       = demux_summary$demux[match(colnames(hb), demux_summary$cell)]
hb$demux_n_snp = demux_summary$demux_n_snp[match(colnames(hb), demux_summary$cell)]
hb$demux_llik  = demux_summary$demux_llik[match(colnames(hb), demux_summary$cell)]
print(paste("Number of cells w/ <50 SNPs", length(which(hb$demux_n_snp < 50))))
print(paste("Number of cells w/ <50 SNPs and are good quality otherwise:", length(which(hb$demux_n_snp < 50 & hb$nFeature_RNA > 200 & hb$nFeature_RNA < 3000))))
print(paste("Number of cells after all filtering:", length(which(hb$demux_n_snp >= 50 & hb$nFeature_RNA > 200 & hb$nFeature_RNA < 3000))))

test = merge(this.demux.wide_1, this.demux.wide_2, by = "BARCODE", suffixes = c("_new", "_old"))
test$agree = test$top_name_new == test$top_name_old
test$new_better = (test$top_value_new > 0.9 & test$second_top_value_new < 0.1) & (test$top_value_old < 0.9 | test$second_top_value_old > 0.1)
test$old_better = (test$top_value_old > 0.9 & test$second_top_value_old < 0.1) & (test$top_value_new < 0.9 | test$second_top_value_new > 0.1)
test$top_value_dif = test$top_value_new - test$top_value_old
test$new_barcode = paste0(this.str, "_", test$BARCODE)
test$label = "Both Agree"
test$label[which(!test$agree)] = "Disagree: Equal Probability"
test$label[which(test$new_better)] = "Disagree: w/ Intergenic Better"
test$label[which(test$old_better)] = "Disagree: w/o Intergenic Better"
ggplot(test, aes(x = top_value_new, y = top_value_dif, color = label)) + geom_point(alpha = 0.5) + theme_bw() + xlab("Best Probability (w/ Intergenic)") + ylab("Best Probability (w/ Intergenic) - Best Probability (w/o Intergenic)")
ggplot(test, aes(x = top_value_old, y = top_value_new)) + geom_point(alpha = 0.5)

test2 = hb@meta.data[test$new_barcode,]
test2$label = test$label
p1 = ggplot(test2, aes(x = nCount_RNA, fill = label)) + geom_histogram(position="identity", alpha = 0.5)
p2 = ggplot(test2, aes(x = nCount_RNA, fill = label)) + geom_density(position="identity", alpha = 0.5)
p3 = ggplot(test2, aes(x = nFeature_RNA, fill = label)) + geom_histogram(position="identity", alpha = 0.5)
p4 = ggplot(test2, aes(x = nFeature_RNA, fill = label)) + geom_density(position="identity", alpha = 0.5)
cowplot::plot_grid(plotlist = list(p1, p3, p2, p4), ncol = 2)

test = hb@meta.data
test$bad_demux = F
test$bad_demux[which(test$demux_top_prob < 0.9 | test$demux_second_top_prob > 0.1)] = T
p1 = ggplot(test, aes(x = nCount_RNA, fill = bad_demux)) + geom_histogram(position="identity", alpha = 0.5)
p2 = ggplot(test, aes(x = nCount_RNA, fill = bad_demux)) + geom_density(position="identity", alpha = 0.5)
cowplot::plot_grid(plotlist = list(p1, p2), ncol = 1)
p3 = ggplot(test, aes(x = nFeature_RNA, fill = bad_demux)) + geom_histogram(position="identity", alpha = 0.5)
p4 = ggplot(test, aes(x = nFeature_RNA, fill = bad_demux)) + geom_density(position="identity", alpha = 0.5)
cowplot::plot_grid(plotlist = list(p3, p4), ncol = 1)
cowplot::plot_grid(plotlist = list(p1, p3, p2, p4), ncol = 2)
# B1 = 200, B2 = 446, B3 = 296, C1 = 298, C2 = 407, C3 = 

# Soup
demux = read.table("~/Downloads/c3_qc_cell_demux.best", sep="\t", header = T)
soup = read.table("~/Downloads/hb_c3_qcc_soup_dbl.tsv", sep="\t", header = T)
soup$best = unlist(lapply(1:nrow(soup), function(x) abs(sort(soup[x,c("cluster0", "cluster1", "cluster2", "cluster3")], decreasing = T)[1]) ))
soup$second_best = unlist(lapply(1:nrow(soup), function(x) abs(sort(soup[x,c("cluster0", "cluster1", "cluster2", "cluster3")], decreasing = T)[2]) ))
soup$prop = 1/(soup$best/soup$second_best)
soup$dif = soup$second_best - soup$best
# soup$col = paste0("cluster", soup$assignment)
soup$col = unlist(lapply(1:nrow(soup), function(x) colnames(soup)[6:10][which.max(soup[x,6:10])] ))
ggplot(soup, aes(x = best, y = prop)) + geom_point(alpha = 0.5)
ggplot(soup, aes(x = best, y = dif)) + geom_point(alpha = 0.5)
demux[, c("soup", "soup_best", "soup_second_best", "soup_prop")] = soup[, c("col", "best", "second_best", "prop")]
demux$prop = 1/(demux$SNG.LLK1 / demux$SNG.LLK2)
demux_big_ovlp = as.data.frame(table( paste0(demux$SNG.1ST, demux$soup) ))
demux_big_ovlp = demux_big_ovlp[which(demux_big_ovlp$Freq >= 5),]
ggplot(demux, aes(x = SNG.1ST, fill = soup)) + geom_bar()

ggplot(demux, aes(x = prop, y = soup_prop)) + geom_point(alpha = 0.5)
ggplot(demux[which(demux$prop > 1.5 & demux$soup_prop > 1.5),], aes(x = SNG.1ST, fill = soup)) + geom_bar()
ggplot(demux, aes(x = SNG.LLK1, y = SNG.LLK2, color = N.SNP, alpha = N.SNP)) + geom_point() + theme_classic() + ggtitle(this.str) + theme(plot.title = element_text(hjust = 0.5)) + scale_color_viridis()

# BB B1 testing
# Does the downsampled have 3138 cells instead of 2532 like it should?
soup.full = read.table("~/Downloads/soup_qc_cell_dbl.tsv", sep="\t", header = T)
soup.down = read.table("~/Downloads/soup_downsampled_dbl.tsv", sep="\t", header = T)
demux.down = read.table("~/Downloads/b1_down_demux.sing2", sep="\t", header = T)
this.demux.wide = reshape(this.demux[, c("BARCODE", "SM_ID", "POSTPRB")], idvar = "BARCODE", timevar = "SM_ID", direction = "wide")
this.demux.wide$top_value = unlist(lapply(1:nrow(this.demux.wide), function(x) max(this.demux.wide[x,2:5]) ))
this.demux.wide$second_top_value = unlist(lapply(1:nrow(this.demux.wide), function(x) sort(this.demux.wide[x,2:5], decreasing = T)[2] ))
this.demux.wide$top_name = unlist(lapply(1:nrow(this.demux.wide), function(x) colnames(this.demux.wide)[2:5][which.max(this.demux.wide[x,2:5])] ))
this.demux.wide$top_name = reshape2::colsplit(this.demux.wide$top_name, "\\.", c('1', '2'))[,2]
this.demux.wide$soup.full = soup.full$assignment
this.demux.wide$soup.down = soup.down$assignment
this.demux.wide = this.demux.wide[which(this.demux.wide$BARCODE %in% bb.meta$cell),]
this.demux.wide$bb = bb.meta$demux[match(bb.meta$cell, this.demux.wide$BARCODE)]
ggplot(this.demux.wide, aes(x = top_name, fill = soup.full)) + geom_bar()
ggplot(this.demux.wide, aes(x = top_name, fill = this.demux.wide$bb)) + geom_bar()
ggplot(this.demux.wide, aes(x = top_name, fill = soup.down)) + geom_bar()

test = rbind(hb_c3_demux2, bb_b1_demux)
p_list[[1]] = ggplot(test, aes(x = RD.TOTL, color = sample, fill = sample)) + geom_histogram(position = "identity", alpha = 0.5) + ggtitle("BB B1 Downsampled Comparison to HB C3 w/ Intergenic")
p_list[[2]] = ggplot(test, aes(x = RD.PASS, color = sample, fill = sample)) + geom_histogram(position = "identity", alpha = 0.5)
p_list[[3]] = ggplot(test, aes(x = N.SNP, color = sample, fill = sample)) + geom_histogram(position = "identity", alpha = 0.5)
cowplot::plot_grid(plotlist = p_list, ncol = 1)

bb_b1_demux = bb_b1_demux[which(bb_b1_demux$BARCODE %in% bb.meta$cell),]
bb_b1_demux$bb = bb.meta$demux[match(bb.meta$cell, bb_b1_demux$BARCODE)]
bb_b1_demux$bb_agree = bb_b1_demux$SNG.1ST == bb_b1_demux$bb
bb_b1_demux = bb_b1_demux[order(bb_b1_demux$bb_agree, decreasing = T),]
ggplot(bb_b1_demux, aes(x = SNG.LLK1, y = SNG.LLK2, color = N.SNP, alpha = N.SNP)) + geom_point() + theme_classic() + ggtitle("BB B1 Downsampled") + theme(plot.title = element_text(hjust = 0.5)) + scale_color_viridis()
ggplot(bb_b1_demux, aes(x = SNG.LLK1, y = N.SNP, color = bb_agree)) + geom_point(alpha = 0.5) + theme_classic() + ggtitle("BB B1 Downsampled") + theme(plot.title = element_text(hjust = 0.5)) + scale_color_viridis_d()
ggplot(bb_b1_demux, aes(x = SNG.LLK1, color = bb_agree, fill = bb_agree)) + geom_histogram(alpha = 0.5, position = "identity") + theme_classic() + ggtitle("BB B1 Downsampled") + theme(plot.title = element_text(hjust = 0.5))
ggplot(bb_b1_demux, aes(x = N.SNP, color = bb_agree, fill = bb_agree)) + geom_histogram(alpha = 0.5, position = "identity") + theme_classic() + ggtitle("BB B1 Downsampled") + theme(plot.title = element_text(hjust = 0.5))

#*******************************************************************************
# Initial clustering ===========================================================
#*******************************************************************************
# Load Data
dir_of_sr_dirs = "~/scratch/brain/bs/PM19/" # Folder where all the individual samples are kept
counts_list = list()
counts_list[["b1"]] = Read10X(paste0(dir_of_sr_dirs, "/PM19_B1_nuc/outs/filtered_feature_bc_matrix/"))
counts_list[["b2"]] = Read10X(paste0(dir_of_sr_dirs, "/PM19_B2_nuc/outs/filtered_feature_bc_matrix/"))
counts_list[["b3"]] = Read10X(paste0(dir_of_sr_dirs, "/PM19_B3_nuc/outs/filtered_feature_bc_matrix/"))
counts_list[["c1"]] = Read10X(paste0(dir_of_sr_dirs, "/PM19_C1_nuc/outs/filtered_feature_bc_matrix/"))
counts_list[["c2"]] = Read10X(paste0(dir_of_sr_dirs, "/PM19_C2_nuc/outs/filtered_feature_bc_matrix/"))
counts_list[["c3"]] = Read10X(paste0(dir_of_sr_dirs, "/PM19_C3_nuc/outs/filtered_feature_bc_matrix/"))

objs = list()
objs[["b1"]] = CreateSeuratObject(counts_list[["b1"]])
objs[["b2"]] = CreateSeuratObject(counts_list[["b2"]])
objs[["b3"]] = CreateSeuratObject(counts_list[["b3"]])
objs[["c1"]] = CreateSeuratObject(counts_list[["c1"]])
objs[["c2"]] = CreateSeuratObject(counts_list[["c2"]])
objs[["c3"]] = CreateSeuratObject(counts_list[["c3"]])

objs[["b1"]]$sample = "b1"; objs[["b1"]]$cond = "BHVE";
objs[["b2"]]$sample = "b2"; objs[["b2"]]$cond = "BHVE";
objs[["b3"]]$sample = "b3"; objs[["b3"]]$cond = "BHVE";
objs[["c1"]]$sample = "c1"; objs[["c1"]]$cond = "CTRL";
objs[["c2"]]$sample = "c2"; objs[["c2"]]$cond = "CTRL";
objs[["c3"]]$sample = "c3"; objs[["c3"]]$cond = "CTRL";

for (s in names(objs)) { objs[[s]] = RenameCells(objs[[s]], paste0(s)) }

# Check the quality
# for (i in names(objs)) {
#   obj = objs[[i]]
#   plot1 = ggplot(data.frame(nCount = obj$nCount_RNA), aes(x = nCount)) + geom_histogram() + theme_classic() + ggtitle("Histogram of nCount")
#   plot2 = VlnPlot(obj, features = "nCount_RNA") + NoLegend() + xlab("") + theme(axis.text.x = element_blank())
#   
#   plot3 = ggplot(data.frame(nFeature = obj$nFeature_RNA), aes(x = nFeature)) + geom_histogram() + theme_classic() + ggtitle("Histogram of nFeature")
#   plot4 = VlnPlot(obj, features = "nFeature_RNA") + NoLegend() + xlab("") + theme(axis.text.x = element_blank())
#   pdf(paste0("~/scratch/brain/results/hb_", i, "quality.pdf"), width = 8, height = 8)
#   print(cowplot::plot_grid(plotlist = list(plot1, plot2, plot3, plot4), ncol = 2))
#   dev.off()
#   print(paste0("# Cells UMI > 300 = ", length(which(obj$nCount_RNA > 300))))
#   print(paste0("# Cells Genes > 300 = ", length(which(obj$nFeature_RNA > 300))))
# }

hb = merge(objs[["b1"]], list(objs[["b2"]],objs[["b3"]],objs[["c1"]],objs[["c2"]],objs[["c3"]]))

gtf = read.delim("~/scratch/m_zebra_ref/gtf_processed.gtf", header = F)
mito.genes = gtf$V10[which(gtf$V1 == "NC_027944.1")]
mito.genes = stringr::str_replace(mito.genes, "_", "-")
hb$pct.mt = colSums(hb@assays$RNA@counts[mito.genes,]) / hb$nCount_RNA

hb = subset(hb, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & pct.mt < 0.05)
print(paste("Number of Cells in hb After Filterning:", ncol(hb)))
hb = NormalizeData(hb, normalization.method = "LogNormalize", scale.factor = 10000)
hb = SCTransform(hb, vars.to.regress = "sample" , verbose = TRUE)
hb@active.assay = "SCT"
hb = RunPCA(hb)
hb = RunUMAP(hb, reduction = "pca", dims = 1:50, min.dist = 0.5, spread = 0.2, n.neighbors = 50, n.epochs = 1000, metric = "euclidean")
hb = FindNeighbors(hb, reduction="umap", dims = 1:2, k.param = 50, n.trees = 500, prune.SNN = 0)
hb = FindClusters(hb, resolution = .30, algorithm = 2)


pdf(paste0("~/scratch/brain/results/hb_cluster.pdf"), width = 8, height = 8)
print(DimPlot(hb, reduction = "umap", label = TRUE, pt.size = 0.7))
dev.off()

pdf(paste0("~/scratch/brain/results/hb_cluster_sample.pdf"), width = 10, height = 6)
print(DimPlot(hb, split.by = "sample", reduction = "umap", label = TRUE, pt.size = 0.7, ncol = 3))
dev.off()


deg = FindAllMarkers(hb, only.pos = T)
write.csv(deg, "~/scratch/brain/results/hb_markers_loose_092122.csv")
write.csv(deg[which(deg$p_val_adj < 0.05),], "~/scratch/brain/results/hb_markers_sig_092122.csv")

# SC3
hb = readRDS("~/research/brain/data/hb_092122.rds")
source('~/research/brain/brain_scripts/all_f.R')
library(SingleCellExperiment)
library(SC3)
# library(scater)
# hb.sce = as.SingleCellExperiment(subset(hb, features = rownames(hb)[which(rowSums(hb@assays$RNA@counts) > 0)]))
# hb.sce = runPCA(hb.sce)
hb.sce = SingleCellExperiment(assays = list(counts = as.matrix(hb@assays$RNA@counts), logcounts = as.matrix(hb@assays$RNA@data)))
rowData(hb.sce)$feature_symbol <- rownames(hb.sce)
rm(hb)
hb.sce = sc3(hb.sce, ks = seq(10, 50, by = 10), gene_filter = F, biology = TRUE, n_cores = 1)
