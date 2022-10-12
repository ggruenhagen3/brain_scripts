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
# chooseR ======================================================================
#*******************************************************************************
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
