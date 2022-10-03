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

#*******************************************************************************
# Demultiplexing ===============================================================
#*******************************************************************************
p_list = list()
cells.to.remove = data.frame()
ind_df = data.frame()
hb$demux = ""
hb$demux_top_prob = 0
hb$demux_second_top_prob = 0
for (this.str in c("b1", "b2", "b3", "c1", "c2", "c3")) {
        this.demux = data.table::fread(paste0("~/research/brain/results/demux/", this.str, "_demux.sing2"), data.table = F)
        this.demux.wide = reshape(this.demux[, c("BARCODE", "SM_ID", "POSTPRB")], idvar = "BARCODE", timevar = "SM_ID", direction = "wide")
        this.demux.wide$top_value = unlist(lapply(1:nrow(this.demux.wide), function(x) max(this.demux.wide[x,2:6]) ))
        this.demux.wide$second_top_value = unlist(lapply(1:nrow(this.demux.wide), function(x) sort(this.demux.wide[x,2:6], decreasing = T)[2] ))
        this.demux.wide$top_name = unlist(lapply(1:nrow(this.demux.wide), function(x) colnames(this.demux.wide)[2:6][which.max(this.demux.wide[x,2:6])] ))
        this.demux.wide$top_name = reshape2::colsplit(this.demux.wide$top_name, "\\.", c('1', '2'))[,2]
        p = ggplot(this.demux.wide, aes(x = top_value, y = second_top_value)) + geom_point(alpha = 0.1) + theme_classic() + ggtitle(this.str) + theme(plot.title = element_text(hjust = 0.5))
        p_list[[this.str]] = ggExtra::ggMarginal(p, type = "histogram")
        this.demux.wide$new_barcode = paste0(this.str, "_", this.demux.wide$BARCODE)
        hb$demux[this.demux.wide$new_barcode] = this.demux.wide$top_name
        hb$demux_top_prob[this.demux.wide$new_barcode] = this.demux.wide$top_value
        hb$demux_second_top_prob[this.demux.wide$new_barcode] = this.demux.wide$second_top_value
        
        # length(which(this.demux.wide$top_value > 0.9 & this.demux.wide$second_top_value < 0.1))
        # length(which(this.demux.wide$top_value < 0.9 | this.demux.wide$second_top_value > 0.1))
        # this.demux.wide.good = this.demux.wide[which(this.demux.wide$top_value > 0.9 & this.demux.wide$second_top_value < 0.1 & this.demux.wide$new_barcode %in% colnames(hb_test)),]
        # cells.to.remove = rbind(cells.to.remove, data.frame(cell = this.demux.wide$new_barcode[which( (this.demux.wide$top_value <= 0.9 | this.demux.wide$second_top_value >= 0.1) & this.demux.wide$new_barcode %in% colnames(hb_test))], sample = this.str))
        # ind_df = rbind(ind_df, data.frame(cell = this.demux.wide.good$new_barcode, ind = this.demux.wide.good$top_name))
}
cowplot::plot_grid(plotlist = p_list, ncol = 3)
hb$demux_good = hb$demux_top_prob > 0.9 & hb$demux_second_top_prob < 0.1
saveRDS(hb, "~/research/brain/data/hb_demux_092922.rds")

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
objs[["b1"]] = CreateSeuratObject(b1_counts)
objs[["b2"]] = CreateSeuratObject(b2_counts)
objs[["b3"]] = CreateSeuratObject(b3_counts)
objs[["c1"]] = CreateSeuratObject(c1_counts)
objs[["c2"]] = CreateSeuratObject(c2_counts)
objs[["c3"]] = CreateSeuratObject(c3_counts)

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

hb = subset(hb, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & pct_mt < 0.05)
print(paste("Number of Cells in hb After Filterning:", ncol(hb)))
hb = NormalizeData(hb, normalization.method = "LogNormalize", scale.factor = 10000)
hb = FindVariableFeatures(hb, selection.method = "vst", nfeatures = 2000)
hb = ScaleData(hb, features = rownames(hb))
hb = RunPCA(hb, features = VariableFeatures(hb))
hb = RunUMAP(hb, reduction = "pca", dims = 1:50)
hb = FindNeighbors(hb, reduction="umap", dims = 1:2)
hb = FindClusters(hb, resolution = .30)


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
