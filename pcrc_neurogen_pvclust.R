# Load Data
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("SeuratObject")
library("parallel")
bb = readRDS(paste0(rna_path, "data/bb_demux_102021.rds"))
Idents(bb) = bb$seurat_clusters

nperm = 10000
print(paste0("Loading Data: ", format(Sys.time(), "%X ")))
neurogen = read.csv("~/scratch/brain/data/conserved_neurogenesis_positive88_zfish_mouse_cichlid.csv")[,3]
pcrc = read.csv("~/scratch/brain/fst/pc_20_rc_20_10kb_bins_25kb_genes_on_lg_11_peak_by_bin.csv")[,2]
pcrc_neurogen_rowSums = rowSums(bb@assays$RNA@counts[c(pcrc, neurogen),] > 0)
data_mat_c = t(bb@assays$RNA@data[names(pcrc_neurogen_rowSums)[which(pcrc_neurogen_rowSums > 0)],])

cl <- makeCluster(24, type = "PSOCK")
fit = pvclust::parPvclust(cl=cl, data = as.matrix(data_mat_c), method.hclust = "average", method.dist = "correlation", nboot = 1000)
saveRDS(fit, "~/scratch/brain/results/pvclust_fit.rds")