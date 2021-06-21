# Load BB
rna_path = "~/scratch/brain/"
source(paste0(rna_path, "brain_scripts/all_f.R"))
library("SeuratObject")
bb = readRDS(paste0(rna_path, "data/bb_cc_04072021.RDS"))
Idents(bb) = bb$seurat_clusters


library(pacman)
p_unload(SeuratDisk)
p_unload(Seurat)
p_load(Seurat)


# Load Correlation Matrices
library("rhdf5")
h5f = H5Fopen("~/scratch/brain/results/test/real_abs_cnr1_cor_bhve.h5")
b_cor = h5f$name
h5closeAll()
h5f = H5Fopen("~/scratch/brain/results/test/real_abs_cnr1_cor_ctrl.h5")
c_cor = h5f$name
h5closeAll()

# Find Genes with at least 5 cells in BHVE and CTRL
mat_b = bb@assays$RNA@counts[,which(bb@assays$RNA@counts["cnr1",] != 0 & bb$cond == "BHVE")]
mat_b[which(mat_b > 1)] = 1
num_cells = data.frame(B = rowSums(mat_b))
mat_c = bb@assays$RNA@counts[,which(bb@assays$RNA@counts["cnr1",] != 0 & bb$cond == "CTRL")]
mat_c[which(mat_c > 1)] = 1
num_cells$C = rowSums(mat_c)
num_cells$gene = rownames(num_cells)
num_cells$num_cells = num_cells$B + num_cells$C
above5 = num_cells$gene[which(num_cells$B >= 5 & num_cells$C >= 5)]

colnames(b_cor) = num_cells$gene
rownames(b_cor) = num_cells$gene
colnames(c_cor) = num_cells$gene
rownames(c_cor) = num_cells$gene

# Find Node Strength
b_ns = rowSums(b_cor, na.rm = T)
c_ns = rowSums(c_cor, na.rm = T)

ns_df = data.frame(B = b_ns, C = c_ns, Dif = b_ns - c_ns, Dif_Abs = abs(b_ns - c_ns), num_cells = num_cells$num_cells)
ns_df$gene = num_cells$gene
sum(ns_df$Dif_Abs[which(ns_df$gene %in% above5)])
ns_df$B_Up = T
ns_df$B_Up[which(ns_df$C > ns_df$B)] = F
ns_df$B_Up = factor(ns_df$B_Up, levels = c(T, F))

# Plots
png("~/scratch/brain/results/cnr1_genes_ns.png", width = 600, height = 600)
print(ggplot(ns_df[which(ns_df$gene %in% above5),], aes(C, B, color = B_Up)) + geom_point() + ylab("BHVE NS") + xlab("CTRL NS"))
dev.off()
png("~/scratch/brain/results/cnr1_genes_ns_cells.png", width = 600, height = 600)
print(ggplot(ns_df[which(ns_df$gene %in% above5),], aes(num_cells, Dif, color = B_Up)) + geom_point() + ylab("NS Dif") + xlab("Number of Cells"))
dev.off()
ns_df2 = reshape2::melt(ns_df[above5, c("B", "C", "gene")])
png("~/scratch/brain/results/cnr1_genes_ns_hist.png", width = 1000, height = 600)
print(ggplot(ns_df2, aes(x = value, fill = variable, color = variable)) + geom_histogram(position="identity", alpha=0.5, bins = 100) + xlab("NS"))
dev.off()
b_cor2 = b_cor[above5, above5]
b_cor2[which(upper.tri(b_cor2, diag = T))] = NA
test = b_cor2[which( !is.na(b_cor2) )]
c_cor2 = c_cor[above5, above5]
c_cor2[which(upper.tri(c_cor2, diag = T))] = NA
test2 = c_cor2[which( !is.na(c_cor2) )]
test3 = data.frame(value = test, variable = "BHVE")
test3 = rbind(test3, data.frame(value = test2, variable = "CTRL"))
png("~/scratch/brain/results/cnr1_hist.png", width = 1000, height = 600)
print(ggplot(test3, aes(x = value, fill = variable, color = variable)) + geom_histogram(position="identity", alpha=0.5, bins = 100) + xlab("Correlation") + xlim(0, 0.2))
dev.off()

# Find p-values for correlations
my_cor_t = function(r, n) (r * sqrt(n - 2))/sqrt(1 - r**2)
my_cor_p = function(t, n) 2*pt(-abs(t), df=n-2)
b_t_mat = my_cor_t(b_cor,   length(which(bb@assays$RNA@counts["cnr1",] != 0 & bb$cond == "BHVE")))
b_p_mat = my_cor_p(b_t_mat, length(which(bb@assays$RNA@counts["cnr1",] != 0 & bb$cond == "BHVE")))
c_t_mat = my_cor_t(c_cor,   length(which(bb@assays$RNA@counts["cnr1",] != 0 & bb$cond == "CTRL")))
c_p_mat = my_cor_p(c_t_mat, length(which(bb@assays$RNA@counts["cnr1",] != 0 & bb$cond == "CTRL")))

b_p_mat = readRDS("~/scratch/brain/data/cnr1_b_p.RDS")
c_p_mat = readRDS("~/scratch/brain/data/cnr1_c_p.RDS")

# Clear Memory
b_t_mat = NULL
c_t_mat = NULL

# Find Significant Edges
b_bon_mat = p.adjust(b_p_mat, method = "bonferroni")
b_bon_mat = matrix(b_bon_mat, nrow = nrow(b_p_mat), byrow = T)
c_bon_mat = p.adjust(c_p_mat, method = "bonferroni")
c_bon_mat = matrix(c_bon_mat, nrow = nrow(c_p_mat), byrow = T)
rownames(b_bon_mat) = colnames(b_bon_mat) = rownames(c_bon_mat) = colnames(c_bon_mat) = colnames(b_p_mat)

min.num.cells = 20
min.bon = 1e-35

# Behave Edges
b_bon_mat[which(upper.tri(b_bon_mat, diag = T))] = NA
b_bon_df = data.frame(which(b_bon_mat < min.bon, arr.ind=TRUE))
b_bon_df$Source = rownames(b_bon_mat)[b_bon_df$row]
b_bon_df$Target = rownames(b_bon_mat)[b_bon_df$col]
b_bon_df$Value = b_cor[which(b_bon_mat < min.bon)]
b_bon_df$Source_B = num_cells$B[match(b_bon_df$Source, num_cells$gene)]
b_bon_df$Source_C = num_cells$C[match(b_bon_df$Source, num_cells$gene)]
b_bon_df$Target_B = num_cells$B[match(b_bon_df$Target, num_cells$gene)]
b_bon_df$Target_C = num_cells$C[match(b_bon_df$Target, num_cells$gene)]
b_bon_df = b_bon_df[which(b_bon_df$Source_B >= min.num.cells & b_bon_df$Source_C >= min.num.cells & b_bon_df$Target_B >= min.num.cells & b_bon_df$Target_C >= min.num.cells),c("Source", "Target", "Value")]
nrow(b_bon_df)
write.csv(b_bon_df, "~/scratch/brain/results/cnr1_b_sig_big.csv", row.names = F)
system(paste0("rclone copy ~/scratch/brain/results/cnr1_b_sig_big.csv dropbox:BioSci-Streelman/George/Brain/bb/results/py_ns/receptors/cnr1_network"))

# Control Edges
c_bon_mat[which(upper.tri(c_bon_mat, diag = T))] = NA
c_bon_df = data.frame(which(c_bon_mat < min.bon, arr.ind=TRUE))
c_bon_df$Source = rownames(c_bon_mat)[c_bon_df$row]
c_bon_df$Target = rownames(c_bon_mat)[c_bon_df$col]
c_bon_df$Value = c_cor[which(c_bon_mat < min.bon)]
c_bon_df$Source_B = num_cells$B[match(c_bon_df$Source, num_cells$gene)]
c_bon_df$Source_C = num_cells$C[match(c_bon_df$Source, num_cells$gene)]
c_bon_df$Target_B = num_cells$B[match(c_bon_df$Target, num_cells$gene)]
c_bon_df$Target_C = num_cells$C[match(c_bon_df$Target, num_cells$gene)]
c_bon_df = c_bon_df[which(c_bon_df$Source_B >= min.num.cells & c_bon_df$Source_C >= min.num.cells & c_bon_df$Target_B >= min.num.cells & c_bon_df$Target_C >= min.num.cells),c("Source", "Target", "Value")]
nrow(c_bon_df)
write.csv(c_bon_df, "~/scratch/brain/results/cnr1_c_sig_big.csv", row.names = F)
system(paste0("rclone copy ~/scratch/brain/results/cnr1_c_sig_big.csv dropbox:BioSci-Streelman/George/Brain/bb/results/py_ns/receptors/cnr1_network"))


# Node Table
ns_df$Id = ns_df$gene
ns_df$Label = ns_df$gene
ns_df = ns_df[,c("Id", "Label", colnames(ns_df)[1:(ncol(ns_df)-2)])]
write.csv(ns_df, "~/scratch/brain/results/cnr1_node_df.csv", row.names = F)
system(paste0("rclone copy ~/scratch/brain/results/cnr1_node_df.csv dropbox:BioSci-Streelman/George/Brain/bb/results/py_ns/receptors/cnr1_network"))

# Zack plot
min.bon = 1e-10
bvc_cor = b_cor - c_cor
bvc_cor = bvc_cor[which(num_cells$B >= 5 & num_cells$C >= 5), which(num_cells$B >= 5 & num_cells$C >= 5)]
bvc_p = r_to_p(b_cor, c_cor, length(which(bb@assays$RNA@counts["cnr1",] > 0 & bb$cond == "BHVE")), length(which(bb@assays$RNA@counts["cnr1",] > 0 & bb$cond == "CTRL")))
bvc_bon_mat = p.adjust(bvc_p2, method = "bonferroni")
bvc_bon_mat = matrix(bvc_bon_mat, nrow = nrow(bvc_p2), byrow = T)
rownames(bvc_bon_mat) = colnames(bvc_bon_mat) = above5
print(paste("Number of Gene-Gene Correlations with Bon < 1e-10: ", length(which(bvc_bon_mat < 1e-10))))

b_cor2 = b_cor[which(num_cells$B >= 5 & num_cells$C >= 5), which(num_cells$B >= 5 & num_cells$C >= 5)]
bvc_bon_mat[which(upper.tri(bvc_bon_mat, diag = T))] = NA
bvc_bon_df = data.frame(which(bvc_bon_mat < min.bon, arr.ind=TRUE))
bvc_bon_df$Source = rownames(bvc_bon_mat)[bvc_bon_df$row]
bvc_bon_df$Target = rownames(bvc_bon_mat)[bvc_bon_df$col]
bvc_bon_df$Value = diag(b_cor2[bvc_bon_df$row, bvc_bon_df$col])
bvc_bon_df$B_Up = diag(b_cor2[bvc_bon_df$row, bvc_bon_df$col]) > diag(c_cor2[bvc_bon_df$row, bvc_bon_df$col])
write.csv(bvc_bon_df, "~/scratch/brain/results/cnr1_b_z2.csv", row.names = F)
system(paste0("rclone copy ~/scratch/brain/results/cnr1_b_z2.csv dropbox:BioSci-Streelman/George/Brain/bb/results/py_ns/receptors/cnr1_network"))

c_cor2 = c_cor[which(num_cells$B >= 5 & num_cells$C >= 5), which(num_cells$B >= 5 & num_cells$C >= 5)]
bvc_bon_mat[which(upper.tri(bvc_bon_mat, diag = T))] = NA
bvc_bon_df = data.frame(which(bvc_bon_mat < min.bon, arr.ind=TRUE))
bvc_bon_df$Source = rownames(bvc_bon_mat)[bvc_bon_df$row]
bvc_bon_df$Target = rownames(bvc_bon_mat)[bvc_bon_df$col]
bvc_bon_df$Value = diag(c_cor2[bvc_bon_df$row, bvc_bon_df$col])
bvc_bon_df$B_Up = diag(b_cor2[bvc_bon_df$row, bvc_bon_df$col]) > diag(c_cor2[bvc_bon_df$row, bvc_bon_df$col])
write.csv(bvc_bon_df, "~/scratch/brain/results/cnr1_c_z2.csv", row.names = F)
system(paste0("rclone copy ~/scratch/brain/results/cnr1_c_z2.csv dropbox:BioSci-Streelman/George/Brain/bb/results/py_ns/receptors/cnr1_network"))
