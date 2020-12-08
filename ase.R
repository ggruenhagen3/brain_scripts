library("RUnit")
library("MBASED")
library("metap")
library("DESeq2")
library("apeglm")
library("EnhancedVolcano")

#=========================================================================================
# New UMD2a Data
#=========================================================================================
# DEG
rna_path = "C:/Users/miles/Downloads/brain/"
SRR904 = read.table(paste(rna_path, "/data/pit_castle_deg/SRR5440904_counter_per_gene.tsv", sep=""), sep="\t", header = F, stringsAsFactors = F)
SRR905 = read.table(paste(rna_path, "/data/pit_castle_deg/SRR5440905_counter_per_gene.tsv", sep=""), sep="\t", header = F, stringsAsFactors = F)
SRR906 = read.table(paste(rna_path, "/data/pit_castle_deg/SRR5440906_counter_per_gene.tsv", sep=""), sep="\t", header = F, stringsAsFactors = F)
SRR907 = read.table(paste(rna_path, "/data/pit_castle_deg/SRR5440907_counter_per_gene.tsv", sep=""), sep="\t", header = F, stringsAsFactors = F)
SRR908 = read.table(paste(rna_path, "/data/pit_castle_deg/SRR5440908_counter_per_gene.tsv", sep=""), sep="\t", header = F, stringsAsFactors = F)
SRR909 = read.table(paste(rna_path, "/data/pit_castle_deg/SRR5440909_counter_per_gene.tsv", sep=""), sep="\t", header = F, stringsAsFactors = F)

SRR904 = SRR904[which(! duplicated(SRR904[,1])),]
SRR905 = SRR905[which(! duplicated(SRR905[,1])),]
SRR906 = SRR906[which(! duplicated(SRR906[,1])),]
SRR907 = SRR907[which(! duplicated(SRR907[,1])),]
SRR908 = SRR908[which(! duplicated(SRR908[,1])),]
SRR909 = SRR909[which(! duplicated(SRR909[,1])),]

genes = SRR904[-c(1:5),1]
mat = as.matrix(cbind(SRR904[-c(1:5),2], SRR905[-c(1:5),2], SRR906[-c(1:5),2], SRR907[-c(1:5),2], SRR908[-c(1:5),2], SRR909[-c(1:5),2]), 
                dimnames=list(genes, c("4", "5", "6", "7", "8", "9")))
mycolData = data.frame(samples=c("4", "5", "6", "7", "8", "9"), 
                       cond=c("pit", "pit", "castle", "castle", "iso", "iso"),
                       isBhve=c("bhve", "bhve", "bhve", "bhve", "ctrl", "ctrl"))

dds = DESeqDataSetFromMatrix(countData = mat,
                             colData = mycolData,
                             design = ~  cond)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="cond_pit_vs_castle")
res <- lfcShrink(dds, coef="cond_pit_vs_castle", type="apeglm")

sig_ind = which(res$padj < 0.05 & res$log2FoldChange > 1)
sig_genes = genes[sig_ind]
res_df = data.frame(gene=genes, logFC=res$log2FoldChange, padj=res$padj)
rownames(res_df) = res_df$gene
EnhancedVolcano(res_df, lab=rownames(res_df), x="logFC", y="padj") + labs(subtitle="Pit v Castle Volcano Plot") + theme(plot.title = element_blank(), plot.caption = element_blank())
sig_genes_hgnc = hgncMzebraInPlace(data.frame(sig_genes), 1, gene_names)
write.table(sig_genes, "C:/Users/miles/Downloads/brain/results/pit_v_castle_deg.txt", quote=F, col.names = F, row.names = F)
write.table(sig_genes_hgnc, "C:/Users/miles/Downloads/brain/results/pit_v_castle_deg_hgnc.txt", quote=F, col.names = F, row.names = F)

# BHVE v CTRL
dds = DESeqDataSetFromMatrix(countData = mat,
                             colData = mycolData,
                             design = ~ isBhve)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="isBhve_ctrl_vs_bhve")
res <- lfcShrink(dds, coef="isBhve_ctrl_vs_bhve", type="apeglm")

sig_ind = which(res$padj < 0.05 & res$log2FoldChange > 1)
sig_genes = genes[sig_ind]
res_df = data.frame(gene=genes, logFC=res$log2FoldChange, padj=res$padj)
rownames(res_df) = res_df$gene
EnhancedVolcano(res_df, lab=rownames(res_df), x="logFC", y="padj") + labs(subtitle="Bhve v Ctrl Volcano Plot") + theme(plot.title = element_blank(), plot.caption = element_blank())
sig_genes_hgnc = hgncMzebraInPlace(data.frame(sig_genes), 1, gene_names)
write.table(sig_genes, "C:/Users/miles/Downloads/brain/results/bhve_v_ctrl_deg.txt", quote=F, col.names = F, row.names = F)
write.table(sig_genes_hgnc, "C:/Users/miles/Downloads/brain/results/bhve_v_ctrl_deg_hgnc.txt", quote=F, col.names = F, row.names = F)

# Sim Pit v Castle
mat_pvc = as.matrix(cbind(SRR904[-c(1:5),2], SRR905[-c(1:5),2], SRR906[-c(1:5),2], SRR907[-c(1:5),2]), 
                    dimnames=list(genes, c("4", "5", "6", "7")))
colData_pvc = data.frame(samples=c("4", "5", "6", "7"), 
                         sim1=c("pit", "castle", "pit", "castle"),
                         sim2=c("castle", "pit", "castle", "pit"))
dds = DESeqDataSetFromMatrix(countData = mat_pvc,
                             colData = colData_pvc,
                             design = ~sim1)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="sim1_pit_vs_castle")
res <- lfcShrink(dds, coef="sim1_pit_vs_castle", type="apeglm")

sig_ind = which(res$padj < 0.05 & res$log2FoldChange > 1)
sig_genes = genes[sig_ind]
res_df = data.frame(gene=genes, logFC=res$log2FoldChange, padj=res$padj)
rownames(res_df) = res_df$gene
EnhancedVolcano(res_df, lab=rownames(res_df), x="logFC", y="padj") + labs(subtitle="Simulated Pit v Castle 1 Volcano Plot") + theme(plot.title = element_blank(), plot.caption = element_blank())
sig_genes_hgnc = hgncMzebraInPlace(data.frame(sig_genes), 1, gene_names)
write.table(sig_genes, "C:/Users/miles/Downloads/brain/results/sim1_pit_v_castle.txt", quote=F, col.names = F, row.names = F)
write.table(sig_genes_hgnc, "C:/Users/miles/Downloads/brain/results/sim1_pit_v_castle_hgnc.txt", quote=F, col.names = F, row.names = F)

dds <- DESeq(dds)
res <- results(dds, name="sim2_pit_vs_castle")
res <- lfcShrink(dds, coef="sim2_pit_vs_castle", type="apeglm")

sig_ind = which(res$padj < 0.05 & res$log2FoldChange > 1)
sig_genes = genes[sig_ind]
res_df = data.frame(gene=genes, logFC=res$log2FoldChange, padj=res$padj)
rownames(res_df) = res_df$gene
EnhancedVolcano(res_df, lab=rownames(res_df), x="logFC", y="padj") + labs(subtitle="Simulated Pit v Castle 2 Volcano Plot") + theme(plot.title = element_blank(), plot.caption = element_blank())
sig_genes_hgnc = hgncMzebraInPlace(data.frame(sig_genes), 1, gene_names)
write.table(sig_genes, "C:/Users/miles/Downloads/brain/results/sim2_pit_v_castle.txt", quote=F, col.names = F, row.names = F)
write.table(sig_genes_hgnc, "C:/Users/miles/Downloads/brain/results/sim2_pit_v_castle_hgnc.txt", quote=F, col.names = F, row.names = F)

# Pit v Iso
mat_pvi = as.matrix(cbind(SRR904[-c(1:5),2], SRR905[-c(1:5),2], SRR908[-c(1:5),2], SRR909[-c(1:5),2]), 
                    dimnames=list(genes, c("4", "5", "8", "9")))
colData_pvi = data.frame(samples=c("4", "5", "8", "9"), 
                         cond=c("pit", "pit", "iso", "iso"))
dds = DESeqDataSetFromMatrix(countData = mat_pvi,
                             colData = colData_pvi,
                             design = ~cond)
dds <- DESeq(dds)
res <- results(dds, name="cond_pit_vs_iso")
res <- lfcShrink(dds, coef="cond_pit_vs_iso", type="apeglm")
sig_ind = which(res$padj < 0.05 & res$log2FoldChange > 1)
sig_genes = genes[sig_ind]
res_df = data.frame(gene=genes, logFC=res$log2FoldChange, padj=res$padj)
rownames(res_df) = res_df$gene
EnhancedVolcano(res_df, lab=rownames(res_df), x="logFC", y="padj") + labs(subtitle="Pit v Isolated Volcano Plot") + theme(plot.title = element_blank(), plot.caption = element_blank())
sig_genes_hgnc = hgncMzebraInPlace(data.frame(sig_genes), 1, gene_names)
write.table(sig_genes, "C:/Users/miles/Downloads/brain/results/pit_v_iso.txt", quote=F, col.names = F, row.names = F)
write.table(sig_genes_hgnc, "C:/Users/miles/Downloads/brain/results/pit_v_iso_hgnc.txt", quote=F, col.names = F, row.names = F)

# Castle v Iso
mat_cvi = as.matrix(cbind(SRR906[-c(1:5),2], SRR907[-c(1:5),2], SRR908[-c(1:5),2], SRR909[-c(1:5),2]), 
                    dimnames=list(genes, c("6", "7", "8", "9")))
colData_cvi = data.frame(samples=c("6", "7", "8", "9"), 
                         cond=c("castle", "castle", "iso", "iso"))
dds = DESeqDataSetFromMatrix(countData = mat_cvi,
                             colData = colData_cvi,
                             design = ~cond)
dds <- DESeq(dds)
res <- results(dds, name="cond_iso_vs_castle")
res <- lfcShrink(dds, coef="cond_iso_vs_castle", type="apeglm")
sig_ind = which(res$padj < 0.05 & res$log2FoldChange > 1)
sig_genes = genes[sig_ind]
res_df = data.frame(gene=genes, logFC=res$log2FoldChange, padj=res$padj)
rownames(res_df) = res_df$gene
EnhancedVolcano(res_df, lab=rownames(res_df), x="logFC", y="padj") + labs(subtitle="Castle v Isolated Volcano Plot") + theme(plot.title = element_blank(), plot.caption = element_blank())
sig_genes_hgnc = hgncMzebraInPlace(data.frame(sig_genes), 1, gene_names)
write.table(sig_genes, "C:/Users/miles/Downloads/brain/results/castle_v_iso.txt", quote=F, col.names = F, row.names = F)
write.table(sig_genes_hgnc, "C:/Users/miles/Downloads/brain/results/castle_v_iso_hgnc.txt", quote=F, col.names = F, row.names = F)

# Dendrogram
mat = matrix(cbind(SRR904[-c(1:5),2], SRR905[-c(1:5),2], SRR906[-c(1:5),2], SRR907[-c(1:5),2], SRR908[-c(1:5),2], SRR909[-c(1:5),2]), ncol = 6, dimnames = list(genes, c("pit", "pit", "castle", "castle", "iso", "iso")))
pit_v_castle_genes = read.table("C:/Users/miles/Downloads/brain/results/pit_v_castle.txt", header=F)
pit_v_castle_genes = as.vector(pit_v_castle_genes$V1)
p = degDend(mat, pit_v_castle_genes, "C:/Users/miles/Downloads/brain/results/pit_v_castle_dend.png", include_samples = c("pit", "castle"))
p = degDend(mat, pit_v_castle_genes, "C:/Users/miles/Downloads/brain/results/pit_v_castle_all_dend.png")

pit_v_iso_genes = read.table("C:/Users/miles/Downloads/brain/results/pit_v_iso.txt", header=F)
pit_v_iso_genes = as.vector(pit_v_iso_genes$V1)
p = degDend(mat, pit_v_iso_genes, "C:/Users/miles/Downloads/brain/results/pit_v_iso_dend.png", include_samples = c("pit", "iso"))

castle_v_iso_genes = read.table("C:/Users/miles/Downloads/brain/results/castle_v_iso.txt", header=F)
castle_v_iso_genes = as.vector(castle_v_iso_genes$V1)
p = degDend(mat, castle_v_iso_genes, "C:/Users/miles/Downloads/brain/results/castle_v_iso_dend.png", include_samples = c("castle", "iso"))


# SNP-level data
SRR904 = read.table(paste(rna_path, "/data/ase/SRR5440904_informative.vcf", sep=""), header = F, stringsAsFactors = F)
SRR905 = read.table(paste(rna_path, "/data/ase/SRR5440905_informative.vcf", sep=""), header = F, stringsAsFactors = F)
SRR906 = read.table(paste(rna_path, "/data/ase/SRR5440906_informative.vcf", sep=""), header = F, stringsAsFactors = F)
SRR907 = read.table(paste(rna_path, "/data/ase/SRR5440907_informative.vcf", sep=""), header = F, stringsAsFactors = F)
SRR908 = read.table(paste(rna_path, "/data/ase/SRR5440908_informative.vcf", sep=""), header = F, stringsAsFactors = F)
SRR909 = read.table(paste(rna_path, "/data/ase/SRR5440909_informative.vcf", sep=""), header = F, stringsAsFactors = F)

SRR904$V6 = as.numeric(as.vector(SRR904$V6))
SRR905$V6 = as.numeric(as.vector(SRR905$V6))
SRR906$V6 = as.numeric(as.vector(SRR906$V6))
SRR907$V6 = as.numeric(as.vector(SRR907$V6))
SRR908$V6 = as.numeric(as.vector(SRR908$V6))
SRR909$V6 = as.numeric(as.vector(SRR909$V6))

SRR904$V7 = as.numeric(as.vector(SRR904$V7))
SRR905$V7 = as.numeric(as.vector(SRR905$V7))
SRR906$V7 = as.numeric(as.vector(SRR906$V7))
SRR907$V7 = as.numeric(as.vector(SRR907$V7))
SRR908$V7 = as.numeric(as.vector(SRR908$V7))
SRR909$V7 = as.numeric(as.vector(SRR909$V7))

SRR904$MC_COUNTS = SRR904$V6
SRR905$MC_COUNTS = SRR905$V6
SRR906$MC_COUNTS = SRR906$V6
SRR907$MC_COUNTS = SRR907$V6
SRR908$MC_COUNTS = SRR908$V6
SRR909$MC_COUNTS = SRR909$V6

SRR904$MC_COUNTS[which(SRR904$V14 == "False")] = SRR904$V7[which(SRR904$V14 == "False")]
SRR905$MC_COUNTS[which(SRR905$V14 == "False")] = SRR905$V7[which(SRR905$V14 == "False")]
SRR906$MC_COUNTS[which(SRR906$V14 == "False")] = SRR906$V7[which(SRR906$V14 == "False")]
SRR907$MC_COUNTS[which(SRR907$V14 == "False")] = SRR907$V7[which(SRR907$V14 == "False")]
SRR908$MC_COUNTS[which(SRR908$V14 == "False")] = SRR908$V7[which(SRR908$V14 == "False")]
SRR909$MC_COUNTS[which(SRR909$V14 == "False")] = SRR909$V7[which(SRR909$V14 == "False")]

SRR904$CV_COUNTS = SRR904$V7
SRR905$CV_COUNTS = SRR905$V7
SRR906$CV_COUNTS = SRR906$V7
SRR907$CV_COUNTS = SRR907$V7
SRR908$CV_COUNTS = SRR908$V7
SRR909$CV_COUNTS = SRR909$V7

SRR904$CV_COUNTS[which(SRR904$V14 == "False")] = SRR904$V6[which(SRR904$V14 == "False")]
SRR905$CV_COUNTS[which(SRR905$V14 == "False")] = SRR905$V6[which(SRR905$V14 == "False")]
SRR906$CV_COUNTS[which(SRR906$V14 == "False")] = SRR906$V6[which(SRR906$V14 == "False")]
SRR907$CV_COUNTS[which(SRR907$V14 == "False")] = SRR907$V6[which(SRR907$V14 == "False")]
SRR908$CV_COUNTS[which(SRR908$V14 == "False")] = SRR908$V6[which(SRR908$V14 == "False")]
SRR909$CV_COUNTS[which(SRR909$V14 == "False")] = SRR909$V6[which(SRR909$V14 == "False")]

SRR904$pos = paste0(SRR904$V1, ":", SRR904$V2, "-", SRR904$V2)
SRR905$pos = paste0(SRR905$V1, ":", SRR905$V2, "-", SRR905$V2)
SRR906$pos = paste0(SRR906$V1, ":", SRR906$V2, "-", SRR906$V2)
SRR907$pos = paste0(SRR907$V1, ":", SRR907$V2, "-", SRR907$V2)
SRR908$pos = paste0(SRR908$V1, ":", SRR908$V2, "-", SRR908$V2)
SRR909$pos = paste0(SRR909$V1, ":", SRR909$V2, "-", SRR909$V2)

pit = inner_join(SRR904, SRR905, by = "pos")
pit_mc = pit$MC_COUNTS.x + pit$MC_COUNTS.y
pit_cv = pit$CV_COUNTS.x + pit$CV_COUNTS.y
names(pit_mc) = pit$pos
names(pit_cv) = pit$pos

castle = inner_join(SRR906, SRR907, by = "pos")
castle_mc = castle$MC_COUNTS.x + castle$MC_COUNTS.y
castle_cv = castle$CV_COUNTS.x + castle$CV_COUNTS.y
names(castle_mc) = castle$pos
names(castle_cv) = castle$pos

iso = inner_join(SRR908, SRR909, by = "pos")
iso_mc = iso$MC_COUNTS.x + iso$MC_COUNTS.y
iso_cv = iso$CV_COUNTS.x + iso$CV_COUNTS.y
names(iso_mc) = iso$pos
names(iso_cv) = iso$pos

pit_v_castle_res = my_MBASED(pit_mc, pit_cv, castle_mc, castle_cv, "pit", "castle", pit$pos, n_boot, isSNP=T)
castle_v_pit_res = my_MBASED(castle_mc, castle_cv, pit_mc, pit_cv, "castle", "pit", pit$pos, n_boot, isSNP=T)
pit_v_castle_pos = pit_v_castle_res[[2]]
castle_v_pit_pos = castle_v_pit_res[[2]]
ovlp_pc_v_cp_pos = pit_v_castle_pos[which(pit_v_castle_pos %in% castle_v_pit_pos)]

pit_v_iso_res = my_MBASED(pit_mc, pit_cv, iso_mc, iso_cv, "pit", "iso", pit$pos, n_boot, isSNP=T)
iso_v_pit_res = my_MBASED(iso_mc, iso_cv, pit_mc, pit_cv, "iso", "pit", pit$pos, n_boot, isSNP=T)
pit_v_iso_pos = pit_v_iso_res[[2]]
iso_v_pit_pos = iso_v_pit_res[[2]]
ovlp_pi_v_ip_pos = pit_v_iso_pos[which(pit_v_iso_pos %in% iso_v_pit_pos)]

castle_v_iso_res = my_MBASED(castle_mc, castle_cv, iso_mc, iso_cv, "castle", "iso", castle$pos, n_boot, isSNP=T)
iso_v_castle_res = my_MBASED(iso_mc, iso_cv, castle_mc, castle_cv, "iso", "csatle", castle$pos, n_boot, isSNP=T)
castle_v_iso_pos = castle_v_iso_res[[2]]
iso_v_castle_pos = iso_v_castle_res[[2]]
ovlp_ci_v_ic_pos = castle_v_iso_pos[which(castle_v_iso_pos %in% iso_v_castle_pos)]

gtf = read.table("C:/Users/miles/Downloads/brain/brain_scripts/full_ens_w_ncbi_gene.gtf", sep="\t", header=F, stringsAsFactors = F)
gtf = gtf[which(gtf[,3] == "gene" & gtf[,1] != "NC_027944.1"),]
gtf_gene_name <- c()
for (i in 1:nrow(gtf)) {
  start <- gregexpr(pattern ='gene_name', gtf$V9[i])[[1]]
  stop  <- gregexpr(pattern =';', substr(gtf$V9[i], start, nchar(gtf$V9[i])))[[1]][1]
  gene_name <- substr(gtf$V9[i], start+10, start+stop-2)
  if (start == -1) {
    gene_name <- substr(gtf$V9[i], start+10, start+stop)
  }
  gtf_gene_name <- c(gtf_gene_name, gene_name)
}
gtf$gene_name <- gtf_gene_name
colnames(gtf) <- c("LG", "source", "type", "start", "stop", "idk", "idk1", "idk2", "info", "gene_name")
gtf = gtf[which(! startsWith(gtf$gene_name, "LOC")),]

pit_v_castle_genes = posToGene(pit_v_castle_pos, gtf)

###################
# Gene Level Data #
###################
rna_path = "C:/Users/miles/Downloads/brain/"
SRR904 = read.table(paste(rna_path, "/data/ase/SRR5440904_RG_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)
SRR905 = read.table(paste(rna_path, "/data/ase/SRR5440905_RG_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)
SRR906 = read.table(paste(rna_path, "/data/ase/SRR5440906_RG_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)
SRR907 = read.table(paste(rna_path, "/data/ase/SRR5440907_RG_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)
SRR908 = read.table(paste(rna_path, "/data/ase/SRR5440908_RG_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)
SRR909 = read.table(paste(rna_path, "/data/ase/SRR5440909_RG_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)

SRR904$GENE = str_replace(SRR904$GENE,"%", " (1 of many)")
SRR905$GENE = str_replace(SRR905$GENE,"%", " (1 of many)")
SRR906$GENE = str_replace(SRR906$GENE,"%", " (1 of many)")
SRR907$GENE = str_replace(SRR907$GENE,"%", " (1 of many)")
SRR908$GENE = str_replace(SRR908$GENE,"%", " (1 of many)")
SRR909$GENE = str_replace(SRR909$GENE,"%", " (1 of many)")

#NCBI
# SRR904 = read.table(paste(rna_path, "/data/ase/SRR5440904_RG_nc_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)
# SRR905 = read.table(paste(rna_path, "/data/ase/SRR5440905_RG_nc_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)
# SRR906 = read.table(paste(rna_path, "/data/ase/SRR5440906_RG_nc_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)
# SRR907 = read.table(paste(rna_path, "/data/ase/SRR5440907_RG_nc_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)
# SRR908 = read.table(paste(rna_path, "/data/ase/SRR5440908_RG_nc_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)
# SRR909 = read.table(paste(rna_path, "/data/ase/SRR5440909_RG_nc_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)


gene_names = SRR904$GENE
# gene_names = str_replace(gene_names,"%", " (1 of many)")
n_boot = 100

# Prepare the Data
pit_mc = SRR904$MC_COUNTS + SRR905$MC_COUNTS
pit_cv = SRR904$CV_COUNTS + SRR905$CV_COUNTS
names(pit_mc) = gene_names
names(pit_cv) = gene_names

castle_mc = SRR906$MC_COUNTS + SRR907$MC_COUNTS
castle_cv = SRR906$CV_COUNTS + SRR907$CV_COUNTS
names(castle_mc) = gene_names
names(castle_cv) = gene_names

iso_mc = SRR908$MC_COUNTS + SRR909$MC_COUNTS
iso_cv = SRR908$CV_COUNTS + SRR909$CV_COUNTS
names(iso_mc) = gene_names
names(iso_cv) = gene_names

# Check for Skews in ASE ratios
hist(log2(SRR904$MC_COUNTS/SRR904$CV_COUNTS), breaks=50)
hist(log2(SRR905$MC_COUNTS/SRR905$CV_COUNTS), breaks=50)
hist(log2(SRR906$MC_COUNTS/SRR906$CV_COUNTS), breaks=50)
hist(log2(SRR907$MC_COUNTS/SRR907$CV_COUNTS), breaks=50)
hist(log2(SRR908$MC_COUNTS/SRR908$CV_COUNTS), breaks=50)
hist(log2(SRR909$MC_COUNTS/SRR909$CV_COUNTS), breaks=50)

pos = which(SRR904$MC_COUNTS > 0 & SRR904$CV_COUNTS > 0)
d = density(log2(SRR904$MC_COUNTS[pos]/SRR904$CV_COUNTS[pos]))
plot(d, main="Pit 1")
pos = which(SRR905$MC_COUNTS > 0 & SRR905$CV_COUNTS > 0)
d = density(log2(SRR905$MC_COUNTS[pos]/SRR905$CV_COUNTS[pos]))
plot(d, main="Pit 2")
pos = which(SRR906$MC_COUNTS > 0 & SRR906$CV_COUNTS > 0)
d = density(log2(SRR906$MC_COUNTS[pos]/SRR906$CV_COUNTS[pos]))
plot(d, main="Castle 1")
pos = which(SRR907$MC_COUNTS > 0 & SRR907$CV_COUNTS > 0)
d = density(log2(SRR907$MC_COUNTS[pos]/SRR907$CV_COUNTS[pos]))
plot(d, main="Castle 2")
pos = which(SRR908$MC_COUNTS > 0 & SRR908$CV_COUNTS > 0)
d = density(log2(SRR908$MC_COUNTS[pos]/SRR908$CV_COUNTS[pos]))
plot(d, main="Isolated 1")
pos = which(SRR909$MC_COUNTS > 0 & SRR909$CV_COUNTS > 0)
d = density(log2(SRR909$MC_COUNTS[pos]/SRR909$CV_COUNTS[pos]))
plot(d, main="Isolated 2")

# Find Discordant ASE
disc_ase = gene_names[which(sign(SRR904$dif) == sign(SRR905$dif) & sign(SRR904$dif) != 0 & 
                            sign(SRR906$dif) == sign(SRR907$dif) & sign(SRR906$dif) != 0 & 
                            sign(SRR906$dif) != sign(SRR904$dif))]
disc_ase_pc = disc_ase
disc_ase_pc_hgnc = sort(hgncMzebraInPlace(data.frame(disc_ase_pc), 1, rownames(tj)))
write.table(disc_ase_pc, "C:/Users/miles/Downloads/brain/results/ase/disc_ASE_pc.txt", quote = F, col.names = F, row.names = F)
write.table(disc_ase_pc_hgnc, "C:/Users/miles/Downloads/brain/results/ase/disc_ASE_pc_hgnc.txt", quote = F, col.names = F, row.names = F)

# Do 1-sampled ASE experiments
pos_all_ind = which(SRR904$MC_COUNTS + SRR904$CV_COUNTS > 0 &
                    SRR905$MC_COUNTS + SRR905$CV_COUNTS > 0 &
                    SRR906$MC_COUNTS + SRR906$CV_COUNTS > 0 &
                    SRR907$MC_COUNTS + SRR907$CV_COUNTS > 0 &
                    SRR908$MC_COUNTS + SRR908$CV_COUNTS > 0 &
                    SRR909$MC_COUNTS + SRR909$CV_COUNTS > 0)

SRR904_1 = my_MBASED_1(SRR904$MC_COUNTS, SRR904$CV_COUNTS, "SRR904 (Pit 1)", "", gene_names, n_boot)
SRR905_1 = my_MBASED_1(SRR905$MC_COUNTS, SRR905$CV_COUNTS, "SRR905 (Pit 2)", "", gene_names, n_boot)
SRR906_1 = my_MBASED_1(SRR906$MC_COUNTS, SRR906$CV_COUNTS, "SRR906 (Castle 1)", "", gene_names, n_boot)
SRR907_1 = my_MBASED_1(SRR907$MC_COUNTS, SRR907$CV_COUNTS, "SRR907 (Castle 2)", "", gene_names, n_boot)
SRR908_1 = my_MBASED_1(SRR908$MC_COUNTS, SRR908$CV_COUNTS, "SRR908 (Isolated 1)", "", gene_names, n_boot)
SRR909_1 = my_MBASED_1(SRR909$MC_COUNTS, SRR909$CV_COUNTS, "SRR909 (Isolated 2)", "", gene_names, n_boot)

SRR904$p = 1; SRR905$p = 1; SRR906$p = 1; SRR907$p = 1; SRR908$p = 1; SRR909$p = 1
SRR904$q = 1; SRR905$q = 1; SRR906$q = 1; SRR907$q = 1; SRR908$q = 1; SRR909$q = 1

SRR904$p[which(SRR904$GENE %in% rownames(assays(SRR904_1[[1]])$pValueASE))] = assays(SRR904_1[[1]])$pValueASE
SRR905$p[which(SRR905$GENE %in% rownames(assays(SRR905_1[[1]])$pValueASE))] = assays(SRR905_1[[1]])$pValueASE
SRR906$p[which(SRR906$GENE %in% rownames(assays(SRR906_1[[1]])$pValueASE))] = assays(SRR906_1[[1]])$pValueASE
SRR907$p[which(SRR907$GENE %in% rownames(assays(SRR907_1[[1]])$pValueASE))] = assays(SRR907_1[[1]])$pValueASE
SRR908$p[which(SRR908$GENE %in% rownames(assays(SRR908_1[[1]])$pValueASE))] = assays(SRR908_1[[1]])$pValueASE
SRR909$p[which(SRR909$GENE %in% rownames(assays(SRR909_1[[1]])$pValueASE))] = assays(SRR909_1[[1]])$pValueASE

SRR904$q[which(SRR904$GENE %in% rownames(assays(SRR904_1[[1]])$pValueASE))] = p.adjust(assays(SRR904_1[[1]])$pValueASE, method="bonferroni")
SRR905$q[which(SRR905$GENE %in% rownames(assays(SRR905_1[[1]])$pValueASE))] = p.adjust(assays(SRR905_1[[1]])$pValueASE, method="bonferroni")
SRR906$q[which(SRR906$GENE %in% rownames(assays(SRR906_1[[1]])$pValueASE))] = p.adjust(assays(SRR906_1[[1]])$pValueASE, method="bonferroni")
SRR907$q[which(SRR907$GENE %in% rownames(assays(SRR907_1[[1]])$pValueASE))] = p.adjust(assays(SRR907_1[[1]])$pValueASE, method="bonferroni")
SRR908$q[which(SRR908$GENE %in% rownames(assays(SRR908_1[[1]])$pValueASE))] = p.adjust(assays(SRR908_1[[1]])$pValueASE, method="bonferroni")
SRR909$q[which(SRR909$GENE %in% rownames(assays(SRR909_1[[1]])$pValueASE))] = p.adjust(assays(SRR909_1[[1]])$pValueASE, method="bonferroni")

SRR904$sig = SRR904$q < 0.05
SRR905$sig = SRR905$q < 0.05
SRR906$sig = SRR906$q < 0.05
SRR907$sig = SRR907$q < 0.05
SRR908$sig = SRR908$q < 0.05
SRR909$sig = SRR909$q < 0.05

SRR904$dif = SRR904$MC_COUNTS - SRR904$CV_COUNTS
SRR905$dif = SRR905$MC_COUNTS - SRR905$CV_COUNTS
SRR906$dif = SRR906$MC_COUNTS - SRR906$CV_COUNTS
SRR907$dif = SRR907$MC_COUNTS - SRR907$CV_COUNTS
SRR908$dif = SRR908$MC_COUNTS - SRR908$CV_COUNTS
SRR909$dif = SRR909$MC_COUNTS - SRR909$CV_COUNTS

SRR904$ase = log2(SRR904$CV_COUNTS / SRR904$MC_COUNTS)
SRR905$ase = log2(SRR905$CV_COUNTS / SRR905$MC_COUNTS)
SRR906$ase = log2(SRR906$CV_COUNTS / SRR906$MC_COUNTS)
SRR907$ase = log2(SRR907$CV_COUNTS / SRR907$MC_COUNTS)
SRR908$ase = log2(SRR908$CV_COUNTS / SRR908$MC_COUNTS)
SRR909$ase = log2(SRR909$CV_COUNTS / SRR909$MC_COUNTS)
df = data.frame(cbind(SRR904$GENE, SRR904$ase, SRR905$ase, SRR906$ase, SRR907$ase, SRR908$ase, SRR909$ase))
df = data.frame(cbind(SRR904, SRR905, SRR906, SRR907, SRR908, SRR909))
df$Digging_1 = as.numeric(as.vector(df$Digging_1))
df$Digging_2 = as.numeric(as.vector(df$Digging_2))
df$Building_1 = as.numeric(as.vector(df$Building_1))
df$Building_2 = as.numeric(as.vector(df$Building_2))
df$Control_1 = as.numeric(as.vector(df$Control_1))
df$Control_2 = as.numeric(as.vector(df$Control_2))
my_ryan = df[which(df[,1] %in% ryan$X),]
my_ryan = my_ryan[match(ryan$X, my_ryan[,1]),]
colnames(my_ryan) = c("GENE", "Digging_1", "Digging_2", "Building_1", "Building_2", "Control_1", "Control_2")
my_ryan$Digging_Mean_ASE = (as.numeric(as.vector(my_ryan$Digging_1)) + as.numeric(as.vector(my_ryan$Digging_2)))/2
my_ryan$Building_Mean_ASE = (as.numeric(as.vector(my_ryan$Building_1)) + as.numeric(as.vector(my_ryan$Building_2)))/2
my_ryan$Control_Mean_ASE = (as.numeric(as.vector(my_ryan$Control_1)) + as.numeric(as.vector(my_ryan$Control_2)))/2

length(which( is.na(my_ryan[,2]) & is.na(my_ryan[,3]) & is.na(my_ryan[,4]) & is.na(my_ryan[,5]) & is.na(my_ryan[,6]) & is.na(my_ryan[,7]) ))
length(which(sign(my_ryan$Digging_Mean_ASE) == sign(ryan$Digging_Mean_ASE) & sign(my_ryan$Building_Mean_ASE) == sign(ryan$Building_Mean_ASE) & sign(my_ryan$Control_Mean_ASE) == sign(ryan$Iso_Mean_ASE)  ))

all_sig_same_dir = SRR904$GENE[which(SRR904$sig & SRR905$sig & SRR906$sig & SRR907$sig & SRR908$sig & SRR909$sig &
                                     sign(SRR905$dif) == sign(SRR904$dif) & 
                                     sign(SRR906$dif) == sign(SRR904$dif) &
                                     sign(SRR907$dif) == sign(SRR904$dif) & 
                                     sign(SRR908$dif) == sign(SRR904$dif) & 
                                     sign(SRR909$dif) == sign(SRR904$dif) )]

all_same_dir = SRR904$GENE[which(sign(SRR905$dif) == sign(SRR904$dif) & SRR905$dif != 0 & SRR904$dif != 0 &
                                 sign(SRR906$dif) == sign(SRR904$dif) & SRR906$dif != 0 & SRR904$dif != 0 &
                                 sign(SRR907$dif) == sign(SRR904$dif) & SRR907$dif != 0 & SRR904$dif != 0 &
                                 sign(SRR908$dif) == sign(SRR904$dif) & SRR908$dif != 0 & SRR904$dif != 0 &
                                 sign(SRR909$dif) == sign(SRR904$dif) & SRR909$dif != 0 & SRR904$dif != 0 )]

# Zack's Method: Combine p-values
# all_p = sapply(1:length(gene_names), function(x) sumlog(c(SRR904$p[x], SRR905$p[x], SRR906$p[x], SRR907$p[x], SRR908$p[x], SRR909$p[x]))$p)
all_p = sapply(1:length(gene_names), function(x) sumz(c(SRR904$p[x], SRR905$p[x], SRR906$p[x], SRR907$p[x], SRR908$p[x], SRR909$p[x]))$p)
all_p[which(SRR904$p == 0 &
            SRR905$p == 0 &
            SRR906$p == 0 &
            SRR907$p == 0 &
            SRR908$p == 0 &
            SRR909$p == 0 )] = 0
all_q = p.adjust(all_p, method = "BH")
agg = gene_names[which(all_q < 0.05 & 
                 sign(SRR905$dif) == sign(SRR904$dif) & 
                 sign(SRR906$dif) == sign(SRR904$dif) &
                 sign(SRR907$dif) == sign(SRR904$dif) & 
                 sign(SRR908$dif) == sign(SRR904$dif) & 
                 sign(SRR909$dif) == sign(SRR904$dif) )]

write.table(all_sig_same_dir, "C:/Users/miles/Downloads/brain/results/ase_all_sig_same_dir_RG.txt", quote = F, col.names = F, row.names = F)
write.table(all_same_dir, "C:/Users/miles/Downloads/brain/results/ase_all_same_dir_RG.txt", quote = F, col.names = F, row.names = F)
write.table(agg, "C:/Users/miles/Downloads/brain/results/ase_agg_sig_same_dir_RG.txt", quote = F, col.names = F, row.names = F)

all_sig_same_dir_hgnc = hgncMzebraInPlace(data.frame(all_sig_same_dir), 1, gene_names)
all_same_dir_hgnc     = hgncMzebraInPlace(data.frame(all_same_dir),     1, gene_names)
agg_hgnc              = hgncMzebraInPlace(data.frame(agg),              1, gene_names)

write.table(all_sig_same_dir_hgnc, "C:/Users/miles/Downloads/brain/results/ase_all_sig_same_dir_hgnc_RG.txt", quote = F, col.names = F, row.names = F)
write.table(all_same_dir_hgnc,     "C:/Users/miles/Downloads/brain/results/ase_all_same_dir_hgnc_RG.txt", quote = F, col.names = F, row.names = F)
write.table(agg_hgnc,              "C:/Users/miles/Downloads/brain/results/ase_agg_sig_same_dir_hgnc_RG.txt", quote = F, col.names = F, row.names = F)


# Do 2-sampled ASE experiments

# Pit v Castle
pit_v_castle_res = my_MBASED(pit_mc, pit_cv, castle_mc, castle_cv, "pit", "castle", gene_names, n_boot)
pit_v_castle_genes = pit_v_castle_res[[2]]
castle_v_pit_res = my_MBASED(castle_mc, castle_cv, pit_mc, pit_cv, "castle", "pit", gene_names, n_boot)
castle_v_pit_genes = castle_v_pit_res[[2]]
ovlp_pc_v_cp = pit_v_castle_genes[which(pit_v_castle_genes %in% castle_v_pit_genes)]

# Pit v Isolated
pit_v_iso_res = my_MBASED(pit_mc, pit_cv, iso_mc, iso_cv, "pit", "iso", gene_names, n_boot)
pit_v_iso_genes = pit_v_iso_res[[2]]
iso_v_pit_res = my_MBASED(iso_mc, iso_cv, pit_mc, pit_cv, "iso", "pit", gene_names, n_boot)
iso_v_pit_genes = iso_v_pit_res[[2]]
ovlp_pi_v_ip = pit_v_iso_genes[which(pit_v_iso_genes %in% iso_v_pit_genes)]

# Castle v Isolated
castle_v_iso_res = my_MBASED(castle_mc, castle_cv, iso_mc, iso_cv, "castle", "iso", gene_names, n_boot)
castle_v_iso_genes = castle_v_iso_res[[2]]
iso_v_castle_res = my_MBASED(iso_mc, iso_cv, castle_mc, castle_cv, "iso", "castle", gene_names, n_boot)
iso_v_castle_genes = iso_v_castle_res[[2]]
ovlp_ci_v_ic = castle_v_iso_genes[which(castle_v_iso_genes %in% iso_v_castle_genes)]

res = data.frame(test=c("pit_v_castle", "castle_v_pit", "pvc_cvp_ovlp", "pit_v_iso", "iso_v_pit", "pvi_ivp_ovlp", "castle_v_iso", "iso_v_castle", "cvi_ivc"),
           num_genes=c(length(pit_v_castle_genes), length(castle_v_pit_genes), length(ovlp_pc_v_cp),
                       length(pit_v_iso_genes), length(iso_v_pit_genes), length(ovlp_pi_v_ip),
                       length(castle_v_iso_genes), length(iso_v_castle_genes), length(ovlp_ci_v_ic)))

write.table(pit_v_castle_genes, "C:/Users/miles/Downloads/brain/results/ase_pit_v_castle_RG.txt", quote = F, col.names = F, row.names = F)
write.table(castle_v_pit_genes, "C:/Users/miles/Downloads/brain/results/ase_castle_v_pit_RG.txt", quote = F, col.names = F, row.names = F)
write.table(pit_v_iso_genes, "C:/Users/miles/Downloads/brain/results/ase_pit_v_iso_RG.txt", quote = F, col.names = F, row.names = F)
write.table(iso_v_pit_genes, "C:/Users/miles/Downloads/brain/results/ase_iso_v_pit_RG.txt", quote = F, col.names = F, row.names = F)
write.table(castle_v_iso_genes, "C:/Users/miles/Downloads/brain/results/ase_castle_v_iso_RG.txt", quote = F, col.names = F, row.names = F)
write.table(iso_v_castle_genes, "C:/Users/miles/Downloads/brain/results/ase_iso_v_castle_RG.txt", quote = F, col.names = F, row.names = F)
write.table(ovlp_pc_v_cp, "C:/Users/miles/Downloads/brain/results/ase_ovlp_pc_v_cp_RG.txt", quote = F, col.names = F, row.names = F)
write.table(ovlp_pi_v_ip, "C:/Users/miles/Downloads/brain/results/ase_ovlp_pi_v_ip_RG.txt", quote = F, col.names = F, row.names = F)
write.table(ovlp_ci_v_ic, "C:/Users/miles/Downloads/brain/results/ase_ovlp_ci_v_ic_RG.txt", quote = F, col.names = F, row.names = F)

pit_v_castle_genes_hgnc = hgncMzebraInPlace(data.frame(pit_v_castle_genes), 1, gene_names)
castle_v_pit_genes_hgnc = hgncMzebraInPlace(data.frame(castle_v_pit_genes), 1, gene_names)
pit_v_iso_genes_hgnc    = hgncMzebraInPlace(data.frame(pit_v_iso_genes),    1, gene_names)
iso_v_pit_genes_hgnc    = hgncMzebraInPlace(data.frame(iso_v_pit_genes),    1, gene_names)
castle_v_iso_genes_hgnc = hgncMzebraInPlace(data.frame(pit_v_iso_genes),    1, gene_names)
iso_v_castle_genes_hgnc = hgncMzebraInPlace(data.frame(iso_v_castle_genes), 1, gene_names)
ovlp_pc_v_cp_hgnc       = hgncMzebraInPlace(data.frame(ovlp_pc_v_cp),       1, gene_names)
ovlp_pi_v_ip_hgnc       = hgncMzebraInPlace(data.frame(ovlp_pi_v_ip),       1, gene_names)
ovlp_ci_v_ic_hgnc       = hgncMzebraInPlace(data.frame(ovlp_ci_v_ic),       1, gene_names)

write.table(pit_v_castle_genes_hgnc, "C:/Users/miles/Downloads/brain/results/ase_pit_v_castle_hgnc.txt", quote = F, col.names = F, row.names = F)
write.table(castle_v_pit_genes_hgnc, "C:/Users/miles/Downloads/brain/results/ase_castle_v_pit_hgnc.txt", quote = F, col.names = F, row.names = F)
write.table(pit_v_iso_genes_hgnc, "C:/Users/miles/Downloads/brain/results/ase_pit_v_iso_hgnc.txt", quote = F, col.names = F, row.names = F)
write.table(iso_v_pit_genes_hgnc, "C:/Users/miles/Downloads/brain/results/ase_iso_v_pit_hgnc.txt", quote = F, col.names = F, row.names = F)
write.table(castle_v_iso_genes_hgnc, "C:/Users/miles/Downloads/brain/results/ase_castle_v_iso_hgnc.txt", quote = F, col.names = F, row.names = F)
write.table(iso_v_castle_genes_hgnc, "C:/Users/miles/Downloads/brain/results/ase_iso_v_castle_hgnc.txt", quote = F, col.names = F, row.names = F)
write.table(ovlp_pc_v_cp_hgnc, "C:/Users/miles/Downloads/brain/results/ase_ovlp_pc_v_cp_hgnc.txt", quote = F, col.names = F, row.names = F)
write.table(ovlp_pi_v_ip_hgnc, "C:/Users/miles/Downloads/brain/results/ase_ovlp_pi_v_ip_hgnc.txt", quote = F, col.names = F, row.names = F)
write.table(ovlp_ci_v_ic_hgnc, "C:/Users/miles/Downloads/brain/results/ase_ovlp_ci_v_ic_hgnc.txt", quote = F, col.names = F, row.names = F)

my_MBASED_1 = function(s1_mc, s1_cv, s1_name, gene_ind, gene_names, n_boot, myIsPhased=T, verbose=T) {
  # Purpose: Run a one sampled MBASED Experiment
  # s1_mc: sample 1 mc counts
  # s1_cv: sample 1 cv counts
  # gene_ind: index of genes to run on (aka subset of gene indexes), set to "" to find pos genes in this sample
  # gene_names: genes the counts are for (aka all genes)
  # n_boot: number of bootstraps in runMBASED
  
  
  # pos_ind = 1:length(gene_names)
  pos_ind = gene_ind
  if (pos_ind == "") {
    pos_ind = which( s1_mc + s1_cv > 0)
  }
  pos_gene = gene_names[pos_ind]
  this_s1_mc = s1_mc[pos_ind]
  this_s1_cv = s1_cv[pos_ind]
  if (verbose) {
    print(paste("Genes Used", length(pos_gene)))
  }
  
  # Create the SummarizedExperiment and run MBASED
  my_granges = GRanges(seqnames = rep("chr1:1-2", length(pos_gene)), aseID=pos_gene)
  # lociAllele1Counts
  s1_exp = SummarizedExperiment(assays=list(
    lociAllele1Counts = matrix( c(this_s1_mc), ncol=1, dimnames = list(pos_gene, s1_name)),
    lociAllele2Counts = matrix( c(this_s1_cv), ncol=1, dimnames = list(pos_gene, s1_name))
  ), rowRanges = my_granges)
  s1 = runMBASED(ASESummarizedExperiment=s1_exp, isPhased = myIsPhased, numSim = n_boot)
  
  # Analyze MBASED Data
  # hist(assays(s1)$majorAlleleFrequencyDifference, main=paste(s1_name, "MAF"), xlab = "Major Allele Frequency")
  # hist(assays(s1)$pValueASE, main=paste(s1_name, "p-value"), xlab = "p-value")
  qvalue = p.adjust(assays(s1)$pValueASE, method="bonferroni")
  s1_genes = pos_gene[which(qvalue < 0.05)]
  
  return(list(s1, s1_genes))
}

my_MBASED = function(s1_mc, s1_cv, s2_mc, s2_cv, s1_name, s2_name, gene_names, n_boot, myIsPhased=T, verbose=T, isSNP=F) {
  # Purpose: Run a two sampled MBASED Experiment
  # s1_mc: sample 1 mc counts
  # s1_cv: sample 1 cv counts
  # s2_mc: sample 2 mc counts
  # s2_cv: sample 2 cv counts
  # s1_name: name of sample 1 (for example "pit")
  # s2_name: name of sample 2 (for example "castle")
  # gene_names: genes the counts are for
  # n_boot: number of bootstraps in runMBASED
  
  # First find non-zero loci bc according to the documentation:
  # "All supplied loci must have total read count (across both alleles) greater than 0 
  # (in each of the two samples, in the case of two-sample analysis)."
  if (isSNP) {
    this_s1_mc = s1_mc[which(names(s1_mc) %in% names(s2_mc))]
    this_s1_cv = s1_cv[which(names(s1_cv) %in% names(s2_cv))]
    this_s2_mc = s2_mc[which(names(s2_mc) %in% names(s1_mc))]
    this_s2_cv = s2_cv[which(names(s2_cv) %in% names(s1_cv))]
    print(paste("SNPs lost from s1:", length(s1_mc) - length(this_s1_mc)))
    print(paste("SNPs lost from s2:", length(s2_mc) - length(this_s2_mc)))
    pos_gene = names(s1_mc)[which(names(s1_mc) %in% names(s2_mc))]
  } else {
    pos_ind = which( s1_mc + s1_cv > 0 & s2_mc + s2_cv > 0 )
    pos_gene = gene_names[pos_ind]
    this_s1_mc = s1_mc[pos_ind]
    this_s1_cv = s1_cv[pos_ind]
    this_s2_mc = s2_mc[pos_ind]
    this_s2_cv = s2_cv[pos_ind]
    if (verbose) {
      print(paste("Genes Used", length(pos_gene)))
    }
  }
 
  
  # Create the SummarizedExperiment and run MBASED
  my_granges = GRanges(seqnames = rep("chr1:1-2", length(pos_gene)), aseID=pos_gene)
  s1_v_s2_exp = SummarizedExperiment(assays=list(
    lociAllele1Counts = matrix( c(this_s1_mc, this_s2_mc), ncol=2, dimnames = list(pos_gene, c(s1_name, s2_name))),
    lociAllele2Counts = matrix( c(this_s1_cv, this_s2_cv), ncol=2, dimnames = list(pos_gene, c(s1_name, s2_name)))
  ), rowRanges = my_granges)
  s1_v_s2 = runMBASED(ASESummarizedExperiment=s1_v_s2_exp, isPhased = myIsPhased, numSim = n_boot)
  
  # Analyze MBASED Data
  hist(assays(s1_v_s2)$majorAlleleFrequencyDifference, main=paste(s1_name, "v", s2_name, "MAF"), xlab = "Major Allele Frequency")
  hist(assays(s1_v_s2)$pValueASE, main=paste(s1_name, "v", s2_name, "p-value"), xlab = "p-value")
  qvalue = p.adjust(assays(s1_v_s2)$pValueASE, method="bonferroni")
  s1_v_s2_genes = pos_gene[which(qvalue < 0.05)]
  
  return(list(s1_v_s2, s1_v_s2_genes))
}

posToGene = function(all_pos, gtf) {
  found_gene = c()
  for (pos in all_pos) {
    stop_1 = gregexpr(pattern = ':', pos)[[1]]
    stop_2 = gregexpr(pattern = '-', pos)[[1]]
    lg = substr(pos, 1, stop_1-1)
    base = substr(pos, stop_1+1, stop_2-1)
    
    this_found = gtf$gene_name[which(gtf$LG == lg & gtf$start+25000 <= base & gtf$stop+25000 >= base)]
    found_gene = c(found_gene, this_found)
  }
  
  return(found_gene)
}

shuffleAlleles = function(s1_mc, s1_cv, s2_mc, s2_cv) {
  all_mc = data.frame(s1_mc, s2_mc)
  ind1 = sample(c(1,2), length(s1_mc), replace = T)
  ind2 = ind1
  ind2 = factor(ind1, levels = c("1", "2"))
  ind2 = plyr::revalue(ind2, replace = c("1" = "2", "2" = "1"))
  new_s1_mc = all_mc[as.matrix(data.frame(1:nrow(all_mc), as.numeric(as.vector(ind1))))]
  new_s2_mc = all_mc[as.matrix(data.frame(1:nrow(all_mc), as.numeric(as.vector(ind2))))]
  
  all_cv = data.frame(s1_cv, s2_cv)
  ind1 = sample(c(1,2), length(s1_cv), replace = T)
  ind2 = ind1
  ind2 = factor(ind1, levels = c("1", "2"))
  ind2 = plyr::revalue(ind2, replace = c("1" = "2", "2" = "1"))
  new_s1_cv = all_cv[as.matrix(data.frame(1:nrow(all_cv), as.numeric(as.vector(ind1))))]
  new_s2_cv = all_cv[as.matrix(data.frame(1:nrow(all_cv), as.numeric(as.vector(ind2))))]
  
  res = data.frame(new_s1_mc, new_s1_cv, new_s2_mc, new_s2_cv)
  
  return(res)
}

#===============#
# Bootstrapping #
#===============#
real_pc = length(pit_v_castle_genes)
real_cp = length(castle_v_pit_genes)
real_ovlp_pc_v_cp = length(ovlp_pc_v_cp)

boot_res = data.frame()
for (n in 1:n_boot) {
  if(n == n_boot) {
    cat(paste(n, "\n"))
  } else if (n %% (n_boot/10) == 0 || n == 1) {
    cat(n)
  } else {
    cat(".")
  }
  
  tryCatch({
    # Pit v Castle
    shuf_res = shuffleAlleles(pit_mc, pit_cv, castle_mc, castle_cv)
    pit_v_castle_res = my_MBASED(shuf_res$new_s1_mc, shuf_res$new_s1_cv, shuf_res$new_s2_mc, shuf_res$new_s2_cv, "pit", "castle", gene_names, n_boot, verbose=F)
    pit_v_castle_genes = pit_v_castle_res[[2]]
    castle_v_pit_res = my_MBASED(shuf_res$new_s2_mc, shuf_res$new_s2_cv, shuf_res$new_s1_mc, shuf_res$new_s1_cv, "castle", "pit", gene_names, n_boot, verbose=F)
    castle_v_pit_genes = castle_v_pit_res[[2]]
    ovlp_pc_v_cp = pit_v_castle_genes[which(pit_v_castle_genes %in% castle_v_pit_genes)]
    
    # # Pit v Isolated
    # pit_v_iso_res = my_MBASED(pit_mc, pit_cv, iso_mc, iso_cv, "pit", "iso", gene_names, n_boot, verbose=F)
    # pit_v_iso_genes = pit_v_iso_res[[2]]
    # iso_v_pit_res = my_MBASED(iso_mc, iso_cv, pit_mc, pit_cv, "iso", "pit", gene_names, n_boot, verbose=F)
    # iso_v_pit_genes = iso_v_pit_res[[2]]
    # ovlp_pi_v_ip = pit_v_iso_genes[which(pit_v_iso_genes %in% iso_v_pit_genes)]
    # 
    # # Castle v Isolated
    # castle_v_iso_res = my_MBASED(castle_mc, castle_cv, iso_mc, iso_cv, "castle", "iso", gene_names, n_boot, verbose=F)
    # castle_v_iso_genes = castle_v_iso_res[[2]]
    # iso_v_castle_res = my_MBASED(iso_mc, iso_cv, castle_mc, castle_cv, "iso", "castle", gene_names, n_boot, verbose=F)
    # iso_v_castle_genes = iso_v_castle_res[[2]]
    # ovlp_ci_v_ic = castle_v_iso_genes[which(castle_v_iso_genes %in% iso_v_castle_genes)]
    # boot_res = rbind(boot_res, t(c(n, ovlp_pc_v_cp, ovlp_pi_v_ip, ovlp_ci_v_ic)))
    boot_res = rbind(boot_res, t(c(n, length(ovlp_pc_v_cp))))
  }, error = function(e) {
    print(paste("Error on boostrap", n))
  })
  
}
# colnames(boot_res) = c("run", "overlap_in_pvc_and_cvp", "overlap_in_pvi_and_ivp", "overlap_in_cvi_and_ivc")
colnames(boot_res) = c("run", "overlap_in_pvc_and_cvp")


boot_res$above = boot_res$overlap_in_pvc_and_cvp > real_ovlp_pc_v_cp
ggplot(boot_res, aes(overlap_in_pvc_and_cvp, alpha=.7, fill=above)) + geom_histogram(alpha=0.5, color = "purple") + geom_vline(aes(xintercept = real_ovlp_pc_v_cp)) + geom_text(aes(x=real_ovlp_pc_v_cp, label="Real Value"), y = Inf, hjust=0, vjust=1, color = "black") + xlab("# of Gene in Overlap Between Pit v Castle and Castle v Pit") + ggtitle("Comparison Between Bootstrap Values and Real Value") + guides(color=F, alpha=F, fill=F)

print(paste("p-value =", length(boot_res$above[which(boot_res$above)]) / length(boot_res$above)))

#=========================================================================================
# Old UMD1 Data
#=========================================================================================
rna_path <- "C:/Users/miles/Downloads/brain/"
data <- read.table(paste(rna_path, "/data/disc_ase.txt", sep=""), header = TRUE)

disc_genes <- c()
for (gene in unique(data$gene)) {
  this_rows <- data[which(data$gene == gene),]
  
  if (gene == "atp1b4") {
    print(this_rows)
  }
  
  if (this_rows$rep_1_ase_ratio[1] > 0 && this_rows$rep_2_ase_ratio[1] > 0 && nrow(this_rows) >= 2) { # both pos
    for (i in 2:nrow(this_rows)) {
      if (this_rows$rep_1_ase_ratio[i] < 0 && this_rows$rep_2_ase_ratio[i] < 0) {
        disc_genes <- c(disc_genes, gene)
      }
    }
  } else if (this_rows$rep_1_ase_ratio[1] < 0 && this_rows$rep_2_ase_ratio[1] < 0 && nrow(this_rows) >= 2) { # both neg
    for (i in 2:nrow(this_rows)) {
      if (this_rows$rep_1_ase_ratio[i] > 0 && this_rows$rep_2_ase_ratio[i] > 0) {
        disc_genes <- c(disc_genes, gene)
      }
    }
  }
  
}

mc_up <- c()
for (gene in unique(data$gene)) {
  this_rows <- data[which(data$gene == gene),]
  build_rows <- this_rows[which(this_rows$condition == "building"),]
  iso_rows <- this_rows[which(this_rows$condition == "isolated"),]
  dig_rows <- this_rows[which(this_rows$condition == "digging"),]
  min_build <- min(build_rows$rep_1_ase_ratio, build_rows$rep_2_ase_ratio)
  
  if (nrow(iso_rows) > 0 && nrow(dig_rows) > 0) { # only both up is considered mc_up
    if (iso_rows$rep_1_ase_ratio[i] < min_build && iso_rows$rep_2_ase_ratio[i] < min_build && dig_rows$rep_1_ase_ratio[i] < min_build && dig_rows$rep_2_ase_ratio[i] < min_build) {
      mc_up <- c(mc_up, gene)
    }
  } else { # either one up, is considered mc_up
    if (nrow(iso_rows) > 0 && iso_rows$rep_1_ase_ratio[i] < min_build && iso_rows$rep_2_ase_ratio[i] < min_build) {
      mc_up <- c(mc_up, gene)
    }
    if (nrow(dig_rows) > 0 && dig_rows$rep_1_ase_ratio[i] < min_build && dig_rows$rep_2_ase_ratio[i] < min_build) {
      mc_up <- c(mc_up, gene)
    }
  }
}
df <- data.frame(gene <- mc_up, bio <- rep("MC_UP", length(mc_up)))
write.table(df, paste(rna_path, "/data/mc_up.txt", sep=""), sep="\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

data = read.csv("C:/Users/miles/Downloads/cichlid_ase_common_genes_all_conditions_filtered_030920.csv", header = T)
test = data[which( sign(data$Digging_Mean_ASE-1) != sign(data$Building_Mean_ASE-1) ),1]


#==============================================================================================
# Single Nuc ASE ==============================================================================
#==============================================================================================
bb = readRDS("~/scratch/brain/data/bb_clustered_102820.rds")
counts = read.table("~/scratch/brain/ase/counts.txt", sep = "\t", header = T, stringsAsFactors=F)

bb_backup = bb

mat_ref = matrix(0L, nrow=nrow(bb_backup), ncol=ncol(bb_backup), dimnames = list(rownames(bb_backup), colnames(bb_backup)))
mat_alt = matrix(0L, nrow=nrow(bb_backup), ncol=ncol(bb_backup), dimnames = list(rownames(bb_backup), colnames(bb_backup)))

for (i in 1:nrow(counts)) {
  if (i%%nrow(counts)/10 == 0) { print(i) }
  gene = counts$GENE[i]
  cell = counts$CELL[i]
  mat_ref[gene, cell] = mat_ref[gene, cell] + counts$REF_COUNT[i]
  mat_alt[gene, cell] = mat_alt[gene, cell] + counts$ALT_COUNT[i]
}
saveRDS(mat_ref, "~/scratch/brain/ase/R/ref_mat.rds")
saveRDS(mat_alt, "~/scratch/brain/ase/R/alt_mat.rds")

# Find the average ASE for the cluster and set the numbers to be the same for every cell in the cluster
mat_clust_ref_alt_15 = matrix(0L, nrow=nrow(bb_backup), ncol=ncol(bb_backup), dimnames = list(rownames(bb_backup), colnames(bb_backup)))
mat_clust_log_ref_alt_15 = matrix(0L, nrow=nrow(bb_backup), ncol=ncol(bb_backup), dimnames = list(rownames(bb_backup), colnames(bb_backup)))
Idents(bb) = bb$seuratclusters15
for (cluster in 0:14) {
  cat(paste0(cluster, "."))
  clust_cells = WhichCells(bb, idents = cluster)
  bhve_cells = colnames(bb)[which(bb$cond == "BHVE")]
  ctrl_cells = colnames(bb)[which(bb$cond == "CTRL")]
  clust_b = clust_cells[which(clust_cells %in% bhve_cells)]
  clust_c = clust_cells[which(clust_cells %in% ctrl_cells)]
  # ase_ref_means = rowMeans(mat_ref[,clust_cells])
  # ase_alt_means = rowMeans(mat_alt[,clust_cells])
  ase_ref_sums_b = rowSums(mat_ref[,clust_b])
  ase_ref_sums_c = rowSums(mat_ref[,clust_c])
  ase_alt_sums_b = rowSums(mat_alt[,clust_b])
  ase_alt_sums_c = rowSums(mat_alt[,clust_c])
  # mat_clust_ref_alt_15[,clust_cells] = matrix( rep(ase_ref_means/ase_alt_means, length(clust_cells)), ncol = length(clust_cells) )
  # mat_clust_ref_alt_15[,clust_cells] = matrix( rep(ase_ref_means/ase_alt_means, length(clust_cells)), ncol = length(clust_cells) )
  mat_clust_log_ref_alt_15[,clust_b] = matrix( rep(log2(ase_ref_sums_b/ase_alt_sums_b), length(clust_b)), ncol = length(clust_b) )
  mat_clust_log_ref_alt_15[,clust_c] = matrix( rep(log2(ase_ref_sums_c/ase_alt_sums_c), length(clust_c)), ncol = length(clust_c) )
}

# bb@assays$RNA@data = mat_clust_log_ref_alt

png_name = "~/scratch/brain/ase/R/log_drd2.png"
png(file = png_name, width = 1000, height = 1000, res = 150)
print(FeaturePlot(bb, "drd2", order = T, pt.size = 1, label =T))
dev.off()
system(paste0("rclone copy ", png_name, " dropbox:BioSci-Streelman/George/Brain/bb/results/sn_ase/"))

bb@assays$RNA@data = bb_backup@assays$RNA@data
png_name = "~/scratch/brain/ase/R/drd2_exp.png"
png(file = png_name, width = 1000, height = 1000, res = 150)
print(FeaturePlot(bb, "drd2", order = T, pt.size = 1, label =T))
dev.off()
system(paste0("rclone copy ", png_name, " dropbox:BioSci-Streelman/George/Brain/bb/results/sn_ase/"))


pos = which(rowSums(mat_ref) > 0 & rowSums(mat_alt) > 0)
d = density(log2( rowSums(mat_ref[pos,]) / rowSums(mat_alt[pos,]) ))
png_name = "~/scratch/brain/ase/R/all_pos.png"
png(file = png_name, width = 1000, height = 1000, res = 150)
print(plot(d, main="All"))
dev.off()
system(paste0("rclone copy ", png_name, " dropbox:BioSci-Streelman/George/Brain/bb/results/sn_ase/"))

for (i in 0:14) {
  clust_cells = WhichCells(bb, idents = cluster)
  pos = which(rowSums(mat_ref[,clust_cells]) > 0 & rowSums(mat_alt[,clust_cells]) > 0)
  d = density(log2( rowSums(mat_ref[pos,clust_cells]) / rowSums(mat_alt[pos,clust_cells]) ))
  png_name = paste0("~/scratch/brain/ase/R/", i, "_pos.png")
  png(file = png_name, width = 1000, height = 1000, res = 150)
  print(plot(d, main=i))
  dev.off()
  system(paste0("rclone copy ", png_name, " dropbox:BioSci-Streelman/George/Brain/bb/results/sn_ase/"))
}

# good_genes = rownames(bb)
# for (i in 0:14) {
#   clust_cells = WhichCells(bb, idents = cluster)
#   ase_ref_sums = rowSums(mat_ref[good_genes, clust_cells])
#   ase_alt_sums = rowSums(mat_alt[good_genes, clust_cells])
#   good_genes = good_genes[which(ase_ref_sums > 5 & ase_alt_sums > 5)]
# }
# length(good_genes)

Idents(bb) = bb$cond
good_genes = rownames(bb)[which(rowSums(mat_ref) > 200 & rowSums(mat_alt) > 200)]
n_boot = 100
# sig_df = data.frame()
# sig_genes = list()
# for (i in 0:14) {
#   for (j in i:14) {
i_cells = WhichCells(bb, idents = "BHVE")
i_ref = rowSums(mat_ref[good_genes, i_cells])
i_alt = rowSums(mat_alt[good_genes, i_cells])
j_cells = WhichCells(bb, idents = "CTRL")
j_ref = rowSums(mat_ref[good_genes, j_cells])
j_alt = rowSums(mat_alt[good_genes, j_cells])

i_j_res = my_MBASED(i_ref, i_alt, j_ref, j_alt, i, j, good_genes, n_boot)
i_j_genes = i_j_res[[2]]
j_i_res = my_MBASED(j_ref, j_alt, i_ref, i_alt, j, i, good_genes, n_boot)
j_i_genes = j_i_res[[2]]
sig_genes = i_j_genes[which(i_j_genes %in% j_i_genes)]
# sig_genes[[paste0(i, "_", j)]] = i_j_genes[which(i_j_genes %in% j_i_genes)]
# sig_df = rbind(sig_df, t(c(i, j, length(i_j_genes), length(j_i_genes), length(sig_genes[[paste0(i, "_", j)]]))))
# colnames(sig_df) = c("cluster_A", "cluster_B", "A_B_genes", "B_A_gnes", "ovlp")
# sig_df = sig_df[which(sig_df$cluster_A != sig_df$cluster_B),]

all_sig_genes = unique(sig_genes)
Idents(bb) = bb$seuratclusters15
for (gene in all_sig_genes) {
  sig_in = unlist(sapply(1:length(sig_genes), function(x) if(gene %in% sig_genes[[x]]) {names(sig_genes)[x]}))
  bb@assays$RNA@data = mat_clust_log_ref_alt_15
  png_name = paste0("~/scratch/brain/ase/R/log_", gene, ".png")
  png(file = png_name, width = 2000, height = 1000, res = 150)
  print(FeaturePlot(bb, gene, order = T, pt.size = 1, label =T, split.by = "cond") + labs(caption = paste0("ASE Ratio in BHVE: ", log2(sum(mat_ref[gene, bhve_cells])/sum(mat_alt[gene, bhve_cells])), ". ASE Ratio in CTRL: ", log2(sum(mat_ref[gene, ctrl_cells])/sum(mat_alt[gene, ctrl_cells])) )))
  dev.off()
  system(paste0("rclone copy ", png_name, " dropbox:BioSci-Streelman/George/Brain/bb/results/sn_ase/genes_bvc/"))
  
  bb@assays$RNA@data = bb_backup@assays$RNA@data
  png_name = paste0("~/scratch/brain/ase/R/exp_", gene, ".png")
  png(file = png_name, width = 1000, height = 1000, res = 150)
  print(FeaturePlot(bb, gene, order = T, pt.size = 1, label =T))
  dev.off()
  system(paste0("rclone copy ", png_name, " dropbox:BioSci-Streelman/George/Brain/bb/results/sn_ase/genes_bvc/"))
}
