library(biomaRt)
library(rentrez)
library(dplyr)
library(stringr)

# Helper Functions
geneCap <- function(gene, gene_names) {
  # Gene the gene name in the right format
  gene_lower <- tolower(gene)
  gene_upper <- toupper(gene)
  gene_title <- str_to_title(gene)
  error <- FALSE
  if (gene_lower %in% gene_names) {
    gene <- gene_lower
  } else if (gene_upper %in% gene_names) {
    gene <- gene_upper
  } else if (gene_title %in% gene_names) {
    gene <- gene_title
  } else {
    error <- TRUE
  }
  
  return(c(gene, error))
}

validGenes <- function(genes, gene_names) {
  valid_genes <- c()
  for (gene in genes) {
    result <- geneCap(gene, gene_names)
    gene <- result[1]
    error <- as.logical(result[2])
    if (! error) {
      valid_genes <- c(valid_genes, gene)
    }
  } # end gene for
  valid_genes <- unique(valid_genes)
  return(valid_genes)
} # end validGenes function

convertToHgnc <- function(genes) {
  human  = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mzebra = useMart(biomart="ensembl", dataset="mzebra_gene_ensembl")
  
  ensembl_genes <- getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = genes , mart = mzebra, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  zfin_genes    <- getLDS(attributes = c("zfin_id_symbol"), filters = "zfin_id_symbol", values = genes , mart = mzebra, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  hgnc_genes    <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = genes , mart = mzebra, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  all_hgnc <- rbind(ensembl_genes, setNames(zfin_genes, names(ensembl_genes)), setNames(hgnc_genes, names(ensembl_genes)))
  all_hgnc = all_hgnc[!duplicated(all_hgnc[,1]),]
  colnames(all_hgnc) <- c("mzebra", "hgnc")
  
  # all_hgnc <- unique(c(ensembl_genes[,2], zfin_genes[,2], hgnc_genes[,2]))
  return(all_hgnc)
}

hgncMzebra <- function(genes, gene_names) {
  pat <- read.table("C:/Users/miles/Downloads/all_research/MZ_treefam_annot_umd2a_ENS_2.bash", header = FALSE, fill = TRUE)
  valid_genes <- validGenes(genes, gene_names)
  all_hgnc <- convertToHgnc(gene_names)
  ind <- match(genes,pat$V2)
  ind <- ind[! is.na(ind)]
  found_names <- as.vector(pat$V7[ind])
  found_names <- found_names[!is.na(found_names)]
  found_names_hgnc <- as.vector(pat$V8[ind])
  found_names_hgnc <- found_names_hgnc[!is.na(found_names_hgnc)]
  
  df1 <- setNames(as.data.frame(found_names_hgnc), c("hgnc"))
  found_names_hgnc <- inner_join(df1, all_hgnc, by = "hgnc")
  
  pseudo_hgnc <- toupper(genes)
  df1 <- setNames(as.data.frame(pseudo_hgnc), c("hgnc"))
  found_mzebra <- inner_join(df1, all_hgnc, by = "hgnc")
  
  found_mzebra <- found_mzebra[,2:1]
  found_names_hgnc <- found_names_hgnc[,2:1]
  good_df <- rbind(all_hgnc, setNames(found_names, names(all_hgnc)), setNames(found_mzebra, names(all_hgnc)), setNames(found_names_hgnc, names(all_hgnc)))
  good_df <- unique(good_df)
  return(good_df)
}

hgncMouse <- function(genes) {
  human  = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse  = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
  
  ensembl_genes <- getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = genes , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(ensembl_genes)
}

hgncMzebraInPlace <- function(df, gene_column, gene_names) {
  bad_genes <- df[,gene_column]
  bad_genes <- unique(bad_genes)
  converter <- hgncMzebra(bad_genes, gene_names)
  df[,gene_column] <- converter[match(df[,gene_column], converter[,1]),2]
  df <- df[which(! is.na(df[,gene_column])),]
  
  return(df)
}

hgncMouseInPlace <- function(df, gene_column) {
  bad_genes <- df[,gene_column]
  bad_genes <- unique(bad_genes)
  converter <- hgncMouse(bad_genes)
  df[,gene_column] <- converter[match(df[,gene_column], converter[,1]),2]
  df <- df[which(! is.na(df[,gene_column])),]
  
  return(df)
}

# THIS IS THE UPDATED/GOOD Function 02/21/2020
hgncGood <- function(genes, gene_names) {
  pat <- read.table("C:/Users/miles/Downloads/all_research/MZ_treefam_annot_umd2a_ENS_2.bash", header = FALSE, fill = TRUE)
  valid_genes <- validGenes(genes, gene_names)
  all_hgnc <- convertToHgnc(gene_names)
  ind <- match(genes,pat$V2)
  ind <- ind[! is.na(ind)]
  found_names <- as.vector(pat$V7[ind])
  found_names <- found_names[!is.na(found_names)]
  found_names_hgnc <- as.vector(pat$V8[ind])
  found_names_hgnc <- found_names_hgnc[!is.na(found_names_hgnc)]
  
  df1 <- setNames(as.data.frame(found_names_hgnc), c("hgnc"))
  found_names_hgnc <- inner_join(df1, all_hgnc, by = "hgnc")$mzebra
  # ind_found_hgnc <- match(found_names_hgnc, all_hgnc$hgnc)
  # ind_found_hgnc <- ind_found_hgnc[! is.na(ind_found_hgnc)]
  # found_names_hgnc <- as.vector(all_hgnc[ind_found_hgnc,1])

  pseudo_hgnc <- toupper(genes)
  df1 <- setNames(as.data.frame(pseudo_hgnc), c("hgnc"))
  found_mzebra <- inner_join(df1, all_hgnc, by = "hgnc")$mzebra
  # ind_hgnc <- match(pseudo_hgnc, all_hgnc$hgnc)
  # ind_hgnc <- ind_hgnc[! is.na(ind_hgnc)]
  # found_mzebra <- as.vector(all_hgnc[ind_hgnc,1])
  
  good <- c(valid_genes, found_names, found_names_hgnc, found_mzebra)
  good <- good[which(good != "")]
  good <- validGenes(good, gene_names)
  good <- unique(good)
  good <- sort(good)
  return(good)
}

# This function is outdated
onlyGood <- function(genes, gene_names) {
  pat <- read.table("C:/Users/miles/Downloads/MZ_treefam_annot_umd2a_ENS_2.bash", header = FALSE, fill = TRUE)
  valid_genes <- validGenes(genes, gene_names)
  ind <- match(genes,pat$V2)
  ind <- ind[! is.na(ind)]
  found_names <- as.vector(pat$V7[ind])
  found_names <- found_names[!is.na(found_names)]
  found_names_hgnc <- as.vector(pat$V8[ind])
  found_names_hgnc <- found_names_hgnc[!is.na(found_names_hgnc)]
  good <- c(valid_genes, found_names, found_names_hgnc)
  good <- unique(good)
  good <- sort(good)
  good <- good[which(good != "")]
  return(good)
}

goodInPlace <- function(data, gene_column, gene_names) {
  # Only keep the rows in the dataframe that have gene names that are in our data
  pat <- read.table("C:/Users/miles/Downloads/all_research/MZ_treefam_annot_umd2a_ENS_2.bash", header = FALSE, fill = TRUE)
  keep_rows <- c()
  for (i in 1:nrow(data)) {
    gene <- as.character(data[i,gene_column])
    result <- geneCap(gene, gene_names)
    gene <- result[1]
    error <- as.logical(result[2])
    ens2_gene <- as.character(pat$V6[which(gene == pat$V2)])[1]
    if (!error) {
      data[i,gene_column] <- gene
      keep_rows <- c(keep_rows, i)
    }
    if (!identical(ens2_gene, character(0)) && ens2_gene %in% gene_names) {
      data[i,gene_column] <- ens2_gene
      keep_rows <- c(keep_rows, i)
    }
  }
  new_data <- data[keep_rows,]
  return(new_data)
}

heatmapComparison <- function(df1, df2, df1_sample, df2_sample, filename) {
  df1_num_clusters <- max(df1$cluster)
  mes_num_clusters <- max(df2$cluster)
  df <- data.frame()
  for (i in 0:df1_num_clusters) {
    for (j in 0:mes_num_clusters) {
      jpool_cluster <- df1[which(df1$cluster == i),]
      mes_cluster <- df2[which(df2$cluster == j),]
      ovlp <- nrow(mes_cluster[which(mes_cluster$gene %in% jpool_cluster$gene),])
      pct <- ovlp / (nrow(jpool_cluster) + nrow(mes_cluster))
      df <- rbind(df, t(c(i, j, ovlp, pct)))
    }
  }
  colnames(df) <- c("jpool_cluster", "mes_cluster", "ovlp", "pct")
  mat <- matrix(df$ovlp,ncol=mes_num_clusters+1,byrow=T)
  rownames(mat) <- 0:df1_num_clusters
  colnames(mat) <- 0:mes_num_clusters
  mat <- t(mat)
  png(paste(rna_path, "/results/", filename, "_ovlp.png", sep=""),  width = 500, height = 750, unit = "px")
  heatmap(mat, Rowv=NA, Colv=NA, revC = TRUE, xlab = paste(df1_sample, "Cluster"), ylab = paste(df2_sample, "Cluster"), main = paste("DEGs in Common b/w", df1_sample, "&", df2_sample,  "Clusters"), scale = "none")
  dev.off()
  for (i in 1:ncol(mat)) {
    max_row <- max(mat[,i])
    mat[which(mat[,i] != max_row),i] <- 0
  }
  png(paste(rna_path, "/results/", filename, "_best_guess.png", sep=""),  width = 500, height = 750, unit = "px")
  heatmap(mat, Rowv=NA, Colv=NA, revC = TRUE, xlab = paste(df1_sample, "Cluster"), ylab = paste(df2_sample, "Cluster"), main = paste("Best Guess b/w", df1_sample, "&", df2_sample), scale = "none")
  dev.off()
  # legend(x="topleft", legend=c("min", "ave", "max"), fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))
  
  pct_mat <- matrix(df$pct,ncol=mes_num_clusters+1,byrow=T)
  rownames(pct_mat) <- 0:df1_num_clusters
  colnames(pct_mat) <- 0:mes_num_clusters
  pct_mat <- t(pct_mat)
  png(paste(rna_path, "/results/", filename, "_pct.png", sep=""),  width = 500, height = 750, unit = "px")
  heatmap(pct_mat, Rowv=NA, Colv=NA, revC = TRUE, xlab = paste(df1_sample, "Cluster"), ylab = paste(df2_sample, "Cluster"), main = paste("% DEGs in Common b/w", df1_sample, "&", df2_sample,  "Clusters"), scale = "none")
  dev.off()
  for (i in 1:ncol(pct_mat)) {
    max_row <- max(pct_mat[,i])
    pct_mat[which(pct_mat[,i] != max_row),i] <- 0
  }
  png(paste(rna_path, "/results/", filename, "_pct_best_guess.png", sep=""),  width = 500, height = 750, unit = "px")
  heatmap(pct_mat, Rowv=NA, Colv=NA, revC = TRUE, xlab = paste(df1_sample, "Cluster"), ylab = paste(df2_sample, "Cluster"), main = paste("% Best Guess b/w", df1_sample, "&", df2_sample), scale = "none")
  dev.off()
  return(mat)
}

########################
# END HELPER FUNCTIONS #
########################
# data <- read.table("C:/Users/miles/Downloads/pnas.1810140115.sd02.txt", sep=",", header = FALSE)
# data <- read.table("C:/Users/miles/Downloads/brain/data/markers/mc_up.txt", sep="\t", header = FALSE)
# data <- readClipboard()

# For Zack
data <- read.table("C:/Users/miles/Downloads/all_markers_B1C1C2MZ_031020.csv", sep=",", header = TRUE, stringsAsFactors = FALSE)
good <- hgncGood(data$gene, rownames(combined@assays$RNA@counts))
new_data <- hgncMzebraInPlace(data, 8, rownames(combined@assays$RNA@counts))
write.table(new_data, file = "C:/Users/miles/Downloads/all_markers_B1C1C2MZ_031020_hgnc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

data <- read.table("C:/Users/miles/Downloads/d_tooth/results/jpool_deg.tsv", sep="\t", header = TRUE, stringsAsFactors = FALSE)
good <- hgncGood(data, rownames(combined@assays$RNA@counts))
new_data <- goodInPlace(data, 1, rownames(tj))

all_hgnc <- convertToHgnc(rownames(brain_combined@assays$RNA@counts))

clust_genes <- lapply(0:data[nrow(data),1], function(x) c())
for (i in 1:nrow(data)) {
  clust <- data[i,1]+1
  gene <- data[i,2]
  clust_genes[[clust]] <- c(clust_genes[[clust]], gene)
}
hgnc_clust_genes <- lapply(0:data[nrow(data),1], function(x) c())
for (i in 1:(data[nrow(data),1]+1)) {
  hgnc_clust_genes[[i]] <- convertToHgnc(clust_genes[[i]])
  write.table(hgnc_clust_genes[[i]], paste("C:/Users/miles/Downloads/d_tooth/results/jpool_toppgene/", i-1, ".txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
# good <- c("myrf", "mpz", "stxbp6l", "mbpa", "adamts14", "ENSMZEG00005026266", "ENSMZEG00005026264", "erbb3b", "ENSMZEG00005014976", "cd59")
df <- data.frame(genes <- good, bio <- rep("ALL_ASE", length(good)))
write.table(df, file = "C:/Users/miles/Downloads/brain/data/markers/rock_up.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

new_data <- new_data[,colnames(new_data)[which(! colnames(new_data) %in% c("X", "X.1"))]]
write.table(good, file = "C:/Users/miles/Downloads/tbud_deg_good_list.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(new_data, file = "C:/Users/miles/Downloads/tbud_deg_good.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Convert DEGs in Place
rna_path <- "C:/Users/miles/Downloads/rna/"
mes_deg <- read.table(paste(rna_path, "/results/clusters/mes_deg_1_5.tsv", sep=""), sep="\t", stringsAsFactors = FALSE, header = TRUE)
hgnc_mes_deg <- hgncMouseInPlace(mes_deg, 7)

epi_deg <- read.table(paste(rna_path, "/results/clusters/epi_deg_1_5.tsv", sep=""), sep="\t", stringsAsFactors = FALSE, header = TRUE)
hgnc_epi_deg <- hgncMouseInPlace(epi_deg, 7)

rna_path <- "C:/Users/miles/Downloads/d_tooth/"
jpool_deg <- read.table(paste(rna_path, "/results/jpool_deg.tsv", sep=""), sep="\t", stringsAsFactors = FALSE, header = TRUE)
hgnc_jpool_deg <- hgncMzebraInPlace(jpool_deg, 2, rownames(jpool))

tj_deg <- read.table(paste(rna_path, "/results/tj_deg.tsv", sep=""), sep="\t", stringsAsFactors = FALSE, header = TRUE)
hgnc_tj_deg <- hgncMzebraInPlace(tj_deg, 2, rownames(jpool))

# JPOOL vs MES Comparison
mat <- heatmapComparison(hgnc_jpool_deg, hgnc_mes_deg, "Jaw", "Mes", "jpool_mes_deg")

# JPOOL vs EPI Comparison
mat <- heatmapComparison(hgnc_jpool_deg, hgnc_epi_deg, "Jaw", "Epi", "jpool_epi_deg")

# TJ vs MES Comparison
mat <- heatmapComparison(hgnc_tj_deg, hgnc_mes_deg, "RT", "Mes", "tj_mes_deg")

# TJ vs EPI Comparison
mat <- heatmapComparison(hgnc_tj_deg, hgnc_epi_deg, "RT", "Epi", "tj_epi_deg")

# OLD #
# data <- as.vector(unique(data[,1]))
# pat <- read.table("C:/Users/miles/Downloads/MZ_treefam_annot_umd2a_ENS_2.bash", header = FALSE, fill = TRUE)
# found_names <- as.vector(pat$V7[which(data %in% pat$V2)])
# found_names <- found_names[!is.na(found_names)]
# data <- c(data, found_names)
# data <- unique(data)
# data <- sort(data)
# data <- data[which(data != "")]
# df <- data.frame(data <- data, bio <- rep("DISC_ASE_MC_UP_BOTH", length(data)))
# write.table(df, file = "C:/Users/miles/Downloads/brain/data/markers/disc_ase_mc_up_both.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 
# gtf <- read.table("C:/Users/miles/Downloads/brain/data/Maylandia_zebra.M_zebra_UMD2a.98.gtf", sep="\t", header = FALSE)
# gtf <- gtf[which(gtf$V3 == "gene"),]
# 
# ensembl = useEnsembl(biomart="ensembl", dataset="mzebra_gene_ensembl")
# 
# "uniprot_gn_symbol"
# "uniprot_gn_id"
# "entrezgene_id"
# 
# ensembl_names   <- unique(getBM(attributes=c('ensembl_gene_id'), filters = 'ensembl_gene_id', values = data, mart = ensembl))
# zfin_names      <- unique(getBM(attributes=c('zfin_id_symbol'), filters = 'zfin_id_symbol', values = data, mart = ensembl))
# hgnc_names      <- unique(getBM(attributes=c('hgnc_symbol'), filters = 'hgnc_symbol', values = data, mart = ensembl))
# external_names  <- unique(getBM(attributes=c('external_gene_name'), filters = 'external_gene_name', values = data, mart = ensembl))
# entrez_names    <- unique(getBM(attributes=c('entrezgene_id'), filters = 'entrezgene_id', values = data, mart = ensembl))
# uniprot_names   <- unique(getBM(attributes=c('uniprot_gn_symbol'), filters = "uniprot_gn_symbol", values = data, mart = ensembl))
# 
# test <- unique(c(tolower(ensembl_names[,1]), tolower(zfin_names[,1]), tolower(hgnc_names[,1]), tolower(external_names[,1]), tolower(entrez_names[,1]), tolower(uniprot_names[,1])))
# test <- sort(test)
# 
# new_genes <- c()
# for (j in 1:nrow(pat)) {
#   pat_gene <- as.character(pat[j,2])
#   # if (j)
#   if (startsWith(pat_gene, "LOC")) {
#     new_gene <- pat_gene
#     pat_lg <- as.character(pat[j, 3])
#     pat_start <- as.numeric(as.character(pat[j, 4]))
#     pat_stop <- as.numeric(as.character(pat[j, 5]))
#     
#     # gtf_rows <- gtf[which(gtf$V1 == pat_lg && abs(pat_start - gtf$v4) < 10000 && abs(pat_stop - gtf$v5) < 10000),]
#     gtf_rows <- gtf[which(gtf$V1 == pat_lg),]
#     gtf_rows$V4 <- abs(pat_start - gtf_rows$V4)
#     gtf_rows$V5 <- abs(pat_stop - gtf_rows$V5)
#     # gtf_rows <- gtf_rows[which(gtf_rows$V4 < 10000 && gtf_rows$V5 < 100000),]
#     gtf_rows <- gtf_rows[which(gtf_rows$V4 < 10000),]
#     gtf_rows <- gtf_rows[which(gtf_rows$V5 < 10000),]
#     
#     if (nrow(gtf_rows) >= 1) {
#       gtf_gene <- substring(as.character(gtf_rows[1,9]),9,26)
#       new_gene <- gtf_gene
#     }
#     
#     # for (i in 1:nrow(gtf)) {
#     #   gtf_gene <- substring(as.character(gtf[i,9]),9,26)
#     #   gtf_lg <- gtf[i,1]
#     #   gtf_start <- gtf[i,4]
#     #   gtf_stop <- gtf[i,5]
#     #   
#     #   if (pat_lg == gtf_lg) {
#     #     if ( abs(pat_start - gtf_start) < 10000 && abs(pat_stop - gtf_stop) < 10000 ) {
#     #       new_gene <- gtf_gene
#     #       break
#     #     }
#     #   }
#     # } # end gtf for
#   } else {
#     new_gene <- pat_gene
#   } # end if LOC
#   new_genes <- c(new_genes, new_gene)
# } # end pat for
# pat$new_gene <- new_genes
# write.table(pat, file = "C:/Users/miles/Downloads/my_MZ_treefam.bash", sep = "\t", row.names = FALSE, quote = FALSE)
