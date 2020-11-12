library(biomaRt)
library(dplyr)

# Helper Functions
getDescription <- function(genes) {
  ensembl = useEnsembl(biomart="ensembl", dataset="mzebra_gene_ensembl")
  
  ensembl_description_all   <- unique(getBM(attributes=c('ensembl_gene_id', 'description'), filters = 'ensembl_gene_id', values = unique(genes), mart = ensembl))
  zfin_description_all      <- unique(getBM(attributes=c('zfin_id_symbol', 'description'), filters = 'zfin_id_symbol', values = unique(genes), mart = ensembl))
  hgnc_description_all      <- unique(getBM(attributes=c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = unique(genes), mart = ensembl))
  
  colnames(ensembl_description_all)   <- c('all_markers', 'all_markers_description')
  colnames(zfin_description_all)      <- c('all_markers', 'all_markers_description_2')
  colnames(hgnc_description_all)      <- c('all_markers', 'all_markers_description_3')
  
  ensembl_description_all = ensembl_description_all[!duplicated(ensembl_description_all$all_markers),]
  zfin_description_all    = zfin_description_all[!duplicated(zfin_description_all$all_markers),]
  hgnc_description_all    = hgnc_description_all[!duplicated(hgnc_description_all$all_markers),]
  
  total <- as.data.frame(genes)
  colnames(total) <- c("all_markers")
  try(total <- left_join(total, ensembl_description_all, by = c("all_markers")), silent = TRUE)
  try(total <- left_join(total, zfin_description_all, by = c("all_markers")), silent = TRUE)
  try(total <- left_join(total, hgnc_description_all, by = c("all_markers")), silent = TRUE)
  
  try(total$all_markers_description[is.na(total$all_markers_description)] <- "", silent = TRUE)
  try(total$all_markers_description_2[is.na(total$all_markers_description_2)] <- "", silent = TRUE)
  try(total$all_markers_description_3[is.na(total$all_markers_description_3)] <- "", silent = TRUE)
  
  all_markers_description  = paste(total$all_markers_description,  total$all_markers_description_2,  total$all_markers_description_3)
  # total$all_markers_description  <- all_markers_description
  # drops <- c("all_markers_description_2","cons_markers_description_2", "all_markers_description_3","cons_markers_description_3")
  # total <- total[ , !(names(total) %in% drops)]
  
  return(all_markers_description)
}

# For Brain
data <- read.table("C:/Users/miles/Downloads/brain/results/c41_conserved_and_all_markers_by_cluster_121319.csv", sep=",", header = TRUE)
names(data)[1]<-"cluster"
data$all_markers_description  <- getDescription(data$all_markers)
data$cons_markers_description <- getDescription(data$conserved_markers)
write.csv(data, "C:/Users/miles/Downloads/brain/results/gene_description.csv")

# For d_tooth
tj.nk.markers$description <- getDescription(tj.nk.markers$gene)
jpool.nk.markers$description <- getDescription(jpool.nk.markers$gene)
lj.markers$description <- getDescription(rownames(lj.markers))
uj.markers$description <- getDescription(rownames(uj.markers))
write.table(tj.nk.markers, "C:/Users/miles/Downloads/d_tooth/results/tj_deg_description.tsv", row.names = FALSE, quote = FALSE, sep="\t")
write.table(jpool.nk.markers, "C:/Users/miles/Downloads/d_tooth/results/jpool_deg_description.tsv", row.names = FALSE, quote = FALSE, sep="\t")
write.table(lj.markers, "C:/Users/miles/Downloads/d_tooth/results/lj_deg_description.tsv", row.names = TRUE, quote = FALSE, sep="\t")
write.table(uj.markers, "C:/Users/miles/Downloads/d_tooth/results/uj_deg_description.tsv", row.names = TRUE, quote = FALSE, sep="\t")


# gene_description <- c()
# ensembl = useEnsembl(biomart="ensembl", dataset="mzebra_gene_ensembl")
# 
# ensembl_description_all   <- unique(getBM(attributes=c('ensembl_gene_id', 'description'), filters = 'ensembl_gene_id', values = unique(data$all_markers), mart = ensembl))
# zfin_description_all      <- unique(getBM(attributes=c('zfin_id_symbol', 'description'), filters = 'zfin_id_symbol', values = unique(data$all_markers), mart = ensembl))
# hgnc_description_all      <- unique(getBM(attributes=c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = unique(data$all_markers), mart = ensembl))
# ensembl_description_cons  <- unique(getBM(attributes=c('ensembl_gene_id', 'description'), filters = 'ensembl_gene_id', values = unique(data$conserved_markers), mart = ensembl))
# zfin_description_cons     <- unique(getBM(attributes=c('zfin_id_symbol', 'description'), filters = 'zfin_id_symbol', values = unique(data$conserved_markers), mart = ensembl))
# hgnc_description_cons     <- unique(getBM(attributes=c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = unique(data$conserved_markers), mart = ensembl))
# 
# colnames(ensembl_description_all)   <- c('all_markers', 'all_markers_description')
# colnames(zfin_description_all)      <- c('all_markers', 'all_markers_description_2')
# colnames(hgnc_description_all)      <- c('all_markers', 'all_markers_description_3')
# colnames(ensembl_description_cons)  <- c('conserved_markers', 'cons_markers_description')
# colnames(zfin_description_cons)     <- c('conserved_markers', 'cons_markers_description_2')
# colnames(hgnc_description_cons)     <- c('conserved_markers', 'cons_markers_description_3')
# 
# total <- left_join(data, ensembl_description_all, by = c("all_markers"))
# total <- left_join(total, zfin_description_all, by = c("all_markers"))
# total <- left_join(total, hgnc_description_all, by = c("all_markers"))
# total <- left_join(total, ensembl_description_cons, by = c("conserved_markers"))
# total <- left_join(total, zfin_description_cons, by = c("conserved_markers"))
# total <- left_join(total, hgnc_description_cons, by = c("conserved_markers"))
# 
# total$all_markers_description[is.na(total$all_markers_description)] <- ""
# total$all_markers_description_2[is.na(total$all_markers_description_2)] <- ""
# total$all_markers_description_3[is.na(total$all_markers_description_3)] <- ""
# total$cons_markers_description[is.na(total$cons_markers_description)] <- ""
# total$cons_markers_description_2[is.na(total$cons_markers_description_2)] <- ""
# total$cons_markers_description_3[is.na(total$cons_markers_description_3)] <- ""
# 
# all_markers_description  = paste(total$all_markers_description,  total$all_markers_description_2,  total$all_markers_description_3)
# cons_markers_description = paste(total$cons_markers_description, total$cons_markers_description_2, total$cons_markers_description_3)
# total$all_markers_description  <- all_markers_description
# total$cons_markers_description <- cons_markers_description

# drops <- c("all_markers_description_2","cons_markers_description_2", "all_markers_description_3","cons_markers_description_3")
# total <- total[ , !(names(total) %in% drops)]
