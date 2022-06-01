#====================================================================================================
# bcftools merge method =============================================================================
#====================================================================================================
# # Pit Castle Stuff
# 47 million variants in bower tangankyan species. 37 million SNPs and 10 million indels. 
# 
# 22 million variants from original pit v castle. 18 million SNPs and 4 million indels.
# 
# 2.8 million SNPs in common between malawi and tanganyka for pit-castle.
# 
# From merged file, 66 million sites were used to calculate FST.

# Pre-process, clean, and select high FST Malawi Pit vs Castle Sites
fst = read.table("~/scratch/brain/fst/tang/fst_pc_malawi_only_on_merged_w_vcf.fst", header = F)
fst_header = scan("~/scratch/brain/fst/tang/malawi_tang_merged_header.vcf", what = character())
fst_header[1] = "CHROM"
fst_header = c(c("CHROM_FST", "POS_FST", "FST"), fst_header)
fst_header[13:27] = substr(fst_header[13:27], 8, 1000L)
colnames(fst) = fst_header
#Clean
fst = fst[complete.cases(fst$FST),]
fst[which(fst$FST < 0),3] = 0
thresh95_clean = quantile(fst$FST, 0.95)
write.table(fst_clean, "~/scratch/brain/fst/tang/fst_pc_malawi_only_on_merged_w_vcf_clean.fst")

library("reshape2")
high_fst = read.table("~/scratch/brain/fst/tang/fst_pc_malawi_only_on_merged_w_vcf_clean_high.fst")

miss_mat = high_fst[9:23] == "./.:.:.:.:."
num_tang_miss = rowSums(miss_mat)
high_fst_max_miss = high_fst[which(num_tang_miss < (23-9)/2 ),]
# high_fst_max_miss = high_fst_max_miss[which(high_fst_max_miss$FST >= 0.90),]

# PCA Method from: https://comppopgenworkshop2019.readthedocs.io/en/latest/contents/03_pca/pca.html
high_fst3 = high_fst_max_miss[1:8]
for (i in 9:51) {
  allele1 = substr(high_fst_max_miss[,i], 0, 1)
  allele2 = substr(high_fst_max_miss[,i], 3, 3)
  
  isRef1 = allele1 == 0
  isHomo = allele1 == allele2
  isMissing1 = allele1 == "."
  
  score = rep(-1, nrow(high_fst_max_miss))
  score[which(isMissing1)] = 9
  score[which(isHomo & isRef1)] = 2
  score[which(isHomo & !isRef1)] = 0
  score[which(! isHomo & ! isMissing1 )] = 1

  high_fst3 = cbind(high_fst3, score)
}
colnames(high_fst3) = colnames(high_fst_max_miss)
res.pca <- prcomp(t(high_fst3[9:51]))

library(factoextra)
sample_df = data.frame(samples = c("SRR9657485", "SRR9657511", "SRR9657512", "SRR9657530", "SRR9657540", "SRR9657558", "SRR9657559", "SRR9665654", "SRR9665657", "SRR9665678", "SRR9665679", "SRR9665706", "SRR9673822", "SRR9673823", "SRR9673911", "AB_all", "AC_all", "CA_all", "CL_all", "CM_all", "CN_all", "CO_all", "CV_all", "DC_all", "DK_all", "FR_all", "GM_all", "LF_all", "LT_all", "MA_all", "MC_all", "ML_all", "MP_all", "MS_all", "MZ_all", "NO_all", "NP_all", "OA_all", "PC_all", "TC_all", "TF_all", "TI_all", "TP_all"),
                       groups = as.factor(c("pit", "pit",      "pit",         "castle",      "pit",      "castle",     "pit",        "pit",        "pit",        "pit",        "pit",        "pit",        "pit",        "pit",        "pit",        "pit",    "pit",     "none",  "castle", "castle", "castle", "none",   'pit',     "pit",    "pit",   "pit",    "none",   'none',   "none",   "castle", "castle", "pit",    "none",   "pit",    "none",   "castle", "pit",    "castle", "none",   "castle",  "castle", "pit",  "pit")),
                       origin = as.factor(c(rep("tang", 15), rep("mlwi", 28))))
png("~/scratch/brain/fst/tang/pca_max_miss_90.png", width = 500, height = 500)
fviz_pca_ind(res.pca,
             col.ind = sample_df$groups, # color by groups
             palette = c("#00AFBB", "gray", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = FALSE
             # max.overlaps = Inf
)
dev.off()

# high_fst2 = high_fst[1:8]
# for (i in 9:51) {
#   print(i)
#   my_alleles = substr(high_fst[,i], 0, 3)
#   my_split = colsplit(my_alleles, pattern = "\\/", names = c(paste0(colnames(high_fst)[i], ".a"), paste0(colnames(high_fst)[i], ".b")))
#   high_fst2 = cbind(high_fst2, my_split)
# }

#====================================================================================================
# Common SNPs method ================================================================================
#====================================================================================================
# Pre-process, clean, and select high FST Malawi Pit vs Castle Sites
fst = read.table("~/scratch/brain/fst/tang/fst_pc_malawi_only_on_common_w_vcf.fst", header = F)
fst_header = scan("~/scratch/brain/fst/tang/malawi_tang_merged_header.vcf", what = character())
fst_header[1] = "CHROM"
vcf_header_junk = c("CHROM", "POS", "ID", "REF", "ALT_tang", "QUAL", "FILTER", "INFO", "FORMAT")
fst_header = c(c("CHROM_FST", "POS_FST", "FST"), fst_header[1:24], vcf_header_junk, fst_header[25:52]) 
# fst_header = c(c("CHROM_FST", "POS_FST", "FST"), fst_header)
fst_header[13:27] = substr(fst_header[13:27], 8, 1000L)
colnames(fst) = fst_header
#Clean
fst = fst[complete.cases(fst$FST),]
fst[which(fst$FST < 0),3] = 0
thresh95_clean = quantile(fst$FST, 0.95)

fst = fst[which(fst$ALT_tang == fst$ALT),] # NEW
fst = fst[,c("CHROM_FST", "POS_FST", "FST", fst_header[13:27], fst_header[37:length(fst_header)])]

high_fst = fst[which(fst$FST > thresh95_clean),]
write.table(fst, "~/scratch/brain/fst/tang/fst_pc_malawi_only_on_common_w_vcf_clean.fst")
write.table(high_fst, "~/scratch/brain/fst/tang/fst_pc_malawi_only_on_common_w_vcf_clean_high.fst")

biallelic = read.table("~/scratch/brain/fst/tang/malawi_tang_common_snp_biallelic.vcf", sep="\t", stringsAsFactors = F, header = F)
high_fst = high_fst[which(paste0(high_fst$CHROM_FST, "_", high_fst$POS_FST) %in% paste0(biallelic[,1], "_", biallelic[,2])),]
write.table(high_fst, "~/scratch/brain/fst/tang/fst_pc_malawi_only_on_common_w_vcf_clean_high.fst")

# Start Here
high_fst = read.table("~/scratch/brain/fst/tang/fst_pc_malawi_only_on_common_w_vcf_clean_high.fst")
backup_high_fst = high_fst
high_fst = backup_high_fst

# PCA Method from: https://comppopgenworkshop2019.readthedocs.io/en/latest/contents/03_pca/pca.html
high_fst = high_fst[which(high_fst$FST >= 0.50),]
scoreSNPs = function(fst_thresh = NULL) {
  if (is.null(fst_thresh))
    high_fst = high_fst
  else
    high_fst = high_fst[which(high_fst$FST >= fst_thresh),]
  
  high_fst3 = high_fst[1:3]
  for (i in 4:ncol(high_fst)) {
    allele1 = substr(high_fst[,i], 0, 1)
    allele2 = substr(high_fst[,i], 3, 3)
    
    isRef1 = allele1 == 0
    isHomo = allele1 == allele2
    isMissing1 = allele1 == "."
    
    score = rep(-1, nrow(high_fst)) # no score should be -1, if any is at the end, something went wrong
    score[which(isMissing1)] = 9
    score[which(!isMissing1 & isHomo & isRef1)] = 2
    score[which(!isMissing1 & isHomo & !isRef1)] = 0
    score[which(!isMissing1 & !isHomo )] = 1
    
    high_fst3 = cbind(high_fst3, score)
  }
  colnames(high_fst3) = colnames(high_fst)
  return(high_fst3)
}


res.pca <- prcomp(t(high_fst3[4:ncol(high_fst3)]))
# high_fst3 = rbind(high_fst3, as.vector(sample_df$groups[match(colnames(high_fst3), sample_df$samples)]))

library(factoextra)
sample_df = data.frame(samples = c("SRR9657485", "SRR9657511", "SRR9657512", "SRR9657530", "SRR9657540", "SRR9657558", "SRR9657559", "SRR9665654", "SRR9665657", "SRR9665678", "SRR9665679", "SRR9665706", "SRR9673822", "SRR9673823", "SRR9673911", "AB_all", "AC_all", "CA_all", "CL_all", "CM_all", "CN_all", "CO_all", "CV_all", "DC_all", "DK_all", "FR_all", "GM_all", "LF_all", "LT_all", "MA_all", "MC_all", "ML_all", "MP_all", "MS_all", "MZ_all", "NO_all", "NP_all", "OA_all", "PC_all", "TC_all", "TF_all", "TI_all", "TP_all"),
                       groups = as.factor(c("pit", "pit",      "pit",         "castle",      "pit",      "castle",     "pit",        "pit",        "pit",        "pit",        "pit",        "pit",        "pit",        "pit",        "pit",        "pit",    "pit",     "none",  "castle", "castle", "castle", "none",   'pit',     "pit",    "pit",   "pit",    "none",   'none',   "none",   "castle", "castle", "pit",    "none",   "pit",    "none",   "castle", "pit",    "castle", "none",   "castle",  "castle", "pit",  "pit")),
                       origin = as.factor(c(rep("tang", 15), rep("mlwi", 28))))
png("~/scratch/brain/fst/tang/pca_common_90.png", width = 500, height = 500)
fviz_pca_ind(res.pca,
             col.ind = sample_df$groups, # color by groups
             palette = c("#00AFBB", "gray", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = FALSE
             # max.overlaps = Inf
)
dev.off()

threshs = seq(from=0.55, to=0.9, by=0.05)
for (thresh in threshs) {
  high_fst3 = scoreSNPs(thresh)
  res.pca <- prcomp(t(high_fst3[4:ncol(high_fst3)]))
  thresh_str = substr(as.character(thresh), 3, 1000L)
  png(paste0("~/scratch/brain/fst/tang/pca_common_", thresh_str, ".png"), width = 500, height = 500)
  print(fviz_pca_ind(res.pca,
               col.ind = sample_df$groups, # color by groups
               palette = c("#00AFBB", "gray", "#FC4E07"),
               addEllipses = TRUE, # Concentration ellipses
               ellipse.type = "confidence",
               legend.title = "Groups",
               repel = FALSE
               # max.overlaps = Inf
  ))
  dev.off()
}

# Find Correlation Between Malawi FST and Tang FST
tang_fst = read.table("~/scratch/brain/fst/tang/fst_pc_tang_unbin.weir.fst", sep = "\t", stringsAsFactors = F, header = T)
colnames(tang_fst) = c("CHROM", "POS", "FST")
high_fst_paste = paste0(high_fst$CHROM_FST, "_", high_fst$POS_FST)
tang_fst_paste = paste0(tang_fst$CHROM, "_", tang_fst$POS)
tang_fst_common = tang_fst[match(high_fst_paste, tang_fst_paste),]
tang_fst_common[which(tang_fst_common$FST < 0),3] = 0
tang_fst_common$FST[which(! complete.cases(tang_fst_common$FST))] = 0
write.csv(tang_fst_common, "~/scratch/brain/fst/tang/fst_pc_tang_unbin_common.weir.fst")

# Plot Correlations
library("ggplot2")
library("viridis")
high_fst3 = high_fst
high_fst3$tang = tang_fst_common$FST
high_fst3 = high_fst3[which(high_fst3$FST >= 0.9),]
colnames(high_fst3)[3] = "malawi"
cor(high_fst3$malawi, high_fst3$tang)
png("~/scratch/brain/fst/tang/tang_common_cor_90.png", width = 400, height = 400)
ggplot(high_fst3, aes(x=malawi, y = tang, color = malawi)) + geom_point() + scale_color_viridis_c()
dev.off()

# Predict Tanginkya based on Malawi
library("parallel")
# high_fst3 = high_fst[which(high_fst$FST >= 0.5),]
numCores = detectCores()
predictPit = function(x) {
  this_row = high_fst3[x,]
  this_row = rbind(this_row, as.vector(sample_df$groups[match(colnames(this_row), sample_df$samples)]))
  this_row = rbind(this_row, as.vector(sample_df$origin[match(colnames(this_row), sample_df$samples)]))
  mlwi_pit_mean = mean( as.numeric(this_row[1,which(this_row[1,] != 9 & this_row[2,] == "pit" & this_row[3,] == "mlwi")]) )
  mlwi_tang_mean = mean( as.numeric(this_row[1,which(this_row[1,] != 9 & this_row[2,] == "castle" & this_row[3,] == "mlwi")]) )
  # print(mlwi_pit_mean)
  # print(mlwi_tang_mean)
  
  # Pit Score
  this_row[4,] = NA
  this_row[4,which(this_row[1,] != 9)] = this_row[1,which(this_row[1,] != 9)]
  this_row[4,] = as.numeric(this_row[4,])
  if (mlwi_pit_mean < mlwi_tang_mean)
    this_row[4,] = 2 - as.numeric(this_row[4,])
  this_row[4,] = as.numeric(this_row[4,]) / 2
  
  # Castle Score
  this_row[5,] = NA 
  this_row[5,which(this_row[1,] != 9)] = this_row[1,which(this_row[1,] != 9)]
  this_row[5,] = as.numeric(this_row[5,])
  if (mlwi_pit_mean >= mlwi_tang_mean)
    this_row[5,] = 2 - as.numeric(this_row[5,])
  this_row[5,] = as.numeric(this_row[5,]) / 2
  
  # tang_pit_score = this_row[1, which(this_row[2,] == "pit" & this_row[3,] == "tang")]
  # tang_castle_score = this_row[1, which(this_row[2,] == "castle" & this_row[3,] == "tang")]
  # if (mlwi_pit_mean < 1)
  #   tang_pit_score = 2 - tang_pit_score
  # else
  #   tang_castle_score = 2 - tang_castle_score
  
  return(list(this_row[4,4:ncol(this_row)], this_row[5,4:ncol(this_row)]))
}

high_fst3_backup = high_fst3
predictions = mclapply(1:nrow(high_fst3), predictPit, mc.cores = numCores)

big_scores = data.frame()
for (thresh in c(0.135, 0.5, 0.6, 0.7, 0.8, 0.9)) {
  high_fst3 = scoreSNPs(thresh)
  sample_df2 = sample_df[match(colnames(high_fst3), sample_df$samples),]
  high_fst3[high_fst3 == 9] = NA
  mlwi_pit_means = rowMeans(high_fst3[,which(as.vector(sample_df2$groups) == "pit" & as.vector(sample_df2$origin) == "mlwi")], na.rm = T)
  mlwi_castle_means = rowMeans(high_fst3[,which(as.vector(sample_df2$groups) == "castle" & as.vector(sample_df2$origin) == "mlwi")], na.rm = T)
  pit_smaller = mlwi_pit_means < mlwi_castle_means
  high_fst_num = high_fst3
  high_fst_num[high_fst_num >= 0] = 1
  scores_df = data.frame()
  for (this_group in c("pit", "castle")) {
    for (this_origin in c("tang", "mlwi")) {
      pit2 = colSums(high_fst3[which(!pit_smaller),which(as.vector(sample_df2$groups) == this_group & as.vector(sample_df2$origin) == this_origin)], na.rm = T)
      pit0 = colSums(high_fst3[which( pit_smaller),which(as.vector(sample_df2$groups) == this_group & as.vector(sample_df2$origin) == this_origin)], na.rm = T)
      
      pit2_num = colSums(high_fst_num[which(!pit_smaller),which(as.vector(sample_df2$groups) == this_group & as.vector(sample_df2$origin) == this_origin)], na.rm = T)
      pit0_num = colSums(high_fst_num[which( pit_smaller),which(as.vector(sample_df2$groups) == this_group & as.vector(sample_df2$origin) == this_origin)], na.rm = T)
      
      pit_all = ((pit2 + (2*pit0_num - pit0)) / (pit2_num + pit0_num)) / 2
      pit2 = (pit2 / pit2_num) / 2
      pit0 = 1-((pit0 / pit0_num) / 2)
      
      castle2 = colSums(high_fst3[which(!castle_smaller),which(as.vector(sample_df2$groups) == this_group & as.vector(sample_df2$origin) == this_origin)], na.rm = T)
      castle0 = colSums(high_fst3[which( castle_smaller),which(as.vector(sample_df2$groups) == this_group & as.vector(sample_df2$origin) == this_origin)], na.rm = T)
      
      castle2_num = colSums(high_fst_num[which(!castle_smaller),which(as.vector(sample_df2$groups) == this_group & as.vector(sample_df2$origin) == this_origin)], na.rm = T)
      castle0_num = colSums(high_fst_num[which( castle_smaller),which(as.vector(sample_df2$groups) == this_group & as.vector(sample_df2$origin) == this_origin)], na.rm = T)
      
      castle_all = ((castle2 + (2*castle0_num - castle0)) / (castle2_num + castle0_num)) / 2
      castle2 = (castle2 / castle2_num) / 2
      castle0 = 1-((castle0 / castle0_num) / 2)
      
      scores_df = rbind(scores_df, t(c(this_origin, this_group, mean(pit2), mean(pit0), mean(pit_all), mean(castle2), mean(castle0), mean(castle_all))))
      # print(paste0(this_origin, ", ", this_group, ": ", mean(pit2), ", ", mean(pit0), ", ", mean(pit_all)))
    }
  }
  colnames(scores_df) = c("lake", "bower", "pit_ref", "pit_alt", "pit_score", "castle_ref", "castle_alt", "castle_score")
  scores_df$thresh = thresh
  scores_df[,c("lake", "bower", "thresh", "pit_score", "castle_score")]
  big_scores = rbind(big_scores, scores_df)
}
print(big_scores[,c("lake", "bower", "thresh", "pit_score", "castle_score")])

big_scores = big_scores[which(big_scores$thresh != 0.90),]
big_scores$combo = paste0(big_scores$lake, "_", big_scores$bower)
big_scores$pit_score = as.numeric(big_scores$pit_score)*100
big_scores$castle_score = as.numeric(big_scores$castle_score)*100

png("~/scratch/brain/fst/tang/pit_scores.png", width = 700, height = 500, res = 90)
ggplot(big_scores, aes(x=combo, y = pit_score, group = thresh, fill = thresh)) + geom_bar(stat = "identity", position=position_dodge2()) + ggtitle("Pit Scores at various FST Thresholds") + xlab("Lake and Bower Type") + ylab("Pit Score") + geom_text(aes(label = round(pit_score, digits = 1)), position = position_dodge(0.9), vjust = -0.2)
dev.off()

png("~/scratch/brain/fst/tang/castle_scores.png", width = 700, height = 500, res = 90)
ggplot(big_scores, aes(x=combo, y = castle_score, group = thresh, fill = thresh)) + geom_bar(stat = "identity", position=position_dodge2()) + ggtitle("Castle Scores at various FST Thresholds") + xlab("Lake and Bower Type") + ylab("Castle Score") + geom_text(aes(label = round(castle_score, digits = 1)), position = position_dodge(0.9), vjust = -0.2)
dev.off()

# pit_scores = colMeans(high_fst3[which(!pit_smaller),which(as.vector(sample_df2$groups) == "pit" & as.vector(sample_df2$origin) == "tang")], na.rm = T)
# pit_scores2 =  (2-colMeans(high_fst3[which(pit_smaller),which(as.vector(sample_df2$groups) == "pit" & as.vector(sample_df2$origin) == "tang")], na.rm = T))/2
# castle_smaller = mlwi_castle_means < mlwi_pit_means
# castle_scores = colMeans(high_fst3[which(!castle_smaller),which(as.vector(sample_df2$groups) == "pit" & as.vector(sample_df2$origin) == "tang")], na.rm = T)


pitis2 = mlwi_pit_means > 1
colMeans(high_fst3[which(pitis2),which(as.vector(sample_df2$groups) == "pit" & as.vector(sample_df2$origin) == "mlwi")], na.rm = T)
castleis2 = mlwi_castle_means > 1
mean( colMeans(high_fst3[which(castleis2),which(as.vector(sample_df2$groups) == "castle" & as.vector(sample_df2$origin) == "mlwi")], na.rm = T) ) / 2

  
pit_lists = list()
castle_lists = list()
for (i in 1:nrow(high_fst3)) {
  res = predictPit(i)
  pit_lists[[i]] = res[[1]]
  castle_lists[[i]] = res[[2]]
}

#*******************************************************************************
# Linkage Disequilibrium =======================================================
#*******************************************************************************
library(dplyr)
library(plyr)
library(data.table)
ld = data.table::fread("~/scratch/brain/fst/sample_vcf/pit_ld_1kb.geno.ld")
ld$BIN_START = round_any(ld$POS1, 10000, f = ceiling)
ld$R2 = ld[,5]
ld$BIN_ID = paste0(ld$CHROM, "_", ld$BIN_START)
test = ld %>% group_by(BIN_START) %>% summarise(across(R2, mean, na.rm = TRUE))
test = ld[,mean(R2), by = list(BIN_ID)]

# ld_complete = ld[which( ! is.na(ld[,5]) ),]
# ld_complete$BIN_START = round_any(ld_complete$POS1, 10000, f = ceiling)
# ld_complete$BIN_START1 = round_any(ld_complete$POS1, 10000, f = ceiling)
# ld_complete$BIN_START2 = round_any(ld_complete$POS2, 10000, f = ceiling)
# ld_complete$BIN_ID = paste0(ld_complete$CHROM, "_", ld_complete$BIN_START)

#*******************************************************************************
# Merge all ====================================================================
#*******************************************************************************
rc = read.delim("rock_castle_10kb.windowed.weir.fst")
rc$BIN_ID = paste0(rc$CHROM, "_", rc$BIN_START)
rc$RANK = order(order(rc$WEIGHTED_FST, decreasing = T))
rc$QUANTILE = 1 - (rc$RANK / nrow(rc))
rc$WEIGHTED_rc[which(rc$WEIGHTED_rc < 0)] = 0
rc$ZFST = (( rc$WEIGHTED_FST - mean(rc$WEIGHTED_FST) ) / sd(rc$WEIGHTED_FST)) + 1
rc$RANKZ = order(order(rc$ZFST, decreasing = T))
rc$QUANTILEZ = 1 - (rc$RANKZ / nrow(rc))
pc = read.delim("pit_castle_10kb.windowed.weir.fst")
pc$BIN_ID = paste0(pc$CHROM, "_", pc$BIN_START)
pc$RANK = order(order(pc$WEIGHTED_FST, decreasing = T))
pc$QUANTILE = 1 - (pc$RANK / nrow(pc))
pc$ZFST = (( pc$WEIGHTED_FST - mean(pc$WEIGHTED_FST) ) / sd(pc$WEIGHTED_FST)) + 1
pc$RANKZ = order(order(pc$ZFST, decreasing = T))
pc$QUANTILEZ = 1 - (pc$RANKZ / nrow(pc))

pcrc_merged = merge(pc, rc, by = "BIN_ID", suffixes = c("_PC", "_RC"))

tang_pc = read.delim("tang_pvc_10kb.windowed.weir.fst")
tang_pc$RANK = order(order(tang_pc$WEIGHTED_FST, decreasing = T))
tang_pc$QUANTILE = 1 - (tang_pc$RANK / nrow(tang_pc))
tang_pc$ZFST = (( tang_pc$WEIGHTED_FST - mean(tang_pc$WEIGHTED_FST) ) / sd(tang_pc$WEIGHTED_FST)) + 1
tang_pc$RANKZ = order(order(tang_pc$ZFST, decreasing = T))
tang_pc$QUANTILEZ = 1 - (tang_pc$RANKZ / nrow(tang_pc))
colnames(tang_pc) = paste0(colnames(tang_pc), "_TANGPC")
tang_pc$BIN_ID = paste0(tang_pc$CHROM_TANGPC, "_", tang_pc$BIN_START_TANGPC)

all_merged = merge(pcrc_merged, tang_pc, by = "BIN_ID")
all_thresh = all_merged[which(all_merged$WEIGHTED_FST_PC >= 0.2 & all_merged$WEIGHTED_FST_RC >= 0.2 & all_merged$WEIGHTED_FST_TANGPC >= 0.2),]
all_thresh_quantile = all_merged[which(all_merged$QUANTILE_PC >= 0.95 & all_merged$QUANTILE_RC >= 0.95 & all_merged$QUANTILE_TANGPC >= 0.95),]
all_thresh_quantilez = all_merged[which(all_merged$QUANTILEZ_PC >= 0.95 & all_merged$QUANTILEZ_RC >= 0.95 & all_merged$QUANTILEZ_TANGPC >= 0.95),]
all_merged$highlight = all_merged$BIN_ID %in% all_thresh_quantilez$BIN_ID
head(all_thresh_quantilez[order(all_thresh_quantilez$QUANTILEZ_TANGPC, decreasing = T),c("BIN_ID", "ZFST_PC", "ZFST_RC", "ZFST_TANGPC")], 10)

write.table(all_thresh_quantilez[,c("CHROM_PC", "BIN_START_PC", "BIN_END_PC")], "~/scratch/brain/fst/all_3.bed", sep = "\t", row.names = F, col.names = F, quote = F)
all_thresh_closest = read.delim("~/scratch/brain/fst/all_3_closest.bed", header = F)
all_thresh_closest$BIN_ID = paste0(all_thresh_closest[,1], "_", all_thresh_closest[,2])
all_thresh_closest$GENE = colsplit(all_thresh_closest[,12], "; gene ", c("1", "2"))[,2]
all_thresh_closest$GENE = colsplit(all_thresh_closest$GENE, "; ", c("1", "2"))[,1]
all_thresh_closest[,13] = abs(all_thresh_closest[,13])
all_thresh_quantilez$DIST = all_thresh_closest[match(all_thresh_quantilez$BIN_ID, all_thresh_closest$BIN_ID), 13]
all_thresh_quantilez$GENE = all_thresh_closest[match(all_thresh_quantilez$BIN_ID, all_thresh_closest$BIN_ID), "GENE"]
all_thresh_quantilez = all_thresh_quantilez[order(all_thresh_quantilez$QUANTILEZ_TANGPC, decreasing = T),]
write.csv(all_thresh_quantilez, "~/scratch/brain/fst/all_3_quantilez_df.csv")
write.csv(all_thresh_quantilez$GENE[which(all_thresh_quantilez$DIST < 25000 & all_thresh_quantilez$GENE != "")], "~/scratch/brain/fst/all_3_quantilez_genes.csv")
write.csv(all_thresh_quantilez$GENE[which(all_thresh_quantilez$DIST < 25000 & all_thresh_quantilez$CHROM_PC == "NC_036790.1" & all_thresh_quantilez$BIN_START_PC > 5950001 & all_thresh_quantilez$BIN_END_PC < 25280000)], "~/scratch/brain/fst/all_3_quantilez_genes_lg11_peak.csv")

p1 = ggplot(all_merged, aes(x = ZFST_PC, y = ZFST_RC,      color = highlight, alpha = highlight)) + geom_point() + scale_color_manual(values = c("gray50", "red")) + scale_alpha_manual(values = c(0.2, 1)) + NoLegend()
p2 = ggplot(all_merged, aes(x = ZFST_PC, y = ZFST_TANGPC, color = highlight, alpha = highlight)) + geom_point() + scale_color_manual(values = c("gray50", "red")) + scale_alpha_manual(values = c(0.2, 1)) + NoLegend()
p3 = ggplot(all_merged, aes(x = ZFST_RC, y = ZFST_TANGPC, color = highlight, alpha = highlight)) + geom_point() + scale_color_manual(values = c("gray50", "red")) + scale_alpha_manual(values = c(0.2, 1)) + NoLegend()

png("~/scratch/brain/fst/all_zfst.png", width = 1200, height = 400, res=100)
print(cowplot::plot_grid(plotlist = list(p1, p2, p3), ncol = 3))
dev.off()

fst = tang_pc
fst$id = 1:nrow(fst)
fst$CHROM <- sub("^NW.*", "unplaced", fst$CHROM)
fst$CHROM <- sub("^NC_027944.1", "unplaced", fst$CHROM)
fst$CHROM <- gsub(".1$", "", fst$CHROM)
fst$highlight = fst$BIN_ID %in% all_thresh_quantilez$BIN_ID
test = sample(1:nrow(fst), 5000)
my_breaks = which(! duplicated(fst$CHROM))
fst = fst[order(fst$highlight),]
image_name <- "~/scratch/brain/fst/tang_pc_highlight_051822.png"
png(image_name, type="cairo", width = 12, height = 4, units = 'in', res=300)
p = ggplot(fst, aes(id, ZFST_TANGPC, color = highlight)) + geom_point(alpha = 0.7, size = 1) + theme_classic() + scale_color_manual(values = rep(brewer.pal(n=8,name="Dark2"), 6)) + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(breaks = my_breaks, labels = unique(fst$CHROM), expand=c(0,0)) + NoLegend() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ggtitle("Pit vs Castle in Tanganyka with Highlight")
print(p)
dev.off()

#*******************************************************************************
# PI ===========================================================================
#*******************************************************************************

# Pit vs Castle
pit_pi = read.delim("~/scratch/brain/fst/sample_vcf/pit_pi_10kb.windowed.pi")
pit_pi$PI[which(pit_pi$PI < 0)] = 0
pit_pi$Z <- (( pit_pi$PI - mean(pit_pi$PI) ) / sd(pit_pi$PI)) + 1
pit_pi$BIN_ID = paste0(pit_pi$CHROM, "_", pit_pi$BIN_START)

castle_pi = read.delim("~/scratch/brain/fst/sample_vcf/castle_pi_10kb.windowed.pi")
castle_pi$PI[which(castle_pi$PI < 0)] = 0
castle_pi$Z <- (( castle_pi$PI - mean(castle_pi$PI) ) / sd(castle_pi$PI)) + 1
castle_pi$BIN_ID = paste0(castle_pi$CHROM, "_", castle_pi$BIN_START)

pvc_pi = merge(pit_pi, castle_pi, by = "BIN_ID", suffixes = c("_PIT", "_CASTLE"))
pvc_pi$Z = pvc_pi$Z_PIT - pvc_pi$Z_CASTLE
pvc_pi$CHROM = pvc_pi$CHROM_PIT
pvc_pi$lg = lgConverter(pvc_pi$CHROM, path_to_info = "~/scratch/m_zebra_ref/M_zebra_UMD2a_assembly_report.txt")
pvc_pi = pvc_pi[which(pvc_pi$lg != "MT"),]
pvc_pi$PI_RANKZ = order(order(pvc_pi$Z, decreasing = T))
pvc_pi$PI_QUANTILEZ = 1 - (pvc_pi$PI_RANKZ / nrow(pvc_pi))

fst = pvc_pi
fst$id = 1:nrow(fst)
fst$CHROM <- sub("^NW.*", "unplaced", fst$CHROM)
fst$CHROM <- sub("^NC_027944.1", "unplaced", fst$CHROM)
fst$CHROM <- gsub(".1$", "", fst$CHROM)
my_breaks = which(! duplicated(fst$lg))
image_name <- "~/scratch/brain/fst/pvc_zpi.png"
png(image_name, type="cairo", width = 12, height = 4, units = 'in', res=300)
p = ggplot(fst, aes(id, Z, color = lg)) + geom_point(alpha = 0.7, size = 1) + theme_classic() + scale_color_manual(values = rep(brewer.pal(n=8,name="Dark2"), 6)) + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(breaks = my_breaks, labels = unique(fst$lg), expand=c(0,0)) + NoLegend() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ggtitle("Pit - Castle ZPI") + ylab("ZPI")
print(p)
dev.off()

pvc_pi_fst = merge(pvc_pi, pc, by = "BIN_ID")
pvc_pi_fst$CHROM = pvc_pi_fst$CHROM_PIT
pvc_pi_fst$CHROM <- sub("^NW.*", "unplaced", pvc_pi_fst$CHROM)
pvc_pi_fst$CHROM <- sub("^NC_027944.1", "unplaced", pvc_pi_fst$CHROM)
pvc_pi_fst$CHROM <- gsub(".1$", "", pvc_pi_fst$CHROM)
pvc_pi_fst$lg = paste0(pvc_pi_fst$lg, " (", pvc_pi_fst$CHROM, ")")
pvc_pi_fst$highlight = "none"
pvc_pi_fst$highlight[which(pvc_pi_fst$QUANTILEZ > 0.95)] = "ZFST +95%"
pvc_pi_fst$highlight[which(pvc_pi_fst$PI_QUANTILEZ > 0.95)] = "ZPI +95%"
pvc_pi_fst$highlight[which(pvc_pi_fst$PI_QUANTILEZ > 0.95 & pvc_pi_fst$QUANTILEZ > 0.95)] = "Both +95%"
pvc_pi_fst$highlight = factor(pvc_pi_fst$highlight, levels = c("none", "ZFST +95%", "ZPI +95%", "Both +95%"))
pvc_pi_fst$lg = factor(pvc_pi_fst$lg, levels = unique(pvc_pi_fst$lg))
pvc_pi_fst = pvc_pi_fst[order(pvc_pi_fst$highlight),]
png("~/scratch/brain/fst/pvc_zpi_fst.png", width = 2000, height = 1600, res = 200)
ggplot(pvc_pi_fst, aes(x = ZFST, y = Z, color = highlight)) + geom_point(alpha = 0.5, size = 0.9) + scale_color_manual(values = c("#ffcdb2", "#e5989b", "#b5838d", "#6d6875")) + ylab("ZPI") + ggtitle("Pit v Castle (ZFST and ZPI)") + facet_wrap(~ lg, ncol = 5) + theme_classic() + guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()

fst = pc
fst$id = 1:nrow(fst)
fst$WEIGHTED_FST[which(fst$WEIGHTED_FST < 0)] = 0
fst$Zfst <- (( fst$WEIGHTED_FST - mean(fst$WEIGHTED_FST) ) / sd(fst$WEIGHTED_FST)) + 1
fst$CHROM <- sub("^NW.*", "unplaced", fst$CHROM)
fst$CHROM <- sub("^NC_027944.1", "unplaced", fst$CHROM)
fst$CHROM <- gsub(".1$", "", fst$CHROM)
my_breaks = which(! duplicated(fst$CHROM))
fst$highlight = pvc_pi_fst$highlight[match(fst$BIN_ID, pvc_pi_fst$BIN_ID)]
fst$highlight[which( is.na(fst$highlight) )] = "none"
fst = fst[order(fst$highlight),]
png("~/scratch/brain/fst/pvc_zpi_fst_fst.png", type="cairo", width = 12, height = 4, units = 'in', res=300)
p = ggplot(fst, aes(id, Zfst, color = highlight)) + geom_point(alpha = 0.7, size = 1) + theme_classic() + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(breaks = my_breaks, labels = unique(fst$CHROM), expand=c(0,0)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ggtitle("Pit v Castle ZFST w/ Highlight") + scale_color_manual(values = c("#ffcdb2", "#e5989b", "#b5838d", "#6d6875")) + guides(colour = guide_legend(override.aes = list(size=3)))
print(p)
dev.off()

# Rock vs Castle
rock_pi = read.delim("~/scratch/brain/fst/sample_vcf/rock_pi_10kb.windowed.pi")
rock_pi$PI[which(rock_pi$PI < 0)] = 0
rock_pi$Z <- (( rock_pi$PI - mean(rock_pi$PI) ) / sd(rock_pi$PI)) + 1
rock_pi$BIN_ID = paste0(rock_pi$CHROM, "_", rock_pi$BIN_START)

rvc_pi = merge(rock_pi, castle_pi, by = "BIN_ID", suffixes = c("_ROCK", "_CASTLE"))
rvc_pi$Z = rvc_pi$Z_ROCK - rvc_pi$Z_CASTLE
rvc_pi$CHROM = rvc_pi$CHROM_ROCK
rvc_pi$lg = lgConverter(rvc_pi$CHROM, path_to_info = "~/scratch/m_zebra_ref/M_zebra_UMD2a_assembly_report.txt")
rvc_pi = rvc_pi[which(rvc_pi$lg != "MT"),]
rvc_pi$PI_RANKZ = order(order(rvc_pi$Z, decreasing = T))
rvc_pi$PI_QUANTILEZ = 1 - (rvc_pi$PI_RANKZ / nrow(rvc_pi))

fst = rvc_pi
fst$id = 1:nrow(fst)
fst$CHROM <- sub("^NW.*", "unplaced", fst$CHROM)
fst$CHROM <- sub("^NC_027944.1", "unplaced", fst$CHROM)
fst$CHROM <- gsub(".1$", "", fst$CHROM)
my_breaks = which(! duplicated(fst$lg))
image_name <- "~/scratch/brain/fst/rvc_zpi.png"
png(image_name, type="cairo", width = 12, height = 4, units = 'in', res=300)
p = ggplot(fst, aes(id, Z, color = lg)) + geom_point(alpha = 0.7, size = 1) + theme_classic() + scale_color_manual(values = rep(brewer.pal(n=8,name="Dark2"), 6)) + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(breaks = my_breaks, labels = unique(fst$lg), expand=c(0,0)) + NoLegend() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ggtitle("Rock - Castle ZPI") + ylab("ZPI")
print(p)
dev.off()

fst = pc
fst$id = 1:nrow(fst)
fst$WEIGHTED_FST[which(fst$WEIGHTED_FST < 0)] = 0
fst$Zfst <- (( fst$WEIGHTED_FST - mean(fst$WEIGHTED_FST) ) / sd(fst$WEIGHTED_FST)) + 1
fst$CHROM <- sub("^NW.*", "unplaced", fst$CHROM)
fst$CHROM <- sub("^NC_027944.1", "unplaced", fst$CHROM)
fst$CHROM <- gsub(".1$", "", fst$CHROM)
my_breaks = which(! duplicated(fst$CHROM))
fst$RVC_ZFST_QUANTILEZ = rc$QUANTILEZ[match(fst$BIN_ID, rc$BIN_ID)]
fst$RVC_PI_QUANTILEZ = rvc_pi$Z[match(fst$BIN_ID, rvc_pi$BIN_ID)]
fst$PVC_PI_QUANTILEZ = pvc_pi$Z[match(fst$BIN_ID, pvc_pi$BIN_ID)]
fst$highlight = "none"
fst$highlight[which(fst$RVC_ZFST_QUANTILEZ > 0.95 & fst$RVC_PI_QUANTILEZ > 0.95 & fst$QUANTILEZ > 0.95 & fst$PVC_PI_QUANTILEZ > 0.95)] = "+95% ZFST&ZPI in PVC and RVC"
fst$highlight = factor(fst$highlight, levels = c("none", "+95% ZFST&ZPI in PVC and RVC"))
fst = fst[order(fst$highlight),]
png("~/scratch/brain/fst/pvc_rvc_zfst_zpi.png", type="cairo", width = 12, height = 4, units = 'in', res=300)
p = ggplot(fst, aes(id, Zfst, color = highlight)) + geom_point(alpha = 0.7, size = 1) + theme_classic() + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(breaks = my_breaks, labels = unique(fst$CHROM), expand=c(0,0)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ggtitle("Pit v Castle ZFST w/ Highlight") + scale_color_manual(values = c("#ffcdb2", "#6d6875")) + guides(colour = guide_legend(override.aes = list(size=3)))
print(p)
dev.off()


#*******************************************************************************
# Tajima's D ===================================================================
#*******************************************************************************

# Pit vs Castle
pit_d = read.delim("~/scratch/brain/fst/sample_vcf/pit_d_10kb.Tajima.D")
colnames(pit_d)[4] = "D"
pit_d$D[which(is.na(pit_d$D))] = 0
pit_d$Z <- (( pit_d$D - mean(pit_d$D) ) / sd(pit_d$D)) + 1
pit_d$BIN_ID = paste0(pit_d$CHROM, "_", pit_d$BIN_START)

castle_d = read.delim("~/scratch/brain/fst/sample_vcf/castle_d_10kb.Tajima.D")
colnames(castle_d)[4] = "D"
castle_d$D[which(is.na(castle_d$D))] = 0
castle_d$Z <- (( castle_d$D - mean(castle_d$D) ) / sd(castle_d$D)) + 1
castle_d$BIN_ID = paste0(castle_d$CHROM, "_", castle_d$BIN_START)

pvc_d = merge(pit_d, castle_d, by = "BIN_ID", suffixes = c("_PIT", "_CASTLE"))
pvc_d$Z = pvc_d$Z_PIT - pvc_d$Z_CASTLE
pvc_d$CHROM = pvc_d$CHROM_PIT
pvc_d$lg = lgConverter(pvc_d$CHROM, path_to_info = "~/scratch/m_zebra_ref/M_zebra_UMD2a_assembly_report.txt")
pvc_d = pvc_d[which(pvc_d$lg != "MT"),]
pvc_d$D_RANKZ = order(order(pvc_d$Z, decreasing = T))
pvc_d$D_QUANTILEZ = 1 - (pvc_d$D_RANKZ / nrow(pvc_d))
pvc_d$DSS_RANKZ = order(order(pvc_d$Z, decreasing = F)) # SS == Selective Sweep
pvc_d$DSS_QUANTILEZ = 1 - (pvc_d$D_RANKZ / nrow(pvc_d))

fst = pvc_d
fst$id = 1:nrow(fst)
fst$CHROM <- sub("^NW.*", "unplaced", fst$CHROM)
fst$CHROM <- sub("^NC_027944.1", "unplaced", fst$CHROM)
fst$CHROM <- gsub(".1$", "", fst$CHROM)
my_breaks = which(! duplicated(fst$lg))
image_name <- "~/scratch/brain/fst/pvc_zd.png"
png(image_name, type="cairo", width = 12, height = 4, units = 'in', res=300)
p = ggplot(fst, aes(id, Z, color = lg)) + geom_point(alpha = 0.7, size = 1) + theme_classic() + scale_color_manual(values = rep(brewer.pal(n=8,name="Dark2"), 6)) + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(breaks = my_breaks, labels = unique(fst$lg), expand=c(0,0)) + NoLegend() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ggtitle("Pit - Castle ZD") + ylab("ZD")
print(p)
dev.off()

# Rock vs Castle
rock_d = read.delim("~/scratch/brain/fst/sample_vcf/rock_d_10kb.Tajima.D")
colnames(rock_d)[4] = "D"
rock_d$D[which(is.na(rock_d$D))] = 0
rock_d$Z <- (( rock_d$D - mean(rock_d$D) ) / sd(rock_d$D)) + 1
rock_d$BIN_ID = paste0(rock_d$CHROM, "_", rock_d$BIN_START)

rvc_d = merge(rock_d, castle_d, by = "BIN_ID", suffixes = c("_ROCK", "_CASTLE"))
rvc_d$Z = rvc_d$Z_ROCK - rvc_d$Z_CASTLE
rvc_d$CHROM = rvc_d$CHROM_ROCK
rvc_d$lg = lgConverter(rvc_d$CHROM, path_to_info = "~/scratch/m_zebra_ref/M_zebra_UMD2a_assembly_report.txt")
rvc_d = rvc_d[which(rvc_d$lg != "MT"),]
rvc_d$D_RANKZ = order(order(rvc_d$Z, decreasing = T))
rvc_d$D_QUANTILEZ = 1 - (rvc_d$D_RANKZ / nrow(rvc_d))
rvc_d$DSS_RANKZ = order(order(rvc_d$Z, decreasing = F)) # SS == Selective Sweep
rvc_d$DSS_QUANTILEZ = 1 - (rvc_d$DSS_RANKZ / nrow(rvc_d))

fst = rvc_d
fst$id = 1:nrow(fst)
fst$CHROM <- sub("^NW.*", "unplaced", fst$CHROM)
fst$CHROM <- sub("^NC_027944.1", "unplaced", fst$CHROM)
fst$CHROM <- gsub(".1$", "", fst$CHROM)
my_breaks = which(! duplicated(fst$lg))
image_name <- "~/scratch/brain/fst/rvc_zd.png"
png(image_name, type="cairo", width = 12, height = 4, units = 'in', res=300)
p = ggplot(fst, aes(id, Z, color = lg)) + geom_point(alpha = 0.7, size = 1) + theme_classic() + scale_color_manual(values = rep(brewer.pal(n=8,name="Dark2"), 6)) + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(breaks = my_breaks, labels = unique(fst$lg), expand=c(0,0)) + NoLegend() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab("") + ggtitle("Rock - Castle ZD") + ylab("ZD")
print(p)
dev.off()

fst = pc
fst$id = 1:nrow(fst)
fst$WEIGHTED_FST[which(fst$WEIGHTED_FST < 0)] = 0
fst$Zfst <- (( fst$WEIGHTED_FST - mean(fst$WEIGHTED_FST) ) / sd(fst$WEIGHTED_FST)) + 1
fst$CHROM <- sub("^NW.*", "unplaced", fst$CHROM)
fst$CHROM <- sub("^NC_027944.1", "unplaced", fst$CHROM)
fst$CHROM <- gsub(".1$", "", fst$CHROM)
my_breaks = which(! duplicated(fst$CHROM))
fst$RVC_ZFST_QUANTILEZ = rc$QUANTILEZ[match(fst$BIN_ID, rc$BIN_ID)]
fst$RVC_PI_QUANTILEZ = rvc_pi$PI_QUANTILEZ[match(fst$BIN_ID, rvc_pi$BIN_ID)]
fst$PVC_PI_QUANTILEZ = pvc_pi$PI_QUANTILEZ[match(fst$BIN_ID, pvc_pi$BIN_ID)]
fst$RVC_D_QUANTILEZ = rvc_d$D_QUANTILEZ[match(fst$BIN_ID, rvc_pi$BIN_ID)]
fst$PVC_D_QUANTILEZ = pvc_d$D_QUANTILEZ[match(fst$BIN_ID, pvc_pi$BIN_ID)]
print(paste0("Number of Bins that are 95th Quantile for ZFST & ZPI in PvC and RvC: ", length(which(fst$RVC_ZFST_QUANTILEZ > 0.95 & fst$RVC_PI_QUANTILEZ > 0.95 & fst$QUANTILEZ > 0.95 & fst$PVC_PI_QUANTILEZ > 0.95))))
print(paste0("Number of Bins that are 95th Quantile for ZFST & ZD in PvC and RvC: ", length(which(fst$RVC_ZFST_QUANTILEZ > 0.95 & fst$RVC_D_QUANTILEZ > 0.95 & fst$QUANTILEZ > 0.95 & fst$PVC_D_QUANTILEZ > 0.95))))
print(paste0("Number of Bins that are 95th Quantile for ZPI & ZD in PvC and RvC: ", length(which(fst$RVC_PI_QUANTILEZ > 0.95 & fst$RVC_D_QUANTILEZ > 0.95 & fst$PVC_PI_QUANTILEZ > 0.95 & fst$PVC_D_QUANTILEZ > 0.95))))
print(paste0("Number of Bins that are 95th Quantile for ZFST & ZPI & ZD in PvC and RvC: ", length(which(fst$RVC_ZFST_QUANTILEZ > 0.95 & fst$RVC_D_QUANTILEZ > 0.95 & fst$RVC_PI_QUANTILEZ > 0.95 & fst$QUANTILEZ > 0.95 & fst$PVC_D_QUANTILEZ > 0.95 & fst$PVC_PI_QUANTILEZ))))
fst$highlight = "none"
fst$highlight[which(fst$RVC_ZFST_QUANTILEZ > 0.95 & fst$RVC_PI_QUANTILEZ > 0.95 & fst$QUANTILEZ > 0.95 & fst$PVC_PI_QUANTILEZ > 0.95)] = "+95% ZFST&ZPI in PVC and RVC"
fst$highlight = factor(fst$highlight, levels = c("none", "+95% ZFST&ZPI in PVC and RVC"))

