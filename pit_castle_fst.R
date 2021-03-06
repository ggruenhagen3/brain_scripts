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
