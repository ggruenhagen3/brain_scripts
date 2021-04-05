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

high_fst = read.table("~/scratch/brain/fst/tang/fst_pc_malawi_only_on_common_w_vcf_clean_high.fst")
backup_high_fst = high_fst

# PCA Method from: https://comppopgenworkshop2019.readthedocs.io/en/latest/contents/03_pca/pca.html
high_fst = high_fst[which(high_fst$FST >= 0.9999),]
high_fst3 = high_fst[1:3]
for (i in 4:ncol(high_fst)) {
  allele1 = substr(high_fst[,i], 0, 1)
  allele2 = substr(high_fst[,i], 3, 3)
  
  isRef1 = allele1 == 0
  isHomo = allele1 == allele2
  isMissing1 = allele1 == "."
  
  score = rep(-1, nrow(high_fst))
  score[which(isMissing1)] = 9
  score[which(isHomo & isRef1)] = 2
  score[which(isHomo & !isRef1)] = 0
  score[which(! isHomo & ! isMissing1 )] = 1
  
  high_fst3 = cbind(high_fst3, score)
}
colnames(high_fst3) = colnames(high_fst)
res.pca <- prcomp(t(high_fst3[4:ncol(high_fst3)]))

library(factoextra)
sample_df = data.frame(samples = c("SRR9657485", "SRR9657511", "SRR9657512", "SRR9657530", "SRR9657540", "SRR9657558", "SRR9657559", "SRR9665654", "SRR9665657", "SRR9665678", "SRR9665679", "SRR9665706", "SRR9673822", "SRR9673823", "SRR9673911", "AB_all", "AC_all", "CA_all", "CL_all", "CM_all", "CN_all", "CO_all", "CV_all", "DC_all", "DK_all", "FR_all", "GM_all", "LF_all", "LT_all", "MA_all", "MC_all", "ML_all", "MP_all", "MS_all", "MZ_all", "NO_all", "NP_all", "OA_all", "PC_all", "TC_all", "TF_all", "TI_all", "TP_all"),
                       groups = as.factor(c("pit", "pit",      "pit",         "castle",      "pit",      "castle",     "pit",        "pit",        "pit",        "pit",        "pit",        "pit",        "pit",        "pit",        "pit",        "pit",    "pit",     "none",  "castle", "castle", "castle", "none",   'pit',     "pit",    "pit",   "pit",    "none",   'none',   "none",   "castle", "castle", "pit",    "none",   "pit",    "none",   "castle", "pit",    "castle", "none",   "castle",  "castle", "pit",  "pit")),
                       origin = as.factor(c(rep("tang", 15), rep("mlwi", 28))))
png("~/scratch/brain/fst/tang/pca_common_100.png", width = 500, height = 500)
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
