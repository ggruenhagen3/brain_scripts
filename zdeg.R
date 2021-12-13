setwd("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/results/deg_lmer_demux/53_clusters")
library(glmmSeq)
library(Seurat)
library(qvalue)
library(tidyverse)
library(tidyr)
library(lme4)
library(lmerTest)
library(edgeR)
library(dplyr)
library(DESeq2)
library(parallel)
library(genefilter)
library(glmGamPoi)
library(scran)

combined <- readRDS("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb_demux_102021.rds")

clusters_3plus_per_subject=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,36)

#DEG analysis
#k=11
for (k in clusters_3plus_per_subject[1]){
  #for (k in 1:1){
  print(paste0("starting cluster ",k))
  #  if(k %in% skip_clusters){
  #    next
  #  }
  data <- subset(x = combined, subset = seuratclusters53 == k)
  
  df_counts <- data.frame(t(as.matrix(data@assays$RNA@counts)))
  head(df_counts[1:5,1:5])
  
  n_genes=ncol(df_counts)
  n_cells=nrow(df_counts)
  
  df_counts_meta = data.frame(rownames(df_counts))
  df_counts_meta$id <- df_counts_meta$rownames.df_counts.
  df_counts_meta$rownames.df_counts. = NULL
  df_counts_meta$trial_id = data$trial_id  
  df_counts_meta$bower_activity_index = as.numeric(data$bower_activity_index)
  df_counts_meta$gsi = as.numeric(data$gsi) 
  df_counts_meta$log_spawn_events = as.numeric(data$log_spawn_events)  
  df_counts_meta$run= as.factor(data$run)
  df_counts_meta$sample = as.factor(data$sample)
  df_counts_meta$pair= as.factor(data$pair)
  
  # cut out genes that were not detected in any cells from dataframe 
  df_counts_no_0 <- df_counts[,which(colSums(df_counts) != 0)]
  
  # store list of detected genes
  n_genes_no_0=ncol(df_counts_no_0)
  
  # add metadata back into the dataframe
  df_counts_no_0 <- cbind(df_counts_no_0,df_counts_meta)
  
  # split based on pair
  df_counts_no_0_split_by_pair <- split(df_counts_no_0, f = df_counts_no_0$pair) #split by pair
  
  #remove genes with 0 counts within each pair  
  for(l in 1:length(df_counts_no_0_split_by_pair)){
    print(paste0("finding good genes for pair ",l))
    temp_pair_l <- data.frame(df_counts_no_0_split_by_pair[[l]])
    temp_pair_l_counts <- temp_pair_l[,1:n_genes_no_0]
    temp_pair_l_counts_no_0 <- temp_pair_l_counts[,which(colSums(temp_pair_l_counts) != 0)]
    #temp <- cbind(good_genes, temp[20407:20410])
    out=data.frame(colnames(temp_pair_l_counts_no_0))
    assign(x=paste0("gene_list_pair_",l),value=get("out"))
  }
  
  # generate list of genes that were detected in all pairs
  good_gene_list <- gene_list_pair_1$colnames.temp_pair_l_counts_no_0. 
  for(m in 2:length(df_counts_no_0_split_by_pair)){
    print(paste0("integrating good genes from pair ",m," of 19"))
    temp_good_gene_list_m <- data.frame(value=get(paste0("gene_list_pair_",m)))
    temp_good_gene_list_m <- temp_good_gene_list_m$colnames.temp_pair_l_counts_no_0.
    good_gene_list <- intersect(good_gene_list,temp_good_gene_list_m)
  }
  
  p <- length(good_gene_list)
  
  # create a count dataframe for genes that were detected in all pairs
  df_counts_no_0_all_pairs <- df_counts_no_0[, good_gene_list]
  head(df_counts_no_0_all_pairs[1:5,1:5])
  
  # create a count matrix for genes that were detected in all pairs  
  count_matrix_final <- as.matrix(df_counts_no_0_all_pairs)
  count_matrix_final <- as.data.frame(t(count_matrix_final))
  
  # add metadata into the dataframe for genes that were detected in all pairs   
  df_counts_no_0_all_pairs <- cbind(df_counts_no_0_all_pairs,df_counts_meta)
  
  #set up coldata for DESeq2
  
  coldata <- df_counts_meta
  coldata$trial_id <- as.factor(coldata$trial_id)
  coldata$bower_activity_index <- as.numeric(coldata$bower_activity_index)
  coldata$log_spawn_events <- as.numeric(coldata$log_spawn_events)
  coldata$gsi <- as.numeric(coldata$gsi)
  coldata$sample <- as.factor(coldata$sample)
  coldata$run <- as.factor(coldata$run)
  coldata$pair <- as.factor(coldata$pair)
  
  cluster_size <- ncol(ncol(count_matrix_final))
  
  dds <- DESeqDataSetFromMatrix(countData = count_matrix_final,
                                colData = coldata,
                                design = ~ run + gsi + log_spawn_events + bower_activity_index)
  
  size_factors <- calculateSumFactors(
    count_matrix_final,
    #sizes = seq(27, 101, 5),
    clusters = NULL,
    ref.clust = NULL,
    max.cluster.size = cluster_size,
    positive = TRUE,
    scaling = NULL,
    min.mean = NULL,
    subset.row = NULL,
    #BPPARAM = SerialParam()
  )                              
  
  sizeFactors(dds) <- size_factors
  
  dds <- estimateDispersions(
    dds,
    fitType = "glmGamPoi",
    #method = "runed",
    #sharingMode = 'maximum',
    useCR = TRUE,
    maxit = 100,
    weightThreshold = 0.01,
    quiet = FALSE,
    modelMatrix = NULL,
    minmu = 1e-06
  )
  
  disp <- as.matrix(mcols(dds))
  disp <- disp[,11]
  names(disp) <- good_gene_list
  
  glmmseq_counts <- data.frame(t(as.matrix(df_counts_no_0_all_pairs)))
  glmmseq_counts <- glmmseq_counts[1:p,]
  glmmseq_counts <- as.data.frame(sapply(glmmseq_counts, as.numeric))
  rownames(glmmseq_counts) <- colnames(df_counts_no_0_all_pairs[,1:p])
  
  results <- glmmSeq(~ bower_activity_index + log_spawn_events + gsi + (1|run/sample/trial_id) + (1|pair/trial_id),
                     #glmmSeq(~ log_spawn_events + bower_activity_index + (1|pair:trial_id) + (1|run:trial_id),
                     id = "trial_id",
                     countdata = glmmseq_counts,
                     metadata = coldata,
                     dispersion = disp,
                     removeDuplicatedMeasures = FALSE,
                     removeSingles=FALSE,
                     progress=TRUE,
                     cores = 3
  )
  
  out_final = data.frame(results@stats)
  hist(out_final$P_bower_activity_index)
  out_final$cluster=k-1
  saveRDS(glmmseq_counts, "glmmseq_counts.rds")
  saveRDS(coldata, "coldata.rds")
  saveRDS(dds, "dds.rds")
  saveRDS(disp, "disp.rds")
  # write.csv(out_final, paste0("bb53_deg_glmmseq_demux_bower_activity_index_log_spawn_events_gsi_control_subjectinsampleinrun_subjectinpair_random_good_genes_by_pair_cluster",k-1,"_102121.csv"))
} 