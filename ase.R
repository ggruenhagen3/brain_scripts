library("RUnit")
library("MBASED")
library("metap")
#=========================================================================================
# New UMD2a Data
#=========================================================================================
rna_path = "C:/Users/miles/Downloads/brain/"
SRR904 = read.table(paste(rna_path, "/data/ase/SRR5440904_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)
SRR905 = read.table(paste(rna_path, "/data/ase/SRR5440905_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)
SRR906 = read.table(paste(rna_path, "/data/ase/SRR5440906_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)
SRR907 = read.table(paste(rna_path, "/data/ase/SRR5440907_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)
SRR908 = read.table(paste(rna_path, "/data/ase/SRR5440908_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)
SRR909 = read.table(paste(rna_path, "/data/ase/SRR5440909_counts.tsv", sep=""), header = TRUE, stringsAsFactors = F)

SRR904$GENE = str_replace(SRR904$GENE,"%", " (1 of many)")
SRR905$GENE = str_replace(SRR905$GENE,"%", " (1 of many)")
SRR906$GENE = str_replace(SRR906$GENE,"%", " (1 of many)")
SRR907$GENE = str_replace(SRR907$GENE,"%", " (1 of many)")
SRR908$GENE = str_replace(SRR908$GENE,"%", " (1 of many)")
SRR909$GENE = str_replace(SRR909$GENE,"%", " (1 of many)")

gene_names = SRR904$GENE
# gene_names = str_replace(gene_names,"%", " (1 of many)")
n_boot = 100

# Prepare the Data
pit_mc = SRR904$MC_COUNTS + SRR904$MC_COUNTS
pit_cv = SRR905$CV_COUNTS + SRR905$CV_COUNTS
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

# Do 1-sampled ASE experiments
pos_all_ind = which(SRR904$MC_COUNTS + SRR904$CV_COUNTS > 0 &
                    SRR905$MC_COUNTS + SRR905$CV_COUNTS > 0 &
                    SRR906$MC_COUNTS + SRR906$CV_COUNTS > 0 &
                    SRR907$MC_COUNTS + SRR907$CV_COUNTS > 0 &
                    SRR908$MC_COUNTS + SRR908$CV_COUNTS > 0 &
                    SRR909$MC_COUNTS + SRR909$CV_COUNTS > 0)
# SRR904_1 = my_MBASED_1(SRR904$MC_COUNTS, SRR904$CV_COUNTS, "SRR904 (Pit 1)", pos_all_ind, gene_names, n_boot)
# SRR905_1 = my_MBASED_1(SRR905$MC_COUNTS, SRR905$CV_COUNTS, "SRR905 (Pit 2)", pos_all_ind, gene_names, n_boot)
# SRR906_1 = my_MBASED_1(SRR906$MC_COUNTS, SRR906$CV_COUNTS, "SRR906 (Castle 1)", pos_all_ind, gene_names, n_boot)
# SRR907_1 = my_MBASED_1(SRR907$MC_COUNTS, SRR907$CV_COUNTS, "SRR907 (Castle 2)", pos_all_ind, gene_names, n_boot)
# SRR908_1 = my_MBASED_1(SRR908$MC_COUNTS, SRR908$CV_COUNTS, "SRR908 (Isolated 1)", pos_all_ind, gene_names, n_boot)
# SRR909_1 = my_MBASED_1(SRR909$MC_COUNTS, SRR909$CV_COUNTS, "SRR909 (Isolated 2)", pos_all_ind, gene_names, n_boot)

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

write.table(all_sig_same_dir, "C:/Users/miles/Downloads/brain/results/ase_all_sig_same_dir.txt", quote = F, col.names = F, row.names = F)
write.table(all_same_dir, "C:/Users/miles/Downloads/brain/results/ase_all_same_dir.txt", quote = F, col.names = F, row.names = F)
write.table(agg, "C:/Users/miles/Downloads/brain/results/ase_agg_sig_same_dir.txt", quote = F, col.names = F, row.names = F)

all_sig_same_dir_hgnc = hgncMzebraInPlace(data.frame(all_sig_same_dir), 1, gene_names)
all_same_dir_hgnc     = hgncMzebraInPlace(data.frame(all_same_dir),     1, gene_names)
agg_hgnc              = hgncMzebraInPlace(data.frame(agg),              1, gene_names)

write.table(all_sig_same_dir_hgnc, "C:/Users/miles/Downloads/brain/results/ase_all_sig_same_dir_hgnc.txt", quote = F, col.names = F, row.names = F)
write.table(all_same_dir_hgnc,     "C:/Users/miles/Downloads/brain/results/ase_all_same_dir_hgnc.txt", quote = F, col.names = F, row.names = F)
write.table(agg_hgnc,              "C:/Users/miles/Downloads/brain/results/ase_agg_sig_same_dir_hgnc.txt", quote = F, col.names = F, row.names = F)


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
  lociAllele1Counts
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

my_MBASED = function(s1_mc, s1_cv, s2_mc, s2_cv, s1_name, s2_name, gene_names, n_boot, myIsPhased=T, verbose=T) {
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
  pos_ind = which( s1_mc + s1_cv > 0 & s2_mc + s2_cv > 0 )
  pos_gene = gene_names[pos_ind]
  this_s1_mc = s1_mc[pos_ind]
  this_s1_cv = s1_cv[pos_ind]
  this_s2_mc = s2_mc[pos_ind]
  this_s2_cv = s2_cv[pos_ind]
  if (verbose) {
    print(paste("Genes Used", length(pos_gene)))
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
ggplot(boot_res, aes(overlap_in_pvc_and_cvp, alpha=.7, fill=above, color = "maroon")) + geom_histogram(alpha=0.5) + geom_vline(aes(xintercept = real_ovlp_pc_v_cp)) + geom_text(aes(x=real_ovlp_pc_v_cp, label="Real Value"), y = Inf, hjust=0, vjust=1, color = "black") + xlab("# of Gene in Overlap Between Pit v Castle and Castle v Pit") + ggtitle("Comparison Between Bootstrap Values and Real Value") + guides(color=F, alpha=F, fill=F)

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
