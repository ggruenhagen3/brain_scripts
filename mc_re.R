args = commandArgs(trailingOnly=TRUE)
fpath = args[1]

cdf = read.delim(fpath, header = F)
cdf$dist = cdf[,ncol(cdf)]
cdf$loc = reshape2::colsplit(cdf[,12], ";", c("1", "2"))[,1]
cdf$loc = reshape2::colsplit(cdf$loc, " ", c("1", "2"))[,2]
cdf$gene_name = reshape2::colsplit(cdf[,12], "; gene_source", c("1", "2"))[,1]
cdf$gene_name = reshape2::colsplit(cdf$gene_name, "gene_name ", c("1", "2"))[,2]
cdf = cdf[which(cdf[,12] != "."),]
cdf$gene = cdf$gene_name = cdf$loc
cdf2 = cdf[which( abs(cdf$dist) < 25000 ),]
cdf2$class = "distal"
cdf2$class[which(cdf2$dist <= 5000 & cdf2$dist > 0)] = "promoter"
cdf2$class[which(cdf2$dist == 0)] = "intragenic"
cdf3 = cdf2[, c("gene", "class")]
write.csv(cdf3, paste0(fpath, ".class.csv"))
