# Mouse Atlas
library(loomR)
library("scrattch.io")
library("Seurat")
library("Matrix")
library("stringr")
library("dplyr")
library("cowplot")
library("ggplot2")
library("biomaRt")
options(stringsAsFactors = FALSE)

####################
# Helper Functions #
####################
convertToHgncObj <- function(obj, organism) {
  # Convert a Seurat object to have all HGNC gene names
  print("Converting Genes Names...")
  genes <- rownames(obj)
  if (organism == "mouse") {
    dataset_name <- "mmusculus_gene_ensembl"
  } else if (organism == "mzebra") {
    dataset_name <- "mzebra_gene_ensembl"
  } else {
    stop("Organism not recognized. Options are mouse or mzebra.")
  }
  human  = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  org = useMart(biomart="ensembl", dataset=dataset_name)
  
  # DF to convert from org to HGNC
  all_hgnc <- getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = genes , mart = org, attributesL = c("external_gene_name"), martL = human, uniqueRows=T)
  
  # Initialize New Matricies
  print("Creating New Matrices...")
  new_counts_matrix <- as.matrix(obj@assays$RNA@counts)
  new_data_matrix   <- as.matrix(obj@assays$RNA@data)
  
  not_i <- 0
  multiple_hgnc <- 0
  bad_multiple_hgnc <- 0
  for (gene in genes) {
    if (gene %in% all_hgnc[,1]) {
      hgnc_gene <- all_hgnc[which(all_hgnc[,1] == gene),2]
      if (length(hgnc_gene) > 1) {
        multiple_hgnc <- multiple_hgnc + 1
        upper_hgnc_gene <- hgnc_gene[which(startsWith(tolower(hgnc_gene), gene))]
        if (length(upper_hgnc_gene) == 1) {
          hgnc_gene <- upper_hgnc_gene
        } else {
          bad_multiple_hgnc <- bad_multiple_hgnc + 1
          hgnc_gene <- hgnc_gene[1]
        } # end bad multiple
      } # end multiple
      
      rownames(new_counts_matrix)[which(rownames(new_counts_matrix) == gene)] <- hgnc_gene
      rownames(new_data_matrix)[which(rownames(new_data_matrix) == gene)] <- hgnc_gene
    } else {
      not_i <- not_i + 1
    } # end gene not an hgnc gene
  } # end gene for
  print(paste("Number of Genes not converted to HGNC:", not_i))
  print(paste("Number of Genes with multple HGNC:", multiple_hgnc))
  print(paste("Number of Genes with multple HGNC and non-ideal circumstance:", bad_multiple_hgnc))
  
  # Merge the duplicated rows
  print("Removing duplicated HGNC rows...")
  dup_genes <- rownames(new_counts_matrix)[which(duplicated(rownames(new_counts_matrix)))]
  dup_ind <- c()
  for (gene in dup_genes) {
    ind <- which(rownames(new_counts_matrix) == gene)
    
    new_counts_row <- new_counts_matrix[ind[1],]
    new_data_row   <- new_data_matrix[ind[1],]
    for (i in 2:length(ind)) {
      new_counts_row <- new_counts_row + new_counts_matrix[ind[i],]
      new_data_row   <- new_data_row   + new_data_matrix[ind[i],]
    }
    new_counts_matrix[ind[1],] <- new_counts_row
    new_data_matrix[ind[1],]   <- new_data_row
    
    # Delete the duplicated rows
    dup_ind <- c(dup_ind, ind[2:length(ind)])
  }
  # Delete all the duplicated rows at once
  new_counts_matrix <- new_counts_matrix[-dup_ind,]
  new_data_matrix   <- new_data_matrix[-dup_ind,]
  
  # Remove mzebra rows
  print("Removing old org rows...")
  ind <- which(rownames(new_counts_matrix) %in% all_hgnc[,2])
  new_counts_matrix <- new_counts_matrix[ind,]
  new_data_matrix   <- new_data_matrix[ind,]
  
  print("Creating New Seurat Object...")
  obj_2 <- CreateSeuratObject(counts = new_counts_matrix, project = obj@project.name)
  obj_2 <- SetAssayData(object = obj_2, slot = 'data', new.data = new_data_matrix)
  obj_2$seurat_clusters <- obj$seurat_clusters
  
  # Add the metadata
  for (col in colnames(obj@meta.data)) {
    obj_2@meta.data[col] <- obj@meta.data[col]
  }
  
  return(obj_2)
}

###############
# Main Script #
###############
rna_path  <- "/nv/hp10/ggruenhagen3/scratch/brain/brain_scripts/"
tome_path <- "/nv/hp10/ggruenhagen3/scratch/brain/data/m_cor_hip/"

tome <- paste(tome_path, "transcrip.tome", sep="")
exons       <- read_tome_dgCMatrix(tome,"data/t_exon")
introns     <- read_tome_dgCMatrix(tome,"data/t_intron")
sample_name <- read_tome_sample_names(tome)
gene_name   <- read_tome_gene_names(tome)

# Plot original TSNE
tsne_coordinates <- read.csv(paste(tome_path, "2d_coordinates.csv", sep=""), header = TRUE)
metadata <- read.csv(paste(tome_path, "sample_annotations.csv", sep=""), header = TRUE)
head(tsne_coordinates)
head(metadata)
original_df <- inner_join(tsne_coordinates, metadata, by = "sample_name")
head(original_df)
# original_df <- cbind(tsne_coordinates[2:3], metadata$cluster_color)
# colnames(original_df) <- c("tsne_1", "tsne_2", "cluster_color")
# png(paste(tome_path, "original_tsne.png", sep=""), width = 4800, height = 3600, unit = "px")
# ggplot(original_df, aes(tsne_1, tsne_2, color=cluster_color)) + geom_point() + theme_classic()
# dev.off()

# for (col in colnames(metadata)[which(! colnames(metadata) %in% c("sample_name", "cell_type_alt_alias_label"))]) {
#   print(col)
#   filename <- paste(tome_path, "/results/original_tsne_", col, ".png", sep="")
#   png(filename, width = 5000, height = 4000, unit = "px", res = 150)
#   p <- ggplot(original_df, aes_(as.name("tsne_1"), as.name("tsne_2"), color=as.name(col))) + geom_point()  + theme(axis.line=element_blank(),
#                                                                                                                                               axis.text.x=element_blank(),
#                                                                                                                                               axis.text.y=element_blank(),
#                                                                                                                                               axis.ticks=element_blank(),
#                                                                                                                                               axis.title.x=element_blank(),
#                                                                                                                                               axis.title.y=element_blank(),
#                                                                                                                                               legend.key = element_rect(fill = "white"),
#                                                                                                                                               panel.background=element_blank(),
#                                                                                                                                               plot.background=element_blank()) + guides(col=guide_legend(ncol=4))
#   print(p)
#   dev.off()
#   # system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/Brain/m_cor_hip/results/", sep=""))
# }
# system(paste("rclone copy ", tome_path, "/results/original_tsne_* dropbox:BioSci-Streelman/George/Brain/m_cor_hip/results/"))

# Redo Clustering in Seurat
# TODO: Look at cell_type_alias_label
colnames(exons) <- sample_name
colnames(introns) <- sample_name
rownames(exons) <- gene_name
rownames(introns) <- gene_name
m_cor_hip <- exons + introns
m_cor_hip <- CreateSeuratObject(counts = m_cor_hip, project = "M_COR_HIP")
m_cor_hip <- NormalizeData(m_cor_hip, normalization.method = "LogNormalize", scale.factor = 100000)
# cat(paste("Number of rows in metadata:", nrow(metadata), "\n"))
# cat(paste("Number of rows in m_cor_hip (no subset):", nrow(metadata), "\n"))

# Add all the metadata
metadata <- metadata %>% slice(match(colnames(m_cor_hip), sample_name))
num_col_obj_metadata <- ncol(m_cor_hip@meta.data)
for (i in 1:ncol(metadata)) {
  m_cor_hip@meta.data[,num_col_obj_metadata + i] <- metadata[,i]
  colnames(m_cor_hip@meta.data)[num_col_obj_metadata + i] <- colnames(metadata)[i]
}
m_cor_hip <- subset(m_cor_hip, subset = nFeature_RNA > 500)

# Exclude cells labelled "Exclude"
Idents(m_cor_hip) <- "class_label"
m_cor_hip <- subset(m_cor_hip, idents = levels(Idents(m_cor_hip))[which(levels(Idents(m_cor_hip)) != "Exclude")])

# Standard Pipeline
m_cor_hip <- FindVariableFeatures(object = m_cor_hip, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
m_cor_hip <- ScaleData(object = m_cor_hip, vars.to.regress = NULL)
m_cor_hip <- RunPCA(m_cor_hip, npcs = 50, verbose = FALSE)
m_cor_hip <- RunTSNE(object = m_cor_hip)
m_cor_hip <- RunUMAP(m_cor_hip, reduction = "pca", dims = 1:20)
m_cor_hip <- FindNeighbors(m_cor_hip, reduction = "umap", dims = 1:2)
m_cor_hip <- FindClusters(m_cor_hip, resolution = 0.1)

m_cor_hip_hgnc <- convertToHgncObj(m_cor_hip, "mouse")
saveRDS(m_cor_hip_hgnc, paste(tome_path, "/m_cor_hip_hgnc.rds", sep=""))

## metadata <- metadata %>% slice(match(colnames(m_cor_hip), sample_name))
## m_cor_hip$class_label <- metadata$class_label[which(metadata$sample_name %in% colnames(m_cor_hip))]

# filename <- paste(tome_path, "new_umap_2.png", sep="")
# png(filename, width = 3600, height = 3600, unit = "px", res = 300)
# Idents(m_cor_hip) <- "seurat_clusters"
# DimPlot(m_cor_hip, reduction = "umap", label = TRUE) + NoLegend()
# dev.off()
# system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/Brain/m_cor_hip/results/", sep=""))
# 
# Idents(m_cor_hip) <- "cell_type_alias_label"
# n = length(levels(Idents(m_cor_hip)))
# hues = seq(15, 375, length = n + 1)
# cols = hcl(h = hues, l = 65, c = 100)[1:n]
# for (i in 1:n) {
#   cell_type <- levels(Idents(m_cor_hip))[i]
#   cell_type_str <- gsub(" ", "_", cell_type)
#   cell_type_str <- gsub("/", "-", cell_type_str)
#   filename <- paste(tome_path, "/results/new_umap/new_umap_", cell_type_str, ".png", sep="")
#   print(filename)
#   png(filename, width = 3600, height = 3600, unit = "px", res = 300)
#   Idents(m_cor_hip) <- "cell_type_alias_label"
#   p <- DimPlot(m_cor_hip[, WhichCells(m_cor_hip, idents = c(cell_type))], reduction = "umap", label = TRUE, cols = cols[i])  + xlim(min(m_cor_hip@reductions$umap@cell.embeddings[,1]), max(m_cor_hip@reductions$umap@cell.embeddings[,1])) + ylim(min(m_cor_hip@reductions$umap@cell.embeddings[,2]), max(m_cor_hip@reductions$umap@cell.embeddings[,2]))
#   print(p)
#   dev.off()
#   system(paste("rclone copy ", filename, " dropbox:BioSci-Streelman/George/Brain/m_cor_hip/results/new_umap/", sep=""))
# }
# 
# png(paste(tome_path, "old_umap.png", sep=""), width = 3600, height = 3600, unit = "px", res = 300)
# Idents(m_cor_hip) <- "class_order"
# DimPlot(m_cor_hip, reduction = "umap", label = TRUE) + NoLegend()
# dev.off()
# 
# png(paste(tome_path, "new_tsne.png", sep=""), width = 3600, height = 3600, unit = "px", res = 300)
# Idents(m_cor_hip) <- "seurat_clusters"
# DimPlot(m_cor_hip, reduction = "tsne", label = TRUE) + NoLegend()
# dev.off()
# 
# png(paste(tome_path, "new_tsne_class.png", sep=""), width = 3600, height = 3600, unit = "px", res = 300)
# Idents(m_cor_hip) <- "class_label"
# DimPlot(m_cor_hip, reduction = "tsne", label = TRUE)
# dev.off()
# 
# png(paste(tome_path, "old_tsne.png", sep=""), width = 3600, height = 3600, unit = "px", res = 300)
# Idents(m_cor_hip) <- "class_order"
# DimPlot(m_cor_hip, reduction = "tsne", label = TRUE) + NoLegend()
# dev.off()
# 
# test_df <- as.data.frame(m_cor_hip@reductions$tsne@cell.embeddings[,1:2])
# test_df$cluster_color <- m_cor_hip$cluster_color
# png(paste(tome_path, "old_tsne_clrs.png", sep=""), width = 3600, height = 3600, unit = "px", res = 300)
# ggplot(test_df, aes(tSNE_1, tSNE_2, color=cluster_color)) + geom_point() + theme_classic() + NoLegend()
# dev.off()
