# Mouse Atlas
library("loomR")
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
keepCommonGenesObj <- function(obj_a, obj_b) {
  # Finds the common genes between two Seurat objects and
  # makes two new Seurat objects with just those genes.
  genes_a <- rownames(obj_a@assays$RNA@counts)
  genes_b <- rownames(obj_b@assays$RNA@counts)
  common <- genes_a[which(genes_a %in% genes_b)]
  print(paste("Found", length(common), "genes in common -", length(common)/length(genes_a) * 100, "of obj_a genes", length(common)/length(genes_b) * 100, "of obj_b genes" ))
  
  ## First Seurat Object
  # Initialize New Matricies
  print("Creating New Matrices (for First Seurat Object)...")
  new_counts_matrix <- as(obj_a@assays$RNA@counts, "sparseMatrix") 
  new_data_matrix   <- as(obj_a@assays$RNA@data, "sparseMatrix")
  
  # Removing Non-Overlapping Genes
  print("Removing Non-Overlapping Genes (for the First Seurat Object)...")
  all_ind_keep <- c()
  all_ind <- 1:length(new_counts_matrix)
  for (gene in common) {
    ind_keep <- which(rownames(new_counts_matrix) == gene)
    all_ind_keep <- c(all_ind_keep, ind_keep)
  }
  all_ind_remove <- all_ind[which(! all_ind %in% all_ind_keep)]
  new_counts_matrix <- new_counts_matrix[-all_ind_remove,]
  new_data_matrix   <- new_data_matrix[-all_ind_remove,]
  
  print("Creating New Seurat Object (for the First Seurat Object)...")
  obj_a_2 <- CreateSeuratObject(counts = new_counts_matrix, project = obj_a@project.name)
  obj_a_2 <- SetAssayData(object = obj_a_2, slot = 'data', new.data = new_data_matrix)
  obj_a_2$seurat_clusters <- obj_a$seurat_clusters
  
  # Add the metadata
  for (col in colnames(obj_a@meta.data)) {
    obj_a_2@meta.data[col] <- obj_a@meta.data[col]
  }
  
  ## Second Seurat Object
  # Initialize New Matricies
  print("Creating New Matrices (for Second Seurat Object)...")
  new_counts_matrix <- as(obj_b@assays$RNA@counts, "sparseMatrix") 
  new_data_matrix   <- as(obj_b@assays$RNA@data, "sparseMatrix")
  
  # Removing Non-Overlapping Genes
  print("Removing Non-Overlapping Genes (for the Second Seurat Object)...")
  all_ind_keep <- c()
  all_ind <- 1:length(new_counts_matrix)
  for (gene in common) {
    ind_keep <- which(rownames(new_counts_matrix) == gene)
    all_ind_keep <- c(all_ind_keep, ind_keep)
  }
  all_ind_remove <- all_ind[which(! all_ind %in% all_ind_keep)]
  new_counts_matrix <- new_counts_matrix[-all_ind_remove,]
  new_data_matrix   <- new_data_matrix[-all_ind_remove,]
  
  print("Creating New Seurat Object (for the Second Seurat Object)...")
  obj_b_2 <- CreateSeuratObject(counts = new_counts_matrix, project = obj_b@project.name)
  obj_b_2 <- SetAssayData(object = obj_b_2, slot = 'data', new.data = new_data_matrix)
  obj_b_2$seurat_clusters <- obj_b$seurat_clusters
  
  # Add the metadata
  for (col in colnames(obj_b@meta.data)) {
    obj_b_2@meta.data[col] <- obj_b@meta.data[col]
  }
  
  return(list(obj_a_2, obj_b_2))
}

convertToHgncObj <- function(obj, organism) {
  # Convert a Seurat object to have all HGNC gene names
  print("Converting Genes Names...")
  genes <- rownames(obj@assays$RNA@counts)
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
  new_counts_matrix <- as(obj@assays$RNA@counts, "sparseMatrix") 
  new_data_matrix   <- as(obj@assays$RNA@data, "sparseMatrix") 
  
  not_i <- 0
  multiple_hgnc <- 0
  bad_multiple_hgnc <- 0
  for (gene in genes) {
    if (gene %in% all_hgnc[,1]) {
      hgnc_gene <- all_hgnc[which(all_hgnc[,1] == gene),2]
      if (length(hgnc_gene) > 1) {
        multiple_hgnc <- multiple_hgnc + 1
        upper_hgnc_gene <- hgnc_gene[which(startsWith(tolower(hgnc_gene), tolower(gene)))]
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
  print(paste("-Number of Genes not converted to HGNC:", not_i))
  print(paste("-Number of Genes with multple HGNC:", multiple_hgnc))
  print(paste("-Number of Genes with multple HGNC and non-ideal circumstance:", bad_multiple_hgnc))
  
  # Remove mzebra rows
  print("Removing old org rows...")
  ind <- which(rownames(new_counts_matrix) %in% all_hgnc[,2])
  new_counts_matrix <- new_counts_matrix[ind,]
  new_data_matrix   <- new_data_matrix[ind,]
  
  # Merge the duplicated rows
  print("Removing duplicated HGNC rows...")
  ptm <- proc.time()
  dup_genes <- unique(rownames(new_counts_matrix)[which(duplicated(rownames(new_counts_matrix)))])
  remove_dup_ind <- c()
  keep_dup_ind <- c()
  j <- 0
  dup_matrix_short_counts <- new_counts_matrix[which(rownames(new_counts_matrix) %in% dup_genes),]
  dup_matrix_short_data   <- new_data_matrix[which(rownames(new_data_matrix) %in% dup_genes),]
  keep_dup_matrix_counts  <- matrix(0L, nrow = length(dup_genes), ncol = ncol(obj@assays$RNA@counts))
  keep_dup_matrix_data    <- matrix(0L, nrow = length(dup_genes), ncol = ncol(obj@assays$RNA@data))
  rownames(keep_dup_matrix_counts) <- dup_genes
  rownames(keep_dup_matrix_data)   <- dup_genes
  for (gene in dup_genes) {
    ind_keep <- which(rownames(dup_matrix_short_counts) == gene)
    keep_dup_matrix_counts[j,] <- colSums(dup_matrix_short_counts[ind_keep,])
    keep_dup_matrix_data[j,]   <- colSums(dup_matrix_short_data[ind_keep,])
    j <- j + 1
  }
  # Delete all the duplicated rows at once
  print(paste("-Number of rows before merging:", nrow(new_counts_matrix)))
  print("-Actually doing the removal now I swear")
  remove_dup_ind <- which(rownames(new_counts_matrix) %in% dup_genes)
  new_counts_matrix <- new_counts_matrix[-remove_dup_ind,]
  new_data_matrix   <- new_data_matrix[-remove_dup_ind,]
  # nrow(new_counts_matrix)
  # Add all the new data at once
  print("-Combining the duplicated rows")
  new_counts_matrix <- rbind(new_counts_matrix, keep_dup_matrix_counts)
  new_data_matrix   <- rbind(new_data_matrix, keep_dup_matrix_data)
  print(paste("-Number of rows after merging:", nrow(new_counts_matrix)))
  proc.time() - ptm
  # new_counts_matrix[keep_dup_ind,] <- keep_dup_matrix_counts
  # new_data_matrix[keep_dup_ind,]   <- keep_dup_matrix_data
  # Delete temporary files
  rm(dup_matrix_short_counts)
  rm(dup_matrix_short_data)
  rm(keep_dup_matrix_counts)
  rm(keep_dup_matrix_data)
  
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
