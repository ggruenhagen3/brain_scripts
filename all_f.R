library("stringr")
library("ggplot2")
library("biomaRt")
library("Seurat")
library("Matrix")
library("reticulate")
library("cowplot")
library("RColorBrewer")
library("dplyr")

####################
# Helper Functions #
####################
convertToMouseObj <- function(obj) {
  # Converts a Mzebra Seurat object to a Mouse Seurat Object.
  
  # Convert a Seurat object to have all HGNC gene names
  print("Converting Genes Names...")
  genes <- rownames(obj@assays$RNA@counts)

  mouse  = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  mzebra = useMart(biomart="ensembl", dataset="mzebra_gene_ensembl")
  
  # DF to convert from org to HGNC
  all_hgnc <- getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = genes , mart = mzebra, attributesL = c("external_gene_name"), martL = mouse, uniqueRows=T)
  
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
  # obj_2$seurat_clusters <- obj$seurat_clusters
  
  # Add the metadata
  for (col in colnames(obj@meta.data)) {
    obj_2@meta.data[col] <- obj@meta.data[col]
  }
  
  return(obj_2)
}

convertMzebraGeneListToMouse <- function(gene_list) {
  mzebra = useEnsembl("ensembl", mirror = "useast", dataset = "mzebra_gene_ensembl")
  mouse  = useEnsembl("ensembl", mirror = "useast", dataset = "mmusculus_gene_ensembl")
  
  # ensembl_genes <- getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = gene_list , mart = mzebra, attributesL = c("external_gene_name"), martL = mouse, uniqueRows=T)
  # zfin_genes    <- getLDS(attributes = c("zfin_id_symbol"), filters = "zfin_id_symbol", values = gene_list , mart = mzebra, attributesL = c("external_gene_name"), martL = mouse, uniqueRows=T)
  # hgnc_genes    <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = gene_list , mart = mzebra, attributesL = c("external_gene_name"), martL = mouse, uniqueRows=T)
  # 
  # all_hgnc <- rbind(ensembl_genes, setNames(zfin_genes, names(ensembl_genes)), setNames(hgnc_genes, names(ensembl_genes)))
  # all_hgnc = all_hgnc[!duplicated(all_hgnc[,1]),]
  # all_mouse <- all_hgnc
  
  # DF to convert from org to HGNC
  all_mouse <- getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = gene_list, mart = mzebra, attributesL = c("external_gene_name"), martL = mouse, uniqueRows=T)
  
  mouse_genes <- unique(all_mouse[,2])
  print(paste0(length(mouse_genes)/length(gene_list) * 100, "% Genes Converted (", length(mouse_genes), "/", length(gene_list), ")"))
  return(all_mouse)
}

convertMzebraDFToMouse <- function(df, gene_column) {
  bad_genes <- df[,gene_column]
  bad_genes <- unique(bad_genes)
  gene_list <- bad_genes

  all_mouse <- convertMzebraGeneListToMouse(bad_genes)
  
  multiple_hgnc <- 0
  bad_multiple_hgnc <- 0
  converter <- data.frame()
  for (gene in bad_genes) {
    mouse_gene <- all_mouse[which(all_mouse[,1] == gene),2]
    new_mouse_gene <- mouse_gene[1]
    if (length(mouse_gene) > 1) {
      multiple_hgnc <- multiple_hgnc + 1
      upper_mouse_gene <- mouse_gene[which(startsWith(tolower(gene), tolower(mouse_gene)))]
      if (length(upper_mouse_gene) == 1) {
        new_mouse_gene <- upper_mouse_gene
      } else {
        bad_multiple_hgnc <- bad_multiple_hgnc + 1
        new_mouse_gene <- mouse_gene[1]
      } # end bad multiple
    } # end multiple
    converter <- rbind(converter, setNames(data.frame(t(c(gene, new_mouse_gene))), c("mzebra", "mouse")) )
  }
  print(paste0("Muliplte Mouse: ", multiple_hgnc))
  print(paste0("Bad Multiple Mouse: ", bad_multiple_hgnc))
  
  df[,gene_column] <- converter[match(df[,gene_column], converter[,1]),2]
  df <- df[which(! is.na(df[,gene_column])),]
  
  return(df)
}

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
  all_ind <- 1:nrow(new_counts_matrix)
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
  # obj_a_2$seurat_clusters <- obj_a$seurat_clusters
  
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
  all_ind <- 1:nrow(new_counts_matrix)
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
  # obj_b_2$seurat_clusters <- obj_b$seurat_clusters
  
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
  # obj_2$seurat_clusters <- obj$seurat_clusters
  
  # Add the metadata
  for (col in colnames(obj@meta.data)) {
    obj_2@meta.data[col] <- obj@meta.data[col]
  }
  
  return(obj_2)
}

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