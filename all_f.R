library("stringr")
library("ggplot2")
library("biomaRt")
library("Seurat")
library("Matrix")
library("reticulate")
library("cowplot")
library("RColorBrewer")
library("dplyr")
library("gplots")
library("viridis")
library("reshape2")

####################
# Helper Functions #
####################
printVectorAsNewVector <- function(vect) {
  str = "newVect = c("
  for (i in 1:length(vect)) {
    element = vect[i]
    if (i == length(vect)) { # if it's the last element, don't put a comma
      str = paste0(str, '"', element, '")') 
    } else {
      str = paste0(str, '"', element, '", ') 
    }
  }
  str = paste0(str, "")
  print(str)
  return(str)
}

bestCombo <- function(obj, targetCluster) {
  # Find the best combo of genes that are exclusive to that cluster
  Idents(obj) <- "seurat_clusters"
  targetCells <- WhichCells(obj, idents = targetCluster)
  targetGenes <- rownames(obj@assays$RNA@counts)[which(rowSums(obj@assays$RNA@counts[,targetCells]) > 100)]
  # maxTarget <- 0
  # minElse   <- ncol(obj)
  geneTargetCells <- list()
  geneOtherCells  <- list()
  for (i in targetGenes) {
    thisSum <- obj@assays$RNA@counts[i,]
    geneTargetCells[[i]] <- names(thisSum)[which(thisSum != 0 & names(thisSum) %in% targetCells)]
    geneOtherCells[[i]]  <- names(thisSum)[which(thisSum != 0 & ! names(thisSum) %in% targetCells)]
  }
  bestScore <- -ncol(obj)
  bestGenes <- c()
  n <- 0
  for (i in targetGenes) {
      if (n %% 1000 == 0) {
        print(n)
        print(bestScore)
      }
    for (j in targetGenes) {
      # print(paste(i,j))
      thisTarget <- length(geneTargetCells[[i]][which(geneTargetCells[[i]] %in% geneTargetCells[[j]])])
      thisOther  <- length(geneOtherCells[[i]][which(geneOtherCells[[i]] %in% geneOtherCells[[j]])])
      score <- thisTarget - thisOther
      if (score > bestScore) {
        bestScore <- score
        bestGenes <- c(i, j)
      }
    }
    n = n+1
  }
  # n <- 0
  # for (i in targetGenes) {
  #   if (n %% 1000 == 0) {
  #     print(n)
  #   }
  #   k <- 0
  #   for (j in targetGenes) {
  #     # print(k)
  #     thisSum <- colSums(obj@assays$RNA@counts[c(i,j),])
  #     thisTarget <- length(thisSum[which(thisSum != 0 & names(thisSum) %in% targetCells)])
  #     thisElse <- length(thisSum[which(thisSum != 0 & ! thisSum %in% targetCells)])
  #     if (thisTarget/thisElse > bestRatio) {
  #       # minElse <- thisElse
  #       # maxTarget <- thisTarget
  #       bestRatio <- thisTarget/thisElse
  #       bestGenes <- c(i, j)
  #     }
  #     k = k + 1
  #   }
  #   n = n + 1
  # }
}

myTotalTrans <- function(obj, slot="data", assay="RNA", features = NULL, cells = NULL) {
  result <- data.frame()
  for (ident in levels(Idents(obj))) {
    print(paste("Averaging Expression for", ident))
    if (identical(features, NULL)) {
      features = rownames(obj)
    }
    
    this_cells <- WhichCells(obj, idents = ident)
    if (! identical(cells, NULL)) {
      this_cells <- this_cells[which(this_cells %in% cells)]
    }
    
    if (length(features) == 1) {
      newRow <- setNames(as.data.frame(sum(GetAssayData(obj, slot=slot, assay=assay)[features,this_cells])), features)
      names(newRow) <- features
      result <- rbind(result, newRow)
    } else {
      result <- rbind(result, setNames(t(as.data.frame(rowSums(GetAssayData(obj, slot=slot, assay=assay)[features,this_cells]))), features)) 
    }
    rownames(result)[length(rownames(result))] <- ident
  }
  result <- as.data.frame(t(result))
  return(result)
}

myAverageExpression <- function(obj, slot="data", assay="RNA", features = NULL, cells=NULL) {
  result <- data.frame()
  for (ident in levels(Idents(obj))) {
    print(paste("Averaging Expression for", ident))
    if (identical(features, NULL)) {
      features = rownames(obj)
    }
    
    this_cells <- WhichCells(obj, idents = ident)
    if (! identical(cells, NULL)) {
      this_cells <- this_cells[which(this_cells %in% cells)]
    }
    
    if (length(features) == 1) {
      newRow <- setNames(as.data.frame(sum(GetAssayData(obj, slot=slot, assay=assay)[features,this_cells])/length(this_cells)), features)
      names(newRow) <- features
      result <- rbind(result, newRow)
    } else {
      result <- rbind(result, setNames(t(as.data.frame(rowSums(GetAssayData(obj, slot=slot, assay=assay)[features,this_cells]))/length(this_cells)), features)) 
    }
    rownames(result)[length(rownames(result))] <- ident
  }
  result <- as.data.frame(t(result))
  return(result)
}

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


convertMouseDataFrameToHgnc = function(mouse_df, gene_column) {
  human  = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse  = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
  m_genes = mouse_df[,gene_column]
  
  hgnc_df = getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = m_genes , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  colnames(hgnc_df) = c("mouse", "hgnc")
  hgnc_df = hgnc_df[which(! is.na(hgnc_df$hgnc)),]
  
  converter = hgnc_df
  mouse_df[,gene_column] <- converter[match(mouse_df[,gene_column], converter[,1]),2]
  mouse_df <- mouse_df[which(! is.na(mouse_df[,gene_column])),]
  
  return(mouse_df)
}

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

hgncMzebra <- function(genes, gene_names, onPACE=F) {
  # genes: genes to find hgnc conversions for
  # gene_names: list of ALL mzebra gene names
  # onPACE: is the script being run on PACE (default = FALSE)
  if (onPACE) {
    pat <- read.table("~/scratch/m_zebra_ref/MZ_treefam_annot_umd2a_ENS_2.bash", header = FALSE, fill = TRUE)
  } else {
    pat <- read.table("C:/Users/miles/Downloads/all_research/MZ_treefam_annot_umd2a_ENS_2.bash", header = FALSE, fill = TRUE)
  }
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

hgncMzebraInPlace <- function(df, gene_column, gene_names, onPACE=F) {
  bad_genes <- df[,gene_column]
  bad_genes <- unique(bad_genes)
  converter <- hgncMzebra(bad_genes, gene_names, onPACE)
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

heatmapComparisonMulti = function(dfs, samples, filename, filepath, labels=F, xlab=F, tri=T) {
  # Input: list of dataframes that are output of Seurat FindAllMarkers
  #        Vector of samples or whatever you want to name those two dataframes
  #        Base file name for png output
  #        Filepath for png output
  clusters = list()
  num_clusters = list()
  all_logFC = c()
  for (i in 1:length(dfs)) {
    clusters[[i]] = unique(as.vector(dfs[[i]]$cluster))
    num_clusters[[i]] = length(clusters[[i]])
    all_logFC = c(all_logFC, dfs[[i]]$avg_logFC)
  }
  
  # Now do Pairwise Comparison of each df's DEGs
  df = data.frame() # big df of all pairwise comparisons
  for (i in 1:length(dfs)) {
    for (i_clust in 1:num_clusters[[i]]) {
      i_clust_df = dfs[[i]][which(dfs[[i]]$cluster == clusters[[i]][i_clust]),]
      
      for (j in 1:length(dfs)) {
        for (j_clust in 1:num_clusters[[j]]) {
          j_clust_df = dfs[[j]][which(dfs[[j]]$cluster == clusters[[j]][j_clust]),]
          ovlp = nrow(j_clust_df[which(j_clust_df$gene %in% i_clust_df$gene),])
          ovlp_same_dir = nrow(j_clust_df[which(j_clust_df$gene %in% i_clust_df$gene & sign(j_clust_df$avg_logFC) == sign(i_clust_df$avg_logFC)),])
          
          pct = (2*ovlp / (nrow(i_clust_df) + nrow(j_clust_df))) * 100
          pct_same_dir = (ovlp_same_dir / (nrow(i_clust_df) + nrow(j_clust_df))) * 100
          
          # Rename the clusters with their sample names to avoid confusion
          sample1_clust = paste0(samples[[i]], " ", clusters[[i]][i_clust])
          sample2_clust = paste0(samples[[j]], " ", clusters[[j]][j_clust])
          df <- rbind(df, t(c(sample1_clust, sample2_clust, ovlp, pct, ovlp_same_dir, pct_same_dir)))
        }
      }
      
    }
    print(paste("Finished Pairwise Comparisons for", samples[[i]]))
  } # end parwise comparison
  
  colnames(df) <- c("df1_cluster", "df2_cluster", "ovlp", "pct", "ovlp_same_dir", "pct_same_dir")
  df$ovlp = as.numeric(as.vector(df$ovlp))
  df$pct = as.numeric(as.vector(df$pct))
  df$ovlp_same_dir = as.numeric(as.vector(df$ovlp_same_dir))
  df$pct_same_dir = as.numeric(as.vector(df$pct_same_dir))
  
  # Extract lower triangle
  if (tri) {
    print("Only Keeping Lower Triangle")
    new_df = data.frame()
    clusters = unique(df$df1_cluster)
    for (i in 1:length(clusters)) {
      i_clust = clusters[i]
      if (i == 1) print(i_clust)
      for (j in 1:i+1) {
        j_clust = clusters[j]
        if (i == 1) print(j_clust)
        new_df = rbind(new_df, df[which(df$df1_cluster == i_clust & df$df2_cluster == j_clust),])
      }
    }
    df = new_df
  }
  
  # Color for text label in heatmap
  df$id = rownames(df)
  df$ovlp_col = df$ovlp > mean(df$ovlp)
  df$pct_col = df$pct > mean(df$pct)
  df$ovlp_same_dir_col = df$ovlp_same_dir > mean(df$ovlp_same_dir)
  df$pct_same_dir_col = df$pct_same_dir > mean(df$pct_same_dir)
  
  # Find Max's
  df$ovlp_best = df$ovlp
  df$pct_best  = df$pct
  df$ovlp_same_dir_best = df$ovlp_same_dir
  df$pct_same_dir_best  = df$pct_same_dir
  for (cluster in unique(df$df1_cluster)) {
    # Find Ovlp Max
    rows = df[which(df$df1_cluster == cluster),]
    max_row = rows$id[which.max(rows$ovlp)]
    df$ovlp_best[which(df$df1_cluster == cluster & df$id != max_row)] = 0
    
    # Find Pct Max
    rows = df[which(df$df1_cluster == cluster),]
    max_row = rows$id[which.max(rows$pct)]
    df$pct_best[which(df$df1_cluster == cluster & df$id != max_row)] = 0
    
    # Find Ovlp Same Direction Max
    rows = df[which(df$df1_cluster == cluster),]
    max_row = rows$id[which.max(rows$ovlp_same_dir)]
    df$ovlp_same_dir_best[which(df$df1_cluster == cluster & df$id != max_row)] = 0
    
    # Find Pct Same Direction Max
    rows = df[which(df$df1_cluster == cluster),]
    max_row = rows$id[which.max(rows$pct_same_dir)]
    df$pct_same_dir_best[which(df$df1_cluster == cluster & df$id != max_row)] = 0
  }
  
  
  if (any(sign(all_logFC) == -1)) {
    print("Up and downregulated DEGs analyzed")
    png1_name = paste(filepath, filename, "_ovlp_same_dir.png", sep="")
    png2_name = paste(filepath, filename, "_best_guess_same_dir.png", sep="")
    png3_name = paste(filepath, filename, "_pct_same_dir.png", sep="")
    png4_name = paste(filepath, filename, "_pct_best_guess_same_dir.png", sep="")
    
    png1_title = paste("DEGs in Common w/ Same Sign b/w Clusters")
    png2_title = paste("Best Guess of DEGs w/ Same Sign")
    png3_title = paste("% DEGs w/ Same Sign in Common Clusters")
    png4_title = paste("% Best Guess of DEGs w/ Same Sign")
    
    df$ovlp = df$ovlp_same_dir
    df$ovlp_col = df$ovlp_same_dir_col
    df$ovlp_best = df$ovlp_same_dir_best
    df$pct = df$pct_same_dir
    df$pct_col = df$pct_same_dir_col
    df$pct_best = df$pct_same_dir_best
  } else {
    print("Only Upregulated Genes")
    png1_name = paste(filepath, filename, "_ovlp.png", sep="")
    png2_name = paste(filepath, filename, "_best_guess.png", sep="")
    png3_name = paste(filepath, filename, "_pct.png", sep="")
    png4_name = paste(filepath, filename, "_pct_best_guess.png", sep="")
    
    png1_title = paste("DEGs in Common b/w Clusters")
    png2_title = paste("Best Guess")
    png3_title = paste("% DEGs in Common b/w Clusters")
    png4_title = paste("% Best Guess")
  }
  
  # Plot 1 - Ovlp
  png(png1_name, width = 250*length(dfs)+50, height = 250*length(dfs), unit = "px", res = 100)
  p = ggplot(df, aes(df1_cluster, df2_cluster, fill=ovlp)) + geom_tile() + scale_fill_viridis(discrete=FALSE) +  ggtitle(png1_title) + guides(color = FALSE) + theme_classic() + theme(line = element_blank()) + coord_fixed() + theme(axis.text.x = element_text(angle = 45))
  if (labels)
    p = p + geom_text(aes(label=ovlp, color=ovlp_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000"))
  if (! xlab)
    p = p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  print(p)
  dev.off()
  print("finished plot 1")
  
  # Plot 2 - Ovlp Best Guess
  png(png2_name, width = 250*length(dfs)+50, height = 250*length(dfs), unit = "px", res = 100)
  p = ggplot(df, aes(df1_cluster, df2_cluster, fill=ovlp_best)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + ggtitle(png2_title) + guides(color = FALSE) + theme_classic() + theme(line = element_blank()) + coord_fixed() + theme(axis.text.x = element_text(angle = 45))
  if (labels)
    p = p + geom_text(data=subset(df, ovlp_same_dir_best > 0), aes(label=ovlp_same_dir_best, color=ovlp_same_dir_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000"))
  if (! xlab)
    p = p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  print(p)
  dev.off()
  print("finished plot 2")
  
  # Plot 3 - Pct
  png(png3_name,  width = 250*length(dfs)+50, height = 250*length(dfs), unit = "px", res = 100)
  p = ggplot(df, aes(df1_cluster, df2_cluster, fill=pct)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + ggtitle(png3_title) + guides(color = FALSE) + theme_classic() + theme(line = element_blank()) + coord_fixed() + theme(axis.text.x = element_text(angle = 45))
  if (labels)
    p = p + geom_text(aes(label=format(round(pct, 1), nsmall = 1), color=pct_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000")) 
  if (! xlab)
    p = p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  print(p)
  dev.off()
  print("finished plot 3")
  
  # Plot 4 - Pct Best Guess
  png(png4_name,  width = 250*length(dfs)+50, height = 250*length(dfs), unit = "px", res = 100)
  p = ggplot(df, aes(df1_cluster, df2_cluster, fill=pct_same_dir_best)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + ggtitle(png4_title) + guides(color = FALSE) + theme_classic() + theme(line = element_blank()) + coord_fixed() + theme(axis.text.x = element_text(angle = 45))
  if (labels)
    p = p + geom_text(data=subset(df, pct_same_dir_best > 0), aes(label=format(round(pct_same_dir_best, 1), nsmall = 1), color=pct_same_dir_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000"))
  if (! xlab)
    p = p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  print(p)
  dev.off()
  print("finished plot 4")
  
  return(df)
}

heatmapComparison <- function(df1, df2, df1_sample, df2_sample, filename, filepath) {
  # Input: 2 dataframes that are output of Seurat FindAllMarkers
  #        The samples or whatever you want to name those two dataframes
  #        Base file name for png output
  #        Filepath for png output
  # Output: 4 png files. A comparison of DEGs in the two input dataframes.
  df1_clusters = unique(as.vector(df1$cluster))
  df2_clusters = unique(as.vector(df2$cluster))
  df1_num_clusters <- length(df1_clusters)
  df2_num_clusters <- length(df2_clusters)
  df <- data.frame()
  for (i in 1:df1_num_clusters) {
    for (j in 1:df2_num_clusters) {
      df1_cluster <- df1[which(df1$cluster == df1_clusters[i]),]
      df2_cluster <- df2[which(df2$cluster == df2_clusters[j]),]
      ovlp <- nrow(df2[which(df2_cluster$gene %in% df1_cluster$gene),])
      ovlp_same_dir = nrow(df2[which(df2_cluster$gene %in% df1_cluster$gene & sign(df2_cluster$avg_logFC) == sign(df1_cluster$avg_logFC)),])
      pct <- (ovlp / (nrow(df1_cluster) + nrow(df2_cluster))) * 100
      pct_same_dir = (ovlp_same_dir / (nrow(df1_cluster) + nrow(df2_cluster))) * 100
      df <- rbind(df, t(c(df1_clusters[i], df2_clusters[j], ovlp, pct, ovlp_same_dir, pct_same_dir)))
    }
  }
  
  colnames(df) <- c("df1_cluster", "df2_cluster", "ovlp", "pct", "ovlp_same_dir", "pct_same_dir")
  df$ovlp = as.numeric(as.vector(df$ovlp))
  df$pct = as.numeric(as.vector(df$pct))
  df$ovlp_same_dir = as.numeric(as.vector(df$ovlp_same_dir))
  df$pct_same_dir = as.numeric(as.vector(df$pct_same_dir))
  
  # Color for text label in heatmap
  df$id = rownames(df)
  df$ovlp_col = df$ovlp > mean(df$ovlp)
  df$pct_col = df$pct > mean(df$pct)
  df$ovlp_same_dir_col = df$ovlp_same_dir > mean(df$ovlp_same_dir)
  df$pct_same_dir_col = df$pct_same_dir > mean(df$pct_same_dir)
  
  # Find Max's
  df$ovlp_best = df$ovlp
  df$pct_best  = df$pct
  df$ovlp_same_dir_best = df$ovlp_same_dir
  df$pct_same_dir_best  = df$pct_same_dir
  for (cluster in unique(df$df1_cluster)) {
    # Find Ovlp Max
    rows = df[which(df$df1_cluster == cluster),]
    max_row = rows$id[which.max(rows$ovlp)]
    df$ovlp_best[which(df$df1_cluster == cluster & df$id != max_row)] = 0
    
    # Find Pct Max
    rows = df[which(df$df1_cluster == cluster),]
    max_row = rows$id[which.max(rows$pct)]
    df$pct_best[which(df$df1_cluster == cluster & df$id != max_row)] = 0
    
    # Find Ovlp Same Direction Max
    rows = df[which(df$df1_cluster == cluster),]
    max_row = rows$id[which.max(rows$ovlp_same_dir)]
    df$ovlp_same_dir_best[which(df$df1_cluster == cluster & df$id != max_row)] = 0
    
    # Find Pct Same Direction Max
    rows = df[which(df$df1_cluster == cluster),]
    max_row = rows$id[which.max(rows$pct_same_dir)]
    df$pct_same_dir_best[which(df$df1_cluster == cluster & df$id != max_row)] = 0
  }
  
  # Plot 1 - Ovlp
  if (any(sign(df1$avg_logFC) == -1) || any(sign(df2$avg_logFC) == -1)) {
    png(paste(filepath, filename, "_ovlp_same_dir.png", sep=""),  width = 500, height = 500, unit = "px", res = 72)
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=ovlp_same_dir)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + geom_text(aes(label=ovlp_same_dir, color=ovlp_same_dir_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000")) + xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("DEGs in Common w/ Same Sign b/w", df1_sample, "&", df2_sample,  "Clusters")) + guides(color = FALSE) + theme_classic() + theme(line = element_blank()))
    dev.off()
  } else {
    png(paste(filepath, filename, "_ovlp.png", sep=""),  width = 500, height = 500, unit = "px", res = 72)
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=ovlp)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + geom_text(aes(label=ovlp, color=ovlp_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000")) + xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("DEGs in Common b/w", df1_sample, "&", df2_sample,  "Clusters")) + guides(color = FALSE) + theme_classic() + theme(line = element_blank()))
    dev.off()
  }
  print("finished 1")

  # Plot 2 - Ovlp Best Guess
  if (any(sign(df1$avg_logFC) == -1) || any(sign(df2$avg_logFC) == -1)) {
    png(paste(filepath, filename, "_best_guess_same_dir.png", sep=""),  width = 500, height = 500, unit = "px")
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=ovlp_same_dir_best)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + geom_text(data=subset(df, ovlp_same_dir_best > 0), aes(label=ovlp_same_dir_best, color=ovlp_same_dir_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000")) + xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("Best Guess of DEGs w/ Same Sign b/w", df1_sample, "&", df2_sample)) + guides(color = FALSE) + theme_classic() + theme(line = element_blank()))
    dev.off()
  } else {
    png(paste(filepath, filename, "_best_guess.png", sep=""),  width = 500, height = 500, unit = "px")
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=ovlp_best)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + geom_text(data=subset(df, ovlp_best > 0), aes(label=ovlp_best, color=ovlp_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000")) + xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("Best Guess b/w", df1_sample, "&", df2_sample)) + guides(color = FALSE) + theme_classic() + theme(line = element_blank()))
    dev.off()
  }
  print("finished 2")

  # Plot 3 - Pct
  if (any(sign(df1$avg_logFC) == -1) || any(sign(df2$avg_logFC) == -1)) {
    png(paste(filepath, filename, "_pct_same_dir.png", sep=""),  width = 500, height = 500, unit = "px")
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=pct_same_dir)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + geom_text(aes(label=format(round(pct_same_dir, 1), nsmall = 1), color=pct_same_dir_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000")) + xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("% DEGs w/ Same Sign in Common b/w", df1_sample, "&", df2_sample,  "Clusters")) + guides(color = FALSE) + theme_classic() + theme(line = element_blank()))
    dev.off()
  } else {
    png(paste(filepath, filename, "_pct.png", sep=""),  width = 500, height = 500, unit = "px")
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=pct)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + geom_text(aes(label=format(round(pct, 1), nsmall = 1), color=pct_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000")) + xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("% DEGs in Common b/w", df1_sample, "&", df2_sample,  "Clusters")) + guides(color = FALSE) + theme_classic() + theme(line = element_blank()))
    dev.off()
  }
  print("finished plot 3")
  
  # Plot 4 - Pct Best Guess
  if (any(sign(df1$avg_logFC) == -1) || any(sign(df2$avg_logFC) == -1)) {
    png(paste(filepath, filename, "_pct_best_guess_same_dir.png", sep=""),  width = 500, height = 500, unit = "px")
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=pct_same_dir_best)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + geom_text(data=subset(df, pct_same_dir_best > 0), aes(label=format(round(pct_same_dir_best, 1), nsmall = 1), color=pct_same_dir_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000")) + xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("% Best Guess of DEGs w/ Same Sign b/w", df1_sample, "&", df2_sample)) + guides(color = FALSE) + theme_classic() + theme(line = element_blank()))
    dev.off()
  } else {
    png(paste(filepath, filename, "_pct_best_guess.png", sep=""),  width = 500, height = 500, unit = "px")
    print(ggplot(df, aes(df1_cluster, df2_cluster, fill=pct_best)) + geom_tile() + scale_fill_viridis(discrete=FALSE) + geom_text(data=subset(df, pct_best > 0), aes(label=format(round(pct_best, 1), nsmall = 1), color=pct_col)) + scale_colour_manual(values=c("#FFFFFF", "#000000")) + xlab(paste(df1_sample, "Cluster")) + ylab(paste(df2_sample, "Cluster")) + ggtitle(paste("% Best Guess b/w", df1_sample, "&", df2_sample)) + guides(color = FALSE) + theme_classic() + theme(line = element_blank()))
    dev.off()
  }
  print("finished plot 4")
  
  return(df)
}


downsample <- function(combined, marker_genes, run) {
  set.seed(run)
  min_trans <- min(combined$nCount_RNA)
  gene_names <- rownames(combined@assays$RNA@counts)
  # new_matrix <- matrix(, nrow = nrow(combined@assays$RNA@counts), ncol = ncol(combined@assays$RNA@counts), dimnames = list(gene_names, colnames(combined@assays$RNA@counts)))
  # new_new_matrix <- matrix(, nrow=nrow(combined@assays$RNA@counts))
  marker_matrix  <- matrix(, nrow=length(marker_genes), ncol = ncol(combined@assays$RNA@counts), dimnames = list(marker_genes, colnames(combined@assays$RNA@counts)))
  i <- 0
  for (cell in colnames(combined@assays$RNA@counts)) {
    # i <- i + 1
    # if (i%%500 == 1) {
    #   print(cell)
    # }
    # start.time <- Sys.time()
    
    trans_names <- rep(gene_names, combined@assays$RNA@counts[,cell])
    ran_trans_names <- sample(trans_names, min_trans)
    ran_trans_names <- ran_trans_names[which(ran_trans_names %in% marker_genes)]
    ran_df <- as.data.frame(table(ran_trans_names))
    zero_gene_names <- marker_genes[which(! marker_genes %in% ran_trans_names)]
    zero_df <- setNames(data.frame(zero_gene_names <- zero_gene_names, Freq <- rep(0, length(zero_gene_names))), c("ran_trans_names", "Freq"))
    ran_df <- rbind(ran_df, zero_df)
    rownames(ran_df) <- ran_df$ran_trans_names
    ran_df <- ran_df[marker_genes,2]
    # new_matrix[,cell] <-as.matrix(ran_df)
    
    
    marker_matrix[,cell] <- as.matrix(ran_df)
    # new_new_matrix <- cbind(new_new_matrix, as.matrix(ran_df))
    
    # end.time <- Sys.time()
    # time.taken <- end.time - start.time
    # print(time.taken)
  }
  
  return(marker_matrix)
}

pickNewCells <- function(combined, num_clusters, num_cells) {
  new_cells <- c()
  for (i in 1:num_cells) {
    ran_cluster <- sample(0:num_clusters, 1)
    this_cells <- names(combined$seurat_clusters[which(combined$seurat_clusters == ran_cluster)])
    new_cells <- c(new_cells, sample(this_cells,1))
  }
  
  return(new_cells)
}

shuffleClusters <- function(combined) {
  # The selection process for a new cluster should be as follows:
  # 1. Pick a random cluster 0-40
  # 2. Pick a random cell from that cluster to be a part of the new cluster
  # This means that the new data set would likely have duplicate cells
  new_cells <- lapply(0:num_clusters, function(x) c())
  num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
  for (i in 0:num_clusters) {
    num_cells <- length(combined$seurat_clusters[which(combined$seurat_clusters == i)])
    new_cells[[i+1]] <- pickNewCells(combined, num_clusters, num_cells)
  }
  
  return(new_cells)
}