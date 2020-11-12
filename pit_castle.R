library("shiny")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("dplyr")
library("ggplot2")

combined <- readRDS("C:/Users/miles/Downloads/brain/brain_scripts/brain_shiny/data/combined.rds")
rna_path <- "C:/Users/miles/Downloads/brain/"
pit_ratio <- read.table(paste(rna_path, "data/cell_pit_castle.tsv", sep=""), header=FALSE, sep="\t", stringsAsFactors=FALSE)
all_cells <- colnames(combined)

found_cells <- c()
old_found <- c()
not_found <- c()
for (cell in pit_ratio$V1) {
  result <- grep(cell, all_cells, value = TRUE)
  if (length(result) > 0) {
    found_cells <- c(found_cells, result[1])
    old_found <- c(old_found, cell)
  } else {
    not_found <- c(not_found, cell)
  }
}

found <- subset(combined, cells = found_cells)
# This is only CTRL cells somehow

pit_ratio_2 <- pit_ratio[which(pit_ratio$V1 %in% old_found),]
matrix2 <- matrix(pit_ratio_2$V2, nrow = 1, ncol = nrow(pit_ratio_2), dimnames = list("pit-ratio", found_cells) )
found[["pit"]] <- CreateAssayObject(counts = matrix2)
FeaturePlot(found, features = "pit-ratio", reduction = "umap", pt.size = 2, label=TRUE, order = TRUE)
