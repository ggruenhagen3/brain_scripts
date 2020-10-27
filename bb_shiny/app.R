#################
# BIG BRAIN APP #
#################
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("ggplot2")
library("shinycssloaders")
library("shinyWidgets")
# options(warn=-1)
combined <- readRDS("data/B1C1C2MZ_combined_031020.rds")
# "C:/Users/miles/Downloads/brain/"
# combined <- readRDS("C:/Users/miles/Downloads/brain/brain_scripts/brain_mz_shiny/data/B1C1C2MZ_combined_031020.rds")

num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
gene_names <- rownames(combined@assays$RNA)
cluster_choices <- c("all", "inhibitory", "excitatory", "non-neuronal", 0:num_clusters)
inhib_clusters <- c(0, 7, 8, 10, 11, 20, 21, 27, 31, 32, 35)
excit_clusters <- c(1, 2, 3, 5, 6, 9, 12, 13, 14, 15, 16, 18, 19, 22, 23, 24, 25, 28, 29, 30, 33, 34, 36, 37, 39)
non_clusters   <- c(17, 26, 38, 40)

# Move some gene names to the front just to present some easy, interesting options
iegs <- c("egr1", "bdnf", "fosb", "CREB1", "npas4a", "socs2", "dusp10")
inhib <- c("gad1b", "gad2", "slc32a1")
gene_names <- c(iegs, inhib, gene_names)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Paint Expression of Cells"),
  
  # Sidebar layout with input and output definitions ----
  # sidebarLayout(
  div(
    # Sidebar panel for inputs ----
    # sidebarPanel(
    wellPanel(
      tags$style(type="text/css", '#leftPanel { width:200px; float:left;}'),
      id = "leftPanel",
      
      selectizeInput(inputId = "gene", label = "Gene Name", choices = NULL, selected = "fosb", multiple = TRUE, options = list(maxOptions = 1000)),
      selectizeInput(inputId = "cluster", label = "Cluster(s) to Split", choices = NULL, multiple = TRUE, options = list(maxOptions = 1000)),
      textInput(inputId = "resolution", label = "Resolution", value = "0.2", width = NULL),
      textInput(inputId = "dims", label = "Dimensions", value = "20", width = NULL),
      textInput(inputId = "npcs", label = "Number of PCs", value = "30", width = NULL),
      radioButtons(inputId = "geneMode", label = "Split Mode", inline = TRUE, choices = list("any", "all")), 
      materialSwitch(inputId = "origClust", label = "Display Original Cluster", value = FALSE, width = "100%", status = "success")
      
    ),
    
    # Main panel for displaying outputs ----
    # mainPanel(
    div(
      style = "flex-grow:1; resize:horizontal; overflow: hidden",
      # mainPanel(
      tabsetPanel(id = "tabs", type = "tabs",
                  tabPanel("BHVE vs CTRL Overlap", value="bc_ovlp_plot", plotOutput("bc_ovlp_plot", width = "100%", height="500px") %>% withSpinner(color="#0dc5c1")),
                  tabPanel("BHVE vs CTRL", value="bcplot", plotOutput("bcplot", width = "100%", height="500px") %>% withSpinner(color="#0dc5c1")),
                  tabPanel("B1 vs C1 Overlap", value="b1c1_ovlp_plot", plotOutput("b1c1_ovlp_plot", width = "100%", height="500px") %>% withSpinner(color="#0dc5c1")),
                  tabPanel("B1 vs C1", value="b1c1plot", plotOutput("b1c1plot", width = "100%", height="500px") %>% withSpinner(color="#0dc5c1")),
                  tabPanel("B1 vs C2", value="b1c2plot", plotOutput("b1c2plot", width = "100%", height="500px") %>% withSpinner(color="#0dc5c1")),
                  tabPanel("C1 vs C2", value="c1c2plot", plotOutput("c1c2plot", width = "100%", height="500px") %>% withSpinner(color="#0dc5c1")),
                  tabPanel("All Samples", value="allplot", plotOutput("allplot", width = "100%", height="500px") %>% withSpinner(color="#0dc5c1")),
                  tabPanel("Split Cluster", value="splClust", plotOutput("splClust", width = "100%", height="500px") %>% withSpinner(color="#0dc5c1")),
                  tabPanel("Paint Split Cluster", value="pntSplClust", plotOutput("pntSplClust", width = "100%", height="500px") %>% withSpinner(color="#0dc5c1")),
                  tabPanel("BarPlot Summary", value="barplot", plotOutput("barplot", width = "100%", height="500px")  %>% withSpinner(color="#0dc5c1")),
                  tabPanel("BarPlot %", value="barplotAvg", plotOutput("barplotAvg", width = "100%", height="500px")  %>% withSpinner(color="#0dc5c1")),
                  tabPanel("Summary", value="summary", verbatimTextOutput("summary") %>% withSpinner(color="#0dc5c1")),
                  tabPanel("Summary %", value="summaryAvg", verbatimTextOutput("summaryAvg") %>% withSpinner(color="#0dc5c1"))
      ),
      downloadButton(outputId = "down", label = "Download the plot"),
      downloadButton(outputId = "degClust", label = "Download DEGs by Cluster"),
      downloadButton(outputId = "degClustCond", label = "Download DEGs by Cluster and Condition")
      # width = 10
      
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  updateSelectizeInput(session, "gene", choices = gene_names, server = TRUE)
  updateSelectizeInput(session, "cluster", choices = cluster_choices, server = TRUE)
  
  findOvlPCells <- function(genes) {
    clean_genes <- c()
    pos_cells <- c()
    for (gene in genes) {
      gene_lower <- tolower(gene)
      gene_upper <- toupper(gene)
      gene_title <- str_to_title(gene)
      if (gene_lower %in% gene_names) {
        gene <- gene_lower
      } else if (gene_upper %in% gene_names) {
        gene <- gene_upper
      } else if (gene_title %in% gene_names) {
        gene <- gene_title
      }
      clean_genes <- c(gene, clean_genes)
      expr <- FetchData(object = combined, vars = gene, slot = "counts")
      pos_cells <- c(colnames(combined[, which(x = expr > 0)]), pos_cells)
    }
    
    # Find the overlap of the positive cells
    num_genes <- length(genes)
    counts <- table(pos_cells)
    ovlp_cells <- rownames(counts[counts >= num_genes])
    combined <- SetIdent(combined, cells=ovlp_cells, value="overlapping_cells")
    combined <- SetIdent(combined, cells=setdiff(WhichCells(combined), ovlp_cells), value="non_overlapping_cells")
    combined$ovlp <- combined@active.ident
    return(combined$ovlp)
  }
  
  createSpltObj <- function(cluster, resolution, dims, npcs, geneMode) {
    combined@active.assay <- "RNA"
    Idents(combined) <- "seurat_clusters"
    split_cells <- c()
    for (this_cluster in cluster) {
      if(this_cluster == "all") {
        new_cells <- colnames(combined)
        print(length(new_cells))
      } else if (this_cluster == "inhibitory") {
        new_cells <- WhichCells(combined, idents = inhib_clusters )
      } else if (this_cluster == "excitatory") {
        new_cells <- WhichCells(combined, idents = excit_clusters)
      } else if (this_cluster == "non-neuronal") {
        new_cells <- WhichCells(combined, idents = non_clusters)
      } else {
        new_cells <- WhichCells(combined, idents = as.numeric(this_cluster))
      }
      split_cells <- unique(c(split_cells, new_cells))
    }
    if (length(cluster) == 0) {
      split_cells <- colnames(combined)
    }
    every_gene_cell <- c()
    if (length(input$gene) > 0) {
      print("splitting by gene(s)")
      for (gene in input$gene) {
        expr <- FetchData(object = combined, vars = gene, slot = "counts")
        gene_cells <- colnames(combined[, which(x = expr > 0)])
        every_gene_cell <- unique(c(every_gene_cell, gene_cells))
        if (geneMode == "all") {
          split_cells <- split_cells[which(split_cells %in% gene_cells)]
        }
      }
      if (geneMode == "any") {
        split_cells <- split_cells[which(split_cells %in% every_gene_cell)]
      }
    }
    
    Idents(combined) <- "seurat_clusters"
    combined$orig.cluster <- combined$seurat_clusters
    obj <- subset(combined, cells = split_cells)
    obj <- FindVariableFeatures(object = obj, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
    obj <- RunPCA(obj, npcs = as.numeric(npcs), verbose = FALSE)
    obj <- RunUMAP(obj, reduction = "pca", dims = 1:as.numeric(dims))
    obj <- FindNeighbors(obj, reduction = "umap", dims = 1:2)
    obj <- FindClusters(obj, resolution = as.numeric(resolution))
    return(obj)
  }
  
  createBarPlot <- function(average) {
    combined$ovlp <- findOvlPCells(input$gene)
    Idents(combined) <- combined$ovlp
    test <- combined[, WhichCells(combined, idents = "overlapping_cells")]
    total_b1 <- c()
    total_b2 <- c()
    total_c1 <- c()
    total_c2 <- c()
    Idents(test) <- test$sample
    Idents(combined) <- combined$sample
    test$cond.cluster <- paste(Idents(test), test$seurat_clusters, sep = "_")
    combined$cond.cluster <- paste(Idents(combined), combined$seurat_clusters, sep = "_")
    Idents(combined) <- combined$cond.cluster
    Idents(test) <- test$cond.cluster
    for (i in 0:num_clusters) {
      b1_cells_in_cluster <- 0
      b2_cells_in_cluster <- 0
      c1_cells_in_cluster <- 0
      c2_cells_in_cluster <- 0
      
      try(b1_cells_in_cluster <- length(WhichCells(test, idents = paste("b1", i, sep="_"))), silent=TRUE)
      try(b2_cells_in_cluster <- length(WhichCells(test, idents = paste("b2", i, sep="_"))), silent=TRUE)
      try(c1_cells_in_cluster <- length(WhichCells(test, idents = paste("c1", i, sep="_"))), silent=TRUE)
      try(c2_cells_in_cluster <- length(WhichCells(test, idents = paste("c2", i, sep="_"))), silent=TRUE)
      
      try(all_b1_cells_in_cluster <- length(WhichCells(combined, idents = paste("b1", i, sep="_"))), silent = TRUE)
      try(all_b2_cells_in_cluster <- length(WhichCells(combined, idents = paste("b2", i, sep="_"))), silent = TRUE)
      try(all_c1_cells_in_cluster <- length(WhichCells(combined, idents = paste("c1", i, sep="_"))), silent = TRUE)
      try(all_c2_cells_in_cluster <- length(WhichCells(combined, idents = paste("c2", i, sep="_"))), silent = TRUE)
      
      if (average == TRUE) {
        total_b1 <- c(total_b1, round( (b1_cells_in_cluster/all_b1_cells_in_cluster) * 100, 2))
        total_b2 <- c(total_b2, round( (b2_cells_in_cluster/all_b2_cells_in_cluster) * 100, 2))
        total_c1 <- c(total_c1, round( (c1_cells_in_cluster/all_c1_cells_in_cluster) * 100, 2))
        total_c2 <- c(total_c2, round( (c2_cells_in_cluster/all_c2_cells_in_cluster) * 100, 2))
      } else {
        total_b1 <- c(total_b1, b1_cells_in_cluster)
        total_b2 <- c(total_b2, b2_cells_in_cluster)
        total_c1 <- c(total_c1, c1_cells_in_cluster)
        total_c2 <- c(total_c2, c2_cells_in_cluster)
      }
      
    }
    df <- data.frame(condition <- c(rep("b1", length(total_b1)), rep("b2", length(total_b2)), rep("c1", length(total_c1)), rep("c2", length(total_c2))),
                     cluster_num <- c(0:num_clusters, 0:num_clusters, 0:num_clusters, 0:num_clusters),
                     value <- c(total_b1, total_b2, total_c1, total_c2))
    colnames(df) <- c("condition", "cluster_num", "value")
    my_title <- paste("Number of Cells Expressing", paste(input$gene, collapse = ' and '), "per Cluster")
    my_ylab <- "Number of Cells"
    if (average == TRUE) {
      my_title <- paste("% Cells Expressing", paste(input$gene, collapse = ' and '), "per Cluster")
      my_ylab <- "% Cells"
    }
    p <- ggplot(df, aes(fill=condition, x=cluster_num, y=value)) +
      geom_bar(position="dodge", stat="identity") +
      theme_minimal() +
      ggtitle(my_title) +
      xlab("Cluster") +
      ylab(my_ylab) +
      scale_x_continuous(breaks = 0:num_clusters) +
      geom_text(aes(label=value), vjust=1.6, color="black", position = position_dodge(0.9), size=3.5)
    theme_minimal()
    p
  }
  
  createSplPlot <- function(cluster, resolution, dims, npcs, geneMode, origClust) {
    obj <- createSpltObj(cluster, resolution, dims, npcs, geneMode)
    if (origClust) {
      print("Displaying two plots")
      p1 <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1.5) + ggtitle("New Clusters")
      Idents(obj) <- obj$orig.cluster 
      p2 <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1.5) + ggtitle("Old Clusters")
      plot_grid(p1,p2)
    } else {
      DimPlot(obj, reduction = "umap", split.by = "cond", label = TRUE, pt.size = 1.5)
    }
  }
  
  createSplFeaturePlot <- function(cluster, resolution, dims, npcs, geneMode) {
    obj <- createSpltObj(cluster, resolution, dims, npcs, geneMode)
    FeaturePlot(obj, features = c(input$gene), split.by = "cond", reduction = "umap", pt.size = 1.5, label=TRUE, order = TRUE)
  }
  
  createPlot <- function(split, samples) {
    combined@active.assay <- "RNA"
    if (split == "cond") {
      FeaturePlot(combined, features = c(input$gene), split.by = "cond", reduction = "umap", pt.size = 1.5, label=TRUE, order = TRUE)
    } else {
      Idents(combined) <- combined$sample
      only_b1_c1 <- combined[,WhichCells(combined,idents = samples)]
      Idents(only_b1_c1) <- only_b1_c1$seurat_clusters
      FeaturePlot(only_b1_c1, features = c(input$gene), split.by = "sample", reduction = "umap", pt.size = 1.5, label=TRUE, order = TRUE)
    }
  }
  
  createOvlpPlot <- function(split) {
    if (length(input$gene) > 1) {
      
      combined$ovlp <- findOvlPCells(input$gene)
      
      Idents(combined) <- combined$ovlp
      obj <- combined
      if (split == "sample") {
        Idents(combined) <- combined$sample
        only_b1_c1 <- combined[,WhichCells(combined,idents = c("b1", "c1"))]
        obj <- only_b1_c1
      }
      
      Idents(obj) <- obj$ovlp
      DimPlot(obj, reduction="umap", group.by = "ident", split.by=split, pt.size=2, order=TRUE)
      
    } else if (length(input$gene) == 1) {
      gene <- input$gene
      gene_lower <- tolower(gene)
      gene_upper <- toupper(gene)
      gene_title <- str_to_title(gene)
      if (gene_lower %in% gene_names) {
        gene <- gene_lower
      } else if (gene_upper %in% gene_names) {
        gene <- gene_upper
      } else if (gene_title %in% gene_names) {
        gene <- gene_title
      } 
      
      Idents(combined) <- combined$seurat_clusters
      obj <- combined
      if (split == "sample") {
        Idents(combined) <- combined$sample
        only_b1_c1 <- combined[,WhichCells(combined,idents = c("b1", "c1"))]
        obj <- only_b1_c1
        Idents(obj) <- obj$seurat_clusters
      }
      obj@active.assay <- "RNA"
      FeaturePlot(obj, features = c(gene), split.by = split, reduction = "umap", pt.size = 1.5, label=TRUE, order = TRUE) 
    }
  }
  
  createSummary <- function(average) {
    combined$ovlp <- findOvlPCells(input$gene)
    Idents(combined) <- combined$ovlp
    
    combined$sample.ovlp <- paste(Idents(combined), combined$sample, sep = "_")
    Idents(combined) <- combined$sample.ovlp
    combined$sample.ovlp.cluster <- paste(Idents(combined), combined$seurat_clusters, sep = "_")
    
    Idents(combined) <- combined$sample
    combined$sample.cluster <- paste(Idents(combined), combined$seurat_clusters, sep = "_")
    
    Idents(combined) <- combined$sample.ovlp
    num_b1_ovlp <- length(WhichCells(combined, idents = "overlapping_cells_b1"))
    num_b2_ovlp <- length(WhichCells(combined, idents = "overlapping_cells_b2"))
    num_c1_ovlp <- length(WhichCells(combined, idents = "overlapping_cells_c1"))
    num_c2_ovlp <- length(WhichCells(combined, idents = "overlapping_cells_c2"))
    
    num_b1 <- length(combined$sample[which(combined$sample == "b1")])
    num_b2 <- length(combined$sample[which(combined$sample == "b2")])
    num_c1 <- length(combined$sample[which(combined$sample == "c1")])
    num_c2 <- length(combined$sample[which(combined$sample == "c2")])
    
    pct_b1_ovlp <- (num_b1_ovlp/num_b1) * 100
    pct_b2_ovlp <- (num_b2_ovlp/num_b2) * 100
    pct_c1_ovlp <- (num_c1_ovlp/num_c1) * 100
    pct_c2_ovlp <- (num_c2_ovlp/num_c2) * 100
    pct_b_ovlp <- (num_b1_ovlp+num_b2_ovlp)/(num_b1+num_b2) * 100
    pct_c_ovlp <- (num_c1_ovlp+num_c2_ovlp)/(num_c1+num_c2) * 100
    
    cat(paste("Total Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", num_b1_ovlp+num_b2_ovlp+num_c1_ovlp, "\n", sep=""))
    cat(paste("BHVE vs CTRL Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", num_b1_ovlp+num_b2_ovlp, " vs ", num_c1_ovlp+num_b2_ovlp, "\n", sep=""))
    cat(paste("% BHVE vs CTRL Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", format(round(pct_b_ovlp, 2), nsmall = 2), "% vs ", format(round(pct_c_ovlp, 2), nsmall = 2), "% \n\n", sep=""))
    
    cat(paste("B1 Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", num_b1_ovlp, "\n", sep=""))
    cat(paste("B2 Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", num_b2_ovlp, "\n", sep=""))
    cat(paste("C1 Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", num_c1_ovlp, "\n", sep=""))
    cat(paste("C2 Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", num_c2_ovlp, "\n\n", sep=""))
    
    cat(paste("% B1 Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", format(round(pct_b1_ovlp, 2), nsmall = 2), "% \n", sep=""))
    cat(paste("% B2 Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", format(round(pct_b2_ovlp, 2), nsmall = 2), "% \n", sep=""))
    cat(paste("% C1 Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", format(round(pct_c1_ovlp, 2), nsmall = 2), "% \n", sep=""))
    cat(paste("% C2 Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", format(round(pct_c2_ovlp, 2), nsmall = 2), "% \n\n", sep=""))
    
    cat(paste("----------------------------------", "\n"))
    Idents(object = combined) <- combined$seurat_clusters
    for (i in 0:num_clusters) {
      num_b1_cells_cluster <- 0
      num_b2_cells_cluster <- 0
      num_c1_cells_cluster <- 0
      num_c2_cells_cluster <- 0
      
      Idents(combined) <- combined$sample.ovlp.cluster
      try(num_b1_cells_cluster <- length(WhichCells(combined, idents = paste("overlapping_cells_b1", i, sep="_"))), silent=TRUE)
      try(num_b2_cells_cluster <- length(WhichCells(combined, idents = paste("overlapping_cells_b2", i, sep="_"))), silent=TRUE)
      try(num_c1_cells_cluster <- length(WhichCells(combined, idents = paste("overlapping_cells_c1", i, sep="_"))), silent=TRUE)
      try(num_c2_cells_cluster <- length(WhichCells(combined, idents = paste("overlapping_cells_c2", i, sep="_"))), silent=TRUE)
      
      if (average == TRUE) {
        Idents(combined) <- combined$sample.cluster
        
        try(num_b1_cells_cluster <- (num_b1_cells_cluster/length(WhichCells(combined, idents = paste("b1", i, sep="_")))) * 100, silent = TRUE)
        try(num_b2_cells_cluster <- (num_b2_cells_cluster/length(WhichCells(combined, idents = paste("b2", i, sep="_")))) * 100, silent = TRUE)
        try(num_c1_cells_cluster <- (num_c1_cells_cluster/length(WhichCells(combined, idents = paste("c1", i, sep="_")))) * 100, silent = TRUE)
        try(num_c2_cells_cluster <- (num_c2_cells_cluster/length(WhichCells(combined, idents = paste("c2", i, sep="_")))) * 100, silent = TRUE)
        
        cat(paste("% B1 Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", format(round(num_b1_cells_cluster, 2), nsmall = 2), "% \n", sep=""))
        cat(paste("% B2 Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", format(round(num_b2_cells_cluster, 2), nsmall = 2), "% \n", sep=""))
        cat(paste("% C1 Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", format(round(num_c1_cells_cluster, 2), nsmall = 2), "% \n", sep=""))
        cat(paste("% C2 Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", format(round(num_c2_cells_cluster, 2), nsmall = 2), "% \n\n", sep=""))
      } else {
        Idents(combined) <- combined$sample.cluster
        cat(paste("B1 Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", num_b1_cells_cluster, "\n", sep=""))
        cat(paste("B2 Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", num_b2_cells_cluster, "\n", sep=""))
        cat(paste("C1 Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", num_c1_cells_cluster, "\n", sep=""))
        cat(paste("C2 Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", num_c2_cells_cluster, "\n\n", sep=""))
      }
      
    }
  }
  
  # Normal Plot
  output$bcplot <- renderPlot({
    createPlot("cond", c())
  })
  
  # Normal Plot (B1 vs C1)
  output$b1c1plot <- renderPlot({
    createPlot("split", c("b1", "c1"))
  })
  
  # Normal Plot (B1 vs C2)
  output$b1c2plot <- renderPlot({
    createPlot("split", c("b1", "c2"))
  })
  
  # Normal Plot (C1 vs C2)
  output$c1c2plot <- renderPlot({
    createPlot("split", c("c1", "c2"))
  })
  
  # Normal Plot (All Samples)
  output$allplot <- renderPlot({
    createPlot("split", c("b1", "b2", "c1", "c2"))
  })
  
  # Plots cells that are overlapping if >1 genes provided (BHVE vs CTRL)
  output$bc_ovlp_plot <- renderPlot({
    createOvlpPlot("cond")
  })
  
  # Plots cells that are overlapping if >1 genes provided (B1 vs C1)
  output$b1c1_ovlp_plot <- renderPlot({
    createOvlpPlot("sample")
  })
  
  output$barplot <- renderPlot({
    createBarPlot(FALSE)
  })
  
  output$barplotAvg <- renderPlot({
    createBarPlot(TRUE)
  })
  
  output$splClust <- renderPlot({
    if (length(input$cluster) > 0 || length(input$gene) > 0) {
      createSplPlot(input$cluster, input$resolution, input$dims, input$npcs, input$geneMode, input$origClust)
    }
  })
  
  output$pntSplClust <- renderPlot({
    if (length(input$cluster) > 0 && length(input$gene) > 0) { 
      createSplFeaturePlot(input$cluster, input$resolution, input$dims, input$npcs, input$geneMode)
    }
  })
  
  output$summary <- renderPrint({
    if (input$gene != "") {
      createSummary(FALSE)
    }
  })
  
  output$summaryAvg <- renderPrint({
    if (input$gene != "") {
      createSummary(TRUE)
    }
  })
  
  # Download Functionality
  output$degClustCond <- downloadHandler(
    filename = function() {
      paste0(paste0(c("degByClustCond", input$cluster, input$gene), collapse = "_"), ".tsv")
    }, 
    content = function(filename) {
      if (input$tabs == "splClust" || input$tabs == "pntSplClust") {
        obj <- createSpltObj(input$cluster, input$resolution, input$dims, input$npcs, input$geneMode)
        deg <- data.frame()
        obj$clust.cond <- paste0(obj$seurat_clusters, obj$cond)
        obj_num_clusters <- as.numeric(tail(levels(obj@meta.data$seurat_clusters), n=1))
        for (i in 0:obj_num_clusters) {
          newRow <- FindMarkers(obj, ident.1 = paste0("CTRL_", i), ident.2 = past0("BHVE_", i))
          newRow$gene <- rownames(newRow)
          newRow$cluster <- i
          deg <- rbind(deg, newRow)
        }
        write.table(deg, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
      }
      
    }
  )
  output$degClust <- downloadHandler(
    filename = function() {
      paste0(paste0(c("degByClust", input$cluster, input$gene), collapse = "_"), ".tsv")
    }, 
    content = function(filename) {
      if (input$tabs == "splClust" || input$tabs == "pntSplClust") {
        obj <- createSpltObj(input$cluster, input$resolution, input$dims, input$npcs, input$geneMode)
        deg <- FindAllMarkers(obj)
        deg$gene <- rownames(deg)
        write.table(deg, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
      }
      
    }
  )
  output$down <- downloadHandler(
    filename = function() {
      if (length(input$gene) > 1) {
        paste(paste(input$gene, collapse = '_'), ".png", sep="")
      } else {
        paste(input$gene, ".png", sep="")
      }
    }, 
    content = function(filename) {
      print(filename)
      print(input$gene)
      if (input$tabs == "bc_ovlp_plot") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createOvlpPlot("cond")
      } else if (input$tabs == "b1c1_ovlp_plot") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createOvlpPlot("sample")
      } else if (input$tabs == "b1c1plot") {
        height <- 500 * length(input$gene)
        png(filename = filename, width = 900, height = height, type="cairo")
        p <- createPlot("sample", c("b1", "c1"))
      } else if (input$tabs == "b1c2plot") {
        height <- 500 * length(input$gene)
        png(filename = filename, width = 900, height = height, type="cairo")
        p <- createPlot("sample", c("b1", "c2"))
      } else if (input$tabs == "c1c2plot") {
        height <- 500 * length(input$gene)
        png(filename = filename, width = 900, height = height, type="cairo")
        p <- createPlot("sample", c("c1", "c2"))
      } else if (input$tabs == "allplot") {
        height <- 500 * length(input$gene)
        png(filename = filename, width = 900, height = height, type="cairo")
        p <- createPlot("sample", c("b1", "b2", "c1", "c2"))
      } else if (input$tabs == "bcplot")  {
        height <- 500 * length(input$gene)
        png(filename = filename, width = 900, height = height, type="cairo")
        p <- createPlot("cond", c())
      } else if (input$tabs == "barplot") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createBarPlot(FALSE)
      } else if (input$tabs == "barplotAvg") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createBarPlot(TRUE)
      } else if (intput$tabs == "splClust") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createSplPlot(input$cluster, input$resolution, input$dims, input$npcs, input$geneMode, input$origClust)
      } else if (intput$tabs == "pntSplClust") {
        png(filename = filename, width = 900, height = 500, type="cairo")
        p <- createSplPlot(input$cluster, input$resolution, input$dims, input$npcs, input$geneMode)
      }
      print(p)
      dev.off()
    }
  )
  
}

shinyApp(ui = ui, server = server)