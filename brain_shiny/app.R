library("shiny")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
# options(warn=-1)
combined <- readRDS("data/combined.rds")
gene_names <- rownames(combined@assays$RNA)

# Move some gene names to the front just to present some easy, interesting options
iegs <- c("egr1", "bdnf", "fosb", "CREB1", "npas4a", "socs2", "dusp10")
inhib <- c("gad1b", "gad2", "slc32a1")
gene_names <- c(iegs, inhib, gene_names)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Paint Expression of Cells"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # textInput(inputId = "gene", label = "Gene Name", value = "fosb", width = NULL,
      #           placeholder = "name of gene")
      selectizeInput(inputId = "gene", label = "Gene Name", choices = NULL, selected = "fosb", multiple = TRUE, options = list(maxOptions = 10))

    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      tabsetPanel(id = "tabs", type = "tabs",
                  tabPanel("BHVE vs CTRL Overlap", value="bc_ovlp_plot", plotOutput("bc_ovlp_plot", width = "100%", height="500px")),
                  tabPanel("BHVE vs CTRL", value="bcplot", plotOutput("bcplot", width = "100%", height="500px")),
                  tabPanel("B1 vs C1 Overlap", value="b1c1_ovlp_plot", plotOutput("b1c1_ovlp_plot", width = "100%", height="500px")),
                  tabPanel("B1 vs C1", value="b1c1plot", plotOutput("b1c1plot", width = "100%", height="500px")),
                  tabPanel("BarPlot Summary", value="barplot", plotOutput("barplot", width = "100%", height="500px")),
                  tabPanel("Summary", value="summary", verbatimTextOutput("summary"))
      ),
      downloadButton(outputId = "down", label = "Download the plot")
      
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, "gene", choices = gene_names, server = TRUE)
  
  createPlot <- function(split) {
    if (split == "cond") {
      FeaturePlot(combined, features = c(input$gene), split.by = "cond", reduction = "umap", pt.size = 2, label=TRUE, order = TRUE)
    } else {
      Idents(combined) <- combined$sample
      only_b1_c1 <- combined[,WhichCells(combined,idents = c("b1", "c1"))]
      Idents(only_b1_c1) <- only_b1_c1$seurat_clusters
      FeaturePlot(only_b1_c1, features = c(input$gene), split.by = "sample", reduction = "umap", pt.size = 2, label=TRUE, order = TRUE)
    }
  }
  
  createOvlpPlot <- function(split) {
    if (length(input$gene) > 1) {
      
      # Get all the genes in the proper format
      genes <- input$gene
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
        expr <- FetchData(object = combined, vars = gene)
        pos_cells <- c(colnames(combined[, which(x = expr > 1)]), pos_cells)
      }
      
      # Find the overlap of the positive cells
      num_genes <- length(genes)
      counts <- table(pos_cells)
      ovlp_cells <- rownames(counts[counts >= num_genes])
      combined <- SetIdent(combined, cells=ovlp_cells, value="overlapping_cells")
      combined <- SetIdent(combined, cells=setdiff(WhichCells(combined), ovlp_cells), value="non_overlapping_cells")
      combined$ovlp <- combined@active.ident
      
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
      
      FeaturePlot(obj, features = c(gene), split.by = split, reduction = "umap", pt.size = 2, label=TRUE, order = TRUE) 
    }
  }
  
  # Normal Plot
  output$bcplot <- renderPlot({
    createPlot("cond")
  })
  
  # Normal Plot (B1 vs C1)
  output$b1c1plot <- renderPlot({
    createPlot("split")
  })
  
  # Plots cells that are overlapping if >1 genes provided (BHVE vs CTRL)
  output$bc_ovlp_plot <- renderPlot({
    createOvlpPlot("cond")
  })
  
  # Plots cells that are overlapping if >1 genes provided (B1 vs C1)
  output$b1c1_ovlp_plot <- renderPlot({
    createOvlpPlot("sample")
  })
  
  output$summary <- renderPrint({
    
    # Get all the genes in the proper format
    genes <- input$gene
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
      expr <- FetchData(object = combined, vars = gene)
      pos_cells <- c(colnames(combined[, which(x = expr > 1)]), pos_cells)
    }
    
    # Find the overlap of the positive cells
    num_genes <- length(genes)
    counts <- table(pos_cells)
    ovlp_cells <- rownames(counts[counts >= num_genes])
    combined <- SetIdent(combined, cells=ovlp_cells, value="overlapping_cells")
    combined <- SetIdent(combined, cells=setdiff(WhichCells(combined), ovlp_cells), value="non_overlapping_cells")
    
    combined$sample.ovlp <- paste(Idents(combined), combined$sample, sep = "_")
    Idents(combined) <- combined$sample.ovlp
    combined$sample.ovlp.cluster <- paste(Idents(combined), combined$seurat_clusters, sep = "_")
    num_b1_ovlp <- length(WhichCells(combined, idents = "overlapping_cells_b1"))
    num_b2_ovlp <- length(WhichCells(combined, idents = "overlapping_cells_b2"))
    num_c1_ovlp <- length(WhichCells(combined, idents = "overlapping_cells_c1"))
    num_b1 <- length(combined$sample[which(combined$sample == "b1")])
    num_b2 <- length(combined$sample[which(combined$sample == "b2")])
    num_c1 <- length(combined$sample[which(combined$sample == "c1")])
    pct_b1_ovlp <- (num_b1_ovlp/num_b1) * 100
    pct_b2_ovlp <- (num_b2_ovlp/num_b2) * 100
    pct_c1_ovlp <- (num_c1_ovlp/num_c1) * 100
    pct_b_ovlp <- (num_b1_ovlp+num_b2_ovlp)/(num_b1+num_b2) * 100
    
    cat(paste("Total Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", length(ovlp_cells), "\n", sep=""))
    cat(paste("BHVE vs CTRL Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", num_b1_ovlp+num_b2_ovlp, " vs ", num_c1_ovlp, "\n", sep=""))
    cat(paste("% BHVE vs CTRL Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", format(round(pct_b_ovlp, 2), nsmall = 2), "% vs ", format(round(pct_c1_ovlp, 2), nsmall = 2), "% \n\n", sep=""))
    
    cat(paste("B1 Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", num_b1_ovlp, "\n", sep=""))
    cat(paste("B2 Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", num_b2_ovlp, "\n", sep=""))
    cat(paste("C1 Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", num_c1_ovlp, "\n\n", sep=""))
    
    cat(paste("% B1 Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", format(round(pct_b1_ovlp, 2), nsmall = 2), "% \n", sep=""))
    cat(paste("% B2 Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", format(round(pct_b2_ovlp, 2), nsmall = 2), "% \n", sep=""))
    cat(paste("% C1 Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", format(round(pct_c1_ovlp, 2), nsmall = 2), "% \n\n", sep=""))
    
    cat(paste("----------------------------------", "\n"))
    Idents(object = combined) <- combined$seurat_clusters
    num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
    Idents(combined) <- combined$sample.ovlp.cluster
    for (i in 1:num_clusters) {
      num_b1_cells_cluster <- 0
      num_b2_cells_cluster <- 0
      num_c1_cells_cluster <- 0
      try(num_b1_cells_cluster <- length(WhichCells(combined, idents = paste("overlapping_cells_b1", i, sep="_"))), silent=TRUE)
      try(num_b2_cells_cluster <- length(WhichCells(combined, idents = paste("overlapping_cells_b2", i, sep="_"))), silent=TRUE)
      try(num_c1_cells_cluster <- length(WhichCells(combined, idents = paste("overlapping_cells_c1", i, sep="_"))), silent=TRUE)
      # tryCatch({
      #   num_b2_cells_cluster <- length(WhichCells(combined, idents = paste("overlapping_cells_b2", i, sep="_")))
      # }, error=function(e)   {   num_b2_cells_cluster <- 0   })
      # tryCatch({
      #   num_c1_cells_cluster <- length(WhichCells(combined, idents = paste("overlapping_cells_c1", i, sep="_")))
      # }, error=function(e)   {   num_c1_cells_cluster <- 0   })

      cat(paste("B1 Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", num_b1_cells_cluster, "\n", sep=""))
      cat(paste("B2 Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", num_b2_cells_cluster, "\n", sep=""))
      cat(paste("C1 Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", num_c1_cells_cluster, "\n\n", sep=""))
    }
  })
  
  # Download Functionality
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
        p <- createPlot("sample")
      } else if (input$tabs == "bcplot")  {
        height <- 500 * length(input$gene)
        png(filename = filename, width = 900, height = height, type="cairo")
        p <- createPlot("cond")
      }
      print(p)
      dev.off()
    }
  )
  
}

shinyApp(ui = ui, server = server)