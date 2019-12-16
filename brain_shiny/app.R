library(shiny)
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
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
      
      # Input: Slider for the number of bins ----
      # sliderInput(inputId = "bins",
      #             label = "Number of bins:",
      #             min = 1,
      #             max = 50,
      #             value = 30)
      
      # textInput(inputId = "gene", label = "Gene Name", value = "fosb", width = NULL,
      #           placeholder = "name of gene")
      selectizeInput(inputId = "gene", label = "Gene Name", choices = NULL, selected = "fosb", multiple = TRUE, options = list(maxOptions = 10))

    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      # plotOutput(outputId = "plot", width = "100%", height="500px"),
      # downloadButton(outputId = "down", label = "Download the plot")
      
      tabsetPanel(type = "tabs",
                  tabPanel("Overlap Plot", plotOutput("ovlp_plot")),
                  tabPanel("Plot", plotOutput("plot")),
                  downloadButton(outputId = "down", label = "Download the plot")
      )
      
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, "gene", choices = gene_names, server = TRUE)
  
  
  createPlot <- function() {
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
      
      DimPlot(combined, reduction="umap", group.by = "ident", split.by="cond", pt.size=2, order=TRUE)
      
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
      
      FeaturePlot(combined, features = c(gene), split.by = "cond", reduction = "umap", pt.size = 2, label=TRUE, order = TRUE) 
    }
  }
  
  # Normal Plot
  output$plot <- renderPlot({
    
   FeaturePlot(combined, features = c(input$gene), split.by = "cond", reduction = "umap", pt.size = 2, label=TRUE, order = TRUE)

    
  })
  
  # Plots cells that are overlapping if >1 genes provided.
  output$ovlp_plot <- renderPlot({
    createPlot()
  })
  
  output$down <- downloadHandler(
    filename = function() {
      if (length(input$gene) > 1) {
        paste(paste(input$gene, collapse = '_'), ".png", sep="")
      } else {
        paste(input$gene, ".png", sep="")
      }
    }, 
    content = function(filename) {
      # plotPNG(createPlot, filename = filename, width = 400, height = 400, res = 72)  
      print(filename)
      print(input$gene)
      png(filename = filename, width = 900, height = 500)
      p <- createPlot()
      print(p)
      dev.off()
    }
  )
  
}

shinyApp(ui = ui, server = server)