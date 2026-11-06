library(shiny)
library(DT)
source("DRI_core.R")

ui <- fluidPage(
  titlePanel("DRI Calculator (v1.0)"),
  p("Upload your bulk RNA-seq expression matrix (CSV format, rows = genes, columns = samples)."),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("expr_file", "Upload Expression Matrix (CSV)", accept = ".csv"),
      textInput("tds_input", "TDS Genes (comma separated)", 
                value = "TG,TSHR,PAX8,SLC5A5,NKX2-1"),
      actionButton("run_btn", "Calculate DRI", class = "btn-primary")
    ),
    
    mainPanel(
      h3("Results"),
      DTOutput("dri_table"),
      downloadButton("download_btn", "Download DRI Scores (CSV)")
    )
  )
)

server <- function(input, output) {
  result_data <- reactiveVal(NULL)
  
  observeEvent(input$run_btn, {
    req(input$expr_file)
    
    # Try-catch block for error handling
    tryCatch({
      # Read the uploaded CSV
      expr <- read.csv(input$expr_file$datapath, row.names = 1, check.names = FALSE)
      
      # Parse the TDS genes input
      tds_genes <- trimws(strsplit(input$tds_input, ",")[[1]])
      
      # Call the core function
      res <- calculate_DRI(expr, tds_genes)
      
      # Store results
      df <- data.frame(Sample = names(res$DRI_score), DRI = as.numeric(res$DRI_score))
      result_data(df)
      
    }, error = function(e) {
      # Show error notification to the user
      showNotification(paste("Error:", e$message), type = "error", duration = 10)
    })
  })
  
  # Render table
  output$dri_table <- renderDT({ 
    req(result_data())
    datatable(result_data(), options = list(pageLength = 10)) 
  })
  
  # Download handler
  output$download_btn <- downloadHandler(
    filename = function() { 
      paste0("DRI_scores_", Sys.Date(), ".csv") 
    },
    content = function(file) { 
      write.csv(result_data(), file, row.names = FALSE) 
    }
  )
}