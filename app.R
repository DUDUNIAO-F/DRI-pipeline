library(shiny)
library(DT)
source("DRI_core.R")

ui <- fluidPage(
  titlePanel("DRI Calculator (v1.0)"),
  p("Upload your bulk RNA-seq expression matrix (CSV format, rows = genes, columns = samples)."),

  sidebarLayout(
    sidebarPanel(
      fileInput("expr_file", "Upload Expression Matrix (CSV)", accept = ".csv"),
      helpText("CSV: first column = gene/tsRNA ID; remaining columns = samples; numeric; no NA/Inf/negative; UTF-8."),
      downloadLink("download_template", "Download example input template (CSV)", class = "btn btn-default btn-sm"),
      tags$details(
        tags$summary("Input format requirements"),
        tags$ul(
          tags$li("First column: unique gene/tsRNA IDs, no duplicates."),
          tags$li("Other columns: one per sample; column names must be non-empty."),
          tags$li("Values: numeric expression (log2CPM/TPM/FPKM); no NA, Inf, or negative values."),
          tags$li("CSV encoding: UTF-8."),
          tags$li("The template already contains TDS genes and DRI genes and can be uploaded directly as a test.")
        )
      ),
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

  output$download_template <- downloadHandler(
    filename = function() "DRI_input_template.csv",
    content = function(file) {
      template_path <- file.path("www", "DRI_input_template.csv")
      if (!file.exists(template_path)) {
        stop("Template file not found: www/DRI_input_template.csv")
      }
      file.copy(template_path, file, overwrite = TRUE)
    }
  )

  observeEvent(input$run_btn, {
    req(input$expr_file)

    tryCatch({
      expr <- read.csv(input$expr_file$datapath, row.names = 1, check.names = FALSE)
      tds_genes <- trimws(strsplit(input$tds_input, ",")[[1]])
      res <- calculate_DRI(expr, tds_genes)
      df <- data.frame(Sample = names(res$DRI_score), DRI = as.numeric(res$DRI_score))
      result_data(df)
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error", duration = 10)
    })
  })

  output$dri_table <- renderDT({
    req(result_data())
    datatable(result_data(), options = list(pageLength = 10))
  })

  output$download_btn <- downloadHandler(
    filename = function() {
      paste0("DRI_scores_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(result_data(), file, row.names = FALSE)
    }
  )
}
shinyApp(ui = ui, server = server)
