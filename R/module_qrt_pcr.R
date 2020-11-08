# Load packages ----
library(DT)
library(dplyr)
library(tidyr)

# Source helper functions ----
# source("helpers_qrt_pcr.R") # automatically sourced

# Setup ----
# delta_h9_list <- load_reference() # from helpers_qrt_pcr.R
data_reference <- load_all_reference_data() # from helpers_qrt_pcr.R
hk_genes_default <- c("ACTB", "GAPDH")

# UI logic ----
uiQPCR <- function(id, label = "qpcr") {
  # `NS(id)` returns a namespace function, which we save as `ns` and will
  # invoke later.
  ns <- NS(id)
  
  tabPanel("QPCR",
    sidebarLayout(
     # input
     sidebarPanel(
       textInput(ns("col_labels"), label = "Column labels", value = "", placeholder = "Paste column labels ..."),
       span(h6(textOutput(ns("col_label_count"))), style="color:gray"),
       textInput(ns("row_labels"), label = "Row labels", value = "", placeholder = "Paste row labels ..."),
       span(h6(textOutput(ns("row_label_count"))), style="color:gray"),
       fileInput(ns("files"), label = "File input", multiple=FALSE, placeholder = "No file selected"),
       span(h6(textOutput(ns("file_count"))), style="color:gray"),
       selectizeInput(ns("genes_housekeeping"), label="Housekeeping genes",
                      # selected = c("ACTB", "GAPDH"), # handled by server
                      choices = NULL, # handled by server
                      multiple = TRUE,
                      options = list(placeholder = "Select houskeeping genes")),
       selectInput(ns("reference_dataset"), "Reference dataset",
                   c("H9 v0 (2018)" = "h9_v0",
                     "H9 v1 (2020)" = "h9_v1"
                     # "RC v0 (2020)" = "rc17_v0"
                     ),
                   selected="h9_v1"),
       actionButton(ns("button_calculate"), "Calculate"),
       actionButton(ns("button_clear_input"), "Clear"),
     ),
     
     
     # output
     mainPanel(
       tabsetPanel(
         tabPanel("Fold Change", br(), dataTableOutput(ns("plate_processed"))),
         tabPanel("Ct", br(), dataTableOutput(ns("plate_raw"))),
         tabPanel("Reference", br(), dataTableOutput(ns("reference"))),
         tabPanel("Log",
            # show / hide this area
            # https://stackoverflow.com/questions/44790028/show-hide-entire-box-element-in-r-shiny
            # https://stackoverflow.com/questions/51333133/using-shinyjs-to-hide-show-ui-elements
            br(),
            h3("Column labels"),
            verbatimTextOutput(ns("col_labels")),
            h3("Row labels"),
            verbatimTextOutput(ns("row_labels")),
            h3("File input"),
            verbatimTextOutput(ns("files"))
         )
       ),
     )
    )
  )
}

# Server logic ----
serverQPCR <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      # -- Notifications -------------------------------------------------------
      notification_warning_nas <- NULL
      notification_warning_zeros <- NULL
      notification_warning_std <- NULL



      # -- Reactive events and values ------------------------------------------
      # Reactive values object for storing all reactive values
      values <- reactiveValues()

      # Column labels - split into vector and store in values
      observeEvent(input$col_labels, {
          values$col_labels <- unlist(strsplit(input$col_labels, split="\\s+")) # unlist cast to vector
      })

      # Row labels - split into vector and store in values
      observeEvent(input$row_labels, {
          values$row_labels <- unlist(strsplit(input$row_labels, split="\\s+")) # unlist cast to vector
      })

      # Housekeeping genes
      # Whenever new hk genes are chosen, update pre-computed hk delta.ct
      # debounce to reduce update freq TODO: adjust this sometime
      genes_housekeeping <- reactive(input$genes_housekeeping)
      genes_housekeeping_d <- genes_housekeeping %>% debounce(2000) 
      observeEvent(genes_housekeeping_d(), {
            # TODO: handle case where none are selected
            # update selected hk genes
            values$hk_genes_selected <- genes_housekeeping_d()

            # update pre-computed delta ct
            values$data_reference_delta_ct <- process_reference_data(values$data_reference_raw,
                                                                     values$hk_genes_selected)
      })

      # Reference dataset
      # Whenever new reference dataset is chosen, update the reactive values,
      # i.e. available hk genes to select, raw ref data and pre-computed hk delta.ct
      observeEvent(input$reference_dataset, {
          # update chosen ref dataset
          values$data_reference_raw <- data_reference[[input$reference_dataset]]

          # update pre-computed hk delta.ct
          values$data_reference_delta_ct <- process_reference_data(values$data_reference_raw,
                                                                   values$hk_genes_selected)

          # get options for which hk genes can be selected from new ref data
          hk_genes <- sort(data_reference[[input$reference_dataset]]$gene)

          # update reactiveValues
          values$hk_genes_selected <- hk_genes_default # Default hk genes

          # update input field
          updateSelectizeInput(session, "genes_housekeeping",
                               selected = values$hk_genes_selected, # Set default
                               choices = hk_genes,
                               server = TRUE)
      })

      # Clear button
      observeEvent(input$button_clear_input, {
        updateTextInput(session, "col_labels", value = "")
        updateTextInput(session, "row_labels", value = "")
      })

      # Calculate button
      plate_processed <- eventReactive(input$button_calculate, {
        # check args, before calling func
        validate(
          need(input$files, label="File input"),
          need(input$col_labels != '', label="Column labels"),
          need(input$row_labels != '', label="Row labels"),
          need((length(values$col_labels) %% 24) == 0, "Number of column labels must be a multiple of 24"),
          need((length(values$row_labels) %% 16) == 0, "Number of row labels must be a multiple of 16")
        )
        # call func
        process_plate(
          in_file=input$files,
          col_names=input$col_labels,
          row_names=input$row_labels,
          reference_data=values$data_reference_delta_ct)
      })



      # -- Output --------------------------------------------------------------

      # Number of col labels vs expected
      # for columns we expect a multiple of 24
      output$col_label_count <- renderText({
          count <- length(values$col_labels)
          expected <- ceiling_to_multiple(count-1, 24)
          return(paste0(count, "/", expected, " labels"))

      })

      # Number of row labels vs expected
      # for rows we expect a multiple of 16
      output$row_label_count <- renderText({
          count <- length(values$row_labels)
          expected <- ceiling_to_multiple(count-1, 16)
          return(paste0(count, "/", expected, " labels"))
      })

      # Number of files vs expected
      # for files (equivalent to 1 plate) we infer the expected number from the
      # number of rows and cols
      output$file_count <- renderText({
          count_files <- length(input$files$name) # count number of file names
          count_cols <- length(values$col_labels)
          count_rows <- length(values$row_labels)
          expected_x <- ceiling_to_multiple(count_cols-1, 24)
          expected_y <- ceiling_to_multiple(count_rows-1, 16)
          # n plates fit in x-axis * n plates fit in y-axis
          expected <- (expected_x / 24) * (expected_y / 16)
          return(paste0(count_files, "/", expected, " files"))
      })
      
      # Output tables
      output$plate_processed <- DT::renderDataTable({ 
        df <- plate_processed()$data_processed
        df <- datatable(
          data = df,
          extensions = c("Buttons"),
          options = list(
            dom = 'Bfrtip',
            pageLength=25,
            buttons = list('copy',
                           'csv',
                           list(extend = 'excel', filename="hippocompute_qpcr_fc", title = NULL),
                           'colvis')
            )
        ) %>%
          formatRound(1:ncol(df), 2)
        
        df })
      
      output$plate_raw <- DT::renderDataTable({ 
        result <- plate_processed()
        
        df <- result$data_raw_ct
        
        mask_zeros <- result$mask_zeros
        mask_nas <- result$mask_nas
        mask_std <- result$mask_std
        mask <- matrix(0L, nrow = dim(df)[1], ncol = dim(df)[2])
        mask[mask_std] <- 1
        mask[mask_zeros] <- 1
        mask[mask_nas] <- 2
        
        # set table cell color according to replacement, see:
        # https://stackoverflow.com/questions/50798941/r-shiny-rendertable-change-cell-colors
        # https://stackoverflow.com/questions/42569302/format-color-of-shiny-datatable-dt-according-to-values-in-a-different-dataset
        df <- datatable(
          data = cbind(df,mask),
          extensions = c("Buttons"),
          options=list(
            dom = 'Bfrtip',
            pageLength = 25,
            columnDefs = list(list(visible=FALSE, targets=c((1+ncol(df)):(ncol(df)+ncol(mask))))),
            buttons = list('copy',
                           'csv',
                           list(extend = 'excel', filename="hippocompute_qpcr_ct", title = NULL),
                           'colvis')
          ),
          selection = "single"
        ) %>%
          formatStyle(
            1:ncol(df),
            valueColumns=(1+ncol(df)):(ncol(df)+ncol(mask)),
            backgroundColor=styleEqual(c(1,2,3), c("khaki", "lightsalmon", "lightcoral"))
          ) %>%
          formatRound(1:ncol(df), 2)
        
        df })
      
      # Reference table
      output$reference <- DT::renderDataTable({
        df_org <- values$data_reference_raw
        
        df <- datatable(
          data=df_org,
          extensions=c("Buttons"),
          options=list(
            dom = 'Bfrtip',
            pageLength = 25,
            buttons = list('copy',
                        'csv',
                        list(extend = 'excel', filename="hippocompute_qpcr_reference", title = NULL),
                        'colvis')
            )
          ) # %>%
          # formatRound(1:ncol(df_org), 2) # from 2 to avoid overwriting gene col
        
        df
      })
      
      # Debug info -- render input as text
      output$col_labels <- renderPrint({
        col_labels <- strsplit(input$col_labels, split="\\s+")
        return(col_labels)
      })
      
      output$row_labels <- renderPrint({
        row_labels <- strsplit(input$row_labels, split="\\s+")
        return(row_labels)
      })
      
      output$files <- renderPrint({ str(input$files) })
    }
  )
}
