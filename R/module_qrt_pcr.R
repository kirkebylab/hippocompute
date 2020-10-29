# Load packages ----
library(DT)
library(dplyr)
library(tidyr)

# Source helper functions ----
# source("helpers_qrt_pcr.R") # automatically sourced

# Setup ----
# delta_h9_list <- load_reference() # from helpers_qrt_pcr.R
data_reference <- load_all_reference_data() # from helpers_qrt_pcr.R

# UI logic ----
tabQrtpcrUI <- function(id, label = "qpcr") {
  # `NS(id)` returns a namespace function, which we save as `ns` and will
  # invoke later.
  ns <- NS(id)
  
  tabPanel("QPCR",
    sidebarLayout(
     # input
     sidebarPanel(
       textInput(ns("col_labels"), label = "Column labels", value = "", placeholder = "Paste column labels ..."),
       textInput(ns("row_labels"), label = "Row labels", value = "", placeholder = "Paste row labels ..."),
       fileInput(ns("data_file"), label = "File input", multiple=FALSE, placeholder = "No file selected"),
       selectInput(ns("reference_dataset"), "Reference",
                   c("H9 v0 (2018)" = "h9_v0"
                     # "RC v0 (2020)" = "rc17_v0"
                     )),
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
            verbatimTextOutput(ns("data_file"))
         )
       ),
     )
    )
  )
}

# Server logic ----
tabQrtpcrServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      notification_warning_nas <- NULL
      notification_warning_zeros <- NULL
      notification_warning_std <- NULL
      
      
      # Reactive events
      observeEvent(input$button_clear_input, {
        updateTextInput(session, "col_labels", value = "")
        updateTextInput(session, "row_labels", value = "")
      })
      
      plate_processed <- eventReactive(input$button_calculate, {
        validate(
          need(input$data_file, label="File input"),
          need(input$col_labels != '', label="Column labels"),
          need(input$row_labels != '', label="Row labels"),
          need(length(unlist(strsplit(input$col_labels, split="\n|\t| "))) > 24, "Too many column labels provided. Max: 24"), # unlist cast to vector
          need(length(unlist(strsplit(input$row_labels, split="\n|\t| "))) > 16, "Too many row labels provided. Max: 16") # unlsit cast to vector))
        )
        process_plate(
          in_file=input$data_file, 
          col_names=input$col_labels, 
          row_names=input$row_labels, 
          reference_data=data_reference[[input$reference_dataset]])
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
        df_org <- data_reference[[input$reference_dataset]] # load the selected ref dataset via key
        df1 <- df_org[1] # col for hk gene ACTB
        df2 <- df_org[2] # col for hk gene GAPDH
        df_comb <- merge(df1, df2, by.x=1, by.y=1)
        row.names(df_comb) <- df_comb[,1] # use gene names as index
        df_comb[,1] <- NULL
        
        df <- datatable(
          data=df_comb,
          extensions=c("Buttons"),
          options=list(
            dom = 'Bfrtip',
            pageLength = 25,
            buttons = list('copy',
                        'csv',
                        list(extend = 'excel', filename="hippocompute_qpcr_reference", title = NULL),
                        'colvis')
            )
          ) %>%
          formatRound(1:ncol(df_comb), 2) # from 2 to avoid overwriting gene col
        
        df
      })
      
      # Debug info -- render input as text
      output$col_labels <- renderPrint({
        col_labels <- strsplit(input$col_labels, split="\n|\t| ")
        return(col_labels)
      })
      
      output$row_labels <- renderPrint({
        row_labels <- strsplit(input$row_labels, split="\n|\t| ")
        return(row_labels)
      })
      
      output$data_file <- renderPrint({ str(input$data_file) })
    }
  )
}
