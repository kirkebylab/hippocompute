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
       # input fields
       h6(textOutput(ns("col_label_count"))),
       textInput(ns("col_labels"), label = NULL, value = "", placeholder = "Paste column labels ..."),
       h6(textOutput(ns("row_label_count"))),
       textInput(ns("row_labels"), label = NULL, value = "", placeholder = "Paste row labels ..."),
       h6(textOutput(ns("file_count"))),
       fileInput(ns("files"), label = NULL, multiple=TRUE, placeholder = "No file selected"),

       # buttons
       fluidRow(
           column(12, align="center", actionButton(ns("button_calculate"), "Calculate",
                                                   class="btn-primary")),
       # actionButton(ns("button_clear_input"), "Clear"),
       ),

       hr(),

       # settings
       selectizeInput(ns("genes_housekeeping"), label="Housekeeping genes",
                      # selected = c("ACTB", "GAPDH"), # handled by server
                      choices = NULL, # handled by server
                      multiple = TRUE,
                      options = list(placeholder = "Select houskeeping genes")),
       selectInput(ns("reference_dataset"), "Reference dataset",
                   # key-value for selecting ref data
                   # the value must match the keys used in load_all_reference_data()
                   c("H9 v0 (2018)" = "h9_v0_2018",
                     "H9 v0 (2020)" = "h9_v0_2020",
                     "H9 v1 (2020)" = "h9_v1_2020",
                     "H9-RC17-avg v1 (2020)" = "h9-rc17_v1_2020",
                     "RC v1 (2020)" = "rc17_v1_2020"
                     ),
                   selected="h9_v1_2020"),
     ),
     
     
     # output
     mainPanel(
       tabsetPanel(id=ns("tabset"), selected="panel_input",
        tabPanel("Input", value="panel_input",
                br(),
                h3("Column labels"),
                verbatimTextOutput(ns("col_labels")),
                h3("Row labels"),
                verbatimTextOutput(ns("row_labels")),
                h3("File input"),
                verbatimTextOutput(ns("files"))
                ),
         tabPanel("Ct", value="panel_ct", br(), uiOutput(ns("plates"))),
         tabPanel("Fold Change", value="panel_fc", br(), dataTableOutput(ns("plate_fold_change"))),
         tabPanel("Reference", value="panel_ref", br(), dataTableOutput(ns("reference")))
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
      ns <- session$ns # session is needed for renderUI to work
      # -- Notifications -------------------------------------------------------
      notification_warning_nas <- NULL
      notification_warning_zeros <- NULL
      notification_warning_std <- NULL


      # -- Reactive events and values ------------------------------------------
      # Reactive values object for storing all reactive values
      values <- reactiveValues()

      values$plates <- NULL # initialize plates to NULL

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
      observeEvent(input$button_calculate, {
        # check args, before calling func
        # TODO: fix validate -- messages are not showing
        validate(
          need(input$files, label="File input"),
          need(input$col_labels != '', label="Column labels"),
          need(input$row_labels != '', label="Row labels"),
          need((length(values$col_labels) %% 24) == 0, "Number of column labels must be a multiple of 24"),
          need((length(values$row_labels) %% 16) == 0, "Number of row labels must be a multiple of 16")
        )

        # call funcs
        tryCatch({
            values$plates <- process_plates(
              files=input$files,
              col_names=input$col_labels,
              row_names=input$row_labels,
              reference_data=values$data_reference_delta_ct)
            values$dt_list <- make_datatables_ct(values$plates$raw)
            updateTabsetPanel(session, "tabset", selected="panel_ct")},
            warning = function(warn){print(warn)
                                  showNotification(paste0(warn),
                                                      duration=0,
                                                      closeButton=TRUE,
                                                      type="warning")},
            error = function(err){print(err)
                                  showNotification(paste0(err),
                                                    duration=0,
                                                    closeButton=TRUE,
                                                    type="err")}
        )
      })


      # -- Rendered Output -----------------------------------------------------

      # Number of col labels vs expected
      # for columns we expect a multiple of 24
      output$col_label_count <- renderText({
          count <- length(values$col_labels)
          expected <- ceiling_to_multiple(count-1, 24)
          return(paste0("Column labels ", count, "/", expected))
      })

      # Number of row labels vs expected
      # for rows we expect a multiple of 16
      output$row_label_count <- renderText({
          count <- length(values$row_labels)
          expected <- ceiling_to_multiple(count-1, 16)
          return(paste0("Row labels ", count, "/", expected))
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
          return(paste0("Files ", count_files, "/", expected))
      })
      
      # Tab - Fold change
      output$plate_fold_change <- DT::renderDataTable({ 
        if (is.null(values$plates)) return()
        df <- values$plates$fold_change

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


      # Tab - Ct / QC tab plates
      output$plates <- renderUI({
        if (is.null(values$dt_list)) return()
        # datatables are calculated elsewhere to enable dynamic rendering of
        # multiple datatables
        func <- function(x) { renderDataTable(x) }
        lapply(names(values$dt_list), function(n){ func(values$dt_list[[n]]) })
      })


      # Tab Reference table
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
          ) %>%
          formatRound(columns=2, digits=2)
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

