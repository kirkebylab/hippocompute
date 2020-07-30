# Load packages ----
# loaded in app.R to avoid reloading packages

# Source helper functions ----
# source("helpers_qrt_pcr.R") # automatically sourced

# Setup ----
delta_h9_list <- load_reference() # from helpers_qrt_pcr.R

tabQrtpcrUI <- function(id, label = "qrt-pcr") {
  # `NS(id)` returns a namespace function, which we save as `ns` and will
  # invoke later.
  ns <- NS(id)
  
  tabPanel("qRT-PCR",
    sidebarLayout(
     # input
     sidebarPanel(
       textInput(ns("col_labels"), label = "Column labels", value = "", placeholder = "Paste column labels ..."),
       textInput(ns("row_labels"), label = "Row labels", value = "", placeholder = "Paste row labels ..."),
       fileInput(ns("data_file"), label = "File input", multiple=FALSE, placeholder = "No file selected"),
       actionButton(ns("button_add_table"), "Add table"),
       actionButton(ns("button_clear_input"), "Clear"),
       uiOutput(ns("clip"), inline=TRUE),
     ),
     
     
     # output
     mainPanel(
       tabsetPanel(
         tabPanel("Processed Ct values", dataTableOutput(ns("plate_processed"))),
         tabPanel("Raw Ct values", dataTableOutput(ns("plate_raw"))),
         tabPanel("Debug",
            p("These boxes reflect the input values and may be used for troubleshooting."),
            # show / hide this area
            # https://stackoverflow.com/questions/44790028/show-hide-entire-box-element-in-r-shiny
            # https://stackoverflow.com/questions/51333133/using-shinyjs-to-hide-show-ui-elements
            column(4,h3("Column labels")),
            column(4,h3("Row labels")),
            column(4,h3("File input")),
            column(4, verbatimTextOutput(ns("col_labels"))),
            column(4, verbatimTextOutput(ns("row_labels"))),
            column(4, verbatimTextOutput(ns("data_file"))),
         )
       ),
     )
    )
  )
}

tabQrtpcrServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      notification_warning_nas <- NULL
      notification_warning_zeros <- NULL
      
      # Reactive events
      observeEvent(input$button_clear_input, {
        updateTextInput(session, "col_labels", value = "")
        updateTextInput(session, "row_labels", value = "")
      })
      
      plate_processed <- eventReactive(input$button_add_table, {
        validate(
          need(input$data_file, label="File input"),
          need(input$col_labels != '', label="Column labels"),
          need(input$row_labels != '', label="Row labels")
        )
        process_plate(
          in_file=input$data_file, 
          col_names=input$col_labels, 
          row_names=input$row_labels, 
          reference_data=delta_h9_list)
      })
      
      # how to use renderui in a module:
      # https://gist.github.com/tbadams45/38f1f56e0f2d7ced3507ef11b0a2fdce
      output$clip <- renderUI({
        # this creates the content of the actual button
        rclipButton(session$ns("clipbtn"), "Copy result", "placeholder", icon("clipboard"))
      })
      observeEvent(input$clipbtn, {
        clipr::write_clip(as.data.frame(plate_processed()$data_processed))
      }) # needed for clip to work
      
      
      # Output tables
      output$plate_processed <- DT::renderDataTable({ 
        df <- plate_processed()$data_processed
        df <- datatable(
          data = df,
          options = list(pageLength=25)
        ) %>%
          formatRound(1:ncol(df), 2)
        
        df })
      
      output$plate_raw <- DT::renderDataTable({ 
        result <- plate_processed()
        
        df <- result$data_raw_ct
        
        mask_zeros <- result$mask_zeros
        mask_nas <- result$mask_nas
        mask <- matrix(0L, nrow = dim(df)[1], ncol = dim(df)[2])
        mask[mask_zeros] <- 1
        mask[mask_nas] <- 2
        
        # set table cell color according to replacement, see:
        # https://stackoverflow.com/questions/50798941/r-shiny-rendertable-change-cell-colors
        # https://stackoverflow.com/questions/42569302/format-color-of-shiny-datatable-dt-according-to-values-in-a-different-dataset
        df <- datatable(
          data = cbind(df,mask), 
          options=list(
            pageLength = 25,
            columnDefs = list(list(visible=FALSE, targets=c((1+ncol(df)):(ncol(df)+ncol(mask)))))
          ),
          selection = "single"
        ) %>%
          formatStyle(
            1:ncol(df),
            valueColumns=(1+ncol(df)):(ncol(df)+ncol(mask)),
            backgroundColor=styleEqual(c(1,2), c("khaki", "lightcoral"))
          ) %>%
          formatRound(1:ncol(df), 2)
        
        df })
      
      
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
