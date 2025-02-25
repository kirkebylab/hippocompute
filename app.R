# Load packages ----
library(shiny)
# N.B. Remember to add dependencies in the DESCRIPTION file

# Source helper functions ----
# All *.R files in the R/ subdir are automatically sourced

# Setup ----
# delta_h9_list <- load_reference()

# User interface ----
ui <- navbarPage(htmlOutput("title"), windowTitle="Hippocompute", 
    # load tabs
    uiQPCR("qpcr"),
    
    tabPanel("Help",
             fluidPage(
               titlePanel("About the App"),
               br(),
               p("Welcome to Hippocompute!"),
               p("Hippocompute is designed for processing qRT-PCR data to easily and efficiently calculate mRNA expression levels.",
               "After running the qRT-PCR in the LightCycler, analyze the data in the software, and export the analyzed file in .txt format."),
               h4("How to use"),
               tags$ol(
                 tags$li("Input column and row labels from your PCR plate, separated by a tab or space."),
                 tags$li("Upload processed qRT-PCR data file."),
                 tags$li("Select the reference dataset for the cell line used."),
                 tags$li("Click 'Calculate' to process the data."),
                 tags$li("View results in the 'Ct' and 'Fold change' tabs to the right.")
               ),
               h4("Additional notes"),
               p("The labels correspond to the primers (gene names) and samples. Column labels are typically the primers, and row labels are the samples.",
               "Revisit your data to determine the correct labels. Replicates should be accounted for.",
               "Primer names must match the reference dataset exactly. Refer to the reference dataset for the correct gene names.",
               "Also, ensure that your input files are in the correct, LightCycler-processed format before uploading."),
               p("Hippocompute was developed by the Tobias Overlund Stannius as a master thesis project at the Kirkeby lab at the University of Copenhagen.", 
               "It is currently (Feb 2025) maintained by Alrik Schörling.",
               p("If you want more information, please visit the GitHub page: ",
                 tags$a(href = "https://github.com/kirkebylab/hippocompute", "Hippocompute GitHub"),
                 ". For example, it is possible to clone the repository and run the app locally.")
               )
             )
    ),
             
    
    # other parameters
    #theme = "bootstrap.css"
    
    #250214 changed the theme (alrik)
    theme = bslib::bs_theme(bootswatch = "flatly"),
    
    tags$head(
      tags$style(HTML("
        .btn-custom {
          height: 34px !important;  /* Adjust button height */
          font-size: 16px !important;  /* Ensure consistent font size */
          padding: 5px 15px !important; /* Uniform padding */
        }
      "))
    )

)

# Server logic ----
server <- function(input, output, session) {
  output$title <- renderText("Hippocompute")
  
    # load tabs
    serverQPCR("qpcr")
}

# Run app ----
shinyApp(ui = ui, server = server)
