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
    
    # other parameters
    theme = "bootstrap.css"
)

# Server logic ----
server <- function(input, output, session) {
    output$title <- renderText("Hippocompute<sup>BETA</sup>")
  
    # load tabs
    serverQPCR("qpcr")
}

# Run app ----
shinyApp(ui = ui, server = server)
