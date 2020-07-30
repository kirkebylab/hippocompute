# Load packages ----
library(shiny)

# Source helper functions ----
# All *.R files in the R/ subdir are automatically sourced

# Setup ----
# delta_h9_list <- load_reference()

# User interface ----
ui <- navbarPage("Hippocompute",
    # load submodules
    tabQrtpcrUI("qrtpcr")
)

# Server logic ----
server <- function(input, output, session) {
    # load submodules
    tabQrtpcrServer("qrtpcr")
}

# Run app ----
shinyApp(ui = ui, server = server)
