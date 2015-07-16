
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
# 
# shiny::runApp(appDir = "/home/htreutle/Code/Java/MetSWATH_GUI")
# 
# shiny::runGist("https://gist.github.com/pssguy/5199135")
# 
# shinyapps::setAccountInfo(name='treutler', token='9E626539E372514577A0015CD5171A60', secret='2aGCt7oyx8VA+WP45/ZDCsIqF5dFIYqa2pAlC5lI')
# shinyapps::deployApp(appDir = "/home/htreutle/Code/Java/MetSWATH_GUI")
# https://treutler.shinyapps.io/MetSWATH_GUI
# 
# htreutle@ipb-halle.de
# 
# 1 vs 2
# >= 30000
# <= -1
# 

library(shiny)

shinyUI(
  ui = navbarPage(title = "MetSWATH GUI", 
    tabPanel(
      ##########################################################################################
      ## tab run
      title = "Run",
      column(width = 4,
        wellPanel(
          helpText("This app does the following...")
        ),
        wellPanel(
          h4("File input"),
          p("Please choose a fragment matrix file"),
          fileInput(
            multiple = FALSE,
            inputId = 'matrixFile', 
            #label = 'Choose fragment matrix file',
            label = NULL,
            accept = c('text/csv', 'text/comma-separated-values,text/plain')
          )
        ),
        conditionalPanel(
          condition = "output.showGUI",
          wellPanel(
            h4("Processed file"),
            verbatimTextOutput("fileInfo"), 
            h4("Filter"),
            fluidRow(
              column(
                width = 6,
                tags$div(title="Please select the first group",
                         selectInput(#width = "100px",
                           inputId = "groupOne", 
                           label = "Group 1",
                           choices = c("")
                         )
                )
              ),
              column(
                width = 6,
                tags$div(title="Please select the second group",
                         selectInput(
                           inputId = "groupTwo", 
                           label = "Group 2", 
                           choices = c("")
                         )
                )
              )
            ),
            #selectInput(inputId = "filter", label = NULL, c("greater", "smaller")),
            textInput(inputId = "filter_average", label = "Average abundance greater"),
            textInput(inputId = "filter_lfc", label = "Log-fold-change more extreme than"),
            textInput(inputId = "filter_ms2_masses", label = "Spectrum comprises mass(es) x, y, z"),
            textInput(inputId = "filter_ms2_ppm", label = "PPM"),
            actionButton(inputId = "applyFilters", label = "Apply filters"),
            h4("Filtered precursors"),
            verbatimTextOutput("filteredPrecursors"),
            tags$head(tags$script(HTML('
              Shiny.addCustomMessageHandler("jsCode", function(message) { eval(message.code); });
            '))),
            actionButton(inputId = "drawPlots", label = "Draw plot"),
            h4("Information"),
            verbatimTextOutput("information"),
            htmlOutput(outputId = "metFragLink")
          )
        )
      ),
      column(width = 7,
        plotOutput(height = 500, 
                   outputId = "clusterPlotDendrogram", 
                   hover = "clusterPlotDendrogram_hover", 
                   click = "clusterPlotDendrogram_click"
                   #dblclick = "plot_dblclick",
                   #brush = "plot_brush"
        ),
        plotOutput(height = 75, 
                   outputId = "clusterPlotHeatmap", 
                   hover = "clusterPlotHeatmap_hover", 
                   click = "clusterPlotHeatmap_click"
        ),
        plotOutput(height = 250, 
                   outputId = "clusterPlotMS2", 
                   hover = "clusterPlotMS2_hover", 
                   click = "clusterPlotMS2_click"
        )
      ),
      column(width = 1,
        plotOutput(#height = 50, 
                   outputId = "clusterPlotHeatmapLegend", 
                   hover = "clusterPlotHeatmapLegend_hover", 
                   click = "clusterPlotHeatmapLegend_click"
        )
      )
    ),
    tabPanel(
      ##########################################################################################
      ## tab about
      title = "About",
      fluidRow(
        column(
          width = 4,
          h4("Authors"),
          p("Hendrik Treutler"),
          br(),
          p("Gerd Balcke"),
          br(),
          p("Steffen Neumann")
        ),
        column(
          width = 4,
          h4("References"),
          p("Submitted...")
        )
      )
    )
  )
)
