
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
# 

library(shiny)
library(shinyBS)
library(shinyjs)

shinyUI(
  ui = navbarPage(title = "MetSWATH GUI", 
    ##########################################################################################
    ##########################################################################################
    ## tab run
    tabPanel(
      title = "Run",
      column(width = 4,
        tabsetPanel(
          id = "runTabs",
          ##############################################################################################
          ## file input
          tabPanel(
            shinyjs::useShinyjs(),
            title = "File input",
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
            ),## well panel
            ##############################################################################################
            ## file info
            wellPanel(
              conditionalPanel(
                condition = "output.showGUI",
                h4("Processed file"),
                bsTooltip(id = "fileInfo", title = "The input file", placement = "bottom", trigger = "hover"),
                verbatimTextOutput("fileInfo")
              ),## conditional panel
              conditionalPanel(
                condition = "!output.showGUI",
                h4("Please import a data file")
              )## conditional panel
            )## well panel
          ),## tab panel
          tabPanel(
            title = "Filter",
            shinyjs::useShinyjs(),
            wellPanel(
              conditionalPanel(
                condition = "output.showGUI",
                ##############################################################################################
                ## filter
                h4("Untargeted"),
                #conditionalPanel(
                #  condition = 'input.analysisType == "HCA"',
                  fluidRow(
                    column(
                      width = 6,
                      tags$div(title="Please select the first group",
                        selectInput(inputId = "groupOne", label = "Group 1", choices = c(""))
                      )
                    ),
                    column(
                      width = 6,
                      tags$div(title="Please select the second group",
                        selectInput(inputId = "groupTwo", label = "Group 2", choices = c(""))
                      )
                    )
                  ),
                #),## conditional
                #conditionalPanel(
                #  condition = 'input.analysisType == "PCA"',
                #  tags$div(title="Please select the set of groups",
                #           selectInput(inputId = "groups", label = "Groups", choices = c(""), multiple = TRUE, selectize = FALSE)
                #  )
                #),## conditional
                
                bsTooltip(id = "filter_average", title = "The average MS1 abundance should be greater than", placement = "bottom", trigger = "hover"),
                textInput(inputId = "filter_average", label = "Average abundance greater"),
                bsTooltip(id = "filter_lfc", title = "The Log-fold-change [ log_2( mean(group two) / mean(group one) ) ] between the average MS1 group abundances should be greater/smaller or equal than", placement = "bottom", trigger = "hover"),
                textInput(inputId = "filter_lfc", label = "Log-fold-change more extreme than"),
                
                h4("Semi-targeted"),
                fluidRow(
                  column(width = 2,
                         h4("")
                  ),
                  column(width = 10,
                         bsTooltip(id = "filter_ms2_masses",  title = "The MS2 spectra should include the following fragment / neutral loss mass(es) (separated by \",\")<p>E.g. 96.969, -162.053 for a compound with a phosphate - fragment (H<sub>2</sub>PO<sub>4</sub>-) and a hexose - neutral loss (C<sub>6</sub>O<sub>5</sub>H<sub>10</sub>)", placement = "bottom", trigger = "hover"),
                         textInput(inputId = "filter_ms2_masses", label = "Spectrum includes mass(es)")
                  )
                ),
                fluidRow(
                  column(width = 2,
                         h4("or")
                  ),
                  column(width = 10,
                         bsTooltip(id = "filter_ms2_masses2", title = "The MS2 spectra should include the following fragment / neutral loss mass(es) (separated by \",\")<p>E.g. 96.969, -162.053 for a compound with a phosphate - fragment (H<sub>2</sub>PO<sub>4</sub>-) and a hexose - neutral loss (C<sub>6</sub>O<sub>5</sub>H<sub>10</sub>)", placement = "bottom", trigger = "hover"),
                         textInput(inputId = "filter_ms2_masses2", label = "Spectrum includes mass(es)")
                  )
                ),
                fluidRow(
                  column(width = 2,
                         h4("or")
                  ),
                  column(width = 10,
                         bsTooltip(id = "filter_ms2_masses3", title = "The MS2 spectra should include the following fragment / neutral loss mass(es) (separated by \",\")<p>E.g. 96.969, -162.053 for a compound with a phosphate - fragment (H<sub>2</sub>PO<sub>4</sub>-) and a hexose - neutral loss (C<sub>6</sub>O<sub>5</sub>H<sub>10</sub>)", placement = "bottom", trigger = "hover"),
                         textInput(inputId = "filter_ms2_masses3", label = "Spectrum includes mass(es)")
                  )
                ),
                
                bsTooltip(id = "filter_ms2_ppm", title = "The MS2 spectra fragment matching allows this error in PPM", placement = "bottom", trigger = "hover"),
                textInput(inputId = "filter_ms2_ppm", label = "PPM"),
                ##############################################################################################
                ## filter buttons
                bsTooltip(id = "applyFilters", title = "Press to determine the set of precursors which fulfill the given filter criteria", placement = "bottom", trigger = "hover"),
                actionButton(inputId = "applyFilters", label = "Apply filters"),
                h4("Filtered precursors"),
                bsTooltip(id = "filteredPrecursors", title = "The number of precursors which fulfill the given filter criteria", placement = "bottom", trigger = "hover"),
                verbatimTextOutput("filteredPrecursors"),
                downloadButton('downloadFilteredPrecursors', 'Download filtered precursors')
              ),## consitionalPanel
              conditionalPanel(
                condition = "!output.showGUI",
                h4("Please import a data file")
              )## conditional panel
            )## well panel
          ),## tab panel
          ##############################################################################################
          ## HCA
          tabPanel(
            title = "HCA", 
            shinyjs::useShinyjs(),
            wellPanel(
              conditionalPanel(
                condition = "output.showGUI && output.filterValid",
                bsTooltip(id = "distanceFunction", title = "The distance function used for clustering", placement = "bottom", trigger = "hover"),
                selectInput(multiple = FALSE, inputId = "distanceFunction", label = "Distance function", choices = c(
                  "Jaccard",
                  "Jaccard (intensity-weighted)",
                  "Jaccard (intensity-weighted map)",
                  "Similarity (intensity-weighted)",
                  "Jaccard (intensity-fragment-count-weighted)",
                  "Similarity (intensity-fragment-count-weighted)",
                  "Jaccard (fragment-count-weighted)",
                  "Manhatten",
                  "NDP"
                )),
                bsTooltip(id = "clusterMethod", title = "The method used for clustering", placement = "bottom", trigger = "hover"),
                selectInput(multiple = FALSE, inputId = "clusterMethod", label = "Cluster method", choices = c(
                  "single", 
                  "complete", 
                  "average", 
                  "mcquitty", 
                  "median", 
                  "centroid", 
                  "ward.D", 
                  "ward.D2"
                )),
                bsTooltip(id = "drawPlots", title = "Press to calculate the distances between all filtered precursors and plot the result", placement = "bottom", trigger = "hover"),
                tags$head(tags$script(HTML('
                  Shiny.addCustomMessageHandler("jsCode", function(message) { eval(message.code); });
                '))),
                actionButton(inputId = "drawHCAplots", label = "Draw hierarchical cluster")
              ),## conditional panel
              conditionalPanel(
                condition = "!output.showGUI",
                h4("Please import a data file")
              ),## conditional panel
              conditionalPanel(
                condition = "output.showGUI && !output.filterValid",
                h4("Please apply a filter")
              )## conditional panel
            )## well panel
          ),## tab panel
          ##############################################################################################
          ## PCA
          tabPanel(
            title = "PCA", 
            shinyjs::useShinyjs(),
            wellPanel(
              conditionalPanel(
                condition = "output.showGUI && output.filterValid",
                checkboxInput(inputId = "pcaUnitVariance", label = "Unit variance", value = TRUE),
                checkboxInput(inputId = "pcaLogTransform", label = "Log transform", value = TRUE),
                fluidRow(
                  column(
                    width = 6,
                    tags$div(title="Please select the first dimension",
                      selectInput(inputId = "pcaDimensionOne", label = "Dimension 1", choices = c("1", "2", "3", "4", "5"), selected = "1")
                    )
                  ),
                  column(
                    width = 6,
                    tags$div(title="Please select the second dimension",
                      selectInput(inputId = "pcaDimensionTwo", label = "Dimension 2", choices = c("1", "2", "3", "4", "5"), selected = "2")
                    )
                  )
                ),
                actionButton(inputId = "drawPCAplots", label = "Draw principal components")
              ),## conditiojal panel
              conditionalPanel(
                condition = "!output.showGUI",
                h4("Please import a data file")
              ),## conditional panel
              conditionalPanel(
                condition = "output.showGUI && !output.filterValid",
                h4("Please apply a filter")
              )## conditional panel
            )## well panel
          ),## tab panel
          ##############################################################################################
          ## search
          tabPanel(
            title = "Search",
            shinyjs::useShinyjs(),
            wellPanel(
              conditionalPanel(
                condition = "output.showGUI",
                h4("TODO search")
              ),## conditional panel
              conditionalPanel(
                condition = "!output.showGUI",
                h4("Please import a data file")
              )## conditional panel
            )## well panel
          )## tab panel
        )## tab set panel
      ),## column
      column(width = 8,
        conditionalPanel(
          condition = 'output.analysisType == "HCA" && output.showHCAplotPanel',
          fluidRow(
            ##############################################################################################
            ## plots
            column(width = 11,
              plotOutput(height = 500, 
                        outputId = "clusterPlotDendrogram", 
                        hover    = "clusterPlotDendrogram_hover", 
                        click    = "clusterPlotDendrogram_click",
                        dblclick = "clusterPlotDendrogram_dblclick",
                        #brush    = "clusterPlotDendrogram_brush"
                        brush = brushOpts(
                          id = "clusterPlotDendrogram_brush",
                          resetOnNew = TRUE,
                          direction = "x",
                          delay = 00,
                          delayType = "debounce"
                        )
              ),
              plotOutput(height = 75, 
                        outputId = "clusterPlotHeatmap",
                        hover    = "clusterPlotHeatmap_hover", 
                        #click = "clusterPlotHeatmap_click"
              )
            ),## column
            column(width = 1,
              plotOutput(outputId = "clusterPlotHeatmapLegend")
            )## column
          )## row
        ),## conditional
        conditionalPanel(
          condition = 'output.analysisType == "PCA" && output.showPCAplotPanel',
          fluidRow(
            ##############################################################################################
            ## plots
            column(width = 6,
              plotOutput(height = 500, 
                        outputId = "pcaPlotScores", 
                        hover    = "pcaPlotScores_hover",
                        click    = "pcaPlotScores_click",
                        dblclick = "pcaPlotScores_dblclick",
                        #brush    = "pcaPlotScores_brush"
                        brush = brushOpts(
                          id = "pcaPlotScores_brush",
                          resetOnNew = TRUE,
                          direction = "xy",
                          delay = 00,
                          delayType = "debounce"
                        )
              )
            ),## column
            column(width = 6,
              plotOutput(height = 500, 
                        outputId = "pcaPlotLoadings", 
                        hover    = "pcaPlotLoadings_hover",
                        click    = "pcaPlotLoadings_click",
                        dblclick = "pcaPlotLoadings_dblclick",
                        #brush    = "pcaPlotScores_brush"
                        brush = brushOpts(
                          id = "pcaPlotLoadings_brush",
                          resetOnNew = TRUE,
                          direction = "xy",
                          delay = 00,
                          delayType = "debounce"
                        )
              )
            )## column
          )## row
        ),## conditional
        conditionalPanel(
          condition = '(output.analysisType == "HCA" && output.showHCAplotPanel) || (output.analysisType == "PCA" && output.showPCAplotPanel)',
          plotOutput(height = 250, 
                    outputId = "clusterPlotMS2",
                    hover    = "clusterPlotMS2_hover",
                    click    = "clusterPlotMS2_click",
                    dblclick = "clusterPlotMS2_dblclick",
                    #brush    = "clusterPlotMS2_brush",
                    brush = brushOpts(
                      id = "clusterPlotMS2_brush",
                      resetOnNew = TRUE,
                      direction = "x",
                      delay = 00,
                      delayType = "debounce"
                    )
          ),
          fluidRow(
            ##############################################################################################
            ## infos
            h4("Information"),
            bsTooltip(id = "information", title = "Information about selected items in the plot", placement = "bottom", trigger = "hover"),
            verbatimTextOutput("information"),
            htmlOutput(outputId = "metFragLink"),
            h4("Tip"),
            bsTooltip(id = "tip", title = "Information about operating options", placement = "bottom", trigger = "hover"),
            verbatimTextOutput("tip"),
            downloadButton('downloadMatrix', 'Download selected precursors')
          )## row
        )## conditional
      )## column
    ),## tab
    ##########################################################################################
    ##########################################################################################
    ## tab about
    tabPanel(
      title = "About",
      ##############################################################################################
      ## intro
      wellPanel(
        helpText("This app does the following...")
      ),
      fluidRow(
        column(
          width = 4,
          h4("Authors"),
          p("Hendrik Treutler"),
          br(),
          p("Gerd Balcke"),
          br(),
          p("Steffen Neumann")
        ),## column
        column(
          width = 4,
          h4("References"),
          p("Submitted to...")
        )## column
      )## row
    )## tab
  )## navBar
)## shinyUI
