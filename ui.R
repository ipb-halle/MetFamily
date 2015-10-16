
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
# 
# shiny::runApp(appDir = "/home/htreutle/Code/Java/MetSWATH_GUI")
# shiny::runApp(appDir = "/home/htreutle/Code/Java/MetSWATH_GUI", host="0.0.0.0")
# 
# shiny::runGist("https://gist.github.com/pssguy/5199135")
# shiny::runGitHub("ksavin/SelectableRows")
# 
# shinyapps::setAccountInfo(name='treutler', token='9E626539E372514577A0015CD5171A60', secret='2aGCt7oyx8VA+WP45/ZDCsIqF5dFIYqa2pAlC5lI')
# shinyapps::deployApp(appDir = "/home/htreutle/Code/Java/MetSWATH_GUI")
# https://treutler.shinyapps.io/MetSWATH_GUI
# 
# htreutle@ipb-halle.de
# 
# sudo cp /home/htreutle/Code/Java/MetSWATH_GUI/ClusteringMS2SpectraGUI.R /vol/R/shiny/srv/shiny-server/MetFam/
# sudo cp /home/htreutle/Code/Java/MetSWATH_GUI/server.R /vol/R/shiny/srv/shiny-server/MetFam/
# sudo cp /home/htreutle/Code/Java/MetSWATH_GUI/ui.R /vol/R/shiny/srv/shiny-server/MetFam/
# 

library(htmltools)
library(shiny)
library(shinyBS)
library(shinyjs)
library(DT)

#source("ClusteringMS2SpectraGUI.R")

shinyUI(
  ui = navbarPage(title = "MetFam", 
    ##########################################################################################
    ##########################################################################################
    ##########################################################################################
    ##########################################################################################
    ## tab run
    tabPanel(
      title = "Run",
      ##############################################################################################
      ##############################################################################################
      ##############################################################################################
      ## side panel
      conditionalPanel(
        condition = "output.showSideBar",
        column(width = 4,
          tabsetPanel(
            id = "runTabs",
            ##############################################################################################
            ##############################################################################################
            ## file input
            tabPanel(
              shinyjs::useShinyjs(),
              title = "File input",
              wellPanel(
                h4("File input"),
                p("Please choose a fragment matrix file"),
                bsTooltip(id = "matrixFile", title = "Press to choose an input file", placement = "bottom", trigger = "hover"),
                fileInput(
                  multiple = FALSE,
                  inputId = 'matrixFile', 
                  label = NULL, #label = 'Choose fragment matrix file',
                  accept = c('text/csv', 'text/comma-separated-values,text/plain')
                )
              ),## well panel
              wellPanel(
                conditionalPanel(
                  condition = "output.showGUI",
                  h4("Processed file"),
                  bsTooltip(id = "fileInfo", title = "The input file", placement = "bottom", trigger = "hover"),
                  verbatimTextOutput(outputId = "fileInfo")
                ),## conditional panel
                conditionalPanel(
                  condition = "!output.showGUI",
                  h4("Please import a data file")
                )## conditional panel
              )## well panel
            ),## tab panel
            ##############################################################################################
            ##############################################################################################
            ## global MS2 filter
            tabPanel(
              title = "MS2 filter",
              shinyjs::useShinyjs(),
              conditionalPanel(
                condition = "output.showGUI",
                wellPanel(
                  ##############################################################################################
                  ## MS2 filter
                  h4("Global MS2 filter"),
                  ##############################################################################################
                  ## MS2 plot
                  fluidRow(
                    column(width = 6,
                           div(style="float:left",
                               h4("Fragment masses")
                           )
                    ),##column
                    column(width = 6,
                           div(style="float:right",
                               bsTooltip(id = "showFragmentPlot", title = "Display fragment masses", placement = "bottom", trigger = "hover"),
                               checkboxInput(inputId = "showFragmentPlot", label = "Show abundant fragments", value = FALSE)
                           )
                    )##column
                  ),##row
                  conditionalPanel(
                    condition = "input.showFragmentPlot",
                    
                    plotOutput(height = 200, 
                               outputId = "fragmentPlot", 
                               #hover    = "fragmentPlot_hover",
                               hover    = hoverOpts(
                                 id = "fragmentPlot_hover",
                                 delay = 50, 
                                 delayType = "debounce"
                               ),
                               click    = "fragmentPlot_click",
                               dblclick = "fragmentPlot_dblclick",
                               #brush    = "fragmentPlot_brush"
                               brush    = brushOpts(
                                 id = "fragmentPlot_brush",
                                 resetOnNew = TRUE,
                                 direction = "x",
                                 delay = 00,
                                 delayType = "debounce"
                               )
                    )
                    ## TODO
                  ),
                  
                  fluidRow(
                    column(width = 11,
                           bsTooltip(id = "globalFilter_ms2_masses1", title = "The MS2 spectra should include the following fragment / neutral loss mass(es) (separated by \",\")<p>E.g. 96.969, -162.053 for a compound with a phosphate - fragment (H<sub>2</sub>PO<sub>4</sub>-) and a hexose - neutral loss (C<sub>6</sub>O<sub>5</sub>H<sub>10</sub>)", placement = "bottom", trigger = "hover"),
                           textInput(inputId = "globalFilter_ms2_masses1", label = "MS2 spectrum includes mass(es) #1")
                    ),
                    column(width = 1,
                           h4("or")
                    )
                  ),
                  fluidRow(
                    column(width = 11,
                           bsTooltip(id = "globalFilter_ms2_masses2", title = "The MS2 spectra should include the following fragment / neutral loss mass(es) (separated by \",\")<p>E.g. 96.969, -162.053 for a compound with a phosphate - fragment (H<sub>2</sub>PO<sub>4</sub>-) and a hexose - neutral loss (C<sub>6</sub>O<sub>5</sub>H<sub>10</sub>)", placement = "bottom", trigger = "hover"),
                           textInput(inputId = "globalFilter_ms2_masses2", label = "MS2 spectrum includes mass(es) #2")
                    ),
                    column(width = 1,
                           h4("or")
                    )
                  ),
                  fluidRow(
                    column(width = 11,
                           bsTooltip(id = "globalFilter_ms2_masses3", title = "The MS2 spectra should include the following fragment / neutral loss mass(es) (separated by \",\")<p>E.g. 96.969, -162.053 for a compound with a phosphate - fragment (H<sub>2</sub>PO<sub>4</sub>-) and a hexose - neutral loss (C<sub>6</sub>O<sub>5</sub>H<sub>10</sub>)", placement = "bottom", trigger = "hover"),
                           textInput(inputId = "globalFilter_ms2_masses3", label = "MS2 spectrum includes mass(es) #3")
                    ),
                    column(width = 1,
                           h4("")
                    )
                  ),
                  bsTooltip(id = "globalFilter_ms2_ppm", title = "The MS2 spectra fragment matching allows this error in PPM", placement = "bottom", trigger = "hover"),
                  textInput(inputId = "globalFilter_ms2_ppm", label = "PPM"),
                  ##############################################################################################
                  ## filter button
                  bsTooltip(id = "applyGlobalMS2filters", title = "Press to determine the global set of precursors which fulfill the given filter criteria", placement = "bottom", trigger = "hover"),
                  actionButton(inputId = "applyGlobalMS2filters", label = "Apply MS2 filter")
                ),## well
                wellPanel(
                  h4("Filtered precursors"),
                  bsTooltip(id = "globalMS2filteredPrecursors", title = "The number of precursors which fulfill the given MS2 filter criteria", placement = "bottom", trigger = "hover"),
                  verbatimTextOutput("globalMS2filteredPrecursors"),
                  conditionalPanel(
                    condition = "output.globalMS2filterValid",
                    bsTooltip(id = "downloadGlobalMS2filteredPrecursors", title = "Download the set of precursors which fulfil the given criteria", placement = "bottom", trigger = "hover"),
                    downloadButton('downloadGlobalMS2filteredPrecursors', 'Download filtered precursors')
                  )
                )## well
              ),## conditional panel
              conditionalPanel(
                condition = "!output.showGUI",
                wellPanel(
                  h4("Please import a data file")
                )
              )## conditional panel
            ),## tab panel
            ##############################################################################################
            ##############################################################################################
            ## HCA
            tabPanel(
              title = "HCA", 
              shinyjs::useShinyjs(),
              conditionalPanel(
                condition = "output.showGUI && output.globalMS2filterValid",
                wellPanel(
                  ##############################################################################################
                  ## HCA group and abundance filter
                  fluidRow(
                    column(width = 6,
                           div(style="float:left",
                               h4("MS1 Abundance filter")
                           )
                    ),##column
                    column(width = 6,
                           div(style="float:right",
                               bsTooltip(id = "showHCAfilterOptions", title = "Display filter settings", placement = "bottom", trigger = "hover"),
                               checkboxInput(inputId = "showHCAfilterOptions", label = "Show filter settings", value = TRUE)
                           )
                    )##column
                  ),##row
                  conditionalPanel(
                    condition = "input.showHCAfilterOptions",
                    fluidRow(
                      column(
                        width = 6,
                        tags$div(
                          title="Please select the first group",
                          radioButtons(inputId = "hcaFilterGroupOne", label = "Group 1", choices = c(""))
                          #selectInput(inputId = "groupOne", label = "Group 1", choices = c(""))
                        )
                      ),
                      column(
                        width = 6,
                        tags$div(
                          title="Please select the second group",
                          radioButtons(inputId = "hcaFilterGroupTwo", label = "Group 2", choices = c(""))
                          #selectInput(inputId = "groupTwo", label = "Group 2", choices = c(""))
                        )
                      )
                    ),
                    bsTooltip(id = "hcaFilter_average", title = "The average MS1 abundance should be greater than", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "hcaFilter_average", label = "MS1 average abundance"),
                    bsTooltip(id = "hcaFilter_lfc", title = "The Log2-fold-change [ log_2( mean(group two) / mean(group one) ) ] between the average MS1 group abundances should be greater/smaller or equal than", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "hcaFilter_lfc", label = "MS1 log2-fold change"),
                    bsTooltip(id = "hcaFilterIncludeIgnoredPrecursors", title = "Include or filter ignored precursors, i.e. precursors which have been annotated as 'ignored'", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "hcaFilterIncludeIgnoredPrecursors", label = "Include ignored precursors", value = FALSE),
                    ##############################################################################################
                    ## filter button
                    bsTooltip(id = "applyHcaFilters", title = "Press to determine the global set of precursors which fulfill the given filter criteria", placement = "bottom", trigger = "hover"),
                    actionButton(inputId = "applyHcaFilters", label = "Apply filter")
                  )## conditional panel
                ),##well panel
                wellPanel(
                  h4("Filtered precursors"),
                  bsTooltip(id = "hcaFilteredPrecursors", title = "The number of precursors which fulfill the given MS2 filter criteria", placement = "bottom", trigger = "hover"),
                  verbatimTextOutput("hcaFilteredPrecursors"),
                  conditionalPanel(
                    condition = "output.hcaFilterValid",
                    bsTooltip(id = "downloadHcaFilteredPrecursors", title = "Download the set of precursors which fulfil the given criteria", placement = "bottom", trigger = "hover"),
                    downloadButton('downloadHcaFilteredPrecursors', 'Download filtered precursors')
                  )##conditional
                ),##well panel
                conditionalPanel(
                  condition = "output.hcaFilterValid",
                  wellPanel(
                    ##############################################################################################
                    ## HCA properties
                    fluidRow(
                      column(width = 6,
                             div(style="float:left",
                                 h4("HCA properties")
                             )
                      ),##column
                      column(width = 6,
                             div(style="float:right",
                                 bsTooltip(id = "showHCAadvancedOptions", title = "Display further HCA settings", placement = "bottom", trigger = "hover"),
                                 checkboxInput(inputId = "showHCAadvancedOptions", label = "Show advanced options", value = FALSE)
                             )
                      )##column
                    ),##row
                    conditionalPanel(
                      condition = "input.showHCAadvancedOptions",
                      bsTooltip(id = "hcaDistanceFunction", title = "The distance function used for clustering", placement = "bottom", trigger = "hover"),
                      selectInput(multiple = FALSE, inputId = "hcaDistanceFunction", label = "Distance function", selected = "Jaccard", choices = c(
                        "Jaccard",
                        "Jaccard (intensity-weighted)",
                        "Jaccard (intensity-weighted map)",
                        "Similarity (intensity-weighted)",
                        "Jaccard (intensity-fragment-count-weighted)",
                        "Similarity (intensity-fragment-count-weighted)",
                        "Jaccard (fragment-count-weighted)",
                        "Manhatten",
                        "NDP"
                      ), selectize = FALSE),
                      bsTooltip(id = "hcaClusterMethod", title = "The method used for clustering", placement = "bottom", trigger = "hover"),
                      selectInput(multiple = FALSE, inputId = "hcaClusterMethod", label = "Cluster method", selected = "ward.D", choices = c(
                        "single", 
                        "complete", 
                        "average", 
                        "mcquitty", 
                        "median", 
                        "centroid", 
                        "ward.D", 
                        "ward.D2"
                      ), selectize = FALSE)
                    ),
                    bsTooltip(id = "drawHCAplots", title = "Display the HCA dendrogram given the set of filtered precursors and HCA settings", placement = "bottom", trigger = "hover"),
                    actionButton(inputId = "drawHCAplots", label = "Draw hierarchical cluster")
                  )## well panel
                )## conditional panel
              ),## conditional panel
              conditionalPanel(
                 condition = "!output.showGUI",
                wellPanel(
                  h4("Please import a data file")
                )## well panel
              ),## conditional panel
              conditionalPanel(
                condition = "output.showGUI && !output.globalMS2filterValid",
                wellPanel(
                  h4("Please apply a valid MS2 filter")
                )## well panel
              )## conditional panel
            ),## tab panel
            ##############################################################################################
            ##############################################################################################
            ## PCA
            tabPanel(
              title = "PCA", 
              shinyjs::useShinyjs(),
              conditionalPanel(
                condition = "output.showGUI && output.globalMS2filterValid",
                wellPanel(
                  ##############################################################################################
                  ## HCA group and abundance filter
                  fluidRow(
                    column(width = 6,
                           div(style="float:left",
                               h4("MS1 abundance filter")
                           )
                    ),##column
                    column(width = 6,
                           div(style="float:right",
                               bsTooltip(id = "showPCAfilterOptions", title = "Display filter settings", placement = "bottom", trigger = "hover"),
                               checkboxInput(inputId = "showPCAfilterOptions", label = "Show filter settings", value = TRUE)
                           )
                    )##column
                  ),##row
                  conditionalPanel(
                    condition = "input.showPCAfilterOptions",
                    tags$div(title="Please select the set of groups",
                             checkboxGroupInput(inputId = "pcaGroups", label = "Groups", choices = c(""))
                             #selectInput(inputId = "groups", label = "Groups", choices = c(""), multiple = TRUE, selectize = FALSE)
                    ),
                    bsTooltip(id = "pcaFilter_average", title = "The average MS1 abundance should be greater than", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "pcaFilter_average", label = "MS1 average abundance"),
                    bsTooltip(id = "pcaFilter_lfc", title = "The Log2-fold-change [ log_2( mean(group two) / mean(group one) ) ] between the average MS1 group abundances should be greater/smaller or equal than", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "pcaFilter_lfc", label = "MS1 log2-fold change"),
                    bsTooltip(id = "pcaFilterIncludeIgnoredPrecursors", title = "Include or filter ignored precursors, i.e. precursors which have been annotated as 'ignored'", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "pcaFilterIncludeIgnoredPrecursors", label = "Include ignored precursors", value = FALSE),
                    ##############################################################################################
                    ## filter button
                    bsTooltip(id = "applyPcaFilters", title = "Press to determine the global set of precursors which fulfill the given filter criteria", placement = "bottom", trigger = "hover"),
                    actionButton(inputId = "applyPcaFilters", label = "Apply filter")
                  )## conditional panel
                ),##well panel
                wellPanel(
                  h4("Filtered precursors"),
                  bsTooltip(id = "pcaFilteredPrecursors", title = "The number of precursors which fulfill the given MS2 filter criteria", placement = "bottom", trigger = "hover"),
                  verbatimTextOutput("pcaFilteredPrecursors"),
                  conditionalPanel(
                    condition = "output.pcaFilterValid",
                    bsTooltip(id = "downloadPcaFilteredPrecursors", title = "Download the set of precursors which fulfil the given criteria", placement = "bottom", trigger = "hover"),
                    downloadButton('downloadPcaFilteredPrecursors', 'Download filtered precursors')
                  )##conditional
                ),##well
                conditionalPanel(
                  condition = "output.pcaFilterValid",
                  wellPanel(
                    ##############################################################################################
                    ## PCA properties
                    fluidRow(
                      column(width = 6,
                             div(style="float:left",
                                 h4("PCA properties")
                             )
                      ),##column
                      column(width = 6,
                             div(style="float:right",
                                 bsTooltip(id = "showPCAadvancedOptions", title = "Display further PCA settings", placement = "bottom", trigger = "hover"),
                                 checkboxInput(inputId = "showPCAadvancedOptions", label = "Show advanced options", value = FALSE)
                             )
                      )##column
                    ),##row
                    conditionalPanel(
                      condition = "input.showPCAadvancedOptions",
                      bsTooltip(id = "pcaScaling", title = "Adjust the scaling of MS1 abundances for PCA", placement = "bottom", trigger = "hover"),
                      selectInput(multiple = FALSE, inputId = "pcaScaling", label = "Scaling", selected = "Pareto", choices = c(
                        "None", 
                        "Mean center", 
                        "Autoscaling (unit variance)",
                        "Pareto"
                        #"Vector normalization", 
                      ), selectize = FALSE),
                      bsTooltip(id = "pcaLogTransform", title = "MS1 abundances for PCA will be log2 transformed", placement = "bottom", trigger = "hover"),
                      checkboxInput(inputId = "pcaLogTransform", label = "Log2 transform", value = FALSE),
                      fluidRow(
                        column(
                          width = 6,
                          tags$div(title="Please select the first dimension",
                            selectInput(inputId = "pcaDimensionOne", label = "Dimension 1", choices = c("1", "2", "3", "4", "5"), selected = "1", selectize = FALSE)
                          )
                        ),
                        column(
                          width = 6,
                          tags$div(title="Please select the second dimension",
                            selectInput(inputId = "pcaDimensionTwo", label = "Dimension 2", choices = c("1", "2", "3", "4", "5"), selected = "2", selectize = FALSE)
                          )
                        )
                      )
                    ),
                    bsTooltip(id = "drawPCAplots", title = "Display the PCA scores and loadings plot given the set of filtered precursors and PCA settings", placement = "bottom", trigger = "hover"),
                    actionButton(inputId = "drawPCAplots", label = "Draw principal components")
                  )##well
                )## conditiojal panel
              ),## conditiojal panel
              conditionalPanel(
                condition = "!output.showGUI",
                wellPanel(
                  h4("Please import a data file")
                )## well panel
              ),## conditional panel
              conditionalPanel(
                condition = "output.showGUI && !output.globalMS2filterValid",
                wellPanel(
                  h4("Please apply a valid MS2 filter")
                )## well panel
              )## conditional panel
            ),## tab panel
            ##############################################################################################
            ##############################################################################################
            ## search
            tabPanel(
              title = "Search",
              shinyjs::useShinyjs(),
              conditionalPanel(
                condition = "output.showGUI && (output.plotHcaShown || output.plotPcaShown)",
                #condition = "output.showGUI",
                wellPanel(
                  h4("Search mode"),
                  radioButtons(inputId = "searchMS1orMS2", label = NULL, choices = c("Precursor mass", "Fragment mass")),
                  conditionalPanel(
                    condition = "input.searchMS1orMS2 == 'Precursor mass'",
                    bsTooltip(id = "searchMS1mass", title = "The average MS1 mass should be similar to the specified mass", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "searchMS1mass", label = "Precursor mass"),
                    bsTooltip(id = "searchMS1massPpm", title = "The specified precursor mass allows this error in PPM", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "searchMS1massPpm", label = "PPM"),
                    bsTooltip(id = "searchMS1includeIgnoredPrecursors", title = "Include or filter ignored precursors, i.e. precursors which have been annotated as 'ignored'", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "searchMS1includeIgnoredPrecursors", label = "Include ignored precursors", value = FALSE),
                    bsTooltip(id = "applySearchMS1", title = "Press to mark the set of precursors which fulfill the given criteria", placement = "bottom", trigger = "hover"),
                    actionButton(inputId = "applySearchMS1", label = "Search")
                  ),## conditional panel
                  conditionalPanel(
                    condition = "input.searchMS1orMS2 == 'Fragment mass'",
                    fluidRow(
                      column(width = 11,
                             bsTooltip(id = "search_ms2_masses1", title = "The MS2 spectra should include the following fragment / neutral loss mass(es) (separated by \",\")<p>E.g. 96.969, -162.053 for a compound with a phosphate - fragment (H<sub>2</sub>PO<sub>4</sub>-) and a hexose - neutral loss (C<sub>6</sub>O<sub>5</sub>H<sub>10</sub>)", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "search_ms2_masses1", label = "MS2 spectrum includes mass(es) #1")
                      ),
                      column(width = 1,
                             h4("or")
                      )
                    ),
                    fluidRow(
                      column(width = 11,
                             bsTooltip(id = "search_ms2_masses2", title = "The MS2 spectra should include the following fragment / neutral loss mass(es) (separated by \",\")<p>E.g. 96.969, -162.053 for a compound with a phosphate - fragment (H<sub>2</sub>PO<sub>4</sub>-) and a hexose - neutral loss (C<sub>6</sub>O<sub>5</sub>H<sub>10</sub>)", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "search_ms2_masses2", label = "MS2 spectrum includes mass(es) #2")
                      ),
                      column(width = 1,
                             h4("or")
                      )
                    ),
                    fluidRow(
                      column(width = 11,
                             bsTooltip(id = "search_ms2_masses3", title = "The MS2 spectra should include the following fragment / neutral loss mass(es) (separated by \",\")<p>E.g. 96.969, -162.053 for a compound with a phosphate - fragment (H<sub>2</sub>PO<sub>4</sub>-) and a hexose - neutral loss (C<sub>6</sub>O<sub>5</sub>H<sub>10</sub>)", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "search_ms2_masses3", label = "MS2 spectrum includes mass(es) #3")
                      ),
                      column(width = 1,
                             h4("")
                      )
                    ),
                    bsTooltip(id = "searchMS2massPpm", title = "The specified precursor mass allows this error in PPM", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "searchMS2massPpm", label = "PPM"),
                    bsTooltip(id = "searchMS2includeIgnoredPrecursors", title = "Include or filter ignored precursors, i.e. precursors which have been annotated as 'ignored'", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "searchMS2includeIgnoredPrecursors", label = "Include ignored precursors", value = FALSE),
                    bsTooltip(id = "applySearchMS2", title = "Press to mark the set of precursors which fulfill the given criteria", placement = "bottom", trigger = "hover"),
                    actionButton(inputId = "applySearchMS2", label = "Search")
                  )## conditional panel
                ),## well panel
                wellPanel(
                  h4("Precursor hits"),
                  bsTooltip(id = "searchInfo", title = "The number of precursors which fulfill the given MS2 filter criteria", placement = "bottom", trigger = "hover"),
                  verbatimTextOutput("searchInfo"),
                  conditionalPanel(
                    condition = "output.searchfilterValid",
                    bsTooltip(id = "downloadSearchPrecursors", title = "Download the set of precursors which fulfil the given criteria", placement = "bottom", trigger = "hover"),
                    downloadButton('downloadSearchPrecursors', 'Download precursors hits')
                  )##conditional
                )##well
              ),## conditional panel
              conditionalPanel(
                #condition = "!(output.analysisType == 'HCA' || output.analysisType == 'PCA')",
                condition = "!output.showGUI",
                wellPanel(
                  h4("Please import a data file")
                )## well
              ),## conditional panel
              conditionalPanel(
                condition = "output.showGUI && !(output.plotHcaShown || output.plotPcaShown)",
                wellPanel(
                  h4("Please search HCA or PCA")
                )## well
              )## conditional panel
            )## tab panel
          )## tab set panel
        )## column
      ),##conditional
      ##############################################################################################
      ##############################################################################################
      ##############################################################################################
      ## plots
      uiOutput("runRightColumn")
    ),## tab
    ##########################################################################################
    ##########################################################################################
    ##########################################################################################
    ##########################################################################################
    ## tab about
    tabPanel(
      title = "About",
      ##############################################################################################
      ## intro
      wellPanel(
        helpText("This app does the following...")
      ),## well panel
      wellPanel(
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
      ),## well panel
      wellPanel(
        h4("R session info"),
        verbatimTextOutput(outputId = "rInfo")
      )## well panel
    )## tab
  )## navBar
)## shinyUI
