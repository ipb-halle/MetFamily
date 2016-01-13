# 
# shiny::runApp(appDir = "/home/htreutle/Code/Java/MetFam")
# shiny::runApp(appDir = "/home/htreutle/Code/Java/MetFam", host="0.0.0.0")
# 
# ############################ shinyapps
# htreutle@ipb-halle.de
# 
# install.packages('devtools')
# devtools::install_github('rstudio/shinyapps')
# library(shinyapps)
# shinyapps::setAccountInfo(name='treutler', token='402AF39746642C3C360E9AED3DDCFF15', secret='21kN1PMYvSvP375O4228XpqzTis75BkxEVggP0nb')
# shinyapps::deployApp('/home/htreutle/Code/Java/MetFam')
# https://treutler.shinyapps.io/MetFam
# 
# ############################ [OLD]
# devtools::install_github('rstudio/rsconnect')
# library(rsconnect)
# shinyapps::setAccountInfo(name='treutler', token='402AF39746642C3C360E9AED3DDCFF15', secret='21kN1PMYvSvP375O4228XpqzTis75BkxEVggP0nb')
# rsconnect::setAccountInfo(name='treutler', token='9E626539E372514577A0015CD5171A60', secret='2aGCt7oyx8VA+WP45/ZDCsIqF5dFIYqa2pAlC5lI')
# deployApp(appDir = "/home/htreutle/Code/Java/MetFam")
# rsconnect::configureApp("APPNAME", size="small")
# 
# shinyapps::setAccountInfo(name='treutler', token='9E626539E372514577A0015CD5171A60', secret='2aGCt7oyx8VA+WP45/ZDCsIqF5dFIYqa2pAlC5lI')
# shinyapps::deployApp(appDir = "/home/htreutle/Code/Java/MetFam")
# https://treutler.shinyapps.io/MetSWATH_GUI
# 
# ############################ deploy locally
# sudo cp /home/htreutle/Code/Java/MetFam/ClusteringMS2SpectraGUI.R /vol/R/shiny/srv/shiny-server/MetFam/
# sudo cp /home/htreutle/Code/Java/MetFam/server.R /vol/R/shiny/srv/shiny-server/MetFam/
# sudo cp /home/htreutle/Code/Java/MetFam/ui.R /vol/R/shiny/srv/shiny-server/MetFam/
# sudo cp /home/htreutle/Code/Java/MetFam/FragmentMatrixFunctions.R /vol/R/shiny/srv/shiny-server/MetFam/
# 
# ############################ debug
# options(warn=1)
# If warn is negative warnings are ignored; if it is zero they are stored and printed after the top–level function has completed; if it is one they are printed as they occur and if it is 2 (or larger) warnings are turned into errors. 
# options(error = recover)
# options(warn = 2, shiny.error = recover)
# options(warn = 2, shiny.error = browser)
# 

library(htmltools)
library(shiny)
library(shinyBS)
library(shinyjs)
library(DT)

#reactiveSvg <- function (outputId){
#  HTML(paste("<div id=\"", outputId, "\" class=\"shiny-network-output\"><svg /></div>", sep=""))
#}

shinyUI(
  ui = navbarPage(title = "MetFam", 
    ##########################################################################################
    ##########################################################################################
    ##########################################################################################
    ##########################################################################################
    ## tab run
    tabPanel(
      ##############################################################################################
      ## enable / disable actionButtons while process is running
      singleton(tags$head(HTML(
        '
          <script type="text/javascript">
            $(document).ready(
              function() {
                // disable start_proc button after a click
                Shiny.addCustomMessageHandler(
                  "disableButton", 
                  function(message) {
                    $("#" + message).attr("disabled", "true");
                  }
                );
                // Enable start_proc button when computation is finished
                Shiny.addCustomMessageHandler(
                  "enableButton", 
                  function(message) {
                    $("#" + message).removeAttr("disabled");
                  }
                );
              }
            )
          </script>
        '
      ))),
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
              title = "Input",
              #selectInput("example", "Choose an example lattice chart (from lattice package documentation):",
              #            choices = c("xyplot", "dotplot", "barchart")),
              #includeHTML("svgfiller.js"),
              #reactiveSvg(outputId = "myImage"),
              #imageOutput(outputId = "myImage"),
              #HTML("<b>hallo2</b><p>"),
              #HTML("<object data='/home/htreutle/Downloads/MetSWATH/svg/tmp2/tooltipLattice.svg' type='image/svg+xml'></object>"),
              #HTML("<img src='/home/htreutle/Downloads/MetSWATH/svg/tmp2/tooltipLattice.svg' alt=''>"),
              #HTML("<img src='img/tooltipLattice.svg' alt=''>"),
              #HTML("<div style='background: url(/home/htreutle/Downloads/MetSWATH/svg/tmp2/tooltipLattice.svg);'>hallo<p><p>holla</div>"),
              #HTML("<svg><path /home/htreutle/Downloads/MetSWATH/svg/tmp2/tooltipLattice.svg /></svg>"),
              #HTML("<iframe height='180' width='100%' src='/home/htreutle/Downloads/MetSWATH/svg/tmp2/tooltipLattice.svg' scrolling='no'><p>SVG-Grafik – hier eine Kopie als PNG</p></iframe>"),
              #HTML("<svg width='96' height='96'> <image xlink:href='/home/htreutle/Downloads/MetSWATH/svg/tmp2/tooltipLattice.svg' src='/home/htreutle/Downloads/MetSWATH/svg/tmp2/tooltipLattice.svg' width='96' height='96'/></svg>"),
              #HTML("<p><b>hallo1</b>"),
              wellPanel(
                bsTooltip(id = "fileInputSelection", title = "The user is able to load a project file or to import external data", placement = "bottom", trigger = "hover"),
                radioButtons(inputId = "fileInputSelection", label = NULL, choices = c("Load project", "Import data", "Example data"), selected = NULL, inline = FALSE),
                shiny::hr(),
                conditionalPanel(
                  condition = 'input.fileInputSelection == "Load project"',
                  h4("Project file input"),
                  p("Please choose a project file"),
                  bsTooltip(id = "matrixFile", title = "Press to choose an input file", placement = "bottom", trigger = "hover"),
                  fileInput(
                    multiple = FALSE,
                    inputId = 'matrixFile', 
                    label = NULL, #label = 'Choose fragment matrix file',
                    accept = c('text/comma-separated-values', 'text/plain', 'text/tab-separated-values')
                  ),
                  bsTooltip(id = "loadProjectData", title = "Press to load the selected project file", placement = "bottom", trigger = "hover"),
                  actionButton(inputId = "loadProjectData", label = "Load project data")
                ),
                conditionalPanel(
                  condition = 'input.fileInputSelection == "Import data"',
                  fluidRow(
                    column(width = 6,
                           div(style="float:left",
                               h4("Parameters for import")
                           )
                    ),##column
                    column(width = 6,
                           div(style="float:right",
                               bsTooltip(id = "showImportParameters", title = "Display parameters for processing during import", placement = "bottom", trigger = "hover"),
                               checkboxInput(inputId = "showImportParameters", label = "Show parameters", value = FALSE)
                           )
                    )##column
                  ),##row
                  conditionalPanel(
                    condition = "input.showImportParameters",
                    h5("Fragment filter"),
                    fluidRow(
                      column(width = 6,
                             bsTooltip(id = "minimumIntensityOfMaximalMS2peak", title = "A MS/MS spectrum is considered iff the MS/MS feature with maximum intensity is greater or equal than this value", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "minimumIntensityOfMaximalMS2peak", label = "Min. MS/MS peak intensity", value = 2000)
                      ),##column
                      column(width = 6,
                             bsTooltip(id = "minimumProportionOfMS2peaks", title = "A MS/MS feature is considered iff the intensity is greater or equal than the maximum intensity times this value", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "minimumProportionOfMS2peaks", label = "MS/MS peak proportion", value = 0.05)
                      )##column
                    ),##row
                    h5("Fragment grouping"),
                    fluidRow(
                      column(width = 6,
                             bsTooltip(id = "mzDeviationAbsolute_grouping", title = "A MS/MS feature is added to a fragment group if the absolute m/z difference is smaller or equal than this value", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "mzDeviationAbsolute_grouping", label = "m/z deviation (abs.)", value = 0.01)
                      ),##column
                      column(width = 6,
                             bsTooltip(id = "mzDeviationInPPM_grouping", title = "A MS/MS feature is added to a fragment group if the absolute m/z difference is smaller or equal than the m/z times this value divided by 1,000,000 (<b>p</b>arts <b>p</b>er <b>m</b>illion)", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "mzDeviationInPPM_grouping", label = "m/z deviation (PPM)", value = 10)
                      )##column
                    ),##row
                    
                    
                    fluidRow(
                      column(width = 6,
                             div(style="float:left",
                                 h4("Advanced parameters")
                             )
                      ),##column
                      column(width = 6,
                             div(style="float:right",
                                 bsTooltip(id = "showImportParametersAdvanced", title = "Display parameters for processing during import", placement = "bottom", trigger = "hover"),
                                 checkboxInput(inputId = "showImportParametersAdvanced", label = "Show parameters", value = FALSE)
                             )
                      )##column
                    ),##row
                    conditionalPanel(
                      condition = "input.showImportParametersAdvanced",
                      h5("MS feature deisotoping"),
                      bsTooltip(id = "doPrecursorDeisotoping", title = "If checked, the set of MS features is deisotoped", placement = "bottom", trigger = "hover"),
                      checkboxInput(inputId = "doPrecursorDeisotoping", label = "MS feature deisotoping", value = TRUE),
                      conditionalPanel(
                        condition = "input.doPrecursorDeisotoping",
                        fluidRow(
                          column(width = 6,
                                 bsTooltip(id = "mzDeviationAbsolute_precursorDeisotoping", title = "A MS feature is considered an +1 isotopic peak if the absolute of the m/z difference to the (putative) monoisotopic peak minus 1.0033548378 (=<sup>13</sup>C - <sup>12</sup>C) is smaller or equal than this value (analog for the +2 isotopic peak)", placement = "bottom", trigger = "hover"),
                                 textInput(inputId = "mzDeviationAbsolute_precursorDeisotoping", label = "m/z deviation (abs.)", value = 0.01)
                          ),##column
                          column(width = 6,
                                 bsTooltip(id = "mzDeviationInPPM_precursorDeisotoping", title = "A MS feature is considered an +1 isotopic peak if the absolute of the m/z difference to the (putative) monoisotopic peak minus 1.0033548378 (=<sup>13</sup>C - <sup>12</sup>C) is smaller or equal than the m/z times this value divided by 1,000,000 (<b>p</b>arts <b>p</b>er <b>m</b>illion, analog for the +2 isotopic peak)", placement = "bottom", trigger = "hover"),
                                 textInput(inputId = "mzDeviationInPPM_precursorDeisotoping", label = "m/z deviation (PPM)", value = 10)
                          )##column
                        ),##row
                        bsTooltip(id = "maximumRtDifference", title = "A MS feature is considered an isotopic peak if the absolute of the retention time difference to the (putative) monoisotopic peak is smaller or equal than this value (in minutes)", placement = "bottom", trigger = "hover"),
                        textInput(inputId = "maximumRtDifference", label = "Retention time difference", value = 0.02)
                      ),##conditional
                      h5("Fragment deisotoping"),
                      bsTooltip(id = "doMs2PeakGroupDeisotoping", title = "If checked, the set of MS/MS features is deisotoped", placement = "bottom", trigger = "hover"),
                      checkboxInput(inputId = "doMs2PeakGroupDeisotoping", label = "Fragment deisotoping", value = TRUE),
                      conditionalPanel(
                        condition = "input.doMs2PeakGroupDeisotoping",
                        fluidRow(
                          column(width = 6,
                                 bsTooltip(id = "mzDeviationAbsolute_ms2PeakGroupDeisotoping", title = "_A MS/MS feature is considered an +1 isotopic peak if the absolute of the m/z difference to the (putative) monoisotopic peak minus 1.0033548378 (=<sup>13</sup>C - <sup>12</sup>C) is smaller or equal than this value", placement = "bottom", trigger = "hover"),
                                 textInput(inputId = "mzDeviationAbsolute_ms2PeakGroupDeisotoping", label = "m/z deviation (abs.)", value = 0.01)
                          ),##column
                          column(width = 6,
                                 bsTooltip(id = "mzDeviationInPPM_ms2PeakGroupDeisotoping", title = "A MS/MS feature is considered an +1 isotopic peak if the absolute of the m/z difference to the (putative) monoisotopic peak minus 1.0033548378 (=<sup>13</sup>C - <sup>12</sup>C) is smaller or equal than the m/z times this value divided by 1,000,000 (<b>p</b>arts <b>p</b>er <b>m</b>illion)", placement = "bottom", trigger = "hover"),
                                 textInput(inputId = "mzDeviationInPPM_ms2PeakGroupDeisotoping", label = "m/z deviation (PPM)", value = 10)
                          )##column
                        )##row
                      )##conditional
                    )
                  ),
                  h4("Data file input"),
                  p("Please choose MS data file"),
                  bsTooltip(id = "ms1DataFile", title = "Press to choose a MS data file", placement = "bottom", trigger = "hover"),
                  fileInput(
                    multiple = FALSE,
                    inputId = 'ms1DataFile', 
                    label = NULL, #label = 'Choose fragment matrix file',
                    accept = c('text/comma-separated-values', 'text/plain', 'text/tab-separated-values')
                  ),
                  p("Please choose MS/MS data file"),
                  bsTooltip(id = "ms2DataFile", title = "Press to choose a MS/MS data file", placement = "bottom", trigger = "hover"),
                  fileInput(
                    multiple = FALSE,
                    inputId = 'ms2DataFile', 
                    label = NULL, #label = 'Choose fragment matrix file',
                    accept = c('text/plain', 'msp')
                  ),
                  bsTooltip(id = "importMs1Ms2Data", title = "Press to import the selected MS data file and MS/MS data file", placement = "bottom", trigger = "hover"),
                  actionButton(inputId = "importMs1Ms2Data", label = "Import MS and MS/MS data")
                ),## conditional
                conditionalPanel(
                  condition = 'input.fileInputSelection == "Example data"',
                  h4("Example data input"),
                  helpText(
                    "The data set used as a showcase in the MetFam publication referenced in the tab 'About'."
                  ),
                  bsTooltip(id = "exampleDataSelection", title = "The user is able to choose the full data set or a reduced data set (only MS features with MS abundance >= 5000)", placement = "bottom", trigger = "hover"),
                  radioButtons(inputId = "exampleDataSelection", label = NULL, choices = c("Example data set (full)", "Example data set (reduced)"), selected = "Example data set (reduced)", inline = FALSE),
                  bsTooltip(id = "loadExampleData", title = "Press to load the example data set", placement = "bottom", trigger = "hover"),
                  actionButton(inputId = "loadExampleData", label = "Load example data")
                )## conditional
              ),##well
              wellPanel(
                h4("Input status"),
                bsTooltip(id = "fileInfo", title = "The current input status", placement = "bottom", trigger = "hover"),
                verbatimTextOutput(outputId = "fileInfo")
              )## well panel
            ),## tab panel
            ##############################################################################################
            ##############################################################################################
            ## global MS2 filter
            tabPanel(
              title = "MS/MS filter",
              shinyjs::useShinyjs(),
              conditionalPanel(
                condition = "output.showGUI",
                wellPanel(
                  ##############################################################################################
                  ## MS2 filter
                  h4("Global MS/MS filter"),
                  ##############################################################################################
                  ## MS2 plot
                  fluidRow(
                    column(width = 6,
                           div(style="float:left",
                               #conditionalPanel(
                              #   condition = "input.showFragmentPlot",
                                 h4("Fragment overview")
                               #)
                           )
                    ),##column
                    column(width = 6,
                           div(style="float:right",
                               bsTooltip(id = "showFragmentPlot", title = "Display m/z and frequency of frequent fragments", placement = "bottom", trigger = "hover"),
                               checkboxInput(inputId = "showFragmentPlot", label = "Show frequent fragments", value = FALSE)
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
                  ),
                  fluidRow(
                    column(width = 11,
                           bsTooltip(id = "globalFilter_ms2_masses1", title = "The MS/MS spectra should include the following fragment / neutral loss mass(es) (separated by \",\")<p>E.g. \"96.969, -162.053\" for a compound with a phosphate - fragment (H<sub>2</sub>PO<sub>4</sub><sup>-</sup>) and a hexose - neutral loss (C<sub>6</sub>O<sub>5</sub>H<sub>10</sub>)", placement = "bottom", trigger = "hover"),
                           textInput(inputId = "globalFilter_ms2_masses1", label = "MS/MS spectrum includes m/z(s) #1")
                    ),
                    column(width = 1,
                           h4("or")
                    )
                  ),
                  fluidRow(
                    column(width = 11,
                           bsTooltip(id = "globalFilter_ms2_masses2", title = "The MS/MS spectra should include the following fragment / neutral loss mass(es) (separated by \",\")<p>E.g. \"96.969, -162.053\" for a compound with a phosphate - fragment (H<sub>2</sub>PO<sub>4</sub><sup>-</sup>) and a hexose - neutral loss (C<sub>6</sub>O<sub>5</sub>H<sub>10</sub>)", placement = "bottom", trigger = "hover"),
                           textInput(inputId = "globalFilter_ms2_masses2", label = "MS/MS spectrum includes m/z(s) #2")
                    ),
                    column(width = 1,
                           h4("or")
                    )
                  ),
                  fluidRow(
                    column(width = 11,
                           bsTooltip(id = "globalFilter_ms2_masses3", title = "The MS/MS spectra should include the following fragment / neutral loss mass(es) (separated by \",\")<p>E.g. \"96.969, -162.053\" for a compound with a phosphate - fragment (H<sub>2</sub>PO<sub>4</sub><sup>-</sup>) and a hexose - neutral loss (C<sub>6</sub>O<sub>5</sub>H<sub>10</sub>)", placement = "bottom", trigger = "hover"),
                           textInput(inputId = "globalFilter_ms2_masses3", label = "MS/MS spectrum includes m/z(s) #3")
                    ),
                    column(width = 1,
                           h4("")
                    )
                  ),
                  bsTooltip(id = "globalFilter_ms2_ppm", title = "The MS/MS feature matching allows this error in PPM (<b>p</b>arts <b>p</b>er <b>m</b>illion)", placement = "bottom", trigger = "hover"),
                  textInput(inputId = "globalFilter_ms2_ppm", label = "PPM"),
                  ##############################################################################################
                  ## filter button
                  bsTooltip(id = "applyGlobalMS2filters", title = "Press to determine the global set of MS features which MS/MS spectra comprise the given MS/MS features", placement = "bottom", trigger = "hover"),
                  actionButton(inputId = "applyGlobalMS2filters", label = "Apply MS/MS filter")
                ),## well
                wellPanel(
                  h4("Filtered MS features"),
                  bsTooltip(id = "globalMS2filteredPrecursors", title = "The number of MS features which MS/MS spectra comprise the given MS/MS features", placement = "bottom", trigger = "hover"),
                  verbatimTextOutput("globalMS2filteredPrecursors"),
                  conditionalPanel(
                    condition = "output.globalMS2filterValid",
                    bsTooltip(id = "downloadGlobalMS2filteredPrecursors", title = "Download a project file which is reduced to the filtered set of MS features", placement = "bottom", trigger = "hover"),
                    downloadButton('downloadGlobalMS2filteredPrecursors', 'Download reduced project file')
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
                    column(width = 7,
                           div(style="float:left",
                               h4("MS abundance filter for HCA")
                           )
                    ),##column
                    column(width = 5,
                           div(style="float:right",
                               bsTooltip(id = "showHCAfilterOptions", title = "Display filter settings", placement = "bottom", trigger = "hover"),
                               checkboxInput(inputId = "showHCAfilterOptions", label = "Show filter settings", value = TRUE)
                           )
                    )##column
                  ),##row
                  conditionalPanel(
                    condition = "input.showHCAfilterOptions",
                    bsTooltip(id = "hcaFilter_average", title = "The average MS abundance should be greater or equal than this value", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "hcaFilter_average", label = "Average MS abundance"),
                    bsTooltip(id = "hcaFilter_lfc", title = "The log<sub>2</sub>-fold change [ log<sub>2</sub>( mean(group two) / mean(group one) ) ] between the average MS abundances should be greater/smaller or equal than this value", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "hcaFilter_lfc", label = "MS log2-fold change"),
                    fluidRow(
                      column(
                        width = 6,
                        tags$div(
                          title="Please select the first replicate group",
                          radioButtons(inputId = "hcaFilterGroupOne", label = "Group 1", choices = c(""))
                          #selectInput(inputId = "groupOne", label = "Group 1", choices = c(""))
                        )
                      ),
                      column(
                        width = 6,
                        tags$div(
                          title="Please select the second replicate group",
                          radioButtons(inputId = "hcaFilterGroupTwo", label = "Group 2", choices = c(""))
                          #selectInput(inputId = "groupTwo", label = "Group 2", choices = c(""))
                        )
                      )
                    ),
                    bsTooltip(id = "hcaFilterIncludeIgnoredPrecursors", title = "Include or filter out ignored MS features, i.e. MS features which have been annotated as \\'Ignore\\'", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "hcaFilterIncludeIgnoredPrecursors", label = "Include ignored MS features", value = FALSE),
                    ##############################################################################################
                    ## filter button
                    bsTooltip(id = "applyHcaFilters", title = "Press to determine the set of MS features which fulfill the given filter criteria", placement = "bottom", trigger = "hover"),
                    actionButton(inputId = "applyHcaFilters", label = "Apply filter")
                  )## conditional panel
                ),##well panel
                wellPanel(
                  h4("Filtered MS features"),
                  bsTooltip(id = "hcaFilteredPrecursors", title = "The number of MS features which fulfill the given filter criteria", placement = "bottom", trigger = "hover"),
                  verbatimTextOutput("hcaFilteredPrecursors"),
                  conditionalPanel(
                    condition = "output.hcaFilterValid",
                    bsTooltip(id = "downloadHcaFilteredPrecursors", title = "Download a project file which is reduced to the filtered set of MS features", placement = "bottom", trigger = "hover"),
                    downloadButton('downloadHcaFilteredPrecursors', 'Download reduced project file')
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
                      selectInput(multiple = FALSE, inputId = "hcaDistanceFunction", label = "Distance function", selected = "Jaccard (intensity-weighted)", choices = c(
                        "Jaccard",
                        "Jaccard (intensity-weighted)",
                        #"Jaccard (intensity-weighted pure)",
                        #"Similarity (intensity-weighted)",
                        #"Jaccard (intensity-fragment-count-weighted)",
                        #"Similarity (intensity-fragment-count-weighted)",
                        "Jaccard (fragment-count-weighted)",
                        #"Manhatten distance",
                        "NDP (Normalized dot product)"
                      ), selectize = FALSE),
                      ## TODO remove cluster method?
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
                    bsTooltip(id = "drawHCAplots", title = "Display the HCA dendrogram given the set of filtered MS features and HCA settings", placement = "bottom", trigger = "hover"),
                    actionButton(inputId = "drawHCAplots", label = "Draw hierarchical cluster"),
                    conditionalPanel(
                      condition = "output.showGUI && output.plotHcaShown",
                      bsTooltip(id = "downloadDistanceMatrix", title = "Download the distance matrix which is the basis of the hierarchical cluster dendrogram", placement = "bottom", trigger = "hover"),
                      downloadButton('downloadDistanceMatrix', 'Download distance matrix')
                    )
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
                  h4("Please apply a valid MS/MS filter")
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
                    column(width = 7,
                           div(style="float:left",
                               h4("MS abundance filter for PCA")
                           )
                    ),##column
                    column(width = 5,
                           div(style="float:right",
                               bsTooltip(id = "showPCAfilterOptions", title = "Display filter settings", placement = "bottom", trigger = "hover"),
                               checkboxInput(inputId = "showPCAfilterOptions", label = "Show filter settings", value = TRUE)
                           )
                    )##column
                  ),##row
                  conditionalPanel(
                    condition = "input.showPCAfilterOptions",
                    bsTooltip(id = "pcaFilter_average", title = "The average MS abundance should be greater or equal than this value", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "pcaFilter_average", label = "Average MS abundance"),
                    bsTooltip(id = "pcaFilter_lfc", title = "The log<sub>2</sub>-fold change [ log<sub>2</sub>( mean(group two) / mean(group one) ) ] between the average MS abundances should be greater/smaller or equal than this value", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "pcaFilter_lfc", label = "MS log2-fold change"),
                    tags$div(title="Please select the set of replicate groups",
                             checkboxGroupInput(inputId = "pcaGroups", label = "Groups", choices = c(""))
                             #selectInput(inputId = "groups", label = "Groups", choices = c(""), multiple = TRUE, selectize = FALSE)
                    ),
                    bsTooltip(id = "pcaFilterIncludeIgnoredPrecursors", title = "Include or filter out ignored MS features, i.e. MS features which have been annotated as \\'Ignore\\'", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "pcaFilterIncludeIgnoredPrecursors", label = "Include ignored MS features", value = FALSE),
                    ##############################################################################################
                    ## filter button
                    bsTooltip(id = "applyPcaFilters", title = "Press to determine the set of MS features which fulfill the given criteria", placement = "bottom", trigger = "hover"),
                    actionButton(inputId = "applyPcaFilters", label = "Apply filter")
                  )## conditional panel
                ),##well panel
                wellPanel(
                  h4("Filtered MS features"),
                  bsTooltip(id = "pcaFilteredPrecursors", title = "The number of MS features which fulfill the given criteria", placement = "bottom", trigger = "hover"),
                  verbatimTextOutput("pcaFilteredPrecursors"),
                  conditionalPanel(
                    condition = "output.pcaFilterValid",
                    bsTooltip(id = "downloadPcaFilteredPrecursors", title = "Download a project file which is reduced to the filtered set of MS features", placement = "bottom", trigger = "hover"),
                    downloadButton('downloadPcaFilteredPrecursors', 'Download reduced project file')
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
                      bsTooltip(id = "pcaScaling", title = "Adjust the scaling of MS abundances for PCA", placement = "bottom", trigger = "hover"),
                      selectInput(multiple = FALSE, inputId = "pcaScaling", label = "Scaling", selected = "Pareto", choices = c(
                        "None", 
                        "Mean center", 
                        "Autoscaling (unit variance)",
                        "Pareto"
                        #"Vector normalization", 
                      ), selectize = FALSE),
                      bsTooltip(id = "pcaLogTransform", title = "MS abundances for PCA will be log<sub>2</sub> transformed", placement = "bottom", trigger = "hover"),
                      checkboxInput(inputId = "pcaLogTransform", label = "Log2 transform", value = FALSE),
                      fluidRow(
                        column(
                          width = 6,
                          tags$div(title="Please select the first principal component",
                            selectInput(inputId = "pcaDimensionOne", label = "Component 1", choices = c("1", "2", "3", "4", "5"), selected = "1", selectize = FALSE)
                          )
                        ),
                        column(
                          width = 6,
                          tags$div(title="Please select the second principal component",
                            selectInput(inputId = "pcaDimensionTwo", label = "Component 2", choices = c("1", "2", "3", "4", "5"), selected = "2", selectize = FALSE)
                          )
                        )
                      )
                    ),
                    bsTooltip(id = "drawPCAplots", title = "Display the PCA scores and loadings plot given the set of filtered MS features and PCA settings", placement = "bottom", trigger = "hover"),
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
                  h4("Please apply a valid MS/MS filter")
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
                  bsTooltip(id = "searchMS1orMS2", title = "Please choose the criterion for selecting MS features", placement = "bottom", trigger = "hover"),
                  radioButtons(inputId = "searchMS1orMS2", label = NULL, choices = c("MS feature m/z", "Fragment m/z")),
                  hr(),
                  conditionalPanel(
                    condition = "input.searchMS1orMS2 == 'MS feature m/z'",
                    bsTooltip(id = "searchMS1mass", title = "The MS feature m/z should be similar to this value", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "searchMS1mass", label = "MS feature m/z"),
                    bsTooltip(id = "searchMS1massPpm", title = "The specified MS feature m/z allows this error in PPM (<b>p</b>arts <b>p</b>er <b>m</b>illion)", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "searchMS1massPpm", label = "PPM"),
                    bsTooltip(id = "searchMS1includeIgnoredPrecursors", title = "Include or filter out ignored MS features, i.e. MS features which have been annotated as \\'Ignore\\'", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "searchMS1includeIgnoredPrecursors", label = "Include ignored MS features", value = FALSE),
                    bsTooltip(id = "applySearchMS1", title = "Press to mark the set of MS features which fulfill the given criteria", placement = "bottom", trigger = "hover"),
                    actionButton(inputId = "applySearchMS1", label = "Search")
                  ),## conditional panel
                  conditionalPanel(
                    condition = "input.searchMS1orMS2 == 'Fragment m/z'",
                    fluidRow(
                      column(width = 11,
                             bsTooltip(id = "search_ms2_masses1", title = "The MS/MS spectra should include the following fragment / neutral loss mass(es) (separated by \",\")<p>E.g. \"96.969, -162.053\" for a compound with a phosphate - fragment (H<sub>2</sub>PO<sub>4</sub><sup>-</sup>) and a hexose - neutral loss (C<sub>6</sub>O<sub>5</sub>H<sub>10</sub>)", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "search_ms2_masses1", label = "MS/MS spectrum includes mass(es) #1")
                      ),
                      column(width = 1,
                             h4("or")
                      )
                    ),
                    fluidRow(
                      column(width = 11,
                             bsTooltip(id = "search_ms2_masses2", title = "The MS/MS spectra should include the following fragment / neutral loss mass(es) (separated by \",\")<p>E.g. \"96.969, -162.053\" for a compound with a phosphate - fragment (H<sub>2</sub>PO<sub>4</sub><sup>-</sup>) and a hexose - neutral loss (C<sub>6</sub>O<sub>5</sub>H<sub>10</sub>)", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "search_ms2_masses2", label = "MS/MS spectrum includes mass(es) #2")
                      ),
                      column(width = 1,
                             h4("or")
                      )
                    ),
                    fluidRow(
                      column(width = 11,
                             bsTooltip(id = "search_ms2_masses3", title = "The MS/MS spectra should include the following fragment / neutral loss mass(es) (separated by \",\")<p>E.g. \"96.969, -162.053\" for a compound with a phosphate - fragment (H<sub>2</sub>PO<sub>4</sub><sup>-</sup>) and a hexose - neutral loss (C<sub>6</sub>O<sub>5</sub>H<sub>10</sub>)", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "search_ms2_masses3", label = "MS/MS spectrum includes mass(es) #3")
                      ),
                      column(width = 1,
                             h4("")
                      )
                    ),
                    bsTooltip(id = "searchMS2massPpm", title = "The specified fragment m/z's allow this error in PPM (<b>p</b>arts <b>p</b>er <b>m</b>illion)", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "searchMS2massPpm", label = "PPM"),
                    bsTooltip(id = "searchMS2includeIgnoredPrecursors", title = "Include or filter out ignored MS features, i.e. MS features which have been annotated as \\'Ignore\\'", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "searchMS2includeIgnoredPrecursors", label = "Include ignored MS features", value = FALSE),
                    bsTooltip(id = "applySearchMS2", title = "Press to mark the set of MS features which fulfill the given criteria", placement = "bottom", trigger = "hover"),
                    actionButton(inputId = "applySearchMS2", label = "Search")
                  ),## conditional panel
                  bsTooltip(id = "clearSearch", title = "Press to clear the selected set of MS feature hits in HCA and PCA", placement = "bottom", trigger = "hover"),
                  actionButton(inputId = "clearSearch", label = "Clear search")
                ),## well panel
                wellPanel(
                  h4("MS feature hits"),
                  bsTooltip(id = "searchInfo", title = "The number of MS features which fulfill the given filter criteria", placement = "bottom", trigger = "hover"),
                  verbatimTextOutput("searchInfo"),
                  conditionalPanel(
                    condition = "output.searchfilterValid & output.filterSearchActive",
                    bsTooltip(id = "downloadSearchPrecursors", title = "Download a project file which is reduced to the searched set of MS features", placement = "bottom", trigger = "hover"),
                    downloadButton('downloadSearchPrecursors', 'Download reduced project file')
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
            ),## tab panel
            tabPanel(
              title = "Export", 
              shinyjs::useShinyjs(),
              conditionalPanel(
                condition = "output.showGUI",
                wellPanel(
                  bsTooltip(id = "downloadAllPrecursors", title = "Download the full project file", placement = "bottom", trigger = "hover"),
                  downloadButton('downloadAllPrecursors', 'Download project file'),
                  conditionalPanel(
                    condition = "output.showGUI && output.showHCAplotPanel",
                    bsTooltip(id = "downloadHcaImage", title = "Download the HCA plots as png image", placement = "bottom", trigger = "hover"),
                    downloadButton('downloadHcaImage', 'Download HCA as image')
                  ),
                  conditionalPanel(
                    condition = "output.showGUI && output.showPCAplotPanel",
                    bsTooltip(id = "downloadPcaImage", title = "Download the PCA plots as png image", placement = "bottom", trigger = "hover"),
                    downloadButton('downloadPcaImage', 'Download PCA as image')
                  )
                  ## TODO annotated msp file
                  #h4("annotated msp file")
                )##well
              ),## conditional panel
              conditionalPanel(
                condition = "!output.showGUI",
                wellPanel(
                  h4("Please import a data file")
                )## well panel
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
        h4(HTML("<b>MetFam 0.9</b>")),
        fluidRow(
          column(width = 9,
            helpText(
              "This web application is designed for the identification of regulated metabolite families. This is possible on the basis of metabolite profiles for a set of MS features as well as one MS/MS spectrum for each MS feature. Group-discriminating MS features are identified using a principal component analysis (PCA) of metabolite profiles and metabolite families are identified using a hierarchical cluster analysis (HCA) of MS/MS spectra. Regulated metabolite families are identified by considering group-discriminating MS features from corporate metabolite families."
            )
          ),##column
          column(width = 3,
            HTML("<a href='http://www.ipb-halle.de/', target='_blank'><img src='logo_ipb_en.png' /></a>")
          )##column
        ),## row
        h4("Reference"),
        h5(HTML("<b>Discovering regulated metabolite families in comprehensive metabolomics studies</b>")),
        br(),
        p(HTML("Hendrik Treutler<sup>1*</sup>, Steffen Neumann<sup>1</sup>, Gerd Balcke<sup>2</sup>")), 
        p(HTML("<sup>1</sup>Leibniz Institute for Plant Biochemistry, Dept. of SEB, Weinberg 3, 06120 Halle, Germany")),
        p(HTML("<sup>2</sup>Leibniz Institute for Plant Biochemistry, Dept. of SZB, Weinberg 3, 06120 Halle, Germany")),
        p(HTML("<sup>*</sup>Corresponding author: Hendrik Treutler <a href='mailto:hendrik.treutler@ipb-halle.de?subject=MetFam%20request'>hendrik.treutler@ipb-halle.de</a>")),
        br(),
        h5(HTML("<b>Abstract</b>")),
        p(HTML("Understanding metabolism is fundamental and the identification and quantification of thousands of metabolites by mass spectrometry in modern metabolomics is a prerequisite for elucidating this area. However, the identification of metabolites is a major bottleneck in traditional approaches hampering advances. Here, we present a novel approach for the untargeted discovery of metabolite families offering a bird's eye view on metabolic regulation in comparative metabolomics.")),
        br(),
        p("Submitted to ...")
      ),## well panel
      wellPanel(
        h4("R session info"),
        verbatimTextOutput(outputId = "rInfo")
      )## well panel
    )## tab
  )## navBar
)## shinyUI
