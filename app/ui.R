
library(htmltools)
library(shiny)
library(shinyBS)
library(shinyjs)
library(DT)

importParameterSetInit <- list(
  minimumIntensityOfMaximalMS2peak = 2000,
  minimumProportionOfMS2peaks = 0.05,
  neutralLossesPrecursorToFragments = TRUE,
  neutralLossesFragmentsToFragments = FALSE,
  mzDeviationAbsolute_grouping = 0.01,
  mzDeviationInPPM_grouping = 10,
  showImportParametersAdvanced = FALSE,
  doPrecursorDeisotoping = TRUE,
  mzDeviationAbsolute_precursorDeisotoping = 0.01,
  mzDeviationInPPM_precursorDeisotoping = 10,
  maximumRtDifference = 0.05,
  doMs2PeakGroupDeisotoping = TRUE,
  mzDeviationAbsolute_ms2PeakGroupDeisotoping = 0.01,
  mzDeviationInPPM_ms2PeakGroupDeisotoping = 10
)

shinyUI(
  ui = navbarPage(title = "MetFamily", 
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
              
              ## TODO 999
              #DT::dataTableOutput("dataTableTest"), ## TODO remove
              
              wellPanel(
                #div(style="display:inline-block",actionButton(inputId = "test", label = "", icon = icon(name = "chevron-up", lib = "font-awesome"))),
                #actionButton(inputId = "test", label = NULL, icon = icon(name = "chevron-up", lib = "font-awesome")),#, width = NULL, ...)
                #downloadButton('downloadReport2', 'Export analysis report'),
                bsTooltip(id = "fileInputSelection", title = "The user is able to load a project file or to import external data", placement = "bottom", trigger = "hover"),
                radioButtons(inputId = "fileInputSelection", label = NULL, choices = c("Import data", "Load project", "Example data"), selected = "Load project", inline = FALSE),
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
                  fluidRow(
                    column(width = 6,
                           bsTooltip(id = "loadProjectData", title = "Press to load the selected project file", placement = "bottom", trigger = "hover"),
                           actionButton(inputId = "loadProjectData", label = "Load project data", class="btn-success", width = "100%")
                    ),##column
                    column(width = 6
                           
                    )##column
                  )##row
                ),
                conditionalPanel(
                  condition = 'input.fileInputSelection == "Import data"',
                  bsTooltip(id = "projectName", title = "Please type the name of the project", placement = "bottom", trigger = "hover"),
                  textInput(inputId = "projectName", label = "Project name", value = paste("MetFamily project (created ", gsub(" ", "_", gsub(":", ".", Sys.time())), ")", sep = "")),
                  bsTooltip(id = "projectDescription", title = "Please type a description of this project as free text", placement = "bottom", trigger = "hover"),
                  tags$style(type="text/css", "textarea {width:100%}"),
                  tags$textarea(id = 'projectDescription', placeholder = 'Comments here', rows = 3, ""),
                  #textInput(inputId = "projectDescription", label = "Project description", value = "", placeholder = "Comments here"),
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
                    fileInput(
                      multiple = FALSE,
                      inputId = 'importParameterFileInput', 
                      label = 'Apply parameters from import parameter file',
                      accept = c('text/comma-separated-values', 'text/plain', 'text/tab-separated-values')
                    ),
                    h5("Spectrum filter"),
                    fluidRow(
                      column(width = 6,
                             bsTooltip(id = "minimumIntensityOfMaximalMS2peak", title = "A MS/MS spectrum is considered iff the MS/MS feature with maximum intensity is greater or equal than this value", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "minimumIntensityOfMaximalMS2peak", label = "Min. spectrum intensity", value = importParameterSetInit$minimumIntensityOfMaximalMS2peak)
                      ),##column
                      column(width = 6,
                             bsTooltip(id = "minimumProportionOfMS2peaks", title = "A MS/MS feature is considered iff the intensity is greater or equal than the maximum intensity times this value", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "minimumProportionOfMS2peaks", label = "MS/MS peak proportion", value = importParameterSetInit$minimumProportionOfMS2peaks)
                      )##column
                    ),##row
                    h5("Neutral losses"),
                    fluidRow(
                      column(width = 6,
                             bsTooltip(id = "neutralLossesPrecursorToFragments", title = "Include neutral losses relative to the precursor ion, i.e. the m/z difference between the m/z of the precursor ion and the m/z of each fragment ion of the corresponding MS/MS spectrum", placement = "bottom", trigger = "hover"),
                             checkboxInput(inputId = "neutralLossesPrecursorToFragments", label = "Fragment vs. precursor", value = importParameterSetInit$neutralLossesPrecursorToFragments)
                      ),##column
                      column(width = 6,
                             bsTooltip(id = "neutralLossesFragmentsToFragments", title = "Include neutral losses amongst fragment ions, i.e. the m/z difference between the m/z of all pairs of fragment ions within each MS/MS spectrum; this involves the incorporation of potentially nonexistent neutral losses and needs more time for processing", placement = "bottom", trigger = "hover"),
                             checkboxInput(inputId = "neutralLossesFragmentsToFragments", label = "Fragment vs. fragment", value = importParameterSetInit$neutralLossesFragmentsToFragments)
                      )##column
                    ),##row
                    h5("Fragment grouping"),
                    fluidRow(
                      column(width = 6,
                             bsTooltip(id = "mzDeviationAbsolute_grouping", title = "A MS/MS feature is added to a fragment group if the absolute m/z difference is smaller or equal than this value", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "mzDeviationAbsolute_grouping", label = "m/z deviation (abs.)", value = importParameterSetInit$mzDeviationAbsolute_grouping)
                      ),##column
                      column(width = 6,
                             bsTooltip(id = "mzDeviationInPPM_grouping", title = "A MS/MS feature is added to a fragment group if the absolute m/z difference is smaller or equal than the m/z times this value divided by 1,000,000 (<b>p</b>arts <b>p</b>er <b>m</b>illion)", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "mzDeviationInPPM_grouping", label = "m/z deviation (PPM)", value = importParameterSetInit$mzDeviationInPPM_grouping)
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
                      h5("MS\u00B9 feature deisotoping"),
                      bsTooltip(id = "doPrecursorDeisotoping", title = "If checked, the set of MS\u00B9 features is deisotoped", placement = "bottom", trigger = "hover"),
                      checkboxInput(inputId = "doPrecursorDeisotoping", label = "MS\u00B9 feature deisotoping", value = importParameterSetInit$doPrecursorDeisotoping),
                      conditionalPanel(
                        condition = "input.doPrecursorDeisotoping",
                        fluidRow(
                          column(width = 6,
                                 bsTooltip(id = "mzDeviationAbsolute_precursorDeisotoping", title = "A MS\u00B9 feature is considered an +1 isotopic peak if the absolute of the m/z difference to the (putative) monoisotopic peak minus 1.0033548378 (=<sup>13</sup>C - <sup>12</sup>C) is smaller or equal than this value (analog for the +2 isotopic peak)", placement = "bottom", trigger = "hover"),
                                 textInput(inputId = "mzDeviationAbsolute_precursorDeisotoping", label = "m/z deviation (abs.)", value = importParameterSetInit$mzDeviationAbsolute_precursorDeisotoping)
                          ),##column
                          column(width = 6,
                                 bsTooltip(id = "mzDeviationInPPM_precursorDeisotoping", title = "A MS\u00B9 feature is considered an +1 isotopic peak if the absolute of the m/z difference to the (putative) monoisotopic peak minus 1.0033548378 (=<sup>13</sup>C - <sup>12</sup>C) is smaller or equal than the m/z times this value divided by 1,000,000 (<b>p</b>arts <b>p</b>er <b>m</b>illion, analog for the +2 isotopic peak)", placement = "bottom", trigger = "hover"),
                                 textInput(inputId = "mzDeviationInPPM_precursorDeisotoping", label = "m/z deviation (PPM)", value = importParameterSetInit$mzDeviationInPPM_precursorDeisotoping)
                          )##column
                        )##row
                      ),##conditional
                      bsTooltip(id = "maximumRtDifference", title = "A MS\u00B9 feature is considered an isotopic peak if the absolute of the retention time difference to the (putative) monoisotopic peak is smaller or equal than this value (in minutes)", placement = "bottom", trigger = "hover"),
                      textInput(inputId = "maximumRtDifference", label = "Retention time difference", value = importParameterSetInit$maximumRtDifference),
                      h5("Fragment deisotoping"),
                      bsTooltip(id = "doMs2PeakGroupDeisotoping", title = "If checked, the set of MS/MS features is deisotoped", placement = "bottom", trigger = "hover"),
                      checkboxInput(inputId = "doMs2PeakGroupDeisotoping", label = "Fragment deisotoping", value = importParameterSetInit$doMs2PeakGroupDeisotoping),
                      conditionalPanel(
                        condition = "input.doMs2PeakGroupDeisotoping",
                        fluidRow(
                          column(width = 6,
                                 bsTooltip(id = "mzDeviationAbsolute_ms2PeakGroupDeisotoping", title = "_A MS/MS feature is considered an +1 isotopic peak if the absolute of the m/z difference to the (putative) monoisotopic peak minus 1.0033548378 (=<sup>13</sup>C - <sup>12</sup>C) is smaller or equal than this value", placement = "bottom", trigger = "hover"),
                                 textInput(inputId = "mzDeviationAbsolute_ms2PeakGroupDeisotoping", label = "m/z deviation (abs.)", value = importParameterSetInit$mzDeviationAbsolute_ms2PeakGroupDeisotoping)
                          ),##column
                          column(width = 6,
                                 bsTooltip(id = "mzDeviationInPPM_ms2PeakGroupDeisotoping", title = "A MS/MS feature is considered an +1 isotopic peak if the absolute of the m/z difference to the (putative) monoisotopic peak minus 1.0033548378 (=<sup>13</sup>C - <sup>12</sup>C) is smaller or equal than the m/z times this value divided by 1,000,000 (<b>p</b>arts <b>p</b>er <b>m</b>illion)", placement = "bottom", trigger = "hover"),
                                 textInput(inputId = "mzDeviationInPPM_ms2PeakGroupDeisotoping", label = "m/z deviation (PPM)", value = importParameterSetInit$mzDeviationInPPM_ms2PeakGroupDeisotoping)
                          )##column
                        )##row
                      )##conditional
                    )
                  ),
                  h4("Data file input"),
                  p("Please choose metabolite profile"),
                  bsTooltip(id = "ms1DataFile", title = "Press to choose a metabolite profile", placement = "bottom", trigger = "hover"),
                  fileInput(
                    multiple = FALSE,
                    inputId = 'ms1DataFile', 
                    label = NULL, #label = 'Choose fragment matrix file',
                    accept = c('text/comma-separated-values', 'text/plain', 'text/tab-separated-values')
                  ),
                  p("Please choose MS/MS library"),
                  bsTooltip(id = "ms2DataFile", title = "Press to choose a MS/MS library", placement = "bottom", trigger = "hover"),
                  fileInput(
                    multiple = FALSE,
                    inputId = 'ms2DataFile', 
                    label = NULL, #label = 'Choose fragment matrix file',
                    accept = c('text/plain', 'msp')
                  ),
                  fluidRow(
                    column(width = 6,
                           bsTooltip(id = "importMs1Ms2Data", title = "Press to import the selected metabolite profile and MS/MS library", placement = "bottom", trigger = "hover"),
                           actionButton(inputId = "importMs1Ms2Data", label = "Import MS\u00B9 and MS/MS data", class="btn-success", width = "100%")
                    ),##column
                    column(width = 6,
                           bsTooltip(id = "importMs2Data", title = "Press to import the selected MS/MS library without a metabolite profile", placement = "bottom", trigger = "hover"),
                           actionButton(inputId = "importMs2Data", label = "Import MS/MS data", class="btn-success", width = "100%")
                    )##column
                  )##row
                ),## conditional
                conditionalPanel(
                  condition = 'input.fileInputSelection == "Example data"',
                  h4("Example data input"),
                  helpText(
                    "The data set used as showcase in the MetFamily publication referenced in the tab 'About'."
                  ),
                  br(),
                  h4("Download original metabolite profile and MS/MS library"),
                  fluidRow(
                    column(width = 6, style="width:50%",
                           div(style="float:left;width:100%",
                               bsTooltip(id = "downloadMsData", title = "Download the original metabolite profile used in the MetFamily publication", placement = "bottom", trigger = "hover"),
                               downloadButton(outputId = "downloadMsData", label = "Download metabolite profile"),
                               tags$style(type='text/css', "#downloadMsData { width:100%}")
                           )
                    ),##column
                    column(width = 6, style="width:50%",
                           div(style="float:right;width:100%",
                               bsTooltip(id = "downloadMsMsData", title = "Download the original MS/MS library used in the MetFamily publication", placement = "bottom", trigger = "hover"),
                               downloadButton(outputId = "downloadMsMsData", label = "Download MS/MS library"),
                               tags$style(type='text/css', "#downloadMsMsData { width:100%}")
                           )
                    )##column
                  ),##row
                  br(),
                  h4("Download generated fragment matrix"),
                  fluidRow(
                    column(width = 6, style="width:50%",
                           bsTooltip(id = "downloadFragmentMatrix", title = "Download the fragment matrix generated from the original metabolite profile and MS/MS library used in the MetFamily publication", placement = "bottom", trigger = "hover"),
                           downloadButton('downloadFragmentMatrix', 'Download fragment matrix'),
                           tags$style(type='text/css', "#downloadFragmentMatrix { width:100%}")
                    ),
                    column(width = 6, style="width:50%"
                           
                    )
                  ),
                  br(),
                  fluidRow(
                    column(width = 12,
                           h4("Download showcase protocol"),
                           bsTooltip(id = "downloadDocShowcaseProtocol", title = "Download the protocol which is the basis of the results of the showcase in the MetFamily publication", placement = "bottom", trigger = "hover"),
                           downloadButton('downloadDocShowcaseProtocol', 'Download showcase protocol')
                    )
                  ),
                  br(),
                  h4("Load example data"),
                  fluidRow(
                    column(width = 6,
                           bsTooltip(id = "loadExampleData", title = "Press to load the example data set", placement = "bottom", trigger = "hover"),
                           actionButton(inputId = "loadExampleData", label = "Load example data", class="btn-success", width = "100%")
                    ),##column
                    column(width = 6
                           
                    )##column
                  )##row
                )## conditional
              ),##well
              wellPanel(
                h4("Input status"),
                bsTooltip(id = "fileInfo", title = "The current input status", placement = "bottom", trigger = "hover"),
                verbatimTextOutput(outputId = "fileInfo")
              ),## well panel
              uiOutput("errorPopupDialog")
            ),## tab panel
            ##############################################################################################
            ##############################################################################################
            ## global MS2 filter
            #navbarMenu("Filter",
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
                               h4("Fragment overview")
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
                  fluidRow(
                    column(width = 6,
                           div(style="float:left;width:100%",
                               bsTooltip(id = "applyGlobalMS2filters", title = "Press to determine the global set of MS\u00B9 features which MS/MS spectra comprise the given MS/MS features", placement = "bottom", trigger = "hover"),
                               actionButton(inputId = "applyGlobalMS2filters", label = "Apply MS/MS filter", class="btn-success", width = "100%")
                           )
                    ),##column
                    column(width = 6,
                           div(style="float:right;width:100%",
                               bsTooltip(id = "clearGlobalMS2filters", title = "Press to clear the global MS/MS filter", placement = "bottom", trigger = "hover"),
                               actionButton(inputId = "clearGlobalMS2filters", label = "Clear MS/MS filter", class="btn-success", width = "100%")
                           )
                    )##column
                  ),##row
                  hr(),
                  h4("Filtered MS\u00B9 features"),
                  bsTooltip(id = "globalMS2filteredPrecursors", title = "The number of MS\u00B9 features which MS/MS spectra comprise the given MS/MS features", placement = "bottom", trigger = "hover"),
                  verbatimTextOutput("globalMS2filteredPrecursors"),
                  conditionalPanel(
                    condition = "output.globalMS2filterValid",
                    bsTooltip(id = "downloadGlobalMS2filteredPrecursors", title = "Download a project file which is reduced to the filtered set of MS\u00B9 features", placement = "bottom", trigger = "hover"),
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
            ## sample selection
            tabPanel(
              title = "Sample filter", 
              shinyjs::useShinyjs(),
              conditionalPanel(
                condition = "output.showGUI",
                wellPanel(
                  ##############################################################################################
                  ## sample table
                  DT::dataTableOutput("sampleTable"),
                  bsTooltip(id = "updateSampleTable", title = "Updates the order and the exclusion status of samples", placement = "bottom", trigger = "hover"),
                  actionButton(inputId = "updateSampleTable", label = "Apply sample order and exclusion status", class="btn-success")
                )##well
              ),## conditiojal panel
              conditionalPanel(
                condition = "!output.showGUI",
                wellPanel(
                  h4("Please import a data file")
                )## well panel
              )## conditional panel
            ),## tab panel
            #),
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
                               h4("MS\u00B9 abundance filter for PCA")
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
                    bsTooltip(id = "pcaFilter_average", title = "The average MS\u00B9 abundance should be greater or equal than this value", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "pcaFilter_average", label = "Average MS\u00B9 abundance"),
                    bsTooltip(id = "pcaFilter_lfc", title = "The log<sub>2</sub>-fold change [ log<sub>2</sub>( mean(group one) / mean(group two) ) ] between the average MS\u00B9 abundances should be greater/smaller or equal than this value", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "pcaFilter_lfc", label = "MS\u00B9 log2-fold change"),
                    
                    fluidRow(style="overflow-y:scroll; 
                                    max-height: 200px;",
                      column(width = 7,
                                 div(style="float:left;
                                            width:100%;",
                                     tags$div(title="Please select the set of sample groups",
                                              checkboxGroupInput(inputId = "pcaGroups", label = "Groups", choices = c(""))
                                     )
                                 )
                      ),##column
                      column(width = 5,
                        fluidRow(style = "vertical-align: top",
                          bsTooltip(id = "filterByPCAgroupSamples", title = "Display and select samples for PCA analysis", placement = "bottom", trigger = "hover"),
                          checkboxInput(inputId = "filterByPCAgroupSamples", label = "Select samples", value = FALSE)
                        ),
                        fluidRow(style = "vertical-align: top",
                          bsTooltip(id = "selectAllPCAGroups", title = "Select all groups", placement = "bottom", trigger = "hover"),
                          actionButton(inputId = "selectAllPCAGroups", label = "Select all", width = "100%"),
                          bsTooltip(id = "selectNoPCAGroups", title = "Deselect all groups", placement = "bottom", trigger = "hover"),
                          actionButton(inputId = "selectNoPCAGroups", label = "Deselect all", width = "100%"),
                          bsTooltip(id = "selectInvertedPCAGroups", title = "Invert the group selection", placement = "bottom", trigger = "hover"),
                          actionButton(inputId = "selectInvertedPCAGroups", label = "Invert", width = "100%")
                        )
                      )##column
                    ),##row
                    conditionalPanel(
                      condition = "input.filterByPCAgroupSamples",
                      fluidRow(style="overflow-y:scroll; 
                                    max-height: 200px;",
                               column(width = 12,
                               div(style="float:left;
                                            width:100%;",
                        tags$div(title="Please select the set of samples",
                                 checkboxGroupInput(inputId = "pcaSamples", label = "Samples", choices = c(""))
                        )))
                      )
                    ),
                    
                    #tags$div(title="Please select the set of sample groups",
                    #         checkboxGroupInput(inputId = "pcaGroups", label = "Groups", choices = c(""))
                    #),
                    
                    bsTooltip(id = "pcaFilterIncludeIgnoredPrecursors", title = "Include or filter out ignored MS\u00B9 features, i.e. MS\u00B9 features which have been annotated as \\'Ignore\\'", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "pcaFilterIncludeIgnoredPrecursors", label = "Include ignored MS\u00B9 features", value = FALSE),
                    ##############################################################################################
                    ## filter button
                    fluidRow(
                      column(width = 6,
                             div(style="float:left;width:100%",
                                 bsTooltip(id = "applyPcaFilters", title = "Press to determine the set of MS\u00B9 features which fulfill the given criteria", placement = "bottom", trigger = "hover"),
                                 actionButton(inputId = "applyPcaFilters", label = "Apply filter", class="btn-success", width = "100%")
                             )
                      ),##column
                      column(width = 6,
                             div(style="float:right;width:100%",
                                 bsTooltip(id = "clearPcaFilters", title = "Press to clear the MS\u00B9 abundance filter for PCA", placement = "bottom", trigger = "hover"),
                                 actionButton(inputId = "clearPcaFilters", label = "Clear filter", class="btn-success", width = "100%")
                             )
                      )##column
                    )##row
                  ),## conditional panel
                  hr(),
                  h4("Filtered MS\u00B9 features"),
                  bsTooltip(id = "pcaFilteredPrecursors", title = "The number of MS\u00B9 features which fulfill the given criteria", placement = "bottom", trigger = "hover"),
                  verbatimTextOutput("pcaFilteredPrecursors"),
                  conditionalPanel(
                    condition = "output.pcaFilterValid",
                    bsTooltip(id = "downloadPcaFilteredPrecursors", title = "Download a project file which is reduced to the filtered set of MS\u00B9 features", placement = "bottom", trigger = "hover"),
                    downloadButton('downloadPcaFilteredPrecursors', 'Download reduced project file')
                  )##conditional
                ),##well
                conditionalPanel(
                  condition = "output.pcaFilterValid",
                  wellPanel(
                    bsTooltip(id = "ms1AnalysisMethod", title = "Please choose the method for the analysis of MS\u00B9 abundances", placement = "bottom", trigger = "hover"),
                    selectInput(multiple = FALSE, inputId = "ms1AnalysisMethod", label = "Method", selected = "PCA", choices = c(
                      "PCA (Principal Component Analysis)", 
                      "sPCA (Sparse Principal Component Analysis)",
                      "PLS-DA (Partial Least Squares Discriminant Analysis)", 
                      "sPLS-DA (Sparse Partial Least Squares Discriminant Analysis)"
                    ), selectize = FALSE),
                    ##############################################################################################
                    ## PCA properties
                    fluidRow(
                      column(width = 6,
                             div(style="float:left",
                                 h4("Method properties")
                             )
                      ),##column
                      column(width = 6,
                             div(style="float:right",
                                 bsTooltip(id = "showPCAadvancedOptions", title = "Display further settings for the selected method", placement = "bottom", trigger = "hover"),
                                 checkboxInput(inputId = "showPCAadvancedOptions", label = "Show advanced options", value = FALSE)
                             )
                      )##column
                    ),##row
                    conditionalPanel(
                      condition = "input.showPCAadvancedOptions",
                      bsTooltip(id = "pcaScaling", title = "Adjust the scaling of MS\u00B9 abundances for the selected method", placement = "bottom", trigger = "hover"),
                      selectInput(multiple = FALSE, inputId = "pcaScaling", label = "Scaling", selected = "Pareto", choices = c(
                        "None", 
                        "Mean center", 
                        "Autoscaling (unit variance)",
                        "Pareto"
                        #"Vector normalization", 
                      ), selectize = FALSE),
                      bsTooltip(id = "pcaLogTransform", title = "MS\u00B9 abundances for the selected method will be log<sub>2</sub> transformed", placement = "bottom", trigger = "hover"),
                      checkboxInput(inputId = "pcaLogTransform", label = "Log2 transformation", value = FALSE),
                      fluidRow(
                        column(
                          width = 6,
                          tags$div(title="Please select the first component",
                            selectInput(inputId = "pcaDimensionOne", label = "Component 1", choices = c("1", "2", "3", "4", "5"), selected = "1", selectize = FALSE)
                          )
                        ),
                        column(
                          width = 6,
                          tags$div(title="Please select the second component",
                            selectInput(inputId = "pcaDimensionTwo", label = "Component 2", choices = c("1", "2", "3", "4", "5"), selected = "2", selectize = FALSE)
                          )
                        )
                      )
                    ),
                    fluidRow(
                      column(width = 6,
                             bsTooltip(id = "drawPCAplots", title = "Display the scores and the loadings plot given the set of filtered MS\u00B9 features and settings of the selected method", placement = "bottom", trigger = "hover"),
                             actionButton(inputId = "drawPCAplots", label = "Perform analysis", class="btn-success", width = "100%")
                      ),##column
                      column(width = 6
                             ## nothing
                      )##column
                    )##row
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
                               h4("MS\u00B9 abundance filter for HCA")
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
                    bsTooltip(id = "hcaFilter_average", title = "The average MS\u00B9 abundance should be greater or equal than this value", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "hcaFilter_average", label = "Average MS\u00B9 abundance"),
                    bsTooltip(id = "hcaFilter_lfc", title = "The log<sub>2</sub>-fold change [ log<sub>2</sub>( mean(group one) / mean(group two) ) ] between the average MS\u00B9 abundances should be greater/smaller or equal than this value", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "hcaFilter_lfc", label = "MS\u00B9 log2-fold change"),
                    fluidRow(style = "overflow-y:scroll; 
                                      max-height: 200px",
                             #border: 1px solid #cccccc;
                      column(
                        width = 6,
                        tags$div(
                          title="Please select the first sample group",
                          radioButtons(inputId = "hcaFilterGroupOne", label = "Group 1", choices = c(""))
                        )
                      ),
                      column(
                        width = 6,
                        tags$div(
                          title="Please select the second sample group",
                          radioButtons(inputId = "hcaFilterGroupTwo", label = "Group 2", choices = c(""))
                        )
                      )
                    ),
                    bsTooltip(id = "hcaFilterIncludeIgnoredPrecursors", title = "Include or filter out ignored MS\u00B9 features, i.e. MS\u00B9 features which have been annotated as \\'Ignore\\'", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "hcaFilterIncludeIgnoredPrecursors", label = "Include ignored MS\u00B9 features", value = FALSE),
                    ##############################################################################################
                    ## filter button
                    fluidRow(
                      column(width = 6, 
                             div(style="float:left;width:100%",
                                 bsTooltip(id = "applyHcaFilters", title = "Press to determine the set of MS\u00B9 features which fulfill the given filter criteria", placement = "bottom", trigger = "hover"),
                                 actionButton(inputId = "applyHcaFilters", label = "Apply filter", class="btn-success", width = "100%")
                             )
                      ),##column
                      column(width = 6, 
                             div(style="float:right;width:100%",
                                 bsTooltip(id = "clearHcaFilters", title = "Press to clear the MS\u00B9 abundance filter for HCA", placement = "bottom", trigger = "hover"),
                                 actionButton(inputId = "clearHcaFilters", label = "Clear filter", class="btn-success", width = "100%")
                             )
                      )##column
                    )##row
                  ),## conditional panel
                # ),##well panel
                # wellPanel(
                  hr(),
                  h4("Filtered MS\u00B9 features"),
                  bsTooltip(id = "hcaFilteredPrecursors", title = "The number of MS\u00B9 features which fulfill the given filter criteria", placement = "bottom", trigger = "hover"),
                  verbatimTextOutput("hcaFilteredPrecursors"),
                  conditionalPanel(
                    condition = "output.hcaFilterValid",
                    bsTooltip(id = "downloadHcaFilteredPrecursors", title = "Download a project file which is reduced to the filtered set of MS\u00B9 features", placement = "bottom", trigger = "hover"),
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
                        "Jaccard (fragment-count-weighted)",
                        "NDP (Normalized dot product)"
                      ), selectize = FALSE)
                      #bsTooltip(id = "hcaClusterMethod", title = "The method used for clustering", placement = "bottom", trigger = "hover"),
                      #selectInput(multiple = FALSE, inputId = "hcaClusterMethod", label = "Cluster method", selected = "ward.D", choices = c(
                      #  "single", 
                      #  "complete", 
                      #  "average", 
                      #  "mcquitty", 
                      #  "median", 
                      #  "centroid", 
                      #  "ward.D", 
                      #  "ward.D2"
                      #), selectize = FALSE)
                    ),
                    fluidRow(
                      column(width = 6,
                             bsTooltip(id = "drawHCAplots", title = "Display the HCA dendrogram given the set of filtered MS\u00B9 features and HCA settings", placement = "bottom", trigger = "hover"),
                             actionButton(inputId = "drawHCAplots", label = "Draw hierarchical cluster", class="btn-success", width = "100%")
                      ),##column
                      column(width = 6
                             
                      )##column
                    ),##row
                    conditionalPanel(
                      condition = "output.showGUI && output.plotHcaShown",
                      br(),
                      bsTooltip(id = "downloadDistanceMatrix", title = "Download the distance matrix of the currently displayed hierarchical cluster dendrogram", placement = "bottom", trigger = "hover"),
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
            ## search
            tabPanel(
              title = "Search",
              shinyjs::useShinyjs(),
              conditionalPanel(
                condition = "output.showGUI && (output.plotHcaShown || output.plotPcaShown)",
                #condition = "output.showGUI",
                wellPanel(
                  h4("Search mode"),
                  bsTooltip(id = "searchMS1orMS2", title = "Please choose the criterion for selecting MS\u00B9 features", placement = "bottom", trigger = "hover"),
                  radioButtons(inputId = "searchMS1orMS2", label = NULL, choices = c("MS1 feature m/z", "Fragment m/z")),
                  hr(),
                  conditionalPanel(
                    condition = "input.searchMS1orMS2 == 'MS1 feature m/z'",
                    fluidRow(
                      column(width = 6,
                             bsTooltip(id = "searchMS1mass", title = "The MS\u00B9 feature m/z should be similar to at least one of the given values (separated by \",\")", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "searchMS1mass", label = "MS\u00B9 feature m/z('s)")
                      ),##column
                      column(width = 6,
                             bsTooltip(id = "searchMS1massPpm", title = "The specified MS\u00B9 feature m/z allows this error in PPM (<b>p</b>arts <b>p</b>er <b>m</b>illion)", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "searchMS1massPpm", label = "PPM")
                      )##column
                    )##row
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
                    textInput(inputId = "searchMS2massPpm", label = "PPM")
                  ),## conditional panel
                  bsTooltip(id = "searchIncludeIgnoredPrecursors", title = "Include or filter out ignored MS\u00B9 features, i.e. MS\u00B9 features which have been annotated as \\'Ignore\\'", placement = "bottom", trigger = "hover"),
                  checkboxInput(inputId = "searchIncludeIgnoredPrecursors", label = "Include ignored MS\u00B9 features", value = FALSE),
                  fluidRow(
                    column(width = 6,
                           div(style="float:left;width:100%",
                               bsTooltip(id = "applySearch", title = "Press to mark the set of MS\u00B9 features which fulfill the given criteria", placement = "bottom", trigger = "hover"),
                               actionButton(inputId = "applySearch", label = "Search", class="btn-success", width = "100%")
                           )
                    ),##column
                    column(width = 6,
                           div(style="float:right;width:100%",
                               bsTooltip(id = "clearSearch", title = "Press to clear the selected set of MS\u00B9 feature hits in HCA and PCA", placement = "bottom", trigger = "hover"),
                               actionButton(inputId = "clearSearch", label = "Clear search", class="btn-success", width = "100%")
                           )
                    )##column
                  ),##row
                  hr(),
                  h4("MS\u00B9 feature hits"),
                  bsTooltip(id = "searchInfo", title = "The number of MS\u00B9 features which fulfill the given filter criteria", placement = "bottom", trigger = "hover"),
                  verbatimTextOutput("searchInfo"),
                  conditionalPanel(
                    condition = "output.searchfilterValid & output.filterSearchActive",
                    bsTooltip(id = "downloadSearchPrecursors", title = "Download a project file which is reduced to the searched set of MS\u00B9 features", placement = "bottom", trigger = "hover"),
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
                  h4("Please perform HCA or PCA prior to search")
                )## well
              )## conditional panel
            ),## tab panel
            tabPanel(
              title = "Classifiers",
              #shinyjs::useShinyjs(),
              conditionalPanel(
                condition = "output.showGUI",
                wellPanel(
                  h4("Classifier selection"),
                  bsTooltip(id = "classifierCount", title = "The number of available classifiers for the semi-automated metabolite family annotation of MS\u00B9 features", placement = "bottom", trigger = "hover"),
                  verbatimTextOutput(outputId = "classifierCount"),
                  DT::dataTableOutput("classifierSelectionTable"),
                  bsTooltip(id = "doAnnotation", title = "Press to automatically annotate metabolite families to spectra ", placement = "bottom", trigger = "hover"),
                  actionButton(inputId = "doAnnotation", label = "Perform scan", class="btn-success", width = "100%")
                )
              ),
              conditionalPanel(
                condition = "!output.showGUI",
                wellPanel(
                  h4("Please import a data file")
                )## well panel
              )## conditional panel
            ),## tab panel
            tabPanel(
              title = "Annotations",
              #shinyjs::useShinyjs(),
              conditionalPanel(
                condition = "output.showGUI",
                wellPanel(
                  h4("Metabolite family selection"),
                  bsTooltip(id = "familyCount2", title = "The number of available metabolite families which were annotated among the MS\u00B9 features", placement = "bottom", trigger = "hover"),
                  verbatimTextOutput(outputId = "familyCount2"),
                  DT::dataTableOutput("familySelectionTable"),
                  h4("Metabolite family properties")
                )
              ),
              conditionalPanel(
                condition = "!output.showGUI",
                wellPanel(
                  h4("Please import a data file")
                )## well panel
              )## conditional panel
            ),## tab panel
            tabPanel(
              title = "Project", 
              shinyjs::useShinyjs(),
              conditionalPanel(
                condition = "output.showGUI",
                wellPanel(
                  h4("MetFamily project"),
                  fluidRow(
                    column(width = 6, div(style="float:left;width:100%",  
                                          h5("Project name")                           
                    )),##column
                    column(width = 6, div(style="float:right;width:100%", 
                                          bsTooltip(id = "projectName2", title = "The name of the project", placement = "bottom", trigger = "hover"),
                                          textInput(inputId = "projectName2", label = NULL, value = "", width = "100%")
                    ))##column
                  ),##row
                  fluidRow(
                    column(width = 6, div(style="float:left;width:100%",  
                                          h5("Project description")
                    )),##column
                    column(width = 6, div(style="float:right;width:100%", 
                                          bsTooltip(id = "projectDescription2", title = "Please update the description of this project as free text", placement = "bottom", trigger = "hover"),
                                          tags$style(type="text/css", "projectDescription2 {width:100%}"), 
                                          tags$textarea(id = 'projectDescription2', placeholder = 'Comments here', rows = 3, ""),
                                          bsTooltip(id = "updateProjectDescription", title = "Press to update the project description", placement = "bottom", trigger = "hover"),
                                          actionButton(inputId = "updateProjectDescription", label = "Update project description", class="btn-success", width = "100%")
                    ))##column
                  ),##row
                  fluidRow(
                    column(width = 6,
                           div(style="float:left",
                               h4("Data import parameters")
                           )
                    ),##column
                    column(width = 6,
                           div(style="float:right",
                               bsTooltip(id = "displayImportParameters", title = "Display parameters used for the initial data import", placement = "bottom", trigger = "hover"),
                               checkboxInput(inputId = "displayImportParameters", label = "Display import parameters", value = FALSE)
                           )
                    )##column
                  ),##row
                  conditionalPanel(
                    condition = "input.displayImportParameters",
                    h5("Fragment filter"),
                    fluidRow(
                      column(width = 6,
                             bsTooltip(id = "minimumIntensityOfMaximalMS2peak2", title = "A MS/MS spectrum is considered iff the MS/MS feature with maximum intensity is greater or equal than this value", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "minimumIntensityOfMaximalMS2peak2", label = "Min. spectrum intensity", value = 2000)
                      ),##column
                      column(width = 6,
                             bsTooltip(id = "minimumProportionOfMS2peaks2", title = "A MS/MS feature is considered iff the intensity is greater or equal than the maximum intensity times this value", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "minimumProportionOfMS2peaks2", label = "MS/MS peak proportion", value = 0.05)
                      )##column
                    ),##row
                    h5("Neutral losses"),
                    fluidRow(
                      column(width = 6,
                             bsTooltip(id = "neutralLossesPrecursorToFragments2", title = "Include neutral losses relative to the precursor ion, i.e. the m/z difference between the m/z of the precursor ion and the m/z of each fragment ion of the corresponding MS/MS spectrum", placement = "bottom", trigger = "hover"),
                             checkboxInput(inputId = "neutralLossesPrecursorToFragments2", label = "Fragment vs. precursor", value = TRUE)
                      ),##column
                      column(width = 6,
                             bsTooltip(id = "neutralLossesFragmentsToFragments2", title = "Include neutral losses amongst fragment ions, i.e. the m/z difference between the m/z's of all pairs of fragment ions within each MS/MS spectrum", placement = "bottom", trigger = "hover"),
                             checkboxInput(inputId = "neutralLossesFragmentsToFragments2", label = "Fragment vs. fragment", value = FALSE)
                      )##column
                    ),##row
                    h5("Fragment grouping"),
                    fluidRow(
                      column(width = 6,
                             bsTooltip(id = "mzDeviationAbsolute_grouping2", title = "A MS/MS feature is added to a fragment group if the absolute m/z difference is smaller or equal than this value", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "mzDeviationAbsolute_grouping2", label = "m/z deviation (abs.)", value = 0.01)
                      ),##column
                      column(width = 6,
                             bsTooltip(id = "mzDeviationInPPM_grouping2", title = "A MS/MS feature is added to a fragment group if the absolute m/z difference is smaller or equal than the m/z times this value divided by 1,000,000 (<b>p</b>arts <b>p</b>er <b>m</b>illion)", placement = "bottom", trigger = "hover"),
                             textInput(inputId = "mzDeviationInPPM_grouping2", label = "m/z deviation (PPM)", value = 10)
                      )##column
                    ),##row
                    h4("Advanced parameters"),
                    h5("MS\u00B9 feature deisotoping"),
                    bsTooltip(id = "doPrecursorDeisotoping2", title = "If checked, the set of MS\u00B9 features is deisotoped", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "doPrecursorDeisotoping2", label = "MS\u00B9 feature deisotoping done", value = TRUE),
                    conditionalPanel(
                      condition = "input.doPrecursorDeisotoping2",
                      fluidRow(
                        column(width = 6,
                               bsTooltip(id = "mzDeviationAbsolute_precursorDeisotoping2", title = "A MS\u00B9 feature is considered an +1 isotopic peak if the absolute of the m/z difference to the (putative) monoisotopic peak minus 1.0033548378 (=<sup>13</sup>C - <sup>12</sup>C) is smaller or equal than this value (analog for the +2 isotopic peak)", placement = "bottom", trigger = "hover"),
                               textInput(inputId = "mzDeviationAbsolute_precursorDeisotoping2", label = "m/z deviation (abs.)", value = 0.01)
                        ),##column
                        column(width = 6,
                               bsTooltip(id = "mzDeviationInPPM_precursorDeisotoping2", title = "A MS\u00B9 feature is considered an +1 isotopic peak if the absolute of the m/z difference to the (putative) monoisotopic peak minus 1.0033548378 (=<sup>13</sup>C - <sup>12</sup>C) is smaller or equal than the m/z times this value divided by 1,000,000 (<b>p</b>arts <b>p</b>er <b>m</b>illion, analog for the +2 isotopic peak)", placement = "bottom", trigger = "hover"),
                               textInput(inputId = "mzDeviationInPPM_precursorDeisotoping2", label = "m/z deviation (PPM)", value = 10)
                        )##column
                      )##row
                    ),##conditional
                    bsTooltip(id = "maximumRtDifference2", title = "A MS\u00B9 feature is considered an isotopic peak if the absolute of the retention time difference to the (putative) monoisotopic peak is smaller or equal than this value (in minutes)", placement = "bottom", trigger = "hover"),
                    textInput(inputId = "maximumRtDifference2", label = "Retention time difference", value = 0.02),
                    h5("Fragment deisotoping"),
                    bsTooltip(id = "doMs2PeakGroupDeisotoping2", title = "If checked, the set of MS/MS features is deisotoped", placement = "bottom", trigger = "hover"),
                    checkboxInput(inputId = "doMs2PeakGroupDeisotoping2", label = "Fragment deisotoping done", value = TRUE),
                    conditionalPanel(
                      condition = "input.doMs2PeakGroupDeisotoping2",
                      fluidRow(
                        column(width = 6,
                               bsTooltip(id = "mzDeviationAbsolute_ms2PeakGroupDeisotoping2", title = "_A MS/MS feature is considered an +1 isotopic peak if the absolute of the m/z difference to the (putative) monoisotopic peak minus 1.0033548378 (=<sup>13</sup>C - <sup>12</sup>C) is smaller or equal than this value", placement = "bottom", trigger = "hover"),
                               textInput(inputId = "mzDeviationAbsolute_ms2PeakGroupDeisotoping2", label = "m/z deviation (abs.)", value = 0.01)
                        ),##column
                        column(width = 6,
                               bsTooltip(id = "mzDeviationInPPM_ms2PeakGroupDeisotoping2", title = "A MS/MS feature is considered an +1 isotopic peak if the absolute of the m/z difference to the (putative) monoisotopic peak minus 1.0033548378 (=<sup>13</sup>C - <sup>12</sup>C) is smaller or equal than the m/z times this value divided by 1,000,000 (<b>p</b>arts <b>p</b>er <b>m</b>illion)", placement = "bottom", trigger = "hover"),
                               textInput(inputId = "mzDeviationInPPM_ms2PeakGroupDeisotoping2", label = "m/z deviation (PPM)", value = 10)
                        )##column
                      )##row
                    )##conditional
                  ),##conditional
                  h4("Export"),
                  fluidRow(
                    column(width = 6, style="width:50%",
                           div(style="float:left;width:100%",
                               bsTooltip(id = "downloadAllPrecursors", title = "Download the full project file", placement = "bottom", trigger = "hover"),
                               downloadButton(outputId = "downloadAllPrecursors", label = "Export project"),
                               tags$style(type='text/css', "#downloadAllPrecursors { width:100%}")
                           )
                    ),##column
                    column(width = 6, style="width:50%",
                           div(style="float:right;width:100%",
                               bsTooltip(id = "downloadImportParameterSet", title = "Download a parameter file with the parameters which have been used for the initial data import", placement = "bottom", trigger = "hover"),
                               downloadButton(outputId = "downloadImportParameterSet", label = "Export import parameter set"),
                               tags$style(type='text/css', "#downloadImportParameterSet { width:100%}")
                           )
                    )##column
                  ),##row
                  br(),
                  fluidRow(
                    column(width = 6, style="width:50%",
                           div(style="float:left;width:100%",
                               bsTooltip(id = "downloadPcaImage", title = "Download the currently displayed PCA plots as image", placement = "bottom", trigger = "hover"),
                               downloadButton('downloadPcaImage', 'Export PCA as image'),
                               tags$style(type='text/css', "#downloadPcaImage { width:100%}")
                           )
                    ),##column
                    column(width = 6, style="width:50%",
                           div(style="float:right;width:100%",
                               bsTooltip(id = "downloadPcaImageType", title = "The user is able to download the PCA plot as image of different types: Portable Network Graphics (*.png) file, Scalable Vector Graphics (*.svg) file, Portable Document Format (*.pdf) file", placement = "bottom", trigger = "hover"),
                               radioButtons(inputId = "downloadPcaImageType", label = NULL, choices = c("png", "svg", "pdf"), selected = "png", inline = TRUE, width = "100%")
                           )
                    )##column
                  ),##row
                  br(),
                  fluidRow(
                    column(width = 6, style="width:50%",
                           div(style="float:left;width:100%",
                               bsTooltip(id = "downloadHcaImage", title = "Download the currently displayed HCA plots as image", placement = "bottom", trigger = "hover"),
                               downloadButton('downloadHcaImage', 'Export HCA as image'),
                               tags$style(type='text/css', "#downloadHcaImage { width:100%}")
                           )
                    ),##column
                    column(width = 6, style="width:50%",
                           div(style="float:right;width:100%",
                               bsTooltip(id = "downloadHcaImageType", title = "The user is able to download the HCA plot as image of different types: Portable Network Graphics (*.png) file, Scalable Vector Graphics (*.svg) file, Portable Document Format (*.pdf) file", placement = "bottom", trigger = "hover"),
                               radioButtons(inputId = "downloadHcaImageType", label = NULL, choices = c("png", "svg", "pdf"), selected = "png", inline = TRUE, width = "100%")
                           )
                    )##column
                  ),##row
                  br(),
                  fluidRow(
                    column(width = 6, style="width:50%",
                           div(style="float:left;width:100%",
                               bsTooltip(id = "downloadReport", title = "Download the currently displayed HCA plots and PCA plots as report", placement = "bottom", trigger = "hover"),
                               downloadButton('downloadReport', 'Export analysis report'),
                               tags$style(type='text/css', "#downloadReport { width:100%}")
                           )
                    ),##column
                    column(width = 6, style="width:50%",
                           div(style="float:right;width:100%"
                               ## nothing here
                           )
                    )##column
                  )##row
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
        h4(HTML("<b>MetFamily 1.0</b>")),
        fluidRow(
          column(width = 9,
                 helpText(
                   "The MetFamily web application is designed for the identification of regulated metabolite families. This is possible on the basis of metabolite profiles for a set of MS\u00B9 features as well as one MS/MS spectrum for each MS\u00B9 feature. Group-discriminating MS\u00B9 features are identified using a principal component analysis (PCA) of metabolite profiles and metabolite families are identified using a hierarchical cluster analysis (HCA) of MS/MS spectra. Regulated metabolite families are identified by considering group-discriminating MS\u00B9 features from corporate metabolite families."
                 )
          ),##column
          column(width = 3,
                 tags$a(imageOutput(outputId = "ipbImage", width = "100%", height = "100%"), href='http://www.ipb-halle.de/en/', target='_blank')
                 #HTML("<a href='http://www.ipb-halle.de/en/', target='_blank'><img src='logo_ipb_en.png' /></a>")
          )##column
        )## row
      ),## well panel
      wellPanel(
        h4(HTML("<b>Published in <a href='http://pubs.acs.org/doi/abs/10.1021/acs.analchem.6b01569', target='_blank'>Analytical Chemistry (ACS Publications)</a>:</b>")),
        br(),
        h5(HTML("<b>Discovering Regulated Metabolite Families in Untargeted Metabolomics Studies</b>")),
        p(HTML("Hendrik Treutler<sup>1</sup>, Hiroshi Tsugawa<sup>2</sup>, Andrea Porzel<sup>3</sup>, Karin Gorzolka<sup>1</sup>, Alain Tissier<sup>4</sup>, Steffen Neumann<sup>1</sup>, and Gerd Ulrich Balcke<sup>4*</sup>")), 
        p(HTML(paste(
          "<FONT SIZE=-1>",
          "<sup>1</sup>Leibniz Institute of Plant Biochemistry, Department of Stress and Developmental Biology, Weinberg 3, D-06120 Halle/Saale, Germany", "<br>",
          "<sup>2</sup>RIKEN Center for Sustainable Resource Science, Yokohama, Kanagawa 230-0045, Japan", "<br>",
          "<sup>3</sup>Leibniz Institute of Plant Biochemistry, Department of Bioorganic Chemistry, Weinberg 3, D-06120 Halle/Saale, Germany", "<br>",
          "<sup>4</sup>Leibniz Institute of Plant Biochemistry, Department of Cell and Metabolic Biology, Weinberg 3, D-06120 Halle/Saale, Germany", "<br>",
          "<sup>*</sup>Corresponding author: Gerd Ulrich Balcke <a href='mailto:Gerd.Balcke@ipb-halle.de?subject=MetFamily%20request'>Gerd.Balcke@ipb-halle.de</a>
          </FONT>", sep=""
        ))),
        br(),
        h5(HTML("<b>Abstract</b>")),
        p(HTML("The identification of metabolites by mass spectrometry constitutes a major bottleneck which considerably limits the throughput of metabolomics studies in biomedical or plant research. Here, we present a novel approach to analyze metabolomics data from untargeted, data-independent LC-MS/MS measurements. By integrated analysis of MS\u00B9 abundances and MS/MS spectra, the identification of regulated metabolite families is achieved. This approach offers a global view on metabolic regulation in comparative metabolomics. We implemented our approach in the web application “MetFamily”, which is freely available at http://msbi.ipb-halle.de/MetFamily/. MetFamily provides a dynamic link between the patterns based on MS\u00B9-signal intensity and the corresponding structural similarity at the MS/MS level. Structurally related metabolites are annotated as metabolite families based on a hierarchical cluster analysis of measured MS/MS spectra. Joint examination with principal component analysis of MS\u00B9 patterns, where this annotation is preserved in the loadings, facilitates the interpretation of comparative metabolomics data at the level of metabolite families. As a proof of concept, we identified two trichome-specific metabolite families from wild-type tomato Solanum habrochaites LA1777 in a fully unsupervised manner and validated our findings based on earlier publications and with NMR.")),
        br(),
        h5(HTML("<b>Cite</b>")),
        p(HTML(paste(
          "Hendrik Treutler, Hiroshi Tsugawa, Andrea Porzel, Karin Gorzolka, Alain Tissier, Steffen Neumann, and Gerd Ulrich U. Balcke.<br>",
          "Discovering Regulated Metabolite Families in Untargeted Metabolomics Studies.<br>",
          "Analytical chemistry, 88(16):8082-8090, August 2016.<br>",
          "<a href='https://dx.doi.org/10.1021/acs.analchem.6b01569', target='_blank'>doi:10.1021/acs.analchem.6b01569</a>", 
          sep = "")))
      ),## well panel
      wellPanel(
        h4(HTML("<b>Documentation</b>")),
        p(HTML("Please find a user guide for the usage of the MetFamily web application below. ")),
        bsTooltip(id = "downloadDocUserGuide", title = "Download the user guide for the MetFamily web application", placement = "bottom", trigger = "hover"),
        downloadButton('downloadDocUserGuide', 'Download user guide'),
        p(HTML("We provide a specification of the input file format. In the MetFamily publication we demonstrated the usage of MetFamily with UPLC-(-)ESI-SWATH-MS/MS data preprocessed with <a href='http://prime.psc.riken.jp/Metabolomics_Software/MS-DIAL/', target='_blank'>MS-DIAL</a>. In the input specification, however, we demonstrate the generation of these files with GC-EI-MS data processed with <a href='https://bioconductor.org/packages/release/bioc/html/xcms.html', target='_blank'>xcms</a> and <a href='https://bioconductor.org/packages/release/bioc/html/CAMERA.html', target='_blank'>CAMERA</a> We also analyzed idMSMS spectra and LC-MS/MS spectra from DDA in MetFamily.")),
        bsTooltip(id = "downloadDocInputSpecification", title = "Download the input specification of the MetFamily web application", placement = "bottom", trigger = "hover"),
        downloadButton('downloadDocInputSpecification', 'Download input specification')
      ),## well panel
      wellPanel(
        h4(HTML("<b>Feedback</b>")),
        p(HTML("The MetFamily web application is designed to support researchers in the interpretation of comparative metabolomics studies at the level of metabolite families. Please help to improve this tool by comments, bug reports, and feature requests. You can contact Dr. Gerd Balcke via <a href='mailto:Gerd.Balcke@ipb-halle.de?subject=MetFamily'>EMail</a> and <a href='https://github.com/Treutler/MetFamily/issues/new', target='_blank'>issues</a> on GitHub."))
      ),## well panel
      wellPanel(
        h4(HTML("<b>Session info</b>")),
        verbatimTextOutput(outputId = "rInfo")
      )## well panel
    )## tab
  )## navBar
)## shinyUI
