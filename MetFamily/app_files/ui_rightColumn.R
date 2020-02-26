
#output$runRightColumn <- renderUI({
  #print(paste("### GUI ### Generate right GUI"))
  #column(width = state$runRightColumnWidth,
  column(width = 8,
         #########################################################################################
         ## show side bar and change plot
         conditionalPanel(
           condition = "(output.showHCAplotPanel & output.analysisType == 'HCA') | (output.showPCAplotPanel & output.analysisType == 'PCA') | (output.showAnnotationplotPanel & output.analysisType == 'Annotation')",
           #condition = "output.showHCAplotPanel | output.showPCAplotPanel",
           fluidRow(
           #  column(width = 6,
           #         div(style="float:right",
           #             bsTooltip(id = "showSideBar", title = "Display or hide the side bar", placement = "bottom", trigger = "hover"),
           #             checkboxInput(inputId = "showSideBar", label = "Display side bar", value = showSideBar)
           #         )
           #  ),##column
           #  column(width = 6,
                    conditionalPanel(
                      condition = "(output.showHCAplotPanel & output.showPCAplotPanel) | (output.showHCAplotPanel & output.showAnnotationplotPanel) | (output.showAnnotationplotPanel & output.showPCAplotPanel)",
                      div(style="float:left",
                          bsTooltip(id = "changePlot", title = "Switch between HCA plots and PCA plots", placement = "bottom", trigger = "hover"),
                          #c("Display HCA", "Display PCA", "Display Annotation")
                          #radioButtons(inputId = "changePlot", label = NULL, choices = plotsToShow, inline = TRUE, selected = plotToShow)
                          radioButtons(inputId = "changePlot", label = NULL, choices = "Display HCA", inline = TRUE, selected = "Display HCA")
                      )
                    )##conditional
           #  )##column
           )##row
         ),##conditional
         ##############################################################################################
         ##############################################################################################
         ## plots
         
         ##############################################################################################
         ## plots controls
         conditionalPanel(## plot properties
           condition = '(output.analysisType == "HCA" & output.showHCAplotPanel) | (output.analysisType == "PCA" & output.showPCAplotPanel)',
           wellPanel(
             fluidRow(
               column(width = 6,
                      div(style="float:left",
                          h4("Plot controls")
                      )
               ),##column
               column(width = 6,
                      div(style="float:right",
                          bsTooltip(id = "showPlotControls", title = "Display control panels for the plots below", placement = "bottom", trigger = "hover"),
                          #checkboxInput(inputId = "showPlotControls", label = "Show plot controls", value = input$showPlotControls)
                          checkboxInput(inputId = "showPlotControls", label = "Show plot controls", value = FALSE)
                      )
               )##column
             ),##row
             fluidRow(
               ################################################
               ## HCA plot controls
               conditionalPanel(## dendrogram properties
                 condition = 'output.analysisType == "HCA" & output.showHCAplotPanel & input.showPlotControls',
                 column(width = 3,
                        h5("Heatmap content"),
                        tags$div(
                          title="Please select the abundance information you would like to display in the heatmap below the dendrogram",
                          #radioButtons(inputId = "heatmapContent", label = NULL, choices = c("Log-fold-change", "Abundance by group", "Abundance by sample"), selected = input$heatmapContent)
                          radioButtons(inputId = "heatmapContent", label = NULL, choices = c("Log-fold-change", "Abundance by group", "Abundance by sample"), selected = "Log-fold-change")
                        )
                 ),
                 column(width = 3,
                        h5("Heatmap ordering"),
                        tags$div(
                          title="Please select the mode of ordering the heatmap rows below the dendrogram",
                          #radioButtons(inputId = "heatmapOrdering", label = NULL, choices = c("Specified order", "MS1 clustering"), selected = input$heatmapOrdering)
                          radioButtons(inputId = "heatmapOrdering", label = NULL, choices = c("Specified order", "MS1 clustering"), selected = "Specified order")
                        )
                 ),
                 column(width = 3,
                        h5("HCA dendrogram"),
                        bsTooltip(id = "showClusterLabels", title = "Display the labels of cluster nodes and MS\u00B9 feature nodes representing the number of characteristic fragments", placement = "bottom", trigger = "hover"),
                        #checkboxInput(inputId = "showClusterLabels", label = "Show node labels", value = input$showClusterLabels)
                        checkboxInput(inputId = "showClusterLabels", label = "Show node labels", value = TRUE)
                 ),
                 column(width = 3,
                        h5("Precursor label"),
                        tags$div(
                          title="Please select the information you would like to display in the labels below the precursors",
                          #radioButtons(inputId = "hcaPrecursorLabels", label = NULL, choices = c("m/z / RT", "Metabolite name", "Metabolite family"), selected = input$hcaPrecursorLabels)
                          radioButtons(inputId = "hcaPrecursorLabels", label = NULL, choices = c("m/z / RT", "Metabolite name", "Metabolite family"), selected = "m/z / RT")
                        )
                 )
               ),## conditional
               ################################################
               ## PCA plot controls
               conditionalPanel(## scores / loadings properties
                 condition = 'output.analysisType == "PCA" & output.showPCAplotPanel & input.showPlotControls',
                 column(width = 3,
                        h5("PCA scores"),
                        bsTooltip(id = "showScoresLabels", title = "Display scores labels", placement = "bottom", trigger = "hover"),
                        #checkboxInput(inputId = "showScoresLabels", label = "Show labels", value = input$showScoresLabels)
                        checkboxInput(inputId = "showScoresLabels", label = "Show labels", value = TRUE)
                 ),
                 column(width = 3,
                        h5("PCA loadings labels"),
                        tags$div(
                          title="Please select the information you would like to display in the precursor labels of the loadings",
                          #radioButtons(inputId = "loadingsLabels", label = NULL, choices = c("None", "m/z / RT", "Metabolite name", "Metabolite family"), selected = input$loadingsLabels)
                          radioButtons(inputId = "loadingsLabels", label = NULL, choices = c("None", "m/z / RT", "Metabolite name", "Metabolite family"), selected = "None")
                        )
                 ),
                 column(width = 3,
                        h5("PCA loadings shown"),
                        tags$div(
                          title="Please select the MS\u00B9 features you would like to display in the loadings plot",
                          bsTooltip(id = "showLoadingsAbundance", title = "Use abundance in MS\u00B9 to scale the size of loadings nodes", placement = "bottom", trigger = "hover"),
                          checkboxGroupInput(inputId = "showLoadingsFeatures", label = NULL, 
                                             choices = c("Annotated", "Not Annotated", "Selected", "Not Selected"), 
                                             selected = c("Annotated", "Not Annotated", "Selected", "Not Selected")
                                             #selected = c(
                                             #  ifelse(test = state_tabPca$showLoadingsFeaturesAnnotated,   yes = "Annotated", no = ""),
                                             #  ifelse(test = state_tabPca$showLoadingsFeaturesUnannotated, yes = "Not Annotated", no = ""),
                                             #  ifelse(test = state_tabPca$showLoadingsFeaturesSelected,    yes = "Selected", no = ""),
                                             #  ifelse(test = state_tabPca$showLoadingsFeaturesUnselected,  yes = "Not Selected", no = "")
                                             #)
                          )
                          #checkboxInput(inputId = "showLoadingsFeaturesAnnotated",   label = "Annotated",     value = input$showLoadingsFeaturesAnnotated),
                          #checkboxInput(inputId = "showLoadingsFeaturesUnannotated", label = "Not Annotated", value = input$showLoadingsFeaturesUnannotated),
                          #checkboxInput(inputId = "showLoadingsFeaturesSelected",    label = "Selected",      value = input$showLoadingsFeaturesSelected),
                          #checkboxInput(inputId = "showLoadingsFeaturesUnselected",  label = "Not Selected",  value = input$showLoadingsFeaturesUnselected)
                        )
                 ),
                 column(width = 3,
                        h5("PCA loadings size"),
                        bsTooltip(id = "showLoadingsAbundance", title = "Use abundance in MS\u00B9 to scale the size of loadings nodes", placement = "bottom", trigger = "hover"),
                        #checkboxInput(inputId = "showLoadingsAbundance", label = "Scale by abundance", value = input$showLoadingsAbundance)
                        checkboxInput(inputId = "showLoadingsAbundance", label = "Scale by abundance", value = FALSE)
                 )
               )## conditional
             )## row
           )## well
         ),## conditional
         ##############################################################################################
         ## plots
         fluidRow(
           #column(width = 12-state$legendColumnWidth,
           column(width = 12-2,
                  ##############################################################################################
                  ## HCA plots
                  conditionalPanel(
                    condition = 'output.analysisType == "HCA" & output.showHCAplotPanel',
                    fluidRow(
                      #plotlyOutput(
                      #   height = state$dendrogramHeatmapHeight, 
                      #   #height = 500, 
                      #   outputId = "plotDendrogram"
                      # )
                      div(style = "position:relative",
                      plotOutput(height = 500, 
                                 outputId = "plotDendrogram", 
                                 #hover    = "plotDendrogram_hover", 
                                 hover    = hoverOpts(
                                   id = "plotDendrogram_hover",
                                   delay = 50, 
                                   delayType = "debounce",
                                   clip = FALSE,
                                   nullOutside = FALSE
                                 ),
                                 click    = "plotDendrogram_click",
                                 dblclick = "plotDendrogram_dblclick",
                                 #brush    = "plotDendrogram_brush"
                                 brush    = brushOpts(
                                   id = "plotDendrogram_brush",
                                   resetOnNew = TRUE,
                                   direction = "x",
                                   delay = 00,
                                   delayType = "debounce"
                                 )
                      ),
                      uiOutput("plotDendrogram_hover_info")
                      )
                    ),## row
                    fluidRow(
                      div(style = "position:relative",
                        uiOutput("ui_plotHeatmap"),
                        uiOutput("plotHeatmap_hover_info")
                      )
                    )## row
                  ),## conditional
                  ##############################################################################################
                  ## PCA plots
                  conditionalPanel(
                    condition = 'output.analysisType == "PCA" & output.showPCAplotPanel',
                    fluidRow(
                      column(width = 6,
                             div(style = "position:relative",
                               plotOutput(height = 500, 
                                          outputId = "plotPcaScores", 
                                          #hover    = "plotPcaScores_hover",
                                          hover    = hoverOpts(
                                            id = "plotPcaScores_hover",
                                            delay = 50, 
                                            delayType = "debounce"
                                          ),
                                          click    = "plotPcaScores_click",
                                          dblclick = "plotPcaScores_dblclick",
                                          #brush    = "plotPcaScores_brush"
                                          brush = brushOpts(
                                            id = "plotPcaScores_brush",
                                            resetOnNew = TRUE,
                                            direction = "xy",
                                            delay = 00,
                                            delayType = "debounce"
                                          )
                               ),
                               uiOutput("plotPcaScores_hover_info")
                             )
                      ),## column
                      column(width = 6,
                             div(style = "position:relative",
                               plotOutput(height = 500, 
                                          outputId = "plotPcaLoadings", 
                                          #hover    = "plotPcaLoadings_hover",
                                          hover    = hoverOpts(
                                            id = "plotPcaLoadings_hover",
                                            delay = 50, 
                                            delayType = "debounce"
                                          ),
                                          click    = "plotPcaLoadings_click",
                                          dblclick = "plotPcaLoadings_dblclick",
                                          #brush    = "plotPcaScores_brush"
                                          brush = brushOpts(
                                            id = "plotPcaLoadings_brush",
                                            resetOnNew = TRUE,
                                            direction = "xy",
                                            delay = 00,
                                            delayType = "debounce"
                                          )
                               ),
                               uiOutput("plotPcaLoadings_hover_info")
                             )
                      )## column
                    )## row
                  )## conditional
           ),##column for plot controls and plots
           #column(width = state$legendColumnWidth,
           column(width = 2,
                  conditionalPanel(
                    condition = 'output.analysisType == "HCA" & output.showHCAplotPanel',
                    splitLayout(
                      style = "border: 1px solid silver;",
                      uiOutput("ui_plotAnnoLegendHCA")
                    )
                  ),## conditional
                  conditionalPanel(
                    condition = 'output.analysisType == "PCA" & output.showPCAplotPanel',
                    splitLayout(
                      style = "border: 1px solid silver;",
                      uiOutput("ui_plotAnnoLegendPCA")
                    )
                  ),## conditional
                  conditionalPanel(
                    condition = '(output.analysisType == "HCA" & output.showHCAplotPanel) | (output.analysisType == "PCA" & output.showPCAplotPanel)',
                    splitLayout(
                      style = "border: 1px solid silver;",
                      plotOutput(outputId = "calcPlotDendrogramLegend", height = 80)
                    )
                  ),## conditional
                  conditionalPanel(
                    condition = 'output.analysisType == "HCA" & output.showHCAplotPanel',
                    splitLayout(
                      style = "border: 1px solid silver;",
                      plotOutput(outputId = "plotHeatmapLegend", height = 150)
                    )
                  ),## conditional
                  conditionalPanel(## loadings properties
                    condition = 'output.analysisType == "PCA" & output.showPCAplotPanel',
                    splitLayout(
                      style = "border: 1px solid silver;",
                      uiOutput("ui_plotScoresGroupsLegend")
                    )
                  ),
                  conditionalPanel(
                    condition = '(output.analysisType == "HCA" & output.showHCAplotPanel) | (output.analysisType == "PCA" & output.showPCAplotPanel)',
                    splitLayout(
                      style = "border: 1px solid silver;",
                      plotOutput(outputId = "plotMS2Legend", height = 80)
                    )
                  ),## conditional
                  conditionalPanel(
                    condition = '(output.analysisType == "HCA" & output.showHCAplotPanel)',
                    splitLayout(
                      style = "border: 1px solid silver;",
                      plotOutput(outputId = "plotFragmentDiscriminativityLegend", height = 100)
                    )
                  )## conditional
           )## column
         ),## row
         #########################################################################################
         ## MS2 plot and info
         conditionalPanel(
           condition = '(output.showHCAplotPanel & output.analysisType == "HCA") | (output.showPCAplotPanel & output.analysisType == "PCA")',
           fluidRow(
             div(style = "position:relative",
             plotOutput(height = 250, 
                        outputId = "plotMS2",
                        #hover    = "plotMS2_hover",
                        hover    = hoverOpts(
                          id = "plotMS2_hover",
                          delay = 50, 
                          delayType = "debounce"
                        ),
                        click    = "plotMS2_click",
                        dblclick = "plotMS2_dblclick",
                        #brush    = "plotMS2_brush",
                        brush = brushOpts(
                          id = "plotMS2_brush",
                          resetOnNew = TRUE,
                          direction = "x",
                          delay = 00,
                          delayType = "debounce"
                        )
             ),
             uiOutput("plotMS2_hover_info")
             ),
             ##############################################################################################
             ## classifier results on dendrogram node
             conditionalPanel(
               condition = 'output.showPutativeAnnotationsTableFromAnalysis',
               DT::dataTableOutput("putativeAnnotationsTableFromAnalysis")
               #putativeAnnotationsTableFromAnalysisRowSelected
             ),
             ##############################################################################################
             ## infos
             wellPanel(
               #h4("Information"),
               #bsTooltip(id = "information", title = "Information about items in the plot", placement = "bottom", trigger = "hover"),
               #verbatimTextOutput("information"),
               h4("Tip"),
               bsTooltip(id = "tip", title = "Information about operating options", placement = "bottom", trigger = "hover"),
               verbatimTextOutput("tip")
             )## well
           )## row
         ),## conditional
         ##############################################################################################
         ## precursor set selection and annotation
         ## change selection
         conditionalPanel(
           condition = '(output.showHCAplotPanel & output.analysisType == "HCA") | (output.showPCAplotPanel & output.analysisType == "PCA")',
           fluidRow(
             wellPanel(
               h4("MS\u00B9 feature selections"),
               bsTooltip(id = "changeSelection", title = "Switch MS\u00B9 feature selection", placement = "bottom", trigger = "hover"),
               #radioButtons(inputId = "changeSelection", label = NULL, choices = c(selectionAnalysisName, selectionFragmentName, selectionSearchName), selected = changeSelectionCurrentSelection, inline = TRUE),
               radioButtons(inputId = "changeSelection", label = NULL, choices = c("Selection by HCA/PCA", "Selection by fragment", "Selection by search"), selected = "Selection by HCA/PCA", inline = TRUE),
               bsTooltip(id = "selectionInfo", title = "The number of MS\u00B9 features in the current selection", placement = "bottom", trigger = "hover"),
               hr(),
               verbatimTextOutput("selectionInfo"),
               conditionalPanel(
                 condition = 'output.precursorSetSelected',
                 tabsetPanel(id = "precursorSelectionTabs",
                             #tabPanel(title = precursorSelectionTabSelection, 
                             tabPanel(title = "Selection", 
                                      wellPanel(
                                        ## selection infos
                                        tags$div(
                                          style="margin-bottom:5px;",
                                          bsTooltip(id = "metFragLink", title = "Press to send the current MS\u00B9 feature as well as the corresponding MS/MS spectrum to MetFrag", placement = "bottom", trigger = "hover"),
                                          htmlOutput(outputId = "metFragLink")
                                        ),
                                        bsTooltip(id = "downloadSelectedPrecursors", title = "Download a project file which is reduced to the selected set of MS\u00B9 features", placement = "bottom", trigger = "hover"),
                                        downloadButton('downloadSelectedPrecursors', 'Download reduced project file'),
                                        bsTooltip(id = "clearSelection", title = "Press to clear this selection", placement = "bottom", trigger = "hover"),
                                        actionButton(inputId = "clearSelection", label = "Clear selection", class="btn-danger")
                                      )## well
                             ),## tab
                             #tabPanel(title = precursorSelectionTabAnnotation, 
                             tabPanel(title = "Annotation", 
                                      wellPanel(
                                        h4("Present annotation(s)"),
                                        fluidRow(
                                          column(
                                            width = 3,
                                            bsTooltip(id = "presentAnnotationValue", title = "The set of present annotations for the set of selected MS\u00B9 features", placement = "bottom", trigger = "hover"),
                                            selectInput(inputId = "presentAnnotationValue", label = NULL, choices = c("[init]"), selectize = FALSE)
                                          ),## column
                                          column(
                                            width = 3,
                                            bsTooltip(id = "setPresentAnnotationPrimary", title = "Sets the selected annotation primary for the set of selected MS\u00B9 features; i.e. this annotation will be used preferentially for coloring in HCA and PCA", placement = "bottom", trigger = "hover"),
                                            actionButton(inputId = "setPresentAnnotationPrimary", label = "Set primary", class="btn-success")
                                          ),## column
                                          column(
                                            width = 6,
                                            bsTooltip(id = "removePresentAnnotation", title = "Removes the selected annotation from the set of selected MS\u00B9 features", placement = "bottom", trigger = "hover"),
                                            actionButton(inputId = "removePresentAnnotation", label = "Remove annotation", class="btn-danger")
                                          )## column
                                        ),##row
                                        fluidRow(
                                          column(
                                            width = 6,
                                            h4("Add new annotation"),
                                            bsTooltip(id = "newAnnotationValue", title = "The name of this annotation", placement = "bottom", trigger = "hover"),
                                            textInput(inputId = "newAnnotationValue", placeholder = 'Metabolite family name here', label = "Type new annotation"),
                                            bsTooltip(id = "newAnnotationColor", title = "The color of this annotation", placement = "bottom", trigger = "hover"),
                                            colourpicker::colourInput(inputId = "newAnnotationColor", label = "Select annotation color", palette = "limited", showColour = "background", allowedCols = c(
                                              "blue",
                                              "red",
                                              "yellow",
                                              "green",
                                              "brown",
                                              "deepskyblue",
                                              "orange",
                                              "deeppink",
                                              "aquamarine",##
                                              "burlywood", 
                                              "cadetblue",
                                              "coral",
                                              "cornflowerblue",
                                              "cyan",##
                                              "darkblue",
                                              "firebrick",
                                              "goldenrod",
                                              "indianred",
                                              "khaki",##
                                              "magenta",
                                              "maroon",
                                              "beige",
                                              "moccasin",
                                              "olivedrab",
                                              "orangered",
                                              "orchid",
                                              "paleturquoise3",##
                                              "rosybrown",
                                              "salmon",
                                              "seagreen3",
                                              "skyblue",
                                              "steelblue"
                                            )),
                                            bsTooltip(id = "submitNewAnnotation", title = "Adds this annotation to the set of selected MS\u00B9 features", placement = "bottom", trigger = "hover"),
                                            actionButton(inputId = "submitNewAnnotation", label = "Add new annotation", class="btn-success")
                                          ),
                                          column(
                                            width = 6,
                                            h4("Add previous annotation"),
                                            bsTooltip(id = "previousAnnotationValue", title = "The set of annotations which have been assigned before", placement = "bottom", trigger = "hover"),
                                            selectInput(inputId = "previousAnnotationValue", label = "Select previous annotation", choices = c("Artifact"), selectize = FALSE),
                                            bsTooltip(id = "submitPreviousAnnotation", title = "Adds this annotation to the set of selected MS\u00B9 features", placement = "bottom", trigger = "hover"),
                                            actionButton(inputId = "submitPreviousAnnotation", label = "Add previous annotation", class="btn-success")
                                          )
                                        )
                                      )## well
                             ),## tab
                             #tabPanel(title = precursorSelectionTabTable, 
                             tabPanel(title = "Table", 
                                      wellPanel(
                                        h4("Selected MS\u00B9 features"),
                                        bsTooltip(id = "updateArtifactsFromCheckboxes", title = "Adds the annotation \\'ignore\\' to the set of checked MS\u00B9 features in the table", placement = "bottom", trigger = "hover"),
                                        actionButton(inputId = "updateArtifactsFromCheckboxes", label = "Apply annotation 'Ignore' to MS\u00B9 features", class="btn-danger"),
                                        DT::dataTableOutput("ms1FeatureTable")
                                      )## well
                             ),## tab
                             #tabPanel(title = precursorSelectionTabSpectrum, 
                             tabPanel(title = "Fragments", 
                                      wellPanel(
                                        bsTooltip(id = "selectedSpectrum", title = "The selected spectrum or frequent fragments / neutral losses", placement = "bottom", trigger = "hover"),
                                        tags$style(type="text/css", "textarea {width:100%}"),
                                        tags$textarea(id = 'selectedSpectrum', placeholder = 'Nothing selected', rows = 10, "")
                                      )## well
                             )## tab
                 )## tab set
               )## conditional
             )## well
           )## row
         ),##conditional
         ##############################################################################################
         ## classifier annotation
         conditionalPanel(
           condition = "(output.showAnnotationplotPanel & output.analysisType == 'Annotation')",
           fluidRow(
             wellPanel(
               h4("Metabolite family selection"),
               fluidRow(
                 column(width = 6, style="width:70%",
                        div(style="float:left;width:100%",
                            bsTooltip(id = "familyCount", title = "The number of metabolite families with one or more potential MS\u00B9 features", placement = "bottom", trigger = "hover"),
                            verbatimTextOutput(outputId = "familyCount"),
                            tags$style(type='text/css', "#familyCount { width:100%}")
                        )
                 ),##column
                 column(width = 6, style="width:30%",
                        div(style="float:right;width:100%",
                            bsTooltip(id = "downloadAllAnnotationResults", title = "Download the annotation results for all putative metabolite families", placement = "bottom", trigger = "hover"),
                            downloadButton('downloadAllAnnotationResults', 'Download annotation results'),
                            tags$style(type='text/css', "#downloadAllAnnotationResults { width:100%}")
                        )
                 )##column
               ),##row
               DT::dataTableOutput("annotationResultTableClass"),
               conditionalPanel(
                 condition = "output.classifierClassSelected",
                 h4("MS\u00B9 feature annotation"),
                 fluidRow(
                   column(width = 6, style="width:70%",
                          div(style="float:left;width:100%",
                              bsTooltip(id = "classToSpectraCount", title = "The number of potential MS\u00B9 feature hits", placement = "bottom", trigger = "hover"),
                              verbatimTextOutput(outputId = "classToSpectraCount"),
                              tags$style(type='text/css', "#classToSpectraCount { width:100%}")
                          )
                   ),##column
                   column(width = 6, style="width:30%",
                          div(style="float:right;width:100%",
                              bsTooltip(id = "downloadMetaboliteFamilyAnnotationResults", title = "Download the annotation results for all putative metabolite families", placement = "bottom", trigger = "hover"),
                              downloadButton('downloadMetaboliteFamilyAnnotationResults', 'Download annotation results'),
                              tags$style(type='text/css', "#downloadMetaboliteFamilyAnnotationResults { width:100%}")
                          )
                   )##column
                 ),##row
                 DT::dataTableOutput("annotationResultTableFeature"),
                 conditionalPanel(
                   condition = "output.classifierClassMS1featureSelected",
                   h4("MS\u00B9 feature spectrum versus Metabolite family"),
                   plotOutput(height = 250, 
                              outputId = "plotMS2vsClass",
                              #hover    = "plotMS2vsClass_hover",
                              #hover    = hoverOpts(
                              #  id = "plotMS2vsClass_hover",
                              #  delay = 50, 
                              #  delayType = "debounce"
                              #),
                              #click    = "plotMS2vsClass_click",
                              dblclick = "plotMS2vsClass_dblclick",
                              ##brush    = "plotMS2vsClass_brush",
                              brush = brushOpts(
                                id = "plotMS2vsClass_brush",
                                resetOnNew = TRUE,
                                direction = "x",
                                delay = 00,
                                delayType = "debounce"
                              )
                   ),
                   ## annotation stuff
                   bsTooltip(id = "newAnnotationValue2", title = "The name of this annotation", placement = "bottom", trigger = "hover"),
                   textInput(inputId = "newAnnotationValue2", placeholder = 'Metabolite family name', label = "Type new annotation"),
                   bsTooltip(id = "newAnnotationColor2", title = "The color of this annotation", placement = "bottom", trigger = "hover"),
                   colourpicker::colourInput(inputId = "newAnnotationColor2", label = "Select annotation color", palette = "limited", showColour = "background", allowedCols = c(
                     "blue",
                     "red",
                     "yellow",
                     "green",
                     "brown",
                     "deepskyblue",
                     "orange",
                     "deeppink",
                     "aquamarine",##
                     "burlywood", 
                     "cadetblue",
                     "coral",
                     "cornflowerblue",
                     "cyan",##
                     "darkblue",
                     "firebrick",
                     "goldenrod",
                     "indianred",
                     "khaki",##
                     "magenta",
                     "maroon",
                     "beige",
                     "moccasin",
                     "olivedrab",
                     "orangered",
                     "orchid",
                     "paleturquoise3",##
                     "rosybrown",
                     "salmon",
                     "seagreen3",
                     "skyblue",
                     "steelblue"
                   )),
                   bsTooltip(id = "confirmAnnotation", title = "Applies the metabolite family annotation to all confirmed MS\u00B9 features", placement = "bottom", trigger = "hover"),
                   actionButton(inputId = "confirmAnnotation", label = "Apply confirmed annotations", class="btn-success")
                 )## cond ms1
               )## cond class
             )## well
           )## row
         )##conditional
  )##column
#})
