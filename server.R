
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)
#source("/home/htreutle/Code/Java/MetSWATH/ClusteringMS2SpectraGUI.R")
source("ClusteringMS2SpectraGUI.R")

#########################################################################################
## global variables

shinyServer(
  func = function(input, output, session) {
    options(shiny.maxRequestSize=500*1024^2) 
    
    #########################################################################################
    ## global variables per user
    lastFile <- NULL
    clusterData <- NULL
    filter <- NULL
    clusterFilter <- NULL
    selectedTreeNode <- NULL
    fragmentsX <- NULL
    fragmentsY <- NULL
    submitButtonValue <- 0
    drawButtonValue <- 0
    
    #########################################################################################
    ## functions
    
    disableActionButton <- function(session, id) {
      session$sendCustomMessage(type="jsCode", list(code = paste(
        "$('#", id, "').prop('disabled', true)", sep=""
      )))
    }
    enableActionButton <- function(session, id) {
      session$sendCustomMessage(type="jsCode", list(code = paste(
        "$('#", id, "').prop('disabled', false)", sep=""
      )))
    }
    
    ##' Parse the input file
    doGetClusterPlot <- function(file){
      print(paste("doGetClusterPlot", is.null(file), is.null(lastFile)))
      
      ## decide whether to update the file content
      useLastFile <- FALSE
      readFile <- FALSE
      if(!is.null(file) & is.null(lastFile))
        readFile <- TRUE
      if(!is.null(file) & !is.null(lastFile)){
        if(file != lastFile)
          readFile <- TRUE
        else
          useLastFile <- TRUE
      }
      if(is.null(file) & !is.null(lastFile))
        useLastFile <- TRUE
      print(paste("doGetClusterPlot???", readFile, useLastFile))
      
      if(readFile){
        lastFile <<- file
        clusterFilter <<- NULL
        selectedTreeNode <- NULL
        fragmentsX <- NULL
        fragmentsY <- NULL
        distances = c(
          "Jaccard", 
          "Jaccard (weighted)", 
          "Similarity (weighted)", 
          "Jaccard (weighted2)", 
          "Similarity (weighted2)", 
          "Jaccard (weighted3)",
          "Manhatten",
          "NDP"
        )
        distance = distances[[6]]
        
        withProgress(message = 'Reading file...', value = 0, {
          #clusterData <<- calcClusterData(file = "/mnt/VOL1/ABT/Alle/Balcke/MetSWATH/data/MS-DIAL/UC Davis/Results/201558139_matrixPrecursorsVersusFragmentsDeisotoped_withoutZerosTest01.txt")
          clusterData <<- readClusterData(file = file, distance = distance, progress = TRUE)
        })
        print(paste("doGetClusterPlot do data", clusterData$minimumMass))
      }
      if(useLastFile){
        ## everything is great
      }
      
      print(paste("doGetClusterPlot:::", is.null(clusterData)))
      
      ## update input values
      switch(as.character(length(clusterData$groups)), 
        "0"={
          groupOne <- NA
          groupTwo <- NA
          selectedOne <- NULL
          selectedTwo <- NULL
        },
        "1"={
          groupOne <- clusterData$groups[[1]]
          groupTwo <- clusterData$groups[[1]]
          selectedOne <- clusterData$groups[[1]]
          selectedTwo <- clusterData$groups[[1]]
        },
        {
          groupOne <- clusterData$groups[[1]]
          groupTwo <- clusterData$groups[[2]]
          selectedOne <- clusterData$groups[[1]]
          selectedTwo <- clusterData$groups[[2]]
        }
      )
      
      updateSelectInput(session = session, inputId = "groupOne", choices = clusterData$groups, selected = selectedOne)
      updateSelectInput(session = session, inputId = "groupTwo", choices = clusterData$groups, selected = selectedTwo)
      updateTextInput(session = session, inputId = "filter_average", value = "0")
      updateTextInput(session = session, inputId = "filter_lfc", value = "0")
      updateTextInput(session = session, inputId = "filter_ms2_masses", value = "")
      updateTextInput(session = session, inputId = "filter_ms2_ppm", value = "20")
      
      #doPerformFiltering()
    }
    
    ##' perform filtering
    doPerformFiltering <- function(){
      
      ####################################################################################
      ## process inputs
      
      ## get inputs
      groupOne <- input$groupOne
      groupTwo <- input$groupTwo
      filter_average <- input$filter_average
      filter_lfc <- input$filter_lfc
      filter_ms2_masses <- input$filter_ms2_masses
      filter_ms2_ppm <- input$filter_ms2_ppm
      print(paste("Observe applyFilters1", groupOne, groupTwo, filter_average, filter_lfc, filter_ms2_masses, filter_ms2_ppm))
      
      ## parse inputs
      if(is.null(filter_average) | nchar(filter_average) == 0)
        filter_average <- NULL
      else
        filter_average <- as.numeric(filter_average)
      
      if(is.null(filter_lfc) | nchar(filter_lfc) == 0)
        filter_lfc <- NULL
      else
        filter_lfc <- as.numeric(filter_lfc)
      
      if(is.null(filter_ms2_masses) | nchar(filter_ms2_masses) == 0)
        filter_ms2_masses <- NULL
      
      if(is.null(filter_ms2_ppm) | nchar(filter_ms2_ppm) == 0)
        filter_ms2_ppm <- NULL
      else
        filter_ms2_ppm <- as.numeric(filter_ms2_ppm)
      
      print(paste("Observe applyFilters2", groupOne, groupTwo, filter_average, filter_lfc, filter_ms2_masses, filter_ms2_ppm))
      print(paste("Observe applyFilters2", is.null(groupOne), is.null(groupTwo), is.null(filter_average), is.null(filter_lfc), is.null(filter_ms2_masses), is.null(filter_ms2_ppm)))
      
      ## check for errors in inputs amd process ms2
      error <- FALSE
      if(!is.null(filter_average))
        error <- error | is.na(filter_average)
      if(!is.null(filter_lfc))
        error <- error | is.na(filter_lfc)
      if(!is.null(filter_ms2_masses)){
        ms2Masses <- strsplit(x = filter_ms2_masses, split = "[,; ]+")[[1]]
        filter_ms2_masses <- vector(mode = "numeric", length = length(ms2Masses))
        for(idx in 1:length(ms2Masses))
          filter_ms2_masses[[idx]] <- as.numeric(ms2Masses[[idx]])
        error <- error | any(is.na(filter_ms2_masses))
      }
      if(!is.null(filter_ms2_ppm))
        error <- error | is.na(filter_ms2_ppm)
      error <- error | (!is.null(filter_ms2_masses) & is.null(filter_ms2_ppm))
      
      print(paste("Observe applyFilters3", error, groupOne, groupTwo, filter_average, filter_lfc, filter_ms2_masses, filter_ms2_ppm))
      if(error){
        output$filteredPrecursors <- renderText({
          print(paste("update output$filteredPrecursors invalid filters", sep = ""))
          paste("There are invalid filter values", sep = "")
        })
        return()
      }
      
      ####################################################################################
      ## do filtering
      filter <- filterData(
        dataList = clusterData, 
        groupOne = groupOne, groupTwo = groupTwo, 
        filter_average = filter_average, filter_lfc = filter_lfc, filter_ms2_masses = filter_ms2_masses, filter_ms2_ppm = filter_ms2_ppm, 
        progress = FALSE
      )
      numberOfPrecursorsFiltered <- length(filter)
      print(paste("Observe applyFilters", numberOfPrecursorsFiltered))
      
      ########################################################################################################
      ## check filter validity
      minimumNumberOfPrecursors <- 6
      maximumNumberOfPrecursors <- 1000000
      if(numberOfPrecursorsFiltered >= minimumNumberOfPrecursors & numberOfPrecursorsFiltered <= maximumNumberOfPrecursors){
        ## filter valid
        print(paste("Observe applyFilters ", minimumNumberOfPrecursors, " <= # <= ", maximumNumberOfPrecursors, sep = ""))
        output$filteredPrecursors <- renderText({
          print(paste("update output$filteredPrecursors ", minimumNumberOfPrecursors, " <= # <= ", maximumNumberOfPrecursors, sep = ""))
          paste("Number of filtered precursors: ", numberOfPrecursorsFiltered, sep = "")
        })
        enableActionButton(session, "drawPlots")
        
        filter <<- list()
        filter$filter <<- filter
        filter$groupOne <<- groupOne
        filter$groupTwo <<- groupTwo
        filter$filter_average <<- filter_average
        filter$filter_lfc <<- filter_lfc
        filter$filter_ms2_masses <<- filter_ms2_masses
        filter$filter_ms2_ppm <<- filter_ms2_ppm
        
      } else {
        ## filter invalid
        filter <<- NULL
        
        disableActionButton(session, "drawPlots")
        
        ## update info
        if(numberOfPrecursorsFiltered == 0){
          print(paste("Observe applyFilters # = 0", sep = ""))
          output$filteredPrecursors <- renderText({
            print(paste("update output$filteredPrecursors # = 0", sep = ""))
            paste("There are no precursors which fulfill the given criteria.", sep = "")
          })
        }
        if(numberOfPrecursorsFiltered > 0 & numberOfPrecursorsFiltered < minimumNumberOfPrecursors){
          print(paste("Observe applyFilters 0 < # < ", minimumNumberOfPrecursors, sep = ""))
          output$filteredPrecursors <- renderText({
            print(paste("update output$filteredPrecursors 0 < # < ", minimumNumberOfPrecursors, sep = ""))
            paste("There are only ", numberOfPrecursorsFiltered, " precursors which fulfill the given criteria. There must be at least more than five precursors to preceed.", sep = "")
          })
        }
        if(numberOfPrecursorsFiltered > maximumNumberOfPrecursors){
          print(paste("Observe applyFilters # > ", maximumNumberOfPrecursors, sep = ""))
          output$filteredPrecursors <- renderText({
            print(paste("update output$filteredPrecursors # > ", maximumNumberOfPrecursors, sep = ""))
            paste("There are ", numberOfPrecursorsFiltered, " precursors which fulfill the given criteria. There must be at most 100 precursors to preceed.", sep = "")
          })
        }
      }
    }
    ##' Calculate plots
    doCalcClusterPlot <- function(){
      print(paste("doCalcClusterPlot", is.null(clusterData), is.null(clusterFilter)))
      if(is.null(clusterData) | is.null(clusterFilter)){
        return()
      }
      
      output$clusterPlotDendrogram <- renderPlot({
        print(paste("output$clusterPlotDendrogram"))
        calcClusterPlotDendrogram(filterList = clusterFilter, nodeIndex = selectedTreeNode)
      })
      output$clusterPlotHeatmap <- renderPlot({
        print(paste("output$clusterPlotHeatmap"))
        calcClusterPlotHeatmap(dataList = clusterData, filterList = clusterFilter)
      })
      output$clusterPlotHeatmapLegend <- renderPlot({
        print(paste("output$clusterPlotHeatmapLegend"))
        calcClusterPlotHeatmapLegend(dataList = clusterData)
      })
      output$clusterPlotMS2 <- renderPlot({
        print(paste("output$clusterPlotMS2"))
        calcClusterPlotMS2(dataList = clusterData, fragmentsX = fragmentsX, fragmentsY = fragmentsY)
      })
    }
    
    #########################################################################################
    ## observer
    
    ##' listen to window resize events
    obsClusterPlotResize <- observe({
      plotWidth  <- session$clientData$output_clusterPlotDendrogram_width
      plotHeight <- session$clientData$output_clusterPlotDendrogram_height
      print(paste("Observe plot resize", plotWidth, plotHeight, is.null(clusterData)))
      
      if(!is.null(clusterData) & !is.null(clusterFilter))
        doCalcClusterPlot()
    })
    
    ##' listen to file input
    obsFile <- observe({
      file <- input$matrixFile$datapath
      fileName <- input$matrixFile$name
      print(paste("Observe file for data", fileName))
      if(!is.null(file)){
        output$fileInfo <- renderText({paste("Processing file ", fileName, "...", sep = "")})
        doGetClusterPlot(file)
        output$fileInfo <- renderText({fileName})
      }
    })
    
    ##' listen to dendrogram clicks
    obsDendrogramClick <- observe({
      posX <- input$clusterPlotDendrogram_click$x
      posY <- input$clusterPlotDendrogram_click$y
      
      plotWidth  <- session$clientData$output_clusterPlotDendrogram_width
      plotHeight <- session$clientData$output_clusterPlotDendrogram_height
      
      print(paste("Observe click", posX, posY, "/", plotWidth, plotHeight))
      
      if(is.null(posX) | is.null(posY))
        return()
      
      ## decide whether the click is close enough to trigger event
      factorX <- plotWidth  / clusterFilter$numberOfPrecursorsFiltered
      factorY <- plotHeight / 1
      
      posX <- posX * factorX
      posY <- posY * factorY
      poiCoordinatesX <- clusterFilter$poiCoordinatesX * factorX
      poiCoordinatesY <- clusterFilter$poiCoordinatesY * factorY
      #posX <- grconvertX(posX, from = "user", to = "npc")
      #posY <- grconvertY(posY, from = "user", to = "npc")
      #poiCoordinatesX <- grconvertX(clusterData$poiCoordinatesX, from = "user", to = "npc")
      #poiCoordinatesY <- grconvertY(clusterData$poiCoordinatesY, from = "user", to = "npc")
      
      distancesX <- poiCoordinatesX - posX
      distancesY <- poiCoordinatesY - posY
      distances <- sqrt(distancesX * distancesX + distancesY * distancesY)
      distanceThreshold <- factorX * 5
      
      minimumIndex <- which.min(distances)
      minimumDistance <- distances[[minimumIndex]]
      minimumLabel <- clusterFilter$poiLabels[[minimumIndex]]
      minimumText <- clusterFilter$poiText[[minimumIndex]]
      
      print(paste("Observe clickEvent?", (minimumDistance <= distanceThreshold)))
      if(minimumDistance > distanceThreshold)
        return()
      
      ## fetch ms2 spectrum
      selectedTreeNode <<- minimumLabel
      resultObj <- getMS2spectrum(dataList = clusterData, filterList = clusterFilter, label = minimumLabel)
      fragmentsX <<- resultObj$fragmentsX
      fragmentsY <<- resultObj$fragmentsY
      
      ## output as message and plots
      output$information <- renderText({
        print(paste("update output$information", resultObj$infoText))
        paste(resultObj$infoText, sep = "; ")
      })
      
      if(!is.null(resultObj$landingPageUrl))
        output$metFragLink <- renderText({
          print(paste("update output$metFragLink", resultObj$landingPageUrl))
          paste("<a href=", gsub(pattern = " ", replacement = "%20", x = resultObj$landingPageUrl)[[1]], " target=\"_blank\">Send to MetFrag!</a>", sep = "")
        })
      else
        output$metFragLink <- renderText({
          print(paste("update output$metFragLink", resultObj$landingPageUrl))
          paste("", sep = "")
        })
      
      output$clusterPlotDendrogram <- renderPlot({
        print(paste("update output$clusterPlotDendrogram"))
        calcClusterPlotDendrogram(filterList = clusterFilter, nodeIndex = minimumLabel)
      })
      
      output$clusterPlotMS2 <- renderPlot({
        print(paste("update output$clusterPlotMS2"))
        calcClusterPlotMS2(
          dataList = clusterData, 
          fragmentsX = resultObj$fragmentsX, 
          fragmentsY = resultObj$fragmentsY
        )
      })
    })
    
    ##' listen to submit button events
    obsApplyFilters <- observe({
      applyFilters <- as.numeric(input$applyFilters)
      
      print(paste("Observe applyFilters", applyFilters))
      
      ####################################################################################
      ## check if button was hit
      if(applyFilters == submitButtonValue)
        return()
      submitButtonValue <<- applyFilters
      
      doPerformFiltering()
    })
    
    ##' listen to draw button events
    obsDraw <- observe({
      drawPlots <- as.numeric(input$drawPlots)
      
      print(paste("Observe drawPlots", drawPlots))
      
      ####################################################################################
      ## check if button was hit
      if(drawPlots == drawButtonValue)
        return()
      drawButtonValue <<- drawPlots
      
      ########################################################################################################
      ## update plots
      
      ## reset selections
      selectedTreeNode <<- NULL
      fragmentsX <<- NULL
      fragmentsY <<- NULL
      
      print(paste("Observe drawPlots2", is.null(filter)))
      
      if(!is.null(filter)){
        ## draw
        print(paste("Observe do drawPlots", sep = ""))
        output$information <- renderText({
          print(paste("update output$information do drawPlots", sep = ""))
          paste("Number of filtered precursors: ", length(filter$filter), sep = "")
        })
        
        clusterFilter <<- calculateCluster(
          dataList = clusterData, filter = filter$filter, 
          groupOne = filter$groupOne, groupTwo = filter$groupTwo, 
          filter_average = filter$filter_average, filter_lfc = filter$filter_lfc, filter_ms2_masses = filter$filter_ms2_masses, filter_ms2_ppm = filter$filter_ms2_ppm
        )
        doCalcClusterPlot()
      } else {
        ## do not draw
        print(paste("Observe do not drawPlots", sep = ""))
        
        ## clear plots
        output$clusterPlotDendrogram <- renderPlot({
          print(paste("reset output$clusterPlotDendrogram"))
          NULL
        })
        output$clusterPlotHeatmap <- renderPlot({
          print(paste("reset output$clusterPlotHeatmap"))
          NULL
        })
        output$clusterPlotHeatmapLegend <- renderPlot({
          print(paste("reset output$clusterPlotHeatmapLegend"))
          NULL
        })
        output$clusterPlotMS2 <- renderPlot({
          print(paste("reset output$clusterPlotMS2"))
          NULL
        })
        
        ## update info
        output$information <- renderText({
          print(paste("update do not drawPlots", sep = ""))
          paste("The given criteria do not result in a valid set of precursors. Please check the set of filters.", sep = "")
        })
      }
    })
    
    ##' suspend listener
    session$onSessionEnded(function() {
      obsClusterPlotResize$suspend()
      obsFile$suspend()
      obsDendrogramClick$suspend()
      #obsGroupSelection$suspend()
      obsApplyFilters$suspend()
      obsDraw$suspend()
    })
    
    #########################################################################################
    ## direct output rendering
    output$fileInfo <- renderText({
      print(paste("init output$fileInfo"))
      if(is.null(lastFile))
        paste("Please select a fragment matrix file in the right panel")
      else
        lastFile
    })
    output$information <- renderText({
      print(paste("init output$information"))
      ""
    })
    
    #########################################################################################
    ## direct output values
    output$showGUI <- reactive({
      print("output$showGUI")
      disableActionButton(session, "drawPlots")
      output$filteredPrecursors <- renderText({
        print(paste("init filteredPrecursors", sep = ""))
        paste("Please perform filtering", sep = "")
      })
      output$information <- renderText({
        print(paste("init information", sep = ""))
        paste("Please perform ploting.", sep = "")
      })
      return(!is.null(input$matrixFile))
    })
    #########################################################################################
    ## output properties
    outputOptions(output, 'showGUI', suspendWhenHidden=FALSE)
  }
)
