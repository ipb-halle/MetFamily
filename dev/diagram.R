library(DiagrammeR)

grViz("digraph{

node [shape = box, fontname = Helvetica]

drawPcaLoadingsPlot -> drawPcaLoadingsPlotImpl -> calcPlotPCAloadings
drawPcaPlots -> {drawPcaScoresPlot, drawPcaLoadingsPlot}

}")


grViz("digraph{
node [shape = box, fontname = Helvetica]
 
doPerformFiltering -> doPerformFiltering_impl -> filterData
 
      
}")