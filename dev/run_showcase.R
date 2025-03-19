library(MetFamily)

# runMetFamily()

project_path <- system.file("extdata/showcase/Project_file_showcase_annotated.csv.gz", package = "MetFamily")
project <- readClusterDataFromProjectFile(file = project_path)

filter_path <- system.file("extdata/testdata/filterObj.Rdata", package = "MetFamily")
load(filter_path) 

pca <- calculatePCA(dataList=project,
                    filterObj=filterObj, 
                    ms1AnalysisMethod="PCA (Principal Component Analysis)", 
                    scaling="None", 
                    logTransform=FALSE)

pcaDimensionOne <- 1
pcaDimensionTwo <- 2

resultObj <- calcPlotPCAscores(
  pcaObj = pca, 
  dataList = project,
  filterObj = filterObj,
  pcaDimensionOne = pcaDimensionOne, 
  pcaDimensionTwo = pcaDimensionTwo, 
  showScoresLabels = FALSE, 
  xInterval = NULL, 
  yInterval = NULL
)

resultObj <- calcPlotPCAloadings(
  pcaObj = pca, 
  dataList = project,
  # filter = filterObj, 
  pcaDimensionOne = pcaDimensionOne, 
  pcaDimensionTwo = pcaDimensionTwo, 
  selectionFragmentPcaLoadingSet = NULL,
  selectionAnalysisPcaLoadingSet = NULL,
  selectionSearchPcaLoadingSet   = NULL,
  xInterval = NULL, 
  yInterval = NULL,
  loadingsLabels = "None", 
  showLoadingsAbundance = FALSE, 
  showLoadingsFeaturesAnnotated   = TRUE,
  showLoadingsFeaturesUnannotated = TRUE,
  showLoadingsFeaturesSelected    = TRUE,
  showLoadingsFeaturesUnselected  = TRUE
)


pcaObj = pca 
dataList = project
filter = filterObj 
pcaDimensionOne = pcaDimensionOne 
pcaDimensionTwo = pcaDimensionTwo 
selectionFragmentPcaLoadingSet = NULL
selectionAnalysisPcaLoadingSet = NULL
selectionSearchPcaLoadingSet   = NULL
xInterval = NULL 
yInterval = NULL
loadingsLabels = "None" 
showLoadingsAbundance = FALSE 
showLoadingsFeaturesAnnotated   = TRUE
showLoadingsFeaturesUnannotated = TRUE
showLoadingsFeaturesSelected    = TRUE
showLoadingsFeaturesUnselected  = TRUE