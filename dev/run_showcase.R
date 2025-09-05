library(MetFamily)

# runMetFamily()


# Using project file ------------------------------------------------------

project_path <- system.file("extdata/showcase/Project_file_showcase_annotated.csv.gz", package = "MetFamily")
project <- readClusterDataFromProjectFile(file = project_path)

# and previously saved filter object
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




# Using separate files ----------------------------------------------------

# MS1 intensities
filePeakMatrixPath <- system.file("extdata/showcase/Metabolite_profile_showcase.txt", package = "MetFamily")

# MS2 spectra
fileSpectra <- system.file("extdata/showcase/MSMS_library_showcase.msp", package = "MetFamily")

# Sirius annotations
fileAnnotation <- system.file("extdata/testdata/canopus/canopus1680.txt", package = "MetFamily")

# We use default import parameters
parameterSet <- parameterSetDefault()
parameterSet$minimumIntensityOfMaximalMS2peak <- 2000
parameterSet$minimumProportionOfMS2peaks <- 0.05


dataList <- projectFromFiles(
  ms1_path = filePeakMatrixPath, 
  ms2_path = fileSpectra, 
  siriusFileColumnName = "ClassyFire class",
  parameterSet = parameterSet,
  annot_path = fileAnnotation
)
