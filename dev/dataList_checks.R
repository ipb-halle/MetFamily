# checks for dataList format

lengths(dataList) %>% tibble::enframe() %>% View


statements <- c(
  # what does it do? Has to do with includeIgnoredPrecursors
  artefacts = !any(dataList$annoArrayIsArtifact)
  ,
  labels = length(dataList$precursorLabels) == dataList$numberOfPrecursors
  ,
  annoLength = length(dataList$annoArrayOfLists) == dataList$numberOfPrecursors
  ,
  shortLists = all(dataList$annoArrayOfLists %>% lengths %in% c(0,1))
  ,
  # dataList$annoArrayOfLists %>% {lengths(.) == 1} %>% sum
  annoLength2 = length(dataList$annoPresentAnnotationsList) == length(dataList$annoPresentColorsList)
  
)

all(statements)
which(!statements)
names(statements)[!statements]

