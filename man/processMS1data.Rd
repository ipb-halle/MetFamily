% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataProcessing.R
\name{processMS1data}
\alias{processMS1data}
\title{Process MS-Dial-like MS1 data.frame}
\usage{
processMS1data(
  sampleNamesToExclude,
  numberOfMS1features,
  precursorLabels,
  grouXXXps,
  metaboliteProfileColumnNames,
  dataColumnIndecesFunctionFromGroupIndex,
  dataColumnsNameFunctionFromGroupIndex,
  dataColumnsNameFunctionFromGroupName,
  dataColumnsNameFunctionFromGroupNames,
  groupNameFunctionFromDataColumnName,
  tagsSector,
  metaboliteProfile,
  progress = FALSE
)
}
\arguments{
\item{sampleNamesToExclude}{}

\item{numberOfMS1features}{}

\item{precursorLabels}{}

\item{grouXXXps}{}

\item{metaboliteProfileColumnNames}{}

\item{dataColumnIndecesFunctionFromGroupIndex}{}

\item{dataColumnsNameFunctionFromGroupIndex}{}

\item{dataColumnsNameFunctionFromGroupName}{}

\item{dataColumnsNameFunctionFromGroupNames}{}

\item{groupNameFunctionFromDataColumnName}{}

\item{tagsSector}{}

\item{metaboliteProfile}{}

\item{progress}{}
}
\description{
Processing of MS-Dial-like MS1 data.frame. Includes calculation 
of MS1 data mean and log-fold-change (LFC) data
}
