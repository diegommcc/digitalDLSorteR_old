% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{preparingToSave}
\alias{preparingToSave}
\title{Prepare \code{DigitalDLSorter} object for saving as RDA file.}
\usage{
preparingToSave(object)
}
\arguments{
\item{object}{\code{\link{DigitalDLSorter}} object with \code{trained.data}
slot.}
}
\description{
Prepare a \code{DigitalDLSorter} object that has a \code{DigitalDLSorterDNN}
object with trained DNN model. \code{keras} models are not able to be stored
natively as R objects (e.g. RData or RDS files). By saving the structure as
JSON character object and weights as list object, it is possible recovering
the model and carrying out predictions.
}
\details{
With this option, the state of optimizer is not saved, only architecture and
weights.

It is possible to save completely the model as HDF5 file with
\code{\link{saveTrainedModelAsH5}} function and to load into
\code{DigitalDLSorter} object with \code{\link{loadTrainedModelFromH5}}
function.

It is also possible to save a \code{DigitalDLSorter} object as RDS file with
\code{saveRDS} function without any type of previous preparation.
}
\seealso{
\code{\link{saveRDS}} \code{\link{saveTrainedModelAsH5}}
}
