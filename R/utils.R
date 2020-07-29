
getProbMatrix <- function(object, type.data) {
  if (class(object) != "DigitalDLSorter") {
    stop("The object provided is not of DigitalDLSorter class")
  } else if (!any(type.data == c("train", "test"))) {
    stop("type.data argument must be 'train' or 'test'")
  }
  return(object@prob.cell.types[[type.data]]@prob.matrix)
}


showProbPlot <- function(object, type.data, set, type.plot = "maxprob") {
  if (class(object) != "DigitalDLSorter") {
    stop("The object provided is not of DigitalDLSorter class")
  } else if (is.null(object@prob.cell.types) | (length(object@prob.cell.types) == 0)) {
    stop("prob.cell.types slot is empty")
  } else if (!any(type.data == c("train", "test"))) {
    stop("type.data argument must be 'train' or 'test'")
  } else if (length(object@prob.cell.types[[type.data]]) == 0) {
    stop("ProbMatrixCellTypes object has not saved plots")
  } else if (set < 1 | set > 6) {
    stop("set argument must be a number from 1 to 6")
  } else if (!any(type.plot == c("violinplot", "boxplot", "linesplot", "maxprob"))) {
    stop("type.plot argument must be one of the next options: violinplot, boxplot, linesplot or maxprob")
  }
  return(object@prob.cell.types[[type.data]]@plots[[set]][[type.plot]])
}


#' Prepare \code{DigitalDLSorter} objects for saving as RDA file.
#'
#' Prepare \code{DigitalDLSorter} objects that have a \code{DigitalDLSorterDNN}
#' object with trained DNN model. \code{keras} models are not able to be stored
#' natively as R objects (e.g. RData or RDS files). By saving the structure as
#' JSON character object and weights as list object, it is possible recovering
#' the model and carrying out perdictions.
#'
#' With this option, the state of optimizer is not saved, only architecture and
#' weights. It is possible to save completely the model as HDF5 file with
#' \code{\link{saveTrainedModelAsH5}} function and to load into \code{DigitalDLSorter}
#' object with \code{\link{loadTrainedModelFromH5}} function.
#'
#' It is also possible to save a \code{DigitalDLSorter} object as RDS file with
#' \code{saveRDS} function without any type of previous preparation.
#'
#' @param object \code{\link{DigitalDLSorter}} object with \code{trained.data}
#' slot.
#'
#' @export
#'
#' @seealso \code{\link{saveRDS}} \code{\link{saveTrainedModelAsH5}}
#'
preparingToSave <- function(object) {
  if (class(object) != "DigitalDLSorter" ||
      class(object) != "DigitalDLSorterDNN") {
    stop("object provided is not a DigitalDLSorter object")
  }
  if (is.null(trained.model(object))) {
    message("Object provided has not a DigitalDLSorterDNN object. It is not necessary ",
            "prepare the object for saving on disk")
    return(object)
  } else if (is.null(trained.model(object)@model)) {
    message("Object provided has not a trained DNN model. It is not necessary ",
            "prepare the object for saving on disk")
    return(object)
  }
  if (class(trained.model(object)@model) == "list") return(object)
  else {
    trained.model.mod <- .saveModelToJSON(trained.model(object))
    trained.model(object) <- trained.model.mod
    return(object)
  }
}
