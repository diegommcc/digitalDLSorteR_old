## meter más funciones aca


getProbMatrix <- function(object, type.data) {
  if (class(object) != "DigitalDLSorter") {
    stop("The object provided is not of DigitalDLSorter class")
  } else if (!any(type.data == c("train", "test"))) {
    stop("type.data argument must be 'train' or 'test'")
  }
  return(object@prob.matrix[[type.data]]@prob.matrix)
}


showProbPlot <- function(object, type.data, set, type.plot = "maxprob") {
  if (class(object) != "DigitalDLSorter") {
    stop("The object provided is not of DigitalDLSorter class")
  } else if (is.null(object@prob.matrix) | (length(object@prob.matrix) == 0)) {
    stop("prob.matrix slot is empty")
  } else if (!any(type.data == c("train", "test"))) {
    stop("type.data argument must be 'train' or 'test'")
  } else if (length(object@prob.matrix[[type.data]]) == 0) {
    stop("ProbMatrixCellTypes object has not saved plots")
  } else if (set < 1 | set > 6) {
    stop("set argument must be a number from 1 to 6")
  } else if (!any(type.plot == c("violinplot", "boxplot", "linesplot", "maxprob"))) {
    stop("type.plot argument must be one of the next options: violinplot, boxplot, linesplot or maxprob")
  }
  return(object@prob.matrix[[type.data]]@plots[[set]][[type.plot]])
}


