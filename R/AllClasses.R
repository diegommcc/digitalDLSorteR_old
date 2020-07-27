#' @importFrom methods setClass setOldClass setClassUnion
#' @import SingleCellExperiment SummarizedExperiment
#' @importClassesFrom splatter ZINBParams
#' @importClassesFrom Matrix dgCMatrix
#' @useDynLib digitalDLSorteR
NULL

setOldClass(Classes = 'package_version')
setOldClass("keras.engine.sequential.Sequential")
setOldClass("keras_training_history")

setClassUnion("MatrixOrNULL", c("matrix", "NULL"))
setClassUnion("ListOrNULL", c("list", "NULL"))
setClassUnion("CharacterOrNULL", c("character", "NULL"))
setClassUnion("SingleCellExperimentOrNULL", c("SingleCellExperiment", "NULL"))
setClassUnion("ZINBParamsOrNULL", c("ZINBParams", "NULL"))
setClassUnion("KerasOrList", c("keras.engine.sequential.Sequential", "list"))


## ProbMatrixCellTypes class ----------------------------------------------------

#' The Class ProbMatrixCellTypes.
#'
#' The ProbMatrixCellTypes class is a data storage class that stores the probability
#' matrix used for the simulation of pseudo-bulk samples. This matrix corresponds
#' with \code{prob.matrix} slot. The rest of slots are additional
#' information generated during the process.
#'
#' As described in Torroja and Sanchez-Cabo, 2019,
#' the probabilities are built by six different methods in order to avoid biases
#' due to the composition of the bulk samples. In \code{plots} slot are stored
#' a representation of these probabilities with the aim of offering a method for
#' monitoring the different sets generated during the process. These plots can
#' be displayed with \code{\link{showProbPlot}} function. See documentation for details.
#'
#' @slot prob.matrix Matrix of probabilities generated for the simulation of
#' bulk samples. Rows correspond with bulk samples which will be generated (\eqn{i}),
#' columns are the cell types present on single-cell data provided (\eqn{j}) and each entry
#' is the proportion of \eqn{i} cell type on \eqn{j} sample.
#' @slot cell.names List with the evolution of the selected metrics during
#' training.
#' @slot set.list Results of the selected metrics on test data.
#' @slot set Matrix with the deconvolution results on test data.
#' Columns are cell types, rows are samples and each entry is the proportion of
#' this cell type on this sample.
#' @slot exclusive.types Optional slot that contains the exclusive cell types
#' on the experiment if they are provided. NULL by default.
#' @slot plots List of lists with the resulting plots generated during the
#' construction of probabilities.
#' @slot type.data Character with the type of data contained: training or test.
#'
#' @references
#' Torroja, C. y Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep Learning algorithm to quantify
#' immune cell populations based on scRNA-Seq data. Frontiers in Genetics 10, 978. doi:
#' \url{10.3389/fgene.2019.00978}
#'
#' @export ProbMatrixCellTypes
#'
ProbMatrixCellTypes <- setClass(
  Class = "ProbMatrixCellTypes",
  slots = c(
    prob.matrix = "MatrixOrNULL",
    cell.names = "MatrixOrNULL",
    set.list = "ListOrNULL",
    set = "CharacterOrNULL",
    exclusive.types = "CharacterOrNULL",
    plots = "ListOrNULL",
    type.data = "CharacterOrNULL"
  )
)

setMethod(
  f = "initialize", signature = "ProbMatrixCellTypes",
  definition = function(
    .Object,
    prob.matrix = NULL,
    cell.names = NULL,
    set.list = NULL,
    set = NULL,
    exclusive.types = NULL,
    plots = NULL,
    type.data = NULL
  ) {
    .Object@prob.matrix <- prob.matrix
    .Object@cell.names <- cell.names
    .Object@set.list <- set.list
    .Object@set <- set
    .Object@exclusive.types = exclusive.types
    .Object@plots <- plots
    .Object@type.data <- type.data
    return(.Object)
  }
)


setValidity(Class = "ProbMatrixCellTypes",
            method = function(object) {
              if (all(object@type.data != c("train", "test"))) {
                return(FALSE)
              } else {
                return(TRUE)
              }
            })

setMethod(f = "show",
          signature = "ProbMatrixCellTypes",
          definition = function(object) {
            # cat("An object of class", class(object), "\n")
            if (is.null(object@prob.matrix)) {
              cat("ProbMatrixCellTypes object empty")
            } else {
              cat(paste0("  Probability matrix for ",
                         object@type.data, ": "))
              cat(paste(dim(object@prob.matrix), c("bulk samples and",
                                                   "cell types"),
                        collapse = " "))
              if (!is.null(object@exclusive.types)) {
                cat("\n    Exclusive types:",
                    paste(object@exclusive.types, collapse = ", "))
              }
            }
          })


## DigitalDLSorterDNN class ----------------------------------------------------

#' The DigitalDLSorterDNN Class.
#'
#' The DigitalDLSorterDNN object stores the trained Deep Neural Network, the
#' training history of selected metrics and the results of prediction on test data.
#'
#' The steps realted with Depp Learning are carried out with \code{keras}
#' package, so the model are stored in a R6 class, system used by the package.
#'
#' @slot model Model resulting from the training of Deep Neural Network. This slot
#' contains a R6 \code{keras.engine.sequential.Sequential} object with the
#' @slot training.history List with the evolution of the selected metrics during
#' training.
#' @slot eval.stats Results of the selected metrics on test data.
#' @slot predict.results Matrix with the deconvolution results on test data.
#' Columns are cell types, rows are samples and each entry is the proportion of
#' this cell type on this sample.
#' @slot cell.types Vector with the cell types to deconvolve.
#' @slot features Vector with features used during training. These features will be
#' used for the folloring predictions.
#'
#' @export DigitalDLSorterDNN
#'
DigitalDLSorterDNN <- setClass(
  Class = "DigitalDLSorterDNN",
  slots = c(
    model = "KerasOrList",
    training.history = "keras_training_history",
    eval.stats = "ListOrNULL",
    predict.results = "matrix",
    cell.types = "character",
    features = "character"
  )
)

setMethod(
  f = "initialize", signature = "DigitalDLSorterDNN",
  definition = function(
    .Object,
    model = NULL,
    training.history = NULL,
    eval.stats = NULL,
    predict.results = NULL,
    cell.types = NULL,
    features = NULL
  ) {
    .Object@model <- model
    .Object@training.history <- training.history
    .Object@eval.stats <- eval.stats
    .Object@predict.results <- predict.results
    .Object@cell.types <- cell.types
    .Object@features <- features
    return(.Object)
  }
)

setMethod(f = "show",
          signature = "DigitalDLSorterDNN",
          definition = function(object) {
            # cat("An object of class", class(object), "\n")
            if (is.null(object@model)) {
              cat("DigitalDLSorterDNN object empty")
            } else {

              cat(paste("Trained model:", object@training.history$params$epochs,
                        "epochs\n"))
              train.metrics <- lapply(object@training.history$metrics,
                                      function(x) x[length(x)])
              cat("  Training metrics (last epoch):\n")
              cat(paste0("    ", names(train.metrics), ": ",
                         lapply(train.metrics, round, 4),
                         collapse = "\n"))
              cat("\n  Evaluation metrics on test data:\n")
              cat(paste0("    ", names(object@eval.stats), ": ",
                         lapply(object@eval.stats, round, 4),
                         collapse = "\n"))
            }
          })

setClassUnion("DigitalDLSorterDNNOrNULL", c("DigitalDLSorterDNN", "NULL"))


## DigitalDLSorter class -------------------------------------------------------

#' The DigitalDLSorter Class.
#'
#' The DigitalDLSorter object is the core of digitalDLSorteR. This object stores
#' the different intermediate data resulting from running pipeline from
#' real single-cell data to the trained Deep Neural Network, including the data
#' on which to carry out the process of devonvolution. In the case that a
#' pre-trained is used, only slots for model and devoncolution data will be used.
#'
#' This object uses other classes to store the different type of data produced
#' during the proccess:
#' \itemize{
#' \item \code{SingleCellExperiment} class for single-cell RNASeq data, using sparse
#' matrix from the \code{Matrix} package (\code{dgCMatrix} class) to store the matrix of
#' counts.
#' \item \code{ZinbParams} class with the estimated params for the simulation of single-cell profiles.
#' \item \code{SummarizedExperiment} class for bulk RNASeq data. In this case, it is possible
#' to load all data in memory or the use of HDF5 files as a backend by
#' \code{DelayedArray} and \code{HDF5Array} packages. See \code{\link{generateBulkSamples}} for details.
#' \item \code{ProbMatrixCellTypes} class for the probability matrices built during the process.
#' See \code{?ProbMatrixCellTypes} for details.
#' \item \code{\link{DigitalDLSorterDNN}} class for store the trained Neural Network. This steps
#' is by \code{keras}, See \code{\link{DigitalDLSorterDNN}} for details.
#' }
#'
#' @slot single.cell.real Real single-cell data stored in a \code{SingleCellExperiment}
#' object. The counts matrix is stored as a \code{dgCMatrix} object for optimize the
#' amount of occupied memory.
#' @slot zinb.params \code{ZinbParams} object with estimated params for the simulation of
#' new single-cell expression profiles.
#' @slot single.cell.sim Simulated single-cell expression profiles.
#' @slot prob.cell.types \code{ProbMatrixCellTypes} class with the probability matrix built
#' for the simulation of bulk RNASeq profiles. These probabilities determine the
#' proportion of single-cell types that will constitute te bulk samples.
#' @slot bulk.sim This slots consists in a list with two elements: train and test simulated bulk RNASeq.
#' This data are stored as a \code{SummarizedExperiment} object. We recommend the
#' use of HDF5 file as a backend due to the large amount of memory they occupy.
#' @slot final.data The final data that will be used for the training and testing the
#' Depp Neural Network. As in the previous slot, it is a list with two items, train
#' and test. The counts matrices are combine with single-cell profiles (if you want),
#' scaled and shuffled for the training.
#' @slot trained.model \code{\link{DigitalDLSorterDNN}} object with the trained model,
#' different metrics obtained during the training and evalutation metrics from
#' the application of the model over test data.
#' @slot deconv.data Optional slot where is possible store the bulk samples for
#' deconvolution. It is a list whose name is the name of the data provided for
#' deconvolution. It is possible store more than one datset to make predictions.
#' It is also possible to carry out the prediction process over
#' files stored in text files. See \code{\link{deconvDigitalDLSorterModel}} for details.
#' @slot deconv.results Slot for store the results from the deconvolution process.
#' It is a list whose name is the name of the data from which they come.
#' @slot project Name of the project
#' @slot version Version of DigitalDLSorteR this object was built under
#'
#' For build a DigitalDLSorter object for train your own model, use
#' \code{\link{CreateDigitalDLSorterObject}} for loading single-cell real data.
#' If you want to deconvolute your data using a pre-trained model, see
#' \code{\link{loadDeconvDataFromFile}} and \code{\link{loadDeconvDataFromSummarizedExperiment}}
#' for details.
#'
#' @exportClass DigitalDLSorter
#' @export DigitalDLSorter
#'
DigitalDLSorter <- setClass(
  Class = "DigitalDLSorter",
  slots = c(
    single.cell.real = "SingleCellExperimentOrNULL",
    zinb.params = "ZINBParamsOrNULL",
    single.cell.sim = "SingleCellExperimentOrNULL",
    prob.cell.types = "ListOrNULL",
    bulk.sim = "ListOrNULL",
    final.data = "ListOrNULL",
    trained.model = "DigitalDLSorterDNNOrNULL",
    deconv.data = "ListOrNULL",
    deconv.results = "ListOrNULL",
    project = "character",
    version = "package_version"
  )
)


setMethod(
  f = "initialize", signature = "DigitalDLSorter",
  definition = function(
    .Object,
    single.cell.real = NULL,
    zinb.params = NULL,
    single.cell.sim = NULL,
    prob.cell.types = NULL,
    bulk.sim = NULL,
    final.data = NULL,
    trained.model = NULL,
    deconv.data = NULL,
    deconv.results = NULL,
    project = "DigitalDLSorterProject",
    version = packageVersion(pkg = "digitalDLSorteR")
  ) {
    .Object@single.cell.real <- single.cell.real
    .Object@zinb.params <- zinb.params
    .Object@single.cell.sim <- single.cell.sim
    .Object@prob.cell.types <- prob.cell.types
    .Object@bulk.sim <- bulk.sim
    .Object@final.data <- final.data
    .Object@trained.model <- trained.model
    .Object@deconv.data <- deconv.data
    .Object@deconv.results <- deconv.results
    .Object@project <- project
    .Object@version <- version
    return(.Object)
  }
)


.sceShow <- function(sce) {
  cat(" ", dim(sce)[1], "features and", dim(sce)[2], "cells\n")
  if (is.null(rownames(sce))) rownames.sce <- "---"
  else rownames.sce <- S4Vectors:::selectSome(rownames(sce), 6)
  if (identical(colnames(sce), character(0))) colnames.sce <- "---"
  else colnames.sce <- S4Vectors:::selectSome(colnames(sce), 6)
  cat("  rownames:", rownames.sce, "\n")
  cat("  colnames:", colnames.sce, "\n")
}

.bulkShow <- function(se) {
  cat("   ", dim(se)[1], "features and", dim(se)[2], "samples\n")
  if (is.null(rowData(se)[[1]])) rownames.se <- "---"
  else rownames.se <- S4Vectors:::selectSome(rowData(se)[[1]], 6)
  if (identical(colnames(se), character(0))) colnames.se <- "---"
  else colnames.se <- S4Vectors:::selectSome(colData(se)[[1]], 6)
  cat("    rownames:", rownames.se, "\n")
  cat("    colnames:", colnames.se, "\n")
}

.finalShow <- function(se) {
  cat("   ", dim(se)[2], "features and", dim(se)[1], "samples: ")
  n.bulk <- sum(grepl("Bulk\\.*", rowData(se)[[1]]))
  n.sc <- abs(n.bulk - dim(se)[1])
  cat(n.bulk, "bulk samples and", n.sc, "single-cell samples\n")
  # if (is.null(rowData(se)[[1]])) rownames.se <- "---"
  # else rownames.se <- S4Vectors:::selectSome(rowData(se)[[1]], 6)
  # if (identical(colnames(se), character(0))) colnames.se <- "---"
  # else colnames.se <- S4Vectors:::selectSome(colData(se)[[1]], 6)
  # cat("    rownames:", rownames.se, "\n")
  # cat("    colnames:", colnames.se, "\n")
}


.zinbModelShow <- function(zinb.model) {
  cat(paste0("ZinbParams object:\n",
             "  ", zinbwave::nSamples(zinb.model), " samples; ",
             "  ", zinbwave::nFeatures(zinb.model), " genes.\n",
             "  ", NCOL(zinbwave::getX_mu(zinb.model)),
             " sample-level covariate(s) (mu); ",
             "  ", NCOL(zinbwave::getX_pi(zinb.model)),
             " sample-level covariate(s) (pi);\n",
             "  ", NCOL(zinbwave::getV_mu(zinb.model)),
             " gene-level covariate(s) (mu); ",
             "  ", NCOL(zinbwave::getV_pi(zinb.model)),
             " gene-level covariate(s) (pi);\n",
             "  ", zinbwave::nFactors(zinb.model), " latent factor(s).\n"))
}


setMethod(f = "show",
          signature = "DigitalDLSorter",
          definition = function(object) {
            cat("An object of class", class(object), "\n")
            if (!is.null(object@single.cell.real)) {
              cat("Real single-cell profiles:\n")
              .sceShow(object@single.cell.real)
            } else {
              cat("Real single-cell profiles:\n")
              .sceShow(DataFrame())
            }
            if (!is.null(object@zinb.params)) {
              .zinbModelShow(object@zinb.params@model)
            }
            if (!is.null(object@single.cell.sim)) {
              cat("Simulated single-cell profiles:\n")
              .sceShow(object@single.cell.sim)
            }
            if (!is.null(object@prob.cell.types)) {
              cat("Probability matrices:\n")
              lapply(X = c("train", "test"), FUN = function(x) {
                if (x %in% names(object@prob.cell.types)) {
                  cat(show(object@prob.cell.types[[x]]), "\n")
                }
              })
            }
            if (!is.null(object@bulk.sim)) {
              cat("Simulated bulk samples:\n")
              lapply(X = c("train", "test"), FUN = function(x) {
                if (x %in% names(object@bulk.sim)) {
                  cat(paste(" ", x, "bulk samples:\n"))
                  .bulkShow(object@bulk.sim[[x]])
                }
              })
            }
            if (!is.null(object@final.data)) {
              cat("Final data samples:\n")
              lapply(X = c("train", "test"), FUN = function(x) {
                if (x %in% names(object@final.data)) {
                  cat(paste(" ", x, "data samples:\n"))
                  .finalShow(object@final.data[[x]])
                }
              })
            }
            if (!is.null(object@trained.model)) {
              cat(show(object@trained.model), "\n")
            }
            cat("Project:", object@project, "\n")
          })
