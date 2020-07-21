# Test file: I want to test if I can add S4 classes as slots in the main class.
# I can do it. Therefore, I can use SCE for the single-cell profiles and even
# SummaryExperiment class for Bulk profiles.

## dependencies: these packages will be introduced in NAMESPACE file
library(SingleCellExperiment)
library(splatter)


# ProbMatrixCellTypes class ----------------------------------------------------

setClassUnion("MatrixOrNULL", c("matrix", "NULL"))
setClassUnion("ListOrNULL", c("list", "NULL"))
setClassUnion("CharacterOrNULL", c("character", "NULL"))
setClassUnion("SingleCellExperimentOrNULL", c("SingleCellExperiment", "NULL"))
setClassUnion("ZINBParamsOrNULL", c("ZINBParams", "NULL"))
setOldClass(Classes = 'package_version')

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


## setValidity <-- para el typedata es necesario por ejemplo
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
                cat("\n    Exclusive types:", paste(object@exclusive.types, collapse = ", "))
              }
            }
          })

# aquí quizás sea recomendable hacer que puedas pasarle tanto un objeto DigitalDLSorter
# como un objeto ProbMatrix

## DigitalDLSorteR
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


getProbMatrix <- function(object, type.data) {
  if (class(object) != "DigitalDLSorter") {
    stop("The object provided is not of DigitalDLSorter class")
  } else if (!any(type.data == c("train", "test"))) {
    stop("type.data argument must be 'train' or 'test'")
  }
  return(object@prob.matrix[[type.data]]@prob.matrix)
}

# getters and setters --> en R no existen los atributos privados, pero igual
# no me interesa generar getters y setters para este objeto, ya que no es algo
# que el usuario debiera tocar

## prob.matrix
setGeneric("prob.matrix", function(object) standardGeneric("prob.matrix"))
setMethod(f = "prob.matrix",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@prob.matrix)

setGeneric("prob.matrix<-", function(object, value) standardGeneric("prob.matrix<-"))
setMethod(f = "prob.matrix<-",
          signature = "ProbMatrixCellTypes",
          definition = function(object, value) {
            object@prob.matrix <- value
            return(object)
          })

## cell.names
setGeneric("cell.names", function(object) standardGeneric("cell.names"))
setMethod(f = "cell.names",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@cell.names)

setGeneric("cell.names<-", function(object, value) standardGeneric("cell.names<-"))
setMethod(f = "cell.names<-",
          signature = "ProbMatrixCellTypes",
          definition = function(object, value) {
            object@cell.names <- value
            return(object)
          })

## set.list
setGeneric("set.list", function(object) standardGeneric("set.list"))
setMethod(f = "set.list",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@set.list)

setGeneric("set.list<-", function(object, value) standardGeneric("set.list<-"))
setMethod(f = "set.list<-",
          signature = "ProbMatrixCellTypes",
          definition = function(object, value) {
            object@set.list <- value
            return(object)
          })

## set
setGeneric("set", function(object) standardGeneric("set"))
setMethod(f = "set",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@set)

setGeneric("set<-", function(object, value) standardGeneric("set<-"))
setMethod(f = "set<-",
          signature = "ProbMatrixCellTypes",
          definition = function(object, value) {
            object@set <- value
            return(object)
          })

## plots
setGeneric("plots", function(object) standardGeneric("plots"))
setMethod(f = "plots",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@plots)

setGeneric("plots<-", function(object, value) standardGeneric("plots<-"))
setMethod(f = "plots<-",
          signature = "ProbMatrixCellTypes",
          definition = function(object, value) {
            object@plots <- value
            return(object)
          })


# DigitalDLSorterDNN class ---------------------------------------------------
setOldClass("keras.engine.sequential.Sequential")
setOldClass("keras_training_history")

DigitalDLSorterDNN <- setClass(
  Class = "DigitalDLSorterDNN",
  slots = c(
    model = "keras.engine.sequential.Sequential",
    training.history = "keras_training_history",
    eval.stats = "ListOrNULL",
    predict.results = "matrix",
    classes = "character",
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
    classes = NULL,
    features = NULL
  ) {
    .Object@model <- model
    .Object@training.history <- training.history
    .Object@eval.stats <- eval.stats
    .Object@predict.results <- predict.results
    .Object@classes <- classes
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

# getters and setters for slots ------------------------------------------------
## model
setGeneric("model", function(object) standardGeneric("model"))
setMethod(f = "model",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@model)

setGeneric("model<-", function(object, value) standardGeneric("model<-"))
setMethod(f = "model<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@model <- value
            return(object)
          })

## training.history
setGeneric("training.history", function(object) standardGeneric("training.history"))
setMethod(f = "training.history",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@training.history)
setGeneric("training.history<-", function(object, value) standardGeneric("training.history<-"))
setMethod(f = "training.history<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@training.history <- value
            return(object)
          })

## eval.stats
setGeneric("eval.stats", function(object) standardGeneric("eval.stats"))
setMethod(f = "eval.stats",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@eval.stats)

setGeneric("eval.stats<-", function(object, value) standardGeneric("eval.stats<-"))
setMethod(f = "eval.stats<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@eval.stats <- value
            return(object)
          })

## predict.results
setGeneric("predict.results", function(object) standardGeneric("predict.results"))
setMethod(f = "predict.results",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@predict.results)

setGeneric("predict.results<-", function(object, value) standardGeneric("predict.results<-"))
setMethod(f = "predict.results<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@predict.results <- value
            return(object)
          })

## classes
setGeneric("classes", function(object) standardGeneric("classes"))
setMethod(f = "classes",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@classes)

setGeneric("classes<-", function(object, value) standardGeneric("classes<-"))
setMethod(f = "classes<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@classes <- value
            return(object)
          })

## features
setGeneric("features", function(object) standardGeneric("features"))
setMethod(f = "features",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@features)

setGeneric("features<-", function(object, value) standardGeneric("features<-"))
setMethod(f = "features<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@features <- value
            return(object)
          })

# DigitalDLSorter class --------------------------------------------------------
DigitalDLSorter <- setClass(
  Class = "DigitalDLSorter",
  slots = c(
    single.cell.real = "SingleCellExperimentOrNULL",
    zinb.params = "ZINBParamsOrNULL",
    single.cell.sim = "SingleCellExperimentOrNULL",
    prob.matrix = "ListOrNULL",
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
    prob.matrix = NULL,
    bulk.sim = NULL,
    final.data = NULL,
    trained.model = NULL,
    deconv.data = NULL,
    deconv.results = NULL,
    project = "DigitalDLSorterProject",
    version = packageVersion(pkg = "digitalDLSorterPackageR")
  ) {
    .Object@single.cell.real <- single.cell.real
    .Object@zinb.params <- zinb.params
    .Object@single.cell.sim <- single.cell.sim
    .Object@prob.matrix <- prob.matrix
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

.checkClasses <- function(value, class1, class2 = "NULL") {
  if (class(value) != class1 & class(value) != class2) {
    stop(paste("Remember that single.cell.real slot must be", class1, "or", class2))
  }
}


# getters and setters for slots ------------------------------------------------
## single.cell.real
setGeneric("single.cell.real", function(object) standardGeneric("single.cell.real"))
setMethod(f = "single.cell.real",
          signature = "DigitalDLSorter",
          definition = function(object) object@single.cell.real)

setGeneric("single.cell.real<-", function(object, value) standardGeneric("single.cell.real<-"))
setMethod(f = "single.cell.real<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@single.cell.real <- value
            return(object)
          })

## zinb.params
setGeneric("zinb.params", function(object) standardGeneric("zinb.params"))
setMethod(f = "zinb.params",
          signature = "DigitalDLSorter",
          definition = function(object) object@zinb.params)
setGeneric("zinb.params<-", function(object, value) standardGeneric("zinb.params<-"))
setMethod(f = "zinb.params<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@zinb.params <- value
            return(object)
          })

## single.cell.sim
setGeneric("single.cell.sim", function(object) standardGeneric("single.cell.sim"))
setMethod(f = "single.cell.sim",
          signature = "DigitalDLSorter",
          definition = function(object) object@single.cell.sim)

setGeneric("single.cell.sim<-", function(object, value) standardGeneric("single.cell.sim<-"))
setMethod(f = "single.cell.sim<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@single.cell.sim <- value
            return(object)
          })

## prob.matrix
setGeneric("prob.matrix", function(object, type.data = "both") standardGeneric("prob.matrix"))
setMethod(f = "prob.matrix",
          signature = "DigitalDLSorter",
          definition = function(object, type.data) {
            if (type.data == "train") object@prob.matrix[["train"]]
            else if (type.data == "test") object@prob.matrix[["test"]]
            else if (type.data == "both") object@prob.matrix
            else stop(paste("No", type.data, "in prob.matrix"))
          })

setGeneric("prob.matrix<-", function(object, value, type.data = "both") standardGeneric("prob.matrix<-"))
setMethod(f = "prob.matrix<-",
          signature = "DigitalDLSorter",
          definition = function(object, value, type.data) {
            if (type.data == "train") object@prob.matrix[["train"]] <- value
            else if (type.data == "test") object@prob.matrix[["test"]] <- value
            else if (type.data == "both") object@prob.matrix <- value
            else stop(paste("No", type.data, "in prob.matrix slot"))
            return(object)
          })

## bulk.sim
setGeneric("bulk.sim", function(object, type.data = "both") standardGeneric("bulk.sim"))
setMethod(f = "bulk.sim",
          signature = "DigitalDLSorter",
          definition = function(object, type.data) {
            if (type.data == "train") object@bulk.sim[["train"]]
            else if (type.data == "test") object@bulk.sim[["test"]]
            else if (type.data == "both") object@bulk.sim
            else stop(paste("No", type.data, "in bulk.sim slot"))
          })

setGeneric("bulk.sim<-", function(object, value, type.data = "both") standardGeneric("bulk.sim<-"))
setMethod(f = "bulk.sim<-",
          signature = "DigitalDLSorter",
          definition = function(object, value, type.data) {
            if (type.data == "train") object@bulk.sim[["train"]] <- value
            else if (type.data == "test") object@bulk.sim[["test"]] <- value
            else if (type.data == "both") object@bulk.sim <- value
            else stop(paste("No", type.data, "in bulk.sim slot"))
            return(object)
          })

## final.data
setGeneric("final.data", function(object, type.data = "both") standardGeneric("final.data"))
setMethod(f = "final.data",
          signature = "DigitalDLSorter",
          definition = function(object, type.data) {
            if (type.data == "train") object@final.data[["train"]]
            else if (type.data == "test") object@final.data[["test"]]
            else if (type.data == "both") object@final.data
            else stop(paste("No", type.data, "in bulk.sim slot"))
          })

setGeneric("final.data<-", function(object, value, type.data = "both") standardGeneric("final.data<-"))
setMethod(f = "final.data<-",
          signature = "DigitalDLSorter",
          definition = function(object, value, type.data) {
            if (type.data == "train") object@final.data[["train"]] <- value
            else if (type.data == "test") object@final.data[["test"]] <- value
            else if (type.data == "both") object@final.data <- value
            else stop(paste("No", type.data, "in bulk.sim slot"))
            return(object)
          })

## trained.model
setGeneric("trained.model", function(object) standardGeneric("trained.model"))
setMethod(f = "trained.model",
          signature = "DigitalDLSorter",
          definition = function(object) object@trained.model)

setGeneric("trained.model<-", function(object, value) standardGeneric("trained.model<-"))
setMethod(f = "trained.model<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@trained.model <- value
            return(object)
          })

## deconv.data
setGeneric("deconv.data", function(object, name.data = NULL) standardGeneric("deconv.data"))
setMethod(f = "deconv.data",
          signature = "DigitalDLSorter",
          definition = function(object, name.data) {
            if (is.null(name.data)) object@deconv.data
            else {
              if (!name.data %in% name(object@deconv.data)) {
                stop("name.data provided does not exists in deconv.data slot")
              }
              return(object@deconv.data[[name.data]])
            }
          })

setGeneric("deconv.data<-", function(object, value, name.data = NULL) standardGeneric("deconv.data<-"))
setMethod(f = "deconv.data<-",
          signature = "DigitalDLSorter",
          definition = function(object, value, name.data) {
            if (is.null(name.data)) object@deconv.data <- value
            else {
              if (!name.data %in% name(object@deconv.data)) {
                stop("name.data provided does not exists in deconv.data slot")
              }
              object@deconv.data[[name.data]] <- value
            }
            return(object)
          })

## deconv.results
setGeneric("deconv.results", function(object, name.data = NULL) standardGeneric("deconv.results"))
setMethod(f = "deconv.results",
          signature = "DigitalDLSorter",
          definition = function(object, name.data) {
            if (is.null(name.data)) object@deconv.results
            else {
              if (!name.data %in% name(object@deconv.results)) {
                stop("name.data provided does not exists in deconv.results slot")
              }
              return(object@deconv.results[[name.data]])
            }
          })

setGeneric("deconv.results<-", function(object, value, name.data = NULL) standardGeneric("deconv.results<-"))
setMethod(f = "deconv.results<-",
          signature = "DigitalDLSorter",
          definition = function(object, value, name.data) {
            if (is.null(name.data)) object@deconv.results <- value
            else {
              if (!name.data %in% name(object@deconv.results)) {
                stop("name.data provided does not exists in deconv.results slot")
              }
              object@deconv.results[[name.data]] <- value
            }
            return(object)
          })

## project
setGeneric("project", function(object) standardGeneric("project"))
setMethod(f = "project",
          signature = "DigitalDLSorter",
          definition = function(object) object@project)

setGeneric("project<-", function(object, value) standardGeneric("project<-"))
setMethod(f = "project<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@project <- value
            return(object)
          })


# show method for display the content of the object ----------------------------
# setGeneric("show", function(obj) standardGeneric("show"))
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
            if (!is.null(object@prob.matrix)) {
              cat("Probability matrices:\n")
              lapply(X = c("train", "test"), FUN = function(x) {
                if (x %in% names(object@prob.matrix)) {
                  cat(show(object@prob.matrix[[x]]), "\n")
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
