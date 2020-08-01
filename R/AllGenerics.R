#' @include AllClasses.R
NULL

## all_generics file

## getters and setters for ProbMatrixCellTypes class ---------------------------

## prob.matrix

#' @title Get and set \code{prob.matrix} slot in a \code{ProbMatrixCellTypes} object.
#'
#' @param object A \code{ProbMatrixCellTypes} object.
#'
#' @rdname prob.matrix
#' @export prob.matrix
#'
setGeneric("prob.matrix", function(object) standardGeneric("prob.matrix"))
setMethod(f = "prob.matrix",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@prob.matrix)


#' @param value \code{matrix} object with cell types as columns and samples as
#' rows.
#'
#' @rdname prob.matrix
#' @export prob.matrix<-
#'
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


## getters and setters for DigitalDLSorterDNN class ----------------------------

## model
#' @title Get and set \code{model} slot in a \code{DigitalDLSorterDNN} object.
#'
#' @param object A \code{DigitalDLSorterDNN} object.
#'
#' @rdname model
#' @export model
#'
setGeneric("model", function(object) standardGeneric("model"))
setMethod(f = "model",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@model)

#' @param value A \code{keras.engine.sequential.Sequential} object with a
#' trained DNN model.
#'
#' @rdname model
#' @export model<-
#'
setGeneric("model<-", function(object, value) standardGeneric("model<-"))
setMethod(f = "model<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@model <- value
            return(object)
          })

## training.history

#' @title Get and set \code{training.history} slot in a \code{DigitalDLSorterDNN}
#' object.
#'
#' @param object A \code{DigitalDLSorterDNN} object.
#'
#' @rdname training.history
#' @export training.history
#'
setGeneric("training.history", function(object) standardGeneric("training.history"))
setMethod(f = "training.history",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@training.history)

#' @param value A \code{keras_training_history} object with training history of
#' DNN model
#' @rdname training.history
#' @export training.history<-
#'
setGeneric("training.history<-", function(object, value) standardGeneric("training.history<-"))
setMethod(f = "training.history<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@training.history <- value
            return(object)
          })

## eval.stats.model

#' @title Get and set \code{eval.stats.model} slot in a \code{DigitalDLSorterDNN}
#' object.
#'
#' @param object A \code{DigitalDLSorterDNN} object.
#'
#' @rdname eval.stats.model
#' @export eval.stats.model
#'
setGeneric("eval.stats.model", function(object) standardGeneric("eval.stats.model"))
setMethod(f = "eval.stats.model",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@eval.stats.model)

#' @param value A \code{list} object with the resulting metrics after prediction
#' on test data with DNN model.
#' @rdname eval.stats.model
#' @export eval.stats.model<-
#'
setGeneric("eval.stats.model<-", function(object, value) standardGeneric("eval.stats.model<-"))
setMethod(f = "eval.stats.model<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@eval.stats.model <- value
            return(object)
          })

## predict.results

#' @title Get and set \code{predict.results} slot in a \code{DigitalDLSorterDNN}
#' object.
#'
#' @param object A \code{DigitalDLSorterDNN} object.
#'
#' @rdname predict.results
#' @export predict.results
#'
setGeneric("predict.results", function(object) standardGeneric("predict.results"))
setMethod(f = "predict.results",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@predict.results)

#' @param value A \code{matrix} object with prediction results on test data.
#' @rdname predict.results
#' @export predict.results<-
#'
setGeneric("predict.results<-", function(object, value) standardGeneric("predict.results<-"))
setMethod(f = "predict.results<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@predict.results <- value
            return(object)
          })

## cell.types

#' @title Get and set \code{cell.types} slot in a \code{DigitalDLSorterDNN}
#' object.
#'
#' @param object A \code{DigitalDLSorterDNN} object.
#'
#' @rdname cell.types
#' @export cell.types
#'
setGeneric("cell.types", function(object) standardGeneric("cell.types"))
setMethod(f = "cell.types",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@cell.types)

#' @param value A \code{vector} with cell types considered by DNN model.
#' @rdname cell.types
#' @export cell.types<-
#'
setGeneric("cell.types<-", function(object, value) standardGeneric("cell.types<-"))
setMethod(f = "cell.types<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@cell.types <- value
            return(object)
          })

## features

#' @title Get and set \code{features} slot in a \code{DigitalDLSorterDNN}
#' object.
#'
#' @param object A \code{DigitalDLSorterDNN} object.
#'
#' @rdname features
#' @export features
#'
setGeneric("features", function(object) standardGeneric("features"))
setMethod(f = "features",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@features)

#' @param value A \code{vector} with features (genes) considered by DNN model.
#' @rdname features
#' @export features<-
#'
setGeneric("features<-", function(object, value) standardGeneric("features<-"))
setMethod(f = "features<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@features <- value
            return(object)
          })

## eval.stats.samples

#' @title Get and set \code{eval.stats.samples} slot in a \code{DigitalDLSorterDNN}
#' object.
#'
#' @param object A \code{DigitalDLSorterDNN} object.
#'
#' @rdname eval.stats.samples
#' @export eval.stats.samples
#'
setGeneric(
  name = "eval.stats.samples",
  def = function(object, metrics = "All") standardGeneric("eval.stats.samples")
)
setMethod(f = "eval.stats.samples",
          signature = "DigitalDLSorterDNN",
          definition = function(object, metrics) {
            if (metrics == "All") object@eval.stats.samples
            else {
              if (!all(metrics %in% names(object@eval.stats.samples)))
                stop("Metric provided is not present in DigitalDLSorterDNN object")
              return(object@eval.stats.samples[[metrics]])
            }
          })

#' @param value A \code{list} with evaluation metrics used for evaluating the
#' performance of the model over each sample from test data.
#' @rdname eval.stats.samples
#' @export eval.stats.samples<-
#'
setGeneric(
  name = "eval.stats.samples<-",
  def = function(object, value, metrics = "All") {
    standardGeneric("eval.stats.samples<-")
  }
)
setMethod(f = "eval.stats.samples<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value, metrics) {
            if (metrics == "All") object@eval.stats.samples <- value
            else {
              if (!all(metrics %in% names(object@eval.stats.samples)))
                stop("Metric provided is not present in DigitalDLSorterDNN object")
              object@eval.stats.samples[[metrics]] <- value
            }
            return(object)
          })


## getters and setters for DigitalDLSorter class -------------------------------

## single.cell.real

#' @title Get and set \code{single.cell.real} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#'
#' @rdname single.cell.real
#' @export single.cell.real
#'
setGeneric("single.cell.real", function(object) standardGeneric("single.cell.real"))
setMethod(f = "single.cell.real",
          signature = "DigitalDLSorter",
          definition = function(object) object@single.cell.real)

#' @param value A \code{SingleCellExperiment} object with real single-cell
#' profiles.
#' @rdname single.cell.real
#' @export single.cell.real<-
#'
setGeneric("single.cell.real<-", function(object, value) standardGeneric("single.cell.real<-"))
setMethod(f = "single.cell.real<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@single.cell.real <- value
            return(object)
          })

## zinb.params

#' @title Get and set \code{zinb.params} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#'
#' @rdname zinb.params
#' @export zinb.params
#'
setGeneric("zinb.params", function(object) standardGeneric("zinb.params"))
setMethod(f = "zinb.params",
          signature = "DigitalDLSorter",
          definition = function(object) object@zinb.params)

#' @param value A \code{ZinbParams} object with ZiNB-WaVE parameters estimated
#' from real single-cell profiles.
#'
#' @rdname zinb.params
#' @export zinb.params<-
#'
setGeneric("zinb.params<-", function(object, value) standardGeneric("zinb.params<-"))
setMethod(f = "zinb.params<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@zinb.params <- value
            return(object)
          })

## single.cell.sim

#' @title Get and set \code{single.cell.sim} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#'
#' @rdname single.cell.sim
#' @export single.cell.sim
#'
setGeneric("single.cell.sim", function(object) standardGeneric("single.cell.sim"))
setMethod(f = "single.cell.sim",
          signature = "DigitalDLSorter",
          definition = function(object) object@single.cell.sim)

#' @param value A \code{SingleCellExperiment} object with real and simulated
#' single-cell profiles.
#' @rdname single.cell.sim
#' @export single.cell.sim<-
#'
setGeneric("single.cell.sim<-", function(object, value) standardGeneric("single.cell.sim<-"))
setMethod(f = "single.cell.sim<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@single.cell.sim <- value
            return(object)
          })

## prob.cell.types

#' @title Get and set \code{prob.cell.types} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#' @param type.data Element of the list. Can be 'train', 'test' or 'both' (the
#' last by default).
#'
#' @rdname prob.cell.types
#' @export prob.cell.types
#'
setGeneric("prob.cell.types", function(object, type.data = "both") standardGeneric("prob.cell.types"))
setMethod(f = "prob.cell.types",
          signature = "DigitalDLSorter",
          definition = function(object, type.data) {
            if (type.data == "train") object@prob.cell.types[["train"]]
            else if (type.data == "test") object@prob.cell.types[["test"]]
            else if (type.data == "both") object@prob.cell.types
            else stop(paste("No", type.data, "in prob.cell.types"))
          })

#' @param value A list with two elements, train and test, each one with a
#' \code{ProbMatrixCellTypes} object.
#' @rdname prob.cell.types
#' @export prob.cell.types<-
#'
setGeneric("prob.cell.types<-", function(object, value, type.data = "both") standardGeneric("prob.cell.types<-"))
setMethod(f = "prob.cell.types<-",
          signature = "DigitalDLSorter",
          definition = function(object, value, type.data) {
            if (type.data == "train") object@prob.cell.types[["train"]] <- value
            else if (type.data == "test") object@prob.cell.types[["test"]] <- value
            else if (type.data == "both") object@prob.cell.types <- value
            else stop(paste("No", type.data, "in prob.cell.types slot"))
            return(object)
          })

## bulk.sim

#' @title Get and set \code{bulk.sim} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#' @param type.data Element of the list. Can be 'train', 'test' or 'both' (the
#' last by default).
#'
#' @rdname bulk.sim
#' @export bulk.sim
#'
setGeneric("bulk.sim", function(object, type.data = "both") standardGeneric("bulk.sim"))
setMethod(f = "bulk.sim",
          signature = "DigitalDLSorter",
          definition = function(object, type.data) {
            if (type.data == "train") object@bulk.sim[["train"]]
            else if (type.data == "test") object@bulk.sim[["test"]]
            else if (type.data == "both") object@bulk.sim
            else stop(paste("No", type.data, "in bulk.sim slot"))
          })

#' @param value A \code{list} with two elements, train and test, each one being a
#' \code{SummarizedExperiment} object with simulated bulk RNA-Seq samples.
#'
#' @rdname bulk.sim
#' @export bulk.sim<-
#'
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

#' @title Get and set \code{final.data} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#' @param type.data Element of the list. Can be 'train', 'test' or 'both' (the
#' last by default).
#'
#' @rdname final.data
#' @export final.data
#'
setGeneric("final.data", function(object, type.data = "both") standardGeneric("final.data"))
setMethod(f = "final.data",
          signature = "DigitalDLSorter",
          definition = function(object, type.data) {
            if (type.data == "train") object@final.data[["train"]]
            else if (type.data == "test") object@final.data[["test"]]
            else if (type.data == "both") object@final.data
            else stop(paste("No", type.data, "in bulk.sim slot"))
          })

#' @param value A \code{list} with two elements, train and test, each one being a
#' \code{SummarizedExperiment} object with simulated bulk RNA-Seq samples prepared for
#' training. This samples have been normalized and shuffled.
#' @rdname final.data
#' @export final.data<-
#'
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

#' @title Get and set \code{trained.model} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#'
#' @rdname trained.model
#' @export trained.model
#'
setGeneric("trained.model", function(object) standardGeneric("trained.model"))
setMethod(f = "trained.model",
          signature = "DigitalDLSorter",
          definition = function(object) object@trained.model)

#' @param value A \code{DigitalDLSorterDNN} object.
#' @rdname trained.model
#' @export trained.model<-
#'
setGeneric("trained.model<-", function(object, value) standardGeneric("trained.model<-"))
setMethod(f = "trained.model<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@trained.model <- value
            return(object)
          })

## deconv.data

#' @title Get and set \code{deconv.data} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#' @param name.data Name of the data. If it is \code{NULL} (by default),
#' all data contained in \code{deconv.data} slot are returned.
#'
#' @rdname deconv.data
#' @export deconv.data
#'
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

#' @param value A \code{list} whose names are the reference of the data stored.
#' @rdname deconv.data
#' @export deconv.data<-
#'
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

#' @title Get and set \code{deconv.results} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#' @param name.data Name of the data. If it is \code{NULL} (by default),
#' all results contained in \code{deconv.results} slot are returned.
#'
#' @rdname deconv.results
#' @export deconv.results
#'
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

#' @param value A \code{list} whose names are the reference of the results stored.
#' @rdname deconv.results
#' @export deconv.results<-
#'
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

#' @title Get and set \code{project} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#'
#' @rdname project
#' @export project
#'
setGeneric("project", function(object) standardGeneric("project"))
setMethod(f = "project",
          signature = "DigitalDLSorter",
          definition = function(object) object@project)

#' @param value A character indicating the name of the project.
#' @rdname project
#' @export project<-
#'
setGeneric("project<-", function(object, value) standardGeneric("project<-"))
setMethod(f = "project<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@project <- value
            return(object)
          })



## generic for save function: save as RDA a keras model
# setGeneric("save", function(object, file, name) standardGeneric("save"))
#
# setMethod("save", "DigitalDLSorterDNN", definition = function(object, file, name) {
#   if ("keras.engine.sequential.Sequential" %in% class(model(object))) {
#     object <- list(object)
#     names(object) <- name
#     object[[1]] <- .saveModelToJSON(object[[1]])
#     base::save(list = names(object), file = file, envir = list2env(object))
#   } else if (class(model(object)) == "list") {
#     object <- list(object)
#     names(object) <- name
#     base::save(list = names(object), file = file, envir = list2env(object))
#   } else {
#     stop("No valid DigitalDLSorterDNN object")
#   }
# })

# setMethod("save", "DigitalDLSorter", definition = function(object, file, name) {
#   if (!is.null(trained.model(object))) {
#     if ("keras.engine.sequential.Sequential" %in% class(model(object@trained.model))) {
#       name.var <- deparse(substitute(object))
#       object <- .saveModelFromJSON(object)
#       assign(x = eval(substitute(name.var)), value = object)
#       base::save(name.var, file = file)
#     } else if (class(model(object)) == "list") {
#       name.var <- deparse(substitute(object))
#       assign(x = eval(substitute(name.var)), value = object)
#       base::save(get(name.var), file = file)
#     } else {
#       stop("No valid DigitalDLSorterDNN object")
#     }
#   } else {
#     base::save(object, file = file)
#   }
#
# })


#' Save \code{DigitalDLSorter} object as RDS file.
#'
#' Save \code{DigitalDLSorter} and \code{DigitalDLSorterDNN} objects as RDS files.
#' We developed this generic with the aim of changing the behavior of the function
#' and saving the structure and weights of DNN model as text objects. This is
#' because \code{keras} models are not able to be stored natively as R objects
#' (e.g. RData or RDS files). By saving the structure as JSON character object and
#' weights as list object, it is possible recovering the model and carrying
#' out perdictions.
#'
#' With this option, the state of optimizer is not saved, only architecture and
#' weights. It is possible to save completely the model as HDF5 file with
#' \code{\link{saveTrainedModelAsH5}} function and to load into \code{DigitalDLSorter}
#' object with \code{\link{loadTrainedModelFromH5}} function.
#'
#' Moreover, if you want to save the object as rda file, it is possible
#' by transforming before the model to an allowed R object with \code{preparingToSave}
#' function.
#
#' @inheritParams saveRDS
#'
#' @export
#'
#' @seealso \code{\link{saveTrainedModelAsH5}} \code{\link{preparingToSave}}
#'
#'
setGeneric("saveRDS", function(
  object,
  file,
  ascii = FALSE,
  version = NULL,
  compress = TRUE,
  refhook = NULL
) {
  standardGeneric("saveRDS")
})

setMethod("saveRDS", "DigitalDLSorterDNN", definition = function(
  object,
  file,
  ascii,
  version,
  compress,
  refhook
) {
  if ("keras.engine.sequential.Sequential" %in% class(model(object))) {
    object <- .saveModelToJSON(object)
    base::saveRDS(
      object = object,
      file = file,
      ascii = ascii,
      version = version,
      compress = compress,
      refhook = refhook
    )
  } else if (class(model(object)) == "list") {
    base::saveRDS(
      object = object,
      file = file,
      ascii = ascii,
      version = version,
      compress = compress,
      refhook = refhook
    )
  } else {
    stop("No valid DigitalDLSorterDNN object")
  }
})

setMethod("saveRDS", "DigitalDLSorter", definition = function(
  object,
  file,
  ascii,
  version,
  compress,
  refhook
) {
  if (!is.null(trained.model(object))) {
    if ("keras.engine.sequential.Sequential" %in% class(trained.model(object)@model)) {
      model.object <- .saveModelToJSON(trained.model(object))
      trained.model(object) <- model.object
    }
  }
  base::saveRDS(
    object = object,
    file = file,
    ascii = ascii,
    version = version,
    compress = compress,
    refhook = refhook
  )
})

setGeneric("barPlotCellTypes", function(
  data,
  colors = NULL,
  color.line = NA,
  x.label = "Bulk samples",
  rm.x.text = FALSE,
  title = "Results of deconvolution",
  legend.title = "Cell types",
  angle = 90,
  ...
) {
  standardGeneric("barPlotCellTypes")
})

setMethod(f = "barPlotCellTypes",
          signature(data = "DigitalDLSorter"),
          definition = function(
            data,
            name.data = NULL,
            colors = NULL,
            color.line = NA,
            x.label = "Bulk samples",
            rm.x.text = FALSE,
            title = "Results of deconvolution",
            legend.title = "Cell types",
            angle = 90
          ) {
            if (is.null(deconv.results(data))) {
              stop("There is not results to show")
            } else if (is.null(name.data)) {
              message("'name.data' not provided. By default, is catch the first results")
              name.data
            }
            plot <- .barPlot(
              data = deconv.results(data)[[name.data]],
              colors = colors,
              color.line = color.line,
              x.label = x.label,
              rm.x.text = rm.x.text,
              title = title,
              legend.title = legend.title,
              angle = angle
            )
            return(plot)
          })

setMethod(f = "barPlotCellTypes",
          signature(data = "ANY"),
          definition = function(
            data,
            colors = NULL,
            color.line = NA,
            x.label = "Bulk samples",
            rm.x.text = FALSE,
            title = "Results of deconvolution",
            legend.title = "Cell types",
            angle = 90
          ) {
            plot <- .barPlot(
              data = data,
              colors = colors,
              color.line = color.line,
              x.label = x.label,
              rm.x.text = rm.x.text,
              title = title,
              legend.title = legend.title,
              angle = angle
            )
            return(plot)
          })
