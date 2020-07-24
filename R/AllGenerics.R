#' @include AllClasses.R

## all_generics file

## getters and setters for ProbMatrixCellTypes class ---------------------------

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


## getters and setters for DigitalDLSorterDNN class ----------------------------

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

## cell.types
setGeneric("cell.types", function(object) standardGeneric("cell.types"))
setMethod(f = "cell.types",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@cell.types)

setGeneric("cell.types<-", function(object, value) standardGeneric("cell.types<-"))
setMethod(f = "cell.types<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@cell.types <- value
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

## getters and setters for DigitalDLSorter class -------------------------------

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
