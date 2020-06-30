# Test file: I want to test if I can add S4 classes as slots in the main class.
# I can do it. Therefore, I can use SCE for the single-cell profiles and even
# SummaryExperiment class for Bulk profiles.

## dependencies: these packages will be introduced in NAMESPACE file
library(SingleCellExperiment)
library(splatter)

setOldClass(Classes = 'package_version')
setClassUnion("SingleCellExperimentOrNULL", c("SingleCellExperiment", "NULL"))
setClassUnion("ZINBParamsOrNULL", c("ZINBParams", "NULL"))

## devuelve un warning porque hay que definir la clase package_version
DigitalDLSorter <- setClass(
  Class = "DigitalDLSorter",
  slots = c(
    single.cell.real = "SingleCellExperimentOrNULL",
    zinb.params = "ZINBParamsOrNULL",
    single.cell.sim = "SingleCellExperimentOrNULL",
    project = "character",
    version = "package_version"
  )
)

# initial constructor (prototype argument is deprecated)
setMethod(
  f = "initialize", signature = "DigitalDLSorter",
  definition = function(.Object,
                        single.cell.real = NULL,
                        zinb.params = NULL,
                        single.cell.sim = NULL,
                        project = "DigitalDLSorterProject",
                        version = packageVersion(pkg = "Seurat")) {
    .Object@single.cell.real <- single.cell.real
    .Object@zinb.params <- zinb.params
    .Object@single.cell.sim <- single.cell.sim
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
  cat("Real single-cell profiles:\n ", dim(sce)[1],
      "features and", dim(sce)[2], "cells\n")
  if (is.null(rownames(sce))) rownames.sce <- "---"
  else rownames.sce <- S4Vectors:::selectSome(rownames(sce), 6)
  if (identical(colnames(sce), character(0))) colnames.sce <- "---"
  else colnames.sce <- S4Vectors:::selectSome(colnames(sce), 6)
  cat("  rownames:", rownames.sce, "\n")
  cat("  colnames:", colnames.sce, "\n")
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
              .sceShow(object@single.cell.real)
            } else {
              .sceShow(DataFrame())
            }
            if (!is.null(object@zinb.params)) {
              .zinbModelShow(object@zinb.params@model)
            }
            if (!is.null(object@single.cell.sim)) {
              .sceShow(object@single.cell.sim)
            }
            cat("Project:", object@project)
          })


