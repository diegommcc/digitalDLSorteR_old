## generateBulkSamples: dar la posibilidad de utilizar un backend H5 o hacerlo normal


# type.data == "test", "train", "both"
library(HDF5Array)

generateBulkSamples <- function(object,
                                threads = 2,
                                type.data = "both",
                                file.backend = NULL,
                                verbose = TRUE) {
  if (class(object) != "DigitalDLSorter") {
    stop("The object provided is not of DigitalDLSorter class")
  } else if (is.null(single.cell.sim(object))) {
    stop("single.cell.sim slot is empty")
  } else if (is.null(prob.matrix(object))) {
    stop("prob.matrix slot is empty")
  } else if (!any(type.data == c("train", "test", "both"))) {
    stop("type.data argument must be one of the next options: train, test or both")
  }
  # else if (object@prob.matrix$train@cell.names) {
  #   stop("prob.matrix is not a valid object")
  # }
  if (!is.null(file.backend)) {
    if (file.exists(file.backend)) {
      stop("file.backend already exists. Please provide a correct file path")
    }
  }
  if (threads <= 0) {
    threads <- 1
  }
  if (verbose) {
    message(paste("=== Set parallel environment to", threads, "threads"))
  }
  sim.counts <- assay(object@single.cell.sim)
  sim.counts <- edgeR::cpm.default(sim.counts)

  if (type.data == "both") {
    bulk.counts <- lapply(X = c("train", "test"),
                          FUN = function(x) {
                            if (verbose) {
                              message(paste("\n=== Generating", x, "bulk samples:"))
                            }
                            .generateBulkProfiles(object = object,
                                                  sim.counts = sim.counts,
                                                  type.data = x,
                                                  file.backend = file.backend,
                                                  threads = threads,
                                                  verbose = verbose)
                            })
    names(bulk.counts) <- c("train", "test")
  } else {
    if (verbose) {
      message(paste("\n=== Generating", type.data, "bulk samples:"))
    }
    bulk.counts <- .generateBulkProfiles(object = object,
                                         sim.counts = sim.counts,
                                         type.data = type.data,
                                         file.backend = file.backend,
                                         threads = threads,
                                         verbose = verbose)
    bulk.counts <- list(type.data = bulk.counts)
  }

  object@bulk.sim <- bulk.counts

  message("\nDONE")
  return(object)
}

setBulks <- function (x, c, i) { # para qué vale i
  return(rowSums(c[, x]))
}

.generateBulkProfiles <- function(object,
                                  sim.counts,
                                  type.data,
                                  file.backend,
                                  threads,
                                  verbose) {
  prob.matrix.names <- object@prob.matrix[[type.data]]@cell.names
  bulk.counts <- pbapply::pbapply(X = prob.matrix.names,
                                  MARGIN = 1,
                                  FUN = setBulks,
                                  c = sim.counts,
                                  cl = threads)
  colnames.bulk <- paste("Bulk", seq(dim(bulk.counts)[2]), sep = "_")
  rownames.bulk <- rownames(bulk.counts)
  if (!is.null(file.backend)) {
    if (verbose) {
      message("\nWriting data on disk:\n")
    }
    bulk.counts <- DelayedArray(bulk.counts)
    bulk.counts <- HDF5Array::writeHDF5Array(bulk.counts,
                                             filepath = file.backend,
                                             name = type.data,
                                             verbose = verbose)
    # en R 4.0 implementaron la posibilidad de guardar dimnames
    # de momento lo meto en un slot
  }
  return(SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = bulk.counts),
    rowData = rownames.bulk,
    colData = colnames.bulk))
}

## generateCombinedScaledShuffledSet function ---------------------------------

generateCombinedScaledShuffledSet <- function(object,
                                              type.data,
                                              file.backend = NULL,
                                              verbose = TRUE) {
  if (class(object) != "DigitalDLSorter") {
    stop("The object provided is not of DigitalDLSorter class")
  } else if (is.null(single.cell.sim(object))) {
    stop("single.cell.sim slot is empty")
  } else if (is.null(prob.matrix(object))) {
    stop("prob.matrix slot is empty")
  } else if (is.null(bulk.sim(object))) {
    stop("bulk.sim slot is empty")
  } else if (!any(type.data == c("train", "test", "both"))) {
    stop("type.data argument must be one of the next options: train, test or both")
  } else if (!any(names(bulk.sim(object)) %in% c("train", "test"))) {
    stop("bulk.sim slot is not correctly built")
  }
  if (type.data == "both") {
    if (!all(names(bulk.sim(object)) %in% c("train", "test")))
      stop("If type.data = 'both', bulk.sim slot must contain train and test data")
  } else {
    if (names(bulk.sim(object)) != type.data)
      stop(paste(type.data, "data is not present in bulk.sim slot"))
  }
  if (is.null(file.backend)) {
    combineBulkSCProfiles <- .combineBulkSCProfilesInMemory
  } else {
    if (file.exists(file.backend)){
      stop("file.backend already exists. Please provide a correct file path")
    }
    combineBulkSCProfiles <- .combineBulkSCProfilesHDF5
  }

  ## 1- Coger los datos: bulk.sim.counts, bulk,probs.matrix, cell.set.list
  ## 2- Combinar bulkCounts y scCounts: intersect entre genes y cbind de ambas matrices
  ## 3- Generar SC probs matrix: genera matrices de probabilidad y lo combina con
  # las probabilidades anteriores
  ## 4- Escalar los datos y escribirlo todo en disco


  ### COMBINE bulkCounts and scCounts
  if (type.data == "both") {
    combined.counts <- lapply(X = c("train", "test"),
                          FUN = function(x) {
                            if (verbose) {
                              message(paste("\n=== Combining, scaling and shuffling", x, "counts"))
                            }
                            combineBulkSCProfiles(object = object,
                                                  type.data = x,
                                                  file.backend = file.backend,
                                                  verbose = verbose)
                          })
    names(combined.counts) <- c("train", "test")
  } else {
    if (verbose) {
      message(paste("\n=== Combining, scaling and shuffling", type.data, "counts"))
    }
    combined.counts <- combineBulkSCProfiles(object = object,
                                             type.data = type.data,
                                             file.backend = file.backend,
                                             verbose = verbose)
    combined.counts <- list(type.data = combined.counts)
  }

  object@final.data <- combined.counts

  message("\nDONE")
  return(object)
}


.setConfigHDF5 <- function(file.backend, name) {
  setHDF5DumpFile(file.backend)
  setHDF5DumpName(name)
}


.combineBulkSCProfilesHDF5 <- function(object = object,
                                       type.data = type.data,
                                       file.backend = file.backend,
                                       verbose = verbose) {
  gene.list <- intersect(rowData(bulk.sim(object)[[type.data]])[[1]],
                         rownames(assay(object@single.cell.sim)))
  counts <- assay(bulk.sim(object)[[type.data]])
  .setConfigHDF5(file.backend = file.backend, name = "tpm")
  counts <- DelayedArray::cbind(counts, DelayedArray(assay(object@single.cell.sim)[,
                                 unlist(object@prob.matrix[[type.data]]@set.list)]))

  # no funciona, es para reestablecer el directorio por defecto
  # .setConfigHDF5(file.backend = getHDF5DumpDir(), name = "tpm")
  sample.names <- c(colData(object@bulk.sim[[type.data]])[[1]],
                    colnames(assay(object@single.cell.sim)[,
                                   unlist(object@prob.matrix[[type.data]]@set.list)]))
  tpsm <- matrix(unlist(sapply(X = names(object@prob.matrix[[type.data]]@set.list),
                               FUN = function (x, l) {
                                 v <- rep(0,length(l))
                                 names(v) <- names(l)
                                 v[x] <- 100
                                 return(rep(v, length(l[[x]])))
                               }, l = object@prob.matrix[[type.data]]@set.list)),
                 ncol = length(object@prob.matrix[[type.data]]@set.list), byrow = T)
  colnames(tpsm) <- names(object@prob.matrix[[type.data]]@set.list)
  tpsm <- tpsm[, colnames(object@prob.matrix[[type.data]]@prob.matrix)]
  rownames(tpsm) <- unlist(object@prob.matrix[[type.data]]@set.list)
  probs.matrix <- rbind(tpsm, object@prob.matrix[[type.data]]@prob.matrix)/100
  rownames(probs.matrix) <- c(rownames(tpsm), colData(object@bulk.sim[[type.data]])[[1]])

  # rownames(assay(DDLSChung.1@bulk.sim$train)) <- gene.list
  ## problema: como en este punto no permiten los dinnames, no puedo hacer subset
  # por el nombre de los genes. De momento voy a continuar sin hacer subset, pero
  # lo mejor va a ser pasarse a R 4.0. Mirar el entorno de conda.

  ## scale counts matrix
  # cpms
  prior <- 1L
  lib.sizes <- DelayedArray::colSums(counts)
  prior.count.scaled <- prior * length(lib.sizes) * lib.sizes / sum(lib.sizes)
  counts <- t(log2((t(counts) + prior.count.scaled) /
                     (lib.sizes +  2 * prior.count.scaled) * 1e+06 ))
  # scale
  c.means <- DelayedArray::colMeans(counts)
  c.sd <- DelayedMatrixStats::colSds(counts)
  counts <- (t(counts) - c.means) / c.sd
  s <- sample(seq(dim(counts)[1]))
  counts <- counts[s, ]
  sample.names <- sample.names[s]
  probs.matrix <- probs.matrix[s, ]
  if (verbose) {
    message("Writing data on disk:\n")
  }
  counts <- HDF5Array::writeHDF5Array(x = counts,
                                      filepath = file.backend,
                                      name = type.data,
                                      verbose = verbose)
  ## remove tpm
  rhdf5::h5delete(file = file.backend, name = "tpm")

  return(SummarizedExperiment(assays = list(sim.counts = counts),
                              rowData = sample.names,
                              colData = gene.list,
                              metadata = list(prob.matrix = probs.matrix)))
}


.combineBulkSCProfilesInMemory <- function(object = object,
                                           type.data = type.data,
                                           file.backend = NULL,
                                           verbose = verbose) {
  gene.list <- intersect(rowData(bulk.sim(object)[[type.data]])[[1]],
                         rownames(assay(object@single.cell.sim)))
  # check if bulk.sim is loaded in memory or in HDF5 file
  if (class(assay(bulk.sim(object)[[type.data]])) == "HDF5Matrix") {
    counts <- matrix(assay(bulk.sim(object)[[type.data]]))
  } else {
    counts <- assay(bulk.sim(object)[[type.data]])
  }
  counts <- cbind(counts,
                  assay(object@single.cell.sim)[, unlist(object@prob.matrix[[type.data]]@set.list)])

  sample.names <- c(colData(object@bulk.sim[[type.data]])[[1]],
                    colnames(assay(object@single.cell.sim)[, unlist(object@prob.matrix[[type.data]]@set.list)]))
  tpsm <- matrix(unlist(sapply(X = names(object@prob.matrix[[type.data]]@set.list),
                               FUN = function (x, l) {
                                 v <- rep(0,length(l))
                                 names(v) <- names(l)
                                 v[x] <- 100
                                 return(rep(v, length(l[[x]])))
                               }, l = object@prob.matrix[[type.data]]@set.list)),
                 ncol = length(object@prob.matrix[[type.data]]@set.list), byrow = T)
  colnames(tpsm) <- names(object@prob.matrix[[type.data]]@set.list)
  tpsm <- tpsm[, colnames(object@prob.matrix[[type.data]]@prob.matrix)]
  rownames(tpsm) <- unlist(object@prob.matrix[[type.data]]@set.list)
  probs.matrix <- rbind(tpsm, object@prob.matrix[[type.data]]@prob.matrix)/100
  rownames(probs.matrix) <- c(rownames(tpsm), colData(object@bulk.sim[[type.data]])[[1]])

  # rownames(assay(DDLSChung.1@bulk.sim$train)) <- gene.list

  ## scale counts matrix
  # cpms
  prior <- 1L
  counts <- edgeR::cpm.default(counts)
  # scaling
  counts <- scale(counts)
  # shuffling
  s <- sample(seq(dim(counts)[2]))
  counts <- t(counts[, s])
  sample.names <- sample.names[s]
  probs.matrix <- probs.matrix[s, ]

  return(SummarizedExperiment(assays = list(sim.counts = counts),
                              rowData = sample.names,
                              colData = gene.list,
                              metadata = list(prob.matrix = probs.matrix)))

}
