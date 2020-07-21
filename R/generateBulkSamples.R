## generateBulkSamples: dar la posibilidad de utilizar un backend H5 o hacerlo normal


# type.data == "test", "train", "both"
library(HDF5Array)

generateBulkSamples <- function(
  object,
  threads = 2,
  type.data = "both",
  file.backend = NULL,
  verbose = TRUE
) {
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
  sim.counts <- edgeR::cpm.default(sim.counts) # sure???

  if (type.data == "both") {
    if (!is.null(object@bulk.sim)) {
      warning("'bulk.sim' slot will be overwritten",
              call. = FALSE, immediate. = TRUE)
    }
    bulk.counts <- lapply(
      X = c("train", "test"),
      FUN = function(x) {
        if (verbose) {
          message(paste("\n=== Generating", x, "bulk samples:"))
        }
        .generateBulkProfiles(
          object = object,
          sim.counts = sim.counts,
          type.data = x,
          file.backend = file.backend,
          threads = threads,
          verbose = verbose
        )
      }
    )
    names(bulk.counts) <- c("train", "test")
    object@bulk.sim <- bulk.counts
  } else {
    if (!is.null(object@bulk.sim) && type.data %in% names(object@bulk.sim)) {
      warning(paste(type.data, "data in 'bulk.sim' slot will be overwritten"),
              call. = FALSE, immediate. = TRUE)
    }
    if (verbose) {
      message(paste("\n=== Generating", type.data, "bulk samples:"))
    }
    bulk.counts <- .generateBulkProfiles(
      object = object,
      sim.counts = sim.counts,
      type.data = type.data,
      file.backend = file.backend,
      threads = threads,
      verbose = verbose
    )
    if (!is.null(object@bulk.counts)) {
      if (type.data %in% names(object@bulk.counts)) {
        object@bulk.sim[[type.data]] <- NULL
      }
      object@bulk.sim <- c(object@bulk.sim, type.data = bulk.counts)
    } else {
      object@bulk.sim[[type.data]] <- bulk.counts
    }
  }

  message("\nDONE")
  return(object)
}

setBulks <- function (x, c, i) {
  return(rowSums(c[, x]))
}

.generateBulkProfiles <- function(
  object,
  sim.counts,
  type.data,
  file.backend,
  threads,
  verbose
) {
  prob.matrix.names <- object@prob.matrix[[type.data]]@cell.names
  bulk.counts <- pbapply::pbapply(
    X = prob.matrix.names,
    MARGIN = 1,
    FUN = setBulks,
    c = sim.counts,
    cl = threads
  )
  colnames.bulk <- paste("Bulk", seq(dim(bulk.counts)[2]), sep = "_")
  rownames.bulk <- rownames(bulk.counts)
  if (!is.null(file.backend)) {
    if (verbose) {
      message("\nWriting data on disk:\n")
    }
    bulk.counts <- DelayedArray::DelayedArray(bulk.counts)
    bulk.counts <- HDF5Array::writeHDF5Array(
      bulk.counts,
      filepath = file.backend,
      name = type.data,
      verbose = verbose
    )
    # en R 4.0 implementaron la posibilidad de guardar dimnames
    # de momento lo meto en un slot
  }
  return(SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = bulk.counts),
    rowData = rownames.bulk,
    colData = colnames.bulk))
}

## generateCombinedScaledShuffledSet function ---------------------------------

generateCombinedScaledShuffledSet <- function(
  object,
  type.data,
  file.backend = NULL,
  verbose = TRUE
) {
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
    if (verbose) {
      message("=== Working in memory")
    }
    combineBulkSCProfiles <- .combineBulkSCProfilesInMemory
  } else {
    if (file.exists(file.backend)){
      stop("file.backend already exists. Please provide a correct file path")
    }
    if (verbose) {
      message("=== Working with HDF5 backend")
    }
    combineBulkSCProfiles <- .combineBulkSCProfilesHDF5
  }

  if (type.data == "both") {
    if (!is.null(object@final.data)) {
      warning("'final.data' slot will be overwritten",
              call. = FALSE, immediate. = TRUE)
    }
    combined.counts <- lapply(
      X = c("train", "test"),
      FUN = function(x) {
        if (verbose) {
          message(paste("\n=== Combining, scaling and shuffling", x, "counts\n"))
        }
        combineBulkSCProfiles(
          object = object,
          type.data = x,
          file.backend = file.backend,
          verbose = verbose
        )
      }
    )
    names(combined.counts) <- c("train", "test")
    object@final.data <- combined.counts
  } else {
    if (verbose) {
      message(paste("\n=== Combining, scaling and shuffling", type.data, "counts\n"))
    }
    combined.counts <- combineBulkSCProfiles(
      object = object,
      type.data = type.data,
      file.backend = file.backend,
      verbose = verbose
    )
    if (!is.null(object@bulk.counts)) {
      if (type.data %in% names(object@bulk.counts)) {
        object@bulk.sim[[type.data]] <- NULL
      }
      object@final.data <- c(object@final.data, type.data = combined.counts)
    } else {
      object@final.data <- list(type.data = combined.counts)
    }
  }

  message("\nDONE")
  return(object)
}


.setConfigHDF5 <- function(file.backend, name) {
  setHDF5DumpFile(file.backend)
  setHDF5DumpName(name)
}

.combineBulkSCProfilesHDF5 <- function(
  object = object,
  type.data = type.data,
  file.backend = file.backend,
  verbose = verbose
) {
  gene.list <- intersect(rowData(bulk.sim(object)[[type.data]])[[1]],
                         rownames(assay(object@single.cell.sim)))
  counts <- assay(bulk.sim(object)[[type.data]])
  # .setConfigHDF5(file.backend = file.backend, name = "tpm")
  counts <- DelayedArray::cbind(DelayedArray(assay(object@single.cell.sim)[,
                                 unlist(object@prob.matrix[[type.data]]@set.list)]),
                                counts)
  sample.names <- c(colnames(assay(object@single.cell.sim)[,
                                   unlist(object@prob.matrix[[type.data]]@set.list)]),
                    colData(object@bulk.sim[[type.data]])[[1]])
  tpsm <- matrix(unlist(sapply(
    X = names(object@prob.matrix[[type.data]]@set.list),
    FUN = function (x, l) {
      v <- rep(0,length(l))
      names(v) <- names(l)
      v[x] <- 100
      return(rep(v, length(l[[x]])))
    }, l = object@prob.matrix[[type.data]]@set.list
  )), ncol = length(object@prob.matrix[[type.data]]@set.list), byrow = T)
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
  lib.sizes <- DelayedArray::colSums(counts)
  prior.count.scaled <- 1L * length(lib.sizes) * lib.sizes / sum(lib.sizes)
  counts <- t(log2((t(counts) + prior.count.scaled) /
                     (lib.sizes +  2 * prior.count.scaled) * 1e+06))
  # scale
  c.means <- DelayedArray::colMeans(counts)
  c.sd <- DelayedMatrixStats::colSds(counts)
  counts <- (t(counts) - c.means) / c.sd
  s <- sample(seq(dim(counts)[1]))
  counts <- counts[s, ]
  sample.names <- sample.names[s]
  probs.matrix <- probs.matrix[s, ]
  if (verbose) {
    message("   Writing data on disk:\n")
  }
  counts <- HDF5Array::writeHDF5Array(
    x = counts,
    filepath = file.backend,
    name = type.data,
    verbose = verbose
  )

  return(SummarizedExperiment::SummarizedExperiment(
    assays = list(sim.counts = counts),
    rowData = sample.names,
    colData = gene.list,
    metadata = list(prob.matrix = probs.matrix)
  ))
}


.combineBulkSCProfilesInMemory <- function(
  object = object,
  type.data = type.data,
  file.backend = NULL,
  verbose = verbose
) {
  gene.list <- intersect(rowData(bulk.sim(object)[[type.data]])[[1]],
                         rownames(assay(object@single.cell.sim)))
  # check if bulk.sim is loaded in memory or in HDF5 file
  if (class(assay(bulk.sim(object)[[type.data]])) == "HDF5Matrix") {
    counts <- matrix(assay(bulk.sim(object)[[type.data]]))
  } else {
    counts <- assay(bulk.sim(object)[[type.data]])
  }
  counts <- cbind(
    assay(object@single.cell.sim)[, unlist(object@prob.matrix[[type.data]]@set.list)],
    counts
  )

  sample.names <- c(colnames(assay(object@single.cell.sim)[,
                             unlist(object@prob.matrix[[type.data]]@set.list)]),
                    colData(object@bulk.sim[[type.data]])[[1]])
  tpsm <- matrix(unlist(sapply(
    X = names(object@prob.matrix[[type.data]]@set.list),
    FUN = function (x, l) {
      v <- rep(0,length(l))
      names(v) <- names(l)
      v[x] <- 100
      return(rep(v, length(l[[x]])))
    }, l = object@prob.matrix[[type.data]]@set.list
  )), ncol = length(object@prob.matrix[[type.data]]@set.list), byrow = T)

  colnames(tpsm) <- names(object@prob.matrix[[type.data]]@set.list)
  tpsm <- tpsm[, colnames(object@prob.matrix[[type.data]]@prob.matrix)]
  rownames(tpsm) <- unlist(object@prob.matrix[[type.data]]@set.list)
  probs.matrix <- rbind(tpsm, object@prob.matrix[[type.data]]@prob.matrix)/100
  rownames(probs.matrix) <- c(rownames(tpsm), colData(object@bulk.sim[[type.data]])[[1]])

  # rownames(assay(DDLSChung.1@bulk.sim$train)) <- gene.list
  ## scale counts matrix
  # cpms
  counts <- edgeR::cpm.default(y = counts, log = TRUE, prior.count = 1)
  # scaling
  counts <- scale(counts)
  # shuffling
  s <- sample(seq(dim(counts)[2]))
  counts <- t(counts[, s])
  sample.names <- sample.names[s]
  probs.matrix <- probs.matrix[s, ]

  return(SummarizedExperiment::SummarizedExperiment(
    assays = list(sim.counts = counts),
    rowData = sample.names,
    colData = gene.list,
    metadata = list(prob.matrix = probs.matrix)
  ))
}
