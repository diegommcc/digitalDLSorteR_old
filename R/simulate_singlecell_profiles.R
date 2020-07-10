## Anotaciones:
## se puede hacer una genérica que funcione tanto para la clase DDLS como para SCE

estimateZinbwaveParams <- function(object,
                                   cell.ID.column,
                                   gene.ID.column,
                                   cell.type.column,
                                   cell.cov.columns,
                                   gene.cov.columns,
                                   set.type = "All",
                                   threads = 2,
                                   verbose = TRUE) {
  if (class(object) != "DigitalDLSorter") {
    stop("The object provided is not of DigitalDLSorter class")
  } else if (is.null(single.cell.real(object))) {
    stop("single.cell.real slot is empty")
  } else if (is.null(cell.cov.columns)) { ## perguntar si exigirlo
    stop("At least one covariate in cell.cov.column is required")
  }
  if (!is.null(zinb.params(object = object))) {
    warning("zinb.params slot already has a ZinbParams object. Note that it will be overwritten\n",
            call. = FALSE, immediate. = TRUE)
  }

  # extract data from SCE to list
  list.data <- .extractDataFromSCE(SCEobject = single.cell.real(object),
                                   filtering = FALSE)
  # check if params are correct
  mapply(function(x, y) {
    .checkColumn(metadata = list.data[[2]],
                 ID.column = x,
                 type.metadata = "cells.metadata",
                 arg = y)},
         c(cell.ID.column, cell.type.column, cell.cov.columns),
         c("cell.ID.column", "cell.type.column",
           rep("cell.cov.columns", length(cell.cov.columns))))

  lapply(gene.cov.columns, function(x) {
    .checkColumn(metadata = list.data[[3]],
                 ID.column = x,
                 type.metadata = "genes.metadata",
                 arg = "gene.cov.columns")
  })
  if (set.type != "All") {
    lapply(set.type, function(x) {
      .checkColumn(metadata = list.data[[2]],
                   ID.column = x,
                   type.metadata = "cells.metadata",
                   arg = "set.type")
    })
  }

  # check if covariates and cell types have at least two levels
  lapply(c(cell.type.column, cell.cov.columns),
         function(x) {
           if (length(unique(list.data[[2]][, x])) < 2) {
             stop(paste(x, "must have 2 or more unique elements"))
           }
         })
  lapply(gene.cov.columns,
         function(x) {
           if (length(unique(list.data[[3]][, x])) < 2) {
             stop(paste(x, "must have 2 or more unique elements"))
           }
         })

  ## preguntar por qué. Al cambiar los colnames, no puedes volver a
  ## ejecutar la función porque hay missing values debido a que el cell.ID
  ## ahora no es el mismo
  # rownames(list.data[[2]]) <- paste(list.data[[2]][, cell.type.column],
  #                                   list.data[[2]][, cell.ID.column],
  #                                   sep = "_")

  # set configuration of parallel computations
  if (threads <= 0) {
    threads <- 1
  }
  if (verbose) {
    message(paste("=== Set parallel environment to", threads, "threads\n"))
  }
  snowParam <- BiocParallel::SnowParam(workers = threads, type = "SOCK")

  if (set.type == "All") {
    list.data[[1]] <- as.matrix(list.data[[1]])
    list.data[[1]] <- list.data[[1]][rowSums(list.data[[1]]) > 0, ]
    # message(dim(counts))
    # message("Estimate parameters for experiment with model matrix")
    formula.cell.model <- as.formula(paste("~", paste(c(cell.cov.columns,
                                                   cell.type.column),
                                                 collapse = "+")))
    if (verbose) {
      message("=== Estimate parameters for the whole experiment\n")
      message(paste("=== Create cell model matrix based on", paste(cell.cov.columns,
                                                                   collapse = ", "),
                    "and", cell.type.column, "columns:"))
      message("\t", formula.cell.model, "\n")
    }
    # no coge siempre un tipo celular menos
    sdm <- model.matrix(formula.cell.model,
                        data = list.data[[2]][match(colnames(list.data[[1]]),
                                                    list.data[[2]][, cell.ID.column]), ])
    sdm.ncol <- ncol(sdm)
    sdm.colnames <- colnames(sdm)
    # message(dim(sdm))
    # message(head(sdm))
  } else {
      message(paste("=== Estimate parameters for", set.type, "from the experiment\n"))
      message(paste("=== Collect counts for", set.type, "cells\n"))
      cell.IDs <- list.data[[2]][which(list.data[[2]][, cell.type.column] == set.type),
                                 cell.ID.column]
      list.data [[1]] <- list.data[[1]][, cell.IDs]
      list.data[[1]] <- as.matrix(list.data[[1]])
      list.data[[1]] <- list.data[[1]][rowSums(list.data[[1]]) > 0,]
      # message(paste(c("Genes","Cells"), dim(list.data[[1]])))
      sdm <- NULL
      sdm.ncol <- 1
      sdm.colnames <- seq(1)
  }
  # covariates for genes
  if (!is.null(gene.cov.columns)) {
    formula.gene.model <- as.formula(paste("~", paste(gene.cov.columns, collapse = "+")))
    if (verbose) {
      message(paste("=== Create gene model matrix with", gene.cov.columns, "covariate(s)"))
      message("\t", formula.gene.model, "\n")
    }
    gdm <- stats::model.matrix(formula.gene.model,
                               data = list.data[[3]][match(rownames(list.data[[1]]),
                                                       list.data[[3]][, gene.ID.column]), ])
  } else {
    if (verbose) {
      message("=== Create gene model matrix without Covariates\n")
    }
    gdm <- stats::model.matrix(~1, data = list.data[[3]][match(rownames(list.data[[1]]),
                                                    list.data[[3]][, gene.ID.column]), ])
  }
  rownames(gdm) <- rownames(list.data[[1]])
  # message(dim(gdm))
  # message(head(gdm))
  start_time <- Sys.time()
  if (verbose) {
    message("=== Run estimation process ",
            paste("(Start time", format(start_time, "%X)"), "\n"))
  }
  zinbParamsObj <- zinbEstimate(ceiling(list.data[[1]])
                                , BPPARAM = snowParam
                                , design.samples = sdm
                                , design.genes = gdm
                                , O_mu = matrix(0, nrow = ncol(list.data[[1]]), ncol = nrow(list.data[[1]])
                                                , dimnames = list(rownames = seq(ncol(list.data[[1]]))
                                                                  ,colnames = rownames(list.data[[1]])
                                                )
                                )
                                , O_pi = matrix(0, nrow = ncol(list.data[[1]]), ncol = nrow(list.data[[1]])
                                                  , dimnames = list(rownames = seq(ncol(list.data[[1]]))
                                                                  ,colnames = rownames(list.data[[1]])
                                                )
                                )
                                , beta_mu = matrix(0, nrow = sdm.ncol, ncol = nrow(list.data[[1]])
                                                   , dimnames = list(rownames = sdm.colnames
                                                                     ,colnames = rownames(list.data[[1]])
                                                   )
                                )
                                , beta_pi = matrix(0, nrow = sdm.ncol, ncol = nrow(list.data[[1]])
                                                   , dimnames = list(rownames = sdm.colnames
                                                                     ,colnames = rownames(list.data[[1]])
                                                   )
                                )
                                , alpha_mu = matrix(0, nrow = 0, ncol = nrow(list.data[[1]])
                                                    , dimnames = list(rownames = NULL
                                                                      ,colnames = rownames(list.data[[1]])
                                                    )
                                )
                                , alpha_pi = matrix(0, nrow = 0, ncol = nrow(list.data[[1]])
                                                    , dimnames = list(rownames = NULL
                                                                      ,colnames = rownames(list.data[[1]])
                                                    )
                                )
                                , verbose = verbose
  )

  # update slots
  sce <- CreateSCEObject(counts = list.data[[1]],
                         cells.metadata = list.data[[2]],
                         genes.metadata = list.data[[3]])
  single.cell.real(object) <- sce
  zinb.params(object) <- zinbParamsObj

  end_time <- Sys.time()
  message("\nDONE\n")
  message(paste("Invested time:", round(end_time - start_time, 2), "mins"))
  return(object)
}

# function for check if slot is correct
.checkSlot <- function(object, slot) {
  ss <- eval(parse(text = paste0(slot, "(",
                   deparse(substitute(object)), ")")))
  if (is.null(ss)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}


simSingleCellProfiles <- function(object,
                                  cell.ID.column,
                                  cell.type.column,
                                  n.cells = 100,
                                  verbose = TRUE) {
  if (is.null(zinb.params(object))) {
    stop(paste("zinb.params slot is empty. To simulate single-cell profiles,",
               "DigitalDLSorter object  must contain the estimated parameters for the data with",
               "estimateZinbwaveParams function"))
  }
  if (is.null(single.cell.real(object))) {
    stop(paste("single.cell.real slot is empty. To simulate single-cell profiles,",
               "DigitalDLSorter object must contain the original data. See ?CreateDigitalDLSorterObject"))
  }
  if (!is.null(single.cell.sim(object = object))) {
    warning("single.cell.sim slot already has a SingleCellExperiment object. Note that it will be overwritten\n",
            call. = FALSE, immediate. = TRUE)
  }

  # extract data
  list.data <- .extractDataFromSCE(SCEobject = single.cell.real(object),
                                   filtering = FALSE)
  zinb.object <- zinb.params(object)
  # check if zinb.params and single.cell.real have the same dimensions
  # en el caso de que set.type sea != All, tendré que eliminar los tipos celulares
  # del cells.metadata. además hay un rowSums > 0 que puede cambiar el número de genes.
  # Lo mejor será quitar esas filas al cargar los datos
  dim.scr <- dim(list.data[[1]])
  if (!zinb.params(object)@nGenes == dim.scr[1] |
      !zinb.params(object)@nCells == dim.scr[2]) {
    stop("zinb.params slot and single.cell.real slot are not compatible,",
         "the number of dimensions are not equal")
  }
  # check if cell.type.column exists in cells.metadata
  .checkColumn(metadata = list.data[[2]],
               ID.column = cell.type.column,
               type.metadata = "cells.metadata",
               arg = "cell.type.column")
  if (n.cells < 5) {
    stop("n.cells must be greater than or equal to 5 cells per cell type")
  }
  # generate metadata for simulated cells
  colnames(list.data[[1]]) <- paste(list.data[[2]][, cell.type.column],
                                    list.data[[2]][, cell.ID.column],
                                    sep = "_")
  list.data[[2]]$simCellName <- paste(list.data[[2]][, cell.type.column],
                                      list.data[[2]][, cell.ID.column],
                                      sep = "_")
  list.data[[2]]$Simulated <- FALSE

  # cell types in model
  cell.set.names <- NULL
  model.cell.types <- grep(pattern = cell.type.column,
                           x = colnames(zinb.object@model@X),
                           value = T)
  cell.type.names <- sub(pattern = cell.type.column,
                         replacement = "",
                         x = model.cell.types)
  names(cell.type.names) <- model.cell.types

  if (verbose) {
    message(paste0("=== Cell Types in Model (", length(cell.type.names), " cell types):"))
    message(paste0("  - ", cell.type.names, collapse = "\n"))
    message()
  }

  for (s in model.cell.types) {
    cell.type.name <- cell.type.names[s]
    # message(paste(s, cell.type.name), "\n")
    cell.index <- rownames(zinb.object@model@X)[which(zinb.object@model@X[, s] == 1)]
    nams <- sample(cell.index, size = n.cells, replace = T)
    if (is.null(cell.set.names)) {
      cell.set.names <- nams
      names(cell.set.names) <- paste(cell.type.name, "_S",
                                     seq(from = 1, to = n.cells), sep = "")
    } else {
      ns <- names(cell.set.names)
      cell.set.names <- c(cell.set.names, nams)
      names(cell.set.names) <- c(ns, paste(cell.type.name, "_S",
              seq(from = length(ns) + 1, to = length(ns) + n.cells), sep = ""))
    }
  }

  # no entiendo esto
  inter.cell.type <- setdiff(levels(factor(list.data[[2]][, cell.type.column])),
                             cell.type.names)
  # cat(inter.cell.type, "\n")
  # To get the intercept cell type the rowSum of all FinalCellType columns should be 0
  cell.index <- rownames(zinb.object@model@X)[rowSums(zinb.object@model@X[,
                         grep(cell.type.column,
                         colnames(zinb.object@model@X), value = T)]) == 0]
  nams <- sample(cell.index, size = n.cells, replace = T)
  ns <- names(cell.set.names)
  cell.set.names <- c(cell.set.names, nams)
  names(cell.set.names) <- c(ns, paste(inter.cell.type, "_S",
                                       seq(from = length(ns) + 1,
                                           to = length(ns) + n.cells),
                                       sep = ""))
  n.cells <- length(cell.set.names)
  mu <- zinbwave::getMu(zinb.object@model) # rows are cells
  pi <- zinbwave::getPi(zinb.object@model) # rows are cells
  theta <- zinbwave::getTheta(zinb.object@model) # for genes

  if (verbose) {
    message("=== Get params from model:")
    message(paste("  - mu:", paste(dim(mu), collapse = ", ")))
    message(paste("  - pi:", paste(dim(pi), collapse = ", ")))
    message(paste("  - Theta:", length(theta)), "\n")
  }

  mu <- mu[cell.set.names, ]
  pi <- pi[cell.set.names, ]

  n <- length(cell.set.names) # as.numeric(nCells)
  J <- zinbwave::nFeatures(zinb.object@model)
  i <- seq(n * J)

  if (verbose) {
    message("=== Simulated Matrix dimensions:")
    message(paste("  - n (cells):", n))
    message(paste("  - J (genes):", J))
    message(paste("  - i (n entries):", length(i)), "\n")
  }

  datanb <- stats::rnbinom(length(i), mu = mu[i], size = theta[ceiling(i/n)])
  data.nb <- matrix(datanb, nrow = n)

  datado <- stats::rbinom(length(i), size=1, prob = pi[i])
  data.dropout <- matrix(datado, nrow = n)

  sim.counts <- data.nb * (1 - data.dropout)
  colnames(sim.counts) <- colnames(mu)
  sim.counts <- t(sim.counts)

  sim.counts <- Matrix::Matrix(sim.counts, dimnames = list(
    rownames = rownames(zinb.object@model@V),
    colnames = names(cell.set.names)))

  sim.cells.metadata <- list.data[[2]][cell.set.names, ]
  sim.cells.metadata$simCellName <- names(cell.set.names)
  rownames(sim.cells.metadata) <- names(cell.set.names)
  sim.cells.metadata$Simulated <- TRUE

  # new sim cells.metadata
  sim.cells.metadata <- rbind(list.data[[2]], sim.cells.metadata)

  # bind simmulated and real counts
  sim.counts <- Matrix::Matrix(as.matrix(sim.counts), sparse = T)
  sim.counts <- cbind(list.data[[1]][rownames(sim.counts), ], sim.counts)

  sim.sce <- CreateSCEObject(counts = sim.counts,
                             cells.metadata = sim.cells.metadata,
                             genes.metadata = list.data[[3]])
  single.cell.sim(object) <- sim.sce

  message("DONE")
  return(object)
}
