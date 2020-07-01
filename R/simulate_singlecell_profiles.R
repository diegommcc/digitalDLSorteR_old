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
    stop("single.cell.real slot si empty")
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
  mapply(function(x, y) {.checkColumn(metadata = list.data[[2]],
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

  ## preguntar por qué
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
      cell.IDs <- list.data[[2]][which(list.data[[2]][, cell.type.column] == set.type), cell.ID.column]
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
    gdm <- model.matrix(formula.gene.model,
                        data = list.data[[3]][match(rownames(list.data[[1]]),
                                                list.data[[3]][, gene.ID.column]), ])
  } else {
    if (verbose) {
      message("=== Create gene model matrix without Covariates\n")
    }
    gdm <- model.matrix(~1, data = list.data[[3]][match(rownames(list.data[[1]]),
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
  zinbParamsAll <- zinbEstimate(ceiling(list.data[[1]])
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
  zinb.params(object) <- zinbParamsAll
  end_time <- Sys.time()
  message("DONE\n")
  message(paste("Invested time:", round(end_time - start_time, 2), "mins"))
  return(object)
}
