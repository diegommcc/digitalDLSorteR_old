## Anotaciones:
## se puede hacer una genérica que funcione tanto para la clase DDLS como para SCE

estimateZinbwaveParams <- function(object,
                                   cell.ID.column,
                                   gene.ID.column,
                                   cell.type.column,
                                   cell.cov.columns,
                                   gene.cov.columns,
                                   set.type = "All",
                                   threads = 2) {
  if (class(object) != "DigitalDLSorter") {
    stop("The object provided is not of DigitalDLSorter class")
  } else if (is.null(single.cell.real(object))) {
    stop("single.cell.real slot si empty")
  }
  # extract data from SCE to list
  list.data <- .extractDataFromSCE(SCEobject = single.cell.real(object),
                                   filtering = FALSE)
  # check if params are correct
  .checkColumn(metadata = list.data[[2]],
               ID.column = cell.ID.column,
               type.metadata = "cells.metadata",
               arg = "cell.ID.column")
  lapply(cell.cov.columns, function(x) {
    .checkColumn(metadata = list.data[[2]],
                 ID.column = x,
                 type.metadata = "cells.metadata",
                 arg = "cell.cov.columns")
  })
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

  rownames(list.data[[2]]) <- paste(list.data[[2]][, cell.type.column],
                                    list.data[[2]][, cell.ID.column],
                                    sep = "_")
  # set configuration of parallel computations
  if (threads <= 0) {
    threads <- 1
  }
  message(paste("=== Set parallel environment to", threads, "threads\n"))
  snowParam <- BiocParallel::SnowParam(workers = threads, type = "SOCK")

  # estimation of parameters for all cell types
  if (set.type == "All") {
    message("=== Estimate parameters for the whole experiment\n")
    list.data[[1]] <- list.data[[1]][rowSums(as.matrix(list.data[[1]])) > 0, ]
    # message(dim(counts))
    # message("Estimate parameters for experiment with model matrix")
    message(paste("=== Create cell model matrix based on:", paste(cell.cov.columns,
                                                         collapse = ", "),
                  "and", cell.type.column))
    formula.cell.model <- as.formula(paste("~", paste(c(cell.cov.columns,
                                                   cell.type.column),
                                                 collapse = "+")))
    message("\t", formula.cell.model)
    sdm <- model.matrix(formula.cell.model,
                        data = list.data[[2]][match(colnames(list.data[[1]]),
                                                    list.data[[2]][, cell.ID.column]), ])
    sdm.ncol <- ncol(sdm)
    sdm.colnames <- colnames(sdm)
    # message(dim(sdm))
    # message(head(sdm))
  } else {
      message(paste("Estimate parameters for", set.type, "from the experiment\n"))
      message(paste("Collect counts for", set.type, "cells"))
      cell.IDs <- list.data[[2]][which(list.data[[2]][, cell.type.column] == set.type), cell.ID.column]
      list.data [[1]] <- list.data[[1]][, cell.IDs]
      # counts <- as.matrix(counts)
      counts <- list.data[[1]][rowSums(as.matrix(list.data[[1]])) > 0,]
      message(paste(c("Genes","Cells"), dim(list.data[[1]])))
      sdm <- NULL
      sdm.ncol <- 1
      sdm.colnames <- seq(1)
  }

  if (!is.null(gene.cov.columns)) {
    message(paste("=== Create gene model matrix with", gene.cov.columns, "covariates"))
    formula.gene.model <- as.formula(paste("~", paste(gene.cov.columns, collapse = "+")))
    gdm <- model.matrix(formula.gene.model,
                        data = list.data[[3]][match(rownames(list.data[[1]]),
                                                list.data[[3]][, gene.ID.column]), ])
  } else {
    message("Create gene model matrix without Covariates\n")
    gdm <- model.matrix(~1, data = list.data[[3]][match(rownames(list.data[[1]]),
                                                    list.data[[3]][, gene.ID.column]), ])
  }
  rownames(gdm) <- rownames(list.data[[1]])
  message(dim(gdm))
  message(head(gdm))

  message("Run estimation process")
  start_time <- Sys.time()
  message(format(start_time, "%X"))
  # zinbParamsAll <- zinbEstimate(ceiling(counts)
  #                               , BPPARAM = snowParam
  #                               , design.samples = sdm
  #                               , design.genes = gdm
  #                               , O_mu = matrix(0, nrow = ncol(counts), ncol = nrow(counts)
  #                                               , dimnames = list(rownames = seq(ncol(counts))
  #                                                                 ,colnames = rownames(counts)
  #                                               )
  #                               )
  #                               , O_pi = matrix(0, nrow = ncol(counts), ncol = nrow(counts)
  #                                                 , dimnames = list(rownames = seq(ncol(counts))
  #                                                                 ,colnames = rownames(counts)
  #                                               )
  #                               )
  #                               , beta_mu = matrix(0, nrow = sdm.ncol, ncol = nrow(counts)
  #                                                  , dimnames = list(rownames = sdm.colnames
  #                                                                    ,colnames = rownames(counts)
  #                                                  )
  #                               )
  #                               , beta_pi = matrix(0, nrow = sdm.ncol, ncol = nrow(counts)
  #                                                  , dimnames = list(rownames = sdm.colnames
  #                                                                    ,colnames = rownames(counts)
  #                                                  )
  #                               )
  #                               , alpha_mu = matrix(0, nrow = 0, ncol = nrow(counts)
  #                                                   , dimnames = list(rownames = NULL
  #                                                                     ,colnames = rownames(counts)
  #                                                   )
  #                               )
  #                               , alpha_pi = matrix(0, nrow = 0, ncol = nrow(counts)
  #                                                   , dimnames = list(rownames = NULL
  #                                                                     ,colnames = rownames(counts)
  #                                                   )
  #                               )
  # )
  end_time <- Sys.time()
  message("DONE\n")
  message(end_time - start_time)
}
