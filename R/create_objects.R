## dependencies: these packages will be introduced in NAMESPACE file
library(Matrix)

# hay funciones que tienen una mejor performance respecto a la lectura de
# ficheros de texto, igual pueden ser una opción, aunque creo que para estos
# ficheros no es necesario

.readTabFiles <- function(file) {
  if (!file.exists(file)) {
    stop("File provided does not exists")
  }
  if (grepl(".tsv", file)) {
    if (grepl(".tsv.gz$", file)) {
      # message("Tab Gzipped format")
      file.obj <- read.delim(file = gzfile(file), sep = "\t",
                             header = T, stringsAsFactors = F)
    } else {
      # message("Tab plain format")
      file.obj <- read.delim(file = file, sep = "\t", header = T,
                             stringsAsFactors = F)
    }
  } else if (grepl(".rds$", file)) {
    # message("RDS R object")
    file.obj <- readRDS(file = file)
  }
  return(file.obj)
}


.readCountsFile <- function(counts.file) {
  if (grepl(".tsv", counts.file) | grepl(".rds$", counts.file)) {
    counts <- .readTabFiles(counts.file)
    return(counts)
  } else if (grepl(".mtx$", counts.file)) {
    if (!file.exists(counts.file)) {
      stop("File provided does not exists")
    }
    message("Sparse Matrix Market format")
    base.dir <- dirname(counts.file)
    if (!file.exists(file.path(base.dir, "genes.tsv"))) {
      stop("No genes.tsv file with the matrix.mtx")
    }
    if (!file.exists(file.path(base.dir, "barcodes.tsv"))) {
      stop("No barcodes.tsv file with the matrix.mtx")
    }
    counts <- Matrix::readMM(counts.file)
    gene.names <- read.delim(file.path(base.dir, "genes.tsv"), header = F,
                             sep = "\t", stringsAsFactors = F)
    rownames(counts) <- gene.names$V1
    cell.names <- read.delim(file.path(base.dir, "barcodes.tsv"), header = F,
                             sep = "\t", stringsAsFactors = F)
    colnames(counts) <- cell.names$V1
    return(list(counts, cell.names, gene.names))
  }
}


CreateSCEObject <- function(counts, cells.metadata, genes.metadata) {
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = cells.metadata,
    rowData = genes.metadata
  )
  return(sce)
}

## modularizaría más la función, pero no estoy seguro si cada vez que paso
## un objeto a una función se duplica en memoria
.processData <- function(counts, cells.metadata, cell.ID.column,
                         genes.metadata, gene.ID.column,
                         min.counts, min.cells) {
  # intersect between cells ----------------------------------------------------
  common.cells <- intersect(colnames(counts), cells.metadata[, cell.ID.column])
  diff <- abs(dim(counts)[2] - length(common.cells))
  disc <- abs(length(cells.metadata[, cell.ID.column]) - length(common.cells))
  if (length(common.cells) < min(dim(counts)[2], dim(cells.metadata)[1]))
    stop(paste("There are", diff, "cells that don't match between counts matrix and metadata"))
  else if (diff != 0){
    warning(paste("There are", diff, "cells that don't match between counts matrix and metadata")) # preguntar
  } else if (disc != 0) {
    message(paste(disc, "cells have been discarded from cells.metadata"))
  }
  cells.metadata <- cells.metadata[cells.metadata[, cell.ID.column] %in%
                                     common.cells, , drop = FALSE]

  # intersect between genes ----------------------------------------------------
  common.genes <- intersect(rownames(counts), genes.metadata[, gene.ID.column])
  diff <- abs(dim(counts)[1] - length(common.genes))
  disc <- abs(length(genes.metadata[, gene.ID.column]) - length(common.genes))
  if (length(common.genes) < min(dim(counts)[1], dim(genes.metadata)[1])) {
    stop(paste("There are", diff, "genes that don't match between counts matrix and metadata"))
  } else if (diff != 0){
    warning(paste("There are", diff, "genes that don't match between counts matrix and metadata"))
  } else if (disc != 0) {
    message(paste(disc, "genes have been discarded from genes.metadata"))
  }
  genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                     common.genes, , drop = FALSE]
  counts <- counts[common.genes, common.cells]

  # filter genes by min.counts and min.cells -----------------------------------
  if (min.counts == 0 & min.cells == 0) {
    return(list(counts, cells.metadata, genes.metadata))
  } else if (min.counts < 0 | min.cells < 0) {
    stop("min.counts and min.cells must be greater than or equal to zero")
  }
  dim.bef <- dim(counts)
  counts <- counts[rowSums(counts > min.counts) >= min.cells, ]
  if (dim(counts)[1] == 0 | dim(counts)[2] == 0) {
    stop(paste("Resulting counts matrix after filtering with min.genes = ",
               min.counts, "and min.cells = ", min.cells,
               "does not have entries"))
  }
  message("=== Filtering genes by min.counts and min.cells:")
  message(paste("Selected genes:",  dim(counts)[1]))
  message(paste("Discarded genes:", dim.bef[1] - dim(counts)[1]))

  genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                     rownames(counts), , drop = FALSE]

  return(list(counts, cells.metadata, genes.metadata))
}



# Con esta función asumo que las sparse matrix están bien construídas
# y el número de elementos en filas y columnas es el mismo que el de los genes y cells
.loadSingleCellData <- function(single.cell, cell.ID.column, gene.ID.column,
                                min.cells, min.counts, real = TRUE) {
  ## check if sce real or simulated
  if (isTRUE(real)) {
    arg <- "single.cell.real"
    if (is.null(single.cell)) {
      stop("single.cell.real cannot be NULL")
    }
  } else {
    arg <- "single.cell.sim"
    if (is.null(single.cell)) {
      return(NULL)
    }
  }
  if (class(single.cell) == "SingleCellExperiment") {
    ## hacer checkeos en el objeto SingleCellExperiment
    return(single.cell)
  } else {
    if (class(single.cell) == "list") {
      if (length(single.cell) == 0) {
        stop("List object is empty")
      } else if (length(single.cell) == 1 & grepl(".mtx", single.cell[[1]])) {
        list.data <- .readCountsFile(single.cell[[1]])
      } else if (length(single.cell) == 3) {
        if (is.null(cell.ID.column) | is.null(gene.ID.column)) {
          stop(paste("cell.ID.column and gene.ID.column are needed if you provide separated files. Please look at the allowed data for",
                     arg, "in ?CreateDigitalDLSorterObject"))
        } else {
          list.data <- .processData(
            counts = .readCountsFile(single.cell[[1]]),
            cells.metadata = .readTabFiles(single.cell[[2]]),
            cell.ID.column = cell.ID.column,
            genes.metadata = .readTabFiles(single.cell[[3]]),
            gene.ID.column = gene.ID.column,
            min.counts = min.counts,
            min.cells = min.cells
          )
        }
      } else {
        stop(paste("Incorrect data elements given. Please look at the allowed data for",
                   arg, "in ?CreateDigitalDLSorterObject"))
      }
    } else if (class(single.cell) == "character" & grepl(".mtx$", single.cell)) {
      ## checkear que la matriz sparse está bien construída
      list.data <- .readCountsFile(single.cell)
    } else {
      stop(paste(arg, "argument is not recognizable. Please look at the allowed data for",
                 arg, "in ?CreateDigitalDLSorterObject"))
    }
    list.data <- .processData(
      counts = list.data[[1]],
      cells.metadata = list.data[[2]],
      cell.ID.column = cell.ID.column,
      genes.metadata = list.data[[3]],
      gene.ID.column = gene.ID.column,
      min.counts = min.counts,
      min.cells = min.cells
    )
  }
  single.cell <- CreateSCEObject(list.data[[1]], list.data[[2]], list.data[[3]])
  return(single.cell)
}


CreateDigitalDLSorterObject <- function(
  single.cell.real,
  cell.ID.column = NULL,
  gene.ID.column = NULL,
  min.cells = 0,
  min.counts = 0,
  project = "DigitalDLSorterProject",
  zinb.params.object = NULL,
  selected.genes = NULL,
  single.cell.sim = NULL
) {
  single.cell.real <- .loadSingleCellData(
    single.cell.real,
    cell.ID.column,
    gene.ID.column,
    min.cells,
    min.counts
  )
  single.cell.sim <- .loadSingleCellData(
    single.cell.sim,
    real = FALSE,
    cell.ID.column,
    gene.ID.column,
    min.cells = 0,
    min.counts = 0
  )

  ddls.object <- new(
    Class = "DigitalDLSorter",
    single.cell.real = single.cell.real,
    zinb.params = zinb.params.object,
    selected.genes = "character",
    single.cell.sim = single.cell.sim,
    project = "DigitalDLSorterProject",
    version = packageVersion(pkg = "digitalDLSorterPackageR")
  )
  return(ddls.object)
}

