.readTabFiles <- function(file) {
  if (!file.exists(file)) {
    stop(paste(file, "file provided does not exists"))
  }
  if (grepl(".tsv", file)) {
    if (grepl(".tsv.gz$", file)) {
      file.obj <- read.delim(file = gzfile(file), sep = "\t",
                             header = T, stringsAsFactors = F)
    } else {
      file.obj <- read.delim(file = file, sep = "\t", header = T,
                             stringsAsFactors = F)
    }
  } else if (grepl(".rds$", file)) {
    file.obj <- readRDS(file = file)
  } else {
    stop("File format is not recognizable. Please look at the allowed data for in ?CreateDigitalDLSorterObject")
  }
  return(file.obj)
}


.readCountsFile <- function(counts.file, gene.column = 1) {
  if (grepl(".tsv", counts.file) | grepl(".rds$", counts.file)) {
    counts <- .readTabFiles(counts.file)
    return(counts)
  } else if (grepl(".mtx$", counts.file)) {
    if (!file.exists(counts.file)) {
      stop(paste(counts.file, "file provided does not exists"))
    }
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
    rownames(counts) <- gene.names[, gene.column]
    cell.names <- read.delim(file.path(base.dir, "barcodes.tsv"), header = F,
                             sep = "\t", stringsAsFactors = F)
    colnames(counts) <- cell.names$V1
    return(counts)
  } else {
    stop("File format is not recognizable. Please look at the allowed data for in ?CreateDigitalDLSorterObject")
  }
}


CreateSCEObject <- function(counts, cells.metadata, genes.metadata) {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = Matrix::Matrix(as.matrix(counts), sparse = TRUE)),
    colData = cells.metadata,
    rowData = genes.metadata
  )
  return(sce)
}


.checkColumn <- function(metadata, ID.column, type.metadata, arg) {
  tryCatch(expr = ID.column <- as.numeric(ID.column),
           error = function(e) invisible(x = NULL),
           warning = function(e) invisible(x = NULL))
  if (class(ID.column) == "numeric") {
    if (!ID.column %in% seq(ncol(metadata))) {
      stop(paste(ID.column, "column number is not present in", type.metadata))
    }
  } else if (class(ID.column) == "character") {
    if (!ID.column %in% colnames(metadata)) {
      stop(paste(ID.column, "column is not present in", type.metadata))
    }
  } else {
    stop(paste(arg, "argument is not recognizable"))
  }
}


# .duplicatedIDs <- function(counts, cells.metadata, cell.ID.column,
#                            genes.metadata, gene.ID.column) {
#   ## duplicated ID cells
#   if (any(duplicated(cells.metadata[, cell.ID.column]))) {
#     warning(paste0("There are duplicated IDs in genes.metadata (",
#                    gene.ID.column, ") column. Making unique"))
#     cells.metadata[, cell.ID.column] <- make.unique(names = cells.metadata[, cell.ID.column])
#   }
#   ## duplicated ID genes
#   if (any(duplicated(genes.metadata[, gene.ID.column]))) {
#     warning(paste0("There are duplicated IDs in genes.metadata (",
#                    gene.ID.column, ") column. Removing duplicated IDs"))
#     genes.metadata <- genes.metadata[!duplicated(genes.metadata[, gene.ID.column]), ]
#   }
# }


.processData <- function(counts, cells.metadata, cell.ID.column,
                         genes.metadata, gene.ID.column,
                         min.counts, min.cells) {
  # check if IDs given exist in metadata
  .checkColumn(metadata = cells.metadata,
               ID.column = cell.ID.column,
               type.metadata = "cells.metadata",
               arg = "cell.ID.column")
  .checkColumn(metadata = genes.metadata,
               ID.column = gene.ID.column,
               type.metadata = "genes.metadata",
               arg = "gene.ID.column")
  # duplicated ID cells --------------------------------------------------------
  if (any(duplicated(cells.metadata[, cell.ID.column]))) {
    warning(paste0("There are duplicated IDs in genes.metadata (column ",
                   gene.ID.column, "). Making unique"))
    cells.metadata[, cell.ID.column] <- make.unique(names = cells.metadata[, cell.ID.column])
  }
  # intersect between cells ----------------------------------------------------
  common.cells <- intersect(colnames(counts), cells.metadata[, cell.ID.column])
  diff <- abs(dim(counts)[2] - length(common.cells))
  disc <- abs(length(cells.metadata[, cell.ID.column]) - length(common.cells))
  if (length(common.cells) < min(dim(counts)[2], dim(cells.metadata)[1])) {
    stop(paste("There are", diff, "cells that don't match between counts matrix and metadata"))
  } else if (diff != 0){
    warning(paste("There are", diff, "cells that don't match between counts matrix and metadata")) # this check includes the last
  } else if (disc != 0) {
    message("=== Intersection between matrix counts and cells.metadata:")
    message(paste("   ", disc, "cells have been discarded from cells.metadata"), "\n")
  }
  cells.metadata <- cells.metadata[cells.metadata[, cell.ID.column] %in%
                                     common.cells, , drop = FALSE]
  # duplicated ID genes --------------------------------------------------------
  if (any(duplicated(genes.metadata[, gene.ID.column]))) {
    message("=== Removing duplicated genes:")
    message(paste0("   There are duplicated IDs in genes.metadata (column ",
                   gene.ID.column, "). Duplicated IDs have been removed"), "\n")
    genes.metadata <- genes.metadata[!duplicated(genes.metadata[, gene.ID.column]), ]
  }
  # intersect between genes ----------------------------------------------------
  common.genes <- intersect(rownames(counts), genes.metadata[, gene.ID.column])
  diff <- abs(dim(counts)[1] - length(common.genes))
  disc <- abs(length(genes.metadata[, gene.ID.column]) - length(common.genes))
  if (length(common.genes) < min(dim(counts)[1], dim(genes.metadata)[1])) {
    stop(paste("There are", diff, "genes that don't match between counts matrix and metadata"))
  } else if (diff != 0){
    warning(paste("There are", diff, "genes that don't match between counts matrix and metadata"))
  } else if (disc != 0) {
    message("=== Intersection between matrix counts and genes.metadata:")
    message(paste("   ", disc, "genes have been discarded from genes.metadata"), "\n")
  }
  genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                     common.genes, , drop = FALSE]
  counts <- counts[common.genes, common.cells]

  # filter genes by min.counts and min.cells -----------------------------------
  filtered.genes <- .filterGenes(counts = counts,
                                 genes.metadata = genes.metadata,
                                 gene.ID.column = gene.ID.column,
                                 min.counts = min.counts,
                                 min.cells = min.cells)

  return(list(filtered.genes[[1]], cells.metadata, filtered.genes[[2]]))
}


.filterGenes <- function(counts, genes.metadata, gene.ID.column,
                         min.counts, min.cells) {
  if (min.counts == 0 & min.cells == 0) {
    return(list(counts, genes.metadata))
  } else if (min.counts < 0 | min.cells < 0) {
    stop("min.counts and min.cells must be greater than or equal to zero")
  }
  dim.bef <- dim(counts)
  counts <- counts[rowSums(as.matrix(counts) > min.counts) >= min.cells, ]
  if (dim(counts)[1] == 0) {
    stop(paste("Resulting counts matrix after filtering with min.genes = ",
               min.counts, "and min.cells = ", min.cells,
               "does not have entries"))
  }
  message("\n=== Filtering features by min.counts and min.cells:")
  message(paste("  - Selected features:",  dim(counts)[1]))
  message(paste("  - Discarded features:", dim.bef[1] - dim(counts)[1]))

  genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                     rownames(counts), , drop = FALSE]

  return(list(counts, genes.metadata))
}


.extractDataFromSCE <- function(SCEobject, cell.ID.column = 1, gene.ID.column = 1,
                                min.counts = 0, min.cells = 0, filtering = TRUE) {
  # extract cells.metadata
  cells.metadata <- SingleCellExperiment::colData(SCEobject)
  if (any(dim(cells.metadata) == 0)) {
    stop("No data provided in colData slot. Metadata about cells is needed, please look ?CreateDigitalDLSorterObject")
  }
  # extract count matrix
  if (length(SummarizedExperiment::assays(SCEobject)) == 0) {
    stop("No data in SingleCellExperiment object provided")
  } else if (length(SummarizedExperiment::assays(SCEobject)) > 1) {
    warning("There are more than one assay, only the first will be used. Remember it must be the original data and not log-transformed data")
  }
  counts <- SummarizedExperiment::assay(SCEobject)
  # extract genes.metadata
  genes.metadata <- SingleCellExperiment::rowData(SCEobject)
  if (any(dim(genes.metadata) == 0)) {
    warning("No data provided in rowData slot. Building a rowData from rownames of counts matrix")
    if (class(gene.ID.column) == "numeric") gene.ID.column <- "gene_names"
    genes.metadata <- DataFrame(gene.ID.column = rownames(counts))
  }

  # check if IDs given exist in genes.metadata. In cells.metadata is not
  # neccesary because the data are provided from an SCE object
  .checkColumn(metadata = genes.metadata,
               ID.column = gene.ID.column,
               type.metadata = "genes.metadata",
               arg = "gene.ID.column")

  # filter genes by min.counts and min.cells only when proccess data
  if (isTRUE(filtering)) {
    filtered.genes <- .filterGenes(counts = counts,
                                   genes.metadata = genes.metadata,
                                   gene.ID.column = gene.ID.column,
                                   min.counts = min.counts,
                                   min.cells = min.cells)
    return(list(filtered.genes[[1]], cells.metadata, filtered.genes[[2]]))
  } else {
    return(list(counts, cells.metadata, genes.metadata))
  }
}


.loadSingleCellData <- function(single.cell, cell.ID.column, gene.ID.column,
                                min.cells, min.counts, real = TRUE) {
  # check if single-cell data are real or simulated
  if (isTRUE(real)) {
    arg <- "single.cell.real"
    if (is.null(single.cell)) {
      stop("single.cell.real cannot be NULL")
    } else if (is.null(cell.ID.column) | is.null(gene.ID.column)) {
      stop("cell.ID.column and gene.ID.column are needed. Please look ?CreateDigitalDLSorterObject")
    }
  } else {
    arg <- "single.cell.sim"
    if (is.null(single.cell)) return(NULL)
  }
  # load data from the allowed sources
  if (class(single.cell) == "SingleCellExperiment") {
    # extract data and filter by min.counts and min.cells
    list.data <- .extractDataFromSCE(SCEobject = single.cell,
                                     cell.ID.column = cell.ID.column,
                                     gene.ID.column = gene.ID.column,
                                     min.counts = min.counts,
                                     min.cells = min.cells)
    single.cell <- CreateSCEObject(counts = list.data[[1]],
                                   cells.metadata = list.data[[2]],
                                   genes.metadata = list.data[[3]])
    return(single.cell)
  } else if (length(single.cell) == 0) {
    stop(paste(arg, "argument is empty"))
  } else if (length(single.cell) == 3) {
    list.data <- list(.readCountsFile(single.cell[[1]]),
                      .readTabFiles(single.cell[[2]]),
                      .readTabFiles(single.cell[[3]]))
  } else {
    stop(paste("Incorrect number of data elements given. Please look at the allowed data for",
               arg, "in ?CreateDigitalDLSorterObject"))
  }
  # process data only for real single-cell profiles from files (not SCE)
  if (isTRUE(real)) {
    list.data <- .processData(counts = list.data[[1]],
                              cells.metadata = list.data[[2]],
                              cell.ID.column = cell.ID.column,
                              genes.metadata = list.data[[3]],
                              gene.ID.column = gene.ID.column,
                              min.counts = min.counts,
                              min.cells = min.cells)
  }
  return(CreateSCEObject(counts = list.data[[1]],
                         cells.metadata = list.data[[2]],
                         genes.metadata = list.data[[3]]))
}



loadRealSCProfiles <- function(
  single.cell.real,
  cell.ID.column = 1,
  gene.ID.column = 1,
  min.cells = 0,
  min.counts = 0,
  project = "DigitalDLSorterProject"
) {
  single.cell.real <- .loadSingleCellData(
    single.cell = single.cell.real,
    cell.ID.column = cell.ID.column,
    gene.ID.column = gene.ID.column,
    min.cells = min.cells,
    min.counts = min.counts
  )
  ddls.object <- new(
    Class = "DigitalDLSorter",
    single.cell.real = single.cell.real,
    project = project,
    version = packageVersion(pkg = "digitalDLSorteR")
  )
  return(ddls.object)
}


