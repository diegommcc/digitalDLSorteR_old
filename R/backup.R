## backup of code to resort in case of problems

# Test file: I want to test if I can add S4 classes as slots in the main class.
# I can do it. Therefore, I can use SCE for the single-cell profiles and even
# SummaryExperiment class for Bulk profiles.

library(SingleCellExperiment)
library(splatter)

DigitalDLSorter <- setClass(
  Class = "DigitalDLSorter",
  slots = c(
    project = "character",
    single.cell.real = "SingleCellExperimentOrNULL",
    zinb.params = "ZINBParamsOrNULL",
    selected.genes = "character",
    single.cell.sim = "SingleCellExperimentOrNULL"
  )
)

setClassUnion("SingleCellExperimentOrNULL", c("SingleCellExperiment", NULL))
setClassUnion("ZINBParamsOrNULL", c("ZINBParams", NULL))


.loadCountsFile <- function(
  countsFile
) {
  if (grepl(".tsv", countsFile)) {
    if (grepl(".tsv.gz", countsFile)) {
      cat("\tTab Gzipped format\n")
      counts <- read.delim(file = gzfile(countsFile), sep = "\t",
                           header = T, stringsAsFactors = F)
    } else {
      cat("\tTab plain format\n")
      counts <- read.delim(file = countsFile, sep = "\t", header = T,
                           stringsAsFactors = F)
    }
  } else if (grepl(".rds", countsFile)) {
    cat("\tRDS R object\n")
    counts <- readRDS(file = countsFile)
  } else if (grepl(".mtx", countsFile)) {
    cat("\tSparse Matrix Market format\n")
    baseDir <- dirname(countsFile)
    if (!file.exists(file.path(baseDir,"genes.tsv"))) {
      stop("No genes.tsv file with the matrix.mtx")
    }
    if (!file.exists(file.path(baseDir,"barcodes.tsv"))) {
      stop("No barcodes.tsv file with the matrix.mtx")
    }
    counts <- readMM(countsFile)
    geneNames <- read.delim(file.path(baseDir,"genes.tsv"),header = F,sep = "\t",stringsAsFactors = F)
    rownames(counts) <- geneNames$V1
    cellNames <- read.delim(file.path(baseDir,"barcodes.tsv"),header = F,sep = "\t",stringsAsFactors = F)
    colnames(counts) <- cellNames$V1


  }
  return(counts)
}


.loadGenesMetadata <- function(
  genesMetadata
) {
  cat("Load Genes Metadata")
  if (grepl(".tsv", genesMetadata)) {
    if (grepl(".tsv.gz", genesMetadata)) {
      cat("\tTab Gzipped format\n")
      genesAnnot <- read.delim(file = gzfile(genesMetadata), sep = "\t",
                               header = T, stringsAsFactors = F)
    } else {
      cat("\tTab plain format\n")
      genesAnnot <- read.delim(file = genesMetadata, sep = "\t",header = T,
                               stringsAsFactors = F)
    }
  } else if (grepl(".rds", genesMetadata)) {
    cat("\tRDS R object\n")
    genesAnnot <- readRDS(file = genesMetadata)
  }
}


.loadCellsMetadata <- function(
  cellsMetadataFile
) {
  cat("Load Cells Metadata\n")
  if (grepl(".tsv",cellsMetadataFile)) {
    if (grepl(".tsv.gz",cellsMetadataFile)) {
      cat("\tTab Gzipped format\n")
      cellsMetadata <- read.delim(file = gzfile(cellsMetadataFile),sep = "\t",header = T,stringsAsFactors = F)
      head(cellsMetadata)
    } else {
      cat("\tTab plain format\n")
      cellsMetadata <- read.delim(file = cellsMetadataFile,sep = "\t",header = T,stringsAsFactors = F)
    }
  } else if (grepl(".rds",cellsMetadataFile)) {
    cat("\tRDS R object\n")
    cellsMetadata <- readRDS(file = cellsMetadataFile)
  }
}


CreateSCEObject <- function(
  counts,
  rowData,
  colData
) {



}

.checkSCEobject <- function(single.cell, real = TRUE) {
  ## check if sce real or simulated
  if (isTRUE(real)) {
    if (is.null(single.cell)) {
      stop("single.cell.real cannot be NULL.")
    }
  } else {
    if (is.null(single.cell)) {
      return(NULL)
    }
  }
  if (class(single.cell.real) == "list") {
    if (length(single.cell.real) == 0) {
      stop("List object is empty.")
    } else if (length(single.cell.real) == 1) {
      list.data <- .loadCountsFile(single.cell.real[[1]])
    } else if (length(single.cell.real) == 3) {
      list.data <- list(.loadCountsFile(single.cell.real[[1]]),
                        .loadCellsMetadata(single.cell.real[[2]]),
                        .loadGenesMetadata(single.cell.real[[3]]))
    }
    SCEObject <- CreateSCEObject(list.data[[1]], list.data[[2]], list.data[[3]])

  } else if (class(single.cell.real) == "SingleCellExperiment") {
    SCEObject <- single.cell.real # esto está mal porque es copiar dos veces el mismo objeto
  } else {
    stop("single.cell.real is not recognizable. Please please look at the allowed data for single.cell.real in ?CreateDigitalDLSorterObject.")
  }
  return(SCEObject)
}


CreateDigitalDLSorterObject <- function(
  single.cell.real,
  project = "DigitalDLSorterProject",
  zimb.params.object = NULL,
  selected.genes = NULL,
  single.cell.sim = NULL
) {
  if (is.null(single.cell.real)) {
    stop("single.cell.real cannot be NULL.")
  } else if (class(single.cell.real) == "list") {
    if (length(single.cell.real) == 0) {
      stop("List object is empty.")
    } else if (length(single.cell.real) == 1) {
      list.data <- .loadCountsFile(single.cell.real[[1]])
    } else if (length(single.cell.real) == 3) {
      list.data <- list(.loadCountsFile(single.cell.real[[1]]),
                        .loadCellsMetadata(single.cell.real[[2]]),
                        .loadGenesMetadata(single.cell.real[[3]]))
    }
    SCEObject <- CreateSCEObject(list.data[[1]], list.data[[2]], list.data[[3]])

  } else if (class(single.cell.real) == "SingleCellExperiment") {
    SCEObject <- single.cell.real # esto está mal porque es copiar dos veces el mismo objeto
  } else {
    stop("single.cell.real is not recognizable. Please please look at the allowed data for single.cell.real in ?CreateDigitalDLSorterObject.")
  }




}



## 25/06/20 --> 00:35

## dependencies: these packages will be introduced in NAMESPACE file
library(Matrix)

# hay funciones que tienen una mejor performance respecto a la lectura de
# ficheros de texto, igual pueden ser una opción, aunque creo que para estos
# ficheros no es necesario

.readTabFiles <- function(filePath) {
  if (!file.exists(filePath)) {
    stop("File provided does not exists")
  }
  if (grepl(".tsv", filePath)) {
    if (grepl(".tsv.gz", filePath)) {
      cat("\tTab Gzipped format\n")
      file <- read.delim(file = gzfile(filePath), sep = "\t",
                         header = T, stringsAsFactors = F)
    } else {
      cat("\tTab plain format\n")
      file <- read.delim(file = filePath, sep = "\t",header = T,
                         stringsAsFactors = F)
    }
  } else if (grepl(".rds", filePath)) {
    cat("\tRDS R object\n")
    file <- readRDS(file = filePath)
  }
  return(file)
}


.readCountsFile <- function(countsFile) {
  if (grepl(".tsv", countsFile) | grepl(".rds", countsFile)) {
    counts <- .readTabFiles(countsFile)
    return(counts)
  } else if (grepl(".mtx", countsFile)) {
    if (!file.exists(countsFile)) {
      stop("File provided does not exists")
    }
    cat("\tSparse Matrix Market format\n")
    baseDir <- dirname(countsFile)
    if (!file.exists(file.path(baseDir, "genes.tsv"))) {
      stop("No genes.tsv file with the matrix.mtx")
    }
    if (!file.exists(file.path(baseDir, "barcodes.tsv"))) {
      stop("No barcodes.tsv file with the matrix.mtx")
    }
    counts <- readMM(countsFile)
    geneNames <- read.delim(file.path(baseDir, "genes.tsv"), header = F,
                            sep = "\t", stringsAsFactors = F)
    rownames(counts) <- geneNames$V1
    cellNames <- read.delim(file.path(baseDir, "barcodes.tsv"), header = F,
                            sep = "\t", stringsAsFactors = F)
    colnames(counts) <- cellNames$V1
    return(list(counts, cellNames, geneNames))
  }
}

## los errores respecto a la construcción del objeto ya están controlados
## por la propia clase. Se pueden hacer, pero no sé hasta qué punto es necesario
## hacerlo. Tampoco sé hasta qué punto estaría bien hacer esta función visible
## para el resto de usuarios. Para nosotros sí.
CreateSCEObject <- function(counts, cellsMetadata, genesMetadata) {
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = cellsMetadata,
    rowData = genesMetadata
  )
  return(sce)
}




.loadSingleCellData <- function(single.cell, real = TRUE) {
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

        list.data <- list(.readCountsFile(single.cell[[1]]),
                          .readTabFiles(single.cell[[2]]),
                          .readTabFiles(single.cell[[3]]))
      } else {
        stop(paste("Incorrect data elements given. Please look at the allowed data for",
                   arg, "in ?CreateDigitalDLSorterObject"))
      }
    } else if (class(single.cell) == "character" & grepl(".mtx", single.cell)) {
      list.data <- .readCountsFile(single.cell)
    } else {
      stop(paste(arg, "argument is not recognizable. Please look at the allowed data for",
                 arg, "in ?CreateDigitalDLSorterObject"))
    }
  }
  single.cell <- CreateSCEObject(list.data[[1]], list.data[[2]], list.data[[3]])
  return(single.cell)
}


CreateDigitalDLSorterObject <- function(
  single.cell.real,
  gene.ID.col = NULL,
  project = "DigitalDLSorterProject",
  zinb.params.object = NULL,
  selected.genes = NULL,
  single.cell.sim = NULL
) {
  single.cell.real <- .loadSingleCellData(single.cell.real)
  single.cell.sim <- .loadSingleCellData(single.cell.sim, real = FALSE)

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



# fragment of .checkCommonItems

# common.genes <- intersect(rownames(counts), genes.metadata[, "external_gene_name"])
# diff <- abs(min(dim(counts)[1], dim(genes.metadata)[1]) - length(common.genes))
# if (length(common.genes) < min(dim(counts)[1], dim(genes.metadata)[1])) {
#   stop(paste("There are", diff, "genes that don't match."))
# } else if (diff != 0){
#   warning(paste("There are", diff, "genes that don't match."))# poner como warning o stop?
# } else {
#   warning(paste(abs(length(genes.metadata[, "external_gene_name"]) - length(common.genes)),
#                 "genes have been filtered."))
# }
# counts <- counts[common.genes, ]
# genes.metadata <- genes.metadata[genes.metadata[, "external_gene_name"] %in% common.genes, ]



.filterGenes <- function(counts, cells.metadata, cell.ID.column, genes.metadata,
                         gene.ID.column, min.counts, min.cells) {
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
  message(paste("Selected genes:",  dim(counts)[1]))
  message(paste("Discarded genes:", dim.bef[1] - dim(counts)[1]))

  genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in% rownames(counts), ]

  return(list(counts, cells.metadata, genes.metadata))
}


## duplicated genes and cells --------------------------------------------------

# remove duplicated genes and cells
# if (anyDuplicated(rownames(counts))) {
#   warning("Non-unique genes present in the counts matrix, making unique")
#   agg.genes <- factor(rownames(counts))
# }
