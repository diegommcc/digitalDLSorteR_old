## predictDigitalDLSorterModel
# Respecto a este paso, creo que hay varias posibilidades.
# * Cargar los datos de predicción en el objeto, lo que creo que no es eficiente
# computacionalmente y que no tiene mucho sentido.
# * Cargarlos como matriz en R, pero seguimos con el mismo problema.
# * Una función que cargue los datos en un fichero H5 con DelayedArray, de
# forma que los datos se guardarán en un objeto de esta clase con backend en un
# fichero H5. Esto realmente es jodido, ya que si son muchos datos, está el problema
# de cargar todo eso en memoria.
# * Tener ya el fichero H5 hecho y linkarlo a un objeto DelayedArray.
# * Pemritir hacer lecturas de ficheros TSV igual que ya está hecho en Python:
#   generar un generador que funcione con ficheros TSV.

# Lo suyo será hacer una función para cargar todo en memoria, otra función para
# generar el fichero H5 con backend y una última función con un generador que funcione
# con ficheros de texto.
# Respecto a lo de cargar todo en memoria, también se puede permitir utiliazr
# el predictor sobre una matriz ya cargada


# si el fichero es demasiado grande, permitir leer por chunks: no creo que sea
# lo mejor, sólo en el caso de que lo quiera cargar en el objeto como HDF5.
# loadDeconvDataHDF5Backend <- function(
#   object,
#   file.path,
#   file.backend,
#   large = FALSE,
#   transposed = TRUE
# ) {
#
# }


## check si el fichero tiene colnames y rownames
loadDeconvDataFromFile <- function(
  object,
  file.path,
  name.data = NULL
) {
  if (class(object) != "DigitalDLSorter") {
    stop("The provided object is not of DigitalDLSorter class")
  }
  se.object <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = .readTabFiles(file = file.path))
  )
  ## generate name for data if is not provided
  if (is.null(name.data)) {
    name.data <- tools::file_path_sans_ext(basename(file.path))
  }
  ## create or not a new list
  if (is.null(object@deconv.data)) list.data <- list()
  else list.data <- object@deconv.data
  ## check if name.data exists
  if (name.data %in% names(list.data)) {
    stop(paste(name.data, "data already exists in deconv.data slot"))
  }
  list.data[[name.data]] <- se.object
  object@deconv.data <- list.data

  return(object)
}


## check si la matriz tiene colnames y rownames
loadDeconvDataFromSummarizedExperiment <- function(
  object,
  se.object,
  name.data = NULL
) {
  if (class(object) != "DigitalDLSorter") {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (class(se.object) != "SummarizedExperiment") {
    stop("The provided object is not of SummarizedExperiment class")
  }
  if (length(assays(se.object)) == 0) {
    stop("assay slot of SummarizedExperiment object is empty")
  } else if (length(assays(se.object)) > 1) {
    warning(paste("There are more than one assays in SummarizedExperiment object,",
                  "only the first assay will be considered. Remember that it is", "
                  recommended that the provided data be of the same nature as",
                  "the data with which the model has been trained (e.g. TPMs)"))
  }
  ## generate name for data if is not provided
  if (is.null(name.data)) {
    name.data <- tools::file_path_sans_ext(basename(file.path))
  }
  ## create or not a new list
  if (is.null(object@deconv.data)) list.data <- list()
  else list.data <- object@deconv.data
  ## check if name.data exists
  if (name.data %in% names(list.data)) {
    stop(paste(name.data, "data already exists in deconv.data slot"))
  }

  list.data[[name.data]] <- se.object

  object@deconv.data <- list.data
  return(object)
}

## función para establecer el fichero desde el que llevar a cabo la predicción
giveFilepathDataForDeconv <- function(
  object,
  file.path,
  name.data,
  chunk.size
) {

}


# pretrained.model = NULLsi si no es nulo
# Si no es NULL, permitir que coja algún modelo preentrenado por nosotros:
# Chung, el de cáncer, etc.

# argumento file.path en el caso de que se quiera entrenar desde un csv en disco.
# de momento voy a obligar al usuario a cargar los datos en memoria
# recordar que los genes deben estar en la misma notación que los datos con los
# que ha sido entrenada la red
deconvDigitalDLSorterModel <- function(
  object,
  name.data = NULL,
  file.path = NULL,
  batch.size = 128,
  # pretrained.model = NULL,
  transpose = FALSE,
  normalize = FALSE,
  verbose = TRUE
) {
  if (class(object) != "DigitalDLSorter") {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (is.null(object@trained.model)) {
    stop("There is not trained model in DigitalDLSorter object")
  } else if (batch.size <= 10) {
    stop("'batch.size' argument must be greater than or equal to 10")
  }
  ## comming soon
  if (!is.null(file.path)) {
    if (!file.exists(file.path)) {
      stop("The provided file does not exists")
    }
    message(paste("=== Deconvolution of data stored on disk in", file.path, "\n"))
    ## set generator that works with CSV files on disk

  } else {
    message("=== Deconvolution of data provided in deconv.data slot of DigitalDLSorter object")
    ## access to data from object. Todo esto lo implemento en una función interna
    if (is.null(name.data)) {
      message("   No name.data provided. Using the first dataset\n")
      name.data <- 1
    }
    deconv.counts <- object@deconv.data[[name.data]]

    if (transpose) {
      deconv.counts <- t(deconv.counts)
    }
    if (is.null(colnames(deconv.counts))) {
      stop("The given matrix does not have column names. You must provide a matrix with feature names in the same notation used in training data")
    }
    if (verbose) {
      message(paste("=== Filtering", sum(filter.features),
                    "features in data that are not present in training data\n"))
      message(paste("=== Setting", sum(fill.features),
                    "features that are not present in training data to zero\n"))
    }
    # this can do it more elegant and probably more efficient
    ## filtering features missing in training data
    filter.features <- colnames(deconv.counts) %in% colData(object@final.data$train)[[1]]
    deconv.counts <- deconv.counts[, filter.features]
    ## set features missing in deconv.data
    fill.features <- !colData(object@final.data$train)[[1]] %in% colnames(deconv.counts)
    m.new <- matrix(0L, ncol = sum(fill.features), nrow = nrow(deconv.counts))
    colnames(m.new) <- colData(object@final.data$train)[[1]][fill.features]
    deconv.counts <- cbind(deconv.counts, m.new)
    deconv.counts <- deconv.counts[, colData(object@final.data$train)[[1]]]

    if (normalize) {
      if (verbose) {
        message("=== Normalizing data\n")
      }
      deconv.counts <- edgeR::cpm.default(deconv.counts)
      deconv.counts <- scale(deconv.counts)
    }
    deconv.generator <- .predictDeconvDataGenerator(
      data = deconv.counts,
      batch.size = 128,
      type.data = "train"
    )
  }
  if (verbose) {
    verbose.model <- 1
    message("=== Predicting cell types present in the provided samples\n")
  } else {
    verbose.model <- 0
  }


  model <- object@trained.model@model
  results <- model %>% predict_generator(
    generator = deconv.generator,
    steps = ceiling(nrow(deconv.counts) / batch.size),
    verbose = verbose.model
  )
  if (!is.null(rownames(deconv.counts))) {
    rownames.deconv <- rownames(deconv.counts)
  } else {
    rownames.deconv <- seq(1, nrow(deconv.counts))
  }
  rownames(results) <- rownames.deconv
  colnames(results) <- object@trained.model@classes

  object@deconv.results[[name.data]] <- results

  return(object)
}

## normalización y transposición mejor con los datos completos en lugar de en el
# generador. Esto solo en el caso de que sea un fichero. Incluso en el HDF5 es mejor
# no hacerlo así, ya que él lo hace por chunks
.predictDeconvDataGenerator <- function(
  data,
  batch.size,
  type.data = "train"
) {
  nb <- 0
  n.samples <- nrow(data)
  n.features <- ncol(object@final.data[[type.data]])
  n.classes <- ncol(object@final.data[[type.data]]@metadata$prob.matrix)
  function() {
    data.index <- seq(nb + 1, nb + batch.size)
    nb <<- nb + batch.size
    if (nb > n.samples) {
      data.index <- data.index[data.index <= n.samples]
      nb <<- 0
    }
    return(list(matrix(data[data.index, ],
                       ncol = n.features,
                       nrow = length(data.index))))
  }
}

.predictDeconvFileGenerator <- function(
  object,
  name.data,
  batch.size,
  transpose,
  normalize
) {

}



