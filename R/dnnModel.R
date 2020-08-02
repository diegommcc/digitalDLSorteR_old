#' @importFrom dplyr %>%
#' @import keras
NULL

# tensorflow::tf$compat$v1$disable_eager_execution()
## this line should be in documentation or vignette, since it is a general option

#' Train DigitalDLSorter DNN model.
#'
#' Train DigitalDLSorter Deep Neural Network model with data store in
#' \code{final.data} slot. Moreover, model is evaluated on test data and
#' prediction results are produced.
#'
#' All the steps related with Deep Neural Network are performed by using
#' \code{keras} package, a API in R for \code{keras} in Python available from CRAN.
#' We recommend use the guide of installation available on \url{https://keras.rstudio.com/}
#' in order to set a custom configuration (type of backend used, CPU or
#' GPU, etc.). In the vignette of \code{digitalDLSorteR} we show the main steps for
#' the installation.
#'
#' Although \code{trainDigitalDLSorterModel} allows to select a custom loss
#' function used during training, we recommend using Kullback-Leibler divergence
#' because its better results. If you want to know more details about the
#' architecture of the DNN and its construction, see Torroja and Sanchez-Cabo,
#' 2019.
#'
#' @param object \code{DigitalDLSorter} object with \code{final.data} slot.
#' @param batch.size Number of samples per gradient update. If unspecified,
#' \code{batch.size} will default to 128.
#' @param numm.epochs Number of epochs to train the model.
#' @param val Boolean that determine if a validation subset is used during
#' training (\code{FALSE} by default).
#' @param freq.val Number between 0.1 and 0.5 that determine the number of
#' samples from training data that will be used as validation subset.
#' @param loss Character indicating loss function selected for training the model
#' (Kullback-Leibler divergence by default).
#' See \code{keras} documentation for more details.
#' @param metrics Vector fo metrics used to evaluate the performance of the model
#' during training and on test data (\code{c("accuracy", "mean_absolute_error",
#' "categorical_accuracy")} by default)
#' @param view.metrics.plots Boolean indicating if show progression plots of loss
#' and metrics during training (\code{TRUE} by default). \code{keras} for R allows
#' to see the progression of the model during training if you are working in RStudio.
#' @param verbose Boolean indicating if show the progression of the model during
#' training. Beisdes, it is shown information about the architecture of the model
#' (\code{TRUE} by default).
#' @return A DigitalDLSorter object with \code{trained.model} slot containing
#' a \code{\link{DigitalDLSorterDNN}} object. For more information about the
#' structure of this class, see \code{\link{?DigitalDLSorterDNN}}.
#'
#' @export
#'
#' @seealso \code{\link{plotTrainingHistory}}, \code{\link{deconvDigitalDLSorterModel}}.
#'
#' @examples
#' DDLSChung <- trainDigitalDLSorterModel(
#'   object = DDLSChung,
#'   batch.size = 128,
#'   num.epochs = 20
#' )
#'
#' @references
#' Torroja, C. y Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep Learning algorithm to quantify
#' immune cell populations based on scRNA-Seq data. Frontiers in Genetics 10, 978. doi:
#' \url{10.3389/fgene.2019.00978}
#'
trainDigitalDLSorterModel <- function(
  object,
  batch.size = 128,
  num.epochs = 20,
  val = FALSE,
  freq.val = 0.1,
  loss = "kullback_leibler_divergence",
  metrics = c("accuracy", "mean_absolute_error",
              "categorical_accuracy"),
  view.metrics.plot = TRUE,
  verbose = TRUE
) {
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (is.null(final.data(object))) {
    stop("'final.data' slot is empty")
  } else if (is.null(prob.cell.types(object))) {
    stop("prob.cell.types slot is empty")
  } else if (num.epochs <= 1) {
    stop("'num.epochs' argument must be greater than or equal to 2")
  } else if (batch.size <= 10) {
    stop("'batch.size' argument must be greater than or equal to 10")
  }
  if (!is.null(trained.model(object))) {
    warning("'trained.model' slot is not empty. For the moment, digitalDLSorteR",
            " does not support for multiple trained models, so the actual model",
            " will be overwritten\n",
            call. = FALSE, immediate. = TRUE)
  }
  if (view.metrics.plot) {
    view.plot <- "auto"
  } else {
    view.plot <- 0
  }
  if (verbose) {
    verbose.model <- 1
  } else {
    verbose.model <- 0
  }
  # number of samples
  n.train <- nrow(object@final.data[["train"]])
  n.test <- nrow(object@final.data[["test"]])

  # set the model. may be we can provide a interface to set a custom model
  model <- keras_model_sequential(name = "DigitalDLSorter") %>%
    layer_dense(units = 200, input_shape = c(ncol(object@final.data$train)),
                name = "Dense1") %>%
    layer_batch_normalization(name = "BatchNormalization1") %>%
    layer_activation(activation = "relu", name = "ActivationReLu1") %>%
    layer_dropout(rate = 0.25, name = "Dropout1") %>%
    layer_dense(units = 200, name = "Dense2") %>%
    layer_batch_normalization(name = "BatchNormalization2") %>%
    layer_activation(activation = "relu", name = "ActivationReLu2") %>%
    layer_dropout(rate = 0.25, name = "Dropout2") %>%
    layer_dense(units = ncol(object@final.data$train@metadata$prob.matrix),
                name = "Dense3") %>%
    layer_batch_normalization(name = "BatchNormalization3") %>%
    layer_activation(activation = "softmax", name = "ActivationSoftmax")

  if (verbose) {
    summary(model)
  }
  # allow set optimizer?
  model %>% compile(
    loss = loss,
    optimizer = optimizer_adam(),
    metrics = metrics
  )
  ## training model
  if (val) {
    if (freq.val < 0.1 || freq.val > 0.5) {
      stop("freq.val must be a float between 0.1 and 0.5")
    }
    # with validation subset. generator divide train set in two subsets
    n.val <- ceiling(n.train * freq.val)
    if (verbose) {
      message(paste("\n=== Training DNN with", n.train - n.val, "samples and evaluating with",
                    n.val, "samples:\n"))
    }
    train.generator <- .trainGenerator(
      object = object,
      type.data = "train",
      batch.size = batch.size,
      shuffle = TRUE,
      min.index = n.val,
      max.index = n.train
    )
    val.generator <- .trainGenerator(
      object = object,
      type.data = type.data,
      batch.size = batch.size,
      shuffle = FALSE,
      min.index = 0,
      max.index = n.val
    )
    history <- model %>% fit_generator(
      generator = train.generator,
      steps_per_epoch = ceiling((n.train - n.val) / batch.size),
      epochs = num.epochs,
      validation_data = val.generator,
      validation_steps = n.val / batch.size,
      verbose = verbose.model,
      view_metrics = view.plot
    )
  } else {
    # without validation subset
    if (verbose) {
      message(paste("\n=== Training DNN with", n.train, "samples:\n"))
    }
    train.generator <- .trainGenerator(
      object = object,
      type.data = "train",
      batch.size = batch.size,
      shuffle = TRUE
    )
    history <- model %>% fit_generator(
      generator = train.generator,
      steps_per_epoch = ceiling(n.train / batch.size),
      epochs = num.epochs,
      verbose = verbose.model,
      view_metrics = view.plot
    )
  }

  if (verbose) {
    message(paste0("\n=== Evaluating DNN in test data (", n.test, " samples)"))
  }
  ## evaluation of the model: set by default, no options?
  test.generator <- .trainGenerator(object = object,
                                    type.data = "test",
                                    batch.size = batch.size,
                                    shuffle = FALSE)
  test.eval <- model %>% evaluate_generator(
    generator = test.generator,
    steps = ceiling(n.test / batch.size)
  )
  ## prediction of test samples
  if (verbose) {
    message(paste0("  - ", names(test.eval), ": ", lapply(test.eval, round, 4),
                   collapse = "\n"))
    message(paste("\n=== Generating prediction results using test data\n"))
  }
  predict.generator <- .predictTestGenerator(
    object = object,
    type.data = "test",
    batch.size = batch.size
  )
  predict.results <- model %>% predict_generator(
    generator = predict.generator,
    steps = ceiling(n.test / batch.size),
    verbose = verbose.model
  )
  rownames(predict.results) <- rowData(object@final.data[["test"]])[[1]]
  colnames(predict.results) <- colnames(object@final.data[["test"]]@metadata[[1]])

  network.object <- new(
    Class = "DigitalDLSorterDNN",
    model = model,
    training.history = history,
    eval.stats.model = test.eval,
    predict.results = predict.results,
    cell.types = colnames(object@final.data[["train"]]@metadata[[1]]),
    features = colData(object@final.data[["test"]])[[1]]
  )
  object@trained.model <- network.object

  message("DONE")
  return(object)
}



.trainGenerator <- function(
  object,
  type.data,
  batch.size,
  shuffle,
  min.index = NULL,
  max.index = NULL
) {
  if (!is.null(min.index) && !is.null(max.index)) {
    n.samples <- length(min.index:max.index)
    nb <- min.index
  } else {
    n.samples <- nrow(assay(object@final.data[[type.data]]))
    min.index <- 1
    nb <- 0
  }
  n.features <- ncol(object@final.data[[type.data]])
  n.classes <- ncol(object@final.data[[type.data]]@metadata$prob.matrix)
  function() {
    data.index <- seq(nb + 1, nb + batch.size)
    nb <<- nb + batch.size
    if (nb > n.samples) {
      data.index <- data.index[data.index <= n.samples]
      fill <- batch.size - length(data.index)
      data.index <- c(data.index, seq(min.index + 1, min.index + fill))
      if (fill <= min.index) nb <<- min.index + 1
      else nb <<- fill
    }
    if (shuffle) {
      shuffling <- sample(1:length(data.index))
      return(list(matrix(assay(object@final.data[[type.data]])[data.index, ],
                         ncol = n.features, nrow = length(data.index))[shuffling, ],
                  matrix(object@final.data[[type.data]]@metadata$prob.matrix[data.index, ],
                         ncol = n.classes, nrow = length(data.index))[shuffling, ]))
    } else {
      return(list(matrix(assay(object@final.data[[type.data]])[data.index, ],
                         ncol = n.features, nrow = length(data.index)),
                  matrix(object@final.data[[type.data]]@metadata$prob.matrix[data.index, ],
                         ncol = n.classes, nrow = length(data.index))))
    }
  }
}

.predictTestGenerator <- function(
  object,
  type.data,
  batch.size
) {
  nb <- 0
  n.samples <- nrow(assay(object@final.data[[type.data]]))
  n.features <- ncol(object@final.data[[type.data]])
  n.classes <- ncol(object@final.data[[type.data]]@metadata$prob.matrix)
  function() {
    data.index <- seq(nb + 1, nb + batch.size)
    nb <<- nb + batch.size
    if (nb > n.samples) {
      data.index <- data.index[data.index <= n.samples]
      nb <<- 0
    }
    return(list(matrix(assay(object@final.data[[type.data]])[data.index, ],
                       ncol = n.features, nrow = length(data.index))))
  }
}


#' Save on disk trained DigitalDLSorter DNN model as HDF5 file.
#'
#' Save on disk the trained model in HDF5 format.
#'
#' @param object \code{DigitalDLSorter} object with \code{trained.model} slot.
#' @param file.path Valid file path where storing model.
#' @param overwrite Overwrite file if it already exists.
#'
#' @export
#'
#' @seealso \code{\link{trainDigitalDLSorterModel}} \code{\link{loadTrainedModel}}.
#'
#' @examples
#'
saveTrainedModelAsH5 <- function(
  object,
  file.path,
  overwrite = FALSE
) {
  if (class(object) != "DigitalDLSorter") {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (is.null(trained.model(object))) {
    stop("'trained.model' slot is empty")
  } else if (is.null(trained.model(object)@model)) {
    stop("There is not model to save on disk. First, train model with trainDigitalDLSorterModel function")
  }
  if (file.exists(file.path)) {
    if (overwrite) {
      message(paste(file.path, "file exists. Since 'overwrite' argument is TRUE, it will be overwritten"))
    } else {
      stop(paste(file.path, "file exists"))
    }
  }
  tryCatch({
    save_model_hdf5(object = object@trained.model@model,
                    filepath = file.path,
                    overwrite = overwrite,
                    include_optimizer = TRUE)
  }, error = function(cond) {
    message(paste("\nProblem during saving", file.path))
    stop(cond)
  })
}


#' Load from HDF5 file a trained DigitalDLSorter DNN model.
#'
#' Load from HDF5 file a trained DigitalDLSorter DNN model and incorporate
#' it in DigitalDLSorter object.
#'
#' @param object \code{DigitalDLSorter} object with \code{trained.model} slot.
#' @param file.path Valid file path where storing model.
#' @param reset.slot Remove \code{trained.slot} if it already exists.
#' A new DigitalDLSorterDNN object will be formed, but it will not contain other
#' slots.
#'
#' @export
#'
#' @seealso \code{\link{trainDigitalDLSorterModel}}, \code{\link{deconvDigitalDLSorterModel}},
#' \code{\link{loadTrainedModel}}
#'
#' @examples
#'
loadTrainedModelFromH5 <- function(
  object,
  file.path,
  reset.slot = TRUE
) {
  if (class(object) != "DigitalDLSorter") {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (!file.exists(file.path)) {
    stop(paste(file.path, "file does not exist. Please provide a valid file path"))
  }
  if (!is.null(trained.model(object))) {
    slot.exists <- TRUE
    message("'trained.model' slot is not empty:")
    if (reset.slot) {
      message("  'reset.slot' is TRUE, 'trained.model' slot will be restart")
    } else {
      message("  'reset.slot' is FALSE, only 'model' slot of DigitalDLSorterDNN object will be overwritten")
    }
  } else {
    slot.exists <- FALSE
  }
  tryCatch({
    loaded.model <- load_model_hdf5(filepath = file.path, compile = TRUE)
  }, error = function(cond) {
    message(paste("\n", file.path, "file provided is not a valid Keras model:"))
    stop(cond)
  })

  if (!slot.exists) {
    model <- new(Class = "DigitalDLSorterDNN",
                 model = loaded.model)
  } else {
    if (reset.slot) {
      model <- new(Class = "DigitalDLSorterDNN",
                   model = loaded.model)
    } else {
      model(object@trained.model) <- loaded.model
      return(object)
    }
  }
  trained.model(object) <- model
  return(object)
}

#' Plot training history of a trained DigitalDLSorter DNN model.
#'
#' Plot training history of a trained DigitalDLSorter DNN model.
#'
#' @param object \code{DigitalDLSorter} object with \code{trained.model} slot.
#' @param title Title of plot.
#' @param metrics Which metrics to plot. If it is equal to \code{NULL} (by default),
#' all metrics will be plotted.
#'
#' @export
#'
#' @seealso \code{\link{trainDigitalDLSorterModel}}, \code{\link{deconvDigitalDLSorterModel}}
#'
#' @examples
#'
plotTrainingHistory <- function(
  object,
  title = "History of metrics during training",
  metrics = NULL
) {
  if (class(object) != "DigitalDLSorter") {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (is.null(trained.model(object))) {
    stop("'trained.model' slot is empty")
  } else if (is.null(trained.model(object)@training.history)) {
    stop("There is not training history on the provided object")
  }
  if (!is.null(metrics)) {
    if (!all(metrics %in% names(object@trained.model@training.history$metrics))) {
      stop("Some of the given metrics are not present in the provided object")
    }
  }
  plot(object@trained.model@training.history,
       metrics = metrics, method = "ggplot2") + ggtitle(title)
}


#' Load data to deconvolute from text file.
#'
#' Load data to deconvolute from text file. Accepted formats are tsv and tas.gz.
#'
#' @param object \code{DigitalDLSorter} object with \code{trained.model} slot.
#' @param file.path File path where data is stored.
#' @param name.data Name with which the data is stored in \code{DigitalDLSorter}
#' object. If \code{name.data} is not provided, basename of file is used.
#'
#' @export
#'
#' @seealso \code{\link{trainDigitalDLSorterModel}}, \code{\link{deconvDigitalDLSorterModel}}
#'
#' @examples
#'
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
## hacer genérica y que funione de forma diferente en función de si es
## SummarizedExperiment o matrix

#' Load data to deconvolute from \code{SummarizedExperiment} object.
#'
#' Load data in \code{DigitalDLSorter} object to deconvolute from
#' \code{SummarizedExperiment} object.
#'
#' @param object \code{DigitalDLSorter} object with \code{trained.model} slot.
#' @param se.object \code{SummarizedExperiment} object.
#' @param name.data Name with which the data is stored in \code{DigitalDLSorter}
#' object.
#'
#' @export
#'
#' @seealso \code{\link{trainDigitalDLSorterModel}}, \code{\link{deconvDigitalDLSorterModel}}
#'
#' @examples
#'
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
    if (is.null(decov.data(object))) {
      name.data <- "deconv_1"
    } else {
      name.data <- paste0("decon_", length(decov.data(object)) + 1)
    }
  }
  ## create or not a new list
  if (is.null(deconv.data(object))) list.data <- list()
  else list.data <- deconv.data(object)
  ## check if name.data exists
  if (name.data %in% names(list.data)) {
    stop(paste(name.data, "data already exists in deconv.data slot"))
  }
  list.data[[name.data]] <- se.object
  deconv.data(object) <- list.data

  return(object)
}

## función para establecer el fichero desde el que llevar a cabo la predicción
# giveFilepathDataForDeconv <- function(
#   object,
#   file.path,
#   name.data,
#   chunk.size
# ) {
# }


#' Deconvolute bulk gene expression samples (RNA-Seq) using a pre-trained
#' DigitalDLSorter model.
#'
#' Deconvolute bulk gene expression samples (RNA-Seq) enumerating and
#' quantifying the proportion of cell types present in a bulk sample. See in
#' Details the available models.
#'
#' This method uses a pre-trained Deep Neural Network model to enumerate
#' and quantify the cell types present in bulk RNA-Seq samples. For the moment,
#' the available models allow to deconvolute the immune infiltration in colorectal
#' (\code{'colon.li'} (Li et al., 2017)) and breast cancer (\code{'breast.chung'}
#' (Chung et al., 2017)).
#'
#' This function is oriented for users that only want to use the method for
#' deconvoluting their bulk RNA-Seq samples. For users that are buliding their
#' own model using the rest of functionalities of \code{digitalDLSoteR} package,
#' see \code{\link{deconvDigitalDLSorterObj}}. The former works with base classes,
#' while the last uses \code{DigitalDLSorter} objects.
#'
#' @param data A matrix bulk gene expression with
#' @param model Pre-trained DNN model to use for deconvoluting process.
#' For the moment, the available
#' models are for deconvoluting RNA-Seq samples from breast cancer ('breast') and
#' colon cancer ('colon') environments.
#' @param batch.size Number of samples per gradient update. If unspecified,
#' \code{batch.size} will default to 128.
#' @param transpose Transpose data for set features (genes) as columns and
#' samples as rows. \code{TRUE} by default.
#' @param normalize Normalize data before deconvolution. \code{TRUE} by default.
#' @param verbose Show messages during the execution.
#' @return A \code{data.frame} with samples (\eqn{i}) as rows and cell types
#' (\eqn{j}) as columns. Each entry represents the proportion of \eqn{j} cell
#' type in \eqn{i} sample.
#'
#' @export
#'
#' @seealso \code{\link{deconvDigitalDLSorterObj}}.
#'
#' @examples
#' DDLSChung <- trainDigitalDLSorterModel(
#'   object = DDLSChung,
#'   batch.size = 128,
#'   num.epochs = 20
#' )
#'
#' @references
#' Chung, W., Eum, H. H., Lee, H. O., Lee, K. M., Lee, H. B., Kim,
#' K. T., et al. (2017). Single-cell RNA-seq enables comprehensive tumour and
#' immune cell profiling in primary breast cancer. Nat. Commun. 8 (1), 15081.
#' doi: \url{10.1038/ncomms15081}.
#'
deconvDigitalDLSorter <- function(
  data,
  model = "breast",
  batch.size = 128,
  transpose = FALSE,
  normalize = FALSE,
  verbose = TRUE
) {
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a matrix or data.frame")
  }
  if (model == "breast") {
    model.dnn <- digitalDLSorteR::breast.chung
    model.dnn <- .loadModelFromJSON(model.dnn)
  } else if (model == "colon") {
    model.dnn <- digitalDLSorteR::colon.li
    model.dnn <- .loadModelFromJSON(model.dnn)
  }

  ## check data --> habría que hacer algún checkeo sobre la matriz, genes,
  ## check si hay duplicados: hacer función para ello y que funcione también
  ## cuando se cargan los datos en loadRealSCProfiles
  results <- .deconvCore(
    deconv.counts = data,
    model = model.dnn,
    batch.size = batch.size,
    transpose = transpose,
    normalize = normalize,
    verbose = verbose
  )
  return(results)
}


#' Deconvolute bulk gene expression samples (RNA-Seq).
#'
#' Deconvolute bulk gene expression samples (RNA-Seq) enumerating and
#' quantifying the proportion of cell types present in a bulk sample. This
#' function needs a \code{DigitalDLSorter} object with a trained DNN model
#' (\code{\link{trained.model}} slot) and data for deconvoluting store in
#' \code{deconv.data} slot.
#'
#' This function is oriented for users that have trained a DNN model using their
#' own data. If you want to use a pre-trained model, see
#' \code{\link{deconvDigitalDLSorter}}.
#'
#' @param object \code{\link{DigitalDLSorter}} object with \code{trained.data}
#' and \code{deconv.data} slots.
#' @param name.data Name of the data store in \code{DigitalDLSorter} object. If
#' it is not provided, the first dataset will be used.
#' @param batch.size Number of samples per gradient update. If unspecified,
#' \code{batch.size} will default to 128.
#' @param transpose Transpose data for set features (genes) as columns and
#' samples as rows. \code{TRUE} by default.
#' @param normalize Normalize data before deconvolution. \code{TRUE} by default.
#' @param verbose Show messages during the execution.
#' @return A \code{data.frame} with samples (\eqn{i}) as rows and cell types
#' (\eqn{j}) as columns. Each entry represents the proportion of \eqn{j} cell
#' type in \eqn{i} sample.
#'
#' @export
#'
#' @seealso \code{\link{trainDigitalDLSorter}} \code{\link{DigitalDLSorter}}
#'
#' @examples
#' DDLSChung <- trainDigitalDLSorterModel(
#'   object = DDLSChung,
#'   batch.size = 128,
#'   num.epochs = 20
#' )
#'
#' @references
#' Torroja, C. y Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep Learning algorithm to quantify
#' immune cell populations based on scRNA-Seq data. Frontiers in Genetics 10, 978. doi:
#' \url{10.3389/fgene.2019.00978}
#'
deconvDigitalDLSorterObj <- function(
  object,
  name.data,
  batch.size = 128,
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

  ## check if model is json format or compiled
  if (is.list(trained.model(object)@model)) {
    model.comp <- .loadModelFromJSON(trained.model(object))
    trained.model(object) <- model.comp
  }

  ## comming soon
  # if (!is.null(file.path)) {
  #   if (!file.exists(file.path)) {
  #     stop("The provided file does not exists")
  #   }
  #   message(paste("=== Deconvolution of data stored on disk in", file.path, "\n"))
  #   ## set generator that works with CSV files on disk
  # } else {
  # message("=== Deconvolution of data provided in deconv.data slot of DigitalDLSorter object")
  ## access to data from object. Todo esto lo implemento en una función interna
  if (missing(name.data)) {
    message("   No name.data provided. Using the first dataset\n")
    name.data <- 1
  }
  deconv.counts <- assay(object@deconv.data[[name.data]])

  ## deconvolution
  results <- .deconvCore(
    deconv.counts = deconv.counts,
    model = trained.model(object),
    batch.size = batch.size,
    transpose = transpose,
    normalize = normalize,
    verbose = verbose
  )
  object@deconv.results[[name.data]] <- results

  return(object)
}


.deconvCore <- function(
  deconv.counts,
  model,
  batch.size,
  transpose,
  normalize,
  verbose
) {
  if (is.null(rownames(deconv.counts))) {
    stop("The given matrix does not have column names. You must provide a matrix",
         " with feature names in the same notation used in training data")
  }
  # this can do it more elegant and probably more efficient
  ## filtering features missing in training data
  filter.features <- rownames(deconv.counts) %in% features(model)
  deconv.counts <- deconv.counts[filter.features, ]
  ## set features missing in deconv.data
  fill.features <- !features(model) %in% rownames(deconv.counts)
  m.new <- matrix(0L, nrow = sum(fill.features), ncol = ncol(deconv.counts))
  rownames(m.new) <- features(model)[fill.features]
  deconv.counts <- rbind(deconv.counts, m.new)
  deconv.counts <- deconv.counts[features(model), ]
  if (verbose) {
    message(paste("=== Filtering", sum(filter.features),
                  "features in data that are not present in training data\n"))
    message(paste("=== Setting", sum(fill.features),
                  "features that are not present in training data to zero\n"))
  }
  if (normalize) {
    if (verbose) {
      message("=== Normalizing data\n")
    }
    deconv.counts <- edgeR::cpm.default(deconv.counts)
    # deconv.counts <- rescale(deconv.counts)
    deconv.counts <- scale(deconv.counts)
  }
  deconv.counts <- t(deconv.counts)
  deconv.generator <- .predictDeconvDataGenerator(
    data = deconv.counts,
    model = model,
    batch.size = batch.size
  )
  if (verbose) {
    verbose.model <- 1
    message("=== Predicting cell types present in the provided samples\n")
  } else {
    verbose.model <- 0
  }
  dnn.model <- model(model)
  results <- dnn.model %>% predict_generator(
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
  colnames(results) <- cell.types(model)

  return(results)
}

.predictDeconvDataGenerator <- function(
  data,
  model,
  batch.size
) {
  nb <- 0
  n.samples <- nrow(data)
  n.features <- length(features(model))
  n.classes <- length(cell.types(model))
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


.loadModelFromJSON <- function(object) {
  model.list <- model(object)
  model.comp <- model_from_json(model.list[[1]])
  model.comp <- set_weights(model.comp, model.list[[2]])
  model(object) <- model.comp

  return(object)
}

.saveModelToJSON <- function(object) {
  model.comp <- model(object)
  model.json <- model_to_json(model.comp)
  weights <- get_weights(model.comp)
  model(object) <- list(model.json, weights)

  return(object)
}

# .predictDeconvFileGenerator <- function(
#   object,
#   name.data,
#   batch.size,
#   transpose,
#   normalize
# ) {
#
# }

