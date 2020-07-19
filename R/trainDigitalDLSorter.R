## digitalDLSorter Keras
# probar a meter en el generador un cat informativo
library(keras)
tensorflow::tf$compat$v1$disable_eager_execution() # no sé qué hacer con esto

trainDigitalDLSorterModel <- function(
  object,
  batch.size = 128,
  num.epochs = 20,
  val = FALSE,
  val.freq = 0.1,
  loss = "kullback_leibler_divergence",
  metrics = c("accuracy", "mean_absolute_error",
              "kullback_leibler_divergence",
              "categorical_accuracy"),
  view.metrics.plot = TRUE,
  verbose = TRUE
) {
  if (class(object) != "DigitalDLSorter") {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (is.null(final.data(object))) {
    stop("'final.data' slot is empty")
  } else if (is.null(prob.matrix(object))) {
    stop("prob.matrix slot is empty")
  } else if (num.epochs <= 1) {
    stop("'num.epochs' argument must be greater than or equal to 2")
  } else if (batch.size <= 10) {
    stop("'batch.size' argument must be greater than or equal to 10")
  }
  # errors with loss and metrics are controled by keras
  # if (!is.null(trained.model(object))) {
  #   warning("'trained.model' slot is not empty. For the moment, digitalDLSorteR does not support for multiple trained models, so the actual model will be overwritten")
  # }

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
    message(paste0("\n=== Evaluating DNN with test data (", n.test, " samples)"))
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
    message(paste("\n=== Generating prediction results on test data\n"))
  }
  predict.generator <- .predictGenerator(
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
    trained.model = model,
    training.history = history,
    eval.stats = test.eval,
    predict.results = predict.results
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

.predictGenerator <- function(
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


