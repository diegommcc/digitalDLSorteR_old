## digitalDLSorter Keras
library(keras)
tensorflow::tf$compat$v1$disable_eager_execution()

trainDigitalDLSorterNetwork <- function(object,
                                        batch.size,
                                        num.epochs,
                                        val,
                                        val.freq = 0.1,
                                        loss = "kullback_leibler_divergence",
                                        plot = TRUE) {



  if (!plot) {
    opt <- options(keras.view_metrics = FALSE)
  }
  model <- keras_model_sequential(name = "DigitalDLSorter") %>%
    layer_dense(units = 200, input_shape = c(dim(object@final.data$train)[2]),
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

  ## dar opciones para elegir?
  model %>% compile(
    loss = loss,
    optimizer = optimizer_adam(),
    metrics = c("accuracy", "mean_absolute_error",
                "mean_absolute_percentage_error",
                "kullback_leibler_divergence",
                "categorical_accuracy")
  )

  if (val) {
    n.val <- ceiling(nrow(object@final.data[["train"]]) * freq.val)
    history <- model %>% fit_generator(generator = .dataGenerator(object = object,
                                                                  type.data = "train",
                                                                  batch.size = batch.size,
                                                                  shuffle = TRUE,
                                                                  min.index = n.val,
                                                                  max.index = nrow(object@final.data[["train"]])),
                            steps_per_epoch = ceiling((nrow(object@final.data[["train"]]) - n.val) / batch.size),
                            epochs = num.epochs,
                            validation_data = .dataGenerator(object = object,
                                                             type.data = type.data,
                                                             batch.size = batch.size,
                                                             shuffle = FALSE,
                                                             min.index = 0,
                                                             max.index = n.val),
                            validation_steps = n.val / bath.size)
  } else {
    # sin validación
    history <- model %>% fit_generator(generator = .dataGenerator(object = object,
                                                                  type.data = type.data,
                                                                  batch.size = batch.size,
                                                                  shuffle = TRUE),
                            steps_per_epoch = ceiling(nrow(object@final.data[["train"]]) / batch.size),
                            epochs = num.epochs)
  }



  history <- model %>% fit_generator(.dataGeneratorShufflePartial(object,
                                                       type.data = type.data,
                                                       batch.size = batch.size),
                          steps_per_epoch = ceiling(nrow(object@final.data[[type.data]]) / batch.size),
                          epochs = num.epochs)

}



.dataGenerator <- function(object,
                           type.data,
                           batch.size,
                           shuffle,
                           min.index = NULL,
                           max.index = NULL) {
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
      # nb <<- 0
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

  opt
}



