## digitalDLSorter Keras


trainDigitalDLSorterNetwork <- function(object,
                                        batch.size,
                                        num.epochs,
                                        loss = "kullback_leibler_divergence") {



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

  model %>% compile(
    loss = loss,
    optimizer = optimizer_adam(),
    metrics = c("accuracy", "metric_mean_absolute_error",
                "metric_mean_absolute_percentage_error",
                "metric_kullback_leibler_divergence",
                "metric_categorical_accuracy")
  )



}










