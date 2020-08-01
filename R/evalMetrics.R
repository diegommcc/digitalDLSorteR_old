#' @importFrom dplyr mutate as_tibble left_join
#' @importFrom tidyr gather

## evaluation metrics
##

# Lista de functiones:
# tidyr: gather
# dplyr: mutate, as_tibble, left_join

color.list <- c(brewer.pal(12, "Paired"),
                brewer.pal(12, "Set3"),
                brewer.pal(8, "Pastel2"),
                colorRampPalette(c("grey20","grey70"))(4))


calculateEvalMetrics <- function(
  object,
  metrics = c("MAE", "MSE")
) {
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (is.null(trained.model(object)) ||
             is.null(trained.model(object)@predict.results)) {
    stop("The provided object does not have a trained model for evaluation")
  } else if (is.null(final.data(object)) ||
             !"test" %in% names(final.data(object))) {
    stop("The provided object does not have the test prediction target for evaluation")
  } else if (is.null(final.data(object, "test")@metadata)) {
    stop("The provided object does not have the test prediction target for evaluation")
  }
  ## validation metrics
  valid.met <- list(MAE = "MAE", MSE = "MSE")
  use.met <- valid.met[names(valid.met) %in% metrics]
  if (length(use.met) == 0) {
    stop("Metrics provided are not valid")
  }
  ## extract information
  testProbsDeconv <- final.data(object, "test")@metadata[[1]]
  predictionsDeconv <- trained.model(object)@predict.results
  ## results test
  tmd <- as_tibble(x = testProbsDeconv)
  tmd <- mutate(tmd, Sample = rownames(testProbsDeconv),
                nMix = factor(rowSums(testProbsDeconv > 0)))
  tmd <- tmd %>% gather(key = CellType, value = "Prob", -Sample, -nMix)
  ## probabilities target test
  pmd <- as_tibble(predictionsDeconv)
  pmd <- mutate(pmd, Sample = rownames(predictionsDeconv))
  pmd <- pmd %>% gather(key = "CellType", value = "Pred", -Sample)
  ## union
  amd <- tmd %>% left_join(pmd, by = c("Sample", "CellType"))
  ## Add bins to Probs
  amd$pBin <- 0
  for (p in seq(from = 0.1, to = 1, by = 0.1)) {
    amd$pBin[amd$Prob <= p & amd$Prob > p - 0.1] <- p
  }
  amd$pBin[amd$Prob == 0] <- 0.1

  ## calculate stats
  amd <- .updateAMD(amd = amd, use.met = use.met)
  eval.stats <- lapply(use.met, function(x) .calculateMetrics(mat = amd, err = x))
  ## update object
  trained.model(object)@eval.stats.samples <- list(amd, eval.stats)

  return(object)
}

## square error
.SqrErr <- function(x) (x$Prob - x$Pred)**2
## proportional square error
.ppSqrErr <- function(x) x$SqrErr / x$pBin

## absolute error
.AbsErr <- function(x) abs(x$Prob - x$Pred)
## proportional absolute error
.ppAbsErr <- function(x) x$AbsErr / x$pBin

## mean error by
.meanErr <- function(x, err, by, name) {
  if (!err %in% colnames(x))
    stop(paste(err, "does not present"))
  if (is(by, "character")) {
    res <- aggregate(x[[err]], FUN = mean, by = list(by = x[[by]]))
    colnames(res) <- c(by, name)
  } else {
    res <- aggregate(x[[err]], FUN = mean, by = by)
    colnames(res) <- c(names(by), name)
  }
  return(res)
}

.statsBlock <- function(x, err, by) {
  if (err == "MAE") {
    err.x <- c("AbsErr", "ppAbsErr")
    name <- c("MAE", "ppMAE")
  } else if (err == "MSE") {
    err.x <- c("SqrErr", "ppSqrErr")
    name <- c("MSE", "ppMSE")
  }
  nor.err <- .meanErr(x = x, err = err.x[1], by = by, name = name[1])
  pp.err <- .meanErr(x = x, err = err.x[2], by = by, name = name[2])
  if (is(by, "character"))
    d.err <- inner_join(x = nor.err, y = pp.err, by = by)
  else
    d.err <- inner_join(x = nor.err, y = pp.err, by = names(by))
  return(d.err)
}

.calculateMetrics <- function(mat, err) {
  pBinNMix <- list(pBinNMix = paste(mat$pBin, mat$nMix, sep = "_"))
  by.stats <- list(Sample = "Sample", pBin = "pBin", nMix = "nMix", pBinNMix = pBinNMix)
  list.stats <- lapply(X = by.stats, function(x) .statsBlock(mat, err = err, by = x))
  return(list.stats)
}

.updateAMD <- function(amd, use.met) {
  for (i in names(use.met)) {
    if (i == "MAE") {
      amd$AbsErr <- .AbsErr(amd)
      amd$ppAbsErr <- .ppAbsErr(amd)
    } else if (i == "MSE") {
      amd$SqrErr <- .SqrErr(amd)
      amd$ppSqrErr <- .ppSqrErr(amd)
    }
  }
  return(amd)
}

## error: AbsErr, ppAbsErr, SqrErr, ppSqrErr
violinError <- function(
  object,
  error = "AbsErr",
  facet.by = "CellType",
  color.by = "nMix",
  colors = color.list,
  ylimit = 1,
  ncol = NULL,
  title = NULL
) {
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (is.null(trained.model(object)) ||
             is.null(trained.model(object)@eval.stats.samples)) {
    stop("The provided object does not have evaluation metrics. Use ",
         "'calculateEvalMetrics'")
  } else if (!is(trained.model(object)@eval.stats.samples[[1]], "tbl_df")) {
    stop("Evaluation metrics are not well built, use 'calculateEvalMetrics'")
  } else if (!error %in% c("AbsErr", "ppAbsErr", "SqrErr", "ppSqrErr")) {
    stop("Error provided is not valid. Errors available are: AbsErr, ",
         "ppAbsErr, SqrErr, ppSqrErr")
  }
  err <- error

  amd <- trained.model(object)@eval.stats.samples[[1]]
  if (is.null(title))
    title <- paste(err, "by Probability Bin")

  plot <- ggplot(amd, aes(x = factor(pBin), y = .data[[err]])) +
    geom_point(size = 0.1, aes(colour = amd[[color.by]]), position = "jitter") +
    geom_violin(trim = TRUE, scale = "width", fill = NA) +
    scale_color_manual(values = colors) +
    # ylim(0, ylimit) +

    ggtitle(title) + xlab("Probability") +  ylab("err") + labs(fill = color.by) +
    facet_wrap( ~ amd[[facet.by]], ncol = ncol) +
    theme(axis.text.x = element_text(size = 8, angle = 45),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.title = element_text(face = "bold"))
  return(plot)
}
