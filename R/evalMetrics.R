#' @importFrom dplyr mutate as_tibble left_join inner_join
#' @importFrom tidyr gather
#' @importFrom RColorBrewer brewer.pal
NULL



color.list <- c(RColorBrewer::brewer.pal(12, "Paired"),
                RColorBrewer::brewer.pal(12, "Set3"),
                RColorBrewer::brewer.pal(8, "Pastel2"),
                "#333333", "#5D5D5D",
                "#888888", "#B3B3B3")

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
  amdf <- amd %>% filter(Prob > 0 & Prob < 1)
  eval.stats <- lapply(use.met, function(x) .calculateMetrics(mat = amd, err = x))
  eval.stats.f <- lapply(use.met, function(x) .calculateMetrics(mat = amdf, err = x))
  ## update object
  trained.model(object)@eval.stats.samples <- list(amd, eval.stats, eval.stats.f)

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

se <- function(x) sqrt(var(x)/length(x))
# res1 <- aggregate(amd[["AbsErr"]], FUN = sd, by = list(by = amd[["CellType"]]))
# res2 <- aggregate(amd[["AbsErr"]], FUN = se, by = list(by = amd[["CellType"]]))
# Reduce(function(x,y) merge(x = x, y = y, by = by),
#        list(mean.res, sd.res, se.red))

## mean error by
.meanErr <- function(x, err, by, name) {
  if (!err %in% colnames(x))
    stop(paste(err, "does not present"))
  if (is(by, "character")) {
    list.res <- lapply(list(mean, sd, se), function(fun) {
      res <- aggregate(x[[err]], FUN = fun, by = list(by = x[[by]]))
      colnames(res) <- c(by, name)
      return(res)
    })
    res <- Reduce(function(x, y) merge(x = x, y = y, by = by),
                  list.res)
    colnames(res) <- c(by, paste0(name, ".mean"),
                       paste0(name, ".sd"),
                       paste0(name, ".se"))
  } else {
    list.res <- lapply(list(mean, sd, se), function(fun) {
      res <- aggregate(x[[err]], FUN = fun, by = by)
      colnames(res) <- c(names(by), name)
      return(res)
    })
    res <- Reduce(function(x, y) merge(x = x, y = y, by = names(by)),
                  list.res)
    colnames(res) <- c(names(by), paste0(name, ".mean"),
                       paste0(name, ".sd"),
                       paste0(name, ".se"))
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
  by.stats <- list(Sample = "Sample", CellType = "CellType",
                   pBin = "pBin", nMix = "nMix", pBinNMix = pBinNMix)
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
  x.by = "pBin",
  facet.by = "nMix",
  color.by = "nMix",
  filter.sc = TRUE,
  colors = color.list,
  ylimit = NULL,
  nrow = NULL,
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
  } else if (!color.by %in% c("nMix", "CellType")) {
    stop("'color.by' provided is not valid. Options available are: nMix, CellType")
  } else if (!x.by %in% c("nMix", "CellType", "pBin")) {
    stop("'x.by' provided is not valid. Options available are: nMix, CellType, pBin")
  }
  amd <- trained.model(object)@eval.stats.samples[[1]]
  if (filter.sc) {
    amd <- amd %>% filter(Prob > 0 & Prob < 1)
  }
  if (length(colors) < length(unique(amd[[color.by]]))) {
    stop("Colors provided are not enought") ## generador de colores
  }
  if (is.null(title))
    title.plot <- paste(error, "by Probability Bin")
  else
    title.plot <- title
  plot <- ggplot(amd, aes(x = factor(.data[[x.by]]), y = .data[[error]]))
  if (!is.null(facet.by)) {
    if (!facet.by %in% c("nMix", "CellType")) {
      stop("'facet.by' provided is not valid. Options available are: nMix, ",
           "CellType or NULL")
    }
    plot <- plot + facet_wrap(as.formula(paste("~", facet.by)),
                              nrow = nrow, ncol = ncol)
  }
  plot <- plot + geom_point(size = 0.1, aes(colour = amd[[color.by]]),
                            position = "jitter") +
    geom_violin(trim = TRUE, scale = "width", fill = NA) +
    scale_color_manual(values = colors, name = color.by) +
    ggtitle(title.plot) + xlab("Probability") +  ylab(error) +
    # facet_wrap(as.formula(paste("~", facet.by)), nrow = nrow, ncol = ncol) +
    guides(colour = guide_legend(override.aes = list(size = 1.5))) +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.title = element_text(face = "bold"))
  if (!is.null(ylimit)) plot <- plot + ylim(0, ylimit)

  return(plot)
}


corrExpPredPlot <- function(
  object,
  facet.by,
  color.by,
  colors = color.list,
  filter.sc = FALSE,
  ncol = NULL,
  nrow = NULL,
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
  } else if (!color.by %in% c("nMix", "CellType")) {
    stop("'color.by' provided is not valid. Options available are: nMix, CellType")
  }
  amd <- trained.model(object)@eval.stats.samples[[1]]
  if (filter.sc) {
    amd <- amd %>% filter(Prob > 0 & Prob < 1)
  }
  if (length(colors) < length(unique(amd[[color.by]]))) {
    stop("Colors provided are not enought") ## generador de colores
  }
  if (is.null(title))
    title.plot <- "Correlation Expected/Predicted"
  else
    title.plot <- title

  plot <- ggplot(amd, aes(x = Prob, y = Pred))
  if (!is.null(facet.by)) {
    if (!facet.by %in% c("nMix", "CellType")) {
      stop("'facet.by' provided is not valid. Options available are: nMix, ",
           "CellType or NULL")
    }
    plot <- plot + facet_wrap(as.formula(paste("~", facet.by)),
                              nrow = nrow, ncol = ncol)
  }
  plot <- plot + geom_point(size = 0.1, aes(colour = amd[[color.by]]),
                            position = "jitter") +
    geom_abline(linetype = "dashed", colour = "gray40") +
    scale_color_manual(values = colors, name = color.by) +
    scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1)) +
    ggtitle(title.plot) + xlab("Expected") + ylab("Predicted") +
    stat_smooth(method = "lm", colour = "darkblue", alpha = 0.8, size = 0.8) +
    stat_cor(method = "pearson", label.x = 0.01, label.y = 0.95, size = 3) +
    guides(colour = guide_legend(override.aes = list(size = 1.5))) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.title = element_text(face = "bold"))

  return(plot)
}


blandAltmanLehPlot <- function(
  object,
  facet.by,
  color.by,
  log.2 = TRUE,
  colors = color.list,
  filter.sc = FALSE,
  density = TRUE,
  ncol = NULL,
  nrow = NULL,
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
  } else if (!color.by %in% c("nMix", "CellType")) {
    stop("'color.by' provided is not valid. Options available are: nMix, CellType")
  }
  amd <- trained.model(object)@eval.stats.samples[[1]]
  if (filter.sc) {
    amd <- amd %>% filter(Prob > 0 & Prob < 1)
  }
  if (log.2) {
    amd <- amd %>% mutate(
      Mean = (log2(Prob + 0.001) + log2(Pred + 0.001)) / 2,
      Diff = log2(Pred + 0.001) - log2(Prob + 0.001)
    )
    add.title <- "(log2 space)"
    x.lab <- "(log2(Pred) + log2(Exp))/2"
    y.lab <- "log2(Pred/Exp)"
  } else {
    amd <- amd %>% mutate(
      Mean = (Prob + Pred) / 2,
      Diff = (Prob - Pred)
    )
    add.title <- "(normal space)"
    x.lab <- "(Pred + Exp)/2"
    y.lab <- "Pred/Exp"
  }
  if (is.null(title))
    title.plot <- paste("Bland-Altman Agreement Plot", add.title)
  else
    title.plot <- title

  plot <- ggplot(amd, aes(x = Mean, y = Diff, colour = amd[[color.by]])) +
    geom_point(size = 0.05) +
    scale_color_manual(values = color.list, name = color.by) +
    scale_x_continuous(labels = c(-10, -7.5, -5, -2.5, 0)) +
    geom_hline(aes(yintercept = mean(Diff)), linetype = "dashed") +
    geom_hline(aes(yintercept = mean(Diff) + 2 * sd(Diff)),
               linetype = "dashed", colour = "red") +
    geom_hline(aes(yintercept = mean(Diff) - 2 * sd(Diff)),
               linetype = "dashed", colour = "red") +
    xlab(x.lab) + ylab(y.lab) +
    facet_wrap(as.formula(paste("~", facet.by))) +
    ggtitle(title.plot) +
    guides(colour = guide_legend(override.aes = list(size = 1.5))) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.title = element_text(face = "bold"))
  if (density)
    plot <- plot + stat_density_2d(colour = "black", alpha = 0.5, linetype = "dashed")

  return(plot)
}


barErrorPlot <- function(
  object,
  error,
  by,
  dispersion = "se",
  filter.sc = TRUE,
  colors = color.list,
  title = NULL,
  angle = 90
) {
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (is.null(trained.model(object)) ||
             is.null(trained.model(object)@eval.stats.samples)) {
    stop("The provided object does not have evaluation metrics. Use ",
         "'calculateEvalMetrics'")
  } else if (!is(trained.model(object)@eval.stats.samples[[1]], "tbl_df")) {
    stop("Evaluation metrics are not well built, use 'calculateEvalMetrics'")
  } else if (!by %in% c("nMix", "CellType")) {
    stop("'by' provided is not valid. Options available are: nMix, CellType")
  } else if (!error %in% c("MAE", "MSE")) {
    stop("Error provided is not valid. Errors available are: AbsErr, ",
         "ppAbsErr, SqrErr, ppSqrErr")
  } else if (!dispersion %in% c("se", "sd")) {
    stop("Dispersion provided is not valid. Options available are: sd (standard",
         " deviation) or se (standard error)")
  }
  if (is.null(title))
    title.plot <- paste0("Bar Error plot by ", by, " (",error, ")")
  else
    title.plot <- title
  ## filter sc data
  if (!filter.sc) index.stats <- 2
  else index.stats <- 3

  data <- trained.model(object)@eval.stats.samples[[index.stats]][[error]][[by]]
  err.mean <- paste0(error, ".mean")
  err.dis <- paste0(error, ".", dispersion)

  plot <- ggplot(data, aes(x = .data[[by]],
                           y = .data[[err.mean]],
                           ymin = .data[[err.mean]] - .data[[err.dis]],
                           ymax = .data[[err.mean]] + .data[[err.dis]])) +
    geom_errorbar(width = 0.2) + geom_point(size = 1.5) +
    xlab(by) + ylab(error) + ggtitle(title.plot) +
    theme(axis.text.x = element_text(size = 8, angle = angle, hjust = 1, vjust = 0.5),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.title = element_text(face = "bold"))

  return(plot)
}

