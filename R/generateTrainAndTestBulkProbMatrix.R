library(ggplot2)




generateTrainAndTestBulkProbMatrix <- function(
  object,
  cell.type.column,
  prob.design,
  train.freq = 2/3,
  exclusive.types = NULL,
  num.bulk.samples = NULL,
  seed = NULL,
  verbose = TRUE
) {
  if (class(object) != "DigitalDLSorter") {
    stop("The object provided is not of DigitalDLSorter class")
  } else if (is.null(single.cell.sim(object))) {
    stop("single.cell.sim slot is empty")
  } else if (!train.freq <= 0.95 | !train.freq >= 0.05) {
    stop("train.seq argument must be less than or equal to 0.95 and greater than or equal to 0.05")
  } else if (!is.data.frame(prob.design)) {
    stop(paste("prob.design must be a data.frame with three column names:",
               "cell.type.column: must be equal to cell.type.column in cells.metadata (colData slot of single.cell.sim)",
               "from: frequency from which the cell type can appear",
               "to: frequency up to which the cell type can appear", sep = "\n   - "))
  }
  if (!is.null(prob.matrix(object)) | !length(prob.matrix(object)) == 0) {
    warning("prob.matrix slot already has the probability matrices. Note that it will be overwritten\n",
            call. = FALSE, immediate. = TRUE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # extract data from SCE to list
  list.data <- .extractDataFromSCE(SCEobject = single.cell.sim(object),
                                   filtering = FALSE)

  # check if cell.type.column is correct
  .checkColumn(metadata = list.data[[2]],
               ID.column = cell.type.column,
               type.metadata = "cells.metadata",
               arg = "cell.type.column")

  # check if prob.design is correctly built
  lapply(X = c(cell.type.column, "from", "to"),
         FUN = function(x) {
           .checkColumn(metadata = prob.design,
                        ID.column = x,
                        type.metadata = "prob.design",
                        arg = "")
         }
  )
  if (any(duplicated(prob.design[, cell.type.column]))) {
    stop(paste("prob.design must not contain duplicated cell types in",
               cell.type.column, "column"))
  } else if (!any(prob.design[, cell.type.column] %in%
              unique(list.data[[2]][, cell.type.column]))) {
    stop("There are some cell types in prob.design that does not appear in cells.metadata. Check that the prob.design matrix is correctly built")
  } else if (any(prob.design$from < 0) | any(prob.design$from > 99)) {
    stop("from column in prob.design must be greater than or equal to 0 and lesser than or equal to 99")
  } else if (any(prob.design$to < 1) | any(prob.design$to > 100)) {
    stop("to column in prob.design must be greater than or equal to 1 and lesser than or equal to 100")
  } else if (any(prob.design$from > prob.design$to)) {
    stop("'from' entries must be lesser than 'to' entries")
  }

  ## set new s.cells if num.bulk.samples is provided
  s.cells <- dim(list.data[[1]])[2]
  if (!is.null(num.bulk.samples)) {
    if (s.cells > num.bulk.samples) {
      stop(paste0("If num.bulk.samples is provided, it must be greater than or equal to the total of simulated cells (",
                  s.cells, " in this case)"))
    }
    s.cells <- ceiling(num.bulk.samples / 17.65)
    message(paste("The number of bulk samples that will be generated has been fixed to", s.cells))
  }

  # split data into training and test sets
  cells <- colnames(list.data[[1]])
  names(cells) <- list.data[[2]][, cell.type.column]

  # train set
  train.set <- sample(cells, size = round(dim(list.data[[1]])[2] * train.freq))
  train.types <- names(train.set)
  train.set.list <- list()
  for (ts in levels(factor(train.types))) {
    train.set.list[[ts]] <- train.set[train.types == ts]
  }

  # test set
  test.set <- cells[!cells %in% train.set]
  test.types <- names(test.set)
  test.set.list <- list()
  for (ts in levels(factor(test.types))) {
    test.set.list[[ts]] <- test.set[test.types == ts]
  }

  if (verbose) {
    message("=== Train Set cells by type:")
    tb <- unlist(lapply(train.set.list, length))
    message(paste0("  - ", names(tb), ": ", tb, collapse = "\n"), "\n")
    message("=== Test Set cells by type:")
    tb <- unlist(lapply(test.set.list, length))
    message(paste0("  - ", names(tb), ": ", tb, collapse = "\n"), "\n")
  }

  prob.list <- apply(X = prob.design,
                      MARGIN = 1,
                      FUN = function (x) {
                        return(seq(from = x['from'], to = x['to']))
                        }
 )
  names(prob.list) <- prob.design[, cell.type.column]

  ## check if there are exclusive types
  if (!is.null(exclusive.types)) {
    if (length(exclusive.types) < 2) {
      stop("'exclusive.types' must be at least 2 different cell types")
    } else if (any(duplicated(exclusive.types))) {
      stop("'exclusive.types' can not contain duplicated elements")
    } else if (!all(exclusive.types %in% unique(list.data[[2]][, "Cell_type"]))) {
      stop("Cell types present in exclusive.types argument must be present in cells.metadata")
    }
    index.ex <- match(exclusive.types, names(prob.list))
  } else {
    index.ex <- NULL
  }

  ## predefined range of cell types
  df <- reshape::melt(prob.list)
  colnames(df) <- c("Perc", "CellType")
  df$CellType <- factor(df$CellType, levels = names(prob.list))
  plot.prob <- .boxPlot(df = df, title = "Predefine Range of cell fractions",
                        y = Perc)

  n.cell.types <- length(unique(train.types))
  functions.list <- list(.generateSet1, .generateSet2, .generateSet3,
                         .generateSet4, .generateSet5, .generateSet6)

  # TRAIN SETS -----------------------------------------------------------------
  train.prob.matrix <- matrix(rep(0, n.cell.types), nrow = 1, byrow = T)
  train.plots <- list()
  nums <- c(1000, 3000, 500, 1000, 3000, 2000)
  n <- 1
  for (fun in functions.list) {
    train.probs <- fun(prob.list = prob.list,
                       prob.matrix = train.prob.matrix,
                       num = nums[n],
                       s.cells = s.cells,
                       n.cell.types = n.cell.types,
                       index.ex = index.ex)
    train.prob.matrix <- rbind(train.prob.matrix, train.probs)
    if (n == 1) train.prob.matrix <- train.prob.matrix[-1,]
    train.plots[[n]] <- .plotsQCSets(probs = train.probs,
                                     prob.matrix = train.prob.matrix,
                                     n = n,
                                     set = "train")
    n <- n + 1
  }

  if (verbose) {
    message("=== Probability matrix for training data:")
    message(paste(c("  - Bulk samples:", "  - Cell types:"),
                  dim(train.prob.matrix),
                  collapse = "\n"), "\n")
  }

  # TEST SETS ------------------------------------------------------------------
  test.prob.matrix <- matrix(rep(0, n.cell.types), nrow = 1, byrow = T)
  test.plots <- list()
  nums <- c(700, 2000, 350, 700, 2000, 1400)
  n <- 1
  for (fun in functions.list) {
    test.probs <- fun(prob.list = prob.list,
                      prob.matrix = test.prob.matrix,
                      num = nums[n],
                      s.cells = s.cells,
                      n.cell.types = n.cell.types,
                      index.ex = index.ex)
    test.prob.matrix <- rbind(test.prob.matrix, test.probs)
    if (n == 1) test.prob.matrix <- test.prob.matrix[-1,]
    test.plots[[n]] <- .plotsQCSets(probs = test.probs,
                                    prob.matrix = test.prob.matrix,
                                    n = n,
                                    set = "test")
    n <- n + 1
  }

  if (verbose) {
    message("=== Probability matrix for test data:")
    message(paste(c("  - Bulk samples:", "  - Cell types:"),
                  dim(test.prob.matrix),
                  collapse = "\n"), "\n")
  }

  # GENERATE PROBS MATRIX NAMES ------------------------------------------------
  setCount <- function (x, setList, sn) {
    names(x) <- sn
    sc <- c()
    for (cType in names(x)) {
      n <- ceiling(x[cType])
      if (n > 0) {
        repl <- ifelse(n > length(setList[[cType]]), TRUE, FALSE)
        sc <- c(sc, sample(setList[[cType]], size = n, replace = repl))
      }
    }
    return(sc[seq(100)])
  }

  train.prob.matrix.names <- t(apply(X = train.prob.matrix,
                                     MARGIN = 1,
                                     FUN = setCount,
                                     setList = train.set.list,
                                     sn = colnames(train.prob.matrix)))

  test.prob.matrix.names <- t(apply(X = test.prob.matrix,
                                    MARGIN = 1,
                                    FUN = setCount,
                                    setList = test.set.list,
                                    sn = colnames(test.prob.matrix)))

  # generate object of ProbMatrixCellTypes class
  train.prob.matrix.object <- new(
    Class = "ProbMatrixCellTypes",
    prob.matrix = train.prob.matrix,
    cell.names = train.prob.matrix.names,
    set.list = train.set.list,
    set = train.set,
    plots = train.plots,
    type.data = "train"
  )

  test.prob.matrix.object <- new(
    Class = "ProbMatrixCellTypes",
    prob.matrix = test.prob.matrix,
    cell.names = test.prob.matrix.names,
    set.list = test.set.list,
    set = test.set,
    plots = test.plots,
    type.data = "test"
  )
  object@prob.matrix <- list(train = train.prob.matrix.object,
                             test = test.prob.matrix.object)

  message("DONE")
  return(object)
}


.violinPlot <- function(df, title, x = CellType, y = Prob) {
  plot <- ggplot(df, aes(x = {{x}}, y = {{y}})) +
    geom_violin() + ggtitle(title)
  return(plot)
}

.boxPlot <- function(df, title, x = CellType, y = Prob) {
  plot <- ggplot(df, aes(x = {{x}}, y = {{y}})) +
    geom_boxplot() + ggtitle(title)
  return(plot)
}

.linesPlot <- function(df, title, x = CellType, y = Prob, group = Sample) {
  plot <- ggplot(df,aes(x = {{x}}, y = {{y}}, group = {{group}})) +
    geom_line(colour = "grey60") + ggtitle(title)
  return(plot)
}

.plotsQCSets <- function(probs, prob.matrix, n, set) {
  title <- paste0("Bulk Probability Dist. Set ", n, " (", set, ")")
  plots.functions <- list(.violinPlot, .boxPlot, .linesPlot)
  df <- reshape::melt(probs)
  colnames(df) <- c("Sample", "CellType", "Prob")
  # first three plots
  plot.list <- lapply(plots.functions, function(f) f(df, title))
  # final plots <-- preguntar por los nombres de los índices
  dummy <- t(apply(prob.matrix, 1, sort, decreasing = T))
  df <- reshape::melt(dummy)
  colnames(df) <- c("Sample", "MaxProb", "Prob")
  df$MaxProb <- factor(df$MaxProb)
  plot.list[[4]] <- .boxPlot(df = df, x = MaxProb, title = title)
  names(plot.list) <- c("violinplot", "boxplot", "linesplot", "maxprob")
  return(plot.list)
}


.cellExcluder <- function(vec, index.ex) {
  sel <- sample(index.ex, length(index.ex) - 1)
  vec[sel] <- 0
  return(list(vec, sel))
}

.setHundredLimit <- function(x, index.ex = NULL) {
  if (sum(x) > 100) {
    while (TRUE) {
      if (is.null(index.ex)) {
        sel <- sample(seq(length(x)), 1)
      } else {
        sel <- sample(seq(length(x))[-index.ex], 1)
      }
      res <- x[sel] - abs(sum(x) - 100)
      if (res >= 0) break
    }
    x[sel] <- res
  } else if (sum(x) < 100) {
    while (TRUE) {
      if (is.null(index.ex)) {
        sel <- sample(seq(length(x)), 1)
      } else {
        sel <- sample(seq(length(x))[-index.ex], 1)
      }
      res <- x[sel] + abs(sum(x) - 100)
      if (res <= 100) break
    }
    x[sel] <- res
  }
  return(x)
}


.adjustHundred <- function(
  x,
  prob.list,
  index.ex = NULL,
  sampling = TRUE
) {
  # w <- rowMedians(as.matrix(prob.design[ind, c("from", "to")]))
  w <- unlist(lapply(prob.list, sample, 1))
  # remove cell types if needed
  if (!is.null(index.ex)) {
    x.list <- .cellExcluder(vec = x, index.ex = index.ex)
    x <- x.list[[1]]
    w[x.list[[2]]] <- 0
  }
  d <- abs(sum(x) - 100)
  if (sum(x) > 100) {
    div.w <- (w / sum(w)) * d
    while (!all(x >= div.w)) {
      index <- which(!x >= div.w)
      w[index] <- 0
      div.w <- (w / sum(w)) * d
    }
    x <- round(x - div.w)
    x <- .setHundredLimit(x = x, index.ex = index.ex)
  } else if (sum(x) < 100) {
    div.w <- (w / sum(w)) * d
    while (!all(100 - div.w > x)) {
      index <- which(!100 - div.w > x)
      w[index] <- 0
      div.w <- (w / sum(w)) * d
    }
    x <- round(x + div.w)
    x <- .setHundredLimit(x = x, index.ex = index.ex)
  }
  return(x)
}


.generateSet1 <- function(
  prob.list,
  prob.matrix,
  num,
  s.cells,
  n.cell.types,
  index.ex
) {
  if (!is.null(index.ex)) {
    sampling <- function(prob.list) {
      x <- .cellExcluder(
        vec = unlist(lapply(X = prob.list, FUN = sample, 1)),
        index.ex = index.ex
      )
      return(x[[1]])
    }
  } else {
    sampling <- function(prob.list) unlist(lapply(X = prob.list, FUN = sample, 1))
  }
  n <- ceiling(num * s.cells/1000)
  while (dim(prob.matrix)[1] <= n) {
    prob.matrix <- rbind(prob.matrix, sampling(prob.list = prob.list))
  }
  prob.matrix <- prob.matrix[-1, ]
  prob.matrix <- round(prob.matrix * 100 / rowSums(prob.matrix))
  prob.matrix <- t(apply(X = prob.matrix, MARGIN = 1, FUN = .setHundredLimit, index.ex))
  return(prob.matrix)
}


.generateSet2 <- function(
  prob.list,
  prob.matrix,
  num,
  s.cells,
  n.cell.types,
  index.ex
) {
  if (!is.null(index.ex)) {
    sampling <- function(p) {
      p <- .cellExcluder(vec = sample(p), index.ex = index.ex)
      return(p[[1]])
    }
  } else {
    sampling <- function(p) sample(p)
  }
  probs <- list()
  n <- ceiling(num * s.cells/1000)
  while (length(probs) < n) {
    p <- rep(0, n.cell.types)
    i <- 1
    while (sum(p) < 100) {
      p[i] <- p[i] + sample(seq(100 - sum(p)), size = 1)
      i <- i + 1
      if (i > n.cell.types) {
        i <- 1
      }
    }
    p <- sampling(p)
    if (sum(p == 0) < n.cell.types) {
      probs[[length(probs) + 1]] <- p
    }
  }
  probs <- matrix(unlist(probs), nrow = n, byrow = T)
  colnames(probs) <- colnames(prob.matrix)
  probs <- round(probs * 100 / rowSums(probs))
  probs <- t(apply(X = probs, MARGIN = 1, FUN = .setHundredLimit, index.ex))
  return(probs)
}


.generateSet3 <- function(
  prob.list,
  prob.matrix,
  num,
  s.cells,
  n.cell.types,
  index.ex
) {
  probs <- list()
  n <- ceiling(num * s.cells/1000)
  while(length(probs) < n) {
    p <- rep(0, n.cell.types)
    names(p) <- names(prob.list)
    i <- 1
    while (sum(p) < 100) {
      dp <- 101
      while (dp > max(prob.list[[i]])) {
        dp <- sample(prob.list[[i]], size = 1)
      }
      p[i] <- dp
      i <- i + 1
      if (i > n.cell.types) i <- 1
    }
    p[1] <- p[1] + 1
    p <- sample(p)
    if (sum(p == 0) < n.cell.types) {
      probs[[length(probs) + 1]] <- p
    }
  }
  probs <- lapply(X = probs, FUN = function (x) return(x[names(prob.list)]))
  probs <- matrix(unlist(probs), nrow = n, byrow = T)
  colnames(probs) <- colnames(prob.matrix)
  probs <- round(probs * 100 / rowSums(probs))
  probs <- t(apply(X = probs, MARGIN = 1,
                   FUN = function(x) {
                     .adjustHundred(
                       x = x,
                       prob.list = prob.list,
                       index.ex = index.ex
                     )
                   }))
  return(probs)
}


.generateSet4 <- function(
  prob.list,
  prob.matrix,
  num,
  s.cells,
  n.cell.types,
  index.ex
) {
  probs <- list()
  n <- ceiling(num * s.cells/1000)
  while(length(probs) < n) {
    p <- rep(0, n.cell.types)
    names(p) <- names(prob.list)
    i <- 1
    while (sum(p) < 100) {
      dp <- sample(prob.list[[i]], size = 1)
      p[i] <- dp
      i <- i + 1
      if (i > n.cell.types) i <- 1
    }
    p[1] <- p[1] + 1
    p <- sample(p)
    if (sum(p == 0) < n.cell.types) {
      probs[[length(probs) + 1]] <- p
    }
  }
  probs <- lapply(X = probs, FUN = sample)
  probs <- matrix(unlist(probs), nrow = n, byrow = T)
  colnames(probs) <- colnames(prob.matrix)
  probs <- round(probs * 100 / rowSums(probs))
  probs <- t(apply(X = probs, MARGIN = 1,
                   FUN = function(x) {
                     .adjustHundred(
                       x = x,
                       prob.list = prob.list,
                       index.ex = index.ex
                     )
                   }))
  return(probs)
}

.generateSet5 <- function(
  prob.list,
  prob.matrix,
  num,
  s.cells,
  n.cell.types,
  index.ex
) {
  n <- ceiling(num * s.cells/1000)
  if (!is.null(index.ex)) {
    generator <- function() {
      gtools::rdirichlet(
        n = 1, alpha = .cellExcluder(rep(1, n.cell.types), index.ex = index.ex)[[1]]
      )
    }
    probs <- t(replicate(n = n, expr = generator(), simplify = TRUE))
  } else {
    probs <- gtools::rdirichlet(n, rep(1, n.cell.types))
  }
  probs <- round(probs * 100)
  probs <- t(apply(X = probs, MARGIN = 1,
                   FUN = function(x) .setHundredLimit(
                     x = x, index.ex = index.ex
                   )))
  colnames(probs) <- colnames(prob.matrix)
  return(probs)
}


.generateSet6 <- function(
  prob.list,
  prob.matrix,
  num,
  s.cells,
  n.cell.types,
  index.ex
) {
  probs <- list()
  n <- ceiling(num * s.cells/1000)
  while (length(probs) < n) {
    probs[[length(probs) + 1]] <- unlist(lapply(X = prob.list, FUN = sample, 1))
  }
  probs <- lapply(X = probs, FUN = function(x) return(round(x * 100 / sum(x))))
  probs <- lapply(X = probs, FUN = sample)
  probs <- lapply(X = probs, FUN = function(x) x[names(prob.list)])
  probs <- matrix(unlist(probs), nrow = n, byrow = T)
  probs <- t(apply(X = probs, 1, FUN = function(x) {
    .adjustHundred(x = x,
                   prob.list = prob.list,
                   index.ex = index.ex
    )
  }))
  colnames(probs) <- colnames(prob.matrix)
  return(probs)
}
