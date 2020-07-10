library(ggplot2)

.adjustHundred <- function (x) {
  d <- abs(sum(x)-100)
  if (sum(x) > 100) {
    x[which(x >= d)[1]] <- x[which(x >= d)[1]] - d
  } else if (sum(x) < 100) {
    x[which(x < (100 - d))[1]] <- x[which(x < (100-d))[1]] + d
  }
  return(x)
}

# generateTrainAndTestBulkProbMatrix: meter opción para seed
generateTrainAndTestBulkProbMatrix <- function(object,
                                               cell.type.column,
                                               prob.design,
                                               train.freq = 2/3,
                                               num.bulk.samples = NULL,
                                               verbose = TRUE) {
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
    stop(paste("prob.design must not contain repeated cell types in",
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

  s.cells <- dim(list.data[[1]])[2]

  if (!is.null(num.bulk.samples)) {
    if (s.cells > num.bulk.samples) {
      stop(paste0("If num.bulk.samples is provided, it must be greater than or equal to the total of simulated cells (",
                  s.cells, " in this case)"))
    }
    s.cells <- round(num.bulk.samples / 17.65) # hay un márgen de error por el número. Preguntar para cambiar
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

  ## predefined range od cell types
  df <- reshape::melt(prob.list)
  colnames(df) <- c("Perc","CellType")
  df$CellType <- factor(df$CellType, levels = names(prob.list))
  plot.prob <- .boxPlot(df = df, title = "Predefine Range of cell fractions",
                        y = Perc)

  n.cell.types <- length(unique(train.types))

  # TRAIN SETS -----------------------------------------------------------------
  train.plots <- list()
  ## train set 1 ---------------------------------------------------------------
  train.prob.matrix <- .generateSet1(prob.list = prob.list,
                                     num = 1000,
                                     s.cells = s.cells,
                                     n.cell.types = n.cell.types)

  ## plots
  train.plots[[1]] <- .plotsQCSets(probs = train.prob.matrix,
                                   prob.matrix = train.prob.matrix,
                                   n = "1",
                                   set = "train")

  ## train set 2 ---------------------------------------------------------------
  train.probs <- .generateSet2(prob.matrix = train.prob.matrix,
                               num = 3000,
                               s.cells = s.cells,
                               n.cell.types = n.cell.types)

  train.prob.matrix <- rbind(train.prob.matrix, train.probs)

  ## plots (antes los plots se hacían antes de actualizar train.prob.matrix)
  train.plots[[2]] <- .plotsQCSets(probs = train.probs,
                                   prob.matrix = train.prob.matrix,
                                   n = "2",
                                   set = "train")

  ## train set 3 ---------------------------------------------------------------

  ## permanece en un bucle infinito porque p no llega nunca a ser mayor que 99
  ## lo he limitado ocn el número de veces que se produce un while.
  ## además, debido a esto, estaba limitado a 9, por lo que lo he peusto a n.cell.types - 1
  ## sin embargo, así no llega nunca al último tipo celular. Creo que ocurre lo mismo por
  ## el if (sum(p == 0) <= 7) --> por eso en los plots generados con la pipeline hay
  ## varios tipos celulares que en el set 3 tienen 0. No sé qué es lo correcto, preguntar
  ## de momento lo he corregido, coge todos los tipos celulares. Con esta corrección,
  ## de hecho, es probable que no sea necesario lo de limitar el número de veces, porque la
  ## condición sum(p) > 99 se cumplirá. En cualquier caso, creo que está bien dejarlo porque
  ## es curar en salud
  train.probs <- .generateSet3(prob.list = prob.list,
                               prob.matrix = train.prob.matrix,
                               num = 500,
                               s.cells = s.cells,
                               n.cell.types = n.cell.types)

  train.prob.matrix <- rbind(train.prob.matrix, train.probs)

  ## plots
  train.plots[[3]] <- .plotsQCSets(probs = train.probs,
                                   prob.matrix = train.prob.matrix,
                                   n = "3",
                                   set = "train")


  ## train set 4 ---------------------------------------------------------------
  train.probs <- .generateSet4(prob.list = prob.list,
                               prob.matrix = train.prob.matrix,
                               num = 1000,
                               s.cells = s.cells,
                               n.cell.types = n.cell.types)

  train.prob.matrix <- rbind(train.prob.matrix, train.probs)

  ## plots
  train.plots[[4]] <- .plotsQCSets(probs = train.probs,
                                   prob.matrix = train.prob.matrix,
                                   n = "4",
                                   set = "train")


  ## train set 5 ---------------------------------------------------------------
  train.probs <- .generateSet5(prob.matrix = train.prob.matrix,
                               num = 3000,
                               s.cells = s.cells,
                               n.cell.types = n.cell.types)

  train.prob.matrix <- rbind(train.prob.matrix, train.probs)

  ## plots
  train.plots[[5]] <- .plotsQCSets(probs = train.probs,
                                   prob.matrix = train.prob.matrix,
                                   n = "5",
                                   set = "train")


  ## train set 6 ---------------------------------------------------------------
  train.probs <- .generateSet6(prob.list = prob.list,
                               prob.matrix = train.prob.matrix,
                               num = 2000,
                               s.cells = s.cells)

  train.prob.matrix <- rbind(train.prob.matrix, train.probs)

  ## plots
  train.plots[[6]] <- .plotsQCSets(probs = train.probs,
                                   prob.matrix = train.prob.matrix,
                                   n = "6",
                                   set = "train")

  if (verbose) {
    message("=== Probability matrix for training data:")
    message(paste(c("  - Bulk samples:", "  - Cell types:"),
                  dim(train.prob.matrix),
                  collapse = "\n"), "\n")
  }

  # TEST SETS ------------------------------------------------------------------
  test.plots <- list()
  ## test set 1 ---------------------------------------------------------------
  test.prob.matrix <- .generateSet1(prob.list = prob.list,
                                    num = 700,
                                    s.cells = s.cells,
                                    n.cell.types = n.cell.types)

  ## plots
  test.plots[[1]] <- .plotsQCSets(probs = test.prob.matrix,
                                  prob.matrix = test.prob.matrix,
                                  n = "1",
                                  set = "test")

  ## test set 2 ---------------------------------------------------------------
  test.probs <- .generateSet2(prob.matrix = test.prob.matrix,
                              num = 2000,
                              s.cells = s.cells,
                              n.cell.types = n.cell.types)

  test.prob.matrix <- rbind(test.prob.matrix, test.probs)

  ## plots (antes los plots se hacían antes de actualizar train.prob.matrix)
  test.plots[[2]] <- .plotsQCSets(probs = test.probs,
                                  prob.matrix = test.prob.matrix,
                                  n = "2",
                                  set = "test")

  ## test set 3 ---------------------------------------------------------------
  test.probs <- .generateSet3(prob.list = prob.list,
                              prob.matrix = test.prob.matrix,
                              num = 350,
                              s.cells = s.cells,
                              n.cell.types = n.cell.types)

  test.prob.matrix <- rbind(test.prob.matrix, test.probs)

  ## plots
  test.plots[[3]] <- .plotsQCSets(probs = test.probs,
                                  prob.matrix = test.prob.matrix,
                                  n = "3",
                                  set = "test")

  ## test set 4 ---------------------------------------------------------------
  test.probs <- .generateSet4(prob.list = prob.list,
                              prob.matrix = test.prob.matrix,
                              num = 700,
                              s.cells = s.cells,
                              n.cell.types = n.cell.types)

  test.prob.matrix <- rbind(test.prob.matrix, test.probs)

  ## plots
  test.plots[[4]] <- .plotsQCSets(probs = test.probs,
                                  prob.matrix = test.prob.matrix,
                                  n = "4",
                                  set = "test")


  ## test set 5 ---------------------------------------------------------------
  test.probs <- .generateSet5(prob.matrix = test.prob.matrix,
                              num = 2000,
                              s.cells = s.cells,
                              n.cell.types = n.cell.types)

  test.prob.matrix <- rbind(test.prob.matrix, test.probs)

  ## plots
  test.plots[[5]] <- .plotsQCSets(probs = test.probs,
                                  prob.matrix = test.prob.matrix,
                                  n = "5",
                                  set = "test")


  ## test set 6 ---------------------------------------------------------------
  test.probs <- .generateSet6(prob.list = prob.list,
                              prob.matrix = test.prob.matrix,
                              num = 1400,
                              s.cells = s.cells)

  test.prob.matrix <- rbind(test.prob.matrix, test.probs)

  ## plots
  test.plots[[6]] <- .plotsQCSets(probs = test.probs,
                                  prob.matrix = test.prob.matrix,
                                  n = "6",
                                  set = "test")

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
    return(sc[seq(100)]) ## por qué 100? --> número de células por bulk sample?
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


.generateSet1 <- function(prob.list, num, s.cells, n.cell.types) {
  n <- ceiling(num * s.cells/1000)
  probs.matrix <- matrix(rep(0, n.cell.types), nrow = 1, byrow = T)
  while (dim(probs.matrix)[1] < n + 1) {
    probs.matrix <- rbind(probs.matrix, unlist(lapply(X = prob.list,
                                                      FUN = sample, 1)))
  }
  probs.matrix <- probs.matrix[-1, ]
  probs.matrix <- round(probs.matrix * 100 / rowSums(probs.matrix))
  probs.matrix <- t(apply(X = probs.matrix, MARGIN = 1, FUN = .adjustHundred))
  return(probs.matrix)
}

.generateSet2 <- function(prob.matrix, num, s.cells, n.cell.types) {
  probs <- list()
  n <- ceiling(num * s.cells/1000)
  while (length(probs) < n) {
    p <- rep(0, n.cell.types)
    i <- 1
    while (sum(p) < 99) {
      p[i] <- p[i] + sample(seq(100 - sum(p)), size = 1)
      i <- i + 1
      if (i > n.cell.types) {
        i <- 1
      }
    }
    p[0] <- p[0] + 1
    if (sum(p == 0) <= 9) { # hardoced? probablemente debería ser n.cell.types y no 9
      p <- sample(p)
      probs[[length(probs) + 1]] <- p
    }
  }
  probs <- matrix(unlist(probs), nrow = n, byrow = T)
  colnames(probs) <- colnames(prob.matrix)
  probs <- round(probs * 100 / rowSums(probs))
  probs <- t(apply(X = probs, MARGIN = 1, FUN = .adjustHundred))
  return(probs)
}



.generateSet3 <- function(prob.list, prob.matrix, num, s.cells, n.cell.types) {
  probs <- list()
  n <- ceiling(num * s.cells/1000)
  while(length(probs) < n) {
    p <- rep(0, n.cell.types)
    names(p) <- names(prob.list)
    i <- 1
    times <- 0
    while (sum(p) < 99 & times <= s.cells * 5) {
      dp <- 101
      while (dp > max(prob.list[[i]])) {
        dp <- sample(prob.list[[i]], size = 1)
      }
      p[i] <- dp
      i <- i + 1
      times <- times + 1
      if (i > n.cell.types) i <- 1 # corrección
    }
    p[1] <- p[1] + 1
    if (sum(p == 0) <= 7) {
      p <- sample(p)
      probs[[length(probs) + 1]] <- p
    }
    # n <- n+1
  }
  probs <- lapply(X = probs, FUN = function (x) return(x[names(prob.list)]))
  probs <- matrix(unlist(probs), nrow = n, byrow = T)
  colnames(probs) <- colnames(prob.matrix)
  probs <- round(probs * 100 / rowSums(probs))
  probs <- t(apply(X = probs, MARGIN = 1, FUN = .adjustHundred))
  return(probs)
}



.generateSet4 <- function(prob.list, prob.matrix, num, s.cells, n.cell.types) {
  probs <- list()
  n <- ceiling(num * s.cells/1000)
  while(length(probs) < n) {
    p <- rep(0, n.cell.types)
    names(p) <- names(prob.list)
    i <- 1
    while (sum(p) < 99) {
      dp <- sample(prob.list[[i]], size = 1)
      p[i] <- dp
      i <- i+1
      if (i > n.cell.types) i <- 1
    }
    p[1] <- p[1] + 1
    if (sum(p == 0) <= 7) { # preguntar --> habrá que cambiarlo
      p <- sample(p)
      probs[[length(probs) + 1]] <- p
    }
    # n <- n+1
  }
  probs <- lapply(X = probs, FUN = sample)
  probs <- matrix(unlist(probs), nrow = n, byrow = T)
  colnames(probs) <- colnames(prob.matrix)
  probs <- round(probs * 100 / rowSums(probs))
  probs <- t(apply(X = probs, MARGIN = 1, FUN = .adjustHundred))
  return(probs)
}


.generateSet5 <- function(prob.matrix, num, s.cells, n.cell.types) {
  n <- ceiling(num * s.cells/1000)
  probs <- gtools::rdirichlet(n, rep(1, n.cell.types))
  probs <- round(probs * 100)
  probs <- t(apply(X = probs, MARGIN = 1, FUN = .adjustHundred))
  colnames(probs) <- colnames(prob.matrix)
  return(probs)
}


.generateSet6 <- function(prob.list, prob.matrix, num, s.cells) {
  probs <- list()
  n <- ceiling(num * s.cells/1000)
  while (length(probs) < n) {
    probs[[length(probs) + 1]] <- unlist(lapply(X = prob.list, FUN = sample, 1))
  }
  probs <- lapply(X = probs, FUN = function(x) return(round(x * 100 / sum(x))))
  probs <- lapply(X = probs, FUN = sample)
  probs <- lapply(X = probs, FUN = .adjustHundred)
  probs <- matrix(unlist(probs), nrow = n, byrow = T)
  colnames(probs) <- colnames(prob.matrix)
  return(probs)
}
