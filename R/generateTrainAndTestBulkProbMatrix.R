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
                                               train.freq = 2/3) {
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
  # split data into training and test sets
  cells <- colnames(list.data[[1]])
  names(cells) <- list.data[[2]][, cell.type.column]
  train.set <- sample(cells, size = round(dim(list.data[[1]])[2] * train.freq))
  train.types <- names(train.set)
  # train.types <- sub("_\\w+", "", train.set, perl = T)
  train.set.list <- list()
  for (ts in levels(factor(train.types))) {
    train.set.list[[ts]] <- train.set[train.types == ts]
  }
  message("Train Set cells by type\n")
  tb <- unlist(lapply(train.set.list, length))
  print(tb)

  test.set <- cells[!cells %in% train.set]
  test.types <- names(test.set)
  # test.types <- sub("_\\w+", "", test.set, perl = T)
  test.set.list <- list()
  for (ts in levels(factor(test.types))) {
    test.set.list[[ts]] <- test.set[test.types == ts]
  }
  message("Test Set cells by type\n")
  tb <- unlist(lapply(test.set.list, length))
  print(tb)

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
  s.cells <- dim(list.data[[1]])[2]
  n.cell.types <- length(unique(train.types))

  # TRAIN SETS -----------------------------------------------------------------

  ## train set 1 ---------------------------------------------------------------
  train.prob.matrix <- .generateSet1(prob.list = prob.list,
                                     num = 1000,
                                     s.cells = s.cells,
                                     n.cell.types = n.cell.types)

  ## plots
  train.plots.1 <- .plotsQCSets(probs = train.prob.matrix,
                                title = "Bulk Probability Dist. Set 1",
                                prob.matrix = train.prob.matrix)

  ## train set 2 ---------------------------------------------------------------
  train.probs <- .generateSet2(prob.matrix = train.prob.matrix,
                               num = 3000,
                               s.cells = s.cells,
                               n.cell.types = n.cell.types)

  train.prob.matrix <- rbind(train.prob.matrix, train.probs)

  ## plots (antes los plots se hacían antes de actualizar train.prob.matrix)
  train.plots.2 <- .plotsQCSets(probs = train.probs,
                                title = "Bulk Probability Dist. Set 2",
                                prob.matrix = train.prob.matrix)

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
  train.plots.3 <- .plotsQCSets(probs = train.probs,
                                title = "Bulk Probability Dist. Set 3",
                                prob.matrix = train.prob.matrix)


  ## train set 4 ---------------------------------------------------------------
  train.probs <- .generateSet4(prob.list = prob.list,
                               prob.matrix = train.prob.matrix,
                               num = 1000,
                               s.cells = s.cells,
                               n.cell.types = n.cell.types)

  train.prob.matrix <- rbind(train.prob.matrix, train.probs)

  ## plots
  train.plots.4 <- .plotsQCSets(probs = train.probs,
                                title = "Bulk Probability Dist. Set 4",
                                prob.matrix = train.prob.matrix)


  ## train set 5 ---------------------------------------------------------------
  train.probs <- .generateSet5(prob.matrix = train.prob.matrix,
                               num = 3000,
                               s.cells = s.cells,
                               n.cell.types = n.cell.types)

  train.prob.matrix <- rbind(train.prob.matrix, train.probs)

  ## plots
  train.plots.5 <- .plotsQCSets(probs = train.probs,
                                title = "Bulk Probability Dist. Set 5",
                                prob.matrix = train.prob.matrix)


  ## train set 6 ---------------------------------------------------------------
  train.probs <- .generateSet6(prob.list = prob.list,
                               prob.matrix = train.prob.matrix,
                               num = 2000,
                               s.cells = s.cells)

  train.prob.matrix <- rbind(train.prob.matrix, train.probs)

  ## plots
  train.plots.6 <- .plotsQCSets(probs = train.probs,
                                title = "Bulk Probability Dist. Set 6",
                                prob.matrix = train.prob.matrix)

  train.plots <- list(train.plots.1, train.plots.2, train.plots.3,
                      train.plots.4, train.plots.5, train.plots.6)

  print(paste(c("bulks","types"),dim(trainProbsMatrix)))

  # TEST SETS ------------------------------------------------------------------

  ## test set 1 ---------------------------------------------------------------
  test.prob.matrix <- .generateSet1(prob.list = prob.list,
                                    num = 700,
                                    s.cells = s.cells,
                                    n.cell.types = n.cell.types)

  ## plots
  test.plots.1 <- .plotsQCSets(probs = test.prob.matrix,
                               title = "Bulk Probability Dist. Set 1",
                               prob.matrix = test.prob.matrix)

  ## test set 2 ---------------------------------------------------------------
  test.probs <- .generateSet2(prob.matrix = test.prob.matrix,
                              num = 2000,
                              s.cells = s.cells,
                              n.cell.types = n.cell.types)

  test.prob.matrix <- rbind(test.prob.matrix, test.probs)

  ## plots (antes los plots se hacían antes de actualizar train.prob.matrix)
  test.plots.2 <- .plotsQCSets(probs = test.probs,
                               title = "Bulk Probability Dist. Set 2",
                               prob.matrix = test.prob.matrix)

  ## test set 3 ---------------------------------------------------------------
  test.probs <- .generateSet3(prob.list = prob.list,
                              prob.matrix = test.prob.matrix,
                              num = 350,
                              s.cells = s.cells,
                              n.cell.types = n.cell.types)

  test.prob.matrix <- rbind(test.prob.matrix, test.probs)

  ## plots
  test.plots.3 <- .plotsQCSets(probs = test.probs,
                                title = "Bulk Probability Dist. Set 3",
                                prob.matrix = test.prob.matrix)


  ## test set 4 ---------------------------------------------------------------
  test.probs <- .generateSet4(prob.list = prob.list,
                              prob.matrix = test.prob.matrix,
                              num = 700,
                              s.cells = s.cells,
                              n.cell.types = n.cell.types)

  test.prob.matrix <- rbind(test.prob.matrix, test.probs)

  ## plots
  test.plots.4 <- .plotsQCSets(probs = test.probs,
                                title = "Bulk Probability Dist. Set 4",
                                prob.matrix = test.prob.matrix)


  ## test set 5 ---------------------------------------------------------------
  test.probs <- .generateSet5(prob.matrix = test.prob.matrix,
                              num = 2000,
                              s.cells = s.cells,
                              n.cell.types = n.cell.types)

  test.prob.matrix <- rbind(test.prob.matrix, test.probs)

  ## plots
  test.plots.5 <- .plotsQCSets(probs = test.probs,
                               title = "Bulk Probability Dist. Set 5",
                               prob.matrix = test.prob.matrix)


  ## test set 6 ---------------------------------------------------------------
  test.probs <- .generateSet6(prob.list = prob.list,
                              prob.matrix = test.prob.matrix,
                              num = 1400,
                              s.cells = s.cells)

  test.prob.matrix <- rbind(test.prob.matrix, test.probs)

  ## plots
  test.plots.6 <- .plotsQCSets(probs = test.probs,
                                title = "Bulk Probability Dist. Set 6",
                                prob.matrix = test.prob.matrix)

  print(paste(c("bulks","types"),dim(testProbsMatrix)))

  test.plots <- list(test.plots.1, test.plots.2, test.plots.3,
                     test.plots.4, test.plots.5, test.plots.6)

  # GENERATE PROBS MATRIX NAMES ------------------------------------------------
  cat(paste(names(train.set.list)),"\n")

  setCount <- function (x, setList, sn) {
    names(x) <- sn
    sc <- c()
    for (cType in names(x)) {
      n <- ceiling(x[cType])
      if (n > 0) {
        # if (cType == "Tumour") {cType <- sample(c("HER2","luminalA","luminalB","TNBC"),1)}
        repl <- ifelse(n > length(setList[[cType]]), TRUE, FALSE)
        sc <- c(sc, sample(setList[[cType]], size = n, replace = repl))
      }
    }
    return(sc[seq(100)]) ## por qué 100? --< número de células por bulk sample?
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

.plotsQCSets <- function(probs, prob.matrix, title) {
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
    probs[[length(probs) + 1]] <- unlist(lapply(X = prob.list,
                                                FUN = sample, 1))
  }
  probs <- lapply(X = probs, FUN = function(x) return(round(x * 100 / sum(x))))
  probs <- lapply(X = probs, FUN = sample)
  probs <- lapply(X = probs, FUN = .adjustHundred)
  probs <- matrix(unlist(probs), nrow = n, byrow = T)
  colnames(probs) <- colnames(prob.matrix)
  return(probs)
}
