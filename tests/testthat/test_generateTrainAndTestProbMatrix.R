context("Generation of cell composition matrix")

## set object with all information needed for generating prob matrix

DDLSChungSmall <- simSingleCellProfiles(
  object = DDLSChungSmall,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_type",
  n.cells = 10,
  verbose = TRUE
)
## set valid prob.matrix
probMatrixValid <- data.frame(
  Cell_type = c("ER+", "HER2+", "ER+ and HER2+", "TNBC",
                "Stromal", "Monocyte", "Tme", "BGC",
                "Bmem", "DC", "Macrophage", "TCD8", "Treg"),
  from = c(rep(30, 4), 1, rep(0, 8)),
  to = c(rep(70, 4), 50, rep(15, 8))
)


## check that object contains all information needed

test_that("Wrong object: single.cell.final missing || Wrong column cell type", {
  ## incorrect object
  DDLSChungSmallBad <- DDLSChungSmall
  single.cell.final(DDLSChungSmallBad) <- NULL
  expect_error(generateTrainAndTestBulkProbMatrix(
    object = DDLSChungSmallBad,
    cell.type.column = "Cell_type",
    prob.design = probMatrixValid,
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "single.cell.final slot is empty")
  ## incorrect column
  expect_error(generateTrainAndTestBulkProbMatrix(
    object = DDLSChungSmall,
    cell.type.column = "non_existent_column",
    prob.design = probMatrixValid,
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "non_existent_column column is not present in cells.metadata")
})



## check if prob.design is correctly built

test_that("Wrong prob.design", {
  ## incorrect cell type
  probMatrixInvalid <- data.frame(
    Cell_type = c("incorrect", "HER2+", "ER+ and HER2+", "TNBC",
                  "Stromal", "Monocyte", "Tme", "BGC",
                  "Bmem", "DC", "Macrophage", "TCD8", "Treg"),
    from = c(rep(30, 4), 1, rep(0, 8)),
    to = c(rep(70, 4), 50, rep(15, 8))
  )
  expect_error(generateTrainAndTestBulkProbMatrix(
    object = DDLSChungSmall,
    cell.type.column = "Cell_type",
    prob.design = probMatrixInvalid,
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "There are some cell types")

  ## new cell types
  probMatrixInvalid <- data.frame(
    Cell_type = c("ER+", "HER2+", "ER+ and HER2+", "TNBC",
                  "Stromal", "Monocyte", "Tme", "BGC",
                  "Bmem", "DC", "Macrophage", "TCD8", "Treg",
                  "new.type"),
    from = c(rep(30, 4), 1, rep(0, 8), 10),
    to = c(rep(70, 4), 50, rep(15, 8), 40)
  )
  expect_error(generateTrainAndTestBulkProbMatrix(
    object = DDLSChungSmall,
    cell.type.column = "Cell_type",
    prob.design = probMatrixInvalid,
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "There are some cell types")

  ## duplicates
  probMatrixInvalid <- data.frame(
    Cell_type = c("ER+", "HER2+", "ER+ and HER2+", "TNBC",
                  "Stromal", "Monocyte", "Tme", "BGC",
                  "Bmem", "DC", "Macrophage", "TCD8", "Treg",
                  "ER+"),
    from = c(rep(30, 4), 1, rep(0, 8), 10),
    to = c(rep(70, 4), 50, rep(15, 8), 40)
  )
  expect_error(generateTrainAndTestBulkProbMatrix(
    object = DDLSChungSmall,
    cell.type.column = "Cell_type",
    prob.design = probMatrixInvalid,
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "prob.design must not contain duplicated cell types")

  ## from greater than to
  probMatrixInvalid <- data.frame(
    Cell_type = c("ER+", "HER2+", "ER+ and HER2+", "TNBC",
                  "Stromal", "Monocyte", "Tme", "BGC",
                  "Bmem", "DC", "Macrophage", "TCD8", "Treg"),
    from = c(rep(30, 4), 1, rep(0, 8)),
    to = c(rep(10, 4), 50, rep(15, 8))
  )
  expect_error(generateTrainAndTestBulkProbMatrix(
    object = DDLSChungSmall,
    cell.type.column = "Cell_type",
    prob.design = probMatrixInvalid,
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "'from' entries must be lesser than 'to' entries")

  ## prob ranges incorrect
  probMatrixInvalid <- data.frame(
    Cell_type = c("ER+", "HER2+", "ER+ and HER2+", "TNBC",
                  "Stromal", "Monocyte", "Tme", "BGC",
                  "Bmem", "DC", "Macrophage", "TCD8", "Treg"),
    from = c(rep(11, 4), 1, rep(0, 8)),
    to = c(rep(90, 4), 50, rep(15, 8))
  )
  expect_error(generateTrainAndTestBulkProbMatrix(
    object = DDLSChungSmall,
    cell.type.column = "Cell_type",
    prob.design = probMatrixInvalid,
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "The sum between")
})


## check proportions arguments

test_that("Wrong proportion arguments", {
  ## incorrect number of elements
  expect_error(generateTrainAndTestBulkProbMatrix(
    object = DDLSChungSmall,
    cell.type.column = "Cell_type",
    prob.design = probMatrixValid,
    proportions.train = c(1, 99),
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "Proportions must be a vector of six elements")

  ## not add 100
  expect_error(generateTrainAndTestBulkProbMatrix(
    object = DDLSChungSmall,
    cell.type.column = "Cell_type",
    prob.design = probMatrixValid,
    proportions.test = c(10, 5, 20, 15, 10, 42),
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "Proportions provided must add up to 100")

  ## not add 100: negative numbers
  expect_error(generateTrainAndTestBulkProbMatrix(
    object = DDLSChungSmall,
    cell.type.column = "Cell_type",
    prob.design = probMatrixValid,
    proportions.test = c(10, 5, 60, 15, 50, -40),
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "Proportions provided must add up to 100")

  ## not add 100: negative numbers
  expect_error(generateTrainAndTestBulkProbMatrix(
    object = DDLSChungSmall,
    cell.type.column = "Cell_type",
    prob.design = probMatrixValid,
    proportions.train = c(0, 5, 20, 15, 10, 50),
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "Proportions can not be equal to or lesser than zero")

})

## check n.cells

test_that(
  "Check number samples and cells: argument control and expected output",
  {
  ## n.cells lesser than n cell types
  expect_error(generateTrainAndTestBulkProbMatrix(
    object = DDLSChungSmall,
    cell.type.column = "Cell_type",
    prob.design = probMatrixValid,
    num.bulk.samples = 200,
    n.cells = 12,
    verbose = TRUE
  ), regexp = "n.cells must be equal to or greater than the number of")

  ## num.bulk.samples lesser than n cell types
  expect_error(generateTrainAndTestBulkProbMatrix(
    object = DDLSChungSmall,
    cell.type.column = "Cell_type",
    prob.design = probMatrixValid,
    num.bulk.samples = 100,
    verbose = TRUE
  ), regexp = "If num.bulk.samples is provided")

  ## dim samples <- 1000 (600 train and 400 test) || n.cells <- 250 ------------
  DDLS.1 <- generateTrainAndTestBulkProbMatrix(
    object = DDLSChungSmall,
    cell.type.column = "Cell_type",
    prob.design = probMatrixValid,
    num.bulk.samples = 1000,
    n.cells = 250,
    verbose = TRUE
  )
  ## number of bulk samples
  # train matrix
  cell.train.matrix <- prob.cell.types(DDLS.1, "train") %>% prob.matrix
  expect_equal(nrow(cell.train.matrix), 600)
  expect_true(all(rowSums(cell.train.matrix) == 100))
  # test matrix
  cell.test.matrix <- prob.cell.types(DDLS.1, "test") %>% prob.matrix
  expect_equal(nrow(cell.test.matrix), 400)
  expect_true(all(rowSums(cell.test.matrix) == 100))

  ## number of cell types
  # train
  n.cell.train <- prob.cell.types(DDLS.1, "train") %>% cell.names
  expect_equal(dim(n.cell.train), c(600, 250))
  # test
  n.cell.test <- prob.cell.types(DDLS.1, "test") %>% cell.names
  expect_equal(dim(n.cell.test), c(400, 250))
  # any shared cell between train and test
  expect_false(any(n.cell.train %in% n.cell.test))


  ## with random numbers -------------------------------------------------------
  set.seed(123)
  n.cells <- ceiling(runif(n = 1, min = 100, max = 500))
  num.bulk.samples <- ceiling(runif(n = 1, min = 200, max = 5000))
  DDLS.2 <- generateTrainAndTestBulkProbMatrix(
    object = DDLSChungSmall,
    cell.type.column = "Cell_type",
    prob.design = probMatrixValid,
    num.bulk.samples = num.bulk.samples,
    n.cells = n.cells,
    verbose = TRUE
  )
  ## number of bulk samples
  # total matrix
  cell.train.matrix <- prob.cell.types(DDLS.2, "train") %>% prob.matrix
  cell.test.matrix <- prob.cell.types(DDLS.2, "test") %>% prob.matrix
  cell.total.matrix <- rbind(cell.train.matrix, cell.test.matrix)
  expect_equal(nrow(cell.total.matrix), num.bulk.samples)
  expect_true(all(rowSums(cell.total.matrix) == 100))

  ## number of cell types
  # train
  n.cell.train <- prob.cell.types(DDLS.2, "train") %>% cell.names
  expect_equal(ncol(n.cell.train), n.cells)
  # test
  n.cell.test <- prob.cell.types(DDLS.2, "test") %>% cell.names
  expect_equal(ncol(n.cell.test), n.cells)
  # any shared cell between train and test
  expect_false(any(n.cell.train %in% n.cell.test))


  ## with random numbers and changing proportions ------------------------------
  n.cells <- ceiling(runif(n = 1, min = 100, max = 500))
  num.bulk.samples <- ceiling(runif(n = 1, min = 200, max = 5000))
  DDLS.3 <- generateTrainAndTestBulkProbMatrix(
    object = DDLSChungSmall,
    cell.type.column = "Cell_type",
    prob.design = probMatrixValid,
    num.bulk.samples = num.bulk.samples,
    proportions.train = c(10, 20, 1, 9, 50, 10),
    proportions.test = c(50, 30, 1, 9, 5, 5),
    n.cells = n.cells,
    verbose = TRUE
  )
  ## number of bulk samples
  # total matrix
  cell.train.matrix <- prob.cell.types(DDLS.3, "train") %>% prob.matrix
  cell.test.matrix <- prob.cell.types(DDLS.3, "test") %>% prob.matrix
  cell.total.matrix <- rbind(cell.train.matrix, cell.test.matrix)
  expect_equal(nrow(cell.total.matrix), num.bulk.samples)
  expect_true(all(rowSums(cell.total.matrix) == 100))

  ## number of cell types
  # train
  n.cell.train <- prob.cell.types(DDLS.3, "train") %>% cell.names
  expect_equal(ncol(n.cell.train), n.cells)
  # test
  n.cell.test <- prob.cell.types(DDLS.3, "test") %>% cell.names
  expect_equal(ncol(n.cell.test), n.cells)
  # any shared cell between train and test
  expect_false(any(n.cell.train %in% n.cell.test))
})

