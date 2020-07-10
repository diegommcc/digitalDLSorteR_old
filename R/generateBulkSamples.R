## generateBulkSamples: dar la posibilidad de utilizar un backend H5 o hacerlo normal


# type.data == "test", "train", "both"


generateBulkSamples <- function(object, threads = 2, type.data) {
  if (class(object) != "DigitalDLSorter") {
    stop("The object provided is not of DigitalDLSorter class")
  } else if (is.null(single.cell.sim(object))) {
    stop("single.cell.sim slot is empty")
  } else if (is.null(prob.matrix(object))) {
    stop("prob.matrix slot is empty")
  }


  simCounts <- edgeR::cpm.default(simCounts)

  ### generate bulk counts
  setBulks <- function (x,c,i) { # para qué vale i
    return(rowSums(c[,x]))
  }

  message("GENERATE BULK COUNTS\n")
  message(paste("Cores",nCores,"\n"))
  bulkCounts <- pbapply(probMatrixNames, 1, FUN = setBulks,
                        c = simCounts, cl = nCores)

  cat(paste0(paste(c("Genes","Samples"),dim(bulkCounts))),"\n")
  cat("DONE\n")
  cat("\n")

    colnames(bulkCounts) <- paste("Bulk",seq(dim(bulkCounts)[2]),sep = "_")
}
