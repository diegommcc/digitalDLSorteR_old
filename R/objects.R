# Test file: I want to test if I can add S4 classes as slots in the main class.
# I can do it. Therefore, I can use SCE for the single-cell profiles and even
# SummaryExperiment class for Bulk profiles.

## dependencies: these packages will be introduced in NAMESPACE file
library(SingleCellExperiment)
library(splatter)

setOldClass(Classes = 'package_version')
setClassUnion("SingleCellExperimentOrNULL", c("SingleCellExperiment", "NULL"))
setClassUnion("ZINBParamsOrNULL", c("ZINBParams", "NULL"))

## devuelve un warning porque hay que definir la clase package_version
DigitalDLSorter <- setClass(
  Class = "DigitalDLSorter",
  slots = c(
    single.cell.real = "SingleCellExperimentOrNULL",
    zinb.params = "ZINBParamsOrNULL",
    selected.genes = "character",
    single.cell.sim = "SingleCellExperimentOrNULL",
    project = "character",
    version = "package_version"
  )
)




