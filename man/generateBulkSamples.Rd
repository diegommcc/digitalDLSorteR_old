% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simBulk.R
\name{generateBulkSamples}
\alias{generateBulkSamples}
\title{Generate training and test simulated bulk RNA-seq samples.}
\usage{
generateBulkSamples(
  object,
  type.data = "both",
  file.backend = NULL,
  threads = 1,
  compression.level = NULL,
  verbose = TRUE
)
}
\arguments{
\item{object}{\code{DigitalDLSorter} object with \code{single.cell.final} and
\code{prob.cell.types} slots.}

\item{type.data}{Type of data to generate among 'train', 'test' or 'both'
(the last by default).}

\item{file.backend}{Valid file path where to save the HDF5 file used as
backend. If it is equal to \code{NULL} (by default), the data are produced
and loaded in memory.}

\item{threads}{Number of threads used during the generation of bulk samples
(2 by default).}

\item{compression.level}{The compression level used if file.backend provided
(6 by default). It is an integer value between 0 (no compression) and 9
(highest and slowest compression).}

\item{verbose}{Show informative messages during the execution.}
}
\value{
A \code{\link{DigitalDLSorter}} object with \code{bulk.sim} slot
containing a list with one or two entries (depending on selected
\code{type.data} argument): 'train' and 'test'. Each entry contains a
\code{SummarizedExperiment} object with simulated bulk samples in
\code{assay} slot, sample names in \code{colData} slot and feature names in
\code{rowData} slot.
}
\description{
Generate training and test bulk profiles using the cell composition matrix
built by \code{\link{generateTrainAndTestBulkProbMatrix}} function. These
samples are generated using the assumption that the expression of gene
\eqn{i} in sample \eqn{j} is given by the sum of the cell type specific
expression \eqn{X_{ijk}} weighted by the proportions of cell type \eqn{k} in
the sample determined by the probability matrix. In practice, as described in
Torroja et al., 2019, these profiles are generated by the summation of 100
cells from different cell types determined by cell composition matrix. The
number of bulk samples is determined by dimensions of cell composition
matrix. See \code{\link{generateTrainAndTestBulkProbMatrix}} for details.
}
\details{
\code{digitalDLSorteR} allows the use of HDF5 files as back-end for the
resulting data using \code{DelayedArray} and \code{HDF5Array} packages in
cases of generating too large bulk expression matrix. This functionality
allows you to work without keeping the data loaded in memory, which will be
of vital importance during some computationally heavy steps such as neural
network training. You must provide a valid file path in \code{file.backend}
argument to store the resulting file with '.h5' extension. The data will be
accessible from R without being loaded into memory. This option slightly
slows down execution times, since subsequent transformations of data will be
carried out by chunks instead of using all data. We recommend this option due
to the large size of the simulated matrices.
}
\examples{
## loading all data in memory
DDLSSmallCompleted <- generateBulkSamples(
  DDLSSmallCompleted,
  threads = 2,
  type.data = "both"
)
\dontrun{
## using HDF5 as backend
DDLSChung <- generateBulkSamples(
  DDLSChung,
  threads = 2,
  type.data = "both",
  file.backend = "DDLSChung.bulk.sim.h5"
)
}

}
\references{
Pagès H, Hickey wcfP, Lun A (2020). DelayedArray: A unified
framework for working transparently with on-disk and in-memory array-like
datasets. R package version 0.14.1.

Pagès H (2020). HDF5Array: HDF5 backend for DelayedArray objects. R package
version 1.16.1.
}
\seealso{
\code{\link{generateTrainAndTestBulkProbMatrix}}
\code{\link{ProbMatrixCellTypes}}
}
