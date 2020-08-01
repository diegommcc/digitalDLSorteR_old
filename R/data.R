#' Pre-trained DigitalDLSorter DNN model for deconvolution of TILs present in
#' breast cancer environment.
#'
#' DigitalDLSorter DNN model built and trained with single-cell data
#' from Chung et al., 2017 (GSE75688). This model allows the enumeration and
#' quantification of immune infiltrated cell types in breast cancer environment.
#' This dataset consists in single-cell profiles from 10 patients from different
#' tumor etiology and stages (see Torroja and Sanchez-Cabo, 2019 for more details).
#' The analysis and characterization of the cells was carried out by the authors
#' of \code{digitalDLSorteR} package.
#'
#' The cell types considered in this model are 13, four of them being different
#' molecular subtypes of breast cancer: ER+, HER2+, ER+ and HER2+, TNBC, Stromal,
#' Monocyte, Tme (memory T cells), BGC (germinal center B cells),
#' Bmem (memory B cells), DC (dendritic cells), Macrophage, TCD8 (CD8+ T cells)
#' and Treg (regulatory T cells).
#'
#' The model consists in 2 hidden layers with 200 neurons per layer trained with
#' 'kullback_leibler_divergence' loss function  batch size equal to 128 and
#' a number of epochs equal to 20.
#'
#' @format A \code{DigitalDLSorterDNN} object with the following slots:
#' \describe{
#'    \item{model}{Trained DNN model.}
#'    \item{training.history}{Evolution of metrics and loss function during training.}
#'    \item{eval.stats}{Metrics and loss results on test data.}
#'    \item{predict.results}{Predictions of cell types on test data.}
#'    \item{cell.types}{Cell types considered by DNN model.}
#'    \item{features}{Features (genes) considered by model.}
#' }
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688}
#'
#' @references
#' Chung, W., Eum, H. H., Lee, H. O., Lee, K. M., Lee, H. B., Kim,
#' K. T., et al. (2017). Single-cell RNA-seq enables comprehensive tumour and
#' immune cell profiling in primary breast cancer. Nat. Commun. 8 (1), 15081.
#' doi: \url{10.1038/ncomms15081}.
#'
#' Torroja, C. y Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep Learning algorithm to quantify
#' immune cell populations based on scRNA-Seq data. Frontiers in Genetics 10, 978. doi:
#' \url{10.3389/fgene.2019.00978}
#'
"breast.chung"




#' Subset from the original dataset used for generating \code{breast.chung}
#' model.
#'
#' 'Toy' subset from the original dataset used for generating \code{breast.chung}
#' model in order to show some examples in vignette and documentation. Data
#' are provided as a \code{SingleCellExperiment} object with counts in
#' \code{assay} slot, cells metadata in \code{colData} slot and genes metadata
#' in \code{rowData} slot.
#'
#' For more information about the complete dataset, see \code{?breast.chung}.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688}
#'
#' @references
#' Chung, W., Eum, H. H., Lee, H. O., Lee, K. M., Lee, H. B., Kim,
#' K. T., et al. (2017). Single-cell RNA-seq enables comprehensive tumour and
#' immune cell profiling in primary breast cancer. Nat. Commun. 8 (1), 15081.
#' doi: \url{10.1038/ncomms15081}.
#'
#' Torroja, C. y Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep Learning algorithm to quantify
#' immune cell populations based on scRNA-Seq data. Frontiers in Genetics 10, 978. doi:
#' \url{10.3389/fgene.2019.00978}
#'
"sc.chung.breast"


#' Colorectal Bulk RNA-Seq samples from TCGA Research Network.
#'
#' @source \url{https://www.cancer.gov/tcga}
#'
#' @references
#' Torroja, C. y Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep Learning algorithm to quantify
#' immune cell populations based on scRNA-Seq data. Frontiers in Genetics 10, 978. doi:
#' \url{10.3389/fgene.2019.00978}
#'
"TCGA.colon"
