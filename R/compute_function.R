#' @title Compute Column Means for Observation Groups
#' @description This function computes the mean of selected columns from 
#'    an HDF5 dataset based on group indices specified in the observation 
#'    column (`obsColName`). It retrieves data from an HDF5 file, calculates
#'    the group-wise means, normalizes the values by dividing by 1000, and 
#'    writes the result to a new dataset in the same HDF5 file. The result 
#'    is saved under a variable name combining the dataset and group names.
#'
#' @param met A `met_obj` object containing metadata. This object must include
#'    information about the file path to the HDF5 file where the data is 
#'    stored.
#' @param dataset A string specifying the dataset path inside the HDF5 file.
#'     It defaults to '/X', but could be customized to other datasets 
#'    (e.g., '/uns/coverage').
#' @param obsColName A string that defines the observation column group 
#'    name within the HDF5 file. This is used to identify the groups (or 
#'    classes) within the dataset.
#'
#' @return This function does not return any value but writes the calculated 
#'    mean values for each group into a new dataset in the HDF5 file.
#' @export
#' @importFrom HDF5Array HDF5Array
#' @importFrom stringi stri_join
#' @importFrom DelayedArray unique which
#' @importFrom DelayedMatrixStats rowMeans2
#' @importFrom rhdf5 h5write
#' @examples
#' ## input data
#' suppressMessages(library(HDF5Array))
#' outputDir <- tempdir()
#' toy <- system.file('extdata', 'toy.h5', package = 'HDF5Array')
#' X <- HDF5Array(toy, 'M1')
#' X <- round(X * 1000, 0)
#' nRows <- dim(X)[1]
#' startN <- seq(from = 10000, by = 1000, length.out = nRows)
#' endN <- startN + 1000
#'
#' ## meta file
#' marker_name <- paste0('chr1_', startN, '_', endN)
#' sample_name <- paste0('sample_', 1:dim(X)[2])
#' group <- c(rep('case', 75), rep('ctrl', 75))
#' vargroup <- c(rep('intron', 3000), rep('exon', 3000), rep('utr', 4000))
#' metH5 <- file.path(outputDir, 'methylTracer_obj_test.h5')
#' unlink(metH5)
#'
#' ## build metH5
#' rhdf5::h5createFile(metH5)
#' rhdf5::h5createGroup(metH5, 'obs')
#' rhdf5::h5createGroup(metH5, 'var')
#' HDF5Array::writeHDF5Array(X, metH5, 'X')
#' rhdf5::h5write(marker_name, metH5, 'var/marker_name')
#' rhdf5::h5write(vargroup, metH5, 'var/vargroup')
#' rhdf5::h5write(sample_name, metH5, 'obs/sample_name')
#' rhdf5::h5write(group, metH5, 'obs/group')
#'
#' ## build methylTracer object
#' met <- build_met_obj(metH5,
#'     sample_name = 'sample_name',
#'     marker_name = 'marker_name'
#' )
#' 
#' ## compute obs group values
#' computeObs(met = met, dataset = '/X', obsColName = 'group')
#' 

computeObs <- function(met, dataset = "/X", obsColName = NULL)
{
  h5_f <- met@seed@filepath
  ## get index of each group
  h5_gn <- HDF5Array::HDF5Array(h5_f, stringi::stri_join("/obs", obsColName, sep = "/"))
  g_ns <- DelayedArray::unique(h5_gn)
  g_i <- list()
  for (i in g_ns)
  {
    idx <- DelayedArray::which(h5_gn == i)
    g_i[[as.character(i)]] <- idx
  }
  ## comupte
  h5_x <- HDF5Array::HDF5Array(h5_f, dataset)
  for (k in seq_along(g_i))
  {
    indices <- g_i[[k]]
    val <- DelayedMatrixStats::rowMeans2(h5_x, cols = indices, na.rm = TRUE)
    # val <- val/1000
    
    nam <- names(g_i)[k]
    d_s <- stringi::stri_join(
      "/var/mean", gsub(".*/", "", dataset),
      nam, sep = "_"
    )
    rhdf5::h5write(val, h5_f, d_s)
    message(stringi::stri_join("Saving results in", d_s, sep = " "))
  }
}


#' @title Compute Column Means for Variable Groups
#'
#' @description This function computes the mean of selected columns 
#'     (variables) from an HDF5 dataset based on specified group indices 
#'     in the regulation column (`var_group`). It retrieves data from an 
#'     HDF5 file, calculates the column-wise means, and writes the result to 
#'     a new dataset in the same HDF5 file. The result is saved under a 
#'     variable name combining the dataset and regulation group names.
#'
#' @param met A `methylTracer` object containing metadata.
#' @param dataset A string specifying the dataset path inside the HDF5 file. 
#'     It defaults to "/X", but can be customized for other datasets such 
#'     as "/coverage".
#' @param varColName A string representing the regulation group within 
#' the HDF5 file. This is used to identify groups (or categories) within 
#' the dataset for computing the column-wise means.
#'
#' @return This function does not return any value but writes the computed mean 
#' values for each group into a new dataset in the HDF5 file. 
#' @export
#' @importFrom HDF5Array HDF5Array
#' @importFrom stringi stri_join
#' @importFrom DelayedArray unique which
#' @importFrom DelayedMatrixStats colMeans2
#' @importFrom rhdf5 h5write
#'
#' @examples
#' ## input data
#' suppressMessages(library(HDF5Array))
#' outputDir <- tempdir()
#' toy <- system.file('extdata', 'toy.h5', package = 'HDF5Array')
#' X <- HDF5Array(toy, 'M1')
#' X <- round(X * 1000, 0)
#' nRows <- dim(X)[1]
#' startN <- seq(from = 10000, by = 1000, length.out = nRows)
#' endN <- startN + 1000
#'
#' ## meta file
#' marker_name <- paste0('chr1_', startN, '_', endN)
#' sample_name <- paste0('sample_', 1:dim(X)[2])
#' group <- c(rep('case', 75), rep('ctrl', 75))
#' vargroup <- c(rep('intron', 3000), rep('exon', 3000), rep('utr', 4000))
#' metH5 <- file.path(outputDir, 'methylTracer_obj_test.h5')
#' unlink(metH5)
#'
#' ## build metH5
#' rhdf5::h5createFile(metH5)
#' rhdf5::h5createGroup(metH5, 'obs')
#' rhdf5::h5createGroup(metH5, 'var')
#' HDF5Array::writeHDF5Array(X, metH5, 'X')
#' rhdf5::h5write(marker_name, metH5, 'var/marker_name')
#' rhdf5::h5write(vargroup, metH5, 'var/vargroup')
#' rhdf5::h5write(sample_name, metH5, 'obs/sample_name')
#' rhdf5::h5write(group, metH5, 'obs/group')
#'
#' ## build methylTracer object
#' met <- build_met_obj(metH5,
#'     sample_name = 'sample_name',
#'     marker_name = 'marker_name'
#' )
#' 
#' 
#' 
#' ## compute var group values
#' computeVar(met = met, dataset = '/X', varColName = 'vargroup')
#' 
computeVar <- function(met, dataset = "/X", varColName = NULL)
{
  h5_f <- met@seed@filepath
  ## get index of each group
  h5_gn <- HDF5Array::HDF5Array(h5_f, stringi::stri_join("/var", varColName, sep = "/"))
  g_ns <- DelayedArray::unique(h5_gn)
  g_i <- list()
  for (i in g_ns)
  {
    idx <- DelayedArray::which(h5_gn == i)
    g_i[[as.character(i)]] <- idx
  }
  ## comupte mean
  h5_x <- HDF5Array::HDF5Array(h5_f, dataset)
  for (k in seq_along(g_i))
  {
    indices <- g_i[[k]]
    val <- DelayedMatrixStats::colMeans2(h5_x, rows = indices, na.rm = TRUE)
    # val <- val/1000
    nam <- names(g_i)[k]
    d_s <- stringi::stri_join(
      "/obs/mean", gsub(".*/", "", dataset),
      nam, sep = "_")
    rhdf5::h5write(val, h5_f, d_s)
    message(stringi::stri_join("Saving results in", d_s, sep = " "))
  }
}
