methylTracer <- setClass(
  Class = "methylTracer", 
  contains = "HDF5Array", 
  slots = c(marker_name = "HDF5Array", sample_name = "HDF5Array"),
  validity = function(object) {
    errors <- character()
    
    ## Check dimensions
    if (length(object@marker_name) !=
        dim(object)[1]) {
      errors <- c(errors, 
                  "marker_name length not equal X now")
    }
    
    if (length(object@sample_name) !=
        dim(object)[2]) {
      errors <- c(errors, 
                  "sample_name length not equal X col")
    }
    
    if (length(errors) ==
        0)
      TRUE else errors
  }
)

setMethod(
  f = "show", signature = "methylTracer", definition = function(object) {
    ## Display basic information about the methylTracer object
    message("methylTracer Object Summary:\n")
    message("Local_path    : ", object@seed@filepath)
    message(
      "Memory use    : ", round(
        utils::object.size(object)/1e+03,
        3
      ),
      " KB"
    )
    ## Display data dimensions
    message("Dimensions:")
    message("X   dims  : ", nrow(object), " x ", ncol(object))
    
    va <- rhdf5::h5ls(object@marker_name@seed@filepath)
    va_c <- sum(va$group == "/var")
    message("var dims  : ", nrow(object), " x ",va_c)
    
    ob <- rhdf5::h5ls(object@sample_name@seed@filepath)
    ob_c <- sum(ob$group == "/obs")
    message("obs dims  : ", ncol(object), " x ", ob_c)
  }
)


#' @title Create a methylTracer object from an HDF5 file
#'
#' @description
#' Construct a \code{methylTracer} object from an HDF5 file that stores
#' methylation data and associated metadata. The HDF5 file is expected to
#' contain a main assay dataset \code{"X"} and accompanying sample and
#' marker identifiers organised in the \code{"/obs"} and \code{"/var"}
#' groups, respectively.
#'
#' @details
#' The typical layout of an HDF5 file for \pkg{methylTracer} is:
#' \itemize{
#'   \item \code{"X"}: a 2D dataset (rows = markers/regions, cols = samples),
#'         storing methylation values as proportions in \eqn{[0, 1]}.
#'   \item \code{"/obs/<sample_name>"}: a 1D dataset containing sample IDs.
#'   \item \code{"/var/<marker_name>"}: a 1D dataset containing marker IDs,
#'         usually of the form \code{"chr_start_end"}.
#' }
#'
#' This function wraps the main assay \code{"X"} and the specified
#' \code{sample_name} and \code{marker_name} datasets into a
#' \code{methylTracer} object using \pkg{HDF5Array} for on-disk storage.
#'
#' @param h5_file Character string giving the path to the HDF5 file
#'   containing methylation data and metadata.
#' @param sample_name Character scalar. Name of the dataset inside the
#'   \code{"/obs"} group that contains sample identifiers
#'   (default: \code{"sample_name"}). The corresponding HDF5 path is
#'   \code{"/obs/<sample_name>"}.
#' @param marker_name Character scalar. Name of the dataset inside the
#'   \code{"/var"} group that contains marker identifiers
#'   (default: \code{"marker_name"}). The corresponding HDF5 path is
#'   \code{"/var/<marker_name>"}.
#'
#' @return
#' A \code{methylTracer} object with:
#' \itemize{
#'   \item the main assay backed by an \code{\link[HDF5Array]{HDF5Array}}
#'         pointing to \code{"X"} in \code{h5_file};
#'   \item \code{marker_name}: an \code{HDF5Array} vector pointing to
#'         \code{"/var/<marker_name>"};
#'   \item \code{sample_name}: an \code{HDF5Array} vector pointing to
#'         \code{"/obs/<sample_name>"}.
#' }
#' If an error occurs (e.g. missing datasets), a diagnostic message is
#' printed and \code{NULL} is returned.
#'
#' @export
#'
#' @importFrom HDF5Array HDF5Array
#' @importFrom rhdf5 h5read h5createFile h5createGroup h5write
#' @importFrom methods show new
#'
#' @examples
#' suppressMessages(library(HDF5Array))
#' toy_h5 <- system.file("extdata", "toy.h5", package = "HDF5Array")
#' X <- HDF5Array::HDF5Array(toy_h5, "M1")
#'
#' num_rows <- dim(X)[1]
#' start_values <- seq(from = 10000, by = 1000, length.out = num_rows)
#' end_values <- start_values + 1000
#' marker_name <- paste0("chr1_", start_values, "_", end_values)
#' sample_name <- paste0("sample_", seq_len(dim(X)[2]))
#' group <- c(rep("A", 50), rep("B", 50), rep("C", 50))
#' regul <- c(rep("exon", 50), rep("intron", 50), rep("TSS", 50))
#'
#' current_dir <- tempdir()
#' met_obj_test <- file.path(current_dir, "methylTracer_obj_test.h5")
#' unlink(met_obj_test)  # Remove existing file if it exists
#'
#' rhdf5::h5createFile(met_obj_test)
#' rhdf5::h5createGroup(met_obj_test, "obs")
#' rhdf5::h5createGroup(met_obj_test, "var")
#'
#' HDF5Array::writeHDF5Array(X, met_obj_test, "X")
#' rhdf5::h5write(marker_name, met_obj_test, "var/marker_name")
#' rhdf5::h5write(regul, met_obj_test, "var/regul")
#' rhdf5::h5write(sample_name, met_obj_test, "obs/sample_name")
#' rhdf5::h5write(group, met_obj_test, "obs/group")
#'
#' met_obj_test <- build_met_obj(
#'   met_obj_test,
#'   sample_name = "sample_name",
#'   marker_name = "marker_name"
#' )
#'
#' met_obj_test
build_met_obj <- function(h5_file = NULL, 
                          sample_name = "sample_name", 
                          marker_name = "marker_name") {
  ## Input validation
  if (is.null(h5_file)) {
    message("h5_file must be provided")
  }
  tryCatch(
    {
      ## Create main data matrix HDF5Array
      met <- HDF5Array::HDF5Array(h5_file, "X")
      marker_vec <- HDF5Array::HDF5Array(
        h5_file,
        paste0("/var/", marker_name)
      )
      sample_vec <- HDF5Array::HDF5Array(
        h5_file, 
        paste0("/obs/", sample_name))
      ## Create methylTracer object
      met <- new("methylTracer", 
                     met, 
                     marker_name = marker_vec, 
                     sample_name = sample_vec)
      return(met)
    }, error = function(e) {
      message("Creating methylTracer object: ", e$message)
    }
  )
}

#' @title Retrieve observations from a methylTracer object
#'
#' @description
#' Extract selected columns (and optionally selected rows) from the
#' \code{"/obs"} group of the HDF5 file underlying a
#' \code{methylTracer} object. The result is returned as a tibble-like
#' data frame with one row per cell/sample.
#'
#' @param met A \code{methylTracer} object containing methylation data.
#'   The path to the HDF5 file is taken from \code{met@seed@filepath}.
#' @param indexList Optional list of length 2 specifying which rows and
#'   columns of the \code{"/obs"} group to read. It should be of the
#'   form \code{list(row_index, col_index)}, where:
#'   \itemize{
#'     \item \code{row_index} is an integer vector of row indices
#'           (1-based). If \code{NULL}, all rows are returned.
#'     \item \code{col_index} is an integer vector of column indices
#'           (1-based) corresponding to datasets in \code{"/obs"}.
#'           If \code{NULL}, all columns are returned.
#'   }
#'   The default (\code{NULL}) returns all rows and all columns in the
#'   \code{"/obs"} group.
#'
#' @return
#' A tibble-like data frame (created via \code{dplyr::bind_cols()})
#' containing the selected observation-level metadata. Each column
#' corresponds to one dataset in \code{"/obs"} (e.g. \code{sample_name},
#' \code{group}, etc.). If no rows are selected or \code{"/obs"} is
#' empty, an empty data frame is returned and a message is printed.
#'
#' @importFrom rhdf5 h5ls h5read h5closeAll
#' @importFrom dplyr bind_cols
#' @export
#'
#' @examples
#' ## build a small HDF5 file with obs/var and create a methylTracer object
#' suppressMessages(library(HDF5Array))
#' outputDir <- tempdir()
#' toy <- system.file("extdata", "toy.h5", package = "HDF5Array")
#' X <- HDF5Array(toy, "M1")
#'
#' nRows <- dim(X)[1]
#' startN <- seq(from = 10000, by = 1000, length.out = nRows)
#' endN <- startN + 1000
#'
#' marker_name <- paste0("chr1_", startN, "_", endN)
#' sample_name <- paste0("sample_", seq_len(dim(X)[2]))
#' group <- c(rep("case", 75), rep("ctrl", 75))
#' vargroup <- c(rep("intron", 3000), rep("exon", 3000), rep("utr", 4000))
#'
#' metH5 <- file.path(outputDir, "methylTracer_obj_test.h5")
#' unlink(metH5)
#'
#' rhdf5::h5createFile(metH5)
#' rhdf5::h5createGroup(metH5, "obs")
#' rhdf5::h5createGroup(metH5, "var")
#' HDF5Array::writeHDF5Array(X, metH5, "X")
#' rhdf5::h5write(marker_name, metH5, "var/marker_name")
#' rhdf5::h5write(vargroup, metH5, "var/vargroup")
#' rhdf5::h5write(sample_name, metH5, "obs/sample_name")
#' rhdf5::h5write(group, metH5, "obs/group")
#'
#' met <- build_met_obj(
#'   metH5,
#'   sample_name = "sample_name",
#'   marker_name = "marker_name"
#' )
#'
#' ## all obs columns, all rows
#' obs_all <- metObs(met)
#'
#' ## only first 5 rows and the first 2 obs columns
#' obs_subset <- metObs(met, indexList = list(1:5, 1:2))

metObs <- function(met, indexList = NULL) {
  ## Basic checks
  if (missing(met) || is.null(met)) {
    stop("'met' must be a non-null methylTracer object.")
  }
  
  hdf5_5mc <- met@seed@filepath
  if (!file.exists(hdf5_5mc)) {
    stop("The HDF5 file referenced by 'met' does not exist: ", hdf5_5mc)
  }
  
  ## List datasets in /obs
  rhdf5::h5closeAll()
  h5_ls <- rhdf5::h5ls(hdf5_5mc)
  h5_ob <- h5_ls[h5_ls$group == "/obs", ]
  h5_ob <- h5_ob[h5_ob$dclass %in% c("FLOAT", "INTEGER", "STRING"), ]
  
  if (nrow(h5_ob) == 0L) {
    message("No datasets found under '/obs' in the HDF5 file.")
    return(dplyr::bind_cols())
  }
  
  ## Determine column indices
  all_cols <- h5_ob$name
  n_cols   <- length(all_cols)
  
  if (is.null(indexList) || length(indexList) < 2L || is.null(indexList[[2]])) {
    col_idx <- seq_len(n_cols)
  } else {
    col_idx <- as.integer(indexList[[2]])
    if (any(col_idx < 1L | col_idx > n_cols)) {
      stop("indexList[[2]] (column indices) out of range. ",
           "There are only ", n_cols, " obs datasets.")
    }
  }
  
  sel_cols <- all_cols[col_idx]
  
  ## Determine row indices
  if (is.null(indexList) || length(indexList) < 1L || is.null(indexList[[1]])) {
    row_idx <- NULL  ## read all rows
  } else {
    row_idx <- as.integer(indexList[[1]])
  }
  
  ## Read selected obs datasets
  ob_list <- vector("list", length(sel_cols))
  names(ob_list) <- sel_cols
  
  for (i in seq_along(sel_cols)) {
    dataset <- paste0("obs/", sel_cols[i])
    if (is.null(row_idx)) {
      vec <- rhdf5::h5read(hdf5_5mc, dataset)
    } else {
      vec <- rhdf5::h5read(hdf5_5mc, dataset, index = list(row_idx))
    }
    ob_list[[i]] <- vec
  }
  
  rhdf5::h5closeAll()
  
  obs <- dplyr::bind_cols(ob_list)
  
  if (nrow(obs) == 0L) {
    message("The 'obs' selection returned an empty result.")
  }
  
  obs
}


#' @title Retrieve variable-level metadata from a methylTracer object
#'
#' @description
#' Extract selected rows and columns from the \code{"/var"} group of the
#' HDF5 file underlying a \code{methylTracer} object. Typical variables
#' include marker IDs, genomic annotations (e.g. region type, gene
#' symbol), etc.
#'
#' @param met A \code{methylTracer} object containing methylation data.
#'   The path to the HDF5 file is taken from \code{met@seed@filepath}.
#' @param indexList Optional list of length 2 specifying which rows and
#'   columns of the \code{"/var"} group to read. It should be of the
#'   form \code{list(row_index, col_index)}, where:
#'   \itemize{
#'     \item \code{row_index} is an integer vector of row indices
#'           (1-based). If \code{NULL}, all rows are returned.
#'     \item \code{col_index} is an integer vector of column indices
#'           (1-based) corresponding to datasets in \code{"/var"}.
#'           If \code{NULL}, all columns are returned.
#'   }
#'   The default (\code{NULL}) returns all rows and all columns in the
#'   \code{"/var"} group.
#'
#' @return
#' A tibble-like data frame (created via \code{dplyr::bind_cols()})
#' containing variable-level metadata from the \code{"/var"} group.
#' Each column corresponds to one dataset in \code{"/var"} (e.g.
#' \code{marker_name}, \code{vargroup}, etc.). If no datasets are found
#' or the selection is empty, an empty data frame is returned and a
#' message is printed.
#'
#' @importFrom rhdf5 h5ls h5read h5closeAll
#' @importFrom dplyr bind_cols
#' @export
#'
#' @examples
#' ## input data
#' suppressMessages(library(HDF5Array))
#' outputDir <- tempdir()
#' toy <- system.file("extdata", "toy.h5", package = "HDF5Array")
#' X <- HDF5Array(toy, "M1")
#'
#' nRows <- dim(X)[1]
#' startN <- seq(from = 10000, by = 1000, length.out = nRows)
#' endN <- startN + 1000
#'
#' ## meta file
#' marker_name <- paste0("chr1_", startN, "_", endN)
#' sample_name <- paste0("sample_", seq_len(dim(X)[2]))
#' group      <- c(rep("case", 75), rep("ctrl", 75))
#' vargroup   <- c(rep("intron", 3000), rep("exon", 3000), rep("utr", 4000))
#'
#' metH5 <- file.path(outputDir, "methylTracer_obj_test.h5")
#' unlink(metH5)
#'
#' ## build metH5
#' rhdf5::h5createFile(metH5)
#' rhdf5::h5createGroup(metH5, "obs")
#' rhdf5::h5createGroup(metH5, "var")
#' HDF5Array::writeHDF5Array(X, metH5, "X")
#' rhdf5::h5write(marker_name, metH5, "var/marker_name")
#' rhdf5::h5write(vargroup,    metH5, "var/vargroup")
#' rhdf5::h5write(sample_name, metH5, "obs/sample_name")
#' rhdf5::h5write(group,       metH5, "obs/group")
#'
#' ## build methylTracer object
#' met <- build_met_obj(
#'   metH5,
#'   sample_name = "sample_name",
#'   marker_name = "marker_name"
#' )
#'
#' ## all variable metadata
#' var_all <- metVar(met)
#'
#' ## only first 5 markers and the first column in /var
#' var_subset <- metVar(met, indexList = list(1:5, 1:1))

metVar <- function(met, indexList = NULL) {
  ## Basic checks
  if (missing(met) || is.null(met)) {
    stop("'met' must be a non-null methylTracer object.")
  }
  
  hdf5_5mc <- met@seed@filepath
  if (!file.exists(hdf5_5mc)) {
    stop("The HDF5 file referenced by 'met' does not exist: ", hdf5_5mc)
  }
  
  ## List datasets in /var
  rhdf5::h5closeAll()
  h5_ls <- rhdf5::h5ls(hdf5_5mc)
  h5_va <- h5_ls[h5_ls$group == "/var", ]
  h5_va <- h5_va[h5_va$dclass %in% c("FLOAT", "INTEGER", "STRING"), ]
  
  if (nrow(h5_va) == 0L) {
    message("No datasets found under '/var' in the HDF5 file.")
    return(dplyr::bind_cols())
  }
  
  ## Determine column indices (/var datasets)
  all_cols <- h5_va$name
  n_cols   <- length(all_cols)
  
  if (is.null(indexList) || length(indexList) < 2L || is.null(indexList[[2]])) {
    col_idx <- seq_len(n_cols)
  } else {
    col_idx <- as.integer(indexList[[2]])
    if (any(col_idx < 1L | col_idx > n_cols)) {
      stop("indexList[[2]] (column indices) out of range. ",
           "There are only ", n_cols, " var datasets.")
    }
  }
  
  sel_cols <- all_cols[col_idx]
  
  ## Determine row indices
  if (is.null(indexList) || length(indexList) < 1L || is.null(indexList[[1]])) {
    row_idx <- NULL  ## read all rows
  } else {
    row_idx <- as.integer(indexList[[1]])
  }
  
  ## Read selected /var datasets
  va_list <- vector("list", length(sel_cols))
  names(va_list) <- sel_cols
  
  for (i in seq_along(sel_cols)) {
    dataset <- paste0("var/", sel_cols[i])
    if (is.null(row_idx)) {
      vec <- rhdf5::h5read(hdf5_5mc, dataset)
    } else {
      vec <- rhdf5::h5read(hdf5_5mc, dataset, index = list(row_idx))
    }
    va_list[[i]] <- vec
  }
  
  rhdf5::h5closeAll()
  
  var_df <- dplyr::bind_cols(va_list)
  
  if (nrow(var_df) == 0L) {
    message("The 'var' selection returned an empty result.")
  }
  
  var_df
}


#' @title Retrieve methylation matrix from a methylTracer object
#'
#' @description
#' Extract a subset (or the entirety) of the main methylation matrix
#' \code{"X"} from the HDF5 file underlying a \code{methylTracer}
#' object. Rows typically correspond to markers/regions and columns to
#' samples.
#'
#' @param met A \code{methylTracer} object containing methylation data.
#' @param indexList Optional list of length 2 specifying which rows and
#'   columns of \code{"X"} to read (see \code{\link{metObs}}).
#' @param as_tibble Logical; if \code{TRUE}, return a tibble with an
#'   extra column \code{marker_name} and one column per sample.
#'   If \code{FALSE} (default), return a numeric matrix.
#'
#' @return
#' If \code{as_tibble = FALSE} (default), a numeric matrix with row and
#' column names. If \code{as_tibble = TRUE}, a tibble where the first
#' column is \code{marker_name} and remaining columns are sample-wise
#' methylation values.
#'
#' @importFrom rhdf5 h5read h5closeAll
#' @importFrom dplyr bind_cols
#' @export
#'
#' @examples
#' \donttest{
#' suppressMessages(library(HDF5Array))
#' outputDir <- tempdir()
#' toy <- system.file("extdata", "toy.h5", package = "HDF5Array")
#' X <- HDF5Array(toy, "M1")  ## toy example matrix in [0, 1]
#'
#' nRows <- dim(X)[1]
#' startN <- seq(from = 10000, by = 1000, length.out = nRows)
#' endN   <- startN + 1000
#'
#' ## meta
#' marker_name <- paste0("chr1_", startN, "_", endN)
#' sample_name <- paste0("sample_", seq_len(dim(X)[2]))
#' group       <- c(rep("case", 75), rep("ctrl", 75))
#' vargroup    <- c(rep("intron", 3000), rep("exon", 3000), rep("utr", 4000))
#'
#' metH5 <- file.path(outputDir, "methylTracer_obj_test.h5")
#' unlink(metH5)
#'
#' rhdf5::h5createFile(metH5)
#' rhdf5::h5createGroup(metH5, "obs")
#' rhdf5::h5createGroup(metH5, "var")
#' HDF5Array::writeHDF5Array(X, metH5, "X")
#' rhdf5::h5write(marker_name, metH5, "var/marker_name")
#' rhdf5::h5write(vargroup,    metH5, "var/vargroup")
#' rhdf5::h5write(sample_name, metH5, "obs/sample_name")
#' rhdf5::h5write(group,       metH5, "obs/group")
#'
#' met <- build_met_obj(
#'   metH5,
#'   sample_name = "sample_name",
#'   marker_name = "marker_name"
#' )
#'
#' ## full matrix
#' X_all <- metX(met)
#'
#' ## first 10 markers, first 5 samples as tibble
#' X_sub <- metX(met, indexList = list(1:10, 1:5), as_tibble = TRUE)
#' }

metX <- function(met, indexList = NULL, as_tibble = TRUE) {
  if (missing(met) || is.null(met)) {
    stop("'met' must be a non-null methylTracer object.")
  }
  
  hdf5_5mc <- met@seed@filepath
  if (!file.exists(hdf5_5mc)) {
    stop("The HDF5 file referenced by 'met' does not exist: ", hdf5_5mc)
  }
  
  x_dim  <- dim(met)
  n_rows <- x_dim[1]
  n_cols <- x_dim[2]
  
  ## rows
  if (is.null(indexList) || length(indexList) < 1L || is.null(indexList[[1]])) {
    row_idx <- seq_len(n_rows)
  } else {
    row_idx <- as.integer(indexList[[1]])
    if (any(row_idx < 1L | row_idx > n_rows)) {
      stop("indexList[[1]] (row indices) out of range. There are ",
           n_rows, " rows in 'X'.")
    }
  }
  
  ## cols
  if (is.null(indexList) || length(indexList) < 2L || is.null(indexList[[2]])) {
    col_idx <- seq_len(n_cols)
  } else {
    col_idx <- as.integer(indexList[[2]])
    if (any(col_idx < 1L | col_idx > n_cols)) {
      stop("indexList[[2]] (column indices) out of range. There are ",
           n_cols, " columns in 'X'.")
    }
  }
  
  rhdf5::h5closeAll()
  X_sub <- rhdf5::h5read(
    file  = hdf5_5mc,
    name  = "X",
    index = list(row_idx, col_idx)
  )
  rhdf5::h5closeAll()
  
  X_mat <- as.matrix(X_sub)
  
  row_names <- as.character(met@marker_name[row_idx])
  col_names <- as.character(met@sample_name[col_idx])
  
  rownames(X_mat) <- row_names
  colnames(X_mat) <- col_names
  
  if (!as_tibble) {
    return(X_mat)
  }
  
  ## 转 tibble：第一列 marker_name，其余每列一个 sample
  df_list <- c(
    list(marker_name = row_names),
    lapply(seq_len(ncol(X_mat)), function(j) X_mat[, j])
  )
  names(df_list) <- c("marker_name", col_names)
  
  dplyr::bind_cols(df_list)
}
