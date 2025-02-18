methylTracer <- setClass(
    Class = "methylTracer", 
    contains = "HDF5Array", 
    slots = c(marker_name = "character", sample_name = "character"),
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
        # Display basic information about the methylTracer object
        message("methylTracer Object Summary:\n")
        message("Local_path    : ", seed(object)@filepath)
        message(
            "Memory use    : ", round(
                  utils::object.size(object)/1e+06,
                  2
              ),
                " MB"
        )
        ## Display data dimensions
        message("Dimensions:")
        message(
            "Main matrix   : ", nrow(object),
                " x ", ncol(object)
        )
        message(
            "Marker names  : ", format(
                  nrow(object),
                  big.mark = ","
              ),
                " entries"
        )
        message(
            "Sample names  : ", format(
                  ncol(object),
                  big.mark = ","
              ),
                " entries"
        )
    }
)


#' @title Create a methylTracer object from HDF5 file
#' @description 
#' This function reads methylation data from an HDF5 file and constructs 
#' a methylTracer object containing the methylation data and associated 
#' metadata. The data should be structured with the sample identifiers 
#' stored in the /obs group and marker identifiers in the /var group.
#'
#' @param h5_file Path to the HDF5 file containing methylation data.
#'                This file should have datasets for both sample and marker
#'                identifiers, as well as the methylation data.
#' @param sample_name Name of the dataset containing sample identifiers 
#'                    in the /obs group. This should match the dataset 
#'                    name where the sample identifiers are stored.
#' @param marker_name Name of the dataset containing marker identifiers 
#'                    in the /var group. This should match the dataset 
#'                    name where the marker identifiers are stored.
#'
#' @return A methylTracer object containing the methylation data and metadata. 
#'         This object will be structured with the sample and marker data 
#'         from the HDF5 file.
#' @export
#' @importFrom HDF5Array HDF5Array
#' @importFrom rhdf5 h5read h5createFile h5createGroup h5write
#' @importFrom methods show new
#'
#' @examples
#' toy_h5 <- system.file('extdata', 'toy.h5', package = 'HDF5Array')
#' X <- HDF5Array(toy_h5, 'M1') # Load data matrix from HDF5 file
#'
#' num_rows <- dim(X)[1]
#' start_values <- seq(from = 10000, by = 1000, length.out = num_rows)
#' end_values <- start_values + 1000
#' marker_name <- paste0('chr1_', start_values, '_', end_values)
#' sample_name <- paste0('sample_', 1:dim(X)[2])
#'
#' current_dir <- tempdir()
#' met_obj_test <- paste0(current_dir, '/methylTracer_obj_test.h5')
#' unlink(met_obj_test) # Remove existing file if it exists
#'
#' rhdf5::h5createFile(met_obj_test)
#' rhdf5::h5createGroup(met_obj_test, 'obs')
#' rhdf5::h5createGroup(met_obj_test, 'var')
#'
#' HDF5Array::writeHDF5Array(X, met_obj_test, 'X')
#' rhdf5::h5write(marker_name, met_obj_test, 'var/marker_name')
#' rhdf5::h5write(sample_name, met_obj_test, 'obs/sample_name')
#'
#' met_obj_test <- build_met_obj(met_obj_test,
#'   sample_name = 'sample_name',
#'   marker_name = 'marker_name'
#' )
#'
#' met_obj_test
#'
#' head(met_obs(met_obj_test))
#'
#' head(met_var(met_obj_test))

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
            met_obj <- HDF5Array(h5_file, "X")

            ## Read marker and sample names as character vectors
            marker_vec <- as.character(
              rhdf5::h5read(h5_file, 
                            paste0("/var/", marker_name)))
            sample_vec <- as.character(
              rhdf5::h5read(h5_file, 
                            paste0("/obs/", sample_name)))

            ## Create methylTracer object
            met_obj <- new("methylTracer", 
                           met_obj, 
                           marker_name = marker_vec, 
                           sample_name = sample_vec)

            return(met_obj)
        }, error = function(e) {
            message("Creating methylTracer object: ", e$message)
        }
    )
}



#' @title Retrieve Observations from methylTracer Object
#'
#' @description This function retrieves specified rows or the entire
#' 'obs' dataset from a `methylTracer` object stored in an HDF5 file.
#' The dataset is read using the `rhdf5::h5read` function, and the
#' result is returned as a `data.frame`.
#'
#' @param met_obj A `methylTracer` object containing methylation
#' data, which has the file path to the HDF5
#' file stored in its `@seed@filepath` slot.
#'
#' @param index Optional; a numeric vector specifying which rows to retrieve
#' from the 'obs' dataset in the HDF5 file. If `NULL` (default),
#' the entire 'obs' dataset is retrieved.
#'
#' @return A `data.frame` containing the observations from the
#' 'obs' dataset. If `index` is specified, only the selected rows 
#' will be returned. If no rows are found or the dataset is empty,
#' an empty `data.frame` is returned.
#'
#' @import rhdf5
#' @export
#' @examples
#' toy_h5 <- system.file('extdata', 'toy.h5', package = 'HDF5Array')
#' X <- HDF5Array(toy_h5, 'M1') # Load data matrix from HDF5 file
#'
#' num_rows <- dim(X)[1]
#' start_values <- seq(from = 10000, by = 1000, length.out = num_rows)
#' end_values <- start_values + 1000
#' marker_name <- paste0('chr1_', start_values, '_', end_values)
#' sample_name <- paste0('sample_', 1:dim(X)[2])
#'
#' current_dir <- tempdir()
#' met_obj_test <- paste0(current_dir, '/methylTracer_obj_test.h5')
#' unlink(met_obj_test) # Remove existing file if it exists
#'
#' rhdf5::h5createFile(met_obj_test)
#' rhdf5::h5createGroup(met_obj_test, 'obs')
#' rhdf5::h5createGroup(met_obj_test, 'var')
#'
#' HDF5Array::writeHDF5Array(X, met_obj_test, 'X')
#' rhdf5::h5write(marker_name, met_obj_test, 'var/marker_name')
#' rhdf5::h5write(sample_name, met_obj_test, 'obs/sample_name')
#'
#' met_obj_test <- build_met_obj(met_obj_test,
#'   sample_name = 'sample_name',
#'   marker_name = 'marker_name'
#' )
#'
#' head(met_obs(met_obj_test))

met_obs <- function(met_obj = NULL, index = NULL) {
    ## Check if met_obj is provided and is of the correct class
    if (is.null(met_obj)) {
        message("met_obj must be a non-null object .")
    }

    ## Retrieve the file path from the met_obj
    hdf5_5mc <- met_obj@seed@filepath

    ## Check if the HDF5 file exists
    if (!file.exists(hdf5_5mc)) {
        message("The specified HDF5 file does not exist.")
    }

    ## Read the 'obs' dataset from the HDF5 file
    obs <- as.data.frame(h5read(hdf5_5mc, "obs", index = index))

    ## Check if the returned data frame is empty
    if (nrow(obs) ==
        0) {
        message("The 'obs' dataset is empty.")
    }

    return(obs)
}

#' @title Retrieve Variables from methylTracer Object
#'
#' @description 
#' This function retrieves specified rows (default is the first 10) 
#' from the 'var' datasets within a `methylTracer` object stored in 
#' an HDF5 file. It reads datasets under the '/var' group and returns 
#' the extracted rows from each dataset as a `data.frame`.
#'
#' @param met_obj A `methylTracer` object containing methylation data, 
#'                which has the file path to the HDF5 file.
#'
#' @return A `data.frame` containing the first few rows (default: 10) 
#'         of each dataset found in the 'var' group of the HDF5 file. 
#'         If no data is retrieved, an empty `data.frame` is returned.
#'
#' @import rhdf5
#' @export
#'
#' @examples
#' toy_h5 <- system.file('extdata', 'toy.h5', package = 'HDF5Array')
#' X <- HDF5Array(toy_h5, 'M1') # Load data matrix from HDF5 file
#'
#' num_rows <- dim(X)[1]
#' start_values <- seq(from = 10000, by = 1000, length.out = num_rows)
#' end_values <- start_values + 1000
#' marker_name <- paste0('chr1_', start_values, '_', end_values)
#' sample_name <- paste0('sample_', 1:dim(X)[2])
#'
#' current_dir <- tempdir()
#' met_obj_test <- paste0(current_dir, '/methylTracer_obj_test.h5')
#' unlink(met_obj_test) # Remove existing file if it exists
#'
#' rhdf5::h5createFile(met_obj_test)
#' rhdf5::h5createGroup(met_obj_test, 'obs')
#' rhdf5::h5createGroup(met_obj_test, 'var')
#'
#' HDF5Array::writeHDF5Array(X, met_obj_test, 'X')
#' rhdf5::h5write(marker_name, met_obj_test, 'var/marker_name')
#' rhdf5::h5write(sample_name, met_obj_test, 'obs/sample_name')
#'
#' met_obj_test <- build_met_obj(
#'   met_obj_test,
#'   sample_name = 'sample_name',
#'   marker_name = 'marker_name'
#' )
#'
#' head(met_var(met_obj_test))

met_var <- function(met_obj = NULL) {
    ## Check if met_obj is provided and is of the correct class
    if (is.null(met_obj)) {
        message("met_obj must be a non-null object.")
    }

    ## Retrieve the file path from the met_obj
    hdf5_5mc <- met_obj@seed@filepath

    ## Check if the HDF5 file exists
    if (!file.exists(hdf5_5mc)) {
        message("The specified HDF5 file does not exist.")
    }
    ## Convert the list to a data frame
    var_df <- as.data.frame(rhdf5::h5read(met_obj@seed@filepath, "/var"))
    ## Check if the resulting data frame is empty
    if (nrow(var_df) ==
        0) {
        message("No data was retrieved from the 'var' datasets.")
    }
    return(var_df)
}
