#' @title Safe HDF5 File Operations Utilities
#' @description
#' Utility functions to ensure proper HDF5 file handle management
#' and prevent file handle leaks through RAII-style resource management.
#'
#' @keywords internal
NULL

#' @title Execute code with HDF5 file safely opened
#' @description
#' Ensures HDF5 file handle is properly closed even if errors occur.
#' Uses on.exit() to guarantee cleanup.
#'
#' @param filepath Path to HDF5 file
#' @param expr Expression to evaluate with file open
#' @param mode File open mode ("H5F_ACC_RDONLY" or "H5F_ACC_RDWR")
#'
#' @return Result of expr
#' @keywords internal
#' @importFrom rhdf5 H5Fopen H5Fclose h5closeAll
with_h5file_safe <- function(filepath, expr, mode = "H5F_ACC_RDONLY") {
  if (!file.exists(filepath)) {
    stop("HDF5 file not found: ", filepath, call. = FALSE)
  }

  h5_handle <- NULL

  tryCatch(
    {
      # Open file handle
      if (mode == "H5F_ACC_RDONLY") {
        h5_handle <- rhdf5::H5Fopen(filepath, flags = "H5F_ACC_RDONLY")
      } else {
        h5_handle <- rhdf5::H5Fopen(filepath, flags = "H5F_ACC_RDWR")
      }

      # Ensure closure on exit (even if error occurs)
      on.exit(
        {
          if (!is.null(h5_handle)) {
            rhdf5::H5Fclose(h5_handle)
          }
        },
        add = TRUE
      )

      # Evaluate expression
      force(expr)
    },
    error = function(e) {
      # Re-throw with context
      stop(
        "Error in HDF5 operation on ",
        filepath,
        ": ",
        e$message,
        call. = FALSE
      )
    }
  )
}

#' @title Check if HDF5 dataset exists safely
#' @description
#' Check dataset existence with proper handle management
#'
#' @param filepath Path to HDF5 file
#' @param dataset_path Path to dataset in HDF5 file
#'
#' @return Logical indicating if dataset exists
#' @keywords internal
#' @importFrom rhdf5 H5Fopen H5Lexists H5Fclose
h5_dataset_exists <- function(filepath, dataset_path) {
  with_h5file_safe(filepath, {
    h5_handle <- rhdf5::H5Fopen(filepath, flags = "H5F_ACC_RDONLY")
    on.exit(rhdf5::H5Fclose(h5_handle), add = TRUE)

    rhdf5::H5Lexists(h5_handle, dataset_path)
  })
}

#' @title Safe HDF5 read operation
#' @description
#' Read from HDF5 with guaranteed handle cleanup
#'
#' @param filepath Path to HDF5 file
#' @param dataset_path Path to dataset
#' @param ... Additional arguments passed to h5read
#'
#' @return Data read from HDF5
#' @keywords internal
#' @importFrom rhdf5 h5read h5closeAll
h5_read_safe <- function(filepath, dataset_path, ...) {
  result <- NULL

  tryCatch(
    {
      result <- rhdf5::h5read(file = filepath, name = dataset_path, ...)

      # Close all handles to prevent leaks
      rhdf5::h5closeAll()

      result
    },
    error = function(e) {
      rhdf5::h5closeAll() # Ensure cleanup even on error
      stop(
        "Failed to read from HDF5 dataset ",
        dataset_path,
        ": ",
        e$message,
        call. = FALSE
      )
    }
  )
}

#' @title Safe HDF5 write operation
#' @description
#' Write to HDF5 with guaranteed handle cleanup
#'
#' @param obj Object to write
#' @param filepath Path to HDF5 file
#' @param dataset_path Path to dataset
#' @param ... Additional arguments passed to h5write
#'
#' @return NULL (invisible)
#' @keywords internal
#' @importFrom rhdf5 h5write h5closeAll
h5_write_safe <- function(obj, filepath, dataset_path, ...) {
  tryCatch(
    {
      rhdf5::h5write(obj = obj, file = filepath, name = dataset_path, ...)

      # Close all handles
      rhdf5::h5closeAll()

      invisible(NULL)
    },
    error = function(e) {
      rhdf5::h5closeAll() # Ensure cleanup even on error
      stop(
        "Failed to write to HDF5 dataset ",
        dataset_path,
        ": ",
        e$message,
        call. = FALSE
      )
    }
  )
}
