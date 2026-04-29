#' Convert a SingleCellExperiment back to a methylTracer object
#'
#' @description
#' Take a \code{SingleCellExperiment} with DMR/marker-level methylation
#' per cell and write it to an HDF5 file in the layout expected by
#' \pkg{methylTracer}, then wrap it as a \code{methylTracer} object via
#' \code{\link{build_met_obj}}.
#'
#' This is essentially the inverse of \code{dmrs_to_sce_memory()}:
#' it uses one assay (default \code{"meth"}) as the main \code{"X"}
#' matrix and stores column/row metadata in the \code{"/obs"} and
#' \code{"/var"} groups, respectively.
#'
#' @param sce A \code{SingleCellExperiment} object. Rows are markers or
#'   DMRs, columns are cells/samples. One assay must contain the
#'   methylation proportions in \eqn{[0, 1]}.
#' @param h5_file Character scalar. Path to the HDF5 file to create.
#'   If the file exists and \code{overwrite = FALSE}, an error is thrown.
#' @param assay_name Character scalar giving the assay in \code{sce} to
#'   use as the main methylation matrix. Default is \code{"meth"}.
#' @param sample_name_col Optional character scalar naming the column
#'   in \code{colData(sce)} to use as \code{sample_name}. If \code{NULL}
#'   (default), the function will use \code{"sample_name"} if present,
#'   otherwise it falls back to \code{colnames(sce)}.
#' @param marker_name_col Optional character scalar naming the column
#'   in \code{rowData(sce)} to use as \code{marker_name}. If \code{NULL}
#'   (default), the function will use \code{"marker_name"} if present;
#'   otherwise, if \code{rowRanges(sce)} is a \code{GRanges}, it will
#'   construct IDs of the form \code{"chr_start_end"}. As a last resort,
#'   it falls back to \code{rownames(sce)}.
#' @param overwrite Logical; whether to overwrite an existing HDF5 file.
#'   Default \code{FALSE}.
#'
#' @return A \code{methylTracer} object pointing to the newly created
#'   HDF5 file.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assay colData rowData rowRanges
#' @importFrom HDF5Array writeHDF5Array HDF5Array
#' @importFrom rhdf5 h5createFile h5createGroup h5write h5ls h5closeAll
#' @export
sce_to_met <- function(
    sce,
    h5_file,
    assay_name      = "meth",
    sample_name_col = NULL,
    marker_name_col = NULL,
    overwrite       = FALSE
) {
  ## basic checks
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("'sce' must be a SingleCellExperiment.")
  }
  if (missing(h5_file) || is.null(h5_file)) {
    stop("'h5_file' must be a non-empty file path.")
  }
  if (file.exists(h5_file)) {
    if (!isTRUE(overwrite)) {
      stop("File already exists: ", h5_file,
           " (set overwrite = TRUE to replace it).")
    } else {
      unlink(h5_file)
    }
  }
  
  if (!assay_name %in% SummarizedExperiment::assayNames(sce)) {
    stop("Assay '", assay_name, "' not found in 'sce'.")
  }
  
  ## extract matrix
  X <- SummarizedExperiment::assay(sce, assay_name)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  
  n_markers <- nrow(X)
  n_cells   <- ncol(X)
  
  ## sample names / colData
  cd <- SummarizedExperiment::colData(sce)
  
  # decide/sample_name vector
  if (is.null(sample_name_col)) {
    if ("sample_name" %in% colnames(cd)) {
      sample_name_col <- "sample_name"
    }
  }
  
  if (!is.null(sample_name_col)) {
    if (!sample_name_col %in% colnames(cd)) {
      stop("Specified 'sample_name_col' not found in colData(sce).")
    }
    sample_name <- as.character(cd[[sample_name_col]])
  } else {
    sample_name <- colnames(sce)
    if (is.null(sample_name)) {
      sample_name <- paste0("cell_", seq_len(n_cells))
    }
    # attach to colData for consistency
    cd$sample_name <- sample_name
  }
  
  # enforce colData rownames = sample_name
  rownames(cd) <- sample_name
  
  ## marker names / rowData
  rd <- SummarizedExperiment::rowData(sce)
  
  if (is.null(marker_name_col)) {
    if ("marker_name" %in% colnames(rd)) {
      marker_name_col <- "marker_name"
    }
  }
  
  if (!is.null(marker_name_col) && marker_name_col %in% colnames(rd)) {
    marker_name <- as.character(rd[[marker_name_col]])
  } else {
    rr <- SummarizedExperiment::rowRanges(sce)
    if (inherits(rr, "GenomicRanges") || inherits(rr, "GRanges")) {
      marker_name <- paste0(
        as.character(GenomeInfoDb::seqnames(rr)), "_",
        BiocGenerics::start(rr), "_",
        BiocGenerics::end(rr)
      )
    } else if (!is.null(rownames(sce))) {
      marker_name <- rownames(sce)
    } else {
      marker_name <- paste0("marker_", seq_len(n_markers))
    }
    rd$marker_name <- marker_name
  }
  
  ## create HDF5 skeleton
  rhdf5::h5createFile(h5_file)
  rhdf5::h5createGroup(h5_file, "obs")
  rhdf5::h5createGroup(h5_file, "var")
  
  ## write main matrix X
  HDF5Array::writeHDF5Array(X, h5_file, "X")
  
  ## write /obs (colData)
  if (ncol(cd) > 0L) {
    for (nm in colnames(cd)) {
      vec <- cd[[nm]]
      # ensure 1D atomic vector
      vec <- as.vector(vec)
      if (length(vec) != n_cells) {
        stop("Column '", nm,
             "' in colData does not have length equal to #cells.")
      }
      rhdf5::h5write(vec, h5_file, file.path("obs", nm))
    }
  } else {
    # sample_name
    rhdf5::h5write(sample_name, h5_file, "obs/sample_name")
  }
  
  ## write /var (rowData)
  if (ncol(rd) > 0L) {
    for (nm in colnames(rd)) {
      vec <- rd[[nm]]
      vec <- as.vector(vec)
      if (length(vec) != n_markers) {
        stop("Column '", nm,
             "' in rowData does not have length equal to #markers.")
      }
      rhdf5::h5write(vec, h5_file, file.path("var", nm))
    }
  } else {
    rhdf5::h5write(marker_name, h5_file, "var/marker_name")
  }
  
  rhdf5::h5closeAll()
  
  ## wrap as methylTracer
  met <- build_met_obj(
    h5_file     = h5_file,
    sample_name = "sample_name",
    marker_name = "marker_name"
  )
  
  met
}
