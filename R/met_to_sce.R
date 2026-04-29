#' @title Convert a methylTracer object to SingleCellExperiment
#'
#' @description
#' Take a \code{methylTracer} object and convert it into a
#' \code{SingleCellExperiment} where rows are markers/windows and
#' columns are cells/samples. The main methylation matrix \code{X}
#' is used as the assay (default name \code{"meth"}), and metadata
#' stored in the HDF5 \code{/obs} and \code{/var} groups are mapped
#' to \code{colData} and \code{rowRanges}/\code{rowData}, respectively.
#'
#' Marker genomic coordinates are inferred from \code{met@marker_name},
#' assuming IDs of the form \code{"chr_start_end"}.
#'
#' @param met A \code{methylTracer} object.
#' @param assay_name Character; name of the assay to store methylation
#'   values. Default \code{"meth"}.
#' @param scale_factor Numeric or \code{NULL}. If not \code{NULL},
#'   the assay is multiplied by this factor (e.g. \code{100} to store
#'   beta values as percentages). Default \code{NULL} (no rescaling).
#' @param verbose Logical; whether to print progress messages.
#'   Default \code{TRUE}.
#'
#' @return
#' A \code{SingleCellExperiment} object with:
#' \itemize{
#'   \item rows = markers/windows from \code{met},
#'   \item columns = cells/samples from \code{met},
#'   \item one assay (by default named \code{"meth"}) containing
#'         methylation proportions in \eqn{[0, 1]} or rescaled by
#'         \code{scale_factor},
#'   \item \code{colData} populated from the HDF5 \code{/obs} group,
#'   \item \code{rowRanges} as a \code{GRanges} built from
#'         \code{marker_name}, with additional \code{/var} fields
#'         stored in \code{mcols(rowRanges)}.
#' }
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom HDF5Array HDF5Array
#' @importFrom stats setNames
#' @export
met_to_sce <- function(
    met,
    assay_name   = "meth",
    scale_factor = NULL,
    verbose      = TRUE
) {
  ## 1. Sanity checks --------------------------------------------------------
  if (is.null(met)) {
    stop("'met' must be a non-null methylTracer object.")
  }
  if (!methods::is(met, "methylTracer")) {
    stop("'met' must be a methylTracer object.")
  }
  
  h5_file <- met@seed@filepath
  if (!file.exists(h5_file)) {
    stop("The HDF5 file referenced by 'met' does not exist: ", h5_file)
  }
  
  if (verbose) {
    message("Converting methylTracer to SingleCellExperiment...")
    message("HDF5 file: ", h5_file)
  }
  
  ## 2. Extract main matrix X (HDF5-backed) ---------------------------------
  X <- HDF5Array::HDF5Array(h5_file, name = "X")
  
  n_markers <- nrow(X)
  n_cells   <- ncol(X)
  
  if (!is.null(scale_factor)) {
    if (!is.numeric(scale_factor) || length(scale_factor) != 1L) {
      stop("'scale_factor' must be a numeric scalar.")
    }
    if (verbose) {
      message("Applying scale_factor = ", scale_factor, " to assay.")
    }
    X <- X * scale_factor
  }
  
  ## 3. Sample names and marker names ---------------------------------------
  sample_name <- as.character(met@sample_name)
  marker_name <- as.character(met@marker_name)
  
  if (length(sample_name) != n_cells) {
    stop("Length of 'met@sample_name' (", length(sample_name),
         ") does not match number of columns in X (", n_cells, ").")
  }
  if (length(marker_name) != n_markers) {
    stop("Length of 'met@marker_name' (", length(marker_name),
         ") does not match number of rows in X (", n_markers, ").")
  }
  
  colnames(X) <- sample_name
  rownames(X) <- marker_name
  
  ## 4. Build colData from /obs ---------------------------------------------
  if (verbose) {
    message("Reading /obs metadata into colData...")
  }
  obs_df <- metObs(met)  ## tibble/data.frame
  col_df <- S4Vectors::DataFrame(obs_df)
  
  ## Ensure rownames(colData) == sample_name
  if (!"sample_name" %in% colnames(col_df)) {
    col_df$sample_name <- sample_name
  }
  rownames(col_df) <- as.character(col_df$sample_name)
  
  ## 5. Build rowRanges and rowData from /var -------------------------------
  if (verbose) {
    message("Reading /var metadata and constructing rowRanges...")
  }
  var_df <- metVar(met)  ## tibble/data.frame
  var_df <- S4Vectors::DataFrame(var_df)
  
  if (!"marker_name" %in% colnames(var_df)) {
    var_df$marker_name <- marker_name
  }
  
  ## Construct GRanges from marker_name "chr_start_end"
  sp <- strsplit(as.character(var_df$marker_name), "_", fixed = TRUE)
  if (length(sp) != n_markers) {
    stop("Failed to parse 'marker_name' into genomic coordinates; ",
         "expected one ID per marker.")
  }
  
  chr   <- vapply(sp, `[`, character(1L), 1L)
  start <- as.integer(vapply(sp, `[`, character(1L), 2L))
  end   <- as.integer(vapply(sp, `[`, character(1L), 3L))
  
  marker_gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges   = IRanges::IRanges(start = start, end = end)
  )
  
  ## Attach additional /var metadata as mcols
  mcols(marker_gr) <- var_df
  
  ## 6. Construct SingleCellExperiment --------------------------------------
  if (verbose) {
    message("Creating SingleCellExperiment object...")
  }
  
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays    = stats::setNames(list(X), assay_name),
    rowRanges = marker_gr,
    colData   = col_df
  )
  
  if (verbose) {
    message("Done. Created SingleCellExperiment with ",
            nrow(sce), " markers x ", ncol(sce), " cells.")
  }
  
  sce
}
