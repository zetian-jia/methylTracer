#' @title Create HDF5 file for methylTracer input
#'
#' @description
#' Create an HDF5 file that stores region-level methylation proportions
#' (0–1), coverage, sample metadata, and marker annotation in a layout
#' compatible with \pkg{methylTracer}.
#'
#' @details
#' The main methylation data matrix is written as a dense double-precision
#' dataset named \code{"X"} at the root of the HDF5 file. Each row
#' corresponds to a genomic marker (region), and each column corresponds
#' to a sample. Values are methylation proportions in the range
#' \eqn{[0, 1]}.
#'
#' Coverage data, if provided, are stored under the dataset
#' \code{"/uns/coverage"} with the same shape as \code{"X"}.
#'
#' Sample metadata are stored under the \code{"/obs"} group, with one
#' dataset per column in \code{sam_info}. Marker annotation columns are
#' stored similarly under the \code{"/var"} group, one dataset per
#' annotation column.
#'
#' @param sam_info Character string. Path to a CSV file containing
#'   metadata for each sample. The file must contain a column named
#'   \code{"sample_name"}, whose values correspond to the column names
#'   of \code{input_file} (excluding the first marker column).
#' @param input_file Character string. Path to the methylation matrix
#'   (typically produced by \code{\link{build_count_matrix}}). The first
#'   column must be \code{marker_name}, and the remaining columns are
#'   sample-wise methylation proportions (0–1).
#' @param input_coverage Optional character string. Path to the coverage
#'   matrix with the same structure as \code{input_file} (first column
#'   \code{marker_name}, remaining columns are coverage values). If
#'   \code{NULL} (default), coverage is not written.
#' @param output_dir Character string. Directory where the HDF5 file
#'   will be saved. Created if it does not exist.
#' @param output_file Character string. File name of the HDF5 file to
#'   create inside \code{output_dir}.
#' @param annotation_file Character string. Path to the genomic annotation
#'   file (e.g. \code{"3.ann.bed"}) produced by
#'   \code{\link{build_count_matrix}}. The file is expected to have a
#'   header and one row per marker.
#' @param chunk_size Integer. Number of rows to read from the input text
#'   file per chunk when streaming into HDF5. Default is \code{1e6}.
#'   Smaller values reduce memory usage but may be slower.
#' @param level Integer. Compression level for the \code{"X"} and
#'   \code{"/uns/coverage"} datasets. A value between 0 (no compression)
#'   and 9 (maximum compression). Default is \code{1}.
#' @param show_progress Logical. Whether to display a text progress bar
#'   while streaming data into the HDF5 file. Default is
#'   \code{getOption("methylTracer.show_progress", TRUE)}.
#'
#' @return
#' Invisibly returns the path to the created HDF5 file. The function is
#' called for its side effects; it writes the following groups/datasets:
#' \itemize{
#'   \item \code{"X"}: dense double matrix of methylation proportions
#'         (rows = markers, cols = samples).
#'   \item \code{"/uns/coverage"}: dense double matrix of coverage
#'         (optional).
#'   \item \code{"/obs/"}: one dataset per column of \code{sam_info}.
#'   \item \code{"/var/"}: one dataset per column of \code{annotation_file}.
#' }
#'
#' @importFrom data.table fread
#' @importFrom rhdf5 h5createFile h5createGroup h5createDataset
#'   h5set_extent h5write H5Sunlimited h5closeAll
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#'
#' @examples
#' ## Minimal toy example (0–1 methylation proportions)
#' output_dir <- tempdir()
#'
#' sam_info_df <- data.frame(
#'   sample_name = c("sample_1", "sample_2"),
#'   group       = c("case", "ctrl"),
#'   stringsAsFactors = FALSE
#' )
#' sam_info <- file.path(output_dir, "sample_info.csv")
#'
#' input_file_df <- data.frame(
#'   marker_name = c("chr1_1000_2000", "chr1_2000_3000", "chr1_3000_4000"),
#'   sample_1    = c(0.10, 0.50, 0.90),
#'   sample_2    = c(0.20, 0.40, 0.80)
#' )
#' input_file <- file.path(output_dir, "methylTracer_1kb.txt")
#'
#' annotation_file_df <- data.frame(
#'   chr         = rep("chr1", 3),
#'   start       = c(1000L, 2000L, 3000L),
#'   end         = c(2000L, 3000L, 4000L),
#'   SYMBOL      = c("gene1", "gene2", "gene3"),
#'   marker_name = c("chr1_1000_2000", "chr1_2000_3000", "chr1_3000_4000"),
#'   stringsAsFactors = FALSE
#' )
#' annotation_file <- file.path(output_dir, "annotation.bed")
#'
#' output_file <- "methylTracer_obj_test.h5"
#'
#' ## Write input files
#' write.csv(sam_info_df, sam_info, row.names = FALSE)
#' write.table(input_file_df, input_file, sep = "\t", row.names = FALSE,
#'             quote = FALSE)
#' write.table(annotation_file_df, annotation_file, sep = "\t",
#'             row.names = FALSE, quote = FALSE)
#'
#' ## Remove existing output if present
#' unlink(file.path(output_dir, output_file), recursive = TRUE)
#'
#' ## Run the function
#' build_h5(
#'   sam_info       = sam_info,
#'   input_file     = input_file,
#'   output_dir     = output_dir,
#'   output_file    = output_file,
#'   annotation_file = annotation_file
#' )
build_h5 <- function(
    sam_info       = NULL,
    input_file     = NULL,
    input_coverage = NULL,
    output_dir     = NULL,
    output_file    = NULL,
    annotation_file = NULL,
    chunk_size     = 1e6,
    level          = 1,
    show_progress  = getOption("methylTracer.show_progress", TRUE)
) {
  infile <- input_file
  
  bu_h5_ch(
    infile          = infile,
    output_dir      = output_dir,
    sam_info        = sam_info,
    annotation_file = annotation_file
  )
  
  hdf5_5mc <- file.path(output_dir, output_file)
  if (file.exists(hdf5_5mc)) {
    unlink(hdf5_5mc)
  }
  
  rhdf5::h5createFile(hdf5_5mc)
  rhdf5::h5createGroup(hdf5_5mc, "obs")
  rhdf5::h5createGroup(hdf5_5mc, "var")
  rhdf5::h5createGroup(hdf5_5mc, "uns")
  
  if (show_progress) message("writing X dataset")
  loop_write_h5(
    chunk_size    = chunk_size,
    infile        = infile,
    hdf5_5mc      = hdf5_5mc,
    dataset       = "X",
    level         = level,
    show_progress = show_progress
  )
  
  if (!is.null(input_coverage)) {
    if (show_progress) message("writing coverage dataset")
    loop_write_h5(
      chunk_size    = chunk_size,
      infile        = input_coverage,
      hdf5_5mc      = hdf5_5mc,
      dataset       = "/uns/coverage",
      level         = level,
      show_progress = show_progress
    )
  }
  
  loop_write_h5_col(
    infile    = sam_info,
    hdf5_5mc  = hdf5_5mc,
    dataset   = "/obs/"
  )
  loop_write_h5_col(
    infile    = annotation_file,
    hdf5_5mc  = hdf5_5mc,
    dataset   = "/var/"
  )
  
  rhdf5::h5closeAll()
  if (show_progress) message("Successfully created HDF5 file: ", hdf5_5mc)
  
  invisible(hdf5_5mc)
}

## build_h5 step-1: sanity checks
bu_h5_ch <- function(
    infile          = NULL,
    output_dir      = NULL,
    sam_info        = NULL,
    annotation_file = NULL
) {
  if (is.null(infile) || !file.exists(infile)) {
    stop("Valid 'input_file' must be provided.")
  }
  if (is.null(sam_info) || !file.exists(sam_info)) {
    stop("Valid 'sam_info' file must be provided.")
  }
  if (is.null(annotation_file) || !file.exists(annotation_file)) {
    stop("Valid 'annotation_file' must be provided.")
  }
  if (is.null(output_dir)) {
    stop("'output_dir' must be provided.")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }
}

## build_h5 step-2: stream numeric matrix into HDF5
loop_write_h5 <- function(
    chunk_size    = 1e5,
    infile        = NULL,
    hdf5_5mc      = NULL,
    dataset       = "X",
    level         = 1,
    show_progress = TRUE
) {
  if (is.null(infile) || is.null(hdf5_5mc)) {
    stop("'infile' and 'hdf5_5mc' must not be NULL.")
  }
  
  ## total number of rows & cols
  rownames_all   <- data.table::fread(infile, select = 1)
  file_total_rows <- nrow(rownames_all)
  header         <- data.table::fread(infile, nrows = 0)
  ncols          <- ncol(header) - 1L
  
  rhdf5::h5createDataset(
    file         = hdf5_5mc,
    dataset      = dataset,
    dims         = c(0L, ncols),
    maxdims      = c(rhdf5::H5Sunlimited(), ncols),
    chunk        = c(chunk_size, ncols),
    storage.mode = "double",
    level        = level
  )
  
  row_start    <- 1L
  current_rows <- 0L
  
  ## progress bar
  if (show_progress) {
    pb <- utils::txtProgressBar(min = 0, max = file_total_rows, style = 3)
  } else {
    pb <- NULL
  }
  
  while (row_start <= file_total_rows) {
    rows_to_read <- min(chunk_size, file_total_rows - row_start + 1L)
    
    chunk_dt <- data.table::fread(
      infile,
      nrows      = rows_to_read,
      skip       = row_start - 1L,
      colClasses = c("character", rep("numeric", ncols)),
      na.strings = c("NA", "NaN", "nan")
    )
    
    chunk <- as.matrix(chunk_dt[, -1, with = FALSE])
    storage.mode(chunk) <- "double"
    
    new_rows <- current_rows + nrow(chunk)
    
    rhdf5::h5set_extent(
      file    = hdf5_5mc,
      dataset = dataset,
      dims    = c(new_rows, ncols)
    )
    
    rhdf5::h5write(
      obj   = chunk,
      file  = hdf5_5mc,
      name  = dataset,
      index = list((current_rows + 1L):new_rows, seq_len(ncols))
    )
    
    current_rows <- new_rows
    row_start    <- row_start + rows_to_read
    
    if (!is.null(pb)) {
      utils::setTxtProgressBar(pb, current_rows)
    }
  }
  
  if (!is.null(pb)) close(pb)
}

## write obs/var columns as separate character datasets
loop_write_h5_col <- function(
    infile   = NULL,
    hdf5_5mc = NULL,
    dataset  = NULL
) {
  if (is.null(infile) || is.null(hdf5_5mc) || is.null(dataset)) {
    stop("'infile', 'hdf5_5mc' and 'dataset' must not be NULL.")
  }
  
  header_dt        <- data.table::fread(infile, nrows = 1)
  sample_col_names <- colnames(header_dt)
  file_total_cols  <- length(sample_col_names)
  
  for (i in seq_len(file_total_cols)) {
    current_col_name <- sample_col_names[i]
    current_dataset  <- paste0(dataset, current_col_name)
    
    col_data <- as.character(
      data.table::fread(infile, select = i)[[1]]
    )
    
    rhdf5::h5write(
      obj  = col_data,
      file = hdf5_5mc,
      name = current_dataset
    )
  }
}
