#' @title Compute QC metrics for a methylTracer object
#'
#' @description
#' Compute basic quality-control (QC) metrics for methylation data stored
#' in a \code{methylTracer} object. Metrics are computed both per cell
#' (column) and per feature/marker (row), and written back into the
#' underlying HDF5 file.
#'
#' @param met A \code{methylTracer} object containing methylation data,
#'   backed by an HDF5 file. The path to the file is taken from
#'   \code{met@seed@filepath}.
#' @param groupname Character scalar giving the name of the dataset in
#'   the HDF5 file to use for QC calculation (default: \code{"X"}).
#'   This should typically be the main methylation matrix, with entries
#'   stored as proportions in \eqn{[0, 1]} and \code{NA} indicating
#'   missing/uncovered sites.
#'
#' @return
#' This function is called for its side effects and does not return a
#' value. It writes the following one-dimensional datasets into the
#' HDF5 file:
#' \itemize{
#'   \item \code{obs/coverage_cells}: number of covered CpG sites
#'         (non-\code{NA}) per cell (column).
#'   \item \code{obs/mean_cell_methylation}: mean methylation proportion
#'         per cell (0–1, excluding \code{NA}).
#'   \item \code{var/coverage_feature}: number of cells (columns)
#'         covering each feature/marker (row).
#'   \item \code{var/mean_feature_methylation}: mean methylation
#'         proportion per feature/marker (0–1, excluding \code{NA}).
#' }
#'
#' Existing datasets with the same names are overwritten.
#'
#' @details
#' Methylation values in \code{groupname} are assumed to be stored as
#' numeric proportions (0–1). Coverage is defined as the number of
#' non-\code{NA} entries. Means are computed with \code{na.rm = TRUE}.
#'
#' @importFrom HDF5Array HDF5Array
#' @importFrom DelayedMatrixStats colCounts colMeans2 rowCounts rowMeans2
#' @importFrom rhdf5 h5ls H5Fopen H5Fclose H5Lexists h5delete h5write h5closeAll
#' @export
#'
#' @examples
#' \donttest{
#' output_dir <- tempdir()
#'
#' ## sample info
#' sam_info_df <- data.frame(
#'   sample_name = paste0("sample_", 1:8),
#'   group       = c(rep("case", 4), rep("ctrl", 4)),
#'   stringsAsFactors = FALSE
#' )
#' sam_info <- file.path(output_dir, "sample_info.csv")
#'
#' ## toy methylation matrix (proportions 0–1)
#' input_file_df <- data.frame(
#'   marker_name = c(
#'     "chr1_1000_2000", "chr1_2000_3000",
#'     "chr1_3000_4000", "chr1_4000_5000"
#'   ),
#'   sample_1 = c(0.10, 0.20, 0.60, 0.90),
#'   sample_2 = c(0.10, 0.20, 0.60, 0.90),
#'   sample_3 = c(0.10, 0.20, 0.60, 0.90),
#'   sample_4 = c(0.10, 0.20, 0.60, 0.90),
#'   sample_5 = c(0.80, 0.90, 0.10, 0.05),
#'   sample_6 = c(0.80, 0.90, 0.10, 0.05),
#'   sample_7 = c(0.80, 0.90, 0.10, 0.05),
#'   sample_8 = c(0.80, 0.90, 0.10, 0.05)
#' )
#' input_file <- file.path(output_dir, "methylTracer_1kb.txt")
#'
#' annotation_file_df <- data.frame(
#'   chr         = rep("chr1", 4),
#'   start       = c(1000L, 2000L, 3000L, 4000L),
#'   end         = c(2000L, 3000L, 4000L, 5000L),
#'   SYMBOL      = c("gene1", "gene2", "gene3", "gene4"),
#'   marker_name = c(
#'     "chr1_1000_2000", "chr1_2000_3000",
#'     "chr1_3000_4000", "chr1_4000_5000"
#'   ),
#'   stringsAsFactors = FALSE
#' )
#' annotation_file <- file.path(output_dir, "annotation.bed")
#'
#' output_file <- "methylTracer_obj_test.h5"
#'
#' unlink(file.path(output_dir, output_file), recursive = TRUE)
#'
#' write.csv(sam_info_df, sam_info, row.names = FALSE)
#' write.table(input_file_df, input_file,
#'             sep = "\t", row.names = FALSE, quote = FALSE)
#' write.table(annotation_file_df, annotation_file,
#'             sep = "\t", row.names = FALSE, quote = FALSE)
#'
#' build_h5(
#'   sam_info       = sam_info,
#'   input_file     = input_file,
#'   output_dir     = output_dir,
#'   output_file    = output_file,
#'   annotation_file = annotation_file
#' )
#'
#' met <- build_met_obj(
#'   file.path(output_dir, output_file),
#'   sample_name = "sample_name",
#'   marker_name = "marker_name"
#' )
#'
#' compute_qc_value(met = met)
#' }
compute_qc_value <- function(met = NULL, groupname = "X") {
  if (is.null(met)) {
    stop("Invalid 'met': must be a non-null methylTracer object.")
  }
  
  h5_path <- met@seed@filepath
  if (!file.exists(h5_path)) {
    stop("The HDF5 file referenced by 'met' does not exist: ", h5_path)
  }
  
  ## check dataset exists
  h5_info <- rhdf5::h5ls(h5_path)
  if (!groupname %in% h5_info$name) {
    stop("Invalid 'groupname': dataset '", groupname,
         "' does not exist in the HDF5 file.")
  }
  
  ## main matrix as HDF5Array
  h5f <- HDF5Array::HDF5Array(h5_path, groupname)
  if (is.null(h5f)) {
    stop("Failed to load HDF5Array from dataset '", groupname, "'.")
  }
  
  x_dim <- dim(h5f)
  n_rows <- x_dim[1]  ## features
  n_cols <- x_dim[2]  ## cells/samples
  
  ## per-cell coverage: #non-NA per column
  coverage_cells <- n_rows - DelayedMatrixStats::colCounts(h5f, value = NA)
  loop_write_h5_vec(
    infile_vec      = coverage_cells,
    sample_col_name = "coverage_cells",
    hdf5_5mc        = h5_path,
    dataset         = "obs/"
  )
  message("Added 'obs/coverage_cells'")
  
  ## mean cell methylation (0–1)
  mean_cell_methylation <- DelayedMatrixStats::colMeans2(h5f, na.rm = TRUE)
  loop_write_h5_vec(
    infile_vec      = mean_cell_methylation,
    sample_col_name = "mean_cell_methylation",
    hdf5_5mc        = h5_path,
    dataset         = "obs/"
  )
  message("Added 'obs/mean_cell_methylation'")
  
  ## per-feature coverage: #non-NA per row
  coverage_feature <- n_cols - DelayedMatrixStats::rowCounts(h5f, value = NA)
  loop_write_h5_vec(
    infile_vec      = coverage_feature,
    sample_col_name = "coverage_feature",
    hdf5_5mc        = h5_path,
    dataset         = "var/"
  )
  message("Added 'var/coverage_feature'")
  
  ## mean feature methylation (0–1)
  mean_feature_methylation <- DelayedMatrixStats::rowMeans2(h5f, na.rm = TRUE)
  loop_write_h5_vec(
    infile_vec      = mean_feature_methylation,
    sample_col_name = "mean_feature_methylation",
    hdf5_5mc        = h5_path,
    dataset         = "var/"
  )
  message("Added 'var/mean_feature_methylation'")
  
  invisible(NULL)
}


loop_write_h5_vec <- function(infile_vec      = NULL,
                              sample_col_name = NULL,
                              hdf5_5mc        = NULL,
                              dataset         = NULL) {
  ## Validate inputs
  if (is.null(infile_vec) || !is.numeric(infile_vec)) {
    stop("'infile_vec' must be a non-null numeric vector.")
  }
  if (is.null(sample_col_name) || length(sample_col_name) != 1L ||
      !is.character(sample_col_name)) {
    stop("'sample_col_name' must be a non-null character scalar.")
  }
  if (is.null(hdf5_5mc) || !file.exists(hdf5_5mc)) {
    stop("'hdf5_5mc' must be a valid HDF5 file path.")
  }
  if (is.null(dataset) || length(dataset) != 1L || !is.character(dataset)) {
    stop("'dataset' must be a non-null character scalar (e.g. 'obs/' or 'var/').")
  }
  
  current_dataset <- paste0(dataset, sample_col_name)
  
  ## check if dataset already exists
  open_h5 <- rhdf5::H5Fopen(hdf5_5mc)
  has_dataset <- rhdf5::H5Lexists(open_h5, current_dataset)
  rhdf5::H5Fclose(open_h5)
  
  ## delete existing dataset if needed
  if (has_dataset) {
    rhdf5::h5delete(hdf5_5mc, current_dataset)
  }
  
  rhdf5::h5closeAll()
  rhdf5::h5write(obj = infile_vec, file = hdf5_5mc, name = current_dataset)
  
  invisible(NULL)
}


#' Filter cells and features in a methylTracer object
#'
#' This function filters cells (columns) and features/markers (rows)
#' based on QC metrics stored in /obs and /var of a methylTracer object,
#' and writes the filtered data into a new HDF5 file. It then rebuilds
#' a new methylTracer object on top of the filtered file.
#'
#' Typical use:
#'   - obs_obj = "coverage_cells"
#'   - var_obj = "coverage_feature"
#' after running compute_qc_value().
#'
#' @param met A \code{methylTracer} object.
#' @param obs_obj Column in \code{metObs(met)} used to filter cells
#'   (e.g. "coverage_cells").
#' @param var_obj Column in \code{metVar(met)} used to filter features
#'   (e.g. "coverage_feature").
#' @param obs_cutoff Numeric; keep cells with \code{obs_obj >= obs_cutoff}.
#'   If \code{NULL}, keep all cells.
#' @param var_cutoff Numeric; keep features with \code{var_obj >= var_cutoff}.
#'   If \code{NULL}, keep all features.
#' @param sample_name Name of the dataset under "/obs" that stores
#'   sample IDs, passed to \code{build_met_obj()}.
#' @param marker_name Name of the dataset under "/var" that stores
#'   marker IDs, passed to \code{build_met_obj()}.
#' @param out_file Optional path to the filtered HDF5 file. If \code{NULL},
#'   a file named \code{filter_<original>.h5} is created in the same
#'   directory as the original HDF5.
#' @param overwrite Logical; overwrite \code{out_file} if it exists.
#' @param copy_coverage Logical; if TRUE, subset and copy "/uns/coverage"
#'   when present; otherwise skip it.
#' @param verbose Logical; print progress messages.
#'
#' @return A new \code{methylTracer} object built on the filtered file.
#' @export
filter_obs_var <- function(
    met,
    obs_obj       = "coverage_cells",
    var_obj       = "coverage_feature",
    obs_cutoff    = NULL,
    var_cutoff    = NULL,
    sample_name   = "sample_name",
    marker_name   = "marker_name",
    out_file      = NULL,
    overwrite     = FALSE,
    copy_coverage = TRUE,
    verbose       = TRUE
) {
  ## basic checks
  if (missing(met) || is.null(met)) {
    stop("'met' must be a non-null methylTracer object.")
  }
  if (!inherits(met, "methylTracer")) {
    stop("'met' must inherit from class 'methylTracer'.")
  }
  
  h5_path <- met@seed@filepath
  if (!file.exists(h5_path)) {
    stop("HDF5 file referenced by 'met' does not exist: ", h5_path)
  }
  
  if (is.null(out_file)) {
    out_file <- file.path(
      dirname(h5_path),
      paste0("filter_", basename(h5_path))
    )
  }
  if (file.exists(out_file)) {
    if (!isTRUE(overwrite)) {
      stop(
        "Output file already exists: ", out_file,
        ". Set 'overwrite = TRUE' to overwrite."
      )
    } else if (verbose) {
      message("Removing existing file: ", out_file)
      unlink(out_file)
    }
  }
  
  ## QC metadata from obs / var
  obs_df <- metObs(met)
  var_df <- metVar(met)
  
  if (!(obs_obj %in% colnames(obs_df))) {
    stop("Column '", obs_obj, "' not found in metObs(met). ",
         "Run 'compute_qc_value()' first or check 'obs_obj'.")
  }
  if (!(var_obj %in% colnames(var_df))) {
    stop("Column '", var_obj, "' not found in metVar(met). ",
         "Run 'compute_qc_value()' first or check 'var_obj'.")
  }
  
  ## build keep indices
  # cells
  if (is.null(obs_cutoff)) {
    obs_keep <- rep(TRUE, nrow(obs_df))
  } else {
    obs_keep <- obs_df[[obs_obj]] >= obs_cutoff
    obs_keep[is.na(obs_keep)] <- FALSE
  }
  obs_index <- which(obs_keep)
  
  # features
  if (is.null(var_cutoff)) {
    var_keep <- rep(TRUE, nrow(var_df))
  } else {
    var_keep <- var_df[[var_obj]] >= var_cutoff
    var_keep[is.na(var_keep)] <- FALSE
  }
  var_index <- which(var_keep)
  
  if (length(obs_index) == 0L) {
    stop("No cells passed the filter (", obs_obj, " >= ", obs_cutoff, ").")
  }
  if (length(var_index) == 0L) {
    stop("No features passed the filter (", var_obj, " >= ", var_cutoff, ").")
  }
  
  if (verbose) {
    message("Keeping ", length(obs_index), " cells and ",
            length(var_index), " features.")
  }
  
  ## subset main matrix X
  if (verbose) message("Filtering dataset 'X' ...")
  X_in  <- HDF5Array::HDF5Array(h5_path, "X")
  X_sub <- X_in[var_index, obs_index, drop = FALSE]
  
  HDF5Array::writeHDF5Array(
    X_sub,
    filepath = out_file,
    name     = "X"
  )
  
  ## subset /uns/coverage if it exists
  h5_ls <- rhdf5::h5ls(h5_path, recursive = TRUE)
  has_coverage <- any(h5_ls$group == "/uns" & h5_ls$name == "coverage")
  
  if (copy_coverage && has_coverage) {
    if (verbose) message("Filtering dataset 'uns/coverage' ...")
    cov_in  <- HDF5Array::HDF5Array(h5_path, "uns/coverage")
    cov_sub <- cov_in[var_index, obs_index, drop = FALSE]
    rhdf5::h5createGroup(out_file, "uns")
    HDF5Array::writeHDF5Array(
      cov_sub,
      filepath = out_file,
      name     = "uns/coverage"
    )
  } else if (copy_coverage && !has_coverage && verbose) {
    message("No '/uns/coverage' dataset found; skipping coverage copy.")
  }
  
  ## write filtered /obs metadata 
  if (verbose) message("Writing filtered '/obs' metadata ...")
  rhdf5::h5createGroup(out_file, "obs")
  
  # 1) sample_name
  sample_vec <- as.vector(met@sample_name[])
  if (length(sample_vec) != nrow(obs_df)) {
    stop("Length of met@sample_name (", length(sample_vec),
         ") does not match nrow(metObs(met)) (", nrow(obs_df), ").")
  }
  rhdf5::h5write(
    obj  = sample_vec[obs_index],
    file = out_file,
    name = paste0("obs/", sample_name)
  )
  
  # 2) other qc
  obs_sub <- obs_df[obs_index, , drop = FALSE]
  for (nm in colnames(obs_sub)) {
    if (nm == sample_name) next
    rhdf5::h5write(
      obj  = obs_sub[[nm]],
      file = out_file,
      name = paste0("obs/", nm)
    )
  }
  
  ## write filtered /var metadata
  if (verbose) message("Writing filtered '/var' metadata ...")
  rhdf5::h5createGroup(out_file, "var")
  
  # 1) marker_name: 
  marker_vec <- as.vector(met@marker_name[])
  if (length(marker_vec) != nrow(var_df)) {
    stop("Length of met@marker_name (", length(marker_vec),
         ") does not match nrow(metVar(met)) (", nrow(var_df), ").")
  }
  rhdf5::h5write(
    obj  = marker_vec[var_index],
    file = out_file,
    name = paste0("var/", marker_name)
  )
  
  # 2) 其他 var 列
  var_sub <- var_df[var_index, , drop = FALSE]
  for (nm in colnames(var_sub)) {
    if (nm == marker_name) next
    rhdf5::h5write(
      obj  = var_sub[[nm]],
      file = out_file,
      name = paste0("var/", nm)
    )
  }
  
  rhdf5::h5closeAll()
  
  ## rebuild methylTracer object
  if (verbose) message("Rebuilding methylTracer object ...")
  met_filtered <- build_met_obj(
    h5_file     = out_file,
    sample_name = sample_name,
    marker_name = marker_name
  )
  
  if (is.null(met_filtered) || !inherits(met_filtered, "methylTracer")) {
    stop("Failed to build methylTracer object from filtered file: ", out_file)
  }
  
  if (verbose) {
    message("Filtered file created: ", met_filtered@seed@filepath)
  }
  
  met_filtered
}



#' @title Impute missing values in a methylTracer object
#'
#' @description
#' Perform simple per-feature imputation of missing methylation values in
#' the main matrix \code{"X"} of a \code{methylTracer} object. Missing
#' entries (\code{NA}) are replaced by the feature-wise mean
#' methylation, as stored in \code{"/var/mean_feature_methylation"}.
#' The imputed matrix is written back to a new HDF5 file and a new
#' \code{methylTracer} object is returned.
#'
#' @param met A \code{methylTracer} object containing methylation data,
#'   backed by an HDF5 file (path taken from \code{met@seed@filepath}).
#' @param method Character scalar specifying the imputation method.
#'   Currently only \code{"mean"} is supported (default).
#' @param chunk_size Integer; number of rows to process per chunk when
#'   reading and writing the matrix \code{"X"} (default: \code{1e5}).
#' @param level Integer scalar giving the compression level (0–9) for
#'   the temporary HDF5 dataset used during imputation (default: \code{1}).
#' @param sample_name Character scalar giving the dataset name in
#'   \code{"/obs"} that contains sample identifiers; passed to
#'   \code{\link{build_met_obj}} (default: \code{"sample_name"}).
#' @param marker_name Character scalar giving the dataset name in
#'   \code{"/var"} that contains marker identifiers; passed to
#'   \code{\link{build_met_obj}} (default: \code{"marker_name"}).
#'
#' @return
#' A new \code{methylTracer} object whose \code{"X"} matrix has had
#' missing values imputed. The new object is backed by a new HDF5 file
#' named \code{"impute_<original.h5>"} stored in the same directory as
#' the original file. The original file is renamed rather than modified
#' in place.
#'
#' @details
#' This function assumes that the main methylation matrix \code{"X"}
#' stores methylation proportions in the range \eqn{[0, 1]}, with
#' \code{NA} indicating missing or uncovered sites, and that the
#' feature-wise mean methylation has been stored in
#' \code{"/var/mean_feature_methylation"} (for example by
#' \code{\link{compute_qc_value}}).
#'
#' If \code{"/var/mean_feature_methylation"} is not present, it will be
#' computed automatically by calling \code{\link{compute_qc_value}}.
#'
#' Imputation is performed in chunks of \code{chunk_size} rows to avoid
#' loading the entire matrix into memory. For each chunk, missing values
#' in a row are replaced by the corresponding entry of
#' \code{mean_feature_methylation}.
#'
#' @importFrom rhdf5 h5createFile h5createDataset h5set_extent
#'   h5write h5ls h5read h5delete h5closeAll H5Sunlimited
#' @importFrom HDF5Array HDF5Array writeHDF5Array
#' @export
#'
#' @examples
#' \donttest{
#' output_dir <- tempdir()
#'
#' sam_info_df <- data.frame(
#'   sample_name = paste0("sample_", 1:8),
#'   group       = c(rep("case", 4), rep("ctrl", 4)),
#'   stringsAsFactors = FALSE
#' )
#' sam_info <- file.path(output_dir, "sample_info.csv")
#'
#' ## toy methylation matrix in [0, 1]
#' input_file_df <- data.frame(
#'   marker_name = c(
#'     "chr1_1000_2000", "chr1_2000_3000",
#'     "chr1_3000_4000", "chr1_4000_5000"
#'   ),
#'   sample_1 = c(0.10, 0.20, 0.60, NA   ),
#'   sample_2 = c(0.10, 0.20, 0.60, 0.90),
#'   sample_3 = c(0.10, NA,   0.60, 0.90),
#'   sample_4 = c(0.10, 0.20, 0.60, 0.90),
#'   sample_5 = c(0.80, 0.90, 0.10, 0.05),
#'   sample_6 = c(0.80, 0.90, 0.10, 0.05),
#'   sample_7 = c(0.80, 0.90, 0.10, 0.05),
#'   sample_8 = c(0.80, 0.90, 0.10, 0.05)
#' )
#' input_file <- file.path(output_dir, "methylTracer_1kb.txt")
#'
#' annotation_file_df <- data.frame(
#'   chr         = rep("chr1", 4),
#'   start       = c(1000L, 2000L, 3000L, 4000L),
#'   end         = c(2000L, 3000L, 4000L, 5000L),
#'   SYMBOL      = c("gene1", "gene2", "gene3", "gene4"),
#'   marker_name = c(
#'     "chr1_1000_2000", "chr1_2000_3000",
#'     "chr1_3000_4000", "chr1_4000_5000"
#'   ),
#'   stringsAsFactors = FALSE
#' )
#' annotation_file <- file.path(output_dir, "annotation.bed")
#'
#' output_file <- "methylTracer_obj_test.h5"
#'
#' unlink(file.path(output_dir, output_file), recursive = TRUE)
#' write.csv(sam_info_df, sam_info, row.names = FALSE)
#' write.table(input_file_df, input_file,
#'   sep = "\t", row.names = FALSE, quote = FALSE
#' )
#' write.table(annotation_file_df, annotation_file,
#'   sep = "\t", row.names = FALSE, quote = FALSE
#' )
#'
#' build_h5(
#'   sam_info        = sam_info,
#'   input_file      = input_file,
#'   output_dir      = output_dir,
#'   output_file     = output_file,
#'   annotation_file = annotation_file
#' )
#'
#' met <- build_met_obj(
#'   file.path(output_dir, output_file),
#'   sample_name = "sample_name",
#'   marker_name = "marker_name"
#' )
#'
#' ## ensure QC metrics exist (including mean_feature_methylation)
#' compute_qc_value(met = met)
#'
#' ## optionally filter first
#' filtered_met <- filter_obs_var(
#'   met         = met,
#'   obs_cutoff  = 1,
#'   var_cutoff  = 1,
#'   sample_name = "sample_name",
#'   marker_name = "marker_name"
#' )
#'
#' ## impute missing values
#' met_imputed <- impute_met_obj(
#'   met         = filtered_met,
#'   method      = "mean",
#'   sample_name = "sample_name",
#'   marker_name = "marker_name"
#' )
#' metX(met)
#' metX(met_imputed)
#' }
impute_met_obj <- function(
    met,
    method      = "mean",
    chunk_size  = 1e5,
    level       = 1,
    sample_name = "sample_name",
    marker_name = "marker_name"
) {
  if (is.null(met)) {
    stop("'met' must be a non-null methylTracer object.")
  }
  if (!identical(method, "mean")) {
    stop("Currently only method = 'mean' is supported.")
  }
  
  hdf5_5mc <- met@seed@filepath
  if (!file.exists(hdf5_5mc)) {
    stop("HDF5 file does not exist: ", hdf5_5mc)
  }
  
  ## ensure mean_feature_methylation exists
  var_list <- rhdf5::h5ls(hdf5_5mc, recursive = TRUE)
  var_datasets <- var_list$name[
    var_list$group == "/var" &
      var_list$dclass != "COMPOUND"
  ]
  
  if (!any(grepl("^mean_feature_methylation$", var_datasets))) {
    message("'var/mean_feature_methylation' not found, computing QC metrics...")
    compute_qc_value(met)
  } else {
    message("'var/mean_feature_methylation' found, using existing values.")
  }
  
  ## dimensions of X
  n_rows <- met@seed@dim[1]
  n_cols <- met@seed@dim[2]
  
  ## create temporary HDF5 file to hold imputed X
  impute_tmp_file <- paste0("tmp_impute_", basename(hdf5_5mc))
  out_dir         <- dirname(hdf5_5mc)
  tmp_h5          <- file.path(out_dir, impute_tmp_file)
  
  if (file.exists(tmp_h5)) {
    unlink(tmp_h5)
  }
  
  rhdf5::h5createFile(tmp_h5)
  met_obj_imputed <- imp_obj(
    tmp_h5          = tmp_h5,
    ncols           = n_cols,
    chunk_size      = chunk_size,
    level           = level,
    file_total_rows = n_rows,
    hdf5_5mc        = hdf5_5mc,
    out_dir         = out_dir,
    sample_name     = sample_name,
    marker_name     = marker_name
  )
  
  met_obj_imputed
}

## impute_met_obj step-1: chunk-wise imputation and file rewrite
imp_obj <- function(
    tmp_h5,
    ncols,
    chunk_size,
    level,
    file_total_rows,
    hdf5_5mc,
    out_dir,
    sample_name,
    marker_name
) {
  if (chunk_size <= 0L) {
    stop("'chunk_size' must be a positive integer.")
  }
  
  row_start    <- 1L
  current_rows <- 0L
  
  ## create extensible dataset "X" in temporary file (double, not integer)
  rhdf5::h5createDataset(
    file         = tmp_h5,
    dataset      = "X",
    dims         = c(0L, ncols),
    maxdims      = c(rhdf5::H5Sunlimited(), ncols),
    chunk        = c(chunk_size, ncols),
    storage.mode = "double",
    level        = level
  )
  
  while (row_start <= file_total_rows) {
    rows_to_read <- min(chunk_size, file_total_rows - row_start + 1L)
    
    ## read chunk from original X
    chunk <- rhdf5::h5read(
      hdf5_5mc,
      "X",
      index = list(
        row_start:(row_start + rows_to_read - 1L),
        seq_len(ncols)
      )
    )
    if (nrow(chunk) == 0L) {
      break
    }
    
    ## read corresponding feature means
    chunk_mean <- rhdf5::h5read(
      hdf5_5mc,
      "var/mean_feature_methylation",
      index = list(row_start:(row_start + rows_to_read - 1L))
    )
    ## chunk_mean is numeric in [0,1]，不再取整
    chunk_mean <- as.numeric(chunk_mean)
    
    ## 使用 C++ 函数填充 NA（注意需要 NumericMatrix 版本，见下文）
    chunk_imputed <- fill_missing_with_mean(chunk, chunk_mean)
    
    ## 扩展临时 X 的维度并写入这个 chunk
    new_rows <- current_rows + nrow(chunk)
    rhdf5::h5set_extent(
      file    = tmp_h5,
      dataset = "X",
      dims    = c(new_rows, ncols)
    )
    rhdf5::h5write(
      obj   = chunk_imputed,
      file  = tmp_h5,
      name  = "X",
      index = list((current_rows + 1L):new_rows, seq_len(ncols))
    )
    
    current_rows <- new_rows
    row_start    <- row_start + rows_to_read
  }
  
  ## 删除原文件中的 X，并用 imputed X 覆盖
  rhdf5::h5delete(hdf5_5mc, "X")
  tmp_h5_X <- HDF5Array::HDF5Array(tmp_h5, "X")
  HDF5Array::writeHDF5Array(tmp_h5_X, hdf5_5mc, "X")
  unlink(tmp_h5)
  
  ## 把原文件整体重命名为 impute_<original.h5>
  impute_file <- file.path(out_dir, paste0("impute_", basename(hdf5_5mc)))
  if (!file.rename(hdf5_5mc, impute_file)) {
    stop("Failed to rename imputed HDF5 file to: ", impute_file)
  }
  
  met_obj_imputed <- build_met_obj(
    h5_file     = impute_file,
    sample_name = sample_name,
    marker_name = marker_name
  )
  message("Imputed file: ", met_obj_imputed@seed@filepath)
  
  met_obj_imputed
}
