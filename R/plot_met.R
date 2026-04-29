
calc_met_stats <- function(
    met = NULL, groupname = "X", 
    diff_group_colname = NULL, case_group = NULL,
    control_group = NULL
) {
  if (is.null(met))               stop("met cannot be NULL")
  if (is.null(diff_group_colname)) stop("diff_group_colname cannot be NULL")
  if (is.null(case_group))        stop("case_group cannot be NULL")
  if (is.null(control_group))     stop("control_group cannot be NULL")
  
  met_h5_file <- met@seed@filepath
  
  ## if /uns/met_stats，delete it
  # suppressWarnings(
  #   rhdf5::h5delete(met_h5_file, "uns/met_stats")
  # )
  safe_h5delete(met_h5_file, "uns/met_stats")
  
  ## group info
  diff_type <- rhdf5::h5read(
    met_h5_file, paste0("/obs/", diff_group_colname)
  )
  if (methods::is(diff_type, "list")) {
    group_vector <- diff_type$categories[diff_type$codes + 1L]
  } else {
    group_vector <- diff_type
  }
  
  sam_info_df <- data.frame(
    sample_name = met@sample_name,
    group       = as.character(group_vector),
    stringsAsFactors = FALSE
  )
  
  calc_me_st(
    met_h5_file        = met_h5_file,
    met                = met,
    sam_info_df        = sam_info_df,
    diff_group_colname = diff_group_colname,
    case_group         = case_group,
    control_group      = control_group
  )
  
  invisible(NULL)
}

calc_me_st <- function(
    met_h5_file,
    met,
    sam_info_df,
    diff_group_colname,
    case_group,
    control_group
) {
  h5_obj <- HDF5Array(met@seed@filepath, name = "X")
  
  row_stat_name <- c("rowMeans2", "rowVars")
  for (st in row_stat_name) {
    calc_group_stat(
      h5_obj      = h5_obj, 
      sam_info_df = sam_info_df, 
      case_group  = case_group, 
      control_group = control_group, 
      met_h5_file = met_h5_file,
      stat_type   = st,
      var_group   = "var"
    )
  }
  
  col_stat_name <- c("colMeans2", "colVars")
  for (st in col_stat_name) {
    calc_sample_stat(
      h5_obj      = h5_obj,
      sam_info_df = sam_info_df, 
      met_h5_file = met_h5_file,
      stat_type   = st,
      var_group   = "obs"
    )
  }
  
  rhdf5::h5closeAll()
  invisible(NULL)
}

calc_group_stat <- function(
    h5_obj, sam_info_df, 
    case_group, control_group,
    met_h5_file, stat_type, var_group
) {
  # comfirm var_group exists
  open_h5 <- rhdf5::H5Fopen(met_h5_file)
  has_group <- rhdf5::H5Lexists(open_h5, var_group)
  rhdf5::H5Fclose(open_h5)
  if (!has_group) {
    rhdf5::h5createGroup(met_h5_file, var_group)
  }
  
  unique_groups <- c(case_group, control_group)
  
  stat_fun <- switch(
    stat_type,
    rowMeans2 = function(...) {
      round(DelayedMatrixStats::rowMeans2(...), 3)
    },
    rowVars = function(...) {
      round(DelayedMatrixStats::rowVars(...), 3)
    },
    stop("Unsupported statistic type: ", stat_type)
  )
  
  for (grp in unique_groups) {
    grp_idx <- which(sam_info_df[["group"]] == grp)
    
    if (!length(grp_idx)) {
      warning("No samples in group: ", grp)
      next
    }
    
    grp_stat <- stat_fun(h5_obj, cols = grp_idx, na.rm = TRUE)
    
    dataset_name <- paste0(var_group, "/", grp, "_", stat_type)
    
    # do not check existence, just try to delete
    # suppressWarnings(
    #   rhdf5::h5delete(met_h5_file, dataset_name)
    # )
    safe_h5delete(met_h5_file, dataset_name)
    
    rhdf5::h5write(
      obj  = grp_stat,
      file = met_h5_file,
      name = dataset_name
    )
  }
  
  invisible(NULL)
}


calc_sample_stat <- function(
    h5_obj, sam_info_df, 
    met_h5_file, stat_type, var_group = "obs"
) {
  # confirm var_group exists
  open_h5 <- rhdf5::H5Fopen(met_h5_file)
  has_group <- rhdf5::H5Lexists(open_h5, var_group)
  rhdf5::H5Fclose(open_h5)
  if (!has_group) {
    rhdf5::h5createGroup(met_h5_file, var_group)
  }
  
  stat_fun <- switch(
    stat_type,
    colMeans2 = function(...) {
      round(DelayedMatrixStats::colMeans2(...), 3)
    },
    colVars = function(...) {
      round(DelayedMatrixStats::colVars(...), 3)
    },
    stop("Unsupported statistic type: ", stat_type)
  )
  
  sample_stat <- stat_fun(h5_obj, na.rm = TRUE)
  
  dataset_name <- paste0(var_group, "/", stat_type)
  
  # suppressWarnings(
  #   rhdf5::h5delete(met_h5_file, dataset_name)
  # )
  safe_h5delete(met_h5_file, dataset_name)
  
  rhdf5::h5write(
    obj  = sample_stat,
    file = met_h5_file,
    name = dataset_name
  )
  
  invisible(NULL)
}


safe_h5delete <- function(file, name) {
  invisible(
    try(
      rhdf5::h5delete(file, name),
      silent = TRUE
    )
  )
}

#' @title Plot histogram of a QC vector from a methylTracer object
#'
#' @description
#' Draw a histogram for a numeric dataset stored in the HDF5 file
#' underlying a \code{methylTracer} object (for example,
#' \code{"obs/coverage_cells"} or \code{"var/mean_feature_methylation"}).
#' Optionally, add a vertical reference line.
#'
#' @param met A \code{methylTracer} object containing methylation data.
#'   The path to the HDF5 file is taken from \code{met@seed@filepath}.
#' @param slot Character scalar giving the HDF5 dataset path to plot,
#'   for example \code{"obs/coverage_cells"} or
#'   \code{"var/coverage_feature"}. Typically these slots are created by
#'   \code{\link{compute_qc_value}}.
#' @param bins Integer giving the number of bins for the histogram
#'   (default: \code{50}).
#' @param vline Optional numeric scalar giving an \eqn{x}-coordinate at
#'   which to draw a vertical reference line. If \code{NULL} (default),
#'   no line is drawn.
#' @param theme A \pkg{ggplot2} theme object used to customize the
#'   appearance of the plot. Default is \code{ggplot2::theme_bw()}.
#'
#' @details
#' The function reads the specified dataset from the HDF5 file, coerces
#' it to numeric, and builds a \pkg{ggplot2} histogram. The
#' \eqn{x}-axis is by default constrained to start at 0 and end slightly
#' above the maximum observed value (10\% padding).
#'
#' This is mainly intended for visualizing QC metrics such as
#' per-cell coverage (\code{"obs/coverage_cells"}) or per-feature
#' coverage (\code{"var/coverage_feature"}).
#'
#' @return
#' A \code{ggplot} object representing the histogram, which can be
#' further modified or printed.
#'
#' @import ggplot2
#' @importFrom rhdf5 h5read h5closeAll
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
#' compute_qc_value(met = met)
#'
#' p <- plot_hist(
#'   met   = met,
#'   slot  = "obs/coverage_cells",
#'   bins  = 50,
#'   vline = 2
#' )
#' p
#' }

plot_hist <- function(met = NULL, slot = NULL, bins = 50, vline = NULL, theme = ggplot2::theme_bw()) {
    # Get the file path of the HDF5 file
    hdf5_5mc <- met@seed@filepath
    # Read data from the specified slot
    tmp_df <- data.frame(tmp_int = as.numeric(h5read(hdf5_5mc, slot)))
    # Calculate the maximum x value for the histogram
    x_max <- (max(tmp_df$tmp_int)/10) +
        max(tmp_df$tmp_int)
    # Create the histogram plot
    p <- ggplot(tmp_df, aes(x = .data$tmp_int)) +
        # Create histogram with specified bins
    geom_histogram(
      bins = bins, fill = "steelblue", color = "white", na.rm = TRUE) +
    theme + 
    geom_vline(
      xintercept = vline, color = "black", linetype = "solid", size = 1) +
        scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(
        expand = c(0.05, 0.05),
        limits = c(0, x_max)
    ) +
        labs(title = "", x = "", y = "")


    return(p)
}


#' @title Plot density of a QC or methylation vector
#'
#' @description
#' Draw a kernel density curve for a numeric dataset stored in the HDF5
#' file underlying a \code{methylTracer} object (for example,
#' \code{"obs/mean_cell_methylation"} or
#' \code{"var/mean_feature_methylation"}). This is particularly useful
#' for visualizing the distribution of methylation proportions in
#' \eqn{[0, 1]}.
#'
#' @param met A \code{methylTracer} object containing methylation data.
#'   The path to the HDF5 file is taken from \code{met@seed@filepath}.
#' @param slot Character scalar giving the HDF5 dataset path to plot,
#'   for example \code{"obs/mean_cell_methylation"} or
#'   \code{"var/mean_feature_methylation"}. The dataset is expected to
#'   be numeric (typically methylation proportions in \eqn{[0, 1]} or
#'   QC-derived summaries).
#' @param bw Numeric scalar giving the bandwidth passed to
#'   \code{ggplot2::geom_density()} (default: \code{0.05}). Smaller
#'   values give a more wiggly curve, larger values a smoother curve.
#' @param vline Optional numeric scalar giving an \eqn{x}-coordinate at
#'   which to draw a vertical reference line. If \code{NULL} (default),
#'   no line is drawn.
#' @param theme A \pkg{ggplot2} theme object used to customize the
#'   appearance of the plot. Default is \code{ggplot2::theme_bw()}.
#'
#' @details
#' This function reads the specified dataset from the HDF5 file,
#' coerces it to numeric, and draws a kernel density estimate using
#' \code{ggplot2::geom_density()}. By default, the \eqn{x}-axis is
#' restricted to \eqn{[0, 1]} to match the usual range of methylation
#' proportions.
#'
#' If you plot a QC vector that is not in \eqn{[0, 1]} (e.g. coverage),
#' you may want to override the default \code{scale_x_continuous()} or
#' rescale the data beforehand.
#'
#' @return
#' A \code{ggplot} object representing the density plot, which can be
#' further modified or printed.
#'
#' @import ggplot2
#' @importFrom rhdf5 h5read h5closeAll
#' @importFrom rlang .data
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
#' compute_qc_value(met = met)
#'
#' ## Example: density of per-cell mean methylation in [0, 1]
#' p <- plot_dens(
#'   met   = met,
#'   slot  = "obs/mean_cell_methylation",
#'   bw    = 0.05,
#'   vline = 0.5
#' )
#' p
#' }
plot_dens <- function(
    met   = NULL,
    slot  = NULL,
    bw    = 0.05,
    vline = NULL,
    theme = ggplot2::theme_bw()
) {
  ## basic checks
  if (is.null(met)) {
    stop("'met' must be a non-null methylTracer object.")
  }
  if (is.null(slot) || !is.character(slot) || length(slot) != 1L) {
    stop("'slot' must be a non-empty character scalar, e.g. 'obs/mean_cell_methylation'.")
  }
  
  hdf5_5mc <- met@seed@filepath
  if (!file.exists(hdf5_5mc)) {
    stop("HDF5 file does not exist: ", hdf5_5mc)
  }
  
  ## read data
  rhdf5::h5closeAll()
  vals <- rhdf5::h5read(hdf5_5mc, slot)
  vals <- as.numeric(vals)
  
  if (!length(vals)) {
    stop("Dataset '", slot, "' is empty.")
  }
  
  tmp_df <- data.frame(value = vals)
  
  p <- ggplot2::ggplot(tmp_df, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_density(
      bw    = bw,
      na.rm = TRUE
    ) +
    theme +
    ggplot2::scale_x_continuous(
      expand = c(0.05, 0.05),
      limits = c(0, 1)
    ) +
    ggplot2::labs(title = NULL, x = NULL, y = NULL)
  
  ## optional vertical line
  if (!is.null(vline)) {
    p <- p + ggplot2::geom_vline(
      xintercept = vline,
      color      = "black",
      linewidth  = 0.5
    )
  }
  
  p
}

#' @title Scatter plot of QC metrics with grouping
#'
#' @description
#' Draw a scatter plot of two sample-level QC vectors stored in the
#' HDF5 file underlying a \code{methylTracer} object, coloured by a
#' grouping variable (e.g. case vs control). A typical use case is
#' plotting per-cell coverage versus mean methylation in \eqn{[0, 1]}.
#'
#' @param met A \code{methylTracer} object containing methylation data.
#'   The path to the HDF5 file is taken from \code{met@seed@filepath}.
#' @param xcol Character scalar giving the name of the dataset under the
#'   \code{"/obs"} group to use for the x-axis, for example
#'   \code{"coverage_cells"}.
#' @param ycol Character scalar giving the name of the dataset under the
#'   \code{"/obs"} group to use for the y-axis, for example
#'   \code{"mean_cell_methylation"} (typically a methylation proportion
#'   in \eqn{[0, 1]}).
#' @param group_col Character scalar giving the name of the dataset
#'   under the \code{"/obs"} group used for colouring points, for
#'   example \code{"group"}.
#' @param theme A \pkg{ggplot2} theme object used to customize the
#'   appearance of the plot. Default is \code{ggplot2::theme_bw()}.
#'
#' @details
#' This function assumes that \code{xcol}, \code{ycol} and
#' \code{group_col} all reside in the \code{"/obs"} group of the HDF5
#' file (i.e. their full paths are \code{"obs/<name>"}). These slots
#' are typically created by \code{\link{compute_qc_value}} (for
#' \code{xcol} / \code{ycol}) and by the user when building the
#' \code{met} object (for \code{group_col}).
#'
#' By default the y-axis is constrained to \eqn{[0, 1]}, which matches
#' the usual range of methylation proportions. If you use a different
#' QC variable on the y-axis, you can override the scale, e.g. via
#' \code{+ ggplot2::scale_y_continuous()}.
#'
#' @return
#' A \code{ggplot} object representing the scatter plot, which can be
#' further modified or printed.
#'
#' @import ggplot2
#' @importFrom rhdf5 h5ls h5read h5closeAll
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
#' compute_qc_value(met = met)
#'
#' p <- plot_scatter(
#'   met       = met,
#'   xcol      = "coverage_cells",
#'   ycol      = "mean_cell_methylation",
#'   group_col = "group"
#' )
#' p
#' }
plot_scatter <- function(
    met,
    xcol,
    ycol,
    group_col,
    theme = ggplot2::theme_bw()
) {
  ## basic checks
  if (missing(met) || is.null(met)) {
    stop("'met' must be a non-null methylTracer object.")
  }
  if (missing(xcol) || is.null(xcol) ||
      !is.character(xcol) || length(xcol) != 1L) {
    stop("'xcol' must be a non-empty character scalar, e.g. 'coverage_cells'.")
  }
  if (missing(ycol) || is.null(ycol) ||
      !is.character(ycol) || length(ycol) != 1L) {
    stop("'ycol' must be a non-empty character scalar, e.g. 'mean_cell_methylation'.")
  }
  if (missing(group_col) || is.null(group_col) ||
      !is.character(group_col) || length(group_col) != 1L) {
    stop("'group_col' must be a non-empty character scalar, e.g. 'group'.")
  }
  
  hdf5_5mc <- met@seed@filepath
  if (!file.exists(hdf5_5mc)) {
    stop("HDF5 file does not exist: ", hdf5_5mc)
  }
  
  rhdf5::h5closeAll()
  h5_info <- rhdf5::h5ls(hdf5_5mc, recursive = TRUE)
  
  ## we expect all three under /obs
  needed <- c(xcol, ycol, group_col)
  obs_rows <- h5_info$group == "/obs"
  obs_names <- h5_info$name[obs_rows]
  
  missing_slots <- setdiff(needed, obs_names)
  if (length(missing_slots) > 0L) {
    stop(
      "The following datasets are missing under '/obs': ",
      paste(missing_slots, collapse = ", "),
      ". If these are QC metrics, please run compute_qc_value() first."
    )
  }
  
  x <- rhdf5::h5read(hdf5_5mc, paste0("obs/", xcol))
  y <- rhdf5::h5read(hdf5_5mc, paste0("obs/", ycol))
  g <- rhdf5::h5read(hdf5_5mc, paste0("obs/", group_col))
  
  scatter_df <- data.frame(
    x = as.numeric(x),
    y = as.numeric(y),
    g = as.factor(g)
  )
  
  p <- ggplot2::ggplot(scatter_df, ggplot2::aes(x = x, y = y, colour = g)) +
    ggplot2::geom_point(size = 1, na.rm = TRUE) +
    theme +
    ggplot2::labs(title = group_col, x = xcol, y = ycol, colour = NULL) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(override.aes = list(size = 4))
    ) +
    ## y in [0, 1] for methylation proportions; override if needed
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0.05, 0.05))
  
  p
}

# Custom theme
met_theme <- theme(
    axis.title.x = ggplot2::element_text(color = "black", size = 12),
    axis.title.y = ggplot2::element_text(color = "black", size = 12),
    axis.text.x = ggplot2::element_text(color = "black", size = 10),
    axis.text.y = ggplot2::element_text(color = "black", size = 10),
    axis.line = ggplot2::element_line(color = "black", linewidth = 0.5),
    axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = ggplot2::unit(0.2, "cm"),
    panel.background = ggplot2::element_rect(fill = "gray90"),
    panel.grid.major = ggplot2::element_line(color = "white", linewidth = 1),
    panel.grid.minor = ggplot2::element_line(color = "white", linewidth = 0.5),
    panel.border = ggplot2::element_rect(color = "white", linewidth = 0, fill = NA)
)




#' @title Plot Differentially Methylated Regions (DMRs)
#'
#' @description
#' Visualise differentially methylated regions (DMRs) over a genomic
#' interval using per-group mean methylation profiles stored in
#' \code{pre_calldmrs_res}. For each group (e.g. case/ctrl) the function
#' draws small “triangles” at CpG positions, where the height encodes the
#' mean methylation proportion in \eqn{[0, 1]}.
#'
#' @param met A \code{methylTracer} object. The underlying HDF5 file is
#'   taken from \code{met@seed@filepath}, and marker coordinates are
#'   parsed from \code{met@marker_name} (format
#'   \code{"chr_start_end"}).
#' @param pre_calldmrs_res A data.frame returned by
#'   \code{\link{pre_calldmrs}}. It must contain at least columns
#'   \code{chr}, \code{pos} (CpG position), and group-wise mean
#'   methylation columns whose names start with \code{"mean"} (e.g.
#'   \code{"mean_case"} / \code{"mean_ctrl"}) in \eqn{[0, 1]}.
#' @param region Character scalar specifying the genomic interval to
#'   plot, in the form \code{"chr:start-end"}, e.g.
#'   \code{"chr1:1000-5000"}. \code{start} and \code{end} should match
#'   values in \code{calldmrs_turbo}.
#' @param shift_up Integer, how many CpG indices upstream of
#'   \code{start} to include for context (default \code{3}).
#' @param shift_down Integer, how many CpG indices downstream of
#'   \code{end} to include for context (default \code{3}).
#' @param gene_col Optional character scalar giving the name of a 1D
#'   dataset under \code{"var/"} (e.g. \code{"SYMBOL"}). The value at
#'   the end marker will be shown in the x-axis label as gene
#'   annotation.
#' @param regul_col Optional character scalar giving the name of a 1D
#'   dataset under \code{"var/"} (e.g. a regulatory category). The value
#'   at the end marker will also be shown in the axis label.
#'
#' @details
#' \itemize{
#'   \item The \code{region} boundaries (\code{start}, \code{end}) are
#'   matched against \code{pre_calldmrs_res$pos} on the specified
#'   chromosome. If either boundary cannot be found an error is raised.
#'
#'   \item The plotted y-axis corresponds to mean methylation
#'   \emph{proportions} in \eqn{[0, 1]} (consistent with storing
#'   methylation as 0–1 values in the package).
#'
#'   \item Upstream/downstream context is controlled by \code{shift_up}
#'   and \code{shift_down} in terms of CpG \emph{indices} rather than
#'   base pairs.
#' }
#'
#' Internally this function calls two helpers:
#' \itemize{
#'   \item \code{plot_dmr_mt()} – builds the triangle polygons per CpG
#'   and group;
#'   \item \code{plot_dmrs_g()} – draws the \pkg{ggplot2} figure.
#' }
#'
#' @return A \code{ggplot} object showing the per-group methylation
#'   profiles over the requested region.
#'
#' @import ggplot2
#' @importFrom rhdf5 h5read h5closeAll
#' @importFrom tidyr starts_with pivot_longer
#' @importFrom rlang .data
#' @importFrom grid unit
#' @export
#' @examples
#' \donttest{
#' ## minimal workflow example (0–1 methylation values)
#' output_dir <- tempdir()
#'
#' ## sample metadata
#' sam_info_df <- data.frame(
#'   sample_name = paste0("sample_", 1:8, ".bed"),
#'   group       = c(rep("case", 4), rep("ctrl", 4)),
#'   stringsAsFactors = FALSE
#' )
#' sam_info <- file.path(output_dir, "sample_info.csv")
#'
#' ## toy methylation matrix (0–1 scale)
#' input_file_df <- data.frame(
#'   marker_name   = c(
#'     "chr1_1000_2000",
#'     "chr1_2000_3000",
#'     "chr1_3000_4000",
#'     "chr1_4000_5000"
#'   ),
#'   sample_1.bed = c(0.10, 0.20, 0.60, 0.90),
#'   sample_2.bed = c(0.10, 0.20, 0.60, 0.90),
#'   sample_3.bed = c(0.10, 0.20, 0.60, 0.90),
#'   sample_4.bed = c(0.10, 0.20, 0.60, 0.90),
#'   sample_5.bed = c(0.80, 0.90, 0.10, 0.05),
#'   sample_6.bed = c(0.80, 0.90, 0.10, 0.05),
#'   sample_7.bed = c(0.80, 0.90, 0.10, 0.05),
#'   sample_8.bed = c(0.80, 0.90, 0.10, 0.05),
#'   check.names = FALSE
#' )
#' input_file <- file.path(output_dir, "methylTracer_1kb.txt")
#'
#' ## annotation (coordinates consistent with marker_name)
#' annotation_file_df <- data.frame(
#'   chr        = rep("chr1", 4),
#'   start      = c(1000, 2000, 3000, 4000),
#'   end        = c(2000, 3000, 4000, 5000),
#'   SYMBOL     = paste0("gene", 1:4),
#'   marker_name = c(
#'     "chr1_1000_2000",
#'     "chr1_2000_3000",
#'     "chr1_3000_4000",
#'     "chr1_4000_5000"
#'   ),
#'   stringsAsFactors = FALSE
#' )
#' annotation_file <- file.path(output_dir, "annotation.bed")
#'
#' output_file <- "methylTracer_obj_test.h5"
#'
#' ## write input files
#' unlink(file.path(output_dir, output_file), recursive = TRUE)
#' write.csv(sam_info_df, sam_info, row.names = FALSE)
#' write.table(input_file_df, input_file, sep = "\t", row.names = FALSE)
#' write.table(annotation_file_df, annotation_file,
#'             sep = "\t", row.names = FALSE)
#'
#' ## build HDF5 and methylTracer object
#' build_h5(
#'   sam_info       = sam_info,
#'   input_file     = input_file,
#'   output_dir     = output_dir,
#'   output_file    = output_file,
#'   annotation_file = annotation_file
#' )
#'
#' met <- build_met_obj(
#'   h5_file     = file.path(output_dir, output_file),
#'   sample_name = "sample_name",
#'   marker_name = "marker_name"
#' )
#'
#' ## pre-compute stats and call DMRs
#' pre_res <- pre_calldmrs(
#'   met           = met,
#'   group_colname = "group",
#'   case_group    = "case",
#'   control_group = "ctrl"
#' )
#'
#' dmr_res <- calldmrs_turbo(
#'   met         = met,
#'   p_threshold = 0.9,
#'   case_group  = "case",
#'   ctrl_group  = "ctrl"
#' )
#'
#' ## pick a region overlapping one of the markers
#' p <- plot_dmrs(
#'   met              = met,
#'   pre_calldmrs_res = pre_res,
#'   region           = "chr1:2000-5000",
#'   shift_up         = 1,
#'   shift_down       = 1,
#'   gene_col         = "SYMBOL"
#' )
#' print(p)
#' }
plot_dmrs <- function(
    met,
    pre_calldmrs_res,
    region,
    shift_up   = 3L,
    shift_down = 3L,
    gene_col   = NULL,
    regul_col  = NULL
) {
  ## basic checks
  if (missing(met) || is.null(met)) {
    stop("'met' must be a non-null methylTracer object.")
  }
  if (missing(pre_calldmrs_res) || is.null(pre_calldmrs_res)) {
    stop("'pre_calldmrs_res' must be provided (output of pre_calldmrs()).")
  }
  if (missing(region) || is.null(region) || !is.character(region) ||
      length(region) != 1L) {
    stop("'region' must be a character scalar like 'chr1:2000-5000'.")
  }
  
  ## parse region: "chr:start-end"
  parts <- strsplit(region, "[:-]")[[1L]]
  if (length(parts) != 3L) {
    stop("Invalid 'region' format. Expected 'chr:start-end'.")
  }
  chr     <- parts[1L]
  start_n <- as.integer(parts[2L])
  end_n   <- as.integer(parts[3L])
  if (is.na(start_n) || is.na(end_n)) {
    stop("'start' and 'end' in 'region' must be integers.")
  }
  
  ## marker_name format: "chr_start_end"
  split_result <- strsplit(as.character(met@marker_name), "_", fixed = TRUE)
  chrtmp <- vapply(split_result, `[`, character(1L), 1L)
  endtmp <- vapply(split_result, `[`, character(1L), 3L)
  met_marker <- paste0(chrtmp, "_", endtmp)
  rm(split_result, chrtmp, endtmp)
  
  ## use end coordinate to pick an annotation row (first match)
  annI <- which(met_marker == paste0(chr, "_", end_n))
  if (length(annI) > 0L) {
    annI <- annI[1L]
  }
  
  ann_gene  <- NA_character_
  ann_regul <- NA_character_
  
  if (length(annI) > 0L) {
    h5file <- met@seed@filepath
    if (!is.null(gene_col)) {
      ann_gene <- rhdf5::h5read(
        h5file,
        paste0("var/", gene_col),
        index = list(annI)
      )
    }
    if (!is.null(regul_col)) {
      ann_regul <- rhdf5::h5read(
        h5file,
        paste0("var/", regul_col),
        index = list(annI)
      )
    }
  }
  
  ## subset stats to this chromosome
  sel <- pre_calldmrs_res[pre_calldmrs_res$chr == chr, , drop = FALSE]
  if (!nrow(sel)) {
    stop("No rows found in 'pre_calldmrs_res' for chromosome ", chr, ".")
  }
  
  start_i <- which(sel$pos == start_n)
  end_i   <- which(sel$pos == end_n)
  if (!length(start_i) || !length(end_i)) {
    stop(
      "Cannot find 'start' or 'end' position in pre_calldmrs_res$pos ",
      "for the requested region."
    )
  }
  ## if multiple matches, use the first
  start_i <- start_i[1L]
  end_i   <- end_i[1L]
  if (start_i > end_i) {
    tmp     <- start_i
    start_i <- end_i
    end_i   <- tmp
  }
  
  triangles_sel <- plot_dmr_mt(
    sel       = sel,
    start_i   = start_i,
    end_i     = end_i,
    shift_up  = shift_up,
    shift_down = shift_down
  )
  
  plot_dmrs_g(
    triangles_sel = triangles_sel,
    chr           = chr,
    start_n       = start_n,
    end_n         = end_n,
    ann_gene      = ann_gene,
    ann_regul     = ann_regul
  )
}


## step 1: build triangle polygons per CpG and group
plot_dmr_mt <- function(
    sel,
    start_i,
    end_i,
    shift_up   = 3L,
    shift_down = 3L
) {
  if (!is.data.frame(sel) || !nrow(sel)) {
    stop("'sel' must be a non-empty data.frame (subset of pre_calldmrs_res).")
  }
  
  n <- nrow(sel)
  start_u <- max(1L, start_i - as.integer(shift_up))
  end_u   <- min(n, end_i + as.integer(shift_down))
  
  sub <- sel[start_u:end_u, , drop = FALSE]
  
  sub_long <- tidyr::pivot_longer(
    sub,
    cols      = tidyr::starts_with("mean"),
    names_to  = "group",
    values_to = "signal"
  )
  
  if (nrow(sub) >= 2L) {
    dif <- abs(sub$pos[2L] - sub$pos[1L])
  } else {
    dif <- 1L
  }
  shif <- if (dif < 50) 1 else dif / 20
  
  sub_long$x1 <- sub_long$pos - shif
  sub_long$y1 <- 0
  sub_long$x2 <- sub_long$pos
  sub_long$y2 <- sub_long$signal
  sub_long$x3 <- sub_long$pos + shif
  sub_long$y3 <- 0
  sub_long$id <- paste0(sub_long$chr, "_", sub_long$pos)
  
  tri <- sub_long[, c("id", "x1", "y1", "x2", "y2", "x3", "y3", "group")]
  
  triangles <- tidyr::pivot_longer(
    tri,
    cols      = tidyr::starts_with("x"),
    names_to  = "point",
    values_to = "x"
  )
  
  triangles$y <- ifelse(
    triangles$point == "x1", triangles$y1,
    ifelse(triangles$point == "x2", triangles$y2, triangles$y3)
  )
  
  triangles_sel <- triangles[, c("id", "x", "y", "group")]
  triangles_sel
}

plot_dmrs_g <- function(
    triangles_sel,
    chr,
    start_n,
    end_n,
    ann_gene  = NA_character_,
    ann_regul = NA_character_
) {
  triangles_sel <- as.data.frame(triangles_sel)
  if (!nrow(triangles_sel)) {
    stop("No points to plot in 'triangles_sel'.")
  }
  
  yl <- max(triangles_sel$y, na.rm = TRUE)
  if (!is.finite(yl) || yl <= 0) {
    yl <- 1
  }
  
  df_region <- data.frame(
    start = start_n,
    end   = end_n,
    color = "region"
  )
  
  label_gene <- if (!is.null(ann_gene)  && length(ann_gene)  &&
                    !is.na(ann_gene[1L]))  paste0("Gene: ", ann_gene[1L]) else ""
  label_reg  <- if (!is.null(ann_regul) && length(ann_regul) &&
                    !is.na(ann_regul[1L])) paste0("Regulatory: ", ann_regul[1L]) else ""
  
  xlab_txt <- paste0(
    chr, ":", start_n, "-", end_n,
    if (label_gene != "")  paste0("\n", label_gene)  else "",
    if (label_reg  != "")  paste0("\n", label_reg)   else ""
  )
  
  p <- ggplot2::ggplot(
    triangles_sel,
    ggplot2::aes(
      x     = .data$x,
      y     = .data$y,
      fill  = .data$group,
      colour = .data$group
    )
  ) +
    ggplot2::geom_polygon(ggplot2::aes(group = .data$id), alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.1) +
    ggplot2::facet_wrap(~ .data$group, strip.position = "left", ncol = 1) +
    ggplot2::xlab(xlab_txt) +
    ggplot2::ylab(
      paste0(
        "Mean methylation proportion\n(range: 0-",
        signif(yl, 3), ")"
      )
    ) +
    ggplot2::ylim(c(0, yl)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position    = "none",
      strip.background   = ggplot2::element_blank(),
      strip.text.y.left  = ggplot2::element_text(angle = 0),
      panel.spacing.y    = grid::unit(0, "line"),
      axis.text.y        = ggplot2::element_blank()
    ) +
    ggplot2::geom_rect(
      data        = df_region,
      inherit.aes = FALSE,
      ggplot2::aes(
        xmin = .data$start,
        xmax = .data$end,
        ymin = 0,
        ymax = yl
      ),
      colour  = "grey50",
      alpha   = 0,
      linetype = "dashed"
    )
  
  p
}





