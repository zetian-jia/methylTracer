plot_met_stats <- function(
    met_obj = NULL, groupname = "X", 
    diff_group_colname = NULL, case_group = NULL,
    control_group = NULL
) {
    ## Process group data and convert to column indices
    met_h5_file <- met_obj@seed@filepath
    ## structure and check for old results
    open_h5 <- rhdf5::H5Fopen(met_h5_file)
    has_old_results <- rhdf5::H5Lexists(open_h5, "uns/met_stats")
    rhdf5::H5Fclose(open_h5)
    if (has_old_results) {
        h5delete(met_h5_file, "uns/met_stats")
    }
    rhdf5::h5closeAll()
    sam_info_df <- met_obj@sample_name
    if (is.null(sam_info_df)) {
        message("Failed to read sample information from HDF5 file")
    }
    diff_type <- rhdf5::h5read(met_h5_file, 
                               paste0("/obs/", diff_group_colname))
    if (methods::is(diff_type, "list")) {
        group_vector <- diff_type$categories[diff_type$codes + 1]
    } else {
        group_vector <- diff_type
    }
    sam_info_df <- data.frame(
      sample_name = met_obj@sample_name, group = as.character(group_vector))

    plot_me_st(met_h5_file=met_h5_file, met_obj=met_obj,
        sam_info_df=sam_info_df, diff_group_colname=diff_group_colname,
        case_group=case_group, control_group=control_group)
}
## plot_met_stats step-1
plot_me_st <- function(
    met_h5_file=met_h5_file,
    met_obj=met_obj,
    sam_info_df=sam_info_df,
    diff_group_colname=diff_group_colname,
    case_group=case_group,
    control_group=control_group
)
{
    open_h5 <- rhdf5::H5Fopen(met_h5_file)
    if (!rhdf5::H5Lexists(open_h5, "var")) {
        rhdf5::H5Fclose(open_h5)
        rhdf5::h5createGroup(met_h5_file, "var")
    }
    rhdf5::H5Fclose(open_h5)
    h5_obj <- HDF5Array(met_obj@seed@filepath, name = "X")
    row_stat_name <- c("rowMeans2", "rowVars")
    for (i in row_stat_name) {
        calc_group_stat(
            h5_obj = h5_obj, 
            sam_info_df = sam_info_df, 
            group_col = diff_group_colname,
            case_group = case_group, 
            control_group = control_group, met_h5_file = met_h5_file,
            stat_type = i, var_group = "var"
        )
    }
    col_stat_name <- c("colMeans2", "colVars")
    for (i in col_stat_name) {
        calc_sample_stat(
            h5_obj = h5_obj, sam_info_df = sam_info_df, 
            met_h5_file = met_h5_file,
            stat_type = i, var_group = "obs"
        )
    }
}




calc_group_stat <- function(
    h5_obj, sam_info_df, group_col = NULL, 
    case_group = NULL, control_group = NULL,
    met_h5_file, stat_type = NULL, var_group = NULL
) {
    open_h5 <- rhdf5::H5Fopen(met_h5_file)
    has_var_group <- rhdf5::H5Lexists(open_h5, var_group)
    rhdf5::H5Fclose(open_h5)
    if (!has_var_group) {
        rhdf5::h5createGroup(met_h5_file, var_group)
    }
    rhdf5::h5closeAll()
    unique_groups <- c(case_group, control_group)
    stat_fun <- switch(
        stat_type, rowMeans2 = function(...) {
            round(
                DelayedMatrixStats::rowMeans2(...)/1000, 3
            )
        },rowVars = function(...) {
            round(
                DelayedMatrixStats::rowVars(...)/1e+06, 3
            )
        }, message("Unsupported statistic type")
    )
    for (group in unique_groups) {
        group_idx <- which(sam_info_df[["group"]] == group)
        ## Calculate statistic for current group
        group_stat <- stat_fun(h5_obj, cols = group_idx, na.rm = TRUE)
        ## Check if dataset exists
        dataset_name <- paste0(var_group, "/", group, "_", stat_type)
        open_h5 <- rhdf5::H5Fopen(met_h5_file)
        has_dataset <- rhdf5::H5Lexists(open_h5, dataset_name)
        rhdf5::H5Fclose(open_h5)
        if (has_dataset) {
            rhdf5::h5delete(met_h5_file, dataset_name)
        }
        rhdf5::h5closeAll()
        rhdf5::h5write(group_stat, file = met_h5_file, name = dataset_name)
    }
    invisible(NULL)
}


calc_sample_stat <- function(
    h5_obj, sam_info_df, 
    met_h5_file, stat_type = NULL, var_group = "obs") {

    open_h5 <- rhdf5::H5Fopen(met_h5_file)
    has_var_group <- rhdf5::H5Lexists(open_h5, var_group)
    rhdf5::H5Fclose(open_h5)

    if (!has_var_group) {
        rhdf5::h5createGroup(met_h5_file, var_group)
    }
    rhdf5::h5closeAll()
    ## Select statistical function
    stat_fun <- switch(
        stat_type, colMeans2 = function(...) {
            round(
                DelayedMatrixStats::colMeans2(...)/1000, 3
            )
        },colVars = function(...) {
            round(
                DelayedMatrixStats::colVars(...)/1e+06,
                3
            )
        }
    )
    group_stat <- stat_fun(h5_obj, na.rm = TRUE)
    dataset_name <- paste0(var_group, "/", stat_type)
    open_h5 <- rhdf5::H5Fopen(met_h5_file)
    has_dataset <- rhdf5::H5Lexists(open_h5, dataset_name)
    rhdf5::H5Fclose(open_h5)
    if (has_dataset) {
        rhdf5::h5delete(met_h5_file, dataset_name)
    }
    rhdf5::h5closeAll()
    rhdf5::h5write(group_stat, file = met_h5_file, name = dataset_name)
    invisible(NULL)
}

#' @title Plot Histogram with Optional Vertical Line
#'
#' @description
#' This function generates a histogram from a specified slot in 
#' the `met_obj` object, with options to customize the number of 
#' bins and add a vertical line at a specified position.
#'
#' @param met_obj A `methylTracer` object containing methylation data.
#' @param slot A character string specifying the slot in the object 
#'             from which the data is read. Example: `'obs/coverage_cells'`. 
#'             Use `methylTracer::compute_qc_values()` to generate slot names.
#' @param bins Integer specifying the number of bins for the histogram. 
#'             Default is `50`.
#' @param vline Numeric value specifying the x-position of a vertical line 
#'              on the histogram. Default is `NULL` (no vertical line).
#'
#' @details
#' The function reads data from the specified slot in `met_obj` and 
#' generates a histogram using `ggplot2`. Optionally, a vertical line 
#' can be added at the specified position (`vline`).
#'
#' @return A `ggplot` object representing the histogram.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' output_dir <- tempdir()
#'
#' sam_info_df <- data.frame(
#'   sample_name = c(
#'     'sample_1.bed', 'sample_2.bed',
#'     'sample_3.bed', 'sample_4.bed',
#'     'sample_5.bed', 'sample_6.bed',
#'     'sample_7.bed', 'sample_8.bed'
#'   ),
#'   group = c(rep('case', 4), rep('ctrl', 4))
#' )
#' sam_info <- file.path(output_dir, 'sample_info.csv')
#'
#' input_file_df <- data.frame(
#'   marker_name = c(
#'     'chr1_1000_2000', 'chr1_2000_3000',
#'     'chr1_3000_4000', 'chr1_4000_5000'
#'   ),
#'   sample_1.bed = c(100, 200, 600, 900),
#'   sample_2.bed = c(100, 200, 600, 900),
#'   sample_3.bed = c(100, 200, 600, 900),
#'   sample_4.bed = c(100, 200, 600, 900),
#'   sample_5.bed = c(800, 900, 100, 50),
#'   sample_6.bed = c(800, 900, 100, 50),
#'   sample_7.bed = c(800, 900, 100, 50),
#'   sample_8.bed = c(800, 900, 100, 50)
#' )
#' input_file <- file.path(output_dir, 'methylTracer_1kb.txt')
#'
#' annotation_file_df <- data.frame(
#'   chr = rep('chr1', 4),
#'   start = c(1000, 2000, 3000, 4000),
#'   end = c(2000, 3000, 4000, 5000),
#'   SYMBOL = c('gene1', 'gene2', 'gene3', 'gene4'),
#'   marker_name = c(
#'     'chr1_1000_2000', 'chr1_2000_3000',
#'     'chr1_3000_4000', 'chr1_4000_5000'
#'   )
#' )
#' annotation_file <- file.path(output_dir, 'annotation.bed')
#'
#' output_file <- 'methylTracer_obj_test.h5'
#'
#' unlink(file.path(output_dir, output_file),
#'   recursive = TRUE
#' )
#'
#' # Write the data to files
#' write.csv(sam_info_df, sam_info,
#'   row.names = FALSE
#' )
#' write.table(input_file_df, input_file,
#'   sep = '\t', row.names = FALSE
#' )
#' write.table(annotation_file_df, annotation_file,
#'   sep = '\t', row.names = FALSE
#' )
#'
#' build_h5(
#'   sam_info = sam_info,
#'   input_file = input_file,
#'   output_dir = output_dir,
#'   output_file = output_file,
#'   annotation_file = annotation_file
#' )
#'
#' met_obj <- build_met_obj(
#'   file.path(output_dir, output_file),
#'   sample_name = 'sample_name',
#'   marker_name = 'marker_name'
#' )
#'
#' compute_qc_value(
#'   met_obj = met_obj
#' )
#'
#' plot_hist(
#'   met_obj = met_obj,
#'   slot = 'obs/coverage_cells',
#'   bins = 50,
#'   vline = 2
#' )

plot_hist <- function(met_obj = NULL, slot = NULL, bins = 50, vline = NULL) {
    # Get the file path of the HDF5 file
    hdf5_5mc <- met_obj@seed@filepath
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
        met_theme + 
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

#' @title Plot Methylation Data (Density Plot)
#'
#' @description 
#' This function creates a density plot (using kernel density estimation)
#' from methylation data stored in a given slot of the `met_obj` object. 
#' It optionally includes a vertical line at a specified position.
#'
#' @param met_obj A methylation object containing the methylation data.
#'        This object should include an HDF5 file with the relevant data.
#' @param slot A character string specifying the specific slot in the HDF5
#'        file to read data from (e.g., 'obs/coverage_cells').
#'        The slot represents the data of interest for the plot.
#' @param bw A numeric value specifying the bandwidth for kernel density 
#'        estimation (KDE). Default is 0.5. This controls the smoothness of 
#'        the density plot.
#' @param vline A numeric value specifying the x-position of a vertical 
#'        line to be drawn on the plot. Default is `NULL`, meaning no line 
#'        will be drawn.
#'
#' @return A `ggplot2` object representing the density plot. Visualizes
#'         the distribution of the values from the specified `slot`.
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' output_dir <- tempdir()
#'
#' sam_info_df <- data.frame(
#'   sample_name = c(
#'     'sample_1.bed', 'sample_2.bed',
#'     'sample_3.bed', 'sample_4.bed',
#'     'sample_5.bed', 'sample_6.bed',
#'     'sample_7.bed', 'sample_8.bed'
#'   ),
#'   group = c(rep('case', 4), rep('ctrl', 4))
#' )
#' sam_info <- file.path(output_dir, 'sample_info.csv')
#'
#' input_file_df <- data.frame(
#'   marker_name = c(
#'     'chr1_1000_2000', 'chr1_2000_3000',
#'     'chr1_3000_4000', 'chr1_4000_5000'
#'   ),
#'   sample_1.bed = c(100, 200, 600, 900),
#'   sample_2.bed = c(100, 200, 600, 900),
#'   sample_3.bed = c(100, 200, 600, 900),
#'   sample_4.bed = c(100, 200, 600, 900),
#'   sample_5.bed = c(800, 900, 100, 50),
#'   sample_6.bed = c(800, 900, 100, 50),
#'   sample_7.bed = c(800, 900, 100, 50),
#'   sample_8.bed = c(800, 900, 100, 50)
#' )
#' input_file <- file.path(output_dir, 'methylTracer_1kb.txt')
#'
#' annotation_file_df <- data.frame(
#'   chr = rep('chr1', 4),
#'   start = c(1000, 2000, 3000, 4000),
#'   end = c(2000, 3000, 4000, 5000),
#'   SYMBOL = c('gene1', 'gene2', 'gene3', 'gene4'),
#'   marker_name = c(
#'     'chr1_1000_2000', 'chr1_2000_3000',
#'     'chr1_3000_4000', 'chr1_4000_5000'
#'   )
#' )
#' annotation_file <- file.path(output_dir, 'annotation.bed')
#'
#' output_file <- 'methylTracer_obj_test.h5'
#'
#' unlink(file.path(output_dir, output_file),
#'   recursive = TRUE
#' )
#'
#' # Write the data to files
#' write.csv(sam_info_df, sam_info,
#'   row.names = FALSE
#' )
#' write.table(input_file_df, input_file,
#'   sep = '\t', row.names = FALSE
#' )
#' write.table(annotation_file_df, annotation_file,
#'   sep = '\t', row.names = FALSE
#' )
#'
#' build_h5(
#'   sam_info = sam_info,
#'   input_file = input_file,
#'   output_dir = output_dir,
#'   output_file = output_file,
#'   annotation_file = annotation_file
#' )
#'
#' met_obj <- build_met_obj(
#'   file.path(output_dir, output_file),
#'   sample_name = 'sample_name',
#'   marker_name = 'marker_name'
#' )
#'
#' compute_qc_value(
#'   met_obj = met_obj
#' )
#'
#' plot_dens(
#'   met_obj = met_obj,
#'   slot = 'obs/coverage_cells',
#'   vline = 2
#' )
plot_dens <- function(met_obj = NULL, slot = NULL, bw = 0.5, vline = NULL) {

    hdf5_5mc <- met_obj@seed@filepath
    tmp_df <- data.frame(tmp_int = as.numeric(h5read(hdf5_5mc, slot)))

    # Kernel Density Estimation
    # (KDE) curve
    p <- ggplot(tmp_df, aes(x = .data$tmp_int)) +
        geom_density(
          color = "steelblue", fill = NA, 
          bw = bw, linewidth = 1.5, na.rm = TRUE) +
        met_theme + scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(
            expand = c(0.05, 0.05),
            limits = c(0, NA)
        ) +
        geom_vline(
          xintercept = vline, color = "black", 
          linetype = "solid", size = 1) +
        labs(title = "", x = "", y = "")

    return(p) 
}

#' @title Plot Scatter Plot with Grouping Variable
#'
#' @description 
#' This function generates a scatter plot from a methylation object, 
#' with customizable axes and a grouping variable for coloring points.
#'
#' @param met_obj A methylation object containing the data.
#' @param xcol The variable to be plotted on the x-axis (e.g., 
#'             `'coverage_cells'`).
#' @param ycol The variable to be plotted on the y-axis (e.g., 
#'             `'mean_cell_methylation'`).
#' @param group_col The variable used for coloring the points (e.g., 
#'                  `'group'`).
#'
#' @return A `ggplot` object representing the scatter plot.
#'
#' @import ggplot2
#' @importFrom rhdf5 h5ls h5read
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @export
#' @examples
#' output_dir <- tempdir()
#'
#' sam_info_df <- data.frame(
#'   sample_name = c(
#'     'sample_1.bed', 'sample_2.bed',
#'     'sample_3.bed', 'sample_4.bed',
#'     'sample_5.bed', 'sample_6.bed',
#'     'sample_7.bed', 'sample_8.bed'
#'   ),
#'   group = c(rep('case', 4), rep('ctrl', 4))
#' )
#' sam_info <- file.path(output_dir, 'sample_info.csv')
#'
#' input_file_df <- data.frame(
#'   marker_name = c(
#'     'chr1_1000_2000', 'chr1_2000_3000',
#'     'chr1_3000_4000', 'chr1_4000_5000'
#'   ),
#'   sample_1.bed = c(100, 200, 600, 900),
#'   sample_2.bed = c(100, 200, 600, 900),
#'   sample_3.bed = c(100, 200, 600, 900),
#'   sample_4.bed = c(100, 200, 600, 900),
#'   sample_5.bed = c(800, 900, 100, 50),
#'   sample_6.bed = c(800, 900, 100, 50),
#'   sample_7.bed = c(800, 900, 100, 50),
#'   sample_8.bed = c(800, 900, 100, 50)
#' )
#' input_file <- file.path(output_dir, 'methylTracer_1kb.txt')
#'
#' annotation_file_df <- data.frame(
#'   chr = rep('chr1', 4),
#'   start = c(1000, 2000, 3000, 4000),
#'   end = c(2000, 3000, 4000, 5000),
#'   SYMBOL = c('gene1', 'gene2', 'gene3', 'gene4'),
#'   marker_name = c(
#'     'chr1_1000_2000', 'chr1_2000_3000',
#'     'chr1_3000_4000', 'chr1_4000_5000'
#'   )
#' )
#' annotation_file <- file.path(output_dir, 'annotation.bed')
#'
#' output_file <- 'methylTracer_obj_test.h5'
#'
#' unlink(file.path(output_dir, output_file),
#'   recursive = TRUE
#' )
#'
#' # Write the data to files
#' write.csv(sam_info_df, sam_info,
#'   row.names = FALSE
#' )
#' write.table(input_file_df, input_file,
#'   sep = '\t', row.names = FALSE
#' )
#' write.table(annotation_file_df, annotation_file,
#'   sep = '\t', row.names = FALSE
#' )
#'
#' build_h5(
#'   sam_info = sam_info,
#'   input_file = input_file,
#'   output_dir = output_dir,
#'   output_file = output_file,
#'   annotation_file = annotation_file
#' )
#'
#' met_obj <- build_met_obj(
#'   file.path(output_dir, output_file),
#'   sample_name = 'sample_name',
#'   marker_name = 'marker_name'
#' )
#'
#' compute_qc_value(
#'   met_obj = met_obj
#' )
#'
#' plot_scatter(
#'   met_obj = met_obj,
#'   xcol = 'coverage_cells',
#'   ycol = 'mean_cell_methylation',
#'   group_col = "group"
#' )
plot_scatter <- function(
    met_obj = NULL, xcol = NULL, 
    ycol = NULL, group_col = NULL) {
    hdf5_5mc <- met_obj@seed@filepath

    # Confirm if the specified variable is in the observations or variables
    group_tmp <- rhdf5::h5ls(hdf5_5mc, recursive = TRUE) %>%
        dplyr::filter(.data$name %in% xcol)

    group_tmp_name <- group_tmp[["group"]]

    # Create a data frame for the scatter plot
    scatter_df <- data.frame(
        x = h5read(hdf5_5mc, paste0(group_tmp_name, "/", xcol)),
        y = h5read(hdf5_5mc, paste0(group_tmp_name, "/", ycol)),
        g = h5read(hdf5_5mc, paste0(group_tmp_name, "/", group_col))
    )  # Read color variable

    # Generate the scatter plot using ggplot2
    p <- ggplot(scatter_df, aes(.data$x, .data$y, color = .data$g)) +
        geom_point(size = 1) +
        met_theme + labs(title = group_col, x = xcol, y = ycol) +
        scale_x_continuous(limits = c(0, NA)) +
        scale_y_continuous(limits = c(0, NA)) +
        guides(color = guide_legend(override.aes = list(size = 6))) +
        theme(
            legend.title = element_blank(), 
            legend.text = element_text(size = 14),
            plot.title = element_text(hjust = 0.5, size = 16),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12)
        )

    return(p)  # Return the plot object
}
# Custom theme
met_theme <- theme(
    axis.title.x = element_text(color = "black", size = 12),
    axis.title.y = element_text(color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    panel.background = element_rect(fill = "gray90"),
    panel.grid.major = element_line(color = "white", linewidth = 1),
    panel.grid.minor = element_line(color = "white", linewidth = 0.5),
    panel.border = element_rect(color = "white", size = 0, fill = NA)
)

#' @title Plot Differentially Methylated Regions (DMRs)
#'
#' @description 
#' This function plots the differentially methylated regions (DMRs) based 
#' on the specified methylation data, DMRs, and gene annotations. The plot 
#' will show the regions of interest with optional shifting, coloring by gene 
#' and regulatory information.
#'
#' @param met_obj A methylation object containing the data, which should 
#' be structured with necessary slots for methylation information.
#' @param pre_calldmrs_res Pre-calculated DMRs results for plotting.
#' @param dmrs_gr GenomicRanges object containing the DMRs.
#' @param region A string specifying the genomic region to plot (e.g., 
#'               "chr1:1000-5000").
#' @param shift_up The number of units to shift the region upward for better 
#' visualization. Default is 3.
#' @param shift_down The number of units to shift the region downward. 
#' Default is 3.
#' @param gene_col The column name that contains gene information for 
#' annotation.
#' @param regul_col The column name that contains regulatory information for 
#' annotation.
#'
#' @details 
#' This function reads in methylation data, DMR information, and gene 
#' annotations to generate a plot. It allows the shifting of the region 
#' of interest up and down for improved visualization and can optionally 
#' display gene and regulatory information.
#'
#' @return A `ggplot2` object representing the DMR plot.
#'
#' @import ggplot2
#' @importFrom rhdf5 h5read
#' @importFrom GenomicRanges GRanges
#' @export
#'
#' @importFrom tidyr starts_with pivot_longer
#'
#' @examples
#' # example code
#' output_dir <- tempdir()
#'
#' sam_info_df <- data.frame(
#'   sample_name = c(
#'     'sample_1.bed', 'sample_2.bed',
#'     'sample_3.bed', 'sample_4.bed',
#'     'sample_5.bed', 'sample_6.bed',
#'     'sample_7.bed', 'sample_8.bed'
#'   ),
#'   group = c(rep('case', 4), rep('ctrl', 4))
#' )
#' sam_info <- file.path(output_dir, 'sample_info.csv')
#'
#' input_file_df <- data.frame(
#'   marker_name = c(
#'     'chr1_1000_2000', 'chr1_2000_3000',
#'     'chr1_3000_4000', 'chr1_4000_5000'
#'   ),
#'   sample_1.bed = c(100, 200, 600, 900),
#'   sample_2.bed = c(100, 200, 600, 900),
#'   sample_3.bed = c(100, 200, 600, 900),
#'   sample_4.bed = c(100, 200, 600, 900),
#'   sample_5.bed = c(800, 900, 100, 50),
#'   sample_6.bed = c(800, 900, 100, 50),
#'   sample_7.bed = c(800, 900, 100, 50),
#'   sample_8.bed = c(800, 900, 100, 50)
#' )
#' input_file <- file.path(output_dir, 'methylTracer_1kb.txt')
#'
#' annotation_file_df <- data.frame(
#'   chr = rep('chr1', 4),
#'   start = c(1000, 2000, 3000, 4000),
#'   end = c(2000, 3000, 4000, 5000),
#'   SYMBOL = c('gene1', 'gene2', 'gene3', 'gene4'),
#'   marker_name = c(
#'     'chr1_1000_2000', 'chr1_2000_3000',
#'     'chr1_3000_4000', 'chr1_4000_5000'
#'   )
#' )
#' annotation_file <- file.path(output_dir, 'annotation.bed')
#'
#' output_file <- 'methylTracer_obj_test.h5'
#'
#' unlink(file.path(output_dir, output_file),
#'   recursive = TRUE
#' )
#'
#' write.csv(sam_info_df, sam_info,
#'   row.names = FALSE
#' )
#' write.table(input_file_df, input_file,
#'   sep = '\t', row.names = FALSE
#' )
#' write.table(annotation_file_df, annotation_file,
#'   sep = '\t', row.names = FALSE
#' )
#'
#' build_h5(
#'   sam_info = sam_info,
#'   input_file = input_file,
#'   output_dir = output_dir,
#'   output_file = output_file,
#'   annotation_file = annotation_file
#' )
#'
#' met_obj <- build_met_obj(
#'   file.path(output_dir, output_file),
#'   sample_name = 'sample_name',
#'   marker_name = 'marker_name'
#' )
#'
#' pre_res <- pre_calldmrs(
#'   met_obj = met_obj,
#'   group_colname = 'group',
#'   case_group = 'case',
#'   control_group = 'ctrl'
#' )
#' dmr_res <- calldmrs_turbo(
#'   met_obj = met_obj,
#'   p_threshold = 0.9,
#'   case_group = "case",
#'   ctrl_group = "ctrl"
#' )
#' plot_dmrs(
#'   met_obj = met_obj,
#'   pre_calldmrs_res = pre_res,
#'   dmrs_gr = dmr_res,
#'   region = "chr1:1000-4000",
#'   shift_up = 1,
#'   shift_down = 1,
#'   gene_col = "SYMBOL"
#' )
plot_dmrs <- function(
    met_obj = NULL, 
    pre_calldmrs_res = NULL, dmrs_gr = NULL, region = NULL, shift_up = 3,
    shift_down = 3, gene_col = "V4", regul_col = NULL
) {
    parts <- strsplit(region, "[:-]")[[1]]
    chr <- parts[1]
    start_n <- as.integer(parts[2])
    end_n <- as.integer(parts[3])
    met_marker <- sub("^([^_]+_[^_]+)_.*$", "\\1", met_obj@marker_name)
    start_n_index <- which(met_marker == paste0(chr, "_", start_n))
    ann_gene <- h5read(
        met_obj@seed@filepath, paste0("/var/", gene_col),
        index = c(start_n_index)
    )
    if (!is.null(regul_col)) {
        ann_regul <- h5read(
            met_obj@seed@filepath, paste0("/var/", regul_col),
            index = c(start_n_index)
        )
    } else {
        ann_regul <- NULL
    }
    start_i <- which(
      pre_calldmrs_res$chr == chr & pre_calldmrs_res$pos == start_n)
    end_i <- which(
      pre_calldmrs_res$chr == chr & pre_calldmrs_res$pos == end_n)
    triangles_sel <- plot_dmr_mt(pre_calldmrs_res=pre_calldmrs_res,
        start_i=start_i, end_i=end_i, shift_up=shift_up,
        shift_down=shift_down
        )
    p <- plot_dmrs_g(triangles_sel=triangles_sel,
        chr=chr, start_n=start_n, end_n=end_n,
        ann_gene=ann_gene, ann_regul=ann_regul)
    return(p)
}
## plot_dmrs step1
plot_dmr_mt <- function(
    pre_calldmrs_res=pre_calldmrs_res,
    start_i=start_i,
    end_i=end_i,
    shift_up=shift_up,
    shift_down=shift_down
)
{
    start_check <- start_i - shift_up
    if (start_check <= 0) {
        start_check <- start_i
    } else {
        start_check <- start_i - shift_up
    }
    triangles <- pre_calldmrs_res[(start_check):(end_i + shift_down), ]
    triangles <- tidyr::pivot_longer(
        triangles, cols = tidyr::starts_with("mean"),
        names_to = "group", values_to = "signal"
    )
    triangles$x1 <- triangles$pos - 1
    triangles$y1 <- 0
    triangles$x2 <- triangles$pos
    triangles$y2 <- triangles$signal
    triangles$x3 <- triangles$pos + 1
    triangles$y3 <- 0
    triangles$id <- paste0(triangles$chr, "_", triangles$pos)

    tri <- triangles[, c("id", "x1", "y1", "x2", "y2", "x3", "y3", "group")]
    triangles <- tri
    triangles <- tidyr::pivot_longer(
        triangles, cols = tidyr::starts_with("x"),
        names_to = "point", values_to = "x"
    )
    triangles$y <- ifelse(
        test = triangles$point == "x1", 
        yes = triangles$y1, 
        no = ifelse(
          test = triangles$point == "x2", 
          yes = triangles$y2, 
          no = triangles$y3)
    )
    triangles_sel <- triangles[, c("id", "x", "y", "group")]
    return(triangles_sel)
}
## plot_dmrs step2
plot_dmrs_g <- function(
    triangles_sel=triangles_sel,
    chr=chr,
    start_n=start_n,
    end_n=end_n,
    ann_gene=ann_gene,
    ann_regul=ann_regul
)
{
    triangles_sel <- triangles_sel %>%
        dplyr::rowwise() %>%
        dplyr::ungroup() %>%
        as.data.frame()
    p <- ggplot(
        data = triangles_sel, 
        mapping = aes(x = .data$x, y = .data$y, 
                      fill = .data$group, color = .data$group)
    )
    p <- p + geom_polygon(aes(group = .data$id),
        alpha = 0.7
    ) +
        geom_hline(yintercept = 0, linewidth = 0.1) +
        facet_wrap(facets = ~.data$group, strip.position = "left", ncol = 1) +
        xlab(label = paste0(as.character(x = chr),
                ":", as.character(x = start_n),
                "-", as.character(x = end_n),
                "\n", ann_gene, "\n", ann_regul
            )
        ) +
        ylab(label = paste0("Methylation signal\n(range: 0-1)")) +
        ylim(c(0, 1)) +
        theme_classic() + theme(
        legend.position = "none", 
        strip.background = element_blank(), 
        strip.text.y.left = element_text(angle = 0),
        panel.spacing.y = unit(x = 0, units = "line"),
        axis.text.y = element_blank()
    )
    df <- data.frame(start = start_n, end = end_n, color = "red")
    p <- p + geom_rect(data = df, inherit.aes = FALSE, 
        aes(xmin = .data$start, 
            xmax = .data$end, ymin = 0, 
            ymax = 1, fill = .data$color),
        color = "black", alpha = 0
    )
    return(p)
}