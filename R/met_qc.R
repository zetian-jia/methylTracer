#' @title QC Metrics for methylTracer Object Methylation Data
#'
#' @description 
#' This function computes quality control (QC) metrics for methylation 
#' data stored in an `methylTracer` object. It calculates the number 
#' of covered CpG sites for each cell, the average methylation level 
#' per cell, the number of cells covering each CpG site, and the 
#' average methylation level for each CpG site.
#'
#' @param met_obj A `methylTracer` object containing methylation data. 
#'                It should have a `seed` slot that contains the data matrix.
#' @param groupname A character string specifying the group name for 
#'                  which the QC values are computed. Default is 'X'.
#'
#' @return This function does not return a value. Instead, it writes 
#'         the computed QC metrics to the associated HDF5 file.
#'         The following datasets are created:
#'         - `obs/coverage_cells`: CpG site coverage per cell.
#'         - `obs/mean_cell_methylation`: Avg. methylation per cell (‰).
#'         - `var/coverage_feature`: Cell coverage per CpG site.
#'         - `var/mean_feature_methylation`: Avg. methylation per CpG site (‰).
#'
#' @details 
#' The methylation levels are recorded in thousandths (‰) for better 
#' numerical precision. The function does not return an object but 
#' updates the corresponding HDF5 datasets.
#'
#' @importFrom HDF5Array HDF5Array writeHDF5Array
#' @import DelayedMatrixStats
#' @import rhdf5
#' @export
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
#' compute_qc_value(met_obj = met_obj)
compute_qc_value <- function(met_obj = NULL, groupname = "X") {
    if (is.null(met_obj)) {
        message("Invalid met_obj: must be a methyTracer object.")
    }
    length1 <- met_obj@seed@dim[1]
    length2 <- met_obj@seed@dim[2] 
    if (!groupname %in% rhdf5::h5ls(met_obj@seed@filepath)$name) {
        message("Invalid groupname: does not exist in the HDF5 file.")
    }
    h5f <- HDF5Array(met_obj@seed@filepath, groupname)
    if (is.null(h5f)) {
        message("Failed to load HDF5 array: h5f is NULL.")
    }
    coverage_cells <- length1 - DelayedMatrixStats::colCounts(h5f, value = NA)
    hdf5_5mc <- met_obj@seed@filepath
    loop_write_h5_vec(
        infile_vec = coverage_cells, 
        sample_col_name = "coverage_cells", 
        hdf5_5mc = hdf5_5mc,
        dataset = "obs/"
    )
    message("add obs/coverage_cells")
    mean_cell_methylation <- DelayedMatrixStats::colMeans2(
      h5f, na.rm = TRUE)/1000
    loop_write_h5_vec(
        infile_vec = mean_cell_methylation, 
        sample_col_name = "mean_cell_methylation",
        hdf5_5mc = hdf5_5mc, dataset = "obs/"
    )
    message("add obs/mean_cell_methylation")
    coverage_feature <- length2 - 
      DelayedMatrixStats::rowCounts(h5f, value = NA)
    loop_write_h5_vec(
        infile_vec = coverage_feature, 
        sample_col_name = "coverage_feature", 
        hdf5_5mc = hdf5_5mc,
        dataset = "var/"
    )
    message("add var/coverage_feature")
    mean_feature_methylation <- DelayedMatrixStats::rowMeans2(
      h5f, na.rm = TRUE)/1000
    loop_write_h5_vec(
        infile_vec = mean_feature_methylation, 
        sample_col_name = "mean_feature_methylation",
        hdf5_5mc = hdf5_5mc, dataset = "var/"
    )
    message("add var/mean_feature_methylation")
}


loop_write_h5_vec <- function(infile_vec = NULL, 
                              sample_col_name = NULL, 
                              hdf5_5mc = NULL, 
                              dataset = NULL) {
    ## Validate inputs
    if (is.null(infile_vec) ||
        !is.numeric(infile_vec)) {
        message("Infile_vec must be a non-null numeric vector.")
    }
    if (is.null(sample_col_name) ||
        !is.character(sample_col_name)) {
        message("Sample_col_name must be a non-null character string.")
    }
    if (is.null(hdf5_5mc) ||
        !file.exists(hdf5_5mc)) {
        message("Hdf5_5mc must be a valid file path.")
    }
    if (is.null(dataset) ||
        !is.character(dataset)) {
        message("Dataset must be a non-null character string.")
    }

    ## Process each column
    current_col_name <- sample_col_name
    current_dataset <- paste0(dataset, current_col_name)

    ## Read entire column at once and convert to character
    col_data <- infile_vec

    dataset_name <- current_dataset
    open_h5 <- rhdf5::H5Fopen(hdf5_5mc)
    has_dataset <- rhdf5::H5Lexists(open_h5, dataset_name)
    rhdf5::H5Fclose(open_h5)

    ## Delete if exists
    if (has_dataset) {
        rhdf5::h5delete(hdf5_5mc, dataset_name)
    }
    rhdf5::h5closeAll()
    rhdf5::h5write(obj = col_data, file = hdf5_5mc, name = current_dataset)

}


#' @title Filter Observations and Variables Based on Cutoff Criteria
#'
#' @description 
#' This function filters observations and variables in a methylation 
#' dataset based on specified cutoff values. It retrieves data from 
#' an HDF5 file and returns a filtered matrix containing only 
#' observations and variables meeting the criteria.
#'
#' @param met_obj A `methylTracer` object containing methylation data. 
#'                It should have a `seed` slot containing the data matrix.
#' @param obs_obj A character string specifying the observation dataset 
#'                to filter (e.g., 'coverage_cells').
#' @param var_obj A character string specifying the variable dataset 
#'                to filter (e.g., 'coverage_feature').
#' @param obs_cutoff A numeric value for filtering observations. 
#'                   Only observations with values ≥ this cutoff are retained.
#' @param var_cutoff A numeric value for filtering variables. 
#'                   Only variables with values ≥ this cutoff are retained.
#' @param sample_name Optional; a character string for sample-based filtering 
#'                    (not implemented).
#' @param marker_name Optional; a character string for marker-based filtering 
#'                    (not implemented).
#'
#' @return A filtered matrix containing only observations and variables 
#'         meeting the specified cutoff criteria.
#'
#' @importFrom HDF5Array HDF5Array writeHDF5Array
#' @import DelayedArray
#' @importFrom methods is
#' @export
#' @examples
#' output_dir <- tempdir()
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
#' compute_qc_value(
#'   met_obj = met_obj
#' )
#'
#' filter_met_obj <- filter_obs_var(
#'   met_obj = met_obj,
#'   obs_cutoff = 0,
#'   var_cutoff = 0,
#'   sample_name = 'sample_name',
#'   marker_name = 'marker_name'
#' )
filter_obs_var <- function(
    met_obj = NULL, obs_obj = "coverage_cells", 
    var_obj = "coverage_feature", 
    obs_cutoff = NULL,
    var_cutoff = NULL, sample_name = NULL, marker_name = NULL
) {
    fi_ch(met_obj=met_obj, obs_obj=obs_obj, var_obj=var_obj,
        obs_cutoff=obs_cutoff, var_cutoff=var_cutoff)

    hdf5_5mc <- met_obj@seed@filepath
    obs_array <- HDF5Array(hdf5_5mc, paste0("obs/", obs_obj))
    var_array <- HDF5Array(hdf5_5mc, paste0("var/", var_obj))
    obs_index <- DelayedArray::which(obs_array >= obs_cutoff)
    var_index <- DelayedArray::which(var_array >= var_cutoff)

    ## Check if indices are valid
    if (length(obs_index) == 0) {
        message("No observations meet the cutoff criteria.")
    }

    X_array <- HDF5Array(hdf5_5mc, "X")
    filtered_X <- X_array[var_index, obs_index]
    ## Create a new HDF5 file and save the data
    filter_file <- paste0("filter_", basename(hdf5_5mc))
    out_dir <- dirname(hdf5_5mc)
    out_file <- file.path(out_dir, filter_file)
    if (file.exists(out_file)) {
        unlink(out_file)
    }
    writeHDF5Array(filtered_X, filepath = out_file, name = "X")
    h5_structure <- h5ls(hdf5_5mc)
    if ("uns/coverage" %in% h5_structure$name) {
        message("Write coverage...")
        h5createGroup(out_file, "uns")
        coverage_array <- HDF5Array(hdf5_5mc, "uns/coverage")
        filtered_coverage <- coverage_array[var_index, obs_index]
        writeHDF5Array(filtered_coverage, 
                       filepath = out_file, name = "uns/coverage")
    }

    me_filter <- fi_lo(met_obj=met_obj, 
                       hdf5_5mc=hdf5_5mc, out_file=out_file,
                       sample_name=sample_name,
                       marker_name=marker_name,
                       obs_index=obs_index,
                       var_index=var_index)

    return(me_filter)
}

## Function to write filtered obs and var data from HDF5Array to a file
writeHDF5Array_index <- function(
    infile_array = NULL, outfile_array = NULL, 
    groupname = NULL, dataset_name = NULL,
    index = NULL
) {
    ## Check if the input array is NULL
    if (is.null(infile_array)) {
        message("Infile_array cannot be NULL.")
    }

    ## Check if the output file is NULL
    if (is.null(outfile_array)) {
        message("Outfile_array cannot be NULL.")
    }

    ## Check if the group name is NULL
    if (is.null(groupname)) {
        message("Groupname cannot be NULL.")
    }

    ## Check if the dataset name is NULL
    if (is.null(dataset_name)) {
        message("Dataset_name cannot be NULL.")
    }

    ## Check if the index is NULL or out of bounds
    if (is.null(index) ||
        any(index < 1) ||
        any(index > length(infile_array))) {
        message("Index must be a valid range within the infile_array.")
    }

    writeHDF5Array(infile_array[index], 
                   outfile_array, 
                   paste0(groupname, "/", dataset_name))
}

## filter_obs_var check step-1
fi_ch <- function(
    met_obj=met_obj,
    obs_obj=obs_obj,
    var_obj=var_obj,
    obs_cutoff=obs_cutoff,
    var_cutoff=var_cutoff
)
{
    if (is.null(met_obj)) {
        message("Met_obj cannot be NULL.")
    }
    if (is.null(obs_obj) ||
        is.null(var_obj)) {
        message("Obs_obj and var_obj must be provided.")
    }
    if (is.null(obs_cutoff) ||
        is.null(var_cutoff)) {
        message("Obs_cutoff and var_cutoff must be provided.")
    }
}
## filter obs_vat loop write step-2
fi_lo <- function(
    met_obj=met_obj,
    hdf5_5mc=hdf5_5mc,
    out_file=out_file,
    sample_name=sample_name,
    marker_name=marker_name,
    obs_index=obs_index,
    var_index=var_index
)
{
    obs_list <- rhdf5::h5ls(hdf5_5mc, recursive = TRUE)
    obs_datasets <- obs_list$name[obs_list$group == "/obs" & 
                                    obs_list$dim == met_obj@seed@dim[2] &
        !is(obs_list$dclass, "COMPOUND")]
    h5createGroup(out_file, "obs")
    total_iterations <- length(obs_datasets)
    for (i in seq_along(obs_datasets)) {
        obs_dataset <- obs_datasets[i]
        writeHDF5Array_index(
            infile_array = HDF5Array(hdf5_5mc, paste0("obs/", obs_dataset)),
            outfile_array = out_file, groupname = "obs", 
            dataset_name = obs_dataset,
            index = obs_index
        )
    }
    h5createGroup(out_file, "var")

    var_list <- rhdf5::h5ls(hdf5_5mc, recursive = TRUE)
    var_datasets <- var_list$name[var_list$group == "/var" & 
                                    var_list$dim == met_obj@seed@dim[1] &
        !methods::is(var_list$dclass, "COMPOUND")]

    total_iterations <- length(var_datasets)
    for (i in seq_along(var_datasets)) {
        var_dataset <- var_datasets[i]
        writeHDF5Array_index(
            infile_array = HDF5Array(hdf5_5mc, paste0("var/", var_dataset)),
            outfile_array = out_file, groupname = "var", 
            dataset_name = var_dataset,
            index = var_index
        )
    }
    met_obj_filtered <- build_met_obj(h5_file = out_file, 
                                      sample_name = sample_name, 
                                      marker_name = marker_name)
    message("Filted file: ", met_obj_filtered@seed@filepath)
    return(met_obj_filtered)
}


#' @title Impute Missing Values in Methylation Data
#'
#' @description 
#' This function imputes missing values in a `methylTracer` object 
#' using the specified method. The imputed values replace missing data 
#' in the HDF5-stored methylation dataset.
#'
#' @param met_obj A `methylTracer` object containing methylation data.
#' @param method The imputation method; options include 'mean' or 'median'. 
#'               Default is 'mean'.
#' @param chunk_size The size of data chunks processed at a time. 
#'                   Default is `1e5`.
#' @param level The compression level for the HDF5 dataset. Default is `1`.
#' @param sample_name A character string representing the sample names.
#' @param marker_name A character string representing the marker names.
#'
#' @return A modified `methylTracer` object with imputed values.
#'
#' @importFrom rhdf5 h5createFile h5createDataset h5set_extent 
#'                    h5write h5ls h5read h5delete h5createGroup
#' @importFrom Rcpp sourceCpp
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
#' filter_met_obj <- filter_obs_var(
#'   met_obj = met_obj,
#'   obs_cutoff = 0,
#'   var_cutoff = 0,
#'   sample_name = 'sample_name',
#'   marker_name = 'marker_name'
#' )
#'
#' met_obj_qc <- impute_met_obj(
#'   met_obj = filter_met_obj,
#'   sample_name = 'sample_name',
#'   marker_name = 'marker_name'
#' )
#'
impute_met_obj <- function(
    met_obj = NULL, method = "mean", chunk_size = 1e+05, 
    level = 1, sample_name = NULL,
    marker_name = NULL
) {
    hdf5_5mc <- met_obj@seed@filepath
    var_list <- rhdf5::h5ls(hdf5_5mc, recursive = TRUE)
    var_datasets <- var_list$name[var_list$group == "/var" & 
                                    var_list$dim == met_obj@seed@dim[1] &
        var_list$dclass != "COMPOUND"]

    if (!any(grepl("mean_feature_methylation", var_datasets))) {
        message("mean_feature_methylation not found, computing...")
        compute_qc_value(met_obj)
    } else {
        message("mean_feature_methylation exist")
    }

    # Create a temporary HDF5 file to store imputation information
    impute_file <- paste0("tmp_impute_", basename(hdf5_5mc))
    out_dir <- dirname(hdf5_5mc)
    tmp_h5 <- file.path(out_dir, impute_file)
    if (file.exists(tmp_h5)) {
        unlink(tmp_h5)
    }
    
    rownames_all <- met_obj@marker_name
    file_total_rows <- met_obj@seed@dim[1]
    ncols <- met_obj@seed@dim[2]
    if (!file.exists(tmp_h5)) {
        h5createFile(tmp_h5)
    }
    met_obj_imputed <- imp_obj(tmp_h5=tmp_h5, ncols=ncols,
        chunk_size=chunk_size, level=level,out_dir=out_dir,
        file_total_rows=file_total_rows, hdf5_5mc=hdf5_5mc,
        sample_name = sample_name, marker_name = marker_name)
    return(met_obj_imputed)
}
## impute_met_obj step-1
imp_obj <- function(
    tmp_h5=tmp_h5, ncols=ncols,
    chunk_size=chunk_size, level=level,
    file_total_rows=file_total_rows,
    hdf5_5mc=hdf5_5mc,out_dir=out_dir,
    sample_name = sample_name, marker_name = marker_name
)
{
    row_start <- 1
    h5createDataset(
        file = tmp_h5, dataset = "X", dims = c(0, ncols),
        maxdims = c(H5Sunlimited(), ncols),
        chunk = c(chunk_size, ncols),
        storage.mode = "integer", level = level
    )
    while (row_start <= file_total_rows) {
        rows_to_read <- min(chunk_size, file_total_rows - row_start + 1)
        chunk <- h5read(
            hdf5_5mc, "X", index = list(
              row_start:(row_start + rows_to_read - 1), seq_len(ncols))
        )
        if (nrow(chunk) == 0)
            break
        chunk_mean <- h5read(
            hdf5_5mc, "var/mean_feature_methylation", 
            index = list(row_start:(row_start + rows_to_read - 1))
        )
        chunk_mean <- round(chunk_mean * 1000, 0)
        chunk_mean <- as.integer(chunk_mean)
        chunk_imputed <- fill_missing_with_mean(chunk, chunk_mean)
        current_rows <- dim(h5read(tmp_h5, "X"))[1]
        new_rows <- current_rows + nrow(chunk)
        h5set_extent(file = tmp_h5, dataset = "X", dims = c(new_rows, ncols))
        rhdf5::h5write(
            obj = chunk_imputed, file = tmp_h5, name = "X", 
            index = list((current_rows + 1):new_rows, seq_len(ncols))
        )
        row_start <- row_start + rows_to_read
    }
    h5delete(hdf5_5mc, "X")
    tmp_h5_X <- HDF5Array(tmp_h5, "X")
    writeHDF5Array(tmp_h5_X, hdf5_5mc, "X")
    unlink(tmp_h5)
    impute_file <- file.path(out_dir, paste0("impute_", basename(hdf5_5mc)))
    system(paste0("mv ", hdf5_5mc, " ", impute_file))
    met_obj_imputed <- build_met_obj(h5_file = impute_file, 
      sample_name = sample_name, marker_name = marker_name)
    return(met_obj_imputed)
}