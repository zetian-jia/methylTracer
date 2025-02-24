#' @title Build Standard HDF5 File for methylTracer Input
#'
#' @description This function creates an HDF5 file containing
#'       methylation data and metadata.
#'
#' @details The main methylation data matrix is stored under the
#'      'X' group as a dense matrix in the HDF5 file.
#'
#' @param sam_info path to file containing metadata for each sample.
#'      The first column must be named 'sample_name'.
#' @param input_file Required; dataframe containing methylation signal.
#' @param input_coverage Optional; coverage dataframe
#'       (default from 'methylTracer::build_count_matrix').
#' @param output_dir Path to the output directory where the HDF5
#'      file will be saved.
#' @param output_file HDF5 filename (without extension) for saving
#'      the methylation data.
#' @param annotation_file Genomic annotation file path, as in
#'      `methylTracer::add_annotation`.
#' @param chunk_size Chunk size for reading data, default is 1e5.
#' @param level Compression level (1-9), default is 1.
#'    Higher values provide better compression but may increase
#'    processing time.
#'
#' @return Creates an HDF5 file with:
#' \itemize{
#'   \item The methylation data matrix ('X' group).
#'   \item Sample information ('obs' group).
#'   \item Marker (genomic) information ('var' group).
#'   \item Genomic annotation data.
#' }
#' @importFrom data.table fread
#' @importFrom rhdf5 h5createFile h5createGroup h5createDataset
#' h5set_extent h5write h5read
#' @importFrom utils txtProgressBar
#' @export
#' @examples
#' output_dir <- tempdir()
#'
#' sam_info_df <- data.frame(
#'   sample_name = c(
#'     'sample_1.bed', 'sample_2.bed',
#'     'sample_3.bed', 'sample_4.bed'
#'   ),
#'   group = c(rep('case', 2), rep('ctrl', 2))
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
#'   sample_3.bed = c(200, 400, 700, 1000),
#'   sample_4.bed = c(300, 400, 700, 1000)
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
#' annotation_file <- file.path(
#'   output_dir,
#'   'annotation.bed'
#' )
#'
#' output_file <- 'methylTracer_obj_test.h5'
#'
#' # Remove the output file if it already exists
#' unlink(file.path(output_dir, output_file), recursive = TRUE)
#'
#' # Write the sample information, input data, and annotation
#'
#' write.csv(sam_info_df, sam_info, row.names = FALSE)
#' write.table(input_file_df, input_file, sep = '\t', row.names = FALSE)
#' write.table(annotation_file_df, annotation_file,
#'   sep = '\t', row.names = FALSE
#' )
#'
#' # Run the function
#' build_h5(
#'   sam_info = sam_info,
#'   input_file = input_file,
#'   output_dir = output_dir,
#'   output_file = output_file,
#'   annotation_file = annotation_file
#' )
#'
build_h5 <- function(
    sam_info = NULL, input_file = NULL, 
    input_coverage = NULL, output_dir = NULL,
    output_file = NULL, 
    annotation_file = NULL, chunk_size = 1e+03,
    level = 1
) {
    infile <- input_file
    bu_h5_ch(infile = infile, output_dir = output_dir, 
        sam_info = sam_info, annotation_file = annotation_file
    ) 
    hdf5_5mc <- file.path(output_dir, output_file)
    if (file.exists(hdf5_5mc)) {
        unlink(hdf5_5mc)
    }
    rhdf5::h5createFile(hdf5_5mc)
    rhdf5::h5createGroup(hdf5_5mc, "obs")
    rhdf5::h5createGroup(hdf5_5mc, "var")
    rhdf5::h5createGroup(hdf5_5mc, "uns")
    
    message("writing X dataset")
    loop_write_h5(
        chunk_size = chunk_size, infile = infile, 
        hdf5_5mc = hdf5_5mc, dataset = "X",
        level = level
    )
    if (!is.null(input_coverage)) {
        message("writing coverage dataset")
        loop_write_h5(
            chunk_size = chunk_size, 
            infile = input_coverage, hdf5_5mc = hdf5_5mc,
            dataset = "/uns/coverage", level = level
        )
    }
    loop_write_h5_col(chunk_size, infile = sam_info, 
                      hdf5_5mc, dataset = "/obs/", level = level)
    loop_write_h5_col(
        chunk_size, infile = annotation_file, hdf5_5mc, 
        dataset = "/var/", level = level
    )
    marker_name <- data.table::fread(sam_info,)
    message("Successfully created HDF5 file: ", hdf5_5mc)
}


## build_h5 step-1
bu_h5_ch <- function(
    infile = infile, output_dir = output_dir, 
    sam_info = sam_info, annotation_file = annotation_file
) {
    if (is.null(infile)) {
        message("Input methylation file must be provided")
    }
    if (is.null(sam_info)) {
        message("Sample information must be provided")
    }
    if (is.null(annotation_file)) {
        message("Annotation file must be provided")
    }
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        message("Created output directory: ", output_dir)
    }
}

## build_h5 step-2
## write matrix to X in chunk loop; all colname chunk rows
## infile:path to matrix file in disk; integr
## hdf5_5mc: path to h5 file in disk
loop_write_h5 <- function(
    chunk_size = 1e+05, infile = NULL, 
    hdf5_5mc = NULL, dataset = "X", level = 1
) 
{
    row_start <- 1
    rownames_all <- data.table::fread(infile, select = 1)
    file_total_rows <- dim(rownames_all)[1]
    header <- data.table::fread(infile, nrows = 0)
    ncols <- ncol(header)-1
    total_chunks <- ceiling(file_total_rows / chunk_size)
    h5createDataset(
        file = hdf5_5mc, dataset = dataset, dims = c(0, ncols),
        maxdims = c(H5Sunlimited(), ncols),
        chunk = c(chunk_size, ncols),
        storage.mode = "integer", level = level
        )
    while (row_start <= file_total_rows) {
        rows_to_read <- min(chunk_size, file_total_rows - row_start + 1)
        chunk <- fread(infile, nrows = rows_to_read, skip = row_start - 1)
        chunk <- as.matrix(chunk[, -1, with = FALSE])
        current_rows <- dim(HDF5Array(hdf5_5mc,"X"))[1]
        new_rows <- current_rows + dim(chunk)[1]
        rhdf5::h5set_extent(file = hdf5_5mc, 
                            dataset = dataset, 
                            dims = c(new_rows, ncols))
        rhdf5::h5write(
            obj = chunk, file = hdf5_5mc, 
            name = dataset, 
            index = list((current_rows + 1):new_rows, seq_len(ncols))
        )
        row_start <- row_start + rows_to_read
        message("Processing chunk ", 
                ceiling(row_start / chunk_size), 
                " of ", 
                total_chunks)
    }
}

## write to slot by column, keep colnames
loop_write_h5_col <- function(
    chunk_size = NULL, infile = NULL, 
    hdf5_5mc = NULL, dataset = NULL, level = NULL
) {

    sample_col <- data.table::fread(infile, nrows = 1)
    sample_col_name <- colnames(sample_col)
    file_total_cols <- length(sample_col_name)
    file_total_rows <- dim(data.table::fread(infile, select = 1))[1]
    for (i in seq_len(file_total_cols)) {
        current_col_name <- sample_col_name[i]
        current_dataset <- paste0(dataset, current_col_name)
        col_data <- as.character(data.table::fread(infile, select = i)[[1]])
        rhdf5::h5write(obj = col_data, 
                       file = hdf5_5mc, name = current_dataset)
    }
}
