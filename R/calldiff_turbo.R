#' @title Call Differential Methylation in Windows or Base in Presto
#' 
#' @description This function calls differential methylation in
#'              windows or base-pairs using the turbo algorithm.
#'
#' @param met_obj Methylation object containing methylation data and
#'               sample information.
#' @param mean_diff_abs Minimum absolute methylation difference
#'                     threshold between groups (default: 0.5).
#' @param groupname Name for the differential analysis (default: 'X').
#' @param diff_group_colname Column name in sample info containing
#'                          group labels.
#' @param chunk_size Number of rows to process at a time (default: 5e4).
#'
#' @details This function uses the turbo algorithm to call differential
#'           methylation in windows or base-pairs. It calculates mean
#'           methylation levels for each group using chunks of data. Then,
#'           it finds significant differences by comparing the absolute
#'           differences between group means for each row. Results are stored
#'           in an HDF5 file with the following structure:
#' \itemize{
#' \item \code{uns/diff_results}: Results of the differential analysis.
#' \item \code{uns/diff_results_idx}: Indices of significant rows.
#' }
#' @importFrom HDF5Array HDF5Array
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate group_split group_by arrange all_of
#' @import DelayedMatrixStats
#' @import presto
#' @import rhdf5
#' @importFrom utils prompt .DollarNames stack relist
#' @importFrom IRanges shift collapse union slice intersect setdiff desc
#' @import BiocParallel
#'
#' @return A data frame containing the results of the differential analysis.
#'
#' @export
#'
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
#' unlink(file.path(output_dir, output_file),
#'   recursive = TRUE
#' )
#' write.csv(sam_info_df, sam_info,
#'   row.names = FALSE)
#' write.table(input_file_df, input_file,
#'   sep = '\t', row.names = FALSE)
#' write.table(annotation_file_df, annotation_file,
#'   sep = '\t', row.names = FALSE)
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
#' # Call differential methylation base presto
#' diff_results <- calldiff_turbo(
#'   met_obj = met_obj,
#'   mean_diff_abs = 0.01,
#'   diff_group_colname = 'group'
#' )
#'
calldiff_turbo <- function(
    met_obj=NULL, mean_diff_abs=0.5, 
    groupname="X", diff_group_colname=NULL,
    chunk_size=50000
) 
{
    infile <- calldiff_turbo_check(
        diff_group_colname=diff_group_colname,
        mean_diff_abs=mean_diff_abs,
        met_obj=met_obj)
    mean_diff_abs <- mean_diff_abs * 1000
    sam_info_df <- data.frame(
        sample_name=met_obj@sample_name, 
        group=as.character(
          rhdf5::h5read(
            met_obj@seed@filepath, paste0("/obs/", diff_group_colname)))
    )
    significant_indices <- sig_ind(
        sam_info_df=sam_info_df,
        met_obj=met_obj,
        mean_diff_abs=mean_diff_abs)
    selected_rows <- Reduce(union, significant_indices)
    BATCH_SIZE <- chunk_size
    ## Initialize HDF5 structure and check for old results
    open_h5 <- rhdf5::H5Fopen(infile)
    has_old_results <- rhdf5::H5Lexists(open_h5, "uns/diff_results")
    has_old_idx <- rhdf5::H5Lexists(open_h5, "uns/diff_results_idx")
    rhdf5::H5Fclose(open_h5)
    ## Delete old results if they exist
    if (has_old_results) {
        h5delete(infile, "uns/diff_results")
    }
    if (has_old_idx) {
        h5delete(infile, "uns/diff_results_idx")
    }
    rhdf5::h5closeAll()
    sample_batch <- rhdf5::h5read(met_obj@seed@filepath, "X", 
        index=list(seq_len(2), NULL))
    diff_results <- pr_result(infile=infile,
        open_h5=open_h5,
        met_obj=met_obj,
        sample_batch=sample_batch,
        sam_info_df=sam_info_df,
        BATCH_SIZE=BATCH_SIZE,
        selected_rows=selected_rows)
    return(diff_results)
}

## calldiff_turbo step-1
calldiff_turbo_check <- function(
    diff_group_colname=diff_group_colname,
    mean_diff_abs=mean_diff_abs,
    met_obj=met_obj
)
{
    if (is.null(diff_group_colname)) {
        message("Group column name cannot be NULL")
    }
    # Validate methylation difference threshold
    if (mean_diff_abs < 0 || mean_diff_abs > 1) {
        message("Methylation difference threshold must be between 0 and 1")
    }
    # Process group data and convert to column indices
    infile <- met_obj@seed@filepath
    if (!file.exists(infile)) {
        message("Input file not found: ", infile)
    }
    sam_info_df <- met_obj@sample_name
    if (is.null(sam_info_df)) {
        message("Failed to read sample information from HDF5 file")
    }
    return(infile)
}
## calldiff_turbo step-2
sig_ind <- function(
    sam_info_df=sam_info_df,
    met_obj=met_obj,
    mean_diff_abs=mean_diff_abs
)
{
    group_indices <- sam_info_df %>%
        dplyr::select(dplyr::all_of("group")) %>%
        dplyr::mutate(row_index = sam_info_df$sample_name) %>%
        dplyr::group_split(!!dplyr::sym("group"))
    group_indices_list <- lapply(
        group_indices, function(group) {
            group[["group"]]
        })
    h5_obj <- HDF5Array(met_obj@seed@filepath, name = "X")
    # Serial processing - process entire matrix at once
    
    significant_indices <- lapply(
        seq_along(group_indices_list),
        function(k) {
            current_group_idx <- which(
                sam_info_df[["group"]] %in% 
                group_indices_list[[k]])
            other_groups_idx <- which(
                sam_info_df[["group"]] %in% 
                unlist(group_indices_list[-k]))
            group1_means <- DelayedMatrixStats::rowMeans2(
                h5_obj, cols = current_group_idx, na.rm = TRUE)
            group2_means <- DelayedMatrixStats::rowMeans2(
                h5_obj, cols = other_groups_idx, na.rm = TRUE)
            which(
                abs(group1_means - group2_means) >=
                mean_diff_abs
            )
        }
    )
    return(significant_indices)
}
## calldiff_turbo step-3
loop_w_presto <- function(
    selected_rows=selected_rows,
    BATCH_SIZE=BATCH_SIZE,
    met_obj=met_obj,
    sam_info_df=sam_info_df,
    result_cols=result_cols,
    infile=infile
)
{
    start_idx <- 1
    while (start_idx <= length(selected_rows)) {
        end_idx <- min(start_idx + BATCH_SIZE - 1, length(selected_rows))
        batch_size <- end_idx - start_idx + 1
        batch_rows <- selected_rows[start_idx:end_idx]
        batch_data <- rhdf5::h5read(
          met_obj@seed@filepath, "X", index = list(batch_rows, NULL))
        colnames(batch_data) <- met_obj@sample_name
        rownames(batch_data) <- met_obj@marker_name[batch_rows]
        batch_res <- presto::wilcoxauc(batch_data, y = sam_info_df[["group"]])
        pb <- utils::txtProgressBar(
            min = 0, max = length(result_cols),style = 3)
        for (i in seq_along(result_cols)) {
            col <- result_cols[i]
            current_rows <- dim(rhdf5::h5read(
                infile, paste0("uns/diff_results/", col)))[1]
            if (is.null(current_rows))
                current_rows <- 0
            new_rows <- current_rows + nrow(batch_res)
            h5set_extent(
                file = infile, dataset = paste0("uns/diff_results/", col),
                dims = c(new_rows)
            )
            rhdf5::h5write(
                obj = as.character(batch_res[[col]]),
                file = infile, name = paste0("uns/diff_results/", col),
                index = list((current_rows + 1):new_rows)
            )
            utils::setTxtProgressBar(pb, i)
        }
        close(pb)
        start_idx <- end_idx + 1
    }
}
## calldiff_turbo step-4
pr_result <- function(
    infile=infile,
    open_h5=open_h5,
    met_obj=met_obj,
    sample_batch=sample_batch,
    sam_info_df=sam_info_df,
    BATCH_SIZE=BATCH_SIZE,
    selected_rows=selected_rows
)
{
    ## Create new group
    open_h5 <- rhdf5::H5Fopen(infile)
    rhdf5::h5createGroup(open_h5, "uns/diff_results")
    rhdf5::H5Fclose(open_h5)
    rhdf5::h5closeAll()
    ## Get column names from sample result
    sample_batch <- rhdf5::h5read(
        met_obj@seed@filepath, "X", 
        index=list(seq_len(2), NULL))
    temp_res <- presto::wilcoxauc(sample_batch, y=sam_info_df[["group"]])
    result_cols <- colnames(temp_res)
    ## Create expandable datasets for each column
    for (col in result_cols) {
        h5createDataset(
            file=infile, dataset=paste0("uns/diff_results/", col),
            dims=c(0),
            maxdims=c(H5Sunlimited()),
            chunk=c(BATCH_SIZE),
            storage.mode="character", size = 100)
    }
    loop_w_presto(
        selected_rows=selected_rows,
        BATCH_SIZE=BATCH_SIZE,
        met_obj=met_obj,
        sam_info_df=sam_info_df,
        result_cols=result_cols,
        infile=infile)
    rhdf5::h5write(
      selected_rows, file=infile, name="uns/diff_results_idx")
    message("Successfully saved new results to /uns/diff_results")
    rhdf5::h5closeAll()
    hdf5_5mc <- met_obj@seed@filepath
    diff_results <- as.data.frame(h5read(hdf5_5mc, "uns/diff_results"))
    return(diff_results)
}
