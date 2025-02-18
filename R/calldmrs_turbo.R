#' @title Compute Differentially Methylated Regions (DMRs) Preparation
#' 
#' @description This function processes methylation data, generates plots,
#'              computes statistical differences between case and control
#'              groups, and returns a dataframe of significant DMRs.
#'
#' @details This function performs the following tasks:
#' - Validates input arguments.
#' - Generates methylation statistics plots for the given groups.
#' - Computes statistical differences between the methylation levels of
#'   the case and control groups.
#' - Filters results based on False Discovery Rate (FDR) and returns
#'   a dataframe with non-NA FDR values.
#'
#' @param met_obj A methylTracer object containing methylation data.
#' @param group_colname Column name containing group information
#'                      (default: 'group').
#' @param case_group Name of the case group in the data (default: 'case').
#' @param control_group Name of the control group in the data
#'                      (default: 'ctrl').
#'
#' @return DMLresult result in case and control adjusted in p.adj.
#'         A dataframe containing non-NA FDR values with integer positions.
#'         If no significant results are found, a message is displayed.
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
pre_calldmrs <- function(
    met_obj = NULL, 
    group_colname = NULL, case_group = NULL, control_group = NULL) {
    # Input validation
    if (is.null(met_obj)) {
        message("met_obj cannot be NULL")
    }
    groupname <- "X"
    plot_met_stats(
        met_obj = met_obj, groupname = groupname, 
        diff_group_colname = group_colname,
        case_group = case_group, control_group = control_group
    )
    stat_df <- compute_Stat(
        met_obj = met_obj, case_group = case_group,
        diff_group_colname = group_colname, 
        control_group = control_group
    )
    stat_df_nona <- stat_df[!is.na(stat_df$fdr),] %>%
        mutate(pos = as.integer(.data$pos))
    if (nrow(stat_df_nona) == 0) {
        message("No significant results found after FDR filtering")
    }
    message("Saving results in /uns/dmrs_results")
    infile <- met_obj@seed@filepath
    open_h5 <- rhdf5::H5Fopen(infile)
    has_old_results <- rhdf5::H5Lexists(open_h5, "uns/dmrs_results")
    rhdf5::H5Fclose(open_h5)
    if (has_old_results) {
        rhdf5::h5delete(infile, "uns/dmrs_results")
    }
    rhdf5::h5closeAll()
    rhdf5::h5write(stat_df_nona, file = infile, name = "uns/dmrs_results")
    return(stat_df_nona)
}


compute_Stat <- function(
    met_obj=NULL, diff_group_colname=NULL, 
    case_group=NULL, control_group=NULL) {

    open_h5 <- rhdf5::H5Fopen(met_obj@seed@filepath)
    case_path <- paste0("var/", case_group, "_rowMeans2")
    control_path <- paste0("var/", control_group, "_rowMeans2")
    if (!rhdf5::H5Lexists(open_h5, case_path)) {
        rhdf5::H5Fclose(open_h5)
        message(sprintf("Group data not found for %s or %s", 
                        case_group, control_group))
    }
    rhdf5::H5Fclose(open_h5)
    mean_1 <- h5read(
      met_obj@seed@filepath, paste0("var/", case_group, "_rowMeans2"))
    mean_2 <- h5read(
      met_obj@seed@filepath, paste0("var/", control_group, "_rowMeans2"))
    var1 <- h5read(
      met_obj@seed@filepath, paste0("var/", case_group, "_rowVars"))
    var2 <- h5read(
      met_obj@seed@filepath, paste0("var/", control_group, "_rowVars"))

    # compute statistics
    result <- computeStatCpp(mean_1, mean_2, var1, var2)
    chr <- sub("_.*", "", met_obj@marker_name)
    pos <- sub("^[^_]+_([^_]+)_.*", "\\1", met_obj@marker_name)
    result <- cbind(chr = chr, pos = pos, result)

    result <- com_sta_out(case_group=case_group,
        control_group=control_group,
        result=result,
        met_obj=met_obj)

    return(result)
}

## compute_Stat step1
com_sta_out <- function(
    case_group=case_group,
    control_group=control_group,
    result=result,
    met_obj=met_obj
)
{
    col_1 <- paste0("mean_", case_group)
    col_2 <- paste0("mean_", control_group)

    colnames(result)[colnames(result) == "mu1"] <- col_1
    colnames(result)[colnames(result) == "mu2"] <- col_2
    # Check and create 'var' group
    open_h5 <- rhdf5::H5Fopen(met_obj@seed@filepath)
    if (!rhdf5::H5Lexists(open_h5, "var")) {
        rhdf5::H5Fclose(open_h5)
        rhdf5::h5createGroup(met_obj@seed@filepath, "var")
    }
    rhdf5::H5Fclose(open_h5)
    # Check and delete old results
    open_h5 <- rhdf5::H5Fopen(met_obj@seed@filepath)
    has_old_results <- rhdf5::H5Lexists(open_h5, "var/stat_results")
    rhdf5::H5Fclose(open_h5)
    if (has_old_results) {
        message("Removing old statistical results...")
        rhdf5::h5delete(met_obj@seed@filepath, "var/stat_results")
    }
    rhdf5::h5closeAll()
    # Write new results
    tryCatch(
        {
            h5write(result, met_obj@seed@filepath, "var/stat_results")
            message("Saving results in /var/stat_results")
        }, error = function(e) {
            message(sprintf("Writing results: %s", e$message))
        }
    )
    return(result)
}


#' @title Call DMRS Turbo
#' 
#' @description This function identifies Differentially Methylated Regions 
#'              (DMRs) based on methylation data using the turbo algorithm.
#'              It calls DMRs by evaluating methylation differences between
#'              case and control groups, applying various thresholds to 
#'              filter the results.
#'
#' @param met_obj A methylTracer object containing methylation data.
#' @param p_threshold p-value threshold for calling DMRs (default: 0.05).
#' @param minlen Minimum length of a DMR, in base pairs (default: 500).
#' @param minCG Minimum number of CG sites within a DMR (default: 3).
#' @param dis_merge Distance for merging DMRs (default: 200).
#' @param pct_sig Percentage of significant sites required for a DMR
#'                to be considered significant (default: 0.1).
#' @param sep Window size for separating potential DMRs (default: 5000).
#' @param case_group The name of the case group in the methylation data
#'                   (default: 'case').
#' @param ctrl_group The name of the control group in the methylation data
#'                   (default: 'ctrl').
#' 
#' @return A GRanges object containing the identified DMRs, which can be 
#'         saved or plotted. The GRanges object will include the genomic 
#'         coordinates and methylation statistics of the DMRs.
#'
#' @details This function computes DMRs by first calculating the methylation
#'          difference between case and control groups. It filters the results
#'          based on various thresholds, including p-value, minimum length, 
#'          and the number of CG sites within a DMR. The final DMRs are merged
#'          based on the specified distance threshold and can be visualized 
#'          or saved for further analysis.
#'
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
calldmrs_turbo <- function(
    met_obj=NULL, 
    p_threshold=5e-02, 
    minlen=50L, 
    minCG=3L, 
    dis_merge=100,
    pct_sig=0.1, 
    sep=5000, 
    case_group=NULL, 
    ctrl_group=NULL
) {
    start_time <- Sys.time()
    DMLresult <- rhdf5::h5read(
      met_obj@seed@filepath, "uns/dmrs_results")
    message("============================== \n")
    total_CG_sites <- dim(DMLresult)[1]
    col_1 <- paste0("mean_", case_group)
    col_2 <- paste0("mean_", ctrl_group)
    colnames(DMLresult)[colnames(DMLresult) ==
        col_1] <- "mu1"
    colnames(DMLresult)[colnames(DMLresult) ==
        col_2] <- "mu2"

    message(sprintf("Total CG sites: %d\n", total_CG_sites))
    dmrs_res <- .Call(
        `_methylTracer_calldmrs_turbo`, 
        DMLresult, p_threshold, 
        minlen, minCG, dis_merge,
        pct_sig, sep)
    message(sprintf("DMRs detected: %d\n", dim(dmrs_res)[1]))
    end_time <- Sys.time()  # End timing
    elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    message(sprintf("Elapsed time: %.2f seconds\n", elapsed_time))
    message("==============================")
    colnames(dmrs_res)[colnames(dmrs_res) ==
        "meanMethy1"] <- col_1
    colnames(dmrs_res)[colnames(dmrs_res) ==
        "meanMethy2"] <- col_2
    dmrs_gr <- GenomicRanges::GRanges(
        seqnames=dmrs_res$chr, 
        ranges=IRanges::IRanges(start=dmrs_res$start, end=dmrs_res$end),
        strand="*", nCG=dmrs_res$nCG, 
        md=dmrs_res$diff.Methy, dmrs_length=dmrs_res$length,
        meanMethy1= dmrs_res[[col_1]], 
        meanMethy2=dmrs_res[[col_2]], 
        areaStat=dmrs_res$areaStat
    )
    return(dmrs_gr)
}
