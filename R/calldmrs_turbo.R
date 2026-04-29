#' @title Prepare differential methylation analysis (DMR pre-computation)
#'
#' @description
#' Compute region-level differential methylation statistics between two
#' groups (case vs control) based on a \code{methylTracer} object. The
#' function ensures that group-wise methylation summaries exist, calls
#' the statistical C++ backend, filters by FDR, and writes the results
#' back into the HDF5 file.
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Validate that \code{met} and group labels are provided.
#'   \item Compute group-wise methylation summaries (means, variances)
#'         via \code{calc_met_stats()} if not already present.
#'   \item Call \code{compute_Stat()} to obtain per-region test
#'         statistics, p-values and FDR.
#'   \item Filter out rows with missing FDR, coerce genomic positions
#'         to integer, and (if any rows remain) write the resulting
#'         table into the HDF5 file under \code{"/uns/dmrs_results"}.
#'   \item In addition, \code{compute_Stat()} writes a more complete
#'         per-region summary table to \code{"/var/stat_results"}.
#' }
#'
#' The function assumes that methylation values in \code{met} are
#' stored as proportions in the range \eqn{[0, 1]}.
#'
#' @param met A \code{methylTracer} object containing region-level
#'   methylation data and a reference to the underlying HDF5 file.
#' @param group_colname Character scalar. Column name in the sample
#'   metadata that encodes the group labels (default: \code{"group"}).
#' @param case_group Character scalar. Label of the case group within
#'   \code{group_colname} (default: \code{"case"}).
#' @param control_group Character scalar. Label of the control group
#'   within \code{group_colname} (default: \code{"ctrl"}).
#'
#' @return
#' A \code{data.frame} containing rows with non-\code{NA} FDR values.
#' Columns include chromosome, genomic position, group-wise means and
#' variances, test statistics, p-values, and FDR (exact column names
#' depend on the C++ backend). If no significant results remain after
#' filtering (i.e. the returned data frame has zero rows), a message
#' is displayed and an empty data frame is returned.
#'
#' As a side effect, the function overwrites (if present) the following
#' datasets in the underlying HDF5 file:
#' \itemize{
#'   \item \code{"/uns/dmrs_results"} – filtered DMR table.
#'   \item \code{"/var/stat_results"} – full per-region statistics
#'         table produced by \code{compute_Stat()}.
#' }
#'
#' @export
#'
#' @examples
#' ## Assuming 'met' is a methylTracer object created with build_met_obj()
#' ## and that 'group' column contains "case" and "ctrl" labels:
#' # pre_res <- pre_calldmrs(
#' #   met            = met,
#' #   group_colname  = "group",
#' #   case_group     = "case",
#' #   control_group  = "ctrl"
#' # )
pre_calldmrs <- function(
  met,
  group_colname = "group",
  case_group = "case",
  control_group = "ctrl"
) {
  ## basic checks
  if (missing(met) || is.null(met)) {
    stop("'met' must be a non-null methylTracer object.")
  }

  infile <- met@seed@filepath
  if (!file.exists(infile)) {
    stop("HDF5 file referenced by 'met' does not exist: ", infile)
  }

  groupname <- "X" ## dataset name used internally by methylTracer

  ## 1. Compute and store per-group means / vars in /var/...
  calc_met_stats(
    met = met,
    groupname = groupname,
    diff_group_colname = group_colname,
    case_group = case_group,
    control_group = control_group
  )

  ## 2. Compute per-marker statistics (mu1, mu2, diff, se, stat, p, fdr)
  stat_df <- compute_Stat(
    met = met,
    case_group = case_group,
    diff_group_colname = group_colname,
    control_group = control_group
  )

  ## Filter out rows with NA FDR
  stat_df_nona <- stat_df[!is.na(stat_df$fdr), , drop = FALSE]
  stat_df_nona <- dplyr::mutate(stat_df_nona, pos = as.integer(.data$pos))

  if (nrow(stat_df_nona) == 0L) {
    message(
      "No significant results found after FDR filtering; ",
      "nothing will be written to /uns/dmrs_results."
    )
    return(stat_df_nona)
  }

  ## in case compute_Stat() failed gracefully
  # if (is.null(stat_df)) {
  #   message("No statistics computed; returning NULL.")
  #   return(invisible(NULL))
  # }

  ## filter out NA FDR rows and coerce position to integer
  # stat_df_nona <- stat_df[!is.na(stat_df$fdr), , drop = FALSE]
  # if (!"pos" %in% colnames(stat_df_nona)) {
  #   stop("Expected column 'pos' in stat_df, but not found.")
  # }
  # stat_df_nona$pos <- as.integer(stat_df_nona$pos)
  #
  # if (nrow(stat_df_nona) == 0) {
  #   message("No significant results found after FDR filtering.")
  #   return(stat_df_nona)
  # }

  ## 3. Overwrite /uns/dmrs_results with filtered results -----------------
  message("Saving filtered results in /uns/dmrs_results")

  # Ensure /uns group exists
  ls_df <- rhdf5::h5ls(infile)
  has_uns_group <- any(
    ls_df$group == "/" & ls_df$name == "uns" & ls_df$otype == "H5I_GROUP"
  )
  if (!has_uns_group) {
    rhdf5::h5createGroup(infile, "uns")
  }

  # Delete any previous dmrs_results
  suppressWarnings(
    try(rhdf5::h5delete(infile, "uns/dmrs_results"), silent = TRUE)
  )

  # Write as character matrix for consistency
  stat_mat <- as.matrix(stat_df_nona)
  storage.mode(stat_mat) <- "character"

  rhdf5::h5write(
    obj = stat_mat,
    file = infile,
    name = "uns/dmrs_results"
  )

  rhdf5::h5closeAll()

  return(stat_df_nona)
}


# compute_Stat <- function(
#     met=NULL, diff_group_colname=NULL,
#     case_group=NULL, control_group=NULL) {
#
#     open_h5 <- rhdf5::H5Fopen(met@seed@filepath)
#     case_path <- paste0("var/", case_group, "_rowMeans2")
#     control_path <- paste0("var/", control_group, "_rowMeans2")
#     if (!rhdf5::H5Lexists(open_h5, case_path)) {
#         rhdf5::H5Fclose(open_h5)
#         message(sprintf("Group data not found for %s or %s",
#                         case_group, control_group))
#     }
#     rhdf5::H5Fclose(open_h5)
#     rhdf5::h5closeAll()
#     mean_1 <- rhdf5::h5read(
#       met@seed@filepath, paste0("var/", case_group, "_rowMeans2"))
#     mean_2 <- rhdf5::h5read(
#       met@seed@filepath, paste0("var/", control_group, "_rowMeans2"))
#     var1 <- rhdf5::h5read(
#       met@seed@filepath, paste0("var/", case_group, "_rowVars"))
#     var2 <- rhdf5::h5read(
#       met@seed@filepath, paste0("var/", control_group, "_rowVars"))
#     nn <- rhdf5::h5read(
#       met@seed@filepath, paste0("obs/", diff_group_colname)
#     )
#     nn_1 <- sum(nn == case_group)
#     nn_2 <- sum(nn == control_group)
#
#     ## compute statistics
#     result <- computeStatCpp(mean_1, mean_2, var1, var2, nn_1, nn_2)
#     split_result <- strsplit(as.character(met@marker_name), "_")
#     chr <- sapply(split_result, function(x) x[1])
#     pos <- sapply(split_result, function(x) x[3])
#     result <- cbind(chr = chr, pos = pos, result)
#     result <- com_sta_out(case_group=case_group,
#         control_group=control_group,
#         result=result,
#         met=met)
#     return(result)
# }

## compute_Stat step1
# com_sta_out <- function(
#     case_group=case_group,
#     control_group=control_group,
#     result=result,
#     met=met
# )
# {
#     col_1 <- paste0("mean_", case_group)
#     col_2 <- paste0("mean_", control_group)
#
#     colnames(result)[colnames(result) == "mu1"] <- col_1
#     colnames(result)[colnames(result) == "mu2"] <- col_2
#     # Check and create 'var' group
#     open_h5 <- rhdf5::H5Fopen(met@seed@filepath)
#     if (!rhdf5::H5Lexists(open_h5, "var")) {
#         rhdf5::H5Fclose(open_h5)
#         rhdf5::h5createGroup(met@seed@filepath, "var")
#     }
#     rhdf5::H5Fclose(open_h5)
#     # Check and delete old results
#     open_h5 <- rhdf5::H5Fopen(met@seed@filepath)
#     has_old_results <- rhdf5::H5Lexists(open_h5, "var/stat_results")
#     rhdf5::H5Fclose(open_h5)
#     if (has_old_results) {
#         message("Removing old statistical results...")
#         rhdf5::h5delete(met@seed@filepath, "var/stat_results")
#     }
#     rhdf5::h5closeAll()
#     # Write new results
#     tryCatch(
#         {
#             h5write(result, met@seed@filepath, "var/stat_results")
#             message("Saving results in /var/stat_results")
#         }, error = function(e) {
#             message(sprintf("Writing results: %s", e$message))
#         }
#     )
#     return(result)
# }

compute_Stat <- function(
  met,
  diff_group_colname,
  case_group,
  control_group
) {
  infile <- met@seed@filepath

  # ============================================================
  # CRITICAL FIX: Use on.exit() to guarantee handle cleanup
  # This prevents file handle leaks even when errors occur
  # ============================================================
  h5_handle <- NULL

  tryCatch(
    {
      # Open handle with guaranteed cleanup
      h5_handle <- rhdf5::H5Fopen(infile, flags = "H5F_ACC_RDONLY")
      on.exit(
        {
          if (!is.null(h5_handle)) {
            rhdf5::H5Fclose(h5_handle)
          }
        },
        add = TRUE
      )

      case_path <- paste0("var/", case_group, "_rowMeans2")
      control_path <- paste0("var/", control_group, "_rowMeans2")

      if (
        !rhdf5::H5Lexists(h5_handle, case_path) ||
          !rhdf5::H5Lexists(h5_handle, control_path)
      ) {
        stop(
          sprintf(
            "Group data not found for '%s' or '%s'; please run calc_met_stats() first.",
            case_group,
            control_group
          ),
          call. = FALSE
        )
      }

      # Close handle before reading (h5read manages its own handles)
      rhdf5::H5Fclose(h5_handle)
      h5_handle <- NULL # Mark as closed

      ## read per-region means and variances with safe wrapper
      mean_1 <- tryCatch(
        {
          rhdf5::h5read(infile, paste0("var/", case_group, "_rowMeans2"))
        },
        finally = {
          rhdf5::h5closeAll()
        }
      )

      mean_2 <- tryCatch(
        {
          rhdf5::h5read(infile, paste0("var/", control_group, "_rowMeans2"))
        },
        finally = {
          rhdf5::h5closeAll()
        }
      )

      var1 <- tryCatch(
        {
          rhdf5::h5read(infile, paste0("var/", case_group, "_rowVars"))
        },
        finally = {
          rhdf5::h5closeAll()
        }
      )

      var2 <- tryCatch(
        {
          rhdf5::h5read(infile, paste0("var/", control_group, "_rowVars"))
        },
        finally = {
          rhdf5::h5closeAll()
        }
      )

      ## sample counts per group
      nn <- tryCatch(
        {
          rhdf5::h5read(infile, paste0("obs/", diff_group_colname))
        },
        finally = {
          rhdf5::h5closeAll()
        }
      )

      nn_1 <- sum(nn == case_group)
      nn_2 <- sum(nn == control_group)

      ## compute statistics in C++ backend
      result <- computeStatCpp(mean_1, mean_2, var1, var2, nn_1, nn_2)

      ## add genomic coordinates derived from marker_name
      split_result <- strsplit(as.character(met@marker_name), "_")
      chr <- vapply(split_result, function(x) x[1], FUN.VALUE = character(1L))
      pos <- vapply(split_result, function(x) x[3], FUN.VALUE = character(1L))

      result <- cbind(chr = chr, pos = pos, result)

      ## post-process + write to HDF5
      result <- com_sta_out(
        case_group = case_group,
        control_group = control_group,
        result = result,
        met = met
      )

      result
    },
    error = function(e) {
      # Ensure cleanup on error
      if (!is.null(h5_handle)) {
        rhdf5::H5Fclose(h5_handle)
      }
      rhdf5::h5closeAll()
      stop("Error in compute_Stat: ", e$message, call. = FALSE)
    }
  )
}

com_sta_out <- function(
  case_group,
  control_group,
  result,
  met
) {
  col_1 <- paste0("mean_", case_group)
  col_2 <- paste0("mean_", control_group)

  colnames(result)[colnames(result) == "mu1"] <- col_1
  colnames(result)[colnames(result) == "mu2"] <- col_2

  infile <- met@seed@filepath

  # ============================================================
  # CRITICAL FIX: Proper handle management with on.exit()
  # ============================================================
  h5_handle <- NULL

  tryCatch(
    {
      ## ensure 'var' group exists
      h5_handle <- rhdf5::H5Fopen(infile, flags = "H5F_ACC_RDWR")
      on.exit(
        {
          if (!is.null(h5_handle)) {
            rhdf5::H5Fclose(h5_handle)
          }
        },
        add = TRUE
      )

      if (!rhdf5::H5Lexists(h5_handle, "var")) {
        rhdf5::H5Fclose(h5_handle)
        h5_handle <- NULL
        rhdf5::h5createGroup(infile, "var")
      } else {
        rhdf5::H5Fclose(h5_handle)
        h5_handle <- NULL
      }

      ## delete any old stat_results
      h5_handle <- rhdf5::H5Fopen(infile, flags = "H5F_ACC_RDONLY")
      has_old_results <- rhdf5::H5Lexists(h5_handle, "var/stat_results")
      rhdf5::H5Fclose(h5_handle)
      h5_handle <- NULL

      if (has_old_results) {
        message("Removing old statistical results in /var/stat_results ...")
        rhdf5::h5delete(infile, "var/stat_results")
      }

      ## write new stat_results with guaranteed cleanup
      tryCatch(
        {
          rhdf5::h5write(result, infile, "var/stat_results")
          message("Saving results in /var/stat_results")
        },
        finally = {
          rhdf5::h5closeAll()
        }
      )

      result
    },
    error = function(e) {
      if (!is.null(h5_handle)) {
        rhdf5::H5Fclose(h5_handle)
      }
      rhdf5::h5closeAll()
      stop("Error in com_sta_out: ", e$message, call. = FALSE)
    }
  )
}

#' @title Call differentially methylated regions (DMRs) using the turbo algorithm
#'
#' @description
#' Identify differentially methylated regions (DMRs) from a
#' \code{methylTracer} object using the C++ \emph{turbo} DMR-calling
#' backend. The function uses per-CpG test statistics and p-values
#' stored in the HDF5 file (under \code{"/uns/dmrs_results"}) to
#' detect contiguous regions enriched for significant CpG sites.
#'
#' @details
#' This function is intended to be run \strong{after}
#' \code{\link{pre_calldmrs}}, which computes per-CpG statistics
#' (means, variances, t-statistics, p-values, FDR) and stores them as a
#' table in \code{"/uns/dmrs_results"} in the underlying HDF5 file.
#'
#' Internally, the function:
#' \enumerate{
#'   \item Reads \code{"/uns/dmrs_results"} from the HDF5 file attached
#'         to \code{met}.
#'   \item Converts the table to a data frame and assigns standard
#'         column names (\code{chr}, \code{pos}, \code{mu1},
#'         \code{mu2}, \code{diff}, \code{diff.se}, \code{stat},
#'         \code{pval}, \code{fdr}).
#'   \item Calls the C++ routine \code{calldmrs_turbo} to group CpGs
#'         into candidate regions based on the p-value threshold and
#'         genomic distance rules, and to compute region-level
#'         statistics:
#'         mean methylation in each group, methylation difference,
#'         and area statistic (sum of per-CpG test statistics).
#'   \item Returns a \link[GenomicRanges]{GRanges} object with DMR
#'         coordinates and associated statistics.
#' }
#'
#' Methylation values are assumed to be stored as proportions in
#' \eqn{[0, 1]}. Region-level mean methylation columns in the return
#' object are renamed to \code{mean_<case_group>} and
#' \code{mean_<ctrl_group>} for convenience.
#'
#' @param met A \code{methylTracer} object containing region-level
#'   methylation data and a reference to the underlying HDF5 file.
#'   The object must have been processed with \code{\link{pre_calldmrs}}
#'   beforehand.
#' @param p_threshold Numeric scalar. Per-CpG p-value threshold used
#'   by the turbo algorithm to define “significant” CpGs
#'   (default: \code{0.05}).
#' @param minlen Integer scalar. Minimum DMR length in base pairs
#'   (default: \code{50}).
#' @param minCG Integer scalar. Minimum number of CpG sites required
#'   for a region to be considered a DMR (default: \code{3}).
#' @param dis_merge Numeric scalar. Maximum distance (in base pairs)
#'   between adjacent significant segments within a region to be
#'   merged into a single DMR (default: \code{100}).
#' @param pct_sig Numeric scalar between 0 and 1. Minimum fraction of
#'   CpG sites within a candidate region that must pass the
#'   \code{p_threshold} to retain the region (default: \code{0.1}).
#' @param sep Numeric scalar. Genomic distance (in base pairs) used to
#'   split the genome into independent windows for bump finding
#'   (default: \code{5000}).
#' @param case_group Character scalar. Label of the case group in the
#'   sample metadata (default: \code{"case"}). Used only to rename
#'   region-level mean methylation columns in the output.
#' @param ctrl_group Character scalar. Label of the control group in
#'   the sample metadata (default: \code{"ctrl"}).
#'
#' @return
#' A \link[GenomicRanges]{GRanges} object with one row per DMR. The
#' metadata columns include:
#' \itemize{
#'   \item \code{nCG}: number of CpG sites in the region.
#'   \item \code{md}: mean methylation difference
#'         (\code{mean_<case_group> - mean_<ctrl_group>}).
#'   \item \code{dmrs_length}: region length in base pairs.
#'   \item \code{mean_<case_group>}: mean methylation proportion (0–1)
#'         in the case group across the region.
#'   \item \code{mean_<ctrl_group>}: mean methylation proportion (0–1)
#'         in the control group across the region.
#'   \item \code{areaStat}: sum of per-CpG test statistics within the
#'         region, reflecting both effect size and consistency.
#' }
#' If no DMRs are detected, an empty \code{GRanges} object is returned.
#'
#' @seealso
#' \code{\link{pre_calldmrs}} for generating the per-CpG statistics
#' that this function consumes.
#'
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
#' ## toy methylation matrix (proportions 0–1)
#' input_file_df <- data.frame(
#'   marker_name = c(
#'     "chr1_1000_2000", "chr1_2000_3000",
#'     "chr1_3000_4000", "chr1_4000_5000"
#'   ),
#'   sample_1    = c(0.10, 0.20, 0.60, 0.90),
#'   sample_2    = c(0.10, 0.20, 0.60, 0.90),
#'   sample_3    = c(0.10, 0.20, 0.60, 0.90),
#'   sample_4    = c(0.10, 0.20, 0.60, 0.90),
#'   sample_5    = c(0.80, 0.90, 0.10, 0.05),
#'   sample_6    = c(0.80, 0.90, 0.10, 0.05),
#'   sample_7    = c(0.80, 0.90, 0.10, 0.05),
#'   sample_8    = c(0.80, 0.90, 0.10, 0.05)
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
#' pre_res <- pre_calldmrs(
#'   met           = met,
#'   group_colname = "group",
#'   case_group    = "case",
#'   control_group = "ctrl"
#' )
#'
#' dmr_res <- calldmrs_turbo(
#'   met         = met,
#'   p_threshold = 0.05,
#'   case_group  = "case",
#'   ctrl_group  = "ctrl"
#' )
#' }
calldmrs_turbo <- function(
  met = NULL,
  p_threshold = 5e-02,
  minlen = 50L,
  minCG = 3L,
  dis_merge = 100,
  pct_sig = 0.1,
  sep = 5000,
  case_group = NULL,
  ctrl_group = NULL
) {
  if (missing(met) || is.null(met)) {
    stop("'met' cannot be NULL. Please provide a methylTracer object.")
  }
  if (is.null(case_group) || is.null(ctrl_group)) {
    stop("'case_group' and 'ctrl_group' must be provided.")
  }

  infile <- met@seed@filepath

  ## check that /uns/dmrs_results exists
  h5f <- rhdf5::H5Fopen(infile)
  has_dmrs <- rhdf5::H5Lexists(h5f, "uns/dmrs_results")
  rhdf5::H5Fclose(h5f)
  rhdf5::h5closeAll()

  if (!has_dmrs) {
    stop(
      "Dataset '/uns/dmrs_results' not found in HDF5 file. ",
      "Please run pre_calldmrs() first."
    )
  }

  start_time <- Sys.time()
  ## read per-CpG stats table (written by pre_calldmrs())
  DMLresult <- rhdf5::h5read(
    met@seed@filepath,
    "uns/dmrs_results"
  )
  message("\n")
  total_CG_sites <- dim(DMLresult)[1]

  ## assign column names (must match computeStatCpp output)
  colN <- c(
    "chr",
    "pos",
    "mu1",
    "mu2",
    "diff",
    "diff.se",
    "stat",
    "pval",
    "fdr"
  )
  if (ncol(DMLresult) != length(colN)) {
    stop(
      "Unexpected number of columns in '/uns/dmrs_results'. ",
      "Got ",
      ncol(DMLresult),
      ", expected ",
      length(colN),
      "."
    )
  }
  colnames(DMLresult) <- colN

  ## convert types
  DMLresult <- as.data.frame(DMLresult)
  DMLresult[, 1] <- as.character(DMLresult[, 1])
  DMLresult[, 2] <- as.integer(DMLresult[, 2])
  DMLresult[, 3:9] <- lapply(DMLresult[, 3:9], as.numeric)

  ## call C++ turbo DMR caller
  message(sprintf("Total CG sites: %d\n", total_CG_sites))
  dmrs_res <- .Call(
    `_methylTracer_calldmrs_turbo`,
    DMLresult,
    p_threshold,
    minlen,
    minCG,
    dis_merge,
    pct_sig,
    sep
  )
  message(sprintf("DMRs detected: %d\n", dim(dmrs_res)[1]))
  end_time <- Sys.time() # End timing
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  message(sprintf("Elapsed time: %.2f seconds\n", elapsed_time))

  ## rename region-level mean methylation columns
  col_1 <- paste0("mean_", case_group)
  col_2 <- paste0("mean_", ctrl_group)
  colnames(dmrs_res)[colnames(dmrs_res) == "meanMethy1"] <- col_1
  colnames(dmrs_res)[colnames(dmrs_res) == "meanMethy2"] <- col_2

  ## construct GRanges
  dmrs_gr <- GenomicRanges::GRanges(
    seqnames = dmrs_res$chr,
    ranges = IRanges::IRanges(start = dmrs_res$start, end = dmrs_res$end),
    strand = "*",
    nCG = dmrs_res$nCG,
    md = dmrs_res$diff.Methy,
    dmrs_length = dmrs_res$length,
    meanMethy1 = dmrs_res[[col_1]],
    meanMethy2 = dmrs_res[[col_2]],
    areaStat = dmrs_res$areaStat
  )
  return(dmrs_gr)
}
