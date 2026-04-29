#' @title Aggregate DMR-level methylation per cell (in-memory) and build SingleCellExperiment
#'
#' @description
#' This function computes, for each DMR and each cell, the mean methylation
#' signal across all windows (markers) overlapping that DMR, by first loading
#' all relevant markers into memory as a dense matrix. This is typically much
#' faster than repeatedly reading small slices from an HDF5-backed matrix,
#' at the cost of higher memory usage. The result is returned as a
#' \code{SingleCellExperiment} object with DMRs as rows and cells as columns.
#'
#' @param met A \code{methylTracer} object containing the per-window
#'   methylation matrix \code{X}. Rows correspond to genomic windows (markers),
#'   columns to cells, and the matrix is stored in an HDF5-backed file.
#' @param dmrs A \code{GRanges} object with DMR coordinates.
#' @param assay_name Character; name of the assay to store DMR-level
#'   methylation values. Default \code{"meth"}.
#' @param group_col Optional character. Name of the column in the HDF5
#'   \code{/obs} group that stores group/cluster labels for cells
#'   (e.g. \code{"Group"}, \code{"cluster"}, \code{"leiden"}). If provided,
#'   this is added to \code{colData}.
#' @param scale_factor Numeric or \code{NULL}. use the default
#'   \code{NULL} (no rescaling).
#' @param verbose Logical; whether to print progress messages. Default TRUE.
#'
#' @details
#' Marker coordinates are inferred from \code{met@marker_name}, assuming
#' they are of the form \code{"chr_start_end"}. Overlaps between markers
#' and DMRs are computed with \code{findOverlaps()}, and only markers that
#' overlap at least one DMR are loaded into memory. For each DMR, the mean
#' methylation per cell is computed as the mean over all overlapping markers,
#' ignoring \code{NA} values.
#'
#' @return
#' A \code{SingleCellExperiment} object with:
#' \itemize{
#'   \item rows = DMRs (\code{rowRanges} set to \code{dmrs});
#'   \item columns = cells (from \code{met@sample_name});
#'   \item assay \code{assay_name} = matrix of DMR × cell mean methylation.
#' }
#'
#' @importFrom GenomicRanges GRanges findOverlaps seqnames
#' @importFrom IRanges IRanges
#' @importFrom HDF5Array HDF5Array
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame queryHits subjectHits
#' @importFrom BiocGenerics start end
#' @importFrom stats setNames
#' @export
dmrs_to_sce_memory <- function(
    met,
    dmrs,
    assay_name   = "meth",
    group_col    = NULL,
    scale_factor = NULL,
    verbose      = TRUE
) {
  ## 1. Sanity checks
  if (is.null(met)) {
    stop("'met' must be a methylTracer object.")
  }
  if (is.null(dmrs)) {
    stop("'dmrs' must be a GRanges object.")
  }
  
  ## 2. Build marker GRanges from met@marker_name, assuming "chr_start_end"
  marker_names <- as.character(met@marker_name)
  sp <- strsplit(marker_names, "_", fixed = TRUE)
  
  chr   <- vapply(sp, `[`, character(1L), 1L)
  start <- as.integer(vapply(sp, `[`, character(1L), 2L))
  end   <- as.integer(vapply(sp, `[`, character(1L), 3L))
  
  marker_gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges   = IRanges::IRanges(start = start, end = end)
  )
  
  ## 3. Find overlaps: markers (query) vs DMRs (subject)
  hits <- GenomicRanges::findOverlaps(marker_gr, dmrs, type = "any")
  if (length(hits) == 0L) {
    stop("No overlaps found between markers and DMRs.")
  }
  
  marker_idx <- S4Vectors::queryHits(hits)   # marker indices in original X
  dmr_idx    <- S4Vectors::subjectHits(hits) # DMR indices
  
  ## Markers that are involved in at least one DMR
  markers_use <- sort(unique(marker_idx))
  
  if (verbose) {
    message("Number of markers overlapping at least one DMR: ",
            length(markers_use))
  }
  
  ## 4. Read the relevant markers into memory as a dense matrix
  h5_file <- met@seed@filepath
  X_h5    <- HDF5Array::HDF5Array(h5_file, name = "X")
  
  if (verbose) {
    message("Reading subset of X into memory (",
            length(markers_use), " markers x ",
            ncol(X_h5), " cells)...")
  }
  
  X_sub <- as.matrix(X_h5[markers_use, , drop = FALSE])
  
  n_cells <- ncol(X_sub)
  n_dmr   <- length(dmrs)
  cell_names <- as.character(met@sample_name)
  colnames(X_sub) <- cell_names
  
  ## 5. Remap marker indices to positions in X_sub
  # marker_idx is in original X; we need positions in X_sub
  pos_in_sub <- match(marker_idx, markers_use)  # same length as marker_idx
  # split by DMR index: for each DMR, which rows in X_sub to use
  dmr_to_rows_sub <- split(pos_in_sub, dmr_idx)
  
  ## 6. Initialize DMR x cell matrix
  dmr_mat <- matrix(NA_real_, nrow = n_dmr, ncol = n_cells)
  colnames(dmr_mat) <- cell_names
  
  if (!is.null(names(dmrs)) && any(nzchar(names(dmrs)))) {
    rownames(dmr_mat) <- names(dmrs)
  } else {
    rownames(dmr_mat) <- paste0(
      as.character(GenomicRanges::seqnames(dmrs)),
      ":",
      BiocGenerics::start(dmrs),
      "-",
      BiocGenerics::end(dmrs)
    )
  }
  
  ## 7. For each DMR, compute mean methylation per cell in-memory
  if (verbose) {
    message("Computing DMR x cell mean methylation in memory...")
  }
  
  for (i in seq_len(n_dmr)) {
    rows_i <- dmr_to_rows_sub[[as.character(i)]]
    if (is.null(rows_i) || length(rows_i) == 0L) {
      next
    }
    sub_mat <- X_sub[rows_i, , drop = FALSE]
    # use base colMeans; X_sub is already a dense matrix
    dmr_mat[i, ] <- colMeans(sub_mat, na.rm = TRUE)
  }
  
  ## 8. Optional rescaling
  if (!is.null(scale_factor)) {
    dmr_mat <- dmr_mat / scale_factor
  }
  
  ## 9. Build colData from /obs (all columns)
  obs_list <- rhdf5::h5read(h5_file, "obs")
  
  col_df <- S4Vectors::DataFrame(obs_list, check.names = FALSE)
  
  # set rownames from sample_name if available
  if ("sample_name" %in% colnames(col_df)) {
    rownames(col_df) <- as.character(col_df$sample_name)
  } else {
    rownames(col_df) <- cell_names
  }
  
  # reorder rows to match met@sample_name
  col_df <- col_df[cell_names, , drop = FALSE]
  
  # add a convenience "cell" column if not present
  if (!"cell" %in% colnames(col_df)) {
    col_df$cell <- cell_names
  }
  
  ## 10. Construct SingleCellExperiment
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays    = setNames(list(dmr_mat), assay_name),
    rowRanges = dmrs,
    colData   = col_df
  )
  
  if (verbose) {
    message("Done. Created SingleCellExperiment with ",
            nrow(sce), " DMRs x ", ncol(sce), " cells.")
  }
  
  sce
}

