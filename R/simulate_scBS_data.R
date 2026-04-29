#' Simulate Single-Cell Bisulfite Sequencing (scBS) Data
#'
#' @description
#' Generate synthetic single-cell bisulfite sequencing data with realistic
#' sparsity patterns and implanted Differentially Methylated Regions (DMRs)
#' for testing and benchmarking methylation analysis pipelines.
#'
#' @details
#' This function simulates scBS data by:
#' \enumerate{
#'   \item Creating a genomic "universe" of CpG sites
#'   \item Defining biological groups (Control vs. Treatment)
#'   \item Implanting DMRs with known methylation differences
#'   \item Adding realistic coverage sparsity (80-95\%)
#'   \item Writing individual cell files to disk
#' }
#'
#' The simulation models key features of real scBS data:
#' \itemize{
#'   \item High sparsity (most sites have no coverage in most cells)
#'   \item Variable coverage depth (Poisson-distributed)
#'   \item Binomial methylation counts
#'   \item Group-specific DMRs for validation
#' }
#'
#' @param n_cells Integer. Number of cells to simulate (default: 100).
#'   Typically 50-1000 for testing, 5000+ for realistic scenarios.
#' @param n_sites Integer. Total number of CpG sites in the "genome"
#'   (default: 1,000,000). Use 1M for testing, 28M for human genome scale.
#' @param sparsity Numeric in [0, 1]. Probability that a site has NO
#'   coverage in a cell (default: 0.85). Typical scBS: 0.80-0.95.
#' @param n_dmrs Integer. Number of DMRs to implant (default: 100).
#' @param dmr_width Integer. Width of each DMR in CpG sites (default: 20).
#' @param group_ratio Numeric. Proportion of cells in control group
#'   (default: 0.5). E.g., 0.5 = 50\% control, 50\% treatment.
#' @param output_dir Character. Directory to save individual cell files
#'   (default: tempdir()). Will be created if it doesn't exist.
#' @param output_format Character. File format: "cov" (4-column coverage)
#'   or "bed" (6-column BED-like) (default: "cov").
#' @param mean_coverage Numeric. Mean sequencing depth per covered site
#'   (default: 10). Higher = more reads.
#' @param background_meth Numeric in [0, 1]. Baseline methylation level
#'   for non-DMR sites (default: 0.5).
#' @param dmr_meth_ctrl Numeric in [0, 1]. Methylation level in DMRs for
#'   control group (default: 0.8).
#' @param dmr_meth_treat Numeric in [0, 1]. Methylation level in DMRs for
#'   treatment group (default: 0.2).
#' @param seed Integer. Random seed for reproducibility (default: 12345).
#' @param verbose Logical. Print progress messages (default: TRUE).
#'
#' @return A list with components:
#' \describe{
#'   \item{metadata}{data.frame with columns: cell_id, file_path, group,
#'     n_covered_sites, mean_coverage}
#'   \item{dmr_regions}{data.frame defining the implanted DMRs: chr, start,
#'     end, dmr_id, true_meth_ctrl, true_meth_treat}
#'   \item{cpg_universe}{data.frame of all CpG sites: chr, pos, site_id}
#'   \item{output_dir}{Path to directory containing simulated files}
#' }
#'
#' @importFrom data.table data.table setkey := fwrite setorder
#' @export
#'
#' @examples
#' \donttest{
#' # Quick test with 20 cells, 10K sites
#' sim <- simulate_scBS_data(
#'   n_cells = 20,
#'   n_sites = 10000,
#'   sparsity = 0.85,
#'   n_dmrs = 10,
#'   output_dir = tempfile("scBS_test")
#' )
#'
#' # Inspect metadata
#' head(sim$metadata)
#'
#' # Check DMR locations
#' head(sim$dmr_regions)
#'
#' # Read one cell file
#' cell_data <- data.table::fread(sim$metadata$file_path[1])
#' head(cell_data)
#' }
simulate_scBS_data <- function(
  n_cells = 100L,
  n_sites = 1000000L,
  sparsity = 0.85,
  n_dmrs = 100L,
  dmr_width = 20L,
  group_ratio = 0.5,
  output_dir = tempdir(),
  output_format = c("cov", "bed"),
  mean_coverage = 10,
  background_meth = 0.5,
  dmr_meth_ctrl = 0.8,
  dmr_meth_treat = 0.2,
  seed = 12345L,
  verbose = TRUE
) {
  # ========================================================================
  # Input Validation
  # ========================================================================

  if (!is.numeric(n_cells) || n_cells < 2) {
    stop("'n_cells' must be an integer >= 2", call. = FALSE)
  }
  if (!is.numeric(n_sites) || n_sites < 1000) {
    stop("'n_sites' must be an integer >= 1000", call. = FALSE)
  }
  if (!is.numeric(sparsity) || sparsity < 0 || sparsity > 1) {
    stop("'sparsity' must be in [0, 1]", call. = FALSE)
  }
  if (!is.numeric(n_dmrs) || n_dmrs < 1) {
    stop("'n_dmrs' must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(group_ratio) || group_ratio <= 0 || group_ratio >= 1) {
    stop("'group_ratio' must be in (0, 1)", call. = FALSE)
  }

  output_format <- match.arg(output_format)

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) message("Created output directory: ", output_dir)
  }

  set.seed(seed)

  # ========================================================================
  # Step 1: Create CpG Universe (Genomic Sites)
  # ========================================================================

  if (verbose) {
    message("[Step 1/5] Creating CpG universe (", n_sites, " sites)...")
  }

  cpg_universe <- create_cpg_universe(
    n_sites = n_sites,
    verbose = verbose
  )

  # ========================================================================
  # Step 2: Define Biological Groups
  # ========================================================================

  if (verbose) {
    message("[Step 2/5] Assigning cells to groups...")
  }

  n_ctrl <- round(n_cells * group_ratio)
  n_treat <- n_cells - n_ctrl

  cell_groups <- c(
    rep("Control", n_ctrl),
    rep("Treatment", n_treat)
  )

  # Shuffle to avoid batch effects
  cell_groups <- sample(cell_groups)

  cell_metadata <- data.table::data.table(
    cell_id = sprintf("cell_%03d", seq_len(n_cells)),
    group = cell_groups
  )

  if (verbose) {
    message("  - Control: ", n_ctrl, " cells")
    message("  - Treatment: ", n_treat, " cells")
  }

  # ========================================================================
  # Step 3: Implant DMRs (Differentially Methylated Regions)
  # ========================================================================

  if (verbose) {
    message("[Step 3/5] Implanting ", n_dmrs, " DMRs...")
  }

  dmr_regions <- implant_dmrs(
    cpg_universe = cpg_universe,
    n_dmrs = n_dmrs,
    dmr_width = dmr_width,
    dmr_meth_ctrl = dmr_meth_ctrl,
    dmr_meth_treat = dmr_meth_treat,
    verbose = verbose
  )

  # ========================================================================
  # Step 4: Simulate Coverage and Methylation for Each Cell
  # ========================================================================

  if (verbose) {
    message("[Step 4/5] Simulating ", n_cells, " cells...")
    pb <- txtProgressBar(min = 0, max = n_cells, style = 3)
  }

  file_paths <- character(n_cells)
  n_covered_sites <- integer(n_cells)
  mean_coverage_vec <- numeric(n_cells)

  for (i in seq_len(n_cells)) {
    cell_id <- cell_metadata$cell_id[i]
    group <- cell_metadata$group[i]

    # Generate cell data
    cell_data <- simulate_single_cell(
      cpg_universe = cpg_universe,
      dmr_regions = dmr_regions,
      group = group,
      sparsity = sparsity,
      mean_coverage = mean_coverage,
      background_meth = background_meth
    )

    # Write to file
    file_path <- write_cell_file(
      cell_data = cell_data,
      cell_id = cell_id,
      output_dir = output_dir,
      output_format = output_format
    )

    file_paths[i] <- file_path
    n_covered_sites[i] <- nrow(cell_data)
    mean_coverage_vec[i] <- mean(cell_data$count_meth + cell_data$count_unmeth)

    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) {
    close(pb)
  }

  # Add file paths to metadata
  cell_metadata[, file_path := file_paths]
  cell_metadata[, n_covered_sites := n_covered_sites]
  cell_metadata[, mean_coverage := round(mean_coverage_vec, 2)]

  # ========================================================================
  # Step 5: Summary Statistics
  # ========================================================================

  if (verbose) {
    message("[Step 5/5] Generating summary...")

    total_possible <- n_cells * n_sites
    total_covered <- sum(n_covered_sites)
    actual_sparsity <- 1 - (total_covered / total_possible)

    message("\n=== Simulation Summary ===")
    message("Cells: ", n_cells)
    message(
      "  - Control: ",
      n_ctrl,
      " (",
      round(n_ctrl / n_cells * 100, 1),
      "%)"
    )
    message(
      "  - Treatment: ",
      n_treat,
      " (",
      round(n_treat / n_cells * 100, 1),
      "%)"
    )
    message("\nCpG Sites: ", format(n_sites, big.mark = ","))
    message("  - DMRs implanted: ", n_dmrs)
    message("  - CpGs in DMRs: ", n_dmrs * dmr_width)
    message(
      "  - Background sites: ",
      format(n_sites - n_dmrs * dmr_width, big.mark = ",")
    )
    message("\nCoverage Statistics:")
    message("  - Target sparsity: ", round(sparsity * 100, 1), "%")
    message("  - Actual sparsity: ", round(actual_sparsity * 100, 1), "%")
    message("  - Mean sites/cell: ", round(mean(n_covered_sites)))
    message("  - Mean depth/site: ", round(mean(mean_coverage_vec), 1), "x")
    message("\nOutput:")
    message("  - Directory: ", output_dir)
    message("  - Format: .", output_format)
    message("  - Files created: ", n_cells)

    total_size <- sum(file.size(file_paths)) / 1e6
    message("  - Total size: ", round(total_size, 2), " MB")
    message("==========================\n")
  }

  # ========================================================================
  # Return Results
  # ========================================================================

  result <- list(
    metadata = as.data.frame(cell_metadata),
    dmr_regions = dmr_regions,
    cpg_universe = cpg_universe,
    output_dir = output_dir,
    simulation_params = list(
      n_cells = n_cells,
      n_sites = n_sites,
      sparsity = sparsity,
      n_dmrs = n_dmrs,
      dmr_width = dmr_width,
      seed = seed
    )
  )

  class(result) <- c("scBS_simulation", "list")

  invisible(result)
}


# ============================================================================
# Helper Functions
# ============================================================================

#' Create CpG Universe (Genomic Coordinates)
#' @keywords internal
create_cpg_universe <- function(n_sites, verbose = TRUE) {
  # Simulate genomic positions across chromosomes
  # Use realistic chromosome distribution (chr1 is largest, etc.)

  chr_names <- paste0("chr", c(1:22, "X"))
  chr_weights <- c(
    seq(22, 1, length.out = 22), # chr1 has most sites
    10 # chrX moderate
  )
  chr_weights <- chr_weights / sum(chr_weights)

  # Assign sites to chromosomes
  chr_assignments <- sample(
    chr_names,
    size = n_sites,
    replace = TRUE,
    prob = chr_weights
  )

  # Generate positions within each chromosome
  # Assume ~250 Mb per chromosome, CpGs every ~100 bp on average
  cpg_universe <- data.table::data.table(
    chr = chr_assignments,
    pos = as.integer(runif(n_sites, min = 1000, max = 250000000))
  )

  # Sort by chromosome and position
  data.table::setorder(cpg_universe, chr, pos)

  # Add unique site IDs
  cpg_universe[, site_id := paste0(chr, "_", pos)]

  # Ensure unique positions (remove duplicates)
  cpg_universe <- unique(cpg_universe, by = "site_id")

  if (verbose) {
    message("  - Generated ", nrow(cpg_universe), " unique CpG sites")
    message("  - Chromosomes: ", length(unique(cpg_universe$chr)))
  }

  as.data.frame(cpg_universe)
}


#' Implant DMRs (Differentially Methylated Regions)
#' @keywords internal
implant_dmrs <- function(
  cpg_universe,
  n_dmrs,
  dmr_width,
  dmr_meth_ctrl,
  dmr_meth_treat,
  verbose = TRUE
) {
  cpg_dt <- data.table::as.data.table(cpg_universe)
  data.table::setorder(cpg_dt, chr, pos)

  # Randomly select DMR start positions
  # Ensure we have enough consecutive sites
  valid_starts <- which(
    seq_len(nrow(cpg_dt)) <= (nrow(cpg_dt) - dmr_width + 1)
  )

  if (length(valid_starts) < n_dmrs) {
    stop(
      "Not enough CpG sites to create ",
      n_dmrs,
      " DMRs of width ",
      dmr_width,
      call. = FALSE
    )
  }

  dmr_start_idx <- sort(sample(valid_starts, n_dmrs, replace = FALSE))

  # Define DMR regions
  dmr_list <- list()

  for (i in seq_len(n_dmrs)) {
    start_idx <- dmr_start_idx[i]
    end_idx <- start_idx + dmr_width - 1

    dmr_sites <- cpg_dt[start_idx:end_idx]

    dmr_list[[i]] <- data.frame(
      chr = dmr_sites$chr[1],
      start = dmr_sites$pos[1],
      end = dmr_sites$pos[nrow(dmr_sites)],
      dmr_id = sprintf("DMR_%03d", i),
      n_sites = dmr_width,
      true_meth_ctrl = dmr_meth_ctrl,
      true_meth_treat = dmr_meth_treat,
      delta_meth = abs(dmr_meth_ctrl - dmr_meth_treat),
      site_ids = I(list(dmr_sites$site_id))
    )
  }

  dmr_regions <- do.call(rbind, dmr_list)

  if (verbose) {
    message("  - DMR width: ", dmr_width, " CpG sites")
    message("  - Control methylation: ", dmr_meth_ctrl)
    message("  - Treatment methylation: ", dmr_meth_treat)
    message("  - Absolute difference: ", abs(dmr_meth_ctrl - dmr_meth_treat))
  }

  dmr_regions
}


#' Simulate a Single Cell's Data
#' @keywords internal
simulate_single_cell <- function(
  cpg_universe,
  dmr_regions,
  group,
  sparsity,
  mean_coverage,
  background_meth
) {
  n_sites <- nrow(cpg_universe)

  # Step 1: Determine which sites have coverage (inverse of sparsity)
  has_coverage <- runif(n_sites) > sparsity
  covered_idx <- which(has_coverage)

  if (length(covered_idx) == 0) {
    # Edge case: no coverage
    return(data.table::data.table(
      chr = character(0),
      pos = integer(0),
      count_meth = integer(0),
      count_unmeth = integer(0)
    ))
  }

  # Step 2: Simulate coverage depth (Poisson-distributed)
  coverage <- rpois(length(covered_idx), lambda = mean_coverage)
  coverage[coverage == 0] <- 1 # At least 1 read if covered

  # Step 3: Determine true methylation level for each site
  true_meth <- rep(background_meth, length(covered_idx))

  # Override with DMR methylation if site is in a DMR
  site_ids_covered <- cpg_universe$site_id[covered_idx]

  for (i in seq_len(nrow(dmr_regions))) {
    dmr_site_ids <- dmr_regions$site_ids[[i]]
    in_dmr <- site_ids_covered %in% dmr_site_ids

    if (any(in_dmr)) {
      if (group == "Control") {
        true_meth[in_dmr] <- dmr_regions$true_meth_ctrl[i]
      } else {
        true_meth[in_dmr] <- dmr_regions$true_meth_treat[i]
      }
    }
  }

  # Step 4: Simulate methylated counts (Binomial)
  count_meth <- rbinom(length(covered_idx), size = coverage, prob = true_meth)
  count_unmeth <- coverage - count_meth

  # Step 5: Create cell data.table
  cell_data <- data.table::data.table(
    chr = cpg_universe$chr[covered_idx],
    pos = cpg_universe$pos[covered_idx],
    count_meth = as.integer(count_meth),
    count_unmeth = as.integer(count_unmeth)
  )

  # Sort by chr and pos
  data.table::setorder(cell_data, chr, pos)

  cell_data
}


#' Write Cell Data to File
#' @keywords internal
write_cell_file <- function(
  cell_data,
  cell_id,
  output_dir,
  output_format
) {
  file_ext <- switch(
    output_format,
    "cov" = ".cov",
    "bed" = ".bed"
  )

  file_path <- file.path(output_dir, paste0(cell_id, file_ext))

  if (output_format == "cov") {
    # 4-column coverage format: chr, pos, count_meth, count_unmeth
    data.table::fwrite(
      cell_data,
      file = file_path,
      sep = "\t",
      col.names = FALSE,
      quote = FALSE,
      scipen = 999
    )
  } else if (output_format == "bed") {
    # 6-column BED-like: chr, start, end, name, count_meth, count_unmeth
    bed_data <- data.table::copy(cell_data)
    bed_data[, `:=`(
      start = pos,
      end = pos + 1L,
      name = "."
    )]

    # Reorder columns: chr, start, end, name, count_meth, count_unmeth
    bed_data <- bed_data[, .(chr, start, end, name, count_meth, count_unmeth)]

    data.table::fwrite(
      bed_data,
      file = file_path,
      sep = "\t",
      col.names = FALSE,
      quote = FALSE,
      scipen = 999
    )
  }

  file_path
}


#' Print method for scBS_simulation objects
#' @export
print.scBS_simulation <- function(x, ...) {
  cat("scBS Simulation Object\n")
  cat("======================\n")
  cat("Cells:", nrow(x$metadata), "\n")
  cat("  - Control:", sum(x$metadata$group == "Control"), "\n")
  cat("  - Treatment:", sum(x$metadata$group == "Treatment"), "\n")
  cat("CpG Sites:", nrow(x$cpg_universe), "\n")
  cat("DMRs Implanted:", nrow(x$dmr_regions), "\n")
  cat("Output Directory:", x$output_dir, "\n")
  cat("\nUse str() for detailed structure\n")
  invisible(x)
}


#' Summary method for scBS_simulation objects
#' @export
summary.scBS_simulation <- function(object, ...) {
  cat("=== scBS Simulation Summary ===\n\n")

  cat("Simulation Parameters:\n")
  cat("  - Cells:", object$simulation_params$n_cells, "\n")
  cat(
    "  - CpG Sites:",
    format(object$simulation_params$n_sites, big.mark = ","),
    "\n"
  )
  cat("  - Sparsity:", object$simulation_params$sparsity, "\n")
  cat("  - DMRs:", object$simulation_params$n_dmrs, "\n")
  cat("  - DMR Width:", object$simulation_params$dmr_width, "sites\n")
  cat("  - Seed:", object$simulation_params$seed, "\n\n")

  cat("Cell Metadata Summary:\n")
  print(summary(object$metadata[, c(
    "group",
    "n_covered_sites",
    "mean_coverage"
  )]))

  cat("\nDMR Summary:\n")
  cat("  - Total DMRs:", nrow(object$dmr_regions), "\n")
  cat(
    "  - Chromosomes with DMRs:",
    length(unique(object$dmr_regions$chr)),
    "\n"
  )
  cat(
    "  - Mean Control Methylation:",
    round(mean(object$dmr_regions$true_meth_ctrl), 3),
    "\n"
  )
  cat(
    "  - Mean Treatment Methylation:",
    round(mean(object$dmr_regions$true_meth_treat), 3),
    "\n"
  )
  cat(
    "  - Mean Absolute Difference:",
    round(mean(object$dmr_regions$delta_meth), 3),
    "\n"
  )

  invisible(object)
}
