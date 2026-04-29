#' @title Build methylation and coverage matrices (SECURE VERSION)
#'
#' @description
#' Process single-sample methylation bedGraph files into region-level
#' methylation and coverage matrices. For each sample, the function
#' (i) sorts and clips chromosomes, (ii) filters sites by coverage,
#' (iii) removes excluded regions (e.g. SNPs, blacklist),
#' (iv) maps per-base counts to the provided genomic annotation, and
#' finally merges all samples into a methylation matrix and a coverage
#' matrix, together with a processed annotation file (`3.ann.bed`).
#'
#' @param sam_info Character string. Path to a CSV file containing
#'   sample metadata. The file must contain a column named
#'   \code{sample_name}; its values must match the methylation file
#'   names in \code{input_dir}.
#' @param input_dir Character string. Directory that contains all
#'   per-sample methylation files (one file per sample).
#' @param annotation_file Character string. Path to a BED-like file
#'   defining genomic regions of interest (0-based), e.g. fixed windows,
#'   promoters, enhancers or CpG sites.
#' @param output_file Character string. Base file name for the merged
#'   methylation matrix to be written in \code{input_dir}. The coverage
#'   matrix will be written to \code{cov.<output_file>} in the same
#'   directory.
#' @param thresholds Numeric scalar. Minimum total coverage
#'   (\eqn{methylated + unmethylated} reads) required at a site
#'   before mapping it to regions.
#' @param bedtools Character string. Path to the \code{bedtools}
#'   executable (version 2.30.0 or higher).
#' @param bedSort Character string. Path to the UCSC \code{bedSort}
#'   executable (e.g. installed via
#'   \code{conda install ucsc-bedsort=469}).
#' @param lib Character string specifying the library type.
#'   One of \code{"TAPS"} (directional conversion; default)
#'   or \code{"WGBS"} (indirect conversion). This controls how
#'   methylated and unmethylated counts are interpreted.
#' @param parallel Integer. Number of threads used by
#'   \pkg{data.table} when reading/writing text files.
#' @param exclude_bed Character vector of paths to BED files
#'   containing regions to be excluded (e.g. SNPs, blacklist).
#'   All regions are concatenated and merged before filtering.
#' @param genome_type Character string. Genome assembly type,
#'   either \code{"human"} (default) or \code{"mouse"}.
#'   Controls which chromosomes are retained.
#' @param limitFiles Integer. Maximum number of file descriptors to
#'   allow when calling \code{paste} (used via \code{ulimit -n}).
#'
#' @importFrom stringi stri_join
#' @importFrom data.table fifelse := setnames fread fwrite .SD
#' @importFrom utils read.csv
#' @export
#'
#' @examples
#' ## See vignette for a full end-to-end example.
build_count_matrix <- function(
  sam_info = NULL,
  input_dir = NULL,
  annotation_file = NULL,
  output_file = NULL,
  thresholds = 1,
  bedtools = NULL,
  bedSort = NULL,
  lib = "TAPS",
  exclude_bed = NULL,
  parallel = 1,
  genome_type = "human",
  limitFiles = 200000
) {
  # ============================================================
  # CRITICAL SECURITY FIX: Input validation before ANY system call
  # ============================================================

  # Validate all file paths and parameters
  validate_inputs(
    sam_info = sam_info,
    input_dir = input_dir,
    annotation_file = annotation_file,
    output_file = output_file,
    bedtools = bedtools,
    bedSort = bedSort,
    exclude_bed = exclude_bed,
    lib = lib,
    genome_type = genome_type
  )

  ## read sample names
  sams <- utils::read.csv(sam_info, stringsAsFactors = FALSE)$sample_name

  # Validate sample names (no shell metacharacters)
  validate_sample_names(sams)

  ## define chromosomes
  if (genome_type == "human") {
    chrom <- c(paste0("chr", 1:22), "chrX")
  } else if (genome_type == "mouse") {
    chrom <- c(paste0("chr", 1:19), "chrX")
  } else {
    stop("'genome_type' must be 'human' or 'mouse'.", call. = FALSE)
  }

  handleAnnExc_secure(
    input_dir,
    parallel,
    annotation_file,
    exclude_bed,
    chrom,
    bedtools,
    bedSort
  )

  handleSam_secure(
    sams,
    input_dir,
    parallel,
    chrom,
    bedtools,
    lib,
    bedSort,
    thresholds
  )

  mergeMeth_secure(bedtools, sams, input_dir, output_file, limitFiles, parallel)
  mergeCov_secure(bedtools, sams, input_dir, output_file, limitFiles, parallel)
  cleanEnv_secure(input_dir)
}

# ============================================================
# SECURITY UTILITIES
# ============================================================

#' Validate all user inputs for security
#' @keywords internal
validate_inputs <- function(
  sam_info,
  input_dir,
  annotation_file,
  output_file,
  bedtools,
  bedSort,
  exclude_bed,
  lib,
  genome_type
) {
  # Check required parameters
  if (is.null(sam_info) || !file.exists(sam_info)) {
    stop("'sam_info' must be a valid file path", call. = FALSE)
  }

  if (is.null(input_dir) || !dir.exists(input_dir)) {
    stop("'input_dir' must be a valid directory", call. = FALSE)
  }

  if (is.null(annotation_file) || !file.exists(annotation_file)) {
    stop("'annotation_file' must be a valid file path", call. = FALSE)
  }

  if (is.null(output_file) || !is.character(output_file)) {
    stop("'output_file' must be a character string", call. = FALSE)
  }

  # Validate executable paths
  if (is.null(bedtools) || !file.exists(bedtools)) {
    stop("'bedtools' executable not found at: ", bedtools, call. = FALSE)
  }

  if (is.null(bedSort) || !file.exists(bedSort)) {
    stop("'bedSort' executable not found at: ", bedSort, call. = FALSE)
  }

  # Check for shell metacharacters in critical paths
  validate_path_safety(sam_info, "sam_info")
  validate_path_safety(input_dir, "input_dir")
  validate_path_safety(annotation_file, "annotation_file")
  validate_path_safety(output_file, "output_file")
  validate_path_safety(bedtools, "bedtools")
  validate_path_safety(bedSort, "bedSort")

  # Validate exclude_bed if provided
  if (!is.null(exclude_bed)) {
    for (i in seq_along(exclude_bed)) {
      if (!file.exists(exclude_bed[i])) {
        stop("exclude_bed file not found: ", exclude_bed[i], call. = FALSE)
      }
      validate_path_safety(exclude_bed[i], paste0("exclude_bed[", i, "]"))
    }
  }

  # Validate lib parameter
  if (!lib %in% c("TAPS", "WGBS")) {
    stop("'lib' must be either 'TAPS' or 'WGBS'", call. = FALSE)
  }

  # Validate genome_type
  if (!genome_type %in% c("human", "mouse")) {
    stop("'genome_type' must be either 'human' or 'mouse'", call. = FALSE)
  }

  invisible(TRUE)
}

#' Check for dangerous shell metacharacters in paths
#' @keywords internal
validate_path_safety <- function(path, param_name) {
  # Dangerous characters that could enable command injection
  dangerous_chars <- c(";", "|", "&", "$", "`", "(", ")", "<", ">", "\n", "\r")

  for (char in dangerous_chars) {
    if (grepl(char, path, fixed = TRUE)) {
      stop(
        "Security error: '",
        param_name,
        "' contains dangerous character '",
        char,
        "'. This could enable command injection.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

#' Validate sample names for safety
#' @keywords internal
validate_sample_names <- function(sams) {
  for (i in seq_along(sams)) {
    validate_path_safety(sams[i], paste0("sample_name[", i, "]"))

    # Additional check: no directory traversal
    if (grepl("\\.\\.", sams[i], fixed = TRUE)) {
      stop(
        "Security error: sample_name '",
        sams[i],
        "' contains '..' (directory traversal attempt)",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

#' Safe wrapper for system2 calls
#' @keywords internal
safe_system2 <- function(command, args, ...) {
  # Quote all arguments to prevent injection
  args_quoted <- shQuote(args)

  # Execute with error checking
  result <- system2(
    command = command,
    args = args_quoted,
    stdout = TRUE,
    stderr = TRUE,
    ...
  )

  # Check exit status
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
    stop(
      "Command failed: ",
      command,
      "\n",
      "Args: ",
      paste(args, collapse = " "),
      "\n",
      "Error: ",
      paste(result, collapse = "\n"),
      call. = FALSE
    )
  }

  invisible(result)
}

# ============================================================
# SECURE IMPLEMENTATIONS
# ============================================================

handleAnnExc_secure <- function(
  input_dir,
  parallel,
  annotation_file,
  exclude_bed,
  chrom,
  bedtools,
  bedSort
) {
  # All paths are now properly quoted using shQuote()

  # Sort annotation
  safe_system2(
    command = bedSort,
    args = c(annotation_file, file.path(input_dir, "1.ann.sort"))
  )

  # Concatenate exclude files using R (avoid shell cat)
  if (!is.null(exclude_bed) && length(exclude_bed) > 0) {
    # Use R's file operations instead of shell 'cat'
    cat_output <- file.path(input_dir, "1.exclude.cat")
    cat_conn <- file(cat_output, "w")
    on.exit(close(cat_conn), add = TRUE)

    for (bed_file in exclude_bed) {
      bed_content <- readLines(bed_file)
      writeLines(bed_content, cat_conn)
    }
    close(cat_conn)
    on.exit() # Remove the close handler

    # Cut first 3 columns using R instead of shell 'cut'
    exclude_data <- data.table::fread(
      cat_output,
      nThread = parallel,
      tmpdir = input_dir,
      verbose = FALSE,
      showProgress = FALSE
    )
    data.table::fwrite(
      exclude_data[, 1:3],
      file.path(input_dir, "2.exclude.3col"),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      sep = "\t",
      scipen = 999,
      nThread = parallel,
      verbose = FALSE,
      showProgress = FALSE
    )

    # Sort exclude
    safe_system2(
      command = bedSort,
      args = c(
        file.path(input_dir, "2.exclude.3col"),
        file.path(input_dir, "3.exclude.sort")
      )
    )

    # Merge exclude regions
    merge_output <- file.path(input_dir, "4.exclude.btmerge")
    safe_system2(
      command = bedtools,
      args = c("merge", "-i", file.path(input_dir, "3.exclude.sort")),
      stdout = merge_output
    )
  }

  # Process annotation (R operations - safe)
  annSort <- data.table::fread(
    file.path(input_dir, "1.ann.sort"),
    nThread = parallel,
    tmpdir = input_dir,
    verbose = FALSE,
    showProgress = FALSE
  )
  annSortC <- annSort[annSort$V1 %in% chrom]

  data.table::fwrite(
    annSortC,
    file.path(input_dir, "3.ann.bed"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t",
    scipen = 999,
    nThread = parallel,
    verbose = FALSE,
    showProgress = FALSE
  )

  data.table::fwrite(
    annSortC[, 1:3],
    file.path(input_dir, "2.ann.3col"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t",
    scipen = 999,
    nThread = parallel,
    verbose = FALSE,
    showProgress = FALSE
  )

  # Process exclude (if exists)
  if (file.exists(file.path(input_dir, "4.exclude.btmerge"))) {
    exc <- data.table::fread(
      file.path(input_dir, "4.exclude.btmerge"),
      nThread = parallel,
      tmpdir = input_dir,
      verbose = FALSE,
      showProgress = FALSE
    )
    excC <- exc[exc$V1 %in% chrom]

    data.table::fwrite(
      excC,
      file.path(input_dir, "5.exclude.bed"),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      sep = "\t",
      scipen = 999,
      nThread = parallel,
      verbose = FALSE,
      showProgress = FALSE
    )
  }
}

handleSam_secure <- function(
  sams,
  input_dir,
  parallel,
  chrom,
  bedtools,
  lib,
  bedSort,
  thresholds
) {
  for (i in seq_along(sams)) {
    sam <- sams[i]
    samF <- file.path(input_dir, sam)

    # Sort sample
    safe_system2(
      command = bedSort,
      args = c(samF, file.path(input_dir, "1.sam.sort"))
    )

    # Filter by coverage threshold (R - safe)
    samFs <- data.table::fread(
      file.path(input_dir, "1.sam.sort"),
      nThread = parallel,
      tmpdir = input_dir,
      verbose = FALSE,
      showProgress = FALSE
    )
    samt <- samFs[(samFs$V5 + samFs$V6) >= thresholds]
    samFsC <- samt[samt$V1 %in% chrom]

    data.table::fwrite(
      samFsC,
      file.path(input_dir, "2.sam.clip"),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      sep = "\t",
      scipen = 999,
      nThread = parallel,
      verbose = FALSE,
      showProgress = FALSE
    )

    # Intersect with exclude regions
    clean_output <- file.path(input_dir, "3.sam.clean")
    safe_system2(
      command = bedtools,
      args = c(
        "intersect",
        "-a",
        file.path(input_dir, "2.sam.clip"),
        "-b",
        file.path(input_dir, "5.exclude.bed"),
        "-v"
      ),
      stdout = clean_output
    )

    # Map to annotation
    map_output <- file.path(input_dir, "4.sam.map")
    safe_system2(
      command = bedtools,
      args = c(
        "map",
        "-a",
        file.path(input_dir, "2.ann.3col"),
        "-b",
        file.path(input_dir, "3.sam.clean"),
        "-c",
        "5,6",
        "-o",
        "sum,sum",
        "-null",
        "nan"
      ),
      stdout = map_output
    )

    # Process results (R - safe)
    if (lib == "WGBS") {
      me <- data.table::fread(
        file.path(input_dir, "4.sam.map"),
        nThread = parallel,
        tmpdir = input_dir,
        verbose = FALSE,
        showProgress = FALSE
      )
    } else if (lib == "TAPS") {
      me <- data.table::fread(
        file.path(input_dir, "4.sam.map"),
        select = c(1, 2, 3, 5, 4),
        nThread = parallel,
        tmpdir = input_dir,
        verbose = FALSE,
        showProgress = FALSE
      )
    }

    data.table::setnames(me, c("V1", "V2", "V3", "V4", "V5")[seq_len(ncol(me))])

    me[,
      beta := data.table::fifelse(
        is.na(.SD$V4) | is.na(.SD$V5) | abs(.SD$V4 + .SD$V5) < 1e-10,
        NaN,
        .SD$V4 / (.SD$V4 + .SD$V5)
      )
    ]

    me[,
      c("cover") := data.table::fifelse(
        is.na(.SD$V4) | is.na(.SD$V5) | abs(.SD$V4 + .SD$V5) < 1e-10,
        NaN,
        .SD$V4 + .SD$V5
      )
    ]

    me$beta <- me$beta
    meC <- me[me$V1 %in% chrom]

    data.table::fwrite(
      meC[, c(1:3, 6)],
      paste0(samF, ".meth"),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      sep = "\t",
      scipen = 999,
      na = "nan",
      nThread = parallel,
      verbose = FALSE,
      showProgress = FALSE
    )

    data.table::fwrite(
      meC[, c(1:3, 7)],
      paste0(samF, ".cov"),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      sep = "\t",
      scipen = 999,
      na = "nan",
      nThread = parallel,
      verbose = FALSE,
      showProgress = FALSE
    )
  }
}

mergeMeth_secure <- function(
  bedtools,
  sams,
  input_dir,
  output_file,
  limitFiles,
  parallel
) {
  # Use R's data operations instead of shell paste
  colN <- data.table::fread(
    file.path(input_dir, "2.ann.3col"),
    nThread = parallel,
    tmpdir = input_dir,
    verbose = FALSE,
    showProgress = FALSE
  )
  colN$marker_name <- paste0(colN$V1, "_", colN$V2, "_", colN$V3)

  # Read all meth files and combine using R (avoid bash/paste)
  meth_list <- vector("list", length(sams) + 1)
  meth_list[[1]] <- colN[, 4, drop = FALSE]

  for (i in seq_along(sams)) {
    fn <- file.path(input_dir, paste0(sams[i], ".meth"))
    meth_data <- data.table::fread(
      fn,
      select = 4,
      nThread = parallel,
      tmpdir = input_dir,
      verbose = FALSE,
      showProgress = FALSE
    )
    meth_list[[i + 1]] <- meth_data
  }

  # Combine all columns
  combined <- do.call(cbind, meth_list)
  colnames(combined) <- c("marker_name", sams)

  # Write output
  data.table::fwrite(
    combined,
    file.path(input_dir, output_file),
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t",
    scipen = 999,
    nThread = parallel,
    verbose = FALSE,
    showProgress = FALSE
  )
}

mergeCov_secure <- function(
  bedtools,
  sams,
  input_dir,
  output_file,
  limitFiles,
  parallel
) {
  # Use R's data operations instead of shell paste
  colN <- data.table::fread(
    file.path(input_dir, "2.ann.3col"),
    nThread = parallel,
    tmpdir = input_dir,
    verbose = FALSE,
    showProgress = FALSE
  )
  colN$marker_name <- paste0(colN$V1, "_", colN$V2, "_", colN$V3)

  # Read all cov files and combine using R (avoid bash/paste)
  cov_list <- vector("list", length(sams) + 1)
  cov_list[[1]] <- colN[, 4, drop = FALSE]

  for (i in seq_along(sams)) {
    fn <- file.path(input_dir, paste0(sams[i], ".cov"))
    cov_data <- data.table::fread(
      fn,
      select = 4,
      nThread = parallel,
      tmpdir = input_dir,
      verbose = FALSE,
      showProgress = FALSE
    )
    cov_list[[i + 1]] <- cov_data
  }

  # Combine all columns
  combined <- do.call(cbind, cov_list)
  colnames(combined) <- c("marker_name", sams)

  # Write output
  data.table::fwrite(
    combined,
    file.path(input_dir, paste0("cov.", output_file)),
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t",
    scipen = 999,
    nThread = parallel,
    verbose = FALSE,
    showProgress = FALSE
  )
}

cleanEnv_secure <- function(input_dir) {
  # Use R's file.remove instead of shell rm for safety
  temp_files <- c(
    file.path(input_dir, "1.ann.sort"),
    file.path(input_dir, "2.ann.3col"),
    file.path(input_dir, "1.exclude.cat"),
    file.path(input_dir, "2.exclude.3col"),
    file.path(input_dir, "3.exclude.sort"),
    file.path(input_dir, "4.exclude.btmerge"),
    file.path(input_dir, "5.exclude.bed"),
    file.path(input_dir, "1.sam.sort"),
    file.path(input_dir, "2.sam.clip"),
    file.path(input_dir, "3.sam.clean"),
    file.path(input_dir, "4.sam.map")
  )

  # Remove files that exist
  for (f in temp_files) {
    if (file.exists(f)) {
      unlink(f, force = TRUE)
      if (file.exists(f)) {
        warning("Failed to delete temporary file: ", f)
      }
    }
  }

  # Remove temp.* files using pattern matching (safe in R)
  temp_pattern_files <- list.files(
    input_dir,
    pattern = "^[12]\\.temp",
    full.names = TRUE
  )

  for (f in temp_pattern_files) {
    unlink(f, force = TRUE)
  }

  invisible(NULL)
}
