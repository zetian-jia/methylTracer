#' @title Build methylation and coverage matrices
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
#' @details
#' Input methylation files are expected to be bedGraph-like text files
#' with at least six columns:
#' \code{chr}, \code{start}, \code{end}, and two columns of
#' methylated and unmethylated counts. After filtering and mapping
#' to the annotation, the function computes, for each region and
#' each sample:
#'
#' \itemize{
#'   \item \strong{methylation proportion} \eqn{beta} in the range
#'         \eqn{[0, 1]}, defined as
#'         \eqn{beta = methylated / (methylated + unmethylated)}.
#'   \item \strong{coverage} as the total read count
#'         \eqn{methylated + unmethylated}.
#' }
#'
#' For each sample, intermediate files \code{<sample>.meth}
#' (region-level beta) and \code{<sample>.cov} (region-level coverage)
#' are written in \code{input_dir}. These are then merged into:
#'
#' \itemize{
#'   \item a methylation matrix: first column \code{marker_name}
#'         (\code{chr_start_end}), subsequent columns one per sample,
#'         values are methylation proportions in \eqn{[0, 1]};
#'   \item a coverage matrix with the same structure
#'         (\code{cov.<output_file>}).
#' }
#'
#' The filtered and chromosome-restricted annotation is written to
#' \code{3.ann.bed} in \code{input_dir}.
#'
#' @return
#' This function is called for its side effects. It writes three main
#' files into \code{input_dir}:
#' \itemize{
#'   \item \code{<output_file>}: merged methylation matrix
#'         (proportions 0–1).
#'   \item \code{cov.<output_file>}: merged coverage matrix.
#'   \item \code{3.ann.bed}: processed genomic annotation BED file.
#' }
#'
#' @importFrom stringi stri_join
#' @importFrom data.table fifelse := setnames fread fwrite .SD
#' @importFrom utils read.csv
#' @export
#'
#' @examples
#' ## See vignette for a full end-to-end example.
build_count_matrix <- function(
    sam_info = NULL, input_dir = NULL, 
    annotation_file = NULL, output_file = NULL,
    thresholds = 1, bedtools = NULL, bedSort = NULL, 
    lib = "TAPS", exclude_bed = NULL, parallel = 1,
    genome_type = "human", limitFiles=200000
) {
  ## read sample names
  sams <- utils::read.csv(sam_info, stringsAsFactors = FALSE)$sample_name
  ## define chromosomes
  if (genome_type == "human") {
    chrom <- c(paste0("chr", 1:22),"chrX")
  } else if (genome_type == "mouse") {
    chrom <- c(paste0("chr", 1:19),"chrX")
  } else {
    message("'genome_type': 'human' or 'mouse'.")
  }
  handleAnnExc(input_dir, parallel, annotation_file, exclude_bed,
               chrom, bedtools, bedSort)
  ## handle samples
  handleSam(sams, input_dir, parallel, chrom, bedtools, lib, bedSort,
            thresholds)
  ## merge meth
  mergeMeth(bedtools, sams, input_dir, output_file, limitFiles, parallel)
  ## merge cov
  mergeCov(bedtools, sams, input_dir, output_file, limitFiles, parallel)
  cleanEnv(input_dir)
}

handleAnnExc <- function(
    input_dir,parallel,annotation_file,exclude_bed,chrom,bedtools,bedSort
)
{
  system2(command = bedSort,
          args = c(annotation_file, paste0(input_dir,"/1.ann.sort")))
  system2(command = "cat",
          args = c(stringi::stri_join(exclude_bed, collapse = " "), ">", 
                   paste0(input_dir,"/1.exclude.cat")))
  system2(command = "cut",
          args = c("-f 1-3", paste0(input_dir,"/1.exclude.cat"), ">", 
                   paste0(input_dir,"/2.exclude.3col")))
  system2(command = bedSort,
          args = c(paste0(input_dir,"/2.exclude.3col"), 
                   paste0(input_dir,"/3.exclude.sort")))
  system2(command = bedtools, 
          args = c("merge -i",
                   paste0(input_dir,"/3.exclude.sort"), ">",
                   paste0(input_dir,"/4.exclude.btmerge")))
  annSort <- data.table::fread(paste0(input_dir,"/1.ann.sort"), 
                               nThread=parallel, tmpdir=input_dir,
                               verbose = FALSE, showProgress = FALSE)
  annSortC <- annSort[annSort$V1 %in% chrom]
  data.table::fwrite(annSortC, paste0(input_dir,"/3.ann.bed"),
                     row.names = FALSE, col.names = FALSE, quote = FALSE,
                     sep = "\t", scipen = 999,
                     nThread=parallel, verbose = FALSE, showProgress = FALSE)
  data.table::fwrite(annSortC[,1:3], paste0(input_dir,"/2.ann.3col"),
                     row.names = FALSE, col.names = FALSE, quote = FALSE,
                     sep = "\t", scipen = 999,
                     nThread=parallel, verbose = FALSE, showProgress = FALSE)
  exc <- data.table::fread(paste0(input_dir,"/4.exclude.btmerge"),
                           nThread=parallel, tmpdir=input_dir,
                           verbose = FALSE, showProgress = FALSE)
  excC <- exc[exc$V1 %in% chrom]
  data.table::fwrite(excC, paste0(input_dir,"/5.exclude.bed"),
                     row.names = FALSE, col.names = FALSE, quote = FALSE,
                     sep = "\t", scipen = 999,
                     nThread=parallel, verbose = FALSE, showProgress = FALSE)
}

handleSam <- function(
    sams,input_dir,parallel,chrom,bedtools, lib, bedSort, thresholds
){
  for (i in seq_along(sams)) {
    sam <- sams[i]
    samF <- file.path(input_dir, sam)
    system2(command = bedSort,
            args = c(samF, paste0(input_dir,"/1.sam.sort")))
    samFs <- data.table::fread(paste0(input_dir,"/1.sam.sort"),
                               nThread=parallel, tmpdir=input_dir,
                               verbose = FALSE, showProgress = FALSE)
    samt <- samFs[(samFs$V5 + samFs$V6) >= thresholds]
    samFsC <- samt[samt$V1 %in% chrom]
    data.table::fwrite(samFsC, paste0(input_dir,"/2.sam.clip"),
                       row.names = FALSE, col.names = FALSE, quote = FALSE,
                       sep = "\t", scipen = 999,
                       nThread=parallel, verbose = FALSE, showProgress = FALSE)
    system2(command = bedtools,  
            args = c("intersect -a", paste0(input_dir,"/2.sam.clip"), "-b",
                     paste0(input_dir,"/5.exclude.bed"), "-v >",
                     paste0(input_dir,"/3.sam.clean")))
    system2(command = bedtools, 
            args = c("map -a", paste0(input_dir,"/2.ann.3col"), 
            "-b", paste0(input_dir,"/3.sam.clean"),
            "-c 5,6 -o sum,sum -null nan >", paste0(input_dir,"/4.sam.map")))
    if (lib == "WGBS") {
      me <- data.table::fread(paste0(input_dir,"/4.sam.map"),
                              nThread=parallel, tmpdir=input_dir,
                              verbose = FALSE, showProgress = FALSE)
    } else if (lib == "TAPS") {
      me <- data.table::fread(paste0(input_dir,"/4.sam.map"),
                              select = c(1, 2, 3, 5, 4),
                              nThread=parallel, tmpdir=input_dir,
                              verbose = FALSE, showProgress = FALSE)
    }
    data.table::setnames(me, 
                         c("V1", "V2", "V3", "V4", "V5")[seq_len(ncol(me))])
    me[, beta := data.table::fifelse(is.na(.SD$V4) | 
                                       is.na(.SD$V5) | .SD$V4 + .SD$V5 == 0, 
                                     NaN, 
                                     .SD$V4 / (.SD$V4 + .SD$V5))]
    # me[, cover := data.table::fifelse(is.na(.SD$V4) | 
    #                                     is.na(.SD$V5) | 
    #                                     .SD$V4 + .SD$V5 == 0, NaN, 
    #                                   .SD$V4 + .SD$V5)]
    me[, c("cover") := data.table::fifelse(is.na(.SD$V4) | is.na(.SD$V5) |
                                               (.SD$V4 + .SD$V5 == 0), NaN, 
                                             .SD$V4 + .SD$V5)]
    # me$beta <- round(me$beta * 1000, 0); meC <- me[me$V1 %in% chrom]
    me$beta <- me$beta; meC <- me[me$V1 %in% chrom]
    data.table::fwrite(meC[,c(1:3, 6)], paste0(samF, ".meth"),
                       row.names = FALSE, col.names = FALSE, quote = FALSE,
                       sep = "\t", scipen = 999, na = "nan",
                       nThread=parallel, verbose = FALSE, showProgress = FALSE)
    data.table::fwrite(meC[,c(1:3, 7)], paste0(samF, ".cov"),
                       row.names = FALSE, col.names = FALSE, quote = FALSE,
                       sep = "\t", scipen = 999, na = "nan",
                       nThread=parallel, verbose = FALSE, showProgress = FALSE)
  }
}

mergeMeth <- function(
    bedtools, sams, input_dir, output_file,limitFiles,parallel
){
  colN <- data.table::fread(paste0(input_dir, "/2.ann.3col"),
                            nThread=parallel, tmpdir=input_dir,
                            verbose = FALSE, showProgress = FALSE)
  colN$marker_name <- paste0(colN$V1, "_", colN$V2, "_", colN$V3)
  data.table::fwrite(colN[,4], paste0(input_dir, "/1.temp.rowname"),
                     row.names = FALSE, col.names = FALSE, quote = FALSE,
                     sep = "\t", scipen = 999,
                     nThread=parallel, verbose = FALSE, showProgress = FALSE)
  meth_files <- character(length(sams))
  for (i in seq_along(sams)) {
    fn <- paste0(input_dir, "/", sams[i], ".meth")
    tmp_fn <- paste0(input_dir, "/2.temp.", sams[i])
    system2(command = "cut", args = c("-f 4", fn, ">", tmp_fn))
    meth_files[i] <- tmp_fn
  }
  pascmd <- paste0("ulimit -n ", as.integer(limitFiles), " && paste ", 
                   paste0(input_dir, "/1.temp.rowname "),
                   stringi::stri_join(meth_files, collapse = " "), " > ",
                   paste0(input_dir, "/2.temp.paste"))
  writeLines(pascmd, con = paste0(input_dir, "/2.temp.paste.sh"))
  system2(command = "bash", args = paste0(input_dir, "/2.temp.paste.sh"))

  header_line <- stringi::stri_join(c("marker_name", sams), collapse = "\t")
  header_file <- paste0(input_dir, "/3.temp.header")
  writeLines(header_line, con = header_file)
  system2(command = "cat", 
          args = c(header_file, 
                   paste0(input_dir, "/2.temp.paste"), 
                    ">", file.path(input_dir, output_file)))
}

mergeCov <- function(
    bedtools, sams, input_dir, output_file, limitFiles, parallel
){
  colN <- data.table::fread(paste0(input_dir, "/2.ann.3col"),
                            nThread=parallel, tmpdir=input_dir,
                            verbose = FALSE, showProgress = FALSE)
  colN$marker_name <- paste0(colN$V1, "_", colN$V2, "_", colN$V3)
  data.table::fwrite(colN[,4], paste0(input_dir, "/1.temp.rowname"),
                     row.names = FALSE, col.names = FALSE, quote = FALSE,
                     sep = "\t", scipen = 999,
                     nThread=parallel, verbose = FALSE, showProgress = FALSE)
  meth_files <- character(length(sams))
  for (i in seq_along(sams)) {
    fn <- paste0(input_dir, "/", sams[i], ".cov")
    tmp_fn <- paste0(input_dir, "/2.temp.", sams[i])
    system2(command = "cut", args = c("-f 4", fn, ">", tmp_fn))
    meth_files[i] <- tmp_fn 
  }
  pascmd <- paste0("ulimit -n ", as.integer(limitFiles), " && paste ", 
                   paste0(input_dir, "/1.temp.rowname "),
                   stringi::stri_join(meth_files, collapse = " "), " > ",
                   paste0(input_dir, "/2.temp.paste"))
  writeLines(pascmd, con = paste0(input_dir, "/2.temp.paste.sh"))
  system2(command = "bash", args = paste0(input_dir, "/2.temp.paste.sh"))
  
  header_line <- stringi::stri_join(c("marker_name", sams), collapse = "\t")
  header_file <- paste0(input_dir, "/3.temp.header")
  writeLines(header_line, con = header_file)
  system2(command = "cat", 
          args = c(header_file, 
                   paste0(input_dir, "/2.temp.paste"), 
                   ">", paste0(input_dir, "/cov.", output_file)))
}

cleanEnv <- function(input_dir){
  system2(command = "rm",
          args = c(paste0(input_dir,"/1.ann.sort"),
                   paste0(input_dir,"/2.ann.3col"),
                   paste0(input_dir,"/1.exclude.cat"),
                   paste0(input_dir,"/2.exclude.3col"),
                   paste0(input_dir,"/3.exclude.sort"),
                   paste0(input_dir,"/4.exclude.btmerge"),
                   paste0(input_dir,"/5.exclude.bed")
          ))
  system2(command = "rm",
          args = c(paste0(input_dir,"/1.sam.sort"),
                   paste0(input_dir,"/2.sam.clip"),
                   paste0(input_dir,"/3.sam.clean"),
                   paste0(input_dir,"/4.sam.map")))
  system2(command = "rm", 
          args = c(paste0(input_dir, "/1.temp.rowname"),
                   paste0(input_dir, "/2.temp.paste"),
                   paste0(input_dir, "/2.temp.paste.sh"),
                   paste0(input_dir, "/3.temp.header")))
  system2(command = "rm", 
          args = c(paste0(input_dir, "/2.temp*")))

}
