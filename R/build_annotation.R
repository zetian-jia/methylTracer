#' @title Build CpG BED File for Base Sites
#' @description 
#' Generates a BED file (0-based, half-open intervals) containing 
#' coordinates of CpG sites based on the hg38 genome. Uses the 
#' BSgenome.Hsapiens.UCSC.hg38 package to extract genomic data and create 
#' a BED file suitable for downstream analysis. Currently supports only 
#' 'CpG' annotations.
#'
#' @param ref_genome A string specifying the annotation, e.g. 
#' 'BSgenome.Hsapiens.UCSC.hg38'.
#' @param output A string for the output file path, e.g. 
#' 'PATH/TO/cpg_sites.bed'.
#' @param genome_type A string specifying the genome type, e.g. 'human' 
#' or 'mouse'.
#'
#' @return 
#' A BED file with CpG site coordinates saved to the specified output 
#' directory.
#'
#' @details 
#' The function builds an annotation for CpG sites using the 
#' BSgenome.Hsapiens.UCSC.hg38 package. Coordinates are formatted into a 
#' BED file.
#'
#' @export
#'
#' @importFrom data.table fwrite as.data.table
#' @importFrom IRanges shift IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom BSgenome getBSgenome
#' @importFrom Biostrings matchPattern
#'
#' @examples
#'
#' ## Example usage (about 3 minutes to run):
#' output_dir <- tempdir()
#' ## build_cpg(ref_genome = "BSgenome.Hsapiens.UCSC.hg38", 
#' ##           output = file.path(output_dir, "cpg_sites.bed"), )
#' 
#' ## $ head cpg_sites.bed
#' ## Output:
#' ## chr1    10468   10469   2       *
#' ## chr1    10470   10471   2       *
#' ## chr1    10483   10484   2       *
#' ## chr1    10488   10489   2       *
build_cpg <- function(
    ref_genome = NULL, output = NULL, genome_type = "human") {
    if (is.null(ref_genome)) {
        message("ref_genome cannot be NULL")
    }
    if (is.null(output)) {
        message("Output filename cannot be NULL")
    }
    ref_genome <- BSgenome::getBSgenome(genome = ref_genome)
    ## https://github.com/CompEpigen/methrix 
    ## https://support.bioconductor.org/p/95239/ 
    if (genome_type == "human") {
        chrs <- paste0("chr", c(1:22, "X", "Y"))
    } else if (genome_type == "mouse") {
        chrs <- paste0("chr", c(1:19, "X", "Y"))
    } else {
    stop("Please use 'human' or 'mouse'.")
    }
    ## Find CpG sites in each chromosome
    message("Finding CpG sites in each chromosome...")
    cgs <- lapply(
        chrs, function(x) {
            start(Biostrings::matchPattern("CG", ref_genome[[x]]))
        }
    )
    ## then shift -1 to get 0-based C position
    cpgr <- do.call(
        c, lapply(
            seq_along(chrs),
            function(x) {
                GenomicRanges::GRanges(
                  seqnames = chrs[x], 
                  ranges = IRanges(start = cgs[[x]], width = 2))
            }
        )
    )
    cpgr <- IRanges::shift(cpgr, shift = -1)
    out <- data.table::as.data.table(cpgr)
    data.table::fwrite(
      out, output, row.names = FALSE, 
      col.names = FALSE, sep = "\t", quote = FALSE)
    message("Successfully wrote CpG sites to ", output)
}
