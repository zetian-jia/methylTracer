#' @title Build a Methylation Count Matrix
#' 
#' @description 
#' Processes a set of .bedGraph files, filters them based on coverage 
#' thresholds, and excludes regions defined in blacklist .bed files. 
#' Uses human or mouse genome as reference, and calculates a methylation 
#' count matrix with bedtools. The result is saved as a data.frame.
#'
#' @param sam_info Path to a sample file with sample or cell annotations. 
#' The first column must be named 'sample_name'.
#' @param input_dir Path to the directory containing .bedGraph files. 
#' Sample names should match those in sam_info. No trailing slash.
#' @param annotation_file A character vector of paths to .bed files defining 
#' genomic windows. The length must match output_file entries.
#' @param output_file A character vector of output file names for methylation 
#' count matrices. The length must match annotation_file.
#' @param thresholds Integer. Coverage threshold for base resolution filtering.
#' @param bedtools Path to bedtools executable.
#' @param lib Character string. Either 'TAPS' (default) or 'WGBS'. 'TAPS' 
#' refers to directional conversion, 'WGBS' refers to indirect conversion.
#' @param parallel Integer. Cores for parallel processing. Default is 
#' 1.
#' @param exclude_bed A character vector of .bed files specifying regions to 
#' exclude.
#' @param genome_type Character string. Either 'human' (default) or 'mouse'.
#'
#' @details
#' This function processes .bedGraph files to calculate a methylation count 
#' matrix, filtering based on coverage thresholds user-specified 
#' regions. The matrix is computed using bedtools and can be based on either 
#' the human or mouse genome. The output is written as a data.frame and saved 
#' to the specified output_file.
#'
#' @return The function writes the methylation count matrix to the specified 
#' output_file.
#' \itemize{
#'    \item \code{columnames}: input filenames
#'    \item \code{rownames}: Unique identifier for each genomic region 
#'    (chr_start_end)
#'    \item \code{values}: Methylation signal (range: 0 to 1000). By default, 
#'    methylTracer saves methylation data scaled by a factor of 1000.
#' }
#'
#' @import utils
#' @importFrom stringi stri_join
#' @export
#' @examples 
#' output_dir <- tempdir()
#'
#' sam_info_df <- data.frame(
#'   sample_name = c('sample_1.bed', 'sample_2.bed', 
#'                   'sample_3.bed', 'sample_4.bed'),
#'   group = c(rep('case', 2), rep('ctrl', 2))
#' )
#' 
#' sam_info <- file.path(output_dir, 'sample_info.csv')
#'
#' input_file_df <- data.frame(
#'   marker_name = c('chr1_1000_2000', 'chr1_2000_3000',
#'                   'chr1_3000_4000', 'chr1_4000_5000'),
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
#'   marker_name = c('chr1_1000_2000', 'chr1_2000_3000',
#'                   'chr1_3000_4000', 'chr1_4000_5000')
#' )
#' 
#' annotation_file <- file.path(output_dir, 'annotation.bed')
#'
#' output_file <- 'methylTracer_obj_test.h5'
#'
#' unlink(file.path(output_dir, output_file), recursive = TRUE)
#'
#' write.csv(sam_info_df, sam_info, row.names = FALSE)
#' write.table(input_file_df, input_file, sep = '\t', 
#'             row.names = FALSE)
#' write.table(annotation_file_df, annotation_file, 
#'              sep = '\t', row.names = FALSE)
#'
#' ## build_count_matrix(sam_info=sam_info, 
#' ##                   input_dir=input_dir,
#' ##                   annotation_file=annotation_file,
#' ##                   output_file=output_file,
#' ##                   thresholds=3,
#' ##                   bedtools="PATH/bedtool")

build_count_matrix <- function(
    sam_info = NULL, input_dir = NULL, 
    annotation_file = NULL, output_file = NULL,
    thresholds = 3, bedtools = NULL, 
    lib = NULL, parallel = 5, exclude_bed = NULL,
    genome_type = "human"
) {
    start_time <- Sys.time()
    sam_info_lines <- readLines(sam_info, n = 1)
    if (!grepl("sample_name", sam_info_lines)) {
        message("'sam_info' must contain a column named 'sample_name'.")
    }
    path_dir <- input_dir
    samples <- utils::read.csv(sam_info, stringsAsFactors = FALSE)$sample_name

    message(sprintf("Total %d cells/samples !", length(samples)))
    ## Define valid chromosomes based on genome type
    valid_chroms <- chrom_pat_fun(genome_type = genome_type)

    ann_sort_fun(
        annotation_file = annotation_file, 
        path_dir = path_dir, valid_chroms = valid_chroms,
        parallel = parallel
    )

    exclude_bed_combined <- exc_sort_fun(
        path_dir = path_dir, 
        genome_type = genome_type, exclude_bed = exclude_bed,
        parallel = parallel
    )

    ## Process each sample file
    message("Starting matrix and annotation building...")
    for (i in seq_along(samples)) {
        sam <- samples[i]
        sample_path <- file.path(path_dir, sam)
        ## Check if sample file exists
        if (!file.exists(sample_path)) {
            message(sprintf("Sample file %s does not exist.", sam))
        }
        processed_sample <- file.path(path_dir, paste0("processed_", sam))
        ## message(paste0('Processing sample:', sam))
        if (lib == "TAPS") {
            system(
                paste0(
                  "cat ", sample_path, " | awk '{if ($5 + $6 >= ", 
                  thresholds, paste0(
                    ") { $4 = int(($6 / ($5 + $6)) * 1000); ", 
                    "printf \"%s\\t%d\\t%d\\t%d\\t%d\\t%d\\n\",",
                    "$1, $2, $2+1, $4, $6, $5 }}' | ", 
                    "sort --buffer-size=5G -k1,1 -k2,2n "
                ),
                  paste0(" --temporary ", path_dir),
                  paste0(" --parallel ", parallel),
                  " | ", bedtools, " intersect -a - -b ", 
                exclude_bed_combined, 
                paste0(
                    " -v -sorted | awk ", 
                    "'{{printf \"%s\\t%d\\t%d\\t%d\\t%d\\t%d\\n\",",
                    "$1, $2, $3, $4, $5, $6}}'"
                ),
                  " > ", processed_sample
              ),
                intern = TRUE, ignore.stderr = TRUE
            )
        } else if (lib == "WGBS") {
            system(
                paste0(
                  "cat ", sample_path, " | awk '{if ($5 + $6 >= ", 
                  thresholds, paste0(
                    ") { $4 = int(($5 / ($5 + $6)) * 1000) ;", 
                    " ", "printf \"%s\\t%d\\t%d\\t%d\\t%d\\t%d\\n\", ",
                    "$1, $2, $2+1, $4, $5, $6 }}'", "|", 
                    " sort --buffer-size=5G  -k1,1 -k2,2n "
                ),
                  paste0(" --temporary ", path_dir),
                  paste0(" --parallel ", parallel),
                  " | ", bedtools, " intersect -a - -b ", 
                exclude_bed_combined, paste0(
                    " -v -sorted | awk '{{printf ", 
                    " \"%s\\t%d\\t%d\\t%d\\t%d\\t%d\\n\",",
                    "$1, $2, $3, $4, $5, $6}}'"
                ),
                  " > ", processed_sample
              ),
                intern = TRUE, ignore.stderr = TRUE
            )
            ## message('WGBS processing completed')
        }

        for (k in seq_along(annotation_file)) {
            map_processed_sample <- file.path(
              path_dir, paste0("map_processed_", output_file[k], sam))
            if (genome_type == "human") {
                annotation_file_k <- annotation_file[k]
                annotation_file_k_filter_sort_k <- file.path(
                  path_dir, paste0("filter_sorted_", 
                                   basename(annotation_file_k)))
                annotation_file_k_filter_sort_k_cut <- file.path(
                  path_dir, paste0("filter_sorted_cut3", 
                                   basename(annotation_file_k)))
                system(
                  paste0(
                    "cut -f 1-3 ", annotation_file_k_filter_sort_k, 
                    " >  ", annotation_file_k_filter_sort_k_cut
                )
              )
                system(
                  paste0(
                    "cat ", processed_sample, paste0(
                      " | awk '$1 ~ /^(chr1$|chr2$|chr3$|chr4$|chr5$|", 
                      "chr6$|chr7$|chr8$|chr9$|chr10$|chr11$|",
                      "chr12$|chr13$|chr14$|chr15$|chr16$|chr17$|chr18$|",
                      "chr19$|chr20$|chr21$|chr22$|chrX$)/ {print $0}' "
                  ),
                    " |  ", bedtools, " map -a ", 
                    annotation_file_k_filter_sort_k_cut,
                    " -b -  -c 5,6 -o sum,sum -null nan ", paste0(
                      " | awk ' {if ($4 != \"nan\" && $5 ", 
                      "!= \"nan\" && ($4 + $5) > 0)",
                      " { $4 = int(($4 / ($4 + $5)) * 1000 ) ; ", 
                      "printf ", "\"%s\\t%d\\t%d\\t%d\\n\", $1, $2, $3, $4 } ",
                      "else { printf", 
                      " \"%s\\t%d\\t%d\\tnan\\n\", $1, $2, $3 } } ' "
                  ),
                    ">", map_processed_sample
                ),
                  intern = TRUE, ignore.stderr = TRUE
              )

                system(paste0("sed -i 's/-2147483647/nan/g' ",
                              map_processed_sample))

                ## Build coverage matrix message('Building matrix...')
                build_coverage_matrix(
                  processed_sample = processed_sample, 
                  annotation_file_k_filter_sort_k_cut = 
                    annotation_file_k_filter_sort_k_cut,
                  map_processed_sample = file.path(
                    path_dir, paste0("map_coverage_processed_", 
                                     output_file[k], sam)),
                  bedtools = bedtools, genome_type = genome_type
              )
            } else if (genome_type == "mouse") {
                ## Similar process for mouse genome
                annotation_file_k <- annotation_file[k]
                annotation_file_k_filter_sort_k <- file.path(
                  path_dir, paste0("filter_sorted_", 
                                   basename(annotation_file_k)))
                annotation_file_k_filter_sort_k_cut <- file.path(
                  path_dir, paste0("filter_sorted_cut3", 
                                   basename(annotation_file_k)))

                system(
                  paste0(
                    "cut -f 1-3 ", annotation_file_k_filter_sort_k, 
                    " >  ", annotation_file_k_filter_sort_k_cut
                )
              )

                # message('Mapping methylation data to mouse annotation
                # regions...')
                system(
                  paste0(
                    "cat ", processed_sample, paste0(
                      " | awk '$1 ~ /^(chr1$|chr2$|chr3$|chr4$|chr5$|", 
                      "chr6$|chr7$|chr8$|chr9$|chr10$|chr11$|chr12$|chr13$|",
                      "chr14$|chr15$|chr16$|chr17$|chr18$|chr19$|", 
                      "chrX$)/ {print $0}' ",
                      " |  "
                  ),
                    bedtools, " map -a ", 
                  annotation_file_k_filter_sort_k_cut, 
                  " -b -  -c 5,6 -o sum,sum -null nan ",
                    paste0(
                      " | awk ' {if ($4 != \"nan\" && $5 !=", 
                      " \"nan\" && ($4 + $5) > 0) ",
                      " { $4 = int(($4 / ($4 + $5)) * 1000 ) ; ", 
                      "printf ", "\"%s\\t%d\\t%d\\t%d\\n\",",
                      " $1, $2, $3, $4 } else ", 
                      "{ printf \"%s\\t%d\\t%d\\tnan\\n\", $1, $2, $3 } } ' "
                  ),
                    " | awk '{if ($4 >= 0) {print $0}}' ", ">", 
                  map_processed_sample
                ),
                  intern = TRUE, ignore.stderr = TRUE
              )

                ## message('Cleaning mouse output data...')
                system(paste0("sed -i 's/-2147483647/nan/g' ", 
                              map_processed_sample))

                # Build coverage matrix for mouse message('Building mouse
                # coverage matrix...')
                build_coverage_matrix(
                  processed_sample = processed_sample, 
                  annotation_file_k_filter_sort_k_cut = 
                    annotation_file_k_filter_sort_k_cut,
                  map_processed_sample = file.path(
                    path_dir, paste0("map_coverage_processed_", 
                                     output_file[k], sam)),
                  bedtools = bedtools, genome_type = genome_type
              )
            }
        }
    }



    message("\n")
    message("Merging all samples for each annotation file...")
    ## Merge individual window files for each annotation file
    pb <- txtProgressBar(
        min = 0, max = length(annotation_file),
        style = 3
    )
    for (k in seq_along(annotation_file)) {
        # message(sprintf('Processing annotation file %d of %d...', k,
        # length(annotation_file)))

        # Merge methylation matrices message('Merging methylation matrices...')
        merge_matrix(
            processed_samples = paste0(path_dir, "/map_processed*"),
            output_file_dir = file.path(path_dir, "tmp.bedGraph"),
            output_file_dir_merge = file.path(path_dir, output_file[k]),
            bedtools = bedtools, samples = samples, path_dir = path_dir
        )

        # Merge coverage matrices message('Merging coverage matrices...')
        merge_matrix(
            processed_samples = paste0(path_dir, "/map_coverage_processed*"),
            output_file_dir = file.path(path_dir, "tmp_coverage.bedGraph"),
            output_file_dir_merge = file.path(
              path_dir, paste0("cov_", output_file[k])),
            bedtools = bedtools, samples = samples, path_dir = path_dir
        )

        # message(sprintf('Completed merging matrices for annotation file %d',
        # k))
        setTxtProgressBar(pb, k)
    }
    close(pb)

    # message('Successfully merged all sample matrices')
    message("\n")
    message(">>>>>>>>>> Count matrix build completed <<<<<<<<<<\n")
    # Clean up temporary files message('Cleaning up temporary files...')
    cleanup_cmd <- paste0(
        "rm ", file.path(path_dir, paste0("map_*")),
        " ", file.path(path_dir, paste0("filter_sorted_cut3*")),
        " ", file.path(path_dir, paste0("tmp*")),
        " ", file.path(path_dir, paste0("process*")),
        " ", exclude_bed_combined
    )

    tryCatch(
        {
            system(cleanup_cmd, intern = TRUE, ignore.stderr = TRUE)
            # message('Successfully cleaned up temporary files')
        }, error = function(e) {
            message("Failed to clean up some temporary files: ", e$message)
        }
    )

    # Calculate and report execution time
    end_time <- Sys.time()
    elapsed_time <- sprintf("%.2f seconds", as.numeric(
      difftime(end_time, start_time, units = "secs")))

    # message('\nCount matrix build completed: ')
    message("Location: ", path_dir)
    message(
        sprintf(
            "1. Matrix file\n   File name: %s", 
            stringi::stri_join(output_file, collapse = " ")
        )
    )
    message(
        sprintf(
            "2. Coverage matrix file\n   File name: %s", stringi::stri_join(
                paste0("cov_", output_file),
                collapse = ", "
            )
        )
    )
    message(
        sprintf(
            "3. Annotation file\n   File name: %s", stringi::stri_join(
                paste0("filter_sorted_", basename(annotation_file)),
                collapse = ", "
            )
        )
    )
    # message('\n========================================')
    message("\n>>>>>>>>>> methylTracer execution completed <<<<<<<<<<")
}

## build_count_matrix step1
chrom_pat_fun <- function(genome_type = NULL) {
    chrom_pattern <- if (genome_type == "human") {
        paste0(
            "^(chr1$|chr2$|chr3$|chr4$|chr5$|chr6$|chr7$|chr8$|chr9$|chr10$|", 
            "chr11$|chr12$|chr13$|chr14$|chr15$|chr16$|chr17$|chr18$|",
            "chr19$|chr20$|chr21$|chr22$|chrX$)"
        )
    } else if (genome_type == "mouse") {
        paste0(
            "^(chr1$|chr2$|chr3$|chr4$|chr5$|chr6$|chr7$|chr8$|chr9$|", 
            "chr10$|chr11$|chr12$|chr13$|chr14$|chr15$|chr16$|chr17$|",
            "chr18$|chr19$|chrX$)"
        )
    } else {
        message("Invalid genome type. Must be 'human' or 'mouse'")
    }
    return(chrom_pattern)
}
## build_count_matrix step2
ann_sort_fun <- function(
    annotation_file = annotation_file, path_dir = path_dir, 
    valid_chroms = valid_chroms,
    parallel
) {
    for (annotation_file_k in annotation_file) {
        annotation_file_k_filter_sort <- file.path(
          path_dir, paste0("filter_sorted_", basename(annotation_file_k)))
        system_result <- system(
            paste0(
                "awk '$1 ~ /^(", valid_chroms, ")/ {print}' ", 
                annotation_file_k,
                "| sort -k1,1 -k2,2n ", paste0(" --temporary ", path_dir),
                paste0(" --parallel ", parallel),
                " > ", annotation_file_k_filter_sort
            ),
            intern = TRUE, ignore.stderr = TRUE
        )
    }
}
## build_count_matrix step3
exc_sort_fun <- function(
    path_dir = path_dir, genome_type = genome_type, 
    exclude_bed = exclude_bed, parallel = parallel
) {
    exclude_bed_combined <- file.path(path_dir, 
                                      paste0("exclude_bed_combined.bed"))
    chrom_pattern <- chrom_pat_fun(genome_type = genome_type)

    system_result <- system(
        paste0(
            "cat ", stringi::stri_join(exclude_bed, collapse = " "),
            "| awk '$1 ~ /^(", chrom_pattern, 
            ")/ {print $0}' |cut -f1-3 | sort -k1,1 -k2,2n",
            paste0(" --temporary ", path_dir),
            paste0(" --parallel ", parallel),
            " > ", exclude_bed_combined
        )
    )
    message("Blacklisted sites are excluded from the analysis!!!")
    return(exclude_bed_combined)
}



build_coverage_matrix <- function(
    processed_sample = NULL, annotation_file_k_filter_sort_k_cut = NULL,
    map_processed_sample = NULL,
    bedtools = NULL, genome_type = NULL
) {
    # Define chromosome patterns based on genome type
    chrom_pattern <- chrom_pat_fun(genome_type = genome_type)
    # Build bedtools command
    cmd <- build_coverage_matrix_step1(
        processed_sample = processed_sample, 
        chrom_pattern = chrom_pattern, bedtools = bedtools,
        ann_fi_k_fi_sort_k_cut = annotation_file_k_filter_sort_k_cut, 
        map_processed_sample = map_processed_sample
    )
    tryCatch(
        {
            system(cmd, intern = TRUE, ignore.stderr = TRUE)
        }, error = function(e) {
            message("Failed to build coverage matrix: ", e$message)
        }
    )
}
## build_coverage_matrix step-1
build_coverage_matrix_step1 <- function(
    processed_sample = NULL, chrom_pattern = NULL, 
    bedtools = NULL, ann_fi_k_fi_sort_k_cut = NULL,
    map_processed_sample = NULL
) {
    cmd <- paste0(
        "cat ", processed_sample, " | awk '$1 ~ /", 
        chrom_pattern, "/ {print $0}' ",
        " | ", bedtools, " map -a ", 
        ann_fi_k_fi_sort_k_cut, 
        " -b - -c 5,6 -o sum,sum -null nan ",
        paste0(
            " | awk '{if ($4 != \"nan\" && $5 != \"nan\" && ($4 + $5) > 0) ", 
            " { $4 = int($4 + $5);printf \"%s\\t%d\\t%d\\t%d\\n\",",
            "$1, $2, $3, $4 } else { printf ", 
            "\"%s\\t%d\\t%d\\tnan\\n\", $1, $2, $3 } }' "
        ),
        "> ", map_processed_sample
    )
    return(cmd)
}
## merge_matrix
merge_matrix <- function(
    processed_samples = NULL, output_file_dir = NULL, 
    output_file_dir_merge = NULL,
    bedtools = NULL, samples = NULL, path_dir = NULL
) {
    merge_matrix_step1(
        bedtools = bedtools, samples = samples, 
        processed_samples = processed_samples,
        output_file_dir = output_file_dir
    )
    tmp1 <- file.path(path_dir, "cut1-3.tmp")
    tmp2 <- file.path(path_dir, "cut4-end.tmp")
    merge_matrix_step2(
        output_file_dir = output_file_dir, tmp1 = tmp1, 
        tmp2 = tmp2, output_file_dir_merge = output_file_dir_merge
    )
}
## merge_matrix step-1
merge_matrix_step1 <- function(
    bedtools = NULL, samples = NULL, processed_samples = NULL, 
    output_file_dir = NULL
) {
    tryCatch(
        {
            system(
                paste0(
                  "LC_ALL=C ", bedtools, " unionbedg -header -names ", 
                  stringi::stri_join(samples, collapse = " "),
                  " -filler nan -i ", processed_samples, " > ", output_file_dir
              ),
                intern = TRUE, ignore.stderr = TRUE
            )
        }, error = function(e) {
            message("Failed to merge coverage matrices: ", e$message)
        }
    )
}
## merge_matrix step-2
merge_matrix_step2 <- function(output_file_dir = NULL, 
                               tmp1 = NULL, tmp2 = NULL, 
                               output_file_dir_merge = NULL) {
    tryCatch(
        {
            system(
                paste0(
                  "cut -f 1-3 ", output_file_dir, paste0(
                    " | awk 'BEGIN {OFS=\"_\"} {print $1, $2, $3}' |", 
                    " sed '1s/.*/marker_name/' > "
                ),
                  tmp1
              ),
                intern = TRUE, ignore.stderr = TRUE
            )
            system(
                paste0("cut -f 4- ", output_file_dir, " > ", tmp2),
                intern = TRUE, ignore.stderr = TRUE
            )
            system(
                paste0("paste -d \"\t\" ", tmp1, " ", tmp2, " > ", 
                       output_file_dir_merge),
                intern = TRUE, ignore.stderr = TRUE
            )
        }, error = function(e) {
            message("Failed to process merged matrix: ", e$message)
        }, finally = {
            unlink(c(tmp1, tmp2))
        }
    )
}
