% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_met.R
\name{plot_dmrs}
\alias{plot_dmrs}
\title{Plot Differentially Methylated Regions (DMRs)}
\usage{
plot_dmrs(
  met_obj = NULL,
  pre_calldmrs_res = NULL,
  dmrs_gr = NULL,
  region = NULL,
  shift_up = 3,
  shift_down = 3,
  gene_col = "V4",
  regul_col = NULL
)
}
\arguments{
\item{met_obj}{A methylation object containing the data, which should
be structured with necessary slots for methylation information.}

\item{pre_calldmrs_res}{Pre-calculated DMRs results for plotting.}

\item{dmrs_gr}{GenomicRanges object containing the DMRs.}

\item{region}{A string specifying the genomic region to plot (e.g.,
"chr1:1000-5000").}

\item{shift_up}{The number of units to shift the region upward for better
visualization. Default is 3.}

\item{shift_down}{The number of units to shift the region downward.
Default is 3.}

\item{gene_col}{The column name that contains gene information for
annotation.}

\item{regul_col}{The column name that contains regulatory information for
annotation.}
}
\value{
A \code{ggplot2} object representing the DMR plot.
}
\description{
This function plots the differentially methylated regions (DMRs) based
on the specified methylation data, DMRs, and gene annotations. The plot
will show the regions of interest with optional shifting, coloring by gene
and regulatory information.
}
\details{
This function reads in methylation data, DMR information, and gene
annotations to generate a plot. It allows the shifting of the region
of interest up and down for improved visualization and can optionally
display gene and regulatory information.
}
\examples{
# example code
output_dir <- tempdir()

sam_info_df <- data.frame(
  sample_name = c(
    'sample_1.bed', 'sample_2.bed',
    'sample_3.bed', 'sample_4.bed',
    'sample_5.bed', 'sample_6.bed',
    'sample_7.bed', 'sample_8.bed'
  ),
  group = c(rep('case', 4), rep('ctrl', 4))
)
sam_info <- file.path(output_dir, 'sample_info.csv')

input_file_df <- data.frame(
  marker_name = c(
    'chr1_1000_2000', 'chr1_2000_3000',
    'chr1_3000_4000', 'chr1_4000_5000'
  ),
  sample_1.bed = c(100, 200, 600, 900),
  sample_2.bed = c(100, 200, 600, 900),
  sample_3.bed = c(100, 200, 600, 900),
  sample_4.bed = c(100, 200, 600, 900),
  sample_5.bed = c(800, 900, 100, 50),
  sample_6.bed = c(800, 900, 100, 50),
  sample_7.bed = c(800, 900, 100, 50),
  sample_8.bed = c(800, 900, 100, 50)
)
input_file <- file.path(output_dir, 'methylTracer_1kb.txt')

annotation_file_df <- data.frame(
  chr = rep('chr1', 4),
  start = c(1000, 2000, 3000, 4000),
  end = c(2000, 3000, 4000, 5000),
  SYMBOL = c('gene1', 'gene2', 'gene3', 'gene4'),
  marker_name = c(
    'chr1_1000_2000', 'chr1_2000_3000',
    'chr1_3000_4000', 'chr1_4000_5000'
  )
)
annotation_file <- file.path(output_dir, 'annotation.bed')

output_file <- 'methylTracer_obj_test.h5'

unlink(file.path(output_dir, output_file),
  recursive = TRUE
)

write.csv(sam_info_df, sam_info,
  row.names = FALSE
)
write.table(input_file_df, input_file,
  sep = '\t', row.names = FALSE
)
write.table(annotation_file_df, annotation_file,
  sep = '\t', row.names = FALSE
)

build_h5(
  sam_info = sam_info,
  input_file = input_file,
  output_dir = output_dir,
  output_file = output_file,
  annotation_file = annotation_file
)

met_obj <- build_met_obj(
  file.path(output_dir, output_file),
  sample_name = 'sample_name',
  marker_name = 'marker_name'
)

pre_res <- pre_calldmrs(
  met_obj = met_obj,
  group_colname = 'group',
  case_group = 'case',
  control_group = 'ctrl'
)
dmr_res <- calldmrs_turbo(
  met_obj = met_obj,
  p_threshold = 0.9,
  case_group = "case",
  ctrl_group = "ctrl"
)
plot_dmrs(
  met_obj = met_obj,
  pre_calldmrs_res = pre_res,
  dmrs_gr = dmr_res,
  region = "chr1:1000-4000",
  shift_up = 1,
  shift_down = 1,
  gene_col = "SYMBOL"
)
}
