% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_met.R
\name{plot_dens}
\alias{plot_dens}
\title{Plot Methylation Data (Density Plot)}
\usage{
plot_dens(met_obj = NULL, slot = NULL, bw = 0.5, vline = NULL)
}
\arguments{
\item{met_obj}{A methylation object containing the methylation data.
This object should include an HDF5 file with the relevant data.}

\item{slot}{A character string specifying the specific slot in the HDF5
file to read data from (e.g., 'obs/coverage_cells').
The slot represents the data of interest for the plot.}

\item{bw}{A numeric value specifying the bandwidth for kernel density
estimation (KDE). Default is 0.5. This controls the smoothness of
the density plot.}

\item{vline}{A numeric value specifying the x-position of a vertical
line to be drawn on the plot. Default is \code{NULL}, meaning no line
will be drawn.}
}
\value{
A \code{ggplot2} object representing the density plot. Visualizes
the distribution of the values from the specified \code{slot}.
}
\description{
This function creates a density plot (using kernel density estimation)
from methylation data stored in a given slot of the \code{met_obj} object.
It optionally includes a vertical line at a specified position.
}
\examples{
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

# Write the data to files
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

compute_qc_value(
  met_obj = met_obj
)

plot_dens(
  met_obj = met_obj,
  slot = 'obs/coverage_cells',
  vline = 2
)
}
