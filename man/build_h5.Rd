% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_h5.R
\name{build_h5}
\alias{build_h5}
\title{Build Standard HDF5 File for methylTracer Input}
\usage{
build_h5(
  sam_info = NULL,
  input_file = NULL,
  input_coverage = NULL,
  output_dir = NULL,
  output_file = NULL,
  annotation_file = NULL,
  chunk_size = 1000,
  level = 1
)
}
\arguments{
\item{sam_info}{path to file containing metadata for each sample.
The first column must be named 'sample_name'.}

\item{input_file}{Required; dataframe containing methylation signal.}

\item{input_coverage}{Optional; coverage dataframe
(default from 'methylTracer::build_count_matrix').}

\item{output_dir}{Path to the output directory where the HDF5
file will be saved.}

\item{output_file}{HDF5 filename (without extension) for saving
the methylation data.}

\item{annotation_file}{Genomic annotation file path, as in
\code{methylTracer::add_annotation}.}

\item{chunk_size}{Chunk size for reading data, default is 1e5.}

\item{level}{Compression level (1-9), default is 1.
Higher values provide better compression but may increase
processing time.}
}
\value{
Creates an HDF5 file with:
\itemize{
\item The methylation data matrix ('X' group).
\item Sample information ('obs' group).
\item Marker (genomic) information ('var' group).
\item Genomic annotation data.
}
}
\description{
This function creates an HDF5 file containing
methylation data and metadata.
}
\details{
The main methylation data matrix is stored under the
'X' group as a dense matrix in the HDF5 file.
}
\examples{
output_dir <- tempdir()

sam_info_df <- data.frame(
  sample_name = c(
    'sample_1.bed', 'sample_2.bed',
    'sample_3.bed', 'sample_4.bed'
  ),
  group = c(rep('case', 2), rep('ctrl', 2))
)
sam_info <- file.path(output_dir, 'sample_info.csv')

input_file_df <- data.frame(
  marker_name = c(
    'chr1_1000_2000', 'chr1_2000_3000',
    'chr1_3000_4000', 'chr1_4000_5000'
  ),
  sample_1.bed = c(100, 200, 600, 900),
  sample_2.bed = c(100, 200, 600, 900),
  sample_3.bed = c(200, 400, 700, 1000),
  sample_4.bed = c(300, 400, 700, 1000)
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
annotation_file <- file.path(
  output_dir,
  'annotation.bed'
)

output_file <- 'methylTracer_obj_test.h5'

# Remove the output file if it already exists
unlink(file.path(output_dir, output_file), recursive = TRUE)

# Write the sample information, input data, and annotation

write.csv(sam_info_df, sam_info, row.names = FALSE)
write.table(input_file_df, input_file, sep = '\t', row.names = FALSE)
write.table(annotation_file_df, annotation_file,
  sep = '\t', row.names = FALSE
)

# Run the function
build_h5(
  sam_info = sam_info,
  input_file = input_file,
  output_dir = output_dir,
  output_file = output_file,
  annotation_file = annotation_file
)

}
