% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_annotation.R
\name{build_cpg}
\alias{build_cpg}
\title{Build CpG BED File for Base Sites}
\usage{
build_cpg(ref_genome = NULL, output = NULL, genome_type = "human")
}
\arguments{
\item{ref_genome}{A string specifying the annotation, e.g.
'BSgenome.Hsapiens.UCSC.hg38'.}

\item{output}{A string for the output file path, e.g.
'PATH/TO/cpg_sites.bed'.}

\item{genome_type}{A string specifying the genome type, e.g. 'human'
or 'mouse'.}
}
\value{
A BED file with CpG site coordinates saved to the specified output
directory.
}
\description{
Generates a BED file (0-based, half-open intervals) containing
coordinates of CpG sites based on the hg38 genome. Uses the
BSgenome.Hsapiens.UCSC.hg38 package to extract genomic data and create
a BED file suitable for downstream analysis. Currently supports only
'CpG' annotations.
}
\details{
The function builds an annotation for CpG sites using the
BSgenome.Hsapiens.UCSC.hg38 package. Coordinates are formatted into a
BED file.
}
\examples{

## Example usage (about 3 minutes to run):
output_dir <- tempdir()
## build_cpg(ref_genome = "BSgenome.Hsapiens.UCSC.hg38", 
##           output = file.path(output_dir, "cpg_sites.bed"), )

## $ head cpg_sites.bed
## Output:
## chr1    10468   10469   2       *
## chr1    10470   10471   2       *
## chr1    10483   10484   2       *
## chr1    10488   10489   2       *
}
