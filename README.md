
# methylTracer

<img src="vignettes/methylTracer.logo.png" width="200" align="right"/>

`methylTracer` is an R package for tile-based, HDF5-backed analysis of
single-cell and large-scale DNA methylation data.

DNA methylation is a key epigenetic modification involved in gene
regulation, development, and disease. With the advent of whole-genome
bisulfite sequencing (WGBS), bisulfite-free chemistries such as
[TAPS](https://www.nature.com/articles/s41587-019-0041-2), and a growing
number of single-cell and multi-omic protocols, researchers can now
profile methylomes at high resolution across many cells. These
technologies, however, generate very large and sparse datasets that are
not well served by traditional bulk-oriented tools.

`methylTracer` addresses this gap by aggregating per-CpG methylation
calls into region-level, coverage-weighted methylation values defined
over arbitrary genomic annotations (e.g. fixed tiles, promoters,
enhancers). The resulting sample-by-region matrix is stored in an
AnnData-style HDF5 file and accessed through a lightweight **met**
object, which keeps only indices and metadata in memory while streaming
methylation values from disk. This design enables analyses of
genome-wide datasets with thousands to tens of thousands of cells on
modest hardware.

On top of this representation, `methylTracer` provides an end-to-end
workflow for de-novo DMR discovery:

- `build_count_matrix()` – construct weighted methylation and coverage
  matrices from per-sample files.  
- `build_h5()` / `build_met_obj()` – create the HDF5 container and the
  on-disk met object.  
- `cluster_met_cells()` – perform PCA, UMAP, and Leiden clustering on
  region-level methylation.  
- `pre_calldmrs()` / `calldmrs_turbo()` – compute region-wise Welch
  t-tests and call de-novo DMRs with a C++-accelerated algorithm.  
- `dmrs_to_sce_memory()` – aggregate DMR-level methylation into a
  [`SingleCellExperiment`](https://www.bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
  object for downstream analysis with standard Bioconductor tools.

Together, these components provide a scalable framework for discovering
and interpreting differential methylation patterns in heterogeneous cell
populations.

## Installation

    if (!requireNamespace("devtools", quietly = TRUE)) {
       install.packages("devtools")
    }

    devtools::install_github("zetian-jia/methylTracer")
