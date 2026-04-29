
if (FALSE) {
  
  set.seed(123)
  
  ncells   <- 100
  nmarkers <- 200
  
  ## 1) 构造 count / meth 矩阵 -----------------------------------------
  # 这里用二项分布保证 0-1 之间的 methylation
  u <- matrix(
    rbinom(nmarkers * ncells, size = 100, prob = 0.5),
    ncol = ncells
  )
  v <- u / 100
  
  # 构造 marker_name: chr1_100_199, chr1_200_299, ...
  starts      <- seq(100, by = 100, length.out = nmarkers)
  ends        <- starts + 99
  marker_name <- paste0("chr1_", starts, "_", ends)
  
  rownames(u) <- marker_name
  rownames(v) <- marker_name
  
  # 构造 sample_name
  sample_name <- paste0("cell_", seq_len(ncells))
  colnames(u) <- sample_name
  colnames(v) <- sample_name
  
  ## 2) 构造 rowRanges & colData --------------------------------------
  rr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", nmarkers),
    ranges   = IRanges::IRanges(start = starts, end = ends)
  )
  mcols(rr)$marker_name <- marker_name
  
  cd <- S4Vectors::DataFrame(
    sample_name = sample_name,
    group       = rep(c("case", "ctrl"), each = ncells / 2)
  )
  rownames(cd) <- cd$sample_name
  
  ## 3) 构造 SingleCellExperiment -------------------------------------
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays    = list(counts = u, meth = v),
    rowRanges = rr,
    colData   = cd
  )
  
  ## 4) 调用 sce_to_met() 写 HDF5，转成 methylTracer -----------------
  h5_file <- tempfile(fileext = ".h5")
  
  met <- sce_to_met(
    sce         = sce,
    h5_file     = h5_file,
    assay_name  = "meth",
    overwrite   = TRUE
  )
  
  met
  
  sce <- met_to_sce(
    met = met
  )
  sce
}