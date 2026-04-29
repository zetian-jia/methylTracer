
if (FALSE) {
  
  ncells <- 100
  u <- matrix(rpois(20000, 5), ncol=ncells)
  v <- u/100
  
  sce <- SingleCellExperiment(assays=list(counts=u, meth=v))
  sce
  
  out_h5 <- file.path(tempdir(), "dmr_level_met.h5")
  met_dmr <- sce_to_met(
    sce      = sce,
    h5_file  = out_h5,
    assay_name = "meth",
    overwrite  = TRUE
  )
  
  met_dmr
}