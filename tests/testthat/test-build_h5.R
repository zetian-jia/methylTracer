
if (FALSE) {
  
  sam_info <- "/home/jiazet/software/R_packages/DEVELOPMENT_r_packages/dev_dmr/methylTracer/methylTracer_personal/tmp_data/TAPS/sam_info.csv"
  input_file <- "/home/jiazet/software/R_packages/DEVELOPMENT_r_packages/dev_dmr/methylTracer/methylTracer_personal/tmp_data/TAPS/met.matrix.CpG.txt"
  input_coverage <- "/home/jiazet/software/R_packages/DEVELOPMENT_r_packages/dev_dmr/methylTracer/methylTracer_personal/tmp_data/TAPS/cov.met.matrix.CpG.txt"
  output_dir <- "/home/jiazet/software/R_packages/DEVELOPMENT_r_packages/dev_dmr/methylTracer/methylTracer_personal/tmp_data/TAPS/"
  output_file <- "test.h5"
  annotation_file <- "/home/jiazet/software/R_packages/DEVELOPMENT_r_packages/dev_dmr/methylTracer/methylTracer_personal/tmp_data/TAPS/3.ann.bed"
  build_h5(
    sam_info = sam_info, 
    input_file = input_file, 
    input_coverage = input_coverage, 
    output_dir = output_dir,
    output_file = output_file, 
    annotation_file = annotation_file, chunk_size = 1e+06,
    level = 0
  ) 
}