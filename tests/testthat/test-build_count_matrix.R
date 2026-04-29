
if (FALSE) {
  sam_info <- "/home/jiazet/software/R_packages/DEVELOPMENT_r_packages/dev_dmr/methylTracer/methylTracer_personal/tmp_data/TAPS/sam_info.csv"
  input_dir <- "/home/jiazet/software/R_packages/DEVELOPMENT_r_packages/dev_dmr/methylTracer/methylTracer_personal/tmp_data/TAPS"
  annotation_file <- "/home/jiazet/software/R_packages/DEVELOPMENT_r_packages/dev_dmr/methylTracer/methylTracer_personal/tmp_data/TAPS/annotation_file.bed"
  exclude_bed <- "/home/jiazet/software/R_packages/DEVELOPMENT_r_packages/dev_dmr/methylTracer/methylTracer_personal/tmp_data/TAPS/exclude_bed"
  parallel = 1
  bedtools <- "/home/jiazet/miniconda3/envs/taps/bin/bedtools"
  lib <- "TAPS"
  output_file <- "met.matrix.CpG.txt"
  thresholds = 1
  genome_type = "human"
  limitFiles=200000
  bedSort="/home/jiazet/miniconda3/envs/taps/bin/bedSort"
  
  build_count_matrix(sam_info=sam_info,
                     input_dir=input_dir,
                     annotation_file=annotation_file,
                     output_file=output_file,
                     thresholds=thresholds,
                     bedtools=bedtools,
                     bedSort=bedSort,
                     lib=lib,
                     parallel=parallel,
                     exclude_bed=exclude_bed,
                     genome_type=genome_type,
                     limitFiles=limitFiles
  )
}