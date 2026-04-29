

if (FALSE) {
  
  output_dir <- "/home/jiazet/software/R_packages/DEVELOPMENT_r_packages/dev_dmr/methylTracer/methylTracer_personal/tmp_data/TAPS"
  output_file <- "test.h5"
  met <- build_met_obj(file.path(output_dir, output_file),
                       sample_name = "sample_name", marker_name = "V5"
  )
  
  pre_res <- pre_calldmrs(
    met = met, group_colname = "Group",
    case_group = "G1", control_group = "G2"
  )
  dmr_res <- calldmrs_turbo(
    met = met, p_threshold = 0.01, minCG=10L,
    case_group = "G1", ctrl_group = "G2"
  )
}