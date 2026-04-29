


if (FALSE) {
  cl <- cluster_met_cells(
    met,
    n_clusters      = 2
  )
  
  pre_res <- pre_calldmrs(
    met          = met,
    group_colname = "cluster",
    case_group    = "C1",
    control_group = "C2"
  )
}