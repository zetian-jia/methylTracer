# Test suite for simulate_scBS_data()

test_that("simulate_scBS_data creates correct number of files", {
  skip_on_cran()

  sim <- simulate_scBS_data(
    n_cells = 10,
    n_sites = 1000,
    sparsity = 0.8,
    n_dmrs = 5,
    output_dir = tempfile("test_sim"),
    verbose = FALSE
  )

  # Check metadata
  expect_equal(nrow(sim$metadata), 10)
  expect_true(all(file.exists(sim$metadata$file_path)))

  # Check groups
  expect_true(all(sim$metadata$group %in% c("Control", "Treatment")))

  # Check DMRs
  expect_equal(nrow(sim$dmr_regions), 5)

  # Cleanup
  unlink(sim$output_dir, recursive = TRUE)
})


test_that("simulate_scBS_data respects sparsity parameter", {
  skip_on_cran()

  # High sparsity test
  sim_sparse <- simulate_scBS_data(
    n_cells = 20,
    n_sites = 10000,
    sparsity = 0.95, # Very sparse
    n_dmrs = 10,
    output_dir = tempfile("test_sparse"),
    verbose = FALSE
  )

  # Calculate actual sparsity
  total_possible <- 20 * 10000
  total_covered <- sum(sim_sparse$metadata$n_covered_sites)
  actual_sparsity <- 1 - (total_covered / total_possible)

  # Should be close to target (within 5%)
  expect_true(abs(actual_sparsity - 0.95) < 0.05)

  # Cleanup
  unlink(sim_sparse$output_dir, recursive = TRUE)
})


test_that("DMRs have correct methylation patterns", {
  skip_on_cran()

  sim <- simulate_scBS_data(
    n_cells = 50,
    n_sites = 5000,
    sparsity = 0.7,
    n_dmrs = 10,
    dmr_width = 10,
    dmr_meth_ctrl = 0.9,
    dmr_meth_treat = 0.1,
    output_dir = tempfile("test_dmr"),
    seed = 42,
    verbose = FALSE
  )

  # Check DMR specifications
  expect_true(all(sim$dmr_regions$true_meth_ctrl == 0.9))
  expect_true(all(sim$dmr_regions$true_meth_treat == 0.1))
  expect_true(all(sim$dmr_regions$n_sites == 10))

  # Read one control and one treatment cell
  ctrl_cell <- sim$metadata[sim$metadata$group == "Control", ][1, ]
  treat_cell <- sim$metadata[sim$metadata$group == "Treatment", ][1, ]

  ctrl_data <- data.table::fread(ctrl_cell$file_path)
  treat_data <- data.table::fread(treat_cell$file_path)

  expect_true(nrow(ctrl_data) > 0)
  expect_true(nrow(treat_data) > 0)

  # Cleanup
  unlink(sim$output_dir, recursive = TRUE)
})


test_that("simulate_scBS_data validates input parameters", {
  expect_error(
    simulate_scBS_data(n_cells = 1, verbose = FALSE),
    "'n_cells' must be an integer >= 2"
  )

  expect_error(
    simulate_scBS_data(n_sites = 100, verbose = FALSE),
    "'n_sites' must be an integer >= 1000"
  )

  expect_error(
    simulate_scBS_data(sparsity = 1.5, verbose = FALSE),
    "'sparsity' must be in"
  )

  expect_error(
    simulate_scBS_data(group_ratio = 1.0, verbose = FALSE),
    "'group_ratio' must be in"
  )
})


test_that("output formats work correctly", {
  skip_on_cran()

  # Test COV format
  sim_cov <- simulate_scBS_data(
    n_cells = 5,
    n_sites = 1000,
    output_dir = tempfile("test_cov"),
    output_format = "cov",
    verbose = FALSE
  )

  # Read one file and check structure
  cov_data <- data.table::fread(sim_cov$metadata$file_path[1])
  expect_equal(ncol(cov_data), 4) # chr, pos, count_meth, count_unmeth

  unlink(sim_cov$output_dir, recursive = TRUE)

  # Test BED format
  sim_bed <- simulate_scBS_data(
    n_cells = 5,
    n_sites = 1000,
    output_dir = tempfile("test_bed"),
    output_format = "bed",
    verbose = FALSE
  )

  bed_data <- data.table::fread(sim_bed$metadata$file_path[1])
  expect_equal(ncol(bed_data), 6) # chr, start, end, name, count_meth, count_unmeth

  unlink(sim_bed$output_dir, recursive = TRUE)
})


test_that("simulation is reproducible with seed", {
  skip_on_cran()

  sim1 <- simulate_scBS_data(
    n_cells = 10,
    n_sites = 1000,
    seed = 999,
    output_dir = tempfile("test_seed1"),
    verbose = FALSE
  )

  sim2 <- simulate_scBS_data(
    n_cells = 10,
    n_sites = 1000,
    seed = 999,
    output_dir = tempfile("test_seed2"),
    verbose = FALSE
  )

  # Same seed should produce same group assignments
  expect_equal(sim1$metadata$group, sim2$metadata$group)

  # Same DMR locations
  expect_equal(sim1$dmr_regions$chr, sim2$dmr_regions$chr)
  expect_equal(sim1$dmr_regions$start, sim2$dmr_regions$start)

  # Cleanup
  unlink(sim1$output_dir, recursive = TRUE)
  unlink(sim2$output_dir, recursive = TRUE)
})


test_that("CpG universe is properly formatted", {
  skip_on_cran()

  sim <- simulate_scBS_data(
    n_cells = 5,
    n_sites = 5000,
    output_dir = tempfile("test_cpg"),
    verbose = FALSE
  )

  # Check CpG universe structure
  expect_true("chr" %in% names(sim$cpg_universe))
  expect_true("pos" %in% names(sim$cpg_universe))
  expect_true("site_id" %in% names(sim$cpg_universe))

  # All positions should be unique
  expect_equal(nrow(sim$cpg_universe), length(unique(sim$cpg_universe$site_id)))

  # Positions should be sorted within each chromosome
  cpg_dt <- data.table::as.data.table(sim$cpg_universe)
  is_sorted <- all(
    cpg_dt[, .SD[, all(diff(pos) >= 0)], by = chr]$V1
  )
  expect_true(is_sorted)

  # Cleanup
  unlink(sim$output_dir, recursive = TRUE)
})


test_that("print and summary methods work", {
  skip_on_cran()

  sim <- simulate_scBS_data(
    n_cells = 10,
    n_sites = 1000,
    output_dir = tempfile("test_print"),
    verbose = FALSE
  )

  # Test print method
  expect_output(print(sim), "scBS Simulation Object")
  expect_output(print(sim), "Cells: 10")

  # Test summary method
  expect_output(summary(sim), "Simulation Parameters")
  expect_output(summary(sim), "DMR Summary")

  # Cleanup
  unlink(sim$output_dir, recursive = TRUE)
})


test_that("coverage statistics are realistic", {
  skip_on_cran()

  mean_cov <- 15

  sim <- simulate_scBS_data(
    n_cells = 30,
    n_sites = 5000,
    mean_coverage = mean_cov,
    output_dir = tempfile("test_coverage"),
    verbose = FALSE
  )

  # Mean coverage should be close to specified value
  expect_true(abs(mean(sim$metadata$mean_coverage) - mean_cov) < 5)

  # All cells should have some coverage
  expect_true(all(sim$metadata$n_covered_sites > 0))

  # Cleanup
  unlink(sim$output_dir, recursive = TRUE)
})
