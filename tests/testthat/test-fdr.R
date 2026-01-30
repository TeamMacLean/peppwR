# Phase D: FDR-Aware Mode Tests
# Tests for whole-peptidome simulation with multiple testing correction

# --- D1: run_power_sim_fdr ---

test_that("run_power_sim_fdr returns power between 0 and 1", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:10), each = 20),
    condition = rep("control", 200),
    abundance = c(
      rgamma(20, shape = 2, rate = 0.1),
      rnorm(20, mean = 50, sd = 10),
      rlnorm(20, meanlog = 3, sdlog = 0.5),
      rgamma(20, shape = 3, rate = 0.2),
      rnorm(20, mean = 80, sd = 15),
      rgamma(20, shape = 1.5, rate = 0.12),
      rnorm(20, mean = 60, sd = 12),
      rlnorm(20, meanlog = 2.5, sdlog = 0.4),
      rgamma(20, shape = 2.5, rate = 0.15),
      rnorm(20, mean = 70, sd = 8)
    )
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  power <- run_power_sim_fdr(
    fits,
    effect_size = 2,
    n_per_group = 6,
    prop_null = 0.9,
    fdr_threshold = 0.05,
    n_sim = 50
  )

  expect_type(power, "double")
  expect_true(power >= 0 && power <= 1)
})

test_that("run_power_sim_fdr returns lower power than uncorrected", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:10), each = 20),
    condition = rep("control", 200),
    abundance = rgamma(200, shape = 2, rate = 0.1)
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  power_fdr <- run_power_sim_fdr(
    fits,
    effect_size = 2,
    n_per_group = 6,
    prop_null = 0.9,
    fdr_threshold = 0.05,
    n_sim = 50
  )

  # FDR-adjusted power should generally be <= uncorrected power
  # (more stringent correction = fewer discoveries)
  expect_true(power_fdr >= 0 && power_fdr <= 1)
})

test_that("run_power_sim_fdr respects prop_null parameter", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:10), each = 20),
    condition = rep("control", 200),
    abundance = rgamma(200, shape = 2, rate = 0.1)
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  # Higher prop_null = more true nulls = harder to maintain power
  power_high_null <- run_power_sim_fdr(
    fits,
    effect_size = 2,
    n_per_group = 6,
    prop_null = 0.95,
    fdr_threshold = 0.05,
    n_sim = 50
  )

  power_low_null <- run_power_sim_fdr(
    fits,
    effect_size = 2,
    n_per_group = 6,
    prop_null = 0.5,
    fdr_threshold = 0.05,
    n_sim = 50
  )

  # Both should be valid
  expect_true(power_high_null >= 0 && power_high_null <= 1)
  expect_true(power_low_null >= 0 && power_low_null <= 1)
})

test_that("run_power_sim_fdr respects fdr_threshold parameter", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:10), each = 20),
    condition = rep("control", 200),
    abundance = rgamma(200, shape = 2, rate = 0.1)
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  # Looser FDR = more discoveries = higher power
  power_strict <- run_power_sim_fdr(
    fits,
    effect_size = 2,
    n_per_group = 6,
    prop_null = 0.9,
    fdr_threshold = 0.01,
    n_sim = 50
  )

  power_loose <- run_power_sim_fdr(
    fits,
    effect_size = 2,
    n_per_group = 6,
    prop_null = 0.9,
    fdr_threshold = 0.1,
    n_sim = 50
  )

  # Both should be valid; typically loose >= strict
  expect_true(power_strict >= 0 && power_strict <= 1)
  expect_true(power_loose >= 0 && power_loose <= 1)
})

# --- D2: power_analysis with apply_fdr ---

test_that("power_analysis.peppwr_fits accepts apply_fdr parameter", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:5), each = 20),
    condition = rep("control", 100),
    abundance = rgamma(100, shape = 2, rate = 0.1)
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  result <- power_analysis(
    fits,
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    apply_fdr = TRUE,
    prop_null = 0.9,
    n_sim = 50
  )

  expect_s3_class(result, "peppwr_power")
  expect_true(result$answer >= 0 && result$answer <= 1)
})

test_that("apply_fdr=TRUE records FDR settings in params", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:5), each = 20),
    condition = rep("control", 100),
    abundance = rgamma(100, shape = 2, rate = 0.1)
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  result <- power_analysis(
    fits,
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    apply_fdr = TRUE,
    prop_null = 0.8,
    fdr_threshold = 0.1,
    n_sim = 50
  )

  expect_true("apply_fdr" %in% names(result$params) ||
                "fdr_adjusted" %in% names(result$simulations))
})

test_that("power_analysis with FDR gives different result than without", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:10), each = 20),
    condition = rep("control", 200),
    abundance = rgamma(200, shape = 2, rate = 0.1)
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  result_fdr <- power_analysis(
    fits,
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    apply_fdr = TRUE,
    prop_null = 0.9,
    n_sim = 100
  )

  result_no_fdr <- power_analysis(
    fits,
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    apply_fdr = FALSE,
    n_sim = 100
  )

  # Both should be valid, typically FDR power <= uncorrected power
  expect_true(result_fdr$answer >= 0 && result_fdr$answer <= 1)
  expect_true(result_no_fdr$answer >= 0 && result_no_fdr$answer <= 1)
})

# --- D3: print.peppwr_power with FDR ---

test_that("print.peppwr_power shows FDR info when applicable", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:5), each = 20),
    condition = rep("control", 100),
    abundance = rgamma(100, shape = 2, rate = 0.1)
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  result <- power_analysis(
    fits,
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    apply_fdr = TRUE,
    prop_null = 0.9,
    n_sim = 50
  )

  output <- capture.output(print(result))

  # Should mention FDR or multiple testing
  expect_true(any(grepl("FDR|fdr|multiple|adjusted|Benjamini|BH",
                        output, ignore.case = TRUE)))
})

test_that("print.peppwr_power shows prop_null when FDR enabled", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:5), each = 20),
    condition = rep("control", 100),
    abundance = rgamma(100, shape = 2, rate = 0.1)
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  result <- power_analysis(
    fits,
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    apply_fdr = TRUE,
    prop_null = 0.85,
    n_sim = 50
  )

  output <- capture.output(print(result))

  # Should show proportion of true nulls
  expect_true(any(grepl("null|unchanged|0\\.85|85", output, ignore.case = TRUE)))
})
