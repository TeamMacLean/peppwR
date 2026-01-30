# Phase B: Empirical Bootstrap Tests
# Tests for simulate_empirical, run_power_sim_empirical, and on_fit_failure="empirical"

# --- B1: simulate_empirical ---

test_that("simulate_empirical returns list with control and treatment", {
  set.seed(42)
  raw_data <- rgamma(50, shape = 2, rate = 0.1)

  result <- simulate_empirical(raw_data, n_per_group = 6, effect_size = 2)

  expect_type(result, "list")
  expect_named(result, c("control", "treatment"))
  expect_length(result$control, 6)
  expect_length(result$treatment, 6)
})

test_that("simulate_empirical applies effect size correctly", {
  set.seed(42)
  raw_data <- rgamma(100, shape = 2, rate = 0.1)

  result <- simulate_empirical(raw_data, n_per_group = 10, effect_size = 2)

  # Treatment should be approximately 2x control on average
  # With enough samples, the means should reflect the effect
  expect_gt(mean(result$treatment), mean(result$control) * 1.5)
})

test_that("simulate_empirical samples from raw_data", {
  set.seed(42)
  raw_data <- c(10, 20, 30, 40, 50)  # Limited unique values

  result <- simulate_empirical(raw_data, n_per_group = 100, effect_size = 1)

  # All control values should be from raw_data
  expect_true(all(result$control %in% raw_data))
})

test_that("simulate_empirical handles small datasets", {
  set.seed(42)
  raw_data <- c(5, 10, 15)  # Only 3 values

  result <- simulate_empirical(raw_data, n_per_group = 10, effect_size = 2)

  expect_length(result$control, 10)
  expect_length(result$treatment, 10)
})

# --- B2: run_power_sim_empirical ---

test_that("run_power_sim_empirical returns power between 0 and 1", {
  set.seed(42)
  raw_data <- rgamma(50, shape = 2, rate = 0.1)

  power <- run_power_sim_empirical(
    raw_data,
    n_per_group = 6,
    effect_size = 2,
    n_sim = 100
  )

  expect_type(power, "double")
  expect_true(power >= 0 && power <= 1)
})

test_that("run_power_sim_empirical returns high power for large effects", {
  set.seed(42)
  raw_data <- rgamma(100, shape = 2, rate = 0.1)

  power <- run_power_sim_empirical(
    raw_data,
    n_per_group = 10,
    effect_size = 5,  # Large effect
    n_sim = 100
  )

  expect_gt(power, 0.8)
})

test_that("run_power_sim_empirical returns low power for small effects", {
  set.seed(42)
  raw_data <- rgamma(100, shape = 2, rate = 0.1)

  power <- run_power_sim_empirical(
    raw_data,
    n_per_group = 3,
    effect_size = 1.1,  # Small effect
    n_sim = 100
  )

  expect_lt(power, 0.5)
})

test_that("run_power_sim_empirical respects alpha parameter", {
  set.seed(42)
  raw_data <- rgamma(50, shape = 2, rate = 0.1)

  power_05 <- run_power_sim_empirical(
    raw_data,
    n_per_group = 6,
    effect_size = 2,
    alpha = 0.05,
    n_sim = 100
  )

  power_01 <- run_power_sim_empirical(
    raw_data,
    n_per_group = 6,
    effect_size = 2,
    alpha = 0.01,
    n_sim = 100
  )

  # Power should be lower with stricter alpha
  expect_lte(power_01, power_05)
})

test_that("run_power_sim_empirical supports different tests", {
  set.seed(42)
  raw_data <- rgamma(50, shape = 2, rate = 0.1)

  power_wilcoxon <- run_power_sim_empirical(
    raw_data,
    n_per_group = 6,
    effect_size = 2,
    test = "wilcoxon",
    n_sim = 50
  )

  power_bootstrap <- run_power_sim_empirical(
    raw_data,
    n_per_group = 6,
    effect_size = 2,
    test = "bootstrap_t",
    n_sim = 50
  )

  expect_type(power_wilcoxon, "double")
  expect_type(power_bootstrap, "double")
})

# --- B3: on_fit_failure = "empirical" integration ---

test_that("power_analysis with on_fit_failure='empirical' works", {
  set.seed(42)

  # Create test data where at least some fits will fail
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:5), each = 20),
    condition = rep("control", 100),
    abundance = c(
      rgamma(20, shape = 2, rate = 0.1),
      rnorm(20, mean = 50, sd = 10),
      # Peptide with potentially problematic data for some distributions
      c(rep(1, 10), rep(2, 10)),  # Discrete values
      rgamma(20, shape = 3, rate = 0.2),
      rnorm(20, mean = 80, sd = 15)
    )
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  # Should work with empirical fallback
  result <- power_analysis(
    fits,
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    on_fit_failure = "empirical",
    n_sim = 50
  )

  expect_s3_class(result, "peppwr_power")
  expect_true(result$answer >= 0 && result$answer <= 1)
})

test_that("empirical mode uses raw data when fit fails", {
  set.seed(42)

  # Create data that will definitely fail distribution fitting
  test_data <- tibble::tibble(
    peptide = rep("pep1", 20),
    condition = rep("control", 20),
    abundance = rep(c(1, 2), 10)  # Only 2 unique values - fitting will likely fail
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  # With empirical fallback, we should still get a result
  result <- power_analysis(
    fits,
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    on_fit_failure = "empirical",
    n_sim = 50
  )

  expect_s3_class(result, "peppwr_power")
})

test_that("power_analysis.peppwr_fits with on_fit_failure='empirical' computes power for failed fits", {
  set.seed(42)

  # Create mix of good and bad data
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:3), each = 20),
    condition = rep("control", 60),
    abundance = c(
      rgamma(20, shape = 2, rate = 0.1),  # Good for gamma
      rep(c(100, 200), 10),               # Bad - will fail
      rlnorm(20, meanlog = 3, sdlog = 0.5) # Good for lognormal
    )
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
    on_fit_failure = "empirical",
    n_sim = 50
  )

  # Should analyze all peptides (none excluded)
  expect_s3_class(result, "peppwr_power")
  # The simulations should track analyzed peptides
  expect_true("n_analyzed" %in% names(result$simulations) ||
                "peptide_power" %in% names(result$simulations))
})
