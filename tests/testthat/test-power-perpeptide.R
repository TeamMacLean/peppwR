# Phase 5: Power Analysis - Per-Peptide Mode
# Tests for power_analysis.peppwr_fits() and on_fit_failure handling

test_that("power_analysis.peppwr_fits returns peppwr_power with per_peptide mode", {
  # Create test data
  test_data <- tibble::tibble(
    peptide = rep(c("pep1", "pep2", "pep3"), each = 20),
    condition = rep("control", 60),
    abundance = c(
      rgamma(20, shape = 2, rate = 0.1),
      rnorm(20, mean = 50, sd = 10),
      rlnorm(20, meanlog = 3, sdlog = 0.5)
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
    n_sim = 50
  )

  expect_s3_class(result, "peppwr_power")
  expect_equal(result$mode, "per_peptide")
})

test_that("power_analysis.peppwr_fits reports power distribution across peptides", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:5), each = 20),
    condition = rep("control", 100),
    abundance = c(
      rgamma(20, shape = 2, rate = 0.1),
      rnorm(20, mean = 50, sd = 10),
      rlnorm(20, meanlog = 3, sdlog = 0.5),
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

  result <- power_analysis(
    fits,
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    n_sim = 50
  )

  # Should have per-peptide power values in simulations
  expect_true("peptide_power" %in% names(result$simulations))
  expect_equal(length(result$simulations$peptide_power), 5)
})

test_that("power_analysis.peppwr_fits answer summarizes across peptides", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:3), each = 20),
    condition = rep("control", 60),
    abundance = c(
      rgamma(20, shape = 2, rate = 0.1),
      rnorm(20, mean = 50, sd = 10),
      rlnorm(20, meanlog = 3, sdlog = 0.5)
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
    target_power = 0.8,
    find = "sample_size",
    n_sim = 50
  )

  # Answer should indicate what proportion of peptides reach target
  expect_true("proportion_powered" %in% names(result$simulations) ||
                is.numeric(result$answer))
})

# --- on_fit_failure handling ---

test_that("on_fit_failure='exclude' skips peptides with failed fits", {
  # Create data where one peptide will likely fail fitting
  test_data <- tibble::tibble(
    peptide = rep(c("pep1", "pep_bad"), each = 20),
    condition = rep("control", 40),
    abundance = c(
      rgamma(20, shape = 2, rate = 0.1),
      rep(c(1, 1000), 10)  # Bimodal, hard to fit
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
    on_fit_failure = "exclude",
    n_sim = 50
  )

  # Should report how many were excluded
  expect_true("n_excluded" %in% names(result$simulations) ||
                "n_analyzed" %in% names(result$simulations))
})

test_that("on_fit_failure='lognormal' uses fallback distribution", {
  test_data <- tibble::tibble(
    peptide = rep("pep1", 20),
    condition = rep("control", 20),
    abundance = rgamma(20, shape = 2, rate = 0.1)
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  # Manually break the fit to test fallback
  fits$fits[[1]] <- tibble::tibble(dist = NA_character_, loglik = NA_real_, aic = NA_real_)
  fits$best[[1]] <- NA_character_

  result <- power_analysis(
    fits,
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    on_fit_failure = "lognormal",
    n_sim = 50
  )

  expect_s3_class(result, "peppwr_power")
})

test_that("on_fit_failure parameter validates input", {
  test_data <- tibble::tibble(
    peptide = rep("pep1", 20),
    condition = rep("control", 20),
    abundance = rgamma(20, shape = 2, rate = 0.1)
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  expect_error(
    power_analysis(
      fits,
      effect_size = 2,
      n_per_group = 6,
      find = "power",
      on_fit_failure = "invalid_option"
    ),
    "on_fit_failure"
  )
})

test_that("per-peptide mode works with find='sample_size'", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:3), each = 25),
    condition = rep("control", 75),
    abundance = c(
      rgamma(25, shape = 2, rate = 0.1),
      rnorm(25, mean = 50, sd = 10),
      rlnorm(25, meanlog = 3, sdlog = 0.5)
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
    target_power = 0.8,
    find = "sample_size",
    n_sim = 50
  )

  expect_s3_class(result, "peppwr_power")
  expect_equal(result$question, "sample_size")
})

# --- find='effect_size' tests ---

test_that("find='effect_size' requires n_per_group", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:3), each = 25),
    condition = rep("control", 75),
    abundance = c(
      rgamma(25, shape = 2, rate = 0.1),
      rnorm(25, mean = 50, sd = 10),
      rlnorm(25, meanlog = 3, sdlog = 0.5)
    )
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  expect_error(
    power_analysis(fits, target_power = 0.8, find = "effect_size"),
    "n_per_group"
  )
})

test_that("find='effect_size' requires target_power", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:3), each = 25),
    condition = rep("control", 75),
    abundance = c(
      rgamma(25, shape = 2, rate = 0.1),
      rnorm(25, mean = 50, sd = 10),
      rlnorm(25, meanlog = 3, sdlog = 0.5)
    )
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  expect_error(
    power_analysis(fits, n_per_group = 6, find = "effect_size"),
    "target_power"
  )
})

test_that("per-peptide mode works with find='effect_size'", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:3), each = 25),
    condition = rep("control", 75),
    abundance = c(
      rgamma(25, shape = 2, rate = 0.1),
      rnorm(25, mean = 50, sd = 10),
      rlnorm(25, meanlog = 3, sdlog = 0.5)
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
    n_per_group = 6,
    target_power = 0.8,
    find = "effect_size",
    n_sim = 50
  )

  expect_s3_class(result, "peppwr_power")
  expect_equal(result$mode, "per_peptide")
  expect_equal(result$question, "effect_size")
})

test_that("find='effect_size' returns effect_curve in simulations", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:3), each = 25),
    condition = rep("control", 75),
    abundance = c(
      rgamma(25, shape = 2, rate = 0.1),
      rnorm(25, mean = 50, sd = 10),
      rlnorm(25, meanlog = 3, sdlog = 0.5)
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
    n_per_group = 6,
    target_power = 0.8,
    find = "effect_size",
    n_sim = 50
  )

  expect_true("effect_curve" %in% names(result$simulations))
  expect_true("effect_size" %in% names(result$simulations$effect_curve))
  expect_true("proportion_powered" %in% names(result$simulations$effect_curve))
})

test_that("find='effect_size' answer is minimum effect achieving threshold", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:3), each = 25),
    condition = rep("control", 75),
    abundance = c(
      rgamma(25, shape = 2, rate = 0.1),
      rnorm(25, mean = 50, sd = 10),
      rlnorm(25, meanlog = 3, sdlog = 0.5)
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
    n_per_group = 6,
    target_power = 0.8,
    find = "effect_size",
    n_sim = 50
  )

  # Answer should be one of the searched effect values
  expect_true(result$answer %in% c(1.5, 2, 3, 5, 10))
})

test_that("proportion_threshold works for find='effect_size'", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:5), each = 25),
    condition = rep("control", 125),
    abundance = c(
      rgamma(25, shape = 2, rate = 0.1),
      rnorm(25, mean = 50, sd = 10),
      rlnorm(25, meanlog = 3, sdlog = 0.5),
      rgamma(25, shape = 3, rate = 0.15),
      rnorm(25, mean = 80, sd = 12)
    )
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  # Lower threshold should need same or smaller effect
  result_50 <- power_analysis(
    fits,
    n_per_group = 6,
    target_power = 0.8,
    find = "effect_size",
    proportion_threshold = 0.5,
    n_sim = 50
  )

  result_80 <- power_analysis(
    fits,
    n_per_group = 6,
    target_power = 0.8,
    find = "effect_size",
    proportion_threshold = 0.8,
    n_sim = 50
  )

  # Higher threshold needs same or larger effect
  expect_true(result_80$answer >= result_50$answer)
})

test_that("proportion_threshold works for find='sample_size'", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:5), each = 25),
    condition = rep("control", 125),
    abundance = c(
      rgamma(25, shape = 2, rate = 0.1),
      rnorm(25, mean = 50, sd = 10),
      rlnorm(25, meanlog = 3, sdlog = 0.5),
      rgamma(25, shape = 3, rate = 0.15),
      rnorm(25, mean = 80, sd = 12)
    )
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  # Lower threshold should need same or smaller N
  result_50 <- power_analysis(
    fits,
    effect_size = 2,
    target_power = 0.8,
    find = "sample_size",
    proportion_threshold = 0.5,
    n_sim = 50
  )

  result_80 <- power_analysis(
    fits,
    effect_size = 2,
    target_power = 0.8,
    find = "sample_size",
    proportion_threshold = 0.8,
    n_sim = 50
  )

  # Higher threshold needs same or larger N
  expect_true(result_80$answer >= result_50$answer)
})
