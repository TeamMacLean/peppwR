# Phase 8: Polish - Summary Methods
# Tests for summary.peppwr_fits and summary.peppwr_power

test_that("summary.peppwr_fits returns summary object", {
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

  s <- summary(fits)

  # Should return a summary object or data frame
expect_true(is.list(s) || is.data.frame(s))
})

test_that("summary.peppwr_fits includes distribution counts", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:10), each = 20),
    condition = rep("control", 200),
    abundance = c(
      rgamma(40, shape = 2, rate = 0.1),
      rnorm(60, mean = 50, sd = 10),
      rlnorm(100, meanlog = 3, sdlog = 0.5)
    )
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  s <- summary(fits)

  # Should have distribution breakdown
  expect_true("best_dist_counts" %in% names(s) ||
                any(grepl("dist|distribution", names(s), ignore.case = TRUE)))
})

test_that("summary.peppwr_fits includes fit statistics", {
  test_data <- tibble::tibble(
    peptide = rep("pep1", 30),
    condition = rep("control", 30),
    abundance = rgamma(30, shape = 2, rate = 0.1)
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  s <- summary(fits)

  # Should have AIC or loglik info
  expect_true(any(grepl("aic|loglik|fit", names(s), ignore.case = TRUE)) ||
                "fit_summary" %in% names(s))
})

test_that("summary.peppwr_power returns summary object", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    n_sim = 100
  )

  s <- summary(result)

  expect_true(is.list(s) || is.data.frame(s))
})

test_that("summary.peppwr_power includes simulation details", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    n_sim = 100
  )

  s <- summary(result)

  # Should include n_sim
  expect_true("n_sim" %in% names(s) ||
                any(grepl("simulation", names(s), ignore.case = TRUE)))
})

test_that("summary.peppwr_power includes confidence interval", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    n_per_group = 10,
    find = "power",
    n_sim = 500
  )

  s <- summary(result)

  # Should have CI for power estimate
  expect_true(any(grepl("ci|conf|interval|lower|upper", names(s), ignore.case = TRUE)))
})

test_that("summary.peppwr_power per-peptide mode includes peptide breakdown", {
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

  s <- summary(result)

  # Should include per-peptide info
  expect_true(any(grepl("peptide|per_peptide", names(s), ignore.case = TRUE)) ||
                "peptide_summary" %in% names(s))
})

test_that("print method exists for summary objects", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    n_sim = 50
  )

  s <- summary(result)

  # Should be printable without error
  expect_output(print(s))
})
