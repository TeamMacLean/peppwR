# Phase 2: Distribution Fitting
# Tests for fit_distributions() and plot.peppwr_fits

test_that("fit_distributions returns peppwr_fits object", {
  # Create minimal test data
  test_data <- tibble::tibble(
    peptide = rep(c("pep1", "pep2"), each = 20),
    condition = rep("control", 40),
    abundance = c(
      rgamma(20, shape = 2, rate = 0.1),
      rnorm(20, mean = 50, sd = 10)
    )
  )

  result <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  expect_s3_class(result, "peppwr_fits")
})

test_that("fit_distributions handles column name arguments", {
  test_data <- tibble::tibble(
    my_id = rep("pep1", 15),
    my_group = rep("A", 15),
    my_value = rgamma(15, shape = 2, rate = 0.1)
  )

  # Should work with custom column names
  result <- fit_distributions(
    test_data,
    id = "my_id",
    group = "my_group",
    value = "my_value"
  )

  expect_s3_class(result, "peppwr_fits")
})

test_that("fit_distributions fits multiple distributions", {
  test_data <- tibble::tibble(
    peptide = rep("pep1", 30),
    condition = rep("control", 30),
    abundance = rgamma(30, shape = 2, rate = 0.1)
  )

  result <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  # Should have fit results with multiple distributions attempted
  expect_true(length(result$fits) > 0)
  # Each fit result should have dist, loglik, aic columns
  expect_true(all(c("dist", "loglik", "aic") %in% names(result$fits[[1]])))
})

test_that("fit_distributions identifies best fit", {
  test_data <- tibble::tibble(
    peptide = rep("pep1", 30),
    condition = rep("control", 30),
    abundance = rgamma(30, shape = 2, rate = 0.1)
  )

  result <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  # Should have best-fit identified

  expect_true(length(result$best) == 1)
  expect_true(is.character(result$best))
})

test_that("fit_distributions stores original call", {
  test_data <- tibble::tibble(
    peptide = rep("pep1", 15),
    condition = rep("control", 15),
    abundance = rnorm(15, 100, 10)
  )

  result <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  expect_true(!is.null(result$call))
  expect_true(is.call(result$call))
})

test_that("plot.peppwr_fits returns ggplot object", {
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

  p <- plot(fits)

  expect_s3_class(p, "ggplot")
})

test_that("plot.peppwr_fits shows distribution summary", {
  test_data <- tibble::tibble(
    peptide = rep(c("pep1", "pep2"), each = 20),
    condition = rep("control", 40),
    abundance = c(
      rgamma(20, shape = 2, rate = 0.1),
      rnorm(20, mean = 50, sd = 10)
    )
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  p <- plot(fits)

  # Plot should have data
  expect_true(nrow(p$data) > 0)
})
