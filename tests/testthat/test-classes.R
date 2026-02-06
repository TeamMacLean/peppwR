# Phase 1: Class Foundation
# Tests for S3 class constructors and print methods

# --- peppwr_fits class ---

test_that("new_peppwr_fits creates valid object", {

# Mock data that would come from fit_distributions()
  mock_data <- tibble::tibble(
    id = c("pep1", "pep2"),
    group = c("A", "A"),
    data = list(
      tibble::tibble(value = rnorm(10, mean = 100, sd = 10)),
      tibble::tibble(value = rnorm(10, mean = 150, sd = 15))
    )
  )

  mock_fits <- list(
    tibble::tibble(dist = "Normal", loglik = -30, aic = 64),
    tibble::tibble(dist = "Gamma", loglik = -28, aic = 60)
  )

  mock_best <- c("Normal", "Gamma")
  mock_call <- quote(fit_distributions(data, id = "id", group = "group", value = "value"))

  result <- new_peppwr_fits(
    data = mock_data,
    fits = mock_fits,
    best = mock_best,
    call = mock_call
  )

  expect_s3_class(result, "peppwr_fits")
  expect_true(all(c("data", "fits", "best", "call") %in% names(result)))
})

test_that("validate_peppwr_fits catches invalid objects", {
  # Missing required components
  bad_obj <- structure(list(data = NULL), class = "peppwr_fits")
  expect_error(validate_peppwr_fits(bad_obj))

  # Wrong class
  expect_error(validate_peppwr_fits(list(a = 1)))
})

test_that("print.peppwr_fits outputs summary", {
  # Create a valid peppwr_fits object
  fits <- new_peppwr_fits(
    data = tibble::tibble(id = "pep1", group = "A", data = list(tibble::tibble(value = 1:10))),
    fits = list(tibble::tibble(dist = "Normal", loglik = -20, aic = 44)),
    best = "Normal",
    call = quote(fit_distributions(x))
  )

  output <- capture.output(print(fits))

  # Should show class name
  expect_true(any(grepl("peppwr_fits", output)))
  # Should show count of peptides
  expect_true(any(grepl("1 peptide", output, ignore.case = TRUE)))
  # Should mention best fit distribution
  expect_true(any(grepl("Normal", output)))
})

# --- peppwr_power class ---

test_that("new_peppwr_power creates valid object", {
  result <- new_peppwr_power(
    mode = "aggregate",
    question = "power",
    answer = 0.82,
    simulations = list(n_sim = 1000, results = rbinom(1000, 1, 0.82)),
    params = list(
      effect_size = 2,
      n_per_group = 6,
      alpha = 0.05,
      test = "wilcoxon"
    ),
    call = quote(power_analysis(effect_size = 2, n_per_group = 6))
  )

  expect_s3_class(result, "peppwr_power")
  expect_true(all(c("mode", "question", "answer", "simulations", "params", "call") %in% names(result)))
  expect_equal(result$mode, "aggregate")
  expect_equal(result$question, "power")
})

test_that("validate_peppwr_power catches invalid objects", {
  # Invalid mode
  bad_obj <- structure(
    list(mode = "invalid", question = "power", answer = 0.8,
         simulations = list(), params = list(), call = NULL),
    class = "peppwr_power"
  )
  expect_error(validate_peppwr_power(bad_obj))

  # Invalid question
  bad_obj2 <- structure(
    list(mode = "aggregate", question = "invalid", answer = 0.8,
         simulations = list(), params = list(), call = NULL),
    class = "peppwr_power"
  )
  expect_error(validate_peppwr_power(bad_obj2))
})

test_that("print.peppwr_power shows clear answer for aggregate mode", {
  power_result <- new_peppwr_power(
    mode = "aggregate",
    question = "power",
    answer = 0.82,
    simulations = list(n_sim = 1000),
    params = list(effect_size = 2, n_per_group = 6, alpha = 0.05, test = "wilcoxon"),
    call = quote(power_analysis(effect_size = 2, n_per_group = 6))
  )

  output <- capture.output(print(power_result))

  # Should clearly state the power
  expect_true(any(grepl("82%|0.82", output)))
  # Should mention sample size
  expect_true(any(grepl("6", output)))
  # Should mention effect size
  expect_true(any(grepl("2", output)))
})

test_that("print.peppwr_power shows clear answer for sample_size question", {
  power_result <- new_peppwr_power(
    mode = "aggregate",
    question = "sample_size",
    answer = 8,
    simulations = list(n_sim = 1000),
    params = list(effect_size = 2, target_power = 0.8, alpha = 0.05, test = "wilcoxon"),
    call = quote(power_analysis(effect_size = 2, target_power = 0.8, find = "sample_size"))
  )

  output <- capture.output(print(power_result))

  # Should clearly state recommended sample size
  expect_true(any(grepl("8", output)))
  expect_true(any(grepl("sample|per group", output, ignore.case = TRUE)))
})

# --- bayes_t print output ---

test_that("print.peppwr_power shows BF threshold for bayes_t test", {
  power_result <- new_peppwr_power(
    mode = "aggregate",
    question = "power",
    answer = 0.65,
    simulations = list(n_sim = 1000),
    params = list(effect_size = 2, n_per_group = 6, alpha = 0.05, test = "bayes_t"),
    call = quote(power_analysis(effect_size = 2, n_per_group = 6, test = "bayes_t"))
  )

  output <- capture.output(print(power_result))

  # Should show BF threshold, not significance level

  expect_true(any(grepl("BF.*3|Bayes|threshold", output, ignore.case = TRUE)))
  # Should NOT show "Significance level" for bayes_t
  expect_false(any(grepl("Significance level", output)))
})

test_that("print.peppwr_power shows significance level for non-Bayesian tests", {
  power_result <- new_peppwr_power(
    mode = "aggregate",
    question = "power",
    answer = 0.75,
    simulations = list(n_sim = 1000),
    params = list(effect_size = 2, n_per_group = 6, alpha = 0.05, test = "wilcoxon"),
    call = quote(power_analysis(effect_size = 2, n_per_group = 6))
  )

  output <- capture.output(print(power_result))

  # Should show significance level for non-Bayesian tests
  expect_true(any(grepl("Significance level|0\\.05", output)))
  # Should NOT show BF threshold
  expect_false(any(grepl("BF.*3", output)))
})

test_that("summary.peppwr_power shows BF threshold for bayes_t test", {
  power_result <- new_peppwr_power(
    mode = "aggregate",
    question = "power",
    answer = 0.65,
    simulations = list(n_sim = 1000),
    params = list(effect_size = 2, n_per_group = 6, alpha = 0.05, test = "bayes_t"),
    call = quote(power_analysis(effect_size = 2, n_per_group = 6, test = "bayes_t"))
  )

  output <- capture.output(print(summary(power_result)))

  # Should show BF threshold in summary
  expect_true(any(grepl("BF.*3|threshold", output, ignore.case = TRUE)))
  # Should NOT show Alpha for bayes_t
  expect_false(any(grepl("Alpha:", output)))
})
