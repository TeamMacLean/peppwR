# Phase 7: Additional Statistical Tests
# Tests for bootstrap_t, bayes_t, and rankprod tests

# --- test_bootstrap_t ---

test_that("test_bootstrap_t returns p-value", {
  control <- rnorm(10, mean = 100, sd = 10)
  treatment <- rnorm(10, mean = 100, sd = 10)

  result <- test_bootstrap_t(control, treatment, n_boot = 500)

  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("test_bootstrap_t detects significant difference", {
  set.seed(42)
  control <- rnorm(15, mean = 100, sd = 10)
  treatment <- rnorm(15, mean = 150, sd = 10)  # Large effect

  result <- test_bootstrap_t(control, treatment, n_boot = 1000)

  expect_true(result < 0.05)
})

test_that("test_bootstrap_t handles small samples", {
  control <- rnorm(5, mean = 100, sd = 10)
  treatment <- rnorm(5, mean = 120, sd = 10)

  result <- test_bootstrap_t(control, treatment, n_boot = 500)

  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("test_bootstrap_t n_boot parameter works", {
  control <- rnorm(10, mean = 100, sd = 10)
  treatment <- rnorm(10, mean = 110, sd = 10)

  # Should run with different n_boot values
  result_500 <- test_bootstrap_t(control, treatment, n_boot = 500)
  result_100 <- test_bootstrap_t(control, treatment, n_boot = 100)

  expect_true(is.numeric(result_500))
  expect_true(is.numeric(result_100))
})

# --- test_bayes_t ---

test_that("test_bayes_t returns Bayes factor", {
  control <- rnorm(10, mean = 100, sd = 10)
  treatment <- rnorm(10, mean = 100, sd = 10)

  result <- test_bayes_t(control, treatment)

  expect_true(is.numeric(result))
  expect_true(result > 0)  # Bayes factors are positive
})

test_that("test_bayes_t gives high BF for large effect", {
  set.seed(42)
  control <- rnorm(20, mean = 100, sd = 10)
  treatment <- rnorm(20, mean = 150, sd = 10)

  result <- test_bayes_t(control, treatment)

  expect_true(result > 3)  # Substantial evidence
})

test_that("test_bayes_t gives low BF for null effect", {
  set.seed(42)
  control <- rnorm(20, mean = 100, sd = 10)
  treatment <- rnorm(20, mean = 100, sd = 10)

  result <- test_bayes_t(control, treatment)

  # BF should favor null or be inconclusive
  expect_true(result < 10)
})

test_that("test_bayes_t can be used with power_analysis", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    n_per_group = 10,
    find = "power",
    test = "bayes_t",
    n_sim = 50
  )

  expect_s3_class(result, "peppwr_power")
  expect_equal(result$params$test, "bayes_t")
})

# --- test_rankprod ---

test_that("test_rankprod returns p-value", {
  control <- rnorm(10, mean = 100, sd = 10)
  treatment <- rnorm(10, mean = 100, sd = 10)

  result <- test_rankprod(control, treatment)

  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("test_rankprod detects significant difference", {
  set.seed(42)
  control <- rnorm(15, mean = 100, sd = 10)
  treatment <- rnorm(15, mean = 160, sd = 10)  # Large effect

  result <- test_rankprod(control, treatment)

  expect_true(result < 0.05)
})

test_that("test_rankprod handles small samples", {
  # Rank products is designed for small samples
  control <- rnorm(4, mean = 100, sd = 10)
  treatment <- rnorm(4, mean = 130, sd = 10)

  result <- test_rankprod(control, treatment)

  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("test_rankprod can be used with power_analysis", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    n_per_group = 8,
    find = "power",
    test = "rankprod",
    n_sim = 50
  )

  expect_s3_class(result, "peppwr_power")
  expect_equal(result$params$test, "rankprod")
})

# --- Integration: all tests work with run_power_sim ---

test_that("run_power_sim works with bootstrap_t test", {
  result <- run_power_sim(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_per_group = 10,
    effect_size = 2,
    alpha = 0.05,
    test = "bootstrap_t",
    n_sim = 50
  )

  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("run_power_sim works with bayes_t test", {
  result <- run_power_sim(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_per_group = 10,
    effect_size = 2,
    alpha = 0.05,  # Used as BF threshold
    test = "bayes_t",
    n_sim = 50
  )

  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("run_power_sim works with rankprod test", {
  result <- run_power_sim(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_per_group = 8,
    effect_size = 2,
    alpha = 0.05,
    test = "rankprod",
    n_sim = 50
  )

  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})
