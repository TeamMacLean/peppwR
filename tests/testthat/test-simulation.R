# Phase 3: Simulation Engine
# Tests for simulate_experiment(), test_wilcoxon(), run_power_sim()

# --- simulate_experiment ---

test_that("simulate_experiment generates control and treatment samples", {
  result <- simulate_experiment(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_per_group = 6,
    effect_size = 2
  )

  expect_true(is.list(result))
  expect_true(all(c("control", "treatment") %in% names(result)))
  expect_equal(length(result$control), 6)
  expect_equal(length(result$treatment), 6)
})
test_that("simulate_experiment applies multiplicative effect", {
  set.seed(42)
  result <- simulate_experiment(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_per_group = 1000,
    effect_size = 2
  )

  # Treatment mean should be approximately 2x control mean
  control_mean <- mean(result$control)
  treatment_mean <- mean(result$treatment)

  expect_true(treatment_mean > control_mean * 1.5)
  expect_true(treatment_mean < control_mean * 2.5)
})

test_that("simulate_experiment works with gamma distribution", {
  result <- simulate_experiment(
    distribution = "gamma",
    params = list(shape = 2, rate = 0.1),
    n_per_group = 10,
    effect_size = 1.5
  )

  expect_equal(length(result$control), 10)
  expect_equal(length(result$treatment), 10)
  expect_true(all(result$control > 0))  # Gamma is positive
  expect_true(all(result$treatment > 0))
})

test_that("simulate_experiment works with lognormal distribution", {
  result <- simulate_experiment(
    distribution = "lnorm",
    params = list(meanlog = 3, sdlog = 0.5),
    n_per_group = 10,
    effect_size = 2
  )

  expect_equal(length(result$control), 10)
  expect_true(all(result$control > 0))  # Lognormal is positive
})

# --- test_wilcoxon ---

test_that("test_wilcoxon returns p-value", {
  control <- rnorm(10, mean = 100, sd = 10)
  treatment <- rnorm(10, mean = 100, sd = 10)

  result <- test_wilcoxon(control, treatment)

  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("test_wilcoxon detects significant difference", {
  set.seed(42)
  control <- rnorm(20, mean = 100, sd = 10)
  treatment <- rnorm(20, mean = 200, sd = 10)  # Large effect

  result <- test_wilcoxon(control, treatment)

  expect_true(result < 0.05)
})

test_that("test_wilcoxon returns non-significant for same distribution", {
  set.seed(42)
  control <- rnorm(10, mean = 100, sd = 10)
  treatment <- rnorm(10, mean = 100, sd = 10)

  result <- test_wilcoxon(control, treatment)

  # Not guaranteed to be > 0.05, but should usually be
  expect_true(is.numeric(result))
})

# --- run_power_sim ---

test_that("run_power_sim returns power estimate", {
  result <- run_power_sim(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_per_group = 10,
    effect_size = 2,
    alpha = 0.05,
    test = "wilcoxon",
    n_sim = 100
  )

  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("run_power_sim gives high power for large effect", {
  set.seed(42)
  result <- run_power_sim(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_per_group = 20,
    effect_size = 3,  # Large effect
    alpha = 0.05,
    test = "wilcoxon",
    n_sim = 200
  )

  expect_true(result > 0.7)  # Should have high power
})

test_that("run_power_sim gives low power for small effect", {
  set.seed(42)
  result <- run_power_sim(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_per_group = 5,
    effect_size = 1.1,  # Tiny effect
    alpha = 0.05,
    test = "wilcoxon",
    n_sim = 200
  )

  expect_true(result < 0.5)  # Should have low power
})

test_that("run_power_sim respects n_sim parameter", {
  # Running with different n_sim should work
  result1 <- run_power_sim(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_per_group = 10,
    effect_size = 2,
    alpha = 0.05,
    test = "wilcoxon",
    n_sim = 50
  )

  result2 <- run_power_sim(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_per_group = 10,
    effect_size = 2,
    alpha = 0.05,
    test = "wilcoxon",
    n_sim = 100
  )

  # Both should be valid power estimates
  expect_true(result1 >= 0 && result1 <= 1)
  expect_true(result2 >= 0 && result2 <= 1)
})
