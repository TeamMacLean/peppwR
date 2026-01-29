# Phase 4: Power Analysis - Aggregate Mode
# Tests for power_analysis() with aggregate mode (no pilot data)

test_that("power_analysis with find='power' returns peppwr_power", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    n_sim = 100
  )

  expect_s3_class(result, "peppwr_power")
  expect_equal(result$mode, "aggregate")
  expect_equal(result$question, "power")
  expect_true(is.numeric(result$answer))
  expect_true(result$answer >= 0 && result$answer <= 1)
})

test_that("power_analysis with find='sample_size' finds required N", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    target_power = 0.8,
    find = "sample_size",
    n_sim = 100
  )

  expect_s3_class(result, "peppwr_power")
  expect_equal(result$question, "sample_size")
  expect_true(is.numeric(result$answer))
  expect_true(result$answer >= 2)  # At least 2 per group
})

test_that("power_analysis with find='effect_size' finds minimum detectable effect", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_per_group = 10,
    target_power = 0.8,
    find = "effect_size",
    n_sim = 100
  )

  expect_s3_class(result, "peppwr_power")
  expect_equal(result$question, "effect_size")
  expect_true(is.numeric(result$answer))
  expect_true(result$answer > 1)  # Effect size > 1 (fold change)
})

test_that("power_analysis stores simulation parameters", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    n_per_group = 6,
    alpha = 0.05,
    test = "wilcoxon",
    find = "power",
    n_sim = 100
  )

  expect_true("params" %in% names(result))
  expect_equal(result$params$effect_size, 2)
  expect_equal(result$params$n_per_group, 6)
  expect_equal(result$params$alpha, 0.05)
  expect_equal(result$params$test, "wilcoxon")
})

test_that("power_analysis uses default alpha of 0.05", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    n_sim = 100
  )

  expect_equal(result$params$alpha, 0.05)
})

test_that("power_analysis uses default test of wilcoxon", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    n_sim = 100
  )

  expect_equal(result$params$test, "wilcoxon")
})

test_that("power_analysis errors with invalid find parameter", {
  expect_error(
    power_analysis(
      distribution = "norm",
      params = list(mean = 100, sd = 10),
      effect_size = 2,
      n_per_group = 6,
      find = "invalid"
    ),
    "find"
  )
})

test_that("power_analysis errors when missing required params for find='power'", {
  # Missing n_per_group
  expect_error(
    power_analysis(
      distribution = "norm",
      params = list(mean = 100, sd = 10),
      effect_size = 2,
      find = "power"
    )
  )
})

test_that("power_analysis errors when missing required params for find='sample_size'", {
  # Missing target_power
  expect_error(
    power_analysis(
      distribution = "norm",
      params = list(mean = 100, sd = 10),
      effect_size = 2,
      find = "sample_size"
    )
  )
})

test_that("power_analysis works with gamma distribution", {
  result <- power_analysis(
    distribution = "gamma",
    params = list(shape = 2, rate = 0.1),
    effect_size = 2,
    n_per_group = 8,
    find = "power",
    n_sim = 100
  )

  expect_s3_class(result, "peppwr_power")
  expect_true(result$answer >= 0 && result$answer <= 1)
})

test_that("power_analysis works with lognormal distribution", {
  result <- power_analysis(
    distribution = "lnorm",
    params = list(meanlog = 3, sdlog = 0.5),
    effect_size = 1.5,
    n_per_group = 10,
    find = "power",
    n_sim = 100
  )

  expect_s3_class(result, "peppwr_power")
})

test_that("power increases with sample size", {
  power_n5 <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 1.5,
    n_per_group = 5,
    find = "power",
    n_sim = 200
  )

  power_n20 <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 1.5,
    n_per_group = 20,
    find = "power",
    n_sim = 200
  )

  expect_true(power_n20$answer > power_n5$answer)
})

test_that("power increases with effect size", {
  power_small <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 1.2,
    n_per_group = 10,
    find = "power",
    n_sim = 200
  )

  power_large <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 3,
    n_per_group = 10,
    find = "power",
    n_sim = 200
  )

  expect_true(power_large$answer > power_small$answer)
})
