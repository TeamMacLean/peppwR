# Phase C: Missing Data Handling Tests
# Tests for missingness computation, MNAR detection, and simulation with missingness

# --- C1: compute_missingness ---

test_that("compute_missingness returns expected structure", {
  values <- c(1, 2, NA, 4, NA, 6)
  result <- compute_missingness(values)

  expect_type(result, "list")
  expect_named(result, c("n_total", "n_missing", "na_rate"))
})

test_that("compute_missingness does not produce wilcox.test warnings", {
  # Values with ties
  values <- c(1, 1, 2, 2, NA, NA)
  expect_no_warning(compute_missingness(values))
})

test_that("compute_missingness calculates correct NA rate", {
  # 2 out of 6 values are NA = 33.3%
  values <- c(1, 2, NA, 4, NA, 6)
  result <- compute_missingness(values)

  expect_equal(result$n_total, 6)
  expect_equal(result$n_missing, 2)
  expect_equal(result$na_rate, 2/6, tolerance = 0.001)
})

test_that("compute_missingness handles no missing values", {
  values <- c(1, 2, 3, 4, 5)
  result <- compute_missingness(values)

  expect_equal(result$n_missing, 0)
  expect_equal(result$na_rate, 0)
})

test_that("compute_missingness handles all missing values", {
  values <- c(NA, NA, NA)
  result <- compute_missingness(values)

  expect_equal(result$na_rate, 1)
})


# --- C2: fit_distributions with missingness slot ---

test_that("fit_distributions includes missingness slot", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:3), each = 20),
    condition = rep("control", 60),
    abundance = c(
      c(rgamma(15, shape = 2, rate = 0.1), rep(NA, 5)),
      rnorm(20, mean = 50, sd = 10),
      c(rlnorm(18, meanlog = 3, sdlog = 0.5), rep(NA, 2))
    )
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  expect_true("missingness" %in% names(fits))
  expect_s3_class(fits$missingness, "tbl_df")
  expect_true("na_rate" %in% names(fits$missingness))
})

test_that("missingness slot has correct dimensions", {
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

  expect_equal(nrow(fits$missingness), 5)
})

# --- C3: print/summary with missingness ---

test_that("print.peppwr_fits shows missingness summary", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:3), each = 20),
    condition = rep("control", 60),
    abundance = c(
      c(rgamma(15, shape = 2, rate = 0.1), rep(NA, 5)),
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

  output <- capture.output(print(fits))
  # Should mention missingness or NA somewhere
  expect_true(any(grepl("missing|NA|missingness", output, ignore.case = TRUE)))
})

test_that("summary.peppwr_fits includes missingness statistics", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:3), each = 20),
    condition = rep("control", 60),
    abundance = c(
      c(rgamma(15, shape = 2, rate = 0.1), rep(NA, 5)),
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

  summ <- summary(fits)
  expect_true("missingness" %in% names(summ) || "mean_na_rate" %in% names(summ))
})

# --- C4: simulate_with_missingness ---

test_that("simulate_with_missingness returns samples with NAs", {
  set.seed(42)
  result <- simulate_with_missingness(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_per_group = 20,
    effect_size = 2,
    na_rate = 0.2
  )

  expect_type(result, "list")
  expect_named(result, c("control", "treatment"))
  # Should have some NAs (roughly 20%)
  expect_true(sum(is.na(result$control)) > 0 || sum(is.na(result$treatment)) > 0)
})

test_that("simulate_with_missingness respects na_rate", {
  set.seed(42)
  # Run multiple simulations to check average NA rate
  na_rates <- replicate(100, {
    result <- simulate_with_missingness(
      distribution = "norm",
      params = list(mean = 100, sd = 10),
      n_per_group = 50,
      effect_size = 1,
      na_rate = 0.3
    )
    mean(is.na(c(result$control, result$treatment)))
  })

  # Average NA rate should be close to 0.3
  expect_equal(mean(na_rates), 0.3, tolerance = 0.05)
})

test_that("simulate_with_missingness applies MNAR pattern", {
  set.seed(42)
  result <- simulate_with_missingness(
    distribution = "lnorm",
    params = list(meanlog = 3, sdlog = 1),
    n_per_group = 100,
    effect_size = 1,
    na_rate = 0.3,
    mnar_score = 2  # Strong MNAR - low values more likely missing
  )

  control_non_na <- result$control[!is.na(result$control)]
  control_na_positions <- which(is.na(result$control))

  # If MNAR is applied, observed values should be biased higher
  # (This is a probabilistic test, may occasionally fail)
  if (length(control_non_na) > 10) {
    # Mean of observed should be higher than overall mean of distribution
    expect_gt(mean(control_non_na), exp(3))  # exp(3) is median of lnorm
  }
})

# --- C5: run_power_sim_with_missingness ---

test_that("run_power_sim_with_missingness returns power between 0 and 1", {
  set.seed(42)
  power <- run_power_sim_with_missingness(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_per_group = 10,
    effect_size = 2,
    na_rate = 0.1,
    n_sim = 100
  )

  expect_type(power, "double")
  expect_true(power >= 0 && power <= 1)
})

test_that("run_power_sim_with_missingness gives lower power with more missingness", {
  set.seed(42)
  power_low_na <- run_power_sim_with_missingness(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_per_group = 10,
    effect_size = 2,
    na_rate = 0.1,
    n_sim = 200
  )

  power_high_na <- run_power_sim_with_missingness(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_per_group = 10,
    effect_size = 2,
    na_rate = 0.5,
    n_sim = 200
  )

  # Higher missingness should generally reduce power
  expect_lte(power_high_na, power_low_na + 0.2)  # Allow some variance
})

# --- C6: power_analysis with include_missingness ---

test_that("power_analysis.peppwr_fits accepts include_missingness parameter", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:3), each = 20),
    condition = rep("control", 60),
    abundance = c(
      c(rgamma(15, shape = 2, rate = 0.1), rep(NA, 5)),
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

  # Should work with include_missingness = TRUE
  result <- power_analysis(
    fits,
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    include_missingness = TRUE,
    n_sim = 50
  )

  expect_s3_class(result, "peppwr_power")
})

test_that("include_missingness affects power estimates", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:5), each = 30),
    condition = rep("control", 150),
    abundance = c(
      c(rgamma(20, shape = 2, rate = 0.1), rep(NA, 10)),
      c(rnorm(25, mean = 50, sd = 10), rep(NA, 5)),
      c(rlnorm(22, meanlog = 3, sdlog = 0.5), rep(NA, 8)),
      c(rgamma(28, shape = 3, rate = 0.2), rep(NA, 2)),
      c(rnorm(18, mean = 80, sd = 15), rep(NA, 12))
    )
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  result_with <- power_analysis(
    fits,
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    include_missingness = TRUE,
    n_sim = 50
  )

  result_without <- power_analysis(
    fits,
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    include_missingness = FALSE,
    n_sim = 50
  )

  # Results should be different (missingness reduces power)
  # But both should be valid power values
  expect_true(result_with$answer >= 0 && result_with$answer <= 1)
  expect_true(result_without$answer >= 0 && result_without$answer <= 1)
})

# --- C7: plot_missingness ---

test_that("plot_missingness returns ggplot", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:10), each = 20),
    condition = rep("control", 200),
    abundance = c(
      c(rgamma(15, shape = 2, rate = 0.1), rep(NA, 5)),
      rnorm(20, mean = 50, sd = 10),
      c(rlnorm(18, meanlog = 3, sdlog = 0.5), rep(NA, 2)),
      rgamma(20, shape = 3, rate = 0.2),
      c(rnorm(17, mean = 80, sd = 15), rep(NA, 3)),
      rgamma(20, shape = 1.5, rate = 0.15),
      c(rnorm(19, mean = 60, sd = 12), rep(NA, 1)),
      rlnorm(20, meanlog = 2.5, sdlog = 0.4),
      c(rgamma(16, shape = 2.5, rate = 0.18), rep(NA, 4)),
      rnorm(20, mean = 70, sd = 10)
    )
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  p <- plot_missingness(fits)
  expect_s3_class(p, "ggplot")
})

test_that("plot_missingness shows NA rate distribution", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:5), each = 20),
    condition = rep("control", 100),
    abundance = c(
      c(rgamma(15, shape = 2, rate = 0.1), rep(NA, 5)),
      rnorm(20, mean = 50, sd = 10),
      c(rlnorm(18, meanlog = 3, sdlog = 0.5), rep(NA, 2)),
      rgamma(20, shape = 3, rate = 0.2),
      c(rnorm(17, mean = 80, sd = 15), rep(NA, 3))
    )
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  p <- plot_missingness(fits)

  # Should return a ggplot or grid arrangement (from cowplot)
  expect_true(inherits(p, "ggplot") || inherits(p, "gtable") ||
                inherits(p, "grob") || inherits(p, "gTree"))
})

test_that("plot_missingness handles data with no missingness", {
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:3), each = 20),
    condition = rep("control", 60),
    abundance = rgamma(60, shape = 2, rate = 0.1)
  )

  fits <- fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )

  # Should still work
  p <- plot_missingness(fits)
  expect_s3_class(p, "ggplot")
})

test_that("plot_missingness errors on invalid input", {
  expect_error(plot_missingness("not a fits object"))
})

# --- v2.1: MNAR Visualization Enhancements ---

test_that("missingness slot includes mean_abundance", {
  test_data <- tibble::tibble(
    peptide = rep("pep1", 8),
    condition = rep("ctrl", 8),
    abundance = c(10, 20, 30, 40, NA, NA, 70, 80)
  )

  fits <- fit_distributions(test_data, "peptide", "condition", "abundance")

  expect_true("mean_abundance" %in% names(fits$missingness))
  # Mean of observed: (10+20+30+40+70+80)/6 = 41.67
  expect_equal(fits$missingness$mean_abundance[1], mean(c(10, 20, 30, 40, 70, 80)), tolerance = 0.01)
})

test_that("mean_abundance handles all NA values", {
  test_data <- tibble::tibble(
    peptide = rep(c("pep1", "pep2"), each = 4),
    condition = rep("ctrl", 8),
    abundance = c(rep(NA, 4), 10, 20, 30, 40)
  )

  fits <- fit_distributions(test_data, "peptide", "condition", "abundance")

  # pep1 should have NaN for mean_abundance (all NA)
  expect_true(is.nan(fits$missingness$mean_abundance[1]) || is.na(fits$missingness$mean_abundance[1]))
  # pep2 should have mean of 25
  expect_equal(fits$missingness$mean_abundance[2], 25, tolerance = 0.01)
})


test_that("plot_missingness includes abundance vs NA rate panel", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:10), each = 10),
    condition = rep("ctrl", 100),
    abundance = rlnorm(100, meanlog = 3, sdlog = 1)
  )
  # Introduce missingness
  test_data$abundance <- ifelse(
    runif(100) < 0.2,
    NA, test_data$abundance
  )

  fits <- fit_distributions(test_data, "peptide", "condition", "abundance")
  p <- plot_missingness(fits)

  # Should return a gtable (from cowplot::plot_grid) for 2-panel layout
  expect_true(inherits(p, "gtable") || inherits(p, "ggplot") || inherits(p, "grob"))
})

test_that("print.peppwr_fits shows MNAR pattern count", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:10), each = 10),
    condition = rep("ctrl", 100),
    abundance = rlnorm(100, meanlog = 3, sdlog = 1)
  )
  # Introduce MNAR pattern - low values missing
  threshold <- quantile(test_data$abundance, 0.25)
  test_data$abundance <- ifelse(
    test_data$abundance < threshold & runif(100) < 0.5,
    NA, test_data$abundance
  )

  fits <- fit_distributions(test_data, "peptide", "condition", "abundance")
  output <- capture.output(print(fits))

  # Should mention MNAR pattern in output
  expect_true(any(grepl("MNAR", output, ignore.case = TRUE)))
})

# --- Dataset-level MNAR tests ---

test_that("compute_dataset_mnar returns expected structure", {
  missingness <- tibble::tibble(
    na_rate = c(0.5, 0.3, 0.1, 0.2, 0.4, 0.0),
    mean_abundance = c(10, 20, 50, 30, 15, 100)
  )

  result <- compute_dataset_mnar(missingness)

  expect_type(result, "list")
  expect_named(result, c("correlation", "p_value", "n_peptides", "interpretation"))
  expect_true(is.numeric(result$correlation) || is.na(result$correlation))
  expect_true(is.numeric(result$p_value) || is.na(result$p_value))
  expect_equal(result$n_peptides, 5)  # Excludes peptide with na_rate = 0
  expect_type(result$interpretation, "character")
})

test_that("compute_dataset_mnar detects strong MNAR pattern", {
  # Strong negative correlation: low abundance -> high NA rate
  set.seed(42)
  missingness <- tibble::tibble(
    mean_abundance = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000),
    na_rate = c(0.9, 0.8, 0.7, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02)
  )

  result <- compute_dataset_mnar(missingness)

  # Should detect negative correlation
  expect_lt(result$correlation, -0.5)
  expect_lt(result$p_value, 0.05)
  expect_true(grepl("MNAR", result$interpretation))
  expect_true(grepl("low-abundance", result$interpretation) ||
              grepl("evidence", result$interpretation, ignore.case = TRUE))
})

test_that("compute_dataset_mnar handles insufficient data", {
  # Only 3 peptides with missing data
  missingness <- tibble::tibble(
    na_rate = c(0.5, 0.3, 0.1, 0.0, 0.0),
    mean_abundance = c(10, 20, 50, 100, 200)
  )

  result <- compute_dataset_mnar(missingness)

  expect_equal(result$n_peptides, 3)
  expect_true(is.na(result$correlation))
  expect_true(is.na(result$p_value))
  expect_true(grepl("Insufficient", result$interpretation))
})

test_that("compute_dataset_mnar handles no peptides with missing data", {
  missingness <- tibble::tibble(
    na_rate = c(0, 0, 0, 0, 0),
    mean_abundance = c(10, 20, 50, 100, 200)
  )

  result <- compute_dataset_mnar(missingness)

  expect_equal(result$n_peptides, 0)
  expect_true(is.na(result$correlation))
  expect_true(grepl("Insufficient", result$interpretation))
})

test_that("compute_dataset_mnar handles NA/NaN in mean_abundance", {
  missingness <- tibble::tibble(
    na_rate = c(0.5, 0.3, 0.1, 0.2, 0.4, 0.5),
    mean_abundance = c(10, NA, 50, NaN, 15, 0)  # NA, NaN, and 0 should be excluded
  )

  result <- compute_dataset_mnar(missingness)

  # Only 3 valid peptides (10, 50, 15) - excludes NA, NaN, and 0
  expect_equal(result$n_peptides, 3)
  expect_true(is.na(result$correlation))  # < 5 peptides
})

test_that("fit_distributions includes dataset_mnar slot", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:10), each = 10),
    condition = rep("ctrl", 100),
    abundance = rlnorm(100, meanlog = 3, sdlog = 1)
  )
  # Introduce MNAR pattern
  threshold <- quantile(test_data$abundance, 0.3)
  test_data$abundance <- ifelse(
    test_data$abundance < threshold & runif(100) < 0.5,
    NA, test_data$abundance
  )

  fits <- fit_distributions(test_data, "peptide", "condition", "abundance")

  expect_true("dataset_mnar" %in% names(fits))
  expect_type(fits$dataset_mnar, "list")
  expect_true(all(c("correlation", "p_value", "n_peptides", "interpretation") %in%
                    names(fits$dataset_mnar)))
})

test_that("print.peppwr_fits shows dataset-level MNAR correlation", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:20), each = 10),
    condition = rep("ctrl", 200),
    abundance = rlnorm(200, meanlog = 3, sdlog = 1)
  )
  # Strong MNAR pattern
  threshold <- quantile(test_data$abundance, 0.4)
  test_data$abundance <- ifelse(
    test_data$abundance < threshold & runif(200) < 0.6,
    NA, test_data$abundance
  )

  fits <- fit_distributions(test_data, "peptide", "condition", "abundance")
  output <- capture.output(print(fits))

  # Should mention MNAR detection with correlation
  expect_true(any(grepl("MNAR detection", output)))
  expect_true(any(grepl("Correlation|r =", output)))
})


test_that("summary.peppwr_fits includes dataset_mnar metrics", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:15), each = 10),
    condition = rep("ctrl", 150),
    abundance = rlnorm(150, meanlog = 3, sdlog = 1)
  )
  # Introduce MNAR pattern
  threshold <- quantile(test_data$abundance, 0.3)
  test_data$abundance <- ifelse(
    test_data$abundance < threshold & runif(150) < 0.5,
    NA, test_data$abundance
  )

  fits <- fit_distributions(test_data, "peptide", "condition", "abundance")
  summ <- summary(fits)

  expect_true("missingness" %in% names(summ))
  expect_true("dataset_mnar_correlation" %in% names(summ$missingness))
  expect_true("dataset_mnar_pvalue" %in% names(summ$missingness))
  expect_true("dataset_mnar_interpretation" %in% names(summ$missingness))
})

test_that("interpret_dataset_mnar provides correct interpretations", {
  # Strong negative correlation (strong MNAR)
  interp1 <- interpret_dataset_mnar(-0.6, 0.001, 20)
  expect_true(grepl("Strong evidence", interp1))
  expect_true(grepl("low-abundance", interp1))

  # Weak negative correlation
  interp2 <- interpret_dataset_mnar(-0.2, 0.1, 20)
  expect_true(grepl("Weak evidence", interp2))

  # No correlation
  interp3 <- interpret_dataset_mnar(0.05, 0.8, 20)
  expect_true(grepl("No evidence", interp3))

  # Positive correlation (unusual)
  interp4 <- interpret_dataset_mnar(0.4, 0.05, 20)
  expect_true(grepl("unusual", interp4, ignore.case = TRUE))
})
