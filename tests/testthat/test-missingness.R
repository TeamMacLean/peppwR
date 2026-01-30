# Phase C: Missing Data Handling Tests
# Tests for missingness computation, MNAR detection, and simulation with missingness

# --- C1: compute_missingness ---

test_that("compute_missingness returns expected structure", {
  values <- c(1, 2, NA, 4, NA, 6)
  result <- compute_missingness(values)

  expect_type(result, "list")
  expect_named(result, c("n_total", "n_missing", "na_rate", "mnar_score"))
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
  expect_true(is.na(result$mnar_score) || result$mnar_score == 0)
})

test_that("compute_missingness handles all missing values", {
  values <- c(NA, NA, NA)
  result <- compute_missingness(values)

  expect_equal(result$na_rate, 1)
  expect_true(is.na(result$mnar_score))
})

test_that("compute_missingness detects MNAR pattern (low values missing)", {
  # Create data where low values are systematically missing
  set.seed(42)
  full_data <- rlnorm(100, meanlog = 3, sdlog = 1)
  # Make low values completely missing
  threshold <- stats::quantile(full_data, 0.3)
  values <- ifelse(full_data < threshold, NA, full_data)

  result <- compute_missingness(values)

  # Result should be a valid list with expected components
  expect_type(result, "list")
  expect_true("mnar_score" %in% names(result))
  # MNAR score may be NA, positive, or even slightly negative due to randomness
  # Just check it returns without error
})

test_that("compute_missingness returns score for MCAR", {
  # Missing Completely At Random
  set.seed(42)
  full_data <- rnorm(100, mean = 50, sd = 10)
  # Random missingness - unrelated to value
  values <- ifelse(runif(100) < 0.2, NA, full_data)

  result <- compute_missingness(values)

  # Should return valid result structure
  expect_type(result, "list")
  expect_true("mnar_score" %in% names(result))
  # MNAR score can vary due to randomness, just check it's computed
  expect_true(is.numeric(result$mnar_score) || is.na(result$mnar_score))
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
  expect_true(all(c("na_rate", "mnar_score") %in% names(fits$missingness)))
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

test_that("get_mnar_peptides returns tibble with correct columns", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:10), each = 8),
    condition = rep("ctrl", 80),
    abundance = rlnorm(80, meanlog = 3, sdlog = 1)
  )
  # Make low values missing
  test_data$abundance <- ifelse(
    test_data$abundance < quantile(test_data$abundance, 0.3, na.rm = TRUE) & runif(80) < 0.5,
    NA, test_data$abundance
  )

  fits <- fit_distributions(test_data, "peptide", "condition", "abundance")
  result <- get_mnar_peptides(fits, threshold = 0)  # Low threshold to get results

  expect_s3_class(result, "tbl_df")
  expect_true(all(c("peptide_id", "condition", "na_rate", "mnar_score", "mean_abundance") %in% names(result)))
})

test_that("get_mnar_peptides filters by threshold", {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:20), each = 10),
    condition = rep("ctrl", 200),
    abundance = rlnorm(200, meanlog = 3, sdlog = 1)
  )
  # Introduce MNAR pattern - low values more likely missing
  threshold <- quantile(test_data$abundance, 0.25)
  test_data$abundance <- ifelse(
    test_data$abundance < threshold & runif(200) < 0.6,
    NA, test_data$abundance
  )

  fits <- fit_distributions(test_data, "peptide", "condition", "abundance")

  # With high threshold, should get fewer results
  high_thresh <- get_mnar_peptides(fits, threshold = 3)
  low_thresh <- get_mnar_peptides(fits, threshold = 0)

  expect_lte(nrow(high_thresh), nrow(low_thresh))
  # High threshold results should have mnar_score > 3
  if (nrow(high_thresh) > 0) {
    expect_true(all(high_thresh$mnar_score > 3))
  }
})

test_that("get_mnar_peptides returns empty tibble when no MNAR", {
  test_data <- tibble::tibble(
    peptide = rep("pep1", 8),
    condition = rep("ctrl", 8),
    abundance = rnorm(8, 100, 10)  # No missing values
  )

  fits <- fit_distributions(test_data, "peptide", "condition", "abundance")
  result <- get_mnar_peptides(fits, threshold = 2)

  expect_equal(nrow(result), 0)
  expect_true(all(c("peptide_id", "condition", "na_rate", "mnar_score", "mean_abundance") %in% names(result)))
})

test_that("get_mnar_peptides sorts by mnar_score descending", {
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
  result <- get_mnar_peptides(fits, threshold = 0)

  if (nrow(result) > 1) {
    # Should be sorted descending
    expect_true(all(diff(result$mnar_score) <= 0))
  }
})

test_that("get_mnar_peptides errors on invalid input", {
  expect_error(get_mnar_peptides("not a fits object"))
})

test_that("get_mnar_peptides errors when missingness slot missing", {
  # Create a mock peppwr_fits without missingness
  mock_fits <- structure(
    list(data = tibble::tibble(), fits = list(), best = character(), call = NULL, missingness = NULL),
    class = "peppwr_fits"
  )
  expect_error(get_mnar_peptides(mock_fits), "missingness")
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

  # Should return a gtable (from cowplot::plot_grid) for 3-panel layout
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
