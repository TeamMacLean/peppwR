# Phase 6: Core Plots
# Tests for plot.peppwr_power (aggregate and per-peptide modes)

test_that("plot.peppwr_power aggregate mode returns ggplot", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    n_sim = 50
  )

  p <- plot(result)

  expect_s3_class(p, "ggplot")
})

test_that("plot.peppwr_power aggregate shows power curve", {
  # Need to compute power across range of N for this plot
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    target_power = 0.8,
    find = "sample_size",
    n_sim = 50
  )

  p <- plot(result)

  expect_s3_class(p, "ggplot")
  # Should have N on x-axis and power on y-axis
  expect_true("x" %in% names(p$mapping) || "n" %in% names(p$data) ||
                "n_per_group" %in% names(p$data))
})

test_that("plot.peppwr_power aggregate includes 80% power line", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    target_power = 0.8,
    find = "sample_size",
    n_sim = 50
  )

  p <- plot(result)

  # Check that there's a horizontal line layer (geom_hline)
  layer_types <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomHline" %in% layer_types || "GeomAbline" %in% layer_types ||
                any(grepl("line", layer_types, ignore.case = TRUE)))
})
test_that("plot.peppwr_power aggregate annotates recommended N", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    target_power = 0.8,
    find = "sample_size",
    n_sim = 50
  )

  p <- plot(result)

  # Should have some text or annotation indicating the recommended N
  # This could be via geom_text, geom_label, geom_vline, or labs
  has_annotation <- any(sapply(p$layers, function(l) {
    grepl("text|label|vline|point", class(l$geom)[1], ignore.case = TRUE)
  })) || !is.null(p$labels$title) || !is.null(p$labels$subtitle)

  expect_true(has_annotation)
})

test_that("plot.peppwr_power per-peptide mode returns ggplot", {
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
    target_power = 0.8,
    find = "sample_size",
    n_sim = 50
  )

  p <- plot(result)

  expect_s3_class(p, "ggplot")
})

test_that("plot.peppwr_power per-peptide shows % peptides at threshold", {
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
    target_power = 0.8,
    find = "sample_size",
    n_sim = 50
  )

  p <- plot(result)

  # Y-axis should represent percentage or proportion
  # Check axis label or data range
  expect_true(
    grepl("percent|proportion|%", p$labels$y, ignore.case = TRUE) ||
      (max(p$data[[2]], na.rm = TRUE) <= 100)  # Assuming y is second column
  )
})

test_that("plot.peppwr_power has informative title", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    n_sim = 50
  )

  p <- plot(result)

  # Should have a title

  expect_true(!is.null(p$labels$title) || !is.null(p$labels$subtitle))
})
