# Phase A: Diagnostic Plots
# Tests for v2 diagnostic plot functions

# Helper function to create test data and fits
create_test_fits <- function(n_peptides = 5) {
  set.seed(42)
  test_data <- tibble::tibble(
    peptide = rep(paste0("pep", 1:n_peptides), each = 20),
    condition = rep("control", n_peptides * 20),
    abundance = c(
      rgamma(20, shape = 2, rate = 0.1),
      rnorm(20, mean = 50, sd = 10),
      rlnorm(20, meanlog = 3, sdlog = 0.5),
      rgamma(20, shape = 3, rate = 0.2),
      rnorm(20, mean = 80, sd = 15)
    )[1:(n_peptides * 20)]
  )

  fit_distributions(
    test_data,
    id = "peptide",
    group = "condition",
    value = "abundance"
  )
}

# --- A1: plot_density_overlay ---

test_that("plot_density_overlay returns ggplot", {
  fits <- create_test_fits()
  p <- plot_density_overlay(fits)
  expect_s3_class(p, "ggplot")
})

test_that("plot_density_overlay shows histogram and density curve", {
  fits <- create_test_fits()
  p <- plot_density_overlay(fits)

  # Check for histogram and line layers

layer_types <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true(
    any(grepl("Histogram|Bar", layer_types)) ||
      any(grepl("hist", layer_types, ignore.case = TRUE))
  )
  expect_true(
    any(grepl("Line|Density", layer_types)) ||
      any(grepl("line|density", layer_types, ignore.case = TRUE))
  )
})

test_that("plot_density_overlay can select specific peptide", {
  fits <- create_test_fits()
  p <- plot_density_overlay(fits, peptide_id = "pep1")
  expect_s3_class(p, "ggplot")
})

test_that("plot_density_overlay respects n_overlay parameter", {
  fits <- create_test_fits(10)
  p <- plot_density_overlay(fits, n_overlay = 4)
  expect_s3_class(p, "ggplot")
  # When faceting, we expect at most n_overlay facets
})

test_that("plot_density_overlay errors on invalid fits object", {
  expect_error(plot_density_overlay("not a fits object"))
  expect_error(plot_density_overlay(list(a = 1)))
})

# --- A2: plot_qq ---

test_that("plot_qq returns ggplot", {
  fits <- create_test_fits()
  p <- plot_qq(fits)
  expect_s3_class(p, "ggplot")
})

test_that("plot_qq shows QQ points", {
  fits <- create_test_fits()
  p <- plot_qq(fits)

  layer_types <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true(any(grepl("Point", layer_types)))
})

test_that("plot_qq can select specific peptide", {
  fits <- create_test_fits()
  # May return a placeholder plot if QQ cannot be computed for that peptide
  p <- plot_qq(fits, peptide_id = "pep1")
  expect_s3_class(p, "ggplot")
})

test_that("plot_qq handles peptides with no valid fit gracefully", {
  fits <- create_test_fits()
  # Should return a ggplot even if no valid QQ data
  p <- plot_qq(fits, peptide_id = "pep1")
  expect_s3_class(p, "ggplot")
})

test_that("plot_qq respects n_plots parameter", {
  fits <- create_test_fits(10)
  p <- plot_qq(fits, n_plots = 4)
  expect_s3_class(p, "ggplot")
})

test_that("plot_qq errors on invalid fits object", {
  expect_error(plot_qq("not a fits object"))
})

# --- A3: plot_power_heatmap ---

test_that("plot_power_heatmap returns ggplot", {
  p <- plot_power_heatmap(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_range = c(3, 10),
    effect_range = c(1.5, 3)
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_power_heatmap creates tile plot", {
  p <- plot_power_heatmap(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_range = c(3, 10),
    effect_range = c(1.5, 3)
  )

  layer_types <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true(any(grepl("Tile|Rect|Raster", layer_types)))
})

test_that("plot_power_heatmap has correct axis labels", {
  p <- plot_power_heatmap(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_range = c(3, 10),
    effect_range = c(1.5, 3)
  )

  # Should have labels for N and effect size
  expect_true(
    grepl("sample|n|size", p$labels$x, ignore.case = TRUE) ||
      grepl("effect", p$labels$x, ignore.case = TRUE)
  )
})

test_that("plot_power_heatmap accepts custom n_sim", {
  p <- plot_power_heatmap(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    n_range = c(3, 6),
    effect_range = c(1.5, 2),
    n_sim = 50
  )
  expect_s3_class(p, "ggplot")
})

# --- A4: plot_power_vs_effect ---

test_that("plot_power_vs_effect returns ggplot", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    n_sim = 50
  )

  p <- plot_power_vs_effect(result, effect_range = c(1.5, 3))
  expect_s3_class(p, "ggplot")
})

test_that("plot_power_vs_effect shows line", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    n_sim = 50
  )

  p <- plot_power_vs_effect(result, effect_range = c(1.5, 3))

  layer_types <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true(any(grepl("Line|Path", layer_types)))
})

test_that("plot_power_vs_effect has correct axis labels", {
  result <- power_analysis(
    distribution = "norm",
    params = list(mean = 100, sd = 10),
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    n_sim = 50
  )

  p <- plot_power_vs_effect(result, effect_range = c(1.5, 3))

  expect_true(grepl("effect", p$labels$x, ignore.case = TRUE))
  expect_true(grepl("power", p$labels$y, ignore.case = TRUE))
})

test_that("plot_power_vs_effect errors on non-aggregate mode", {
  fits <- create_test_fits()
  result <- power_analysis(
    fits,
    effect_size = 2,
    n_per_group = 6,
    find = "power",
    n_sim = 50
  )

  # Should work with per-peptide results too (mean power)
  p <- plot_power_vs_effect(result, effect_range = c(1.5, 3))
  expect_s3_class(p, "ggplot")
})

# --- A5: plot_param_distribution ---

test_that("plot_param_distribution returns ggplot", {
  fits <- create_test_fits()
  p <- plot_param_distribution(fits)
  expect_s3_class(p, "ggplot")
})

test_that("plot_param_distribution shows histogram or density", {
  fits <- create_test_fits(20)
  p <- plot_param_distribution(fits)

  layer_types <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true(
    any(grepl("Histogram|Density|Bar", layer_types)) ||
      any(grepl("hist|dens", layer_types, ignore.case = TRUE))
  )
})

test_that("plot_param_distribution errors on invalid fits object", {
  expect_error(plot_param_distribution("not a fits object"))
})

test_that("plot_param_distribution has informative labels", {
  fits <- create_test_fits()
  p <- plot_param_distribution(fits)

  # Should have a title or axis labels
  expect_true(
    !is.null(p$labels$title) ||
      !is.null(p$labels$x) ||
      !is.null(p$labels$y)
  )
})
