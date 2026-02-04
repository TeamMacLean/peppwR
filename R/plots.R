# Plot methods for peppwR classes

#' Plot method for peppwr_fits
#'
#' Creates a bar chart showing the count of best-fit distributions
#'
#' @param x A peppwr_fits object
#' @param ... Additional arguments (ignored)
#'
#' @return A ggplot object
#' @export
#' @method plot peppwr_fits
plot.peppwr_fits <- function(x, ...) {
  validate_peppwr_fits(x)

  # Create data for plot: count of best distributions by AIC and LogLik
  fits_df <- x$data

  ll <- fits_df |>
    dplyr::mutate(best = purrr::map(.data$fits, \(f) f$dist[which.max(f$loglik)])) |>
    tidyr::unnest(cols = best) |>
    dplyr::count(best) |>
    dplyr::mutate(metric = "LogLikelihood")

  la <- fits_df |>
    dplyr::mutate(best = purrr::map(.data$fits, \(f) f$dist[which.min(f$aic)])) |>
    tidyr::unnest(cols = best) |>
    dplyr::count(best) |>
    dplyr::mutate(metric = "AIC")

  plot_data <- dplyr::bind_rows(ll, la)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = best, y = n, fill = best)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_wrap(~ metric, ncol = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Distribution",
      y = "Count",
      fill = "Distribution",
      title = "Best Fitted Distributions"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot method for peppwr_power
#'
#' Creates power curves or % peptides at threshold plots
#'
#' @param x A peppwr_power object
#' @param ... Additional arguments (ignored)
#'
#' @return A ggplot object
#' @export
#' @method plot peppwr_power
plot.peppwr_power <- function(x, ...) {
  validate_peppwr_power(x)

  if (x$mode == "aggregate") {
    plot_power_aggregate(x)
  } else {
    plot_power_perpeptide(x)
  }
}

#' Plot power curve for aggregate mode
#'
#' @param x A peppwr_power object
#' @return A ggplot object
#' @keywords internal
plot_power_aggregate <- function(x) {
  # Get power curve data
  if (x$question == "sample_size" && "power_curve" %in% names(x$simulations)) {
    plot_data <- x$simulations$power_curve
    target_power <- x$params$target_power
    recommended_n <- x$answer

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = n_per_group, y = power)) +
      ggplot2::geom_line(linewidth = 1, color = "steelblue") +
      ggplot2::geom_point(size = 2, color = "steelblue") +
      ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
      ggplot2::geom_vline(xintercept = recommended_n, linetype = "dotted", color = "darkgreen") +
      ggplot2::geom_point(
        data = data.frame(n_per_group = recommended_n, power = target_power),
        ggplot2::aes(x = n_per_group, y = power),
        size = 4, color = "darkgreen"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Sample Size (per group)",
        y = "Power",
        title = paste0("Recommended sample size: N=", recommended_n, " per group"),
        subtitle = paste0("Effect size: ", x$params$effect_size, "-fold, Alpha: ", x$params$alpha)
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      ggplot2::annotate(
        "text",
        x = recommended_n,
        y = target_power + 0.05,
        label = paste0("N=", recommended_n),
        color = "darkgreen"
      )
  } else if (x$question == "effect_size" && "effect_curve" %in% names(x$simulations)) {
    plot_data <- x$simulations$effect_curve
    target_power <- x$params$target_power
    min_effect <- x$answer

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = effect_size, y = power)) +
      ggplot2::geom_line(linewidth = 1, color = "steelblue") +
      ggplot2::geom_point(size = 2, color = "steelblue") +
      ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "red") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Effect Size (fold change)",
        y = "Power",
        title = paste0("Minimum detectable effect: ", min_effect, "-fold"),
        subtitle = paste0("N=", x$params$n_per_group, " per group, Alpha: ", x$params$alpha)
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format())
  } else {
    # Simple single-point result
    plot_data <- data.frame(
      n = x$params$n_per_group,
      power = x$answer
    )

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = n, y = power)) +
      ggplot2::geom_point(size = 4, color = "steelblue") +
      ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Sample Size (per group)",
        y = "Power",
        title = paste0("Power: ", round(x$answer * 100), "%"),
        subtitle = paste0("N=", x$params$n_per_group, ", Effect: ", x$params$effect_size, "-fold")
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format())
  }

  p
}

#' Plot % peptides at threshold for per-peptide mode
#'
#' @param x A peppwr_power object
#' @return A ggplot object
#' @keywords internal
plot_power_perpeptide <- function(x) {
  if (x$question == "effect_size" && "effect_curve" %in% names(x$simulations)) {
    plot_data <- x$simulations$effect_curve
    target_power <- x$params$target_power
    min_effect <- x$answer
    proportion_threshold <- x$params$proportion_threshold %||% 0.5

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = effect_size, y = proportion_powered * 100)) +
      ggplot2::geom_line(linewidth = 1, color = "steelblue") +
      ggplot2::geom_point(size = 2, color = "steelblue") +
      ggplot2::geom_vline(xintercept = min_effect, linetype = "dotted", color = "darkgreen") +
      ggplot2::geom_point(
        data = data.frame(effect_size = min_effect, proportion_powered = proportion_threshold),
        ggplot2::aes(x = effect_size, y = proportion_powered * 100),
        size = 4, color = "darkgreen"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Effect Size (fold change)",
        y = paste0("% Peptides at ", target_power * 100, "% Power"),
        title = paste0("Minimum detectable effect: ", min_effect, "-fold"),
        subtitle = paste0(
          "N=", x$params$n_per_group, " per group | ",
          round(proportion_threshold * 100), "% of peptides must reach ",
          round(target_power * 100), "% power"
        )
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 100)) +
      ggplot2::annotate(
        "text",
        x = min_effect,
        y = proportion_threshold * 100 + 5,
        label = paste0(min_effect, "-fold"),
        color = "darkgreen"
      )
  } else if (x$question == "sample_size" && "power_curve" %in% names(x$simulations)) {
    plot_data <- x$simulations$power_curve
    target_power <- x$params$target_power
    proportion_threshold <- x$params$proportion_threshold %||% 0.5

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = n_per_group, y = proportion_powered * 100)) +
      ggplot2::geom_line(linewidth = 1, color = "steelblue") +
      ggplot2::geom_point(size = 2, color = "steelblue") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Sample Size (per group)",
        y = paste0("% Peptides at ", round(target_power * 100), "% Power"),
        title = paste0("Peptide Power Analysis (target: ", target_power * 100, "% power)"),
        subtitle = paste0(
          "Effect size: ", x$params$effect_size, "-fold | ",
          round(proportion_threshold * 100), "% of peptides must reach ",
          round(target_power * 100), "% power"
        )
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 100))
  } else {
    # For find="power", show histogram of per-peptide power values
    if ("peptide_power" %in% names(x$simulations)) {
      plot_data <- data.frame(power = x$simulations$peptide_power)
      plot_data <- plot_data[!is.na(plot_data$power), , drop = FALSE]

      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = power)) +
        ggplot2::geom_histogram(bins = 20, fill = "steelblue", color = "white") +
        ggplot2::geom_vline(xintercept = 0.8, linetype = "dashed", color = "red") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          x = "Power",
          y = "Number of Peptides",
          title = paste0("Per-Peptide Power Distribution (mean: ", round(x$answer * 100), "%)"),
          subtitle = paste0("N=", x$params$n_per_group, ", Effect: ", x$params$effect_size, "-fold")
        ) +
        ggplot2::scale_x_continuous(limits = c(0, 1), labels = scales::percent_format())
    } else {
      # Fallback
      plot_data <- data.frame(power = x$answer)
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = 1, y = power)) +
        ggplot2::geom_point(size = 4, color = "steelblue") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste0("Mean power: ", round(x$answer * 100), "%")
        )
    }
  }

  p
}

fill_scale <- function(name = "name"){
  cols <- RColorBrewer::brewer.pal(length(avail_dists()), "Set3")
  names(cols) <- d2n(avail_dists())
  ggplot2::scale_fill_manual(name = name, values = cols)
}


plot_best <- function(fits_df, fit_col="fits") {
  #TODO check and test df has fit_col 'fits'

  if (! fit_col %in% colnames(fits_df)){
    stop(paste("dataframe does not have a column callled", fit_col, ". Use fit_col = <col name> to specify the nested column containg fit results."))
  }
 ll <- fits_df |>
    dplyr::mutate(best = purrr::map(.data[[fit_col]], \(x) x$dist[which.max(x$loglik)] )) |>
    tidyr::unnest(cols = best) |>
    dplyr::count( best ) |>
    dplyr::mutate(metric = "LogLikelihood")


  la <- fits_df |>
    dplyr::mutate(best = purrr::map(.data[[fit_col]], \(x) x$dist[which.min(x$aic)] )) |>
    tidyr::unnest(cols = best) |>
    dplyr::count( best ) |>
    dplyr::mutate(metric = "AIC")



  dplyr::bind_rows(ll, la) |>
    dplyr::mutate(x = 1) |>
    ggplot2::ggplot() +
    ggplot2::aes(fill=best, y=n, x=x ) +
    ggplot2::geom_bar(position = "stack", stat="identity") +
    ggplot2::facet_wrap(~ metric, ncol = 1) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    fill_scale() +
    ggplot2::labs(x = ggplot2::element_blank(),
                  y = "Phosphopeptide Count",
                  fill = "Model") +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()
    ) +
    ggplot2::xlim(0,2) +
    ggplot2::ggtitle("Best Fitted Models by Metric")

}


plot_failed <- function(fits_df, fit_col = "fits"){

  if (! fit_col %in% colnames(fits_df)){
    stop(paste("dataframe does not have a column callled", fit_col, ". Use fit_col = <col name> to specify the nested column containg fit results."))
  }
  fits_df |>
    tidyr::unnest(cols = .data[[fit_col]]) |>
    dplyr::filter(is.na(aic)) |>
      dplyr::count(dist) |>
    dplyr::mutate(x = 1) |>
    ggplot2::ggplot() +
    ggplot2::aes(fill=dist, y=n, x=x ) +
    ggplot2::geom_bar(position = "stack", stat="identity") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    fill_scale() +
    ggplot2::labs(x = ggplot2::element_blank(),
                  y = "Failed Models Count",
                  fill = "Failed Model") +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()
    ) +
    ggplot2::xlim(0,2) +
    ggplot2::ggtitle("Count of Failed Model Types")


}

fake_leg <- function(name="Dist") {
  p <- data.frame(
    dist = d2n(avail_dists()),
    count = rep(1, length(d2n(avail_dists())) ),
    x = rep(1, length(d2n(avail_dists())) )
  ) |>
    ggplot2::ggplot() +
    ggplot2::aes(count, fill = dist) +
    ggplot2::geom_bar() +
    fill_scale(name=name) +
    ggplot2::theme(legend.position = "bottom")
    cowplot::get_legend(p)
}

evaldist <- function(fits_df, fit_col="fits"){
  a <- plot_best(fits_df, fit_col = fit_col)
  b <- plot_failed(fits_df, fit_col=fit_col)

  p <- cowplot::plot_grid(
    a + ggplot2::theme(legend.position = "none"),
    b + ggplot2::theme(legend.position = "none"),
    rel_heights = c(2, 1),
    nrow = 2
  )


  lgd <-fake_leg() #cowplot::get_legend(a + ggplot2::theme(legend.position = "bottom"))
  cowplot::plot_grid(p, lgd, nrow = 2, rel_heights = c(1, .15))
}

# ============================================================================
# Phase A: Diagnostic Plots (v2)
# ============================================================================

#' Plot density overlay: observed histogram with fitted density curve
#'
#' @param fits A peppwr_fits object
#' @param peptide_id Specific peptide to plot (NULL for multiple)
#' @param n_overlay Number of peptides to overlay when peptide_id is NULL
#'
#' @return A ggplot object
#' @export
plot_density_overlay <- function(fits, peptide_id = NULL, n_overlay = 6) {
  validate_peppwr_fits(fits)

  # Get nested data with fits
  data_df <- fits$data

  if (!is.null(peptide_id)) {
    # Single peptide
    id_col <- names(data_df)[1]  # First column is ID
    data_df <- data_df[data_df[[id_col]] == peptide_id, ]
    if (nrow(data_df) == 0) {
      stop("Peptide ID '", peptide_id, "' not found in fits object")
    }
  } else {
    # Sample up to n_overlay peptides
    if (nrow(data_df) > n_overlay) {
      data_df <- data_df[sample(nrow(data_df), n_overlay), ]
    }
  }

  id_col <- names(data_df)[1]
  value_col <- names(data_df$data[[1]])[1]

  # Build plot data: unnest raw values
  plot_data <- data_df |>
    dplyr::select(tidyr::all_of(id_col), data) |>
    tidyr::unnest(cols = data)
  names(plot_data)[2] <- "value"

  # Get density curves from best fits
  density_data <- purrr::map2_dfr(
    seq_len(nrow(data_df)),
    data_df[[id_col]],
    function(i, pep_id) {
      raw_vals <- data_df$data[[i]][[1]]
      raw_vals <- raw_vals[!is.na(raw_vals)]

      if (length(raw_vals) < 3) return(NULL)

      fit_df <- data_df$fits[[i]]
      best_dist_idx <- which(fits$data[[id_col]] == pep_id)
      if (length(best_dist_idx) == 0) {
        return(NULL)
      }
      best_dist <- fits$best[best_dist_idx[1]]  # Take first match

      if (length(best_dist) == 0 || is.na(best_dist)) {
        return(NULL)
      }

      # Generate x range for density
      min_val <- min(raw_vals)
      max_val <- max(raw_vals)
      if (!is.finite(min_val) || !is.finite(max_val)) return(NULL)

      x_range <- seq(min_val, max_val, length.out = 100)

      # Get density function and compute
      y_vals <- compute_fitted_density(raw_vals, best_dist, x_range)

      if (is.null(y_vals)) return(NULL)

      result <- tibble::tibble(
        peptide_id = pep_id,
        x = x_range,
        density = y_vals
      )
      names(result)[1] <- id_col
      result
    }
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 30,
      fill = "steelblue",
      alpha = 0.5,
      color = "white"
    ) +
    ggplot2::theme_minimal()

  if (nrow(density_data) > 0) {
    p <- p + ggplot2::geom_line(
      data = density_data,
      ggplot2::aes(x = x, y = density),
      color = "darkred",
      linewidth = 1
    )
  }

  if (is.null(peptide_id) && nrow(data_df) > 1) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", id_col)), scales = "free")
  }

  p + ggplot2::labs(
    x = "Abundance",
    y = "Density",
    title = "Observed Data with Fitted Density"
  )
}

#' Compute fitted density values
#'
#' @param raw_vals Raw observed values
#' @param dist_name Distribution name
#' @param x_range X values to compute density at
#'
#' @return Density values or NULL on failure
#' @keywords internal
compute_fitted_density <- function(raw_vals, dist_name, x_range) {
  tryCatch({
    # Map display names back to short names
    short_name <- names(which(d2n(avail_dists()) == dist_name))
    if (length(short_name) == 0) short_name <- dist_name

    # Fit the distribution using univariateML
    fit <- univariateML::model_select(raw_vals, models = short_name, na.rm = TRUE)
    dfunc <- get_dfunc(short_name)
    params <- as.list(fit)
    do.call(dfunc, c(list(x = x_range), params))
  }, error = function(e) {
    NULL
  })
}

#' Get density function for a distribution
#'
#' @param distribution Distribution name
#' @return Density function
#' @keywords internal
get_dfunc <- function(distribution) {
  switch(distribution,
    "norm" = stats::dnorm,
    "gamma" = stats::dgamma,
    "lnorm" = stats::dlnorm,
    "invgamma" = function(x, shape, rate) {
      (rate^shape / gamma(shape)) * x^(-shape - 1) * exp(-rate / x)
    },
    "invgauss" = function(x, mean, shape) {
      sqrt(shape / (2 * pi * x^3)) * exp(-shape * (x - mean)^2 / (2 * mean^2 * x))
    },
    stats::dnorm  # Fallback
  )
}

#' Plot QQ plots for goodness-of-fit
#'
#' @param fits A peppwr_fits object
#' @param peptide_id Specific peptide to plot (NULL for multiple)
#' @param n_plots Number of peptides to plot when peptide_id is NULL
#'
#' @return A ggplot object
#' @export
plot_qq <- function(fits, peptide_id = NULL, n_plots = 6) {
  validate_peppwr_fits(fits)

  data_df <- fits$data

  if (!is.null(peptide_id)) {
    id_col <- names(data_df)[1]
    data_df <- data_df[data_df[[id_col]] == peptide_id, ]
    if (nrow(data_df) == 0) {
      stop("Peptide ID '", peptide_id, "' not found in fits object")
    }
  } else {
    if (nrow(data_df) > n_plots) {
      data_df <- data_df[sample(nrow(data_df), n_plots), ]
    }
  }

  id_col <- names(data_df)[1]

  # Generate QQ data for each peptide
  qq_data <- purrr::map2_dfr(
    seq_len(nrow(data_df)),
    data_df[[id_col]],
    function(i, pep_id) {
      raw_vals <- data_df$data[[i]][[1]]
      raw_vals <- raw_vals[!is.na(raw_vals)]
      best_dist_idx <- which(fits$data[[id_col]] == pep_id)
      if (length(best_dist_idx) == 0) {
        return(NULL)
      }
      best_dist <- fits$best[best_dist_idx[1]]  # Take first match

      if (length(best_dist) == 0 || is.na(best_dist) || length(raw_vals) < 3) {
        return(NULL)
      }

      # Compute theoretical quantiles
      qq_result <- compute_qq_points(raw_vals, best_dist)
      if (is.null(qq_result)) return(NULL)

      result <- tibble::tibble(
        peptide_id = pep_id,
        theoretical = qq_result$theoretical,
        sample = qq_result$sample
      )
      names(result)[1] <- id_col
      result
    }
  )

  if (nrow(qq_data) == 0) {
    # Return a simple placeholder plot when no QQ data is available
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = "No valid QQ data available for selected peptides")
    )
  }

  p <- ggplot2::ggplot(qq_data, ggplot2::aes(x = theoretical, y = sample)) +
    ggplot2::geom_point(alpha = 0.6, color = "steelblue") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Theoretical Quantiles",
      y = "Sample Quantiles",
      title = "QQ Plot: Observed vs Fitted Distribution"
    )

  if (is.null(peptide_id) && length(unique(qq_data[[id_col]])) > 1) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", id_col)), scales = "free")
  }

  p
}

#' Compute QQ points for a distribution fit
#'
#' @param raw_vals Raw observed values
#' @param dist_name Distribution name
#'
#' @return List with theoretical and sample quantiles, or NULL
#' @keywords internal
compute_qq_points <- function(raw_vals, dist_name) {
  tryCatch({
    short_name <- names(which(d2n(avail_dists()) == dist_name))
    if (length(short_name) == 0) short_name <- dist_name

    fit <- univariateML::model_select(raw_vals, models = short_name, na.rm = TRUE)
    qfunc <- get_qfunc(short_name)
    params <- as.list(fit)

    n <- length(raw_vals)
    p <- stats::ppoints(n)
    theoretical <- do.call(qfunc, c(list(p = p), params))
    sample <- sort(raw_vals)

    list(theoretical = theoretical, sample = sample)
  }, error = function(e) {
    NULL
  })
}

#' Get quantile function for a distribution
#'
#' @param distribution Distribution name
#' @return Quantile function
#' @keywords internal
get_qfunc <- function(distribution) {
  switch(distribution,
    "norm" = stats::qnorm,
    "gamma" = stats::qgamma,
    "lnorm" = stats::qlnorm,
    "invgamma" = function(p, shape, rate) rate / stats::qgamma(1 - p, shape = shape, rate = 1),
    "invgauss" = function(p, mean, shape) {
      # Approximate via Newton-Raphson or use median approximation
      mean * (1 + stats::qnorm(p) * sqrt(mean / shape))
    },
    stats::qnorm  # Fallback
  )
}

#' Plot power heatmap: N x effect size grid
#'
#' @param distribution Distribution name
#' @param params Distribution parameters
#' @param n_range Range of sample sizes (vector of length 2)
#' @param effect_range Range of effect sizes (vector of length 2)
#' @param n_steps Number of grid points per dimension
#' @param n_sim Number of simulations per grid cell
#' @param test Statistical test to use
#' @param alpha Significance level
#'
#' @return A ggplot object
#' @export
plot_power_heatmap <- function(distribution, params, n_range, effect_range,
                                n_steps = 6, n_sim = 100, test = "wilcoxon",
                                alpha = 0.05) {
  # Generate grid
  n_vals <- seq(n_range[1], n_range[2], length.out = n_steps)
  n_vals <- unique(round(n_vals))
  effect_vals <- seq(effect_range[1], effect_range[2], length.out = n_steps)

  grid <- expand.grid(n_per_group = n_vals, effect_size = effect_vals)

  # Compute power for each cell
  grid$power <- purrr::map2_dbl(
    grid$n_per_group,
    grid$effect_size,
    function(n, eff) {
      run_power_sim(
        distribution = distribution,
        params = params,
        n_per_group = n,
        effect_size = eff,
        alpha = alpha,
        test = test,
        n_sim = n_sim
      )
    }
  )

  ggplot2::ggplot(grid, ggplot2::aes(x = n_per_group, y = effect_size, fill = power)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "white",
      mid = "steelblue",
      high = "darkblue",
      midpoint = 0.5,
      limits = c(0, 1),
      labels = scales::percent_format()
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Sample Size (per group)",
      y = "Effect Size (fold change)",
      fill = "Power",
      title = "Power Heatmap",
      subtitle = paste("Distribution:", distribution, "| Test:", test)
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::percent(power, accuracy = 1)),
      color = ifelse(grid$power > 0.5, "white", "black"),
      size = 3
    )
}

#' Plot power vs effect size at fixed N
#'
#' @param power_result A peppwr_power object
#' @param effect_range Range of effect sizes to explore
#' @param n_steps Number of effect size values to compute
#' @param n_sim Number of simulations per point
#'
#' @return A ggplot object
#' @export
plot_power_vs_effect <- function(power_result, effect_range, n_steps = 10, n_sim = NULL) {
  validate_peppwr_power(power_result)

  # Get n_per_group from the result
  n_per_group <- power_result$params$n_per_group
  if (is.null(n_per_group)) {
    n_per_group <- power_result$answer  # If question was sample_size
  }

  # Get simulation parameters
  if (is.null(n_sim)) {
    n_sim <- power_result$simulations$n_sim %||% 100
  }

  # Generate effect size range
  effect_vals <- seq(effect_range[1], effect_range[2], length.out = n_steps)

  # Compute power at each effect size
  if (power_result$mode == "aggregate") {
    distribution <- power_result$params$distribution
    params <- power_result$params$params

    power_curve <- purrr::map_dbl(effect_vals, function(eff) {
      run_power_sim(
        distribution = distribution,
        params = params,
        n_per_group = n_per_group,
        effect_size = eff,
        alpha = power_result$params$alpha %||% 0.05,
        test = power_result$params$test %||% "wilcoxon",
        n_sim = n_sim
      )
    })
  } else {
    # Per-peptide mode: compute mean power across peptides
    power_curve <- purrr::map_dbl(effect_vals, function(eff) {
      # Re-run power analysis with different effect size
      # For simplicity, use aggregate mode based on median params
      0.5  # Placeholder - will be computed properly
    })

    # Actually recompute using the fits data
    fits_data <- power_result$simulations$fits_data
    if (!is.null(fits_data)) {
      power_curve <- purrr::map_dbl(effect_vals, function(eff) {
        pep_powers <- compute_perpeptide_power(
          fits_data,
          n_per_group = n_per_group,
          effect_size = eff,
          alpha = power_result$params$alpha %||% 0.05,
          test = power_result$params$test %||% "wilcoxon",
          n_sim = n_sim
        )
        mean(pep_powers, na.rm = TRUE)
      })
    }
  }

  plot_data <- tibble::tibble(
    effect_size = effect_vals,
    power = power_curve
  )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = effect_size, y = power)) +
    ggplot2::geom_line(linewidth = 1, color = "steelblue") +
    ggplot2::geom_point(size = 2, color = "steelblue") +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    ggplot2::theme_minimal() +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
    ggplot2::labs(
      x = "Effect Size (fold change)",
      y = "Power",
      title = paste0("Power vs Effect Size (N=", n_per_group, " per group)"),
      subtitle = paste("Mode:", power_result$mode)
    )
}

#' Compute per-peptide power (helper for plot_power_vs_effect)
#'
#' @keywords internal
compute_perpeptide_power <- function(fits_data, n_per_group, effect_size, alpha, test, n_sim) {
  # Simplified: return vector of powers
  purrr::map_dbl(seq_len(nrow(fits_data)), function(i) {
    fit_df <- fits_data$fits[[i]]
    valid_fits <- fit_df[!is.na(fit_df$aic), ]
    if (nrow(valid_fits) == 0) return(NA_real_)

    best_dist <- valid_fits$dist[which.min(valid_fits$aic)]
    short_name <- names(which(d2n(avail_dists()) == best_dist))
    if (length(short_name) == 0) short_name <- "norm"

    raw_vals <- fits_data$data[[i]][[1]]
    raw_vals <- raw_vals[!is.na(raw_vals)]

    tryCatch({
      fit <- univariateML::model_select(raw_vals, models = short_name, na.rm = TRUE)
      params <- as.list(fit)

      run_power_sim(
        distribution = short_name,
        params = params,
        n_per_group = n_per_group,
        effect_size = effect_size,
        alpha = alpha,
        test = test,
        n_sim = n_sim
      )
    }, error = function(e) NA_real_)
  })
}

#' Plot distribution of fitted parameters across peptidome
#'
#' @param fits A peppwr_fits object
#'
#' @return A ggplot object
#' @export
plot_param_distribution <- function(fits) {
  validate_peppwr_fits(fits)

  # Extract AIC values from all fits as a summary statistic
  # (Parameter extraction varies by distribution, so we use AIC as proxy)
  id_col <- names(fits$data)[1]

  # Collect best AIC per peptide and distribution info
  param_data <- purrr::map2_dfr(
    seq_len(length(fits$best)),
    fits$best,
    function(i, best_dist) {
      if (is.na(best_dist)) return(NULL)

      fit_df <- fits$fits[[i]]
      valid_fits <- fit_df[!is.na(fit_df$aic), ]
      if (nrow(valid_fits) == 0) return(NULL)

      best_row <- valid_fits[which.min(valid_fits$aic), ]

      tibble::tibble(
        best_dist = best_dist,
        aic = best_row$aic,
        loglik = best_row$loglik
      )
    }
  )

  if (nrow(param_data) == 0) {
    stop("No valid fits to plot")
  }

  # Create faceted plot showing AIC distribution by best distribution
  p1 <- ggplot2::ggplot(param_data, ggplot2::aes(x = aic, fill = best_dist)) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.7) +
    ggplot2::facet_wrap(~ best_dist, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "AIC",
      y = "Count",
      fill = "Distribution",
      title = "Distribution of Fit Quality (AIC) Across Peptidome",
      subtitle = paste(nrow(param_data), "peptides with valid fits")
    ) +
    ggplot2::theme(legend.position = "none")

  p1
}

#' Null-coalescing operator
#' @param x Value to check for NULL
#' @param y Default value if x is NULL
#' @return x if not NULL, otherwise y
#' @name null_coalesce
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

# Declare global variables to avoid R CMD check notes
utils::globalVariables(c("value", "density", "theoretical", "best_dist",
                         "na_rate", "mnar_score", "mean_abundance"))

#' Plot missingness statistics
#'
#' Shows distribution of NA rates, MNAR scores, and relationship between
#' abundance and missingness across peptides.
#'
#' @param fits A peppwr_fits object
#'
#' @return A ggplot object or gtable (combined panels)
#' @export
plot_missingness <- function(fits) {
  validate_peppwr_fits(fits)

  if (is.null(fits$missingness)) {
    stop("Fits object does not contain missingness data")
  }


  miss_data <- fits$missingness

  # Panel 1: NA rate distribution
  p1 <- ggplot2::ggplot(miss_data, ggplot2::aes(x = na_rate)) +
    ggplot2::geom_histogram(
      bins = 20,
      fill = "steelblue",
      color = "white",
      alpha = 0.7
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "NA Rate",
      y = "Number of Peptides",
      title = "Distribution of Missing Data Rates"
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(),
      limits = c(0, max(miss_data$na_rate, 0.1) * 1.1)
    )

  # MNAR score plot (only if there's variation)
  mnar_scores <- miss_data$mnar_score[!is.na(miss_data$mnar_score)]

  if (length(mnar_scores) > 1 && stats::var(mnar_scores) > 0.01) {
    # Panel 2: MNAR score histogram
    p2 <- ggplot2::ggplot(miss_data[!is.na(miss_data$mnar_score), ],
                          ggplot2::aes(x = mnar_score)) +
      ggplot2::geom_histogram(
        bins = 20,
        fill = "darkorange",
        color = "white",
        alpha = 0.7
      ) +
      ggplot2::geom_vline(xintercept = 2, linetype = "dashed", color = "red") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "MNAR Score",
        y = "Number of Peptides",
        subtitle = "Score > 2 suggests MNAR pattern (red line)"
      ) +
      ggplot2::annotate(
        "text",
        x = 2.2,
        y = Inf,
        label = "MNAR threshold",
        hjust = 0,
        vjust = 1.5,
        color = "red",
        size = 3
      )

    # Panel 3: Abundance vs NA rate scatter (only peptides with missing data)
    miss_subset <- miss_data[miss_data$na_rate > 0 & !is.na(miss_data$mean_abundance), ]

    if (nrow(miss_subset) > 0) {
      p3 <- ggplot2::ggplot(miss_subset,
                            ggplot2::aes(x = mean_abundance, y = na_rate)) +
        ggplot2::geom_point(ggplot2::aes(color = mnar_score), alpha = 0.6, size = 2) +
        ggplot2::scale_color_gradient2(
          low = "steelblue", mid = "gray", high = "red",
          midpoint = 0, name = "MNAR\nScore"
        ) +
        ggplot2::scale_x_log10() +
        ggplot2::scale_y_continuous(labels = scales::percent_format()) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          x = "Mean Abundance (log scale)",
          y = "NA Rate",
          subtitle = "Negative trend suggests MNAR pattern"
        )

      # Add trend line if enough points
      if (nrow(miss_subset) >= 5) {
        p3 <- p3 + ggplot2::geom_smooth(
          method = "loess",
          se = FALSE,
          color = "black",
          linetype = "dashed",
          formula = y ~ x
        )
      }

      # Combine 3 panels
      cowplot::plot_grid(p1, p2, p3, nrow = 3, rel_heights = c(1, 1, 1.2))
    } else {
      # No peptides with missing data - just show 2 panels
      cowplot::plot_grid(p1, p2, nrow = 2, rel_heights = c(1, 1))
    }
  } else {
    # Just show NA rate plot
    p1 +
      ggplot2::labs(subtitle = paste0(
        sum(miss_data$na_rate > 0), " peptides with missing data"
      ))
  }
}


