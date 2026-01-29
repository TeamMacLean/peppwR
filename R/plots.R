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
  if (x$question == "sample_size" && "power_curve" %in% names(x$simulations)) {
    plot_data <- x$simulations$power_curve
    target_power <- x$params$target_power

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = n_per_group, y = proportion_powered * 100)) +
      ggplot2::geom_line(linewidth = 1, color = "steelblue") +
      ggplot2::geom_point(size = 2, color = "steelblue") +
      ggplot2::geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Sample Size (per group)",
        y = "% Peptides at Target Power",
        title = paste0("Peptide Power Analysis (target: ", target_power * 100, "% power)"),
        subtitle = paste0("Effect size: ", x$params$effect_size, "-fold")
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


