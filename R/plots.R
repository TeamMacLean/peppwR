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


