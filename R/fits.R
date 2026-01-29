# Distribution fitting functions for peppwR

#' Fit distributions to peptide abundance data
#'
#' @param data A data frame containing peptide abundance data
#' @param id Column name for peptide identifier
#' @param group Column name for group/condition
#' @param value Column name for abundance values
#' @param distributions Which distributions to fit: "all" or a character vector
#'
#' @return A peppwr_fits object
#' @export
fit_distributions <- function(data, id, group, value, distributions = "all") {
  the_call <- match.call()

  # Nest data by peptide and group
  nested <- tidyr::nest(data, .by = tidyr::all_of(c(id, group)), data = tidyr::all_of(value))

  # Fit distributions to each peptide
  nested <- dplyr::mutate(
    nested,
    fits = purrr::map(.data$data, do_fits)
  )

  # Extract best fit for each peptide (by AIC)
  best_fits <- purrr::map_chr(nested$fits, function(fit_df) {
    valid_fits <- fit_df[!is.na(fit_df$aic), ]
    if (nrow(valid_fits) == 0) {
      return(NA_character_)
    }
    valid_fits$dist[which.min(valid_fits$aic)]
  })

  new_peppwr_fits(
    data = nested,
    fits = nested$fits,
    best = best_fits,
    call = the_call
  )
}

single_fit <- function(df, dist){
  # Extract numeric vector from data frame (df[[1]] gives vector, df[,1] gives tibble)
  x <- df[[1]]

  if (dist %in% c("nbinom")){
    result <- tryCatch(
      {fitdistrplus::fitdist(x, dist)},
      error = \(e){list(dist = dist) }
    )
  } else {
    result <- tryCatch(
      {
        univariateML::model_select(x, models = dist, na.rm = TRUE)
      },
      error = \(e){ list(dist = dist) }
    )
  }
  result
}


avail_dists <- function() {
  c("gamma", "norm",
    "snorm", "invgamma",
    "invgauss", "lnorm",
    "lgamma", "pareto",
    "nbinom")
}

d2n <- function(tag){
  v <- c("Gamma", "Normal", "Skew Normal", "InvGamma", "Inverse Gaussian",
           "Lognormal", "Log Gamma", "Pareto", "Negative Binomial")
  names(v) <- avail_dists()
  v[tag]
}

do_fits <- function(df) {
  dists <- avail_dists()

  fits <- lapply(dists, \(d) single_fit(df, d))
  names(fits) <- dists
  dplyr::bind_rows( lapply(fits, squash_fits) )

}

squash_fits <- function(fit){

  if (length(fit) == 1){
    return( tibble::tibble(
      dist = d2n(fit$dist),
      loglik = as.numeric(NA),
      aic = as.numeric(NA)
    ) )
  }
  if (methods::is(fit, "fitdist")) {
      parse_fitdist(fit)
  } else {
      parse_univariateML(fit)
    }

}

parse_fitdist <- function(fit){
  tibble::tibble(
    dist = fit$distname,
    loglik = fit$loglik,
    aic = fit$aic
  )
}

parse_univariateML <- function(fit){
  tibble::tibble(
    dist = attr(fit, "model"),
    loglik = attr(fit, "logLik"),
    aic = stats::AIC(fit)
  )
}
