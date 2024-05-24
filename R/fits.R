#' Fit a Single Distribution
#'
#' @description This fits a specified distribution to the first column of a data frame.
#'
#' @param df The data frame containing the data to fit.
#' @param dist The distribution name or a list of distribution models to try fitting. One of "nbinom", "gamma", "snorm", "norm",
#'        "invgamma", "invgauss", "lnorm", "lgamma",  and "pareto"
#'
#'@return A list containing the fit result or the distribution name in case there is an error.
#'
#' @export
single_fit <- function(df, dist){
  if (dist %in% c("nbinom")){
    result <- tryCatch(
      {fitdistrplus::fitdist(df[,1], dist)},
      error = \(x){result <- list(dist = dist) }
    )
  } else {
    result <- tryCatch(
      {
        univariateML::model_select(df[,1], models=dist, na.rm=TRUE)
        },
      error = \(x){ result <- list(dist = dist) }
    )
  }

}

#' Available Distributions
#'
#' @description This provides a vector of available distributions for fitting.
#'
#' @return A character vector of available distribution names.
#'
#'
#'
#' @export
avail_dists <- function() {
  c("gamma", "norm",
    "snorm", "invgamma",
    "invgauss", "lnorm",
    "lgamma", "pareto",
    "nbinom")
}


#' Distribution Name to Normalized Name
#'
#' @description This converts distribution tags to names that humans can read and understand.
#'
#' @param tag A character string representing the distribution tag.
#'
#' @return A character string with the human-readable distribution name.
#'
#'
#' @export
d2n <- function(tag){
  v <- c("Gamma", "Normal", "Skew Normal", "InvGamma", "Inverse Gaussian",
           "Lognormal", "Log Gamma", "Pareto", "Negative Binomial")
  names(v) <- avail_dists()
  v[tag]
}


#' Fit Multiple Distributions
#'
#' @description This is used to fits multiple distributions to the first column of a data frame.
#'
#' @param df A data frame with the data to fit.
#'
#' @return A data frame with the fit results for each distribution.
#'
#'
#' @export
do_fits <- function(df) {
  dists <- avail_dists()

  fits <- lapply(dists, \(d) single_fit(df, d))
  names(fits) <- dists
  dplyr::bind_rows( lapply(fits, squash_fits) )

}


#' Squash Fit Results
#'
#' @description This processes fit results into a tidy format.
#'
#' @param fit A list containing the fit result.
#'
#' @return A tibble with the distribution name, log-likelihood, and Akaike Information Criterion (AIC).
#'
#' @export
squash_fits <- function(fit){

  if (length(fit) == 1){
    return( tibble::tibble(
      dist = d2n(fit$dist),
      loglik = as.numeric(NA),
      aic = as.numeric(NA)
    ) )
  }
  if (is(fit,"fitdist")) {
      parse_fitdist(fit)
  } else {
      parse_univariateML(fit)
    }

}



#' Parse Fit Results from fitdistrplus
#'
#' @description This converts fit results from the `fitdistrplus` package into a tidy format.
#'
#' @param fit An object returned by `fitdistrplus::fitdist`.
#'
#' @return A tibble with the distribution name, log-likelihood, and AIC.
#'
#'
#' @export
parse_fitdist <- function(fit){
  tibble::tibble(
    dist = fit$distname,
    loglik = fit$loglik,
    aic = fit$aic
  )
}


#' Parse Fit Results from univariateML
#'
#' @description This converts fit results from the `univariateML` package into a tidy format.
#'
#' @param fit An object returned by `univariateML::model_select`.
#'
#' @return A tibble with the distribution name, log-likelihood, and AIC.
#'
#' @export
parse_univariateML <- function(fit){
  tibble::tibble(
    dist = attr(fit, "model"),
    loglik = attr(fit, "logLik"),
    aic = AIC(fit)
  )
}
