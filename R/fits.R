
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
  if (is(fit,"fitdist")) {
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
    aic = AIC(fit)
  )
}
