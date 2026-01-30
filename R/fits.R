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

  # Compute missingness statistics for each peptide
  missingness <- purrr::map_dfr(nested$data, function(df) {
    values <- df[[1]]
    miss_stats <- compute_missingness(values)
    tibble::tibble(
      n_total = miss_stats$n_total,
      n_missing = miss_stats$n_missing,
      na_rate = miss_stats$na_rate,
      mnar_score = miss_stats$mnar_score,
      mnar_pvalue = miss_stats$mnar_pvalue
    )
  })

  # Add peptide IDs to missingness table
  missingness <- dplyr::bind_cols(
    tibble::tibble(peptide_idx = seq_len(nrow(nested))),
    missingness
  )

  new_peppwr_fits(
    data = nested,
    fits = nested$fits,
    best = best_fits,
    call = the_call,
    missingness = missingness
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

# ============================================================================
# Phase C: Missing Data Handling
# ============================================================================

#' Compute missingness statistics for a vector of values
#'
#' Calculates NA rate and MNAR (Missing Not At Random) score.
#' MNAR detection uses the observation that under MCAR, the mean rank
#' of observed values should be (n+1)/2. Under MNAR with low values
#' missing, the mean rank will be higher.
#'
#' @param values Numeric vector (may contain NAs)
#'
#' @return List with:
#'   - n_total: Total number of values
#'   - n_missing: Number of NA values
#'   - na_rate: Proportion missing (0-1)
#'   - mnar_score: Normalized deviation from expected mean rank
#'     (positive = low values more likely missing)
#'   - mnar_pvalue: P-value for MNAR pattern (Wilcoxon one-sample test)
#'
#' @export
compute_missingness <- function(values) {
  n_total <- length(values)
  n_missing <- sum(is.na(values))
  na_rate <- n_missing / n_total

  # MNAR detection
  observed_values <- values[!is.na(values)]
  n_obs <- length(observed_values)

  if (n_obs < 2 || n_missing == 0) {
    # Can't compute MNAR score without both observed and missing
    mnar_score <- NA_real_
    mnar_pvalue <- NA_real_
  } else {
    # Under MCAR, the expected mean rank of observed values is (n_total + 1) / 2
    # Under MNAR (low missing), observed values will have higher ranks
    all_ranks <- rank(c(observed_values, rep(NA, n_missing)), na.last = "keep")
    observed_ranks <- all_ranks[1:n_obs]

    expected_mean_rank <- (n_total + 1) / 2
    observed_mean_rank <- mean(observed_ranks, na.rm = TRUE)

    # Standard error of mean rank under MCAR
    se_rank <- sqrt(n_total * (n_total + 1) / (12 * n_obs))

    mnar_score <- (observed_mean_rank - expected_mean_rank) / se_rank

    # P-value using Wilcoxon signed-rank test (one-sample)
    mnar_pvalue <- tryCatch({
      test <- stats::wilcox.test(
        observed_ranks,
        mu = expected_mean_rank,
        alternative = "greater"
      )
      test$p.value
    }, error = function(e) NA_real_)
  }

  list(
    n_total = n_total,
    n_missing = n_missing,
    na_rate = na_rate,
    mnar_score = mnar_score,
    mnar_pvalue = mnar_pvalue
  )
}
