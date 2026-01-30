# Distribution fitting functions for peppwR

#' Get distributions for a preset
#'
#' @param preset One of "continuous" (default), "counts", or "all"
#' @return Character vector of distribution names
#' @export
get_distribution_preset <- function(preset = "continuous") {
  continuous_dists <- c("gamma", "norm", "snorm", "invgamma",
                        "invgauss", "lnorm", "lgamma", "pareto")
  count_dists <- c("nbinom")

  switch(preset,
    "continuous" = continuous_dists,
    "counts" = c(continuous_dists, count_dists),
    "all" = c(continuous_dists, count_dists),
    stop("Unknown preset: ", preset, ". Use 'continuous', 'counts', or 'all'")
  )
}

#' Check if data appears to be count data (non-negative integers)
#'
#' @param x Numeric vector
#' @return TRUE if data looks like counts, FALSE otherwise
#' @export
is_count_data <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(FALSE)
  all(x >= 0) && all(x == floor(x))
}

#' Fit distributions to peptide abundance data
#'
#' @param data A data frame containing peptide abundance data
#' @param id Column name for peptide identifier
#' @param group Column name for group/condition
#' @param value Column name for abundance values
#' @param distributions Which distributions to fit: "continuous" (default),
#'   "counts", "all", or a character vector of distribution names
#'
#' @return A peppwr_fits object
#' @export
fit_distributions <- function(data, id, group, value, distributions = "continuous") {
  the_call <- match.call()

  # Handle distributions parameter
  if (length(distributions) == 1 && distributions %in% c("continuous", "counts", "all")) {
    dist_list <- get_distribution_preset(distributions)

    # Auto-detection messaging for "all"
    if (distributions == "all") {
      # Sample some values to check data type
      sample_values <- data[[value]][1:min(100, nrow(data))]
      if (!is_count_data(sample_values)) {
        message("i Skipping nbinom: data appears continuous (non-integer values detected)")
        dist_list <- setdiff(dist_list, "nbinom")
      }
    }
  } else {
    # User provided explicit list
    dist_list <- distributions
  }

  # Nest data by peptide and group
  nested <- tidyr::nest(data, .by = tidyr::all_of(c(id, group)), data = tidyr::all_of(value))

  # Fit distributions to each peptide
  nested <- dplyr::mutate(
    nested,
    fits = purrr::map(.data$data, \(df) do_fits(df, dists = dist_list))
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
      mean_abundance = mean(values, na.rm = TRUE)
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


avail_dists <- function(include_counts = FALSE) {
  if (include_counts) {
    get_distribution_preset("all")
  } else {
    get_distribution_preset("continuous")
  }
}

d2n <- function(tag){
  # Full name mapping for all distributions (including counts)
  all_dists <- c("gamma", "norm", "snorm", "invgamma",
                 "invgauss", "lnorm", "lgamma", "pareto", "nbinom")
  v <- c("Gamma", "Normal", "Skew Normal", "InvGamma", "Inverse Gaussian",
           "Lognormal", "Log Gamma", "Pareto", "Negative Binomial")
  names(v) <- all_dists
  v[tag]
}

do_fits <- function(df, dists = NULL) {
  if (is.null(dists)) {
    dists <- get_distribution_preset("continuous")
  }

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
    dist = d2n(fit$distname),
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

#' Get peptides showing MNAR pattern
#'
#' Returns peptides where the MNAR score exceeds a threshold, indicating
#' that low values are systematically missing (Missing Not At Random).
#'
#' @param fits A peppwr_fits object
#' @param threshold MNAR score threshold (default 2, corresponds to ~95% confidence)
#'
#' @return A tibble with columns:
#'   - peptide_id: Peptide identifier
#'   - condition: Group/condition
#'   - na_rate: Proportion of missing values
#'   - mnar_score: Z-score for MNAR pattern
#'   - mean_abundance: Mean of observed values
#'   Sorted by mnar_score descending. Returns empty tibble if no peptides
#'   exceed threshold.
#'
#' @export
get_mnar_peptides <- function(fits, threshold = 2) {
  validate_peppwr_fits(fits)

  if (is.null(fits$missingness)) {
    stop("Fits object does not contain missingness data")
  }

  # Get peptide IDs from fits$data
  id_col <- names(fits$data)[1]
  group_col <- names(fits$data)[2]

  # Build result tibble
  result <- tibble::tibble(
    peptide_id = fits$data[[id_col]],
    condition = fits$data[[group_col]],
    na_rate = fits$missingness$na_rate,
    mnar_score = fits$missingness$mnar_score,
    mean_abundance = fits$missingness$mean_abundance
  )

  # Filter to threshold and sort
  result <- result[!is.na(result$mnar_score) & result$mnar_score > threshold, ]
  result <- result[order(-result$mnar_score), ]

  result
}

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
#'   - mnar_score: Z-score measuring MNAR pattern. Positive values indicate
#'     low values are more likely to be missing. Values > 2 suggest MNAR.
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
  }

  list(
    n_total = n_total,
    n_missing = n_missing,
    na_rate = na_rate,
    mnar_score = mnar_score
  )
}
