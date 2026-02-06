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
#' Fits candidate distributions to each peptide's abundance values and selects
#' the best fit by AIC. Also computes missingness statistics including
#' dataset-level MNAR detection.
#'
#' @section Missingness Tracking:
#' The returned object includes:
#' - Per-peptide NA rates and MNAR scores (in `$missingness`)
#' - Dataset-level MNAR correlation (in `$dataset_mnar`)
#'
#' The dataset-level MNAR metric correlates log(mean_abundance) with NA rate

#' across peptides. A negative correlation indicates low-abundance peptides
#' have more missing values - typical of detection-limit-driven MNAR.
#'
#' Print the result to see both metrics. For small sample sizes (N < 15),
#' the dataset-level correlation is more reliable than per-peptide scores.
#'
#' @param data A data frame containing peptide abundance data
#' @param id Column name for peptide identifier
#' @param group Column name for group/condition
#' @param value Column name for abundance values
#' @param distributions Which distributions to fit: "continuous" (default),
#'   "counts", "all", or a character vector of distribution names
#'
#' @return A peppwr_fits object containing:
#'   - `$data`: Nested tibble with original data and fit results
#'   - `$best`: Best-fitting distribution for each peptide
#'   - `$missingness`: Per-peptide missingness statistics
#'   - `$dataset_mnar`: Dataset-level MNAR correlation and interpretation
#'
#' @seealso [compute_dataset_mnar()] for dataset-level MNAR details,
#'   [plot_missingness()] to visualize missingness patterns
#'
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
      mean_abundance = mean(values, na.rm = TRUE)
    )
  })

  # Add peptide IDs to missingness table
  missingness <- dplyr::bind_cols(
    tibble::tibble(peptide_idx = seq_len(nrow(nested))),
    missingness
  )

  # Compute dataset-level MNAR metric
  dataset_mnar <- compute_dataset_mnar(missingness)

  new_peppwr_fits(
    data = nested,
    fits = nested$fits,
    best = best_fits,
    call = the_call,
    missingness = missingness,
    dataset_mnar = dataset_mnar
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

#' Compute missingness statistics for a vector of values
#'
#' Calculates the number and proportion of missing (NA) values.
#'
#' @param values Numeric vector (may contain NAs)
#'
#' @return List with:
#'   - n_total: Total number of values
#'   - n_missing: Number of NA values
#'   - na_rate: Proportion missing (0-1)
#'
#' @seealso [compute_dataset_mnar()] for dataset-level MNAR detection
#'
#' @export
compute_missingness <- function(values) {
  n_total <- length(values)
  n_missing <- sum(is.na(values))
  na_rate <- n_missing / n_total

  list(
    n_total = n_total,
    n_missing = n_missing,
    na_rate = na_rate
  )
}

#' Compute dataset-level MNAR metric
#'
#' Calculates Spearman correlation between log(mean_abundance) and NA rate
#' across peptides to detect whether low-abundance peptides have more missing
#' values than high-abundance peptides.
#'
#' @section MNAR Detection:
#' MNAR (Missing Not At Random) in mass spectrometry typically manifests as
#' low-abundance peptides having higher rates of missing values due to
#' detection limits. This function detects this pattern by correlating
#' mean abundance with NA rate across all peptides.
#'
#' A negative correlation indicates that low-abundance peptides have more
#' missing values - the hallmark of detection-limit-driven MNAR.
#'
#' @section Interpretation:
#' | Correlation (r) | Interpretation |
#' |-----------------|----------------|
#' | r > -0.1 | No evidence of MNAR |
#' | -0.3 < r < -0.1 | Weak evidence |
#' | -0.5 < r < -0.3 | Moderate evidence |
#' | r < -0.5 | Strong evidence of MNAR |
#'
#' @param missingness A tibble with columns na_rate and mean_abundance
#'   (as produced by fit_distributions)
#'
#' @return A list with:
#'   - correlation: Spearman correlation coefficient (negative = MNAR pattern)
#'   - p_value: p-value for the correlation
#'   - n_peptides: Number of peptides with missing data used in calculation
#'   - interpretation: Human-readable interpretation string
#'
#' @export
compute_dataset_mnar <- function(missingness) {
  # Filter to peptides with missing data and valid mean abundance
  subset <- missingness[
    missingness$na_rate > 0 &
    !is.na(missingness$mean_abundance) &
    !is.nan(missingness$mean_abundance) &
    missingness$mean_abundance > 0,
  ]

  n_peptides <- nrow(subset)

  # Need at least 5 peptides for meaningful correlation

  if (n_peptides < 5) {
    return(list(
      correlation = NA_real_,
      p_value = NA_real_,
      n_peptides = n_peptides,
      interpretation = "Insufficient peptides with missing data (< 5)"
    ))
  }

  # Compute Spearman correlation between log(mean_abundance) and na_rate
  cor_test <- stats::cor.test(
    log(subset$mean_abundance),
    subset$na_rate,
    method = "spearman",
    exact = FALSE
  )

  correlation <- as.numeric(cor_test$estimate)
  p_value <- cor_test$p.value

  # Generate interpretation
  interpretation <- interpret_dataset_mnar(correlation, p_value, n_peptides)

  list(
    correlation = correlation,
    p_value = p_value,
    n_peptides = n_peptides,
    interpretation = interpretation
  )
}

#' Interpret dataset-level MNAR result
#'
#' @param correlation Spearman correlation coefficient
#' @param p_value p-value for the correlation
#' @param n_peptides Number of peptides in calculation
#'
#' @return Interpretation string
#' @keywords internal
interpret_dataset_mnar <- function(correlation, p_value, n_peptides) {
  if (is.na(correlation)) {
    return("Unable to compute correlation")
  }

  # Determine strength based on correlation magnitude
  abs_cor <- abs(correlation)
  if (abs_cor < 0.1) {
    strength <- "No evidence"
  } else if (abs_cor < 0.3) {
    strength <- "Weak evidence"
  } else if (abs_cor < 0.5) {
    strength <- "Moderate evidence"
  } else {
    strength <- "Strong evidence"
  }

  # Direction
  if (correlation < -0.1) {
    direction <- "low-abundance peptides have more missing values"
  } else if (correlation > 0.1) {
    direction <- "high-abundance peptides have more missing values (unusual)"
  } else {
    direction <- "no clear abundance-missingness relationship"
  }

  # Significance
  if (!is.na(p_value) && p_value < 0.05) {
    sig <- sprintf("(p = %.2g)", p_value)
  } else if (!is.na(p_value)) {
    sig <- sprintf("(p = %.2g, not significant)", p_value)
  } else {
    sig <- ""
  }

  paste(strength, "of MNAR:", direction, sig)
}
