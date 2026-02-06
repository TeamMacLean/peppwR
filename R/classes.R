# S3 class constructors and validators for peppwR

# --- peppwr_fits class ---

#' Create a new peppwr_fits object
#'
#' @param data Original data (nested tibble)
#' @param fits Fit results per peptide (list of tibbles with dist, loglik, aic)
#' @param best Best-fitting distribution per peptide (character vector)
#' @param call Original function call
#' @param missingness Tibble with missingness statistics per peptide (optional)
#'
#' @return A peppwr_fits object
#' @export
new_peppwr_fits <- function(data, fits, best, call, missingness = NULL) {
  structure(
    list(
      data = data,
      fits = fits,
      best = best,
      call = call,
      missingness = missingness
    ),
    class = "peppwr_fits"
  )
}

#' Validate a peppwr_fits object
#'
#' @param x Object to validate
#'
#' @return The validated object (invisibly)
#' @export
validate_peppwr_fits <- function(x) {
  if (!inherits(x, "peppwr_fits")) {
    stop("Object must be of class 'peppwr_fits'", call. = FALSE)
  }

  required <- c("data", "fits", "best", "call")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop(
      "peppwr_fits object missing required components: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(x)
}

#' Print method for peppwr_fits
#'
#' @param x A peppwr_fits object
#' @param ... Additional arguments (ignored)
#'
#' @return The object invisibly
#' @export
#' @method print peppwr_fits
print.peppwr_fits <- function(x, ...) {
  validate_peppwr_fits(x)

  n_peptides <- length(x$best)
  peptide_word <- if (n_peptides == 1) "peptide" else "peptides"

  cat("peppwr_fits object\n")
  cat("------------------\n")
  cat(sprintf("%d %s fitted\n\n", n_peptides, peptide_word))

  # Count best distributions
  best_counts <- table(x$best)
  cat("Best fit distribution counts:\n")
  for (dist_name in names(best_counts)) {
    cat(sprintf("  %s: %d\n", dist_name, best_counts[dist_name]))
  }

  # Count failures (NA in best)
  n_failed <- sum(is.na(x$best))
  if (n_failed > 0) {
    cat(sprintf("\nFailed fits: %d (%.1f%%)\n", n_failed, 100 * n_failed / n_peptides))
  }

  # Missingness summary
  if (!is.null(x$missingness) && nrow(x$missingness) > 0) {
    total_na <- sum(x$missingness$n_missing)
    total_obs <- sum(x$missingness$n_total)
    peps_with_na <- sum(x$missingness$na_rate > 0)

    if (total_na > 0) {
      cat(sprintf("\nMissingness: %d/%d values NA (%.1f%%)\n",
                  total_na, total_obs, 100 * total_na / total_obs))
      cat(sprintf("Peptides with missing data: %d\n", peps_with_na))

      # MNAR summary
      mnar_scores <- x$missingness$mnar_score[!is.na(x$missingness$mnar_score)]
      if (length(mnar_scores) > 0 && peps_with_na > 0) {
        n_mnar <- sum(mnar_scores > 2, na.rm = TRUE)  # MNAR score > 2 suggests MNAR
        cat(sprintf("  MNAR pattern: %d of %d peptides with missing data (score > 2)\n",
                    n_mnar, peps_with_na))
      }
    }
  }

  invisible(x)
}


# --- peppwr_power class ---

#' Create a new peppwr_power object
#'
#' @param mode Either "aggregate" or "per_peptide"
#' @param question What was solved for: "power", "sample_size", or "effect_size"
#' @param answer The computed answer
#' @param simulations List of simulation details
#' @param params List of input parameters used
#' @param call Original function call
#'
#' @return A peppwr_power object
#' @export
new_peppwr_power <- function(mode, question, answer, simulations, params, call) {
  structure(
    list(
      mode = mode,
      question = question,
      answer = answer,
      simulations = simulations,
      params = params,
      call = call
    ),
    class = "peppwr_power"
  )
}

#' Validate a peppwr_power object
#'
#' @param x Object to validate
#'
#' @return The validated object (invisibly)
#' @export
validate_peppwr_power <- function(x) {
  if (!inherits(x, "peppwr_power")) {
    stop("Object must be of class 'peppwr_power'", call. = FALSE)
  }

  required <- c("mode", "question", "answer", "simulations", "params", "call")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop(
      "peppwr_power object missing required components: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  valid_modes <- c("aggregate", "per_peptide")
  if (!x$mode %in% valid_modes) {
    stop(
      "Invalid mode '", x$mode, "'. Must be one of: ",
      paste(valid_modes, collapse = ", "),
      call. = FALSE
    )
  }

  valid_questions <- c("power", "sample_size", "effect_size")
  if (!x$question %in% valid_questions) {
    stop(
      "Invalid question '", x$question, "'. Must be one of: ",
      paste(valid_questions, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(x)
}

#' Print method for peppwr_power
#'
#' @param x A peppwr_power object
#' @param ... Additional arguments (ignored)
#'
#' @return The object invisibly
#' @export
#' @method print peppwr_power
print.peppwr_power <- function(x, ...) {
  validate_peppwr_power(x)

  cat("peppwr_power analysis\n")
  cat("---------------------\n")
  cat(sprintf("Mode: %s\n\n", x$mode))

  if (x$question == "power") {
    power_pct <- round(x$answer * 100)
    cat(sprintf("Power: %d%%\n", as.integer(power_pct)))
    if (!is.null(x$params$n_per_group)) {
      cat(sprintf("Sample size: %d per group\n", x$params$n_per_group))
    }
    if (!is.null(x$params$effect_size)) {
      cat(sprintf("Effect size: %.2f-fold\n", x$params$effect_size))
    }
  } else if (x$question == "sample_size") {
    cat(sprintf("Recommended sample size: N=%d per group\n", x$answer))
    if (!is.null(x$params$target_power)) {
      cat(sprintf("Target power: %.0f%%\n", x$params$target_power * 100))
    }
    if (!is.null(x$params$effect_size)) {
      cat(sprintf("Effect size: %.2f-fold\n", x$params$effect_size))
    }
  } else if (x$question == "effect_size") {
    cat(sprintf("Minimum detectable effect: %.2f-fold\n", x$answer))
    if (!is.null(x$params$n_per_group)) {
      cat(sprintf("Sample size: %d per group\n", x$params$n_per_group))
    }
    if (!is.null(x$params$target_power)) {
      cat(sprintf("Target power: %.0f%%\n", x$params$target_power * 100))
    }
  }

  # Show test-appropriate threshold
  if (!is.null(x$params$test)) {
    cat(sprintf("\nStatistical test: %s\n", x$params$test))
    if (x$params$test == "bayes_t") {
      cat("Decision threshold: BF > 3 (substantial evidence)\n")
    } else if (!is.null(x$params$alpha)) {
      cat(sprintf("Significance level: %.2f\n", x$params$alpha))
    }
  } else if (!is.null(x$params$alpha)) {
    cat(sprintf("\nSignificance level: %.2f\n", x$params$alpha))
  }

  # FDR information
  if (isTRUE(x$params$apply_fdr) || isTRUE(x$simulations$fdr_adjusted)) {
    cat("\nFDR-adjusted analysis (Benjamini-Hochberg)\n")
    if (!is.null(x$params$prop_null)) {
      cat(sprintf("Proportion true nulls: %.0f%%\n", x$params$prop_null * 100))
    }
    if (!is.null(x$params$fdr_threshold)) {
      cat(sprintf("FDR threshold: %.0f%%\n", x$params$fdr_threshold * 100))
    }
  }

  invisible(x)
}

# --- Summary methods ---

#' Summary method for peppwr_fits
#'
#' @param object A peppwr_fits object
#' @param ... Additional arguments (ignored)
#'
#' @return A list with summary statistics
#' @export
#' @method summary peppwr_fits
summary.peppwr_fits <- function(object, ...) {
  validate_peppwr_fits(object)

  n_peptides <- length(object$best)
  n_failed <- sum(is.na(object$best))

  # Count best distributions
  best_dist_counts <- table(object$best, useNA = "no")

  # Collect AIC and loglik statistics
  all_aic <- numeric(0)
  all_loglik <- numeric(0)
  for (fit_df in object$fits) {
    valid <- !is.na(fit_df$aic)
    all_aic <- c(all_aic, fit_df$aic[valid])
    all_loglik <- c(all_loglik, fit_df$loglik[valid])
  }

  fit_summary <- list(
    aic_range = if (length(all_aic) > 0) range(all_aic) else c(NA, NA),
    aic_median = if (length(all_aic) > 0) stats::median(all_aic) else NA,
    loglik_range = if (length(all_loglik) > 0) range(all_loglik) else c(NA, NA),
    loglik_median = if (length(all_loglik) > 0) stats::median(all_loglik) else NA
  )

  result <- list(
    n_peptides = n_peptides,
    n_failed = n_failed,
    best_dist_counts = best_dist_counts,
    fit_summary = fit_summary
  )

  # Missingness summary
  if (!is.null(object$missingness) && nrow(object$missingness) > 0) {
    na_rates <- object$missingness$na_rate
    mnar_scores <- object$missingness$mnar_score

    result$missingness <- list(
      total_missing = sum(object$missingness$n_missing),
      total_values = sum(object$missingness$n_total),
      mean_na_rate = mean(na_rates, na.rm = TRUE),
      median_na_rate = stats::median(na_rates, na.rm = TRUE),
      max_na_rate = max(na_rates, na.rm = TRUE),
      n_peptides_with_na = sum(na_rates > 0),
      mean_mnar_score = mean(mnar_scores, na.rm = TRUE),
      n_potential_mnar = sum(mnar_scores > 2, na.rm = TRUE)
    )
  }

  class(result) <- "summary.peppwr_fits"
  result
}

#' Print method for summary.peppwr_fits
#'
#' @param x A summary.peppwr_fits object
#' @param ... Additional arguments (ignored)
#'
#' @return The object invisibly
#' @export
#' @method print summary.peppwr_fits
print.summary.peppwr_fits <- function(x, ...) {
  cat("Summary of peppwr_fits\n")
  cat("======================\n\n")

  cat(sprintf("Peptides fitted: %d\n", x$n_peptides))
  cat(sprintf("Failed fits: %d (%.1f%%)\n\n", x$n_failed,
              100 * x$n_failed / x$n_peptides))

  cat("Best distribution counts:\n")
  for (dist_name in names(x$best_dist_counts)) {
    cat(sprintf("  %s: %d\n", dist_name, x$best_dist_counts[dist_name]))
  }

  cat("\nFit statistics:\n")
  cat(sprintf("  AIC range: [%.1f, %.1f]\n",
              x$fit_summary$aic_range[1], x$fit_summary$aic_range[2]))
  cat(sprintf("  AIC median: %.1f\n", x$fit_summary$aic_median))
  cat(sprintf("  LogLik range: [%.1f, %.1f]\n",
              x$fit_summary$loglik_range[1], x$fit_summary$loglik_range[2]))
  cat(sprintf("  LogLik median: %.1f\n", x$fit_summary$loglik_median))

  invisible(x)
}

#' Summary method for peppwr_power
#'
#' @param object A peppwr_power object
#' @param ... Additional arguments (ignored)
#'
#' @return A list with summary statistics
#' @export
#' @method summary peppwr_power
summary.peppwr_power <- function(object, ...) {
  validate_peppwr_power(object)

  result <- list(
    mode = object$mode,
    question = object$question,
    answer = object$answer,
    n_sim = object$simulations$n_sim,
    params = object$params
  )

  # Add confidence interval for power estimate (using binomial approximation)
  if (object$question == "power" && !is.null(object$simulations$n_sim)) {
    n <- object$simulations$n_sim
    p <- object$answer
    se <- sqrt(p * (1 - p) / n)
    result$ci_lower <- max(0, p - 1.96 * se)
    result$ci_upper <- min(1, p + 1.96 * se)
  }

  # Add per-peptide summary if available
  if (object$mode == "per_peptide") {
    if ("peptide_power" %in% names(object$simulations)) {
      pep_power <- object$simulations$peptide_power
      pep_power <- pep_power[!is.na(pep_power)]
      result$peptide_summary <- list(
        n_analyzed = length(pep_power),
        mean_power = mean(pep_power),
        median_power = stats::median(pep_power),
        min_power = min(pep_power),
        max_power = max(pep_power)
      )
    }
    if ("n_analyzed" %in% names(object$simulations)) {
      result$n_analyzed <- object$simulations$n_analyzed
    }
    if ("n_excluded" %in% names(object$simulations)) {
      result$n_excluded <- object$simulations$n_excluded
    }
  }

  class(result) <- "summary.peppwr_power"
  result
}

#' Print method for summary.peppwr_power
#'
#' @param x A summary.peppwr_power object
#' @param ... Additional arguments (ignored)
#'
#' @return The object invisibly
#' @export
#' @method print summary.peppwr_power
print.summary.peppwr_power <- function(x, ...) {
  cat("Summary of peppwr_power\n")
  cat("=======================\n\n")

  cat(sprintf("Mode: %s\n", x$mode))
  cat(sprintf("Question: %s\n", x$question))
  cat(sprintf("Answer: %s\n\n", format_answer(x$answer, x$question)))

  cat("Simulation parameters:\n")
  cat(sprintf("  Number of simulations: %d\n", x$n_sim))
  if (!is.null(x$params$effect_size)) {
    cat(sprintf("  Effect size: %.2f-fold\n", x$params$effect_size))
  }
  if (!is.null(x$params$n_per_group)) {
    cat(sprintf("  Sample size: %d per group\n", x$params$n_per_group))
  }
  if (!is.null(x$params$test)) {
    cat(sprintf("  Test: %s\n", x$params$test))
    if (x$params$test == "bayes_t") {
      cat("  Decision threshold: BF > 3\n")
    } else if (!is.null(x$params$alpha)) {
      cat(sprintf("  Alpha: %.2f\n", x$params$alpha))
    }
  } else if (!is.null(x$params$alpha)) {
    cat(sprintf("  Alpha: %.2f\n", x$params$alpha))
  }

  # Confidence interval
  if (!is.null(x$ci_lower) && !is.null(x$ci_upper)) {
    cat(sprintf("\n95%% CI for power: [%.1f%%, %.1f%%]\n",
                x$ci_lower * 100, x$ci_upper * 100))
  }

  # Per-peptide summary
  if (!is.null(x$peptide_summary)) {
    cat("\nPer-peptide summary:\n")
    cat(sprintf("  Peptides analyzed: %d\n", x$peptide_summary$n_analyzed))
    cat(sprintf("  Mean power: %.1f%%\n", x$peptide_summary$mean_power * 100))
    cat(sprintf("  Median power: %.1f%%\n", x$peptide_summary$median_power * 100))
    cat(sprintf("  Range: [%.1f%%, %.1f%%]\n",
                x$peptide_summary$min_power * 100, x$peptide_summary$max_power * 100))
  }

  invisible(x)
}

#' Format answer based on question type
#'
#' @param answer The answer value
#' @param question The question type
#' @return Formatted string
#' @keywords internal
format_answer <- function(answer, question) {
  switch(question,
    "power" = sprintf("%.1f%%", answer * 100),
    "sample_size" = sprintf("N=%d per group", answer),
    "effect_size" = sprintf("%.2f-fold", answer),
    as.character(answer)
  )
}
