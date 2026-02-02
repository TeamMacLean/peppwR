# Power analysis functions for peppwR

#' Power analysis for peptide abundance data
#'
#' @param distribution Distribution name (character) or peppwr_fits object for per-peptide mode
#' @param ... Additional arguments
#'
#' @return A peppwr_power object
#' @export
power_analysis <- function(distribution, ...) {
  UseMethod("power_analysis")
}

#' Power analysis with specified distribution (aggregate mode)
#'
#' @param distribution Distribution name (e.g., "norm", "gamma", "lnorm")
#' @param params List of distribution parameters
#' @param effect_size Fold change to detect
#' @param n_per_group Sample size per group (required for find="power")
#' @param target_power Target power (required for find="sample_size" or find="effect_size")
#' @param alpha Significance level (default 0.05)
#' @param test Statistical test to use (default "wilcoxon")
#' @param find What to solve for: "power", "sample_size", or "effect_size"
#' @param n_sim Number of simulations (default 1000)
#' @param ... Additional arguments (ignored)
#'
#' @return A peppwr_power object
#' @export
#' @method power_analysis character
power_analysis.character <- function(distribution, params, effect_size = NULL, n_per_group = NULL,
                                     target_power = NULL, alpha = 0.05, test = "wilcoxon",
                                     find = "power", n_sim = 1000, ...) {
  the_call <- match.call()

  # Validate find parameter
  valid_finds <- c("power", "sample_size", "effect_size")
  if (!find %in% valid_finds) {
    stop("find must be one of: ", paste(valid_finds, collapse = ", "))
  }

  # Validate required parameters based on find
  if (find == "power") {
    if (is.null(n_per_group)) {
      stop("n_per_group is required when find='power'")
    }
    if (is.null(effect_size)) {
      stop("effect_size is required when find='power'")
    }
  } else if (find == "sample_size") {
    if (is.null(target_power)) {
      stop("target_power is required when find='sample_size'")
    }
    if (is.null(effect_size)) {
      stop("effect_size is required when find='sample_size'")
    }
  } else if (find == "effect_size") {
    if (is.null(n_per_group)) {
      stop("n_per_group is required when find='effect_size'")
    }
    if (is.null(target_power)) {
      stop("target_power is required when find='effect_size'")
    }
  }

  # Perform the analysis based on what we're solving for
  if (find == "power") {
    answer <- run_power_sim(distribution, params, n_per_group, effect_size,
                            alpha, test, n_sim)
    simulations <- list(n_sim = n_sim, power = answer)
  } else if (find == "sample_size") {
    result <- find_sample_size(distribution, params, effect_size, target_power,
                               alpha, test, n_sim)
    answer <- result$n
    n_per_group <- result$n
    simulations <- list(n_sim = n_sim, power_curve = result$power_curve)
  } else if (find == "effect_size") {
    result <- find_effect_size(distribution, params, n_per_group, target_power,
                               alpha, test, n_sim)
    answer <- result$effect_size
    effect_size <- result$effect_size
    simulations <- list(n_sim = n_sim, effect_curve = result$effect_curve)
  }

  # Store parameters
  stored_params <- list(
    distribution = distribution,
    params = params,
    effect_size = effect_size,
    n_per_group = n_per_group,
    target_power = target_power,
    alpha = alpha,
    test = test
  )

  new_peppwr_power(
    mode = "aggregate",
    question = find,
    answer = answer,
    simulations = simulations,
    params = stored_params,
    call = the_call
  )
}

#' Default power analysis method (aggregate mode with defaults)
#'
#' @rdname power_analysis.character
#' @export
#' @method power_analysis default
power_analysis.default <- function(distribution, ...) {
  # If distribution is a string, use the character method
  if (is.character(distribution)) {
    power_analysis.character(distribution, ...)
  } else {
    stop("Unknown input type for power_analysis")
  }
}

#' Power analysis for per-peptide mode using fitted distributions
#'
#' @param distribution A peppwr_fits object from fit_distributions()
#' @param effect_size Fold change to detect
#' @param n_per_group Sample size per group (required for find="power")
#' @param target_power Target power (required for find="sample_size")
#' @param alpha Significance level (default 0.05)
#' @param test Statistical test to use (default "wilcoxon")
#' @param find What to solve for: "power" or "sample_size"
#' @param n_sim Number of simulations per peptide (default 1000)
#' @param on_fit_failure How to handle failed fits: "exclude", "empirical", or "lognormal"
#' @param proportion_threshold Proportion of peptides that must reach target_power (default 0.5)
#' @param include_missingness If TRUE, incorporate peptide-specific NA rates into simulations
#' @param apply_fdr If TRUE, use FDR-aware simulation with Benjamini-Hochberg correction
#' @param prop_null Proportion of true null peptides (default 0.9 = 90% unchanged)
#' @param fdr_threshold FDR threshold for calling discoveries (default 0.05)
#' @param ... Additional arguments (ignored)
#'
#' @return A peppwr_power object
#' @export
#' @method power_analysis peppwr_fits
power_analysis.peppwr_fits <- function(distribution, effect_size = NULL, n_per_group = NULL,
                                       target_power = NULL, alpha = 0.05, test = "wilcoxon",
                                       find = "power", n_sim = 1000,
                                       on_fit_failure = "exclude",
                                       proportion_threshold = 0.5,
                                       include_missingness = FALSE,
                                       apply_fdr = FALSE, prop_null = 0.9,
                                       fdr_threshold = 0.05, ...) {
  the_call <- match.call()
  fits <- distribution

  # Extract missingness info if needed
  missingness_data <- if (include_missingness && !is.null(fits$missingness)) {
    fits$missingness
  } else {
    NULL
  }

  # FDR mode - use whole-peptidome simulation
  if (apply_fdr && find == "power") {
    if (is.null(n_per_group)) {
      stop("n_per_group is required when find='power'")
    }
    if (is.null(effect_size)) {
      stop("effect_size is required when find='power'")
    }

    fdr_power <- run_power_sim_fdr(
      fits, effect_size, n_per_group,
      prop_null = prop_null, fdr_threshold = fdr_threshold,
      alpha = alpha, test = test, n_sim = n_sim
    )

    return(new_peppwr_power(
      mode = "per_peptide",
      question = "power",
      answer = fdr_power,
      simulations = list(
        n_sim = n_sim,
        fdr_adjusted = TRUE,
        prop_null = prop_null,
        fdr_threshold = fdr_threshold
      ),
      params = list(
        effect_size = effect_size,
        n_per_group = n_per_group,
        alpha = alpha,
        test = test,
        apply_fdr = TRUE,
        prop_null = prop_null,
        fdr_threshold = fdr_threshold
      ),
      call = the_call
    ))
  }

  # Validate on_fit_failure
  valid_failures <- c("exclude", "empirical", "lognormal")
  if (!on_fit_failure %in% valid_failures) {
    stop("on_fit_failure must be one of: ", paste(valid_failures, collapse = ", "))
  }

  # Validate find parameter
  valid_finds <- c("power", "sample_size", "effect_size")
  if (!find %in% valid_finds) {
    stop("find must be one of: ", paste(valid_finds, collapse = ", "))
  }

  # Validate required parameters
  if (find == "power") {
    if (is.null(n_per_group)) {
      stop("n_per_group is required when find='power'")
    }
    if (is.null(effect_size)) {
      stop("effect_size is required when find='power'")
    }
  } else if (find == "sample_size") {
    if (is.null(target_power)) {
      stop("target_power is required when find='sample_size'")
    }
    if (is.null(effect_size)) {
      stop("effect_size is required when find='sample_size'")
    }
  } else if (find == "effect_size") {
    if (is.null(n_per_group)) {
      stop("n_per_group is required when find='effect_size'")
    }
    if (is.null(target_power)) {
      stop("target_power is required when find='effect_size'")
    }
  }

  # Get peptide-specific power estimates
  n_peptides <- length(fits$best)
  peptide_power <- numeric(n_peptides)
  n_excluded <- 0
  n_analyzed <- 0

  for (i in seq_len(n_peptides)) {
    best_dist <- fits$best[i]
    fit_df <- fits$fits[[i]]
    raw_data <- fits$data$data[[i]][[1]]

    # Handle fit failure
    if (is.na(best_dist)) {
      if (on_fit_failure == "exclude") {
        peptide_power[i] <- NA
        n_excluded <- n_excluded + 1
        next
      } else if (on_fit_failure == "lognormal") {
        # Use lognormal with moment-matched parameters
        best_dist <- "Lognormal"
        raw_data <- raw_data[raw_data > 0]  # Lognormal needs positive values
        if (length(raw_data) < 2) {
          peptide_power[i] <- NA
          n_excluded <- n_excluded + 1
          next
        }
        m <- mean(raw_data)
        s <- stats::sd(raw_data)
        meanlog <- log(m^2 / sqrt(s^2 + m^2))
        sdlog <- sqrt(log(1 + s^2 / m^2))
        params <- list(meanlog = meanlog, sdlog = sdlog)
      } else if (on_fit_failure == "empirical") {
        # Bootstrap from empirical data
        if (length(raw_data) >= 2) {
          if (find == "power") {
            peptide_power[i] <- run_power_sim_empirical(
              raw_data, n_per_group, effect_size, alpha, test, n_sim
            )
          }
          n_analyzed <- n_analyzed + 1
        } else {
          peptide_power[i] <- NA
          n_excluded <- n_excluded + 1
        }
        next
      }
    } else {
      # Get parameters from the best-fitting distribution
      params <- get_dist_params(best_dist, fit_df, raw_data)
    }

    # Map distribution name to R function name
    dist_rfunc <- map_dist_name(best_dist)

    if (is.null(dist_rfunc) || is.null(params)) {
      peptide_power[i] <- NA
      n_excluded <- n_excluded + 1
      next
    }

    n_analyzed <- n_analyzed + 1

    # Run power simulation for this peptide
    if (find == "power") {
      if (!is.null(missingness_data) && missingness_data$na_rate[i] > 0) {
        # Use missingness-aware simulation
        peptide_power[i] <- run_power_sim_with_missingness(
          dist_rfunc, params, n_per_group, effect_size,
          na_rate = missingness_data$na_rate[i],
          mnar_score = missingness_data$mnar_score[i] %||% 0,
          alpha = alpha, test = test, n_sim = n_sim
        )
      } else {
        peptide_power[i] <- run_power_sim(dist_rfunc, params, n_per_group, effect_size,
                                          alpha, test, n_sim)
      }
    }
  }

  # Helper for null coalescing
  `%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

  # For find = "sample_size", search over N values
  if (find == "sample_size") {
    n_values <- c(3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30)
    proportion_powered <- sapply(n_values, function(n) {
      peps <- sapply(seq_len(n_peptides), function(i) {
        best_dist <- fits$best[i]
        fit_df <- fits$fits[[i]]
        raw_data <- fits$data$data[[i]][[1]]

        if (is.na(best_dist)) {
          if (on_fit_failure == "lognormal") {
            raw_data <- raw_data[raw_data > 0]
            if (length(raw_data) < 2) return(NA)
            m <- mean(raw_data)
            s <- stats::sd(raw_data)
            meanlog <- log(m^2 / sqrt(s^2 + m^2))
            sdlog <- sqrt(log(1 + s^2 / m^2))
            params <- list(meanlog = meanlog, sdlog = sdlog)
            dist_rfunc <- "lnorm"
          } else {
            return(NA)
          }
        } else {
          params <- get_dist_params(best_dist, fit_df, raw_data)
          dist_rfunc <- map_dist_name(best_dist)
          if (is.null(dist_rfunc) || is.null(params)) return(NA)
        }
        run_power_sim(dist_rfunc, params, n, effect_size, alpha, test, n_sim)
      })
      mean(peps >= target_power, na.rm = TRUE)
    })

    power_curve <- data.frame(n_per_group = n_values, proportion_powered = proportion_powered)

    # Find N where proportion of peptides reaches threshold
    above_target <- proportion_powered >= proportion_threshold
    if (any(above_target)) {
      answer <- n_values[which(above_target)[1]]
    } else {
      answer <- max(n_values)
    }

    simulations <- list(
      n_sim = n_sim,
      power_curve = power_curve,
      proportion_powered = proportion_powered,
      n_analyzed = n_analyzed,
      n_excluded = n_excluded
    )
  } else if (find == "effect_size") {
    effect_values <- c(1.5, 2, 3, 5, 10)
    proportion_powered <- sapply(effect_values, function(es) {
      peps <- sapply(seq_len(n_peptides), function(i) {
        best_dist <- fits$best[i]
        fit_df <- fits$fits[[i]]
        raw_data <- fits$data$data[[i]][[1]]

        if (is.na(best_dist)) {
          if (on_fit_failure == "lognormal") {
            raw_data <- raw_data[raw_data > 0]
            if (length(raw_data) < 2) return(NA)
            m <- mean(raw_data)
            s <- stats::sd(raw_data)
            meanlog <- log(m^2 / sqrt(s^2 + m^2))
            sdlog <- sqrt(log(1 + s^2 / m^2))
            params <- list(meanlog = meanlog, sdlog = sdlog)
            dist_rfunc <- "lnorm"
          } else {
            return(NA)
          }
        } else {
          params <- get_dist_params(best_dist, fit_df, raw_data)
          dist_rfunc <- map_dist_name(best_dist)
          if (is.null(dist_rfunc) || is.null(params)) return(NA)
        }
        run_power_sim(dist_rfunc, params, n_per_group, es, alpha, test, n_sim)
      })
      mean(peps >= target_power, na.rm = TRUE)
    })

    effect_curve <- data.frame(effect_size = effect_values, proportion_powered = proportion_powered)

    # Find minimum effect where proportion reaches threshold
    above_target <- proportion_powered >= proportion_threshold
    if (any(above_target)) {
      answer <- effect_values[which(above_target)[1]]
    } else {
      answer <- max(effect_values)
    }

    simulations <- list(
      n_sim = n_sim,
      effect_curve = effect_curve,
      proportion_powered = proportion_powered,
      n_analyzed = sum(!is.na(fits$best)),
      n_excluded = sum(is.na(fits$best))
    )
  } else {
    # find = "power"
    answer <- mean(peptide_power, na.rm = TRUE)
    simulations <- list(
      n_sim = n_sim,
      peptide_power = peptide_power,
      n_analyzed = n_analyzed,
      n_excluded = n_excluded
    )
  }

  # Store parameters
  stored_params <- list(
    effect_size = effect_size,
    n_per_group = n_per_group,
    target_power = target_power,
    alpha = alpha,
    test = test,
    proportion_threshold = proportion_threshold
  )

  new_peppwr_power(
    mode = "per_peptide",
    question = find,
    answer = answer,
    simulations = simulations,
    params = stored_params,
    call = the_call
  )
}

#' Get distribution parameters from fitted results
#'
#' @param dist_name Distribution name
#' @param fit_df Fit results data frame
#' @param raw_data Raw data values
#'
#' @return List of parameters
#' @keywords internal
get_dist_params <- function(dist_name, fit_df, raw_data) {
  # Use moment-matched parameters based on distribution type
  raw_data <- raw_data[!is.na(raw_data)]
  if (length(raw_data) < 2) return(NULL)

  m <- mean(raw_data)
  s <- stats::sd(raw_data)
  if (s <= 0) s <- 0.001 * abs(m)

  switch(dist_name,
    "Normal" = list(mean = m, sd = s),
    "Gamma" = {
      if (m <= 0) m <- 0.001
      shape <- (m / s)^2
      rate <- m / s^2
      list(shape = shape, rate = rate)
    },
    "Lognormal" = {
      raw_data <- raw_data[raw_data > 0]
      if (length(raw_data) < 2) return(NULL)
      m <- mean(raw_data)
      s <- stats::sd(raw_data)
      meanlog <- log(m^2 / sqrt(s^2 + m^2))
      sdlog <- sqrt(log(1 + s^2 / m^2))
      list(meanlog = meanlog, sdlog = sdlog)
    },
    "Skew Normal" = list(mean = m, sd = s),  # Approximate with normal
    "InvGamma" = {
      if (m <= 0) m <- 0.001
      shape <- 2 + m^2 / s^2
      rate <- m * (shape - 1)
      list(shape = shape, rate = rate)
    },
    "Inverse Gaussian" = {
      if (m <= 0) m <- 0.001
      list(mean = m, shape = m^3 / s^2)
    },
    "Log Gamma" = {
      raw_data <- raw_data[raw_data > 0]
      if (length(raw_data) < 2) return(NULL)
      m <- mean(raw_data)
      s <- stats::sd(raw_data)
      meanlog <- log(m^2 / sqrt(s^2 + m^2))
      sdlog <- sqrt(log(1 + s^2 / m^2))
      list(meanlog = meanlog, sdlog = sdlog)  # Approximate with lognormal
    },
    "Pareto" = {
      raw_data <- raw_data[raw_data > 0]
      if (length(raw_data) < 2) return(NULL)
      xm <- min(raw_data)
      alpha_par <- length(raw_data) / sum(log(raw_data / xm))
      list(xmin = xm, alpha = alpha_par)
    },
    "Negative Binomial" = {
      if (m <= 0) m <- 1
      size <- m^2 / max(s^2 - m, 0.001)
      list(size = size, mu = m)
    },
    NULL
  )
}

#' Map distribution name to R function prefix
#'
#' @param dist_name Distribution name
#' @return R function name for random generation
#' @keywords internal
map_dist_name <- function(dist_name) {
  switch(dist_name,
    "Normal" = "norm",
    "Gamma" = "gamma",
    "Lognormal" = "lnorm",
    "Skew Normal" = "norm",  # Use normal as approximation
    "InvGamma" = "invgamma",
    "Inverse Gaussian" = "invgauss",
    "Log Gamma" = "lnorm",  # Use lognormal as approximation
    "Pareto" = NULL,  # Not commonly supported in base R
    "Negative Binomial" = NULL,  # Discrete distribution
    NULL
  )
}

#' Find required sample size for target power
#'
#' @param distribution Distribution name
#' @param params Distribution parameters
#' @param effect_size Effect size to detect
#' @param target_power Target power level
#' @param alpha Significance level
#' @param test Statistical test
#' @param n_sim Number of simulations per sample size
#'
#' @return List with n and power_curve
#' @keywords internal
find_sample_size <- function(distribution, params, effect_size, target_power,
                             alpha, test, n_sim) {
  # Search over sample sizes
  n_values <- c(3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50)

  power_values <- sapply(n_values, function(n) {
    run_power_sim(distribution, params, n, effect_size, alpha, test, n_sim)
  })

  power_curve <- data.frame(n_per_group = n_values, power = power_values)

  # Find the smallest N that achieves target power
  above_target <- power_values >= target_power
  if (any(above_target)) {
    n <- n_values[which(above_target)[1]]
  } else {
    # If we never reach target, return the largest N tested
    n <- max(n_values)
  }

  list(n = n, power_curve = power_curve)
}

#' Find minimum detectable effect size
#'
#' @param distribution Distribution name
#' @param params Distribution parameters
#' @param n_per_group Sample size per group
#' @param target_power Target power level
#' @param alpha Significance level
#' @param test Statistical test
#' @param n_sim Number of simulations per effect size
#'
#' @return List with effect_size and effect_curve
#' @keywords internal
find_effect_size <- function(distribution, params, n_per_group, target_power,
                             alpha, test, n_sim) {
  # Search over effect sizes
  effect_values <- c(1.1, 1.2, 1.3, 1.5, 1.75, 2, 2.5, 3, 4, 5)

  power_values <- sapply(effect_values, function(es) {
    run_power_sim(distribution, params, n_per_group, es, alpha, test, n_sim)
  })

  effect_curve <- data.frame(effect_size = effect_values, power = power_values)

  # Find the smallest effect size that achieves target power
  above_target <- power_values >= target_power
  if (any(above_target)) {
    es <- effect_values[which(above_target)[1]]
  } else {
    # If we never reach target, return the largest effect tested
    es <- max(effect_values)
  }

  list(effect_size = es, effect_curve = effect_curve)
}
