# Simulation engine for peppwR

#' Simulate an experiment with control and treatment groups
#'
#' @param distribution Distribution name (e.g., "norm", "gamma", "lnorm")
#' @param params List of distribution parameters
#' @param n_per_group Number of samples per group
#' @param effect_size Fold change for treatment (multiplicative effect)
#'
#' @return List with control and treatment vectors
#' @export
simulate_experiment <- function(distribution, params, n_per_group, effect_size) {
  # Generate random samples from the specified distribution
  rfunc <- get_rfunc(distribution)

  # Generate control samples
  control <- do.call(rfunc, c(list(n = n_per_group), params))

  # Generate treatment samples (apply multiplicative effect)
  treatment <- do.call(rfunc, c(list(n = n_per_group), params)) * effect_size

  list(control = control, treatment = treatment)
}

#' Get the random generation function for a distribution
#'
#' @param distribution Distribution name
#' @return Random generation function
#' @keywords internal
get_rfunc <- function(distribution) {
  switch(distribution,
    "norm" = stats::rnorm,
    "gamma" = stats::rgamma,
    "lnorm" = stats::rlnorm,
    "invgamma" = function(n, shape, rate) 1 / stats::rgamma(n, shape = shape, rate = rate),
    "invgauss" = function(n, mean, shape) {
      # Inverse Gaussian sampling using Chhikara & Folks method
      mu <- mean
      lambda <- shape
      y <- stats::rnorm(n)^2
      x <- mu + (mu^2 * y) / (2 * lambda) - (mu / (2 * lambda)) * sqrt(4 * mu * lambda * y + mu^2 * y^2)
      u <- stats::runif(n)
      ifelse(u <= mu / (mu + x), x, mu^2 / x)
    },
    stop("Unknown distribution: ", distribution)
  )
}

#' Wilcoxon rank-sum test
#'
#' @param control Control group values
#' @param treatment Treatment group values
#'
#' @return p-value from the test
#' @export
test_wilcoxon <- function(control, treatment) {
  result <- stats::wilcox.test(control, treatment)
  result$p.value
}

#' Run power simulation
#'
#' @param distribution Distribution name
#' @param params List of distribution parameters
#' @param n_per_group Number of samples per group
#' @param effect_size Fold change for treatment
#' @param alpha Significance level
#' @param test Statistical test to use ("wilcoxon", "bootstrap_t", "bayes_t", "rankprod")
#' @param n_sim Number of simulations
#'
#' @return Power estimate (proportion of significant results)
#' @export
run_power_sim <- function(distribution, params, n_per_group, effect_size,
                          alpha = 0.05, test = "wilcoxon", n_sim = 1000) {
  # Get the test function
  test_func <- get_test_func(test)

  # Run simulations
  results <- replicate(n_sim, {
    samples <- simulate_experiment(distribution, params, n_per_group, effect_size)
    p_or_bf <- test_func(samples$control, samples$treatment)

    # For Bayes factor tests, check if BF > threshold (use 3 as default)
    if (test == "bayes_t") {
      return(p_or_bf > 3)  # BF > 3 indicates substantial evidence
    } else {
      return(p_or_bf < alpha)
    }
  })

  mean(results)
}

#' Get the test function for a given test name
#'
#' @param test Test name
#' @return Test function
#' @keywords internal
get_test_func <- function(test) {
  switch(test,
    "wilcoxon" = test_wilcoxon,
    "bootstrap_t" = test_bootstrap_t,
    "bayes_t" = test_bayes_t,
    "rankprod" = test_rankprod,
    stop("Unknown test: ", test)
  )
}

# ============================================================================
# Phase B: Empirical Bootstrap
# ============================================================================

#' Simulate experiment using empirical bootstrap
#'
#' Resamples from observed data instead of using parametric distributions.
#' Useful when distribution fitting fails or for non-parametric analysis.
#'
#' @param raw_data Numeric vector of observed values
#' @param n_per_group Number of samples per group
#' @param effect_size Fold change for treatment (multiplicative effect)
#'
#' @return List with control and treatment vectors
#' @export
simulate_empirical <- function(raw_data, n_per_group, effect_size) {
  # Remove NA values
  raw_data <- raw_data[!is.na(raw_data)]

  if (length(raw_data) < 2) {
    stop("raw_data must have at least 2 non-NA values")
  }

  # Bootstrap resample for control group
  control <- sample(raw_data, size = n_per_group, replace = TRUE)

  # Bootstrap resample for treatment group, then apply effect
  treatment <- sample(raw_data, size = n_per_group, replace = TRUE) * effect_size

  list(control = control, treatment = treatment)
}

#' Run power simulation using empirical bootstrap
#'
#' Estimates power by repeatedly resampling from observed data
#' rather than sampling from fitted distributions.
#'
#' @param raw_data Numeric vector of observed values
#' @param n_per_group Number of samples per group
#' @param effect_size Fold change for treatment
#' @param alpha Significance level
#' @param test Statistical test to use
#' @param n_sim Number of simulations
#'
#' @return Power estimate (proportion of significant results)
#' @export
run_power_sim_empirical <- function(raw_data, n_per_group, effect_size,
                                     alpha = 0.05, test = "wilcoxon",
                                     n_sim = 1000) {
  # Get the test function
  test_func <- get_test_func(test)

  # Run simulations
  results <- replicate(n_sim, {
    samples <- simulate_empirical(raw_data, n_per_group, effect_size)
    p_or_bf <- test_func(samples$control, samples$treatment)

    # For Bayes factor tests, check if BF > threshold
    if (test == "bayes_t") {
      return(p_or_bf > 3)
    } else {
      return(p_or_bf < alpha)
    }
  })

  mean(results)
}

# ============================================================================
# Phase C: Simulation with Missingness
# ============================================================================

#' Simulate experiment with realistic missingness
#'
#' Generates control and treatment samples from a distribution,
#' then introduces missing values according to the specified rate and
#' MNAR pattern.
#'
#' @param distribution Distribution name (e.g., "norm", "gamma", "lnorm")
#' @param params List of distribution parameters
#' @param n_per_group Number of samples per group
#' @param effect_size Fold change for treatment (multiplicative effect)
#' @param na_rate Proportion of values to make NA (0-1)
#' @param mnar_score MNAR intensity: 0 = MCAR, positive = low values more
#'   likely to be missing. Typical values: 0-3.
#'
#' @return List with control and treatment vectors (may contain NAs)
#' @export
simulate_with_missingness <- function(distribution, params, n_per_group,
                                       effect_size, na_rate = 0,
                                       mnar_score = 0) {
  # First generate complete data
  samples <- simulate_experiment(distribution, params, n_per_group, effect_size)

  if (na_rate <= 0) {
    return(samples)
  }

  # Apply missingness
  samples$control <- apply_missingness(samples$control, na_rate, mnar_score)
  samples$treatment <- apply_missingness(samples$treatment, na_rate, mnar_score)

  samples
}

#' Apply missingness to a vector
#'
#' @param x Numeric vector
#' @param na_rate Proportion to make NA
#' @param mnar_score MNAR intensity (0 = MCAR)
#'
#' @return Vector with some values replaced by NA
#' @keywords internal
apply_missingness <- function(x, na_rate, mnar_score = 0) {
  n <- length(x)
  n_missing <- round(n * na_rate)

  if (n_missing == 0) return(x)
  if (n_missing >= n) return(rep(NA_real_, n))

  if (abs(mnar_score) < 0.01) {
    # MCAR: random selection
    missing_idx <- sample(n, n_missing)
  } else {
    # MNAR: low values more likely to be missing when mnar_score > 0
    # Use ranks to determine missingness probability
    ranks <- rank(x)
    # Transform ranks to probabilities using logistic function
    # Higher mnar_score = stronger bias toward missing low values
    probs <- 1 / (1 + exp(mnar_score * (ranks - stats::median(ranks)) / (n / 4)))
    probs <- probs / sum(probs)  # Normalize

    missing_idx <- sample(n, n_missing, prob = probs, replace = FALSE)
  }

  x[missing_idx] <- NA
  x
}

#' Run power simulation with missingness
#'
#' Estimates power accounting for realistic missing data patterns.
#'
#' @param distribution Distribution name
#' @param params List of distribution parameters
#' @param n_per_group Number of samples per group
#' @param effect_size Fold change for treatment
#' @param na_rate Proportion of values that will be NA
#' @param mnar_score MNAR intensity (0 = MCAR)
#' @param alpha Significance level
#' @param test Statistical test to use
#' @param n_sim Number of simulations
#'
#' @return Power estimate (proportion of significant results)
#' @export
run_power_sim_with_missingness <- function(distribution, params, n_per_group,
                                            effect_size, na_rate = 0,
                                            mnar_score = 0, alpha = 0.05,
                                            test = "wilcoxon", n_sim = 1000) {
  test_func <- get_test_func(test)

  results <- replicate(n_sim, {
    samples <- simulate_with_missingness(
      distribution, params, n_per_group, effect_size, na_rate, mnar_score
    )

    # Remove NAs for testing
    control <- samples$control[!is.na(samples$control)]
    treatment <- samples$treatment[!is.na(samples$treatment)]

    # Need at least 2 values in each group
    if (length(control) < 2 || length(treatment) < 2) {
      return(FALSE)
    }

    p_or_bf <- test_func(control, treatment)

    if (test == "bayes_t") {
      return(p_or_bf > 3)
    } else {
      return(p_or_bf < alpha)
    }
  })

  mean(results)
}

# ============================================================================
# Phase D: FDR-Aware Power Simulation
# ============================================================================

#' Run FDR-aware power simulation for whole peptidome
#'
#' Simulates an entire peptidome experiment with a mix of true nulls and
#' true alternatives, then applies Benjamini-Hochberg correction to compute
#' FDR-adjusted power.
#'
#' @param fits A peppwr_fits object
#' @param effect_size Fold change for treatment
#' @param n_per_group Number of samples per group
#' @param prop_null Proportion of peptides that are true nulls (no effect).
#'   Default 0.9 (90% unchanged).
#' @param fdr_threshold FDR threshold for calling discoveries. Default 0.05.
#' @param alpha Nominal significance level (used for simulation). Default 0.05.
#' @param test Statistical test to use. Default "wilcoxon".
#' @param n_sim Number of simulation iterations. Default 1000.
#'
#' @return FDR-adjusted power estimate (proportion of true alternatives detected)
#' @export
run_power_sim_fdr <- function(fits, effect_size, n_per_group,
                               prop_null = 0.9, fdr_threshold = 0.05,
                               alpha = 0.05, test = "wilcoxon",
                               n_sim = 1000) {
  # Get peptide-specific distribution info
  n_peptides <- length(fits$best)
  test_func <- get_test_func(test)

  # Prepare distribution parameters for each peptide
  peptide_info <- lapply(seq_len(n_peptides), function(i) {
    best_dist <- fits$best[i]
    fit_df <- fits$fits[[i]]
    raw_data <- fits$data$data[[i]][[1]]

    if (is.na(best_dist)) {
      return(NULL)  # Skip peptides with failed fits
    }

    params <- get_dist_params(best_dist, fit_df, raw_data)
    dist_name <- map_dist_name(best_dist)

    if (is.null(dist_name) || is.null(params)) {
      return(NULL)
    }

    list(dist = dist_name, params = params)
  })

  # Filter out NULL entries (failed fits)
  valid_idx <- which(!sapply(peptide_info, is.null))
  peptide_info <- peptide_info[valid_idx]
  n_valid <- length(peptide_info)

  if (n_valid < 2) {
    warning("Need at least 2 peptides with valid fits for FDR simulation")
    return(NA_real_)
  }

  # Run simulations
  powers <- replicate(n_sim, {
    # Assign peptides to null (no effect) or alternative (has effect)
    is_null <- stats::runif(n_valid) < prop_null
    n_alt <- sum(!is_null)

    if (n_alt == 0) {
      # No true alternatives in this simulation
      return(0)
    }

    # Simulate each peptide and collect p-values
    p_values <- sapply(seq_len(n_valid), function(i) {
      info <- peptide_info[[i]]

      # Effect size: 1 (no effect) for nulls, actual effect for alternatives
      eff <- if (is_null[i]) 1 else effect_size

      # Simulate experiment
      samples <- tryCatch({
        simulate_experiment(info$dist, info$params, n_per_group, eff)
      }, error = function(e) NULL)

      if (is.null(samples)) return(NA)

      # Get p-value (bayes_t is blocked at power_analysis level for FDR mode)
      tryCatch({
        test_func(samples$control, samples$treatment)
      }, error = function(e) NA)
    })

    # Apply BH correction
    p_adj <- stats::p.adjust(p_values, method = "BH")

    # Count discoveries among true alternatives
    discoveries <- sum(p_adj[!is_null] < fdr_threshold, na.rm = TRUE)

    # Power = proportion of true alternatives discovered
    discoveries / n_alt
  })

  mean(powers, na.rm = TRUE)
}
