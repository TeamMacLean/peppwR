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
