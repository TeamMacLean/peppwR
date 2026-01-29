# Statistical tests for peppwR

#' Bootstrap-t test
#'
#' Performs a bootstrap-t test comparing two groups
#'
#' @param control Control group values
#' @param treatment Treatment group values
#' @param n_boot Number of bootstrap iterations (default 1000)
#'
#' @return p-value from the test
#' @export
test_bootstrap_t <- function(control, treatment, n_boot = 1000) {
  # Observed t-statistic
  obs_t <- compute_t_stat(control, treatment)

  # Pool the data under the null hypothesis
  pooled <- c(control, treatment)
  n_control <- length(control)
  n_treatment <- length(treatment)

  # Bootstrap distribution of t under the null
  boot_t <- replicate(n_boot, {
    # Resample from pooled data
    boot_control <- sample(pooled, n_control, replace = TRUE)
    boot_treatment <- sample(pooled, n_treatment, replace = TRUE)
    compute_t_stat(boot_control, boot_treatment)
  })

  # Two-sided p-value
  p_value <- mean(abs(boot_t) >= abs(obs_t))
  max(p_value, 1 / n_boot)  # Ensure non-zero p-value
}

#' Compute t-statistic for two groups
#'
#' @param x First group
#' @param y Second group
#' @return t-statistic
#' @keywords internal
compute_t_stat <- function(x, y) {
  n_x <- length(x)
  n_y <- length(y)

  mean_diff <- mean(x) - mean(y)
  var_x <- stats::var(x)
  var_y <- stats::var(y)

  # Pooled standard error
  se <- sqrt(var_x / n_x + var_y / n_y)

  if (se == 0) return(0)
  mean_diff / se
}

#' Bayes factor t-test
#'
#' Computes Bayes factor for difference between two groups
#'
#' @param control Control group values
#' @param treatment Treatment group values
#'
#' @return Bayes factor (BF10: evidence for alternative over null)
#' @export
test_bayes_t <- function(control, treatment) {
  # Use JZS Bayes factor approximation
  # Based on Rouder et al. (2009)

  n_x <- length(control)
  n_y <- length(treatment)
  n <- n_x + n_y

  # Effect size (Cohen's d)
  pooled_sd <- sqrt(((n_x - 1) * stats::var(control) + (n_y - 1) * stats::var(treatment)) / (n - 2))
  if (pooled_sd == 0) pooled_sd <- 0.001

  d <- (mean(treatment) - mean(control)) / pooled_sd

  # t-statistic
  t_stat <- d * sqrt(n_x * n_y / n)

  # Degrees of freedom
df <- n - 2

  # JZS Bayes factor approximation
  # Using the formula from Rouder et al. (2009)
  bf <- bf_jzs(t_stat, n_x, n_y)

  bf
}

#' JZS Bayes factor approximation
#'
#' @param t t-statistic
#' @param n1 Sample size group 1
#' @param n2 Sample size group 2
#' @return Bayes factor
#' @keywords internal
bf_jzs <- function(t, n1, n2) {
  # Approximate JZS Bayes factor using BIC approximation
  n <- n1 + n2
  df <- n - 2

  # BIC-based approximation
  r2 <- t^2 / (t^2 + df)
  bf <- sqrt((n + 1) / n) * ((1 + t^2 / df)^(-(n + 1) / 2))
  bf <- 1 / bf  # BF10

  # Ensure reasonable bounds
  bf <- max(bf, 0.001)
  bf <- min(bf, 1000)

  bf
}

#' Rank Products test
#'
#' Performs rank products test comparing two groups
#' Simplified implementation for two-group comparison
#'
#' @param control Control group values
#' @param treatment Treatment group values
#' @param n_perm Number of permutations for p-value estimation (default 1000)
#'
#' @return p-value from the test
#' @export
test_rankprod <- function(control, treatment, n_perm = 1000) {
  n_c <- length(control)
  n_t <- length(treatment)

  # Compute observed rank product
  # For up-regulated: rank treatment high values
  # For down-regulated: rank control high values

  all_vals <- c(control, treatment)
  ranks <- rank(-all_vals)  # Negative for descending order

  control_ranks <- ranks[1:n_c]
  treatment_ranks <- ranks[(n_c + 1):(n_c + n_t)]

  # Geometric mean of ranks for treatment
  obs_rp_up <- exp(mean(log(treatment_ranks)))
  obs_rp_down <- exp(mean(log(n_c + n_t + 1 - treatment_ranks)))

  # Use the more extreme one
  obs_rp <- min(obs_rp_up, obs_rp_down)

  # Permutation test
  perm_rp <- replicate(n_perm, {
    perm <- sample(n_c + n_t)
    perm_treatment_ranks <- ranks[perm[(n_c + 1):(n_c + n_t)]]
    rp_up <- exp(mean(log(perm_treatment_ranks)))
    rp_down <- exp(mean(log(n_c + n_t + 1 - perm_treatment_ranks)))
    min(rp_up, rp_down)
  })

  # p-value
  p_value <- mean(perm_rp <= obs_rp)
  max(p_value, 1 / n_perm)  # Ensure non-zero
}
