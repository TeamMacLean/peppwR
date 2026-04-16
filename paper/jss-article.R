# ======================================================================
# R Script for: peppwR: Simulation-Based Power Analysis for
#                Phosphoproteomics Experiments in R
#
# Journal of Statistical Software - JSS Article Code
#
# This script contains all R code from the manuscript for full
# reproducibility. No caching is used - all computations are performed.
# Total runtime: approximately 2+ hours for case studies.
#
# Author: Dan MacLean
# Date: 2026-04-14
# ======================================================================

# ----------------------------------------------------------------------
# INITIAL SETUP
# ----------------------------------------------------------------------

options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE, digits = 7)

# Load required libraries
library("peppwR")
library("dplyr")
library("ggplot2")

cat("peppwR JSS Article - Full Reproducible Script\n")
cat("============================================\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ----------------------------------------------------------------------
# SECTION 3: AGGREGATE MODE EXAMPLES
# ----------------------------------------------------------------------

cat("SECTION 3: AGGREGATE MODE\n")
cat("-------------------------\n")


# Sample Size Determination
cat("3.1 Sample Size Determination\n")
set.seed(123)

result <- power_analysis(
  distribution = "gamma",
  params = list(shape = 2, rate = 0.05),
  effect_size = 2,
  target_power = 0.8,
  find = "sample_size",
  n_sim = 5000
)

print(result)

# Generate sample size plot
pdf("jss-article-sample-size-plot.pdf", width = 7, height = 5)
plot(result)
dev.off()
cat("Sample size plot saved as: jss-article-sample-size-plot.pdf\n\n")

# Power Estimation
cat("3.2 Power Estimation\n")
set.seed(123)

result <- power_analysis(
  distribution = "gamma",
  params = list(shape = 2, rate = 0.05),
  effect_size = 2,
  n_per_group = 6,
  find = "power",
  n_sim = 5000
)

print(result)
cat("\n")

# Minimum Detectable Effect Size
cat("3.3 Minimum Detectable Effect Size\n")
set.seed(123)

result <- power_analysis(
  distribution = "gamma",
  params = list(shape = 2, rate = 0.05),
  n_per_group = 6,
  target_power = 0.8,
  find = "effect_size",
  n_sim = 5000
)

print(result)

# Generate effect size plot
pdf("jss-article-effect-size-plot.pdf", width = 7, height = 5)
plot(result)
dev.off()
cat("Effect size plot saved as: jss-article-effect-size-plot.pdf\n\n")

# ----------------------------------------------------------------------
# SECTION 4: PER-PEPTIDE MODE EXAMPLES
# ----------------------------------------------------------------------

cat("SECTION 4: PER-PEPTIDE MODE\n")
cat("---------------------------\n")

# Generate heterogeneous pilot data
cat("4.1 Generating Pilot Data\n")
set.seed(123)

n_peptides <- 100
n_per_group <- 4

peptide_params <- tibble::tibble(
  peptide_id = paste0("pep_", sprintf("%04d", 1:n_peptides)),
  shape = runif(n_peptides, 1.5, 5),
  rate = runif(n_peptides, 0.01, 0.1)
)

pilot_data <- peptide_params |>
  dplyr::rowwise() |>
  dplyr::mutate(
    data = list(tibble::tibble(
      condition = rep(c("control", "treatment"), each = n_per_group),
      replicate = rep(1:n_per_group, 2),
      abundance = rgamma(n_per_group * 2, shape = shape, rate = rate)
    ))
  ) |>
  dplyr::ungroup() |>
  dplyr::select(peptide_id, data) |>
  tidyr::unnest(data)

cat("Pilot data generated for", n_distinct(pilot_data$peptide_id), "peptides\n")

# Fit distributions
cat("4.2 Distribution Fitting\n")
fits <- fit_distributions(
  pilot_data,
  id = "peptide_id",
  group = "condition",
  value = "abundance",
  distributions = "continuous"
)

print(fits)
cat("\n")

# Per-peptide power analysis
cat("4.3 Per-Peptide Power Analysis\n")
set.seed(123)

result_pp <- power_analysis(
  fits,
  effect_size = 2,
  n_per_group = 6,
  find = "power",
  n_sim = 500
)

print(result_pp)

# Generate per-peptide power plot
pdf("jss-article-per-peptide-power-plot.pdf", width = 6, height = 4)
plot(result_pp)
dev.off()
cat("Per-peptide power plot saved as: jss-article-per-peptide-power-plot.pdf\n\n")

# Per-peptide sample size determination
cat("4.4 Per-Peptide Sample Size Determination\n")
set.seed(123)

result_n <- power_analysis(
  fits,
  effect_size = 2,
  target_power = 0.8,
  find = "sample_size",
  n_sim = 500
)

print(result_n)

# Generate per-peptide sample size plot
pdf("jss-article-per-peptide-sample-size-plot.pdf", width = 6, height = 4)
plot(result_n)
dev.off()
cat("Per-peptide sample size plot saved as: jss-article-per-peptide-sample-size-plot.pdf\n\n")

# ----------------------------------------------------------------------
# SECTION 5.1: DDA CASE STUDY
# ----------------------------------------------------------------------

cat("SECTION 5.1: DDA CASE STUDY\n")
cat("----------------------------\n")

# Load libraries again for clarity
library(peppwR)
library(dplyr)
library(ggplot2)

# Load the DDA experiment data from peppwR package
cat("5.1.1 Loading DDA Data\n")
dda_path <- system.file("extdata", "dda_data.csv", package = "peppwR")
dda <- read.csv(dda_path)

# Filter to early vs late timepoints and format for peppwR
pilot <- dda |>
  filter(timepoints %in% c(0, 600)) |>
  transmute(
    peptide_id = new_annotation,
    condition = paste0("t", timepoints),
    abundance = timepoints_values
  )

# Dataset summary
cat("Unique peptides:", n_distinct(pilot$peptide_id), "\n")
cat("Observations per condition:\n")
print(pilot |> count(condition))

# Fit distributions to each peptide
cat("\n5.1.2 DDA Distribution Fitting\n")
cat("This may take several minutes...\n")
fits <- fit_distributions(pilot, "peptide_id", "condition", "abundance")
print(fits)

# Statistical test comparison
cat("\n5.1.3 Statistical Test Comparison\n")
cat("This may take 30+ minutes for bootstrap test...\n")
set.seed(123)

# Compare three statistical tests at N=3
power_wilcox <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                               find = "power", test = "wilcoxon", n_sim = 500)

power_boot <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                             find = "power", test = "bootstrap_t", n_sim = 500)

power_bayes <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                              find = "power", test = "bayes_t", n_sim = 500)

# Create comparison summary
comparison <- tibble(
  Test = c("Wilcoxon rank-sum", "Bootstrap-t", "Bayes factor"),
  `Median Power` = c(
    median(power_wilcox$simulations$peptide_power, na.rm = TRUE),
    median(power_boot$simulations$peptide_power, na.rm = TRUE),
    median(power_bayes$simulations$peptide_power, na.rm = TRUE)
  ),
  `Peptides >80% Power` = c(
    mean(power_wilcox$simulations$peptide_power > 0.8, na.rm = TRUE) * 100,
    mean(power_boot$simulations$peptide_power > 0.8, na.rm = TRUE) * 100,
    mean(power_bayes$simulations$peptide_power > 0.8, na.rm = TRUE) * 100
  )
)

print(comparison)

# DDA Power Assessment
cat("\n5.1.4 DDA Power Assessment\n")
set.seed(123)

power_current <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                find = "power", test = "bayes_t", n_sim = 500)
print(power_current)

# Generate DDA power plot
pdf("jss-article-dda-power-plot.pdf", width = 6, height = 4)
plot(power_current)
dev.off()
cat("DDA power plot saved as: jss-article-dda-power-plot.pdf\n")

# DDA Sample Size Determination
cat("\n5.1.5 DDA Sample Size Determination\n")
set.seed(123)

# Question 2: Required sample size for adequate power
sample_size_dda <- power_analysis(fits, effect_size = 2, target_power = 0.8,
                                  find = "sample_size", test = "bayes_t", n_sim = 500)
print(sample_size_dda)

# ----------------------------------------------------------------------
# SECTION 5.2: PRM CASE STUDY
# ----------------------------------------------------------------------

cat("\nSECTION 5.2: PRM CASE STUDY\n")
cat("----------------------------\n")

# Load the PRM experiment data from peppwR package
cat("5.2.1 Loading PRM Data\n")
prm_path <- system.file("extdata", "prm_data.csv", package = "peppwR")
prm <- read.csv(prm_path)

# Filter to early (0) vs late (6) timepoints
# Average technical replicates within each biological replicate
pilot <- prm |>
  filter(timepoint %in% c(0, 6)) |>
  group_by(peptide_modified_sequence, genotype, timepoint, bio_rep) |>
  summarise(abundance = mean(total_area, na.rm = TRUE), .groups = "drop") |>
  transmute(
    peptide_id = peptide_modified_sequence,
    condition = paste0("t", timepoint),
    abundance = abundance
  )

# Handle NaN values from averaging (when both tech reps are NA)
pilot <- pilot |>
  mutate(abundance = ifelse(is.nan(abundance), NA, abundance))

cat("PRM data loaded and processed\n")

# Fit distributions - missingness is tracked automatically
cat("\n5.2.2 PRM Distribution Fitting and Missingness Analysis\n")
fits <- fit_distributions(pilot, "peptide_id", "condition", "abundance")

# Display summary including missingness information
print(fits)

# Visualize missingness patterns
pdf("jss-article-prm-missingness-plot.pdf", width = 8, height = 6)
plot_missingness(fits)
dev.off()
cat("PRM missingness plot saved as: jss-article-prm-missingness-plot.pdf\n")

# Statistical test comparison
cat("\n5.2.3 PRM Statistical Test Comparison\n")
# Compare three statistical tests at N=3 with 2-fold effect
power_wilcox <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                               find = "power", test = "wilcoxon", n_sim = 100)

power_boot <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                             find = "power", test = "bootstrap_t", n_sim = 100)

power_bayes <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                              find = "power", test = "bayes_t", n_sim = 100)

# Create comparison table
comparison <- data.frame(
  Test = c("Wilcoxon rank-sum", "Bootstrap-t", "Bayes factor"),
  Median_Power = c(
    median(power_wilcox$simulations$peptide_power, na.rm = TRUE),
    median(power_boot$simulations$peptide_power, na.rm = TRUE),
    median(power_bayes$simulations$peptide_power, na.rm = TRUE)
  ),
  Pct_Above_50 = c(
    mean(power_wilcox$simulations$peptide_power > 0.5, na.rm = TRUE) * 100,
    mean(power_boot$simulations$peptide_power > 0.5, na.rm = TRUE) * 100,
    mean(power_bayes$simulations$peptide_power > 0.5, na.rm = TRUE) * 100
  ),
  Pct_Above_80 = c(
    mean(power_wilcox$simulations$peptide_power > 0.8, na.rm = TRUE) * 100,
    mean(power_boot$simulations$peptide_power > 0.8, na.rm = TRUE) * 100,
    mean(power_bayes$simulations$peptide_power > 0.8, na.rm = TRUE) * 100
  )
)

print(comparison)

# Impact of missingness on power estimates
cat("\n5.2.4 Impact of Missingness on Power\n")
# Power without accounting for missingness (optimistic)
power_no_miss <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                find = "power", test = "bayes_t",
                                include_missingness = FALSE, n_sim = 100)

# Power accounting for missingness (realistic)
power_with_miss <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                  find = "power", test = "bayes_t",
                                  include_missingness = TRUE, n_sim = 100)

cat("Median power WITHOUT missingness:",
    round(median(power_no_miss$simulations$peptide_power, na.rm = TRUE), 3), "\n")
cat("Median power WITH missingness:   ",
    round(median(power_with_miss$simulations$peptide_power, na.rm = TRUE), 3), "\n")

# ----------------------------------------------------------------------
# SCRIPT COMPLETION
# ----------------------------------------------------------------------

cat("\n")
cat("======================================================================\n")
cat("SCRIPT COMPLETED SUCCESSFULLY\n")
cat("======================================================================\n")
cat("End time:", format(Sys.time()), "\n")

# List generated files
cat("\nGenerated files:\n")
generated_files <- c(
  "jss-article-sample-size-plot.pdf",
  "jss-article-effect-size-plot.pdf",
  "jss-article-per-peptide-power-plot.pdf",
  "jss-article-per-peptide-sample-size-plot.pdf",
  "jss-article-dda-power-plot.pdf",
  "jss-article-prm-missingness-plot.pdf"
)

for(file in generated_files) {
  if(file.exists(file)) {
    cat("✓", file, "\n")
  } else {
    cat("✗", file, "(not found)\n")
  }
}


