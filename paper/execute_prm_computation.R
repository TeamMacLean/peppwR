# PRM Case Study Computation Script
# This script executes all computations for the PRM case study section

library(peppwR)
library(dplyr)
library(ggplot2)

cat("=== PRM CASE STUDY COMPUTATION ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# Load the PRM experiment data
cat("Loading PRM data...\n")
prm <- read.csv("../sample_data/prm_data.csv")

# Check for missing values in raw data
cat("Raw PRM Data Summary:\n")
cat("Total observations:", nrow(prm), "\n")
cat("Missing values:", sum(is.na(prm$total_area)), "\n")
cat("Missing rate:", round(mean(is.na(prm$total_area)) * 100, 1), "%\n\n")

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

cat("After averaging technical replicates:\n")
cat("Unique peptides:", n_distinct(pilot$peptide_id), "\n")
cat("Missing values:", sum(is.na(pilot$abundance)), "\n")
cat("Missing rate:", round(mean(is.na(pilot$abundance)) * 100, 1), "%\n\n")

# Fit distributions with missingness analysis
cat("Starting distribution fitting and missingness analysis...\n")
start_fit_time <- Sys.time()
fits <- fit_distributions(pilot, "peptide_id", "condition", "abundance")
end_fit_time <- Sys.time()
cat("Distribution fitting completed in:", round(as.numeric(end_fit_time - start_fit_time, units = "mins"), 2), "minutes\n")

print(fits)

# Generate missingness plot
cat("\n=== GENERATING PRM MISSINGNESS PLOT ===\n")
pdf("jss-article-prm-missingness-plot.pdf", width = 8, height = 4)
plot_missingness(fits)
dev.off()
cat("PRM missingness plot saved as: jss-article-prm-missingness-plot.pdf\n")

# Statistical test comparison
cat("\n=== STATISTICAL TEST COMPARISON ===\n")
set.seed(123)

cat("Running Wilcoxon test...\n")
start_wilcox <- Sys.time()
power_wilcox <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                               find = "power", test = "wilcoxon", n_sim = 100)
end_wilcox <- Sys.time()
cat("Wilcoxon completed in:", round(as.numeric(end_wilcox - start_wilcox, units = "mins"), 2), "minutes\n")

cat("Running Bootstrap-t test...\n")
start_boot <- Sys.time()
power_boot <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                             find = "power", test = "bootstrap_t", n_sim = 100)
end_boot <- Sys.time()
cat("Bootstrap-t completed in:", round(as.numeric(end_boot - start_boot, units = "mins"), 2), "minutes\n")

cat("Running Bayes factor test...\n")
start_bayes <- Sys.time()
power_bayes <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                              find = "power", test = "bayes_t", n_sim = 100)
end_bayes <- Sys.time()
cat("Bayes factor completed in:", round(as.numeric(end_bayes - start_bayes, units = "mins"), 2), "minutes\n")

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

cat("\nTest Comparison Results:\n")
print(comparison)

# Impact of missingness
cat("\n=== IMPACT OF MISSINGNESS ON POWER ===\n")
cat("Running power analysis without missingness...\n")
start_no_miss <- Sys.time()
power_no_miss <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                find = "power", test = "bayes_t",
                                include_missingness = FALSE, n_sim = 100)
end_no_miss <- Sys.time()
cat("Completed in:", round(as.numeric(end_no_miss - start_no_miss, units = "mins"), 2), "minutes\n")

cat("Running power analysis with missingness...\n")
start_with_miss <- Sys.time()
power_with_miss <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                  find = "power", test = "bayes_t",
                                  include_missingness = TRUE, n_sim = 100)
end_with_miss <- Sys.time()
cat("Completed in:", round(as.numeric(end_with_miss - start_with_miss, units = "mins"), 2), "minutes\n")

cat("Median power WITHOUT missingness:",
    round(median(power_no_miss$simulations$peptide_power, na.rm = TRUE), 3), "\n")
cat("Median power WITH missingness:   ",
    round(median(power_with_miss$simulations$peptide_power, na.rm = TRUE), 3), "\n")

# FDR-aware analysis
cat("\n=== FDR-AWARE POWER ANALYSIS ===\n")
cat("Running standard power analysis...\n")
start_nominal <- Sys.time()
power_nominal <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                               find = "power", test = "wilcoxon",
                               apply_fdr = FALSE, n_sim = 100)
end_nominal <- Sys.time()
cat("Standard power analysis completed in:", round(as.numeric(end_nominal - start_nominal, units = "mins"), 2), "minutes\n")

cat("Running FDR-aware power analysis...\n")
start_fdr <- Sys.time()

# Try FDR analysis with error handling
tryCatch({
  power_fdr <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                             find = "power", test = "wilcoxon",
                             apply_fdr = TRUE, prop_null = 0.8,
                             fdr_threshold = 0.05, n_sim = 100)
  end_fdr <- Sys.time()
  cat("FDR-aware power analysis completed in:", round(as.numeric(end_fdr - start_fdr, units = "mins"), 2), "minutes\n")

  # Extract power values with safe handling
  nominal_power <- power_nominal$simulations$peptide_power
  fdr_power <- power_fdr$simulations$peptide_power

  cat("Nominal power (no FDR correction):\n")
  cat("  Median power:", round(median(nominal_power, na.rm = TRUE), 3), "\n")
  cat("  % peptides > 80% power:",
      round(mean(nominal_power > 0.8, na.rm = TRUE) * 100, 1), "%\n")

  cat("\nFDR-aware power (BH correction, 80% true nulls):\n")
  cat("  Median power:", round(median(fdr_power, na.rm = TRUE), 3), "\n")
  cat("  % peptides > 80% power:",
      round(mean(fdr_power > 0.8, na.rm = TRUE) * 100, 1), "%\n")

}, error = function(e) {
  cat("\nFDR-aware analysis encountered an error (likely due to limited data):\n")
  cat("Error:", e$message, "\n")
  cat("Proceeding with nominal power results only.\n")

  # Just report nominal power
  nominal_power <- power_nominal$simulations$peptide_power
  power_fdr <- NULL  # Set to NULL to handle in save statement

  cat("Nominal power (no FDR correction):\n")
  cat("  Median power:", round(median(nominal_power, na.rm = TRUE), 3), "\n")
  cat("  % peptides > 80% power:",
      round(mean(nominal_power > 0.8, na.rm = TRUE) * 100, 1), "%\n")
  cat("\nFDR-aware analysis skipped due to computational limitations.\n")
})

# Save results for integration
if(exists("power_fdr") && !is.null(power_fdr)) {
  save(fits, power_wilcox, power_boot, power_bayes, power_no_miss, power_with_miss,
       power_nominal, power_fdr, comparison,
       file = "prm_computation_results.RData")
} else {
  save(fits, power_wilcox, power_boot, power_bayes, power_no_miss, power_with_miss,
       power_nominal, comparison,
       file = "prm_computation_results.RData")
}

total_time <- Sys.time()
cat("\n=== PRM COMPUTATION COMPLETED ===\n")
cat("Total computation time: TBD\n")
cat("End time:", format(Sys.time()), "\n")
cat("Results saved to: prm_computation_results.RData\n")
cat("Figure saved to: jss-article-prm-missingness-plot.pdf\n")