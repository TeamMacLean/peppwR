# DDA Case Study Computation Script
# This script executes all computations for the DDA case study section

library(peppwR)
library(dplyr)
library(ggplot2)

cat("=== DDA CASE STUDY COMPUTATION ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# Load the DDA experiment data
cat("Loading DDA data...\n")
dda <- read.csv("../sample_data/dda_data.csv")

# Filter to early vs late timepoints and format for peppwR
pilot <- dda |>
  filter(timepoints %in% c(0, 600)) |>
  transmute(
    peptide_id = new_annotation,
    condition = paste0("t", timepoints),
    abundance = timepoints_values
  )

# Dataset summary
cat("DDA Dataset Summary:\n")
cat("Unique peptides:", n_distinct(pilot$peptide_id), "\n")
pilot |> count(condition) |> print()

# Fit distributions to each peptide
cat("\nStarting distribution fitting for", n_distinct(pilot$peptide_id), "peptides...\n")
start_fit_time <- Sys.time()
fits <- fit_distributions(pilot, "peptide_id", "condition", "abundance")
end_fit_time <- Sys.time()
cat("Distribution fitting completed in:", round(as.numeric(end_fit_time - start_fit_time, units = "mins"), 2), "minutes\n")

print(fits)

# Statistical test comparison
cat("\n=== STATISTICAL TEST COMPARISON ===\n")
set.seed(123)

cat("Running Wilcoxon test...\n")
start_wilcox <- Sys.time()
power_wilcox <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                               find = "power", test = "wilcoxon", n_sim = 500)
end_wilcox <- Sys.time()
cat("Wilcoxon completed in:", round(as.numeric(end_wilcox - start_wilcox, units = "mins"), 2), "minutes\n")

cat("Running Bootstrap-t test...\n")
start_boot <- Sys.time()
power_boot <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                             find = "power", test = "bootstrap_t", n_sim = 500)
end_boot <- Sys.time()
cat("Bootstrap-t completed in:", round(as.numeric(end_boot - start_boot, units = "mins"), 2), "minutes\n")

cat("Running Bayes factor test...\n")
start_bayes <- Sys.time()
power_bayes <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                              find = "power", test = "bayes_t", n_sim = 500)
end_bayes <- Sys.time()
cat("Bayes factor completed in:", round(as.numeric(end_bayes - start_bayes, units = "mins"), 2), "minutes\n")

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

cat("\nTest Comparison Results:\n")
print(comparison)

# Current power assessment
cat("\n=== CURRENT POWER ASSESSMENT ===\n")
set.seed(123)
start_current <- Sys.time()
power_current <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                find = "power", test = "bayes_t", n_sim = 500)
end_current <- Sys.time()
cat("Current power analysis completed in:", round(as.numeric(end_current - start_current, units = "mins"), 2), "minutes\n")

print(power_current)

# Generate DDA power plot
cat("\n=== GENERATING DDA POWER PLOT ===\n")
pdf("jss-article-dda-power-plot.pdf", width = 6, height = 4)
plot(power_current)
dev.off()
cat("DDA power plot saved as: jss-article-dda-power-plot.pdf\n")

# Sample size determination
cat("\n=== SAMPLE SIZE DETERMINATION ===\n")
set.seed(123)
start_sample <- Sys.time()
sample_size_dda <- power_analysis(fits, effect_size = 2, target_power = 0.8,
                                  find = "sample_size", test = "bayes_t", n_sim = 500)
end_sample <- Sys.time()
cat("Sample size analysis completed in:", round(as.numeric(end_sample - start_sample, units = "mins"), 2), "minutes\n")

print(sample_size_dda)

# Save results for integration
save(fits, power_wilcox, power_boot, power_bayes, power_current, sample_size_dda,
     file = "dda_computation_results.RData")

total_time <- Sys.time()
cat("\n=== DDA COMPUTATION COMPLETED ===\n")
cat("Total computation time: TBD\n")
cat("End time:", format(Sys.time()), "\n")
cat("Results saved to: dda_computation_results.RData\n")
cat("Figure saved to: jss-article-dda-power-plot.pdf\n")