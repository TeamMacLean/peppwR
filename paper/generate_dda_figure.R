#!/usr/bin/env Rscript
# Generate DDA case study multipanel figure for JSS article

library(peppwR)
library(dplyr)
library(ggplot2)
library(cowplot)

set.seed(42)

# Load and prepare DDA data (extracted from vignette)
dda <- read.csv("../sample_data/dda_data.csv")

pilot <- dda |>
  filter(timepoints %in% c(0, 600)) |>
  transmute(
    peptide_id = new_annotation,
    condition = paste0("t", timepoints),
    abundance = timepoints_values
  )

cat("Processing", n_distinct(pilot$peptide_id), "peptides\n")

# Fit distributions
fits <- fit_distributions(pilot, "peptide_id", "condition", "abundance")

# Panel A: Distribution fitting counts
panel_a <- plot(fits) +
  labs(title = "Distribution Fitting Results") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel B: Statistical test comparison (simplified version with fewer simulations)
# Note: Using n_sim = 50 for faster generation, real analysis used more
power_wilcox <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                               find = "power", test = "wilcoxon", n_sim = 50)

power_bayes <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                              find = "power", test = "bayes_t", n_sim = 50)

# Create comparison data for plotting
test_comparison <- data.frame(
  Test = c("Wilcoxon", "Bayes Factor"),
  Median_Power = c(
    median(power_wilcox$simulations$peptide_power, na.rm = TRUE),
    median(power_bayes$simulations$peptide_power, na.rm = TRUE)
  )
)

panel_b <- ggplot(test_comparison, aes(x = Test, y = Median_Power)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_text(aes(label = round(Median_Power, 3)), vjust = -0.5) +
  labs(title = "Test Performance Comparison",
       x = "Statistical Test",
       y = "Median Power (N=3, 2-fold effect)") +
  theme_minimal() +
  ylim(0, max(test_comparison$Median_Power) * 1.2)

# Panel C: Power heatmap with data-derived parameters
log_abundances <- log(pilot$abundance[pilot$abundance > 0])
derived_params <- list(
  meanlog = median(log_abundances),
  sdlog = mad(log_abundances, constant = 1)
)

panel_c <- plot_power_heatmap(
  distribution = "lnorm",
  params = derived_params,
  n_range = c(3, 12),
  effect_range = c(1.2, 3),
  test = "wilcoxon"
) + labs(title = "Experimental Design Guidance")

# Panel D: Question 1 - Current Power (find = "power")
power_current <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                find = "power", test = "bayes_t", n_sim = 50)
panel_d <- plot(power_current) +
  labs(title = "Q1: Current Power (N=3, 2-fold)") +
  theme_minimal()

# Panel E: Question 2 - Sample Size (find = "sample_size")
sample_size <- power_analysis(fits, effect_size = 2, target_power = 0.8,
                              find = "sample_size", test = "bayes_t", n_sim = 50)
panel_e <- plot(sample_size) +
  labs(title = "Q2: Sample Size for 80% Power") +
  theme_minimal()

# Panel F: Question 3 - Minimum Effect (find = "effect_size")
min_effect <- power_analysis(fits, n_per_group = 3, target_power = 0.8,
                             find = "effect_size", test = "bayes_t", n_sim = 50)
panel_f <- plot(min_effect) +
  labs(title = "Q3: Minimum Detectable Effect") +
  theme_minimal()

# Combine panels in 2x3 grid (3 across, 2 down)
final_plot <- plot_grid(
  panel_a, panel_b, panel_c,
  panel_d, panel_e, panel_f,
  labels = c("A", "B", "C", "D", "E", "F"),
  nrow = 2, ncol = 3
)

# Save figure
ggsave("figures/dda_case_study.pdf", final_plot,
       width = 18, height = 12, device = "pdf")
ggsave("figures/dda_case_study.png", final_plot,
       width = 18, height = 12, dpi = 300)

cat("DDA case study 6-panel figure saved to figures/\n")