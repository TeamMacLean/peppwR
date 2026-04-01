#!/usr/bin/env Rscript
# Generate PRM case study multipanel figure for JSS article

library(peppwR)
library(dplyr)
library(ggplot2)
library(cowplot)

set.seed(42)

# Load and prepare PRM data (extracted from vignette)
prm <- read.csv("../sample_data/prm_data.csv")

# Filter and average technical replicates
pilot <- prm |>
  filter(timepoint %in% c(0, 6)) |>
  group_by(peptide_modified_sequence, genotype, timepoint, bio_rep) |>
  summarise(abundance = mean(total_area, na.rm = TRUE), .groups = "drop") |>
  transmute(
    peptide_id = peptide_modified_sequence,
    condition = paste0("t", timepoint),
    abundance = abundance
  ) |>
  mutate(abundance = ifelse(is.nan(abundance), NA, abundance))

cat("Processing", n_distinct(pilot$peptide_id), "peptides\n")
cat("Missing rate:", round(mean(is.na(pilot$abundance)) * 100, 1), "%\n")

# Fit distributions (includes missingness analysis)
fits <- fit_distributions(pilot, "peptide_id", "condition", "abundance")

# Panel A: MNAR correlation plot from plot_missingness()
missingness_plots <- plot_missingness(fits)
# Check if it's a single plot or list, handle accordingly
if (is.list(missingness_plots) && length(missingness_plots) > 1) {
  panel_a <- missingness_plots[[2]] +
    labs(title = "MNAR Pattern Detection") +
    theme_minimal()
} else {
  panel_a <- missingness_plots +
    labs(title = "MNAR Pattern Detection") +
    theme_minimal()
}

# Panel B: Distribution fitting counts
panel_b <- plot(fits) +
  labs(title = "Distribution Fitting Results") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel C: Power comparison with/without missingness
# Simplified version with fewer simulations for faster generation
power_no_miss <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                 find = "power", test = "bayes_t",
                                 include_missingness = FALSE, n_sim = 50)

power_with_miss <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                  find = "power", test = "bayes_t",
                                  include_missingness = TRUE, n_sim = 50)

# Create comparison data
missingness_comparison <- data.frame(
  Condition = c("Without Missingness", "With Missingness"),
  Median_Power = c(
    median(power_no_miss$simulations$peptide_power, na.rm = TRUE),
    median(power_with_miss$simulations$peptide_power, na.rm = TRUE)
  )
)

panel_c <- ggplot(missingness_comparison, aes(x = Condition, y = Median_Power)) +
  geom_col(fill = c("lightblue", "darkblue"), alpha = 0.7) +
  geom_text(aes(label = round(Median_Power, 3)), vjust = -0.5) +
  labs(title = "Impact of Missingness Modeling",
       x = "Analysis Type",
       y = "Median Power (N=3, 2-fold effect)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  ylim(0, max(missingness_comparison$Median_Power) * 1.2)

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
ggsave("figures/prm_case_study.pdf", final_plot,
       width = 18, height = 12, device = "pdf")
ggsave("figures/prm_case_study.png", final_plot,
       width = 18, height = 12, dpi = 300)

cat("PRM case study 6-panel figure saved to figures/\n")