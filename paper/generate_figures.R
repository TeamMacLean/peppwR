#!/usr/bin/env Rscript
#
# generate_figures.R
#
# This script generates manuscript figures using real sample data included
# with the peppwR package.
# #
# Output files:
#   - dda_power_analysis.pdf/.png: DDA phosphoproteomics power analysis
#   - prm_missingness_analysis.pdf/.png: PRM targeted proteomics analysis with missingness
#
# ==============================================================================
# SETUP AND DATA LOADING
# ==============================================================================

# Load required packages
suppressPackageStartupMessages({
  library(peppwR)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

# Set reproducible seed
set.seed(1549)

# Loading sample data from package

# Load sample data files
# Try multiple paths for flexibility in different environments

# Path 1: Check if we're in the paper directory (JSS submission context)
if (file.exists("../sample_data/dda_data.csv") && file.exists("../sample_data/prm_data.csv")) {
  dda <- read.csv("../sample_data/dda_data.csv")
  prm <- read.csv("../sample_data/prm_data.csv")
} else if (file.exists("sample_data/dda_data.csv") && file.exists("sample_data/prm_data.csv")) {
  # Path 2: Check if we're in the package root directory
  dda <- read.csv("sample_data/dda_data.csv")
  prm <- read.csv("sample_data/prm_data.csv")
} else {
  # Path 3: Try using system.file for installed package
  pkg_root <- system.file(package = "peppwR")
  if (pkg_root != "") {
    # For installed package, data would be in inst/extdata
    dda_path <- system.file("extdata", "dda_data.csv", package = "peppwR")
    prm_path <- system.file("extdata", "prm_data.csv", package = "peppwR")

    if (file.exists(dda_path) && file.exists(prm_path)) {
      dda <- read.csv(dda_path)
      prm <- read.csv(prm_path)
    } else {
      stop("Sample data files not found in package installation.")
    }
  } else {
    stop("Sample data files not found. Please ensure script is run from package directory or peppwR package is installed.")
  }
}

# Prepare DDA pilot data (subset to t=0 vs t=600 timepoints)
dda_pilot <- dda %>%
  filter(timepoints %in% c(0, 600)) %>%
  transmute(
    peptide_id = new_annotation,
    condition = paste0("t", timepoints),
    abundance = timepoints_values
  )

# Prepare PRM pilot data (subset to t=0 vs t=6 timepoints)
prm_pilot <- prm %>%
  filter(timepoint %in% c(0, 6)) %>%
  group_by(peptide_modified_sequence, genotype, timepoint, bio_rep) %>%
  summarise(abundance = mean(total_area, na.rm = TRUE), .groups = "drop") %>%
  transmute(
    peptide_id = peptide_modified_sequence,
    condition = paste0("t", timepoint),
    abundance = ifelse(is.nan(abundance), NA, abundance)
  )

# Data loaded successfully

# Fit distributions using peppwR functions
dda_fits <- fit_distributions(dda_pilot, "peptide_id", "condition", "abundance")
prm_fits <- fit_distributions(prm_pilot, "peptide_id", "condition", "abundance")

# ==============================================================================
# DDA PHOSPHOPROTEOMICS POWER ANALYSIS
# ==============================================================================

# Creating DDA Power Analysis

# Panel A: Power heatmap using aggregate mode with DDA-derived parameters
# Extract lognormal parameters from actual DDA data for realistic simulations
dda_log_vals <- log(dda_pilot$abundance[dda_pilot$abundance > 0])
dda_lnorm_params <- list(
  meanlog = mean(dda_log_vals, na.rm = TRUE),
  sdlog = sd(dda_log_vals, na.rm = TRUE)
)

# Define parameter space for heatmap
effect_sizes <- c(1.5, 2, 2.5, 3)
sample_sizes <- c(4, 6, 8, 10, 12)

# Generate power heatmap data using aggregate mode
heatmap_results <- list()
for(i in seq_along(effect_sizes)) {
  for(j in seq_along(sample_sizes)) {
    power_result <- power_analysis(
      distribution = "lnorm",
      params = dda_lnorm_params,
      effect_size = effect_sizes[i],
      n_per_group = sample_sizes[j],
      find = "power",
      n_sim = 300
    )
    heatmap_results[[length(heatmap_results) + 1]] <- data.frame(
      effect_size = effect_sizes[i],
      n_per_group = sample_sizes[j],
      power = power_result$answer * 100
    )
  }
}
heatmap_data <- do.call(rbind, heatmap_results)

# Create power heatmap panel
panel_1a <- ggplot(heatmap_data, aes(x = factor(effect_size), y = factor(n_per_group),
                                     fill = power)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(power)), size = 3, color = "white") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 50, limits = c(0, 100)) +
  labs(
    x = "Effect size (fold-change)",
    y = "Sample size per group",
    title = "A) Power heatmap (DDA parameters)",
    fill = "Power (%)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    panel.grid = element_blank()
  )

# Panel B: Sample size requirement curves using per-peptide analysis
sample_sizes_curve <- 3:12

# Calculate percentage of peptides achieving 80% power for different sample sizes
wilcox_results <- sapply(sample_sizes_curve, function(n) {
  power_res <- power_analysis(dda_fits, effect_size = 2, n_per_group = n,
                             find = "power", test = "wilcoxon", n_sim = 200)
  mean(power_res$simulations$peptide_power >= 0.8, na.rm = TRUE) * 100
})

bayes_results <- sapply(sample_sizes_curve, function(n) {
  power_res <- power_analysis(dda_fits, effect_size = 2, n_per_group = n,
                             find = "power", test = "bayes_t", n_sim = 200)
  mean(power_res$simulations$peptide_power >= 0.8, na.rm = TRUE) * 100
})

curve_data <- data.frame(
  n_per_group = rep(sample_sizes_curve, 2),
  power = c(wilcox_results, bayes_results),
  test = rep(c("Wilcoxon", "Bayes factor"), each = length(sample_sizes_curve))
)

# Create sample size curves panel
panel_1b <- ggplot(curve_data, aes(x = n_per_group, y = power, color = test, linetype = test)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 80, linetype = "dotted", color = "gray40") +
  scale_color_manual(values = c("Wilcoxon" = "#D73027", "Bayes factor" = "#4575B4")) +
  scale_x_continuous(breaks = seq(3, 12, 2)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    x = "Sample size per group",
    y = "% peptides with power ≥ 80%",
    title = "B) Sample size requirements (DDA)",
    color = "Test", linetype = "Test"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    legend.position = c(0.7, 0.3),
    panel.grid.minor = element_blank()
  )

# Panel C: FDR correction impact

# Per-peptide analysis without FDR correction
regular_power <- power_analysis(dda_fits, effect_size = 2, n_per_group = 8,
                               find = "power", test = "wilcoxon", n_sim = 200)

# FDR-aware analysis with Benjamini-Hochberg correction
fdr_power <- power_analysis(dda_fits, effect_size = 2, n_per_group = 8,
                           find = "power", test = "wilcoxon", n_sim = 200,
                           apply_fdr = TRUE, prop_null = 0.9)

fdr_data <- data.frame(
  Analysis = c("Per-peptide", "FDR-adjusted"),
  Power = c(
    median(regular_power$simulations$peptide_power, na.rm = TRUE) * 100,
    median(fdr_power$simulations$peptide_power, na.rm = TRUE) * 100
  )
)

# Create FDR impact panel
panel_1c <- ggplot(fdr_data, aes(x = Analysis, y = Power, fill = Analysis)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = paste0(round(Power), "%")), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("#5AAE61", "#F46D43")) +
  scale_y_continuous(limits = c(0, max(fdr_data$Power) * 1.2)) +
  labs(
    x = "Analysis type",
    y = "Median power (%)",
    title = "C) FDR correction impact (DDA)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# Combine and save DDA power analysis figure
dda_power_analysis <- plot_grid(panel_1a, panel_1b, panel_1c, nrow = 1, align = "hv")
ggsave("dda_power_analysis.pdf", dda_power_analysis, width = 12, height = 4, device = "pdf")
ggsave("dda_power_analysis.png", dda_power_analysis, width = 12, height = 4, dpi = 300)

# ==============================================================================
# PRM TARGETED PROTEOMICS WITH MISSINGNESS ANALYSIS
# ==============================================================================

# Creating PRM missingness analysis

# Panel A: MNAR (Missing Not At Random) detection
# Analyze relationship between abundance and missingness rate
miss_detailed <- prm_pilot %>%
  group_by(peptide_id) %>%
  summarise(
    miss_rate = mean(is.na(abundance)) * 100,
    mean_abundance = mean(abundance, na.rm = TRUE),
    n_obs = sum(!is.na(abundance)),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean_abundance) & mean_abundance > 0)

# Calculate correlation with statistical test
cor_test_result <- cor.test(log10(miss_detailed$mean_abundance),
                           miss_detailed$miss_rate, use = "complete.obs")

# Create MNAR detection panel
panel_2a <- ggplot(miss_detailed, aes(x = log10(mean_abundance), y = miss_rate)) +
  geom_point(alpha = 0.6, color = "#762A83", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "#2166AC", linewidth = 1) +
  labs(
    x = "Mean abundance (log₁₀)",
    y = "Missingness rate (%)",
    title = "A) MNAR detection (PRM)",
    subtitle = bquote(rho == .(round(cor_test_result$estimate, 3)) ~
                     ", p" == .(format.pval(cor_test_result$p.value, digits = 3)))
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

# Panel B: Impact of missingness on power analysis

# Power analysis excluding missing data (complete cases only)
prm_complete <- prm_pilot %>% filter(!is.na(abundance))
prm_fits_complete <- fit_distributions(prm_complete, "peptide_id", "condition", "abundance")

power_no_missing <- power_analysis(prm_fits_complete, effect_size = 2, n_per_group = 6,
                                  find = "power", test = "wilcoxon", n_sim = 200)

# Power analysis including missingness patterns
power_with_missing <- power_analysis(prm_fits, effect_size = 2, n_per_group = 6,
                                    find = "power", test = "wilcoxon", n_sim = 200)

missing_impact_data <- data.frame(
  Analysis = c("Complete cases only", "With missingness"),
  Power = c(
    median(power_no_missing$simulations$peptide_power, na.rm = TRUE) * 100,
    median(power_with_missing$simulations$peptide_power, na.rm = TRUE) * 100
  )
)

# Calculate power reduction percentage
power_reduction <- (missing_impact_data$Power[1] - missing_impact_data$Power[2]) /
                   missing_impact_data$Power[1] * 100

# Create missingness impact panel
panel_2b <- ggplot(missing_impact_data, aes(x = Analysis, y = Power, fill = Analysis)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = paste0(round(Power), "%")), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_y_continuous(limits = c(0, max(missing_impact_data$Power) * 1.2)) +
  labs(
    x = "Analysis approach",
    y = "Median power (%)",
    title = "B) Missingness impact on power",
    subtitle = paste0("Power reduction: ", round(power_reduction, 1), "%")
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# Panel C: Distribution fitting results for PRM data
prm_dist_summary <- table(prm_fits$best) %>%
  as.data.frame() %>%
  rename(Distribution = Var1, Count = Freq) %>%
  arrange(desc(Count)) %>%
  mutate(
    Percentage = Count / sum(Count) * 100,
    Distribution = factor(Distribution, levels = rev(Distribution))
  )

n_fitted <- sum(prm_dist_summary$Count)

# Create distribution fitting panel
panel_2c <- ggplot(prm_dist_summary, aes(x = Distribution, y = Count, fill = Distribution)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            hjust = -0.1, size = 3) +
  coord_flip() +
  labs(
    x = "Best-fit distribution",
    y = "Number of peptides",
    title = "C) Distribution fitting (PRM)",
    subtitle = paste0("n = ", n_fitted, " peptides fitted")
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray60"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# Combine and save PRM missingness analysis figure
prm_missingness_analysis <- plot_grid(panel_2a, panel_2b, panel_2c, nrow = 1, align = "hv")
ggsave("prm_missingness_analysis.pdf", prm_missingness_analysis, width = 12, height = 4, device = "pdf")
ggsave("prm_missingness_analysis.png", prm_missingness_analysis, width = 12, height = 4, dpi = 300)

