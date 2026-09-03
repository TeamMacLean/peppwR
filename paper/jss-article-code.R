# <<preliminaries, echo=FALSE, results=hide>>=
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE, digits = 7)
library("peppwR")
library("dplyr")
library("ggplot2")

# Cache control: set to TRUE to use fast cached results
USE_CACHED_RESULTS <- TRUE

# <<library-load>>=
library(peppwR)

# <<sample-size-demo>>=
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

# <<sample-size-plot, fig=TRUE, include=FALSE>>=
plot(result)

# <<power-estimation>>=
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

# <<effect-size-demo>>=
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

# <<effect-size-plot, fig=TRUE, include=FALSE>>=
plot(result)

# <<generate-pilot-data>>=
set.seed(123)

# Generate heterogeneous pilot data
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

# <<fit-distributions>>=
fits <- fit_distributions(
  pilot_data,
  id = "peptide_id",
  group = "condition",
  value = "abundance",
  distributions = "continuous"
)

print(fits)

# <<per-peptide-power>>=
set.seed(123)

result_pp <- power_analysis(
  fits,
  effect_size = 2,
  n_per_group = 6,
  find = "power",
  n_sim = 500
)

print(result_pp)

# <<per-peptide-power-plot, fig=TRUE, include=FALSE, width=6, height=4>>=
plot(result_pp)

# <<per-peptide-sample-size>>=
set.seed(123)

result_n <- power_analysis(
  fits,
  effect_size = 2,
  target_power = 0.8,
  find = "sample_size",
  n_sim = 500
)

print(result_n)

# <<per-peptide-sample-size-plot, fig=TRUE, include=FALSE, width=6, height=4>>=
plot(result_n)

# <<load-dda-data>>=
library(peppwR)
library(dplyr)
library(ggplot2)

# Load the DDA experiment data
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
cat("Unique peptides:", n_distinct(pilot$peptide_id), "\n")
cat("Observations per condition:\n")
pilot |> count(condition)

# <<dda-distribution-fitting-display, eval=FALSE>>=
# Fit distributions to each peptide
fits <- fit_distributions(pilot, "peptide_id", "condition", "abundance")

# <<dda-distribution-fitting-execute, echo=FALSE>>=
# Actually execute (with caching)
if(USE_CACHED_RESULTS && file.exists("dda_computation_results.RData") && !exists("fits")) {
  load("dda_computation_results.RData")
} else if(!exists("fits")) {
  fits <- fit_distributions(pilot, "peptide_id", "condition", "abundance")
}

# <<dda-distribution-fitting-output>>=
print(fits)

# <<dda-test-comparison-display, eval=FALSE>>=
set.seed(123)

# Compare three statistical tests at N=3
power_wilcox <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                               find = "power", test = "wilcoxon", n_sim = 500)

power_boot <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                             find = "power", test = "bootstrap_t", n_sim = 500)

# <<dda-test-comparison-execute, echo=FALSE>>=
# Actually execute (with caching)
set.seed(123)
if(!USE_CACHED_RESULTS || !exists("power_wilcox")) {
  power_wilcox <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                 find = "power", test = "wilcoxon", n_sim = 500)
  power_boot <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                               find = "power", test = "bootstrap_t", n_sim = 500)
}

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

# <<dda-current-power>>=
set.seed(123)

power_current <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                find = "power", test = "bayes_t", n_sim = 500)
print(power_current)

# <<dda-power-plot, fig=TRUE, include=FALSE, width=6, height=4>>=
plot(power_current)

# <<dda-sample-size>>=
set.seed(123)

# Question 2: Required sample size for adequate power
sample_size_dda <- power_analysis(fits, effect_size = 2, target_power = 0.8,
                                  find = "sample_size", test = "bayes_t", n_sim = 500)
print(sample_size_dda)

# <<prm-data-loading>>=
# Load the PRM experiment data
prm <- read.csv("../sample_data/prm_data.csv")

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

# <<prm-missingness-analysis-display, eval=FALSE>>=
# Fit distributions - missingness is tracked automatically
fits <- fit_distributions(pilot, "peptide_id", "condition", "abundance")

# <<prm-missingness-analysis-execute, echo=FALSE>>=
# Actually execute (with caching)
if(FALSE && file.exists("prm_computation_results.RData") && !exists("fits")) {
  load("prm_computation_results.RData")
} else  {
  fits <- fit_distributions(pilot, "peptide_id", "condition", "abundance")
}

# <<prm-missingness-analysis-output>>=
# Display summary including missingness information
print(fits)

# <<jss-article-prm-missingness-plot, echo=FALSE, fig=TRUE>>=
# Visualize missingness patterns
plot_missingness(fits)

# <<prm-test-comparison-display, eval=FALSE>>=
# Compare three statistical tests at N=3 with 2-fold effect
power_wilcox <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                               find = "power", test = "wilcoxon", n_sim = 100)

power_boot <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                             find = "power", test = "bootstrap_t", n_sim = 100)

power_bayes <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                              find = "power", test = "bayes_t", n_sim = 100)

# <<prm-test-comparison-execute, echo=FALSE>>=
# Actually execute (with caching)
if(!USE_CACHED_RESULTS || !exists("power_wilcox") || !exists("power_boot") || !exists("power_bayes")) {
  power_wilcox <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                 find = "power", test = "wilcoxon", n_sim = 100)
  power_boot <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                               find = "power", test = "bootstrap_t", n_sim = 100)
  power_bayes <- power_analysis(fits, effect_size = 2, n_per_group = 3,
                                find = "power", test = "bayes_t", n_sim = 100)
}

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

# <<prm-missingness-impact>>=
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

