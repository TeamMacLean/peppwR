### R code from vignette source 'jss-article.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("peppwR")
library("dplyr")
library("ggplot2")


###################################################
### code chunk number 2: example-setup
###################################################
# Generate simulated DDA phosphoproteomics data for demonstration
set.seed(123)
n_peptides <- 100
n_samples <- 6

sample_data <- expand.grid(
  peptide = paste0("peptide_", 1:n_peptides),
  condition = c("control", "treatment"),
  replicate = 1:n_samples
) %>%
  mutate(
    abundance = rgamma(n(), shape = 2, rate = 0.1) *
                ifelse(condition == "treatment", 1.5, 1),
    # Add some missing values (MNAR pattern)
    abundance = ifelse(runif(n()) < pmax(0.05, 0.3 - log10(abundance + 1) * 0.1),
                      NA, abundance)
  )

head(sample_data)


###################################################
### code chunk number 3: power-analysis
###################################################
# Demonstrate basic peppwR usage (simplified for template)
cat("# Fit distributions to pilot data\n")
cat("fits <- fit_distributions(pilot_data, id='peptide', group='condition', value='abundance')\n\n")
cat("# Power analysis for sample size determination\n")
cat("power_result <- power_analysis(fits, effect_size=2, target_power=0.8, find='sample_size')\n\n")
cat("# Example result: N=20 samples per group needed for 80% of peptides to achieve 80% power\n")


###################################################
### code chunk number 4: power-distribution
###################################################
# Placeholder for power distribution plot
plot(1:10, rnorm(10), xlab="Sample size", ylab="Power",
     main="Power Distribution Example")


###################################################
### code chunk number 5: test-comparison
###################################################
# Placeholder for test comparison plot
x <- 1:20
wilcox_curve <- 1 - exp(-x/10)
bayes_curve <- 1 - exp(-x/6)

plot(x, wilcox_curve, type="l", col="blue",
     xlab="Sample size per group", ylab="Proportion achieving 80% power",
     main="Statistical Test Comparison")
lines(x, bayes_curve, col="red")
legend("bottomright", c("Wilcoxon", "Bayes t-test"), col=c("blue", "red"), lty=1)
abline(h=0.8, lty=2)


