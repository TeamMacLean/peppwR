# Utility functions for peppwR

# Global variables used in NSE (non-standard evaluation)
# This avoids R CMD check NOTEs about undefined global variables
utils::globalVariables(c(
  ".data",
  "aic",
  "best",
  "count",
  "data",
  "dist",
  "effect_size",
  "fits",
  "n",
  "n_per_group",
  "power",
  "proportion_powered",
  "x"
))
