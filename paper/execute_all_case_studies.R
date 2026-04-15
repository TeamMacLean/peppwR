# Master Script for All Case Study Computations
# Executes both DDA and PRM case study computations sequentially

cat("======================================================\n")
cat("CASE STUDY COMPUTATION EXECUTION\n")
cat("======================================================\n")
cat("Start time:", format(Sys.time()), "\n")
cat("Working directory:", getwd(), "\n\n")

total_start_time <- Sys.time()

# Execute DDA computation
cat("PHASE 1: DDA CASE STUDY COMPUTATION\n")
cat("----------------------------------------------------\n")
dda_start_time <- Sys.time()

tryCatch({
  source("execute_dda_computation.R")
  dda_end_time <- Sys.time()
  dda_duration <- as.numeric(dda_end_time - dda_start_time, units = "mins")
  cat("DDA computation successful! Duration:", round(dda_duration, 2), "minutes\n\n")
}, error = function(e) {
  cat("ERROR in DDA computation:", e$message, "\n\n")
  stop("DDA computation failed")
})

# Execute PRM computation
cat("PHASE 2: PRM CASE STUDY COMPUTATION\n")
cat("----------------------------------------------------\n")
prm_start_time <- Sys.time()

tryCatch({
  source("execute_prm_computation.R")
  prm_end_time <- Sys.time()
  prm_duration <- as.numeric(prm_end_time - prm_start_time, units = "mins")
  cat("PRM computation successful! Duration:", round(prm_duration, 2), "minutes\n\n")
}, error = function(e) {
  cat("ERROR in PRM computation:", e$message, "\n\n")
  stop("PRM computation failed")
})

# Summary
total_end_time <- Sys.time()
total_duration <- as.numeric(total_end_time - total_start_time, units = "mins")

cat("======================================================\n")
cat("COMPUTATION SUMMARY\n")
cat("======================================================\n")
cat("DDA computation time:", round(dda_duration, 2), "minutes\n")
cat("PRM computation time:", round(prm_duration, 2), "minutes\n")
cat("Total computation time:", round(total_duration, 2), "minutes\n")
cat("End time:", format(Sys.time()), "\n\n")

cat("Generated files:\n")
cat("- jss-article-dda-power-plot.pdf\n")
cat("- jss-article-prm-missingness-plot.pdf\n")
cat("- dda_computation_results.RData\n")
cat("- prm_computation_results.RData\n\n")

cat("All case study computations completed successfully!\n")
cat("======================================================\n")