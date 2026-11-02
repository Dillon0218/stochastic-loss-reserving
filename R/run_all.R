# ==============================================================================
# run_all.R -- run the full pipeline end to end.
#   Rscript R/run_all.R
# Expects the six CAS files in data/raw/ (see scripts/download_data.sh).
# ==============================================================================

steps <- c("R/01_prepare_data.R",
           "R/02_deterministic.R",
           "R/03_stochastic.R",
           "R/04_backtest.R",
           "R/05_capital.R",
           "R/06_figures.R")

# Dotted names for the loop state: the step scripts are sourced into this same
# environment and freely use short names like `tri`, `res`, `d`. Anything not
# dotted would be silently clobbered mid-loop.
for (.step in steps) {
  cat(sprintf("\n\n################ %s ################\n", .step))
  .t0 <- Sys.time()
  source(.step, echo = FALSE)
  cat(sprintf("\n[%s completed in %.1f s]\n", .step,
              as.numeric(difftime(Sys.time(), .t0, units = "secs"))))
  # Clear each step's objects so a later script cannot accidentally read a
  # stale in-memory object instead of the .rds it is supposed to load.
  rm(list = setdiff(ls(all.names = TRUE), c("steps", ".step", ".t0")))
}
cat("\nPipeline complete. Tables in outputs/tables/, figures in outputs/figures/.\n")
