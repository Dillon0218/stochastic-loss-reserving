# ==============================================================================
# 03_stochastic.R
#
# Full predictive distributions of the lag-10 outcome for every triangle, from
# the two standard stochastic reserving models:
#
#   Mack (1993)            distribution-free chain ladder: gives the mean and the
#                          mean squared error of prediction (MSEP). As in Mack's
#                          paper and in Meyers (2015), a lognormal with the same
#                          mean and standard error turns that into a distribution.
#                          Run on paid AND incurred triangles.
#   Bootstrap ODP          England & Verrall (2002): overdispersed Poisson GLM on
#                          incremental paid losses, Pearson residuals resampled
#                          10,000 times, gamma process error. Paid only: the ODP
#                          variance is phi x mean, so it needs non-negative
#                          incrementals, and incurred incrementals routinely go
#                          negative when case reserves are released.
#
# For each model the script records the predicted mean, standard deviation, the
# 99.5th percentile, and -- the input to the backtest -- the percentile at which
# the ACTUAL outcome landed.
#
# Runs on the screened sample and on Meyers' (2015) 200 triangles.
# ==============================================================================

source("R/00_common.R")
suppressPackageStartupMessages(library(ChainLadder))

tri <- readRDS(file.path(DIR_PROC, "triangles.rds"))
cat(sprintf("Running Mack (paid, incurred) and bootstrap ODP (paid, %s sims) on %d triangles\n",
            format(N_BOOT, big.mark = ","), length(tri)))

#' Mack on one upper triangle -> total mean, s.e., percentile of the outcome.
#' A link ratio from a zero or negative base is undefined, so those cells get
#' zero weight. The screened sample has none (rule 2 in 01_prepare_data.R);
#' 17 of Meyers' 200 triangles do, and without this ChainLadder returns NaN.
fit_mack <- function(up, outcome) {
  w  <- matrix(1, nrow(up), ncol(up))
  w[!is.na(up) & up <= 0] <- 0
  m  <- MackChainLadder(up, weights = w, est.sigma = "log-linear")
  mu <- sum(m$FullTriangle[, N_AY])
  se <- as.numeric(m$Total.Mack.S.E)
  if (!is.finite(se) || se <= 0) return(list(ok = FALSE, reason = "Mack standard error not finite"))
  list(ok = TRUE, mean = mu, sd = se,
       q995 = qnt_lnorm(0.995, mu, se),
       pct  = pct_lnorm(outcome, mu, se))
}

#' Bootstrap ODP on one paid upper triangle.
#' Run as ChainLadder ships it, as Meyers did, so the replication compares like
#' with like. The share of pseudo-triangles that went non-positive is recorded
#' so 04_backtest.R can check the verdict without the unstable ones.
#' One silent failure is caught here: if a development column nets to exactly
#' zero (link ratio 1.000) while cells in it move, a fitted incremental of 0
#' gives an infinite Pearson residual, the scale parameter becomes infinite,
#' the gamma draws come back NA, and ChainLadder fills them with 0 -- every
#' simulated reserve is zero. That is a failed fit, not a prediction.
fit_odp <- function(up, outcome, seed) {
  b    <- BootChainLadder(up, R = N_BOOT, process.distr = "gamma", seed = seed)
  sims <- sum(latest_diag(up)) + b$IBNR.Totals
  if (!all(is.finite(sims)) || sd(sims) == 0)
    return(list(ok = FALSE, reason = "bootstrap scale parameter infinite (all simulated reserves zero)"))
  list(ok = TRUE, mean = mean(sims), sd = sd(sims),
       q995 = unname(quantile(sims, 0.995)),
       pct  = pct_sim(outcome, sims),
       invalid_share = pseudo_invalid_share(b$simClaims))
}

# ChainLadder flags development periods with no variation (common in the last
# few lags of small triangles) as an "Information:" warning. Muffle exactly
# that message; anything else still surfaces.
quiet_info <- function(expr) withCallingHandlers(expr, warning = function(w) {
  if (grepl("^Information: essentially no variation", conditionMessage(w)))
    invokeRestart("muffleWarning")
})

t0 <- Sys.time()
res <- quiet_info(rbindlist(lapply(seq_along(tri), function(k) {
  t <- tri[[k]]
  # Seed from the identity of the triangle, not its position in the list, so
  # each triangle's simulation is reproducible whatever else is in the run.
  seed <- t$grcode * 10L + match(t$line, LINES$line)
  out <- list()
  for (basis in c("paid", "incurred")) {
    up <- upper_only(t[[basis]])
    outcome <- sum(t[[basis]][, N_AY])
    r <- tryCatch(fit_mack(up, outcome), error = function(e) list(ok = FALSE, reason = conditionMessage(e)))
    out[[length(out) + 1L]] <- c(list(model = paste("Mack", basis), outcome = outcome), r)
    if (basis == "paid") {
      r <- tryCatch(suppressWarnings(fit_odp(up, outcome, seed)),
                    error = function(e) list(ok = FALSE, reason = conditionMessage(e)))
      out[[length(out) + 1L]] <- c(list(model = "ODP paid", outcome = outcome), r)
    }
  }
  if (k %% 50 == 0) cat(sprintf("  %d / %d triangles (%.0f s)\n", k, length(tri),
                                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  d <- rbindlist(out, fill = TRUE)
  d[, `:=`(line = t$line, grcode = t$grcode, grname = t$grname,
           selected = t$selected, meyers200 = t$meyers200)]
  d
}), fill = TRUE))
cat(sprintf("Done in %.0f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

setcolorder(res, c("line", "grcode", "grname", "selected", "meyers200", "model"))
res[, model := factor(model, c("Mack incurred", "Mack paid", "ODP paid"))]

# Fits that did not produce a usable distribution are reported, then dropped.
# (Meyers hit the same thing: ChainLadder returns NA for the Mack standard error
# when a whole development column has no variation.)
fails <- res[ok == FALSE | is.na(ok)]
cat(sprintf("\nFits without a usable distribution: %d\n", nrow(fails)))
if (nrow(fails)) print(fails[, .(line, grcode, model, selected, meyers200, reason)])
save_table(fails[, .(line, grcode, grname, model, selected, meyers200, reason)], "03_failed_fits")

res <- res[ok == TRUE]
res[, exceed995 := outcome > q995]
save_table(res, "03_predictive_by_triangle")
saveRDS(res, file.path(DIR_PROC, "predictive.rds"))

# Width of the two paid models side by side, same screened triangles
w <- dcast(res[selected == TRUE & model != "Mack incurred"], line + grcode ~ model, value.var = "sd")
w <- w[complete.cases(w)]
cat(sprintf("\nScreened sample, ODP sd / Mack sd on the same paid triangle: median %.2f (IQR %.2f-%.2f)\n",
            median(w$`ODP paid` / w$`Mack paid`),
            quantile(w$`ODP paid` / w$`Mack paid`, 0.25),
            quantile(w$`ODP paid` / w$`Mack paid`, 0.75)))
