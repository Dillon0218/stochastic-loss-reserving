# ==============================================================================
# 04_backtest.R
#
# Did the predicted distributions contain what actually happened?
#
# For each triangle, 03_stochastic.R recorded the percentile p at which the
# actual lag-10 outcome fell in the model's predictive distribution. If the
# model is right, p is uniform on (0, 1) across triangles: an outcome is as
# likely to land at the 10th percentile as at the 90th. Departures have a
# recognisable shape (Meyers 2015, Figure 1):
#
#   too many p near 0 AND near 1   distribution too narrow (light tails)
#   too many p near 0              model biased high (it over-reserves)
#   too many p near 1              model biased low
#
# Tests: Kolmogorov-Smirnov against U(0,1), shown on a p-p plot with the 5%
# band 1.36/sqrt(n). Plus the number that matters for capital: how often the
# outcome exceeded the model's own 99.5th percentile, against the 0.5% the
# model promises.
#
# Then the same calculation on Meyers' 200 triangles, against his published
# KS statistics, as a check on this implementation.
# ==============================================================================

source("R/00_common.R")

res <- readRDS(file.path(DIR_PROC, "predictive.rds"))
main <- res[selected == TRUE]
cat(sprintf("Screened sample: %d triangles x 3 models\n", uniqueN(main[, .(line, grcode)])))

# ------------------------------------------------------------------------------
# 1. UNIFORMITY TESTS BY MODEL AND LINE
# ------------------------------------------------------------------------------
uniformity <- function(p, exceed) {
  n  <- length(p)
  ks <- suppressWarnings(ks.test(p, "punif"))
  list(n = n,
       ks_d = unname(ks$statistic), ks_p = ks$p.value,
       ks_crit_5pct = 1.36 / sqrt(n),
       mean_pct = mean(p),
       share_below_5 = mean(p < 0.05), share_above_95 = mean(p > 0.95),
       n_exceed_995 = sum(exceed), exp_exceed_995 = 0.005 * n,
       # one-sided binomial: is the 99.5th percentile breached more than it should be?
       exceed_p = binom.test(sum(exceed), n, 0.005, alternative = "greater")$p.value)
}

by_line <- main[, uniformity(pct, exceed995), by = .(model, line)]
pooled  <- main[, uniformity(pct, exceed995), by = .(model)][, line := "All"]
bt <- rbind(by_line, pooled, use.names = TRUE)
bt[, line := factor(line, c(LINES$line, "All"))]
setorder(bt, model, line)

cat("\n=========== UNIFORMITY OF OUTCOME PERCENTILES (screened sample) ===========\n")
print(bt[line == "All", .(model, n, ks_d = round(ks_d, 3), ks_p = signif(ks_p, 2),
                          crit = round(ks_crit_5pct, 3), mean_pct = round(mean_pct, 3),
                          below_5 = round(share_below_5, 3), above_95 = round(share_above_95, 3),
                          n_exceed_995, exp_exceed_995)])
cat("\nKS D, then KS p-value, by line:\n")
print(dcast(bt, line ~ model, value.var = "ks_d")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
print(dcast(bt, line ~ model, value.var = "ks_p")[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 2) else x)])
save_table(bt, "04_uniformity_tests")

# p-p plot coordinates: sorted percentiles against their expected values
pp <- rbind(main[, .(line = line, pct = pct), by = model],
            main[, .(line = "All", pct = pct), by = model])
pp <- pp[order(model, line, pct)]
pp[, `:=`(expected = seq_len(.N) / (.N + 1), n = .N), by = .(model, line)]
save_table(pp, "04_pp_points")

# Histogram counts in deciles, the other way Meyers shows the same thing
dec <- main[, .(n = .N), by = .(model, decile = pmin(10L, floor(pct * 10) + 1L))]
dec <- dcast(dec, decile ~ model, value.var = "n", fill = 0)
cat("\nOutcomes by predicted-percentile decile (all lines; about 17 per decile if uniform):\n")
print(dec)
save_table(dec, "04_percentile_deciles")

# ------------------------------------------------------------------------------
# 2. TAIL CALIBRATION: how often was the model's own 99.5th percentile breached?
# This is a VaR backtest in the market-risk sense, one observation per
# triangle: a 99.5% quantile should be exceeded by 0.5% of outcomes. The
# Clopper-Pearson interval shows how much the rate itself is pinned down by
# 172 triangles.
# ------------------------------------------------------------------------------
tail_cal <- main[, {
  k  <- sum(exceed995)
  ci <- binom.test(k, .N, 0.005)$conf.int
  .(n = .N, n_exceed_995 = k, exceed_995 = k / .N,
    exceed_lo = ci[1], exceed_hi = ci[2],
    times_nominal = (k / .N) / 0.005,
    exceed_95 = mean(pct > 0.95))
}, by = model]
cat("\nTail calibration:\n")
print(tail_cal[, lapply(.SD, function(x) if (is.double(x)) round(x, 4) else x)])
save_table(tail_cal, "04_tail_calibration")

# ------------------------------------------------------------------------------
# 3. REPLICATION OF MEYERS (2015) ON HIS 200 TRIANGLES
# His KS statistic is max |p_(i) - i/n| in percent (see ks_meyers()). Values
# transcribed from Figures 2-4 of the monograph.
# ------------------------------------------------------------------------------
pub <- data.table(
  model = rep(c("Mack incurred", "Mack paid", "ODP paid"), each = 5),
  line  = rep(c("CA", "PA", "WC", "OL", "All"), 3),
  meyers_d = c(16.4, 14.7, 25.0, 14.2, 15.4,
               20.3, 48.3, 32.4, 10.6, 24.9,
               19.1, 48.3, 29.1, 12.1, 25.6))

m200 <- res[meyers200 == TRUE]
rep_line <- m200[, .(n = .N, this_repo_d = 100 * ks_meyers(pct)), by = .(model, line)]
rep_all  <- m200[, .(n = .N, this_repo_d = 100 * ks_meyers(pct), line = "All"), by = model]
replic <- merge(pub, rbind(rep_line, rep_all, use.names = TRUE), by = c("model", "line"))
replic[, crit_5pct := 136 / sqrt(n)]
replic[, `:=`(meyers_rejects = meyers_d > 136 / sqrt(ifelse(line == "All", 200, 50)),
              this_repo_rejects = this_repo_d > crit_5pct)]
replic[, line := factor(line, c("CA", "PA", "WC", "OL", "All"))]
setorder(replic, model, line)
cat("\n=========== REPLICATION: MEYERS (2015) 200 TRIANGLES ===========\n")
print(replic[, .(model, line, n, meyers_d, this_repo_d = round(this_repo_d, 1),
                 meyers_rejects, this_repo_rejects)])
cat(sprintf("Same reject/accept verdict in %d of %d cells\n",
            sum(replic$meyers_rejects == replic$this_repo_rejects), nrow(replic)))
save_table(replic, "04_meyers_replication")

# ------------------------------------------------------------------------------
# 4. DO RESERVE SURPRISES MOVE TOGETHER ACROSS LINES OF ONE COMPANY?
# Input for the copula in 05_capital.R. For every company with two or more
# lines in the screened sample, pair the outcome percentiles of its lines and
# measure rank dependence across companies. Kendall's tau converts to the
# correlation parameter of a Gaussian or t copula via rho = sin(pi * tau / 2).
# This sees company-specific co-movement (one claims department, one reserving
# philosophy) but not a shock that hits every insurer in the same year, so it
# is a floor for the dependence a capital model should carry, not an estimate
# of all of it.
# ------------------------------------------------------------------------------
pairs_for <- function(d) {
  d <- d[order(line)]
  if (nrow(d) < 2) return(NULL)
  cmb <- combn(nrow(d), 2)
  data.table(pair = paste(d$line[cmb[1, ]], d$line[cmb[2, ]], sep = "-"),
             u1 = d$pct[cmb[1, ]], u2 = d$pct[cmb[2, ]])
}
dep <- rbindlist(lapply(c("ODP paid", "Mack incurred"), function(mdl) {
  pr <- main[model == mdl, pairs_for(.SD), by = grcode]
  tau <- cor(pr$u1, pr$u2, method = "kendall")
  # bootstrap over companies, since one company contributes several pairs
  set.seed(42)
  comps <- unique(pr$grcode)
  taus <- replicate(2000, {
    s <- sample(comps, replace = TRUE)
    b <- pr[data.table(grcode = s), on = "grcode", allow.cartesian = TRUE]
    cor(b$u1, b$u2, method = "kendall")
  })
  data.table(model = mdl, n_companies = length(comps), n_pairs = nrow(pr),
             kendall_tau = tau, tau_lo = quantile(taus, 0.025), tau_hi = quantile(taus, 0.975),
             rho = sin(pi * tau / 2),
             rho_lo = sin(pi * quantile(taus, 0.025) / 2),
             rho_hi = sin(pi * quantile(taus, 0.975) / 2))
}))
cat("\n=========== CROSS-LINE DEPENDENCE OF OUTCOME PERCENTILES ===========\n")
print(dep[, lapply(.SD, function(x) if (is.double(x)) round(x, 3) else x)])
save_table(dep, "04_cross_line_dependence")

# ------------------------------------------------------------------------------
# 5. WHY THE PAID MODELS OVER-RESERVE: A DIAGNOSTIC KNOWN AT YEAR END 1997
# Chain ladder on paid losses assumes each accident year pays out on the same
# schedule as the older ones. If claims are being settled faster, the recent
# years have already paid more of their ultimate than the history implies, and
# multiplying them by the old development factors overshoots.
#
# The classic tell (Berquist & Sherman 1977) is the paid-to-incurred ratio at a
# given age drifting upward across accident years. For lags 1-8:
#   shift_j = (paid / incurred on the 1997 diagonal at lag j)
#             - (mean paid / incurred of the earlier accident years at lag j)
# and the P/I shift is the average of shift_j. Upper triangle only -- nothing
# here uses the outcomes, so an actuary could have computed it in 1998.
# Positive = more of what has been reported is already paid than used to be:
# faster settlement, or case reserves set weaker than before.
# ------------------------------------------------------------------------------
tri <- readRDS(file.path(DIR_PROC, "triangles.rds"))
tri <- tri[sapply(tri, `[[`, "selected")]
pi_shift <- rbindlist(lapply(tri, function(t) {
  r <- upper_only(t$paid) / upper_only(t$incurred)
  s <- sapply(1:8, function(j) {
    k <- which(!is.na(r[, j])); last <- max(k)
    r[last, j] - mean(r[setdiff(k, last), j])
  })
  data.table(line = t$line, grcode = t$grcode, pi_shift = mean(s))
}))
det <- readRDS(file.path(DIR_PROC, "deterministic.rds"))
diag <- merge(pi_shift,
              dcast(det[method == "Chain ladder"], line + grcode ~ basis, value.var = "err_reserve"),
              by = c("line", "grcode"))
setnames(diag, c("paid", "incurred"), c("cl_paid_error", "cl_incurred_error"))
diag[, tercile := cut(pi_shift, quantile(pi_shift, 0:3 / 3), include.lowest = TRUE,
                      labels = c("Lowest third", "Middle third", "Highest third"))]
save_table(diag, "04_pi_shift_by_triangle")

cat(sprintf("\n=========== PAID-TO-INCURRED SHIFT vs PAID CHAIN LADDER ERROR ===========\n"))
cat(sprintf("Spearman correlation, P/I shift vs paid CL reserve error: %.2f (n = %d)\n",
            cor(diag$pi_shift, diag$cl_paid_error, method = "spearman"), nrow(diag)))

by_t <- merge(main, diag[, .(line, grcode, tercile)],
              by = c("line", "grcode"))
terc <- by_t[, {
  ks <- suppressWarnings(ks.test(pct, "punif"))
  .(n = .N, mean_pct = mean(pct), ks_d = unname(ks$statistic), ks_p = ks$p.value)
}, by = .(tercile, model)]
terc <- merge(terc,
              diag[, .(pi_shift_median = median(pi_shift),
                       cl_paid_median_error = median(cl_paid_error),
                       cl_incurred_median_error = median(cl_incurred_error)), by = tercile],
              by = "tercile")
setorder(terc, tercile, model)
print(terc[, lapply(.SD, function(x) if (is.double(x)) signif(x, 3) else x)])
save_table(terc, "04_pi_shift_terciles")

# ------------------------------------------------------------------------------
# 5b. ONE TRIANGLE AS THE EXAMPLE
# Of the five triangles with the highest P/I shift, the one with the largest
# reserve: the speed-up is visible by accident year in the paid-to-incurred
# ratio at lag 1 (known at year end 1997) and, with hindsight, in the share of
# the lag-10 amount paid in the first year.
# ------------------------------------------------------------------------------
top5 <- diag[order(-pi_shift)][1:5]
top5[, cl_paid_reserve := sapply(paste(line, grcode, sep = "_"), function(k) {
  up <- upper_only(tri[[k]]$paid); sum(cl_ultimate(up)) - sum(latest_diag(up)) })]
ex_key <- top5[which.max(cl_paid_reserve), paste(line, grcode, sep = "_")]
ex <- tri[[ex_key]]
P <- upper_only(ex$paid); I <- upper_only(ex$incurred)
ex_ay <- data.table(ay = FIRST_AY:(FIRST_AY + N_AY - 1L),
                    paid_to_incurred_lag1 = P[, 1] / I[, 1],
                    first_year_share_of_lag10_paid = ex$paid[, 1] / ex$paid[, N_AY])
ex_sum <- data.table(
  line = ex$line, grcode = ex$grcode, grname = ex$grname,
  pi_shift = diag[line == ex$line & grcode == ex$grcode, pi_shift],
  pi_shift_rank = diag[, frank(-pi_shift)][diag$line == ex$line & diag$grcode == ex$grcode],
  n_triangles = nrow(diag),
  cl_paid_reserve = sum(cl_ultimate(P)) - sum(latest_diag(P)),
  cl_incurred_reserve = sum(cl_ultimate(I)) - sum(latest_diag(P)),
  actual_unpaid = sum(ex$paid[, N_AY]) - sum(latest_diag(P)),
  mack_paid_pct = main[model == "Mack paid" & line == ex$line & grcode == ex$grcode, pct],
  mack_incurred_pct = main[model == "Mack incurred" & line == ex$line & grcode == ex$grcode, pct])
cat(sprintf("\nExample triangle: %s %s\n", ex$grname, ex$line))
print(ex_sum)
print(ex_ay[, lapply(.SD, function(x) if (is.double(x)) round(x, 3) else x)])
save_table(ex_ay, "04_example_triangle_by_ay")
save_table(ex_sum, "04_example_triangle_summary")

# ------------------------------------------------------------------------------
# 6. ROBUSTNESS: THE BOOTSTRAP'S OWN INSTABILITY
# On triangles with small early payments (typical of other liability, product
# liability and med mal), many bootstrap pseudo-triangles go non-positive and
# their projections explode, stretching the simulated right tail. Re-run the
# verdict on the triangles where no pseudo-triangle did, with Mack paid on the
# same triangles alongside, to show the conclusion does not depend on them.
# ------------------------------------------------------------------------------
stab <- main[model == "ODP paid", .(line, grcode, stable = invalid_share == 0)]
rob <- merge(main[model != "Mack incurred"], stab, by = c("line", "grcode"))
robust <- rob[, c(.(subset = fifelse(stable, "Stable bootstrap", "Some invalid pseudo-triangles")),
                  uniformity(pct, exceed995)[c("n", "ks_d", "ks_p", "mean_pct", "n_exceed_995")]),
              by = .(model, stable)][, stable := NULL]
robust <- rbind(robust, main[model != "Mack incurred",
                             c(.(subset = "All screened"),
                               uniformity(pct, exceed995)[c("n", "ks_d", "ks_p", "mean_pct", "n_exceed_995")]),
                             by = model])
setorder(robust, model, subset)
cat("\n=========== ROBUSTNESS: BOOTSTRAP STABILITY ===========\n")
cat(sprintf("Screened triangles with any invalid pseudo-triangle: %d of %d (by line: %s)\n",
            sum(!stab$stable), nrow(stab),
            stab[stable == FALSE, .N, by = line][, paste(line, N, sep = " ", collapse = ", ")]))
print(robust[, lapply(.SD, function(x) if (is.double(x)) signif(x, 3) else x)])
save_table(robust, "04_bootstrap_stability")
