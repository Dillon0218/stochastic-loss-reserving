# ==============================================================================
# 05_capital.R
#
# Reserve risk capital for one multi-line insurer at 99.5%, VaR and ES, with
# lines aggregated through Gaussian and Student-t copulas.
#
# Company, chosen mechanically:
#   1. three or more lines in the screened sample;
#   2. a stable bootstrap on every line: no pseudo-triangle went non-positive
#      in the 10,000-simulation run of 03_stochastic.R. A capital model is not
#      signed off on a bootstrap that throws exploding simulations, and the ES
#      would be measuring those explosions (Federal Insurance's product
#      liability: 58% of pseudo-triangles invalid);
#   3. of those, the least concentrated reserve (smallest largest-line share of
#      the chain ladder paid reserve). Aggregation only has something to say
#      when no single line dominates; for State Farm, personal auto is 87%.
#
# Marginals: bootstrap ODP on paid, 100,000 simulations per line (ten
# independent runs of 10,000 pooled, which is the same thing as one run of
# 100,000 because every simulation resamples independently).
#
# Aggregation by reordering (Iman & Conover 1982, the way most capital models
# do it): draw 100,000 points from the copula, and pair each line's simulations
# in the rank order of the copula draws. Each line keeps exactly its own
# simulated values in every scenario -- only which values occur together
# changes -- so standalone capital is identical across scenarios and, at
# correlation 1, aggregate capital equals the sum of standalone exactly.
#
# Two views of the same simulations
#   Ultimate view  -- all future payments to lag 10. Capital = VaR(reserve) -
#                     mean(reserve). What this project's backtest scores.
#   One-year view  -- Solvency II: next year's payments plus the reserve the
#                     actuary would re-estimate at the end of next year, given
#                     the simulated next diagonal ("actuary in the box";
#                     ChainLadder's NYCost). Only one year of surprises can
#                     emerge, so it is narrower. Merz-Wuthrich, the analytic
#                     one-year counterpart of Mack, is shown as a check.
#
# Capital here is the amount above the mean reserve; no discounting, no risk
# margin, no allowance for the model error the backtest found (see README).
# ==============================================================================

source("R/00_common.R")
suppressPackageStartupMessages(library(ChainLadder))

tri   <- readRDS(file.path(DIR_PROC, "triangles.rds"))
multi <- fread(file.path(DIR_TAB, "01_multiline_companies.csv"))

P        <- 0.995
N_BATCH  <- 10L                  # x N_BOOT = 100,000 simulations per line
N_COP    <- N_BATCH * N_BOOT     # one copula draw per simulation (reordering)
T_DF     <- 4                    # Student-t degrees of freedom
RHO_GRID <- seq(0, 1, by = 0.05)
RHO_SHOW <- c(0, 0.25, 0.5, 0.75, 1)

# ------------------------------------------------------------------------------
# 1. CHOOSE THE COMPANY
# ------------------------------------------------------------------------------
odp  <- readRDS(file.path(DIR_PROC, "predictive.rds"))[selected == TRUE & model == "ODP paid"]
cand <- multi[n_lines >= 3]
mix <- rbindlist(lapply(cand$grcode, function(g) {
  ks <- names(tri)[sapply(tri, function(t) t$grcode == g && t$selected)]
  r  <- sapply(ks, function(k) { up <- upper_only(tri[[k]]$paid)
                                 sum(cl_ultimate(up)) - sum(latest_diag(up)) })
  inv <- odp[grcode == g, invalid_share]
  data.table(grcode = g, grname = tri[[ks[1]]]$grname, n_lines = length(ks),
             lines = paste(sapply(ks, function(k) tri[[k]]$line), collapse = " "),
             cl_paid_reserve = sum(r), largest_line_share = max(r) / sum(r),
             worst_invalid_share = if (length(inv) == length(ks)) max(inv) else NA_real_)
}))
mix[, stable := !is.na(worst_invalid_share) & worst_invalid_share == 0]
setorder(mix, -stable, largest_line_share)
cat("Companies with three or more lines in the screened sample:\n")
print(mix[, .(grcode, grname, lines, largest_line_share = round(largest_line_share, 3),
              worst_invalid_share = round(worst_invalid_share, 3), stable)])
save_table(mix, "05_company_choice")

G <- mix[stable == TRUE]$grcode[1]
keys <- names(tri)[sapply(tri, function(t) t$grcode == G && t$selected)]
keys <- keys[order(match(sapply(keys, function(k) tri[[k]]$line), LINES$line))]
lines <- sapply(keys, function(k) tri[[k]]$line)
cat(sprintf("\nCompany: %s (NAIC group %d), lines %s\n",
            tri[[keys[1]]]$grname, G, paste(lines, collapse = ", ")))

# ------------------------------------------------------------------------------
# 2. STANDALONE DISTRIBUTIONS PER LINE
# ------------------------------------------------------------------------------
quiet_info <- function(expr) withCallingHandlers(expr, warning = function(w) {
  if (grepl("^Information: essentially no variation", conditionMessage(w)))
    invokeRestart("muffleWarning")
})

sims <- list()
standalone <- rbindlist(lapply(seq_along(keys), function(i) {
  t  <- tri[[keys[i]]]
  up <- upper_only(t$paid)
  reserve <- one_year <- numeric(0)
  for (b in seq_len(N_BATCH)) {
    bt <- BootChainLadder(up, R = N_BOOT, process.distr = "gamma",
                          seed = 20261101L + 100L * i + b)
    reserve  <- c(reserve, bt$IBNR.Totals)
    one_year <- c(one_year, apply(bt$NYCost.ByOrigin, 3, sum))
  }
  sims[[t$line]] <<- list(reserve = reserve, one_year = one_year)

  m   <- quiet_info(MackChainLadder(up, est.sigma = "log-linear"))
  mw  <- quiet_info(CDR(m))
  ult <- var_es(reserve, P); one <- var_es(one_year, P)
  actual <- sum(t$paid[, N_AY]) - sum(latest_diag(up))
  data.table(
    line = t$line,
    cl_reserve = sum(cl_ultimate(up)) - sum(latest_diag(up)),
    mean = mean(reserve), sd = sd(reserve), cv = sd(reserve) / mean(reserve),
    var995 = ult[["VaR"]], es995 = ult[["ES"]],
    cap_var = ult[["VaR"]] - mean(reserve), cap_es = ult[["ES"]] - mean(reserve),
    one_year_mean = mean(one_year), one_year_sd = sd(one_year),
    one_year_cap_var = one[["VaR"]] - mean(one_year),
    one_year_cap_es  = one[["ES"]] - mean(one_year),
    mack_se = as.numeric(m$Total.Mack.S.E),
    merz_wuthrich_se = mw["Total", "CDR(1)S.E."],
    actual_unpaid = actual,
    actual_pct = pct_sim(actual, reserve),
    # Monte Carlo precision: VaR from two independent halves of the simulations
    var_half_diff = abs(diff(c(var_es(reserve[1:(length(reserve) / 2)], P)[["VaR"]],
                               var_es(reserve[-(1:(length(reserve) / 2))], P)[["VaR"]]))) / ult[["VaR"]])
}))
standalone[, `:=`(one_year_to_ult_var = one_year_cap_var / cap_var,
                  one_year_to_ult_sd  = one_year_sd / sd,
                  mw_to_mack          = merz_wuthrich_se / mack_se)]

cat("\n=========== STANDALONE RESERVE RISK, ULTIMATE VIEW (USD 000) ===========\n")
print(standalone[, .(line, cl_reserve = round(cl_reserve), mean = round(mean), cv = round(cv, 3),
                     var995 = round(var995), es995 = round(es995),
                     cap_var = round(cap_var), cap_es = round(cap_es),
                     actual_unpaid = round(actual_unpaid), actual_pct = round(actual_pct, 4))])
cat("\nOne-year vs ultimate view:\n")
print(standalone[, .(line, cap_var = round(cap_var), one_year_cap_var = round(one_year_cap_var),
                     one_year_to_ult_var = round(one_year_to_ult_var, 3),
                     one_year_to_ult_sd = round(one_year_to_ult_sd, 3),
                     mw_to_mack = round(mw_to_mack, 3))])
cat(sprintf("\nMonte Carlo check: standalone VaR from two halves of the simulations differs by at most %.1f%%\n",
            100 * max(standalone$var_half_diff)))
save_table(standalone, "05_standalone")

# ------------------------------------------------------------------------------
# 3. AGGREGATION ACROSS LINES
# Every line's simulated distribution is kept exactly as simulated (the
# marginal); the copula only decides which simulated outcomes happen together.
# ------------------------------------------------------------------------------
set.seed(20261101L)
d   <- length(lines)
f   <- rnorm(N_COP)
e   <- matrix(rnorm(N_COP * d), N_COP, d)
chi <- rchisq(N_COP, T_DF)

aggregate_view <- function(view) {
  marg <- lapply(lines, function(l) sort(sims[[l]][[view]]))
  st   <- sapply(marg, function(x) var_es(x, P) - mean(x))     # standalone, fixed
  rbindlist(lapply(c("Gaussian", "t"), function(cop) {
    rbindlist(lapply(RHO_GRID, function(rho) {
      u <- copula_u(f, e, rho, df = if (cop == "t") T_DF else Inf, chi = chi)
      # reorder: the simulation with the k-th smallest value of line j goes
      # where the copula put its k-th smallest draw for line j
      x <- sapply(seq_len(d), function(j) marg[[j]][rank(u[, j], ties.method = "first")])
      s <- rowSums(x)
      ag <- var_es(s, P) - mean(s)
      data.table(view = view, copula = cop, rho = rho,
                 mean = mean(s),
                 cap_var = ag[["VaR"]], cap_es = ag[["ES"]],
                 sum_standalone_var = sum(st["VaR", ]), sum_standalone_es = sum(st["ES", ]))
    }))
  }))
}

agg <- rbind(aggregate_view("reserve"), aggregate_view("one_year"))
agg[, `:=`(div_var = 1 - cap_var / sum_standalone_var,
           div_es  = 1 - cap_es / sum_standalone_es)]
agg[, view := fifelse(view == "reserve", "Ultimate", "One-year")]

cat("\n=========== AGGREGATE RESERVE RISK CAPITAL AT 99.5% (USD 000) ===========\n")
print(agg[rho %in% RHO_SHOW, .(view, copula, rho, cap_var = round(cap_var), cap_es = round(cap_es),
                              div_var = round(div_var, 3), div_es = round(div_es, 3))])

# t minus Gaussian: what tail dependence adds at the same correlation
td <- dcast(agg[view == "Ultimate"], rho ~ copula, value.var = "cap_var")
td[, t_over_gauss := t / Gaussian - 1]
cat("\nUltimate view, VaR capital: t(4) copula relative to Gaussian at the same rho:\n")
print(td[rho %in% RHO_SHOW, .(rho, Gaussian = round(Gaussian), t = round(t),
                              t_over_gauss = round(t_over_gauss, 3))])
save_table(agg, "05_aggregation")

# Where the realised run-off landed in the aggregate distribution, as a check
# on the whole chain: the actual total unpaid against the Gaussian rho = 0.25
# aggregate (the median of the simulated total is reported alongside).
u  <- copula_u(f, e, 0.25)
s  <- rowSums(sapply(seq_len(d), function(j)
  sort(sims[[lines[j]]]$reserve)[rank(u[, j], ties.method = "first")]))
outcome <- data.table(actual_total_unpaid = sum(standalone$actual_unpaid),
                      median_simulated = median(s), mean_simulated = mean(s),
                      actual_pct = pct_sim(sum(standalone$actual_unpaid), s))
cat("\nActual run-off against the aggregate distribution:\n")
print(outcome)
save_table(outcome, "05_actual_vs_aggregate")

saveRDS(list(grcode = G, lines = lines, sims = sims), file.path(DIR_PROC, "capital_sims.rds"))
