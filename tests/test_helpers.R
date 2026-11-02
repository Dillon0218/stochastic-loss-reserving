# ==============================================================================
# tests/test_helpers.R -- checks on the hand-written helpers in R/00_common.R
#   Rscript tests/test_helpers.R      (from the repository root)
#
# The package functions (MackChainLadder, BootChainLadder) are ChainLadder's
# responsibility. What this repo wrote itself -- the chain ladder and BF
# projections, the percentile and KS helpers, VaR/ES, the copula sampler, the
# bootstrap validity check -- is checked against known answers or an
# independent implementation.
# ==============================================================================

source("R/00_common.R")
suppressPackageStartupMessages(library(ChainLadder))
set.seed(1)

fails <- 0
check <- function(label, got, want, tol = 1e-6) {
  ok <- isTRUE(abs(got - want) < tol)
  cat(sprintf("%-58s %12.6f  (expected %10.4f)  %s\n",
              label, got, want, if (ok) "PASS" else "FAIL"))
  if (!ok) fails <<- fails + 1
}

# --- chain ladder vs the ChainLadder package ----------------------------------
for (nm in c("RAA", "GenIns")) {
  tri <- unclass(get(nm)); tri <- matrix(tri, nrow(tri), dimnames = dimnames(tri))
  ref <- predict(chainladder(tri))[, ncol(tri)]
  check(sprintf("%s: chain ladder ultimates vs ChainLadder (max rel diff)", nm),
        max(abs(cl_ultimate(tri) - ref) / ref), 0, tol = 1e-10)
}

# --- Bornhuetter-Ferguson ------------------------------------------------------
tri  <- unclass(RAA); tri <- matrix(tri, nrow(tri))
prem <- rep(20000, nrow(tri))
ult  <- cl_ultimate(tri)
# with each year's a priori equal to its own CL loss ratio, BF must equal CL
check("BF = CL when ELR is each year's CL loss ratio",
      max(abs(bf_ultimate(tri, prem, ult / prem) - ult)), 0, tol = 1e-8)
# the fully developed year has nothing left to develop
check("BF on the oldest year = latest diagonal",
      bf_ultimate(tri, prem, 0.7)[1], latest_diag(tri)[1], tol = 1e-10)

# --- lognormal moment matching and percentiles ---------------------------------
lp <- lnorm_par(1000, 250)
check("lognormal: mean recovered",  exp(lp$meanlog + lp$sdlog^2 / 2), 1000, tol = 1e-8)
check("lognormal: sd recovered",
      sqrt(exp(lp$sdlog^2) - 1) * exp(lp$meanlog + lp$sdlog^2 / 2), 250, tol = 1e-8)
check("lognormal: percentile of the median = 0.5", pct_lnorm(exp(lp$meanlog), 1000, 250), 0.5)
check("lognormal: quantile and percentile invert",
      pct_lnorm(qnt_lnorm(0.995, 1000, 250), 1000, 250), 0.995)
check("simulated percentile: mid-rank on ties", pct_sim(2, c(1, 2, 2, 3)), 0.5)

# --- KS statistic ---------------------------------------------------------------
n <- 200
check("Meyers KS: perfectly spread percentiles -> 0", ks_meyers((1:n) / n), 0)
p <- runif(n)
# Meyers' D drops the (i-1)/n side, so it sits between ks.test's D - 1/n and D
d_ks <- unname(ks.test(p, "punif")$statistic)
check("Meyers KS between ks.test D - 1/n and D",
      as.numeric(ks_meyers(p) <= d_ks + 1e-12 && ks_meyers(p) >= d_ks - 1 / n - 1e-12), 1)

# --- VaR and ES against the normal distribution ----------------------------------
x <- rnorm(2e6, 100, 10)
z <- qnorm(0.995)
ve <- var_es(x, 0.995)
check("VaR 99.5% of N(100, 10^2)", ve[["VaR"]], 100 + 10 * z, tol = 0.15)
check("ES 99.5% of N(100, 10^2)",  ve[["ES"]],  100 + 10 * dnorm(z) / 0.005, tol = 0.15)

# --- copula sampler ----------------------------------------------------------------
m  <- 40000
f  <- rnorm(m); e <- matrix(rnorm(m * 2), m, 2); chi <- rchisq(m, 4)
ug <- copula_u(f, e, 0.5)
ut <- copula_u(f, e, 0.5, df = 4, chi = chi)
# for any elliptical copula, Kendall's tau = (2 / pi) * asin(rho)
tau_want <- 2 / pi * asin(0.5)
check("Gaussian copula: Kendall tau = 2/pi asin(rho)",
      cor(ug[1:8000, 1], ug[1:8000, 2], method = "kendall"), tau_want, tol = 0.02)
check("t copula: Kendall tau = 2/pi asin(rho)",
      cor(ut[1:8000, 1], ut[1:8000, 2], method = "kendall"), tau_want, tol = 0.02)
check("t copula: margins uniform (KS p > 0.01)",
      as.numeric(ks.test(ut[, 1], "punif")$p.value > 0.01), 1)
# tail dependence: joint exceedance of the 99th percentile is likelier under t
joint <- function(u) mean(u[, 1] > 0.99 & u[, 2] > 0.99) / 0.01
check("t copula has more upper-tail dependence than Gaussian",
      as.numeric(joint(ut) > joint(ug)), 1)
check("rho = 1 gives identical uniforms (comonotonic)",
      max(abs(copula_u(f, e, 1)[, 1] - copula_u(f, e, 1)[, 2])), 0)

# --- bootstrap validity check --------------------------------------------------------
sim <- array(NA_real_, c(3, 3, 4))
for (k in 1:4) { sim[1, , k] <- c(10, 5, 1); sim[2, 1:2, k] <- c(8, 4); sim[3, 1, k] <- 6 }
sim[2, 2, 3] <- -9                               # cumulative 8 - 9 = -1 in simulation 3
check("pseudo_invalid_share flags one bad simulation in four",
      pseudo_invalid_share(sim), 0.25)

cat(sprintf("\n%s\n", if (fails == 0) "All checks passed." else sprintf("%d check(s) FAILED.", fails)))
if (fails > 0) quit(status = 1)
