# ==============================================================================
# 00_common.R -- shared paths, line metadata, triangle and risk-measure helpers
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# --- paths --------------------------------------------------------------------
# All scripts are meant to be run from the repository root:
#   Rscript R/01_prepare_data.R
PROJ <- getwd()
if (!dir.exists(file.path(PROJ, "R")) && dir.exists("../R")) PROJ <- normalizePath("..")
stopifnot(dir.exists(file.path(PROJ, "R")))

DIR_RAW  <- file.path(PROJ, "data", "raw")
DIR_REF  <- file.path(PROJ, "data", "reference")
DIR_PROC <- file.path(PROJ, "data", "processed")
DIR_FIG  <- file.path(PROJ, "outputs", "figures")
DIR_TAB  <- file.path(PROJ, "outputs", "tables")
for (d in c(DIR_PROC, DIR_FIG, DIR_TAB)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# --- the six Schedule P lines -------------------------------------------------
# Two-letter codes follow Meyers (2015) for the four lines he used.
LINES <- data.table(
  line = c("PA", "CA", "WC", "OL", "PL", "MM"),
  file = c("ppauto_pos", "comauto_pos", "wkcomp_pos", "othliab_pos",
           "prodliab_pos", "medmal_pos"),
  name = c("Private passenger auto", "Commercial auto", "Workers' compensation",
           "Other liability", "Product liability", "Medical malpractice"))

N_AY      <- 10L     # accident years 1988-1997, development lags 1-10
FIRST_AY  <- 1988L
EVAL_YEAR <- 1997L   # the upper triangle is what was known at year end 1997

# Screening limits (see 01_prepare_data.R for why these values)
CV_PREM_MAX     <- 0.40
CV_NET_DIR_MAX  <- 0.125

N_BOOT <- 10000L     # bootstrap ODP simulations per triangle in the backtest

#' Write a table to outputs/tables.
save_table <- function(dt, name) {
  fwrite(dt, file.path(DIR_TAB, paste0(name, ".csv")))
  invisible(dt)
}

# --- triangle helpers ---------------------------------------------------------
# A triangle is a 10 x 10 numeric matrix, rows = accident year, columns =
# development lag. "Full" matrices hold the upper AND lower triangle (the CAS
# database has both); "upper" matrices have the unknown future set to NA.

#' Pivot one company-line from long format into a full 10 x 10 matrix.
full_matrix <- function(d, col) {
  m <- matrix(NA_real_, N_AY, N_AY,
              dimnames = list(origin = FIRST_AY:(FIRST_AY + N_AY - 1L), dev = 1:N_AY))
  m[cbind(d$ay - FIRST_AY + 1L, d$lag)] <- d[[col]]
  m
}

#' Blank out everything after the year-end-1997 diagonal.
upper_only <- function(m) {
  i <- row(m); j <- col(m)
  m[i + j > N_AY + 1L] <- NA
  m
}

#' Latest diagonal of an upper triangle (cumulative amount at year end 1997).
latest_diag <- function(tri) {
  n <- nrow(tri)
  tri[cbind(1:n, n:1)]
}

#' Volume-weighted chain ladder link ratios f_j = sum_i C[i,j+1] / sum_i C[i,j],
#' summed over the accident years where both cells are known.
cl_factors <- function(tri) {
  n <- ncol(tri)
  sapply(1:(n - 1L), function(j) {
    k <- 1:(n - j)                      # rows with lag j+1 observed
    sum(tri[k, j + 1L]) / sum(tri[k, j])
  })
}

#' Cumulative development factor to lag 10 for each accident year.
cdf_to_ult <- function(f) {
  n <- length(f) + 1L
  # accident year i sits at lag n-i+1; it still needs factors (n-i+1) .. (n-1)
  sapply(1:n, function(i) {
    lag_now <- n - i + 1L
    if (lag_now >= n) 1 else prod(f[lag_now:(n - 1L)])
  })
}

#' Chain ladder ultimates (to lag 10, no tail factor).
cl_ultimate <- function(tri) latest_diag(tri) * cdf_to_ult(cl_factors(tri))

#' Bornhuetter-Ferguson ultimates: latest + (1 - 1/CDF) * ELR * premium.
bf_ultimate <- function(tri, premium, elr) {
  cdf <- cdf_to_ult(cl_factors(tri))
  latest_diag(tri) + (1 - 1 / cdf) * elr * premium
}

# --- predictive distributions and risk measures -------------------------------

#' Lognormal parameters with a given mean and standard deviation.
lnorm_par <- function(mean, sd) {
  s2 <- log(1 + (sd / mean)^2)
  list(meanlog = log(mean) - s2 / 2, sdlog = sqrt(s2))
}

#' Percentile of an outcome under a lognormal with matched mean and sd. This is
#' how Meyers (2015), following Mack, turns Mack's mean and standard error into
#' a full predictive distribution.
pct_lnorm <- function(x, mean, sd) {
  p <- lnorm_par(mean, sd)
  plnorm(x, p$meanlog, p$sdlog)
}

qnt_lnorm <- function(q, mean, sd) {
  p <- lnorm_par(mean, sd)
  qlnorm(q, p$meanlog, p$sdlog)
}

#' Percentile of an outcome within a simulated distribution (mid-rank on ties).
pct_sim <- function(x, sims) mean(sims < x) + 0.5 * mean(sims == x)

#' Share of bootstrap pseudo-triangles with a non-positive cumulative value.
#' A large negative residual resampled into a small early cell can push a
#' pseudo cumulative amount to zero or below; the re-fitted link ratios are
#' then undefined or enormous and that simulation's projection explodes.
#' sim_inc is BootChainLadder's simClaims: origin x dev x simulation,
#' incremental, upper triangle only.
pseudo_invalid_share <- function(sim_inc) {
  cum <- sim_inc
  for (j in 2:dim(cum)[2]) cum[, j, ] <- cum[, j - 1L, ] + cum[, j, ]
  mean(apply(cum <= 0, 3, any, na.rm = TRUE))
}

#' Meyers' (2015) version of the KS statistic: max |p_(i) - i/n|. The textbook
#' statistic (ks.test) also checks p_(i) - (i-1)/n, so the two differ by at most
#' 1/n. Used only to compare like with like against his published values.
ks_meyers <- function(p) {
  n <- length(p)
  max(abs(sort(p) - (1:n) / n))
}

#' Value at Risk and Expected Shortfall of a simulated loss distribution.
#' ES is the mean of the worst (1 - p) share of outcomes.
var_es <- function(x, p = 0.995) {
  xs <- sort(x)
  n  <- length(xs)
  k  <- ceiling(p * n)
  c(VaR = xs[k], ES = mean(xs[k:n]))
}

#' Inverse empirical CDF: map uniforms onto a simulated marginal distribution.
inv_ecdf <- function(sims, u) {
  xs <- sort(sims)
  xs[pmax(1L, ceiling(u * length(xs)))]
}

#' Uniforms from a Gaussian (df = Inf) or Student-t copula with exchangeable
#' correlation rho, written out rather than via the `copula` package so the
#' construction is visible. One-factor form:
#'   x_j = sqrt(rho) * f + sqrt(1 - rho) * e_j   gives corr(x_i, x_j) = rho,
#' and dividing every x_j by the same sqrt(chi2_df / df) makes x multivariate t
#' (the shared divisor is what creates tail dependence). The normals and the
#' chi-square are arguments so every rho can reuse the same draws: common
#' random numbers, so a difference between two rho values is not noise.
#'   f: n-vector N(0,1); e: n x d matrix N(0,1); chi: n-vector chi2(df)
copula_u <- function(f, e, rho, df = Inf, chi = NULL) {
  x <- sqrt(rho) * f + sqrt(1 - rho) * e
  if (is.finite(df)) pt(x / sqrt(chi / df), df) else pnorm(x)
}
