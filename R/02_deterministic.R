# ==============================================================================
# 02_deterministic.R
#
# Chain ladder and Bornhuetter-Ferguson point estimates on the screened sample,
# scored against what each company actually went on to pay or report.
#
#   Chain ladder (CL)  volume-weighted link ratios, projected to lag 10.
#   BF                 latest + (1 - 1/CDF) x ELR x net earned premium. The
#                      a priori ELR is the premium-weighted CL ultimate loss
#                      ratio of the five most mature accident years (1988-92),
#                      i.e. "premium and early experience". It is not adjusted
#                      for rate changes or trend, which Schedule P does not
#                      report; see the README limitations.
#
# Both methods run on the paid and on the incurred triangle. Each estimate is
# compared with the outcome on its own basis (paid with paid at lag 10, incurred
# with incurred at lag 10). Errors are expressed as a share of the actual unpaid
# amount at year end 1997 -- the reserve the company really needed -- so paid
# and incurred errors sit on the same scale.
# ==============================================================================

source("R/00_common.R")
suppressPackageStartupMessages(library(ChainLadder))

tri <- readRDS(file.path(DIR_PROC, "triangles.rds"))
tri <- tri[sapply(tri, `[[`, "selected")]
cat(sprintf("Screened sample: %d triangles\n", length(tri)))

MATURE <- 1:5   # accident years 1988-92: at least six years of development

one_triangle <- function(t) {
  rbindlist(lapply(c("paid", "incurred"), function(basis) {
    full <- t[[basis]]
    up   <- upper_only(full)
    ult_cl <- cl_ultimate(up)
    elr    <- sum(ult_cl[MATURE]) / sum(t$premium[MATURE])
    ult_bf <- bf_ultimate(up, t$premium, elr)
    actual <- full[, N_AY]
    # the reserve that was really needed: paid at lag 10 less paid to date
    unpaid <- sum(t$paid[, N_AY]) - sum(latest_diag(upper_only(t$paid)))
    data.table(line = t$line, grcode = t$grcode, basis = basis,
               method = c("Chain ladder", "Bornhuetter-Ferguson"),
               elr = elr,
               est_ult = c(sum(ult_cl), sum(ult_bf)),
               act_ult = sum(actual),
               actual_unpaid = unpaid,
               # most recent accident year (1997, one year of development)
               est_ay97 = c(ult_cl[N_AY], ult_bf[N_AY]),
               act_ay97 = actual[N_AY])
  }))
}

det <- rbindlist(lapply(tri, one_triangle))
det[, err_reserve := (est_ult - act_ult) / actual_unpaid]
det[, err_ay97    := (est_ay97 - act_ay97) / act_ay97]

# Cross-check the hand-written chain ladder against the ChainLadder package on
# every triangle: same ultimates to rounding, or stop.
cl_diff <- max(sapply(tri, function(t) {
  up  <- upper_only(t$paid)
  ref <- predict(chainladder(up))[, N_AY]
  max(abs(cl_ultimate(up) - ref) / ref)
}))
cat(sprintf("Max relative difference vs ChainLadder::chainladder ultimates: %.2e\n", cl_diff))
stopifnot(cl_diff < 1e-10)

# Medians, not means: a company with a 50K reserve and a 200% miss would
# otherwise outweigh forty well-behaved ones. The dollar-weighted error is the
# portfolio view -- what the sample as a whole was over- or under-reserved by.
err_stats <- function(d) d[, .(n = .N,
                               median_error     = median(err_reserve),
                               median_abs_error = median(abs(err_reserve)),
                               dollar_weighted_error = sum(est_ult - act_ult) / sum(actual_unpaid),
                               share_over       = mean(err_reserve > 0))]
rnd <- function(d) d[, lapply(.SD, function(x) if (is.double(x)) round(x, 3) else x)]

cat("\n=========== RESERVE ERROR vs ACTUAL RUN-OFF (share of actual unpaid) ===========\n")
summ <- det[, err_stats(.SD), by = .(basis, method)]
print(rnd(summ))
save_table(summ, "02_reserve_error_summary")

by_line <- det[, err_stats(.SD), by = .(line, basis, method)]
by_line[, line := factor(line, LINES$line)]
setorder(by_line, line, basis, method)
cat("\nMedian error by line:\n")
print(dcast(rnd(by_line), line ~ basis + method, value.var = "median_error"))
save_table(by_line, "02_reserve_error_by_line")

cat("\n=========== MOST RECENT ACCIDENT YEAR (1997) ===========\n")
ay97 <- det[, .(median_error = median(err_ay97),
                median_abs_error = median(abs(err_ay97))), by = .(basis, method)]
print(rnd(ay97))
save_table(ay97, "02_ay1997_error_summary")

# ------------------------------------------------------------------------------
# Why BF misses where it does: the a priori ELR comes from 1988-92, so BF is
# only right if later accident years run at the same loss ratio. Hindsight loss
# ratios (incurred at lag 10 / net earned premium, pooled over the sample) show
# whether they did. Diagnostic only -- nothing here feeds back into a model.
# ------------------------------------------------------------------------------
lr <- rbindlist(lapply(tri, function(t)
  data.table(line = t$line, ay = FIRST_AY:(FIRST_AY + N_AY - 1L),
             prem = t$premium, ult_inc = t$incurred[, N_AY])))
lr_shift <- lr[, .(lr_1988_92 = sum(ult_inc[ay <= 1992]) / sum(prem[ay <= 1992]),
                   lr_1993_97 = sum(ult_inc[ay >= 1993]) / sum(prem[ay >= 1993])), by = line]
lr_shift <- merge(lr_shift,
                  by_line[basis == "paid", .(line, method, median_error)] |>
                    dcast(line ~ method, value.var = "median_error"),
                  by = "line")
setnames(lr_shift, c("Bornhuetter-Ferguson", "Chain ladder"),
         c("bf_paid_median_error", "cl_paid_median_error"))
lr_shift[, line := factor(line, LINES$line)]
setorder(lr_shift, line)
cat("\nHindsight loss ratio shift vs BF paid error:\n")
print(rnd(lr_shift))
save_table(lr_shift, "02_loss_ratio_shift")

save_table(det, "02_deterministic_by_triangle")
saveRDS(det, file.path(DIR_PROC, "deterministic.rds"))
