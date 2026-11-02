# ==============================================================================
# 01_prepare_data.R
#
# Load the six CAS Schedule P files, build paid and incurred triangles for every
# company-line, and screen out the ones the models should not be judged on.
#
# Definitions
#   paid      = Schedule P Part 3 cumulative paid loss and DCC, net of reinsurance
#   incurred  = Part 2 incurred minus Part 4 bulk + IBNR, i.e. reported (case)
#               incurred. This is the definition in Meyers (2015, footnote 6);
#               chain-laddering a number that already contains the company's own
#               IBNR estimate would double count.
#   premium   = net earned premium by accident year
#   outcome   = the lag-10 value of each accident year, summed. "Ultimate" in
#               this project always means development lag 10; no tail factor.
#
# Screening (applied in this order, every drop counted in 01_screening.csv)
#   1. Net and direct earned premium positive in all ten accident years.
#   2. Paid and incurred positive in every cell of the upper triangle, so every
#      link ratio is defined, and the paid triangle still develops after 1997
#      (positive chain ladder paid reserve). This removes group 38997, whose
#      losses are all paid in the first year -- Meyers dropped it by hand.
#   3. CV of net earned premium across accident years below 0.40.
#   4. CV of the net-to-direct premium ratio below 0.125.
# Rules 3 and 4 are Meyers' (2015, Appendix A) two signals of an operational
# change: a swing in volume, or a change in reinsurance. He tuned the limits by
# line to get 50 insurers in each of four lines; 0.40 / 0.125 are his limits
# for commercial auto and other liability, rounded (0.399 and 0.390 on premium;
# personal auto used 0.45, workers' comp 0.95 / 0.20), applied here uniformly
# to all six lines and fixed before any model was run.
#
# Meyers' own 200 triangles (Table 12 of the monograph) are flagged separately
# and used in 04_backtest.R to check this implementation against his results.
# ==============================================================================

source("R/00_common.R")

# ------------------------------------------------------------------------------
# 1. LOAD
# ------------------------------------------------------------------------------
read_line <- function(line, file) {
  d <- fread(file.path(DIR_RAW, paste0(file, ".csv")))
  # strip the line-of-business suffix (_B, _C, _D, _h1, _R1, _F2)
  setnames(d, sub("_[A-Za-z0-9]+$", "", names(d)))
  d[, .(line = line, grcode = GRCODE, grname = trimws(GRNAME),
        ay = AccidentYear, lag = DevelopmentLag,
        paid = as.numeric(CumPaidLoss),
        incurred = as.numeric(IncurLoss - BulkLoss),
        prem_net = as.numeric(EarnedPremNet),
        prem_dir = as.numeric(EarnedPremDIR),
        single = Single)]
}

long <- rbindlist(Map(read_line, LINES$line, LINES$file))
long[, upper := ay + lag - 1L <= EVAL_YEAR]
setkey(long, line, grcode, ay, lag)

cat(sprintf("Loaded %s rows: %d company-line combinations across %d lines\n",
            format(nrow(long), big.mark = ","),
            uniqueN(long[, .(line, grcode)]), uniqueN(long$line)))
stopifnot(long[, .N, by = .(line, grcode)][, all(N == N_AY^2)])

# ------------------------------------------------------------------------------
# 2. SCREENING STATISTICS PER COMPANY-LINE
# ------------------------------------------------------------------------------
stats <- long[, {
  p  <- prem_net[lag == 1L]; pd <- prem_dir[lag == 1L]
  r  <- p / pd
  up <- upper_only(full_matrix(.SD, "paid"))
  .(grname    = grname[1],
    prem_min  = min(p), prem_dir_min = min(pd),
    prem_tot  = sum(p),
    cv_prem   = sd(p) / mean(p),
    cv_netdir = sd(r) / mean(r),
    paid_upper_min = min(paid[upper]),
    inc_upper_min  = min(incurred[upper]),
    # computed only where every cell is positive, otherwise link ratios blow up
    cl_paid_reserve = if (min(paid[upper]) > 0) sum(cl_ultimate(up)) - sum(latest_diag(up))
                      else NA_real_)
}, by = .(line, grcode)]

stats[, `:=`(
  r1_premium  = prem_min > 0 & prem_dir_min > 0,
  r2_positive = paid_upper_min > 0 & inc_upper_min > 0 &
                !is.na(cl_paid_reserve) & cl_paid_reserve > 0,
  r3_volume   = is.finite(cv_prem) & cv_prem < CV_PREM_MAX,
  r4_reins    = is.finite(cv_netdir) & cv_netdir < CV_NET_DIR_MAX)]
stats[, selected := r1_premium & r2_positive & r3_volume & r4_reins]

meyers <- fread(file.path(DIR_REF, "meyers2015_groups.csv"))
stats[, meyers200 := paste(line, grcode) %in% paste(meyers$line, meyers$grcode)]
stopifnot(sum(stats$meyers200) == 200)

waterfall <- stats[, .(
  in_database     = .N,
  premium_present = sum(r1_premium),
  triangle_usable = sum(r1_premium & r2_positive),
  stable_volume   = sum(r1_premium & r2_positive & r3_volume),
  stable_reins    = sum(selected),
  meyers_2015     = sum(meyers200)), by = line]
waterfall <- merge(LINES[, .(line, name)], waterfall, by = "line", sort = FALSE)
waterfall <- rbind(waterfall,
                   waterfall[, c(.(line = "All", name = "All lines"),
                                 lapply(.SD, sum)), .SDcols = -c("line", "name")])
cat("\n=========== SCREENING WATERFALL (company-lines remaining) ===========\n")
print(waterfall)
save_table(waterfall, "01_screening")

# How much of the industry the selected sample covers, by premium
cover <- stats[r1_premium == TRUE, .(
  prem_all = sum(prem_tot), prem_selected = sum(prem_tot[selected])), by = line]
cover[, share := prem_selected / prem_all]
cat("\nShare of 1988-97 net earned premium (among companies with premium) kept:\n")
print(cover[, .(line, share = round(share, 3))])

# ------------------------------------------------------------------------------
# 3. TRIANGLES
# Keep every company-line that is in either the screened sample or Meyers' 200.
# ------------------------------------------------------------------------------
keep <- stats[selected | meyers200]
tri <- lapply(seq_len(nrow(keep)), function(k) {
  d <- long[.(keep$line[k], keep$grcode[k])]
  list(line = keep$line[k], grcode = keep$grcode[k], grname = keep$grname[k],
       selected = keep$selected[k], meyers200 = keep$meyers200[k],
       premium  = d[lag == 1L, prem_net],
       paid     = full_matrix(d, "paid"),
       incurred = full_matrix(d, "incurred"))
})
names(tri) <- paste(keep$line, keep$grcode, sep = "_")

cat(sprintf("\nTriangles kept: %d screened + %d Meyers (%d in both) = %d company-lines\n",
            sum(keep$selected), sum(keep$meyers200),
            sum(keep$selected & keep$meyers200), length(tri)))

# Companies with several lines in the screened sample (candidates for the
# multi-line capital exercise in 05_capital.R)
multi <- stats[selected == TRUE, .(n_lines = .N, lines = paste(line, collapse = " "),
                                   prem_tot = sum(prem_tot)), by = .(grcode)]
multi <- merge(multi, unique(stats[, .(grcode, grname)], by = "grcode"), by = "grcode")
setorder(multi, -n_lines, -prem_tot)
cat("\nCompanies with the most lines in the screened sample:\n")
print(head(multi, 8))
save_table(multi, "01_multiline_companies")

save_table(stats, "01_company_line_stats")
saveRDS(tri, file.path(DIR_PROC, "triangles.rds"))
