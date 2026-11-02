# ==============================================================================
# 06_figures.R
#
# README figures. Reads only the CSVs written by steps 01-05, so it never
# re-fits anything and cannot drift from the reported numbers.
# ==============================================================================

source("R/00_common.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

INK   <- "#1F6F6F"   # primary series
AMBER <- "#D08A3E"   # secondary series
GREY  <- "#6B7280"   # reference lines, annotation
BAND  <- "#E5E7EB"   # KS band
# teal vs amber: OKLab Delta E x100 17.5 (protan) / 25.3 (deutan) -- above the 8 target

theme_set(theme_minimal(base_size = 11) +
          theme(panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold", size = 12),
                plot.subtitle = element_text(colour = GREY, size = 9.5),
                plot.caption = element_text(colour = GREY, size = 8.5, hjust = 0),
                strip.text = element_text(face = "bold", size = 10),
                legend.position = "top"))

MODEL_LAB <- c("Mack incurred" = "Mack, incurred triangles",
               "Mack paid"     = "Mack, paid triangles",
               "ODP paid"      = "Bootstrap ODP, paid triangles")

# ------------------------------------------------------------------------------
# HERO: p-p plots and percentile histograms, all six lines pooled
# ------------------------------------------------------------------------------
pp  <- fread(file.path(DIR_TAB, "04_pp_points.csv"))[line == "All"]
bt  <- fread(file.path(DIR_TAB, "04_uniformity_tests.csv"))[line == "All"]
dec <- melt(fread(file.path(DIR_TAB, "04_percentile_deciles.csv")), id.vars = "decile",
            variable.name = "model", value.name = "n")

pp[, model := factor(MODEL_LAB[model], MODEL_LAB)]
bt[, model := factor(MODEL_LAB[model], MODEL_LAB)]
dec[, model := factor(MODEL_LAB[as.character(model)], MODEL_LAB)]

n_all <- bt$n[1]
crit  <- 1.36 / sqrt(n_all)
band  <- data.table(x = seq(0, 1, length.out = 201))
band[, `:=`(lo = pmax(0, x - crit), hi = pmin(1, x + crit))]

bt[, verdict := sprintf("n = %d, KS D = %.3f, p = %s\n%s",
                        n, ks_d, formatC(ks_p, format = "g", digits = 2),
                        fifelse(mean_pct < 0.45, sprintf("mean percentile %.2f: biased high", mean_pct),
                                sprintf("too narrow: %.0f%% outside 5th-95th", 100 * (share_below_5 + share_above_95))))]

p_pp <- ggplot(pp, aes(expected, pct)) +
  geom_ribbon(data = band, aes(x = x, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = BAND) +
  annotate("segment", x = 0, y = 0, xend = 1, yend = 1, colour = GREY, linewidth = 0.4) +
  geom_line(colour = INK, linewidth = 0.9) +
  geom_text(data = bt, aes(x = 0.03, y = 0.97, label = verdict), hjust = 0, vjust = 1,
            size = 3, colour = "grey15", lineheight = 1.05) +
  facet_wrap(~ model, nrow = 1) +
  scale_x_continuous(labels = scales::percent, breaks = c(0, 0.5, 1)) +
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1)) +
  coord_equal() +
  labs(x = "Expected percentile (uniform)", y = "Actual outcome's percentile, sorted",
       title = "Where the actual run-off landed in each model's predicted distribution",
       subtitle = sprintf("%d company-line triangles across six Schedule P lines, accident years 1988-97, scored against the lower triangle.\nA calibrated model stays on the diagonal; the grey band is the 5%% Kolmogorov-Smirnov region.", n_all))

p_hist <- ggplot(dec, aes(decile, n)) +
  geom_col(fill = INK, width = 0.82) +
  geom_hline(data = bt, aes(yintercept = n / 10), colour = GREY, linewidth = 0.4) +
  facet_wrap(~ model, nrow = 1) +
  scale_x_continuous(breaks = c(1, 5.5, 10), labels = c("0-10%", "", "90-100%")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Decile of the predicted distribution the outcome fell in",
       y = "Triangles",
       caption = "Grey line: the count per decile if the model were calibrated. Mack: lognormal with Mack's mean and standard error. Bootstrap ODP: 10,000 simulations.") +
  theme(panel.grid.major.x = element_blank(), strip.text = element_blank())

ggsave(file.path(DIR_FIG, "00_hero_pp_plots.png"), p_pp / p_hist + plot_layout(heights = c(1.15, 0.6)),
       width = 10, height = 6.6, dpi = 160, bg = "white")

# ------------------------------------------------------------------------------
# P/I SHIFT: the diagnostic an actuary could have run in 1998
# ------------------------------------------------------------------------------
terc <- fread(file.path(DIR_TAB, "04_pi_shift_terciles.csv"))
tt <- unique(terc[, .(tercile, pi_shift_median, cl_paid_median_error, cl_incurred_median_error)])
tt <- melt(tt, id.vars = c("tercile", "pi_shift_median"), variable.name = "basis", value.name = "err")
tt[, basis := fifelse(basis == "cl_paid_median_error", "Paid chain ladder", "Incurred chain ladder")]
tt[, tercile := factor(sprintf("%s\n(median shift %+.1f pts)", tercile, 100 * pi_shift_median),
                       unique(sprintf("%s\n(median shift %+.1f pts)", tercile, 100 * pi_shift_median)))]
tt[, basis := factor(basis, c("Paid chain ladder", "Incurred chain ladder"))]

p_pi <- ggplot(tt, aes(tercile, err, fill = basis)) +
  geom_hline(yintercept = 0, colour = GREY, linewidth = 0.4) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  geom_text(aes(label = sprintf("%+.0f%%", 100 * err),
                vjust = ifelse(err >= 0, -0.45, 1.35)),
            position = position_dodge(width = 0.78), size = 3.2, colour = "grey15") +
  scale_fill_manual(values = c("Paid chain ladder" = AMBER, "Incurred chain ladder" = INK), name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.12, 0.12))) +
  labs(x = "Triangles grouped by the shift in paid-to-incurred ratio on the 1997 diagonal",
       y = "Median reserve error\n(share of actual unpaid)",
       title = "Paid chain ladder over-reserves where claims were being paid faster",
       subtitle = "The shift uses only the upper triangle: it was computable at year end 1997, before any outcome was known") +
  theme(panel.grid.major.x = element_blank())
ggsave(file.path(DIR_FIG, "01_pi_shift.png"), p_pi, width = 8.5, height = 4.6, dpi = 160, bg = "white")

# ------------------------------------------------------------------------------
# BF: error tracks the move in the loss ratio its a priori was taken from
# ------------------------------------------------------------------------------
lrs <- fread(file.path(DIR_TAB, "02_loss_ratio_shift.csv"))
lrs[, lr_change := lr_1993_97 - lr_1988_92]
p_bf <- ggplot(lrs, aes(100 * lr_change, bf_paid_median_error)) +
  geom_hline(yintercept = 0, colour = GREY, linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = GREY, linewidth = 0.4) +
  geom_point(colour = INK, size = 3) +
  geom_text(aes(label = line), nudge_x = 1.4, hjust = 0, size = 3.4, colour = "grey15") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Change in hindsight loss ratio, AY 1993-97 vs 1988-92 (points)",
       y = "BF paid: median reserve error",
       title = "Bornhuetter-Ferguson is only as good as its a priori loss ratio",
       subtitle = "ELR from accident years 1988-92. Where later years ran at a lower loss ratio (WC, PA), BF over-reserved;\nwhere they ran higher (MM, PL), it under-reserved")
ggsave(file.path(DIR_FIG, "02_bf_vs_loss_ratio.png"), p_bf, width = 7.5, height = 4.6, dpi = 160, bg = "white")

# ------------------------------------------------------------------------------
# CAPITAL: aggregate 99.5% capital against the correlation assumption
# ------------------------------------------------------------------------------
agg <- fread(file.path(DIR_TAB, "05_aggregation.csv"))[view == "Ultimate"]
cc  <- fread(file.path(DIR_TAB, "05_company_choice.csv"))[stable == TRUE][1]
long_cap <- rbind(agg[, .(copula, rho, measure = "99.5% VaR", cap = cap_var, sum_st = sum_standalone_var)],
                  agg[, .(copula, rho, measure = "99.5% Expected Shortfall", cap = cap_es, sum_st = sum_standalone_es)])
long_cap[, measure := factor(measure, c("99.5% VaR", "99.5% Expected Shortfall"))]
long_cap[, copula := factor(fifelse(copula == "t", "Student-t copula (4 df)", "Gaussian copula"),
                            c("Gaussian copula", "Student-t copula (4 df)"))]
ref <- long_cap[, .(sum_st = mean(sum_st)), by = measure]
lab <- long_cap[rho == 0.25]

p_cap <- ggplot(long_cap, aes(rho, cap / 1000, colour = copula)) +
  geom_hline(data = ref, aes(yintercept = sum_st / 1000), colour = GREY, linewidth = 0.4) +
  geom_text(data = ref, aes(x = 0, y = sum_st / 1000, label = "sum of standalone capital"),
            inherit.aes = FALSE, hjust = 0, vjust = -0.5, size = 3, colour = GREY) +
  geom_line(linewidth = 0.9) +
  geom_point(data = lab, size = 2.4) +
  facet_wrap(~ measure, scales = "free_y") +
  scale_colour_manual(values = c("Gaussian copula" = INK, "Student-t copula (4 df)" = AMBER), name = NULL) +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(labels = scales::dollar_format(suffix = "M"), expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = "Correlation between lines (same for every pair)",
       y = "Reserve risk capital above the mean",
       title = sprintf("%s: %s aggregated, ultimate view", cc$grname, gsub(" ", ", ", cc$lines)),
       subtitle = "Bootstrap ODP marginals, 100,000 simulations per line, paired by copula rank (Iman-Conover). Dots mark correlation 0.25",
       caption = "Schedule P amounts are in USD thousands; shown here in USD millions.")
ggsave(file.path(DIR_FIG, "03_capital_vs_correlation.png"), p_cap, width = 9, height = 4.6, dpi = 160, bg = "white")

cat("Figures written to outputs/figures/\n")
