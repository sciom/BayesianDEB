## =============================================================
## Real-data DEBtox case study: Cd toxicity to E. andrei
## Data: Van Gestel et al. (1991) via EGrowth (Mathieu 2018)
## =============================================================

library(BayesianDEB)
library(posterior)
library(ggplot2)

## --- 1. Download and prepare data --------------------------------

curves <- read.delim(
  "https://raw.githubusercontent.com/JeromeMathieuEcology/EGrowth/master/curves.txt"
)

# Natural soil Cd series: gr0119 (0), gr0120 (10), gr0121 (32),
#                         gr0122 (100), gr0123 (320 mg Cd/kg)
ids <- c("gr0119", "gr0120", "gr0121", "gr0122", "gr0123")
conc_map <- c(gr0119 = 0, gr0120 = 10, gr0121 = 32,
              gr0122 = 100, gr0123 = 320)

vg <- curves[curves$CURVE_ID %in% ids, ]

# Convert wet mass (mg) to structural length (cm)
d_V <- 1.05  # g/cm^3
vg$length <- ((vg$bm / 1000) / d_V)^(1/3)
vg$concentration <- conc_map[vg$CURVE_ID]

# Format: each curve = one "group" (group means, not individuals)
# Shift t=0 to small positive value (Stan ODE requires ts > t0)
vg$time[vg$time == 0] <- 0.01

debtox_df <- data.frame(
  id = as.integer(factor(vg$CURVE_ID, levels = ids)),
  time = vg$time,
  length = vg$length,
  concentration = vg$concentration
)

cat("=== Van Gestel 1991 Cd data ===\n")
cat("Concentrations:", sort(unique(debtox_df$concentration)), "mg Cd/kg\n")
cat("Time points per group:", table(debtox_df$id)[1], "\n")
cat("Length range:", round(range(debtox_df$length), 3), "cm\n\n")

# Raw data plot
p_raw <- ggplot(debtox_df, aes(time, length,
                                colour = factor(concentration))) +
  geom_point(size = 2) + geom_line(alpha = 0.5) +
  theme_bw() +
  labs(x = "Time (days)", y = "Structural length (cm)",
       colour = "Cd (mg/kg)",
       title = "E. andrei growth under Cd exposure (Van Gestel 1991)")
ggsave("fig_debtox_rawdata.pdf", p_raw, width = 6, height = 4)

## --- 2. Fit DEBtox model -----------------------------------------

conc_named <- setNames(unique(debtox_df$concentration),
                       as.character(unique(debtox_df$concentration)))

dat <- bdeb_data(growth = debtox_df,
                 concentration = conc_named, f_food = 1.0)

mod <- bdeb_tox(dat, stress = "assimilation",
  priors = list(
    z_w = prior_lognormal(mu = 3.0, sigma = 1.0),
    b_w = prior_lognormal(mu = -4.0, sigma = 1.5)
  )
)

cat("Fitting DEBtox model...\n")
fit <- bdeb_fit(mod, chains = 4, iter_sampling = 2000,
                adapt_delta = 0.95, max_treedepth = 12,
                parallel_chains = 4,
                refresh = 200, seed = 77,
                init = 0.5)

## --- 3. Diagnostics ----------------------------------------------

diag <- bdeb_diagnose(fit)
cat("\n=== Diagnostics ===\n")
cat("Divergences:", diag$n_divergent, "\n")
cat("Max treedepth:", diag$n_max_treedepth, "\n")

## --- 4. EC50 and NEC ---------------------------------------------

ec50 <- bdeb_ec50(fit, prob = 0.90)
cat("\n=== EC50 / NEC ===\n")
cat(sprintf("EC50: %.1f [%.1f, %.1f] mg Cd/kg\n",
            ec50$summary$median, ec50$summary$q5, ec50$summary$q95))
cat(sprintf("NEC:  %.1f [%.1f, %.1f] mg Cd/kg\n",
            median(ec50$NEC), quantile(ec50$NEC, 0.05),
            quantile(ec50$NEC, 0.95)))

## --- 5. Posterior summary ----------------------------------------

tox_pars <- c("p_Am", "p_M", "kappa", "v", "E_G",
              "sigma_L", "k_d", "z_w", "b_w")
smry <- bdeb_summary(fit, pars = tox_pars, prob = 0.90)
cat("\n=== Posterior summary ===\n")
for (i in seq_len(nrow(smry))) {
  cat(sprintf("%-10s  %.4f  [%.4f, %.4f]\n",
              smry$variable[i], smry$median[i], smry[[5]][i], smry[[6]][i]))
}

cat("\n=== DONE ===\n")
