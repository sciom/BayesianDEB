## =============================================================
## Real data case study: Eisenia fetida growth
## Data source: EGrowth database (Mathieu 2018)
##   - Neuhauser, Hartenstein & Kaplan (1980) activated sludge
##   - Helling et al. (2000) Cu ecotoxicology
## =============================================================

library(BayesianDEB)
library(ggplot2)

## --- 1. Download data from EGrowth (Zenodo/GitHub) ---

curves <- read.delim(
  "https://raw.githubusercontent.com/JeromeMathieuEcology/EGrowth/master/curves.txt"
)

## --- 2. Neuhauser 1980: E. fetida growth on activated sludge ---
## 37 time points over 250 days, group means (n=20 worms), 25 C
## Body mass in mg -> convert to structural length

neuh <- curves[curves$CURVE_ID == "gr0226", c("time", "bm")]
names(neuh) <- c("time", "wet_mass_mg")

# Convert wet mass to structural length:
#   W = d_V * V + d_E * E  ~  d_V * V  (for structural mass dominance)
#   V = W / d_V,  L = V^(1/3)
#   For E. fetida: d_V ~ 1.05 g/cm^3, delta_M ~ 0.24
#   Structural volume: V = (W_mg / 1000) / d_V
#   Structural length: L = V^(1/3)
d_V <- 1.05  # g/cm^3, structural density
neuh$V <- (neuh$wet_mass_mg / 1000) / d_V  # cm^3
neuh$length <- neuh$V^(1/3)  # structural length in cm

# Prepare for BayesianDEB (single "individual" = group mean)
eisenia_real <- data.frame(
  id = 1,
  time = neuh$time,
  length = neuh$length
)

cat("=== Neuhauser 1980 growth data ===\n")
cat("Time points:", nrow(eisenia_real), "\n")
cat("Time range:", range(eisenia_real$time), "days\n")
cat("Length range:", round(range(eisenia_real$length), 3), "cm\n")

# Plot raw data
p_raw <- ggplot(eisenia_real, aes(time, length)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.3) +
  theme_bw() +
  labs(x = "Time (days)", y = "Structural length (cm)",
       title = "E. fetida growth (Neuhauser et al. 1980)")
print(p_raw)
ggsave("fig_eisenia_raw.pdf", p_raw, width = 6, height = 4)

## --- 3. Fit individual DEB model ---

dat <- bdeb_data(growth = eisenia_real, f_food = 1.0)

mod <- bdeb_model(dat, type = "individual",
  priors = list(
    p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
    p_M     = prior_lognormal(mu = -1.0, sigma = 0.5),
    kappa   = prior_beta(a = 3, b = 2),
    v       = prior_lognormal(mu = -1.5, sigma = 0.5),
    E_G     = prior_lognormal(mu = 6.0, sigma = 0.5),
    sigma_L = prior_halfnormal(sigma = 0.05)
  ),
  temperature = list(T = 298.15, T_ref = 293.15, T_A = 8000)  # 25C
)

fit <- bdeb_fit(mod, chains = 4, iter_sampling = 2000,
                adapt_delta = 0.95, max_treedepth = 12, seed = 42)

## --- 4. Diagnostics ---

diag <- bdeb_diagnose(fit)
cat("\n=== Diagnostics ===\n")
cat("Divergences:", diag$n_divergent, "\n")
cat("Max treedepth:", diag$n_max_treedepth, "\n")

# Trace plots
p_trace <- plot(fit, type = "trace",
                pars = c("p_Am", "p_M", "kappa", "sigma_L"))
ggsave("fig_trace.pdf", p_trace, width = 8, height = 6)

# Posterior densities
p_post <- plot(fit, type = "posterior",
               pars = c("p_Am", "p_M", "kappa", "v", "E_G"))
ggsave("fig_posterior.pdf", p_post, width = 8, height = 6)

# Pairs
p_pairs <- plot(fit, type = "pairs",
                pars = c("p_Am", "p_M", "kappa", "E_G"))
ggsave("fig_pairs.pdf", p_pairs, width = 8, height = 8)

## --- 5. Posterior summary ---

cat("\n=== Posterior summary ===\n")
s <- bdeb_summary(fit, pars = c("p_Am", "p_M", "kappa", "v", "E_G",
                                 "sigma_L"), prob = 0.90)
print(as.data.frame(s), digits = 3)

## --- 6. Derived quantities ---

der <- bdeb_derived(fit, quantities = c("L_m", "L_inf", "k_M",
                                         "g", "growth_rate"), f = 1.0)
cat("\n=== Derived quantities ===\n")
print(as.data.frame(posterior::summarise_draws(der, "mean", "sd",
  "q5" = ~quantile(.x, 0.05), "q95" = ~quantile(.x, 0.95))),
  digits = 3)

## --- 7. Trajectory plot ---

p_traj <- plot(fit, type = "trajectory", n_draws = 200, seed = 1)
ggsave("fig_trajectory.pdf", p_traj, width = 6, height = 4)

## --- 8. PPC ---

ppc <- bdeb_ppc(fit, type = "growth")
p_ppc <- plot(ppc, n_draws = 200)
ggsave("fig_ppc.pdf", p_ppc, width = 6, height = 4)

## --- 9. Posterior correlation matrix ---

draws <- posterior::as_draws_matrix(fit$fit$draws())
pars_draws <- draws[, c("p_Am", "p_M", "kappa", "v", "E_G", "sigma_L")]
cormat <- cor(pars_draws, method = "spearman")
cat("\n=== Posterior correlation (Spearman) ===\n")
print(round(cormat, 2))

## --- 10. LOO-CV ---

loo_result <- bdeb_loo(fit)
cat("\n=== LOO-CV ===\n")
print(loo_result)

## --- 11. Session info ---

bdeb_session_info(fit)

cat("\n=== DONE ===\n")
