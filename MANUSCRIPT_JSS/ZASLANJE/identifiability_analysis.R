## =============================================================
## Identifiability analysis for BayesianDEB
## Prior-to-posterior contraction as a function of data size
## =============================================================
##
## Quantifies how much the data narrow the posterior relative to
## the prior, for different numbers of observations.
##
## Output:
##   fig_contraction.pdf    - contraction vs N_obs figure
##   (console)              - LaTeX table for manuscript
## =============================================================

library(BayesianDEB)
library(posterior)
library(ggplot2)

## --- True parameters (E. fetida AmP values) ----------------------

TRUE_PARS <- list(
  p_Am = 5.0, p_M = 0.5, kappa = 0.75, v = 0.2,
  E_G = 400, E0 = 1.0, L0 = 0.1, sigma_L = 0.015
)

## --- Priors (same as manuscript individual growth example) --------

fit_priors <- list(
  p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
  p_M     = prior_lognormal(mu = -1.0, sigma = 0.5),
  kappa   = prior_beta(a = 3, b = 2),
  v       = prior_lognormal(mu = -1.5, sigma = 0.5),
  E_G     = prior_lognormal(mu = 6.0, sigma = 0.5),
  sigma_L = prior_halfnormal(sigma = 0.05)
)

## --- Compute prior 90% CI widths --------------------------------

set.seed(42)
N_prior <- 50000
prior_samp <- data.frame(
  p_Am    = rlnorm(N_prior, 1.5, 0.5),
  p_M     = rlnorm(N_prior, -1.0, 0.5),
  kappa   = rbeta(N_prior, 3, 2),
  v       = rlnorm(N_prior, -1.5, 0.5),
  E_G     = rlnorm(N_prior, 6.0, 0.5),
  sigma_L = abs(rnorm(N_prior, 0, 0.05))
)
# Derived quantities from prior
prior_samp$L_inf <- with(prior_samp, kappa * p_Am / p_M)
prior_samp$growth_rate <- with(prior_samp, {
  k_M <- p_M / E_G
  g   <- E_G * v / (kappa * p_Am)
  k_M * g / (3 * (1 + g))  # f = 1
})

all_pars <- c("p_Am", "p_M", "kappa", "v", "E_G", "sigma_L",
              "L_inf", "growth_rate")
prior_width <- sapply(prior_samp[all_pars], function(x)
  unname(diff(quantile(x, c(0.05, 0.95)))))

cat("=== Prior 90% CI widths ===\n")
print(round(prior_width, 4))

## --- Simulate base trajectory ------------------------------------

sim <- deb_simulate(t_max = 84, p_Am = 5, p_M = 0.5, kappa = 0.75,
                    v = 0.2, E_G = 400, E0 = 1, L0 = 0.1, f = 1, dt = 0.05)

## --- Time-point configurations -----------------------------------

t_configs <- list(
  "5"  = seq(0, 84, length.out = 5),
  "9"  = seq(0, 84, length.out = 9),
  "13" = seq(0, 84, by = 7),
  "25" = seq(0, 84, length.out = 25),
  "37" = seq(0, 84, length.out = 37)
)

## --- Fit each configuration --------------------------------------

cat("\n=== Fitting models ===\n")

results <- list()
for (name in names(t_configs)) {
  t_obs <- t_configs[[name]]
  L_true <- approx(sim$time, sim$L, xout = t_obs)$y

  set.seed(123)
  L_obs <- rnorm(length(t_obs), L_true, TRUE_PARS$sigma_L)

  df <- data.frame(id = 1, time = t_obs, length = L_obs)
  dat <- bdeb_data(growth = df, f_food = 1.0)
  mod <- bdeb_model(dat, type = "individual", priors = fit_priors)

  fit <- bdeb_fit(mod, chains = 4, iter_sampling = 2000,
                  adapt_delta = 0.95, parallel_chains = 4,
                  refresh = 0, seed = 42)

  # Extract posterior as plain data frame
  draws <- as.data.frame(as_draws_matrix(fit$fit$draws()))

  # Compute derived quantities directly from posterior draws
  draws$L_inf <- draws$kappa * draws$p_Am / draws$p_M
  draws$growth_rate <- with(draws, {
    k_M <- p_M / E_G
    g   <- E_G * v / (kappa * p_Am)
    k_M * g / (3 * (1 + g))
  })

  # Posterior 90% CI widths
  post_width <- sapply(draws[, all_pars, drop = FALSE], function(x)
    unname(diff(quantile(x, c(0.05, 0.95)))))

  # Contraction: 1 - posterior_width / prior_width
  contraction <- 1 - post_width / prior_width
  contraction <- pmax(contraction, 0)  # floor at 0

  results[[name]] <- list(
    n_obs = length(t_obs),
    post_width = post_width,
    contraction = contraction
  )

  cat(sprintf("N_obs = %2d: done\n", length(t_obs)))
}

## --- Contraction table -------------------------------------------

cat("\n=== Contraction ratios (1 = fully data-determined, 0 = prior-dominated) ===\n\n")

cat(sprintf("%-15s", "Parameter"))
for (name in names(results)) {
  cat(sprintf(" N=%2d", results[[name]]$n_obs))
}
cat("\n")
cat(paste(rep("-", 15 + 6 * length(results)), collapse = ""), "\n")

for (p in all_pars) {
  cat(sprintf("%-15s", p))
  for (name in names(results)) {
    cat(sprintf(" %4.2f", results[[name]]$contraction[p]))
  }
  cat("\n")
}

## --- LaTeX table -------------------------------------------------

cat("\n=== LaTeX table rows ===\n\n")

tex_names <- c(
  p_Am = "$\\pAm$", p_M = "$\\pM$", kappa = "$\\kappa$",
  v = "$v$", E_G = "$\\EG$", sigma_L = "$\\sigma_L$",
  L_inf = "$L_\\infty$", growth_rate = "$\\dot{r}_B$"
)

n_configs <- c("5", "13", "37")
for (p in all_pars) {
  vals <- sapply(n_configs, function(n) results[[n]]$contraction[p])
  cat(sprintf("%s & %s \\\\\n",
              tex_names[p],
              paste(sprintf("%.2f", vals), collapse = " & ")))
}

## --- Contraction figure ------------------------------------------

cont_df <- do.call(rbind, lapply(names(results), function(name) {
  data.frame(
    n_obs = results[[name]]$n_obs,
    parameter = all_pars,
    contraction = results[[name]]$contraction[all_pars],
    stringsAsFactors = FALSE
  )
}))

# Labels for facets/legend
par_labels <- c(
  p_Am = "p[Am]", p_M = "p[M]", kappa = "kappa",
  v = "v", E_G = "E[G]", sigma_L = "sigma[L]",
  L_inf = "L[infinity]", growth_rate = "dot(r)[B]"
)
cont_df$par_label <- factor(par_labels[cont_df$parameter],
                            levels = par_labels)

# Split into primary and derived
cont_df$type <- ifelse(cont_df$parameter %in% c("L_inf", "growth_rate"),
                       "Derived", "Primary")

p_cont <- ggplot(cont_df, aes(x = n_obs, y = contraction,
                               colour = par_label, linetype = type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey40") +
  annotate("text", x = 36, y = 0.53, label = "c = 0.5", size = 3,
           colour = "grey40", hjust = 1) +
  scale_colour_manual(
    values = c("p[Am]" = "#E41A1C", "p[M]" = "#377EB8",
               "kappa" = "#4DAF4A", "v" = "#984EA3",
               "E[G]" = "#FF7F00", "sigma[L]" = "#A65628",
               "L[infinity]" = "#000000", "dot(r)[B]" = "#666666"),
    labels = scales::parse_format()
  ) +
  scale_linetype_manual(values = c("Primary" = "solid", "Derived" = "dashed")) +
  scale_x_continuous(breaks = c(5, 9, 13, 25, 37)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 10) +
  labs(x = "Number of observations",
       y = "Contraction  c",
       colour = "Parameter", linetype = "Type") +
  theme(legend.position = "right",
        legend.text.align = 0)

ggsave("fig_contraction.pdf", p_cont, width = 7, height = 4)
cat("\nContraction figure saved to fig_contraction.pdf\n")

cat("=== DONE ===\n")
