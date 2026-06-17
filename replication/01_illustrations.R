# =============================================================
# Section 5: Illustrations
# Reproduces all figures and tables in Section 5 of the manuscript:
#   5.1 Individual growth (simulated E. fetida)
#   5.2 Growth with reproduction (simulated F. candida)
#   5.3 Hierarchical model (simulated 21 individuals)
#   5.4 DEBtox analysis (Cd toxicity, eisenia_cd dataset)
# =============================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

.script_dir <- (function() {
	args <- commandArgs(trailingOnly = FALSE)
	m <- grep("^--file=", args, value = TRUE)
	if (length(m)) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
	ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
	if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
	getwd()
})()
source(file.path(.script_dir, "00_setup.R"))

# Archived fits (loaded unless BDEB_RECOMPUTE=true); see 00_setup.R.
# Note: the archived `fit_cd` is the converged cadmium DEBtox fit used for
# Table 4 and Figure 2.  With BDEB_RECOMPUTE=true the eisenia_cd group-mean
# data are refit here; `init = 0.5` is required for this stiff TKTD ODE to
# converge (matching the manuscript code).
cache01 <- load_cache(file.path(OUTPUTS_DIR, "01_illustrations_fits.rds"))


# -------------------------------------------------------------
# 5.1 Individual growth — simulated E. fetida
# -------------------------------------------------------------

cli::cli_h1("Section 5.1: Individual growth (simulated)")

data("eisenia_growth")
df1 <- eisenia_growth[eisenia_growth$id == 5, ]
dat <- bdeb_data(growth = df1, f_food = 1.0)

mod <- bdeb_model(dat, type = "individual",
  priors = list(
    p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
    kappa   = prior_beta(a = 3, b = 2),
    sigma_L = prior_halfnormal(sigma = 0.05)))

# Code listing on manuscript p. 7: prior_species("Eisenia_fetida")
# is illustrated here -- inclusion required by reviewer comment.
cat("\n--- prior_species('Eisenia_fetida') ---\n")
print(prior_species("Eisenia_fetida"))

fit <- cache01$fit %||% bdeb_fit(mod,
                chains        = ITER_CFG$chains,
                iter_warmup   = ITER_CFG$warmup,
                iter_sampling = ITER_CFG$sampling,
                adapt_delta   = 0.9, seed = 42, refresh = 0)

print(bdeb_diagnose(fit))

cat("\n--- Section 5.1: Derived quantities (Table 3) ---\n")
der <- bdeb_derived(fit,
  quantities = c("L_m", "L_inf", "k_M", "g", "growth_rate"),
  f = 1.0)
print(as.data.frame(summarise_draws(der, "mean", "sd",
  "q5"  = ~quantile(.x, 0.05),
  "q95" = ~quantile(.x, 0.95))), digits = 3)

# LOO comparison: Gaussian vs lognormal observation
mod_ln <- bdeb_model(dat, type = "individual",
  observation = list(growth = obs_lognormal()),
  priors = list(
    p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
    kappa   = prior_beta(a = 3, b = 2),
    sigma_L = prior_halfnormal(sigma = 0.05)))

fit_ln <- cache01$fit_ln %||% bdeb_fit(mod_ln,
                   chains        = ITER_CFG$chains,
                   iter_warmup   = ITER_CFG$warmup,
                   iter_sampling = ITER_CFG$sampling,
                   adapt_delta   = 0.9, seed = 42, refresh = 0)

cat("\n--- Section 5.1: LOO comparison ---\n")
# Note: with only n = 13 observations, loo() emits a Pareto-k warning
# (a few k > 0.7); the importance-sampling ELPD is therefore indicative.
# The qualitative conclusion (models indistinguishable, |dELPD| < SE) is
# robust across reruns -- see manuscript Section 5.1.
print(loo::loo_compare(list(gaussian  = bdeb_loo(fit),
                            lognormal = bdeb_loo(fit_ln))))

# Manuscript listing on p. 18: switching to a Student-t observation
# model is a one-liner (no Stan recompilation).  Construction shown
# here to satisfy the reviewer's audit; not fitted in this section.
mod_t <- bdeb_model(dat, type = "individual",
  observation = list(growth = obs_student_t(nu = 5)),
  priors = list(
    p_Am    = prior_lognormal(mu = 1.5, sigma = 0.5),
    kappa   = prior_beta(a = 3, b = 2),
    sigma_L = prior_halfnormal(sigma = 0.05)))
cat("\n--- Section 5.1: Student-t observation model spec ---\n")
print(mod_t)


# -------------------------------------------------------------
# 5.2 Growth with reproduction — simulated F. candida
# -------------------------------------------------------------

cli::cli_h1("Section 5.2: Growth with reproduction (simulated)")

growth_fc <- data.frame(id = 1,
  time = seq(0, 28, by = 2),
  length = c(0.048, 0.053, 0.058, 0.063, 0.067, 0.071,
             0.075, 0.078, 0.081, 0.084, 0.086, 0.088,
             0.089, 0.090, 0.091))
repro_fc <- data.frame(id = 1, t_start = 0, t_end = 28, count = 118)

dat_gr <- bdeb_data(growth = growth_fc, reproduction = repro_fc,
                    f_food = 1.0)

mod_gr <- bdeb_model(dat_gr, type = "growth_repro",
  priors = prior_species("Folsomia_candida", type = "growth_repro"),
  observation = list(growth = obs_normal(),
                     reproduction = obs_negbinom()))

fit_gr <- cache01$fit_gr %||% bdeb_fit(mod_gr,
                   chains        = ITER_CFG$chains,
                   iter_warmup   = ITER_CFG$warmup,
                   iter_sampling = ITER_CFG$sampling,
                   adapt_delta   = 0.9, seed = 42, refresh = 0)

cat("\n--- Section 5.2: Diagnostics ---\n")
print(bdeb_diagnose(fit_gr))

cat("\n--- Section 5.2: Posterior summary ---\n")
print(as.data.frame(summary(fit_gr,
  pars = c("p_Am", "p_M", "kappa", "k_J", "k_R", "phi_R", "sigma_L"),
  prob = 0.90)), digits = 3)


# -------------------------------------------------------------
# 5.3 Hierarchical model — simulated, 21 individuals
# -------------------------------------------------------------

cli::cli_h1("Section 5.3: Hierarchical model")

dat_all <- bdeb_data(growth = eisenia_growth, f_food = 1.0)

mod_h <- bdeb_model(dat_all, type = "hierarchical",
  priors = list(
    mu_log_p_Am    = prior_normal(mu = 1.5, sigma = 0.5),
    sigma_log_p_Am = prior_exponential(rate = 2)))

fit_h <- cache01$fit_h %||% bdeb_fit(mod_h,
                  chains        = ITER_CFG$chains,
                  iter_warmup   = max(ITER_CFG$warmup, 600L),
                  iter_sampling = ITER_CFG$sampling,
                  adapt_delta   = if (BDEB_MODE == "full") 0.95 else 0.99,
                  max_treedepth = if (BDEB_MODE == "full") 10L else 8L,
                  threads_per_chain = if (BDEB_MODE == "full") 4 else 2,
                  seed = 123, refresh = 0)

cat("\n--- Section 5.3: Hierarchical diagnostics ---\n")
print(bdeb_diagnose(fit_h))


# -------------------------------------------------------------
# 5.4 DEBtox analysis — Cd toxicity (eisenia_cd dataset)
# -------------------------------------------------------------

cli::cli_h1("Section 5.4: DEBtox -- Cd toxicity")

data("eisenia_cd")
conc_map <- setNames(unique(eisenia_cd$concentration),
                     as.character(unique(eisenia_cd$concentration)))
dat_cd <- bdeb_data(growth = eisenia_cd, concentration = conc_map,
                    f_food = 1.0)

# Figure: raw data (fig:debtox_rawdata)
p_debtox_raw <- ggplot(eisenia_cd,
  aes(time, length, colour = factor(concentration))) +
  geom_point(size = 2) + geom_line() +
  theme_bw(base_size = 12) +
  labs(x = "Time (days)", y = "Structural length (cm)",
       colour = "Cd (mg/kg)")
save_fig(p_debtox_raw, "fig_debtox_rawdata.pdf", width = 7, height = 4)

# Code listing on manuscript p. 7: species-specific DEBtox priors.
# Included here so every manuscript listing appears in the replication
# (reviewer comment, 2026-05 round).
cat("\n--- prior_species('Daphnia_magna', type = 'debtox') ---\n")
print(prior_species("Daphnia_magna", type = "debtox"))

mod_cd <- bdeb_tox(dat_cd, stress = "assimilation",
  priors = list(
    z_w = prior_lognormal(mu = 3.0, sigma = 1.0),
    b_w = prior_lognormal(mu = -4.0, sigma = 1.5)))

fit_cd <- cache01$fit_cd %||% bdeb_fit(mod_cd,
                   chains        = ITER_CFG$chains,
                   iter_warmup   = max(ITER_CFG$warmup, 600L),
                   iter_sampling = ITER_CFG$sampling,
                   adapt_delta   = if (BDEB_MODE == "full") 0.95 else 0.99,
                   max_treedepth = if (BDEB_MODE == "full") 10L else 8L,
                   init = 0.5,   # essential for this stiff TKTD ODE to converge
                   seed = 42, refresh = 0)

cat("\n--- Section 5.4: DEBtox posterior (Table 4) ---\n")
print(as.data.frame(summary(fit_cd,
  pars = c("p_Am", "p_M", "kappa", "k_d", "z_w", "b_w", "sigma_L"),
  prob = 0.90)), digits = 3)

ec <- bdeb_ec50(fit_cd, prob = 0.90)
cat("\n--- Section 5.4: EC50 / NEC ---\n")
print(ec$summary, digits = 3)

# Figure: dose-response (fig:debtox_doseresponse)
p_dr <- plot_dose_response(fit_cd, n_draws = 100, seed = 1)
save_fig(p_dr, "fig_debtox_doseresponse.pdf", width = 6, height = 4)


# -------------------------------------------------------------
# Persist key fits for downstream scripts
# -------------------------------------------------------------

# Persist only when refitting; the shipped archive is authoritative and
# must not be overwritten by a cache-mode rerun.
if (BDEB_RECOMPUTE) {
	saveRDS(list(fit = fit, fit_ln = fit_ln, fit_gr = fit_gr,
	             fit_h = fit_h, fit_cd = fit_cd),
	        file.path(OUTPUTS_DIR, "01_illustrations_fits.rds"))
}

cli::cli_alert_success("Section 5 illustrations: done.")
