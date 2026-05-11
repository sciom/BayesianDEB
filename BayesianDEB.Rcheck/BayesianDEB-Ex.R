pkgname <- "BayesianDEB"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "BayesianDEB-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('BayesianDEB')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("arrhenius")
### * arrhenius

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: arrhenius
### Title: Arrhenius Temperature Correction
### Aliases: arrhenius

### ** Examples

# Correction at 25 C relative to 20 C reference
arrhenius(298.15, T_ref = 293.15, T_A = 8000)  # ~ 1.74

# No correction at reference temperature
arrhenius(293.15)  # exactly 1



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("arrhenius", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bdeb_data")
### * bdeb_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bdeb_data
### Title: Prepare Data for BDEB Models
### Aliases: bdeb_data

### ** Examples

# Simple growth data
df <- data.frame(
  id = rep(1, 10),
  time = seq(0, 45, by = 5),
  length = c(0.1, 0.15, 0.22, 0.30, 0.38, 0.44, 0.49, 0.52, 0.54, 0.55)
)
dat <- bdeb_data(growth = df)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bdeb_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bdeb_derived")
### * bdeb_derived

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bdeb_derived
### Title: Compute Derived Biological Quantities from the Posterior
### Aliases: bdeb_derived bdeb_derived.default bdeb_derived.bdeb_fit

### ** Examples

## Not run: 
##D data(eisenia_growth)
##D dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
##D fit <- bdeb_fit(bdeb_model(dat, type = "individual"))
##D bdeb_derived(fit, quantities = c("L_inf", "k_M"))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bdeb_derived", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bdeb_diagnose")
### * bdeb_diagnose

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bdeb_diagnose
### Title: MCMC Convergence Diagnostics
### Aliases: bdeb_diagnose

### ** Examples

## Not run: 
##D data(eisenia_growth)
##D dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
##D fit <- bdeb_fit(bdeb_model(dat, type = "individual"),
##D                 chains = 2, iter_warmup = 200, iter_sampling = 200)
##D d <- bdeb_diagnose(fit)
##D summary(d)
##D plot(d, type = "rhat")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bdeb_diagnose", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bdeb_ec50")
### * bdeb_ec50

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bdeb_ec50
### Title: Extract EC50 and NEC from a DEBtox Fit
### Aliases: bdeb_ec50

### ** Examples

## Not run: 
##D data(debtox_growth)
##D dat <- bdeb_data(growth = debtox_growth,
##D                  concentration = c("1" = 0, "11" = 20,
##D                                   "21" = 80, "31" = 200))
##D fit <- bdeb_fit(bdeb_tox(dat, stress = "assimilation"))
##D bdeb_ec50(fit, prob = 0.95)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bdeb_ec50", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bdeb_fit")
### * bdeb_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bdeb_fit
### Title: Fit a BDEB Model via Hamiltonian Monte Carlo
### Aliases: bdeb_fit

### ** Examples

# bdeb_fit() requires the external CmdStan toolchain (Suggests:
# cmdstanr) and a fresh fit takes > 30 seconds (Stan compilation
# plus MCMC).  The example is wrapped in \donttest{} and gated on
# cmdstanr availability so that R CMD check skips it on CRAN's
# toolchain-free workers but power users can run it after
# cmdstanr::install_cmdstan().
## No test: 
if (requireNamespace("cmdstanr", quietly = TRUE) &&
    nzchar(tryCatch(cmdstanr::cmdstan_path(),
                    error = function(e) ""))) {
  data(eisenia_growth)
  dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
  mod <- bdeb_model(dat, type = "individual")
  fit <- bdeb_fit(mod, chains = 1, iter_warmup = 200,
                  iter_sampling = 200, refresh = 0)
  print(fit)
}
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bdeb_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bdeb_loo")
### * bdeb_loo

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bdeb_loo
### Title: LOO Cross-Validation for Model Comparison
### Aliases: bdeb_loo

### ** Examples

## Not run: 
##D data(eisenia_growth)
##D dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
##D fit <- bdeb_fit(bdeb_model(dat, type = "individual"))
##D bdeb_loo(fit)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bdeb_loo", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bdeb_model")
### * bdeb_model

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bdeb_model
### Title: Specify a BDEB Model
### Aliases: bdeb_model

### ** Examples

df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
dat <- bdeb_data(growth = df)
mod <- bdeb_model(dat, type = "individual")
print(mod)
summary(mod)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bdeb_model", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bdeb_ppc")
### * bdeb_ppc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bdeb_ppc
### Title: Posterior Predictive Checks for BDEB Models
### Aliases: bdeb_ppc

### ** Examples

## Not run: 
##D data(eisenia_growth)
##D dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
##D fit <- bdeb_fit(bdeb_model(dat, type = "individual"))
##D ppc <- bdeb_ppc(fit, type = "growth")
##D plot(ppc)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bdeb_ppc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bdeb_prior-methods")
### * bdeb_prior-methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bdeb_prior-methods
### Title: Methods for 'bdeb_prior' Objects
### Aliases: bdeb_prior-methods print.bdeb_prior summary.bdeb_prior
###   print.summary.bdeb_prior plot.bdeb_prior

### ** Examples

p <- prior_lognormal(mu = 1.5, sigma = 0.5)
print(p)
summary(p)
plot(p)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bdeb_prior-methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bdeb_prior_predictive")
### * bdeb_prior_predictive

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bdeb_prior_predictive
### Title: Prior Predictive Check
### Aliases: bdeb_prior_predictive

### ** Examples

dat <- bdeb_data(growth = data.frame(id = 1, time = 0:10,
  length = seq(0.1, 0.5, length.out = 11)))
mod <- bdeb_model(dat, type = "individual")
pp <- bdeb_prior_predictive(mod, n_draws = 100, seed = 1)
plot(pp, n_draws = 30)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bdeb_prior_predictive", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bdeb_session_info")
### * bdeb_session_info

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bdeb_session_info
### Title: Reproducibility Report
### Aliases: bdeb_session_info

### ** Examples

bdeb_session_info()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bdeb_session_info", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bdeb_summary")
### * bdeb_summary

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bdeb_summary
### Title: Posterior Summary for BDEB Parameters (deprecated)
### Aliases: bdeb_summary
### Keywords: internal

### ** Examples

## Not run: 
##D fit <- bdeb_fit(mod)
##D bdeb_summary(fit)  # equivalent to summary(fit); deprecated
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bdeb_summary", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bdeb_tox")
### * bdeb_tox

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bdeb_tox
### Title: DEBtox Model Specification
### Aliases: bdeb_tox

### ** Examples

# R-side specification only (no Stan sampling)
data(debtox_growth)
# one replicate per concentration avoids aggregation warning
dt <- debtox_growth[debtox_growth$id %in% c(1, 11, 21, 31), ]
dat <- bdeb_data(growth = dt, concentration = c(0, 20, 80, 200))
mod <- bdeb_tox(dat, stress = "assimilation")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bdeb_tox", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("coef.bdeb_fit")
### * coef.bdeb_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: coef.bdeb_fit
### Title: Extract Point Estimates from a BDEB Fit
### Aliases: coef.bdeb_fit

### ** Examples

## Not run: 
##D fit <- bdeb_fit(mod)
##D coef(fit)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("coef.bdeb_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("confint.bdeb_fit")
### * confint.bdeb_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: confint.bdeb_fit
### Title: Posterior Credible Intervals for BDEB Model Parameters
### Aliases: confint.bdeb_fit

### ** Examples

## Not run: 
##D fit <- bdeb_fit(mod)
##D confint(fit, level = 0.90)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("confint.bdeb_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("deb_fluxes")
### * deb_fluxes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: deb_fluxes
### Title: Compute DEB Energy Fluxes
### Aliases: deb_fluxes

### ** Examples

# Energy fluxes for a 1 cm^3 organism with typical Eisenia parameters
deb_fluxes(E = 1, V = 1, f = 1, p_Am = 5, p_M = 0.5,
  kappa = 0.75, v = 0.2, E_G = 4400)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("deb_fluxes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("deb_simulate")
### * deb_simulate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: deb_simulate
### Title: Simulate DEB Growth Trajectory
### Aliases: deb_simulate

### ** Examples

# Simulate E. fetida growth for 84 days
traj <- deb_simulate(t_max = 84, p_Am = 5, p_M = 0.5,
  kappa = 0.75, v = 0.2, E_G = 400, E0 = 1, L0 = 0.1)
plot(traj$time, traj$L, type = "l", xlab = "Days", ylab = "L (cm)")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("deb_simulate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("debtox_growth")
### * debtox_growth

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: debtox_growth
### Title: Simulated DEBtox Growth Data
### Aliases: debtox_growth
### Keywords: datasets

### ** Examples

data(debtox_growth)
head(debtox_growth)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("debtox_growth", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("debtox_simulate")
### * debtox_simulate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: debtox_simulate
### Title: Simulate DEBtox Growth Under Toxicant Exposure
### Aliases: debtox_simulate

### ** Examples

traj <- debtox_simulate(t_max = 42, p_Am = 5, p_M = 0.5,
  kappa = 0.75, v = 0.2, E_G = 400, E0 = 1, L0 = 0.1,
  k_d = 0.3, z_w = 15, b_w = 0.003, C_w = 80)
plot(traj$time, traj$L, type = "l")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("debtox_simulate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("eisenia_cd")
### * eisenia_cd

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: eisenia_cd
### Title: Eisenia andrei Cadmium Toxicity Data (Van Gestel 1991)
### Aliases: eisenia_cd
### Keywords: datasets

### ** Examples

data(eisenia_cd)
plot(eisenia_cd$time, eisenia_cd$length,
     col = as.factor(eisenia_cd$concentration),
     xlab = "Days", ylab = "Structural length (cm)")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("eisenia_cd", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("eisenia_growth")
### * eisenia_growth

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: eisenia_growth
### Title: Simulated Eisenia fetida Growth Data
### Aliases: eisenia_growth
### Keywords: datasets

### ** Examples

data(eisenia_growth)
head(eisenia_growth)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("eisenia_growth", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("eisenia_neuhauser")
### * eisenia_neuhauser

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: eisenia_neuhauser
### Title: Eisenia fetida Growth Data (Neuhauser 1980)
### Aliases: eisenia_neuhauser
### Keywords: datasets

### ** Examples

data(eisenia_neuhauser)
plot(eisenia_neuhauser$time, eisenia_neuhauser$length,
     xlab = "Days", ylab = "Structural length (cm)")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("eisenia_neuhauser", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fitted.bdeb_fit")
### * fitted.bdeb_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fitted.bdeb_fit
### Title: Fitted Values from a BDEB Fit
### Aliases: fitted.bdeb_fit

### ** Examples

## Not run: 
##D fit <- bdeb_fit(mod)
##D fitted(fit)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fitted.bdeb_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("folsomia_repro")
### * folsomia_repro

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: folsomia_repro
### Title: Simulated Folsomia candida Reproduction Data
### Aliases: folsomia_repro
### Keywords: datasets

### ** Examples

data(folsomia_repro)
head(folsomia_repro)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("folsomia_repro", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("logLik.bdeb_fit")
### * logLik.bdeb_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: logLik.bdeb_fit
### Title: Log-Likelihood of a BDEB Fit
### Aliases: logLik.bdeb_fit

### ** Examples

## Not run: 
##D fit <- bdeb_fit(mod)
##D logLik(fit)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("logLik.bdeb_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("nobs.bdeb_fit")
### * nobs.bdeb_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: nobs.bdeb_fit
### Title: Number of Observations Used in a BDEB Fit
### Aliases: nobs.bdeb_fit

### ** Examples

## Not run: 
##D fit <- bdeb_fit(mod)
##D nobs(fit)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("nobs.bdeb_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("observation_models")
### * observation_models

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: observation_models
### Title: Observation Model Specifications
### Aliases: observation_models obs_normal obs_lognormal obs_student_t
###   obs_poisson obs_negbinom

### ** Examples

obs_normal()
obs_lognormal()
obs_negbinom()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("observation_models", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.bdeb_data")
### * plot.bdeb_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.bdeb_data
### Title: Plot a BDEB Data Object
### Aliases: plot.bdeb_data

### ** Examples

df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
dat <- bdeb_data(growth = df)
plot(dat)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.bdeb_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.bdeb_diagnostics")
### * plot.bdeb_diagnostics

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.bdeb_diagnostics
### Title: Plot Convergence Diagnostics
### Aliases: plot.bdeb_diagnostics

### ** Examples

## Not run: 
##D fit <- bdeb_fit(mod)
##D plot(bdeb_diagnose(fit), type = "ess")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.bdeb_diagnostics", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.bdeb_fit")
### * plot.bdeb_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.bdeb_fit
### Title: Plot a BDEB Fit
### Aliases: plot.bdeb_fit

### ** Examples

## Not run: 
##D data(eisenia_growth)
##D dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
##D fit <- bdeb_fit(bdeb_model(dat, type = "individual"))
##D plot(fit, type = "trace")
##D plot(fit, type = "trajectory", n_draws = 50)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.bdeb_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.bdeb_model")
### * plot.bdeb_model

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.bdeb_model
### Title: Plot a BDEB Model Specification
### Aliases: plot.bdeb_model

### ** Examples

df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
dat <- bdeb_data(growth = df)
mod <- bdeb_model(dat, type = "individual")
plot(mod, pars = c("p_Am", "kappa"))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.bdeb_model", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.bdeb_ppc")
### * plot.bdeb_ppc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.bdeb_ppc
### Title: Plot Posterior Predictive Checks
### Aliases: plot.bdeb_ppc

### ** Examples

## Not run: 
##D fit <- bdeb_fit(mod)
##D plot(bdeb_ppc(fit, type = "growth"))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.bdeb_ppc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.bdeb_prediction")
### * plot.bdeb_prediction

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.bdeb_prediction
### Title: Plot Posterior Predictive Trajectories
### Aliases: plot.bdeb_prediction

### ** Examples

## Not run: 
##D fit <- bdeb_fit(mod)
##D plot(predict(fit), n_draws = 50)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.bdeb_prediction", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_dose_response")
### * plot_dose_response

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_dose_response
### Title: Plot DEBtox Dose-Response
### Aliases: plot_dose_response

### ** Examples

## Not run: 
##D data(debtox_growth)
##D dat <- bdeb_data(growth = debtox_growth,
##D                  concentration = c("1" = 0, "11" = 20,
##D                                   "21" = 80, "31" = 200))
##D fit <- bdeb_fit(bdeb_tox(dat, stress = "assimilation"))
##D plot_dose_response(fit, n_draws = 50)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_dose_response", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict.bdeb_fit")
### * predict.bdeb_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict.bdeb_fit
### Title: Predict from a BDEB Fit
### Aliases: predict.bdeb_fit bdeb_predict

### ** Examples

## Not run: 
##D data(eisenia_growth)
##D dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
##D fit <- bdeb_fit(bdeb_model(dat, type = "individual"))
##D pred <- predict(fit, n_draws = 200)
##D summary(pred); plot(pred)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict.bdeb_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.bdeb_diagnostics")
### * print.bdeb_diagnostics

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.bdeb_diagnostics
### Title: Print a BDEB Diagnostics Report
### Aliases: print.bdeb_diagnostics

### ** Examples

## Not run: 
##D fit <- bdeb_fit(mod)
##D print(bdeb_diagnose(fit))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.bdeb_diagnostics", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.bdeb_prediction")
### * print.bdeb_prediction

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.bdeb_prediction
### Title: Print a BDEB Prediction
### Aliases: print.bdeb_prediction

### ** Examples

## Not run: 
##D fit <- bdeb_fit(mod)
##D print(predict(fit))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.bdeb_prediction", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("prior_default")
### * prior_default

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: prior_default
### Title: Default Priors for DEB Parameters
### Aliases: prior_default

### ** Examples

prior_default("individual")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("prior_default", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("prior_species")
### * prior_species

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: prior_species
### Title: Species-Specific Priors from the AmP Collection
### Aliases: prior_species

### ** Examples

# Use AmP-calibrated priors for E. fetida
prior_species("Eisenia_fetida")

# Combine with model specification (R-side, no Stan sampling)
data(eisenia_growth)
dat <- bdeb_data(growth = eisenia_growth[eisenia_growth$id == 1, ])
mod <- bdeb_model(dat, type = "individual",
  priors = prior_species("Eisenia_fetida"))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("prior_species", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("priors")
### * priors

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: priors
### Title: Prior Distribution Specifications for BDEB Models
### Aliases: priors prior_lognormal prior_normal prior_beta
###   prior_halfnormal prior_halfcauchy prior_exponential

### ** Examples

# Log-normal prior on assimilation rate: median ~ exp(1.5) ~ 4.5
prior_lognormal(mu = 1.5, sigma = 0.5)

# Beta(2,2) prior on kappa — symmetric, favouring 0.5
prior_beta(a = 2, b = 2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("priors", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("repro_to_intervals")
### * repro_to_intervals

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: repro_to_intervals
### Title: Convert Cumulative Reproduction to Intervals
### Aliases: repro_to_intervals

### ** Examples

cumul <- data.frame(
  id = rep(1, 5),
  time = c(0, 7, 14, 21, 28),
  cumulative = c(0, 10, 30, 60, 100)
)
repro_to_intervals(cumul)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("repro_to_intervals", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("residuals.bdeb_fit")
### * residuals.bdeb_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: residuals.bdeb_fit
### Title: Residuals from a BDEB Fit
### Aliases: residuals.bdeb_fit

### ** Examples

## Not run: 
##D fit <- bdeb_fit(mod)
##D residuals(fit)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("residuals.bdeb_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.bdeb_data")
### * summary.bdeb_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.bdeb_data
### Title: Summary of a BDEB Data Object
### Aliases: summary.bdeb_data print.summary.bdeb_data

### ** Examples

df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
dat <- bdeb_data(growth = df)
summary(dat)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.bdeb_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.bdeb_diagnostics")
### * summary.bdeb_diagnostics

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.bdeb_diagnostics
### Title: Compact Summary of a BDEB Diagnostics Report
### Aliases: summary.bdeb_diagnostics

### ** Examples

## Not run: 
##D fit <- bdeb_fit(mod)
##D summary(bdeb_diagnose(fit))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.bdeb_diagnostics", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.bdeb_fit")
### * summary.bdeb_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.bdeb_fit
### Title: Posterior Summary for a BDEB Fit
### Aliases: summary.bdeb_fit

### ** Examples

## Not run: 
##D fit <- bdeb_fit(mod)
##D summary(fit, pars = c("p_Am", "kappa"), prob = 0.95)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.bdeb_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.bdeb_prediction")
### * summary.bdeb_prediction

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.bdeb_prediction
### Title: Tabulate Posterior Predictive Quantiles
### Aliases: summary.bdeb_prediction

### ** Examples

## Not run: 
##D fit <- bdeb_fit(mod)
##D summary(predict(fit), prob = 0.95)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.bdeb_prediction", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("vcov.bdeb_fit")
### * vcov.bdeb_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: vcov.bdeb_fit
### Title: Posterior Covariance Matrix of BDEB Model Parameters
### Aliases: vcov.bdeb_fit

### ** Examples

## Not run: 
##D fit <- bdeb_fit(mod)
##D vcov(fit)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("vcov.bdeb_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
