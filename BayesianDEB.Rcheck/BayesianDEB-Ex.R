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
nameEx("bdeb_fit")
### * bdeb_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bdeb_fit
### Title: Fit a BDEB Model via Hamiltonian Monte Carlo
### Aliases: bdeb_fit

### ** Examples

## Not run: 
##D dat <- bdeb_data(growth = growth_df)
##D mod <- bdeb_model(dat, type = "individual")
##D fit <- bdeb_fit(mod, chains = 2, iter_sampling = 500)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bdeb_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bdeb_prior_predictive")
### * bdeb_prior_predictive

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bdeb_prior_predictive
### Title: Prior Predictive Check
### Aliases: bdeb_prior_predictive

### ** Examples

## Not run: 
##D dat <- bdeb_data(growth = data.frame(id = 1, time = 0:10,
##D   length = seq(0.1, 0.5, length.out = 11)))
##D mod <- bdeb_model(dat, type = "individual")
##D pp <- bdeb_prior_predictive(mod, n_draws = 200)
##D plot(pp)
## End(Not run)



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
nameEx("bdeb_tox")
### * bdeb_tox

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bdeb_tox
### Title: DEBtox Model Specification
### Aliases: bdeb_tox

### ** Examples

## Not run: 
##D conc <- c("ctrl" = 0, "low" = 5, "mid" = 20, "high" = 100)
##D dat <- bdeb_data(growth = growth_df, concentration = conc)
##D mod <- bdeb_tox(dat, stress = "assimilation")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bdeb_tox", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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

# Combine with model specification
## Not run: 
##D mod <- bdeb_model(dat, type = "individual",
##D   priors = prior_species("Eisenia_fetida"))
## End(Not run)



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
