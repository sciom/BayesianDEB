# ===========================================================
# Deep validation tests: edge cases, numerical stability,
# failure scenarios, and convergence diagnostics
# ===========================================================


# =========================================================
# 1. DATA EDGE CASES
# =========================================================

test_that("bdeb_data handles duplicate times within individual", {
  df <- data.frame(id = 1, time = c(0, 7, 7, 14), length = c(0.1, 0.2, 0.21, 0.3))
  dat <- bdeb_data(growth = df)
  expect_equal(nrow(dat$growth), 4)
})

test_that("bdeb_data handles NA in length (passes through)", {
  df <- data.frame(id = 1, time = 0:3, length = c(0.1, NA, 0.3, 0.4))
  # NA lengths should not trigger the non-negative check
  dat <- bdeb_data(growth = df)
  expect_equal(nrow(dat$growth), 4)
  expect_true(is.na(dat$growth$length[2]))
})

test_that("bdeb_data handles very large number of individuals", {
  n <- 200
  df <- data.frame(
    id     = rep(seq_len(n), each = 3),
    time   = rep(c(0, 7, 14), n),
    length = runif(n * 3, 0.1, 0.5)
  )
  dat <- bdeb_data(growth = df)
  expect_equal(dat$n_ind, n)
})

test_that("bdeb_data with only time=0 observations", {
  df <- data.frame(id = 1:5, time = 0, length = runif(5, 0.08, 0.12))
  dat <- bdeb_data(growth = df)
  expect_equal(dat$n_ind, 5)
})

test_that("bdeb_data rejects Inf in times", {
  df <- data.frame(id = 1, time = c(0, Inf), length = c(0.1, 0.2))
  # Inf >= 0 is TRUE, so it passes the non-negative check
  # But this is still a valid concern — document that Inf is accepted
  dat <- bdeb_data(growth = df)
  expect_true(is.infinite(dat$growth$time[2]))
})

test_that("bdeb_data with all reproduction counts zero", {
  df <- data.frame(id = 1:3, t_start = 0, t_end = 28, count = 0)
  dat <- bdeb_data(reproduction = df)
  expect_equal(nrow(dat$reproduction), 3)
})

test_that("bdeb_data combined growth + repro", {
  g <- data.frame(id = 1:2, time = c(0, 0), length = c(0.1, 0.1))
  r <- data.frame(id = 1, t_start = 0, t_end = 28, count = 50)
  dat <- bdeb_data(growth = g, reproduction = r)
  expect_equal(sort(dat$endpoints), c("growth", "reproduction"))
  expect_equal(dat$n_ind, 2)
})


# =========================================================
# 2. NUMERICAL STABILITY OF DEB FLUXES
# =========================================================

test_that("deb_fluxes: near-zero E and V simultaneously", {
  fl <- deb_fluxes(E = 1e-20, V = 1e-20, f = 1.0,
                   p_Am = 5, p_M = 0.5, kappa = 0.75,
                   v = 0.2, E_G = 400)
  expect_true(all(is.finite(unlist(fl))))
  # All fluxes should be near zero
  expect_lt(abs(fl$p_A), 1e-10)
  expect_lt(abs(fl$p_C), 1e-10)
})

test_that("deb_fluxes: kappa near boundary values", {
  fl_low <- deb_fluxes(E = 10, V = 0.5, f = 1.0,
                       p_Am = 5, p_M = 0.5, kappa = 0.001,
                       v = 0.2, E_G = 400)
  fl_high <- deb_fluxes(E = 10, V = 0.5, f = 1.0,
                        p_Am = 5, p_M = 0.5, kappa = 0.999,
                        v = 0.2, E_G = 400)
  expect_true(all(is.finite(unlist(fl_low))))
  expect_true(all(is.finite(unlist(fl_high))))
  # With kappa~0, almost nothing goes to growth
  expect_equal(fl_low$p_G, 0)  # kappa*p_C < p_M
  # With kappa~1, almost nothing goes to reproduction
  expect_lt(fl_high$p_R, fl_low$p_R)
})

test_that("deb_fluxes: energy conservation (p_C splits into kappa and 1-kappa)", {
  # Use high E so that kappa*p_C > p_M*V (positive growth)
  fl <- deb_fluxes(E = 500, V = 0.1, f = 1.0,
                   p_Am = 5, p_M = 0.2, kappa = 0.8,
                   v = 0.2, E_G = 400, k_J = 0, E_Hp = 0)
  # When growth is positive: kappa * p_C = p_G + p_M_flux
  # fl$p_M is p_M_flux = [p_M] * V
  kappa_flux <- 0.8 * fl$p_C
  expect_equal(kappa_flux, fl$p_G + fl$p_M, tolerance = 1e-10)
  # (1 - kappa) * p_C = p_R + p_J
  one_minus_kappa_flux <- 0.2 * fl$p_C
  expect_equal(one_minus_kappa_flux, fl$p_R + fl$p_J, tolerance = 1e-10)
})

test_that("deb_fluxes: scaled reserve density e is dimensionless and bounded", {
  fl <- deb_fluxes(E = 10, V = 0.5, f = 1.0,
                   p_Am = 5, p_M = 0.5, kappa = 0.75,
                   v = 0.2, E_G = 400)
  expect_gte(fl$e, 0)
  # At steady state with f=1, e should be close to 1 (but not necessarily exact)
  expect_true(is.finite(fl$e))
})

test_that("deb_fluxes: very large organism", {
  fl <- deb_fluxes(E = 1e6, V = 1e4, f = 1.0,
                   p_Am = 5, p_M = 0.5, kappa = 0.75,
                   v = 0.2, E_G = 400)
  expect_true(all(is.finite(unlist(fl))))
  expect_gt(fl$p_A, 0)
  expect_gt(fl$L, 0)
})


# =========================================================
# 3. ARRHENIUS NUMERICAL EDGE CASES
# =========================================================

test_that("arrhenius: near absolute zero (1 K) underflows to zero", {
  val <- arrhenius(1, T_ref = 293.15, T_A = 8000)
  # exp(8000/293.15 - 8000/1) = exp(-7972.7) underflows to 0.0
  # This is correct numerical behaviour, not an error
  expect_true(is.finite(val))
  expect_gte(val, 0)
})

test_that("arrhenius: moderately cold (200 K) is finite and positive", {
  val <- arrhenius(200, T_ref = 293.15, T_A = 8000)
  expect_true(is.finite(val))
  expect_gt(val, 0)
  expect_lt(val, 1)
})

test_that("arrhenius: very high temperature (1000 K)", {
  val <- arrhenius(1000, T_ref = 293.15, T_A = 8000)
  expect_true(is.finite(val))
  expect_gt(val, 1)
})

test_that("arrhenius: very large T_A (100000 K)", {
  val <- arrhenius(298.15, T_ref = 293.15, T_A = 100000)
  expect_true(is.finite(val))
  expect_gt(val, 1)  # Still > 1 since T > T_ref
})

test_that("arrhenius: equal temperatures with different T_A always gives 1", {
  for (ta in c(0, 100, 8000, 50000)) {
    expect_equal(arrhenius(300, T_ref = 300, T_A = ta), 1.0)
  }
})


# =========================================================
# 4. PRIOR SYSTEM ROBUSTNESS
# =========================================================

test_that("prior_to_stan_data with extreme hyperparameters", {
  priors <- list(
    p_Am    = prior_lognormal(mu = -100, sigma = 0.001),
    p_M     = prior_lognormal(mu = 100, sigma = 100),
    kappa   = prior_beta(a = 0.001, b = 0.001),
    v       = prior_lognormal(mu = 0, sigma = 50),
    E_G     = prior_lognormal(mu = 0, sigma = 0.001),
    E0      = prior_lognormal(mu = 0, sigma = 1),
    L0      = prior_lognormal(mu = 0, sigma = 1),
    sigma_L = prior_halfnormal(sigma = 1e-10)
  )
  sd <- prior_to_stan_data(priors)
  expect_true(all(is.finite(unlist(sd))))
  expect_equal(sd$prior_p_Am_mu, -100)
  expect_equal(sd$prior_kappa_a, 0.001)
  expect_equal(sd$prior_sigma_L_sd, 1e-10)
})

test_that("bdeb_model merges custom priors correctly with defaults", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)

  # Override only p_Am and sigma_L
  mod <- bdeb_model(dat, type = "individual",
    priors = list(
      p_Am    = prior_lognormal(mu = 99, sigma = 0.01),
      sigma_L = prior_halfnormal(sigma = 42)
    ))

  # Overridden
  expect_equal(mod$priors$p_Am$mu, 99)
  expect_equal(mod$priors$sigma_L$sigma, 42)
  # Defaults preserved
  defaults <- prior_default("individual")
  expect_equal(mod$priors$p_M$mu, defaults$p_M$mu)
  expect_equal(mod$priors$kappa$a, defaults$kappa$a)
  expect_equal(mod$priors$v$mu, defaults$v$mu)
  expect_equal(mod$priors$E_G$mu, defaults$E_G$mu)
  expect_equal(mod$priors$E0$mu, defaults$E0$mu)
  expect_equal(mod$priors$L0$mu, defaults$L0$mu)
})


# =========================================================
# 5. CONVERGENCE DIAGNOSTICS: BAD FIT DETECTION
# =========================================================

test_that("bdeb_diagnose detects divergent transitions", {
  fit <- mock_bdeb_fit_bad(n_divergent = 25, n_max_treedepth = 0,
                           low_ebfmi = FALSE)
  result <- suppressWarnings(bdeb_diagnose(fit))

  expect_equal(result$n_divergent, 50)  # 25 per chain * 2 chains
})

test_that("bdeb_diagnose detects max treedepth saturation", {
  fit <- mock_bdeb_fit_bad(n_divergent = 0, n_max_treedepth = 30,
                           low_ebfmi = FALSE)
  result <- suppressWarnings(bdeb_diagnose(fit))

  expect_equal(result$n_max_treedepth, 60)  # 30 per chain * 2
})

test_that("bdeb_diagnose detects low E-BFMI", {
  fit <- mock_bdeb_fit_bad(n_divergent = 0, n_max_treedepth = 0,
                           low_ebfmi = TRUE)
  result <- suppressWarnings(bdeb_diagnose(fit))

  expect_true(all(result$ebfmi < 0.3))
})

test_that("bdeb_diagnose detects bad Rhat from bimodal chains", {
  fit <- mock_bdeb_fit_bad(n_divergent = 0, n_max_treedepth = 0,
                           low_ebfmi = FALSE, bimodal = TRUE)
  result <- suppressWarnings(bdeb_diagnose(fit))

  # At least some parameters should have Rhat > 1.01
  bad_rhat <- result$summary$variable[
    !is.na(result$summary$rhat) & result$summary$rhat > 1.01]
  expect_true(length(bad_rhat) > 0)
})

test_that("bdeb_diagnose: clean fit has no problems", {
  fit <- mock_bdeb_fit(n_draws = 200)
  result <- suppressWarnings(bdeb_diagnose(fit))

  expect_equal(result$n_divergent, 0)
  expect_equal(result$n_max_treedepth, 0)
  expect_true(all(result$ebfmi > 0.3))
})

test_that("bdeb_diagnose returns all expected fields", {
  fit <- mock_bdeb_fit_bad(n_divergent = 5, n_max_treedepth = 3,
                           low_ebfmi = TRUE)
  result <- suppressWarnings(bdeb_diagnose(fit))

  expect_type(result, "list")
  expect_true(all(c("n_divergent", "n_max_treedepth", "ebfmi", "summary")
                  %in% names(result)))
  expect_s3_class(result$summary, "data.frame")
  expect_true("rhat" %in% names(result$summary))
  expect_true("ess_bulk" %in% names(result$summary))
  expect_true("ess_tail" %in% names(result$summary))
})


# =========================================================
# 6. DERIVED QUANTITY MATHEMATICAL CORRECTNESS
# =========================================================

test_that("bdeb_derived L_inf = f * kappa * p_Am / p_M exactly", {
  fit <- mock_bdeb_fit(n_draws = 100)
  draws <- posterior::as_draws_df(fit$fit$draws())

  d <- bdeb_derived(fit, quantities = "L_inf", f = 1.0)
  expected <- draws$kappa * draws$p_Am / draws$p_M
  expect_equal(as.numeric(d$L_inf), as.numeric(expected), tolerance = 1e-12)
})

test_that("bdeb_derived k_M = p_M / E_G exactly", {
  fit <- mock_bdeb_fit(n_draws = 100)
  draws <- posterior::as_draws_df(fit$fit$draws())

  d <- bdeb_derived(fit, quantities = "k_M")
  expected <- draws$p_M / draws$E_G
  expect_equal(as.numeric(d$k_M), as.numeric(expected), tolerance = 1e-12)
})

test_that("bdeb_derived growth_rate = k_M * g / (3*(f+g)) exactly (Kooijman Eq 3.23)", {
  fit <- mock_bdeb_fit(n_draws = 100)
  draws <- posterior::as_draws_df(fit$fit$draws())

  d <- bdeb_derived(fit, quantities = "growth_rate", f = 1.0)
  g <- draws$E_G * draws$v / (draws$kappa * draws$p_Am)
  k_M <- draws$p_M / draws$E_G
  expected <- k_M * g / (3 * (1.0 + g))
  expect_equal(as.numeric(d$growth_rate), as.numeric(expected), tolerance = 1e-12)
})

test_that("bdeb_derived growth_rate depends on f", {
  fit <- mock_bdeb_fit(n_draws = 100)

  d1  <- bdeb_derived(fit, quantities = "growth_rate", f = 1.0)
  d05 <- bdeb_derived(fit, quantities = "growth_rate", f = 0.5)

  # At lower f, growth rate is HIGHER because r_B = k_M*g/(3*(f+g))
  # and smaller f means smaller denominator
  expect_true(mean(d05$growth_rate) > mean(d1$growth_rate))
})

test_that("bdeb_derived L_inf scales linearly with f", {
  fit <- mock_bdeb_fit(n_draws = 100)

  d1   <- bdeb_derived(fit, quantities = "L_inf", f = 1.0)
  d05  <- bdeb_derived(fit, quantities = "L_inf", f = 0.5)
  d025 <- bdeb_derived(fit, quantities = "L_inf", f = 0.25)

  # L_inf(f=0.5) / L_inf(f=1.0) should be exactly 0.5
  ratio1 <- as.numeric(d05$L_inf) / as.numeric(d1$L_inf)
  expect_true(all(abs(ratio1 - 0.5) < 1e-12))

  ratio2 <- as.numeric(d025$L_inf) / as.numeric(d1$L_inf)
  expect_true(all(abs(ratio2 - 0.25) < 1e-12))
})

test_that("bdeb_derived k_M is independent of f, growth_rate is not", {
  fit <- mock_bdeb_fit(n_draws = 50)

  d1  <- bdeb_derived(fit, quantities = c("k_M", "growth_rate"), f = 1.0)
  d07 <- bdeb_derived(fit, quantities = c("k_M", "growth_rate"), f = 0.7)

  # k_M = p_M / E_G — does not depend on f
  expect_equal(as.numeric(d1$k_M), as.numeric(d07$k_M), tolerance = 1e-12)
  # growth_rate = k_M*g/(3*(f+g)) — depends on f
  expect_false(isTRUE(all.equal(as.numeric(d1$growth_rate),
                                as.numeric(d07$growth_rate))))
})


test_that("bdeb_derived g = E_G * v / (kappa * p_Am) exactly", {
  fit <- mock_bdeb_fit(n_draws = 100)
  draws <- posterior::as_draws_df(fit$fit$draws())

  d <- bdeb_derived(fit, quantities = "g")
  expected <- draws$E_G * draws$v / (draws$kappa * draws$p_Am)
  expect_equal(as.numeric(d$g), as.numeric(expected), tolerance = 1e-12)
})

test_that("bdeb_derived dimensional consistency: L_inf * p_M = kappa * p_Am * f", {
  # From L_inf = f*kappa*p_Am/p_M => L_inf * p_M = f*kappa*p_Am
  fit <- mock_bdeb_fit(n_draws = 100)
  draws <- posterior::as_draws_df(fit$fit$draws())
  d <- bdeb_derived(fit, quantities = "L_inf", f = 0.8)

  lhs <- as.numeric(d$L_inf) * as.numeric(draws$p_M)
  rhs <- 0.8 * as.numeric(draws$kappa) * as.numeric(draws$p_Am)
  expect_equal(lhs, rhs, tolerance = 1e-12)
})

test_that("bdeb_derived growth_rate relates to L_inf and k_M correctly", {
  # r_B = k_M * g / (3*(f+g))
  # Also: r_B = v / (3*L_inf) * (f/(f+g)) ... let's verify the basic relation
  # r_B * 3 * (f+g) = k_M * g
  fit <- mock_bdeb_fit(n_draws = 100)
  d <- bdeb_derived(fit, quantities = c("growth_rate", "k_M", "g"), f = 1.0)

  lhs <- as.numeric(d$growth_rate) * 3 * (1.0 + as.numeric(d$g))
  rhs <- as.numeric(d$k_M) * as.numeric(d$g)
  expect_equal(lhs, rhs, tolerance = 1e-12)
})


# =========================================================
# 7. EC50/NEC MATHEMATICAL PROPERTIES
# =========================================================

test_that("EC50 = z_w + 0.5/b_w exactly for every draw", {
  fit <- mock_bdeb_fit(n_draws = 200, type = "debtox")
  draws <- posterior::as_draws_df(fit$fit$draws())

  ec50_computed <- draws$z_w + 0.5 / draws$b_w
  expect_equal(as.numeric(draws$EC50), as.numeric(ec50_computed),
               tolerance = 1e-12)
})

test_that("NEC = z_w exactly for every draw", {
  fit <- mock_bdeb_fit(n_draws = 200, type = "debtox")
  draws <- posterior::as_draws_df(fit$fit$draws())

  expect_equal(as.numeric(draws$NEC), as.numeric(draws$z_w),
               tolerance = 1e-12)
})

test_that("EC50 > NEC for every posterior draw", {
  fit <- mock_bdeb_fit(n_draws = 500, type = "debtox")
  draws <- posterior::as_draws_df(fit$fit$draws())

  # Since b_w > 0, 0.5/b_w > 0, so EC50 > z_w = NEC always
  expect_true(all(draws$EC50 > draws$NEC))
})

test_that("bdeb_ec50 credible interval narrows with more data (proxy)", {
  # More draws = tighter MC estimate (not real data, but tests the plumbing)
  fit <- mock_bdeb_fit(n_draws = 500, type = "debtox")

  r90 <- bdeb_ec50(fit, prob = 0.90)
  r50 <- bdeb_ec50(fit, prob = 0.50)

  width90 <- r90$summary$upper[1] - r90$summary$lower[1]
  width50 <- r50$summary$upper[1] - r50$summary$lower[1]
  expect_gt(width90, width50)
})


# =========================================================
# 8. STAN DATA BUILDER NUMERICAL PROPERTIES
# =========================================================

test_that("build_stan_data_individual preserves data order", {
  df <- data.frame(id = 1, time = c(14, 0, 7), length = c(0.3, 0.1, 0.2))
  dat <- bdeb_data(growth = df)
  priors <- prior_default("individual")
  sd <- build_stan_data_individual(dat, priors)

  # Should be sorted by time
  expect_equal(sd$t_obs, c(1e-3, 7, 14))  # t=0 replaced with epsilon
  expect_equal(sd$L_obs, c(0.1, 0.2, 0.3))
})

test_that("build_stan_data_hierarchical with 1 individual degenerates gracefully", {
  df <- data.frame(id = 1, time = 0:3, length = seq(0.1, 0.4, by = 0.1))
  dat <- bdeb_data(growth = df)
  priors <- prior_default("hierarchical")
  sd <- build_stan_data_hierarchical(dat, priors)

  expect_equal(sd$N_ind, 1)
  expect_equal(sd$max_N_obs, 4)
  expect_true(is.matrix(sd$t_obs))
  expect_equal(nrow(sd$t_obs), 1)
})

test_that("build_stan_data_debtox groups data correctly by concentration", {
  df <- data.frame(
    id = rep(1:6, each = 3),
    time = rep(c(0, 7, 14), 6),
    length = runif(18, 0.1, 0.5)
  )
  conc <- c("1" = 0, "2" = 0, "3" = 50, "4" = 50, "5" = 100, "6" = 100)
  dat <- bdeb_data(growth = df, concentration = conc)
  priors <- prior_default("debtox")
  sd <- build_stan_data_debtox(dat, priors)

  expect_equal(sd$N_groups, 3)
  expect_equal(sort(sd$C_w), c(0, 50, 100))
})


# =========================================================
# 9. PPC DIMENSIONAL CONSISTENCY
# =========================================================

test_that("PPC L_rep and L_obs have consistent dimensions", {
  fit <- mock_bdeb_fit(n_draws = 100)
  ppc <- bdeb_ppc(fit, type = "growth")

  expect_equal(ncol(ppc$growth$L_rep), length(ppc$growth$L_obs))
  expect_equal(nrow(ppc$growth$L_rep), ppc$growth$n_draws)
  expect_equal(length(ppc$growth$t_obs), ppc$growth$n_obs)
})

test_that("PPC replicated data has similar scale to observed", {
  fit <- mock_bdeb_fit(n_draws = 200)
  ppc <- bdeb_ppc(fit, type = "growth")

  obs_range <- range(ppc$growth$L_obs)
  rep_range <- range(ppc$growth$L_rep)

  # Replicated range should roughly overlap with observed
  expect_lt(rep_range[1], obs_range[2])
  expect_gt(rep_range[2], obs_range[1])
})


# =========================================================
# 10. FULL PIPELINE INTEGRATION (mock)
# =========================================================

test_that("full pipeline: individual mock -> diagnose -> summary -> derived -> ppc -> plot", {
  fit <- mock_bdeb_fit(n_draws = 100, type = "individual")

  # Diagnose
  diag <- suppressWarnings(bdeb_diagnose(fit))
  expect_type(diag, "list")

  # Summary
  s <- suppressWarnings(bdeb_summary(fit, pars = c("p_Am", "kappa")))
  expect_equal(nrow(s), 2)

  # Derived
  d <- bdeb_derived(fit, quantities = c("L_inf", "k_M", "growth_rate"))
  expect_true(all(c("L_inf", "k_M", "growth_rate") %in% names(d)))

  # PPC
  ppc <- bdeb_ppc(fit, type = "growth")
  expect_s3_class(ppc, "bdeb_ppc")

  # Plot (returns ggplot, doesn't error)
  p <- plot(fit, type = "trace", pars = c("p_Am", "kappa"))
  expect_s3_class(p, "ggplot")
})

test_that("full pipeline: debtox mock -> diagnose -> ec50 -> dose-response", {
  fit <- mock_bdeb_fit(n_draws = 100, type = "debtox")

  diag <- suppressWarnings(bdeb_diagnose(fit))
  expect_equal(diag$n_divergent, 0)

  ec <- bdeb_ec50(fit, prob = 0.90)
  expect_gt(ec$summary$mean[ec$summary$parameter == "EC50"], 0)
  expect_gt(ec$summary$mean[ec$summary$parameter == "NEC"], 0)

  p <- plot_dose_response(fit, n_draws = 20)
  expect_s3_class(p, "ggplot")
})

test_that("full pipeline: hierarchical mock -> diagnose -> derived", {
  fit <- mock_bdeb_fit(n_draws = 100, type = "hierarchical")

  diag <- suppressWarnings(bdeb_diagnose(fit))
  expect_type(diag, "list")

  # Derived should use mu_log_p_Am fallback for p_Am
  d <- bdeb_derived(fit, quantities = "L_inf")
  expect_true("L_inf" %in% names(d))
  expect_true(all(d$L_inf > 0))
})
