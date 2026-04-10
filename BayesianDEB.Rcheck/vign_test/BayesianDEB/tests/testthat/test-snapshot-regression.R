# ===========================================================
# C. Snapshot / regression tests
# Fixed inputs -> fixed outputs. Guards against silent API drift.
# ===========================================================

# --- Fixed dataset ---

test_that("eisenia_growth has expected dimensions and columns", {
  data(eisenia_growth, package = "BayesianDEB", envir = environment())
  expect_equal(nrow(eisenia_growth), 273)
  expect_equal(ncol(eisenia_growth), 3)
  expect_named(eisenia_growth, c("id", "time", "length"))
  expect_equal(length(unique(eisenia_growth$id)), 21)
  expect_equal(sort(unique(eisenia_growth$time)),
               seq(0, 84, by = 7))
})

test_that("folsomia_repro has expected dimensions and columns", {
  data(folsomia_repro, package = "BayesianDEB", envir = environment())
  expect_equal(nrow(folsomia_repro), 30)
  expect_named(folsomia_repro,
               c("id", "concentration", "t_start", "t_end", "count"))
  expect_equal(sort(unique(folsomia_repro$concentration)),
               c(0, 10, 50, 100, 500))
})

test_that("debtox_growth has expected dimensions and columns", {
  data(debtox_growth, package = "BayesianDEB", envir = environment())
  expect_equal(nrow(debtox_growth), 280)
  expect_named(debtox_growth, c("id", "concentration", "time", "length"))
  expect_equal(sort(unique(debtox_growth$concentration)),
               c(0, 20, 80, 200))
  expect_equal(length(unique(debtox_growth$id)), 40)
})


# --- Fixed bdeb_data output ---

test_that("bdeb_data on eisenia id=1 produces stable structure", {
  data(eisenia_growth, package = "BayesianDEB", envir = environment())
  df1 <- eisenia_growth[eisenia_growth$id == 1, ]
  dat <- bdeb_data(growth = df1, f_food = 1.0)

  expect_equal(dat$n_ind, 1)
  expect_equal(dat$endpoints, "growth")
  expect_equal(dat$f_food, 1.0)
  expect_equal(nrow(dat$growth), 13)
  expect_equal(dat$growth$time, seq(0, 84, by = 7))
})


# --- Fixed prior_default output ---

test_that("prior_default('individual') is stable", {
  d <- prior_default("individual")
  expect_equal(length(d), 8)
  expect_named(d, c("p_Am", "p_M", "kappa", "v", "E_G",
                     "E0", "L0", "sigma_L"))
  expect_equal(d$p_Am$family, "lognormal")
  expect_equal(d$p_Am$mu, 1.0)
  expect_equal(d$p_Am$sigma, 1.0)
  expect_equal(d$kappa$family, "beta")
  expect_equal(d$kappa$a, 2)
  expect_equal(d$kappa$b, 2)
  expect_equal(d$sigma_L$family, "halfnormal")
  expect_equal(d$sigma_L$sigma, 0.1)
})

test_that("prior_default('debtox') adds exactly k_d, z_w, b_w, k_R, phi_R", {
  d <- prior_default("debtox")
  extra <- setdiff(names(d), names(prior_default("individual")))
  expect_equal(sort(extra), c("b_w", "k_R", "k_d", "phi_R", "z_w"))
})

test_that("prior_default('hierarchical') adds exactly mu/sigma_log_p_Am", {
  d <- prior_default("hierarchical")
  extra <- setdiff(names(d), names(prior_default("individual")))
  expect_equal(sort(extra), c("mu_log_p_Am", "sigma_log_p_Am"))
})


# --- Fixed mock posterior: deterministic outputs ---

test_that("mock_bdeb_fit produces deterministic draws (seed = 42)", {
  fit1 <- mock_bdeb_fit(n_draws = 50, type = "individual")
  fit2 <- mock_bdeb_fit(n_draws = 50, type = "individual")

  d1 <- posterior::as_draws_df(fit1$fit$draws())
  d2 <- posterior::as_draws_df(fit2$fit$draws())

  # Same seed => identical draws
  expect_equal(as.numeric(d1$p_Am), as.numeric(d2$p_Am))
  expect_equal(as.numeric(d1$kappa), as.numeric(d2$kappa))
})

test_that("bdeb_summary on fixed mock gives reproducible means", {
  fit <- mock_bdeb_fit(n_draws = 200, type = "individual")
  s <- suppressWarnings(bdeb_summary(fit, pars = c("p_Am", "kappa")))
  s_df <- as.data.frame(s)

  # p_Am: mock draws are rlnorm(., 1.5, 0.1) => mean ~ exp(1.5 + 0.01/2) ~ 4.5
  expect_gt(s_df$mean[s_df$variable == "p_Am"], 3.5)
  expect_lt(s_df$mean[s_df$variable == "p_Am"], 5.5)

  # kappa: mock draws are rbeta(., 8, 3) => mean ~ 8/11 ~ 0.727
  expect_gt(s_df$mean[s_df$variable == "kappa"], 0.6)
  expect_lt(s_df$mean[s_df$variable == "kappa"], 0.85)
})

test_that("bdeb_derived on fixed mock gives reproducible L_inf", {
  fit <- mock_bdeb_fit(n_draws = 200, type = "individual")
  d <- bdeb_derived(fit, quantities = "L_inf", f = 1.0)

  # L_inf = kappa * p_Am / p_M
  # kappa ~ 0.73, p_Am ~ 4.5, p_M ~ 0.50 => L_inf ~ 6.6
  mean_L <- mean(d$L_inf)
  expect_gt(mean_L, 4)
  expect_lt(mean_L, 10)
})

test_that("bdeb_ec50 on fixed mock gives reproducible range", {
  fit <- mock_bdeb_fit(n_draws = 200, type = "debtox")
  ec <- bdeb_ec50(fit, prob = 0.90)

  # z_w ~ exp(2.7) ~ 14.9, b_w ~ exp(-5.5) ~ 0.004
  # EC50 ~ 14.9 + 0.5/0.004 ~ 140
  expect_gt(ec$summary$mean[1], 50)
  expect_lt(ec$summary$mean[1], 300)

  # NEC ~ 14.9
  expect_gt(ec$summary$mean[2], 5)
  expect_lt(ec$summary$mean[2], 30)
})


# --- Stan data structure regression ---

test_that("build_stan_data_individual has stable field names", {
  df <- data.frame(id = 1, time = c(0, 7, 14), length = c(0.1, 0.2, 0.3))
  dat <- bdeb_data(growth = df)
  sd <- build_stan_data_individual(dat, prior_default("individual"))

  expected_fields <- c("N_obs", "t_obs", "L_obs", "f_food",
                       "has_temperature", "T_obs", "T_ref", "T_A",
                       "prior_p_Am_mu", "prior_p_Am_sd",
                       "prior_kappa_a", "prior_kappa_b",
                       "prior_sigma_L_sd")
  for (fld in expected_fields) {
    expect_true(fld %in% names(sd), info = paste("missing:", fld))
  }
})

test_that("bdeb_model('individual') stan_data structure is stable", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, "individual")

  expect_equal(mod$stan_data$N_obs, 6)
  expect_equal(mod$stan_data$f_food, 1.0)
  expect_true(is.numeric(mod$stan_data$prior_p_Am_mu))
  expect_true(is.numeric(mod$stan_data$prior_kappa_a))
})


# --- arrhenius: known numerical values ---

test_that("arrhenius at 25C/20C/8000 gives known value", {
  val <- arrhenius(298.15, T_ref = 293.15, T_A = 8000)
  # exp(8000/293.15 - 8000/298.15) = exp(27.295 - 26.829) = exp(0.466) ~ 1.594
  expect_equal(val, exp(8000/293.15 - 8000/298.15), tolerance = 1e-10)
  expect_gt(val, 1.5)
  expect_lt(val, 1.7)
})


# --- deb_fluxes: known numerical values ---

test_that("deb_fluxes with known inputs gives known L", {
  fl <- deb_fluxes(E = 10, V = 0.125, f = 1.0,
                   p_Am = 5, p_M = 0.5, kappa = 0.75,
                   v = 0.2, E_G = 400)
  # V = 0.125 => L = 0.5
  expect_equal(fl$L, 0.5)
  # p_A = 1.0 * 5 * 0.5^2 = 1.25
  expect_equal(fl$p_A, 1.25)
  # p_M = 0.5 * 0.125 = 0.0625
  expect_equal(fl$p_M, 0.0625)
})
