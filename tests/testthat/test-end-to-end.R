# ===========================================================
# End-to-end integration tests
# Full pipeline without Stan sampling, plus optional
# skip_if_not_installed("cmdstanr") mini-fit test
# ===========================================================


# --- E2E: bdeb_data -> bdeb_model -> stan_data validation ---

test_that("E2E individual: data -> model -> stan_data is complete", {
  data(eisenia_growth, package = "BayesianDEB", envir = environment())
  df1 <- eisenia_growth[eisenia_growth$id == 1, ]

  dat <- bdeb_data(growth = df1, f_food = 1.0)
  mod <- bdeb_model(dat, type = "individual",
    priors = list(p_Am = prior_lognormal(mu = 1.5, sigma = 0.5)))

  sd <- mod$stan_data
  expect_equal(sd$N_obs, 13)
  expect_equal(sd$f_food, 1.0)
  expect_equal(sd$obs_growth, 1L)  # normal default
  expect_equal(sd$has_temperature, 0L)
  expect_equal(sd$prior_p_Am_mu, 1.5)  # custom
  expect_equal(sd$prior_p_Am_sd, 0.5)  # custom
  expect_equal(sd$prior_kappa_a, 2)    # default
  expect_true(all(sd$t_obs >= 0))
  expect_true(all(sd$L_obs > 0))
})

test_that("E2E hierarchical: data -> model -> stan_data shapes", {
  data(eisenia_growth, package = "BayesianDEB", envir = environment())

  dat <- bdeb_data(growth = eisenia_growth, f_food = 1.0)
  mod <- bdeb_model(dat, type = "hierarchical")

  sd <- mod$stan_data
  expect_equal(sd$N_ind, 21)
  expect_equal(sd$max_N_obs, 13)
  expect_true(is.matrix(sd$t_obs))
  expect_equal(nrow(sd$t_obs), 21)
  expect_equal(ncol(sd$t_obs), 13)
  expect_true(is.matrix(sd$L_obs))
  expect_equal(sd$obs_growth, 1L)
  expect_equal(sd$has_temperature, 0L)
  expect_true("prior_mu_log_p_Am_mu" %in% names(sd))
  expect_true("prior_sigma_log_p_Am_rate" %in% names(sd))
})

test_that("E2E growth_repro: data -> model -> stan_data indices", {
  g <- data.frame(id = 1, time = seq(0, 56, by = 7),
                  length = c(0.10, 0.18, 0.28, 0.37, 0.44,
                             0.49, 0.52, 0.54, 0.55))
  r <- data.frame(id = 1, t_start = c(28, 35, 42, 49),
                  t_end = c(35, 42, 49, 56), count = c(5, 12, 18, 15))

  dat <- bdeb_data(growth = g, reproduction = r, f_food = 1.0)
  mod <- bdeb_model(dat, type = "growth_repro",
    observation = list(growth = obs_lognormal(),
                       reproduction = obs_poisson()))

  sd <- mod$stan_data
  expect_equal(sd$N_L, 9)
  expect_equal(sd$N_R, 4)
  expect_equal(sd$obs_growth, 2L)  # lognormal
  expect_equal(sd$obs_repro, 2L)   # poisson
  expect_true(sd$N_times >= max(sd$N_L, sd$N_R))
  # Index bounds
  expect_true(all(sd$idx_L >= 1 & sd$idx_L <= sd$N_times))
  expect_true(all(sd$idx_R_start >= 1 & sd$idx_R_start <= sd$N_times))
  expect_true(all(sd$idx_R_end >= 1 & sd$idx_R_end <= sd$N_times))
  expect_true(all(sd$idx_R_end > sd$idx_R_start))
})

test_that("E2E debtox: data -> model -> stan_data groups", {
  data(debtox_growth, package = "BayesianDEB", envir = environment())

  conc_levels <- unique(debtox_growth$concentration)
  conc_map <- setNames(conc_levels, as.character(conc_levels))

  dat <- bdeb_data(growth = debtox_growth, concentration = conc_map)
  mod <- bdeb_model(dat, type = "debtox",
    observation = list(growth = obs_student_t(nu = 5),
                       reproduction = obs_negbinom()))

  sd <- mod$stan_data
  expect_equal(sd$N_groups, 4)
  expect_equal(sort(sd$C_w), c(0, 20, 80, 200))
  expect_true(is.matrix(sd$t_obs))
  expect_true(is.matrix(sd$L_obs))
  expect_equal(sd$obs_growth, 3L)  # student_t
  expect_equal(sd$obs_nu, 5)
  expect_equal(sd$obs_repro, 1L)   # negbinom
  expect_equal(sd$has_repro, 0L)   # no repro data provided
  expect_true("prior_k_d_mu" %in% names(sd))
  expect_true("prior_z_w_mu" %in% names(sd))
})


# --- E2E: temperature correction flows through ---

test_that("E2E temperature: data -> model -> stan_data at 25C", {
  data(eisenia_growth, package = "BayesianDEB", envir = environment())
  df1 <- eisenia_growth[eisenia_growth$id == 1, ]

  dat <- bdeb_data(growth = df1)
  temp <- list(T_obs = 298.15, T_ref = 293.15, T_A = 8000)
  mod <- bdeb_model(dat, "individual", temperature = temp)

  sd <- mod$stan_data
  expect_equal(sd$has_temperature, 1L)
  expect_equal(sd$T_obs, 298.15)
  expect_equal(sd$T_ref, 293.15)
  expect_equal(sd$T_A, 8000)
})


# --- E2E: observation model flows through ---

test_that("E2E obs model: student_t(nu=3) flows into stan_data", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, "individual",
    observation = list(growth = obs_student_t(nu = 3)))

  expect_equal(mod$stan_data$obs_growth, 3L)
  expect_equal(mod$stan_data$obs_nu, 3)
})


# --- E2E: full mock pipeline per model type ---

test_that("E2E mock individual: fit -> diagnose -> summary -> derived -> ppc -> plot", {
  fit <- mock_bdeb_fit(n_draws = 100, type = "individual")

  diag <- suppressWarnings(bdeb_diagnose(fit))
  expect_equal(diag$n_divergent, 0)

  s <- suppressWarnings(bdeb_summary(fit, pars = c("p_Am", "kappa")))
  expect_equal(nrow(s), 2)

  d <- bdeb_derived(fit, quantities = c("L_m", "L_inf", "k_M", "g", "growth_rate"))
  expect_true(all(c("L_m", "L_inf", "k_M", "g", "growth_rate") %in% names(d)))

  ppc <- bdeb_ppc(fit, type = "growth")
  expect_s3_class(ppc, "bdeb_ppc")

  p1 <- plot(fit, type = "trace", pars = c("p_Am", "kappa"))
  expect_s3_class(p1, "ggplot")

  p2 <- plot(fit, type = "trajectory", n_draws = 10)
  expect_s3_class(p2, "ggplot")
})

test_that("E2E mock hierarchical: fit -> diagnose -> derived -> trajectory", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "hierarchical")

  diag <- suppressWarnings(bdeb_diagnose(fit))
  expect_type(diag, "list")

  d <- bdeb_derived(fit, quantities = c("L_inf", "g"))
  expect_true(all(d$L_inf > 0))

  p <- plot(fit, type = "trajectory", n_draws = 5)
  expect_s3_class(p, "ggplot")
  expect_true("FacetWrap" %in% class(p$facet))
})

test_that("E2E mock debtox: fit -> diagnose -> ec50 -> dose-response -> trajectory", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "debtox")

  diag <- suppressWarnings(bdeb_diagnose(fit))
  expect_equal(diag$n_divergent, 0)

  ec <- bdeb_ec50(fit, prob = 0.90)
  expect_true(all(ec$draws > 0))
  expect_gt(ec$summary$mean[1], ec$summary$mean[2])  # EC50 > NEC

  p1 <- plot_dose_response(fit, n_draws = 5)
  expect_s3_class(p1, "ggplot")

  p2 <- plot(fit, type = "trajectory", n_draws = 5)
  expect_s3_class(p2, "ggplot")
  expect_true("FacetWrap" %in% class(p2$facet))
})


# --- E2E: optional CmdStan compilation test ---

test_that("Stan models compile (requires cmdstanr + CmdStan)", {
  skip_if_not_installed("cmdstanr")
  skip_if(
    is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
    "CmdStan not installed"
  )

  models <- c("bdeb_individual_growth", "bdeb_growth_repro",
              "bdeb_hierarchical_growth", "bdeb_debtox")

  tmp <- tempdir()
  for (m in models) {
    path <- stan_file(m)
    # Compile to temp dir to avoid polluting inst/stan with binaries
    tmp_stan <- file.path(tmp, basename(path))
    file.copy(path, tmp_stan, overwrite = TRUE)
    mod <- cmdstanr::cmdstan_model(tmp_stan, compile = TRUE, quiet = TRUE)
    expect_true(file.exists(mod$exe_file()), info = m)
  }
})
