# ===========================================================
# A. API contract tests
# Explicit tests for documented behaviour at user-facing boundaries
# ===========================================================


# --- bdeb_model: type vs data contracts ---

test_that("individual model with >1 individual errors immediately", {
  df <- data.frame(id = rep(1:3, each = 5), time = rep(0:4, 3),
                   length = runif(15, 0.1, 0.5))
  dat <- bdeb_data(growth = df)
  expect_error(bdeb_model(dat, type = "individual"),
               "requires exactly 1 individual")
})

test_that("growth_repro without reproduction data errors", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  expect_error(bdeb_model(dat, type = "growth_repro"), "reproduction")
})

test_that("debtox without concentration errors", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  expect_error(bdeb_model(dat, type = "debtox"), "concentration")
})

test_that("invalid model type errors with available options", {
  df <- data.frame(id = 1, time = 0:3, length = 1:4 * 0.1)
  dat <- bdeb_data(growth = df)
  expect_error(bdeb_model(dat, type = "banana"), "should be one of")
})


# --- bdeb_model: temperature contract ---

test_that("temperature = NULL produces has_temperature = 0 in stan_data", {
  df <- data.frame(id = 1, time = 0:3, length = 1:4 * 0.1)
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, "individual", temperature = NULL)
  expect_equal(mod$stan_data$has_temperature, 0L)
})

test_that("temperature with all three fields produces has_temperature = 1", {
  df <- data.frame(id = 1, time = 0:3, length = 1:4 * 0.1)
  dat <- bdeb_data(growth = df)
  temp <- list(T_obs = 298.15, T_ref = 293.15, T_A = 8000)
  mod <- bdeb_model(dat, "individual", temperature = temp)
  expect_equal(mod$stan_data$has_temperature, 1L)
  expect_equal(mod$stan_data$T_obs, 298.15)
  expect_equal(mod$stan_data$T_ref, 293.15)
  expect_equal(mod$stan_data$T_A, 8000)
})

test_that("temperature missing T_ref errors with field name", {
  df <- data.frame(id = 1, time = 0:3, length = 1:4 * 0.1)
  dat <- bdeb_data(growth = df)
  expect_error(
    bdeb_model(dat, "individual", temperature = list(T_obs = 298, T_A = 8000)),
    "T_ref"
  )
})

test_that("temperature missing T_A errors with field name", {
  df <- data.frame(id = 1, time = 0:3, length = 1:4 * 0.1)
  dat <- bdeb_data(growth = df)
  expect_error(
    bdeb_model(dat, "individual", temperature = list(T_obs = 298, T_ref = 293)),
    "T_A"
  )
})


# --- bdeb_data: input type contracts ---

test_that("bdeb_data growth must be a data.frame", {
  expect_error(bdeb_data(growth = list(id = 1, time = 0, length = 0.1)))
  expect_error(bdeb_data(growth = matrix(1:9, 3, 3)))
})

test_that("bdeb_data with no arguments errors", {
  expect_error(bdeb_data(), "At least one")
})


# --- Input validation contracts ---

test_that("arrhenius rejects T <= 0", {
  expect_error(arrhenius(0), "positive")
  expect_error(arrhenius(-100), "positive")
})

test_that("arrhenius rejects T_ref <= 0", {
  expect_error(arrhenius(300, T_ref = 0), "positive")
})

test_that("arrhenius rejects T_A < 0", {
  expect_error(arrhenius(300, T_A = -1), "non-negative")
})

test_that("prior_lognormal rejects sigma <= 0", {
  expect_error(prior_lognormal(sigma = 0), "positive")
  expect_error(prior_lognormal(sigma = -1), "positive")
})

test_that("prior_beta rejects a <= 0 or b <= 0", {
  expect_error(prior_beta(a = 0), "positive")
  expect_error(prior_beta(b = -1), "positive")
})

test_that("prior_exponential rejects rate <= 0", {
  expect_error(prior_exponential(rate = 0), "positive")
})

test_that("prior_halfnormal rejects sigma <= 0", {
  expect_error(prior_halfnormal(sigma = -0.1), "positive")
})

test_that("obs_student_t rejects nu < 1", {
  expect_no_error(obs_student_t(nu = 1))
  expect_error(obs_student_t(nu = 0), "nu.*>= 1")
  expect_error(obs_student_t(nu = -1), "nu.*>= 1")
})

test_that("bdeb_data rejects f_food outside [0,1]", {
  df <- data.frame(id = 1, time = 0:3, length = 1:4 * 0.1)
  expect_error(bdeb_data(growth = df, f_food = -0.1), "\\[0, 1\\]")
  expect_error(bdeb_data(growth = df, f_food = 1.5), "\\[0, 1\\]")
})

test_that("bdeb_fit rejects invalid sampling parameters", {
  df <- data.frame(id = 1, time = 0:3, length = 1:4 * 0.1)
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, "individual")

  expect_error(bdeb_fit(mod, chains = 0), "chains")
  expect_error(bdeb_fit(mod, iter_sampling = 0), "iter_sampling")
  expect_error(bdeb_fit(mod, adapt_delta = 0), "adapt_delta")
  expect_error(bdeb_fit(mod, adapt_delta = 1), "adapt_delta")
  expect_error(bdeb_fit(mod, max_treedepth = 0), "max_treedepth")
  expect_error(bdeb_fit(mod, threads_per_chain = 0), "threads_per_chain")
})


# --- Prior override contract ---

test_that("custom priors override exactly those specified", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)

  custom <- list(p_Am = prior_lognormal(mu = 99, sigma = 0.01))
  mod <- bdeb_model(dat, "individual", priors = custom)

  # Overridden
  expect_equal(mod$priors$p_Am$mu, 99)
  expect_equal(mod$priors$p_Am$sigma, 0.01)

  # All others from default
  defaults <- prior_default("individual")
  for (nm in setdiff(names(defaults), "p_Am")) {
    expect_equal(mod$priors[[nm]]$family, defaults[[nm]]$family,
                 info = paste("prior", nm, "should be default"))
  }
})


# --- bdeb_fit: input contracts ---

test_that("bdeb_fit rejects non-bdeb_model", {
  expect_error(bdeb_fit(list()), "bdeb_model")
  expect_error(bdeb_fit(data.frame()), "bdeb_model")
  expect_error(bdeb_fit(42), "bdeb_model")
})


# --- bdeb_ppc: model type contracts ---

test_that("bdeb_ppc on individual model returns growth PPC", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "individual")
  ppc <- bdeb_ppc(fit, type = "growth")
  expect_s3_class(ppc, "bdeb_ppc")
  expect_false(is.null(ppc$growth))
})

test_that("bdeb_ppc on hierarchical errors with informative message", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "hierarchical")
  expect_error(bdeb_ppc(fit), "not available.*hierarchical")
})

test_that("bdeb_ppc on debtox errors with informative message", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "debtox")
  expect_error(bdeb_ppc(fit), "not available.*debtox")
})


# --- bdeb_ec50: model type contract ---

test_that("bdeb_ec50 on non-debtox model errors", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "individual")
  expect_error(bdeb_ec50(fit), "DEBtox")
})


# --- bdeb_diagnose, summary, bdeb_derived: class contracts ---

test_that("bdeb_diagnose rejects wrong class", {
  expect_error(bdeb_diagnose(list()), "bdeb_fit")
})

test_that("summary.bdeb_fit rejects wrong class on direct dispatch", {
  expect_error(summary.bdeb_fit(list()), "bdeb_fit")
})

test_that("bdeb_derived rejects wrong class", {
  expect_error(bdeb_derived(list()), "bdeb_fit")
})


# --- Observation model defaults ---

test_that("default observation model for growth is normal", {
  df <- data.frame(id = 1, time = 0:3, length = 1:4 * 0.1)
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, "individual")
  expect_equal(mod$observation$growth$family, "normal")
})

test_that("default observation model for reproduction is negbinom", {
  df <- data.frame(id = 1, time = 0:3, length = 1:4 * 0.1)
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, "individual")
  expect_equal(mod$observation$reproduction$family, "negbinom")
})
