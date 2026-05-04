# ===========================================================
# Tests: input validation across the public API (B10 audit).
# ===========================================================

# --- bdeb_data ---

test_that("bdeb_data rejects non-data-frame growth/reproduction", {
  expect_error(bdeb_data(growth = "oops"), "data frame")
  expect_error(bdeb_data(reproduction = matrix(1, 2, 2)), "data frame")
})

test_that("bdeb_data rejects non-numeric / out-of-range f_food", {
  df <- data.frame(id = 1, time = 0:2, length = c(0.1, 0.2, 0.3))
  expect_error(bdeb_data(growth = df, f_food = "1"), "f_food")
  expect_error(bdeb_data(growth = df, f_food = 1.5), "f_food")
  expect_error(bdeb_data(growth = df, f_food = -0.1), "f_food")
  expect_error(bdeb_data(growth = df, f_food = NA_real_), "f_food")
})


# --- bdeb_fit ---

test_that("bdeb_fit rejects bad numeric arguments", {
  fake_model <- structure(list(type = "individual"), class = "bdeb_model")
  expect_error(bdeb_fit(fake_model, chains = 0), "chains")
  expect_error(bdeb_fit(fake_model, iter_warmup = -1), "iter_warmup")
  expect_error(bdeb_fit(fake_model, iter_sampling = 0), "iter_sampling")
  expect_error(bdeb_fit(fake_model, adapt_delta = 1.0), "adapt_delta")
  expect_error(bdeb_fit(fake_model, adapt_delta = 0.0), "adapt_delta")
  expect_error(bdeb_fit(fake_model, max_treedepth = 0), "max_treedepth")
  expect_error(bdeb_fit(fake_model, threads_per_chain = 0),
               "threads_per_chain")
})


# --- summary.bdeb_fit ---

test_that("summary.bdeb_fit validates prob and pars", {
  fit <- mock_bdeb_fit(n_draws = 40)
  expect_error(summary(fit, prob = 0), "prob")
  expect_error(summary(fit, prob = 1), "prob")
  expect_error(summary(fit, prob = -0.1), "prob")
  expect_error(summary(fit, prob = NA_real_), "prob")
  expect_error(summary(fit, pars = 1:3), "character")
})


# --- confint.bdeb_fit ---

test_that("confint.bdeb_fit validates level", {
  fit <- mock_bdeb_fit(n_draws = 40)
  expect_error(confint(fit, level = 0), "level")
  expect_error(confint(fit, level = 1), "level")
  expect_error(confint(fit, level = c(0.5, 0.9)), "level")
})


# --- bdeb_diagnose ---

test_that("bdeb_diagnose validates fit class and pars", {
  expect_error(bdeb_diagnose("not a fit"), "bdeb_fit")
  fit <- mock_bdeb_fit(n_draws = 40)
  expect_error(bdeb_diagnose(fit, pars = 1:3), "character")
})


# --- bdeb_predict / predict.bdeb_fit ---

test_that("bdeb_predict validates n_draws / dt / seed", {
  fit <- mock_bdeb_fit(n_draws = 40)
  expect_error(bdeb_predict(fit, n_draws = 0), "n_draws")
  expect_error(bdeb_predict(fit, n_draws = -1), "n_draws")
  expect_error(bdeb_predict(fit, n_draws = NA_real_), "n_draws")
  expect_error(bdeb_predict(fit, dt = 0), "dt")
  expect_error(bdeb_predict(fit, dt = -1), "dt")
  expect_error(bdeb_predict(fit, seed = "abc"), "seed")
})


# --- bdeb_prior_predictive ---

test_that("bdeb_prior_predictive validates n_draws / dt", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, type = "individual")
  expect_error(bdeb_prior_predictive(mod, n_draws = 0), "n_draws")
  expect_error(bdeb_prior_predictive(mod, dt = 0), "dt")
  expect_error(bdeb_prior_predictive(mod, seed = "abc"), "seed")
})


# --- summary.bdeb_prediction ---

test_that("summary.bdeb_prediction validates prob", {
  obj <- structure(list(
    L_hat = matrix(rnorm(20), nrow = 5),
    t = 1:4, n_draws = 5, model_type = "individual"
  ), class = "bdeb_prediction")
  expect_error(summary(obj, prob = 0), "prob")
  expect_error(summary(obj, prob = 1.5), "prob")
})


# --- bdeb_ec50 ---

test_that("bdeb_ec50 validates fit class, type, prob, verbose", {
  expect_error(bdeb_ec50("oops"), "bdeb_fit")
  fit_indiv <- mock_bdeb_fit(n_draws = 40, type = "individual")
  expect_error(bdeb_ec50(fit_indiv), "DEBtox")
  fit_tox <- mock_bdeb_fit(n_draws = 40, type = "debtox")
  expect_error(bdeb_ec50(fit_tox, prob = 0), "prob")
  expect_error(bdeb_ec50(fit_tox, prob = 1.5), "prob")
  expect_error(bdeb_ec50(fit_tox, verbose = "yes"), "verbose")
  expect_error(bdeb_ec50(fit_tox, verbose = NA), "verbose")
})


# --- plot_dose_response ---

test_that("plot_dose_response validates numeric args", {
  fit_tox <- mock_bdeb_fit(n_draws = 40, type = "debtox")
  expect_error(plot_dose_response(fit_tox, n_draws = 0), "n_draws")
  expect_error(plot_dose_response(fit_tox, n_conc = 1), "n_conc")
  expect_error(plot_dose_response(fit_tox, dt = -1), "dt")
  expect_error(plot_dose_response(fit_tox, t_end = 0), "t_end")
})


# --- deb_fluxes ---

test_that("deb_fluxes validates argument types and ranges", {
  base <- list(E = 1, V = 1, f = 1, p_Am = 5, p_M = 0.5,
               kappa = 0.75, v = 0.2, E_G = 4400)
  expect_error(do.call(deb_fluxes,
                       modifyList(base, list(E = "x"))),
               "scalar")
  expect_error(do.call(deb_fluxes,
                       modifyList(base, list(E = -1))),
               "non-negative")
  expect_error(do.call(deb_fluxes,
                       modifyList(base, list(f = 1.5))),
               "f")
  expect_error(do.call(deb_fluxes,
                       modifyList(base, list(kappa = 1.5))),
               "kappa")
  expect_error(do.call(deb_fluxes,
                       modifyList(base, list(kappa = 0))),
               "kappa")
})


# --- arrhenius ---

test_that("arrhenius validates temperature inputs", {
  expect_error(arrhenius("hot"), "temp")
  expect_error(arrhenius(-10), "temp")
  expect_error(arrhenius(298.15, T_ref = 0), "T_ref")
  expect_error(arrhenius(298.15, T_A = -100), "T_A")
})


# --- prior constructors ---

test_that("prior constructors validate parameters", {
  expect_error(prior_lognormal(sigma = 0), "sigma")
  expect_error(prior_lognormal(sigma = -1), "sigma")
  expect_error(prior_normal(sigma = 0), "sigma")
  expect_error(prior_beta(a = 0, b = 1), "a")
  expect_error(prior_beta(a = 1, b = -1), "b")
  expect_error(prior_halfnormal(sigma = 0), "sigma")
  expect_error(prior_halfcauchy(sigma = 0), "sigma")
  expect_error(prior_exponential(rate = 0), "rate")
})


# --- bdeb_loo ---

test_that("bdeb_loo errors clearly on unsupported model types", {
  fit <- mock_bdeb_fit(n_draws = 40, type = "hierarchical")
  expect_error(bdeb_loo(fit), "hierarchical")
})


# --- bdeb_ppc ---

test_that("bdeb_ppc rejects unsupported model types", {
  fit <- mock_bdeb_fit(n_draws = 40, type = "hierarchical")
  expect_error(bdeb_ppc(fit), "not available")
})
