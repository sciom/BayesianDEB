# ===========================================================
# Tests: model_spec.R  (bdeb_model, observation models)
# ===========================================================

# --- bdeb_model: happy path ---

test_that("bdeb_model creates individual model", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, type = "individual")

  expect_s3_class(mod, "bdeb_model")
  expect_equal(mod$type, "individual")
  expect_equal(mod$stan_model_name, "bdeb_individual_growth")
  expect_type(mod$stan_data, "list")
  expect_true("N_obs" %in% names(mod$stan_data))
})

test_that("bdeb_model fills default priors", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, type = "individual", priors = list())

  # Should have all default priors
  expect_s3_class(mod$priors$p_Am, "bdeb_prior")
  expect_s3_class(mod$priors$kappa, "bdeb_prior")
  expect_s3_class(mod$priors$sigma_L, "bdeb_prior")
})

test_that("bdeb_model rejects mismatched type/data", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)

  expect_error(bdeb_model(dat, type = "growth_repro"), "reproduction")
  expect_error(bdeb_model(dat, type = "debtox"), "concentration")
})

test_that("bdeb_model hierarchical builds correct Stan data", {
  df <- data.frame(
    id = rep(1:3, each = 6),
    time = rep(0:5, 3),
    length = runif(18, 0.1, 0.5)
  )
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, type = "hierarchical")

  expect_equal(mod$stan_data$N_ind, 3)
  expect_equal(mod$stan_data$max_N_obs, 6)
  expect_true(is.matrix(mod$stan_data$t_obs))
})

test_that("print.bdeb_model works", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, type = "individual")
  expect_invisible(print(mod))
})

test_that("observation model constructors work", {
  expect_s3_class(obs_normal(), "bdeb_obs")
  expect_s3_class(obs_lognormal(), "bdeb_obs")
  expect_s3_class(obs_student_t(nu = 3), "bdeb_obs")
  expect_s3_class(obs_poisson(), "bdeb_obs")
  expect_s3_class(obs_negbinom(), "bdeb_obs")

  expect_equal(obs_student_t(nu = 7)$nu, 7)
})


# --- bdeb_model: edge cases ---

test_that("bdeb_model rejects non-bdeb_data input", {
  expect_error(bdeb_model(data.frame(x = 1)), "bdeb_data")
  expect_error(bdeb_model(list(a = 1)), "bdeb_data")
  expect_error(bdeb_model(NULL), "bdeb_data")
})

test_that("bdeb_model rejects invalid type", {
  df <- data.frame(id = 1, time = 0:3, length = seq(0.1, 0.4, by = 0.1))
  dat <- bdeb_data(growth = df)
  expect_error(bdeb_model(dat, type = "invalid_type"), "should be one of")
})

test_that("bdeb_model errors for individual with multiple IDs", {
  df <- data.frame(
    id = rep(1:3, each = 5),
    time = rep(0:4, 3),
    length = runif(15, 0.1, 0.5)
  )
  dat <- bdeb_data(growth = df)
  expect_error(
    bdeb_model(dat, type = "individual"),
    "requires exactly 1 individual"
  )
})

test_that("bdeb_model growth_repro rejects multiple individuals", {
  g <- data.frame(id = rep(1:3, each = 5), time = rep(0:4, 3),
                  length = runif(15, 0.1, 0.5))
  r <- data.frame(id = 1:3, t_start = 0, t_end = 28, count = c(10, 20, 30))
  dat <- bdeb_data(growth = g, reproduction = r)
  expect_error(
    bdeb_model(dat, type = "growth_repro"),
    "single-individual model"
  )
})

test_that("bdeb_model custom priors override defaults", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)

  custom_p_Am <- prior_lognormal(mu = 99, sigma = 0.01)
  mod <- bdeb_model(dat, type = "individual",
                    priors = list(p_Am = custom_p_Am))

  expect_equal(mod$priors$p_Am$mu, 99)
  expect_equal(mod$priors$p_Am$sigma, 0.01)
  # Other priors should still be defaults
  defaults <- prior_default("individual")
  expect_equal(mod$priors$kappa$a, defaults$kappa$a)
  expect_equal(mod$priors$p_M$mu, defaults$p_M$mu)
})

test_that("bdeb_model passes extra priors without error", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  # extra_par is not a standard DEB prior — should not break anything
  mod <- bdeb_model(dat, type = "individual",
                    priors = list(extra_par = prior_normal(0, 1)))
  expect_true("extra_par" %in% names(mod$priors))
})

test_that("bdeb_model sets default observation models", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, type = "individual")

  expect_s3_class(mod$observation$growth, "bdeb_obs")
  expect_equal(mod$observation$growth$family, "normal")
  expect_s3_class(mod$observation$reproduction, "bdeb_obs")
  expect_equal(mod$observation$reproduction$family, "negbinom")
})

test_that("bdeb_model stores temperature correction", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  temp <- list(T_obs = 298.15, T_ref = 293.15, T_A = 8000)
  mod <- bdeb_model(dat, type = "individual", temperature = temp)

  expect_equal(mod$temperature$T_obs, 298.15)
  expect_equal(mod$temperature$T_ref, 293.15)
  expect_equal(mod$temperature$T_A, 8000)
})

test_that("bdeb_model individual stan_data has correct N_obs", {
  n <- 20
  df <- data.frame(id = 1, time = seq(0, 95, length.out = n),
                   length = seq(0.1, 0.5, length.out = n))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, type = "individual")
  expect_equal(mod$stan_data$N_obs, n)
})


# --- bdeb_model: Stan model name mapping ---

test_that("bdeb_model maps type to correct Stan model name", {
  g <- data.frame(id = 1, time = 0:3, length = seq(0.1, 0.4, by = 0.1))
  r <- data.frame(id = 1, t_start = 0, t_end = 28, count = 50)
  dat_g <- bdeb_data(growth = g)
  dat_gr <- bdeb_data(growth = g, reproduction = r)

  conc <- c("1" = 0, "2" = 10)
  g2 <- data.frame(id = rep(1:2, each = 4), time = rep(0:3, 2),
                   length = runif(8, 0.1, 0.5))
  dat_tox <- bdeb_data(growth = g2, concentration = conc)

  expect_equal(bdeb_model(dat_g, "individual")$stan_model_name,
               "bdeb_individual_growth")
  expect_equal(bdeb_model(dat_gr, "growth_repro")$stan_model_name,
               "bdeb_growth_repro")
  expect_equal(bdeb_model(dat_g, "hierarchical")$stan_model_name,
               "bdeb_hierarchical_growth")
  expect_equal(bdeb_model(dat_tox, "debtox")$stan_model_name,
               "bdeb_debtox")
})


# --- Observation model constructors: edge cases ---

test_that("obs_normal has correct family", {
  o <- obs_normal()
  expect_equal(o$family, "normal")
})

test_that("obs_lognormal has correct family", {
  o <- obs_lognormal()
  expect_equal(o$family, "lognormal")
})

test_that("obs_student_t default nu is 5", {
  o <- obs_student_t()
  expect_equal(o$nu, 5)
})

test_that("obs_student_t custom nu", {
  o <- obs_student_t(nu = 3)
  expect_equal(o$nu, 3)
})

test_that("obs_poisson has correct family", {
  o <- obs_poisson()
  expect_equal(o$family, "poisson")
})

test_that("obs_negbinom has correct family", {
  o <- obs_negbinom()
  expect_equal(o$family, "negbinom")
})


# --- Growth+repro Stan data: edge cases ---

test_that("bdeb_model growth_repro builds combined time grid", {
  g <- data.frame(id = 1, time = c(0, 7, 14, 21, 28),
                  length = c(0.1, 0.2, 0.3, 0.4, 0.45))
  r <- data.frame(id = 1, t_start = c(14, 21), t_end = c(21, 28),
                  count = c(5, 10))
  dat <- bdeb_data(growth = g, reproduction = r)
  mod <- bdeb_model(dat, type = "growth_repro")

  expect_true("N_times" %in% names(mod$stan_data))
  expect_true("t_all" %in% names(mod$stan_data))
  expect_true("idx_L" %in% names(mod$stan_data))
  expect_true("idx_R_start" %in% names(mod$stan_data))
  expect_true("idx_R_end" %in% names(mod$stan_data))
  # All observation indices should be within t_all range
  expect_true(all(mod$stan_data$idx_L <= mod$stan_data$N_times))
  expect_true(all(mod$stan_data$idx_R_end <= mod$stan_data$N_times))
})


# --- DEBtox Stan data: edge cases ---

# --- Temperature correction ---

test_that("bdeb_model with temperature adds has_temperature=1 to stan_data", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  temp <- list(T_obs = 298.15, T_ref = 293.15, T_A = 8000)
  mod <- bdeb_model(dat, type = "individual", temperature = temp)

  expect_equal(mod$stan_data$has_temperature, 1L)
  expect_equal(mod$stan_data$T_obs, 298.15)
  expect_equal(mod$stan_data$T_ref, 293.15)
  expect_equal(mod$stan_data$T_A, 8000)
})

test_that("bdeb_model without temperature sets has_temperature=0", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, type = "individual")

  expect_equal(mod$stan_data$has_temperature, 0L)
  # T_obs defaults don't matter (c_T will be 1.0)
  expect_true("T_obs" %in% names(mod$stan_data))
})

test_that("bdeb_model rejects incomplete temperature specification", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)

  expect_error(bdeb_model(dat, "individual", temperature = list(T_obs = 298)),
               "T_ref")
  expect_error(bdeb_model(dat, "individual", temperature = list(T_obs = 298, T_ref = 293)),
               "T_A")
})

test_that("temperature_to_stan_data produces correct c_T=1 placeholder", {
  sd <- temperature_to_stan_data(NULL)
  expect_equal(sd$has_temperature, 0L)
  # With T_A=0, exp(0/T_ref - 0/T_obs) = exp(0) = 1
  expect_equal(sd$T_A, 0.0)
})

test_that("temperature works for all model types", {
  temp <- list(T_obs = 298.15, T_ref = 293.15, T_A = 8000)

  # Individual
  df <- data.frame(id = 1, time = 0:3, length = seq(0.1, 0.4, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, "individual", temperature = temp)
  expect_equal(mod$stan_data$has_temperature, 1L)

  # Hierarchical
  mod_h <- bdeb_model(dat, "hierarchical", temperature = temp)
  expect_equal(mod_h$stan_data$has_temperature, 1L)

  # Growth+repro
  r <- data.frame(id = 1, t_start = 0, t_end = 28, count = 10)
  dat_gr <- bdeb_data(growth = df, reproduction = r)
  mod_gr <- bdeb_model(dat_gr, "growth_repro", temperature = temp)
  expect_equal(mod_gr$stan_data$has_temperature, 1L)

  # Debtox
  df2 <- data.frame(id = rep(1:2, each = 4), time = rep(0:3, 2),
                    length = runif(8, 0.1, 0.5))
  conc <- c("1" = 0, "2" = 50)
  dat_t <- bdeb_data(growth = df2, concentration = conc)
  mod_t <- bdeb_model(dat_t, "debtox", temperature = temp)
  expect_equal(mod_t$stan_data$has_temperature, 1L)
})


test_that("bdeb_model debtox with named vector concentrations", {
  df <- data.frame(
    id = rep(1:4, each = 5),
    time = rep(0:4, 4),
    length = runif(20, 0.1, 0.5)
  )
  conc <- c("1" = 0, "2" = 10, "3" = 50, "4" = 200)
  dat <- bdeb_data(growth = df, concentration = conc)
  mod <- bdeb_model(dat, type = "debtox")

  expect_equal(mod$stan_data$N_groups, 4)
  expect_equal(sort(mod$stan_data$C_w), c(0, 10, 50, 200))
  expect_true(is.matrix(mod$stan_data$t_obs))
  expect_true(is.matrix(mod$stan_data$L_obs))
})
