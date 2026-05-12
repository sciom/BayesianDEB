# ===========================================================
# Tests: print/summary/plot methods on bdeb_data, bdeb_model,
# bdeb_prior (B1, B2, B3).
# ===========================================================

# --- bdeb_data: summary + plot ---

test_that("summary.bdeb_data returns the documented structure", {
  df <- data.frame(id = c(1, 1, 2, 2),
                   time = c(0, 7, 0, 14),
                   length = c(0.1, 0.2, 0.1, 0.3))
  dat <- bdeb_data(growth = df, f_food = 0.7)
  s <- summary(dat)

  expect_s3_class(s, "summary.bdeb_data")
  expect_equal(s$n_obs_growth, 4)
  expect_equal(s$n_individuals, 2)
  expect_true(s$has_growth)
  expect_false(s$has_reproduction)
  expect_equal(s$f_food, 0.7)
  expect_equal(s$time_range, c(0, 14))
})

test_that("summary.bdeb_data captures concentration levels", {
  g <- data.frame(id = 1:4, time = rep(0, 4), length = rep(0.1, 4))
  conc <- c("1" = 0, "2" = 20, "3" = 80, "4" = 200)
  dat <- bdeb_data(growth = g, concentration = conc)
  s <- summary(dat)

  expect_equal(s$conc_levels, c(0, 20, 80, 200))
})

test_that("print.summary.bdeb_data is invisible", {
  df <- data.frame(id = 1, time = 0:3, length = (1:4) * 0.1)
  dat <- bdeb_data(growth = df)
  s <- summary(dat)
  expect_invisible(print(s))
})

test_that("plot.bdeb_data returns a ggplot for growth and reproduction", {
  g <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = g)
  p_g <- plot(dat)
  expect_s3_class(p_g, "ggplot")

  r <- data.frame(id = 1, t_start = c(0, 7), t_end = c(7, 14),
                  count = c(10, 20))
  dat_r <- bdeb_data(reproduction = r)
  p_r <- plot(dat_r, endpoint = "reproduction")
  expect_s3_class(p_r, "ggplot")
})

test_that("plot.bdeb_data errors when endpoint is missing", {
  r <- data.frame(id = 1, t_start = c(0, 7), t_end = c(7, 14),
                  count = c(10, 20))
  dat <- bdeb_data(reproduction = r)
  expect_error(plot(dat, endpoint = "growth"), "growth")
})


# --- bdeb_model: summary + plot ---

test_that("summary.bdeb_model exposes structure and priors", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, type = "individual")

  s <- summary(mod)
  expect_s3_class(s, "summary.bdeb_model")
  expect_equal(s$type, "individual")
  expect_equal(s$stan_model_name, "bdeb_individual_growth")
  expect_true("p_Am" %in% s$par_names)
  expect_invisible(print(s))
})

test_that("plot.bdeb_model returns a ggplot", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, type = "individual")

  p <- plot(mod, pars = c("p_Am", "kappa"))
  expect_s3_class(p, "ggplot")
})

test_that("plot.bdeb_model rejects unknown parameters", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, type = "individual")
  expect_error(plot(mod, pars = "no_such_param"), "Unknown")
})


# --- bdeb_prior: print / summary / plot ---

test_that("print.bdeb_prior emits and returns invisibly", {
  p <- prior_lognormal(mu = 1.0, sigma = 0.5)
  expect_invisible(print(p))
})

test_that("summary.bdeb_prior reports theoretical statistics", {
  ln <- prior_lognormal(mu = 0, sigma = 1)
  s_ln <- summary(ln)
  expect_s3_class(s_ln, "summary.bdeb_prior")
  expect_equal(s_ln$family, "lognormal")
  expect_equal(s_ln$mean, exp(0.5), tolerance = 1e-10)

  bt <- prior_beta(a = 2, b = 5)
  s_bt <- summary(bt)
  expect_equal(s_bt$mean, 2 / 7, tolerance = 1e-10)

  hn <- prior_halfnormal(sigma = 1)
  s_hn <- summary(hn)
  expect_equal(s_hn$mean, sqrt(2 / pi), tolerance = 1e-10)

  ex <- prior_exponential(rate = 2)
  s_ex <- summary(ex)
  expect_equal(s_ex$mean, 0.5, tolerance = 1e-10)
})

test_that("summary.bdeb_prior leaves halfcauchy mean as NA", {
  hc <- prior_halfcauchy(sigma = 1)
  s_hc <- summary(hc)
  expect_true(is.na(s_hc$mean))
  expect_true(is.finite(s_hc$q50))
})

test_that("plot.bdeb_prior returns a ggplot for every supported family", {
  fams <- list(
    prior_lognormal(0, 1),
    prior_normal(0, 1),
    prior_beta(2, 2),
    prior_halfnormal(1),
    prior_halfcauchy(1),
    prior_exponential(1)
  )
  for (p in fams) {
    expect_s3_class(plot(p), "ggplot")
  }
})

test_that("print.summary.bdeb_prior is invisible", {
  p <- prior_lognormal(1, 0.5)
  s <- summary(p)
  expect_invisible(print(s))
})
