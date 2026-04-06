# ===========================================================
# Tests: debtox.R  (bdeb_tox, bdeb_ec50, plot_dose_response)
# ===========================================================

# --- bdeb_tox ---

test_that("bdeb_tox rejects missing concentration", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  expect_error(bdeb_tox(dat), "concentration")
})

test_that("bdeb_tox creates debtox model", {
  df <- data.frame(
    id = rep(1:2, each = 5), time = rep(0:4, 2),
    length = runif(10, 0.1, 0.5)
  )
  conc <- c("1" = 0, "2" = 50)
  dat <- bdeb_data(growth = df, concentration = conc)

  mod <- bdeb_tox(dat, stress = "assimilation")
  expect_s3_class(mod, "bdeb_model")
  expect_equal(mod$type, "debtox")
  expect_equal(mod$stan_model_name, "bdeb_debtox")
})

test_that("bdeb_tox falls back to assimilation for unsupported stress modes", {
  df <- data.frame(
    id = rep(1:2, each = 5), time = rep(0:4, 2),
    length = runif(10, 0.1, 0.5)
  )
  conc <- c("1" = 0, "2" = 50)
  dat <- bdeb_data(growth = df, concentration = conc)

  # cli_alert_warning is not a base R warning, so just verify the model
  # is created successfully (stress silently falls back to assimilation)
  mod <- suppressMessages(bdeb_tox(dat, stress = "maintenance"))
  expect_s3_class(mod, "bdeb_model")
  expect_equal(mod$type, "debtox")
})

test_that("bdeb_tox accepts custom priors", {
  df <- data.frame(
    id = rep(1:2, each = 5), time = rep(0:4, 2),
    length = runif(10, 0.1, 0.5)
  )
  conc <- c("1" = 0, "2" = 50)
  dat <- bdeb_data(growth = df, concentration = conc)

  mod <- bdeb_tox(dat,
    priors = list(z_w = prior_lognormal(mu = 3.0, sigma = 0.5)))
  expect_equal(mod$priors$z_w$mu, 3.0)
})


# --- bdeb_ec50 ---

test_that("bdeb_ec50 rejects non-bdeb_fit", {
  expect_error(bdeb_ec50(list()), "bdeb_fit")
  expect_error(bdeb_ec50(NULL), "bdeb_fit")
})

test_that("bdeb_ec50 rejects non-debtox fit", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "individual")
  expect_error(bdeb_ec50(fit), "DEBtox")
})

test_that("bdeb_ec50 works with debtox mock", {
  fit <- mock_bdeb_fit(n_draws = 100, type = "debtox")
  result <- bdeb_ec50(fit, prob = 0.90)

  expect_type(result, "list")
  expect_true("draws" %in% names(result))
  expect_true("summary" %in% names(result))
  expect_true("NEC" %in% names(result))

  # Summary should have EC50 and NEC rows
  expect_equal(nrow(result$summary), 2)
  expect_true("EC50" %in% result$summary$parameter)
  expect_true("NEC" %in% result$summary$parameter)
})

test_that("bdeb_ec50 EC50 > NEC", {
  fit <- mock_bdeb_fit(n_draws = 100, type = "debtox")
  result <- bdeb_ec50(fit, prob = 0.90)

  ec50_mean <- result$summary$mean[result$summary$parameter == "EC50"]
  nec_mean  <- result$summary$mean[result$summary$parameter == "NEC"]
  expect_gt(ec50_mean, nec_mean)
})

test_that("bdeb_ec50 returns correct number of draws", {
  fit <- mock_bdeb_fit(n_draws = 80, type = "debtox")
  result <- bdeb_ec50(fit)
  expect_length(result$draws, 80)
})

test_that("bdeb_ec50 respects prob argument", {
  fit <- mock_bdeb_fit(n_draws = 200, type = "debtox")

  r90 <- bdeb_ec50(fit, prob = 0.90)
  r50 <- bdeb_ec50(fit, prob = 0.50)

  # 90% interval wider than 50%
  width90 <- r90$summary$upper[1] - r90$summary$lower[1]
  width50 <- r50$summary$upper[1] - r50$summary$lower[1]
  expect_gt(width90, width50)
})

test_that("bdeb_ec50 summary has correct columns", {
  fit <- mock_bdeb_fit(n_draws = 100, type = "debtox")
  result <- bdeb_ec50(fit)

  expect_true(all(c("parameter", "mean", "median", "sd", "lower", "upper")
                  %in% names(result$summary)))
})

test_that("bdeb_ec50 all draws are positive", {
  fit <- mock_bdeb_fit(n_draws = 100, type = "debtox")
  result <- bdeb_ec50(fit)
  expect_true(all(result$draws > 0))
})


# --- plot_dose_response ---

test_that("plot_dose_response rejects non-debtox fit", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "individual")
  expect_error(plot_dose_response(fit), "DEBtox")
})

test_that("plot_dose_response rejects non-bdeb_fit", {
  expect_error(plot_dose_response(list()), "DEBtox|bdeb_fit")
})

test_that("plot_dose_response returns ggplot for debtox mock", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "debtox")
  p <- plot_dose_response(fit, n_draws = 20)
  expect_s3_class(p, "ggplot")
})
