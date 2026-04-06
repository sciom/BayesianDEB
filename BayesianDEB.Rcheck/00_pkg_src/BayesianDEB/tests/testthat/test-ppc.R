# ===========================================================
# Tests: ppc.R  (bdeb_ppc, print.bdeb_ppc, predict.bdeb_fit, bdeb_predict)
# ===========================================================

# --- bdeb_ppc ---

test_that("bdeb_ppc rejects non-bdeb_fit", {
  expect_error(bdeb_ppc(list()), "bdeb_fit")
  expect_error(bdeb_ppc(NULL), "bdeb_fit")
  expect_error(bdeb_ppc("fit"), "bdeb_fit")
})

test_that("bdeb_ppc rejects hierarchical model with clear message", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "hierarchical")
  expect_error(bdeb_ppc(fit), "not available.*hierarchical")
})

test_that("bdeb_ppc rejects debtox model with clear message", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "debtox")
  expect_error(bdeb_ppc(fit), "not available.*debtox")
})

test_that("bdeb_ppc extracts growth PPC from mock", {
  fit <- mock_bdeb_fit(n_draws = 80)
  ppc <- bdeb_ppc(fit, type = "growth")

  expect_s3_class(ppc, "bdeb_ppc")
  expect_equal(ppc$type, "growth")
  expect_false(is.null(ppc$growth))
  expect_true(is.matrix(ppc$growth$L_rep))
  expect_equal(ppc$growth$n_draws, 80)
  expect_equal(ppc$growth$n_obs, length(ppc$growth$L_obs))
})

test_that("bdeb_ppc L_rep dimensions match", {
  fit <- mock_bdeb_fit(n_draws = 60)
  ppc <- bdeb_ppc(fit, type = "growth")

  # L_rep rows = n_draws, columns = n_obs
  expect_equal(nrow(ppc$growth$L_rep), 60)
  expect_equal(ncol(ppc$growth$L_rep), ppc$growth$n_obs)
})

test_that("bdeb_ppc includes t_obs for individual model", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "individual")
  ppc <- bdeb_ppc(fit, type = "growth")

  expect_false(is.null(ppc$growth$t_obs))
  expect_equal(length(ppc$growth$t_obs), ppc$growth$n_obs)
})

test_that("print.bdeb_ppc works", {
  fit <- mock_bdeb_fit(n_draws = 50)
  ppc <- bdeb_ppc(fit, type = "growth")
  expect_invisible(print(ppc))
})


# --- bdeb_prior_predictive ---

test_that("bdeb_prior_predictive rejects non-bdeb_model", {
  expect_error(bdeb_prior_predictive(list()), "bdeb_model")
})

test_that("bdeb_prior_predictive returns correct structure", {
  df <- data.frame(id = 1, time = seq(0, 42, by = 7),
                   length = seq(0.1, 0.4, length.out = 7))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, "individual")

  pp <- bdeb_prior_predictive(mod, n_draws = 50, seed = 42)

  expect_s3_class(pp, "bdeb_prior_predictive")
  expect_equal(nrow(pp$L), 50)
  expect_equal(length(pp$L_inf), 50)
  expect_true(all(pp$L_inf > 0))
  expect_true(all(pp$t >= 0))
  expect_lte(max(pp$t), 42)
  expect_gte(max(pp$t), 35)  # close to t_max
})

test_that("bdeb_prior_predictive is reproducible with seed", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, "individual")

  pp1 <- bdeb_prior_predictive(mod, n_draws = 20, seed = 99)
  pp2 <- bdeb_prior_predictive(mod, n_draws = 20, seed = 99)
  expect_equal(pp1$L_inf, pp2$L_inf)
  expect_equal(pp1$L, pp2$L)
})

test_that("print.bdeb_prior_predictive works", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, "individual")
  pp <- bdeb_prior_predictive(mod, n_draws = 20, seed = 1)
  expect_invisible(print(pp))
})

test_that("plot.bdeb_prior_predictive returns ggplot", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, "individual")
  pp <- bdeb_prior_predictive(mod, n_draws = 20, seed = 1)
  p <- plot(pp, n_draws = 10)
  expect_s3_class(p, "ggplot")
})


# --- bdeb_predict / predict.bdeb_fit ---

test_that("bdeb_predict rejects non-bdeb_fit", {
  expect_error(bdeb_predict(list()), "bdeb_fit")
  expect_error(bdeb_predict(NULL), "bdeb_fit")
})

test_that("bdeb_predict returns prediction object", {
  fit <- mock_bdeb_fit(n_draws = 80)
  pred <- bdeb_predict(fit, n_draws = 40)

  expect_s3_class(pred, "bdeb_prediction")
  expect_equal(pred$n_draws, 40)
  expect_true(is.matrix(pred$L_hat))
  expect_equal(nrow(pred$L_hat), 40)
  expect_equal(pred$model_type, "individual")
})

test_that("predict.bdeb_fit dispatches correctly", {
  fit <- mock_bdeb_fit(n_draws = 80)
  pred <- predict(fit, n_draws = 30)

  expect_s3_class(pred, "bdeb_prediction")
  expect_equal(pred$n_draws, 30)
})

test_that("bdeb_predict respects n_draws limit", {
  fit <- mock_bdeb_fit(n_draws = 50)

  # Request more than available
  pred <- bdeb_predict(fit, n_draws = 200)
  expect_lte(pred$n_draws, 50)

  # Request fewer
  pred2 <- bdeb_predict(fit, n_draws = 10)
  expect_equal(pred2$n_draws, 10)
})

test_that("bdeb_predict with newdata requires t_predict", {
  fit <- mock_bdeb_fit(n_draws = 50)
  expect_error(bdeb_predict(fit, newdata = list(f_food = 0.5)), "t_predict")
})

test_that("bdeb_predict with newdata forward-simulates", {
  fit <- mock_bdeb_fit(n_draws = 50)
  pred <- bdeb_predict(fit,
    newdata = list(t_predict = seq(0, 100, by = 10), f_food = 0.7),
    n_draws = 10)

  expect_s3_class(pred, "bdeb_prediction")
  expect_equal(length(pred$t), 11)
  expect_equal(nrow(pred$L_hat), 10)
  expect_equal(ncol(pred$L_hat), 11)
  # All predicted lengths should be positive
  expect_true(all(pred$L_hat > 0, na.rm = TRUE))
})

test_that("bdeb_predict newdata uses different f_food", {
  fit <- mock_bdeb_fit(n_draws = 50)
  t_pred <- seq(0, 60, by = 10)

  pred_f1  <- bdeb_predict(fit, newdata = list(t_predict = t_pred, f_food = 1.0),
                           n_draws = 10)
  pred_f05 <- bdeb_predict(fit, newdata = list(t_predict = t_pred, f_food = 0.5),
                           n_draws = 10)

  # Lower food -> smaller final length
  mean_f1  <- mean(pred_f1$L_hat[, ncol(pred_f1$L_hat)])
  mean_f05 <- mean(pred_f05$L_hat[, ncol(pred_f05$L_hat)])
  expect_gt(mean_f1, mean_f05)
})

test_that("bdeb_predict newdata with debtox and concentration", {
  fit <- mock_bdeb_fit(n_draws = 30, type = "debtox")
  pred <- bdeb_predict(fit,
    newdata = list(t_predict = seq(0, 42, by = 7), concentration = 100),
    n_draws = 5)

  expect_s3_class(pred, "bdeb_prediction")
  expect_equal(pred$model_type, "debtox")
  expect_true(all(pred$L_hat > 0, na.rm = TRUE))
})

test_that("plot.bdeb_prediction returns ggplot", {
  fit <- mock_bdeb_fit(n_draws = 50)
  pred <- bdeb_predict(fit,
    newdata = list(t_predict = seq(0, 50, by = 10)),
    n_draws = 10)
  p <- plot(pred, n_draws = 5)
  expect_s3_class(p, "ggplot")
})

test_that("bdeb_predict includes time vector", {
  fit <- mock_bdeb_fit(n_draws = 50)
  pred <- bdeb_predict(fit, n_draws = 20)

  expect_false(is.null(pred$t))
  expect_true(length(pred$t) > 0)
})
