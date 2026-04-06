# ===========================================================
# Tests: ppc.R  (bdeb_ppc, print.bdeb_ppc, predict.bdeb_fit, bdeb_predict)
# ===========================================================

# --- bdeb_ppc ---

test_that("bdeb_ppc rejects non-bdeb_fit", {
  expect_error(bdeb_ppc(list()), "bdeb_fit")
  expect_error(bdeb_ppc(NULL), "bdeb_fit")
  expect_error(bdeb_ppc("fit"), "bdeb_fit")
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

test_that("bdeb_predict includes time vector", {
  fit <- mock_bdeb_fit(n_draws = 50)
  pred <- bdeb_predict(fit, n_draws = 20)

  expect_false(is.null(pred$t))
  expect_true(length(pred$t) > 0)
})
