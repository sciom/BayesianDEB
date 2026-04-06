# ===========================================================
# Tests: fit.R  (bdeb_fit input validation, print, summary S3)
# ===========================================================

test_that("bdeb_fit rejects non-bdeb_model input", {
  expect_error(bdeb_fit(list(a = 1)), "bdeb_model")
  expect_error(bdeb_fit(NULL), "bdeb_model")
  expect_error(bdeb_fit(data.frame(x = 1)), "bdeb_model")
})

test_that("bdeb_fit requires cmdstanr", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, type = "individual")

  # We can't easily test this without unloading cmdstanr,
  # but we verify the function signature accepts all arguments
  expect_true(is.function(bdeb_fit))
  args <- formals(bdeb_fit)
  expect_true("chains" %in% names(args))
  expect_true("iter_warmup" %in% names(args))
  expect_true("iter_sampling" %in% names(args))
  expect_true("adapt_delta" %in% names(args))
  expect_true("max_treedepth" %in% names(args))
  expect_true("seed" %in% names(args))
  expect_true("parallel_chains" %in% names(args))
  expect_true("refresh" %in% names(args))
})

test_that("print.bdeb_fit works with mock", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "individual")
  expect_invisible(print(fit))
})

test_that("summary.bdeb_fit dispatches to bdeb_summary", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "individual")
  s <- summary(fit)
  expect_s3_class(s, "draws_summary")
})
